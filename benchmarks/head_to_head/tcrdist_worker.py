#!/usr/bin/env python3
"""Standalone tcrdist3 worker invoked by ``benchmarks/head_to_head/run.py``.

Runs in an isolated Python environment (e.g. micromamba + bioconda ``parasail-python``
or ``Dockerfile.tcrdist``) so ``tcrdist3`` and ``parasail`` stay completely decoupled
from ``tcrsift``'s core dependencies.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import re
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse as sp
from tcrdist.paths import path_to_db
from tcrdist.repertoire import TCRrep

_DB_DF = pd.read_csv(Path(path_to_db) / "alphabeta_gammadelta_db.tsv", sep="\t")
_VALID_HUMAN_GENES = set(_DB_DF[_DB_DF["organism"] == "human"]["id"].astype(str))


def _norm_v(gene: str) -> str:
    g = str(gene or "").strip()
    if not g or g.lower() in ("nan", "none"):
        return ""
    # Normalize Cell Ranger slashless dual-locus TRAV names (e.g. TRAV14DV4 -> TRAV14/DV4)
    g = re.sub(r"(TRAV[0-9\-]+)(DV[0-9]+)", r"\1/\2", g)
    if g in _VALID_HUMAN_GENES:
        return g
    base = g.split("*")[0]
    if f"{base}*01" in _VALID_HUMAN_GENES:
        return f"{base}*01"
    for cand in sorted(_VALID_HUMAN_GENES):
        if cand.startswith(f"{base}*"):
            return cand
    return ""


def _modal_species(target_vals: np.ndarray) -> str:
    """Consensus (modal) species across nearest-distance reference hits."""
    counts: Counter[str] = Counter()
    for val in target_vals:
        for tok in str(val).split(";"):
            t = tok.strip()
            if t:
                counts[t] += 1
    if not counts:
        return ""
    # Deterministic alphabetical tie-break matching tcrsift's consensus tie-break
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def _eval_sparse(
    csr_mat: sp.csr_matrix,
    q_indices: np.ndarray,
    q_truth: np.ndarray,
    r_vals: np.ndarray,
    radii: list[int],
    mode: str,
    n_total_q: int,
    cohort: str,
    prefix: str,
    elapsed: float,
) -> list[dict]:
    out = []
    for rad in radii:
        n_m = 0
        n_c = 0
        for row_i in range(csr_mat.shape[0]):
            start, end = csr_mat.indptr[row_i], csr_mat.indptr[row_i + 1]
            if start == end:
                continue
            cols = csr_mat.indices[start:end]
            raw_d = csr_mat.data[start:end]
            dists = np.where(raw_d < 0, 0, raw_d)
            mask = dists <= rad
            if not np.any(mask):
                continue
            hits = cols[mask]
            hit_d = dists[mask]
            n_m += 1
            qi = int(q_indices[row_i])
            truth = str(q_truth[qi])
            min_d = hit_d.min()
            nearest_hits = hits[hit_d == min_d]
            if mode == "species":
                modal_sp = _modal_species(r_vals[nearest_hits])
                if truth.lower() in modal_sp.lower():
                    n_c += 1
            else:
                near_eps = {
                    tok.strip()
                    for v in r_vals[nearest_hits]
                    for tok in str(v).split(";")
                    if tok.strip()
                }
                if truth in near_eps:
                    n_c += 1
        out.append(
            {
                "cohort": cohort,
                "method": f"{prefix}_r{rad}",
                "n_query_clones": n_total_q,
                "n_matched": n_m,
                "n_cognate": n_c,
                "cognate_precision": round(n_c / n_m, 4) if n_m else None,
                "cognate_recall": round(n_c / n_total_q, 4) if n_total_q else None,
                "elapsed_sec": round(elapsed, 3),
            }
        )
    return out


def _chunked_rect_sparse(
    tr_q: TCRrep,
    tr_r: TCRrep,
    max_rad: int,
    paired: bool,
    chunk: int = 500,
) -> sp.csr_matrix:
    rows_l, cols_l, vals_l = [], [], []
    n_q = len(tr_q.clone_df)
    n_r = len(tr_r.clone_df)
    for s0 in range(0, n_q, chunk):
        sub_q = tr_q.clone_df.iloc[s0 : s0 + chunk].reset_index(drop=True)
        tr_q.compute_rect_distances(df=sub_q, df2=tr_r.clone_df)
        dmat = (tr_q.rw_alpha + tr_q.rw_beta) if paired else tr_q.rw_beta
        ri, ci = np.where(dmat <= max_rad)
        if len(ri):
            dv = dmat[ri, ci]
            rows_l.append(ri + s0)
            cols_l.append(ci)
            # Encode exact 0 distance as -1 so CSR matrix retains zero-distance entries
            vals_l.append(np.where(dv == 0, -1, dv))
    if rows_l:
        return sp.csr_matrix(
            (np.concatenate(vals_l), (np.concatenate(rows_l), np.concatenate(cols_l))),
            shape=(n_q, n_r),
        )
    return sp.csr_matrix((n_q, n_r))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--query-tsv", type=Path, required=True)
    ap.add_argument("--ref-tsv", type=Path, required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--mode", choices=("species", "epitope"), required=True)
    ap.add_argument("--radii-ab", default="24,48")
    ap.add_argument("--radii-b", default="0")
    ap.add_argument("--out-json", type=Path, required=True)
    args = ap.parse_args()

    warnings.filterwarnings("ignore")
    ver = importlib.metadata.version("tcrdist3")
    radii_ab = [int(x) for x in args.radii_ab.split(",") if x.strip()]
    radii_b = [int(x) for x in args.radii_b.split(",") if x.strip()]

    q_raw = pd.read_csv(args.query_tsv, sep="\t", keep_default_na=False)
    r_raw = pd.read_csv(args.ref_tsv, sep="\t", keep_default_na=False)

    q_df = pd.DataFrame(
        {
            "cdr3_a_aa": q_raw["CDR3_alpha"].astype(str),
            "v_a_gene": q_raw["v_alpha"].map(_norm_v),
            "j_a_gene": "TRAJ1*01",
            "cdr3_b_aa": q_raw["CDR3_beta"].astype(str),
            "v_b_gene": q_raw["v_beta"].map(_norm_v),
            "j_b_gene": "TRBJ1-1*01",
            "count": 1,
            "q_idx": np.arange(len(q_raw)),
        }
    )
    r_df = pd.DataFrame(
        {
            "cdr3_a_aa": r_raw["cdr3_a"].astype(str),
            "v_a_gene": r_raw["v_a"].map(_norm_v),
            "j_a_gene": "TRAJ1*01",
            "cdr3_b_aa": r_raw["cdr3_b"].astype(str),
            "v_b_gene": r_raw["v_b"].map(_norm_v),
            "j_b_gene": "TRBJ1-1*01",
            "count": 1,
            "target_val": r_raw["target_val"].astype(str),
        }
    )

    rows: list[dict] = []

    if radii_ab:
        q_paired = q_df[
            (q_df["cdr3_a_aa"].str.len() >= 6)
            & (q_df["cdr3_b_aa"].str.len() >= 6)
            & (q_df["v_a_gene"] != "")
            & (q_df["v_b_gene"] != "")
        ].reset_index(drop=True)
        r_paired_raw = r_df[
            (r_df["cdr3_a_aa"].str.len() >= 6)
            & (r_df["cdr3_b_aa"].str.len() >= 6)
            & (r_df["v_a_gene"] != "")
            & (r_df["v_b_gene"] != "")
        ]
        r_paired = r_paired_raw.groupby(
            ["cdr3_a_aa", "v_a_gene", "j_a_gene", "cdr3_b_aa", "v_b_gene", "j_b_gene"],
            as_index=False,
        ).agg(
            count=("count", "sum"),
            target_val=("target_val", lambda s: ";".join(x for v in s for x in str(v).split(";") if x)),
        )
        t0 = time.perf_counter()
        tr_q = TCRrep(
            cell_df=q_paired,
            organism="human",
            chains=["alpha", "beta"],
            db_file="alphabeta_gammadelta_db.tsv",
            deduplicate=False,
            compute_distances=False,
            cpus=1,
        )
        tr_r = TCRrep(
            cell_df=r_paired,
            organism="human",
            chains=["alpha", "beta"],
            db_file="alphabeta_gammadelta_db.tsv",
            deduplicate=False,
            compute_distances=False,
            cpus=1,
        )
        for obj in (tr_q, tr_r):
            obj.clone_df = obj.clone_df.dropna(
                subset=["cdr1_a_aa", "cdr2_a_aa", "pmhc_a_aa", "cdr1_b_aa", "cdr2_b_aa", "pmhc_b_aa"]
            ).reset_index(drop=True)
        csr_ab = _chunked_rect_sparse(tr_q, tr_r, max(radii_ab), paired=True)
        elapsed_ab = time.perf_counter() - t0
        rows.extend(
            _eval_sparse(
                csr_ab,
                tr_q.clone_df["q_idx"].to_numpy(),
                q_raw["truth"].to_numpy(),
                tr_r.clone_df["target_val"].to_numpy(),
                radii_ab,
                args.mode,
                len(q_raw),
                args.cohort,
                f"tcrdist3_v{ver}_paired",
                elapsed_ab,
            )
        )

    if radii_b:
        q_beta = q_df[
            (q_df["cdr3_b_aa"].str.len() >= 6) & (q_df["v_b_gene"] != "")
        ].reset_index(drop=True)
        r_beta_raw = r_df[
            (r_df["cdr3_b_aa"].str.len() >= 6) & (r_df["v_b_gene"] != "")
        ]
        r_beta = r_beta_raw.groupby(
            ["cdr3_b_aa", "v_b_gene", "j_b_gene"],
            as_index=False,
        ).agg(
            count=("count", "sum"),
            target_val=("target_val", lambda s: ";".join(x for v in s for x in str(v).split(";") if x)),
        )
        t0 = time.perf_counter()
        tr_qb = TCRrep(
            cell_df=q_beta,
            organism="human",
            chains=["beta"],
            db_file="alphabeta_gammadelta_db.tsv",
            deduplicate=False,
            compute_distances=False,
            cpus=1,
        )
        tr_rb = TCRrep(
            cell_df=r_beta,
            organism="human",
            chains=["beta"],
            db_file="alphabeta_gammadelta_db.tsv",
            deduplicate=False,
            compute_distances=False,
            cpus=1,
        )
        for obj in (tr_qb, tr_rb):
            obj.clone_df = obj.clone_df.dropna(
                subset=["cdr1_b_aa", "cdr2_b_aa", "pmhc_b_aa"]
            ).reset_index(drop=True)
        csr_b = _chunked_rect_sparse(tr_qb, tr_rb, max(radii_b), paired=False)
        elapsed_b = time.perf_counter() - t0
        rows.extend(
            _eval_sparse(
                csr_b,
                tr_qb.clone_df["q_idx"].to_numpy(),
                q_raw["truth"].to_numpy(),
                tr_rb.clone_df["target_val"].to_numpy(),
                radii_b,
                args.mode,
                len(q_raw),
                args.cohort,
                f"tcrdist3_v{ver}_trb",
                elapsed_b,
            )
        )

    args.out_json.write_text(json.dumps({"tcrdist3_version": ver, "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
