#!/usr/bin/env python3
"""Run the public paired blood-tumor (`blood-tumor`), post-ACT (`T30`), and CD4+ TIL (`--tcell-type cd4`) benchmark.

    python benchmarks/public_blood_cd4/download.py --data-dir CACHE
    python benchmarks/public_blood_cd4/run.py --data-dir CACHE --work-dir CACHE/runs

Exercises ``tcrsift`` across four open-access human single-cell cohorts:

1. Pulliam et al., Cell Reports Medicine 2024 (PMID:38340723; GSE227708):
   10 demultiplexed peripheral blood libraries and 2 matched Merkel cell
   carcinoma (`MCPyV`-driven) tumor libraries (`>=20` cells across 5 hashed 10x
   runs / 4 patients `Z1210`, `Z1249`, `Z1308`, `Z1493`), linked to pMHC
   multimer-verified `MCPyV` antitumor vs `CMV`/`EBV`/`Flu` bystander viral TCRs
   in Table S4. Evaluates unmodified `tcrsift prioritize --context blood-tumor`
   (`NeoTCR_PBL`), `--context blood`, blood-to-tumor clonal overlap (`in_tumor`),
   and `--exclude-known-viral` (`strict_ab` vs `ab_with_partial`).
2. Chiffelle et al., Immunity 2024 / Barras et al., Science Immunology 2024
   (PMID:39321802; GSE222448):
   Paired baseline pre-ACT (`T0`, 8 biopsies) and Day-30 post-ACT (`T30`, 6
   biopsies) melanoma tumors across Patients 1, 5, 7, 8, 10, and 13, linked to
   autologous tumor-reactive (`TR`) vs non-reactive (`NTR`) TCRs and infusion
   product expansion (`Clonally expanded in ACTP`) in Table S3.
3. Veatch et al., Cancer Cell 2022 (PMID:35413271; GSE198265):
   Human melanoma scRNA/AbSeq + TCR profiling across Patients X197, X198, X205,
   and X422 with experimentally verified autologous neoantigen-specific (`CD4_neo`)
   and melanoma self-antigen (`CD4_self`) CD4+ TIL clonotypes.
4. Oh et al., Cell 2020 (PMID:32497499; GSE149652):
   Whole-transcriptome 10x scRNA-seq + paired alpha/beta TCR sequencing of human
   bladder cancer CD4+ TILs and matched adjacent non-malignant bladder CD4+ T
   cells across 7 patients (13 libraries), evaluating unmodified
   `tcrsift prioritize --context solid-tumor --tcell-type cd4` (`NeoTCR4` 40/40
   genes), `pre_sorted=CD4` coreceptor-dropout rescue, and `--exclude-treg`
   (`FOXP3+` Treg exclusion).
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent


def is_current(out: Path) -> bool:
    from tcrsift import __version__

    cfg = out / "prioritization.json"
    return cfg.exists() and json.loads(cfg.read_text()).get("version") == __version__


def auroc(scores, positive) -> float:
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
    ok = np.isfinite(scores)
    scores, positive = scores[ok], positive[ok]
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if not n_pos or not n_neg:
        return float("nan")
    ranks = rankdata(scores)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auroc_ci(scores, positive, draws: int = 2000, seed: int = 0) -> tuple[float, float]:
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
    ok = np.isfinite(scores)
    pos, neg = scores[ok & positive], scores[ok & ~positive]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(draws):
        p = rng.choice(pos, len(pos))
        n = rng.choice(neg, len(neg))
        vals.append(auroc(np.concatenate([p, n]), np.arange(len(p) + len(n)) < len(p)))
    low, high = np.quantile(vals, [0.025, 0.975])
    return float(low), float(high)


def eval_auroc_block(
    df: pd.DataFrame,
    pos_mask: np.ndarray,
    signatures: list[str],
    *,
    col_fmt: str = "signature_{}_percentile_max",
    size_col: str = "size_percentile",
    **tags,
) -> list[dict]:
    rows = []
    pos_mask = np.asarray(pos_mask, dtype=bool)
    for sig in [*signatures, "clone_size"]:
        col = size_col if sig == "clone_size" else col_fmt.format(sig)
        if col not in df.columns:
            continue
        vals = df[col].to_numpy(dtype=float)
        ok = np.isfinite(vals)
        est = auroc(vals[ok], pos_mask[ok])
        low, high = auroc_ci(vals[ok], pos_mask[ok])
        rows.append(
            {
                **tags,
                "signature": sig,
                "n_positive": int(pos_mask[ok].sum()),
                "n_negative": int((~pos_mask[ok]).sum()),
                "auroc": round(est, 4) if np.isfinite(est) else np.nan,
                "auroc_low": round(low, 4) if np.isfinite(low) else np.nan,
                "auroc_high": round(high, 4) if np.isfinite(high) else np.nan,
            }
        )
    return rows


def collapse_by_beta(table: pd.DataFrame, any_cols: list[str], label_col: str = "label", untested_val: str = "untested") -> pd.DataFrame:
    """Collapse beta-labeled clones across alpha pairings (preserving untested clones)."""
    labelled = table[table[label_col].ne(untested_val)].sort_values(
        ["cell_count", "CDR3_alpha"], ascending=[False, True], kind="mergesort"
    )
    if labelled.empty:
        return table.copy()
    grouped = labelled.groupby(["donor", "CDR3_beta"], sort=False)
    collapsed = grouped.head(1).drop(columns=any_cols).merge(
        grouped[any_cols].any().reset_index(), on=["donor", "CDR3_beta"]
    )
    if "selection_rank" in labelled.columns and "selected_signature" in labelled.columns:
        sel = labelled[labelled["selection_rank"].notna()].sort_values(
            ["selection_rank", "CDR3_alpha"], ascending=[True, True], kind="mergesort"
        )
        if not sel.empty:
            best_sel = sel.drop_duplicates(["donor", "CDR3_beta"])[
                ["donor", "CDR3_beta", "selection_rank", "selected_signature"]
            ].rename(columns={"selection_rank": "_best_rank", "selected_signature": "_best_sig"})
            collapsed = collapsed.merge(best_sel, on=["donor", "CDR3_beta"], how="left")
            collapsed["selection_rank"] = collapsed["selection_rank"].combine_first(collapsed.pop("_best_rank"))
            collapsed["selected_signature"] = collapsed["selected_signature"].combine_first(collapsed.pop("_best_sig"))
    return pd.concat([table[table[label_col].eq(untested_val)], collapsed], ignore_index=True)


# ---------------------------------------------------- 1. Pulliam 2024 Blood + Tumor


def demux_pulliam_all(data_dir: Path, work_dir: Path):
    """Demultiplex each individual Blood and Tumor hashtag tube into its own 10x sample."""
    root = data_dir / "pulliam"
    meta = pd.read_csv(root / "meta_data.csv.gz", encoding="utf-8-sig")
    blood_samples, tumor_samples, demux_counts = [], [], {}

    for run_dir in sorted((root / "runs").iterdir()):
        if not run_dir.is_dir():
            continue
        run = run_dir.name
        tubes = meta[meta["Run_id"].eq(run)]
        features = pd.read_csv(run_dir / "features.tsv.gz", sep="\t", header=None)
        barcodes = pd.read_csv(run_dir / "barcodes.tsv.gz", header=None)[0].to_numpy()
        matrix = scipy.io.mmread(run_dir / "matrix.mtx.gz").tocsr()
        tags = features.index[features[0].str.startswith("Hashtag")].to_numpy()
        hto = matrix[tags].toarray().T
        ordered = np.sort(hto, axis=1)
        clear = (ordered[:, -1] >= 10) & (ordered[:, -1] >= 3 * np.maximum(ordered[:, -2], 1))
        top_num = features.loc[tags[hto.argmax(axis=1)], 0].str[len("Hashtag") :].astype(int).to_numpy()
        assigned_ht = np.where(clear, top_num, 0)

        genes_mask = features[2].eq("Gene Expression").to_numpy()
        clean_features = features[genes_mask].replace(r"^GRCh38_", "", regex=True)
        contigs = pd.read_csv(run_dir / "all_contig_annotations.csv.gz", low_memory=False)
        called = contigs["is_cell"].astype(str).str.lower().eq("true")
        confident = contigs["high_confidence"].astype(str).str.lower().eq("true")
        valid_contigs = contigs[called & confident]

        demux_counts[run] = {}
        for _, trow in tubes.iterrows():
            ht = int(trow["Hashtag"])
            tis = str(trow["tissue"]).strip()
            pat = str(trow["Patient_id"]).strip()
            tube_id = str(trow["Tube"]).strip()
            note = str(trow.get("Note") or "").strip()
            cell_mask = assigned_ht == ht
            n_cells = int(cell_mask.sum())
            demux_counts[run][f"H{ht}_{tube_id}_{tis}"] = n_cells
            if n_cells < 20:
                continue
            sname = f"{run}_H{ht}_{tube_id}_{tis.lower()}"
            sdir = work_dir / "pulliam_demux" / sname
            gex_dir = sdir / "gex" / "filtered_feature_bc_matrix"
            if not (gex_dir / "matrix.mtx.gz").exists():
                gex_dir.mkdir(parents=True, exist_ok=True)
                (sdir / "vdj").mkdir(parents=True, exist_ok=True)
                with gzip.open(gex_dir / "matrix.mtx.gz", "wb") as handle:
                    scipy.io.mmwrite(handle, matrix[genes_mask][:, cell_mask].tocoo(), field="integer")
                clean_features.to_csv(gex_dir / "features.tsv.gz", sep="\t", header=False, index=False)
                pd.Series(barcodes[cell_mask]).to_csv(gex_dir / "barcodes.tsv.gz", header=False, index=False)
                valid_contigs[valid_contigs["barcode"].isin(barcodes[cell_mask])].to_csv(
                    sdir / "vdj" / "filtered_contig_annotations.csv", index=False
                )
            rec = {"sample": sname, "patient_id": pat, "tissue": tis.lower(), "note": note, "dir": sdir}
            if tis.lower() == "blood":
                blood_samples.append(rec)
            elif tis.lower() == "tumor":
                tumor_samples.append(rec)

    return blood_samples, tumor_samples, demux_counts


def pulliam_s4_labels(data_dir: Path) -> pd.DataFrame:
    root = data_dir / "pulliam"
    meta = pd.read_csv(root / "meta_data.csv.gz", encoding="utf-8-sig")
    plain = lambda s: s.astype(str).str.lower().str.replace(r"[^a-z0-9]", "", regex=True)  # noqa: E731
    patient_map = dict(zip(plain(meta["Tube"]), meta["Patient_id"]))
    raw = pd.read_excel(root / "table_s4.xlsx", header=1)
    df = pd.DataFrame(
        {
            "donor": plain(raw["subject"]).map(patient_map),
            "CDR3_beta": raw["cdr3_b_aa"].astype(str).str.strip(),
            "CDR3_alpha_s4": raw["cdr3_a_aa"].astype(str).str.strip(),
            "virus": raw["virus"].astype(str).str.strip(),
        }
    ).dropna(subset=["donor"])
    df = df[df["virus"].isin(["MCPyV", "CMV", "EBV", "Flu"])].copy()
    df["label"] = np.where(df["virus"].eq("MCPyV"), "MCPyV", "control_virus")
    # If a beta has conflicting labels within a donor, drop it cleanly
    n_lbl = df.groupby(["donor", "CDR3_beta"])["label"].transform("nunique")
    return df[n_lbl == 1].drop_duplicates(["donor", "CDR3_beta"])[["donor", "CDR3_beta", "label", "virus"]]


def run_pulliam_blood(data_dir: Path, work_dir: Path, cache_dir: Path | None):
    from tcrsift.cli import main as tcrsift_main
    from tcrsift.loader import _pivot_vdj_by_barcode, load_cellranger_vdj

    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []
    blood_samples, tumor_samples, demux_counts = demux_pulliam_all(data_dir, work_dir)
    labels = pulliam_s4_labels(data_dir)

    # Collect tumor TCR beta & alpha-beta repertoire per patient from matched Tumor hashtags
    tum_betas_by_donor: dict[str, set[str]] = {}
    tum_abs_by_donor: dict[str, set[str]] = {}
    for ts in tumor_samples:
        vdf = load_cellranger_vdj(ts["dir"] / "vdj", ts["sample"], verbose=False)
        piv = _pivot_vdj_by_barcode(vdf)
        pat = ts["patient_id"]
        tum_betas_by_donor.setdefault(pat, set()).update(
            piv["CDR3_beta"].dropna().astype(str).loc[lambda s: s != ""].tolist()
        )
        paired = piv[
            piv["CDR3_alpha"].notna()
            & (piv["CDR3_alpha"] != "")
            & piv["CDR3_beta"].notna()
            & (piv["CDR3_beta"] != "")
        ]
        tum_abs_by_donor.setdefault(pat, set()).update(paired["CDR3ab"].dropna().astype(str).tolist())

    blood_sheet = work_dir / "pulliam_blood_sheet.csv"
    lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    for bs in blood_samples:
        lines.append(f"{bs['sample']},{bs['patient_id']},blood,{bs['dir'] / 'vdj'},{bs['dir'] / 'gex'}")
    blood_sheet.write_text("\n".join(lines) + "\n")

    runs = [
        ("blood_tumor_default", ["--context", "blood-tumor"]),
        ("blood_tumor_min1", ["--context", "blood-tumor", "--min-cells", "1"]),
        ("blood_default", ["--context", "blood"]),
        ("blood_tumor_excl_viral_strict", ["--context", "blood-tumor", "--exclude-known-viral", "--database-match", "strict_ab"]),
        ("blood_tumor_excl_viral_partial", ["--context", "blood-tumor", "--exclude-known-viral", "--database-match", "ab_with_partial"]),
    ]
    auroc_rows, shortlist_rows, facts = [], [], {}

    for rname, extra in runs:
        out = work_dir / f"pulliam_{rname}"
        if not is_current(out):
            tcrsift_main(["prioritize", str(blood_sheet), "-o", str(out), *extra, *cache_args])
        cfg = json.loads((out / "prioritization.json").read_text())
        sigs = cfg["resolved_signatures"]
        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones.groupby("donor")["cell_count"].rank(pct=True)
        clones["in_matched_tumor_ab"] = [
            ab in tum_abs_by_donor.get(d, set()) for d, ab in zip(clones["donor"], clones["CDR3ab"])
        ]
        clones["in_matched_tumor_beta"] = [
            b in tum_betas_by_donor.get(d, set()) for d, b in zip(clones["donor"], clones["CDR3_beta"])
        ]
        clones = clones.merge(labels, on=["donor", "CDR3_beta"], how="left")
        clones["label"] = clones["label"].fillna("untested")

        any_cols = [
            "eligible_for_review",
            "selected_for_review",
            "in_matched_tumor_ab",
            "in_matched_tumor_beta",
            "known_viral_match",
            *[f"signature_{s}_passes_cutoff_max" for s in sigs],
        ]
        collapsed = collapse_by_beta(clones, any_cols)

        if rname == "blood_tumor_min1":
            # Evaluate AUROCs across all scored Table S4 clones (>=1 cell) AND expanded (>=2 cells)
            s4_all = collapsed[collapsed["label"].isin(["MCPyV", "control_virus"])].copy()
            auroc_rows += eval_auroc_block(
                s4_all,
                s4_all["label"].eq("MCPyV").to_numpy(),
                sigs,
                dataset="pulliam_blood",
                subset="all_clones_ge1",
                contrast="MCPyV_vs_CMV_EBV_Flu",
            )
            s4_ge2 = s4_all[s4_all["cell_count"] >= 2].copy()
            auroc_rows += eval_auroc_block(
                s4_ge2,
                s4_ge2["label"].eq("MCPyV").to_numpy(),
                sigs,
                dataset="pulliam_blood",
                subset="expanded_ge2",
                contrast="MCPyV_vs_CMV_EBV_Flu",
            )
            # MCPyV vs unlabeled >=2-cell blood clones in donors with MCPyV labels
            pos_donors = collapsed.loc[collapsed["label"].eq("MCPyV"), "donor"].unique()
            exp_unlab = collapsed[
                (collapsed["cell_count"] >= 2)
                & collapsed["donor"].isin(pos_donors)
                & collapsed["label"].isin(["MCPyV", "untested"])
            ].copy()
            auroc_rows += eval_auroc_block(
                exp_unlab,
                exp_unlab["label"].eq("MCPyV").to_numpy(),
                sigs,
                dataset="pulliam_blood",
                subset="expanded_ge2",
                contrast="MCPyV_vs_untested_blood",
            )
            # Blood-to-tumor clonal overlap (`in_matched_tumor_ab`) restricted to donors
            # that actually have sequenced matched tumor hashtag libraries (Z1249 & Z1308)
            tum_donor_clones = clones[clones["donor"].isin(tum_abs_by_donor)].copy()
            auroc_rows += eval_auroc_block(
                tum_donor_clones,
                tum_donor_clones["in_matched_tumor_ab"].to_numpy(dtype=bool),
                sigs,
                dataset="pulliam_blood",
                subset="matched_tumor_donors_ge1",
                contrast="tumor_infiltrating_ab_vs_blood_only",
            )
            tum_donor_ge2 = tum_donor_clones[tum_donor_clones["cell_count"] >= 2].copy()
            auroc_rows += eval_auroc_block(
                tum_donor_ge2,
                tum_donor_ge2["in_matched_tumor_ab"].to_numpy(dtype=bool),
                sigs,
                dataset="pulliam_blood",
                subset="matched_tumor_donors_ge2",
                contrast="tumor_infiltrating_ab_vs_blood_only",
            )

        sel = collapsed["selected_for_review"].to_numpy(dtype=bool)
        sel_raw = clones["selected_for_review"].to_numpy(dtype=bool)
        shortlist_rows.append(
            {
                "cohort": "pulliam_blood",
                "run": rname,
                "retained_cells": int(cfg["cell_filters"][-1]["after"]),
                "total_clones": len(clones),
                "expanded_ge2": int((clones["cell_count"] >= 2).sum()),
                "shortlisted": int(sel_raw.sum()),
                "mcpyv_total_ge2": int(((collapsed["label"] == "MCPyV") & (collapsed["cell_count"] >= 2)).sum()),
                "mcpyv_shortlisted": int((sel & collapsed["label"].eq("MCPyV").to_numpy()).sum()),
                "control_viral_total_ge2": int(
                    ((collapsed["label"] == "control_virus") & (collapsed["cell_count"] >= 2)).sum()
                ),
                "control_viral_shortlisted": int((sel & collapsed["label"].eq("control_virus").to_numpy()).sum()),
                "tumor_infiltrating_ab_total_ge2": int(
                    (clones["in_matched_tumor_ab"] & (clones["cell_count"] >= 2)).sum()
                ),
                "tumor_infiltrating_ab_shortlisted": int(
                    (sel_raw & clones["in_matched_tumor_ab"].to_numpy(dtype=bool)).sum()
                ),
            }
        )
        facts[rname] = shortlist_rows[-1]

    return auroc_rows, shortlist_rows, {"demux_cells": demux_counts, "runs": facts}


# ---------------------------------------------------- 2. Chiffelle 2024 T0 vs T30 ACT


def chiffelle_s3_labels(data_dir: Path) -> pd.DataFrame:
    raw = pd.read_excel(data_dir / "chiffelle" / "table_s3.xlsx", header=2)
    raw["Patient"] = raw["Patient"].ffill()
    raw = raw[raw["Tumor-reactivity"].isin(["yes", "no"]) & raw.iloc[:, 1].notna()].copy()
    df = pd.DataFrame(
        {
            "donor": "P" + raw["Patient"].astype(int).astype(str),
            "CDR3_beta": raw.iloc[:, 1].astype(str).str.split("_").str[1].str[:-1],
            "label": np.where(raw["Tumor-reactivity"].eq("yes"), "TR", "NTR"),
            "in_actp": raw["Clonally expanded in ACTP"].notna(),
            "in_t0_s3": raw["Clonally expanded in T0"].notna(),
        }
    )
    n_lbl = df.groupby(["donor", "CDR3_beta"])["label"].transform("nunique")
    df = df[n_lbl == 1]
    return (
        df.groupby(["donor", "CDR3_beta"], as_index=False)
        .agg(label=("label", "first"), in_actp=("in_actp", "any"), in_t0_s3=("in_t0_s3", "any"))
    )


def run_chiffelle_act(data_dir: Path, work_dir: Path, cache_dir: Path | None):
    from tcrsift.cli import main as tcrsift_main

    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []
    sroot = data_dir / "chiffelle" / "samples"
    labels = chiffelle_s3_labels(data_dir)
    all_samples = sorted(p for p in sroot.iterdir() if p.is_dir())
    t0_samples = [p for p in all_samples if "_T0" in p.name]
    t30_samples = [p for p in all_samples if "_T30" in p.name]

    auroc_rows, shortlist_rows, facts = [], [], {}
    clones_by_tp = {}

    for tp_name, chosen in (("T0_pre_ACT", t0_samples), ("T30_post_ACT", t30_samples)):
        sheet = work_dir / f"chiffelle_{tp_name}.csv"
        lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
        for s in chosen:
            lines.append(f"{s.name},{s.name.split('_')[0]},tumor,{s / 'vdj'},{s / 'gex'}")
        sheet.write_text("\n".join(lines) + "\n")
        out = work_dir / f"chiffelle_{tp_name}"
        if not is_current(out):
            tcrsift_main(["prioritize", str(sheet), "-o", str(out), "--context", "solid-tumor", *cache_args])
        cfg = json.loads((out / "prioritization.json").read_text())
        sigs = cfg["resolved_signatures"]
        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones.groupby("donor")["cell_count"].rank(pct=True)
        clones = clones.merge(labels, on=["donor", "CDR3_beta"], how="left")
        clones["label"] = clones["label"].fillna("untested")
        clones["in_actp"] = clones["in_actp"].eq(True)

        any_cols = [
            "eligible_for_review",
            "selected_for_review",
            "in_actp",
            *[f"signature_{s}_passes_cutoff_max" for s in sigs],
        ]
        collapsed = collapse_by_beta(clones, any_cols)
        clones_by_tp[tp_name] = collapsed

        tested = collapsed[collapsed["label"].isin(["TR", "NTR"])].copy()
        auroc_rows += eval_auroc_block(
            tested,
            tested["label"].eq("TR").to_numpy(),
            sigs,
            dataset="chiffelle_act_6patients",
            subset=tp_name,
            contrast="TR_vs_NTR",
        )
        tested_actp = tested[tested["in_actp"]].copy()
        if (tested_actp["label"] == "TR").sum() >= 2 and (tested_actp["label"] == "NTR").sum() >= 2:
            auroc_rows += eval_auroc_block(
                tested_actp,
                tested_actp["label"].eq("TR").to_numpy(),
                sigs,
                dataset="chiffelle_act_6patients",
                subset=f"{tp_name}_ACTP_expanded",
                contrast="TR_vs_NTR_in_ACTP",
            )

        sel = collapsed["selected_for_review"].to_numpy(dtype=bool)
        tr_mask = collapsed["label"].eq("TR").to_numpy()
        ntr_mask = collapsed["label"].eq("NTR").to_numpy()
        shortlist_rows.append(
            {
                "cohort": "chiffelle_act",
                "run": tp_name,
                "retained_cells": int(cfg["cell_filters"][-1]["after"]),
                "total_clones": len(clones),
                "expanded_ge2": int((clones["cell_count"] >= 2).sum()),
                "shortlisted": int(clones["selected_for_review"].sum()),
                "tr_scored": int(tr_mask.sum()),
                "tr_shortlisted": int((sel & tr_mask).sum()),
                "ntr_scored": int(ntr_mask.sum()),
                "ntr_shortlisted": int((sel & ntr_mask).sum()),
                "tested_shortlist_precision": round(
                    int((sel & tr_mask).sum()) / max(int((sel & (tr_mask | ntr_mask)).sum()), 1), 4
                ),
            }
        )
        facts[tp_name] = shortlist_rows[-1]

    # Persistence of T0 clones into T30
    t0_b = set(zip(clones_by_tp["T0_pre_ACT"]["donor"], clones_by_tp["T0_pre_ACT"]["CDR3_beta"]))
    t30_df = clones_by_tp["T30_post_ACT"]
    persisted = [pair in t0_b for pair in zip(t30_df["donor"], t30_df["CDR3_beta"])]
    facts["t30_persisted_from_t0_tr"] = int((np.asarray(persisted) & t30_df["label"].eq("TR").to_numpy()).sum())
    facts["t30_persisted_from_t0_ntr"] = int((np.asarray(persisted) & t30_df["label"].eq("NTR").to_numpy()).sum())
    facts["t30_newly_detected_tr"] = int((~np.asarray(persisted) & t30_df["label"].eq("TR").to_numpy()).sum())
    facts["t30_newly_detected_ntr"] = int((~np.asarray(persisted) & t30_df["label"].eq("NTR").to_numpy()).sum())
    return auroc_rows, shortlist_rows, facts


# ---------------------------------------------------- 3. Veatch 2022 Melanoma CD4+ TILs


def run_veatch_cd4(data_dir: Path):
    from tcrsift.prioritize import (
        PRIORITY_CONTEXT_SIGNATURES,
        REGISTRY,
        resolve_signatures,
        select_round_robin,
        signature_pass_mask,
    )
    from tcrsift.signature_background import calibrate_signature_background
    from tcrsift.signature_methods import score_by_name

    vroot = data_dir / "veatch"
    sigs = resolve_signatures("solid-tumor", "cd4", species="human")

    adatas = []
    panel_coverage: dict[str, str] = {}
    for pat in ("X197", "X198", "X205", "X422"):
        pdir = vroot / pat
        vdf = pd.read_csv(next(pdir.glob("*VDJ*")))
        vdf = vdf[
            vdf["TCRVb?"].eq(True)
            & vdf["Cell_Type_Experimental"].isin(["T_CD4_memory", "T_CD4_naive"])
            & vdf["TCR_Beta_Delta_CDR3_Translation_Dominant"].notna()
        ].copy()
        if "sample" in vdf.columns:
            vdf = vdf[~vdf["sample"].isin(["Multiplet", "Undetermined"])].copy()
            # In X198, restrict to tumor harvest (`8_26_main_harvest`) for TIL scoring;
            # keep X422 tumorA and tumorB as distinct tumor replicates!
            if pat == "X198":
                vdf = vdf[vdf["sample"].eq("8_26_main_harvest")].copy()
                vdf["sample_id"] = "X198_tumor"
            else:
                vdf["sample_id"] = pat + "_" + vdf["sample"].astype(str)
        else:
            vdf["sample_id"] = f"{pat}_tumor"

        gfile = next(pdir.glob("*MolsPerCell*"))
        gdf = pd.read_csv(gfile)
        mrna_cols = [c for c in gdf.columns if c != "Cell_Index" and not c.endswith("|pAbO")]
        abseq_cols = [c for c in gdf.columns if c.endswith("|pAbO")]
        mrna_genes = [c.split("|")[0].strip() for c in mrna_cols]

        if not panel_coverage:
            mrna_upper = {g.upper() for g in mrna_genes}
            for sname in sigs:
                all_g = REGISTRY[sname].all_genes
                n_pres = sum(1 for g in all_g if g.upper() in mrna_upper)
                panel_coverage[sname] = f"{n_pres}/{len(all_g)}"

        joined = vdf.merge(gdf, on="Cell_Index", how="inner")
        X = sp.csr_matrix(joined[mrna_cols].to_numpy(dtype=float))
        fox_col = next((c for c, g in zip(mrna_cols, mrna_genes) if g.upper() == "FOXP3"), None)
        fox_pos = (joined[fox_col].to_numpy(dtype=float) > 0) if fox_col else np.zeros(len(joined), dtype=bool)

        is_neo = (
            joined["CD4_neo"].eq(1.0).to_numpy()
            if "CD4_neo" in joined.columns
            else joined["Neoantigen"].notna().to_numpy()
        )
        if "CD4_self" in joined.columns:
            is_self = joined["CD4_self"].eq(1.0).to_numpy()
        elif "self" in joined.columns:
            # Exclude X197 CD8 self annotations (`magea3_cd8`, `tyro_cd8`)
            sval = joined["self"].fillna("").astype(str)
            is_self = (sval != "") & (~sval.str.contains("cd8", case=False)).to_numpy()
        else:
            is_self = np.zeros(len(joined), dtype=bool)

        label_arr = np.where(is_neo, "CD4_neo", np.where(is_self, "CD4_self", "untested"))
        # Extract PD-1 and CD39 AbSeq log1p protein if present
        pd1_ab = next((c for c in abseq_cols if "PDCD1|" in c), None)
        cd39_ab = next((c for c in abseq_cols if "ENTPD1|" in c), None)

        obs = pd.DataFrame(
            {
                "donor": pat,
                "sample": joined["sample_id"].to_numpy(),
                "lineage": "cd4",
                "CDR3ab": joined["TCR_Beta_Delta_CDR3_Translation_Dominant"].astype(str).str.strip().to_numpy(),
                "cell_label": label_arr,
                "foxp3_pos": fox_pos.astype(float),
                "abseq_pd1_log1p": np.log1p(joined[pd1_ab].to_numpy(dtype=float)) if pd1_ab else np.nan,
                "abseq_cd39_log1p": np.log1p(joined[cd39_ab].to_numpy(dtype=float)) if cd39_ab else np.nan,
            },
            index=[f"{pat}_{cid}" for cid in joined["Cell_Index"]],
        )
        part = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=mrna_genes))
        part.var_names_make_unique()
        sc.pp.normalize_total(part, target_sum=10_000)
        sc.pp.log1p(part)
        for sname in sigs:
            part.obs[f"signature_{sname}"] = score_by_name(
                part,
                sname,
                log1p=False,
                groups=part.obs["sample"],
                on_missing="ignore",
                species="human",
            ).to_numpy(dtype=float)
        adatas.append(part)

    full = ad.concat(adatas, join="inner")
    stratum = ["donor", "sample", "lineage"]
    sig_cols = [f"signature_{s}" for s in sigs]
    per_cell = full.obs[[*stratum, "CDR3ab", *sig_cols]].copy()
    grp = full.obs.groupby([*stratum, "CDR3ab"], observed=True)
    scores = grp[sig_cols].mean().join(grp.size().rename("cells")).reset_index()
    scores["frequency"] = scores["cells"] / scores.groupby(stratum, observed=True)["cells"].transform("sum")
    for col in sig_cols:
        scores[f"{col}_percentile"] = scores.groupby(stratum, observed=True)[col].rank(method="average", pct=True)
    scores, _ = calibrate_signature_background(scores, per_cell, sigs, quantile=0.99)
    for sname in sigs:
        scores[f"signature_{sname}_passes_cutoff"] = signature_pass_mask(scores, sname, evidence="either")

    # Aggregate to donor-level clone table (matching X422 tumorA + tumorB)
    c_grp = full.obs.groupby(["donor", "CDR3ab"], observed=True)
    clones = c_grp.agg(
        cell_count=("sample", "size"),
        foxp3_pos_fraction=("foxp3_pos", "mean"),
        label=(
            "cell_label",
            lambda s: "CD4_neo"
            if (s == "CD4_neo").any()
            else ("CD4_self" if (s == "CD4_self").any() else "untested"),
        ),
        abseq_pd1_log1p=("abseq_pd1_log1p", "mean"),
        abseq_cd39_log1p=("abseq_cd39_log1p", "mean"),
    ).reset_index()
    max_cols = [f"signature_{s}_percentile" for s in sigs] + [f"signature_{s}_passes_cutoff" for s in sigs]
    sig_max = (
        scores.groupby(["donor", "CDR3ab"], observed=True)[[*max_cols, "frequency"]]
        .max()
        .rename(columns={"frequency": "max_frequency", **{c: f"{c}_max" for c in max_cols}})
        .reset_index()
    )
    clones = clones.merge(sig_max, on=["donor", "CDR3ab"], how="left")
    clones["size_percentile"] = clones.groupby("donor")["cell_count"].rank(pct=True)
    clones["is_treg_like"] = clones["foxp3_pos_fraction"] >= 0.50
    clones["any_sig_pass"] = clones[[f"signature_{s}_passes_cutoff_max" for s in sigs]].any(axis=1)

    exp = clones[clones["cell_count"] >= 2].copy()
    auroc_rows = []
    # Contrast A: Verified CD4_neo vs unlabeled >=2-cell intratumoral CD4 clones
    sub_neo_unlab = exp[exp["label"].isin(["CD4_neo", "untested"])].copy()
    auroc_rows += eval_auroc_block(
        sub_neo_unlab,
        sub_neo_unlab["label"].eq("CD4_neo").to_numpy(),
        sigs,
        dataset="veatch_melanoma_cd4",
        subset="expanded_ge2",
        contrast="CD4_neo_vs_untested_cd4",
    )
    # Contrast B: Verified CD4_neo vs non-Treg (`foxp3 < 0.5`) unlabeled >=2-cell CD4 clones
    sub_neo_notreg = exp[
        exp["label"].eq("CD4_neo") | (exp["label"].eq("untested") & ~exp["is_treg_like"])
    ].copy()
    auroc_rows += eval_auroc_block(
        sub_neo_notreg,
        sub_neo_notreg["label"].eq("CD4_neo").to_numpy(),
        sigs,
        dataset="veatch_melanoma_cd4",
        subset="expanded_ge2_no_treg",
        contrast="CD4_neo_vs_non_treg_untested",
    )
    # Contrast C: X205 CD4_neo vs CD4_self (both tumor-antigen-specific CD4+ TILs in the same tumor!)
    sub_x205 = exp[(exp["donor"] == "X205") & exp["label"].isin(["CD4_neo", "CD4_self"])].copy()
    if (sub_x205["label"] == "CD4_neo").sum() >= 2 and (sub_x205["label"] == "CD4_self").sum() >= 2:
        auroc_rows += eval_auroc_block(
            sub_x205,
            sub_x205["label"].eq("CD4_neo").to_numpy(),
            sigs,
            dataset="veatch_melanoma_cd4",
            subset="X205_expanded_ge2",
            contrast="CD4_neo_vs_CD4_self",
        )

    # Round-robin selection with and without `--exclude-treg`
    priority_sigs = [s for s in PRIORITY_CONTEXT_SIGNATURES["solid-tumor"] if s in sigs]
    shortlist_rows = []
    for mode_name, excl_treg in (("default_keep_treg", False), ("exclude_treg_0.50", True)):
        c_in = clones.copy()
        c_in["eligible_for_review"] = (
            (c_in["cell_count"] >= 2)
            & (c_in["max_frequency"] >= 0.001)
            & c_in["any_sig_pass"]
            & (~c_in["is_treg_like"] if excl_treg else True)
        )
        picked = select_round_robin(
            c_in,
            scores,
            sigs,
            max_clones=200,
            evidence="either",
            priority_signatures=priority_sigs,
        )
        sel = picked["selected_for_review"].to_numpy(dtype=bool)
        shortlist_rows.append(
            {
                "cohort": "veatch_melanoma_cd4",
                "run": mode_name,
                "retained_cells": int(full.n_obs),
                "total_clones": len(picked),
                "expanded_ge2": int((picked["cell_count"] >= 2).sum()),
                "shortlisted": int(sel.sum()),
                "cd4_neo_ge2": int(((picked["label"] == "CD4_neo") & (picked["cell_count"] >= 2)).sum()),
                "cd4_neo_shortlisted": int((sel & picked["label"].eq("CD4_neo").to_numpy()).sum()),
                "cd4_self_ge2": int(((picked["label"] == "CD4_self") & (picked["cell_count"] >= 2)).sum()),
                "cd4_self_shortlisted": int((sel & picked["label"].eq("CD4_self").to_numpy()).sum()),
                "treg_ge2": int((picked["is_treg_like"] & (picked["cell_count"] >= 2)).sum()),
                "treg_shortlisted": int((sel & picked["is_treg_like"].to_numpy()).sum()),
            }
        )

    facts = {
        "targeted_mrna_genes": int(full.n_vars),
        "signature_panel_coverage": panel_coverage,
        "retained_cd4_cells": int(full.n_obs),
        "total_clones": len(clones),
        "expanded_clones_ge2": len(exp),
        "cd4_neo_clones_ge2": int((exp["label"] == "CD4_neo").sum()),
        "cd4_self_clones_ge2": int((exp["label"] == "CD4_self").sum()),
        "treg_clones_ge2": int(exp["is_treg_like"].sum()),
        "runs": {r["run"]: r for r in shortlist_rows},
    }
    return auroc_rows, shortlist_rows, facts


# ---------------------------------------------------- 4. Oh 2020 Bladder Cancer CD4+ TILs


def run_oh_cd4(data_dir: Path, work_dir: Path, cache_dir: Path | None):
    from tcrsift.cli import main as tcrsift_main
    from tcrsift.loader import _pivot_vdj_by_barcode, load_cellranger_vdj

    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []
    sroot = data_dir / "oh" / "samples"
    tumor_dirs = sorted(p for p in sroot.iterdir() if p.is_dir() and p.name.endswith("_tumor"))
    normal_dirs = sorted(p for p in sroot.iterdir() if p.is_dir() and p.name.endswith("_normal"))

    # Map author per-cell CD4 subtype (`cell_types` in cellinfo.csv.gz) + normal-bladder presence per (donor, CDR3ab)
    normal_clones_by_donor: dict[str, set[str]] = {}
    for nd in normal_dirs:
        donor = nd.name[: -len("_normal")]
        vdf = load_cellranger_vdj(nd / "vdj", nd.name, verbose=False)
        piv = _pivot_vdj_by_barcode(vdf)
        paired = piv[
            piv["CDR3_alpha"].notna()
            & (piv["CDR3_alpha"] != "")
            & piv["CDR3_beta"].notna()
            & (piv["CDR3_beta"] != "")
        ]
        normal_clones_by_donor.setdefault(donor, set()).update(paired["CDR3ab"].dropna().astype(str).tolist())

    subtype_rows = []
    for td in tumor_dirs:
        donor = td.name[: -len("_tumor")]
        ci = pd.read_csv(td / "cellinfo.csv.gz")
        ci["barcode"] = ci["index"].astype(str).str.split("-").str[0] + "-1"
        vdf = load_cellranger_vdj(td / "vdj", td.name, verbose=False)
        piv = _pivot_vdj_by_barcode(vdf).reset_index()
        m = piv.merge(ci[["barcode", "cell_types"]], on="barcode", how="inner")
        m = m[m["CDR3ab"].notna() & (m["CDR3ab"] != "")]
        for cdr3ab, grp in m.groupby("CDR3ab"):
            vc = grp["cell_types"].value_counts()
            n_tot = len(grp)
            n_eff = int(grp["cell_types"].isin(["CD4CXCL13", "CD4GZMB", "CD4GZMK"]).sum())
            n_cxcl13 = int((grp["cell_types"] == "CD4CXCL13").sum())
            n_cyto = int(grp["cell_types"].isin(["CD4GZMB", "CD4GZMK"]).sum())
            n_treg = int(grp["cell_types"].isin(["CD4IL2RAHI", "CD4IL2RALO"]).sum())
            n_treg_hi = int((grp["cell_types"] == "CD4IL2RAHI").sum())
            subtype_rows.append(
                {
                    "donor": donor,
                    "CDR3ab": cdr3ab,
                    "dom_subtype": str(vc.index[0]),
                    "eff_cxcl13_cyto_frac": n_eff / n_tot,
                    "cxcl13_frac": n_cxcl13 / n_tot,
                    "cyto_frac": n_cyto / n_tot,
                    "treg_hi_frac": n_treg_hi / n_tot,
                    "treg_any_frac": n_treg / n_tot,
                }
            )
    sub_df = pd.DataFrame(subtype_rows)

    sheet_pre = work_dir / "oh_tumor_presorted.csv"
    sheet_raw = work_dir / "oh_tumor_unsorted.csv"
    lines_p = ["sample,patient_id,tissue,pre_sorted,vdj_dir,gex_dir"]
    lines_r = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    for td in tumor_dirs:
        donor = td.name[: -len("_tumor")]
        lines_p.append(f"{td.name},{donor},tumor,CD4,{td / 'vdj'},{td / 'gex'}")
        lines_r.append(f"{td.name},{donor},tumor,{td / 'vdj'},{td / 'gex'}")
    sheet_pre.write_text("\n".join(lines_p) + "\n")
    sheet_raw.write_text("\n".join(lines_r) + "\n")

    runs = [
        ("presorted_default", sheet_pre, []),
        ("presorted_exclude_treg_0.50", sheet_pre, ["--exclude-treg", "--treg-foxp3-fraction", "0.5"]),
        ("presorted_exclude_treg_0.25", sheet_pre, ["--exclude-treg", "--treg-foxp3-fraction", "0.25"]),
        ("no_presorted_cd4_mrna_only", sheet_raw, []),
    ]
    auroc_rows, shortlist_rows, facts = [], [], {}

    for rname, sh, extra in runs:
        out = work_dir / f"oh_{rname}"
        if not is_current(out):
            tcrsift_main(
                [
                    "prioritize",
                    str(sh),
                    "-o",
                    str(out),
                    "--context",
                    "solid-tumor",
                    "--tcell-type",
                    "cd4",
                    *extra,
                    *cache_args,
                ]
            )
        cfg = json.loads((out / "prioritization.json").read_text())
        sigs = cfg["resolved_signatures"]
        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones.groupby("donor")["cell_count"].rank(pct=True)
        clones = clones.merge(sub_df, on=["donor", "CDR3ab"], how="left")
        clones["in_matched_normal"] = [
            ab in normal_clones_by_donor.get(d, set()) for d, ab in zip(clones["donor"], clones["CDR3ab"])
        ]

        if rname == "presorted_default":
            exp = clones[clones["cell_count"] >= 2].copy()
            is_eff = exp["dom_subtype"].isin(["CD4CXCL13", "CD4GZMB", "CD4GZMK"])
            is_cxcl13 = exp["dom_subtype"].eq("CD4CXCL13")
            is_treg = exp["dom_subtype"].eq("CD4IL2RAHI")
            is_cm_th17 = exp["dom_subtype"].isin(["CD4CM", "CD4TH17"])

            # 1. Effector (`CD4CXCL13` + `CD4GZMB`/`CD4GZMK`) vs resting bystander (`CD4CM` + `CD4TH17`)
            sub1 = exp[is_eff | is_cm_th17].copy()
            auroc_rows += eval_auroc_block(
                sub1,
                is_eff.loc[sub1.index].to_numpy(),
                sigs,
                dataset="oh_bladder_cd4",
                subset="expanded_ge2",
                contrast="CD4CXCL13_GZMB_GZMK_vs_CD4CM_TH17",
            )
            # 2. `CD4CXCL13` specifically vs all non-Treg `CD4CM`/`CD4TH17` clones
            sub2 = exp[is_cxcl13 | is_cm_th17].copy()
            auroc_rows += eval_auroc_block(
                sub2,
                is_cxcl13.loc[sub2.index].to_numpy(),
                sigs,
                dataset="oh_bladder_cd4",
                subset="expanded_ge2",
                contrast="CD4CXCL13_vs_CD4CM_TH17",
            )
            # 3. Effector (`CD4CXCL13` + `CD4GZMB`/`CD4GZMK`) vs `CD4IL2RAHI` Treg clones (quantifies Treg overlap!)
            sub3 = exp[is_eff | is_treg].copy()
            auroc_rows += eval_auroc_block(
                sub3,
                is_eff.loc[sub3.index].to_numpy(),
                sigs,
                dataset="oh_bladder_cd4",
                subset="expanded_ge2",
                contrast="CD4CXCL13_GZMB_GZMK_vs_CD4IL2RAHI_Treg",
            )
            # 4. Tumor-exclusive non-Treg expanded clones vs matched-normal-shared bystander clones
            donors_with_norm = set(normal_clones_by_donor)
            sub4 = exp[exp["donor"].isin(donors_with_norm) & ~is_treg].copy()
            auroc_rows += eval_auroc_block(
                sub4,
                (~sub4["in_matched_normal"]).to_numpy(),
                sigs,
                dataset="oh_bladder_cd4",
                subset="paired_normal_patients_non_treg_ge2",
                contrast="tumor_exclusive_vs_normal_shared",
            )

        sel = clones["selected_for_review"].to_numpy(dtype=bool)
        dom = clones["dom_subtype"].fillna("unmapped")
        shortlist_rows.append(
            {
                "cohort": "oh_bladder_cd4",
                "run": rname,
                "retained_cells": int(cfg["cell_filters"][-1]["after"]),
                "total_clones": len(clones),
                "expanded_ge2": int((clones["cell_count"] >= 2).sum()),
                "shortlisted": int(sel.sum()),
                "eff_cxcl13_cyto_shortlisted": int(
                    (sel & dom.isin(["CD4CXCL13", "CD4GZMB", "CD4GZMK"]).to_numpy()).sum()
                ),
                "cxcl13_shortlisted": int((sel & dom.eq("CD4CXCL13").to_numpy()).sum()),
                "cyto_gzmb_gzmk_shortlisted": int((sel & dom.isin(["CD4GZMB", "CD4GZMK"]).to_numpy()).sum()),
                "treg_il2rahi_shortlisted": int((sel & dom.eq("CD4IL2RAHI").to_numpy()).sum()),
                "cm_th17_shortlisted": int((sel & dom.isin(["CD4CM", "CD4TH17"]).to_numpy()).sum()),
                "normal_shared_shortlisted": int((sel & clones["in_matched_normal"].to_numpy(dtype=bool)).sum()),
            }
        )
        facts[rname] = shortlist_rows[-1]

    return auroc_rows, shortlist_rows, {"runs": facts}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional tcrsift cache directory containing vdjdb",
    )
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    from tcrsift import __version__
    from tcrsift.datacache import ensure_cached

    vdjdb_path = ensure_cached("vdjdb", data_dir=args.cache_dir, auto_download=False)
    vdjdb_release = "none"
    if vdjdb_path is not None:
        subs = sorted(p.name for p in Path(vdjdb_path).glob("vdjdb-*") if p.is_dir())
        vdjdb_release = subs[-1] if subs else Path(vdjdb_path).name

    t0 = time.perf_counter()
    print("1/4 Running Pulliam 2024 paired Blood + Tumor MCC benchmark (GSE227708)...")
    pul_aurocs, pul_shortlists, pul_facts = run_pulliam_blood(args.data_dir, args.work_dir, args.cache_dir)

    print("2/4 Running Chiffelle/Barras 2024 pre-ACT T0 vs post-ACT T30 benchmark (GSE222448)...")
    chi_aurocs, chi_shortlists, chi_facts = run_chiffelle_act(args.data_dir, args.work_dir, args.cache_dir)

    print("3/4 Running Veatch 2022 melanoma neoantigen CD4+ TIL benchmark (GSE198265)...")
    vea_aurocs, vea_shortlists, vea_facts = run_veatch_cd4(args.data_dir)

    print("4/4 Running Oh 2020 whole-transcriptome bladder CD4+ TIL benchmark (GSE149652)...")
    oh_aurocs, oh_shortlists, oh_facts = run_oh_cd4(args.data_dir, args.work_dir, args.cache_dir)

    pd.DataFrame(pul_aurocs + chi_aurocs + vea_aurocs + oh_aurocs).to_csv(
        results_dir / "auroc_blood_cd4.csv", index=False
    )
    pd.DataFrame(pul_shortlists).to_csv(results_dir / "shortlist_pulliam_blood.csv", index=False)
    pd.DataFrame(chi_shortlists).to_csv(results_dir / "shortlist_chiffelle_act.csv", index=False)
    pd.DataFrame(vea_shortlists).to_csv(results_dir / "shortlist_veatch_cd4.csv", index=False)
    pd.DataFrame(oh_shortlists).to_csv(results_dir / "shortlist_oh_cd4.csv", index=False)

    summary = {
        "tcrsift_version": __version__,
        "vdjdb_release": vdjdb_release,
        "elapsed_seconds": round(time.perf_counter() - t0, 1),
        "ground_truth_provenance": {
            "pulliam_blood": {
                "study": "Pulliam et al., Cell Reports Medicine 2024 (PMID:38340723)",
                "accession": "GSE227708 + Table S4 (mmc5.xlsx)",
                "ground_truth": "DNA-barcoded MCPyV T-antigen multimer+ antitumor CD8+ T cells vs CMV/EBV/Influenza viral control CD8+ T cells across 10 demultiplexed peripheral blood libraries + 2 matched tumor libraries (>=20 cells across 5 hashed 10x runs) in 4 MCC patients",
            },
            "chiffelle_act": {
                "study": "Chiffelle et al., Immunity 2024 / Barras et al., Science Immunology 2024 (PMID:39321802)",
                "accession": "GSE222448 + Table S3 (mmc3.xlsx)",
                "ground_truth": "Autologous melanoma co-culture verified tumor-reactive (TR) vs non-reactive (NTR) TCRs across paired baseline T0 (8 samples) and Day-30 post-ACT T30 (6 samples) biopsies in Patients 1, 5, 7, 8, 10, and 13",
            },
            "veatch_cd4": {
                "study": "Veatch et al., Cancer Cell 2022 (PMID:35413271)",
                "accession": "GSE198265 (GSM5942454-GSM5942471)",
                "ground_truth": "Autologous neoantigen-specific (CD4_neo: BRAF V600E, ERAP1, SNRPB, HTATSF1, ITGB5, B7/C8/F3) and melanoma self-antigen (CD4_self) CD4+ TIL clonotypes across 4 melanoma patients (397-gene targeted BD Rhapsody mRNA + AbSeq)",
            },
            "oh_cd4": {
                "study": "Oh et al., Cell 2020 (PMID:32497499)",
                "accession": "GSE149652 (GSM4506643-GSM4506662)",
                "ground_truth": "Whole-transcriptome 10x scRNA-seq + paired alpha/beta TCR across 7 human bladder tumors and 6 matched normal bladder CD4+ T-cell libraries with author-curated CD4CXCL13, cytotoxic CD4GZMB/CD4GZMK, FOXP3+ CD4IL2RAHI Treg, and CD4CM states",
            },
        },
        "pulliam_blood": pul_facts,
        "chiffelle_act": chi_facts,
        "veatch_cd4": vea_facts,
        "oh_cd4": oh_facts,
    }
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Done in {summary['elapsed_seconds']}s -> {results_dir}")


if __name__ == "__main__":
    main()
