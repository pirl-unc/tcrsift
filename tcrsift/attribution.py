# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Partial-information clonotype attribution (#176).

Recovers signal that integer clone counting discards — alpha-dropout cells and
droplet doublets — while keeping every output keyed on complete paired
alpha-beta clones. Each cell distributes a total weight of 1.0 (or 0.0 when it
can't be attributed) across complete clones:

1. Complete single-chain cell -> weight 1.0 to its own ``(a1, b1)`` clone.
2. Alpha-dropout cell (valid beta, no valid alpha) -> weight 1.0 spread across
   complete clones sharing that beta, proportional to each clone's complete-cell
   abundance (``beta_sharing``). No matching complete clone -> dropped.
3. Dual-alpha cell (two alphas, one beta):
   - Recurrence-merge: if the sorted ``(a1, a2, b)`` triple recurs in
     ``>= dual_alpha_min_cells`` cells, treat it as one biological clone
     (allelic inclusion) and merge the ``(a1, b)`` / ``(a2, b)`` clonotypes;
     the cell gets weight 1.0 on the merged clone.
   - Else (singleton doublet): split weight across ``(a1, b)`` / ``(a2, b)``
     proportional to abundance (``doublet_split``).
4. Dual-beta cell -> split across ``(a, b1)`` / ``(a, b2)`` proportionally.

Two passes: pass 1 computes complete-clone abundances (priors) from unambiguous
complete cells; pass 2 distributes ambiguous cells using those priors. An
``em_iterations`` hook is reserved; v1 does a single distribution pass.

The public entry point :func:`attribute_cells` returns a per-(cell, clone) long
table of weights and a ``merge_map`` mapping merged constituent clone IDs to
their canonical ID.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

LONG_COLUMNS = ["cell_barcode", "CDR3ab", "sample", "weight", "kind"]


def _make_key(alpha: str, beta: str) -> str:
    """Clone ID string, matching ``aggregate_clonotypes``' ``alpha_beta``."""
    return f"{alpha}_{beta}"


def _resolve_merge_components(
    edges: list[tuple[str, str]],
    prior_of: dict[str, int],
) -> dict[str, str]:
    """Build a transitive ``merge_map`` over connected dual-alpha pairs."""
    parent: dict[str, str] = {}

    def _find(x: str) -> str:
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    for a, b in edges:
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[rb] = ra

    comps: dict[str, list[str]] = {}
    for node in parent:
        comps.setdefault(_find(node), []).append(node)

    merge_map: dict[str, str] = {}
    for members in comps.values():
        # Pick highest complete-cell prior; lexical ascending tie-break.
        canon = min(members, key=lambda k: (-prior_of.get(k, 0), k))
        for m in members:
            merge_map[m] = canon
    return merge_map


def _attribute_single_pool(
    df: pd.DataFrame,
    config,
    min_umi: int,
    cell_indices: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, str], dict[str, dict]]:
    """Run attribution within a single donor pool (preserving global ``cell_indices``)."""
    n = len(df)
    if n == 0:
        return pd.DataFrame(columns=[*LONG_COLUMNS, "cell_idx"]), {}, {}

    bc = np.asarray(df.index.astype(str))
    sample = (
        np.asarray(df["sample"].fillna("").astype(str))
        if "sample" in df.columns
        else np.array([""] * n)
    )

    def _str(col: str) -> pd.Series:
        if col in df.columns:
            return df[col].astype("string")
        return pd.Series(pd.array([pd.NA] * n, dtype="string"), index=df.index)

    def _cdr3_valid(s: pd.Series) -> np.ndarray:
        present = s.fillna("").str.strip()
        lowered = present.str.lower()
        return (
            (present.str.len() > 0) & (~lowered.isin({"nan", "none", "<na>"}))
        ).to_numpy(dtype=bool)

    def _primary_str(primary_col: str, slot_col: str) -> pd.Series:
        if primary_col in df.columns and slot_col in df.columns:
            p = _str(primary_col)
            s = _str(slot_col)
            p_empty = p.fillna("").str.strip().eq("")
            return p.where(~p_empty, s)
        if primary_col in df.columns:
            return _str(primary_col)
        return _str(slot_col)

    def _series_valid(s: pd.Series, umi_col: str) -> np.ndarray:
        v = _cdr3_valid(s)
        if umi_col in df.columns:
            v = v & (
                pd.to_numeric(df[umi_col], errors="coerce").fillna(0).to_numpy(dtype=float)
                >= min_umi
            )
        return v

    s_a1 = _primary_str("CDR3_alpha", "TRA_1_cdr3")
    s_a2 = _str("TRA_2_cdr3")
    s_b1 = _primary_str("CDR3_beta", "TRB_1_cdr3")
    s_b2 = _str("TRB_2_cdr3")
    A1 = s_a1.fillna("").str.strip().to_numpy(dtype=object)
    A2 = s_a2.fillna("").str.strip().to_numpy(dtype=object)
    B1 = s_b1.fillna("").str.strip().to_numpy(dtype=object)
    B2 = s_b2.fillna("").str.strip().to_numpy(dtype=object)
    va1 = _series_valid(s_a1, "TRA_1_umis")
    va2 = _series_valid(s_a2, "TRA_2_umis")
    vb1 = _series_valid(s_b1, "TRB_1_umis")
    vb2 = _series_valid(s_b2, "TRB_2_umis")

    # Collapse homotypic AA contigs within a cell (same AA, distinct NT/contig)
    # so synonymous dual contigs aren't misclassified as dual-alpha/dual-beta.
    va2 = va2 & ~(va1 & (A1 == A2))
    vb2 = vb2 & ~(vb1 & (B1 == B2))

    # Promote a valid secondary slot into an invalid primary slot so cells where
    # slot 1 failed validity/UMI while slot 2 passed are never silently dropped.
    swap_a = (~va1) & va2
    if swap_a.any():
        A1 = np.where(swap_a, A2, A1)
        va1 = va1 | swap_a
        va2 = va2 & (~swap_a)
    swap_b = (~vb1) & vb2
    if swap_b.any():
        B1 = np.where(swap_b, B2, B1)
        vb1 = vb1 | swap_b
        vb2 = vb2 & (~swap_b)

    # --- Pass 1: priors from unambiguous complete-single cells -------------
    complete_single = va1 & vb1 & (~va2) & (~vb2)
    cs_keys = np.array(
        [_make_key(a, b) for a, b in zip(A1[complete_single], B1[complete_single])]
    )
    priors = (
        pd.Series(cs_keys).value_counts() if len(cs_keys)
        else pd.Series(dtype=int)
    )
    prior_of = priors.to_dict()

    # --- Dual-alpha recurrence-merge detection ----------------------------
    dual_alpha = va1 & va2 & vb1 & (~vb2)
    merge_map: dict[str, str] = {}
    merge_stats: dict[str, dict] = {}
    if config.dual_alpha_merge and dual_alpha.any():
        beta_paired = vb1 & (va1 | va2)
        beta_cells: dict[str, int] = {}
        beta_alphas: dict[str, set] = {}
        for i in np.where(beta_paired)[0]:
            b = B1[i]
            beta_cells[b] = beta_cells.get(b, 0) + 1
            partners = beta_alphas.setdefault(b, set())
            if va1[i]:
                partners.add(A1[i])
            if va2[i]:
                partners.add(A2[i])

        min_fraction = getattr(config, "dual_alpha_min_fraction", 0.0)
        max_partners = getattr(config, "dual_alpha_max_beta_partners", 0)

        da_idx = np.where(dual_alpha)[0]
        triples = []
        for i in da_idx:
            a_lo, a_hi = sorted((A1[i], A2[i]))
            triples.append((a_lo, a_hi, B1[i]))
        triple_counts = pd.Series(triples).value_counts()
        n_rejected = 0
        accepted_edges: list[tuple[str, str]] = []
        accepted_meta: list[tuple[str, str, str, int, int]] = []
        for (a_lo, a_hi, b), count in triple_counts.items():
            if count < config.dual_alpha_min_cells:
                continue
            denom = beta_cells.get(b, count)
            if min_fraction > 0 and denom and (count / denom) < min_fraction:
                n_rejected += 1
                continue
            if max_partners and len(beta_alphas.get(b, ())) > max_partners:
                n_rejected += 1
                continue
            k_lo, k_hi = _make_key(a_lo, b), _make_key(a_hi, b)
            accepted_edges.append((k_lo, k_hi))
            accepted_meta.append((k_lo, k_hi, b, int(count), int(denom)))
        if accepted_edges:
            merge_map = _resolve_merge_components(accepted_edges, prior_of)
            for k_lo, _k_hi, b, count, denom in accepted_meta:
                canon = merge_map[k_lo]
                prev_n = merge_stats.get(canon, {}).get("dual_alpha_n_dual_cells", 0)
                total_n = prev_n + count
                merge_stats[canon] = {
                    "dual_alpha_cofraction": round(total_n / denom, 4) if denom else None,
                    "dual_alpha_beta_partners": len(beta_alphas.get(b, ())),
                    "dual_alpha_n_dual_cells": int(total_n),
                }
        if n_rejected:
            logger.info(
                "  Attribution: rejected %d recurrent dual-alpha triple(s) as "
                "public/promiscuous-beta coincidences (not allelic inclusion, #198)",
                n_rejected,
            )

    def _canon(key: str) -> str:
        return merge_map.get(key, key)

    # Canonical-clone priors (merged constituents pooled).
    canon_prior: dict[str, float] = {}
    for key, p in prior_of.items():
        canon_prior[_canon(key)] = canon_prior.get(_canon(key), 0.0) + float(p)

    uniform = config.share_weighting == "uniform"

    def _eff_prior(key: str) -> float:
        if uniform:
            return 1.0 if key in canon_prior else 0.0
        return float(canon_prior.get(key, 0.0))

    rows_idx: list = []
    rows_bc: list = []
    rows_key: list = []
    rows_sample: list = []
    rows_w: list = []
    rows_kind: list = []

    def _emit(idx, keys, weights, kind):
        rows_idx.extend(cell_indices[idx])
        rows_bc.extend(bc[idx])
        rows_key.extend(keys)
        rows_sample.extend(sample[idx])
        rows_w.extend(weights)
        rows_kind.extend([kind] * len(idx))

    # --- Category 1: complete single -> weight 1.0 ------------------------
    cs_idx = np.where(complete_single)[0]
    if len(cs_idx):
        keys = [_canon(_make_key(A1[i], B1[i])) for i in cs_idx]
        _emit(cs_idx, keys, np.ones(len(cs_idx)), "complete")

    def _split_fixed(mask, target_fns, kind_split, kind_primary, do_split):
        idx = np.where(mask)[0]
        if not len(idx):
            return
        if not do_split:
            keys = [_canon(target_fns[0](i)) for i in idx]
            _emit(idx, keys, np.ones(len(idx)), kind_primary)
            return
        target_keys = [[_canon(fn(i)) for i in idx] for fn in target_fns]
        target_pr = [
            np.array([_eff_prior(k) for k in tk]) for tk in target_keys
        ]
        totals = np.sum(target_pr, axis=0)
        m = len(target_fns)
        for j in range(m):
            w = np.where(totals > 0, target_pr[j] / np.where(totals > 0, totals, 1), 1.0 / m)
            _emit(idx, target_keys[j], w, kind_split)

    # --- Category 2: dual-alpha (merge or split) --------------------------
    if dual_alpha.any():
        merged_mask = np.zeros(n, dtype=bool)
        if merge_map:
            da_idx = np.where(dual_alpha)[0]
            for i in da_idx:
                k1 = _make_key(A1[i], B1[i])
                k2 = _make_key(A2[i], B1[i])
                if (k1 in merge_map or k2 in merge_map) and _canon(k1) == _canon(k2):
                    merged_mask[i] = True
        m_idx = np.where(merged_mask)[0]
        if len(m_idx):
            keys = [_canon(_make_key(A1[i], B1[i])) for i in m_idx]
            _emit(m_idx, keys, np.ones(len(m_idx)), "dual_alpha_merge")
        _split_fixed(
            dual_alpha & (~merged_mask),
            [lambda i: _make_key(A1[i], B1[i]), lambda i: _make_key(A2[i], B1[i])],
            "dual_alpha_split",
            "primary",
            config.dual_alpha_split,
        )

    # --- Category 3: dual-beta -> split -----------------------------------
    dual_beta = va1 & vb1 & vb2 & (~va2)
    _split_fixed(
        dual_beta,
        [lambda i: _make_key(A1[i], B1[i]), lambda i: _make_key(A1[i], B2[i])],
        "dual_beta_split",
        "primary",
        config.dual_beta_split,
    )

    # --- Category 4: quad (two alphas, two betas) -> split ----------------
    quad = va1 & va2 & vb1 & vb2
    if quad.any():
        a_list = [A1] + ([A2] if config.dual_alpha_split else [])
        b_list = [B1] + ([B2] if config.dual_beta_split else [])
        quad_fns = [
            (lambda i, aa=aa, bb=bb: _make_key(aa[i], bb[i]))
            for aa in a_list
            for bb in b_list
        ]
        _split_fixed(quad, quad_fns, "quad_split", "primary", len(quad_fns) > 1)

    # --- Category 5: alpha-dropout -> beta-share --------------------------
    n_dropped_beta = 0
    alpha_dropout = (~va1) & (~va2) & (vb1 | vb2)
    if alpha_dropout.any() and config.beta_sharing:
        clone_keys = list(canon_prior.keys())
        if clone_keys:
            clone_beta = pd.DataFrame({
                "clone": clone_keys,
                "beta": [k.split("_", 1)[1] if "_" in k else "" for k in clone_keys],
                "prior": [_eff_prior(k) for k in clone_keys],
            })
            ad_idx = np.where(alpha_dropout)[0]
            cell_betas = []
            for i in ad_idx:
                betas = set()
                if vb1[i]:
                    betas.add(B1[i])
                if vb2[i]:
                    betas.add(B2[i])
                for b in betas:
                    cell_betas.append((i, b))
            cb = pd.DataFrame(cell_betas, columns=["i", "beta"])
            cand = cb.merge(clone_beta, on="beta", how="inner")
            attributed = set(cand["i"].unique())
            n_dropped_beta = len(ad_idx) - len(attributed)
            if len(cand):
                tot = cand.groupby("i")["prior"].transform("sum")
                cnt = cand.groupby("i")["prior"].transform("size")
                cand["weight"] = np.where(tot > 0, cand["prior"] / tot.where(tot > 0, 1), 1.0 / cnt)
                _emit(
                    cand["i"].to_numpy(),
                    cand["clone"].tolist(),
                    cand["weight"].to_numpy(),
                    "beta_share",
                )
        else:
            n_dropped_beta = int(alpha_dropout.sum())
    elif alpha_dropout.any():
        n_dropped_beta = int(alpha_dropout.sum())

    if n_dropped_beta:
        logger.info(
            "  Attribution: %d alpha-dropout cells had no matching complete "
            "clone (dropped)", n_dropped_beta,
        )

    long_table = pd.DataFrame({
        "cell_idx": rows_idx,
        "cell_barcode": rows_bc,
        "CDR3ab": rows_key,
        "sample": rows_sample,
        "weight": rows_w,
        "kind": rows_kind,
    })
    if len(long_table):
        # Group by unique positional `cell_idx` (NOT `cell_barcode` alone, since
        # 10x barcodes legitimately repeat across concatenated samples).
        long_table = (
            long_table.groupby(
                ["cell_idx", "cell_barcode", "CDR3ab", "sample", "kind"],
                as_index=False,
                observed=True,
            )["weight"].sum()
        )
        long_table = long_table[long_table["weight"] > 0].reset_index(drop=True)

    return long_table, merge_map, merge_stats


def attribute_cells(
    df: pd.DataFrame,
    config,
    group_by: str = "CDR3ab",
    min_umi: int = 2,
    *,
    include_cell_idx: bool = False,
) -> tuple[pd.DataFrame, dict[str, str], dict[str, dict]]:
    """Distribute each cell's weight across complete paired clones."""
    out_cols = [*LONG_COLUMNS, "cell_idx"] if include_cell_idx else list(LONG_COLUMNS)
    if group_by != "CDR3ab" or len(df) == 0:
        return pd.DataFrame(columns=out_cols), {}, {}

    all_indices = np.arange(len(df), dtype=np.int64)
    if "patient_id" in df.columns:
        pid = df["patient_id"]
        bucket_keys = pd.Series(
            np.where(pid.isna(), "__missing_patient_id__", pid.astype(str)),
            index=df.index,
        )
        if bucket_keys.nunique() > 1:
            tables: list[pd.DataFrame] = []
            merged_all: dict[str, str] = {}
            stats_all: dict[str, dict] = {}
            for _bk, pos in pd.Series(all_indices).groupby(bucket_keys.to_numpy(), sort=False):
                pos_arr = pos.to_numpy(dtype=np.int64)
                sub_lt, sub_mm, sub_ms = _attribute_single_pool(
                    df.iloc[pos_arr], config, min_umi, pos_arr
                )
                if len(sub_lt):
                    tables.append(sub_lt)
                merged_all.update(sub_mm)
                stats_all.update(sub_ms)
            long_table = (
                pd.concat(tables, ignore_index=True)
                if tables
                else pd.DataFrame(columns=[*LONG_COLUMNS, "cell_idx"])
            )
            return long_table[out_cols].reset_index(drop=True), merged_all, stats_all

    long_table, merge_map, merge_stats = _attribute_single_pool(
        df, config, min_umi, all_indices
    )
    return long_table[out_cols].reset_index(drop=True), merge_map, merge_stats
