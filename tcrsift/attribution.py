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


def attribute_cells(
    df: pd.DataFrame,
    config,
    group_by: str = "CDR3ab",
    min_umi: int = 2,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Distribute each cell's weight across complete paired clones.

    Parameters
    ----------
    df : pd.DataFrame
        Cell table (``adata.obs``), indexed by cell barcode. Uses
        ``CDR3_alpha``/``CDR3_beta`` (primary chains), ``TRA_2_cdr3`` /
        ``TRB_2_cdr3`` (secondary), and the matching ``*_umis`` columns.
    config : AttributionConfig
        Attribution settings.
    group_by : str
        Only ``"CDR3ab"`` (paired) is supported; other values return empty.
    min_umi : int
        Minimum per-chain UMI for a chain to count as present.

    Returns
    -------
    (long_table, merge_map)
        ``long_table`` columns: ``cell_barcode, CDR3ab, sample, weight, kind``;
        one row per (cell, attributed clone), weights per cell summing to 1.0
        (attributed) or 0.0 (dropped). ``merge_map`` maps each merged
        constituent clone ID to its canonical (surviving) clone ID.
    """
    if group_by != "CDR3ab":
        return pd.DataFrame(columns=LONG_COLUMNS), {}

    n = len(df)
    bc = np.asarray(df.index.astype(str))
    sample = (
        np.asarray(df["sample"].astype(str)) if "sample" in df.columns
        else np.array([""] * n)
    )

    def _str(col: str) -> pd.Series:
        if col in df.columns:
            return df[col].astype("string")
        return pd.Series(pd.array([pd.NA] * n, dtype="string"), index=df.index)

    def _cdr3_valid(s: pd.Series) -> np.ndarray:
        return (
            s.notna() & (s.str.len() > 0) & (s.str.lower() != "nan")
        ).to_numpy()

    def _chain_valid(cdr3_col: str, umi_col: str) -> np.ndarray:
        v = _cdr3_valid(_str(cdr3_col))
        if umi_col in df.columns:
            v = v & (df[umi_col].fillna(0).to_numpy(dtype=float) >= min_umi)
        return v

    A1 = _str("CDR3_alpha").fillna("").to_numpy()
    A2 = _str("TRA_2_cdr3").fillna("").to_numpy()
    B1 = _str("CDR3_beta").fillna("").to_numpy()
    B2 = _str("TRB_2_cdr3").fillna("").to_numpy()
    va1 = _chain_valid("CDR3_alpha", "TRA_1_umis")
    va2 = _chain_valid("TRA_2_cdr3", "TRA_2_umis")
    vb1 = _chain_valid("CDR3_beta", "TRB_1_umis")
    vb2 = _chain_valid("TRB_2_cdr3", "TRB_2_umis")

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
    if config.dual_alpha_merge and dual_alpha.any():
        da_idx = np.where(dual_alpha)[0]
        triples = []
        for i in da_idx:
            a_lo, a_hi = sorted((A1[i], A2[i]))
            triples.append((a_lo, a_hi, B1[i]))
        triple_counts = pd.Series(triples).value_counts()
        for (a_lo, a_hi, b), count in triple_counts.items():
            if count < config.dual_alpha_min_cells:
                continue
            k_lo, k_hi = _make_key(a_lo, b), _make_key(a_hi, b)
            # Canonical = higher-abundance constituent (lexical tie-break) so
            # the merged ID stays a real observed clonotype string.
            p_lo, p_hi = prior_of.get(k_lo, 0), prior_of.get(k_hi, 0)
            if (p_lo, k_hi) >= (p_hi, k_lo):
                canon = k_lo
            else:
                canon = k_hi
            merge_map[k_lo] = canon
            merge_map[k_hi] = canon

    def _canon(key: str) -> str:
        return merge_map.get(key, key)

    # Canonical-clone priors (merged constituents pooled).
    canon_prior: dict[str, float] = {}
    for key, p in prior_of.items():
        canon_prior[_canon(key)] = canon_prior.get(_canon(key), 0.0) + float(p)

    uniform = config.share_weighting == "uniform"

    def _eff_prior(key: str) -> float:
        return 1.0 if uniform else float(canon_prior.get(key, 0.0))

    rows_bc: list = []
    rows_key: list = []
    rows_sample: list = []
    rows_w: list = []
    rows_kind: list = []

    def _emit(idx, keys, weights, kind):
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

    def _split_fixed(mask, target_fns, kind_split, kind_primary):
        """Split each masked cell across fixed candidate clones by prior.

        ``target_fns`` is a list of callables i -> clone-key. With
        ``doublet_split`` off, all weight goes to the first (primary) target.
        """
        idx = np.where(mask)[0]
        if not len(idx):
            return
        if not config.doublet_split:
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
            # Proportional where any prior is positive; uniform fallback when
            # every candidate has zero prior.
            w = np.where(totals > 0, target_pr[j] / np.where(totals > 0, totals, 1), 1.0 / m)
            _emit(idx, target_keys[j], w, kind_split)

    # --- Category 2: dual-alpha (merge or split) --------------------------
    if dual_alpha.any():
        merged_mask = np.zeros(n, dtype=bool)
        if merge_map:
            da_idx = np.where(dual_alpha)[0]
            for i in da_idx:
                if _make_key(A1[i], B1[i]) in merge_map:
                    merged_mask[i] = True
        # Merged dual-alpha cells: weight 1.0 on the canonical merged clone.
        m_idx = np.where(merged_mask)[0]
        if len(m_idx):
            keys = [_canon(_make_key(A1[i], B1[i])) for i in m_idx]
            _emit(m_idx, keys, np.ones(len(m_idx)), "dual_alpha_merge")
        # Non-merged dual-alpha: split (or primary).
        _split_fixed(
            dual_alpha & (~merged_mask),
            [lambda i: _make_key(A1[i], B1[i]), lambda i: _make_key(A2[i], B1[i])],
            "dual_alpha_split",
            "primary",
        )

    # --- Category 3: dual-beta -> split -----------------------------------
    dual_beta = va1 & vb1 & vb2 & (~va2)
    _split_fixed(
        dual_beta,
        [lambda i: _make_key(A1[i], B1[i]), lambda i: _make_key(A1[i], B2[i])],
        "dual_beta_split",
        "primary",
    )

    # --- Category 4: quad (two alphas, two betas) -> split ----------------
    quad = va1 & va2 & vb1 & vb2
    _split_fixed(
        quad,
        [
            lambda i: _make_key(A1[i], B1[i]),
            lambda i: _make_key(A2[i], B1[i]),
            lambda i: _make_key(A1[i], B2[i]),
            lambda i: _make_key(A2[i], B2[i]),
        ],
        "quad_split",
        "primary",
    )

    # --- Category 5: alpha-dropout -> beta-share --------------------------
    n_dropped_beta = 0
    alpha_dropout = (~va1) & (~va2) & (vb1 | vb2)
    if alpha_dropout.any() and config.beta_sharing:
        # Clone-by-beta lookup: each canonical complete clone and its beta.
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
        "cell_barcode": rows_bc,
        "CDR3ab": rows_key,
        "sample": rows_sample,
        "weight": rows_w,
        "kind": rows_kind,
    })
    if len(long_table):
        # Collapse any within-cell duplicate clone targets (e.g. two split
        # candidates that canonicalize to the same merged clone) by summing.
        long_table = (
            long_table.groupby(
                ["cell_barcode", "CDR3ab", "sample", "kind"], as_index=False, observed=True
            )["weight"].sum()
        )
        # Drop zero-weight attributions — a proportional split sends 0 to a
        # zero-prior candidate (e.g. a doublet's partner clone that was never
        # observed complete), which would otherwise surface as a phantom clone
        # with cell_count 0.0.
        long_table = long_table[long_table["weight"] > 0]
        long_table = long_table[LONG_COLUMNS].reset_index(drop=True)

    return long_table, merge_map
