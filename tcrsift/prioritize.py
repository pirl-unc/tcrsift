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

"""Explicit context presets and round-robin clone selection.

Presets are starting points, not tissue classifiers or specificity predictions.
See docs/user-guide/til-signatures.md for their rationale and limitations.
"""

from __future__ import annotations

from collections import deque

import numpy as np
import pandas as pd

from .signature_methods import NEOANTIGEN_SIGNATURES, SIGNATURES

REGISTRY = {**SIGNATURES, **NEOANTIGEN_SIGNATURES}
CORE_SIGNATURES = (
    "Differentiated", "AntigenExperienced", "Cytolytic",
    "AcuteActivation", "Proliferation", "AIM",
)
CONTEXT_SIGNATURES = {
    "generic": CORE_SIGNATURES,
    "blood": (*CORE_SIGNATURES, "CirculatingMemory"),
    "blood-tumor": (*CORE_SIGNATURES, "CirculatingMemory", "NeoTCR_PBL"),
    "solid-tumor": (*CORE_SIGNATURES, "TumorReactive", "MANAscore", "NeoTCR8", "NeoTCR4"),
    # Avoid the epithelial-residency component of TumorReactive in effusions.
    "mpe": (*CORE_SIGNATURES, "CirculatingMemory", "MANAscore", "NeoTCR8", "NeoTCR4"),
    # No leukemia-specific predictor is shipped. Use broad state coverage.
    "heme": (*CORE_SIGNATURES, "CirculatingMemory"),
}
SIGNATURE_LINEAGE = {
    "NeoTCR4": "cd4", "NeoTCR8": "cd8", "NeoTCR_PBL": "cd8", "MANAscore": "cd8",
}
QC_DEFAULTS = {
    "min_genes": 250, "max_genes": 15000,
    "min_counts": 500, "max_counts": 100000,
    "min_mito_pct": 0.0, "max_mito_pct": 8.0,
}


def resolve_signatures(context, tcell_type, signatures=None, exclude_signatures=None):
    """Resolve names once; reject incompatible explicit lineage selections."""
    lookup = {name.lower(): name for name in REGISTRY}

    def canonical(names):
        result = []
        for name in names:
            if name.lower() not in lookup:
                raise ValueError(f"Unknown signature {name!r}; choose from {', '.join(REGISTRY)}")
            result.append(lookup[name.lower()])
        return list(dict.fromkeys(result))

    names = canonical(signatures if signatures is not None else CONTEXT_SIGNATURES[context])
    excluded = set(canonical(exclude_signatures or []))
    names = [name for name in names if name not in excluded]
    incompatible = [n for n in names if tcell_type != "both"
                    and SIGNATURE_LINEAGE.get(n, tcell_type) != tcell_type]
    if signatures is not None and incompatible:
        raise ValueError(f"Signatures {incompatible} are incompatible with --tcell-type {tcell_type}")
    names = [name for name in names if name not in incompatible]
    if not names:
        raise ValueError("No signatures remain after context, lineage, and exclusion settings")
    return names


def parse_expression_limits(values):
    """Parse repeatable GENE=VALUE thresholds without evaluating expressions."""
    limits = []
    for value in values or []:
        try:
            gene, number = value.split("=", 1)
            threshold = float(number)
            if not gene.strip() or not np.isfinite(threshold) or threshold < 0:
                raise ValueError
        except ValueError as exc:
            raise ValueError(f"Expected GENE=VALUE with a finite nonnegative value; got {value!r}") from exc
        limits.append((gene.strip().upper(), threshold))
    return limits


def filter_cells(adata, args):
    """Apply one shared cell population to abundance, frequencies, and scores."""
    from .signature_methods import expression_frame_from_adata

    obs = adata.obs
    masks = {}
    masks["lineage"] = (obs.is_CD8 | obs.is_CD4 if args.tcell_type == "both"
                        else obs[f"is_{args.tcell_type.upper()}"])
    if args.min_cd3 > 0:
        if "CD3" not in obs:
            raise ValueError("CD3 expression is required for --min-cd3; set it to 0 to disable")
        masks["cd3"] = pd.to_numeric(obs.CD3, errors="coerce").ge(args.min_cd3)
    for chain in ("alpha", "beta"):
        col = f"CDR3_{chain}"
        if col not in obs:
            raise ValueError(f"Paired VDJ data are required: missing {col}")
        masks[f"{chain}_cdr3"] = obs[col].astype("string").str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+", na=False)
    for metric, threshold in (("umis", args.min_vdj_umis), ("reads", args.min_vdj_reads)):
        if threshold == 0:
            continue
        for chain in ("TRA", "TRB"):
            col = f"{chain}_1_{metric}"
            if col not in obs:
                raise ValueError(f"Missing {col}, required for --min-vdj-{metric}")
            masks[col] = pd.to_numeric(obs[col], errors="coerce").ge(threshold)
    limits = [("min", *item) for item in parse_expression_limits(args.min_expression)]
    limits += [("max", *item) for item in parse_expression_limits(args.max_expression)]
    if limits:
        # CP10K uses all measured genes, not only the requested marker panel.
        expr = expression_frame_from_adata(adata, [gene for _, gene, _ in limits])
        totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
        normalized = np.log1p(expr.to_numpy() * (10000 / np.maximum(totals, 1))[:, None])
        expr = pd.DataFrame(normalized, columns=expr.columns)
        for direction, gene, threshold in limits:
            values = expr[gene].to_numpy()
            masks[f"{direction}_expression:{gene}={threshold}"] = (
                values >= threshold if direction == "min" else values <= threshold
            )
    keep = np.ones(adata.n_obs, dtype=bool)
    counts = []
    for label, mask in masks.items():
        before = int(keep.sum())
        keep &= np.asarray(mask, dtype=bool)
        counts.append({"filter": label, "before": before, "after": int(keep.sum())})
    counts_by_sample = (
        obs.assign(retained=keep).groupby("sample", observed=True).retained
        .agg(["size", "sum"]).rename(columns={"size": "after_gex_qc", "sum": "retained"})
        .reset_index().to_dict("records")
    )
    if not keep.any():
        raise ValueError("No cells remain after lineage, paired-chain, and expression filters")
    result = adata[keep].copy()
    result.obs["lineage"] = np.where(result.obs.is_CD8, "cd8", "cd4")
    return result, counts, counts_by_sample


def gene_filter_mask(clones, include, exclude, segment):
    """Match exact, allele-insensitive V/J calls in either chain."""
    def base(value):
        return str(value).strip().upper().split("*", 1)[0]

    include, exclude = set(map(base, include or [])), set(map(base, exclude or []))
    allowed = pd.Series(True, index=clones.index)
    if not include and not exclude:
        return allowed
    columns = [f"{chain}_{segment}_gene" for chain in ("alpha", "beta")]
    if any(col not in clones for col in columns):
        raise ValueError(f"Both chain {segment.upper()}-gene calls are required for gene filters")
    genes = clones[columns].apply(lambda col: col.fillna("").map(base))
    # Missing calls cannot establish that a clone passes an exclusion filter.
    allowed &= genes.ne("").all(axis=1)
    if include:
        allowed &= genes.isin(include).any(axis=1)
    if exclude:
        allowed &= ~genes.isin(exclude).any(axis=1)
    return allowed


def signature_pass_mask(scores, signature, quantile=0.9, min_score=0.0):
    """Qualify one list's informative head before clone exclusions or budgeting.

    The percentile is relative to all scored clones in the same stratum.
    A positive score and a non-flat list avoid filling with zero/constant
    programs. These are transparent selection heuristics, not significance tests.
    """
    col = f"signature_{signature}"
    values = scores[col].where(np.isfinite(scores[col]))
    groups = [scores[key] for key in ("donor", "sample", "lineage")]
    grouped = values.groupby(groups, observed=True)
    low, high = grouped.transform("min"), grouped.transform("max")
    varying = ~np.isclose(low, high, rtol=1e-9, atol=1e-12)
    return (values.gt(min_score) & scores[f"{col}_percentile"].ge(quantile) & varying).fillna(False)


def select_round_robin(clones, scores, signatures, max_clones=200, quantile=0.9, min_score=0.0):
    """Take one unseen clone per (signature, sample, lineage) list each round.

    One budget covers the entire run. Clone identity includes the donor;
    shared sequences in different patients remain separate candidate rows.
    Only the informative head of each list participates. Exhausted lists
    contribute no more turns; other lists continue without quota limits.
    """
    result = clones.copy()
    result["selection_rank"] = pd.Series(pd.NA, index=result.index, dtype="Int64")
    for col in ("selected_signature", "selected_sample", "selected_lineage"):
        result[col] = ""
    eligible = result[result.eligible_for_review]
    by_clone = dict(zip(zip(eligible.donor, eligible.CDR3ab), eligible.index))
    lists = []
    for signature in signatures:
        col = f"signature_{signature}"
        passing = scores[signature_pass_mask(scores, signature, quantile, min_score)]
        for (donor, sample, lineage), group in passing.groupby(
            ["donor", "sample", "lineage"], sort=True, observed=True,
        ):
            ranked = group.sort_values(
                [col, "frequency", "cells", "CDR3ab"], ascending=[False, False, False, True],
            )
            queue = deque((donor, clone) for clone in ranked.CDR3ab if (donor, clone) in by_clone)
            if queue:
                lists.append((signature, sample, lineage, queue))
    selected = set()
    limit = max_clones if max_clones else len(eligible)
    while len(selected) < limit:
        before = len(selected)
        for signature, sample, lineage, queue in lists:
            while queue and queue[0] in selected:
                queue.popleft()
            if not queue:
                continue
            clone = queue.popleft()
            selected.add(clone)
            idx = by_clone[clone]
            result.loc[idx, ["selection_rank", "selected_signature", "selected_sample", "selected_lineage"]] = [
                len(selected), signature, sample, lineage,
            ]
            if len(selected) == limit:
                break
        if len(selected) == before:
            break
    result["selected_for_review"] = result.selection_rank.notna()
    return result.sort_values(["selection_rank", "donor", "CDR3ab"], na_position="last")
