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

import logging
from collections import Counter, deque
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .signature_methods import NEOANTIGEN_SIGNATURES, SIGNATURES
from .species import normalize_species

logger = logging.getLogger(__name__)

REGISTRY = {**SIGNATURES, **NEOANTIGEN_SIGNATURES}
CORE_SIGNATURES = (
    "Differentiated", "AntigenExperienced", "Cytolytic",
    "AcuteActivation", "Proliferation", "AIM",
)
CONTEXT_SIGNATURES = {
    "generic": CORE_SIGNATURES,
    "blood": ("CirculatingMemory", *CORE_SIGNATURES),
    "blood-tumor": ("NeoTCR_PBL", "CirculatingMemory", *CORE_SIGNATURES),
    "solid-tumor": ("TumorReactive", "MANAscore", "NeoTCR8", "NeoTCR4", *CORE_SIGNATURES),
    # Avoid the epithelial-residency component of TumorReactive in effusions.
    "mpe": ("MANAscore", "NeoTCR8", "NeoTCR4", "CirculatingMemory", *CORE_SIGNATURES),
    # No leukemia-specific predictor is shipped. Use broad state coverage.
    "heme": ("CirculatingMemory", *CORE_SIGNATURES),
}
# Human neoantigen-reactivity signatures (NeoTCR4/8/PBL, MANAscore) contain
# primate-restricted genes and human lncRNAs, so default mouse presets omit them.
CONTEXT_SIGNATURES_MOUSE = {
    "generic": CORE_SIGNATURES,
    "blood": ("CirculatingMemory", *CORE_SIGNATURES),
    "blood-tumor": ("CirculatingMemory", *CORE_SIGNATURES),
    "solid-tumor": ("TumorReactive", *CORE_SIGNATURES),
    "mpe": ("CirculatingMemory", *CORE_SIGNATURES),
    "heme": ("CirculatingMemory", *CORE_SIGNATURES),
}
# Dedicated tumor/neoantigen-reactivity programs receive Tier-1 round-robin
# priority under round_robin_mode="tiered" (default when --signatures is unset)
# before broader state panels backfill any remaining budget.
PRIORITY_CONTEXT_SIGNATURES = {
    "generic": (),
    "blood": (),
    "blood-tumor": ("NeoTCR_PBL",),
    "solid-tumor": ("TumorReactive", "MANAscore", "NeoTCR8", "NeoTCR4"),
    "mpe": ("MANAscore", "NeoTCR8", "NeoTCR4"),
    "heme": (),
}
SIGNATURE_LINEAGE = {
    "NeoTCR4": "cd4", "NeoTCR8": "cd8", "NeoTCR_PBL": "cd8", "MANAscore": "cd8",
}
QC_DEFAULTS = {
    "min_genes": 250, "max_genes": 15000,
    "min_counts": 500, "max_counts": 100000,
    "min_mito_pct": 0.0, "max_mito_pct": 8.0,
}


@dataclass
class PrioritizeConfig:
    """Typed configuration for :func:`tcrsift.til_prioritize.run_til_prioritize`."""

    sample_sheet: Path
    output_dir: Path
    context: str = "generic"
    tcell_type: str = "cd8"
    species: str = "auto"
    signatures: list[str] | None = None
    exclude_signatures: list[str] | None = None
    max_clones: int = 200
    min_cells: int = 2
    min_frequency: float = 0.001
    signature_quantile: float = 0.9
    signature_cutoff: str = "background"
    signature_evidence: str | None = None
    signature_high_quantile: float = 0.9
    min_signature_high_cells: int = 5
    min_signature_high_fraction: float = 0.0
    signature_subset_pvalue: float = 0.01
    signature_background_quantile: float = 0.99
    signature_background_draws: int = 2000
    signature_background_seed: int = 0
    min_signature_score: float | None = None
    min_signature_support: int = 1
    vdjdb: Path | None = None
    iedb: Path | None = None
    iedb_epitope: Path | None = None
    cedar: Path | None = None
    cache_dir: Path | None = None
    auto_download: bool = False
    database_match: str = "strict_ab"
    round_robin_mode: str = "tiered"
    exclude_known_viral: bool = False
    exclude_known_mart1: bool = True
    exclude_trav12_2: bool = False
    exclude_public_quantile: float | None = None
    exclude_treg: bool = False
    treg_foxp3_fraction: float = 0.5
    call_multimers: bool = False
    multimer_min_umi: float = 10.0
    multimer_min_nc_ratio: float = 5.0
    multimer_min_concordance: float = 0.5
    multimer_min_specific_cells: int = 2
    min_genes: int = 250
    max_genes: int = 15000
    min_counts: int = 500
    max_counts: int = 100000
    min_mito_pct: float = 0.0
    max_mito_pct: float = 8.0
    min_vdj_umis: int = 2
    min_vdj_reads: int = 0
    min_cd3: int = 3
    min_alpha_cdr3_length: int = 0
    min_beta_cdr3_length: int = 0
    min_expression: list[str] = field(default_factory=list)
    max_expression: list[str] = field(default_factory=list)
    include_v_gene: list[str] = field(default_factory=list)
    exclude_v_gene: list[str] = field(default_factory=list)
    include_j_gene: list[str] = field(default_factory=list)
    exclude_j_gene: list[str] = field(default_factory=list)
    verbose: bool = False

    def __post_init__(self) -> None:
        self.sample_sheet = Path(self.sample_sheet)
        self.output_dir = Path(self.output_dir)
        for attr in ("vdjdb", "iedb", "iedb_epitope", "cedar", "cache_dir"):
            val = getattr(self, attr)
            if isinstance(val, str) and not val.strip():
                setattr(self, attr, None)
            elif val is not None and not isinstance(val, Path):
                setattr(self, attr, Path(val))
        for list_attr in (
            "min_expression",
            "max_expression",
            "include_v_gene",
            "exclude_v_gene",
            "include_j_gene",
            "exclude_j_gene",
        ):
            if getattr(self, list_attr) is None:
                setattr(self, list_attr, [])

    @classmethod
    def from_namespace(cls, args: Any) -> PrioritizeConfig:
        """Build a :class:`PrioritizeConfig` from an ``argparse.Namespace`` or similar object."""
        if isinstance(args, cls):
            return args
        kwargs: dict[str, Any] = {}
        for f in fields(cls):
            if hasattr(args, f.name):
                kwargs[f.name] = getattr(args, f.name)
        return cls(**kwargs)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def resolve_signatures(
    context,
    tcell_type,
    signatures=None,
    exclude_signatures=None,
    *,
    species: str = "human",
):
    """Resolve names once; reject incompatible explicit lineage selections."""
    lookup = {name.lower(): name for name in REGISTRY}

    def canonical(names):
        result = []
        for name in names:
            if name.lower() not in lookup:
                raise ValueError(f"Unknown signature {name!r}; choose from {', '.join(REGISTRY)}")
            result.append(lookup[name.lower()])
        return list(dict.fromkeys(result))

    ctx_map = (
        CONTEXT_SIGNATURES_MOUSE
        if normalize_species(species) == "mouse"
        else CONTEXT_SIGNATURES
    )
    names = canonical(signatures if signatures is not None else ctx_map[context])
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
        # Lineage selection is expected to remove many cells; any other single
        # gate taking more than half usually means a threshold does not suit
        # this dataset (a raw-count floor on a shallow library, for example).
        removed = before - int(keep.sum())
        if label != "lineage" and before and removed > before / 2:
            logger.warning(
                "Cell filter '%s' removed %s of %s remaining cells (%.0f%%). Check that "
                "its threshold suits this dataset before relying on the shortlist.",
                label, f"{removed:,}", f"{before:,}", 100 * removed / before,
            )
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


def signature_pass_mask(scores, signature, quantile=0.9, min_score=None, cutoff="background",
                        evidence="mean", min_high_cells=5, min_high_fraction=0.0, max_subset_p=0.01):
    """Qualify mean or high-cell subset evidence within the same stratum."""
    if cutoff not in ("background", "legacy"):
        raise ValueError(f"Unknown signature cutoff {cutoff!r}")
    if evidence not in ("mean", "subset", "either"):
        raise ValueError(f"Unknown signature evidence {evidence!r}")
    if cutoff == "legacy" and evidence != "mean":
        raise ValueError("Subset evidence requires --signature-cutoff background")
    if evidence == "either":
        return signature_pass_mask(scores, signature, quantile, min_score, cutoff) | signature_pass_mask(
            scores, signature, quantile, min_score, cutoff, "subset",
            min_high_cells, min_high_fraction, max_subset_p,
        )
    col = f"signature_{signature}"
    values = scores[col].where(np.isfinite(scores[col]))
    if evidence == "subset":
        tail = scores[f"{col}_subset_tail_probability"]
        passing = (tail.between(0, max_subset_p) & scores[f"{col}_high_cells"].ge(min_high_cells)
                   & scores[f"{col}_high_fraction"].ge(min_high_fraction)
                   & scores[f"{col}_subset_percentile"].ge(quantile) & values.notna())
        if min_score is not None:
            passing &= values.gt(min_score)
        return passing.fillna(False)
    if cutoff == "background":
        floor = scores[f"{col}_noise_floor"]
        above_floor = values.gt(floor) & ~np.isclose(values, floor, rtol=1e-9, atol=1e-12)
    else:  # legacy; other values were rejected above
        above_floor = values.gt(0 if min_score is None else min_score)
    if min_score is not None:
        above_floor &= values.gt(min_score)
    groups = [scores[key] for key in ("donor", "sample", "lineage")]
    grouped = values.groupby(groups, observed=True)
    low, high = grouped.transform("min"), grouped.transform("max")
    varying = ~np.isclose(low, high, rtol=1e-9, atol=1e-12)
    return (above_floor & scores[f"{col}_percentile"].ge(quantile) & varying).fillna(False)


def select_round_robin(clones, scores, signatures, max_clones=200, quantile=0.9,
                       min_score=None, cutoff="background", evidence="mean",
                       min_high_cells=5, min_high_fraction=0.0, max_subset_p=0.01,
                       *, priority_signatures=None):
    """Take one unseen clone per (signature, sample, lineage) list each round.

    One budget covers the entire run. Clone identity includes the donor;
    shared sequences in different patients remain separate candidate rows.
    Only the qualifying head of each list participates. Exhausted lists
    contribute no more turns; other lists continue without quota limits.
    When ``priority_signatures`` is non-empty, those dedicated context
    signatures exhaust their qualifying lists first (Tier 1) before the
    remaining signatures backfill any leftover budget (Tier 2).
    When the budget ends inside a round, each remaining turn goes to the
    signature, then the sample, with the fewest turns in that round.
    """
    scores = scores.reset_index(drop=True)
    result = clones.copy()
    result["selection_rank"] = pd.Series(pd.NA, index=result.index, dtype="Int64")
    for col in ("selected_signature", "selected_sample", "selected_lineage", "selected_evidence"):
        result[col] = ""
    for col in ("selected_cells", "selected_mean_score", "selected_high_cells", "selected_high_fraction"):
        result[col] = np.nan
    eligible = result[result.eligible_for_review]
    by_clone = dict(zip(zip(eligible.donor, eligible.CDR3ab), eligible.index))
    lists = []
    for signature in signatures:
        col = f"signature_{signature}"
        for route in (("mean", "subset") if evidence == "either" else (evidence,)):
            passing = scores[signature_pass_mask(
                scores, signature, quantile, min_score, cutoff, route,
                min_high_cells, min_high_fraction, max_subset_p,
            )]
            for (donor, sample, lineage), group in passing.groupby(
                ["donor", "sample", "lineage"], sort=True, observed=True,
            ):
                sort_cols = ([col, "frequency", "cells", "CDR3ab"] if route == "mean" else
                             [f"{col}_subset_tail_probability", f"{col}_high_fraction", "cells", "CDR3ab"])
                ranked = group.sort_values(sort_cols, ascending=[route == "subset", False, False, True])
                queue = deque((donor, row.CDR3ab, row_id) for row_id, row in ranked.iterrows()
                              if (donor, row.CDR3ab) in by_clone)
                if queue:
                    lists.append((signature, sample, lineage, route, queue))
    selected = set()
    limit = max_clones or float("inf")

    def unseen(entry):
        """Drop clones another list already supplied; True if any remain."""
        queue = entry[-1]
        while queue and queue[0][:2] in selected:
            queue.popleft()
        return bool(queue)

    def take(entry):
        signature, sample, lineage, route, queue = entry
        donor, cdr3, row_id = queue.popleft()
        clone = (donor, cdr3)
        selected.add(clone)
        idx = by_clone[clone]
        result.loc[idx, ["selection_rank", "selected_signature", "selected_sample", "selected_lineage", "selected_evidence"]] = [
            len(selected), signature, sample, lineage, route,
        ]
        row = scores.loc[row_id]
        result.loc[idx, ["selected_cells", "selected_mean_score", "selected_high_cells", "selected_high_fraction"]] = [
            row.cells, row[f"signature_{signature}"], row.get(f"signature_{signature}_high_cells", np.nan),
            row.get(f"signature_{signature}_high_fraction", np.nan),
        ]

    def _run_pool(pool_lists):
        while len(selected) < limit:
            before = len(selected)
            live = [entry for entry in pool_lists if unseen(entry)]
            if limit - len(selected) >= len(live):
                for entry in live:
                    if unseen(entry):
                        take(entry)
            else:
                turns, pending = Counter(), live
                while len(selected) < limit:
                    pending = [entry for entry in pending if unseen(entry)]
                    if not pending:
                        break
                    entry = min(pending, key=lambda item: (
                        turns[item[0]], turns[(item[-1][0][0], item[1], item[2])]))
                    pending = [item for item in pending if item is not entry]
                    turns[entry[0]] += 1
                    turns[(entry[-1][0][0], entry[1], entry[2])] += 1
                    take(entry)
            if len(selected) == before:
                break

    priority_set = set(priority_signatures or ())
    pools = (
        ([e for e in lists if e[0] in priority_set], [e for e in lists if e[0] not in priority_set])
        if priority_set
        else (lists,)
    )
    for pool in pools:
        if pool and len(selected) < limit:
            _run_pool(pool)

    result["selected_for_review"] = result.selection_rank.notna()
    for col in ("selected_cells", "selected_high_cells"):
        result[col] = result[col].astype("Int64")
    return result.sort_values(["selection_rank", "donor", "CDR3ab"], na_position="last")
