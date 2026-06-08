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

"""Composable in-silico filter layer (#149).

A **stackable filter layer** that composes *on top of* the assay/enrichment
groups (the sorts) and the selection rules. Each layer is a per-clone
percentile predicate on an existing feature — Ppost (#143) or a GEX
signature score (#141/#142/#144) — so a group can be progressively refined:

    CTYneg
    CTYneg · Ppost^low
    CTYneg · Ppost^low · GEXnaive^low
    CTYneg · Ppost^low · GEXnaive^low · GEXresponse^high

Applying the full stack yields one **named composite** subgroup (default
label ``IS`` for *in-silico*), so each assay group gains a filtered twin:
``CTYneg`` and ``CTYneg+IS`` — N assay groups → 2N groups. These twins feed
the cross-group Jaccard-overlap analysis (:func:`insilico_overlap_long`) so
filtered subsets are directly comparable across assays.

Scoring conventions (from the pilot):

- Each predicate is a **lower-is-better percentile rank** in ``[0, 1]``
  (0 = best). ``direction="high"`` flips the rank so high scores are best.
- Tie handling keeps ``Ppost == 0`` a valid *best* value and groups all
  equal scores together (pandas ``rank`` with ``method="average"`` /
  ``"min"``) — never dropped.
- A clone **passes** a predicate when its percentile rank ≤ the predicate's
  ``percentile`` (the best fraction). Missing scores never pass.
- Ranks are computed **within each assay group** when a group column is
  given, so the filter refines each group against its own distribution.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_LABEL = "IS"
DEFAULT_GROUP_SEP = "+"


@dataclass(frozen=True)
class FilterPredicate:
    """One in-silico filter dimension.

    ``score`` is a per-clone column; ``direction`` is ``"low"`` (low score
    is better, e.g. ``Ppost^low``) or ``"high"`` (high is better, e.g.
    ``GEXresponse^high``); ``percentile`` is the best fraction kept
    (``0.5`` → best half). ``rank_method`` is the pandas tie method.
    """

    score: str
    direction: str = "low"
    percentile: float = 0.5
    rank_method: str = "average"

    def __post_init__(self):
        if self.direction not in ("low", "high"):
            raise ValueError(
                f"FilterPredicate.direction must be 'low'/'high', got "
                f"{self.direction!r}"
            )
        if not 0.0 < self.percentile <= 1.0:
            raise ValueError(
                f"FilterPredicate.percentile must be in (0, 1], got "
                f"{self.percentile!r}"
            )

    @property
    def lower_is_better(self) -> bool:
        return self.direction == "low"


def predicates_from_config(cfg: dict) -> list[FilterPredicate]:
    """Parse the ``insilico_filter.predicates`` list into predicates.

    Each entry is a dict ``{score, direction?, percentile?, rank_method?}``.
    """
    out: list[FilterPredicate] = []
    for entry in cfg.get("predicates", []) or []:
        if "score" not in entry:
            raise ValueError(f"insilico_filter predicate missing 'score': {entry!r}")
        out.append(
            FilterPredicate(
                score=entry["score"],
                direction=entry.get("direction", "low"),
                percentile=float(entry.get("percentile", 0.5)),
                rank_method=entry.get("rank_method", "average"),
            )
        )
    return out


# Canonical PRISM criteria (#158) — the SINGLE source of truth. Both the
# FilterPredicate form (annotate_tcrs.PRISM_DEFAULT_PREDICATES) and the dict
# form (selection.PRISM_TERMS) derive from this so they cannot drift apart
# (their drift is exactly how the two PRISM engines diverged): low α/β Ppost,
# high antigen-response, low naive.
PRISM_PREDICATES: list[FilterPredicate] = [
    FilterPredicate("ppost_alpha", "low"),
    FilterPredicate("ppost_beta", "low"),
    FilterPredicate("antigen_response_score", "high"),
    FilterPredicate("naive_score", "low"),
]


def percentile_rank(
    values,
    *,
    lower_is_better: bool = True,
    method: str = "average",
) -> pd.Series:
    """Per-clone percentile rank in ``[0, 1]`` (0 = best). NaN-safe.

    Lower values rank best when ``lower_is_better`` (e.g. Ppost). Ties are
    grouped via pandas ``rank(method=...)`` so equal scores — including a
    cluster of ``Ppost == 0`` — share a rank rather than being dropped or
    arbitrarily ordered. NaN stays NaN (and never passes a predicate).
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce")
    ranks = s.rank(method=method, ascending=lower_is_better, na_option="keep")
    n = int(ranks.notna().sum())
    if n <= 1:
        # 0 (best) for the single/non-NaN value(s); NaN stays NaN.
        return pd.Series(np.where(s.notna(), 0.0, np.nan), index=s.index)
    return (ranks - 1.0) / (n - 1.0)


def _grouped_percentile_rank(
    df: pd.DataFrame, pred: FilterPredicate, group_col: str | None,
) -> pd.Series:
    if pred.score not in df.columns:
        raise ValueError(
            f"insilico_filter: score column {pred.score!r} not in frame "
            f"(have {list(df.columns)})"
        )
    if group_col is not None and group_col in df.columns:
        return df.groupby(group_col, observed=True)[pred.score].transform(
            lambda s: percentile_rank(
                s, lower_is_better=pred.lower_is_better, method=pred.rank_method,
            )
        )
    return percentile_rank(
        df[pred.score], lower_is_better=pred.lower_is_better,
        method=pred.rank_method,
    )


def insilico_mask(
    df: pd.DataFrame,
    predicates: list[FilterPredicate],
    *,
    group_col: str | None = None,
) -> pd.Series:
    """Boolean Series: clones passing **all** predicates (composite members).

    A clone passes a predicate when its percentile rank ≤ the predicate's
    ``percentile``. Clones with a missing score for any predicate are
    excluded (a NaN rank never satisfies ``≤``). When ``group_col`` is set,
    ranks are computed within each group so the filter refines each assay
    group against its own distribution.
    """
    mask = pd.Series(True, index=df.index)
    if not predicates:
        return mask
    for pred in predicates:
        pr = _grouped_percentile_rank(df, pred, group_col)
        mask &= pr <= pred.percentile
    mask = mask.fillna(False)
    return mask


def average_percentile_rank(
    df: pd.DataFrame,
    predicates: list[FilterPredicate],
    *,
    group_col: str | None = None,
    weights: list[float] | None = None,
    require_complete: bool = True,
) -> pd.Series:
    """(Weighted) mean of the per-dimension percentile ranks (0 = best).

    The single row-wise composite engine for in-silico ranking and PRISM. Each
    predicate contributes a lower-is-better percentile rank in [0, 1]; the means
    are averaged (optionally weighted by ``weights``).

    ``require_complete`` (default True) makes a row missing ANY dimension get
    NaN — i.e. don't score a clone you can't fully rank. Set False for a
    NaN-aware mean over the dimensions that ARE present.
    """
    if not predicates:
        return pd.Series(np.nan, index=df.index)
    mat = pd.concat(
        [_grouped_percentile_rank(df, pred, group_col) for pred in predicates],
        axis=1,
    )
    if weights is None:
        return mat.mean(axis=1, skipna=not require_complete)
    w = np.asarray(weights, dtype=float)
    if len(w) != len(predicates):
        raise ValueError("average_percentile_rank: weights length != predicates length")
    if require_complete:
        return (mat * w).sum(axis=1, skipna=False) / w.sum()
    present = mat.notna()
    num = (mat.fillna(0.0) * w).sum(axis=1)
    den = (present * w).sum(axis=1)
    return num / den.where(den > 0)


def apply_insilico_filter(
    df: pd.DataFrame,
    predicates: list[FilterPredicate],
    *,
    group_col: str | None = None,
    label: str = DEFAULT_LABEL,
    pass_col: str = "insilico_pass",
    rank_col: str = "insilico_avg_rank",
    composite_col: str = "insilico_group",
    group_sep: str = DEFAULT_GROUP_SEP,
) -> pd.DataFrame:
    """Annotate ``df`` with in-silico pass flag, avg rank, and composite name.

    Returns a copy with three added columns:

    - ``pass_col`` — bool, passes all predicates (composite membership).
    - ``rank_col`` — average percentile rank (0 = best), for top-k ranking.
    - ``composite_col`` — for passing rows, ``"<group>+<label>"`` (or just
      ``label`` when no ``group_col``); NaN for non-passing rows.

    Does not mutate the input.
    """
    out = df.copy()
    passed = insilico_mask(out, predicates, group_col=group_col)
    out[pass_col] = passed
    out[rank_col] = average_percentile_rank(out, predicates, group_col=group_col)
    if group_col is not None and group_col in out.columns:
        composite = out[group_col].astype(str) + group_sep + label
    else:
        composite = pd.Series(label, index=out.index)
    out[composite_col] = composite.where(passed, other=np.nan)
    n_pass = int(passed.sum())
    logger.info(
        "insilico_filter: %d/%d clones pass all %d predicate(s) → '%s'",
        n_pass, len(out), len(predicates), label,
    )
    return out


def expand_insilico_twins(
    df: pd.DataFrame,
    predicates: list[FilterPredicate],
    *,
    group_col: str,
    label: str = DEFAULT_LABEL,
    out_group_col: str = "overlap_group",
    group_sep: str = DEFAULT_GROUP_SEP,
) -> pd.DataFrame:
    """Long table where each assay group gains a filtered ``+label`` twin.

    Every input row keeps its base group in ``out_group_col``; rows passing
    the in-silico filter are **duplicated** with ``out_group_col`` set to
    ``"<group>+<label>"``. N base groups → up to 2N groups, ready for
    :func:`tcrsift.clonotype.compute_sample_overlap_matrices` with
    ``sample_col=out_group_col``.
    """
    annotated = apply_insilico_filter(
        df, predicates, group_col=group_col, label=label,
        composite_col="_composite", group_sep=group_sep,
    )
    base = annotated.copy()
    base[out_group_col] = base[group_col].astype(str)
    twins = annotated[annotated["insilico_pass"]].copy()
    twins[out_group_col] = twins["_composite"]
    combined = pd.concat([base, twins], ignore_index=True)
    return combined.drop(columns=["_composite"])


def insilico_overlap_long(
    df: pd.DataFrame,
    predicates: list[FilterPredicate],
    *,
    group_col: str,
    clone_col: str = "CDR3ab",
    cells_col: str = "cells",
    label: str = DEFAULT_LABEL,
    group_sep: str = DEFAULT_GROUP_SEP,
) -> dict[str, pd.DataFrame]:
    """Jaccard-overlap matrices across base groups **and** their ``+label`` twins.

    Expands twins via :func:`expand_insilico_twins`, then computes the
    standard clone-set and cell-weighted Jaccard matrices over the combined
    group axis. The result lets you ask whether ``AIMpos+IS`` overlaps
    ``tetpos+IS`` more than the unfiltered versions do.
    """
    from .clonotype import compute_sample_overlap_matrices

    expanded = expand_insilico_twins(
        df, predicates, group_col=group_col, label=label,
        out_group_col="overlap_group", group_sep=group_sep,
    )
    return compute_sample_overlap_matrices(
        expanded, clone_col=clone_col, sample_col="overlap_group",
        cells_col=cells_col,
    )
