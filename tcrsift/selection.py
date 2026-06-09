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

"""Shared clone-selection / per-sample-evidence primitives (#125).

This module is the single importable home for the per-``(clone, sample)``
and per-``(clone, method)`` building blocks that report, selection, and
cohort code all depend on. The primitives themselves already live in
:mod:`tcrsift.filter` and :mod:`tcrsift.clonotype`; gathering them here
makes the shared dependency explicit at the library level instead of
buried inside subcommand internals, and gives downstream Python
consumers one place to import from.

It is intentionally **not** wired into the CLI — these are functions that
existing commands (``filter``, ``assemble``, ``run``, and the future
``cohort``) and external callers import.

Foundation increment of the #125 CLI-first roadmap.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd

from .clonotype import build_clone_method_long, build_clone_sample_long
from .filter import DEFAULT_THRESHOLD_TIERS, per_sample_tier
from .insilico_filter import PRISM_PREDICATES as _PRISM_PREDICATES
from .validation import TCRsiftValidationError

logger = logging.getLogger(__name__)

# Strict→permissive walk order shared by the per-sample tier helpers.
DEFAULT_TIER_ORDER: tuple[str, ...] = (
    "tier1", "tier2", "tier3", "tier4", "tier5",
)

__all__ = [
    "DEFAULT_THRESHOLD_TIERS",
    "DEFAULT_TIER_ORDER",
    "build_clone_sample_long",
    "build_clone_method_long",
    "per_sample_tier",
    "attach_per_sample_tiers",
    "attach_method_tiers",
    "build_selection_rules",
    "select_from_clone_sample_long",
    "extract_per_method_evidence",
    "build_pdf_annotations",
    "select_specificity_candidates",
    "pivot_per_sample_tiers",
    "percentile_rank_score",
    "select_freq_prism_per_condition",
    "prism_candidates",
    "freq_prism_grid",
    "PRISM_TERMS",
    "DEFAULT_TIE_BREAK",
]

# Tier label -> numeric rank (lower = stronger enrichment). Used to
# evaluate "tier3 or better" predicates in the selection language.
_TIER_RANK: dict[str, int] = {
    "tier1": 1, "tier2": 2, "tier3": 3, "tier4": 4, "tier5": 5,
}
_UNRANKED = 99  # tier None / unknown — never meets a threshold.

# Default predicate threshold for the "tier3+" exclusion rules (#122).
_DEFAULT_EXCLUDE_TIER = "tier3"


def _tier_rank(tier: object) -> int:
    if tier is None or (isinstance(tier, float) and pd.isna(tier)):
        return _UNRANKED
    return _TIER_RANK.get(str(tier), _UNRANKED)


def _meets_tier(tier: object, threshold: object) -> bool:
    """True when ``tier`` is at least as strong as ``threshold``."""
    return _tier_rank(tier) <= _tier_rank(threshold)


def attach_per_sample_tiers(
    clone_sample_long: pd.DataFrame,
    *,
    tier_defs: dict[str, dict] | None = None,
    tier_order: tuple[str, ...] = DEFAULT_TIER_ORDER,
    cells_col: str = "cells",
    freq_col: str = "frequency",
    out_col: str = "per_sample_tier",
) -> pd.DataFrame:
    """Return a copy of a clone-sample-long frame with a per-sample tier.

    Labels each ``(clone, sample)`` row with its abundance tier via
    :func:`tcrsift.filter.per_sample_tier`, using the row's within-sample
    ``cells`` and ``frequency``. The global clonotype ``tier`` summarises
    a clone's peak signal across all samples; this is the per-sample
    equivalent that report code needs.

    Parameters
    ----------
    clone_sample_long
        Output of :func:`build_clone_sample_long` (one row per
        ``(clone, sample)`` with ``cells`` and ``frequency`` columns).
    tier_defs
        Tier definition dict; defaults to
        :data:`tcrsift.filter.DEFAULT_THRESHOLD_TIERS` (abundance-only).
    tier_order
        Strict→permissive walk order.
    cells_col, freq_col
        Source column names for per-sample cell count and frequency.
    out_col
        Name of the tier column to add.

    Returns
    -------
    pd.DataFrame
        A copy with ``out_col`` added (values are tier labels like
        ``"tier1"`` or ``None`` for rows below the most permissive tier).
    """
    df = clone_sample_long.copy()
    if df.empty:
        df[out_col] = pd.Series([], dtype=object)
        return df
    missing = {cells_col, freq_col} - set(df.columns)
    if missing:
        raise KeyError(
            f"clone_sample_long is missing required column(s) {sorted(missing)}; "
            "expected the output of build_clone_sample_long."
        )
    df[out_col] = [
        per_sample_tier(cells, freq, tier_defs=tier_defs, tier_order=tier_order)
        for cells, freq in zip(df[cells_col], df[freq_col])
    ]
    return df


def attach_method_tiers(
    clone_method_long: pd.DataFrame,
    *,
    tier_defs: dict[str, dict] | None = None,
    tier_order: tuple[str, ...] = DEFAULT_TIER_ORDER,
    cells_col: str = "cells_in_method",
    freq_col: str = "max_freq_in_method",
    out_col: str = "tier",
) -> pd.DataFrame:
    """Return a copy of a clone-method-long frame with a per-method tier.

    Per-method analogue of :func:`attach_per_sample_tiers`, operating on
    the output of :func:`build_clone_method_long` (one row per
    ``(clone, method)`` with ``cells_in_method`` / ``max_freq_in_method``).
    The per-method tier is the input the selection language uses to
    evaluate ``include_tier`` / ``exclude_tier3plus`` predicates.
    """
    df = clone_method_long.copy()
    if df.empty:
        df[out_col] = pd.Series([], dtype=object)
        return df
    missing = {cells_col, freq_col} - set(df.columns)
    if missing:
        raise KeyError(
            f"clone_method_long is missing required column(s) {sorted(missing)}; "
            "expected the output of build_clone_method_long."
        )
    df[out_col] = [
        per_sample_tier(cells, freq, tier_defs=tier_defs, tier_order=tier_order)
        for cells, freq in zip(df[cells_col], df[freq_col])
    ]
    return df


# Columns emitted by build_selection_rules, in order.
_SELECTION_COLUMNS: tuple[str, ...] = (
    "selection_rule", "rank_within_rule", "ranking_metric",
    "ranking_value", "global_rank",
)


def _empty_selection(clone_col: str) -> pd.DataFrame:
    return pd.DataFrame(columns=[clone_col, *_SELECTION_COLUMNS])


def percentile_rank_score(
    clone_scores: pd.DataFrame,
    terms: list[dict],
    *,
    clone_col: str = "CDR3ab",
    group_col: str | None = None,
    require_complete: bool = False,
) -> dict:
    """Weighted mean-percentile composite per clone — the PRISM math (#193).

    ``terms`` is a list of ``{"col", "direction": "low"|"high", "weight"}``.
    Each term contributes a ``[0, 1]`` percentile rank (0 = best;
    ``direction="high"`` means high values are best — e.g. an antigen-response
    GEX score), averaged (optionally weighted). Returns ``{clone: composite}``
    where a **low** composite = a strong multi-criterion candidate (private +
    antigen-responding + non-naive, when given those terms). Percentile ranks
    are computed within ``group_col`` (e.g. per condition) when provided.

    Missing-data policy (explicit, #193 audit): with ``require_complete=False``
    (default) a clone missing one term still ranks on its remaining terms
    (NaN-aware) — used by the config-driven ``select`` ``rank_by`` route, and
    NaN never corrupts ``_rank_candidates`` ordering (it coerces to +inf). With
    ``require_complete=True`` a clone missing ANY term gets a NaN composite and
    is therefore never PRISM-picked — "don't pick a clone we can't fully score",
    which is the behavior the validated per-condition recipe relies on.

    This lets a ``tcrsift select`` route ``rank_by`` reproduce PRISM exactly by
    combining Ppost with GEX/RNA score terms — a deliberate, config-driven
    choice vs. the validated freq + low-Ppost route (#146, GEX excluded).
    """
    from .insilico_filter import percentile_rank

    if clone_scores is None or len(clone_scores) == 0 or not terms:
        return {}
    weighted = []
    for t in terms:
        col = t.get("col")
        if not col or col not in clone_scores.columns:
            continue
        lower = str(t.get("direction", "low")).lower() == "low"
        w = float(t.get("weight", 1.0))
        if group_col and group_col in clone_scores.columns:
            pr = clone_scores.groupby(group_col, observed=True)[col].transform(
                lambda s, _l=lower: percentile_rank(s, lower_is_better=_l)
            )
        else:
            pr = percentile_rank(clone_scores[col], lower_is_better=lower)
        weighted.append((pr, w))
    if not weighted:
        return {}
    prs = np.vstack([np.asarray(pr, dtype=float) for pr, _ in weighted])
    wts = np.array([w for _, w in weighted], dtype=float)
    if require_complete:
        # Any missing term -> NaN composite (clone not fully scorable).
        num = (prs * wts[:, None]).sum(axis=0)
        composite = num / wts.sum()
    else:
        # NaN-aware weighted mean over the terms that ARE present; a clone
        # missing every term gets NaN (no information).
        present = ~np.isnan(prs)
        num = np.nansum(prs * wts[:, None], axis=0)
        den = (present * wts[:, None]).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            composite = np.where(den > 0, num / den, np.nan)
    return dict(zip(clone_scores[clone_col].astype(str), composite))


def _composite_for(metric, clone_scores, clone_col):
    """Per-clone composite dict when ``rank_by`` is a percentile spec, else None."""
    if isinstance(metric, dict) and metric.get("percentile_rank"):
        return percentile_rank_score(
            clone_scores, metric["percentile_rank"],
            clone_col=clone_col, group_col=metric.get("group_col"),
        )
    return None


def _rank_candidates(cands: list[dict], composite) -> str | None:
    """Sort candidate dicts in place. With a ``composite`` map → ascending by
    composite (low = best); else descending by the frequency ``value``. Returns
    the metric label to record (``"percentile_rank"`` for a composite, else
    None so the caller keeps its configured metric)."""
    if composite is not None:
        for c in cands:
            v = composite.get(str(c["clone"]))
            # A missing OR NaN composite sorts last — Python's list.sort does
            # not sink NaN, so it must be coerced to +inf to avoid corrupting
            # the order of the well-scored candidates around it.
            c["value"] = float(v) if v is not None and np.isfinite(v) else float("inf")
        # Secondary key = clone id so ties are reproducible across runs/versions
        # (#221), not an artifact of input order.
        cands.sort(key=lambda r: (r["value"], str(r["clone"])))
        return "percentile_rank"
    cands.sort(key=lambda r: (-float(r["value"]), str(r["clone"])))
    return None


def build_selection_rules(
    clone_method_long: pd.DataFrame,
    config: dict,
    *,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
    tier_col: str = "tier",
    freq_col: str = "max_freq_in_method",
    exclude_clones: set | None = None,
    clone_scores: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Apply a config-driven multi-rule selection language (#122/#125).

    Consumes a per-``(clone, method)`` table carrying a per-method
    abundance ``tier`` (see :func:`attach_method_tiers`) and per-method
    frequency, and assigns each selected clone to a rule:

    * ``shared`` — clones whose strongest per-method tier is in
      ``include_tiers``, ranked by max frequency across methods.
    * ``private_<method>`` — top ``top_n`` clones that meet
      ``include_tier`` in one method and are **not** ``tier3+`` (or the
      configured ``exclude_tier``) in any other method.
    * ``method_pair_<name>`` — top ``top_n`` clones that meet
      ``require_tier_in_all_members`` in **every** member of a method
      pair and are not ``tier3+`` outside the pair, ranked by mean
      frequency across the pair.

    Routes are evaluated in ``global_rank.rule_order`` (default
    ``[shared, method_pair, private]``) and a clone is assigned to the first
    rule that claims it — no clone appears twice. ``global_rank`` is a
    1-based interleave of the blocks in that order.

    Returns a DataFrame with ``clone_col`` plus ``selection_rule``,
    ``rank_within_rule``, ``ranking_metric``, ``ranking_value``,
    ``global_rank``. Empty input or no matches yields an empty frame with
    those columns.
    """
    rules_cfg = config.get("rules", {}) or {}
    global_cfg = config.get("global_rank", {}) or {}
    rule_order = global_cfg.get(
        "rule_order", ["shared", "method_pair", "private"]
    )

    if clone_method_long.empty or not rules_cfg:
        return _empty_selection(clone_col)

    # Per-clone method->tier and method->frequency maps.
    tier_by: dict[object, dict] = {}
    freq_by: dict[object, dict] = {}
    for clone, g in clone_method_long.groupby(
        clone_col, sort=False, observed=True,
    ):
        tier_by[clone] = dict(zip(g[method_col], g[tier_col]))
        freq_by[clone] = dict(zip(g[method_col], g[freq_col]))
    # Drop excluded clones up front (e.g. public-DB viral bystanders) so
    # no rule can select them (#122 exclude_viral).
    exclude_clones = exclude_clones or set()
    all_clones = [c for c in tier_by if c not in exclude_clones]

    selected: set = set()
    # rule_name -> list of dicts (clone, selection_rule, ranking_metric,
    # ranking_value), already ranked within the rule.
    rule_rows: dict[str, list[dict]] = {}

    def _max_other_tier_rank(clone, exclude_methods: set) -> int:
        ranks = [
            _tier_rank(t)
            for m, t in tier_by[clone].items()
            if m not in exclude_methods
        ]
        return min(ranks) if ranks else _UNRANKED

    for block in rule_order:
        cfg = rules_cfg.get(block)
        if not cfg:
            continue

        if block == "shared":
            include = set(cfg.get("include_tiers", []))
            metric = cfg.get("rank_by", "max_frequency")
            # Optional bounds so the shared block doesn't dominate the
            # output: a frequency floor and a top-N cap (#122 follow-up).
            min_freq = cfg.get("min_frequency")
            top_n = cfg.get("top_n")  # None = uncapped (back-compat)
            rows = []
            for c in all_clones:
                if c in selected:
                    continue
                best_rank = min(
                    (_tier_rank(t) for t in tier_by[c].values()),
                    default=_UNRANKED,
                )
                best_label = next(
                    (t for t in tier_by[c].values()
                     if _tier_rank(t) == best_rank),
                    None,
                )
                if best_label in include:
                    val = max(freq_by[c].values(), default=0.0)
                    if min_freq is not None and val < min_freq:
                        continue
                    rows.append({"clone": c, "value": float(val)})
            metric_label = _rank_candidates(
                rows, _composite_for(metric, clone_scores, clone_col)
            ) or metric
            if top_n is not None:
                rows = rows[:top_n]
            rule_rows["shared"] = [
                {"clone": r["clone"], "selection_rule": "shared",
                 "ranking_metric": metric_label, "ranking_value": r["value"]}
                for r in rows
            ]
            selected.update(r["clone"] for r in rows)

        elif block == "private":
            include_tier = cfg.get("include_tier", "tier3")
            exclude_tier = (
                _DEFAULT_EXCLUDE_TIER
                if cfg.get("exclude_tier3plus_in_other_methods", True)
                else None
            )
            top_n = cfg.get("top_n", 3)
            metric = cfg.get("rank_by", "frequency_in_method")
            methods = cfg.get("apply_to_methods")
            if methods is None:
                methods = sorted(
                    {m for maps in tier_by.values() for m in maps}
                )
            for method in methods:
                cands = []
                for c in all_clones:
                    if c in selected:
                        continue
                    if not _meets_tier(tier_by[c].get(method), include_tier):
                        continue
                    if exclude_tier is not None:
                        other = _max_other_tier_rank(c, {method})
                        if other <= _tier_rank(exclude_tier):
                            continue
                    cands.append(
                        {"clone": c, "value": float(freq_by[c].get(method, 0.0))}
                    )
                metric_label = _rank_candidates(
                    cands, _composite_for(metric, clone_scores, clone_col)
                ) or metric
                chosen = cands[:top_n]
                rule_rows[f"private_{method}"] = [
                    {"clone": r["clone"],
                     "selection_rule": f"private_{method}",
                     "ranking_metric": metric_label, "ranking_value": r["value"]}
                    for r in chosen
                ]
                selected.update(r["clone"] for r in chosen)

        elif block == "method_pair":
            pairs = cfg.get("pairs", {})
            require = cfg.get("require_tier_in_all_members", "tier3")
            exclude_tier = (
                _DEFAULT_EXCLUDE_TIER
                if cfg.get("exclude_tier3plus_outside_pair", True)
                else None
            )
            top_n = cfg.get("top_n", 3)
            metric = cfg.get("rank_by", "mean_freq_across_pair")
            for name, members in pairs.items():
                members = list(members)
                member_set = set(members)
                cands = []
                for c in all_clones:
                    if c in selected:
                        continue
                    if not all(
                        _meets_tier(tier_by[c].get(m), require) for m in members
                    ):
                        continue
                    if exclude_tier is not None:
                        other = _max_other_tier_rank(c, member_set)
                        if other <= _tier_rank(exclude_tier):
                            continue
                    mean_freq = sum(
                        float(freq_by[c].get(m, 0.0)) for m in members
                    ) / max(len(members), 1)
                    cands.append({"clone": c, "value": mean_freq})
                metric_label = _rank_candidates(
                    cands, _composite_for(metric, clone_scores, clone_col)
                ) or metric
                chosen = cands[:top_n]
                rule_rows[f"method_pair_{name}"] = [
                    {"clone": r["clone"],
                     "selection_rule": f"method_pair_{name}",
                     "ranking_metric": metric_label, "ranking_value": r["value"]}
                    for r in chosen
                ]
                selected.update(r["clone"] for r in chosen)

    # Assemble in block order, interleaving sub-rules (private_*, pair_*)
    # in config order, with a 1-based global rank.
    ordered_rule_names: list[str] = []
    for block in rule_order:
        if block == "shared":
            if "shared" in rule_rows:
                ordered_rule_names.append("shared")
        elif block == "private":
            cfg = rules_cfg.get("private") or {}
            methods = cfg.get("apply_to_methods") or sorted(
                {m for maps in tier_by.values() for m in maps}
            )
            ordered_rule_names += [
                f"private_{m}" for m in methods if f"private_{m}" in rule_rows
            ]
        elif block == "method_pair":
            cfg = rules_cfg.get("method_pair") or {}
            ordered_rule_names += [
                f"method_pair_{n}" for n in cfg.get("pairs", {})
                if f"method_pair_{n}" in rule_rows
            ]

    out_rows: list[dict] = []
    global_rank = 0
    for rule_name in ordered_rule_names:
        for within, row in enumerate(rule_rows.get(rule_name, []), start=1):
            global_rank += 1
            out_rows.append({
                clone_col: row["clone"],
                "selection_rule": row["selection_rule"],
                "rank_within_rule": within,
                "ranking_metric": row["ranking_metric"],
                "ranking_value": row["ranking_value"],
                "global_rank": global_rank,
            })

    if not out_rows:
        return _empty_selection(clone_col)
    return pd.DataFrame(out_rows, columns=[clone_col, *_SELECTION_COLUMNS])


def select_from_clone_sample_long(
    clone_sample_long: pd.DataFrame,
    config: dict,
    *,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
    exclude_clones: set | None = None,
    clone_scores: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """End-to-end selection from a clone-sample-long table.

    Convenience wrapper that runs the full chain so callers (``tcrsift
    run`` and downstream Python) don't re-stitch it:

    ``clone_sample_long`` → :func:`build_clone_method_long` →
    :func:`attach_method_tiers` → :func:`build_selection_rules`.

    Returns the selection frame (one row per selected clone with
    ``selection_rule`` / ``rank_within_rule`` / ``ranking_metric`` /
    ``ranking_value`` / ``global_rank``). Empty when there's no method
    axis, no rules configured, or nothing matches.
    """
    if clone_sample_long.empty or not (config or {}).get("rules"):
        return _empty_selection(clone_col)
    # No method axis → no per-method tiers → nothing the rule language
    # can act on (the rules are all defined over methods/method-pairs).
    if method_col not in clone_sample_long.columns:
        return _empty_selection(clone_col)
    cml = build_clone_method_long(
        clone_sample_long, clone_col=clone_col, method_col=method_col
    )
    if cml.empty:
        return _empty_selection(clone_col)
    cml = attach_method_tiers(cml)
    return build_selection_rules(
        cml, config, clone_col=clone_col, exclude_clones=exclude_clones,
        clone_scores=clone_scores,
    )


def extract_per_method_evidence(
    clone_method_long: pd.DataFrame,
    clone: object,
    *,
    methods: list | None = None,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
) -> pd.DataFrame:
    """Per-method evidence rows for one clone (#123/#125).

    Returns the ``(method, cells_in_method, max_freq_in_method, tier)``
    rows for ``clone`` from a clone-method-long table (attach per-method
    tiers first via :func:`attach_method_tiers`). Optionally restricted
    to ``methods`` and ordered by that list.
    """
    sub = clone_method_long[clone_method_long[clone_col] == clone]
    if methods is not None:
        sub = sub[sub[method_col].isin(methods)]
        order = {m: i for i, m in enumerate(methods)}
        sub = sub.assign(
            __o=sub[method_col].map(order)
        ).sort_values("__o").drop(columns="__o")
    cols = [method_col] + [
        c for c in ("cells_in_method", "max_freq_in_method", "tier")
        if c in sub.columns
    ]
    return sub[cols].reset_index(drop=True)


def build_pdf_annotations(
    clone_method_long: pd.DataFrame,
    *,
    selection_df: pd.DataFrame | None = None,
    methods: list | None = None,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
) -> dict:
    """Build per-clone annotation line-blocks for the sequence PDF.

    Returns ``{clone: [line, ...]}``. Each block lists the clone's
    per-method evidence (cells, within-method frequency, tier) and, when
    ``selection_df`` (output of :func:`build_selection_rules`) is given,
    a header with the clone's selection rule and global rank.
    """
    rule_by: dict = {}
    if selection_df is not None and not selection_df.empty:
        for _, r in selection_df.iterrows():
            rule_by[r[clone_col]] = (r["selection_rule"], r["global_rank"])

    out: dict = {}
    for clone, ev in clone_method_long.groupby(
        clone_col, sort=False, observed=True,
    ):
        lines: list[str] = []
        if clone in rule_by:
            rule, rank = rule_by[clone]
            lines.append(f"selection: {rule} (global rank {int(rank)})")
        ev = extract_per_method_evidence(
            clone_method_long, clone, methods=methods,
            clone_col=clone_col, method_col=method_col,
        )
        for _, r in ev.iterrows():
            cells = r.get("cells_in_method")
            freq = r.get("max_freq_in_method")
            tier = r.get("tier")
            cells_s = "" if cells is None or pd.isna(cells) else f"{int(cells)} cells"
            freq_s = "" if freq is None or pd.isna(freq) else f"{100 * float(freq):.1f}%"
            tier_s = "-" if tier is None or (isinstance(tier, float) and pd.isna(tier)) else str(tier)
            lines.append(
                f"  {r[method_col]}: {cells_s}, {freq_s}  [{tier_s}]".replace(", ,", ",")
            )
        out[clone] = lines
    return out


def select_specificity_candidates(
    clonotypes: pd.DataFrame,
    *,
    freq_col: str = "max_frequency_per_method",
    min_frequency: float = 0.001,
    ppost_col: str = "ppost_alpha",
    percentile: float = 25.0,
    abs_log10_ppost: float | None = None,
    alpha_col: str = "CDR3_alpha",
    beta_col: str = "CDR3_beta",
) -> pd.DataFrame:
    """Flag private / rare-precursor clones — candidate high-specificity (#146).

    The validated B1-2/B1-3 recipe (building on #143 Pgen/Ppost): among complete
    αβ clones above a frequency gate, rank by **low Ppost(α)** (rare precursor =
    candidate private receptor) and flag the bottom ``percentile`` (default
    bottom quartile). RNA terms were tested and reproduce in only one donor, so
    they are deliberately **not** part of this gate — only the sequence axis
    (Ppost) is used.

    Adds:
    - ``specificity_gated`` — passes the frequency + complete-αβ gate.
    - ``specificity_ppost_rank`` — 1 = lowest Ppost(α) among gated clones.
    - ``specificity_candidate`` — gated AND in the bottom ``percentile`` of
      Ppost(α) (and, if ``abs_log10_ppost`` is set, also ``log10 Ppost ≤`` it;
      the spec's absolute guide is ≈ −7).

    Caveat (documented intent): this prioritizes private / rare-precursor
    candidates as a *proxy* for high specificity / low cross-reactivity. It is
    **not** validated against true specificity or avidity — confirmation needs a
    day-0 precursor baseline and a functional-avidity assay. It ranks
    candidates; it does not measure specificity. Returns a copy.
    """
    out = clonotypes.copy()
    if ppost_col not in out.columns:
        raise ValueError(
            f"select_specificity_candidates: missing {ppost_col!r}. "
            "Run add_pgen_ppost / annotate_tcrs first."
        )

    def _valid(col):
        if col not in out.columns:
            return pd.Series(False, index=out.index)
        s = out[col].astype("string")
        return s.notna() & (s.str.len() > 0) & (s.str.lower() != "nan")

    complete = _valid(alpha_col) & _valid(beta_col)
    # Frequency gate. Prefer the per-method column; fall back to max_frequency
    # so the route still works on runs without an enrichment_method axis (a
    # missing per-method column would otherwise gate out every clone silently).
    if freq_col in out.columns:
        freq = out[freq_col]
    elif "max_frequency" in out.columns:
        logger.warning(
            "select_specificity_candidates: %r absent; gating on 'max_frequency'.",
            freq_col,
        )
        freq = out["max_frequency"]
    else:
        logger.warning(
            "select_specificity_candidates: no frequency column (%r / "
            "'max_frequency'); frequency gate disabled (all complete clones pass).",
            freq_col,
        )
        freq = pd.Series(float("inf"), index=out.index)
    gated = complete & (freq.fillna(0) > min_frequency) & out[ppost_col].notna()
    out["specificity_gated"] = gated

    # Rank gated clones by ascending Ppost(alpha) (lowest = most private).
    ppost = out[ppost_col]
    out["specificity_ppost_rank"] = (
        ppost.where(gated).rank(method="min", ascending=True).astype("Int64")
    )

    candidate = pd.Series(False, index=out.index)
    if gated.any():
        cutoff = ppost[gated].quantile(percentile / 100.0)
        candidate = gated & (ppost <= cutoff)
        if abs_log10_ppost is not None:
            candidate = candidate & (np.log10(ppost.clip(lower=1e-300)) <= abs_log10_ppost)
    out["specificity_candidate"] = candidate
    return out


def pivot_per_sample_tiers(
    clone_sample_long: pd.DataFrame,
    *,
    clone_col: str = "CDR3ab",
    out_prefix: str = "tier_in_",
) -> pd.DataFrame:
    """Wide per-sample tier table for `tcrsift filter --emit-per-sample-tier` (#122).

    Labels each ``(clone, sample)`` row with its within-sample abundance tier
    (:func:`attach_per_sample_tiers`) and pivots to one ``tier_in_{sample}``
    column per sample, keyed by ``clone_col``. The selection language consumes
    these to evaluate per-method ``include_tier`` / ``exclude_tier3plus``
    predicates without recomputing. Returns a frame with ``clone_col`` + one
    column per sample (empty frame with just ``clone_col`` on empty input).
    """
    if clone_sample_long.empty or "sample" not in clone_sample_long.columns:
        return pd.DataFrame(columns=[clone_col])
    tiered = attach_per_sample_tiers(clone_sample_long)
    wide = tiered.pivot_table(
        index=clone_col, columns="sample", values="per_sample_tier", aggfunc="first"
    )
    wide.columns = [f"{out_prefix}{c}" for c in wide.columns]
    return wide.reset_index()


# Default PRISM terms (tcrsift column names), matching the validated
# selection_analysis recipe: low private Ppost + high antigen-response GEX +
# low naive GEX. Derived from the single source of truth
# (insilico_filter.PRISM_PREDICATES) so the dict form here and the
# FilterPredicate form in annotate_tcrs can't drift. Override `score_terms`
# to change the criteria.
PRISM_TERMS: list[dict] = [
    {"col": _p.score, "direction": _p.direction, "weight": 1.0}
    for _p in _PRISM_PREDICATES
]


# Principled, deterministic, NON-PRISM tie-break (#221/#223). Among clones tied
# on the primary key (frequency, or an exact PRISM score), prefer — in order:
#   single_alpha — a single-α clone over a merged dual-α (allelic-inclusion) one
#   umi          — higher per-clone UMI coverage (better-sequenced)
#   cd8_purity   — higher phenotype purity (consistent CD8/CD4 across the clone)
#   cdr3ab       — lexical CDR3ab, full determinism (always terminates the chain)
# Keys whose column is absent are skipped; cdr3ab guarantees a reproducible order
# regardless of input row order or version (fixing #221 even with no metadata).
DEFAULT_TIE_BREAK: tuple[str, ...] = ("single_alpha", "umi", "cd8_purity", "cdr3ab")


def _tie_break_sort(cand: pd.DataFrame, tie_break, clone_col: str):
    """Return ``(sort_cols, ascending)`` for the configured tie-break, adding
    private ``_tb_*`` derived columns to ``cand``. Always ends on ``clone_col``."""
    cols: list[str] = []
    asc: list[bool] = []
    for key in tie_break or ():
        if key == "single_alpha" and "merged_alpha_partners" in cand.columns:
            p = cand["merged_alpha_partners"].astype(str).str.strip().str.lower()
            cand["_tb_dual"] = cand["merged_alpha_partners"].notna() & ~p.isin(("", "nan"))
            cols.append("_tb_dual")
            asc.append(True)  # single (False) before dual (True)
        elif key == "umi":
            umi_cols = [c for c in ("TRA_1_umis", "TRB_1_umis") if c in cand.columns]
            if umi_cols:
                cand["_tb_umi"] = (
                    cand[umi_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
                )
                cols.append("_tb_umi")
                asc.append(False)  # higher coverage first
        elif key in ("cd8_purity", "cd4_purity", "purity") and "Tcell_type_purity" in cand.columns:
            cand["_tb_purity"] = pd.to_numeric(
                cand["Tcell_type_purity"], errors="coerce"
            ).fillna(0.0)
            cols.append("_tb_purity")
            asc.append(False)  # higher purity first
        elif key == "cdr3ab" and clone_col not in cols:
            cols.append(clone_col)
            asc.append(True)
    if clone_col not in cols:  # always terminate deterministically
        cols.append(clone_col)
        asc.append(True)
    return cols, asc


def select_freq_prism_per_condition(
    feat: pd.DataFrame,
    *,
    condition_col: str,
    freq_col: str = "frequency",
    clone_col: str = "CDR3ab",
    score_terms: list[dict] | None = None,
    gate: float = 0.001,
    top_freq: int = 10,
    top_prism: int = 5,
    rescue_target: int = 0,
    rescue_rank_col: str | None = None,
    tie_break: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Per-condition ``top-K freq ∪ top-K PRISM`` selection (#193).

    Faithful reproduction of the validated downstream recipe (selection_analysis
    ``make_plots.py``), so it can be retired without changing the science. For
    each condition: gate to ``freq > gate``; compute PRISM = mean percentile
    rank of ``score_terms`` **within that condition's candidate set** (low =
    best — exactly the downstream's per-condition ranking); then select the
    ``top_freq`` clones by frequency UNION the ``top_prism`` clones by PRISM.

    ``score_terms`` defaults to :data:`PRISM_TERMS` (low ppost_alpha/beta, high
    antigen_response_score, low naive_score). Optional coverage rescue: when a
    condition has fewer than ``rescue_target`` gated clones, the best
    sub-threshold clones by ``rescue_rank_col`` (descending) are added back.

    Note: the ppost terms are ~0.9 |corr| with CDR3 length (longer CDR3 ->
    lower Pgen/Ppost), inherent to any generation-probability score, so the
    PRISM ppost axes partly rank on length, not just rarity. Raises (when
    ``top_prism > 0``) if a PRISM score column is missing or entirely NaN, since
    PRISM would otherwise pick nothing.

    Both route sorts use a deterministic, principled, NON-PRISM ``tie_break``
    (default :data:`DEFAULT_TIE_BREAK` — single-α > dual-α, then UMI coverage,
    then phenotype purity, then CDR3ab) so boundary ties among equally-frequent
    (or equally-scored) clones are reproducible and quality-ranked rather than
    an arbitrary artifact of input row order (#221/#223). Tie-break keys whose
    column is absent are skipped; the CDR3ab fallback always guarantees
    determinism. Pass ``[]`` for CDR3ab-only.

    Returns one row per (clone, condition) with ``selection_route`` (``freq`` /
    ``prism`` / ``both``), ``rank_within_route``, ``prism_score``, the
    condition, the frequency, and ``rescued`` (True when the clone cleared no
    gate and was added back by coverage rescue — a low-confidence pick). A
    clone selected in multiple conditions appears once per condition.
    """
    terms = score_terms or PRISM_TERMS
    if feat.empty or condition_col not in feat.columns:
        return pd.DataFrame(
            columns=[clone_col, condition_col, "selection_route",
                     "rank_within_route", "prism_score", freq_col, "rescued"]
        )
    # Guard: PRISM requested but its score columns are unusable. With
    # require_complete scoring, an all-NaN or missing term means NO clone is
    # fully scorable, so the PRISM route would silently pick nothing (e.g. when
    # ppost was never populated). Fail loudly instead of degenerating.
    if top_prism > 0:
        term_cols = [t.get("col") for t in terms if t.get("col")]
        unusable = [
            c for c in term_cols
            if c not in feat.columns or not feat[c].notna().any()
        ]
        if unusable:
            raise TCRsiftValidationError(
                f"PRISM requested (top_prism={top_prism}) but its score "
                f"column(s) {unusable} are missing or entirely NaN — every "
                "clone is unscorable, so PRISM would pick nothing.",
                hint="Populate the PRISM score columns first (e.g. "
                "annotate_tcrs.add_pgen_ppost / add_gex_signature_scores for "
                "ppost_alpha/ppost_beta/antigen_response_score/naive_score), "
                "pass score_terms for the columns you do have, or set "
                "top_prism=0 to use the frequency route alone.",
            )
    out_rows: list[dict] = []
    for cond, grp in feat.groupby(condition_col, observed=True):
        cand = grp[grp[freq_col].fillna(0) > gate].copy()
        rescued_clones: set[str] = set()
        if rescue_target and len(cand) < rescue_target and rescue_rank_col in grp.columns:
            need = rescue_target - len(cand)
            sub = grp[grp[freq_col].fillna(0) <= gate].copy()
            # Deterministic: break rescue-rank ties by clone id so which thin
            # clones are rescued is reproducible, not input-row-order dependent.
            sub = sub.sort_values(
                [rescue_rank_col, clone_col], ascending=[False, True],
            ).head(need)
            rescued_clones = set(sub[clone_col].astype(str))
            cand = pd.concat([cand, sub], ignore_index=True)
        if cand.empty:
            continue
        # Collapse to one row per clone within the condition. The input may be
        # per-(clone, sample) (e.g. condition_col='method' over clone_sample_long),
        # so a clone present in several samples of one condition would otherwise
        # occupy multiple top-K slots — undercounting distinct clones, skewing the
        # within-condition percentile ranks, and producing wrong rank_within_route.
        # Keep each clone's best (highest-frequency) row.
        # Deterministic dedup: keep each clone's highest-frequency row, breaking
        # ties by clone id so which row survives is reproducible (not row-order).
        cand = (
            cand.sort_values([freq_col, clone_col], ascending=[False, True])
            .drop_duplicates(subset=[clone_col], keep="first")
        )
        # PRISM percentile ranks computed WITHIN this condition's candidates.
        # require_complete=True: a clone missing any term gets a NaN PRISM score
        # and sorts last, so it is not PRISM-picked — "don't pick a clone we
        # can't fully score", matching the validated recipe (a clone scored on a
        # subset of criteria is not comparable to fully-scored ones).
        composite = percentile_rank_score(
            cand, terms, clone_col=clone_col, require_complete=True,
        )
        cand["prism_score"] = cand[clone_col].astype(str).map(composite)
        # Deterministic, principled NON-PRISM tie-break (#221/#223) applied to
        # BOTH routes so tied boundary picks are reproducible and meaningful.
        tb_cols, tb_asc = _tie_break_sort(
            cand, DEFAULT_TIE_BREAK if tie_break is None else tie_break, clone_col,
        )
        top_f = list(
            cand.sort_values([freq_col, *tb_cols], ascending=[False, *tb_asc])[
                clone_col
            ].head(top_freq)
        )
        top_p = list(
            cand.sort_values(["prism_score", *tb_cols], ascending=[True, *tb_asc])[
                clone_col
            ].head(top_prism)
        )
        freq_set, prism_set = set(top_f), set(top_p)
        f_rank = {c: i + 1 for i, c in enumerate(top_f)}
        p_rank = {c: i + 1 for i, c in enumerate(top_p)}
        # Emit in deterministic "order of selected clones": the freq route in
        # frequency-rank order first (incl. 'both'), then prism-only clones in
        # PRISM-rank order. top_f/top_p are already rank-ordered lists; iterating
        # the set union directly would be order-nondeterministic across runs and
        # churn selection_native.csv row order with no semantic change (#230).
        ordered = top_f + [c for c in top_p if c not in freq_set]
        for c in ordered:
            in_f, in_p = c in freq_set, c in prism_set
            route = "both" if (in_f and in_p) else ("freq" if in_f else "prism")
            row = cand[cand[clone_col] == c].iloc[0]
            out_rows.append({
                clone_col: c,
                condition_col: cond,
                "selection_route": route,
                "rank_within_route": f_rank.get(c) if in_f else p_rank.get(c),
                "prism_score": float(row["prism_score"]) if pd.notna(row["prism_score"]) else None,
                freq_col: float(row[freq_col]) if pd.notna(row[freq_col]) else 0.0,
                # True when this clone cleared no condition gate and was added
                # back by coverage rescue (#228) — a low-confidence pick.
                "rescued": str(c) in rescued_clones,
            })
    return pd.DataFrame(
        out_rows,
        columns=[clone_col, condition_col, "selection_route",
                 "rank_within_route", "prism_score", freq_col, "rescued"],
    )


def prism_candidates(
    feat: pd.DataFrame,
    *,
    condition_col: str,
    freq_col: str = "frequency",
    clone_col: str = "CDR3ab",
    score_terms: list[dict] | None = None,
    gate: float = 0.001,
    top_freq: int = 10,
    top_prism: int = 5,
    tie_break: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Per-condition gated-candidate frame for the freq×PRISM selection plots.

    Mirrors :func:`select_freq_prism_per_condition`'s gating, per-clone dedup,
    and within-condition PRISM percentile scoring, but emits EVERY gated
    candidate (``frequency > gate``) — not just the picks — each tagged with the
    ``selection_route`` it would receive (``freq`` / ``prism`` / ``both`` /
    ``unselected``). This is the 2-D selection space (#248) and the PRISM-vs-
    background view (#249): the per-candidate PRISM score only exists at select
    time, so it's surfaced here for plotting.

    Returns one row per (clone, condition) with ``clone_col``, ``condition_col``,
    ``freq_col``, ``prism_score`` (low = better; NaN when not fully scorable),
    and ``selection_route``. Empty when ``feat`` is empty or the condition
    column is absent.
    """
    terms = score_terms or PRISM_TERMS
    if feat.empty or condition_col not in feat.columns:
        return pd.DataFrame(
            columns=[clone_col, condition_col, freq_col, "prism_score",
                     "selection_route"]
        )
    out_rows: list[dict] = []
    for cond, grp in feat.groupby(condition_col, observed=True):
        cand = grp[grp[freq_col].fillna(0) > gate].copy()
        if cand.empty:
            continue
        cand = (
            cand.sort_values([freq_col, clone_col], ascending=[False, True])
            .drop_duplicates(subset=[clone_col], keep="first")
        )
        composite = percentile_rank_score(
            cand, terms, clone_col=clone_col, require_complete=True,
        )
        cand["prism_score"] = cand[clone_col].astype(str).map(composite)
        tb_cols, tb_asc = _tie_break_sort(
            cand, DEFAULT_TIE_BREAK if tie_break is None else tie_break, clone_col,
        )
        freq_set = set(
            cand.sort_values([freq_col, *tb_cols], ascending=[False, *tb_asc])[
                clone_col
            ].head(top_freq)
        )
        prism_set = set(
            cand.sort_values(["prism_score", *tb_cols], ascending=[True, *tb_asc])[
                clone_col
            ].head(top_prism)
        ) if top_prism > 0 else set()
        for _, row in cand.iterrows():
            c = row[clone_col]
            in_f, in_p = c in freq_set, c in prism_set
            route = (
                "both" if (in_f and in_p)
                else "freq" if in_f
                else "prism" if in_p
                else "unselected"
            )
            out_rows.append({
                clone_col: c,
                condition_col: cond,
                freq_col: float(row[freq_col]) if pd.notna(row[freq_col]) else 0.0,
                "prism_score": float(row["prism_score"]) if pd.notna(row["prism_score"]) else None,
                "selection_route": route,
            })
    return pd.DataFrame(
        out_rows,
        columns=[clone_col, condition_col, freq_col, "prism_score",
                 "selection_route"],
    )


def freq_prism_grid(
    feat: pd.DataFrame,
    *,
    condition_col: str,
    freq_values: Sequence[int] = (0, 5, 10, 15, 20),
    prism_values: Sequence[int] = (0, 5, 10, 15, 20),
    freq_col: str = "frequency",
    clone_col: str = "CDR3ab",
    score_terms: list[dict] | None = None,
    gate: float = 0.001,
) -> pd.DataFrame:
    """Sweep ``(top_freq, top_prism)`` and count total distinct clones (#207).

    For each cell of the ``freq_values`` × ``prism_values`` grid, runs
    :func:`select_freq_prism_per_condition` with that ``(top_freq, top_prism)``
    and records the number of *distinct* clones selected across all conditions
    (the union — a clone picked in two conditions counts once). This is the
    trade-off the downstream ``make_plots.py`` heatmap shows when picking how
    many clones to take by frequency vs. by PRISM.

    A ``0`` in either axis means "take none by that route", so the row/column at
    0 isolates the pure-PRISM / pure-frequency yield. Returns a tidy DataFrame
    with columns ``top_freq``, ``top_prism``, ``n_clones`` (one row per cell);
    pivot to ``top_prism`` × ``top_freq`` for a heatmap via
    :func:`tcrsift.plots.plot_freq_prism_grid`.
    """
    rows: list[dict] = []
    for nf in freq_values:
        for npr in prism_values:
            sel = select_freq_prism_per_condition(
                feat,
                condition_col=condition_col,
                freq_col=freq_col,
                clone_col=clone_col,
                score_terms=score_terms,
                gate=gate,
                top_freq=int(nf),
                top_prism=int(npr),
            )
            n = sel[clone_col].nunique() if not sel.empty else 0
            rows.append({"top_freq": int(nf), "top_prism": int(npr), "n_clones": int(n)})
    return pd.DataFrame(rows, columns=["top_freq", "top_prism", "n_clones"])
