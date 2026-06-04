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

import pandas as pd

from .clonotype import build_clone_method_long, build_clone_sample_long
from .filter import DEFAULT_THRESHOLD_TIERS, per_sample_tier

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


def build_selection_rules(
    clone_method_long: pd.DataFrame,
    config: dict,
    *,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
    tier_col: str = "tier",
    freq_col: str = "max_freq_in_method",
    exclude_clones: set | None = None,
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
            rows.sort(key=lambda r: r["value"], reverse=True)
            if top_n is not None:
                rows = rows[:top_n]
            rule_rows["shared"] = [
                {"clone": r["clone"], "selection_rule": "shared",
                 "ranking_metric": metric, "ranking_value": r["value"]}
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
                cands.sort(key=lambda r: r["value"], reverse=True)
                chosen = cands[:top_n]
                rule_rows[f"private_{method}"] = [
                    {"clone": r["clone"],
                     "selection_rule": f"private_{method}",
                     "ranking_metric": metric, "ranking_value": r["value"]}
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
                cands.sort(key=lambda r: r["value"], reverse=True)
                chosen = cands[:top_n]
                rule_rows[f"method_pair_{name}"] = [
                    {"clone": r["clone"],
                     "selection_rule": f"method_pair_{name}",
                     "ranking_metric": metric, "ranking_value": r["value"]}
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
    freq = (
        out[freq_col] if freq_col in out.columns
        else pd.Series(0.0, index=out.index)
    )
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
            import numpy as np

            candidate = candidate & (np.log10(ppost.clip(lower=1e-300)) <= abs_log10_ppost)
    out["specificity_candidate"] = candidate
    return out
