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
    "build_selection_routes",
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


# Columns emitted by build_selection_routes, in order.
_SELECTION_COLUMNS: tuple[str, ...] = (
    "selection_route", "rank_within_route", "ranking_metric",
    "ranking_value", "global_rank",
)


def _empty_selection(clone_col: str) -> pd.DataFrame:
    return pd.DataFrame(columns=[clone_col, *_SELECTION_COLUMNS])


def build_selection_routes(
    clone_method_long: pd.DataFrame,
    config: dict,
    *,
    clone_col: str = "CDR3ab",
    method_col: str = "method",
    tier_col: str = "tier",
    freq_col: str = "max_freq_in_method",
) -> pd.DataFrame:
    """Apply a config-driven multi-route selection language (#122/#125).

    Consumes a per-``(clone, method)`` table carrying a per-method
    abundance ``tier`` (see :func:`attach_method_tiers`) and per-method
    frequency, and assigns each selected clone to a route:

    * ``shared`` — clones whose strongest per-method tier is in
      ``include_tiers``, ranked by max frequency across methods.
    * ``private_<method>`` — top ``top_n`` clones that meet
      ``include_tier`` in one method and are **not** ``tier3+`` (or the
      configured ``exclude_tier``) in any other method.
    * ``cty_pair_<name>`` — top ``top_n`` clones that meet
      ``require_tier_in_all_members`` in **every** member of a method
      pair and are not ``tier3+`` outside the pair, ranked by mean
      frequency across the pair.

    Routes are evaluated in ``global_rank.block_order`` (default
    ``[shared, cty_pair, private]``) and a clone is assigned to the first
    route that claims it — no clone appears twice. ``global_rank`` is a
    1-based interleave of the blocks in that order.

    Returns a DataFrame with ``clone_col`` plus ``selection_route``,
    ``rank_within_route``, ``ranking_metric``, ``ranking_value``,
    ``global_rank``. Empty input or no matches yields an empty frame with
    those columns.
    """
    routes_cfg = config.get("routes", {}) or {}
    global_cfg = config.get("global_rank", {}) or {}
    block_order = global_cfg.get(
        "block_order", ["shared", "cty_pair", "private"]
    )

    if clone_method_long.empty or not routes_cfg:
        return _empty_selection(clone_col)

    # Per-clone method->tier and method->frequency maps.
    tier_by: dict[object, dict] = {}
    freq_by: dict[object, dict] = {}
    for clone, g in clone_method_long.groupby(
        clone_col, sort=False, observed=True,
    ):
        tier_by[clone] = dict(zip(g[method_col], g[tier_col]))
        freq_by[clone] = dict(zip(g[method_col], g[freq_col]))
    all_clones = list(tier_by)

    selected: set = set()
    # route_name -> list of dicts (clone, selection_route, ranking_metric,
    # ranking_value), already ranked within the route.
    block_rows: dict[str, list[dict]] = {}

    def _max_other_tier_rank(clone, exclude_methods: set) -> int:
        ranks = [
            _tier_rank(t)
            for m, t in tier_by[clone].items()
            if m not in exclude_methods
        ]
        return min(ranks) if ranks else _UNRANKED

    for block in block_order:
        cfg = routes_cfg.get(block)
        if not cfg:
            continue

        if block == "shared":
            include = set(cfg.get("include_tiers", []))
            metric = cfg.get("rank_by", "max_frequency")
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
                    rows.append({"clone": c, "value": float(val)})
            rows.sort(key=lambda r: r["value"], reverse=True)
            block_rows["shared"] = [
                {"clone": r["clone"], "selection_route": "shared",
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
                block_rows[f"private_{method}"] = [
                    {"clone": r["clone"],
                     "selection_route": f"private_{method}",
                     "ranking_metric": metric, "ranking_value": r["value"]}
                    for r in chosen
                ]
                selected.update(r["clone"] for r in chosen)

        elif block == "cty_pair":
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
                block_rows[f"cty_pair_{name}"] = [
                    {"clone": r["clone"],
                     "selection_route": f"cty_pair_{name}",
                     "ranking_metric": metric, "ranking_value": r["value"]}
                    for r in chosen
                ]
                selected.update(r["clone"] for r in chosen)

    # Assemble in block order, interleaving sub-routes (private_*, pair_*)
    # in config order, with a 1-based global rank.
    ordered_route_names: list[str] = []
    for block in block_order:
        if block == "shared":
            if "shared" in block_rows:
                ordered_route_names.append("shared")
        elif block == "private":
            cfg = routes_cfg.get("private") or {}
            methods = cfg.get("apply_to_methods") or sorted(
                {m for maps in tier_by.values() for m in maps}
            )
            ordered_route_names += [
                f"private_{m}" for m in methods if f"private_{m}" in block_rows
            ]
        elif block == "cty_pair":
            cfg = routes_cfg.get("cty_pair") or {}
            ordered_route_names += [
                f"cty_pair_{n}" for n in cfg.get("pairs", {})
                if f"cty_pair_{n}" in block_rows
            ]

    out_rows: list[dict] = []
    global_rank = 0
    for route_name in ordered_route_names:
        for within, row in enumerate(block_rows.get(route_name, []), start=1):
            global_rank += 1
            out_rows.append({
                clone_col: row["clone"],
                "selection_route": row["selection_route"],
                "rank_within_route": within,
                "ranking_metric": row["ranking_metric"],
                "ranking_value": row["ranking_value"],
                "global_rank": global_rank,
            })

    if not out_rows:
        return _empty_selection(clone_col)
    return pd.DataFrame(out_rows, columns=[clone_col, *_SELECTION_COLUMNS])
