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
]


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
