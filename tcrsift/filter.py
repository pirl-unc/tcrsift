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
Clonotype filtering for TCRsift.

Implements tiered filtering to identify antigen-specific TCR clones.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .model import fit_frequency_logistic_model
from .validation import (
    TCRsiftValidationError,
    safe_percentage,
    validate_clonotype_df,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


# Default tier definitions for threshold-based filtering
#
# As of 1.1.0 (#99) the default tiers encode **abundance only** — the
# specificity dimension (how broadly a clone is distributed across
# antigen conditions) is no longer baked into the tier label. Mixing
# the two made the tier opaque ("tier3" could mean a borderline-cells
# clone OR a strongly-enriched clone demoted by ``max_conditions``)
# and impossible to interpret per-sample.
#
# Tier criteria (abundance):
#   - min_cells: minimum cells sharing the clonotype (global cell_count)
#   - min_frequency: minimum within-sample frequency (max_frequency)
#
# Specificity is now an explicit, separate downstream filter — use
# ``n_conditions <= K`` on the assembled frame, or pass
# ``max_conditions=K`` to :func:`filter_clonotypes_threshold`. The
# tier label answers "how strong is the enrichment?"; ``n_conditions``
# answers "how narrowly does it appear?". Keeping them orthogonal
# makes the tier interpretable per-sample (via
# :func:`per_sample_tier`) and easier to audit downstream.
#
# Back-compat: :func:`clone_clears_tier` still honors ``max_conditions``
# when present in a *user-supplied* tier dict — only the bundled
# defaults dropped the key.
#
DEFAULT_THRESHOLD_TIERS = {
    "tier1": {  # Strictest abundance criteria
        "min_cells": 10,
        "min_frequency": 0.01,  # 1%
    },
    "tier2": {
        "min_cells": 5,
        "min_frequency": 0.005,  # 0.5%
    },
    "tier3": {
        "min_cells": 3,
        "min_frequency": 0.001,  # 0.1%
    },
    "tier4": {
        "min_cells": 2,
        "min_frequency": 0.0005,  # 0.05%
    },
    "tier5": {  # Most permissive
        "min_cells": 2,
        "min_frequency": 0.0,
    },
}

# Default FDR values for tiers
DEFAULT_FDR_TIERS = [0.0001, 0.001, 0.01, 0.1, 0.15]


def clone_clears_tier(
    cells: int | float,
    max_frequency: float,
    tier_def: dict,
    *,
    n_conditions: int | None = None,
) -> bool:
    """Test whether a clone's stats clear a tier definition (#85).

    Returns True iff:
    - ``cells >= tier_def['min_cells']``
    - ``max_frequency >= tier_def['min_frequency']``
    - ``n_conditions <= tier_def['max_conditions']`` (skipped when None)

    Downstream callers were reimplementing this logic against
    :data:`DEFAULT_THRESHOLD_TIERS`. Centralizing it keeps the
    threshold definitions in lockstep across the per-method-recovery
    plot, the selection-route heatmap, and the funnel.
    """
    try:
        if float(cells) < tier_def.get("min_cells", 0):
            return False
        if float(max_frequency) < tier_def.get("min_frequency", 0.0):
            return False
    except (TypeError, ValueError):
        return False
    if n_conditions is not None:
        max_c = tier_def.get("max_conditions")
        if max_c is not None and n_conditions > max_c:
            return False
    return True


def strictest_tier_met(
    cells: int | float,
    max_frequency: float,
    tier_defs: dict[str, dict] | None = None,
    *,
    n_conditions: int | None = None,
    tier_order: tuple[str, ...] = ("tier1", "tier2", "tier3", "tier4", "tier5"),
) -> str | None:
    """Return the strictest tier the clone qualifies for, or None (#85).

    Walks ``tier_order`` from strict→permissive and returns the first
    label whose tier definition the clone clears via
    :func:`clone_clears_tier`. ``None`` when the clone fails even the
    most permissive tier.
    """
    tier_defs = tier_defs if tier_defs is not None else DEFAULT_THRESHOLD_TIERS
    for name in tier_order:
        if name not in tier_defs:
            continue
        if clone_clears_tier(
            cells, max_frequency, tier_defs[name], n_conditions=n_conditions
        ):
            return name
    return None


def per_sample_tier(
    cells: int | float,
    frequency: float,
    tier_defs: dict[str, dict] | None = None,
    *,
    tier_order: tuple[str, ...] = ("tier1", "tier2", "tier3", "tier4", "tier5"),
) -> str | None:
    """Return the abundance tier for a single (clone, sample) row (#99).

    Thin wrapper over :func:`strictest_tier_met` aimed at downstream
    report code (per-method PDFs, per-sample heatmaps) that needs to
    label a clone's enrichment within ONE sample rather than across
    the whole cohort.

    The global ``tier`` column on a clonotype frame summarises the
    clone's peak signal across all samples; downstream callers showing
    per-sample breakdowns previously had to re-implement the threshold
    table to get a per-sample equivalent. This helper centralises that
    so the tier definitions stay in one place.

    Parameters
    ----------
    cells
        Number of cells supporting this clone in this sample.
    frequency
        Within-sample frequency of this clone (cells / sample-total).
    tier_defs
        Tier definition dict; defaults to :data:`DEFAULT_THRESHOLD_TIERS`
        (abundance-only since #99). Passing a custom dict that includes
        ``max_conditions`` is supported but unused here — per-sample
        tiers ignore specificity by definition.
    tier_order
        Walk order from strictest to most permissive. Override only
        if you've defined custom tier names.

    Returns
    -------
    str | None
        Tier label (e.g. ``"tier1"``) or ``None`` when the row clears
        no tier.
    """
    return strictest_tier_met(
        cells, frequency, tier_defs=tier_defs, tier_order=tier_order
    )


# Named filter modes (#15). Pre-canned compositions of the donor/method
# knobs that match common multi-donor / multi-method study designs.
# `fdr` is the existing default and is documented here for completeness;
# the named-mode resolver below ignores it and lets normal dispatch run.
#
# `shared-high-freq` defaults: K=2, F=0.01 — calibrated against the canonical
# 2-donor x 7-method MART-1 study where the FDR-tier filter's tier 1
# happened to produce exactly this set (issue #15).
#
# `cross-donor-public` requires ≥2 donors and ≥1 method within each, with
# a frequency floor; per #15 it warns rather than errors when the cohort
# isn't multi-donor, since blocking valid analyses is worse than a flag.
VALID_FILTER_MODES = ("fdr", "shared-high-freq", "cross-donor-public")

# FDR scope (#26): controls the null distribution the FDR-tier filter
# computes against.
VALID_FDR_SCOPES = ("auto", "global", "per-donor", "per-sample")


def resolve_fdr_scope(
    requested: str,
    n_donors: int,
    donors_share_antigen: bool,
) -> str:
    """Resolve `auto` to a concrete scope based on cohort metadata.

    Auto-resolution: ``per-donor`` when ``n_donors > 1`` and the sheet
    does not set ``donors_share_antigen``, else ``global``. Explicit
    values (``global``/``per-donor``/``per-sample``) pass through
    unchanged. Raises on unknown values.
    """
    if requested not in VALID_FDR_SCOPES:
        raise TCRsiftValidationError(
            f"Invalid fdr_scope: {requested!r}",
            hint=f"Valid options: {VALID_FDR_SCOPES}",
        )
    if requested != "auto":
        return requested
    if n_donors > 1 and not donors_share_antigen:
        return "per-donor"
    return "global"


def split_clonotypes_by_donor(
    clonotypes: pd.DataFrame,
    long_df: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Partition the clonotype table into per-donor subsets.

    Each donor gets every clone present in any of that donor's samples
    (so a public clone appears under both donors). When ``long_df`` lacks
    a ``donor`` column, returns ``{"all": clonotypes}`` so callers can
    treat per-donor scope uniformly.

    Used by ``--fdr-scope per-donor`` to run the FDR filter against each
    donor's own pooled-across-methods null.
    """
    if "donor" not in long_df.columns:
        return {"all": clonotypes.copy()}
    out: dict[str, pd.DataFrame] = {}
    for donor, group in long_df.groupby("donor", observed=True):
        donor_clones = set(group["CDR3ab"].astype(str).unique())
        out[str(donor)] = clonotypes[
            clonotypes["CDR3ab"].astype(str).isin(donor_clones)
        ].copy()
    return out


def split_clonotypes_by_sample(
    clonotypes: pd.DataFrame,
    long_df: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Partition the clonotype table into per-sample subsets.

    Used by ``--fdr-scope per-sample`` — each sample gets its own FDR null.
    """
    if "sample" not in long_df.columns:
        return {"all": clonotypes.copy()}
    out: dict[str, pd.DataFrame] = {}
    for sample, group in long_df.groupby("sample", observed=True):
        sample_clones = set(group["CDR3ab"].astype(str).unique())
        out[str(sample)] = clonotypes[
            clonotypes["CDR3ab"].astype(str).isin(sample_clones)
        ].copy()
    return out

NAMED_FILTER_MODE_DEFAULTS: dict[str, dict] = {
    "shared-high-freq": {
        "min_methods_per_donor": 2,
        "min_frequency_per_method": 0.01,
    },
    "cross-donor-public": {
        "min_donors": 2,
        "min_methods_per_donor": 1,
        "min_frequency_per_method": 0.005,
    },
}


def resolve_filter_mode_kwargs(
    mode: str,
    user_kwargs: dict | None = None,
) -> dict:
    """Map a named filter mode to a kwargs dict for `filter_clonotypes`.

    Returns an empty dict for `fdr` (the default mode — no override). For
    named modes, applies the mode's preset thresholds, then layers any
    explicit user_kwargs on top so users can override individual knobs
    while keeping the rest of the preset.

    Raises TCRsiftValidationError on an unknown mode.
    """
    if mode not in VALID_FILTER_MODES:
        raise TCRsiftValidationError(
            f"Invalid filter_mode: '{mode}'",
            hint=f"Valid modes: {VALID_FILTER_MODES}",
        )
    resolved = dict(NAMED_FILTER_MODE_DEFAULTS.get(mode, {}))
    if user_kwargs:
        # Strip None / 0 / 0.0 from user_kwargs so unset CLI flags don't
        # clobber preset values.
        for k, v in user_kwargs.items():
            if v not in (None, 0, 0.0):
                resolved[k] = v
    return resolved


def _get_condition_count_info(clonotypes: pd.DataFrame) -> tuple[pd.Series | None, str | None]:
    """Return the best available condition-count series and its column name."""
    for col in ("n_conditions", "n_antigens", "n_samples"):
        if col in clonotypes.columns:
            return clonotypes[col], col
    return None, None


def _assign_logistic_tiers_from_thresholds(
    clonotypes: pd.DataFrame,
    fdr_to_threshold: dict[float, float],
) -> pd.DataFrame:
    """Assign logistic tiers from precomputed FDR thresholds."""
    df = clonotypes.copy()
    df["tier"] = None

    sorted_fdrs = sorted(fdr_to_threshold, reverse=True)  # Highest FDR (lowest tier) first

    for i, fdr in enumerate(sorted_fdrs):
        tier_name = f"tier{len(sorted_fdrs) - i}"
        threshold = fdr_to_threshold[fdr]
        mask = df["max_frequency"] >= threshold
        df.loc[mask, "tier"] = tier_name
        df.loc[mask, "fdr_threshold"] = fdr

    return df


def filter_clonotypes_threshold(
    clonotypes: pd.DataFrame,
    min_cells: int = 2,
    min_frequency: float = 0.0,
    max_conditions: int | None = None,
    require_complete: bool = True,
    tcell_type: str | None = None,
    exclude_viral: bool = False,
    min_donors: int = 0,
    min_methods_per_donor: int = 0,
    min_cells_per_method: int = 0,
    min_frequency_per_method: float = 0.0,
    min_timepoints: int = 0,
    min_timepoints_per_donor: int = 0,
    min_apcs: int = 0,
    min_apcs_per_donor: int = 0,
    min_til_cells_per_donor: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Filter clonotypes using simple threshold criteria.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    min_cells : int
        Minimum cell count per clone
    min_frequency : float
        Minimum frequency in any condition
    max_conditions : int or None
        Maximum number of conditions a clone can appear in. ``None`` (default)
        disables the filter. Uses antigen/condition counts when available,
        otherwise sample count is used as a proxy.
    require_complete : bool
        Require both alpha and beta chains
    tcell_type : str, optional
        Filter to specific T cell type ("cd8" or "cd4")
    exclude_viral : bool
        Exclude clones flagged as viral
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        Filtered clonotypes
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_filtering=True)
    validate_numeric_param(min_cells, "min_cells", min_value=0)
    validate_numeric_param(min_frequency, "min_frequency", min_value=0, max_value=1)
    if max_conditions is not None:
        validate_numeric_param(max_conditions, "max_conditions", min_value=1)

    if tcell_type is not None:
        valid_types = ["cd8", "cd4"]
        if tcell_type.lower() not in valid_types:
            raise TCRsiftValidationError(
                f"Invalid tcell_type: '{tcell_type}'",
                hint=f"Valid options are: {valid_types}, or None for no filtering",
            )

    df = clonotypes.copy()
    initial_count = len(df)

    if verbose:
        logger.info(f"Filtering {initial_count:,} clonotypes with threshold method")

    # Cell count filter
    if min_cells > 0:
        before = len(df)
        df = df[df["cell_count"] >= min_cells]
        if verbose:
            logger.info(
                f"  min_cells >= {min_cells}: {before:,} -> {len(df):,} ({before - len(df):,} removed)"
            )

    # Frequency filter
    if min_frequency > 0 and "max_frequency" in df.columns:
        before = len(df)
        df = df[df["max_frequency"] >= min_frequency]
        if verbose:
            logger.info(
                f"  min_frequency >= {min_frequency}: {before:,} -> {len(df):,} ({before - len(df):,} removed)"
            )

    # Condition specificity filter
    condition_counts, condition_count_col = _get_condition_count_info(df)
    if max_conditions is not None and condition_counts is not None:
        before = len(df)
        df = df[condition_counts <= max_conditions]
        if verbose:
            logger.info(
                f"  max_conditions <= {max_conditions} using {condition_count_col}: "
                f"{before:,} -> {len(df):,} ({before - len(df):,} removed)"
            )

    # Donor / method / timepoint / APC / TIL axis filters (#15, #9, #20).
    # Each is a no-op when its underlying column isn't on the clonotype
    # table — preserves backwards compat for designs that don't supply
    # the relevant axis. Previously fully silent; now emits a logger.warning
    # so a user passing --min-donors 2 against a sheet without patient_id
    # discovers the mismatch instead of seeing all clones pass through.
    axis_filters: list[tuple[float, str, str]] = [
        (min_donors, "n_donors", "min_donors"),
        (min_methods_per_donor, "max_methods_per_donor", "min_methods_per_donor"),
        (min_cells_per_method, "max_cells_per_method", "min_cells_per_method"),
        (min_frequency_per_method, "max_frequency_per_method", "min_frequency_per_method"),
        (min_timepoints, "n_timepoints", "min_timepoints"),
        (min_timepoints_per_donor, "max_timepoints_per_donor", "min_timepoints_per_donor"),
        (min_apcs, "n_apcs", "min_apcs"),
        (min_apcs_per_donor, "max_apcs_per_donor", "min_apcs_per_donor"),
        (min_til_cells_per_donor, "max_til_cells_per_donor", "min_til_cells_per_donor"),
    ]
    for value, column, label in axis_filters:
        if value <= 0:
            continue
        if column not in df.columns:
            logger.warning(
                "%s=%s requested but column %r is not on the clonotype "
                "table; skipping. Check that the relevant sample-sheet axis "
                "(patient_id / enrichment_method / timepoint / apc_type / "
                "TIL data with patient_id) is populated.",
                label, value, column,
            )
            continue
        before = len(df)
        df = df[df[column] >= value]
        if verbose:
            logger.info(
                f"  {label} >= {value}: {before:,} -> {len(df):,} "
                f"({before - len(df):,} removed)"
            )

    # Complete TCR filter
    if require_complete:
        before = len(df)
        has_alpha = df["CDR3_alpha"].notna() & (df["CDR3_alpha"] != "")
        has_beta = df["CDR3_beta"].notna() & (df["CDR3_beta"] != "")
        df = df[has_alpha & has_beta]
        if verbose:
            logger.info(
                f"  require_complete TCR: {before:,} -> {len(df):,} ({before - len(df):,} removed)"
            )

    # T cell type filter
    if tcell_type and "Tcell_type_consensus" in df.columns:
        before = len(df)
        if tcell_type.lower() == "cd8":
            df = df[df["Tcell_type_consensus"].str.contains("CD8", na=False)]
        elif tcell_type.lower() == "cd4":
            df = df[df["Tcell_type_consensus"].str.contains("CD4", na=False)]
        if verbose:
            logger.info(
                f"  tcell_type={tcell_type}: {before:,} -> {len(df):,} ({before - len(df):,} removed)"
            )

    # Viral exclusion
    if exclude_viral and "is_viral" in df.columns:
        before = len(df)
        n_viral = df["is_viral"].sum()
        df = df[~df["is_viral"]]
        if verbose:
            logger.info(
                f"  exclude_viral: {before:,} -> {len(df):,} ({n_viral:,} viral clones removed)"
            )

    if verbose:
        pct_retained = safe_percentage(len(df), initial_count, default=0.0)
        logger.info(
            f"  Final: {initial_count:,} -> {len(df):,} clonotypes ({pct_retained:.1f}% retained)"
        )

    if len(df) == 0:
        raise TCRsiftValidationError(
            "No clonotypes remain after filtering",
            hint=f"Try relaxing filter criteria. Current: min_cells={min_cells}, "
            f"min_frequency={min_frequency}, tcell_type={tcell_type}",
        )

    return df


def assign_tiers_threshold(
    clonotypes: pd.DataFrame,
    tier_definitions: dict | None = None,
    tcell_type: str | None = None,
    exclude_viral: bool = False,
) -> pd.DataFrame:
    """
    Assign quality tiers to clonotypes using threshold-based method.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    tier_definitions : dict, optional
        Custom tier definitions (default: DEFAULT_THRESHOLD_TIERS)
    tcell_type : str, optional
        Filter to specific T cell type
    exclude_viral : bool
        Exclude viral clones from all tiers

    Returns
    -------
    pd.DataFrame
        Clonotypes with 'tier' column added
    """
    if tier_definitions is None:
        tier_definitions = DEFAULT_THRESHOLD_TIERS

    df = clonotypes.copy()
    df["tier"] = None

    # Apply T cell type filter if specified
    if tcell_type and "Tcell_type_consensus" in df.columns:
        if tcell_type.lower() == "cd8":
            type_mask = df["Tcell_type_consensus"].str.contains("CD8", na=False)
        elif tcell_type.lower() == "cd4":
            type_mask = df["Tcell_type_consensus"].str.contains("CD4", na=False)
        else:
            type_mask = pd.Series(True, index=df.index)
    else:
        type_mask = pd.Series(True, index=df.index)

    # Apply viral exclusion if specified
    if exclude_viral and "is_viral" in df.columns:
        viral_mask = ~df["is_viral"]
    else:
        viral_mask = pd.Series(True, index=df.index)

    # Assign from most-permissive to strictest so the strictest tier a clone
    # qualifies for wins. Order by actual threshold strictness — (min_cells,
    # min_frequency) ascending — NOT by sorted tier NAME: lexical name order
    # breaks for custom tier labels (e.g. "tier10" sorts before "tier2", or
    # non-"tierN" names), silently mis-ranking the override precedence. For the
    # default tier1..tier5 this yields the identical order.
    def _strictness(name: str) -> tuple[float, float]:
        d = tier_definitions[name]
        return (d.get("min_cells", 0), d.get("min_frequency", 0.0))

    for tier_name in sorted(tier_definitions.keys(), key=_strictness):
        tier_def = tier_definitions[tier_name]

        cell_mask = df["cell_count"] >= tier_def["min_cells"]

        if "max_frequency" in df.columns:
            freq_mask = df["max_frequency"] >= tier_def["min_frequency"]
        else:
            freq_mask = pd.Series(True, index=df.index)

        # max_conditions is optional in the tier definition (#99):
        # the bundled defaults dropped it, but user-supplied tier
        # dicts may still include it. Skip the specificity check when
        # absent — that's the new abundance-only semantics.
        max_cond = tier_def.get("max_conditions")
        condition_counts, _ = _get_condition_count_info(df)
        if max_cond is not None and condition_counts is not None:
            cond_mask = condition_counts <= max_cond
        else:
            cond_mask = pd.Series(True, index=df.index)

        tier_mask = cell_mask & freq_mask & cond_mask & type_mask & viral_mask
        df.loc[tier_mask, "tier"] = tier_name

    # Log tier distribution
    tier_counts = df["tier"].value_counts()
    logger.info(f"Tier distribution:\n{tier_counts.to_string()}")

    return df


def filter_clonotypes_logistic(
    clonotypes: pd.DataFrame,
    fdr_tiers: list | None = None,
    min_freq_threshold: float = 0.09,
    default_freq_threshold: float = 0.5,
    only_avoid_viral: bool = True,
) -> pd.DataFrame:
    """
    Filter clonotypes using logistic regression model.

    This method fits a logistic model to predict clone quality based on
    frequency and assigns tiers based on FDR thresholds.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame with max_frequency column
    fdr_tiers : list, optional
        FDR values for tier assignment (default: DEFAULT_FDR_TIERS)
    min_freq_threshold : float
        Minimum frequency to consider for model fitting
    default_freq_threshold : float
        Fallback threshold if model fitting fails
    only_avoid_viral : bool
        If True, model target is non-viral; if False, target is single-culture specific

    Returns
    -------
    pd.DataFrame
        Clonotypes with tier assignments and threshold information
    """
    if fdr_tiers is None:
        fdr_tiers = DEFAULT_FDR_TIERS

    df = clonotypes.copy()

    if "max_frequency" not in df.columns:
        logger.warning("max_frequency column not found, falling back to threshold method")
        return assign_tiers_threshold(df)

    # Prepare model target
    target_above_min = (df["max_frequency"] > min_freq_threshold).values

    if only_avoid_viral and "is_viral" in df.columns:
        target = target_above_min & (~df["is_viral"]).values
    else:
        condition_counts, _ = _get_condition_count_info(df)
        if condition_counts is not None:
            # Single-condition specificity
            target = target_above_min & (condition_counts == 1).values
        else:
            target = target_above_min

    # Fit logistic regression
    try:
        model = fit_frequency_logistic_model(df["max_frequency"].values, target.astype(float))
        weight = float(model.params[0])
    except Exception as e:
        logger.warning(
            f"Model fitting failed: {e}. Falling back to default frequency threshold."
        )
        fallback_thresholds = {float(fdr): float(default_freq_threshold) for fdr in fdr_tiers}
        result_df = _assign_logistic_tiers_from_thresholds(df, fallback_thresholds)
        result_df.attrs["logistic_model_weight"] = 0.0
        result_df.attrs["fdr_to_threshold"] = fallback_thresholds
        result_df.attrs["logistic_fallback_reason"] = "fit_failed"
        return result_df

    if weight < 0:
        logger.warning(
            "Data too noisy for adaptive thresholds, falling back to default frequency threshold"
        )
        fallback_thresholds = {float(fdr): float(default_freq_threshold) for fdr in fdr_tiers}
        result_df = _assign_logistic_tiers_from_thresholds(df, fallback_thresholds)
        result_df.attrs["logistic_model_weight"] = weight
        result_df.attrs["fdr_to_threshold"] = fallback_thresholds
        result_df.attrs["logistic_fallback_reason"] = "negative_weight"
        return result_df

    # Calculate thresholds for each FDR level
    x_range = np.linspace(df["max_frequency"].min(), df["max_frequency"].max(), 10000)
    y_pred = model.predict(x_range)

    fdr_to_threshold = {}
    for fdr in fdr_tiers:
        y_target = 1.0 - fdr
        threshold_idx = np.argmin(np.abs(y_pred - y_target))
        fdr_to_threshold[fdr] = max(min_freq_threshold, x_range[threshold_idx])

    df = _assign_logistic_tiers_from_thresholds(df, fdr_to_threshold)

    # Store model info
    df.attrs["logistic_model_weight"] = weight
    df.attrs["fdr_to_threshold"] = fdr_to_threshold

    return df


def pick_strictest_populated_tier(clonotypes: pd.DataFrame) -> str:
    """Return the strictest *populated* tier label, falling back to ``"*"``.

    Walks ``tier1`` → ``tier5`` and returns the first one that has any rows.
    Unexpanded cohorts (no clone >=1% frequency) have no ``tier1`` clones, so a
    hardcoded ``tier1`` target yields an empty method-recovery table and the
    panel is silently dropped (#167). Returns ``"*"`` (all clones) when there is
    no ``tier`` column or no tier is populated.
    """
    if "tier" not in clonotypes.columns:
        return "*"
    for tier in ("tier1", "tier2", "tier3", "tier4", "tier5"):
        if (clonotypes["tier"] == tier).any():
            return tier
    return "*"


def resolve_fdr_tiers_for_method(
    method: str,
    configured_tiers: list | None,
    default_tiers: list | None,
) -> tuple[list | None, bool]:
    """Resolve which ``fdr_tiers`` to pass on, given the active tiering method.

    ``fdr_tiers`` only drives the logistic (FDR) tiering path. Under any other
    method it is inert (#170). Returns ``(tiers_to_pass, warn_inert)``:

    - logistic: pass the configured tiers through; never warn.
    - non-logistic: pass ``None`` (so it can't trigger a redundant downstream
      warning), and set ``warn_inert=True`` only when the user *customized*
      ``fdr_tiers`` away from the default — that's the footgun worth flagging;
      a default list left untouched under threshold tiering is the normal path
      and stays silent.
    """
    if method == "logistic":
        return configured_tiers, False
    customized = sorted(configured_tiers or []) != sorted(default_tiers or [])
    return None, customized


def filter_clonotypes(
    clonotypes: pd.DataFrame,
    method: str = "threshold",
    tcell_type: str = "cd8",
    min_cells: int = 2,
    min_frequency: float = 0.0,
    require_complete: bool = True,
    exclude_viral: bool = False,
    fdr_tiers: list | None = None,
    tier_definitions: dict | None = None,
    min_donors: int = 0,
    min_methods_per_donor: int = 0,
    min_cells_per_method: int = 0,
    min_frequency_per_method: float = 0.0,
    min_timepoints: int = 0,
    min_timepoints_per_donor: int = 0,
    min_apcs: int = 0,
    min_apcs_per_donor: int = 0,
    min_til_cells_per_donor: int = 0,
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Main filtering function that dispatches to appropriate method.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    method : str
        Filtering method: "threshold" or "logistic"
    tcell_type : str
        T cell type filter: "cd8", "cd4", or "both"
    min_cells : int
        Minimum cells per clone
    min_frequency : float
        Minimum frequency
    require_complete : bool
        Require complete TCR
    exclude_viral : bool
        Exclude viral clones
    fdr_tiers : list, optional
        FDR tiers for logistic method
    tier_definitions : dict, optional
        Tier definitions for threshold method
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Filtered and tiered clonotypes
    """
    # Validate method
    valid_methods = ["threshold", "logistic"]
    if method not in valid_methods:
        raise TCRsiftValidationError(
            f"Invalid filtering method: '{method}'",
            hint=f"Valid options are: {valid_methods}",
        )

    # Validate tcell_type
    valid_tcell_types = ["cd8", "cd4", "both"]
    if tcell_type.lower() not in valid_tcell_types:
        raise TCRsiftValidationError(
            f"Invalid tcell_type: '{tcell_type}'",
            hint=f"Valid options are: {valid_tcell_types}",
        )

    # fdr_tiers only drives the logistic (FDR) tiering path. Under threshold
    # tiering it is silently inert — the tiers come from fixed abundance cutoffs
    # (DEFAULT_THRESHOLD_TIERS). Warn rather than ignore so a config that sets
    # fdr_tiers + method='threshold' isn't a silent footgun (#170).
    if fdr_tiers is not None and method != "logistic":
        logger.warning(
            "fdr_tiers was provided but method=%r — fdr_tiers only applies to "
            "method='logistic'. It is being ignored; tiers come from fixed "
            "abundance thresholds. Set method='logistic' for FDR-based tiering.",
            method,
        )

    logger.info(f"Filtering clonotypes using {method} method")

    # Basic filtering first
    df = filter_clonotypes_threshold(
        clonotypes,
        min_cells=min_cells,
        min_frequency=min_frequency,
        require_complete=require_complete,
        tcell_type=tcell_type if tcell_type != "both" else None,
        exclude_viral=exclude_viral,
        min_donors=min_donors,
        min_methods_per_donor=min_methods_per_donor,
        min_cells_per_method=min_cells_per_method,
        min_frequency_per_method=min_frequency_per_method,
        min_timepoints=min_timepoints,
        min_timepoints_per_donor=min_timepoints_per_donor,
        min_apcs=min_apcs,
        min_apcs_per_donor=min_apcs_per_donor,
        min_til_cells_per_donor=min_til_cells_per_donor,
        verbose=verbose,
    )

    # Assign tiers
    if verbose:
        logger.info("Assigning confidence tiers...")

    if method == "logistic":
        df = filter_clonotypes_logistic(df, fdr_tiers=fdr_tiers)
    else:
        df = assign_tiers_threshold(
            df,
            tier_definitions=tier_definitions,
            tcell_type=tcell_type if tcell_type != "both" else None,
            exclude_viral=exclude_viral,
        )

    # Log tier distribution
    if verbose and "tier" in df.columns:
        tier_counts = df["tier"].value_counts().sort_index()
        logger.info("  Tier distribution:")
        for tier, count in tier_counts.items():
            if tier is not None:
                pct = count / len(df) * 100
                logger.info(f"    {tier}: {count:,} ({pct:.1f}%)")

    return df


def split_by_tier(clonotypes: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Split clonotypes DataFrame by tier.

    Returns
    -------
    dict
        Mapping from tier name to DataFrame
    """
    if "tier" not in clonotypes.columns:
        raise ValueError("Clonotypes must have 'tier' column. Run filter_clonotypes first.")

    result = {}
    for tier in clonotypes["tier"].unique():
        if tier is not None:
            result[tier] = clonotypes[clonotypes["tier"] == tier].copy()

    return result


def bucket_by_donor_sharing(
    clonotypes: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Split clonotypes into private-to-donor and public-across-donors buckets.

    "private_to_donor" = ``n_donors == 1`` (or unset / NaN treated as 1).
    "public_across_donors" = ``n_donors >= 2``.

    Returns an empty dict if the ``n_donors`` column isn't on the table; the
    caller should treat that as "no bucketing applicable for this design"
    rather than as an error. Empty buckets are dropped from the result so
    downstream callers can skip writing zero-row CSVs for irrelevant
    designs (single-donor cohorts, etc).
    """
    if "n_donors" not in clonotypes.columns:
        return {}

    # to_numeric first so an object-dtype column doesn't trigger the
    # fillna object-downcast FutureWarning; coerce non-numeric -> NaN -> 1.
    n_donors = (
        pd.to_numeric(clonotypes["n_donors"], errors="coerce").fillna(1).astype(int)
    )
    buckets = {
        "private_to_donor": clonotypes[n_donors <= 1].copy(),
        "public_across_donors": clonotypes[n_donors >= 2].copy(),
    }
    return {name: df for name, df in buckets.items() if len(df) > 0}


def get_filter_summary(clonotypes: pd.DataFrame) -> dict:
    """
    Get summary of filtering results.

    Returns
    -------
    dict
        Summary statistics by tier
    """
    if "tier" not in clonotypes.columns:
        return {"total_clonotypes": len(clonotypes)}

    summary = {
        "total_clonotypes": len(clonotypes),
        "tier_counts": clonotypes["tier"].value_counts().to_dict(),
    }

    for tier in clonotypes["tier"].unique():
        if tier is not None:
            tier_df = clonotypes[clonotypes["tier"] == tier]
            summary[f"{tier}_cells"] = tier_df["cell_count"].sum()
            summary[f"{tier}_median_freq"] = (
                tier_df["max_frequency"].median() if "max_frequency" in tier_df.columns else None
            )

    return summary
