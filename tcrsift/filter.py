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

import logging
from typing import Optional

import numpy as np
import pandas as pd

from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


# Default tier definitions for threshold-based filtering
DEFAULT_THRESHOLD_TIERS = {
    "tier1": {  # Highest confidence
        "min_cells": 10,
        "min_frequency": 0.01,  # 1%
        "max_conditions": 2,
    },
    "tier2": {
        "min_cells": 5,
        "min_frequency": 0.005,  # 0.5%
        "max_conditions": 3,
    },
    "tier3": {
        "min_cells": 3,
        "min_frequency": 0.001,  # 0.1%
        "max_conditions": 5,
    },
    "tier4": {
        "min_cells": 2,
        "min_frequency": 0.0005,  # 0.05%
        "max_conditions": 10,
    },
    "tier5": {  # Lowest confidence
        "min_cells": 2,
        "min_frequency": 0.0,
        "max_conditions": 999,
    },
}

# Default FDR values for tiers
DEFAULT_FDR_TIERS = [0.0001, 0.001, 0.01, 0.1, 0.15]


def filter_clonotypes_threshold(
    clonotypes: pd.DataFrame,
    min_cells: int = 2,
    min_frequency: float = 0.0,
    max_conditions: int = 999,
    require_complete: bool = True,
    tcell_type: Optional[str] = None,
    exclude_viral: bool = False,
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
    max_conditions : int
        Maximum number of conditions clone can appear in
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
            logger.info(f"  min_cells >= {min_cells}: {before:,} -> {len(df):,} ({before - len(df):,} removed)")

    # Frequency filter
    if min_frequency > 0 and "max_frequency" in df.columns:
        before = len(df)
        df = df[df["max_frequency"] >= min_frequency]
        if verbose:
            logger.info(f"  min_frequency >= {min_frequency}: {before:,} -> {len(df):,} ({before - len(df):,} removed)")

    # Condition specificity filter
    if max_conditions < 999 and "n_samples" in df.columns:
        before = len(df)
        df = df[df["n_samples"] <= max_conditions]
        if verbose:
            logger.info(f"  max_conditions <= {max_conditions}: {before:,} -> {len(df):,} ({before - len(df):,} removed)")

    # Complete TCR filter
    if require_complete:
        before = len(df)
        has_alpha = df["CDR3_alpha"].notna() & (df["CDR3_alpha"] != "")
        has_beta = df["CDR3_beta"].notna() & (df["CDR3_beta"] != "")
        df = df[has_alpha & has_beta]
        if verbose:
            logger.info(f"  require_complete TCR: {before:,} -> {len(df):,} ({before - len(df):,} removed)")

    # T cell type filter
    if tcell_type and "Tcell_type_consensus" in df.columns:
        before = len(df)
        if tcell_type.lower() == "cd8":
            df = df[df["Tcell_type_consensus"].str.contains("CD8", na=False)]
        elif tcell_type.lower() == "cd4":
            df = df[df["Tcell_type_consensus"].str.contains("CD4", na=False)]
        if verbose:
            logger.info(f"  tcell_type={tcell_type}: {before:,} -> {len(df):,} ({before - len(df):,} removed)")

    # Viral exclusion
    if exclude_viral and "is_viral" in df.columns:
        before = len(df)
        n_viral = df["is_viral"].sum()
        df = df[~df["is_viral"]]
        if verbose:
            logger.info(f"  exclude_viral: {before:,} -> {len(df):,} ({n_viral:,} viral clones removed)")

    if verbose:
        logger.info(f"  Final: {initial_count:,} -> {len(df):,} clonotypes ({len(df)/initial_count*100:.1f}% retained)")

    if len(df) == 0:
        raise TCRsiftValidationError(
            "No clonotypes remain after filtering",
            hint=f"Try relaxing filter criteria. Current: min_cells={min_cells}, "
            f"min_frequency={min_frequency}, tcell_type={tcell_type}",
        )

    return df


def assign_tiers_threshold(
    clonotypes: pd.DataFrame,
    tier_definitions: Optional[dict] = None,
    tcell_type: Optional[str] = None,
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

    # Assign tiers from highest to lowest (so higher tiers override)
    for tier_name in reversed(sorted(tier_definitions.keys())):
        tier_def = tier_definitions[tier_name]

        cell_mask = df["cell_count"] >= tier_def["min_cells"]

        if "max_frequency" in df.columns:
            freq_mask = df["max_frequency"] >= tier_def["min_frequency"]
        else:
            freq_mask = pd.Series(True, index=df.index)

        if "n_samples" in df.columns:
            cond_mask = df["n_samples"] <= tier_def["max_conditions"]
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
    fdr_tiers: Optional[list] = None,
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
    try:
        import statsmodels.api as sm
    except ImportError:
        logger.warning("statsmodels not installed, falling back to threshold method")
        return assign_tiers_threshold(clonotypes)

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
    elif "n_samples" in df.columns:
        # Single-culture specificity
        target = target_above_min & (df["n_samples"] == 1).values
    else:
        target = target_above_min

    # Fit logistic regression
    try:
        model = sm.Logit(target.astype(float), df["max_frequency"].values)
        result = model.fit(disp=False)
        weight = result.params[0]
    except Exception as e:
        logger.warning(f"Model fitting failed: {e}. Using default thresholds.")
        df["tier"] = "tier5"
        df.loc[df["max_frequency"] >= default_freq_threshold, "tier"] = "tier1"
        return df

    if weight < 0:
        logger.warning("Data too noisy for adaptive thresholds, using defaults")
        df["tier"] = "tier5"
        df.loc[df["max_frequency"] >= default_freq_threshold, "tier"] = "tier1"
        return df

    # Calculate thresholds for each FDR level
    x_range = np.linspace(df["max_frequency"].min(), df["max_frequency"].max(), 10000)
    y_pred = result.predict(x_range)

    fdr_to_threshold = {}
    for fdr in fdr_tiers:
        y_target = 1.0 - fdr
        threshold_idx = np.argmin(np.abs(y_pred - y_target))
        fdr_to_threshold[fdr] = max(min_freq_threshold, x_range[threshold_idx])

    # Assign tiers based on thresholds
    df["tier"] = None
    sorted_fdrs = sorted(fdr_tiers, reverse=True)  # Highest FDR (lowest tier) first

    for i, fdr in enumerate(sorted_fdrs):
        tier_name = f"tier{len(sorted_fdrs) - i}"
        threshold = fdr_to_threshold[fdr]
        df.loc[df["max_frequency"] >= threshold, "tier"] = tier_name
        df.loc[df["max_frequency"] >= threshold, "fdr_threshold"] = fdr

    # Store model info
    df.attrs["logistic_model_weight"] = weight
    df.attrs["fdr_to_threshold"] = fdr_to_threshold

    return df


def filter_clonotypes(
    clonotypes: pd.DataFrame,
    method: str = "threshold",
    tcell_type: str = "cd8",
    min_cells: int = 2,
    min_frequency: float = 0.0,
    require_complete: bool = True,
    exclude_viral: bool = False,
    fdr_tiers: Optional[list] = None,
    tier_definitions: Optional[dict] = None,
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

    logger.info(f"Filtering clonotypes using {method} method")

    # Basic filtering first
    df = filter_clonotypes_threshold(
        clonotypes,
        min_cells=min_cells,
        min_frequency=min_frequency,
        require_complete=require_complete,
        tcell_type=tcell_type if tcell_type != "both" else None,
        exclude_viral=exclude_viral,
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
            summary[f"{tier}_median_freq"] = tier_df["max_frequency"].median() if "max_frequency" in tier_df.columns else None

    return summary
