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
Multi-experiment unification for TCRsift.

Functions for merging clonotype data from multiple experiments into a
unified table with source tracking and combined statistics.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

from .validation import TCRsiftValidationError

logger = logging.getLogger(__name__)


def merge_experiments(
    experiments: list[tuple[pd.DataFrame, str]],
    *,
    key_cols: list[str] | None = None,
    count_cols: list[str] | None = None,
    gex_cols: list[str] | None = None,
    add_occurrence_flags: bool = True,
    add_combined_stats: bool = True,
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Merge clonotype data from multiple experiments into a unified table.

    Creates a single DataFrame with prefixed columns from each experiment,
    occurrence flags (occurs_in_*), and combined statistics.

    Parameters
    ----------
    experiments : list of (DataFrame, name) tuples
        List of (clonotype_df, experiment_name) pairs
    key_cols : list of str, optional
        Columns to use as merge keys (default: CDR3_pair, CDR3_alpha, CDR3_beta)
    count_cols : list of str, optional
        Count columns to combine across sources (default: auto-detect)
    gex_cols : list of str, optional
        Gene expression columns to combine (default: auto-detect)
    add_occurrence_flags : bool
        Add 'occurs_in_*' boolean columns for each source
    add_combined_stats : bool
        Add combined statistics across all sources
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Unified clonotype table with prefixed columns and combined stats

    Examples
    --------
    >>> unified = merge_experiments([
    ...     (til_df, "TIL"),
    ...     (culture_df, "Culture"),
    ...     (amplify_df, "Amplify"),
    ... ])
    """
    if not experiments:
        raise TCRsiftValidationError(
            "No experiments provided to merge",
            hint="Provide a list of (DataFrame, name) tuples",
        )

    if key_cols is None:
        key_cols = ["CDR3_pair", "CDR3_alpha", "CDR3_beta"]

    if verbose:
        logger.info(f"Merging {len(experiments)} experiments")
        for df, name in experiments:
            logger.info(f"  {name}: {len(df):,} clonotypes")

    # Collect all unique clonotypes
    all_clonotypes = set()
    for df, name in experiments:
        if "CDR3_pair" not in df.columns:
            raise TCRsiftValidationError(
                f"Experiment '{name}' missing CDR3_pair column",
                hint="Ensure all DataFrames have standardized clonotype columns",
            )
        all_clonotypes.update(df["CDR3_pair"].dropna())

    if verbose:
        logger.info(f"  Total unique clonotypes: {len(all_clonotypes):,}")

    # Add prefixes to each DataFrame
    prefixed_dfs = []
    source_names = []
    sources_with_counts = []
    sources_with_gex = []

    for df, name in experiments:
        df_copy = df.copy()

        # Track which sources have count/gex columns
        if "total_cells.count" in df_copy.columns:
            sources_with_counts.append(name)
        if any("gex." in c and ".sum" in c for c in df_copy.columns):
            sources_with_gex.append(name)

        # Rename columns with prefix (except key columns)
        rename_dict = {
            col: f"{name}.{col}"
            for col in df_copy.columns
            if col not in key_cols
        }
        df_copy = df_copy.rename(columns=rename_dict)
        prefixed_dfs.append(df_copy)
        source_names.append(name)

    # Merge all DataFrames
    if verbose:
        logger.info("  Merging DataFrames...")

    merged = prefixed_dfs[0]
    for df in prefixed_dfs[1:]:
        merged = pd.merge(merged, df, on=key_cols, how="outer")

    if verbose:
        logger.info(f"  Merged result: {len(merged):,} rows, {len(merged.columns)} columns")

    # Add basic metadata
    merged["has_CDR3_alpha"] = merged["CDR3_alpha"].fillna("").str.len() > 0
    merged["has_CDR3_beta"] = merged["CDR3_beta"].fillna("").str.len() > 0
    merged["is_paired"] = merged["has_CDR3_alpha"] & merged["has_CDR3_beta"]

    # Add occurrence flags
    if add_occurrence_flags:
        if verbose:
            logger.info("  Adding occurrence flags...")
        for df, name in experiments:
            clone_set = set(df["CDR3_pair"].dropna())
            col_name = f"occurs_in_{name}"
            merged[col_name] = merged["CDR3_pair"].isin(clone_set)
            if verbose:
                n_true = merged[col_name].sum()
                logger.info(f"    {col_name}: {n_true:,} True")

    # Add combined statistics
    if add_combined_stats:
        if verbose:
            logger.info("  Computing combined statistics...")

        # Combine count columns
        count_col_names = count_cols or ["total_cells", "CD4_only", "CD8_only"]
        for count_col in count_col_names:
            input_cols = [f"{source}.{count_col}.count" for source in sources_with_counts]
            input_cols = [c for c in input_cols if c in merged.columns]
            if input_cols:
                output_col = f"combined.{count_col}.count"
                merged[output_col] = merged[input_cols].fillna(0).sum(axis=1)
                if verbose:
                    total = merged[output_col].sum()
                    logger.info(f"    {output_col}: {total:,.0f} total")

        # Combine GEX columns
        if sources_with_gex:
            gex_col_names = gex_cols or ["gex.CD3", "gex.CD4", "gex.CD8"]
            for gex_col in gex_col_names:
                input_cols = [f"{source}.{gex_col}.sum" for source in sources_with_gex]
                input_cols = [c for c in input_cols if c in merged.columns]
                if input_cols:
                    output_col = f"combined.{gex_col}.sum"
                    merged[output_col] = merged[input_cols].fillna(0).sum(axis=1)

        # Compute CD4/CD8 fractions
        if "combined.total_cells.count" in merged.columns:
            total_cells = merged["combined.total_cells.count"]
            for cell_type in ["CD4_only", "CD8_only"]:
                count_col = f"combined.{cell_type}.count"
                if count_col in merged.columns:
                    frac_col = f"combined.{cell_type}.frac"
                    merged[frac_col] = merged[count_col] / total_cells.replace(0, np.nan)

    return merged


def add_phenotype_confidence(
    df: pd.DataFrame,
    *,
    cd4_sum_col: str = "combined.gex.CD4.sum",
    cd8_sum_col: str = "combined.gex.CD8.sum",
    ratio_threshold: float = 10.0,
    til_evidence_cols: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Add phenotype confidence columns based on combined evidence.

    Computes Confident_CD4, Confident_CD8, Likely_CD4, Likely_CD8 based on:
    - Gene expression ratios across all experiments
    - TIL occurrence patterns
    - Culture type evidence

    Parameters
    ----------
    df : pd.DataFrame
        Unified clonotype DataFrame
    cd4_sum_col : str
        Column with combined CD4 expression
    cd8_sum_col : str
        Column with combined CD8 expression
    ratio_threshold : float
        Ratio threshold for confident classification
    til_evidence_cols : list of str, optional
        Columns indicating TIL CD4/CD8 evidence
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        DataFrame with confidence columns added
    """
    if verbose:
        logger.info("Adding phenotype confidence columns...")

    cd4_sum = df.get(cd4_sum_col, pd.Series(0, index=df.index)).fillna(0)
    cd8_sum = df.get(cd8_sum_col, pd.Series(0, index=df.index)).fillna(0)

    # Confident classifications based on expression ratio
    df["Confident_CD4"] = (
        ((cd4_sum > 0) & (cd8_sum == 0)) |
        (cd4_sum > (ratio_threshold * (1 + cd8_sum)))
    )
    df["Confident_CD8"] = (
        ((cd8_sum > 0) & (cd4_sum == 0)) |
        (cd8_sum > (ratio_threshold * (1 + cd4_sum)))
    )

    # TIL evidence
    til_cd8_mask = pd.Series(False, index=df.index)
    til_cd4_mask = pd.Series(False, index=df.index)

    if til_evidence_cols:
        for col in til_evidence_cols:
            if col in df.columns:
                if "CD8" in col:
                    til_cd8_mask |= df[col].fillna(False)
                elif "CD4" in col:
                    til_cd4_mask |= df[col].fillna(False)
    else:
        # Auto-detect TIL columns
        for col in df.columns:
            if "occurs_in_TIL" in col:
                if "CD8" in col:
                    til_cd8_mask |= df[col].fillna(False)
                elif "CD4" in col:
                    til_cd4_mask |= df[col].fillna(False)

    # Likely classifications
    df["Likely_CD8"] = df["Confident_CD8"] | til_cd8_mask
    df["Likely_CD4"] = df["Confident_CD4"] | til_cd4_mask | (
        ~df["Likely_CD8"] & (cd4_sum > 0)
    )

    if verbose:
        logger.info(f"  Confident CD4+: {df['Confident_CD4'].sum():,}")
        logger.info(f"  Confident CD8+: {df['Confident_CD8'].sum():,}")
        logger.info(f"  Likely CD4+: {df['Likely_CD4'].sum():,}")
        logger.info(f"  Likely CD8+: {df['Likely_CD8'].sum():,}")

    return df


def compute_condition_statistics(
    df: pd.DataFrame,
    conditions: list[str],
    *,
    source_prefix: str | None = None,
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Compute per-condition statistics for antigen specificity analysis.

    Parameters
    ----------
    df : pd.DataFrame
        Unified clonotype DataFrame
    conditions : list of str
        Condition names (e.g., ["pool1", "pool2", "pool3"])
    source_prefix : str, optional
        Prefix for condition columns (e.g., "culture" -> "culture.condition_pool1.count")
    verbose : bool
        Print progress
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        DataFrame with condition statistics added
    """
    if verbose:
        logger.info(f"Computing statistics for {len(conditions)} conditions")

    # Find condition columns for each condition
    for condition in conditions:
        # Find fraction columns for this condition
        frac_pattern = f"condition_{condition}" if not source_prefix else f"{source_prefix}.condition_{condition}"
        frac_cols = [c for c in df.columns if frac_pattern in c and c.endswith(".frac")]

        if not frac_cols:
            if verbose:
                logger.info(f"  No fraction columns found for {condition}")
            continue

        # Compute sum and mean of fractions across sources
        data = df[frac_cols]
        df[f"combined.{condition}.frac.sum"] = data.sum(axis=1, skipna=True)
        df[f"combined.{condition}.frac.mean"] = data.fillna(0).mean(axis=1)

        # Standard deviation and CV
        df[f"combined.{condition}.frac.std"] = data.fillna(0).std(axis=1)
        mean_vals = data.fillna(0).mean(axis=1)
        std_vals = data.fillna(0).std(axis=1)
        cv = np.where(mean_vals > 0, std_vals / mean_vals, np.nan)
        df[f"combined.{condition}.frac.cv"] = cv

        # Consistency (1 / (1 + CV))
        df[f"combined.{condition}.frac.consistency"] = np.nan_to_num(1 / (1 + cv), nan=0.0)

        if verbose:
            n_nonzero = (df[f"combined.{condition}.frac.sum"] > 0).sum()
            logger.info(f"  {condition}: {n_nonzero:,} clonotypes with evidence")

    return df


def find_top_condition(
    df: pd.DataFrame,
    conditions: list[str],
    *,
    stat: str = "frac.mean",
    require_clear_winner: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Find the top condition for each clonotype based on a statistic.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with condition statistics
    conditions : list of str
        Condition names
    stat : str
        Statistic to use for ranking (default: "frac.mean")
    require_clear_winner : bool
        Only assign top condition if it's strictly higher than second place
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        DataFrame with top_condition column added
    """
    stat_cols = [f"combined.{c}.{stat}" for c in conditions]
    stat_cols = [c for c in stat_cols if c in df.columns]

    if not stat_cols:
        if verbose:
            logger.warning(f"No columns found matching pattern 'combined.*.{stat}'")
        return df

    data = df[stat_cols].copy()

    def get_top(row):
        sorted_vals = row.sort_values(ascending=False, na_position="last")
        if len(sorted_vals.dropna()) < 1:
            return None
        if require_clear_winner and len(sorted_vals.dropna()) >= 2:
            if sorted_vals.iloc[0] <= sorted_vals.iloc[1]:
                return None
        return sorted_vals.index[0].replace("combined.", "").replace(f".{stat}", "")

    df["top_condition"] = data.apply(get_top, axis=1)
    df["top_condition.value"] = data.max(axis=1)

    if verbose:
        n_assigned = df["top_condition"].notna().sum()
        logger.info(f"  Top condition assigned to {n_assigned:,} clonotypes")
        if n_assigned > 0:
            top_counts = df["top_condition"].value_counts()
            for cond, count in top_counts.head(5).items():
                logger.info(f"    {cond}: {count:,}")

    return df


def get_unify_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Get summary statistics for a unified clonotype table.

    Returns
    -------
    dict
        Summary statistics
    """
    summary = {
        "total_clonotypes": len(df),
        "paired_clonotypes": df.get("is_paired", pd.Series(False)).sum(),
    }

    # Source occurrence counts
    occurs_cols = [c for c in df.columns if c.startswith("occurs_in_")]
    for col in occurs_cols:
        source = col.replace("occurs_in_", "")
        summary[f"from_{source}"] = df[col].sum()

    # Combined counts
    if "combined.total_cells.count" in df.columns:
        summary["total_cells"] = df["combined.total_cells.count"].sum()

    # Phenotype counts
    for col in ["Confident_CD4", "Confident_CD8", "Likely_CD4", "Likely_CD8"]:
        if col in df.columns:
            summary[col] = df[col].sum()

    return summary
