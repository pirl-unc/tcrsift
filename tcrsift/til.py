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
TIL (Tumor-Infiltrating Lymphocyte) matching for TCRsift.

Identifies culture-validated TCRs in TIL samples.
Supports loading TIL data from multiple formats:
- h5ad: Pre-processed AnnData files
- CSV: Simple CSV with CDR3_alpha, CDR3_beta columns
- vdj_dir: CellRanger VDJ output directory
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .sample_sheet import Sample, SampleSheet
from .validation import TCRsiftValidationError

logger = logging.getLogger(__name__)


def load_til_data(
    source_type: str,
    path: str | Path,
    sample_name: str | None = None,
) -> pd.DataFrame:
    """
    Load TIL data from various formats into a standard DataFrame.

    Parameters
    ----------
    source_type : str
        Type of source: "h5ad", "csv", or "vdj_dir"
    path : str or Path
        Path to the data source
    sample_name : str, optional
        Name to assign to this sample (used in output)

    Returns
    -------
    pd.DataFrame
        DataFrame with at least CDR3_alpha, CDR3_beta, and sample columns.
        Each row represents one cell.

    Raises
    ------
    TCRsiftValidationError
        If the file doesn't exist or has invalid format
    """
    path = Path(path)
    sample_name = sample_name or path.stem

    if source_type == "h5ad":
        return _load_til_from_h5ad(path, sample_name)
    elif source_type == "csv":
        return _load_til_from_csv(path, sample_name)
    elif source_type == "vdj_dir":
        return _load_til_from_vdj_dir(path, sample_name)
    else:
        raise TCRsiftValidationError(
            f"Unknown TIL source type: {source_type}",
            hint="Valid types are: h5ad, csv, vdj_dir",
        )


def _load_til_from_h5ad(path: Path, sample_name: str) -> pd.DataFrame:
    """Load TIL data from an h5ad file."""
    if not path.exists():
        raise TCRsiftValidationError(
            f"TIL h5ad file not found: {path}",
            hint="Check that the path is correct.",
        )

    adata = ad.read_h5ad(path)
    df = adata.obs.copy()

    # Ensure sample column
    if "sample" not in df.columns:
        df["sample"] = sample_name

    # Validate CDR3 columns exist
    if "CDR3_alpha" not in df.columns and "CDR3_beta" not in df.columns:
        raise TCRsiftValidationError(
            f"TIL h5ad file missing CDR3 columns: {path}",
            hint="Expected CDR3_alpha and/or CDR3_beta columns in .obs",
        )

    logger.info(f"Loaded {len(df)} TIL cells from h5ad: {path}")
    return df


def _load_til_from_csv(path: Path, sample_name: str) -> pd.DataFrame:
    """Load TIL data from a CSV file."""
    if not path.exists():
        raise TCRsiftValidationError(
            f"TIL CSV file not found: {path}",
            hint="Check that the path is correct.",
        )

    df = pd.read_csv(path)

    # Ensure sample column
    if "sample" not in df.columns:
        df["sample"] = sample_name

    # Check for CDR3 columns (try various common names)
    cdr3_alpha_candidates = ["CDR3_alpha", "cdr3_alpha", "CDR3a", "cdr3a", "TRA_cdr3", "cdr3_TRA"]
    cdr3_beta_candidates = ["CDR3_beta", "cdr3_beta", "CDR3b", "cdr3b", "TRB_cdr3", "cdr3_TRB"]

    # Find alpha column
    alpha_col = None
    for col in cdr3_alpha_candidates:
        if col in df.columns:
            alpha_col = col
            break

    # Find beta column
    beta_col = None
    for col in cdr3_beta_candidates:
        if col in df.columns:
            beta_col = col
            break

    if alpha_col is None and beta_col is None:
        raise TCRsiftValidationError(
            f"TIL CSV file missing CDR3 columns: {path}",
            hint=f"Expected CDR3_alpha and/or CDR3_beta columns. Found: {list(df.columns)}",
        )

    # Rename to standard names
    if alpha_col and alpha_col != "CDR3_alpha":
        df = df.rename(columns={alpha_col: "CDR3_alpha"})
    if beta_col and beta_col != "CDR3_beta":
        df = df.rename(columns={beta_col: "CDR3_beta"})

    logger.info(f"Loaded {len(df)} TIL cells from CSV: {path}")
    return df


def _load_til_from_vdj_dir(path: Path, sample_name: str) -> pd.DataFrame:
    """Load TIL data from a CellRanger VDJ directory."""
    # Import here to avoid circular imports
    from .loader import load_cellranger_vdj, pivot_vdj_to_cells

    if not path.exists():
        raise TCRsiftValidationError(
            f"TIL VDJ directory not found: {path}",
            hint="Check that the path is correct.",
        )

    # Load VDJ data
    vdj_df = load_cellranger_vdj(path, sample_name, verbose=False)

    # Pivot to per-cell format
    df = pivot_vdj_to_cells(vdj_df, sample_name)

    logger.info(f"Loaded {len(df)} TIL cells from VDJ directory: {path}")
    return df


def load_til_samples(
    samples: list[Sample] | SampleSheet,
) -> dict[str, pd.DataFrame]:
    """
    Load TIL data from multiple samples.

    Parameters
    ----------
    samples : list[Sample] or SampleSheet
        Samples to load (only TIL samples will be processed)

    Returns
    -------
    dict[str, pd.DataFrame]
        Dictionary mapping sample names to DataFrames
    """
    if isinstance(samples, SampleSheet):
        samples = samples.get_til_samples()
    else:
        samples = [s for s in samples if s.is_til()]

    til_data = {}
    for sample in samples:
        source = sample.get_til_data_source()
        if source is None:
            logger.warning(f"TIL sample '{sample.sample}' has no data source, skipping")
            continue

        source_type, source_path = source
        df = load_til_data(source_type, source_path, sample.sample)

        # Add tissue info if available
        if sample.tissue and "tissue" not in df.columns:
            df["tissue"] = sample.tissue

        til_data[sample.sample] = df

    return til_data


def match_til(
    culture_clonotypes: pd.DataFrame,
    til_data: ad.AnnData | pd.DataFrame | dict[str, pd.DataFrame],
    match_by: str = "CDR3ab",
    min_til_cells: int = 1,
) -> pd.DataFrame:
    """
    Match culture-validated clonotypes against TIL data.

    Supports single TIL sample (AnnData or DataFrame) or multiple TIL samples
    (dict mapping sample names to DataFrames).

    Parameters
    ----------
    culture_clonotypes : pd.DataFrame
        Filtered clonotypes from culture experiments
    til_data : AnnData, DataFrame, or dict[str, DataFrame]
        TIL data with TCR information. Can be:
        - AnnData: Single TIL sample (uses .obs)
        - DataFrame: Single TIL sample with CDR3_alpha, CDR3_beta columns
        - dict: Multiple TIL samples mapping name -> DataFrame
    match_by : str
        Matching strategy: "CDR3ab" or "CDR3b_only"
    min_til_cells : int
        Minimum TIL cells to count as present (per sample for multi-sample)

    Returns
    -------
    pd.DataFrame
        Culture clonotypes with TIL match information:
        - til_match: Found in any TIL sample
        - til_samples: Comma-separated list of matching TIL samples
        - til_cell_count: Total cells across all TIL samples
        - til_frequency: Combined frequency
        - til_cell_count.{sample}: Cells in specific sample (multi-sample only)
        - til_frequency.{sample}: Frequency in specific sample (multi-sample only)
    """
    # Normalize til_data to dict format
    if isinstance(til_data, ad.AnnData):
        til_dict = {"TIL": til_data.obs.copy()}
    elif isinstance(til_data, pd.DataFrame):
        til_dict = {"TIL": til_data.copy()}
    elif isinstance(til_data, dict):
        til_dict = {k: v.copy() for k, v in til_data.items()}
    else:
        raise TypeError(f"til_data must be AnnData, DataFrame, or dict, got {type(til_data)}")

    n_samples = len(til_dict)
    logger.info(f"Matching {len(culture_clonotypes)} culture clonotypes against {n_samples} TIL sample(s)")

    df = culture_clonotypes.copy()

    # Process each TIL sample
    sample_stats = {}
    for sample_name, til_df in til_dict.items():
        # Build clone_id
        if match_by == "CDR3ab":
            til_df["clone_id"] = (
                til_df.get("CDR3_alpha", pd.Series("", index=til_df.index)).fillna("") +
                "_" +
                til_df.get("CDR3_beta", pd.Series("", index=til_df.index)).fillna("")
            )
        else:
            til_df["clone_id"] = til_df.get("CDR3_beta", pd.Series("", index=til_df.index)).fillna("")

        # Count cells per clone
        clone_counts = til_df["clone_id"].value_counts().to_dict()
        total_til = len(til_df[til_df["clone_id"] != "_"])

        sample_stats[sample_name] = {
            "clone_counts": clone_counts,
            "total_til": total_til,
        }

    # Initialize combined columns
    df["til_match"] = False
    df["til_samples"] = ""
    df["til_cell_count"] = 0
    df["til_frequency"] = 0.0

    # Initialize per-sample columns if multiple samples
    if n_samples > 1:
        for sample_name in til_dict.keys():
            df[f"til_cell_count.{sample_name}"] = 0
            df[f"til_frequency.{sample_name}"] = 0.0

    # Total TIL cells across all samples (for combined frequency)
    total_til_all = sum(s["total_til"] for s in sample_stats.values())

    # Match each culture clone
    for idx, row in df.iterrows():
        if match_by == "CDR3ab":
            clone_id = row.get("clone_id", "")
        else:
            clone_id = row.get("CDR3_beta", "")

        if not clone_id:
            continue

        matching_samples = []
        total_count = 0

        for sample_name, stats in sample_stats.items():
            if clone_id in stats["clone_counts"]:
                count = stats["clone_counts"][clone_id]
                if count >= min_til_cells:
                    matching_samples.append(sample_name)
                    total_count += count

                    # Per-sample stats
                    if n_samples > 1:
                        df.loc[idx, f"til_cell_count.{sample_name}"] = count
                        df.loc[idx, f"til_frequency.{sample_name}"] = (
                            count / stats["total_til"] if stats["total_til"] > 0 else 0
                        )

        if matching_samples:
            df.loc[idx, "til_match"] = True
            df.loc[idx, "til_samples"] = ",".join(matching_samples)
            df.loc[idx, "til_cell_count"] = total_count
            df.loc[idx, "til_frequency"] = total_count / total_til_all if total_til_all > 0 else 0

    n_matches = df["til_match"].sum()
    logger.info(f"Found {n_matches} culture clonotypes present in TILs ({n_matches/len(df)*100:.1f}%)")

    if n_samples > 1:
        for sample_name in til_dict.keys():
            sample_matches = (df[f"til_cell_count.{sample_name}"] > 0).sum()
            logger.info(f"  - {sample_name}: {sample_matches} matches")

    return df


def get_til_enrichment(
    matched_clonotypes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculate enrichment statistics for TIL-matched clonotypes.

    Parameters
    ----------
    matched_clonotypes : pd.DataFrame
        Clonotypes with TIL match information

    Returns
    -------
    pd.DataFrame
        Enrichment statistics per clone
    """
    df = matched_clonotypes.copy()

    if "til_match" not in df.columns:
        raise ValueError("Clonotypes must have TIL match information. Run match_til first.")

    # Only calculate for matched clones
    matched = df[df["til_match"]].copy()

    if len(matched) == 0:
        logger.warning("No TIL-matched clonotypes found")
        return matched

    # Enrichment: TIL frequency vs culture frequency
    if "max_frequency" in matched.columns and "til_frequency" in matched.columns:
        matched["til_enrichment"] = np.log2(
            (matched["til_frequency"] + 1e-6) / (matched["max_frequency"] + 1e-6)
        )
    else:
        matched["til_enrichment"] = 0

    return matched


def get_til_summary(
    matched_clonotypes: pd.DataFrame,
) -> dict:
    """
    Get summary of TIL matching results.

    Returns
    -------
    dict
        Summary statistics
    """
    if "til_match" not in matched_clonotypes.columns:
        return {"error": "No TIL match information"}

    matched = matched_clonotypes[matched_clonotypes["til_match"]]

    summary = {
        "total_culture_clones": len(matched_clonotypes),
        "til_matched_clones": len(matched),
        "til_recovery_rate": len(matched) / len(matched_clonotypes) if len(matched_clonotypes) > 0 else 0,
        "total_til_cells_matched": matched["til_cell_count"].sum(),
        "median_til_frequency": matched["til_frequency"].median() if len(matched) > 0 else 0,
    }

    # By tier if available
    if "tier" in matched_clonotypes.columns:
        tier_recovery = {}
        for tier in matched_clonotypes["tier"].unique():
            if tier is not None:
                tier_df = matched_clonotypes[matched_clonotypes["tier"] == tier]
                tier_matched = tier_df["til_match"].sum()
                tier_recovery[tier] = tier_matched / len(tier_df) if len(tier_df) > 0 else 0
        summary["recovery_by_tier"] = tier_recovery

    # By antigen if available
    if "antigens" in matched.columns:
        antigen_recovery = matched.groupby("antigens")["til_cell_count"].sum().to_dict()
        summary["til_cells_by_antigen"] = antigen_recovery

    return summary


def identify_til_specific_clones(
    til_data: ad.AnnData,
    culture_clonotypes: pd.DataFrame | None = None,
    min_cells: int = 2,
) -> pd.DataFrame:
    """
    Identify clones that are abundant in TILs but not in culture.

    These could be tumor-reactive TCRs not captured in the culture system.

    Parameters
    ----------
    til_data : ad.AnnData
        TIL data
    culture_clonotypes : pd.DataFrame, optional
        Culture clonotypes to exclude
    min_cells : int
        Minimum cells in TIL to consider

    Returns
    -------
    pd.DataFrame
        TIL-specific clones
    """
    til_df = til_data.obs.copy()

    # Build clone identifier
    til_df["clone_id"] = (
        til_df.get("CDR3_alpha", pd.Series("", index=til_df.index)).fillna("") +
        "_" +
        til_df.get("CDR3_beta", pd.Series("", index=til_df.index)).fillna("")
    )

    # Aggregate TIL clones
    til_clones = til_df.groupby("clone_id").agg({
        "sample": "first",
    }).reset_index()

    til_clones["til_cell_count"] = til_df.groupby("clone_id").size().values
    til_clones = til_clones[til_clones["til_cell_count"] >= min_cells]

    # Extract CDR3 sequences
    til_clones["CDR3_alpha"] = til_clones["clone_id"].str.split("_").str[0]
    til_clones["CDR3_beta"] = til_clones["clone_id"].str.split("_").str[1]

    # Filter out culture clones if provided
    if culture_clonotypes is not None:
        culture_ids = set(culture_clonotypes["clone_id"].values)
        til_clones = til_clones[~til_clones["clone_id"].isin(culture_ids)]
        logger.info(f"Found {len(til_clones)} TIL-specific clones not in culture")
    else:
        logger.info(f"Found {len(til_clones)} expanded TIL clones")

    return til_clones
