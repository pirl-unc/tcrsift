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
"""

import logging
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def match_til(
    culture_clonotypes: pd.DataFrame,
    til_data: ad.AnnData,
    match_by: str = "CDR3ab",
    min_til_cells: int = 1,
) -> pd.DataFrame:
    """
    Match culture-validated clonotypes against TIL data.

    Parameters
    ----------
    culture_clonotypes : pd.DataFrame
        Filtered clonotypes from culture experiments
    til_data : ad.AnnData
        TIL data with TCR information
    match_by : str
        Matching strategy: "CDR3ab" or "CDR3b_only"
    min_til_cells : int
        Minimum TIL cells to count as present

    Returns
    -------
    pd.DataFrame
        Culture clonotypes with TIL match information
    """
    logger.info(f"Matching {len(culture_clonotypes)} culture clonotypes against TIL data")

    df = culture_clonotypes.copy()

    # Build TIL clone lookup
    til_df = til_data.obs.copy()

    if match_by == "CDR3ab":
        # Build CDR3ab identifier in TIL data
        til_df["clone_id"] = (
            til_df.get("CDR3_alpha", pd.Series("", index=til_df.index)).fillna("") +
            "_" +
            til_df.get("CDR3_beta", pd.Series("", index=til_df.index)).fillna("")
        )
    else:
        til_df["clone_id"] = til_df.get("CDR3_beta", pd.Series("", index=til_df.index)).fillna("")

    # Count TIL cells per clone
    til_clone_counts = til_df["clone_id"].value_counts().to_dict()

    # Initialize match columns
    df["til_match"] = False
    df["til_cell_count"] = 0
    df["til_frequency"] = 0.0

    # Total TIL cells with complete TCR
    total_til = len(til_df[til_df["clone_id"] != "_"])

    # Match each culture clone
    for idx, row in df.iterrows():
        if match_by == "CDR3ab":
            clone_id = row.get("clone_id", "")
        else:
            clone_id = row.get("CDR3_beta", "")

        if clone_id and clone_id in til_clone_counts:
            til_count = til_clone_counts[clone_id]
            if til_count >= min_til_cells:
                df.loc[idx, "til_match"] = True
                df.loc[idx, "til_cell_count"] = til_count
                df.loc[idx, "til_frequency"] = til_count / total_til if total_til > 0 else 0

    n_matches = df["til_match"].sum()
    logger.info(f"Found {n_matches} culture clonotypes present in TILs ({n_matches/len(df)*100:.1f}%)")

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
    culture_clonotypes: Optional[pd.DataFrame] = None,
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
