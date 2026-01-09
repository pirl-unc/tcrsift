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
T cell phenotyping for TCRsift.

Classifies cells as CD4+ or CD8+ based on gene expression markers.
"""
from __future__ import annotations

import logging

import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)


# T cell type categories in order of confidence
TCELL_TYPE_CATEGORIES = [
    "Confident CD8+",
    "Confident CD4+",
    "Likely CD8+",
    "Likely CD4+",
    "Unknown",
]


def classify_tcell_type(
    cd4_expr: float,
    cd8_expr: float,
    cd4_cd8_ratio: float = 3.0,
) -> str:
    """
    Classify a single cell as CD4+ or CD8+ based on expression.

    Parameters
    ----------
    cd4_expr : float
        CD4 gene expression (raw counts)
    cd8_expr : float
        CD8A + CD8B combined expression (raw counts)
    cd4_cd8_ratio : float
        Ratio threshold for confident classification

    Returns
    -------
    str
        T cell type classification
    """
    # Avoid division by zero
    cd4_safe = cd4_expr + 1
    cd8_safe = cd8_expr + 1

    ratio_cd8_over_cd4 = cd8_safe / cd4_safe
    ratio_cd4_over_cd8 = cd4_safe / cd8_safe

    if ratio_cd8_over_cd4 > cd4_cd8_ratio:
        return "Confident CD8+"
    elif ratio_cd4_over_cd8 > cd4_cd8_ratio:
        return "Confident CD4+"
    elif cd8_expr > 0 and cd4_expr == 0:
        return "Likely CD8+"
    elif cd4_expr > 0 and cd8_expr == 0:
        return "Likely CD4+"
    else:
        return "Unknown"


def phenotype_cells(
    adata: ad.AnnData,
    cd4_cd8_ratio: float = 3.0,
    min_cd3_reads: int = 10,
) -> ad.AnnData:
    """
    Add T cell phenotype classification to AnnData object.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with T cell marker expression in obs (CD3, CD4, CD8)
    cd4_cd8_ratio : float
        Ratio threshold for confident CD4/CD8 classification
    min_cd3_reads : int
        Minimum CD3 reads to consider a valid T cell

    Returns
    -------
    ad.AnnData
        AnnData with Tcell_type column added to obs
    """
    logger.info(f"Phenotyping {len(adata)} cells with CD4/CD8 ratio threshold {cd4_cd8_ratio}")

    # Check for required columns
    required_cols = ["CD4", "CD8"]
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        # Try to compute CD8 from CD8A + CD8B
        if "CD8" in missing and "CD8A" in adata.obs.columns and "CD8B" in adata.obs.columns:
            adata.obs["CD8"] = adata.obs["CD8A"] + adata.obs["CD8B"]
            missing.remove("CD8")
        if missing:
            raise ValueError(f"Missing required columns for phenotyping: {missing}")

    # Classify each cell
    tcell_types = []
    for _, row in adata.obs.iterrows():
        cd4_expr = row.get("CD4", 0) or 0
        cd8_expr = row.get("CD8", 0) or 0
        tcell_type = classify_tcell_type(cd4_expr, cd8_expr, cd4_cd8_ratio)
        tcell_types.append(tcell_type)

    adata.obs["Tcell_type"] = pd.Categorical(tcell_types, categories=TCELL_TYPE_CATEGORIES)

    # Add convenience columns
    adata.obs["is_CD8"] = adata.obs["Tcell_type"].str.contains("CD8")
    adata.obs["is_CD4"] = adata.obs["Tcell_type"].str.contains("CD4")
    adata.obs["is_confident"] = adata.obs["Tcell_type"].str.startswith("Confident")

    # CD3 filter
    if "CD3" in adata.obs.columns:
        adata.obs["filter:min_cd3"] = adata.obs["CD3"] >= min_cd3_reads
    else:
        adata.obs["filter:min_cd3"] = True

    # Combined filter for confident T cells with complete TCR
    has_tcr = adata.obs.get("has_both_chains", pd.Series(True, index=adata.obs_names))
    adata.obs["confident_and_complete"] = (
        adata.obs["is_confident"]
        & has_tcr
        & adata.obs["filter:min_cd3"]
    )

    # Log summary
    type_counts = adata.obs["Tcell_type"].value_counts()
    logger.info(f"Phenotyping results:\n{type_counts.to_string()}")

    return adata


def validate_phenotype_vs_expected(
    adata: ad.AnnData,
) -> list[str]:
    """
    Compare observed phenotypes to expected based on sample metadata.

    Returns list of warning messages for mismatches.
    """
    warnings = []

    if "expected_tcell_type" not in adata.obs.columns:
        return warnings

    for sample in adata.obs["sample"].unique():
        sample_mask = adata.obs["sample"] == sample
        sample_data = adata.obs[sample_mask]

        expected = sample_data["expected_tcell_type"].iloc[0]
        if pd.isna(expected) or expected == "mixed":
            continue

        # Count observed types
        observed_cd8 = sample_data["is_CD8"].sum()
        observed_cd4 = sample_data["is_CD4"].sum()
        total = len(sample_data)

        if total == 0:
            continue

        cd8_pct = observed_cd8 / total * 100
        cd4_pct = observed_cd4 / total * 100

        if expected == "CD8" and cd4_pct > 50:
            warnings.append(
                f"Sample '{sample}': Expected CD8+ but found {cd4_pct:.1f}% CD4+ cells. "
                "This may indicate antigen cross-presentation or incorrect metadata."
            )
        elif expected == "CD4" and cd8_pct > 50:
            warnings.append(
                f"Sample '{sample}': Expected CD4+ but found {cd8_pct:.1f}% CD8+ cells. "
                "This may indicate MHC-I restricted epitopes or incorrect metadata."
            )

    return warnings


def filter_by_tcell_type(
    adata: ad.AnnData,
    tcell_type: str = "cd8",
) -> ad.AnnData:
    """
    Filter AnnData to only include cells of specified T cell type.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with Tcell_type column
    tcell_type : str
        Type to keep: "cd8", "cd4", or "both"

    Returns
    -------
    ad.AnnData
        Filtered AnnData
    """
    if "Tcell_type" not in adata.obs.columns:
        raise ValueError("AnnData must have Tcell_type column. Run phenotype_cells first.")

    if tcell_type.lower() == "cd8":
        mask = adata.obs["is_CD8"]
        logger.info(f"Filtering to CD8+ cells: {mask.sum()} of {len(adata)}")
    elif tcell_type.lower() == "cd4":
        mask = adata.obs["is_CD4"]
        logger.info(f"Filtering to CD4+ cells: {mask.sum()} of {len(adata)}")
    elif tcell_type.lower() == "both":
        mask = adata.obs["is_CD8"] | adata.obs["is_CD4"]
        logger.info(f"Keeping both CD4+ and CD8+ cells: {mask.sum()} of {len(adata)}")
    else:
        raise ValueError(f"Invalid tcell_type: {tcell_type}. Use 'cd8', 'cd4', or 'both'")

    return adata[mask].copy()


def get_phenotype_summary(adata: ad.AnnData) -> pd.DataFrame:
    """
    Get summary of T cell phenotypes by sample.

    Returns
    -------
    pd.DataFrame
        Summary with counts and percentages per sample
    """
    if "Tcell_type" not in adata.obs.columns:
        raise ValueError("AnnData must have Tcell_type column. Run phenotype_cells first.")

    summary_data = []
    for sample in adata.obs["sample"].unique():
        sample_data = adata.obs[adata.obs["sample"] == sample]
        total = len(sample_data)

        row = {
            "sample": sample,
            "total_cells": total,
        }

        for tcell_type in TCELL_TYPE_CATEGORIES:
            count = (sample_data["Tcell_type"] == tcell_type).sum()
            row[tcell_type] = count
            row[f"{tcell_type}_pct"] = count / total * 100 if total > 0 else 0

        # Add expected type if available
        if "expected_tcell_type" in sample_data.columns:
            row["expected_type"] = sample_data["expected_tcell_type"].iloc[0]

        summary_data.append(row)

    return pd.DataFrame(summary_data)
