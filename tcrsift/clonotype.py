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
Clonotype aggregation for TCRsift.

Groups cells by TCR CDR3 sequences to identify clonal populations.
"""

import logging

import anndata as ad
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from .validation import (
    TCRsiftValidationError,
    validate_anndata,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


def aggregate_clonotypes(
    adata: ad.AnnData,
    group_by: str = "CDR3ab",
    min_umi: int = 2,
    handle_doublets: str = "flag",
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Aggregate cells into clonotypes based on CDR3 sequences.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with VDJ and phenotype information
    group_by : str
        How to group clones: "CDR3ab" (alpha+beta) or "CDR3b_only" (beta only)
    min_umi : int
        Minimum UMI count for a chain to be considered
    handle_doublets : str
        How to handle cells with multiple chains: "flag", "remove", "keep-primary"
    verbose : bool
        Print detailed progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per unique clonotype
    """
    # Validate inputs
    adata = validate_anndata(adata, "input AnnData", min_cells=1)
    validate_numeric_param(min_umi, "min_umi", min_value=0)

    valid_group_by = ["CDR3ab", "CDR3b_only"]
    if group_by not in valid_group_by:
        raise TCRsiftValidationError(
            f"Invalid group_by: '{group_by}'",
            hint=f"Valid options are: {valid_group_by}",
        )

    valid_doublet_handling = ["flag", "remove", "keep-primary"]
    if handle_doublets not in valid_doublet_handling:
        raise TCRsiftValidationError(
            f"Invalid handle_doublets: '{handle_doublets}'",
            hint=f"Valid options are: {valid_doublet_handling}",
        )

    # Check for required columns
    if group_by == "CDR3ab":
        required = ["CDR3_alpha", "CDR3_beta"]
    else:
        required = ["CDR3_beta"]

    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        available = [c for c in adata.obs.columns if "CDR3" in c or "cdr3" in c.lower()]
        raise TCRsiftValidationError(
            f"Missing required CDR3 columns for {group_by} grouping: {missing}",
            hint=f"Available CDR3-related columns: {available}. "
            "Make sure VDJ data was loaded correctly.",
        )

    logger.info(f"Aggregating clonotypes by {group_by} from {len(adata):,} cells")

    df = adata.obs.copy()

    # Handle doublets
    if "multi_chain" in df.columns:
        n_doublets = df["multi_chain"].sum()
        if n_doublets > 0:
            if verbose:
                logger.info(f"  Found {n_doublets:,} cells with multiple chains ({n_doublets/len(df)*100:.1f}%)")

            if handle_doublets == "remove":
                df = df[~df["multi_chain"]]
                if verbose:
                    logger.info(f"  Removed doublets, {len(df):,} cells remaining")
            elif handle_doublets == "flag":
                # Keep but flag
                df["is_doublet"] = df["multi_chain"]
                if verbose:
                    logger.info("  Flagging doublets (keeping all cells)")
            # keep-primary: use primary chain (already sorted by UMI)
            elif verbose:
                logger.info("  Using primary chain for doublets")

    # Apply UMI filter
    if "TRA_1_umis" in df.columns and min_umi > 0:
        df["TRA_pass_umi"] = df["TRA_1_umis"].fillna(0) >= min_umi
        df["TRB_pass_umi"] = df["TRB_1_umis"].fillna(0) >= min_umi
    else:
        df["TRA_pass_umi"] = True
        df["TRB_pass_umi"] = True

    # Build clone identifier
    if group_by == "CDR3ab":
        # Require both chains for complete clone
        df["clone_id"] = df["CDR3_alpha"].fillna("") + "_" + df["CDR3_beta"].fillna("")
        df["is_complete_clone"] = (
            df["CDR3_alpha"].notna()
            & (df["CDR3_alpha"] != "")
            & df["CDR3_beta"].notna()
            & (df["CDR3_beta"] != "")
            & df["TRA_pass_umi"]
            & df["TRB_pass_umi"]
        )
    elif group_by == "CDR3b_only":
        df["clone_id"] = df["CDR3_beta"].fillna("")
        df["is_complete_clone"] = (
            df["CDR3_beta"].notna()
            & (df["CDR3_beta"] != "")
            & df["TRB_pass_umi"]
        )
    else:
        raise ValueError(f"Invalid group_by: {group_by}. Use 'CDR3ab' or 'CDR3b_only'")

    # Filter to complete clones
    df_complete = df[df["is_complete_clone"]].copy()
    if verbose:
        logger.info(f"  Found {len(df_complete):,} cells with complete TCR ({len(df_complete)/len(df)*100:.1f}%)")

    if len(df_complete) == 0:
        raise TCRsiftValidationError(
            "No complete clones found after filtering",
            hint=f"Check that cells have valid CDR3 sequences. "
            f"Grouping by: {group_by}, min_umi: {min_umi}. "
            f"Total cells: {len(df)}, cells with CDR3_beta: {df['CDR3_beta'].notna().sum() if 'CDR3_beta' in df.columns else 'N/A'}",
        )

    # Count unique clones for progress bar
    n_unique_clones = df_complete["clone_id"].nunique()
    if verbose:
        logger.info(f"  Aggregating {n_unique_clones:,} unique clone IDs...")

    # Aggregate by clone
    clonotypes = _aggregate_clone_data(df_complete, group_by, show_progress=show_progress)

    if verbose:
        logger.info(f"  Found {len(clonotypes):,} unique clonotypes")
        # Log clone size distribution
        if len(clonotypes) > 0:
            n_singletons = (clonotypes["cell_count"] == 1).sum()
            n_expanded = (clonotypes["cell_count"] > 1).sum()
            max_size = clonotypes["cell_count"].max()
            logger.info(f"    Singletons: {n_singletons:,}, Expanded clones: {n_expanded:,}, Max clone size: {max_size:,}")

    return clonotypes


def _aggregate_clone_data(df: pd.DataFrame, group_by: str, show_progress: bool = True) -> pd.DataFrame:
    """Aggregate cell-level data to clone-level."""

    clone_data = []

    # Create iterator with optional progress bar
    grouped = df.groupby("clone_id")
    if show_progress:
        grouped = tqdm(
            grouped,
            desc="Aggregating clones",
            unit="clone",
            total=df["clone_id"].nunique(),
        )

    for clone_id, clone_df in grouped:
        if clone_id == "_" or clone_id == "":
            continue

        record = {
            "clone_id": clone_id,
            "cell_count": len(clone_df),
            "cell_barcodes": ";".join(clone_df.index.tolist()),
        }

        # CDR3 sequences
        if group_by == "CDR3ab":
            parts = clone_id.split("_")
            record["CDR3_alpha"] = parts[0] if len(parts) > 0 else ""
            record["CDR3_beta"] = parts[1] if len(parts) > 1 else ""
        else:
            record["CDR3_alpha"] = clone_df["CDR3_alpha"].mode().iloc[0] if "CDR3_alpha" in clone_df.columns else ""
            record["CDR3_beta"] = clone_id

        # Sample and condition information
        record["samples"] = ";".join(clone_df["sample"].unique())
        record["n_samples"] = clone_df["sample"].nunique()

        # Antigen information if available
        if "antigen_description" in clone_df.columns:
            antigens = clone_df["antigen_description"].dropna().unique()
            record["antigens"] = ";".join(str(a) for a in antigens)
            record["n_antigens"] = len(antigens)

        # Source information
        if "source" in clone_df.columns:
            record["sources"] = ";".join(clone_df["source"].unique())

        # T cell type consensus
        if "Tcell_type" in clone_df.columns:
            type_counts = clone_df["Tcell_type"].value_counts()
            record["Tcell_type_consensus"] = type_counts.index[0]
            record["Tcell_type_purity"] = type_counts.iloc[0] / len(clone_df)

            # CD4/CD8 counts
            record["n_CD8"] = clone_df["is_CD8"].sum() if "is_CD8" in clone_df.columns else 0
            record["n_CD4"] = clone_df["is_CD4"].sum() if "is_CD4" in clone_df.columns else 0

        # Gene usage
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            for gene in ["v_gene", "j_gene", "c_gene"]:
                col = f"{prefix}_{gene}"
                if col in clone_df.columns:
                    mode_val = clone_df[col].mode()
                    record[f"{chain}_{gene}"] = mode_val.iloc[0] if len(mode_val) > 0 else None

        # VDJ sequences if available
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            vdj_col = f"{prefix}_vdj_aa"
            if vdj_col in clone_df.columns:
                mode_val = clone_df[vdj_col].mode()
                record[f"VDJ_{chain}_aa"] = mode_val.iloc[0] if len(mode_val) > 0 else None

            vdj_nt_col = f"{prefix}_vdj_nt"
            if vdj_nt_col in clone_df.columns:
                mode_val = clone_df[vdj_nt_col].mode()
                record[f"VDJ_{chain}_nt"] = mode_val.iloc[0] if len(mode_val) > 0 else None

        # Quality metrics
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            for metric in ["umis", "reads"]:
                col = f"{prefix}_{metric}"
                if col in clone_df.columns:
                    record[f"{chain}_{metric}_mean"] = clone_df[col].mean()
                    record[f"{chain}_{metric}_sum"] = clone_df[col].sum()

        # Contig IDs
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            contig_col = f"{prefix}_contig_id"
            if contig_col in clone_df.columns:
                record[f"{chain}_contig_ids"] = ";".join(clone_df[contig_col].dropna().astype(str).tolist())

        # Doublet information
        if "is_doublet" in clone_df.columns:
            record["n_doublet_cells"] = clone_df["is_doublet"].sum()

        # Calculate clone frequency within each sample
        sample_freqs = []
        for sample in clone_df["sample"].unique():
            sample_total = df[df["sample"] == sample]["is_complete_clone"].sum()
            sample_count = (clone_df["sample"] == sample).sum()
            if sample_total > 0:
                freq = sample_count / sample_total
                sample_freqs.append(freq)
        record["max_frequency"] = max(sample_freqs) if sample_freqs else 0
        record["mean_frequency"] = np.mean(sample_freqs) if sample_freqs else 0

        clone_data.append(record)

    return pd.DataFrame(clone_data)


def calculate_clone_frequencies(
    clonotypes: pd.DataFrame,
    adata: ad.AnnData,
) -> pd.DataFrame:
    """
    Calculate detailed frequency information for each clone.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame from aggregate_clonotypes
    adata : ad.AnnData
        Original AnnData with cell-level data

    Returns
    -------
    pd.DataFrame
        Clonotypes with additional frequency columns
    """
    df = adata.obs.copy()

    # Calculate total complete TCRs per sample
    sample_totals = df.groupby("sample")["is_complete_clone"].sum().to_dict()

    freq_data = []
    for _, clone_row in clonotypes.iterrows():
        clone_id = clone_row["clone_id"]
        clone_cells = df[df["clone_id"] == clone_id]

        sample_freqs = {}
        for sample in clone_cells["sample"].unique():
            sample_count = (clone_cells["sample"] == sample).sum()
            sample_total = sample_totals.get(sample, 1)
            sample_freqs[sample] = sample_count / sample_total if sample_total > 0 else 0

        freq_data.append({
            "clone_id": clone_id,
            "sample_frequencies": sample_freqs,
            "max_frequency": max(sample_freqs.values()) if sample_freqs else 0,
            "n_conditions_present": len(sample_freqs),
        })

    freq_df = pd.DataFrame(freq_data)

    # Merge back
    clonotypes = clonotypes.merge(
        freq_df[["clone_id", "sample_frequencies", "n_conditions_present"]],
        on="clone_id",
        how="left",
        suffixes=("", "_new"),
    )

    return clonotypes


def get_clonotype_summary(clonotypes: pd.DataFrame) -> dict:
    """
    Get summary statistics for clonotypes.

    Returns
    -------
    dict
        Summary statistics
    """
    return {
        "n_clonotypes": len(clonotypes),
        "n_cells": clonotypes["cell_count"].sum(),
        "median_clone_size": clonotypes["cell_count"].median(),
        "max_clone_size": clonotypes["cell_count"].max(),
        "n_singletons": (clonotypes["cell_count"] == 1).sum(),
        "n_expanded": (clonotypes["cell_count"] > 1).sum(),
        "n_multi_sample": (clonotypes["n_samples"] > 1).sum() if "n_samples" in clonotypes.columns else 0,
    }


def export_clonotypes_airr(clonotypes: pd.DataFrame, output_path: str):
    """
    Export clonotypes in AIRR format.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    output_path : str
        Output file path (.tsv)
    """
    # Map to AIRR standard columns
    airr_mapping = {
        "clone_id": "clone_id",
        "CDR3_alpha": "junction_aa_tra",
        "CDR3_beta": "junction_aa_trb",
        "alpha_v_gene": "v_call_tra",
        "alpha_j_gene": "j_call_tra",
        "beta_v_gene": "v_call_trb",
        "beta_d_gene": "d_call_trb",
        "beta_j_gene": "j_call_trb",
        "cell_count": "clone_count",
    }

    airr_df = pd.DataFrame()
    for src_col, dst_col in airr_mapping.items():
        if src_col in clonotypes.columns:
            airr_df[dst_col] = clonotypes[src_col]

    # Add required AIRR fields with defaults
    if "sequence_id" not in airr_df.columns:
        airr_df["sequence_id"] = [f"clone_{i}" for i in range(len(airr_df))]
    if "productive" not in airr_df.columns:
        airr_df["productive"] = "T"

    airr_df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Exported {len(airr_df)} clonotypes to AIRR format: {output_path}")
