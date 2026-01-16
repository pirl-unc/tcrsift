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
Data loading functions for TCRsift.

Handles loading CellRanger VDJ and GEX outputs into unified data structures.
"""
from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from tqdm.auto import tqdm

from .sample_sheet import Sample, load_sample_sheet
from .validation import (
    TCRsiftValidationError,
    validate_cellranger_gex_dir,
    validate_cellranger_vdj_dir,
    validate_file_exists,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


# Standard column names for VDJ annotations
VDJ_COLUMNS = {
    "barcode": "barcode",
    "contig_id": "contig_id",
    "chain": "chain",
    "v_gene": "v_gene",
    "d_gene": "d_gene",
    "j_gene": "j_gene",
    "c_gene": "c_gene",
    "cdr3": "cdr3",
    "cdr3_nt": "cdr3_nt",
    "reads": "reads",
    "umis": "umis",
    "productive": "productive",
    "full_length": "full_length",
}

# VDJ segment columns for full sequence assembly
VDJ_SEGMENT_COLS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
VDJ_SEGMENT_NT_COLS = [c + "_nt" for c in VDJ_SEGMENT_COLS]


def load_cellranger_vdj(
    vdj_dir: str | Path,
    sample_name: str,
    annotations_filename: str = "filtered_contig_annotations.csv",
    clonotypes_filename: str = "clonotypes.csv",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load CellRanger VDJ output files.

    Parameters
    ----------
    vdj_dir : str or Path
        Path to CellRanger VDJ output directory
    sample_name : str
        Name to assign to this sample
    annotations_filename : str
        Name of the contig annotations CSV file
    clonotypes_filename : str
        Name of the clonotypes CSV file
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        DataFrame with VDJ annotations for all cells
    """
    vdj_dir = Path(vdj_dir)

    # Validate directory exists and has expected files
    try:
        vdj_dir = validate_cellranger_vdj_dir(vdj_dir)
    except TCRsiftValidationError:
        # Re-raise with more context
        raise TCRsiftValidationError(
            f"Invalid CellRanger VDJ directory for sample '{sample_name}': {vdj_dir}",
            hint="Make sure this is the 'outs' directory from 'cellranger vdj'. "
            "It should contain 'filtered_contig_annotations.csv' or 'all_contig_annotations.csv'.",
        )

    annotations_path = vdj_dir / annotations_filename
    clonotypes_path = vdj_dir / clonotypes_filename

    if not annotations_path.exists():
        # Try alternative filename
        alt_path = vdj_dir / "all_contig_annotations.csv"
        if alt_path.exists():
            annotations_path = alt_path
            if verbose:
                logger.info("  Using all_contig_annotations.csv (filtered not found)")
        else:
            raise TCRsiftValidationError(
                f"VDJ annotations file not found: {annotations_path}",
                hint=f"Expected one of: {annotations_filename}, all_contig_annotations.csv "
                f"in directory: {vdj_dir}",
            )

    logger.info(f"Loading VDJ annotations from {annotations_path}")
    df = pd.read_csv(annotations_path)

    # Validate the VDJ data
    if len(df) == 0:
        raise TCRsiftValidationError(
            f"VDJ annotations file is empty: {annotations_path}",
            hint="Check that CellRanger VDJ ran successfully. "
            "The file should contain contig annotations.",
        )

    required_cols = ["barcode", "chain"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise TCRsiftValidationError(
            f"VDJ annotations missing required columns: {missing_cols}",
            hint=f"Available columns: {list(df.columns)[:15]}. "
            "This doesn't look like a CellRanger VDJ output file.",
        )

    # Log summary statistics
    n_contigs = len(df)
    n_cells = df["barcode"].nunique()
    n_productive = df["productive"].sum() if "productive" in df.columns else "unknown"
    if verbose:
        logger.info(f"  Loaded {n_contigs:,} contigs from {n_cells:,} cells ({n_productive} productive)")

    # Validate chain types
    valid_chains = {"TRA", "TRB", "TRD", "TRG", "IGH", "IGK", "IGL", "Multi"}
    invalid_chains = set(df["chain"].unique()) - valid_chains
    if invalid_chains:
        logger.warning(f"  Unexpected chain types found: {invalid_chains}")

    # Add sample information
    df["sample"] = sample_name
    df["vdj_dir"] = str(vdj_dir)

    # Load clonotypes if available
    if clonotypes_path.exists():
        logger.info(f"Loading clonotypes from {clonotypes_path}")
        df_clonotypes = pd.read_csv(clonotypes_path)
        # Clonotypes contain MAIT/NKT evidence
        if "mait_evidence" in df_clonotypes.columns or "inkt_evidence" in df_clonotypes.columns:
            # Merge clonotype info
            clonotype_cols = ["clonotype_id"]
            if "mait_evidence" in df_clonotypes.columns:
                clonotype_cols.append("mait_evidence")
            if "inkt_evidence" in df_clonotypes.columns:
                clonotype_cols.append("inkt_evidence")
            df_clonotypes_subset = df_clonotypes[clonotype_cols].copy()
            df_clonotypes_subset = df_clonotypes_subset.rename(columns={"clonotype_id": "raw_clonotype_id"})
            df = df.merge(df_clonotypes_subset, on="raw_clonotype_id", how="left")

    # Combine VDJ segments into full sequence if available
    if all(col in df.columns for col in VDJ_SEGMENT_COLS):
        df["vdj_aa"] = df[VDJ_SEGMENT_COLS].fillna("").agg("".join, axis=1)
    if all(col in df.columns for col in VDJ_SEGMENT_NT_COLS):
        df["vdj_nt"] = df[VDJ_SEGMENT_NT_COLS].fillna("").agg("".join, axis=1)

    return df


def load_cellranger_gex(
    gex_dir: str | Path,
    sample_name: str,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Load CellRanger gene expression output.

    Parameters
    ----------
    gex_dir : str or Path
        Path to CellRanger count output directory
    sample_name : str
        Name to assign to this sample
    min_genes : int
        Minimum genes detected per cell
    max_genes : int
        Maximum genes detected per cell
    min_counts : int
        Minimum UMI counts per cell
    max_counts : int
        Maximum UMI counts per cell
    max_mito_pct : float
        Maximum mitochondrial percentage
    min_mito_pct : float
        Minimum mitochondrial percentage
    verbose : bool
        Print progress information

    Returns
    -------
    ad.AnnData
        AnnData object with gene expression data
    """
    # Validate numeric parameters
    validate_numeric_param(min_genes, "min_genes", min_value=0)
    validate_numeric_param(max_genes, "max_genes", min_value=1)
    validate_numeric_param(min_counts, "min_counts", min_value=0)
    validate_numeric_param(max_counts, "max_counts", min_value=1)
    validate_numeric_param(min_mito_pct, "min_mito_pct", min_value=0, max_value=100)
    validate_numeric_param(max_mito_pct, "max_mito_pct", min_value=0, max_value=100)

    if min_genes > max_genes:
        raise TCRsiftValidationError(
            f"min_genes ({min_genes}) cannot be greater than max_genes ({max_genes})",
            hint="Check your QC filter parameters.",
        )
    if min_counts > max_counts:
        raise TCRsiftValidationError(
            f"min_counts ({min_counts}) cannot be greater than max_counts ({max_counts})",
            hint="Check your QC filter parameters.",
        )
    if min_mito_pct > max_mito_pct:
        raise TCRsiftValidationError(
            f"min_mito_pct ({min_mito_pct}) cannot be greater than max_mito_pct ({max_mito_pct})",
            hint="Check your QC filter parameters.",
        )

    gex_dir = Path(gex_dir)

    # Validate directory
    try:
        gex_dir = validate_cellranger_gex_dir(gex_dir)
    except TCRsiftValidationError:
        raise TCRsiftValidationError(
            f"Invalid CellRanger GEX directory for sample '{sample_name}': {gex_dir}",
            hint="Make sure this is the 'outs' directory from 'cellranger count'. "
            "It should contain 'filtered_feature_bc_matrix' or 'filtered_feature_bc_matrix.h5'.",
        )

    # Try standard CellRanger output locations
    matrix_dir = gex_dir / "filtered_feature_bc_matrix"
    if not matrix_dir.exists():
        matrix_dir = gex_dir / "outs" / "filtered_feature_bc_matrix"
    if not matrix_dir.exists():
        # Try h5 file
        h5_path = gex_dir / "filtered_feature_bc_matrix.h5"
        if not h5_path.exists():
            h5_path = gex_dir / "outs" / "filtered_feature_bc_matrix.h5"
        if h5_path.exists():
            logger.info(f"Loading GEX from h5 file: {h5_path}")
            adata = sc.read_10x_h5(str(h5_path))
        else:
            available = [f.name for f in gex_dir.iterdir()][:15]
            raise TCRsiftValidationError(
                f"Gene expression matrix not found in: {gex_dir}",
                hint=f"Expected 'filtered_feature_bc_matrix' directory or 'filtered_feature_bc_matrix.h5'. "
                f"Available files/directories: {available}",
            )
    else:
        logger.info(f"Loading GEX from matrix directory: {matrix_dir}")
        adata = sc.read_10x_mtx(str(matrix_dir), var_names="gene_ids")

    # Validate loaded data
    if adata.n_obs == 0:
        raise TCRsiftValidationError(
            f"Gene expression matrix contains no cells: {gex_dir}",
            hint="Check that CellRanger count ran successfully.",
        )

    if verbose:
        logger.info(f"  Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Add sample information
    adata.obs["sample"] = sample_name
    adata.obs["gex_dir"] = str(gex_dir)

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.contains("^ENSG.*MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Rename QC columns for consistency
    if "pct_counts_mt" in adata.obs.columns:
        adata.obs["percent_mt"] = adata.obs["pct_counts_mt"]
    if "n_genes_by_counts" in adata.obs.columns:
        adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    if "total_counts" in adata.obs.columns:
        adata.obs["n_counts"] = adata.obs["total_counts"]

    # Add QC filter flags
    adata.obs["filter:min_genes"] = adata.obs["n_genes"] >= min_genes
    adata.obs["filter:max_genes"] = adata.obs["n_genes"] <= max_genes
    adata.obs["filter:min_counts"] = adata.obs["n_counts"] >= min_counts
    adata.obs["filter:max_counts"] = adata.obs["n_counts"] <= max_counts
    adata.obs["filter:min_mito"] = adata.obs["percent_mt"] >= min_mito_pct
    adata.obs["filter:max_mito"] = adata.obs["percent_mt"] <= max_mito_pct
    adata.obs["filter:pass_qc"] = (
        adata.obs["filter:min_genes"]
        & adata.obs["filter:max_genes"]
        & adata.obs["filter:min_counts"]
        & adata.obs["filter:max_counts"]
        & adata.obs["filter:min_mito"]
        & adata.obs["filter:max_mito"]
    )

    return adata


def _extract_tcell_markers(adata: ad.AnnData) -> pd.DataFrame:
    """
    Extract T cell marker gene expression from AnnData.

    Returns DataFrame with CD3, CD4, CD8 expression per cell.
    """
    # Gene name mappings (ENSEMBL IDs to gene symbols)
    gene_mappings = {
        "ENSG00000167286": "CD3D",
        "ENSG00000198851": "CD3E",
        "ENSG00000160654": "CD3G",
        "ENSG00000010610": "CD4",
        "ENSG00000153563": "CD8A",
        "ENSG00000172116": "CD8B",
    }

    markers = ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"]
    marker_df = pd.DataFrame(index=adata.obs_names)

    for marker in markers:
        # Try direct gene name match
        if marker in adata.var_names:
            marker_df[marker] = adata[:, marker].X.toarray().flatten()
        else:
            # Try ENSEMBL ID match
            found = False
            for ensembl_id, gene_name in gene_mappings.items():
                if gene_name == marker:
                    matching_vars = [v for v in adata.var_names if v.startswith(ensembl_id)]
                    if matching_vars:
                        marker_df[marker] = adata[:, matching_vars[0]].X.toarray().flatten()
                        found = True
                        break
            if not found:
                marker_df[marker] = 0
                logger.warning(f"T cell marker {marker} not found in gene expression data")

    return marker_df


def combine_gex_and_vdj(
    adata: ad.AnnData,
    vdj_df: pd.DataFrame,
    sample_name: str,
) -> ad.AnnData:
    """
    Combine gene expression and VDJ data for a single sample.

    Parameters
    ----------
    adata : ad.AnnData
        Gene expression data
    vdj_df : pd.DataFrame
        VDJ annotations
    sample_name : str
        Sample name

    Returns
    -------
    ad.AnnData
        Combined AnnData with VDJ info in obs
    """
    # Extract T cell markers
    marker_df = _extract_tcell_markers(adata)
    for col in marker_df.columns:
        adata.obs[col] = marker_df[col].values

    # Calculate combined markers
    adata.obs["CD3"] = adata.obs["CD3D"] + adata.obs["CD3E"] + adata.obs["CD3G"]
    adata.obs["CD8"] = adata.obs["CD8A"] + adata.obs["CD8B"]

    # Pivot VDJ data to get one row per barcode with chain info
    if len(vdj_df) > 0:
        vdj_pivoted = _pivot_vdj_by_barcode(vdj_df)

        # Join on barcode
        # Make sure barcodes match format (sometimes have -1 suffix)
        adata_barcodes = set(adata.obs_names)
        vdj_barcodes = set(vdj_pivoted.index)

        # Try to match barcodes
        if not adata_barcodes.intersection(vdj_barcodes):
            # Try stripping -1 suffix from GEX barcodes
            adata.obs_names = [b.split("-")[0] for b in adata.obs_names]

        # Add VDJ columns to adata.obs
        for col in vdj_pivoted.columns:
            adata.obs[col] = vdj_pivoted[col].reindex(adata.obs_names).values

    return adata


def _pivot_vdj_by_barcode(vdj_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot VDJ data to get one row per barcode with TRA/TRB chain info.

    Handles doublets by keeping track of multiple chains per barcode.
    """
    # Sort by UMI count to prioritize high-quality chains
    vdj_df = vdj_df.sort_values(["barcode", "chain", "umis", "reads"], ascending=[True, True, False, False])

    # Create entry ID for each chain per barcode (1, 2, etc.)
    vdj_df["entry_id"] = vdj_df.groupby(["barcode", "chain"]).cumcount() + 1

    # Columns to pivot - include all CDR/FWR segments if available
    pivot_cols = ["cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "umis", "reads", "contig_id"]

    # Add VDJ segment columns (fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4)
    # Skip columns already in pivot_cols to avoid duplicates
    for seg_col in VDJ_SEGMENT_COLS:
        if seg_col in vdj_df.columns and seg_col not in pivot_cols:
            pivot_cols.append(seg_col)
    for seg_col in VDJ_SEGMENT_NT_COLS:
        if seg_col in vdj_df.columns and seg_col not in pivot_cols:
            pivot_cols.append(seg_col)

    # Add combined VDJ sequences
    if "vdj_aa" in vdj_df.columns:
        pivot_cols.append("vdj_aa")
    if "vdj_nt" in vdj_df.columns:
        pivot_cols.append("vdj_nt")

    # Filter to only include first 2 entries per chain (to handle doublets)
    vdj_df = vdj_df[vdj_df["entry_id"] <= 2]

    # Pivot
    pivot_df = vdj_df.pivot_table(
        index="barcode",
        columns=["chain", "entry_id"],
        values=pivot_cols,
        aggfunc="first",
    )

    # Flatten column names
    pivot_df.columns = [f"{chain}_{entry}_{col}" for col, chain, entry in pivot_df.columns]

    # Add chain count and doublet flags
    for chain in ["TRA", "TRB"]:
        umi_cols = [c for c in pivot_df.columns if c.startswith(f"{chain}_") and c.endswith("_umis")]
        pivot_df[f"{chain}_count"] = pivot_df[umi_cols].notna().sum(axis=1)
        pivot_df[f"has_{chain}"] = pivot_df[f"{chain}_count"] > 0
        pivot_df[f"multi_{chain}"] = pivot_df[f"{chain}_count"] > 1

    pivot_df["multi_chain"] = pivot_df["multi_TRA"] | pivot_df["multi_TRB"]
    pivot_df["has_both_chains"] = pivot_df["has_TRA"] & pivot_df["has_TRB"]

    # Create combined CDR3ab identifier
    pivot_df["CDR3_alpha"] = pivot_df.get("TRA_1_cdr3", pd.Series(index=pivot_df.index))
    pivot_df["CDR3_beta"] = pivot_df.get("TRB_1_cdr3", pd.Series(index=pivot_df.index))
    pivot_df["CDR3ab"] = pivot_df["CDR3_alpha"].fillna("") + "_" + pivot_df["CDR3_beta"].fillna("")

    return pivot_df


def load_sample(
    sample: Sample,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
) -> ad.AnnData:
    """
    Load all data for a single sample.

    Parameters
    ----------
    sample : Sample
        Sample object with paths and metadata
    min_genes, max_genes, min_counts, max_counts, max_mito_pct, min_mito_pct
        QC filter parameters for GEX data

    Returns
    -------
    ad.AnnData
        Combined AnnData with GEX and VDJ data
    """
    adata = None
    vdj_df = None

    # Load GEX if available
    if sample.gex_dir:
        adata = load_cellranger_gex(
            sample.gex_dir,
            sample.sample,
            min_genes=min_genes,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            max_mito_pct=max_mito_pct,
            min_mito_pct=min_mito_pct,
        )

    # Load VDJ if available
    if sample.vdj_dir:
        vdj_df = load_cellranger_vdj(sample.vdj_dir, sample.sample)

    # Combine or create from VDJ only
    if adata is not None and vdj_df is not None:
        adata = combine_gex_and_vdj(adata, vdj_df, sample.sample)
    elif vdj_df is not None:
        # Create minimal AnnData from VDJ data
        vdj_pivoted = _pivot_vdj_by_barcode(vdj_df)
        adata = ad.AnnData(obs=vdj_pivoted)
        adata.obs["sample"] = sample.sample

    # Add sample metadata
    if adata is not None:
        adata.obs["antigen_type"] = sample.antigen_type
        adata.obs["antigen_description"] = sample.antigen_description
        adata.obs["source"] = sample.source
        adata.obs["expected_tcell_type"] = sample.get_expected_tcell_type()

    return adata


def load_samples(
    sample_sheet_path: str | Path,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
    verbose: bool = True,
    show_progress: bool = True,
) -> ad.AnnData:
    """
    Load all samples from a sample sheet into a single AnnData object.

    Parameters
    ----------
    sample_sheet_path : str or Path
        Path to sample sheet (CSV or YAML)
    min_genes, max_genes, min_counts, max_counts, max_mito_pct, min_mito_pct
        QC filter parameters
    verbose : bool
        Print detailed progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    ad.AnnData
        Combined AnnData with all samples
    """
    # Validate sample sheet path
    sample_sheet_path = validate_file_exists(sample_sheet_path, "sample sheet")

    sample_sheet = load_sample_sheet(sample_sheet_path)

    if len(sample_sheet) == 0:
        raise TCRsiftValidationError(
            f"Sample sheet is empty: {sample_sheet_path}",
            hint="Add sample entries to the sample sheet.",
        )

    logger.info(f"Loading {len(sample_sheet)} samples from {sample_sheet_path}")

    # Pre-validate all sample paths to fail fast
    validation_errors = []
    for i, sample in enumerate(sample_sheet):
        if sample.vdj_dir and not Path(sample.vdj_dir).exists():
            validation_errors.append(f"Sample '{sample.sample}': VDJ directory not found: {sample.vdj_dir}")
        if sample.gex_dir and not Path(sample.gex_dir).exists():
            validation_errors.append(f"Sample '{sample.sample}': GEX directory not found: {sample.gex_dir}")

    if validation_errors:
        raise TCRsiftValidationError(
            f"Sample sheet validation failed with {len(validation_errors)} error(s):\n" +
            "\n".join(f"  - {e}" for e in validation_errors[:5]),
            hint="Check that all paths in the sample sheet are correct and accessible.",
        )

    adatas = []
    total_cells = 0

    # Create progress bar iterator
    sample_iter = sample_sheet
    if show_progress:
        sample_iter = tqdm(
            sample_sheet,
            desc="Loading samples",
            unit="sample",
            disable=not show_progress,
        )

    for sample in sample_iter:
        if show_progress:
            sample_iter.set_postfix(sample=sample.sample[:20])

        if verbose:
            logger.info(f"Loading sample: {sample.sample}")

        try:
            adata = load_sample(
                sample,
                min_genes=min_genes,
                max_genes=max_genes,
                min_counts=min_counts,
                max_counts=max_counts,
                max_mito_pct=max_mito_pct,
                min_mito_pct=min_mito_pct,
            )
            if adata is not None:
                adatas.append(adata)
                total_cells += adata.n_obs
                if verbose:
                    logger.info(f"  Sample {sample.sample}: {adata.n_obs:,} cells")
        except TCRsiftValidationError:
            raise
        except Exception as e:
            raise TCRsiftValidationError(
                f"Failed to load sample '{sample.sample}': {e}",
                hint="Check that the CellRanger output directories are valid and complete. "
                f"VDJ: {sample.vdj_dir}, GEX: {sample.gex_dir}",
            ) from e

    if not adatas:
        raise TCRsiftValidationError(
            "No samples loaded successfully",
            hint="Check that at least one sample has valid VDJ or GEX data.",
        )

    # Concatenate all samples
    logger.info(f"Concatenating {len(adatas)} samples ({total_cells:,} total cells)")

    if show_progress:
        # Show progress for concatenation (can be slow for large datasets)
        with tqdm(total=1, desc="Concatenating samples", unit="step") as pbar:
            combined = ad.concat(
                adatas,
                join="outer",
                label="sample",
                keys=[a.obs["sample"].iloc[0] for a in adatas],
            )
            pbar.update(1)
    else:
        combined = ad.concat(
            adatas,
            join="outer",
            label="sample",
            keys=[a.obs["sample"].iloc[0] for a in adatas],
        )

    # Store sample sheet as uns
    combined.uns["sample_sheet"] = sample_sheet.to_dataframe().to_dict()

    logger.info(f"Successfully loaded {combined.n_obs:,} cells from {len(adatas)} samples")

    return combined
