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
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from anndata.experimental import concat_on_disk
from scipy.sparse import csr_matrix, issparse
from tqdm.auto import tqdm

from .genes import TCELL_MARKERS, find_column_for_gene
from .sample_sheet import Sample, SampleSheet, load_sample_sheet
from .validation import (
    TCRsiftValidationError,
    validate_cellranger_gex_dir,
    validate_cellranger_vdj_dir,
    validate_file_exists,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


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
    if "productive" in df.columns:
        productive = df["productive"]
        if productive.dtype == bool:
            n_productive = int(productive.sum())
        else:
            n_productive = int((productive.astype(str) == "True").sum())
    else:
        n_productive = "unknown"
    if verbose:
        logger.info(
            f"  Loaded {n_contigs:,} contigs from {n_cells:,} cells ({n_productive} productive)"
        )

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
            df_clonotypes_subset = df_clonotypes_subset.rename(
                columns={"clonotype_id": "raw_clonotype_id"}
            )
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
    min_mito_pct: float = 0.0,
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
        Minimum mitochondrial percentage (a FLOOR: cells below this are
        dropped). Default 0 = no floor; non-standard, opt-in only (#168).
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

    # Calculate QC metrics - detect mitochondrial genes
    # Handle both loading paths:
    # - h5 path: var_names are gene symbols (MT-ND1, etc.)
    # - mtx path with var_names="gene_ids": var_names are ENSEMBL IDs, symbols in var['gene_symbols']
    if "gene_symbols" in adata.var.columns:
        adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    else:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
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

    # The min_mito_pct FLOOR discards *low*-mito cells, the inverse of standard
    # scRNA-seq QC (where low mito is good). Default is 0 (no floor); when a
    # positive floor is in effect, report its drop explicitly so it can't
    # silently shrink the cohort unnoticed (#168). Loud (warning) regardless of
    # verbose, since a floor is an unusual choice.
    if min_mito_pct > 0:
        n_below_floor = int((~adata.obs["filter:min_mito"]).sum())
        if n_below_floor:
            logger.warning(
                f"  min_mito_pct={min_mito_pct} is a FLOOR: dropping "
                f"{n_below_floor:,}/{adata.n_obs:,} cells with <{min_mito_pct}% "
                "mitochondrial content (low mito is usually good — set "
                "min_mito_pct=0 to keep them)."
            )

    adata.obs["filter:pass_qc"] = (
        adata.obs["filter:min_genes"]
        & adata.obs["filter:max_genes"]
        & adata.obs["filter:min_counts"]
        & adata.obs["filter:max_counts"]
        & adata.obs["filter:min_mito"]
        & adata.obs["filter:max_mito"]
    )

    # Apply the QC mask. Prior to issue #39 these flags were advisory only —
    # the column was computed but cells failing QC were not dropped, so the
    # user-facing --min-mito/--max-mito/--min-genes/etc. parameters were
    # silently no-ops.
    n_before = adata.n_obs
    adata = adata[adata.obs["filter:pass_qc"]].copy()
    if verbose:
        logger.info(
            f"  QC: {n_before:,} -> {adata.n_obs:,} cells pass "
            f"(dropped {n_before - adata.n_obs:,})"
        )
    # Catch "silently wrong" config mistakes (e.g. --min-mito 1.0 typo for
    # 10) where the filter eats most of the sample and downstream artifacts
    # just look weird rather than failing clearly (#41). 50% is a reasonable
    # line — peripheral blood routinely passes 90%+, tumor sits around 70%,
    # below half is almost always a config or dataset-mismatch issue.
    pass_rate = adata.n_obs / n_before if n_before else 0.0
    if pass_rate < 0.5:
        logger.warning(
            f"  QC dropped {(1 - pass_rate) * 100:.0f}% of cells "
            f"({n_before - adata.n_obs:,}/{n_before:,}) — "
            "check --min-genes/--max-genes/--min-counts/--max-counts/"
            "--min-mito/--max-mito if this is unexpected."
        )

    return adata


def _extract_tcell_markers(adata: ad.AnnData) -> pd.DataFrame:
    """
    Extract T cell marker gene expression from AnnData.

    Returns DataFrame with CD3, CD4, CD8 expression per cell.
    Uses version-robust ENSEMBL ID matching from genes.py.

    Resolves all present markers via a single column slice on `adata`
    rather than N per-marker slices — one anndata-view creation and one
    sparse-to-dense conversion instead of N (issue #35, finding #2).
    """
    var_names = list(adata.var_names)

    present_cols: list[str] = []
    present_symbols: list[str] = []
    missing_symbols: list[str] = []
    for gene in TCELL_MARKERS:
        col = find_column_for_gene(gene, var_names)
        if col is not None:
            present_cols.append(col)
            present_symbols.append(gene.symbol)
        else:
            missing_symbols.append(gene.symbol)

    for sym in missing_symbols:
        logger.warning(f"T cell marker {sym} not found in gene expression data")

    if present_cols:
        block = adata[:, present_cols].X
        if hasattr(block, "toarray"):
            block = block.toarray()
        # asarray is defensive: covers the h5py-dataset path under backed
        # mode where .X is neither sparse nor an ndarray.
        block = np.asarray(block)
        marker_df = pd.DataFrame(
            block, index=adata.obs_names, columns=present_symbols
        )
        # Drop the local handle regardless of whether pandas shares or
        # copies its backing. If shared, harmless (marker_df still holds
        # the array); if copied (e.g. later re-blocked when we add
        # missing-marker columns below), this frees the source so peak
        # memory stays at ~n_obs*n_present, not 2x.
        del block
    else:
        marker_df = pd.DataFrame(index=adata.obs_names)
    for sym in missing_symbols:
        marker_df[sym] = 0

    return marker_df


def combine_gex_and_vdj(
    adata: ad.AnnData,
    vdj_df: pd.DataFrame,
) -> ad.AnnData:
    """
    Combine gene expression and VDJ data for a single sample.

    Parameters
    ----------
    adata : ad.AnnData
        Gene expression data
    vdj_df : pd.DataFrame
        VDJ annotations

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

        # Build barcode mapping (handles gem group suffixes like -1, -2, -3)
        # Map each GEX barcode to matching VDJ barcode without modifying adata
        vdj_barcodes = set(vdj_pivoted.index)

        # Helper to strip gem group suffix (e.g., -1, -2, -3)
        def _strip_suffix(bc: str) -> str:
            return bc.rsplit("-", 1)[0] if "-" in bc else bc

        # Start with direct matches where possible
        barcode_to_vdj = {bc: bc for bc in adata.obs_names if bc in vdj_barcodes}

        # Build base-barcode mapping for VDJ barcodes
        vdj_by_base: dict[str, list[str]] = {}
        for bc in vdj_barcodes:
            base = _strip_suffix(bc)
            vdj_by_base.setdefault(base, []).append(bc)

        # Try base matching for any unmapped GEX barcodes
        for bc in adata.obs_names:
            if bc in barcode_to_vdj:
                continue
            bc_base = _strip_suffix(bc)
            if bc_base in vdj_by_base:
                candidates = vdj_by_base[bc_base]
                if bc in candidates:
                    barcode_to_vdj[bc] = bc
                else:
                    barcode_to_vdj[bc] = candidates[0]

        # Add VDJ columns to adata.obs using the mapping
        mapped_barcodes = [barcode_to_vdj.get(bc) for bc in adata.obs_names]
        for col in vdj_pivoted.columns:
            src = vdj_pivoted[col]
            reindexed = src.reindex(mapped_barcodes)
            values = reindexed.values
            # reindex upcasts bool/int to object when missing barcodes introduce NaN,
            # which then breaks h5ad serialization. Restore the original dtype with
            # the natural "absent" fill value.
            missing = pd.isna(values)
            if src.dtype == bool:
                values = np.where(missing, False, values).astype(bool)
            elif pd.api.types.is_integer_dtype(src.dtype):
                values = np.where(missing, 0, values).astype(src.dtype)
            adata.obs[col] = values

    return adata


def _pivot_vdj_by_barcode(vdj_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot VDJ data to get one row per barcode with TRA/TRB chain info.

    Handles doublets by keeping track of multiple chains per barcode.
    """
    # Sort by UMI count to prioritize high-quality chains
    vdj_df = vdj_df.sort_values(
        ["barcode", "chain", "umis", "reads"], ascending=[True, True, False, False]
    )

    # Create entry ID for each chain per barcode (1, 2, etc.)
    vdj_df["entry_id"] = vdj_df.groupby(
        ["barcode", "chain"], observed=True,
    ).cumcount() + 1

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
        umi_cols = [
            c for c in pivot_df.columns if c.startswith(f"{chain}_") and c.endswith("_umis")
        ]
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
    min_mito_pct: float = 0.0,
) -> ad.AnnData | None:
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
    ad.AnnData or None
        Combined AnnData with GEX and VDJ data, or None if sample has neither
        gex_dir nor vdj_dir.
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
        adata = combine_gex_and_vdj(adata, vdj_df)
    elif vdj_df is not None:
        # Create minimal AnnData from VDJ data
        vdj_pivoted = _pivot_vdj_by_barcode(vdj_df)
        adata = ad.AnnData(obs=vdj_pivoted)
        adata.obs["sample"] = sample.sample

    # Add sample metadata. Skip fields the user didn't supply: an all-None
    # object column survives concat but breaks anndata.write_h5ad, and
    # downstream readers already guard with `if col in adata.obs.columns`.
    if adata is not None:
        metadata = [
            ("antigen_type", sample.antigen_type),
            ("antigen_description", sample.antigen_description),
            ("antigen_name", sample.antigen_name),
            ("antigen_sequence", sample.antigen_sequence),
            ("epitope_sequence", sample.epitope_sequence),
            ("mhc_allele", sample.mhc_allele),
            ("antigen_names", sample.antigen_names),
            ("antigen_sequences", sample.antigen_sequences),
            ("epitope_sequences", sample.epitope_sequences),
            ("source", sample.source),
            ("patient_id", sample.patient_id),
            ("enrichment_method", sample.enrichment_method),
            ("timepoint", sample.timepoint),
            ("apc_type", sample.apc_type),
            ("tissue", sample.tissue),
            ("expected_tcell_type", sample.get_expected_tcell_type()),
        ]
        for col, val in metadata:
            if val is not None:
                adata.obs[col] = val

    return adata


def _ensure_x(adata: ad.AnnData) -> None:
    """Normalize `adata.X` so the spilled h5ad merges via `concat_on_disk`.

    `concat_on_disk` requires every input store to have an X group, and
    along the obs axis only supports CSR sparse or dense. Two cases:

    1. VDJ-only AnnData has no expression matrix at all. Synthesize a
       zero-column CSR placeholder. The placeholder is stripped post-merge
       in `load_samples` when `combined.n_vars == 0`, preserving the
       contract that VDJ-only loads have `X is None`.

    2. Any non-CSR sparse format raises `NotImplementedError: Concat of
       following not supported: [...]` from `concat_on_disk`. Convert to
       CSR. The original report was CellRanger's CSC output (#37), but
       the same shape applies to COO, `csc_array` (scipy 1.11+ sparse
       arrays), etc.
    """
    if adata.X is None:
        adata.X = csr_matrix((adata.n_obs, 0))
    elif issparse(adata.X) and getattr(adata.X, "format", "csr") != "csr":
        adata.X = adata.X.tocsr()


def load_samples(
    sample_sheet_path: str | Path | SampleSheet,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 0.0,
    verbose: bool = True,
    show_progress: bool = True,
    tmpdir: str | Path | None = None,
) -> ad.AnnData:
    """
    Load all samples from a sample sheet into a single AnnData object.

    Parameters
    ----------
    sample_sheet_path : str or Path or SampleSheet
        Path to sample sheet (CSV or YAML), or a SampleSheet instance.
    min_genes, max_genes, min_counts, max_counts, max_mito_pct, min_mito_pct
        QC filter parameters.
    verbose : bool
        Print detailed progress information.
    show_progress : bool
        Show progress bar.
    tmpdir : str or Path or None
        Parent directory for the spill tempdir (see Notes). Defaults to the
        system temp location ($TMPDIR / /tmp). Pass an explicit disk-backed
        directory when the system temp is on tmpfs.

    Returns
    -------
    ad.AnnData
        Combined AnnData with all samples. `X` is None when every input was
        VDJ-only, otherwise it carries the outer-joined expression matrix.

    Notes
    -----
    Memory: each per-sample AnnData is spilled to a tempfile h5ad after
    load and the in-memory copy is freed, then merged with
    `anndata.experimental.concat_on_disk`, which streams inputs and output
    in bounded chunks (~400 MB default). The merged file is read back into
    memory once. Peak memory ≈ max(one sample, output), versus the prior
    in-memory `ad.concat` peak of ~2 × Σ(samples).

    Temp disk: spilled per-sample h5ads plus the merged output, so ~2× the
    total sparse dataset size. Redirect with `tmpdir=` when $TMPDIR is on
    tmpfs.

    VDJ-only samples: see `_ensure_x` — a zero-column placeholder is
    synthesized to keep the merge path uniform, then stripped after the
    merge so `combined.X` is None for VDJ-only loads.
    """
    # Load sample sheet (path or object)
    if isinstance(sample_sheet_path, SampleSheet):
        sample_sheet = sample_sheet_path
        sample_sheet_label = "<SampleSheet>"
    else:
        sample_sheet_path = validate_file_exists(sample_sheet_path, "sample sheet")
        sample_sheet = load_sample_sheet(sample_sheet_path)
        sample_sheet_label = str(sample_sheet_path)

    if len(sample_sheet) == 0:
        raise TCRsiftValidationError(
            f"Sample sheet is empty: {sample_sheet_label}",
            hint="Add sample entries to the sample sheet.",
        )

    logger.info(f"Loading {len(sample_sheet)} samples from {sample_sheet_label}")

    # Pre-validate all sample paths to fail fast
    validation_errors = []
    for sample in sample_sheet:
        if sample.vdj_dir and not Path(sample.vdj_dir).exists():
            validation_errors.append(
                f"Sample '{sample.sample}': VDJ directory not found: {sample.vdj_dir}"
            )
        if sample.gex_dir and not Path(sample.gex_dir).exists():
            validation_errors.append(
                f"Sample '{sample.sample}': GEX directory not found: {sample.gex_dir}"
            )

    if validation_errors:
        raise TCRsiftValidationError(
            f"Sample sheet validation failed with {len(validation_errors)} error(s):\n"
            + "\n".join(f"  - {e}" for e in validation_errors[:5]),
            hint="Check that all paths in the sample sheet are correct and accessible.",
        )

    sample_keys: list[str] = []
    spill_paths: list[Path] = []
    total_cells = 0

    with tempfile.TemporaryDirectory(prefix="tcrsift_load_", dir=tmpdir) as spill_dir_str:
        spill_dir = Path(spill_dir_str)
        sample_iter = tqdm(
            sample_sheet,
            desc="Loading samples",
            unit="sample",
            disable=not show_progress,
        )

        for sample in sample_iter:
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
            except TCRsiftValidationError:
                raise
            except Exception as e:
                raise TCRsiftValidationError(
                    f"Failed to load sample '{sample.sample}': {e}",
                    hint="Check that the CellRanger output directories are valid and complete. "
                    f"VDJ: {sample.vdj_dir}, GEX: {sample.gex_dir}",
                ) from e

            if adata is None:
                continue

            _ensure_x(adata)
            sample_keys.append(adata.obs["sample"].iloc[0])
            total_cells += adata.n_obs

            # Number files by the load order of successful samples (skipped
            # samples don't leave gaps), so the spill dir is easy to inspect.
            path = spill_dir / f"sample_{len(spill_paths):04d}.h5ad"
            adata.write_h5ad(path)
            spill_paths.append(path)

            if verbose:
                logger.info(f"  Sample {sample.sample}: {adata.n_obs:,} cells")
            del adata  # release the in-memory copy now that it's on disk

        if not sample_keys:
            raise TCRsiftValidationError(
                "No samples loaded successfully",
                hint="Check that at least one sample has valid VDJ or GEX data.",
            )

        logger.info(f"Concatenating {len(sample_keys)} samples ({total_cells:,} total cells)")
        out_path = spill_dir / "combined.h5ad"
        with tqdm(
            total=1, desc="Concatenating samples", unit="step", disable=not show_progress
        ) as pbar:
            concat_on_disk(
                spill_paths, out_path,
                join="outer", label="sample", keys=sample_keys,
            )
            combined = ad.read_h5ad(out_path)
            pbar.update(1)

    # Tempdir is cleaned up here; `combined` is fully in-memory.
    # Inverse of `_ensure_x`: an all-VDJ-only sheet produces n_vars == 0,
    # in which case the synthesized empty placeholders should look like the
    # original no-X state to downstream callers.
    if combined.n_vars == 0:
        combined.X = None

    # Store sample sheet as uns as a JSON string. The previous form,
    # `to_dataframe().to_dict()`, produced `{col: {row_idx: val}}` whose
    # int row-index keys break h5ad serialization. Other structured forms
    # (DataFrame directly, `to_dict(orient="list")`) hit anndata/h5py
    # mixed-None-object-column limitations. JSON is unambiguous about
    # None and stores as a single string. Read back with
    # `pd.read_json(io.StringIO(adata.uns["sample_sheet"]), orient="records")`.
    combined.uns["sample_sheet"] = sample_sheet.to_dataframe().to_json(
        orient="records"
    )

    logger.info(f"Successfully loaded {combined.n_obs:,} cells from {len(sample_keys)} samples")

    return combined
