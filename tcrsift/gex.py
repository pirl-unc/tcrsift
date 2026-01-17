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
Flexible gene expression augmentation for TCRsift.

Functions for adding gene expression data to TCR/clonotype DataFrames,
with support for custom gene lists and groupings.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

from .validation import TCRsiftValidationError, validate_file_exists

logger = logging.getLogger(__name__)

# Default gene lists for T cell analysis
DEFAULT_GENE_LIST = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    "GZMA", "GZMB", "PRF1", "IFNG",
    "PDCD1", "CD274", "PDCD1LG2",  # PD-1/PD-L1
    "CTLA4", "LAG3", "HAVCR2", "TIGIT",  # Checkpoint
    "FOXP3", "IL2RA",  # Treg
]

DEFAULT_GENE_GROUPS = {
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
}


def augment_with_gex(
    df: pd.DataFrame,
    gex_path: str | Path,
    *,
    barcode_col: str = "barcode",
    gene_list: list[str] | None = None,
    gene_groups: dict[str, list[str]] | None = None,
    col_prefix: str = "gex",
    include_qc: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Augment a DataFrame with gene expression data from a 10x HDF5 file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with barcode column
    gex_path : str or Path
        Path to 10x filtered_feature_bc_matrix.h5 file
    barcode_col : str
        Name of barcode column in df
    gene_list : list of str, optional
        Genes to extract (default: T cell markers)
    gene_groups : dict, optional
        Gene groups to compute (e.g., {"CD8": ["CD8A", "CD8B"]})
    col_prefix : str
        Prefix for new columns (default: "gex")
    include_qc : bool
        Include QC metrics (n_reads, n_genes, pct_mito)
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        DataFrame with gene expression columns added
    """
    import scanpy as sc
    from scipy.sparse import issparse

    gex_path = validate_file_exists(Path(gex_path), "GEX file")

    if gene_list is None:
        gene_list = DEFAULT_GENE_LIST.copy()
    if gene_groups is None:
        gene_groups = DEFAULT_GENE_GROUPS.copy()

    # Add genes from groups to gene_list
    all_genes = set(gene_list)
    for genes in gene_groups.values():
        all_genes.update(genes)
    gene_list = list(all_genes)

    if barcode_col not in df.columns:
        raise TCRsiftValidationError(
            f"Barcode column '{barcode_col}' not found in DataFrame",
            hint=f"Available columns: {list(df.columns)[:10]}",
        )

    if verbose:
        logger.info(f"Loading gene expression from {gex_path}")

    # Load expression matrix
    adata = sc.read_10x_h5(str(gex_path))
    gene_names = adata.var_names

    # Handle duplicated gene symbols
    duplicated = gene_names[gene_names.duplicated()].unique().tolist()
    gene_to_best_idx: dict[str, int] = {}
    if duplicated:
        if verbose:
            logger.info(f"  Found {len(duplicated)} duplicated gene symbols")
        for symbol in duplicated:
            locs = np.where(gene_names == symbol)[0]
            means = np.asarray(adata.X[:, locs].mean(axis=0)).ravel()
            gene_to_best_idx[symbol] = locs[np.argmax(means)]

    # Compute per-cell QC metrics
    if issparse(adata.X):
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel().astype(int)
        n_genes_detected = np.asarray((adata.X > 0).sum(axis=1)).ravel().astype(int)
    else:
        total_counts = adata.X.sum(axis=1).astype(int)
        n_genes_detected = (adata.X > 0).sum(axis=1).astype(int)

    # Mitochondrial percentage
    mito_mask = gene_names.str.upper().str.startswith("MT-")
    if mito_mask.any():
        mito_counts = (
            adata.X[:, mito_mask].sum(axis=1).A.ravel()
            if issparse(adata.X)
            else adata.X[:, mito_mask].sum(axis=1)
        )
        pct_mito = mito_counts / np.maximum(total_counts, 1) * 100
    else:
        pct_mito = np.full(len(total_counts), np.nan)

    # Map genes to indices
    available_genes = []
    missing_genes = []
    gene_to_idx: dict[str, int] = {}
    for gene in gene_list:
        if gene in duplicated:
            gene_to_idx[gene] = gene_to_best_idx[gene]
            available_genes.append(gene)
        elif gene in gene_names:
            gene_to_idx[gene] = gene_names.get_loc(gene)
            available_genes.append(gene)
        else:
            missing_genes.append(gene)

    if missing_genes and verbose:
        logger.info(f"  {len(missing_genes)} genes not found: {missing_genes[:5]}...")

    if not available_genes:
        raise TCRsiftValidationError(
            "None of the requested genes were found in the expression matrix",
            hint=f"Requested: {gene_list[:5]}... Available: {list(gene_names[:5])}...",
        )

    # Initialize output columns
    for gene in available_genes:
        df[f"{col_prefix}.{gene}"] = np.nan

    if include_qc:
        df[f"{col_prefix}.n_reads"] = np.nan
        df[f"{col_prefix}.n_genes"] = np.nan
        df[f"{col_prefix}.pct_mito"] = np.nan

    # Barcode mapping
    bc_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}

    # Populate values
    missing_barcodes = 0
    matched_barcodes = 0
    for i, bc in enumerate(df[barcode_col]):
        if bc not in bc_to_idx:
            missing_barcodes += 1
            continue

        idx = bc_to_idx[bc]
        matched_barcodes += 1

        # QC metrics
        if include_qc:
            df.at[i, f"{col_prefix}.n_reads"] = total_counts[idx]
            df.at[i, f"{col_prefix}.n_genes"] = n_genes_detected[idx]
            df.at[i, f"{col_prefix}.pct_mito"] = pct_mito[idx]

        # Gene expression
        for gene in available_genes:
            v = adata.X[idx, gene_to_idx[gene]]
            if hasattr(v, "toarray"):
                v = v.toarray()[0, 0]
            df.at[i, f"{col_prefix}.{gene}"] = v

    # Compute gene groups
    if verbose:
        logger.info("  Computing gene group signatures...")
    for group_name, genes in gene_groups.items():
        cols = [f"{col_prefix}.{g}" for g in genes if f"{col_prefix}.{g}" in df.columns]
        if cols:
            df[f"{col_prefix}.{group_name}"] = df[cols].mean(axis=1)
            if verbose:
                logger.info(f"    {group_name}: mean of {len(cols)} genes")

    if verbose:
        logger.info(f"  Matched {matched_barcodes:,}/{len(df):,} barcodes")
        logger.info(f"  Added {len(available_genes)} gene columns + {len(gene_groups)} group columns")

    return df


def aggregate_gex_by_clonotype(
    df: pd.DataFrame,
    group_col: str = "CDR3_pair",
    *,
    gex_prefix: str = "gex",
    operations: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Aggregate gene expression values by clonotype.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with per-cell GEX columns
    group_col : str
        Column to group by (default: CDR3_pair)
    gex_prefix : str
        Prefix of GEX columns (default: "gex")
    operations : list of str
        Aggregation operations (default: ["sum", "mean"])
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        Aggregated data with per-clonotype statistics
    """
    if operations is None:
        operations = ["sum", "mean"]

    # Find GEX columns
    gex_cols = [c for c in df.columns if c.startswith(f"{gex_prefix}.")]
    qc_cols = [c for c in gex_cols if any(x in c for x in ["n_reads", "n_genes", "pct_mito"])]
    expr_cols = [c for c in gex_cols if c not in qc_cols]

    if verbose:
        logger.info(f"Aggregating {len(expr_cols)} expression columns by {group_col}")

    result_data = []
    for clonotype, group in df.groupby(group_col):
        row = {group_col: clonotype, "total_cells.count": len(group)}

        for col in expr_cols:
            values = group[col].dropna()
            if len(values) > 0:
                if "sum" in operations:
                    row[f"{col}.sum"] = values.sum()
                if "mean" in operations:
                    row[f"{col}.mean"] = values.mean()
                if "max" in operations:
                    row[f"{col}.max"] = values.max()
                if "min" in operations:
                    row[f"{col}.min"] = values.min()

        result_data.append(row)

    result = pd.DataFrame(result_data)

    if verbose:
        logger.info(f"  Aggregated to {len(result):,} clonotypes")

    return result


def compute_cd4_cd8_counts(
    df: pd.DataFrame,
    group_col: str = "CDR3_pair",
    *,
    gex_prefix: str = "gex",
    cd4_col: str | None = None,
    cd8_col: str | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compute CD4-only and CD8-only cell counts per clonotype.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with per-cell data
    group_col : str
        Column to group by
    gex_prefix : str
        GEX column prefix
    cd4_col : str, optional
        CD4 expression column (auto-detected if not specified)
    cd8_col : str, optional
        CD8 expression column (auto-detected if not specified)
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        Clonotype counts with CD4_only.count and CD8_only.count columns
    """
    # Auto-detect columns
    if cd4_col is None:
        candidates = [f"{gex_prefix}.CD4", "CD4", "gex.CD4"]
        for c in candidates:
            if c in df.columns:
                cd4_col = c
                break
        if cd4_col is None:
            raise TCRsiftValidationError(
                "Could not find CD4 expression column",
                hint=f"Available columns: {list(df.columns)[:10]}",
            )

    if cd8_col is None:
        candidates = [f"{gex_prefix}.CD8", "CD8", "gex.CD8"]
        for c in candidates:
            if c in df.columns:
                cd8_col = c
                break
        if cd8_col is None:
            raise TCRsiftValidationError(
                "Could not find CD8 expression column",
                hint=f"Available columns: {list(df.columns)[:10]}",
            )

    if verbose:
        logger.info(f"Computing CD4/CD8 counts using {cd4_col} and {cd8_col}")

    result_data = []
    for clonotype, group in df.groupby(group_col):
        cd4_vals = group[cd4_col].fillna(0)
        cd8_vals = group[cd8_col].fillna(0)

        cd4_only = ((cd4_vals > 0) & (cd8_vals == 0)).sum()
        cd8_only = ((cd8_vals > 0) & (cd4_vals == 0)).sum()
        total = len(group)

        result_data.append({
            group_col: clonotype,
            "total_cells.count": total,
            "CD4_only.count": cd4_only,
            "CD8_only.count": cd8_only,
        })

    result = pd.DataFrame(result_data)

    if verbose:
        logger.info(f"  {result['CD4_only.count'].sum():,} CD4-only cells")
        logger.info(f"  {result['CD8_only.count'].sum():,} CD8-only cells")

    return result
