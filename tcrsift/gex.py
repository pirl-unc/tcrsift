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

from .signatures import MARKER_PANEL_HGNC, NAIVE_STEM_GENES_HGNC
from .validation import TCRsiftValidationError, validate_file_exists

logger = logging.getLogger(__name__)

# Default gene lists for T cell analysis (includes full T_CELL_SIGNATURES support)
DEFAULT_GENE_LIST = list(
    dict.fromkeys(
        [
            "CD3D",
            "CD3E",
            "CD3G",
            "CD4",
            "CD8A",
            "CD8B",
            "GZMA",
            "GZMB",
            "PRF1",
            "IFNG",
            "PDCD1",
            "CD274",
            "PDCD1LG2",  # PD-1/PD-L1
            "CTLA4",
            "LAG3",
            "HAVCR2",
            "TIGIT",  # Checkpoint
            "FOXP3",
            "IL2RA",  # Treg
            *MARKER_PANEL_HGNC,
            *NAIVE_STEM_GENES_HGNC,
        ]
    )
)

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
    from .genes import adata_mito_mask, resolve_adata_var_map

    if gex_path.suffix.lower() == ".h5ad":
        adata = sc.read_h5ad(str(gex_path))
    else:
        adata = sc.read_10x_h5(str(gex_path))
    gene_names = adata.var_names

    # Handle duplicated gene symbols for informational logging
    duplicated = gene_names[gene_names.duplicated()].unique().tolist()
    if duplicated and verbose:
        logger.info(f"  Found {len(duplicated)} duplicated gene symbols")

    # Compute per-cell QC metrics
    if issparse(adata.X):
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel().astype(int)
        n_genes_detected = np.asarray((adata.X > 0).sum(axis=1)).ravel().astype(int)
    else:
        total_counts = adata.X.sum(axis=1).astype(int)
        n_genes_detected = (adata.X > 0).sum(axis=1).astype(int)

    # Mitochondrial percentage (supports human MT-*, mouse mt-*, and Ensembl IDs)
    mito_mask = adata_mito_mask(adata)
    if mito_mask.any():
        mito_counts = (
            np.asarray(adata.X[:, mito_mask].sum(axis=1)).ravel()
            if issparse(adata.X)
            else adata.X[:, mito_mask].sum(axis=1)
        )
        pct_mito = mito_counts / np.maximum(total_counts, 1) * 100
    else:
        pct_mito = np.full(len(total_counts), np.nan)

    # Map requested genes (with HGNC/MGI/Ensembl + duplicate & synonym resolution) to column indices
    sym_upper, sym_to_var = resolve_adata_var_map(adata, gene_list)
    var_to_idx = {str(vn): idx for idx, vn in enumerate(adata.var_names)}
    available_genes = []
    missing_genes = []
    gene_to_idx: dict[str, int] = {}
    for gene in gene_list:
        g_up = str(gene).upper()
        mapped_var = sym_to_var.get(g_up)
        if mapped_var is not None and mapped_var in var_to_idx:
            gene_to_idx[gene] = var_to_idx[mapped_var]
            available_genes.append(gene)
            continue
        locs = np.where(sym_upper == g_up)[0]
        if locs.size == 0:
            missing_genes.append(gene)
        elif locs.size == 1 or adata.n_obs == 0:
            gene_to_idx[gene] = int(locs[0])
            available_genes.append(gene)
        else:
            sub_x = adata.X[:, locs]
            if issparse(sub_x):
                means = np.asarray(sub_x.mean(axis=0)).ravel()
            else:
                means = np.nanmean(np.asarray(sub_x, dtype=float), axis=0)
            means = np.where(np.isfinite(means), means, -np.inf)
            gene_to_idx[gene] = int(locs[int(np.argmax(means))])
            available_genes.append(gene)

    if missing_genes and verbose:
        logger.info(f"  {len(missing_genes)} genes not found: {missing_genes[:5]}...")

    if not available_genes:
        raise TCRsiftValidationError(
            "None of the requested genes were found in the expression matrix",
            hint=f"Requested: {gene_list[:5]}... Available: {list(gene_names[:5])}...",
        )

    # Copy input DataFrame so callers' frames are never mutated in place
    df = df.copy()

    # Barcode mapping
    bc_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}

    # Populate values positionally (safe when df.index contains duplicates)
    n_rows = len(df)
    adata_indices = np.full(n_rows, -1, dtype=np.int64)
    for pos, bc in enumerate(df[barcode_col].to_numpy()):
        if bc in bc_to_idx:
            adata_indices[pos] = bc_to_idx[bc]
    matched_mask = adata_indices >= 0
    matched_barcodes = int(matched_mask.sum())
    matched_adata_idx = adata_indices[matched_mask]

    if include_qc:
        n_reads_arr = np.full(n_rows, np.nan, dtype=float)
        n_genes_arr = np.full(n_rows, np.nan, dtype=float)
        pct_mito_arr = np.full(n_rows, np.nan, dtype=float)
        if matched_barcodes > 0:
            n_reads_arr[matched_mask] = total_counts[matched_adata_idx]
            n_genes_arr[matched_mask] = n_genes_detected[matched_adata_idx]
            pct_mito_arr[matched_mask] = pct_mito[matched_adata_idx]
        df[f"{col_prefix}.n_reads"] = n_reads_arr
        df[f"{col_prefix}.n_genes"] = n_genes_arr
        df[f"{col_prefix}.pct_mito"] = pct_mito_arr

    if matched_barcodes > 0:
        gene_cols_idx = [gene_to_idx[g] for g in available_genes]
        sub_mat = adata.X[matched_adata_idx, :][:, gene_cols_idx]
        if hasattr(sub_mat, "toarray"):
            sub_mat = sub_mat.toarray()
        sub_mat = np.asarray(sub_mat, dtype=float)
        for j, gene in enumerate(available_genes):
            vals = np.full(n_rows, np.nan, dtype=float)
            vals[matched_mask] = sub_mat[:, j]
            df[f"{col_prefix}.{gene}"] = vals
    else:
        for gene in available_genes:
            df[f"{col_prefix}.{gene}"] = np.nan

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
        logger.info(
            f"  Added {len(available_genes)} gene columns + {len(gene_groups)} group columns"
        )

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
    for clonotype, group in df.groupby(group_col, observed=True):
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


def _build_prefixed_col_lookup(columns, gex_prefix: str) -> dict[str, str]:
    """Case-insensitive + GENE_SYNONYMS-aware lookup for ``{gex_prefix}.<gene>`` columns."""
    from .genes import GENE_SYNONYMS

    prefix = f"{gex_prefix}."
    lookup: dict[str, str] = {}
    for col in columns:
        if col.startswith(prefix):
            lookup.setdefault(col[len(prefix) :].upper(), col)
    for group in GENE_SYNONYMS:
        upper_group = [m.upper() for m in group]
        found = next((lookup[m] for m in upper_group if m in lookup), None)
        if found is not None:
            for m in upper_group:
                lookup.setdefault(m, found)
    return lookup


def _score_gene_block(
    frame: pd.DataFrame,
    cols: list[str],
    *,
    zscore: bool,
    apply_log1p: bool = True,
) -> np.ndarray:
    """Per-cell score for one gene block: mean over genes of ``log1p(expr)``
    (or raw/already-logged ``expr`` when ``apply_log1p=False``), optionally
    z-scored per gene across the cells in ``frame`` first.

    The single per-cell scoring primitive shared by the per-clonotype and
    per-cell signature scorers (#313) — ``zscore=False`` reproduces the
    historical log1p-mean; ``zscore=True`` gives the z-scored program score
    (``mean_z(genes)``) the atlas colorings and vs-background panels use.
    """
    vals = frame[cols].fillna(0).to_numpy(dtype=float)
    if apply_log1p:
        vals = np.log1p(np.clip(vals, a_min=0.0, a_max=None))
    if zscore:
        mu = vals.mean(axis=0)
        sd = vals.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        vals = (vals - mu) / sd
    return vals.mean(axis=1)


def compute_signature_scores_per_clonotype(
    df: pd.DataFrame,
    *,
    signatures: dict[str, tuple[str, ...]] | None = None,
    group_col: str = "CDR3_pair",
    gex_prefix: str = "gex",
    cd8_only: bool | None = None,
    cd8_col: str | None = None,
    zscore: bool = False,
    species: str = "human",
    verbose: bool = True,
) -> pd.DataFrame:
    """Compute per-clonotype signature scores (#74).

    For each named signature (gene set), compute the per-cell
    ``log1p(expression)`` mean across the gene set (or the z-scored program
    score when ``zscore=True``, #313), then aggregate per clonotype as the
    **cell-count-weighted mean** (i.e. each cell contributes equally; cells
    aren't first summarized per sample).

    ``cd8_only`` defaults to ``None``, which behaves as ``True`` (CD8+ only)
    for backward compatibility **but drops CD4 cells** — pass it explicitly.
    For CD4-inclusive scoring (or any non-CD8 analysis) pass
    ``cd8_only=False``; the per-cell scorer
    (:func:`compute_signature_scores_per_cell`) is CD4-inclusive by default.
    ``zscore=True`` standardizes each gene across the (post-CD8-filter) cells
    before averaging, so a program score isn't dominated by its
    highest-expression gene.

    Parameters
    ----------
    df : pd.DataFrame
        Per-cell frame with GEX columns named ``{gex_prefix}.GENE``
        (the shape :func:`augment_with_gex` writes).
    signatures : dict[str, tuple[str, ...]] | None
        Mapping of signature name → gene-symbol tuple. Defaults to
        :data:`tcrsift.signatures.T_CELL_SIGNATURES`.
    group_col : str
        Column identifying the clonotype (default ``CDR3_pair``).
    gex_prefix : str
        Prefix on the GEX columns to look up (default ``gex``).
    cd8_only : bool
        When True (default) and a CD8 column is present, restrict to
        cells with ``CD8 > 0`` before aggregating. Catches mixed CD4/CD8
        clones where the relevant signature lives in the CD8 subset.
    cd8_col : str | None
        Override CD8 column name. Auto-detected from ``{gex_prefix}.CD8A``,
        ``{gex_prefix}.CD8``, etc. when None.
    verbose : bool
        Log progress.

    Returns
    -------
    pd.DataFrame
        One row per clonotype with columns ``[{group_col}] +
        ["signature_{name}" for name in signatures]``. Signatures whose
        gene set has no overlap with the available GEX columns yield
        NaN for every clone (with a single warning log line).
    """
    if signatures is None:
        from .signatures import MOUSE_T_CELL_SIGNATURES, T_CELL_SIGNATURES
        from .species import normalize_species

        signatures = (
            MOUSE_T_CELL_SIGNATURES
            if normalize_species(species) == "mouse"
            else T_CELL_SIGNATURES
        )

    if group_col not in df.columns:
        raise ValueError(
            f"compute_signature_scores_per_clonotype: missing {group_col!r} column"
        )

    if cd8_only is None:
        cd8_only = True
        if verbose:
            logger.warning(
                "compute_signature_scores_per_clonotype: cd8_only defaulted to "
                "True, which DROPS CD4 cells — pass cd8_only explicitly "
                "(cd8_only=False for CD4-inclusive scoring)."
            )

    col_lookup = _build_prefixed_col_lookup(df.columns, gex_prefix)
    sub = df
    if cd8_only:
        if cd8_col is None:
            for c in (
                f"{gex_prefix}.CD8A",
                f"{gex_prefix}.CD8",
                col_lookup.get("CD8A"),
                col_lookup.get("CD8"),
                "CD8A",
                "CD8",
                "Cd8a",
            ):
                if c and c in df.columns:
                    cd8_col = c
                    break
        if cd8_col is not None and cd8_col in df.columns:
            sub = df[df[cd8_col].fillna(0) > 0]
            if verbose:
                logger.info(
                    f"compute_signature_scores_per_clonotype: restricting to "
                    f"{len(sub):,}/{len(df):,} CD8+ cells via {cd8_col!r}"
                )

    score_cols: dict[str, str] = {}
    for sig_name, genes in signatures.items():
        gex_cols = list(
            dict.fromkeys(
                f"{gex_prefix}.{g}"
                if f"{gex_prefix}.{g}" in sub.columns
                else col_lookup[str(g).upper()]
                for g in genes
                if f"{gex_prefix}.{g}" in sub.columns or str(g).upper() in col_lookup
            )
        )
        if not gex_cols:
            if verbose:
                logger.warning(
                    f"compute_signature_scores_per_clonotype: signature "
                    f"{sig_name!r} — none of {list(genes)} found under "
                    f"prefix {gex_prefix!r}; emitting NaN"
                )
            sub = sub.assign(**{f"_signature_{sig_name}": np.nan})
            score_cols[sig_name] = f"_signature_{sig_name}"
            continue
        per_cell = _score_gene_block(sub, gex_cols, zscore=zscore)
        sub = sub.assign(**{f"_signature_{sig_name}": per_cell})
        score_cols[sig_name] = f"_signature_{sig_name}"

    if len(sub) == 0:
        # No rows after filtering — emit NaN-valued result keyed by all
        # clonotypes in the original frame.
        result = (
            df[[group_col]]
            .drop_duplicates()
            .assign(**{f"signature_{s}": np.nan for s in signatures})
        )
        return result.reset_index(drop=True)

    grouped = sub.groupby(group_col, observed=True)[
        list(score_cols.values())
    ].mean().reset_index()

    grouped = grouped.rename(
        columns={src: f"signature_{name}" for name, src in score_cols.items()}
    )

    if verbose:
        logger.info(
            f"compute_signature_scores_per_clonotype: scored "
            f"{len(grouped):,} clonotypes across {len(signatures)} signatures"
        )
    return grouped


def compute_signature_scores_per_cell(
    adata,
    *,
    signatures: dict[str, tuple[str, ...]] | None = None,
    zscore: bool = True,
    layer: str | None = None,
    log1p: bool | None = None,
    key_prefix: str = "signature_",
    on_missing: str = "warn",
    species: str | None = None,
    verbose: bool = True,
):
    """Per-CELL signature scores written to ``adata.obs`` (#313).

    The per-cell sibling of :func:`compute_signature_scores_per_clonotype`:
    the scores that colour the single-cell UMAP and drive the
    signature-vs-background panels. For each named gene set, extract the
    symbol-keyed expression (via the shared Ensembl-aware resolver) and score
    each cell as the mean over present genes of z-scored ``log1p`` expression
    (``zscore=True``, the default) or plain ``log1p``-mean (``zscore=False``),
    writing it to ``adata.obs[f"{key_prefix}{name}"]``.

    When ``log1p=None`` (default), automatically avoids double-logging matrices
    already marked ``adata.uns["log1p"]``, matrices with negative (scaled)
    values, or non-integer normalized expression matrices, and applies CP10K
    normalization before ``log1p`` on raw count matrices. Pass ``log1p=True``
    or ``log1p=False`` to force or disable ``np.log1p`` explicitly.

    Unlike the per-clonotype scorer this is CD4-INCLUSIVE (no CD8 gate) — it
    scores every cell; restrict upstream for a lineage subset. For signed /
    published registry signatures (MANAscore, NeoTCR) use
    :func:`tcrsift.signature_methods.score_by_name` instead.

    Both scorers share the one per-cell primitive (:func:`_score_gene_block`)
    and default to the same :data:`tcrsift.signatures.T_CELL_SIGNATURES`, so a
    gene set is defined once and scored identically per-cell and per-clone.

    Returns the same ``adata`` (mutated in place); signatures with no gene
    overlap get an all-NaN column and a warning.
    """
    from .signature_methods import expression_frame_from_adata

    if signatures is None:
        from .signatures import MOUSE_T_CELL_SIGNATURES, T_CELL_SIGNATURES
        from .species import resolve_species

        resolved_species = resolve_species(
            species or getattr(adata, "uns", {}).get("species", "auto"),
            adata=adata,
            warn_mismatch=False,
        )
        signatures = (
            MOUSE_T_CELL_SIGNATURES
            if resolved_species == "mouse"
            else T_CELL_SIGNATURES
        )

    all_genes = sorted({g for genes in signatures.values() for g in genes})
    expr = expression_frame_from_adata(
        adata, all_genes, layer=layer, on_missing=on_missing
    )

    # Decide whether to CP10K-normalize / log1p based on matrix scale
    if log1p is not None:
        do_log1p = bool(log1p)
    else:
        src_mat = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
        if expr.empty:
            do_log1p = False
        else:
            arr = expr.to_numpy(dtype=float)
            finite = arr[np.isfinite(arr)]
            has_neg = bool(finite.size and finite.min() < 0)
            already_logged = bool(layer is None and "log1p" in getattr(adata, "uns", {}))
            if has_neg or already_logged:
                do_log1p = False
            else:
                # Non-integer bounded values without uns["log1p"] are already normalized/logged
                x_max = float(finite.max()) if finite.size else 0.0
                is_int_like = bool(
                    finite.size and np.allclose(finite, np.round(finite), atol=1e-6)
                )
                if not is_int_like and x_max <= 20.0:
                    do_log1p = False
                else:
                    # Raw counts: gate CP10K normalization on the full matrix (not
                    # just the queried signature subset, which may be low-count TFs/checkpoints).
                    mat_max = (
                        float(src_mat.max())
                        if (src_mat.shape[0] and src_mat.shape[1])
                        else x_max
                    )
                    if mat_max > 25.0 or is_int_like or layer == "counts":
                        totals = np.asarray(src_mat.sum(axis=1)).ravel().astype(float)
                        scale = np.where(totals > 0, 1e4 / totals, 0.0)
                        expr = expr.mul(scale, axis=0)
                    do_log1p = True

    present = set(expr.columns)
    n_scored = 0
    for name, genes in signatures.items():
        cols = [g for g in genes if g in present]
        col_name = f"{key_prefix}{name}"
        if not cols:
            if verbose:
                logger.warning(
                    "compute_signature_scores_per_cell: signature %r — none of "
                    "%s present; emitting NaN", name, list(genes),
                )
            adata.obs[col_name] = np.nan
            continue
        adata.obs[col_name] = _score_gene_block(
            expr, cols, zscore=zscore, apply_log1p=do_log1p
        )
        n_scored += 1
    if verbose:
        logger.info(
            "compute_signature_scores_per_cell: scored %d/%d signatures over "
            "%d cells", n_scored, len(signatures), adata.n_obs,
        )
    return adata


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
    # Auto-detect columns (case-insensitive + mouse MGI/synonym aware)
    col_lookup = _build_prefixed_col_lookup(df.columns, gex_prefix)
    if cd4_col is None:
        candidates = [
            f"{gex_prefix}.CD4",
            col_lookup.get("CD4"),
            "CD4",
            "Cd4",
            "gex.CD4",
            "gex.Cd4",
        ]
        for c in candidates:
            if c and c in df.columns:
                cd4_col = c
                break
        if cd4_col is None:
            raise TCRsiftValidationError(
                "Could not find CD4 expression column",
                hint=f"Available columns: {list(df.columns)[:10]}",
            )

    if cd8_col is None:
        candidates = [
            f"{gex_prefix}.CD8",
            col_lookup.get("CD8"),
            col_lookup.get("CD8A"),
            col_lookup.get("CD8B"),
            "CD8",
            "CD8A",
            "Cd8a",
            "Cd8b1",
            "gex.CD8",
        ]
        for c in candidates:
            if c and c in df.columns:
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
    for clonotype, group in df.groupby(group_col, observed=True):
        cd4_vals = group[cd4_col].fillna(0)
        cd8_vals = group[cd8_col].fillna(0)

        cd4_only = ((cd4_vals > 0) & (cd8_vals == 0)).sum()
        cd8_only = ((cd8_vals > 0) & (cd4_vals == 0)).sum()
        total = len(group)

        result_data.append(
            {
                group_col: clonotype,
                "total_cells.count": total,
                "CD4_only.count": cd4_only,
                "CD8_only.count": cd8_only,
            }
        )

    result = (
        pd.DataFrame(result_data)
        if result_data
        else pd.DataFrame(
            columns=[group_col, "total_cells.count", "CD4_only.count", "CD8_only.count"]
        )
    )

    if verbose:
        logger.info(f"  {result['CD4_only.count'].sum():,} CD4-only cells")
        logger.info(f"  {result['CD8_only.count'].sum():,} CD8-only cells")

    return result


# =============================================================================
# Externally-chosen clone annotation (#302)
# =============================================================================


def _normalize_cdr3_series(s: pd.Series) -> pd.Series:
    """Normalize a CDR3 column for key matching: NaN→"", stripped, upper.

    Keeping this in one place is the whole point of :func:`annotate_chosen` —
    a chosen-clone list assembled in one experiment and a scored table from
    another must agree on whitespace/case before set membership means
    anything.
    """
    return s.fillna("").astype(str).str.strip().str.upper()


def _normalize_cdr3_value(value) -> str:
    """Scalar counterpart of :func:`_normalize_cdr3_series` (NaN/None→"")."""
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value).strip().upper()


def _extract_chosen_chains(chosen, alpha_col: str, beta_col: str, match: str):
    """Pull ``(alphas, betas)`` lists out of a heterogeneous chosen input.

    Accepts a :class:`pandas.DataFrame` (reads ``alpha_col`` / ``beta_col``),
    an iterable of ``(alpha, beta)`` tuples, or — for ``match="beta"`` — an
    iterable of bare β strings. ``alphas`` entries are ``None`` when only β is
    available/needed.
    """
    if isinstance(chosen, pd.DataFrame):
        if beta_col not in chosen.columns:
            raise TCRsiftValidationError(
                f"chosen DataFrame missing β column {beta_col!r}",
                hint=f"Available columns: {list(chosen.columns)[:10]}",
            )
        betas = chosen[beta_col].tolist()
        if match == "pair":
            if alpha_col not in chosen.columns:
                raise TCRsiftValidationError(
                    f"chosen DataFrame missing α column {alpha_col!r} "
                    "(required for match='pair')",
                    hint=f"Available columns: {list(chosen.columns)[:10]}",
                )
            alphas = chosen[alpha_col].tolist()
        else:
            alphas = [None] * len(betas)
        return alphas, betas

    items = list(chosen)
    if match == "beta":
        # A bare β string or the β slot of an (α, β) tuple both work.
        betas = [it[1] if isinstance(it, (tuple, list)) else it for it in items]
        return [None] * len(betas), betas

    alphas, betas = [], []
    for it in items:
        if not isinstance(it, (tuple, list)) or len(it) < 2:
            raise TCRsiftValidationError(
                "match='pair' needs (alpha, beta) pairs",
                hint=f"Got element {it!r}; pass a DataFrame or (α, β) tuples.",
            )
        alphas.append(it[0])
        betas.append(it[1])
    return alphas, betas


def annotate_chosen(
    table: pd.DataFrame,
    chosen,
    *,
    alpha_col: str = "CDR3_alpha",
    beta_col: str = "CDR3_beta",
    match: str = "pair",
    name: str = "chosen",
) -> pd.Series:
    """Flag rows of ``table`` that appear in an externally-chosen clone list (#302).

    Centralizes the αβ-pair / CDR3β key matching that otherwise gets
    re-implemented (slightly differently) every time a clone list selected in
    one experiment is contextualized against a clone table scored in another
    (e.g. peptide-culture picks projected onto the TIL repertoire). The
    returned mask is the single primitive behind "chosen vs background",
    "chosen vs all", and "chosen + overlap pool" views — they differ only in
    which table you pass and which other masks you combine it with.

    Matching is exact on normalized CDR3 strings (NaN→"", stripped,
    upper-cased on both sides), never fuzzy — selection lists are curated, so
    a near-miss is a data-entry problem to surface, not silently absorb. Empty
    keys never match.

    Parameters
    ----------
    table : pd.DataFrame
        Scored / harmonized clone table to annotate. Must contain ``beta_col``
        (and ``alpha_col`` when ``match="pair"``). One clone per row, as
        produced by :func:`tcrsift.til_select.build_harmonized_table`.
    chosen : pd.DataFrame | iterable
        The externally-chosen clones. A DataFrame (uses ``alpha_col`` /
        ``beta_col``), an iterable of ``(alpha, beta)`` tuples, or — for
        ``match="beta"`` — an iterable of bare β strings.
    alpha_col, beta_col : str
        CDR3α / CDR3β column names (shared by ``table`` and a chosen
        DataFrame). Defaults match the harmonized-table convention.
    match : {"pair", "beta"}
        ``"pair"`` requires both α and β to match (exact αβ pair).
        ``"beta"`` matches on CDR3β alone — the looser key for when chosen α
        chains are unavailable or unreliable.
    name : str
        Name of the returned Series.

    Returns
    -------
    pd.Series
        Boolean Series aligned to ``table.index``; True where the row's
        clone is in ``chosen``.
    """
    if match not in {"pair", "beta"}:
        raise TCRsiftValidationError(
            f"Invalid match: {match!r}",
            hint="Valid options are: 'pair', 'beta'",
        )
    if beta_col not in table.columns:
        raise TCRsiftValidationError(
            f"table missing β column {beta_col!r}",
            hint=f"Available columns: {list(table.columns)[:10]}",
        )
    if match == "pair" and alpha_col not in table.columns:
        raise TCRsiftValidationError(
            f"table missing α column {alpha_col!r} (required for match='pair')",
            hint=f"Available columns: {list(table.columns)[:10]}",
        )

    chosen_alphas, chosen_betas = _extract_chosen_chains(
        chosen, alpha_col, beta_col, match
    )

    tbl_beta = _normalize_cdr3_series(table[beta_col])
    if match == "beta":
        keyset = {_normalize_cdr3_value(b) for b in chosen_betas}
        keyset.discard("")
        mask = tbl_beta.isin(keyset)
    else:
        keyset = {
            (_normalize_cdr3_value(a), _normalize_cdr3_value(b))
            for a, b in zip(chosen_alphas, chosen_betas)
        }
        keyset = {(a, b) for (a, b) in keyset if a and b}
        tbl_alpha = _normalize_cdr3_series(table[alpha_col])
        pairs = pd.Series(zip(tbl_alpha, tbl_beta), index=table.index)
        mask = pairs.isin(keyset)

    return pd.Series(
        np.asarray(mask, dtype=bool), index=table.index, name=name
    )
