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

"""Single-cell embedding for the cell-level atlas (#311).

The cell-embedding step that runs *before* the clonotype aggregation: it turns
a per-cell AnnData into an integrated 2-D atlas (UMAP + Leiden clusters) the
annotator (#312) and atlas plots (#315) consume. Structure is computed on
**analytic Pearson residuals** of an informative-gene panel — a variance-
stabilizing transform that clusters on genuine biology rather than depth /
ambient noise — then PCA → optional Harmony batch integration → neighbors →
UMAP → Leiden. Signature scoring and labeling stay on log1p CP10K
(:mod:`tcrsift.gex`); this step only computes structure.

Heavy, optional dependencies (scanpy is already core; harmonypy + igraph +
leidenalg are the ``atlas`` extra) are imported lazily so a plain
``import tcrsift`` stays light.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np

logger = logging.getLogger(__name__)


def embed_cells(adata, config=None):
    """Variance-stabilized single-cell embedding (#311).

    Pearson residuals on an informative-gene panel → PCA → optional Harmony
    (over ``config.batch_key``) → neighbors → UMAP → Leiden. Returns a **copy**
    of ``adata`` with ``obsm["X_pca"]`` (and ``obsm["X_pca_harmony"]`` when a
    batch key is given), ``obsm["X_umap"]``, ``obs["leiden"]``, and the panel in
    ``uns["embed_panel"]``. Signature scores are NOT computed here.

    ``adata.X`` must be raw counts (or provide an ``adata.layers["counts"]``) —
    Pearson residuals are defined on counts. ``config`` is an
    :class:`tcrsift.config.EmbedConfig` (defaults used when None).

    Harmony is called through ``harmonypy.run_harmony`` directly rather than the
    scanpy wrapper: recent harmonypy returns ``Z_corr`` as ``(cells × features)``
    while the scanpy wrapper still assumes the legacy ``(features × cells)`` and
    transposes, corrupting ``X_pca_harmony`` — so we assign with an explicit
    orientation check instead.
    """
    import scanpy as sc

    from .signature_methods import is_receptor_gene

    if config is None:
        from .config import EmbedConfig

        config = EmbedConfig()

    adata = adata.copy()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Candidate pool for the clustering panel: drop TCR/Ig receptor transcripts
    # (so structure isn't driven by clonotype) and mitochondrial genes, then
    # drop genes seen in too few cells.
    var = adata.var_names
    keep = np.ones(adata.n_vars, dtype=bool)
    if config.exclude_receptor_genes:
        keep &= ~np.array([is_receptor_gene(g) for g in var])
    keep &= ~np.asarray(var.str.upper().str.startswith("MT-"))
    pool = adata[:, keep].copy()
    pool.X = pool.layers["counts"].copy()
    sc.pp.filter_genes(pool, min_cells=config.min_cells)

    if config.informative_genes:
        want = {str(g).upper() for g in config.informative_genes}
        panel = [g for g in pool.var_names if str(g).upper() in want]
        if not panel:
            raise ValueError(
                "embed_cells: none of informative_genes are present in the "
                "matrix after filtering"
            )
    else:
        # Highly-variable genes on a log-normalized copy (stable seurat flavor).
        lognorm = pool.copy()
        sc.pp.normalize_total(lognorm, target_sum=1e4)
        sc.pp.log1p(lognorm)
        n_top = min(config.n_top_genes, lognorm.n_vars)
        sc.pp.highly_variable_genes(lognorm, n_top_genes=n_top, flavor="seurat")
        panel = list(lognorm.var_names[lognorm.var["highly_variable"].to_numpy()])

    logger.info("embed_cells: clustering on %d informative genes", len(panel))

    # Embed on the Pearson residuals of the panel counts.
    emb = pool[:, panel].copy()
    emb.X = emb.layers["counts"].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.experimental.pp.normalize_pearson_residuals(emb)
    emb.X = np.nan_to_num(np.asarray(emb.X), nan=0.0, posinf=0.0, neginf=0.0)

    n_comps = min(config.n_pcs, emb.n_obs - 1, emb.n_vars - 1)
    sc.pp.pca(emb, n_comps=n_comps, random_state=config.seed)

    rep = "X_pca"
    if config.batch_key:
        if config.batch_key not in emb.obs.columns:
            raise ValueError(
                f"embed_cells: batch_key {config.batch_key!r} not in adata.obs"
            )
        import harmonypy

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ho = harmonypy.run_harmony(
                emb.obsm["X_pca"], emb.obs, [config.batch_key],
                max_iter_harmony=20,
            )
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] != emb.n_obs:  # orientation guard (see docstring)
            Z = Z.T
        emb.obsm["X_pca_harmony"] = Z
        rep = "X_pca_harmony"

    n_neighbors = min(config.n_neighbors, max(2, emb.n_obs - 1))
    sc.pp.neighbors(
        emb, use_rep=rep, n_neighbors=n_neighbors, random_state=config.seed
    )
    sc.tl.umap(emb, random_state=config.seed)
    sc.tl.leiden(
        emb, resolution=config.leiden_resolution, random_state=config.seed,
        flavor="igraph", n_iterations=2, directed=False,
    )

    # Write structure back to the full-gene adata (same cells, same order).
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    if "X_pca_harmony" in emb.obsm:
        adata.obsm["X_pca_harmony"] = emb.obsm["X_pca_harmony"]
    adata.obsm["X_umap"] = emb.obsm["X_umap"]
    adata.obs["leiden"] = emb.obs["leiden"].to_numpy()
    adata.uns["embed_panel"] = list(panel)
    return adata
