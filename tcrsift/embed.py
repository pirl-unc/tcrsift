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
from dataclasses import replace

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


def refine_cluster(
    adata,
    mask,
    *,
    relabel_fn,
    informative_genes=None,
    resolution=1.0,
    label_col="phenotype",
    final_cluster_col=None,
    min_cells=20,
    random_state=0,
    embed_config=None,
):
    """Recovery / disaggregation pass: re-embed ONE compartment on its own
    Pearson residuals, sub-cluster it, and relabel each sub-cluster via a
    caller callback (#352).

    ``annotate_clusters`` labels clusters over the GLOBAL embedding, whose axes
    are dominated by whole-atlas structure; within-compartment biology — rare
    pDC/cDC1 co-embedded in myeloid, genuine γδ vs ambient-αβ CD8 sharing one
    cytotoxic cluster, tumor sub-states — does not resolve there. This
    primitive isolates a compartment (``mask``), re-embeds it *alone* via
    :func:`embed_cells` (Pearson residuals → PCA → UMAP → Leiden at
    ``resolution`` over ``informative_genes``), then asks ``relabel_fn`` to name
    each resulting sub-cluster. It is the shared engine for hand-rolled passes
    like γδ/αβ-CD8 disaggregation, DC recovery, and tumor sub-state splitting —
    the discriminating logic lives entirely in ``relabel_fn``.

    Parameters
    ----------
    adata
        The annotated atlas. ``.X`` must be raw counts (or provide
        ``.layers["counts"]`` — see :func:`embed_cells`). Modified in place:
        ``adata.obs[label_col]`` is overwritten for the masked cells only; the
        rest of the atlas (labels, embedding) is untouched.
    mask
        Selects the compartment to refine. One of: a boolean mask (array/list)
        over ``adata.obs`` rows; a positional-index list; a label string matched
        against ``adata.obs[label_col]``; or a predicate ``adata -> boolean``.
    relabel_fn
        ``relabel_fn(sub_adata) -> label`` — called once per sub-cluster with a
        copy of the AnnData slice for that sub-cluster (raw counts + the local
        re-embedding in ``.obsm`` / ``.obs["leiden"]``). Compute whatever
        aggregate signal you need (mean ``TRDC``, αβ-contig fraction, dominant
        program, …) and return the new label string; a ``(label, final_cluster)``
        pair to also stamp ``final_cluster_col``; or ``None`` to leave that
        sub-cluster's existing labels untouched. A ``None`` label inside the pair
        likewise leaves the label alone while still stamping ``final_cluster``.
    informative_genes
        Clustering panel for the local re-embedding (``None`` → highly-variable
        genes picked within the compartment). Defaults to ``embed_config``'s.
    resolution
        Leiden resolution for the local re-embedding.
    label_col
        ``obs`` column read for a string ``mask`` and written with new labels
        (default ``"phenotype"``). Created (object dtype, ``None``-filled) if
        absent and ``mask`` is not a string.
    final_cluster_col
        Optional ``obs`` column for the callback's ``final_cluster`` tag.
    min_cells
        Skip re-embedding (labels unchanged) if the compartment is smaller than
        this — too few cells to resolve sub-structure.
    random_state
        Seed threaded to the local embedding.
    embed_config
        Optional base :class:`~tcrsift.config.EmbedConfig`; ``enabled``,
        ``informative_genes``, ``leiden_resolution``, ``seed`` and ``batch_key``
        are overridden (a single compartment is not batch-integrated).

    Returns
    -------
    The same ``adata``, mutated in place.
    """
    import pandas as pd

    # --- resolve the compartment mask -------------------------------------
    if callable(mask):
        sel = np.asarray(mask(adata))
    elif isinstance(mask, str):
        if label_col not in adata.obs.columns:
            raise ValueError(
                f"refine_cluster: label_col {label_col!r} not in adata.obs; "
                "cannot match a string mask"
            )
        sel = (adata.obs[label_col].astype(object) == mask).to_numpy()
    else:
        sel = np.asarray(mask)
    if sel.dtype != bool:  # a positional-index list
        idx = np.zeros(adata.n_obs, dtype=bool)
        idx[sel] = True
        sel = idx
    if sel.shape[0] != adata.n_obs:
        raise ValueError(
            f"refine_cluster: mask length {sel.shape[0]} != n_obs {adata.n_obs}"
        )

    n_sel = int(sel.sum())
    if n_sel < min_cells:
        logger.warning(
            "refine_cluster: compartment has %d cells (< min_cells=%d); "
            "skipping re-embedding, labels unchanged",
            n_sel, min_cells,
        )
        return adata

    # --- re-embed the compartment alone -----------------------------------
    from .config import EmbedConfig

    if embed_config is None:
        embed_config = EmbedConfig()
    cfg = replace(
        embed_config,
        enabled=True,
        informative_genes=(
            informative_genes
            if informative_genes is not None
            else embed_config.informative_genes
        ),
        leiden_resolution=resolution,
        batch_key=None,
        seed=random_state,
    )
    sub = embed_cells(adata[sel], cfg)

    # --- relabel each sub-cluster via the callback ------------------------
    if label_col not in adata.obs.columns:
        adata.obs[label_col] = pd.Series(
            [None] * adata.n_obs, index=adata.obs_names, dtype=object
        )
    labels = adata.obs[label_col].astype(object)
    finals = None
    if final_cluster_col is not None:
        existing = adata.obs.get(final_cluster_col)
        finals = (
            existing.astype(object)
            if existing is not None
            else pd.Series(
                [None] * adata.n_obs, index=adata.obs_names, dtype=object
            )
        )

    leiden = sub.obs["leiden"].to_numpy()
    n_relabeled = 0
    for sub_cl in pd.unique(leiden):
        sc_mask = leiden == sub_cl
        cell_names = sub.obs_names[sc_mask]
        result = relabel_fn(sub[sc_mask].copy())
        if result is None:
            continue
        new_label, final_tag = result if isinstance(result, tuple) else (result, None)
        if new_label is not None:
            labels.loc[cell_names] = new_label
            n_relabeled += len(cell_names)
        if finals is not None and final_tag is not None:
            finals.loc[cell_names] = final_tag

    adata.obs[label_col] = labels
    if finals is not None:
        adata.obs[final_cluster_col] = finals
    logger.info(
        "refine_cluster: re-embedded %d cells into %d sub-clusters; "
        "relabeled %d",
        n_sel, len(pd.unique(leiden)), n_relabeled,
    )
    return adata
