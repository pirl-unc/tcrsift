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


def _reembed_leiden(
    adata,
    informative_genes,
    *,
    n_pcs,
    n_neighbors,
    resolution,
    n_top_genes,
    random_state,
):
    """Deterministic sub-cluster partition for :func:`refine_cluster`.

    counts → ``normalize_pearson_residuals`` → ``nan_to_num`` → PCA
    (``n_comps=min(n_pcs, n_obs-1, n_vars-1)``) → neighbors → Leiden
    (``flavor="igraph"``, ``n_iterations=2``, ``directed=False``), everything on
    ``random_state``. When ``informative_genes`` is given the panel is used
    **verbatim** (var-order intersection, no MT/receptor/min-cells curation) so a
    hand-rolled pass running the identical recipe reproduces the partition
    byte-for-byte; ``None`` falls back to seurat highly-variable genes on a
    receptor/MT-excluded pool. No UMAP — the partition is all that is needed.

    Returns the Leiden labels (object ndarray) aligned to ``adata``'s cells.
    """
    import scanpy as sc

    def _counts(a):
        return a.layers["counts"] if "counts" in a.layers else a.X

    if informative_genes:
        want = {str(g).upper() for g in informative_genes}
        panel = [g for g in adata.var_names if str(g).upper() in want]
        if not panel:
            raise ValueError(
                "refine_cluster: none of informative_genes are present in the "
                "compartment"
            )
    else:
        from .signature_methods import is_receptor_gene

        pool = adata.copy()
        counts = _counts(pool)
        pool.X = counts.copy() if hasattr(counts, "copy") else np.asarray(counts)
        keep = ~np.array([is_receptor_gene(g) for g in pool.var_names])
        keep &= ~np.asarray(pool.var_names.str.upper().str.startswith("MT-"))
        pool = pool[:, keep].copy()
        sc.pp.normalize_total(pool, target_sum=1e4)
        sc.pp.log1p(pool)
        n_top = min(n_top_genes, pool.n_vars)
        sc.pp.highly_variable_genes(pool, n_top_genes=n_top, flavor="seurat")
        panel = list(pool.var_names[pool.var["highly_variable"].to_numpy()])

    emb = adata[:, panel].copy()
    counts = _counts(emb)
    emb.X = counts.copy() if hasattr(counts, "copy") else np.asarray(counts)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.experimental.pp.normalize_pearson_residuals(emb)
    emb.X = np.nan_to_num(np.asarray(emb.X), nan=0.0, posinf=0.0, neginf=0.0)

    n_comps = min(n_pcs, emb.n_obs - 1, emb.n_vars - 1)
    sc.pp.pca(emb, n_comps=n_comps, random_state=random_state)
    k = min(n_neighbors, max(2, emb.n_obs - 1))
    sc.pp.neighbors(emb, use_rep="X_pca", n_neighbors=k, random_state=random_state)
    sc.tl.leiden(
        emb, resolution=resolution, random_state=random_state,
        flavor="igraph", n_iterations=2, directed=False,
    )
    return emb.obs["leiden"].to_numpy()


def refine_cluster(
    adata,
    mask,
    *,
    relabel_fn,
    informative_genes=None,
    resolution=1.0,
    n_pcs=30,
    n_neighbors=15,
    n_top_genes=2000,
    label_col="phenotype",
    final_cluster_col=None,
    min_cells=20,
    random_state=0,
):
    """Recovery / disaggregation pass: re-embed ONE compartment on its own
    Pearson residuals, sub-cluster it, and relabel each sub-cluster via a
    caller callback (#352).

    ``annotate_clusters`` labels clusters over the GLOBAL embedding, whose axes
    are dominated by whole-atlas structure; within-compartment biology — rare
    pDC/cDC1 co-embedded in myeloid, genuine γδ vs ambient-αβ CD8 sharing one
    cytotoxic cluster, tumor sub-states — does not resolve there. This
    primitive isolates a compartment (``mask``), re-embeds it *alone* and asks
    ``relabel_fn`` to name each resulting sub-cluster. It is the shared engine
    for hand-rolled passes like γδ/αβ-CD8 disaggregation (see
    :func:`gd_cd8_relabel`), DC recovery, and tumor sub-state splitting — the
    discriminating logic lives entirely in ``relabel_fn``.

    **The re-embed is your recipe, exposed:** counts → Pearson residuals →
    ``nan_to_num`` → ``pca(n_comps=min(n_pcs, n_obs-1, n_vars-1))`` →
    ``neighbors(n_neighbors)`` → ``leiden(resolution, flavor="igraph",
    n_iterations=2, directed=False)``, all on ``random_state``, no UMAP. When
    ``informative_genes`` is provided the clustering panel is used **verbatim**
    (no MT/receptor/min-cells curation), so a hand-rolled pass running the same
    steps/seeds reproduces the sub-partition byte-for-byte and can be deleted in
    favor of this primitive.

    Parameters
    ----------
    adata
        The annotated atlas. Raw counts must live in ``.layers["counts"]`` (or
        in ``.X`` when there is no counts layer) — the re-embed needs counts for
        Pearson residuals. ``.X`` is what ``relabel_fn`` sees for gene levels, so
        keep it on the scale your callback's thresholds assume (h37: log1p
        CP10K). Modified in place: ``obs[label_col]`` (and ``final_cluster_col``)
        are overwritten for the masked cells only; the rest of the atlas is
        untouched.
    mask
        Selects the compartment to refine. One of: a boolean mask (array/list)
        over ``adata.obs`` rows; a positional-index list; a label string matched
        against ``adata.obs[label_col]``; or a predicate ``adata -> boolean``.
    relabel_fn
        ``relabel_fn(sub_adata) -> label`` — called once per sub-cluster with a
        copy of the ORIGINAL-scale AnnData slice for that sub-cluster (``.X`` as
        in ``adata``, ``obs`` intact, plus ``obs["leiden"]`` = the local
        sub-cluster id). Compute whatever aggregate signal you need (mean
        ``TRDC``, αβ-contig fraction, dominant program, …) and return the new
        label string; a ``(label, final_cluster)`` pair to also stamp
        ``final_cluster_col``; or ``None`` to leave that sub-cluster's existing
        labels untouched. A ``None`` label inside the pair likewise leaves the
        label alone while still stamping ``final_cluster``.
    informative_genes
        Clustering panel, used **verbatim** for a byte-faithful re-embed
        (``None`` → seurat highly-variable genes on a receptor/MT-excluded pool,
        ``n_top_genes`` of them).
    resolution, n_pcs, n_neighbors, n_top_genes
        Re-embed knobs (defaults ``1.0`` / ``30`` / ``15`` / ``2000``); see the
        recipe above.
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
        Seed threaded to PCA, neighbors and Leiden.

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

    # --- re-embed the compartment alone (caller's recipe, verbatim panel) --
    compartment = adata[sel].copy()
    leiden = _reembed_leiden(
        compartment,
        informative_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        resolution=resolution,
        n_top_genes=n_top_genes,
        random_state=random_state,
    )
    compartment.obs["leiden"] = leiden

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

    n_relabeled = 0
    for sub_cl in pd.unique(leiden):
        sc_mask = leiden == sub_cl
        cell_names = compartment.obs_names[sc_mask]
        result = relabel_fn(compartment[sc_mask].copy())
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


def gd_cd8_relabel(
    *,
    gene="TRDC",
    ab_col="has_ab_contig",
    trdc_min=0.9,
    ab_max=0.25,
    gd_label="gd T",
    cd8_label="CD8 effector/cytotoxic",
    gd_cluster="gdt",
    cd8_cluster="gdt_cd8",
    ambiguous_band=(0.25, 0.55),
):
    """Build a :func:`refine_cluster` ``relabel_fn`` that splits a mixed
    γδ / ambient-αβ CD8 sub-cluster (#351).

    Each re-embedded sub-cluster is called γδ **iff both** conditions hold: its
    mean ``gene`` level (default ``TRDC``) is ``>= trdc_min`` **and** its mean
    αβ-contig capture (``obs[ab_col]``) is ``< ab_max`` — the h37 rule
    (``TRDC >= 0.9`` on log1p CP10K, αβ ``< 0.25``). Aggregating per sub-cluster
    makes it robust to the ~46% per-cell ``TRDC`` dropout that defeats a per-cell
    threshold. Returns ``(label, final_cluster)`` — γδ →
    ``(gd_label, gd_cluster)``, else ``(cd8_label, cd8_cluster)``.

    ``gene`` is read from ``sub_adata.X``, which must be on the scale the
    thresholds assume (h37: log1p CP10K). ``ab_col`` is a caller-populated
    boolean obs column — materialize it however your VDJ capture is encoded
    (h37, from a ``"α/β"`` paired-CDR3 string: ``(s != "/") & (s != "") &
    s.str.contains("[A-Z]")``); this factory does not hard-code any CDR3 format.
    When a sub-cluster's mean capture lands in ``ambiguous_band`` (between the γδ
    and αβ regimes) a warning is logged — the call is made, not failed.
    """
    lo, hi = ambiguous_band

    def relabel(sub):
        if ab_col not in sub.obs.columns:
            raise KeyError(
                f"gd_cd8_relabel: obs[{ab_col!r}] not found — populate the "
                "αβ-contig-capture boolean column before refine_cluster"
            )
        col = sub[:, gene].X
        col = col.toarray() if hasattr(col, "toarray") else np.asarray(col)
        trdc_mean = float(col.mean()) if col.size else 0.0
        ab_mean = float(np.asarray(sub.obs[ab_col], dtype=float).mean())
        is_gd = (trdc_mean >= trdc_min) and (ab_mean < ab_max)
        if lo < ab_mean < hi:
            logger.warning(
                "gd_cd8_relabel: sub-cluster mean αβ-capture %.2f is between the "
                "γδ (<%.2f) and αβ (>=%.2f) regimes; calling it %s (TRDC=%.2f)",
                ab_mean, ab_max, hi, gd_label if is_gd else cd8_label, trdc_mean,
            )
        return (gd_label, gd_cluster) if is_gd else (cd8_label, cd8_cluster)

    return relabel
