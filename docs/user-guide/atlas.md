# Single-Cell Atlas Workflow

Alongside clonotype prioritization, TCRsift provides a three-step single-cell
gene-expression workflow (quality control, embedding, and cluster annotation)
for building cell-level atlases from raw-count `AnnData` objects:

```
cell_qc_funnel  →  embed_cells  →  annotate_cells / annotate_clusters
  (coverage,         (Pearson-residual        (cluster-level cell-type
   mito, and          PCA → Harmony →          and T/B-state scoring)
   doublet filters)   UMAP → Leiden)
```

These functions are available in Python and through the
[`tcrsift cells qc`, `cells embed`, and `cells annotate`](cli.md#tcrsift-cells)
CLI subcommands. All three steps expect raw UMI counts in `adata.X` (or
`adata.layers["counts"]`), such as the `AnnData` returned by
`tcrsift.load_sample`.

## End-to-end Python example

```python
from tcrsift import annotate_cells, cell_qc_funnel, embed_cells
from tcrsift.config import EmbedConfig, LoadConfig

# 1. Per-cell QC and doublet filtering.
#    Passing LoadConfig() enables full PBMC coverage and mito bounds;
#    waterfall_df records how many cells each step removed.
adata, waterfall_df = cell_qc_funnel(adata, config=LoadConfig())
print(waterfall_df)  # columns: [step, reason, removed, remaining]

# 2. Pearson-residual PCA, Harmony batch integration, UMAP, and Leiden clustering.
atlas = embed_cells(adata, EmbedConfig(batch_key="sample"))

# 3. Cluster-level cell-type and T/B-state assignment
#    (writes obs["phenotype_base"] and obs["phenotype"]).
annotate_cells(atlas)
print(atlas.obs["phenotype"].value_counts())
```

## Choosing PBMC vs solid-tumor references

By default, `annotate_cells` scores clusters against `CELL_TYPE_SIGNATURES`,
`T_STATE_SIGNATURES`, and `B_STATE_SIGNATURES` from `tcrsift.signatures`, and
`cell_qc_funnel` checks cross-lineage doublets using `DEFAULT_LINEAGE_PROGRAMS`.
**Both defaults are designed for blood and PBMC datasets.** Because
`CELL_TYPE_SIGNATURES` omits non-hematopoietic solid-tissue lineages
(mesothelial, osteoclast, skeletal muscle, adipocyte, Schwann/neural), applying
the default reference to a solid-tumor sample will misclassify tumor or stromal
clusters (often as `Fibroblast`).

For solid tumors, pass the solid-tumor doublet and cell-type registries
explicitly. Custom `reference`, `t_state_reference`, and `b_state_reference`
dictionaries are forwarded to both signature scoring and cluster assignment:

```python
from tcrsift.qc import SOLID_TUMOR_LINEAGE_PROGRAMS
from tcrsift.signatures import SOLID_TUMOR_CELL_TYPE_SIGNATURES

# Solid-tumor doublet rules (merges osteoclasts into myeloid and keeps
# fibroblast collagen to catch tumor-immune doublets):
adata, wf = cell_qc_funnel(
    adata,
    config=LoadConfig(),
    lineage_sets=SOLID_TUMOR_LINEAGE_PROGRAMS,
)

# Solid-tumor cell-type reference:
annotate_cells(atlas, reference=SOLID_TUMOR_CELL_TYPE_SIGNATURES)
```

`PBMC_CELL_TYPE_SIGNATURES` is exported as an explicit alias for
`CELL_TYPE_SIGNATURES`.

**Mouse (`Mus musculus`) datasets:** `cell_qc_funnel`, `embed_cells`, and
`annotate_cells` auto-detect mouse `AnnData` inputs (from `mt-*` mitochondrial
prefixes, MGI Title-case symbols such as `Cd8b1`, or `ENSMUSG*` Ensembl IDs; override
with `species="mouse"` or `--species mouse`) and switch automatically to the
murine marker registries (`MOUSE_CELL_TYPE_SIGNATURES`,
`MOUSE_SOLID_TUMOR_CELL_TYPE_SIGNATURES`, `MOUSE_LINEAGE_PROGRAMS`, and
`MOUSE_GATES`, including NK1.1/`Klrb1c`, `Ncr1`, `Foxp3`, and `Trdc` gates).

## Identifying malignant clusters with `tumor_override`

Neither cell-type registry includes a generic malignant signature, because
malignant transcriptional programs vary widely by cancer type and collagen-rich
tumor clusters otherwise score as fibroblasts. Instead, `tumor_override` labels
malignant clusters from cancer-testis antigen (CTA) gene co-expression and
optional lineage transcription factors:

- **Primary rule:** relabel a cluster `Tumor` when at least a target fraction
  of cells each express ≥2 distinct genes from the marker panel.
- **Rescue rule:** relabel `Tumor` when a tissue lineage transcription factor
  is elevated **and** ≥1 panel marker is expressed in a target fraction of
  cells (recovering lower-coverage malignant clusters).

By default, the marker panel combines `oncoref.CTA_gene_names()` with
`oncoref.cta_clinical_target_gene_names()` (retaining clinically targeted CTA
genes such as `CTAG1A`/`CTAG1B` [encoding NY-ESO-1], `CTAG2` [encoding LAGE-1],
and `MAGEA11` alongside strictly testis-restricted genes). Only the rescue
transcription factors need to be specified per tumor histology:

```python
import oncoref
from tcrsift.annotate_cells import tumor_override

# Osteosarcoma: default CTA panel + bone/mesenchymal lineage TFs
osteosarcoma = tumor_override(lineage_tfs=["RUNX2", "SATB2"])
annotate_cells(
    atlas,
    reference=SOLID_TUMOR_CELL_TYPE_SIGNATURES,
    overrides=[osteosarcoma],
)

# Melanoma
melanoma = tumor_override(lineage_tfs=["MITF", "SOX10"])

# Primary CTA rule only (no lineage-TF rescue)
cta_only = tumor_override()

# Custom stricter panel (testis-restricted CTAs only)
strict_melanoma = tumor_override(
    oncoref.CTA_testis_restricted_gene_names(),
    lineage_tfs=["MITF"],
)
```

`tumor_override` returns a `MarkerCountOverride` instance; construct
`MarkerCountOverride` directly when customizing `min_expr`, `count_col`, or
cluster fractions.

## Sub-clustering compartments with `refine_cluster`

Global Leiden clustering prioritizes broad lineage differences and often merges
finer populations within a compartment, such as pDCs and cDC1s inside myeloid
clusters, or γδ T cells co-clustering with cytotoxic αβ CD8 T cells.
`refine_cluster` isolates a subset of cells, re-computes Pearson-residual PCA
and Leiden clusters within that subset alone, and passes each sub-cluster to a
callback (`relabel_fn`) that assigns updated labels.

The re-clustering pipeline runs `counts → normalize_pearson_residuals →
nan_to_num → pca(n_comps=min(n_pcs, n_obs-1, n_vars-1)) →
neighbors(n_neighbors) → leiden(resolution, flavor="igraph", n_iterations=2,
directed=False)` with a fixed `random_state` and no UMAP step. When
`informative_genes` is supplied, that gene list is used without further
filtering.

TCRsift includes `gd_cd8_relabel` for separating γδ T cells from αβ CD8 T
cells:

```python
from tcrsift import gd_cd8_relabel, refine_cluster

# Populate a boolean obs column indicating paired/productive αβ VDJ capture:
s = atlas.obs["CDR3a/b"].astype(str)
atlas.obs["has_ab_contig"] = (
    (s != "/") & (s != "") & s.str.contains("[A-Z]", regex=True)
)

refine_cluster(
    atlas,                          # modified in place
    "CD8 T",                        # label in label_col, bool mask, index, or callable
    relabel_fn=gd_cd8_relabel(),    # mean TRDC >= 0.9 (log1p CP10K) and αβ capture < 0.25
    informative_genes=panel,
    resolution=0.5,
    n_pcs=30,
    min_cells=51,                   # skip compartments smaller than min_cells
    label_col="phenotype_base",
    final_cluster_col="final_cluster",
)
```

`gd_cd8_relabel` classifies a sub-cluster as γδ T (`("gd T", "gdt")`) only when
its sub-cluster mean `TRDC` expression in `.X` (on `log1p(CP10K)`) meets
`trdc_min` **and** its mean `has_ab_contig` fraction falls below `ab_max`;
averaging across the sub-cluster compensates for single-cell `TRDC` dropout.
Otherwise it returns `("CD8 effector/cytotoxic", "gdt_cd8")`, logging a warning
when αβ capture falls in an intermediate range. Gene symbol, `obs` column,
thresholds, and output labels are all configurable.

Custom `relabel_fn` callbacks receive a copy of one sub-cluster (`adata` with
original `.X`, `obs`, and local `obs["leiden"]`) and return a new label string,
a `(label, final_cluster)` tuple, or `None` to leave existing labels unchanged.

## Parameter and QC details

- **Threshold precedence in `cell_qc_funnel`.** Without a `config` argument,
  only `min_genes=250`, `min_counts=500`, and `max_mito_pct=8` are active
  (`max_genes`, `max_counts`, and `min_mito_pct` default to `None`). Passing a
  `LoadConfig` or top-level `TCRsiftConfig` supplies full PBMC bounds, and
  explicit keyword arguments override `config` values (passing `None` disables
  any individual bound). Any single step that removes more than `warn_frac`
  (default 10%) of entering cells logs a warning.
- **Mitochondrial percentage calculation.** If `adata.obs["pct_counts_mt"]` is
  already present, `cell_qc_funnel` leaves that column unchanged; however, the
  mitochondrial filter itself always evaluates mitochondrial fraction directly
  from raw counts.
- **Stimulated PBMC cultures.** For Ficoll-isolated PBMC cultures without solid
  stroma, pass `allowed_types=tcrsift.signatures.PBMC_CULTURE_TYPES` to
  `annotate_cells` so activated monocyte-derived DCs or macrophages are not
  assigned stromal labels.
