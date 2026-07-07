# Single-Cell Atlas Path

Alongside the TCR-selection pipeline, TCRsift ships a gene-expression **atlas
path** that QCs, embeds, and annotates an `AnnData` of single cells:

```
cell_qc_funnel  →  embed_cells  →  annotate_cells / annotate_clusters
  (per-cell QC       (Pearson-residual        (per-cluster cell-type +
   + doublet          PCA → Harmony →           T/B-state argmax with
   funnel)            UMAP → Leiden)            biology-aware gates)
```

These are library functions; on the command line they are the
[`tcrsift cells qc` / `cells embed` / `cells annotate`](cli.md#tcrsift-cells)
subcommands. Everything below runs on a raw-count `AnnData` — for example the
object returned by `tcrsift.load_sample`.

## End to end

```python
import scanpy as sc
from tcrsift import cell_qc_funnel, embed_cells, annotate_cells
from tcrsift.config import LoadConfig, EmbedConfig

# adata.X is RAW COUNTS (or provide adata.layers["counts"]).

# 1) Per-cell QC + doublet funnel. Bare defaults are gentle (only the
#    near-universal low-quality floors); the opinionated PBMC preset is one
#    explicit config away. A returned waterfall makes the cull visible.
adata, waterfall = cell_qc_funnel(adata, config=LoadConfig())
print(waterfall)  # [step, reason, removed, remaining] per QC/doublet class

# 2) Variance-stabilized embedding → Leiden clusters.
atlas = embed_cells(adata, EmbedConfig(batch_key="sample"))

# 3) Annotate: per-cluster cell-type + T/B-state argmax, written to
#    obs["phenotype_base"] and obs["phenotype"].
annotate_cells(atlas)
print(atlas.obs["phenotype"].value_counts())
```

## The default registries are blood/PBMC-oriented

`annotate_cells` argmaxes each cluster over `CELL_TYPE_SIGNATURES`,
`T_STATE_SIGNATURES`, and `B_STATE_SIGNATURES` from `tcrsift.signatures`. **These
defaults are calibrated for blood / PBMC data.** `CELL_TYPE_SIGNATURES` covers
the immune lineages and common stroma, but omits solid-tissue types a tumor
atlas needs (mesothelial, osteoclast, skeletal muscle, adipocyte, Schwann/nerve),
and its shared-type gene lists are PBMC-tuned. Running a solid-tumor atlas against
the default gives confidently wrong calls (e.g. tumor cells argmax'd to
Fibroblast) with no signal that the reference is blood-shaped.

Override `reference=` (and, for tissue-specific states, `t_state_reference=` /
`b_state_reference=`) — they are threaded to **both** the scorer and the
per-cluster reader, so a caller-defined type or state is consistently scored and
used rather than silently dropped:

```python
from tcrsift.signatures import SOLID_TUMOR_CELL_TYPE_SIGNATURES
from tcrsift.qc import SOLID_TUMOR_LINEAGE_PROGRAMS

# QC doublet gate for solid tumor (folds osteoclasts into myeloid, keeps
# fibroblast as the tumor+immune doublet handle):
adata, wf = cell_qc_funnel(adata, config=LoadConfig(),
                           lineage_sets=SOLID_TUMOR_LINEAGE_PROGRAMS)

# Annotation with the solid-tumor cell-type registry:
annotate_cells(atlas, reference=SOLID_TUMOR_CELL_TYPE_SIGNATURES)
```

`PBMC_CELL_TYPE_SIGNATURES` is an explicit alias for the default, so an override
reads symmetrically against `SOLID_TUMOR_CELL_TYPE_SIGNATURES`.

## Typing malignant cells

Malignant cells are **not** a signature in either registry — a tumor cluster
often shares collagen with fibroblasts and loses the argmax. Type it on positive
marker-count evidence instead, with `MarkerCountOverride` (a caller-supplied
antigen/marker panel; TCRsift only counts + relabels):

```python
from tcrsift.annotate_cells import MarkerCountOverride

tumor = MarkerCountOverride(
    "Tumor", gene_set=MY_CTA_PANEL, min_distinct=2, min_cluster_frac=0.4,
    # Rescue a low-coverage cluster: lineage-TF high AND >=1 marker in >=10% of
    # cells (rescue_min_distinct=1 gates on the broad any-marker fraction).
    rescue=(["RUNX2", "SATB2"], 1.0, 0.1, 1),
)
annotate_cells(atlas, reference=SOLID_TUMOR_CELL_TYPE_SIGNATURES, overrides=[tumor])
```

## Notes

- **Config**: `cell_qc_funnel(config=...)` accepts a top-level `TCRsiftConfig`
  (unwrapped to its `.load` section) or a flat `LoadConfig`; explicit kwargs win
  over `config` win over the built-in defaults. Passing `None` for any bound
  disables it.
- **`pct_counts_mt`**: the funnel leaves a caller's pre-existing `pct_counts_mt`
  (e.g. from `scanpy.pp.calculate_qc_metrics`) untouched and only computes it when
  absent; the mito gate runs off its own arithmetic either way.
- **Culture restriction**: for a Ficoll-isolated PBMC / peptide-stimulated
  culture, pass `allowed_types=tcrsift.signatures.PBMC_CULTURE_TYPES` so a
  cultured moDC/macrophage cluster can't win a stromal signature it shares an
  activation program with.
