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

## Typing malignant cells (any cancer type)

Malignant cells are **not** a signature in either registry — a tumor cluster
often shares collagen with fibroblasts and loses the argmax, and there is no
reliable pan-cancer tumor expression program. Type it on positive marker-count
evidence instead. `tumor_override` wires the h37-style rule in one line:

- **primary** — relabel `Tumor` if a fraction of cells each express ≥2 distinct
  markers from your panel; **or**
- **rescue** — a lineage-TF is high **and** ≥1 marker shows in a fraction of
  cells (catches low-coverage tumor whose sparse markers dropped below the
  primary bar).

The marker panel defaults to **oncoref's** curated pan-cancer CTA set —
`oncoref.CTA_gene_names()` unioned with `oncoref.cta_clinical_target_gene_names()`
(canonical clinically-expressed CTAs like NY-ESO-1/CTAG2 and MAGEA11 that
testis-restriction filtering drops but are real tumor markers). CTAs are shared
across solid tumors, so the panel is reused across cancers and only the
lineage-TF rescue is tissue-specific:

```python
from tcrsift.annotate_cells import tumor_override

osteosarcoma = tumor_override(lineage_tfs=["RUNX2", "SATB2"])  # panel auto-filled
annotate_cells(atlas, reference=SOLID_TUMOR_CELL_TYPE_SIGNATURES,
               overrides=[osteosarcoma])
```

Adapting to another cancer is usually just swapping the rescue TFs (omit
`lineage_tfs` for the primary bar alone); override `gene_set` when you want a
tighter or tuned oncoref panel:

```python
tumor_override(lineage_tfs=["MITF", "SOX10"])                    # melanoma
tumor_override()                                                 # no rescue signal
# a tighter, testis-restricted panel:
import oncoref
tumor_override(oncoref.CTA_testis_restricted_gene_names(), lineage_tfs=["MITF"])
```

`tumor_override` returns a `MarkerCountOverride`, so you can still reach for the
class directly when you need finer control (custom `min_expr`, `count_col`, …).

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
