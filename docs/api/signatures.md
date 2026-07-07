# Signatures API

Canonical T-cell gene-set signatures used across the per-sample
signature scatter (`tcrsift.plots`), the TIL-selection scoring
(`tcrsift.til_select`), and the signature-based candidate shortlist
(`tcrsift.candidate`). Pulled out of `til_select.py` so they can drive
non-TIL selections too — antigen-response screens, exhaustion-state
phenotyping, healthy-donor panels.

## Overview

All gene symbols are HGNC (human). Gene sets grouped by intent:

| Constant | Genes | Use |
| --- | --- | --- |
| `EFFECTOR_GENES_HGNC` | `IFNG, GZMB, PRF1, GNLY, NKG7` | Cytotoxic-effector panel. (`ACTIVATION_GENES_HGNC` is a deprecated alias — #142: this is effector differentiation, not immediate-early activation.) |
| `NAIVE_STEM_GENES_HGNC` | `TCF7, LEF1, CCR7, SELL, IL7R, CD27, CD28` | Naïve / stem-memory program — the "down" pole of the effector−naïve differentiation contrast (#141). |
| `ANTIGEN_RESPONSE_GENES_HGNC` | `TNFRSF9, MKI67` | AIM-assay marker (4-1BB / CD137) + Ki-67. Runs *inverse* to in-vitro clonal expansion at snapshot (#142). |
| `CYTOLYTIC_GENES_HGNC` | `PRF1, GZMB` | Canonical cytotoxic effector readout (Caushi 2021, Krishna 2021, Hanada 2022) |
| `EXHAUSTION_GENES_HGNC` | `PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4` | Canonical exhausted-T-cell surface markers |
| `TUMOR_REACTIVE_GENES_HGNC` | `CXCL13, ENTPD1` | TIL-resident tumor-reactive phenotype (Workel 2019, Duhen 2018) |

`T_CELL_SIGNATURES` is a snake-case-name → tuple dict for convenient iteration.

## Cell-type / state registries (atlas annotator)

The per-cell annotator (`tcrsift.annotate_cells`) argmaxes each Leiden cluster
over these registries. The defaults are **blood/PBMC-calibrated**; a solid-tumor
atlas should override `reference=` (see the [Single-Cell Atlas
Path](../user-guide/atlas.md) guide).

| Constant | What |
| --- | --- |
| `CELL_TYPE_SIGNATURES` | Default cell-type registry (immune + common stroma), blood/PBMC-calibrated. |
| `PBMC_CELL_TYPE_SIGNATURES` | Explicit alias for the default, naming its PBMC orientation. |
| `SOLID_TUMOR_CELL_TYPE_SIGNATURES` | Default **plus** solid-tissue lineages (mesothelial, osteoclast, skeletal muscle, adipocyte, Schwann/nerve). Parallels `tcrsift.qc.SOLID_TUMOR_LINEAGE_PROGRAMS`. Malignant cells stay `MarkerCountOverride` territory. |
| `T_STATE_SIGNATURES` / `B_STATE_SIGNATURES` | T-cell / B-cell sub-state registries (injectable via `t_state_reference=` / `b_state_reference=`). |
| `PBMC_CULTURE_TYPES` | Type allow-list for a PBMC/culture (`annotate_clusters(allowed_types=...)`). |

## Usage

```python
from tcrsift import signatures
from tcrsift.plots import plot_clone_freq_vs_signature_per_sample

# Use a focal signature for the scatter.
plot_clone_freq_vs_signature_per_sample(
    adata, clonotypes,
    gene_ids=signatures.ANTIGEN_RESPONSE_GENES_HGNC,
    signature_label="antigen-response",
    output_dir="figs/",
)

# Iterate over all five.
for name, genes in signatures.T_CELL_SIGNATURES.items():
    print(name, genes)
```

::: tcrsift.signatures
    options:
      members:
        - EFFECTOR_GENES_HGNC
        - ACTIVATION_GENES_HGNC
        - NAIVE_STEM_GENES_HGNC
        - ANTIGEN_RESPONSE_GENES_HGNC
        - CYTOLYTIC_GENES_HGNC
        - EXHAUSTION_GENES_HGNC
        - TUMOR_REACTIVE_GENES_HGNC
        - T_CELL_SIGNATURES
