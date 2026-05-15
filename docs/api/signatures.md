# Signatures API

Canonical T-cell gene-set signatures used across the per-sample
signature scatter (`tcrsift.plots`), the TIL-selection scoring
(`tcrsift.til_select`), and the signature-based candidate shortlist
(`tcrsift.candidate`). Pulled out of `til_select.py` so they can drive
non-TIL selections too — antigen-response screens, exhaustion-state
phenotyping, healthy-donor panels.

## Overview

All gene symbols are HGNC (human). Five signatures grouped by intent:

| Constant | Genes | Use |
| --- | --- | --- |
| `ACTIVATION_GENES_HGNC` | `IFNG, GZMB, PRF1, GNLY, NKG7` | Broad effector activation panel |
| `ANTIGEN_RESPONSE_GENES_HGNC` | `TNFRSF9, MKI67` | AIM-assay activation marker (4-1BB / CD137) + Ki-67 proliferation |
| `CYTOLYTIC_GENES_HGNC` | `PRF1, GZMB` | Canonical cytotoxic effector readout (Caushi 2021, Krishna 2021, Hanada 2022) |
| `EXHAUSTION_GENES_HGNC` | `PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4` | Canonical exhausted-T-cell surface markers |
| `TUMOR_REACTIVE_GENES_HGNC` | `CXCL13, ENTPD1` | TIL-resident tumor-reactive phenotype (Workel 2019, Duhen 2018) |

`T_CELL_SIGNATURES` is a snake-case-name → tuple dict for convenient iteration.

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
        - ACTIVATION_GENES_HGNC
        - ANTIGEN_RESPONSE_GENES_HGNC
        - CYTOLYTIC_GENES_HGNC
        - EXHAUSTION_GENES_HGNC
        - TUMOR_REACTIVE_GENES_HGNC
        - T_CELL_SIGNATURES
