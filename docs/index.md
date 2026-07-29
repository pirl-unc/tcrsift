# TCRsift

Prioritize T-cell receptor clonotypes from paired single-cell VDJ and gene
expression data.

## Overview

TCRsift combines observed clonal expansion, T-cell phenotype and expression
state, public-database annotation, sequence publicness, and optional TIL
matching. These are prioritization signals: antigen specificity still requires
experimental validation.

TCRsift takes standard 10x Genomics Cell Ranger outputs as input:

- **VDJ output** from `cellranger vdj`, containing TCR contigs and clonotypes
- **Matching GEX output** from `cellranger count`, when the workflow uses
  CD4/CD8 phenotyping or gene-expression signatures

Using those inputs, it can:

- prioritize expanded clonotypes with biology-aware filtering
- score T-cell phenotype and published expression signatures
- annotate known database matches and sequence publicness
- find culture-enriched TCRs in matched tumor samples
- assemble full-length TCR sequences for selected candidates

VDJ-only workflows are also supported when expression-based analyses are not
needed. See [Input requirements](getting-started/sample-sheets.md) for the
expected files and sample-sheet fields.

## Quick example

Install TCRsift and create a minimal sample sheet pointing to the Cell Ranger
output directories:

```bash
pip install tcrsift
```

```yaml title="samples.yaml"
samples:
  - sample: "Patient1_Culture"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
```

Run the pipeline:

```bash
tcrsift run --sample-sheet samples.yaml --output-dir results/
```

The result directory contains per-cell data, clonotype tables, plots, and a
record of the resolved configuration. Database annotation and sequence assembly
outputs appear when their optional inputs are provided. `source` defaults to
`culture`. TIL-only analyses use `source: "til"` and the
[multi-sample TIL workflow](user-guide/til-signatures.md), not this
culture-oriented `run` command.

See the [Quick Start](getting-started/quickstart.md) for a complete example or
the [Python API](user-guide/pipeline.md) for step-by-step control.

## Choose a workflow

| Starting data | Start here |
| --- | --- |
| Antigen-stimulated culture, optionally with matched TIL | [Quick Start](getting-started/quickstart.md) |
| Multiple TIL VDJ + GEX samples | [Multi-sample TIL Prioritization](user-guide/til-signatures.md) |
| A single-cell atlas that needs QC, embedding, and cell typing | [Single-Cell Atlas Path](user-guide/atlas.md) |

## Installation

The quick example uses the core PyPI package. See
[Installation](getting-started/installation.md) for optional dependencies and
source installation.

## Next Steps

- [Quick Start Guide](getting-started/quickstart.md)
- [Multi-sample TIL Prioritization](user-guide/til-signatures.md)
- [Sample Sheet Format](getting-started/sample-sheets.md)
- [Pipeline Overview](user-guide/pipeline.md)
- [API Reference](api/sample_sheet.md)
