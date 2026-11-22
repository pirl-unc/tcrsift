# TCRsift

Prioritize T-cell receptor clonotypes from paired 10x Genomics single-cell VDJ
and gene-expression data in **human** (*Homo sapiens*) and **mouse** (*Mus
musculus*) studies.

## Overview

TCRsift ranks candidate TCR clonotypes using clonal abundance, CD4/CD8
phenotype, gene-expression signatures, public-database matches (VDJdb, IEDB,
CEDAR), sequence generation probability, and optional overlap with matched
tumor-infiltrating lymphocytes (TIL). Output tables and full-length receptor
constructs are intended to narrow candidates for functional testing; TCRsift
does not predict epitope specificity from sequence or expression alone.

Standard 10x Genomics Cell Ranger directories serve as input:

- **VDJ output** (`cellranger vdj`): productive α/β contig annotations, CDR3
  sequences, and optional FASTA contigs for leader and constant-region assembly.
- **Gene-expression output** (`cellranger count`): feature-barcode count
  matrices for CD4/CD8 phenotyping, signature scoring, and cell-type annotation
  (optional for VDJ-only workflows).

For an overview of the data flow and selection algorithms, see
[How TCRsift works](how-it-works.md). For empirical evaluations on public
human cohorts, see [Benchmark on public TIL data](benchmarks.md) and
[Benchmark on infectious disease and blood](benchmarks-infectious.md).

## Choose a workflow

| Study design | Command | Guide |
| --- | --- | --- |
| Antigen-stimulated T-cell culture (with optional matched TIL) | `tcrsift run` | [Quick Start](getting-started/quickstart.md), [Culture Pipeline](user-guide/pipeline.md) |
| Fresh TIL, blood, pleural effusion, or marrow (paired VDJ + GEX) | `tcrsift prioritize --context ...` | [Clone Prioritization](user-guide/til-signatures.md) |
| Longitudinal TIL samples in the legacy v2 file layout | `tcrsift til-select` | [CLI Reference](user-guide/cli.md#tcrsift-til-select) |
| Single-cell GEX QC, batch-integrated embedding, and cell typing | `tcrsift cells` | [Single-Cell Atlas](user-guide/atlas.md) |

## Quick example

Install TCRsift from PyPI:

```bash
pip install tcrsift
```

Create a sample sheet pointing to your Cell Ranger output directories:

```yaml title="samples.yaml"
samples:
  - sample: "Patient1_Culture"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
```

Run the culture expansion pipeline:

```bash
tcrsift run --sample-sheet samples.yaml --output-dir results/
```

Or prioritize clones by expression state in a fresh tumor sample sheet:

```bash
tcrsift prioritize til_samples.yaml -o candidates/ --context solid-tumor
```

See [Installation](getting-started/installation.md) for optional extras,
[Sample Sheets](getting-started/sample-sheets.md) for all metadata fields, and
the [CLI Reference](user-guide/cli.md) or [API Reference](api/sample_sheet.md)
for individual commands and Python functions.
