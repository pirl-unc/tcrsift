# TCRsift

Prioritize T-cell receptor clonotypes from paired 10x Genomics single-cell VDJ
and gene-expression data in human (*Homo sapiens*) and mouse (*Mus
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
[How TCRsift works](how-it-works.md). For evaluations across 13 public human and
mouse datasets (15 cohort evaluations) and head-to-head comparisons against `stitchr`,
`DNAChisel`, `scirpy`, `muon` DSB, `OLGA`, `soNNia`, `Scanpy score_genes`, and
`pyucell` UCell, see the
[Benchmark overview](benchmarks-overview.md) and the cohort pages on
[public CD8 TIL data](benchmarks.md),
[infectious disease and blood](benchmarks-infectious.md),
[paired blood-tumor, post-ACT, and CD4 TILs](benchmarks-blood-cd4.md), and
[murine models (`--species mouse`)](benchmarks-mouse.md).

## Choose a workflow

| Study design | Command | Guide |
| --- | --- | --- |
| Fresh TIL, blood, pleural effusion, or marrow (paired VDJ + GEX) | `tcrsift prioritize samples.yaml -o out/ --context ... [--assemble]` | [Clone Prioritization](user-guide/til-signatures.md) |
| Antigen-stimulated T-cell culture or short restimulation (with optional matched TIL) | `tcrsift run -s samples.yaml -o out/ [--method stimulation]` | [Quick Start](getting-started/quickstart.md), [Culture Pipeline](user-guide/pipeline.md) |
| Standalone full-length TCR construct assembly from a clone CSV | `tcrsift assemble -i clones.csv -o seqs.csv` | [CLI Reference](user-guide/cli.md#tcrsift-assemble) |
| Single-cell GEX QC, batch-integrated embedding, and cell typing | `tcrsift cells {qc,embed,annotate}` | [Single-Cell Atlas](user-guide/atlas.md) |
| Longitudinal TIL samples in the legacy v2 file layout | `tcrsift til-select` | [CLI Reference](user-guide/cli.md#tcrsift-til-select) |

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

Or prioritize clones by expression state in a fresh tumor sample sheet and
assemble full-length `β-T2A-α` constructs in one command:

```bash
tcrsift prioritize til_samples.yaml -o candidates/ --context solid-tumor --assemble
```

See [Installation](getting-started/installation.md) for optional extras,
[Sample Sheets](getting-started/sample-sheets.md) for all metadata fields, and
the [CLI Reference](user-guide/cli.md) or [API Reference](api/sample_sheet.md)
for individual commands and Python functions.
