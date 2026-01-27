# Pipeline Overview

TCRsift processes single-cell TCR data through a series of steps, each building on the previous.

## Pipeline Steps

```
┌─────────────────┐
│  Sample Sheet   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Load Data     │  ← CellRanger VDJ + GEX
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Phenotype     │  ← CD4/CD8 classification
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Clonotype     │  ← Aggregate by CDR3
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    Filter       │  ← Tiered selection
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Annotate      │  ← VDJdb, IEDB, CEDAR
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Match TIL     │  ← Optional TIL matching
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Assemble      │  ← Full-length sequences
└─────────────────┘
```

## 1. Load Data

Parses CellRanger VDJ and GEX outputs:

- **VDJ**: `filtered_contig_annotations.csv` for TCR sequences
- **GEX**: `filtered_feature_bc_matrix/` for gene expression

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

The output is an AnnData object with:

- Cell barcodes as observations
- Gene expression in `X`
- VDJ annotations in `obs` columns

## 2. Phenotype Cells

Classifies cells as CD4+ or CD8+ based on gene expression:

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --ratio 3.0
```

Classification categories:

| Category | Criteria |
|----------|----------|
| Confident CD8+ | CD8/CD4 ratio > threshold |
| Confident CD4+ | CD4/CD8 ratio > threshold |
| Likely CD8+ | CD8 > 0 and CD4 = 0 |
| Likely CD4+ | CD4 > 0 and CD8 = 0 |
| Unknown | Similar expression levels |

## 3. Aggregate Clonotypes

Groups cells by CDR3 sequences:

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

Grouping options:

- `CDR3ab`: Match by both alpha and beta chains (stricter)
- `CDR3b_only`: Match by beta chain only (more permissive)

## 4. Filter Clonotypes

Applies tiered filtering to prioritize antigen-specific clones:

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

See [Filtering Strategies](filtering.md) for detailed options.

## 5. Annotate Clonotypes

Matches against public TCR databases:

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb
```

Annotations include:

- Known epitope specificity
- Viral vs tumor antigens
- Database source

## 6. Match TIL (Optional)

For tumor studies, match culture clonotypes against TIL:

```bash
tcrsift match-til -i annotated.csv --til-csv til_clonotypes.csv -o matched.csv
```

This identifies clones that:

- Were expanded in culture AND present in tumor
- Are TIL-specific (not in culture)

TIL samples are excluded from culture aggregation/filtering and only used for matching.

## 7. Assemble Full Sequences

Build full-length TCR sequences:

```bash
tcrsift assemble -i annotated.csv -o sequences.csv \
    --include-constant \
    --linker T2A \
    --export-fasta sequences.fasta
```

Output includes:

- Leader peptide (from contigs)
- Variable region (VDJ)
- Constant region (from Ensembl)
- Single-chain construct (beta-T2A-alpha)

## Running the Complete Pipeline

Use `tcrsift run` to execute all steps:

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb \
    --tcell-type cd8 \
    --method threshold \
    # report generation is enabled by default; use --no-report to disable
```

This creates a complete output directory with all intermediate files and a summary report.
