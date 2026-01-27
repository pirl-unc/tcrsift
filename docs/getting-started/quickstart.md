# Quick Start

This guide walks you through a typical TCRsift workflow.

## Step 1: Prepare Your Sample Sheet

Create a YAML file describing your samples:

```yaml title="samples.yaml"
samples:
  - sample: "Patient1_CMV"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_description: "CMV pp65"
    source: "culture"

  - sample: "Patient1_TIL"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
```

## Step 2: Run the Pipeline

### Command Line (Recommended)

Run the complete pipeline with a single command:

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb
```

This will:

1. Load all samples
2. Phenotype cells as CD4+ or CD8+
3. Aggregate clonotypes
4. Apply tiered filtering
5. Annotate with VDJdb
6. Generate a summary report

### Python API

```python
import tcrsift

# Load sample sheet
sample_sheet = tcrsift.load_sample_sheet("samples.yaml")

# Load all samples into AnnData
adata = tcrsift.load_samples(sample_sheet)

# Phenotype cells
adata = tcrsift.phenotype_cells(adata, cd4_cd8_ratio=3.0)

# Aggregate clonotypes
clonotypes = tcrsift.aggregate_clonotypes(adata, group_by="CDR3ab")

# Filter clonotypes (default: CD8+ with threshold method)
filtered = tcrsift.filter_clonotypes(
    clonotypes,
    method="threshold",
    tcell_type="cd8",
)

# Annotate with public databases
annotated = tcrsift.annotate_clonotypes(
    filtered,
    vdjdb_path="/path/to/vdjdb",
    exclude_viral=True,
)

# Save results
annotated.to_csv("results/annotated_clonotypes.csv", index=False)
```

## Step 3: Explore Results

### Output Files

The pipeline creates:

```
results/
├── data/
│   ├── loaded.h5ad           # Raw loaded data
│   ├── phenotyped.h5ad       # With CD4/CD8 classification
│   ├── clonotypes.csv        # All clonotypes
│   ├── filtered_tier1.csv    # Highest confidence clones
│   ├── filtered_tier2.csv
│   ├── filtered_tier3.csv
│   ├── filtered_tier4.csv
│   ├── filtered_tier5.csv
│   ├── annotated.csv         # With database annotations (if provided)
│   ├── til_matched.csv       # TIL matching results (if TIL samples provided)
│   └── full_sequences.csv    # Assembled sequences (if assembly enabled)
├── plots/
│   ├── qc.pdf
│   ├── phenotype.pdf
│   ├── clonotypes.pdf
│   └── tcrsift_report.pdf    # Summary report (if enabled)
└── config.yaml               # Resolved config used for the run
```

### Key Columns

The output CSV files contain:

| Column | Description |
|--------|-------------|
| `CDR3ab` | Unique identifier (CDR3_alpha_CDR3_beta) |
| `CDR3_alpha` | Alpha chain CDR3 sequence |
| `CDR3_beta` | Beta chain CDR3 sequence |
| `cell_count` | Number of cells |
| `max_frequency` | Maximum frequency |
| `tier` | Quality tier (1 = best) |
| `Tcell_type_consensus` | CD4+ or CD8+ |
| `db_match` | Matched in public database |
| `is_viral` | Known viral specificity |

## Next Steps

- [Sample Sheet Format](sample-sheets.md) - Detailed sample sheet options
- [Pipeline Overview](../user-guide/pipeline.md) - Understanding each step
- [Filtering Strategies](../user-guide/filtering.md) - Customizing filters
