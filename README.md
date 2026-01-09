# TCRsift

[![Tests](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml)
[![Documentation](https://github.com/pirl-unc/tcrsift/actions/workflows/docs.yml/badge.svg)](https://pirl-unc.github.io/tcrsift/)
[![Coverage Status](https://coveralls.io/repos/github/pirl-unc/tcrsift/badge.svg)](https://coveralls.io/github/pirl-unc/tcrsift)
[![PyPI version](https://badge.fury.io/py/tcrsift.svg)](https://pypi.org/project/tcrsift/)

T-cell receptor selection for TCR-T studies from antigen-specific culture and scRNA/VDJ sequencing.

TCRsift is a comprehensive pipeline for identifying antigen-specific T cell receptor clones from single-cell sequencing data. It supports loading CellRanger outputs, CD4/CD8 phenotyping, clonotype aggregation, tiered filtering, annotation with public TCR databases, TIL matching, and full-length TCR sequence assembly.

## Installation

```bash
pip install tcrsift
```

Or install from source:

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
pip install -e .
```

### Optional Dependencies

For full functionality, install optional dependencies:

```bash
# For PDF report generation
pip install reportlab pdfkit
brew install wkhtmltopdf  # macOS

# For constant region sequences from Ensembl
pip install pyensembl
pyensembl install --release 93 --species human
```

## Quick Start

### Command Line Interface

Run the complete pipeline:

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb \
    --report
```

Or run individual steps:

```bash
# Load data
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad

# Phenotype cells
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad

# Aggregate clonotypes
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv

# Filter clonotypes
tcrsift filter -i clonotypes.csv -o filtered/

# Annotate with public databases
tcrsift annotate -i filtered/tier1.csv -o annotated.csv --vdjdb /path/to/vdjdb

# Assemble full-length sequences
tcrsift assemble -i annotated.csv -o full_sequences.csv --include-constant
```

### Python API

```python
import tcrsift

# Load sample sheet
sample_sheet = tcrsift.load_sample_sheet("samples.yaml")

# Load samples
adata = tcrsift.load_samples(sample_sheet)

# Phenotype cells
adata = tcrsift.phenotype_cells(adata)

# Aggregate clonotypes
clonotypes = tcrsift.aggregate_clonotypes(adata)

# Filter clonotypes
filtered = tcrsift.filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")

# Annotate with VDJdb
annotated = tcrsift.annotate_clonotypes(filtered, vdjdb_path="/path/to/vdjdb")

# Assemble full sequences
assembled = tcrsift.assemble_full_sequences(annotated, include_constant=True)
```

## Sample Sheet Format

TCRsift accepts sample sheets in CSV or YAML format.

### YAML Format

```yaml
samples:
  - sample: "Patient1_Culture"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_description: "CMV pp65"
    culture_days: 14
    source: "culture"

  - sample: "Patient1_TIL"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
    tissue: "tumor"
```

### CSV Format

```csv
sample,vdj_dir,gex_dir,antigen_type,antigen_description,source
Patient1_Culture,/data/patient1/vdj,/data/patient1/gex,short_peptide,CMV pp65,culture
Patient1_TIL,/data/patient1_til/vdj,,,til
```

### Sample Sheet Fields

| Field | Required | Description |
|-------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `vdj_dir` | Yes* | Path to CellRanger VDJ output |
| `gex_dir` | Yes* | Path to CellRanger GEX output |
| `antigen_type` | No | Type of antigen (see below) |
| `antigen_description` | No | Description of the antigen |
| `culture_days` | No | Duration of culture |
| `source` | No | Sample source: `culture`, `til`, `tetramer`, `sct` |
| `tcell_type_expected` | No | Expected T cell type: `CD4`, `CD8`, `mixed` |
| `pre_sorted` | No | Pre-sorting: `CD4`, `CD8` |
| `mhc_blocking` | No | MHC blocking: `MHC-I`, `MHC-II` |

*At least one of `vdj_dir` or `gex_dir` is required.

### Antigen Types

TCRsift uses antigen type to set biology-aware defaults:

| Antigen Type | Expected T Cell | Description |
|--------------|-----------------|-------------|
| `short_peptide` | CD8 | 8-11aa peptides, direct MHC-I binding |
| `long_peptide` | mixed | 15-25+aa peptides, requires processing |
| `minigene` | mixed | Minigene constructs |
| `whole_protein` | mixed | Full protein antigens |
| `tetramer_mhc1` | CD8 | MHC-I tetramer selection |
| `tetramer_mhc2` | CD4 | MHC-II tetramer selection |
| `sct_mhc1` | CD8 | Single-cell tetramer, MHC-I |
| `sct_mhc2` | CD4 | Single-cell tetramer, MHC-II |

## Pipeline Steps

### 1. Load Data (`tcrsift load`)

Loads CellRanger VDJ and GEX outputs for all samples.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

### 2. Phenotype Cells (`tcrsift phenotype`)

Classifies cells as CD4+ or CD8+ based on gene expression.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --ratio 3.0
```

Categories:
- **Confident CD8+**: CD8/CD4 ratio > threshold
- **Confident CD4+**: CD4/CD8 ratio > threshold
- **Likely CD8+**: CD8 > 0 and CD4 = 0
- **Likely CD4+**: CD4 > 0 and CD8 = 0
- **Unknown**: Similar expression levels

### 3. Aggregate Clonotypes (`tcrsift clonotype`)

Groups cells by CDR3 sequences into clonotypes.

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

Options:
- `--group-by CDR3ab`: Match by both alpha and beta chains (default)
- `--group-by CDR3b_only`: Match by beta chain only
- `--min-umi 2`: Minimum UMI count per chain
- `--handle-doublets flag|remove|keep-primary`: Handle cells with multiple chains

### 4. Filter Clonotypes (`tcrsift filter`)

Applies tiered filtering to identify antigen-specific clones.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

#### Threshold Method (default)

Uses configurable thresholds for cell count, frequency, and condition specificity:

| Tier | Min Cells | Min Frequency | Max Conditions |
|------|-----------|---------------|----------------|
| Tier 1 | 10 | 1% | 2 |
| Tier 2 | 5 | 0.5% | 3 |
| Tier 3 | 3 | 0.1% | 5 |
| Tier 4 | 2 | 0.05% | 10 |
| Tier 5 | 2 | 0% | unlimited |

#### Logistic Method

Uses logistic regression to adaptively set thresholds based on FDR:

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method logistic --fdr-tiers 0.0001,0.001,0.01,0.1,0.15
```

### 5. Annotate Clonotypes (`tcrsift annotate`)

Matches clonotypes against public TCR databases.

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb \
    --cedar /path/to/cedar
```

Options:
- `--match-by CDR3ab|CDR3b_only`: Matching strategy
- `--exclude-viral`: Remove clones matching viral epitopes
- `--flag-only`: Flag viral clones but don't remove

### 6. Match TIL (`tcrsift match-til`)

Matches culture clonotypes against TIL (tumor-infiltrating lymphocyte) data.

```bash
tcrsift match-til -i annotated.csv --til-data til_clonotypes.csv -o til_matched.csv
```

### 7. Assemble Full Sequences (`tcrsift assemble`)

Builds full-length TCR sequences with leader peptides and constant regions.

```bash
tcrsift assemble -i annotated.csv -o full_sequences.csv \
    --include-constant \
    --contigs-dir /path/to/contigs \
    --linker T2A
```

Options:
- `--include-constant`: Add constant regions from Ensembl
- `--include-leader`: Extract leader peptides from contigs
- `--linker T2A`: Add T2A linker for single-chain constructs
- `--export-fasta sequences.fasta`: Export to FASTA format

## API Reference

### Sample Sheet

```python
from tcrsift import Sample, SampleSheet, load_sample_sheet, validate_sample_sheet

# Load sample sheet
sample_sheet = load_sample_sheet("samples.yaml")

# Access samples
for sample in sample_sheet:
    print(sample.sample, sample.get_expected_tcell_type())

# Get specific sample types
culture_samples = sample_sheet.get_culture_samples()
til_samples = sample_sheet.get_til_samples()
tetramer_samples = sample_sheet.get_tetramer_samples()

# Validate
warnings = validate_sample_sheet(sample_sheet)
```

### Data Loading

```python
from tcrsift import load_cellranger_vdj, load_cellranger_gex, load_sample, load_samples

# Load individual outputs
vdj_df = load_cellranger_vdj("/path/to/vdj")
adata = load_cellranger_gex("/path/to/gex")

# Load single sample
adata = load_sample(sample)

# Load all samples
adata = load_samples(sample_sheet)
```

### Phenotyping

```python
from tcrsift import phenotype_cells, classify_tcell_type, filter_by_tcell_type, get_phenotype_summary

# Classify cells
adata = phenotype_cells(adata, cd4_cd8_ratio=3.0)

# Filter to specific type
cd8_cells = filter_by_tcell_type(adata, tcell_type="cd8")

# Get summary
summary = get_phenotype_summary(adata)
```

### Clonotyping

```python
from tcrsift import aggregate_clonotypes, get_clonotype_summary, export_clonotypes_airr

# Aggregate clonotypes
clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab", min_umi=2)

# Get summary
summary = get_clonotype_summary(clonotypes)

# Export AIRR format
export_clonotypes_airr(clonotypes, "clonotypes.airr.tsv")
```

### Filtering

```python
from tcrsift import filter_clonotypes, filter_clonotypes_threshold, assign_tiers_threshold, split_by_tier, get_filter_summary

# Filter with default settings
filtered = filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")

# Custom tier definitions
custom_tiers = {
    "tier1": {"min_cells": 20, "min_frequency": 0.02, "max_conditions": 1},
    "tier2": {"min_cells": 10, "min_frequency": 0.01, "max_conditions": 2},
}
filtered = filter_clonotypes(clonotypes, tier_definitions=custom_tiers)

# Split by tier
tier_dfs = split_by_tier(filtered)
for tier, df in tier_dfs.items():
    df.to_csv(f"{tier}.csv")
```

### Annotation

```python
from tcrsift import load_vdjdb, load_iedb, load_cedar, annotate_clonotypes, get_annotation_summary

# Load databases
vdjdb = load_vdjdb("/path/to/vdjdb")
iedb = load_iedb("/path/to/iedb")

# Annotate clonotypes
annotated = annotate_clonotypes(
    clonotypes,
    vdjdb_path="/path/to/vdjdb",
    iedb_path="/path/to/iedb",
    match_by="CDR3ab",
    exclude_viral=True,
)

# Get summary
summary = get_annotation_summary(annotated)
```

### TIL Matching

```python
from tcrsift import match_til, get_til_summary, identify_til_specific_clones

# Match against TIL data
matched = match_til(clonotypes, til_clonotypes)

# Get recovery statistics
summary = get_til_summary(matched)

# Find TIL-specific clones
til_specific = identify_til_specific_clones(matched)
```

### Sequence Assembly

```python
from tcrsift import assemble_full_sequences, translate_dna, validate_sequences, export_fasta

# Assemble full sequences
assembled = assemble_full_sequences(
    clonotypes,
    contigs_dir="/path/to/contigs",
    include_leader=True,
    include_constant=True,
    linker="T2A",
)

# Validate sequences
warnings = validate_sequences(assembled)

# Export FASTA
export_fasta(assembled, "sequences.fasta", sequence_col="single_chain_aa")
```

### Utilities

```python
from tcrsift import tcr_name

# Generate mnemonic name for TCR
name = tcr_name("CASSLGQAYEQYF", "CAVSDGGSQGNLIF")
```

## Output Files

### Clonotypes CSV

| Column | Description |
|--------|-------------|
| `clone_id` | Unique clonotype identifier (CDR3a_CDR3b) |
| `CDR3_alpha` | Alpha chain CDR3 sequence |
| `CDR3_beta` | Beta chain CDR3 sequence |
| `cell_count` | Number of cells in clonotype |
| `samples` | Samples where clone is found |
| `n_samples` | Number of samples |
| `max_frequency` | Maximum frequency in any sample |
| `Tcell_type_consensus` | Consensus T cell type |
| `tier` | Quality tier (after filtering) |
| `db_match` | Whether matched in public database |
| `is_viral` | Whether specificity is viral |

### AIRR Format

TCRsift can export clonotypes in [AIRR format](https://docs.airr-community.org/) for interoperability:

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --airr clonotypes.airr.tsv
```

## Documentation

Full documentation is available at [https://pirl-unc.github.io/tcrsift/](https://pirl-unc.github.io/tcrsift/)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

## Citation

If you use TCRsift in your research, please cite:

```bibtex
@software{tcrsift,
  author = {Rubinsteyn, Alex},
  title = {TCRsift: T-cell receptor selection from antigen-specific culture},
  url = {https://github.com/pirl-unc/tcrsift},
  year = {2024}
}
```
