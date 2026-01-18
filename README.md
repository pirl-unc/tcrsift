# TCRsift

[![Tests](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml)
[![Documentation](https://github.com/pirl-unc/tcrsift/actions/workflows/docs.yml/badge.svg)](https://pirl-unc.github.io/tcrsift/)
[![Coverage Status](https://coveralls.io/repos/github/pirl-unc/tcrsift/badge.svg)](https://coveralls.io/github/pirl-unc/tcrsift)
[![PyPI version](https://badge.fury.io/py/tcrsift.svg)](https://pypi.org/project/tcrsift/)

Select antigen-specific TCRs from single-cell sequencing data.

```bash
pip install tcrsift
tcrsift run --sample-sheet samples.yaml -o results/
```

## Contents

- [Architecture](#architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Sample Sheet Format](#sample-sheet-format)
- [Core Pipeline Steps](#core-pipeline-steps)
- [Supplementary Tools](#supplementary-tools)
- [Workflows](#workflows)
- [API Reference](#api-reference)
- [Output Files](#output-files)

## Architecture

```
CellRanger VDJ + GEX  -->  tcrsift run  -->  clonotypes.csv
                              |
        load -> phenotype -> clonotype -> filter -> annotate -> assemble

Supplementary: load-sct, annotate-gex, match-til, unify
```

### Key Data Structures

| Stage | Format | Description |
|-------|--------|-------------|
| Load → Phenotype | `AnnData` | Per-cell data with expression matrix + VDJ annotations in `.obs` |
| Clonotype → Assemble | `DataFrame` | Per-clonotype data with aggregated statistics |

The transition from AnnData to DataFrame happens at the clonotype aggregation step. After that point, all operations work on clonotype-level DataFrames.

### CellRanger Requirements

TCRsift expects standard 10x Genomics CellRanger output directories:

**VDJ Directory** (from `cellranger vdj`):
```
vdj_outs/
├── filtered_contig_annotations.csv   # Required: contig info per cell
├── clonotypes.csv                    # Optional: CellRanger clonotype calls
├── consensus_annotations.csv         # Optional: for sequence assembly
└── filtered_contig.fasta             # Optional: for native leader extraction
```

Required columns in `filtered_contig_annotations.csv`:
- `barcode`, `chain` (TRA/TRB)
- `cdr3` (amino acid sequence)
- `v_gene`, `j_gene`, `c_gene`
- `umis`, `reads`
- `productive`, `full_length`

**GEX Directory** (from `cellranger count`):
```
gex_outs/
├── filtered_feature_bc_matrix.h5     # Preferred: HDF5 format
└── filtered_feature_bc_matrix/       # Alternative: MTX directory
    ├── matrix.mtx.gz
    ├── features.tsv.gz
    └── barcodes.tsv.gz
```

The GEX matrix must contain T cell marker genes (CD3D, CD3E, CD3G, CD4, CD8A, CD8B) for phenotyping. Gene names or ENSEMBL IDs are both supported.

**Barcode Matching**: CellRanger VDJ and GEX may use different barcode suffixes (e.g., `ACGT-1` vs `ACGT-2`). TCRsift strips suffixes and matches on the core barcode.

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

```bash
# For PDF report generation
pip install reportlab pdfkit
brew install wkhtmltopdf  # macOS

# For constant region sequences from Ensembl
pip install pyensembl
pyensembl install --release 93 --species human
```

## Quick Start

### Run the Complete Pipeline

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb
```

This runs: load → phenotype → clonotype → filter → annotate → assemble

### With Configuration File

```bash
# Generate example config with all defaults
tcrsift generate-config -o my_config.yaml

# Edit config, then run
tcrsift run --config my_config.yaml --sample-sheet samples.yaml -o results/
```

### Python API

```python
import tcrsift

# Load samples from sample sheet
adata = tcrsift.load_samples("samples.yaml")

# Phenotype cells (CD4/CD8 classification)
adata = tcrsift.phenotype_cells(adata)

# Aggregate to clonotypes
clonotypes = tcrsift.aggregate_clonotypes(adata)

# Filter by expansion
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
  # Culture sample with peptide stimulation
  - sample: "Patient1_Culture"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"
    mhc_allele: "HLA-A*02:01"
    source: "culture"

  # TIL sample (no antigen info needed)
  - sample: "Patient1_TIL"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
    tissue: "tumor"
```

### CSV Format

```csv
sample,vdj_dir,gex_dir,antigen_type,source
Patient1_Culture,/data/patient1/vdj,/data/patient1/gex,short_peptide,culture
Patient1_TIL,/data/patient1_til/vdj,,,til
```

### Required Fields

| Field | Required | Description |
|-------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `vdj_dir` | Yes* | Path to CellRanger VDJ output |
| `gex_dir` | No | Path to CellRanger GEX output |
| `source` | No | Sample type: `culture`, `til`, `tetramer`, `sct` |

*At least one of `vdj_dir` or `gex_dir` is required.

### Antigen Types

| Antigen Type | Expected T Cell | Description |
|--------------|-----------------|-------------|
| `short_peptide` | CD8 | 8-11aa peptides (direct MHC-I binding) |
| `long_peptide` | mixed | 15-25+aa (requires processing) |
| `whole_protein` | mixed | Full protein antigens |
| `tetramer_mhc1` | CD8 | MHC-I tetramer selection |
| `tetramer_mhc2` | CD4 | MHC-II tetramer selection |
| `sct` | CD8 | Single-chain trimer (pMHC-I fusion) |

## Core Pipeline Steps

### 1. Load Data

Loads CellRanger VDJ and GEX outputs, extracts T cell markers (CD3, CD4, CD8), and combines into a unified AnnData object.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

**What happens:**
1. Reads `filtered_contig_annotations.csv` for VDJ data
2. Reads `filtered_feature_bc_matrix.h5` for gene expression
3. Matches barcodes between VDJ and GEX
4. Extracts CD3D/E/G, CD4, CD8A/B expression per cell
5. Pivots VDJ to get one row per cell with TRA/TRB info

### 2. Phenotype Cells

Classifies each cell as CD4+ or CD8+ based on gene expression ratios.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --ratio 3.0
```

**Classification logic:**
- **Confident CD8+**: `(CD8A + CD8B + 1) / (CD4 + 1) > ratio` (default: 3.0)
- **Confident CD4+**: `(CD4 + 1) / (CD8A + CD8B + 1) > ratio`
- **Likely CD8+**: CD8 > 0 and CD4 = 0 (any CD8 without CD4)
- **Likely CD4+**: CD4 > 0 and CD8 = 0 (any CD4 without CD8)
- **Unknown**: Similar expression or both near zero

### 3. Aggregate Clonotypes

Groups cells by CDR3 sequences into clonotypes with aggregated statistics.

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

**Grouping options:**
- `CDR3ab`: Match by both alpha and beta CDR3 (strict pairing)
- `CDR3b_only`: Match by beta chain only (allows alpha variation)

**Output columns:**
- `clone_id`: Unique identifier (CDR3α_CDR3β)
- `cell_count`: Number of cells with this TCR
- `frequency`: Proportion of total cells
- `Tcell_type_consensus`: Most common phenotype
- `samples`: Which samples contain this clone

### 4. Filter Clonotypes

Applies tiered filtering to prioritize expanded clones.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

**Tier thresholds (default):**

| Tier | Min Cells | Min Frequency | Max Conditions |
|------|-----------|---------------|----------------|
| 1 | 10 | 1% | 2 |
| 2 | 5 | 0.5% | 3 |
| 3 | 3 | 0.1% | 5 |
| 4 | 2 | 0.05% | 10 |
| 5 | 2 | 0% | unlimited |

### 5. Annotate Clonotypes

Matches against public TCR databases to identify known specificities.

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb
```

**Supported databases:**
- **VDJdb**: Curated TCR-epitope pairs
- **IEDB**: Immune Epitope Database
- **CEDAR**: Cancer Epitope Database and Analysis Resource

**Viral flagging:** Clones matching CMV, EBV, HIV, Influenza, etc. are flagged as `is_viral=True` for review.

### 6. Assemble Full Sequences

Builds full-length TCR sequences with leader peptides and constant regions.

```bash
tcrsift assemble -i annotated.csv -o full_sequences.csv \
    --alpha-leader CD28 --beta-leader CD8A --include-constant
```

**Sequence structure:**
```
[Leader] + [V(D)J variable region] + [Constant region]

Single-chain construct:
[Beta full] + [T2A linker] + [Alpha full]
```

**Leader options:** CD8A, CD28, IgK, TRAC, TRBC, or `from_contig` (extract native)

## Supplementary Tools

These tools handle data outside the standard CellRanger workflow.

### Load SCT Data

Loads TCR data from SCT (single-cell TCR) platform Excel files.

```bash
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate
```

**When to use:** You have SCT platform data (pMHC tetramer with paired TCR sequencing) that wasn't processed through CellRanger.

**Quality filters applied:**
- `high_quality`: SNR ≥ 2.0, reads ≥ 10 per chain, mutation match
- `chosen`: Stricter criteria (SNR ≥ 3.4, reads ≥ 50, comPACT match)

### TIL Matching (Automatic)

When you include TIL samples in your sample sheet with `source: til`, the `run` command **automatically** detects them and adds TIL matching columns to culture clonotypes.

```yaml
# samples.yaml - TIL samples are auto-detected
samples:
  - sample: "Culture_Pool1"
    vdj_dir: "/data/culture/vdj"
    source: "culture"
  - sample: "Patient1_TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"  # This sample will be used for TIL matching
```

```bash
tcrsift run --sample-sheet samples.yaml -o results/
# TIL matching happens automatically - no extra flags needed!

# You can also explicitly specify TIL samples (overrides auto-detection):
tcrsift run --sample-sheet samples.yaml -o results/ --til-samples Patient1_TIL
```

**Output columns added:**
- `til_match` (bool): Clone found in TIL
- `til_cell_count`: Number of TIL cells with this TCR
- `til_frequency`: Frequency in TIL repertoire

**Why TIL matching matters:** Clones that appear in both antigen-stimulated culture AND tumor tissue provide orthogonal evidence of tumor-reactivity.

### Match TIL (Cross-Run)

Use `match-til` only when TIL data was processed in a **separate** pipeline run.

```bash
tcrsift match-til \
    -i culture_clonotypes.csv \
    --til-data til_processed.h5ad \
    -o matched.csv
```

**When to use:**
- TIL from a different patient or experiment
- Retrospective matching against archived TIL data
- TIL processed with different parameters

### Annotate with Gene Expression (`annotate-gex`)

Adds gene expression data from a 10x HDF5 file to TCR DataFrames.

**`annotate` vs `annotate-gex`:**
| Command | Data Source | Purpose |
|---------|-------------|---------|
| `annotate` | Public databases (VDJdb, IEDB) | Label clonotypes with known epitope specificities |
| `annotate-gex` | 10x HDF5 expression file | Add per-cell gene expression values |

**When GEX data is available:**
- **Standard pipeline**: If `gex_dir` is in your sample sheet, GEX is loaded automatically at the `load` step and used for CD4/CD8 phenotyping
- **VDJ-only workflows**: Use `annotate-gex` to add expression from a separate HDF5 file

**When to use `annotate-gex`:**
- You loaded VDJ-only data (no `gex_dir` in sample sheet)
- You have a separate 10x HDF5 file with expression data
- You want genes beyond the default CD3/CD4/CD8 markers

```bash
# Add per-cell expression from HDF5 file
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file filtered_feature_bc_matrix.h5 \
    -o cells_with_gex.csv

# Add GEX and aggregate to clonotype level
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file filtered_feature_bc_matrix.h5 \
    --aggregate \
    --cd4-cd8-counts \
    -o clonotype_gex.csv

# Custom gene list
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file matrix.h5 \
    --genes "GZMA,GZMB,PRF1,IFNG,TNF" \
    -o cytotoxicity_markers.csv
```

**Output columns:**
- `gex.{GENE}`: Expression per cell
- `gex.{GENE}.sum`, `gex.{GENE}.mean`: Aggregated per clonotype (with `--aggregate`)
- `gex.n_reads`, `gex.n_genes`, `gex.pct_mito`: QC metrics per cell
- `CD4_only.count`, `CD8_only.count`: Cells with exclusive expression (with `--cd4-cd8-counts`)

**Note:** For most workflows, use `gex_dir` in your sample sheet and the standard pipeline will handle GEX automatically during loading.

**Python API:**

```python
from tcrsift import augment_with_gex, aggregate_gex_by_clonotype, compute_cd4_cd8_counts

cells_df = augment_with_gex(cells_df, "filtered_feature_bc_matrix.h5")
clonotype_gex = aggregate_gex_by_clonotype(cells_df, group_col="CDR3_pair")
cd4_cd8 = compute_cd4_cd8_counts(cells_df, group_col="CDR3_pair")
```

### Unify Multiple Experiments

Merges clonotype data from multiple **independent** pipeline runs into a unified table.

```bash
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv sct_clonotypes.csv \
    -o unified.csv
```

**When to use:** You have results from multiple independent runs and want to compare or combine them.

**`run` vs `unify`:**
| Scenario | Use |
|----------|-----|
| One patient, culture + TIL in same sample sheet | `run` (TIL auto-detected) |
| One patient, culture + TIL processed separately | `match-til` |
| Multiple patients or experiments | `unify` |
| Comparing results across different data sources | `unify` |

**Output includes:**
- Prefixed columns from each source (e.g., `TIL.cell_count`, `Culture.cell_count`)
- Occurrence flags (`occurs_in_TIL`, `occurs_in_Culture`)
- Combined statistics (`combined.total_cells.count`)
- Phenotype confidence based on combined evidence

## Workflows

### Standard Single-Experiment Analysis

```bash
tcrsift run --sample-sheet samples.yaml -o results/ --vdjdb /path/to/vdjdb
```

### Culture + TIL Together

When TIL and culture samples are in the same sample sheet, TIL matching happens automatically:

```yaml
# samples.yaml
samples:
  - sample: "Culture_Pool1"
    vdj_dir: "/data/culture/vdj"
    gex_dir: "/data/culture/gex"
    source: "culture"
  - sample: "TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"  # Auto-detected for TIL matching
```

```bash
tcrsift run --sample-sheet samples.yaml -o results/
# No --til-samples flag needed - auto-detected from source: til
```

This automatically matches culture clones against TIL samples and adds `til_match`, `til_cell_count`, `til_frequency` columns.

### Multi-Source Unification

When combining data from different sources processed separately:

```bash
# Process each source
tcrsift run --sample-sheet til_samples.yaml -o til_results/
tcrsift run --sample-sheet culture_samples.yaml -o culture_results/
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate

# Unify
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv sct_clonotypes.csv \
    -o unified_clonotypes.csv
```

## API Reference

### Data Loading

```python
from tcrsift import load_samples, load_cellranger_vdj, load_cellranger_gex

# Load all samples from sample sheet
adata = load_samples("samples.yaml")

# Load individual CellRanger outputs
vdj_df = load_cellranger_vdj("/path/to/vdj", sample_name="S1")
adata = load_cellranger_gex("/path/to/gex", sample_name="S1")
```

### Phenotyping

```python
from tcrsift import phenotype_cells, filter_by_tcell_type, get_phenotype_summary

# Classify cells
adata = phenotype_cells(adata, cd4_cd8_ratio=3.0)

# Filter to CD8+ only
cd8_cells = filter_by_tcell_type(adata, tcell_type="cd8")

# Get summary by sample
summary = get_phenotype_summary(adata)
```

### Clonotyping

```python
from tcrsift import aggregate_clonotypes, get_clonotype_summary

# Aggregate by CDR3 pair
clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab", min_umi=2)

# Get summary
summary = get_clonotype_summary(clonotypes)
```

### Filtering

```python
from tcrsift import filter_clonotypes, split_by_tier

# Filter with default tiers
filtered = filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")

# Split into separate DataFrames by tier
tier_dfs = split_by_tier(filtered)
```

### Annotation

```python
from tcrsift import annotate_clonotypes, load_vdjdb

# Load database
vdjdb = load_vdjdb("/path/to/vdjdb")

# Annotate clonotypes
annotated = annotate_clonotypes(
    clonotypes,
    vdjdb_path="/path/to/vdjdb",
    match_by="CDR3ab",
    exclude_viral=True,
)
```

### TIL Matching

```python
from tcrsift import match_til, get_til_summary

# Match culture clones against TIL data
matched = match_til(culture_clonotypes, til_adata, match_by="CDR3ab")

# Get recovery statistics
summary = get_til_summary(matched)
```

### SCT Data

```python
from tcrsift import load_sct, aggregate_sct, get_sct_specificities

# Load and filter SCT data
df = load_sct("sct_data.xlsx", min_snr=2.0, min_reads_per_chain=10)
hq = df[df.high_quality]

# Aggregate to clonotypes
clonotypes = aggregate_sct(df)

# Get specificity mapping
specificities = get_sct_specificities(clonotypes)
```

### Multi-Experiment Unification

```python
from tcrsift import merge_experiments, add_phenotype_confidence

# Prepare experiments
experiments = [
    (til_clonotypes, "TIL"),
    (culture_clonotypes, "Culture"),
]

# Merge with occurrence flags and combined stats
unified = merge_experiments(experiments, add_occurrence_flags=True)

# Add phenotype confidence
unified = add_phenotype_confidence(unified, ratio_threshold=10.0)
```

### Sequence Assembly

```python
from tcrsift import assemble_full_sequences, export_fasta

# Assemble with leaders and constant regions
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader="CD28",
    beta_leader="CD8A",
    include_constant=True,
    linker="T2A",
)

# Export FASTA
export_fasta(assembled, "sequences.fasta", sequence_col="single_chain_aa")
```

## Output Files

### clonotypes.csv

| Column | Description |
|--------|-------------|
| `clone_id` | Unique identifier (CDR3α_CDR3β) |
| `CDR3_alpha` | Alpha chain CDR3 sequence |
| `CDR3_beta` | Beta chain CDR3 sequence |
| `cell_count` | Number of cells |
| `frequency` | Proportion of total cells |
| `Tcell_type_consensus` | Consensus T cell type |
| `tier` | Quality tier (1-5) |
| `db_match` | Matched in public database |
| `is_viral` | Viral specificity flag |

### full_sequences.csv

| Column | Description |
|--------|-------------|
| `clone_id` | Clonotype identifier |
| `alpha_full_aa` | Full alpha chain (leader + VDJ + constant) |
| `beta_full_aa` | Full beta chain |
| `single_chain_aa` | Beta-2A-Alpha construct |
| `single_chain_nt` | DNA sequence |

## Documentation

Full documentation: [https://pirl-unc.github.io/tcrsift/](https://pirl-unc.github.io/tcrsift/)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

## Citation

```bibtex
@software{tcrsift,
  author = {Rubinsteyn, Alex},
  title = {TCRsift: T-cell receptor selection from antigen-specific culture},
  url = {https://github.com/pirl-unc/tcrsift},
  year = {2024}
}
```
