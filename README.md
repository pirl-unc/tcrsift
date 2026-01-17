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
    --vdjdb /path/to/vdjdb
```

With a configuration file (recommended for reproducibility):

```bash
# Generate example config
tcrsift generate-config -o my_config.yaml

# Edit config, then run
tcrsift run --config my_config.yaml --sample-sheet samples.yaml -o results/
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
  # Minimal peptide culture (antigen == epitope)
  - sample: "Patient1_CMV"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"  # same as antigen for minimal peptides
    mhc_allele: "HLA-A*02:01"
    culture_days: 14
    source: "culture"

  # Whole protein culture (antigen >> epitope)
  - sample: "Patient1_Protein"
    vdj_dir: "/data/patient1_protein/vdj"
    gex_dir: "/data/patient1_protein/gex"
    antigen_type: "whole_protein"
    antigen_name: "PRAME"
    # epitope_sequence unknown - will be processed by APCs
    culture_days: 21
    source: "culture"

  # Peptide pool stimulation
  - sample: "Patient1_Pool"
    vdj_dir: "/data/patient1_pool/vdj"
    gex_dir: "/data/patient1_pool/gex"
    antigen_type: "peptide_pool"
    antigen_names:  # required when >1 antigen
      - "KRAS_G12D"
      - "TP53_R175H"
      - "PIK3CA_H1047R"
    antigen_sequences:  # optional but helpful
      - "GADGVGKSAL"
      - "HMTEVVRHC"
      - "ARHGGWTTKM"
    culture_days: 14
    source: "culture"

  # SCT selection (epitope is known from the construct)
  - sample: "Patient1_SCT"
    vdj_dir: "/data/patient1_sct/vdj"
    antigen_type: "sct"
    antigen_name: "PRAME"  # source protein
    epitope_sequence: "SLLQHLIGL"  # what's in the SCT
    mhc_allele: "HLA-A*02:01"
    source: "sct"

  # TIL sample (no antigen info needed)
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
| `antigen_name` | No | Name of source antigen (e.g., "PRAME", "CMV pp65") |
| `antigen_sequence` | No | Sequence of source antigen (may be long) |
| `epitope_sequence` | No | Minimal peptide AA sequence that binds MHC |
| `mhc_allele` | No | MHC restriction (e.g., "HLA-A*02:01") |
| `antigen_names` | No | List of source antigen names (for pools/libraries) |
| `antigen_sequences` | No | List of source antigen sequences (for pools/libraries) |
| `epitope_sequences` | No | List of minimal epitope sequences (for pools, if known) |

*At least one of `vdj_dir` or `gex_dir` is required.

**Antigen vs Epitope:**

- **Antigen** = what you gave to APCs (whole protein, long peptide, minigene, mRNA, etc.)
- **Epitope** = the minimal peptide that sits in the MHC groove (8-11aa for MHC-I, 13-25aa for MHC-II)

For minimal peptide stimulation, antigen == epitope. For whole proteins, the antigen is much larger than the processed epitope. Only the epitope can have an MHC restriction.

**Field usage:**

- **Single antigen**: Use `antigen_name` and optionally `antigen_sequence`. If you know the minimal epitope, add `epitope_sequence` and `mhc_allele`.

- **Tetramer/SCT**: The epitope is known (it's in the tetramer). Provide `epitope_sequence`, `mhc_allele`, and `antigen_name` to describe the source (e.g., "HLA-A*02:01/SLLQHLIGL from PRAME").

- **Pools/libraries**: Use `antigen_names` (required if >1 antigen) and optionally `antigen_sequences`. Add `epitope_sequences` if you know the minimal epitopes.

**Naming rules**: At least a name or sequence must be provided. If only a sequence is given, it becomes the name. If multiple sequences are given without names, that's ambiguous.

### Antigen Types

TCRsift uses antigen type to set biology-aware defaults:

| Antigen Type | Expected T Cell | Description |
|--------------|-----------------|-------------|
| `short_peptide` | CD8 | 8-11aa peptides, direct MHC-I binding |
| `long_peptide` | mixed | 15-25+aa peptides, requires processing |
| `peptide_pool` | mixed | Pool of peptides for stimulation |
| `minigene` | mixed | Single minigene construct |
| `minigene_library` | mixed | Library of multiple minigene constructs |
| `whole_protein` | mixed | Full protein antigens |
| `mrna` | mixed | mRNA encoding one or more antigens |
| `tetramer_mhc1` | CD8 | MHC-I tetramer selection (single antigen) |
| `tetramer_mhc2` | CD4 | MHC-II tetramer selection (single antigen) |
| `sct` | CD8 | Single-chain trimer (pMHC-I: alpha-B2M-peptide fusion) |

## Pipeline Steps

### 1. Load Data (`tcrsift load`)

Loads CellRanger VDJ and GEX outputs for all samples.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

**How it works:**

1. **VDJ Loading**: Parses `filtered_contig_annotations.csv` from CellRanger VDJ output to extract:
   - Cell barcodes and contig sequences
   - V(D)J gene assignments (TRAV, TRAJ, TRBV, TRBJ, etc.)
   - CDR3 amino acid and nucleotide sequences
   - UMI counts per chain

2. **GEX Loading**: Uses scanpy to read `filtered_feature_bc_matrix.h5` containing:
   - Sparse expression matrix (cells × genes)
   - Cell barcodes (must match VDJ barcodes after suffix adjustment)
   - Gene metadata

3. **Barcode Matching**: CellRanger outputs may have different barcode suffixes between VDJ and GEX. TCRsift strips suffixes (e.g., `-1`) and matches cells by their core barcode.

4. **Quality Filtering**: Cells are filtered based on:
   - Minimum/maximum genes detected
   - Minimum/maximum UMI counts
   - Maximum mitochondrial percentage (indicative of dying cells)

### 2. Phenotype Cells (`tcrsift phenotype`)

Classifies cells as CD4+ or CD8+ based on gene expression ratios.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --ratio 3.0
```

**How it works:**

1. **Expression Extraction**: Reads CD4 and CD8 (CD8A+CD8B) gene expression from the GEX matrix.

2. **Ratio Calculation**: For each cell, computes:
   - CD8 signal = CD8A + CD8B expression
   - Ratio = CD8 / (CD4 + 1) and CD4 / (CD8 + 1)

3. **Classification Logic**:
   - **Confident CD8+**: CD8/CD4 ratio > threshold (default: 3.0)
   - **Confident CD4+**: CD4/CD8 ratio > threshold
   - **Likely CD8+**: CD8 > 0 and CD4 = 0 (any CD8 without CD4)
   - **Likely CD4+**: CD4 > 0 and CD8 = 0 (any CD4 without CD8)
   - **Unknown**: Similar expression levels or both near zero

4. **Confidence Levels**: For downstream unification, assigns:
   - `Confident_CD4/CD8`: High ratio evidence
   - `Likely_CD4/CD8`: Some evidence but not definitive

### 3. Aggregate Clonotypes (`tcrsift clonotype`)

Groups cells by CDR3 sequences into clonotypes with aggregated statistics.

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

**How it works:**

1. **Clone Key Construction**: Creates a unique identifier for each clonotype:
   - `CDR3ab`: `{CDR3_alpha}/{CDR3_beta}` (strict pairing)
   - `CDR3b_only`: Just CDR3_beta (allows alpha chain variation)

2. **Grouping**: Cells with identical clone keys are grouped together.

3. **Aggregation Statistics** computed per clonotype:
   - `cell_count`: Total cells with this TCR
   - `frequency`: Proportion of total cells
   - `samples`: Which samples contain this clone
   - `Tcell_type_consensus`: Most common phenotype among cells
   - V/J gene consensus (most frequent assignment)
   - UMI statistics (min, max, median)

4. **Doublet Handling**: Cells with multiple productive chains:
   - `flag`: Mark as doublet but keep
   - `remove`: Exclude from analysis
   - `keep-primary`: Keep the most abundant chain pair

**Options:**
- `--group-by CDR3ab`: Match by both alpha and beta chains (default)
- `--group-by CDR3b_only`: Match by beta chain only
- `--min-umi 2`: Minimum UMI count per chain
- `--handle-doublets flag|remove|keep-primary`: Handle cells with multiple chains

### 4. Filter Clonotypes (`tcrsift filter`)

Applies tiered filtering to identify antigen-specific clones.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

**How it works:**

1. **T Cell Type Filtering**: Optionally restricts to CD4+ or CD8+ clonotypes based on consensus phenotype.

2. **Tier Assignment**: Each clonotype is assigned to the highest tier whose criteria it meets.

3. **Threshold Method (default)**:
   Uses configurable thresholds based on three metrics:
   - **Cell count**: More cells = more confident
   - **Frequency**: Higher frequency = more likely antigen-specific
   - **Condition specificity**: Clones in fewer conditions are more specific

| Tier | Min Cells | Min Frequency | Max Conditions |
|------|-----------|---------------|----------------|
| Tier 1 | 10 | 1% | 2 |
| Tier 2 | 5 | 0.5% | 3 |
| Tier 3 | 3 | 0.1% | 5 |
| Tier 4 | 2 | 0.05% | 10 |
| Tier 5 | 2 | 0% | unlimited |

4. **Logistic Method** (advanced):
   Fits a logistic regression model to predict whether a clone is likely antigen-specific based on:
   - Log cell count
   - Log frequency
   - Number of conditions

   FDR-based tiers are assigned using model-predicted probabilities.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method logistic --fdr-tiers 0.0001,0.001,0.01,0.1,0.15
```

### 5. Annotate Clonotypes (`tcrsift annotate`)

Matches clonotypes against public TCR databases to identify known specificities.

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb \
    --cedar /path/to/cedar
```

**How it works:**

1. **Database Loading**: Parses VDJdb, IEDB, and CEDAR formats into a unified structure with:
   - CDR3α and CDR3β sequences
   - Epitope/antigen information
   - Species (for viral flagging)

2. **Matching Strategy**:
   - `CDR3ab`: Exact match on both chains (most specific)
   - `CDR3b_only`: Match on beta chain only (more permissive)

3. **Viral Flagging**: Clonotypes matching viral epitopes (CMV, EBV, HIV, Influenza, etc.) are flagged as `is_viral=True`.

4. **Exclusion Options**:
   - `--exclude-viral`: Remove viral-specific clones from output
   - `--flag-only`: Keep viral clones but mark them for review

**Options:**
- `--match-by CDR3ab|CDR3b_only`: Matching strategy
- `--exclude-viral`: Remove clones matching viral epitopes
- `--flag-only`: Flag viral clones but don't remove

### 6. Match TIL (`tcrsift match-til`)

Matches culture clonotypes against TIL (tumor-infiltrating lymphocyte) data.

```bash
tcrsift match-til -i annotated.csv --til-data til_clonotypes.csv -o til_matched.csv
```

**How it works:**

1. **TIL Reference Building**: Creates a lookup of CDR3 sequences found in TIL samples.

2. **Matching**: Each culture clonotype is checked against the TIL reference:
   - `til_match`: Boolean indicating presence in TIL
   - `til_frequency`: Frequency of the clone in TIL (if present)
   - `til_cell_count`: Number of TIL cells with this clone

3. **Enrichment Analysis**: Calculates whether culture-enriched clones are over-represented in TIL vs. expected by chance.

### 7. Load Amplify Data (`tcrsift load-amplify`)

Loads TCR data from the Amplify single-cell platform (Excel format).

```bash
tcrsift load-amplify -i amplify_data.xlsx -o amplify_clonotypes.csv --aggregate
```

**How it works:**

1. **Excel Parsing**: Reads the "Cell" sheet from Amplify output containing:
   - CDR3α/β sequences
   - Signal-to-noise ratio (SNR)
   - Read counts per chain
   - Mutation/specificity calls (PE and APC channels)

2. **Quality Filtering**:
   - `high_quality`: SNR ≥ 2.0, reads ≥ 10 per chain, mutation match
   - `chosen`: Stricter criteria (SNR ≥ 3.4, reads ≥ 50, comPACT match)

3. **Standardization**: Adds standardized columns (`CDR3_pair`, `CDR3_alpha`, `CDR3_beta`) for compatibility with other TCRsift functions.

4. **Aggregation** (optional): Groups by CDR3 pair and computes:
   - Min/median/max for numeric columns
   - Any/all for boolean columns
   - Consistent mutation assignment

### 8. Unify Experiments (`tcrsift unify`)

Merges clonotype data from multiple experiments into a unified table.

```bash
tcrsift unify -i til_clonotypes.csv culture_clonotypes.csv amplify_clonotypes.csv \
    -o unified.csv
```

**How it works:**

1. **Prefix Addition**: Each input file's columns are prefixed with the experiment name (derived from filename), except for key columns (`CDR3_pair`, `CDR3_alpha`, `CDR3_beta`).

2. **Outer Merge**: All clonotypes are combined with outer join, preserving clonotypes unique to each experiment.

3. **Occurrence Flags**: Boolean columns `occurs_in_*` indicate which experiments contain each clonotype:
   ```
   occurs_in_TIL, occurs_in_Culture, occurs_in_Amplify
   ```

4. **Combined Statistics**: Aggregates metrics across experiments:
   - `combined.total_cells.count`: Sum of cell counts
   - `combined.gex.CD4.sum`: Sum of CD4 expression
   - `combined.CD4_only.frac`: Fraction of CD4-only cells

5. **Phenotype Confidence**: Uses combined evidence for more confident CD4/CD8 calls:
   - `Confident_CD4/CD8`: Strong expression ratio (default: 10×)
   - `Likely_CD4/CD8`: Any evidence from TIL sorting or expression

### 9. GEX Augmentation (Python API)

Adds gene expression data to TCR DataFrames for downstream analysis.

```python
from tcrsift import augment_with_gex, aggregate_gex_by_clonotype

# Add per-cell expression
cells = augment_with_gex(cells_df, "filtered_feature_bc_matrix.h5")

# Aggregate to clonotype level
clonotypes = aggregate_gex_by_clonotype(cells, operations=["sum", "mean"])
```

**How it works:**

1. **HDF5 Loading**: Uses scanpy to read the 10x-formatted expression matrix.

2. **Gene Selection**: Extracts expression for specified genes (default: T cell markers like CD3D, CD4, CD8A, etc.) and QC metrics.

3. **Gene Groups**: Computes aggregate scores for gene groups:
   - `gex.CD3` = mean(CD3D, CD3E, CD3G)
   - `gex.CD8` = mean(CD8A, CD8B)

4. **QC Metrics** per cell:
   - `gex.n_reads`: Total UMI counts
   - `gex.n_genes`: Genes detected (>0)
   - `gex.pct_mito`: Mitochondrial percentage

5. **Clonotype Aggregation**: When aggregating by clonotype:
   - `gex.CD4.sum`: Total CD4 expression across cells
   - `gex.CD4.mean`: Average CD4 expression per cell
   - `total_cells.count`: Number of cells in clonotype

### 10. Assemble Full Sequences (`tcrsift assemble`)

Builds full-length TCR sequences with leader peptides and constant regions.

```bash
# Default: CD28 on alpha, CD8A on beta (distinct leaders)
tcrsift assemble -i annotated.csv -o full_sequences.csv --include-constant

# Custom leaders
tcrsift assemble -i annotated.csv -o full_sequences.csv \
    --alpha-leader CD8A --beta-leader CD28 --linker P2A

# No leaders
tcrsift assemble -i annotated.csv -o full_sequences.csv --no-leaders

# Extract native leaders from contig FASTAs
tcrsift assemble -i annotated.csv -o full_sequences.csv \
    --leaders-from-contigs --contigs-dir /path/to/contigs
```

Options:
- `--alpha-leader`: Leader for alpha chain (default: CD28)
- `--beta-leader`: Leader for beta chain (default: CD8A)
- `--no-leaders`: Disable leaders on both chains
- `--leaders-from-contigs`: Extract native leaders from FASTA (requires `--contigs-dir`)
- `--include-constant`: Add constant regions from Ensembl (TRAC, TRBC1, TRBC2)
- `--linker T2A`: Linker for single-chain constructs (T2A, P2A, E2A, F2A)
- `--fasta sequences.fasta`: Export to FASTA format

**Available Leader Sequences:**

| Leader | Source | Species | Use Case |
|--------|--------|---------|----------|
| CD8A | CD8A signal peptide | Human | Common choice for TCR expression |
| CD28 | CD28 signal peptide | Human | Alternative signal peptide |
| IgK | IgGκ light chain | Mouse | High secretion efficiency |
| TRAC | TCR alpha constant | Human | Native-like expression |
| TRBC | TCR beta constant | Human | Native-like expression |

By default, TCRsift uses **distinct leaders** on each chain (CD28 on alpha, CD8A on beta) to facilitate identification in downstream applications.

**How it works:**

1. **VDJ Sequence Extraction**: Gets variable region sequences from input clonotypes:
   - From `VDJ_alpha_aa`/`VDJ_beta_aa` columns if present
   - Or from contig FASTA files via barcode/contig ID lookup

2. **Leader Selection**: For each chain, the leader sequence is:
   - None (no leader)
   - A named signal peptide (CD8A, CD28, IgK, TRAC, TRBC)
   - Extracted from the contig sequence (native leader)

3. **Constant Region Lookup**: If `include_constant=True`:
   - TRAC for alpha chain
   - TRBC1 or TRBC2 for beta chain (based on C gene annotation)
   - Sequences fetched from Ensembl via pyensembl or built-in sequences

4. **Full Chain Assembly**: For each chain:
   ```
   [Leader] + [V(D)J variable region] + [Constant region]
   ```

5. **Single-Chain Construct**: Joins chains with 2A peptide linker:
   ```
   [Beta full chain] + [2A linker] + [Alpha full chain]
   ```

   Note: Beta is placed first because it's typically more important for antigen recognition and the first ORF after the promoter is expressed at higher levels.

6. **Nucleotide Back-Translation**: DNA sequences are assembled using codon tables, preserving original codons where possible.

## Typical Workflows

### Single Experiment Analysis

For analyzing a single antigen-specific culture:

```bash
tcrsift run --sample-sheet samples.yaml -o results/ --vdjdb /path/to/vdjdb
```

This runs the full pipeline: load → phenotype → clonotype → filter → annotate → assemble.

### Multi-Experiment Unification

For combining data from multiple sources (TIL, culture, Amplify):

```bash
# Step 1: Process each experiment separately
tcrsift run --sample-sheet til_samples.yaml -o til_results/
tcrsift run --sample-sheet culture_samples.yaml -o culture_results/
tcrsift load-amplify -i amplify.xlsx -o amplify_clonotypes.csv --aggregate

# Step 2: Unify all experiments
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv amplify_clonotypes.csv \
    -o unified_clonotypes.csv
```

### With GEX Augmentation (Python)

For detailed gene expression analysis:

```python
import tcrsift

# Load and process
cells = tcrsift.load_cellranger_vdj("vdj_dir")
cells = tcrsift.augment_with_gex(cells, "gex_dir/filtered_feature_bc_matrix.h5")

# Aggregate with expression data
clonotypes = tcrsift.aggregate_gex_by_clonotype(cells)

# Compute CD4/CD8 counts
cd4_cd8 = tcrsift.compute_cd4_cd8_counts(cells)
```

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

### Amplify Data

```python
from tcrsift import load_amplify, aggregate_amplify, get_amplify_specificities

# Load Amplify data
df = load_amplify("amplify_data.xlsx", min_snr=2.0, min_reads_per_chain=10)

# Filter to high-quality cells
hq = df[df.high_quality]

# Aggregate to clonotypes
clonotypes = aggregate_amplify(df)

# Get specificity mapping
specificities = get_amplify_specificities(clonotypes)
# Returns: {"CAVXXX/CASSYYY": "KRAS_G12D", ...}
```

### GEX Augmentation

```python
from tcrsift import augment_with_gex, aggregate_gex_by_clonotype, compute_cd4_cd8_counts

# Add gene expression to per-cell data
cells = augment_with_gex(
    cells_df,
    "filtered_feature_bc_matrix.h5",
    gene_list=["CD3D", "CD4", "CD8A", "GZMA"],  # optional custom genes
    gene_groups={"CD3": ["CD3D", "CD3E", "CD3G"]},  # optional groups
)

# Aggregate expression by clonotype
clonotypes = aggregate_gex_by_clonotype(cells, operations=["sum", "mean"])

# Compute CD4-only and CD8-only cell counts
cd4_cd8 = compute_cd4_cd8_counts(cells)
```

### Multi-Experiment Unification

```python
from tcrsift import merge_experiments, add_phenotype_confidence, get_unify_summary

# Prepare experiment list
experiments = [
    (til_clonotypes, "TIL"),
    (culture_clonotypes, "Culture"),
    (amplify_clonotypes, "Amplify"),
]

# Merge with prefixed columns
unified = merge_experiments(experiments, add_occurrence_flags=True, add_combined_stats=True)

# Add phenotype confidence based on combined evidence
unified = add_phenotype_confidence(unified, ratio_threshold=10.0)

# Get summary
summary = get_unify_summary(unified)
# Returns: {"total_clonotypes": 1234, "from_TIL": 567, "Confident_CD8": 890, ...}
```

### Sequence Assembly

```python
from tcrsift import assemble_full_sequences, translate_dna, validate_sequences, export_fasta

# Assemble with default leaders (CD28 on alpha, CD8A on beta)
assembled = assemble_full_sequences(clonotypes, include_constant=True)

# Custom leaders
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader="CD8A",    # or None, "from_contig"
    beta_leader="CD28",     # or None, "from_contig"
    include_constant=True,
    linker="P2A",
)

# Extract native leaders from contigs
assembled = assemble_full_sequences(
    clonotypes,
    contigs_dir="/path/to/contigs",
    alpha_leader="from_contig",
    beta_leader="from_contig",
)

# No leaders
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader=None,
    beta_leader=None,
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
