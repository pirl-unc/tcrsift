# TCRsift

[![Tests](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml/badge.svg)](https://github.com/pirl-unc/tcrsift/actions/workflows/tests.yml)
[![Documentation](https://github.com/pirl-unc/tcrsift/actions/workflows/docs.yml/badge.svg)](https://pirl-unc.github.io/tcrsift/)
[![Coverage Status](https://coveralls.io/repos/github/pirl-unc/tcrsift/badge.svg)](https://coveralls.io/github/pirl-unc/tcrsift)
[![PyPI version](https://badge.fury.io/py/tcrsift.svg)](https://pypi.org/project/tcrsift/)

Prioritize T-cell receptor clonotypes from paired 10x Genomics single-cell VDJ
and gene-expression data for **human** (*Homo sapiens*) and **mouse** (*Mus
musculus*) studies. Standard Cell Ranger directories serve as input:
`cellranger vdj` supplies TCR contigs and clonotypes, and `cellranger count`
supplies gene expression for CD4/CD8 phenotyping, signature scoring, and
cell-type annotation (VDJ-only workflows are also supported).

TCRsift combines clonal abundance, T-cell lineage, gene-expression state
signatures, public TCR database matches (VDJdb, IEDB, CEDAR), sequence
generation probability (`Pgen`/`Ppost`), and optional TIL overlap to produce
auditable shortlists and synthesis-ready full-length constructs. None of these
signals predicts epitope specificity on its own; functional testing is still
required.

- **Documentation:** [pirl-unc.github.io/tcrsift](https://pirl-unc.github.io/tcrsift/)
- **Method overview:** [How TCRsift works](https://pirl-unc.github.io/tcrsift/how-it-works/)
- **Empirical evaluation (solid-tumor CD8 TIL, 5 cohorts / 33 patients):** [Benchmark on public TIL data](https://pirl-unc.github.io/tcrsift/benchmarks/)
- **Empirical evaluation (infectious disease & blood, 3 cohorts):** [Benchmark on infectious disease and blood](https://pirl-unc.github.io/tcrsift/benchmarks-infectious/)
- **Empirical evaluation (paired blood-tumor, post-ACT & CD4+ TIL, 4 cohorts):** [Benchmark on paired blood-tumor, post-ACT, and CD4+ TILs](https://pirl-unc.github.io/tcrsift/benchmarks-blood-cd4/)
- **Empirical evaluation (murine `--species mouse`, 3 cohorts / 19 libraries):** [Benchmark on murine models](https://pirl-unc.github.io/tcrsift/benchmarks-mouse/)

## Quick example

```bash
pip install tcrsift
```

Create `samples.yaml` pointing to your Cell Ranger output directories:

```yaml
samples:
  - sample: "Patient1_Culture"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
```

Run the antigen-stimulated culture pipeline:

```bash
tcrsift run --sample-sheet samples.yaml --output-dir results/
```

Or rank clones by expression state in fresh TIL, blood, pleural effusion, or
marrow:

```bash
tcrsift prioritize til_samples.yaml -o candidates/ --context solid-tumor
```

See [Cell Ranger Requirements](#cell-ranger-requirements) for expected input
files and the [Clone Prioritization guide](https://pirl-unc.github.io/tcrsift/user-guide/til-signatures/)
for context presets and cutoff rules.

## Contents

- [Quick example](#quick-example)
- [Choose a Workflow](#choose-a-workflow)
- [Architecture](#architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Sample Sheet Format](#sample-sheet-format)
- [Core Pipeline Steps](#core-pipeline-steps)
- [Supplementary Tools](#supplementary-tools)
- [Workflows](#workflows)
- [API Reference](#api-reference)
- [Output Files](#output-files)

## Choose a Workflow

| Study design | Command |
| --- | --- |
| Antigen-stimulated culture, optionally with matched TIL | `tcrsift run` |
| Fresh TIL, blood, MPE, or hematologic samples (paired VDJ + GEX) | `tcrsift prioritize --context ...` ([guide](https://pirl-unc.github.io/tcrsift/user-guide/til-signatures/)) |
| Longitudinal TIL samples in the legacy v2 per-timepoint layout | `tcrsift til-select` |
| Single-cell GEX QC, Harmony/UMAP/Leiden embedding, and cell typing | `tcrsift cells` ([guide](https://pirl-unc.github.io/tcrsift/user-guide/atlas/)) |

## Architecture

```
Cell Ranger VDJ + GEX  -->  tcrsift run  -->  clonotypes.csv + full_sequences.csv
                                 |
          load -> phenotype -> clonotype -> filter -> annotate -> assemble

Additional workflows: prioritize / til-prioritize, cells, data, pgen,
                      load-sct, annotate-gex, match-til, til-clonotype, til-select, unify
```

### Key Data Structures

| Stage | Object | Contents |
| --- | --- | --- |
| Load → Phenotype | `AnnData` | Per-cell expression counts in `.X` + paired VDJ annotations in `.obs` |
| Clonotype → Assemble | `DataFrame` | Per-clonotype counts, frequencies, scores, annotations, and sequences |

### Cell Ranger Requirements

**VDJ directory** (from `cellranger vdj`):

```
vdj_outs/
├── filtered_contig_annotations.csv   # Required (or all_contig_annotations.csv)
├── clonotypes.csv                    # Optional: Cell Ranger clonotype calls
├── consensus_annotations.csv         # Optional: consensus V(D)J sequences
└── filtered_contig.fasta             # Optional: donor leader & J->C contig reads
```

Columns read from `filtered_contig_annotations.csv`:

- `barcode`, `chain` (`TRA`/`TRB`), `cdr3` (amino-acid sequence)
- `v_gene`, `d_gene`, `j_gene`, `c_gene`
- `umis`, `reads`, `contig_id`
- `productive`: non-productive contigs are excluded before selecting each
  cell's highest-UMI chain

**GEX directory** (from `cellranger count`):

```
gex_outs/
├── filtered_feature_bc_matrix/       # Preferred when present
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
└── filtered_feature_bc_matrix.h5     # Used when matrix directory is absent
```

Phenotyping reads `CD3D`, `CD3E`, `CD3G`, `CD4`, `CD8A`, and `CD8B` in human
samples or `Cd3d`, `Cd3e`, `Cd3g`, `Cd4`, `Cd8a`, and `Cd8b1` in mouse samples
(matched by HGNC/MGI gene symbol or `ENSG*`/`ENSMUSG*` Ensembl ID).

**Barcode matching:** Barcodes are matched exactly between VDJ and GEX. Only
when a sample has zero exact matches (for example, `ACGT-1` vs `ACGT-2` library
suffixes) does TCRsift fall back to core barcode prefixes, skipping ambiguous
core barcodes and warning if nothing matches.

## Installation

```bash
pip install tcrsift
```

From source:

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
pip install -e .
```

### Optional Dependencies

```bash
pip install "tcrsift[full]"     # All optional extras below
pip install "tcrsift[tcrpeg]"   # Deep autoregressive TCRpeg/PyTorch sequence probability
pip install "tcrsift[atlas]"    # Harmony + Leiden clustering (`tcrsift cells embed`)
pip install "tcrsift[plots]"    # UpSet set-overlap plots
pip install "tcrsift[reports]"  # Legacy pdfkit HTML-to-PDF helper
```

Sequence assembly, Excel support, `reportlab` PDF sequence reports, `oncoref`
CTA panels, and bundled NumPy k-mer `Pgen`/`Ppost` probability models are
included in the base package.

## Quick Start

### Command Line

```bash
# Optional one-time cache download for VDJdb, IEDB, CEDAR, and IMGT V- and J-region references:
tcrsift data download

tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/
```

### With a YAML Configuration File

```bash
tcrsift generate-config -o my_config.yaml
tcrsift run --config my_config.yaml --sample-sheet samples.yaml -o results/
```

### Python API

```python
import tcrsift

adata = tcrsift.load_samples("samples.yaml")
adata = tcrsift.phenotype_cells(adata)
clonotypes = tcrsift.aggregate_clonotypes(adata)
filtered = tcrsift.filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")
annotated = tcrsift.annotate_clonotypes(filtered, vdjdb_path="/path/to/vdjdb")
assembled = tcrsift.assemble_full_sequences(annotated, include_constant=True)
```

## Sample Sheet Format

Sample sheets can be written in YAML or CSV.

### YAML

```yaml
samples:
  - sample: "Patient1_Culture"
    patient_id: "Patient1"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"
    mhc_allele: "HLA-A*02:01"
    source: "culture"

  - sample: "Patient1_TIL"
    patient_id: "Patient1"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
    tissue: "tumor"
```

### CSV

```csv
sample,patient_id,vdj_dir,gex_dir,antigen_type,source
Patient1_Culture,Patient1,/data/patient1/vdj,/data/patient1/gex,short_peptide,culture
Patient1_TIL,Patient1,/data/patient1_til/vdj,,,til
```

### Key Fields

| Field | Required | Description |
| --- | --- | --- |
| `sample` | Yes | Unique sample identifier |
| `vdj_dir` | Yes* | Cell Ranger VDJ output directory |
| `gex_dir` | No* | Cell Ranger GEX output directory (or `gex_path` for `.h5` / `.h5ad`) |
| `source` | No | `culture` (default), `control`, `til`, `tetramer`, or `sct` |
| `patient_id` | No | Donor or patient ID; enables per-donor tiering and sharing tables |
| `enrichment_method` | No | Sorting or stimulation label (e.g. `AIMpos`, `tetpos`); enables per-method outputs |
| `timepoint` | No | Longitudinal timepoint label (e.g. `D7`, `D14`, `pre`, `post`) |
| `apc_type` | No | Antigen-presenting cell label (e.g. `mDC`, `B-LCL`) |

\*At least one data path (`vdj_dir`, `gex_dir`/`gex_path`, `sct_path`, `til_csv`, or
`til_h5ad`) is required; `prioritize` requires `vdj_dir` plus `gex_dir` (or `gex_path`).

### Multi-donor and multi-method studies

Populating `patient_id` and `enrichment_method` enables:

- Per-`(donor, method)` clonotype summaries (`n_methods_per_donor`, `max_methods_per_donor`)
- Ranked CSVs per donor and method under `data/filtered_by_method/`
- Pairwise and UpSet method-overlap tables and plots, plus tier-recovery charts
- Per-donor tiering by default (`--fdr-scope auto` resolves to `per-donor`)

When multiple donors share the target antigen and MHC restriction and should be
ranked against a single pooled distribution, set `donors_share_antigen: true`
at the root of the YAML sheet (or pass `--fdr-scope global`).

### Antigen Types

| `antigen_type` | Default T-cell lineage | Description |
| --- | --- | --- |
| `short_peptide` | `CD8` | Minimal 8–11 aa MHC-I peptides |
| `long_peptide` | `mixed` | 15–25+ aa peptides requiring processing |
| `peptide_pool` | `mixed` | Multi-peptide stimulation pool |
| `minigene`, `minigene_library` | `mixed` | Single or pooled minigene constructs |
| `whole_protein`, `mrna` | `mixed` | Full-length protein or mRNA |
| `tetramer_mhc1`, `sct` | `CD8` | MHC-I tetramer or single-chain trimer sort |
| `tetramer_mhc2` | `CD4` | MHC-II tetramer sort |

## Core Pipeline Steps

### 1. Load Data

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

Reads productive VDJ contigs and GEX matrices, pairs cells by barcode, extracts
`CD3D/E/G`, `CD4`, and `CD8A/B` counts, and stores primary α/β chains per cell.

### 2. Phenotype Cells

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --cd4-cd8-ratio 3.0
```

- **Confident CD8+:** `(CD8A + CD8B + 1) / (CD4 + 1) > ratio` (default `3.0`)
- **Confident CD4+:** `(CD4 + 1) / (CD8A + CD8B + 1) > ratio`
- **Likely CD8+ / CD4+:** non-zero expression of one lineage and zero of the other
- **Unknown:** comparable expression of both or zero of both

### 3. Aggregate Clonotypes

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv --group-by CDR3ab
```

Groups cells by paired `CDR3ab` (or `CDR3b_only`) and computes `cell_count`,
`frequency`, `max_frequency`, `Tcell_type_consensus`, and `samples` (`tcrsift
run` also attaches bundled k-mer `Pgen`/`Ppost` log-probabilities automatically;
run `tcrsift pgen annotate` when executing steps individually).

### 4. Filter Clonotypes

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

Default abundance tiers:

| Tier | Min cells | Min frequency |
| --- | --- | --- |
| `tier1` | 10 | 1% |
| `tier2` | 5 | 0.5% |
| `tier3` | 3 | 0.1% |
| `tier4` | 2 | 0.05% |
| `tier5` | 2 | 0% |

#### Tiering scope (`--fdr-scope`)

| Value | Behavior | Tier files written |
| --- | --- | --- |
| `auto` (default) | `per-donor` when `n_donors > 1` and `donors_share_antigen` is false; otherwise `global` | Depends on resolved scope |
| `per-donor` | Ranks and tiers each donor independently | `filtered_tier{1..5}_<donor>.csv` |
| `per-sample` | Ranks and tiers each sample independently | `filtered_tier{1..5}_<sample>.csv` |
| `global` | Pools all samples into one ranking | `filtered_tier{1..5}.csv` |

#### Method overlap and recovery (`enrichment_method`)

When two or more samples are loaded (or `--emit-clone-sample-long always` is
passed), `tcrsift run` writes `data/clone_sample_long.csv` with one row per
`(clone, sample)` pair and per-sample abundance tiers. When `enrichment_method`
is also populated on the sample sheet, `tcrsift run` writes:

```
data/filtered_by_method/<donor>__<m>.csv    # Top-N clones per (donor, method)
data/method_overlap_<donor>.csv             # Pairwise method overlap (--method-overlap-similarity)
plots/method_overlap_<donor>.png            # Pairwise overlap heatmap
data/method_overlap_sets.csv                # Multi-way intersection counts
plots/method_overlap_upset.png              # UpSet / intersection plot
data/method_recovery.csv                    # Per-(donor, method) recovery of strictest populated tier
plots/method_recovery.png                   # Paired-by-donor recovery bar chart
```

### 5. Annotate Clonotypes and Manage Reference Databases

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb
```

Matches clonotypes against **VDJdb**, **IEDB**, and **CEDAR**, setting
`db_match`, separating official encoding gene symbols (`db_antigen_gene`, e.g.
`CTAG1A`/`CTAG1B` vs `CTAG2`, `MLANA`/`Mlana`, `DCT`/`Dct`, `CEACAM5`; `None`
for viral or foreign proteins such as chicken OVA) from canonical protein labels
(`db_protein_canonical`, e.g. `NY-ESO-1`, `LAGE-1`, `MART-1`, `TRP-2`, `CEA`),
and flagging viral epitopes with `is_viral=True`.

When explicit paths are omitted, `run`, `annotate`, `prioritize`,
`til-prioritize`, and `til-select` automatically check the managed cache
(`--cache-dir` → `$TCRSIFT_DATA_DIR` → `$XDG_CACHE_HOME/tcrsift` →
`~/.cache/tcrsift`), and `--auto-download` fetches missing databases on demand:

```bash
tcrsift data list                             # Show cached databases, sizes, and paths
tcrsift data download                         # Fetch VDJdb, IEDB, CEDAR, and human/mouse IMGT germlines
tcrsift data download --group germline_mouse  # Fetch mouse IMGT V-REGION, J-REGION, and leader references
tcrsift data download --db iedb               # Fetch a single database (--force refreshes)
tcrsift data clear                            # Clear cached files (--db NAME for one)
```

### 6. Assemble Full-Length Sequences

```bash
tcrsift assemble -i annotated.csv -o full_sequences.csv \
    --contigs-dir /path/to/vdj \
    --alpha-leader CD28 --beta-leader CD8A \
    --include-constant --single-chain --linker T2A
```

Builds `[Leader] + [V(D)J] + [Constant]` chains and `[Beta] + [2A linker] +
[Alpha]` single-chain constructs using `CD8A`, `CD28`, `mCD8B`, `mCD28`,
`H2Kb`, `IgK`, `TRAC`, `TRBC`, or native `from_contig` leaders and packaged
canonical human (`TRAC`/`TRBC1`/`TRBC2`) or mouse (`Trac`/`Trbc1`/`Trbc2`)
constant regions (`--species {auto,human,mouse}`). For chimeric human-VDJ /
murine-constant TCR-T constructs, pass `--species human --constant-species
mouse`. Pass `--allow-canonical-fallback` when contig FASTAs are not available.

## Supplementary Tools

### Expression-State Clone Prioritization (`prioritize` / `til-prioritize`)

Select a run-wide shortlist (default: CD8, `--min-cd3 3`, up to 200 clones) by
round-robin across context-specific expression signatures calibrated against
size-matched random cell groups:

```bash
tcrsift prioritize samples.yaml -o candidates/ --context solid-tumor
```

Contexts: `generic`, `blood`, `blood-tumor`, `solid-tumor` (default for
`til-prioritize`), `mpe`, and `heme` (supports `--species {auto,human,mouse}`).
Known MART-1 (`MLANA` / `Mlana` protein product) matches are excluded by default
and viral matches are retained unless `--exclude-known-viral` is set; reference
databases are resolved from `--vdjdb`/`--iedb`/`--cedar` or the `tcrsift data`
cache (`--auto-download` fetches missing tables automatically).
Writes `candidate_clones.csv`, `all_scored_clones.csv`,
`clone_sample_scores.csv`, `signature_background.csv`, and
`prioritization.json`. See the
[Clone Prioritization guide](https://pirl-unc.github.io/tcrsift/user-guide/til-signatures/).

### Single-Cell Atlas Commands (`cells`)

```bash
tcrsift cells qc       -i raw.h5ad   -o qc.h5ad    --min-genes 250 --max-mito 8
tcrsift cells embed    -i qc.h5ad    -o embed.h5ad --batch-key sample --resolution 1.0
tcrsift cells annotate -i embed.h5ad -o atlas.h5ad --solid-tumor
```

See the [Single-Cell Atlas guide](https://pirl-unc.github.io/tcrsift/user-guide/atlas/).

### Automatic and Cross-Run TIL Matching

When a sample sheet contains both `source: culture` and `source: til` entries,
`tcrsift run` automatically excludes TIL cells from culture aggregation and
adds `til_match`, `til_cell_count`, and `til_frequency` to culture clonotypes.

When culture and TIL were processed in separate runs, use `match-til`:

```bash
tcrsift match-til -i culture_clonotypes.csv --til-h5ad til_processed.h5ad -o matched.csv

tcrsift match-til -i culture_clonotypes.csv -o matched.csv \
    --til-sample T1=csv:/path/to/til_t1.csv \
    --til-sample T2=h5ad:/path/to/til_t2.h5ad \
    --til-sample T3=vdj:/path/to/til_t3_vdj_outs
```

### TIL-Only Aggregation (`til-clonotype`) and Legacy Timepoints (`til-select`)

Harmonize TIL-only counts and frequencies across samples:

```bash
tcrsift til-clonotype -o til_clonotypes.csv \
    --til-sample T1=csv:/path/to/til_t1.csv \
    --til-sample T2=h5ad:/path/to/til_t2.h5ad
```

For longitudinal TIL cohorts stored in the legacy v2 per-timepoint layout
(`consensus_annotations.<TP>.csv`, `clonotypes.<TP>.csv`,
`filtered_contig_annotations.<TP>.csv`,
`sample_filtered_feature_bc_matrix.<TP>.h5`):

```bash
tcrsift til-select \
  --data-dir /path/to/til_timepoints \
  --vdjdb /path/to/vdjdb.txt \
  --iedb /path/to/iedb_tcr_full_v3.tsv \
  --cedar /path/to/cedar_tcr_full_v3.tsv \
  --rank-by marker_score_z_mean
```

### Load SCT Platform Workbooks (`load-sct`)

```bash
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate
```

Applies `high_quality` (SNR ≥ 2.0, reads ≥ 10/chain, mutation match) and
`chosen` (SNR ≥ 3.4, reads ≥ 50, comPACT match) flags.

### Augment VDJ Tables with Standalone GEX (`annotate-gex`)

Attach per-cell or clonotype-aggregated expression from a standalone 10x `.h5`
matrix when `gex_dir` was omitted during `load`:

```bash
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file filtered_feature_bc_matrix.h5 \
    --aggregate \
    --cd4-cd8-counts \
    -o clonotype_gex.csv
```

### Merge Independent Experiments (`unify`) and Pronounceable Names (`mnemonic`)

```bash
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv sct_clonotypes.csv \
    -o unified.csv

tcrsift mnemonic -i clonotypes.csv -o clonotypes_named.csv
```

| Scenario | Command |
| --- | --- |
| Culture + TIL on the same sample sheet | `tcrsift run` (auto-matched) |
| Culture + TIL processed in separate runs | `tcrsift match-til` |
| Fresh TIL, blood, MPE, or heme VDJ + GEX prioritization | `tcrsift prioritize` / `til-prioritize` |
| Longitudinal TIL in legacy v2 per-timepoint files | `tcrsift til-select` |
| Combining clonotype tables across independent runs or platforms | `tcrsift unify` |

## Workflows

### Culture + Matched TIL

```yaml
samples:
  - sample: "Culture_Pool1"
    vdj_dir: "/data/culture/vdj"
    gex_dir: "/data/culture/gex"
    source: "culture"
  - sample: "TIL"
    vdj_dir: "/data/til/vdj"
    source: "til"
```

```bash
tcrsift run --sample-sheet samples.yaml -o results/
```

### Multi-Platform Unification

```bash
tcrsift til-clonotype -o til_results/clonotypes.csv --sample-sheet til_samples.yaml
tcrsift run --sample-sheet culture_samples.yaml -o culture_results/
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate

tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv sct_clonotypes.csv \
    -o unified_clonotypes.csv
```

## API Reference

```python
from tcrsift import (
    add_phenotype_confidence,
    aggregate_clonotypes,
    aggregate_sct,
    annotate_clonotypes,
    assemble_full_sequences,
    export_fasta,
    filter_by_tcell_type,
    filter_clonotypes,
    get_clonotype_summary,
    get_phenotype_summary,
    get_sct_specificities,
    get_til_summary,
    load_cellranger_gex,
    load_cellranger_vdj,
    load_samples,
    load_sct,
    load_vdjdb,
    match_til,
    merge_experiments,
    phenotype_cells,
    split_by_tier,
)

# Load and phenotype
adata = load_samples("samples.yaml")
adata = phenotype_cells(adata, cd4_cd8_ratio=3.0)
cd8_cells = filter_by_tcell_type(adata, tcell_type="cd8")

# Clonotype, tier, annotate, match TIL, and assemble
clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab", min_umi=2)
filtered = filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")
tier_dfs = split_by_tier(filtered)
annotated = annotate_clonotypes(filtered, vdjdb_path="/path/to/vdjdb", exclude_viral=True)
assembled = assemble_full_sequences(
    annotated,
    alpha_leader="CD28",
    beta_leader="CD8A",
    include_constant=True,
    linker="T2A",
)
export_fasta(assembled, "sequences.fasta", sequence_col="single_chain_aa")
```

## Output Files

### `data/clonotypes.csv`

| Column | Description |
| --- | --- |
| `CDR3ab` | Paired clonotype identifier (`CDR3_alpha` + `_` + `CDR3_beta`) |
| `CDR3_alpha`, `CDR3_beta` | Chain CDR3 amino-acid sequences |
| `cell_count`, `frequency`, `max_frequency` | Total cells and within-sample frequencies |
| `Tcell_type_consensus` | Majority cell phenotype |
| `ppost_alpha`, `ppost_beta`, `pgen_alpha`, `pgen_beta` | Natural-log repertoire and generation probabilities |
| `tier` | Strictest confidence tier (`tier1`–`tier5`) |
| `db_match`, `is_viral` | Public-database match and viral-epitope flag |
| `n_donors`, `n_methods_per_donor`, `max_methods_per_donor` | Populated when `patient_id` and `enrichment_method` are set |

### `data/full_sequences.csv`

| Column | Description |
| --- | --- |
| `CDR3ab` | Clonotype identifier |
| `full_alpha_aa`, `full_beta_aa` | Full-length chain proteins (leader + V(D)J + constant) |
| `single_chain_aa`, `single_chain_nt` | Assembled β-2A-α protein and DNA construct |

## Contributing and License

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines. Licensed
under the Apache License 2.0 ([LICENSE](LICENSE)).

## Citation

```bibtex
@software{tcrsift,
  author = {Rubinsteyn, Alex},
  title = {TCRsift: T-cell receptor clonotype prioritization from single-cell VDJ and gene-expression data},
  url = {https://github.com/pirl-unc/tcrsift},
  year = {2024}
}
```
