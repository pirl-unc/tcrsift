# CLI Reference

TCRsift provides a command-line interface with subcommands for each pipeline step.

## Global Options

```bash
tcrsift --version      # Show version
tcrsift -v             # Show version (short form)
tcrsift --help         # Show help
```

## Commands

### `tcrsift run`

Run the complete pipeline.

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb \
    --cedar /path/to/cedar \
    --tcell-type cd8 \
    --method threshold
```

| Option | Description |
|--------|-------------|
| `--sample-sheet`, `-s` | Path to sample sheet (required) |
| `--output-dir`, `-o` | Output directory (required) |
| `--config`, `-c` | YAML configuration file |
| `--vdjdb` | Path to VDJdb |
| `--iedb` | Path to IEDB |
| `--cedar` | Path to CEDAR |
| `--til-samples` | Comma-separated TIL sample names from the sample sheet |
| `--til-sample` | Repeatable direct TIL sample spec: `NAME=TYPE:PATH` (TYPE: `csv`, `h5ad`, `vdj`) |
| `--til-match-by` | TIL matching mode: `CDR3ab` or `CDR3b_only` |
| `--tcell-type` | Filter type: `cd8`, `cd4`, `both` |
| `--method` | Filter method: `threshold`, `logistic` |
| `--no-report` | Skip generating summary report |
| `--skip-plots` | Skip generating plots |

Use either `--til-samples` or `--til-sample`, not both.

---

### `tcrsift load`

Load CellRanger data.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

| Option | Description |
|--------|-------------|
| `--sample-sheet`, `-s` | Path to sample sheet (required) |
| `-o`, `--output` | Output .h5ad file (required) |

---

### `tcrsift cells`

The single-cell gene-expression path — QC → embed → annotate — as three
`.h5ad`-in/`.h5ad`-out subcommands wrapping `cell_qc_funnel`, `embed_cells`, and
`annotate_cells`. See the [Single-Cell Atlas Path](atlas.md) guide for the
library equivalents and the solid-tumor overrides.

```bash
tcrsift cells qc       -i raw.h5ad   -o qc.h5ad    --min-genes 250 --max-mito 8
tcrsift cells embed    -i qc.h5ad    -o embed.h5ad --batch-key sample --resolution 1.0
tcrsift cells annotate -i embed.h5ad -o atlas.h5ad --solid-tumor
```

**`cells qc`** — per-cell QC + doublet funnel.

| Option | Description |
|--------|-------------|
| `-i`, `--input` / `-o`, `--output` | Input (raw counts) / output .h5ad (required) |
| `--min-genes` / `--max-genes` / `--min-counts` / `--max-counts` | Coverage bounds (`N` or `none` to disable; default: the funnel's) |
| `--min-mito` / `--max-mito` | Mito-percent floor / ceiling (`N` or `none`) |
| `--pbmc-preset` | Use the opinionated PBMC `LoadConfig()` thresholds as the baseline |
| `--solid-tumor` | Use `SOLID_TUMOR_LINEAGE_PROGRAMS` for the cross-lineage doublet gate |
| `--waterfall-csv` | Also write the QC waterfall to a CSV |

**`cells embed`** — Pearson-residual PCA → Harmony → UMAP → Leiden.

| Option | Description |
|--------|-------------|
| `-i`, `--input` / `-o`, `--output` | Input (raw counts) / output .h5ad (required) |
| `--batch-key` | obs column to integrate over with Harmony (default: none) |
| `--resolution` / `--n-pcs` / `--n-neighbors` / `--n-top-genes` / `--seed` | Clustering knobs |

**`cells annotate`** — per-cluster cell-type + T/B-state labels.

| Option | Description |
|--------|-------------|
| `-i`, `--input` (has `obs['leiden']`) / `-o`, `--output` | Input / output .h5ad (required) |
| `--leiden-col` | Cluster column in obs (default: `leiden`) |
| `--solid-tumor` | Use `SOLID_TUMOR_CELL_TYPE_SIGNATURES` instead of the PBMC default |
| `--pbmc-culture` | Restrict the argmax to `PBMC_CULTURE_TYPES` |
| `--min-subtype-frac` / `--no-suffix` | Distinguishing-gene-suffix controls |

Requires the `atlas` extra (`pip install 'tcrsift[atlas]'`) for `embed`.

---

### `tcrsift phenotype`

Classify cells as CD4+ or CD8+.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --cd4-cd8-ratio 3.0 --min-cd3-reads 10
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input .h5ad file (required) |
| `-o`, `--output` | Output .h5ad file (required) |
| `--cd4-cd8-ratio` | CD4/CD8 ratio threshold (default: 3.0) |
| `--min-cd3-reads` | Minimum CD3 reads (default: 10) |
| `--plot-phenotype` | Save phenotype plot |

---

### `tcrsift clonotype`

Aggregate cells into clonotypes.

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv \
    --group-by CDR3ab \
    --min-umi 2 \
    --handle-doublets flag
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input .h5ad file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--group-by` | Grouping: `CDR3ab` or `CDR3b_only` |
| `--min-umi` | Minimum UMI count (default: 2) |
| `--handle-doublets` | `flag`, `remove`, `keep-primary` |
| `--airr` | Export AIRR format to this path |
| `--plot-clonotypes` | Save clonotype plot |

---

### `tcrsift filter`

Filter and tier clonotypes.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ \
    --method threshold \
    --tcell-type cd8 \
    --min-cells 2 \
    --exclude-viral
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output directory (required) |
| `--method` | `threshold` or `logistic` |
| `--tcell-type` | `cd8`, `cd4`, or `both` |
| `--min-cells` | Minimum cell count (default: 2) |
| `--min-frequency` | Minimum frequency |
| `--exclude-viral` | Exclude viral clones |
| `--fdr-tiers` | FDR values for logistic method |
| `--plot-filter` | Save filter plot |

---

### `tcrsift annotate`

Annotate with public databases.

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb \
    --cedar /path/to/cedar \
    --match-by CDR3ab \
    --exclude-viral
```

If no database paths are provided, `annotate` still succeeds and writes default
annotation columns (`db_match=False`, `is_viral=False`, etc.).

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--vdjdb` | Path to VDJdb |
| `--iedb` | Path to IEDB |
| `--cedar` | Path to CEDAR |
| `--match-by` | `CDR3ab` or `CDR3b_only` |
| `--exclude-viral` | Remove viral matches |
| `--flag-only` | Flag but don't remove viral |

---

### `tcrsift match-til`

Match against TIL data. Supports multiple input formats and multi-sample matching.

**Basic usage (single TIL source):**

```bash
# From h5ad file
tcrsift match-til -i clonotypes.csv --til-h5ad til.h5ad -o matched.csv

# From CSV file
tcrsift match-til -i clonotypes.csv --til-csv til.csv -o matched.csv

# From CellRanger VDJ directory
tcrsift match-til -i clonotypes.csv --til-vdj-dir /path/to/vdj -o matched.csv
```

**Multi-sample TIL matching (via sample sheet):**

```bash
tcrsift match-til -i clonotypes.csv -s til_samples.yaml -o matched.csv
```

**Multi-sample TIL matching (direct specs, no sample sheet):**

```bash
tcrsift match-til -i clonotypes.csv -o matched.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad \
  --til-sample T3=vdj:/path/to/til_t3_vdj_outs
```

Sample sheet format for multiple TIL samples:
```yaml
samples:
  - sample: "TIL_Sample1"
    source: til
    vdj_dir: "/path/to/vdj"
  - sample: "TIL_Sample2"
    source: til
    til_csv: "/path/to/til.csv"
```

**Required inputs:**

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Culture clonotypes CSV (required) |
| `-o`, `--output` | Output CSV file (required) |

**TIL data source (provide one):**

| Option | Description |
|--------|-------------|
| `-s`, `--sample-sheet` | Sample sheet with TIL samples (YAML or CSV) |
| `--til-h5ad` | Single TIL h5ad file |
| `--til-csv` | Single TIL CSV file (must have CDR3_alpha/CDR3_beta) |
| `--til-vdj-dir` | Single CellRanger VDJ output directory |
| `--til-sample` | Repeatable direct sample spec: `NAME=TYPE:PATH` (TYPE: `csv`, `h5ad`, `vdj`) |

**Matching options:**

| Option | Description |
|--------|-------------|
| `--match-by` | `CDR3ab` (default) or `CDR3b_only` |
| `--min-til-cells` | Min TIL cells to count as match (default: 1) |

**Output columns:**

For all inputs:
- `til_match`: Boolean indicating if clone was found in TIL
- `til_samples`: Comma-separated list of matching TIL samples
- `til_cell_count`: Total cells across all TIL samples
- `til_frequency`: Combined frequency

For multi-sample inputs, per-sample columns are added:
- `til_cell_count.{sample}`: Cells in specific TIL sample
- `til_frequency.{sample}`: Frequency in specific TIL sample

---

### `tcrsift til-clonotype`

Aggregate TIL-only data into clonotype counts/frequencies across one or more TIL samples.

```bash
# Single TIL source
tcrsift til-clonotype -o til_clonotypes.csv --til-csv til.csv

# Multi-sample direct specs (no sample sheet)
tcrsift til-clonotype -o til_clonotypes.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad \
  --til-sample T3=vdj:/path/to/til_t3_vdj_outs
```

| Option | Description |
|--------|-------------|
| `-o`, `--output` | Output CSV file (required) |
| `-s`, `--sample-sheet` | TIL sample sheet (YAML/CSV) |
| `--til-h5ad` | Single TIL h5ad file |
| `--til-csv` | Single TIL CSV file |
| `--til-vdj-dir` | Single TIL VDJ directory |
| `--til-sample` | Repeatable direct sample spec: `NAME=TYPE:PATH` |
| `--match-by` | `CDR3ab` (default) or `CDR3b_only` |
| `--min-cells` | Minimum total TIL cells per clonotype (default: 1) |

---

### `tcrsift til-select`

Select promising TIL clonotypes from two or more ordered 10x VDJ+GEX tumor
samples/timepoints.

This command:

- harmonize clonotypes across timepoints,
- score per-clone marker expression from GEX,
- apply a CD8 + minimum-cell + known-virus base mask,
- add cytolytic, antigen-response, AIM, exhaustion, and frequency-change
  branches,
- export v2-style subset tables and reports.

The final union is a review shortlist, despite the legacy
`is_candidate_tumor_reactive` column name. It does not prove tumor-antigen
specificity. Frequency-increase branches are meaningful only when the input
order is a genuine longitudinal series.

```bash
# Auto-discover inputs from data directory:
tcrsift til-select --data-dir data --rank-by marker_score_z_mean

# Explicit mappings (v2-compatible):
tcrsift til-select \
  --samples \
    T1=consensus_annotations.T1.csv,clonotypes.T1.csv \
    T2=consensus_annotations.T2.csv,clonotypes.T2.csv \
    T3=consensus_annotations.T3.csv,clonotypes.T3.csv
```

Per-timepoint required files for marker scoring:
- `filtered_contig_annotations.<TP>.csv`
- `sample_filtered_feature_bc_matrix.<TP>.h5`

| Option | Description |
|--------|-------------|
| `--data-dir` | Directory with per-timepoint files (default: `./data`) |
| `--samples`, `--inputs` | Explicit mappings `LABEL=CONSENSUS,CLONOTYPES` |
| `--config` | YAML mapping with `timepoints:` block |
| `--vdjdb`, `--iedb`, `--cedar` | Optional public DB annotation sources |
| `--match-by` | DB matching mode: `CDR3ab` or `CDR3b_only` (default: `CDR3b_only`) |
| `--marker-genes` | Marker panel for GEX scoring |
| `--cytolytic-genes`, `--antigen-response-genes` | Compact score branches |
| `--activation-genes`, `--exhaustion-genes` | AIM and chronic-antigen score branches |
| `--min-cells-per-clone` | Base selection minimum total cells |
| `--min-cd8-cp10k` | Base selection CD8 threshold |
| `--max-cd4-to-cd8-ratio` | Base selection CD4/CD8 ratio ceiling |
| `--increase-ratio-nonzero-min` | Increasing branch ratio threshold |
| `--immunogenic-percentile` | Percentile cutoff for immunogenic branch |
| `--rank-by` | Top-k ranking metric |
| `--fig-dir` | Output directory for plots/subsets (default: `./figures`) |
| `--out-table` | Harmonized CSV output (default: `./abTCR_harmonized.csv`) |
| `--out-selected-report` | Selected-clone PDF report path |

Main outputs include:
- `abTCR_harmonized.csv`
- `figures/abTCR_master_table.csv`
- `figures/abTCR_annotated.csv`
- `figures/selection_masks.csv`
- `figures/subset_*.csv`
- `figures/selection_funnel.png`
- `figures/selected_clones_report.pdf`

v2 CSV compatibility mode (default):
- `til-select` emits v2-compatible CSV column sets and ordering.
- With identical inputs/options, CSV outputs are expected to match legacy `v2/harmonize_abtcr_timepoints.py`.
- Plot/PDF artifacts are generated but may differ byte-for-byte across environments.

For standard CellRanger directories and published NeoTCR/MANAscore methods, see
[Multi-sample TIL Prioritization](til-signatures.md).

---

### `tcrsift assemble`

Assemble full-length sequences.

```bash
tcrsift assemble -i annotated.csv -o sequences.csv \
    --alpha-leader CD28 \
    --beta-leader CD8A \
    --include-constant \
    --linker T2A \
    --fasta sequences.fasta
```

**Required inputs:**

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output CSV file (required) |

**Conditionally required:**

| Option | Description |
|--------|-------------|
| `--contigs-dir` | Directory with CellRanger contig FASTAs. Required when using `--leaders-from-contigs` or `--alpha/beta-leader=from_contig` |

**Leader peptide options:**

| Option | Description |
|--------|-------------|
| `--alpha-leader` | Alpha chain leader: `CD8A`, `CD28`, `IgK`, `TRAC`, `TRBC`, `from_contig`, `none` (default: CD28) |
| `--beta-leader` | Beta chain leader (default: CD8A) |
| `--no-leaders` | Disable leaders on both chains |
| `--leaders-from-contigs` | Extract native leaders from CellRanger FASTAs (requires `--contigs-dir`) |

**Sequence options:**

| Option | Description |
|--------|-------------|
| `--include-constant` | Add constant regions |
| `--constant-source` | `canonical` (or legacy alias `ensembl`) for packaged references; `from-data` uses `*_constant_aa/nt` columns |
| `--single-chain` | Generate single-chain constructs (beta-linker-alpha) |
| `--linker` | Linker for single-chain (default: T2A) |

**Output options:**

| Option | Description |
|--------|-------------|
| `--fasta` | Export to FASTA file |
| `--airr` | Export to AIRR format |

---

### `tcrsift annotate-gex`

Annotate TCR data with gene expression.

```bash
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file filtered_feature_bc_matrix.h5 \
    -o cells_with_gex.csv \
    --aggregate \
    --cd4-cd8-counts
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--gex-file` | 10x filtered_feature_bc_matrix.h5 file (required) |
| `--barcode-col` | Column with cell barcodes (default: barcode) |
| `--genes` | Comma-separated list of genes |
| `--aggregate` | Aggregate expression by clonotype |
| `--cd4-cd8-counts` | Compute CD4/CD8 cell counts per clonotype |

---

### `tcrsift mnemonic`

Generate mnemonic names for TCRs.

```bash
tcrsift mnemonic -i clonotypes.csv -o clonotypes_with_names.csv
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--cdr3-col` | Column with CDR3 sequences (auto-detected) |
| `--name-col` | Output column name for mnemonics |

---

### `tcrsift load-sct`

Load TCR data from SCT platform Excel files.

```bash
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input Excel file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--aggregate` | Aggregate to clonotype level |
| `--min-snr` | Minimum SNR (default: 2.0) |
| `--min-reads` | Minimum reads per chain (default: 10) |

---

### `tcrsift unify`

Merge clonotype data from multiple experiments.

```bash
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv \
    -o unified.csv
```

| Option | Description |
|--------|-------------|
| `-i`, `--inputs` | Input CSV files (multiple, required) |
| `-o`, `--output` | Output CSV file (required) |

---

### `tcrsift generate-config`

Generate an example configuration file.

```bash
tcrsift generate-config -o my_config.yaml
```

| Option | Description |
|--------|-------------|
| `-o`, `--output` | Optional YAML path; without it, print to stdout |
| `--force` | Allow replacement of an existing output file |
