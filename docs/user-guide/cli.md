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
| `--tcell-type` | Filter type: `cd8`, `cd4`, `both` |
| `--method` | Filter method: `threshold`, `logistic` |
| `--no-report` | Skip generating summary report |
| `--skip-plots` | Skip generating plots |

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
| `--constant-source` | `ensembl` or `from-data` (use `*_constant_aa/nt` columns if present) |
| `--single-chain` | Generate single-chain constructs |
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
| `-i`, `--input` | Input CSV files (multiple, required) |
| `-o`, `--output` | Output CSV file (required) |

---

### `tcrsift generate-config`

Generate an example configuration file.

```bash
tcrsift generate-config -o my_config.yaml
```

| Option | Description |
|--------|-------------|
| `-o`, `--output` | Output YAML file (required) |
