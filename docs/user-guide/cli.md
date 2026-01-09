# CLI Reference

TCRsift provides a command-line interface with subcommands for each pipeline step.

## Global Options

```bash
tcrsift --version      # Show version
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
    --method threshold \
    --report
```

| Option | Description |
|--------|-------------|
| `--sample-sheet` | Path to sample sheet (required) |
| `--output-dir`, `-o` | Output directory (required) |
| `--vdjdb` | Path to VDJdb |
| `--iedb` | Path to IEDB |
| `--cedar` | Path to CEDAR |
| `--tcell-type` | Filter type: `cd8`, `cd4`, `both` |
| `--method` | Filter method: `threshold`, `logistic` |
| `--report` | Generate summary report |

---

### `tcrsift load`

Load CellRanger data.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

| Option | Description |
|--------|-------------|
| `--sample-sheet` | Path to sample sheet (required) |
| `-o`, `--output` | Output .h5ad file (required) |

---

### `tcrsift phenotype`

Classify cells as CD4+ or CD8+.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --ratio 3.0 --min-cd3 10
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input .h5ad file (required) |
| `-o`, `--output` | Output .h5ad file (required) |
| `--ratio` | CD4/CD8 ratio threshold (default: 3.0) |
| `--min-cd3` | Minimum CD3 reads (default: 10) |
| `--plot` | Save phenotype plot |

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
| `--plot` | Save clonotype plot |

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
| `--plot` | Save filter plot |

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

Match against TIL data.

```bash
tcrsift match-til -i annotated.csv --til-data til.csv -o matched.csv
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Culture clonotypes (required) |
| `--til-data` | TIL clonotypes (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--match-by` | `CDR3ab` or `CDR3b_only` |

---

### `tcrsift assemble`

Assemble full-length sequences.

```bash
tcrsift assemble -i annotated.csv -o sequences.csv \
    --contigs-dir /path/to/contigs \
    --include-leader \
    --include-constant \
    --linker T2A \
    --export-fasta sequences.fasta
```

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input CSV file (required) |
| `-o`, `--output` | Output CSV file (required) |
| `--contigs-dir` | Directory with contig FASTAs |
| `--include-leader` | Extract leader peptides |
| `--include-constant` | Add constant regions |
| `--constant-source` | `ensembl` or `from-data` |
| `--linker` | Linker for single-chain (default: T2A) |
| `--export-fasta` | Export to FASTA file |

---

### `tcrsift mnemonic`

Generate mnemonic names for TCRs.

```bash
tcrsift mnemonic --cdr3b CASSLGQAYEQYF --cdr3a CAVSDGGSQGNLIF
```

| Option | Description |
|--------|-------------|
| `--cdr3b` | Beta chain CDR3 sequence |
| `--cdr3a` | Alpha chain CDR3 sequence |
