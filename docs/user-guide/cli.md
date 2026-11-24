# CLI Reference

Command-line interface for TCRsift. Pass `--help` to any subcommand for its
complete flag list.

## Global options

```bash
tcrsift --version      # Print version (short form: -v)
tcrsift --help         # Tiered command overview (Primary workflows / Exploration / Modular steps)
```

`tcrsift --help` groups all 25 commands into three sections:

1. **Primary workflows:** `prioritize`, `run`, `assemble`, `data`
2. **Exploration & annotation:** `cells`, `annotate`, `report`
3. **Modular pipeline steps & legacy commands:** `load`, `phenotype`, `clonotype`, `filter`, `select`, `pgen`, `log-pgen`, `match-til`, `til-clonotype`, `til-prioritize`, `til-select`, `annotate-gex`, `load-sct`, `audit-alleles`, `cohort`, `unify`, `mnemonic`, `generate-config`

---

## Primary workflows

### `tcrsift run`

Run the end-to-end culture prioritization pipeline (`load → phenotype →
clonotype → filter → annotate → match-til → assemble`).

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
| --- | --- |
| `-s`, `--sample-sheet` | Sample sheet in YAML or CSV format (required) |
| `-o`, `--output-dir` | Output directory (required) |
| `-c`, `--config` | YAML configuration file |
| `--vdjdb`, `--iedb`, `--cedar` | Paths to public TCR databases (auto-discovered from `tcrsift data` cache if omitted) |
| `--cache-dir` | Override managed reference-database cache directory |
| `--auto-download` | Automatically download missing reference databases into the cache |
| `--species` | Organism (`auto`, `human`, or `mouse`; default `auto`) |
| `--til-samples` | Comma-separated TIL sample names from the sample sheet |
| `--til-sample` | Repeatable external TIL source: `NAME=TYPE:PATH` (`TYPE` in `csv`, `h5ad`, `vdj`) |
| `--til-match-by` | TIL matching key: `CDR3ab` (default) or `CDR3b_only` |
| `--tcell-type` | Lineage filter: `cd8`, `cd4`, or `both` |
| `--method` | Tiering method: `threshold` (default), `logistic`, or `stimulation` |
| `--fdr-scope` | Tiering scope: `auto` (default), `global`, `per-donor`, or `per-sample` |
| `--call-multimers` | Extract 10x Feature Barcode pMHC multimers and attach clonal concordance columns |
| `--germline-vdj-fallback` | Rescue 5'-truncated contigs and bare `(v_gene, CDR3, j_gene)` rows from cached IMGT V+J germlines |
| `--auto-download-germline` | Automatically fetch missing IMGT V-REGION and J-REGION references when `--germline-vdj-fallback` is set |
| `--no-pgen-ppost` | Skip default k-mer `Pgen`/`Ppost` scoring during clonotyping |
| `--no-report` | Skip PDF summary report generation |
| `--skip-plots` | Skip PNG plot generation |

Pass either `--til-samples` or `--til-sample`, not both.

---

### `tcrsift prioritize` / `tcrsift til-prioritize`

Select a run-wide shortlist of clonotypes from one or more paired Cell Ranger
VDJ + GEX samples by round-robin across gene-expression signature lists.
Defaults: `cd8`, `--min-cd3 3`, up to 200 clones total, known MART-1 exclusion
on, viral exclusion off (reference databases resolve from `--vdjdb`/`--iedb`/`--cedar`
or the `tcrsift data` cache; add `--auto-download` to fetch missing databases).
`til-prioritize` is an alias with `--context solid-tumor` as its default.

```bash
tcrsift prioritize samples.yaml -o candidates/ --context solid-tumor
```

| Option | Description |
| --- | --- |
| `--context` | Preset: `generic` (default for `prioritize`), `blood`, `blood-tumor`, `solid-tumor` (default for `til-prioritize`), `mpe`, `heme` |
| `--tcell-type` | `cd8` (default), `cd4`, or `both` |
| `--signatures NAME ...` | Explicit signature list replacing the context preset |
| `--exclude-signatures NAME ...` | Signatures to drop from the preset |
| `--max-clones N` | Run-wide shortlist cap (default `200`; `0` removes the cap while keeping cutoffs) |
| `--signature-cutoff` | `background` (default size-matched resampling) or `legacy` (`score > 0`) |
| `--signature-evidence` | `either` (default for `background`), `mean` (default for `legacy`), or `subset` |
| `--signature-quantile` | Within-stratum percentile rank floor (default `0.90`; `0` disables) |
| `--signature-background-quantile` | Random-group quantile for the mean floor (default `0.99`) |
| `--min-signature-score` | Optional raw-score minimum (default unset in `background` mode, `0` in `legacy`) |
| `--min-signature-support` | Minimum distinct signatures passing all cutoffs (each within a single stratum; default `1`) |
| `--min-cd3` | Minimum summed raw `CD3D/E/G` counts per cell (default `3`; `0` disables) |
| `--min-cells`, `--min-frequency` | Clone abundance floors (defaults `2`, `0.001`) |
| `--vdjdb`, `--iedb`, `--cedar`, `--iedb-epitope` | Reference databases (auto-discovered from `tcrsift data` cache if omitted) |
| `--cache-dir`, `--auto-download` | Cache directory override and on-demand database download |
| `--species` | `auto` (default), `human`, or `mouse` |
| `--round-robin-mode` | `tiered` (default when `--signatures` is unset; schedules dedicated tumor/neoantigen signatures ahead of general state panels) or `flat` |
| `--database-match` | Chain-matching policy for database annotations and exclusions: `strict_ab` (default for paired α/β `prioritize`), `ab_with_partial`, or `b_only` |
| `--exclude-known-viral` | Exclude clones matching viral entries in VDJdb/IEDB/CEDAR |
| `--no-exclude-known-mart1` | Retain known MART-1 database matches |
| `--exclude-treg`, `--treg-foxp3-fraction` | Opt-in exclusion of `FOXP3`/`Foxp3`+ Treg-like clones (audited in `risk_flags` as `treg_FOXP3`; default fraction threshold `0.5`) |
| `--call-multimers` | Call 10x Feature Barcode pMHC multimers and attach clonal specificity/concordance columns |
| `--assemble` | Assemble validated full-length `β-2A-α` TCR constructs for shortlisted candidates (`candidate_constructs.csv` + `candidate_constructs.fasta`; auto-detects sample `*contig*.fasta` leaders) |
| `--germline-vdj-fallback` | Reconstruct 5'-truncated or bare `(v_gene, CDR3, j_gene)` chains from cached IMGT V-REGION + J-REGION references when `--assemble` is enabled |
| `--single-chain-linker` | 2A linker for `--assemble` constructs: `T2A` (default), `P2A`, `E2A`, `F2A`, or `none` |

See [Clone Prioritization](til-signatures.md) for all QC, VDJ, expression, V/J
gene, CDR3 length, and subset-evidence options, and `tcrsift prioritize --help`
for the full flag list.

---

### `tcrsift cells`

Single-cell GEX atlas subcommands (`qc → embed → annotate`) operating on
`.h5ad` files. See [Single-Cell Atlas Workflow](atlas.md) for Python equivalents
and tumor-typing overrides.

```bash
tcrsift cells qc       -i raw.h5ad   -o qc.h5ad    --min-genes 250 --max-mito 8
tcrsift cells embed    -i qc.h5ad    -o embed.h5ad --batch-key sample --resolution 1.0
tcrsift cells annotate -i embed.h5ad -o atlas.h5ad --solid-tumor
```

**`tcrsift cells qc`**: coverage, mitochondrial, and doublet filtering:

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Raw-count input `.h5ad` and filtered output `.h5ad` (required) |
| `--min-genes`, `--max-genes`, `--min-counts`, `--max-counts` | Per-cell coverage bounds (`N` or `none` to disable) |
| `--min-mito`, `--max-mito` | Mitochondrial percentage floor and ceiling (`N` or `none`) |
| `--pbmc-preset` | Start from full `LoadConfig()` PBMC bounds |
| `--solid-tumor` | Use `SOLID_TUMOR_LINEAGE_PROGRAMS` for cross-lineage doublet detection |
| `--waterfall-csv` | Write step-by-step cell removal counts to CSV |

**`tcrsift cells embed`**: Pearson-residual PCA, Harmony, UMAP, and Leiden
(requires `pip install "tcrsift[atlas]"`):

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input (raw counts) and embedded output `.h5ad` (required) |
| `--batch-key` | `obs` column integrated by Harmony (default: none) |
| `--resolution`, `--n-pcs`, `--n-neighbors`, `--n-top-genes`, `--seed` | PCA, neighbor-graph, and Leiden parameters |

**`tcrsift cells annotate`**: cluster-level cell-type and T/B-state assignment:

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Clustered input (`obs["leiden"]`) and annotated output `.h5ad` (required) |
| `--leiden-col` | Cluster column in `obs` (default `leiden`) |
| `--solid-tumor` | Score `SOLID_TUMOR_CELL_TYPE_SIGNATURES` instead of the PBMC default |
| `--pbmc-culture` | Restrict candidate cell types to `PBMC_CULTURE_TYPES` |
| `--min-subtype-frac`, `--no-suffix` | Controls for distinguishing-gene state suffixes |

---

## Step-by-step pipeline commands

### `tcrsift load`

Load Cell Ranger VDJ and GEX directories from a sample sheet into one `.h5ad`
file.

```bash
tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
```

| Option | Description |
| --- | --- |
| `-s`, `--sample-sheet` | Sample sheet path (required) |
| `-o`, `--output` | Output `.h5ad` path (required) |

---

### `tcrsift phenotype`

Classify cells into CD4+ and CD8+ categories and apply the CD3 count floor.

```bash
tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad --cd4-cd8-ratio 3.0 --min-cd3-reads 3
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input and output `.h5ad` paths (required) |
| `--cd4-cd8-ratio` | Expression ratio threshold for confident calls (default `3.0`) |
| `--min-cd3-reads` | Minimum summed `CD3D/E/G` counts (default `3`; `0` disables) |
| `--plot-phenotype` | Generate phenotype summary plots (written to `--output-dir`) |

---

### `tcrsift clonotype`

Group cells into clonotypes by CDR3 sequence.

```bash
tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv \
    --group-by CDR3ab \
    --min-umi 2 \
    --handle-doublets flag
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input `.h5ad` and output `.csv` paths (required) |
| `--group-by` | `CDR3ab` (default) or `CDR3b_only` |
| `--min-umi` | Minimum UMI count per chain (default `2`) |
| `--handle-doublets` | `flag`, `remove`, or `keep-primary` |
| `--airr` | Optional path for AIRR-formatted TSV export |
| `--plot-clonotypes` | Generate clonotype abundance plots (written to `--output-dir`) |

---

### `tcrsift filter`

Assign clonotypes to confidence tiers (`tier1`–`tier5`).

```bash
tcrsift filter -i clonotypes.csv -o filtered/ \
    --method threshold \
    --tcell-type cd8 \
    --min-cells 2 \
    --exclude-viral
```

| Option | Description |
| --- | --- |
| `-i`, `--input` | Input clonotype CSV (required) |
| `-o`, `--output` | Output directory for tier CSVs (required) |
| `--method` | `threshold` (default), `logistic`, or `stimulation` |
| `--tcell-type` | `cd8`, `cd4`, or `both` |
| `--min-cells`, `--min-frequency` | Global abundance floors (default `--min-cells 2`) |
| `--exclude-viral` | Drop rows with `is_viral=True` |
| `--fdr-tiers` | Comma-separated FDR cutoffs for `--method logistic` |
| `--stim-fdr-threshold`, `--stim-log2-fc-threshold`, `--stim-delta-threshold` | Expansion FDR/fold-change and paired transcriptional delta cutoffs for `--method stimulation` (defaults `0.1`, `1.0`, `0.5`) |
| `--plot-filter` | Generate tier distribution plots (written to `--output`) |

---

### `tcrsift annotate`

Match clonotypes against VDJdb, IEDB, and CEDAR. Succeeds with default columns
(`db_match=False`, `is_viral=False`) when no database is supplied or cached.

```bash
tcrsift annotate -i filtered/tier1.csv -o annotated.csv \
    --vdjdb /path/to/vdjdb \
    --iedb /path/to/iedb \
    --cedar /path/to/cedar \
    --match-by CDR3ab \
    --exclude-viral
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input and output CSV paths (required) |
| `--vdjdb`, `--iedb`, `--cedar` | Database file/directory paths (auto-discovered from `tcrsift data` cache if omitted) |
| `--cache-dir`, `--auto-download` | Cache directory override and on-demand download |
| `--species` | Query TCR species: `auto` (default), `human`, or `mouse` |
| `--host-species-filter` | `all` (default; marks cross-species matches with `_cross`) or `same` (drop heterologous host entries) |
| `--database-match` | Chain-matching policy: `ab_with_partial` (default), `strict_ab` (require paired αβ on both query and database entry), or `b_only` |
| `--match-mode` | CDR3 comparison: `exact` (default) or `levenshtein` |
| `--max-distance` | Maximum edit distance when `--match-mode levenshtein` is set (default `1`) |
| `--match-by` | Legacy alias for `--database-match`: `CDR3ab` (`ab_with_partial`) or `CDR3b_only` (`b_only`) |
| `--exclude-viral` | Drop clonotypes matching viral epitopes |
| `--flag-only` | Add annotation columns without dropping viral rows |

---

### `tcrsift assemble`

Build full-length α, β, and single-chain (`β-linker-α`) protein and DNA
constructs.

```bash
tcrsift assemble -i annotated.csv -o sequences.csv \
    --contigs-dir /path/to/vdj \
    --alpha-leader CD28 \
    --beta-leader CD8A \
    --include-constant \
    --single-chain \
    --linker T2A \
    --fasta sequences.fasta
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input and output CSV paths (required) |
| `--contigs-dir` | Directory of Cell Ranger contig FASTAs (`*.fasta`, `*.fasta.gz`) or pre-CR6 `all_contig_annotations.json(.gz)` files (required for contig verification unless `--allow-canonical-fallback` is set, and required for `--leaders-from-contigs` / `from_contig` leaders) |
| `--allow-canonical-fallback` | Permit canonical constant-region assembly when Cell Ranger contig FASTAs are unavailable |
| `--alpha-leader`, `--beta-leader` | `CD8A`, `CD28`, `mCD8B`, `mCD28`, `H2Kb`, `IgK`, `TRAC`, `TRBC`, `from_contig`, or `none` (defaults: α `CD28`, β `CD8A`) |
| `--no-leaders` | Omit signal peptides on both chains |
| `--leaders-from-contigs` | Extract native V-gene leaders from Cell Ranger FASTAs |
| `--include-constant` | Append TRAC/TRBC constant regions |
| `--constant-source` | `canonical` (packaged human/mouse references; legacy alias `ensembl`) or `from-data` |
| `--species` | VDJ organism: `auto` (default), `human`, or `mouse` |
| `--constant-species`, `--codon-species` | Override constant-region or codon-optimization species (e.g. `--constant-species mouse` for chimeric human-VDJ / murine-constant constructs) |
| `--single-chain` | Emit β-linker-α single-chain construct (`single_chain_aa`, `single_chain_nt`) |
| `--linker` | 2A peptide: `T2A` (default), `P2A`, `E2A`, or `F2A` |
| `--germline-vdj-fallback` | Reconstruct missing 5' V-region prefixes (`partial_VDJ_*`) or bare `(v_gene, CDR3, j_gene)` chains from cached IMGT V-REGION + J-REGION references (`construct_contig_verified=False`) |
| `--auto-download-germline` | Automatically fetch missing IMGT V-REGION and J-REGION references into the cache when `--germline-vdj-fallback` is set |
| `--fasta`, `--airr` | Optional FASTA or AIRR output paths |

---

## Reference cache, sequence probability, and supplementary commands

### `tcrsift data`

Inspect, download, or clear cached public reference databases (VDJdb, IEDB,
CEDAR, and human/mouse IMGT `TRAV`/`TRBV` V-REGION, `TRAJ`/`TRBJ` J-REGION, and
leader references). Cache location resolves in order: `--cache-dir`,
`$TCRSIFT_DATA_DIR`, `$XDG_CACHE_HOME/tcrsift`, `~/.cache/tcrsift`.

```bash
tcrsift data list                             # Show database status, sizes, and paths
tcrsift data download                         # Download all auto-fetchable databases
tcrsift data download --group germline        # Download human IMGT V-REGION + J-REGION FASTAs
tcrsift data download --group germline_mouse  # Download mouse IMGT V-REGION, J-REGION + leader FASTAs
tcrsift data download --db iedb               # Download one database (--force refreshes)
tcrsift data clear --db iedb                  # Remove one database (omit --db for all)
```

---

### `tcrsift pgen`

Score clonotype CSVs with sequence generation (`Pgen`) and observed-repertoire
publicness (`Ppost`) log-probabilities, or manage custom TCRpeg/k-mer models.

```bash
tcrsift pgen annotate clones.csv -o clones_prob.csv --chain both --backend kmer
tcrsift pgen status
```

---

### `tcrsift match-til`

Annotate culture clonotypes with presence, cell counts, and frequencies across
one or more TIL datasets.

```bash
# Single TIL input (.h5ad, .csv, or Cell Ranger vdj directory)
tcrsift match-til -i clonotypes.csv --til-h5ad til.h5ad -o matched.csv
tcrsift match-til -i clonotypes.csv --til-csv til.csv -o matched.csv
tcrsift match-til -i clonotypes.csv --til-vdj-dir /path/to/vdj -o matched.csv

# Multiple TIL samples via sample sheet or direct NAME=TYPE:PATH specs
tcrsift match-til -i clonotypes.csv -s til_samples.yaml -o matched.csv
tcrsift match-til -i clonotypes.csv -o matched.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad \
  --til-sample T3=vdj:/path/to/til_t3_vdj_outs
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input culture CSV and output CSV (required) |
| `-s`, `--sample-sheet` | YAML/CSV sample sheet containing `source: til` entries |
| `--til-h5ad`, `--til-csv`, `--til-vdj-dir` | Single TIL source |
| `--til-sample` | Repeatable `NAME=TYPE:PATH` (`csv`, `h5ad`, `vdj`) |
| `--match-by` | `CDR3ab` (default) or `CDR3b_only` |
| `--min-til-cells` | Minimum TIL cells required to flag `til_match=True` (default `1`) |

Adds `til_match`, `til_samples`, `til_cell_count`, and `til_frequency` (plus
`til_cell_count.{sample}` and `til_frequency.{sample}` when multiple TIL samples
are given).

---

### `tcrsift til-clonotype`

Harmonize one or more TIL-only inputs into a single clonotype count and
frequency table.

```bash
tcrsift til-clonotype -o til_clonotypes.csv \
  --til-sample T1=csv:/path/to/til_t1.csv \
  --til-sample T2=h5ad:/path/to/til_t2.h5ad
```

Accepts `-s/--sample-sheet`, `--til-h5ad`, `--til-csv`, `--til-vdj-dir`, or
repeatable `--til-sample` specs, plus `--match-by` (`CDR3ab` or `CDR3b_only`)
and `--min-cells` (default `1`).

---

### `tcrsift til-select`

Prioritize TIL clonotypes across two or more longitudinal timepoints stored in
the legacy v2 per-timepoint file layout (`consensus_annotations.<TP>.csv`,
`clonotypes.<TP>.csv`, `filtered_contig_annotations.<TP>.csv`, and
`sample_filtered_feature_bc_matrix.<TP>.h5`). For standard Cell Ranger
directories, use [Clone Prioritization](til-signatures.md) (`tcrsift
prioritize --context solid-tumor`) instead.

`til-select`:

- harmonizes clonotypes across ordered timepoints,
- scores per-clone marker expression from GEX,
- applies a CD8, minimum-cell, and non-viral base mask,
- unions cytolytic, antigen-response, AIM, exhaustion, and frequency-increase
  selection branches, and
- writes v2-compatible CSV tables (`abTCR_harmonized.csv`,
  `figures/abTCR_master_table.csv`, `figures/selection_masks.csv`,
  `figures/subset_*.csv`) and PDF/PNG reports.

Frequency-increase branches are meaningful only when inputs form a genuine
time series, and the legacy `is_candidate_tumor_reactive` column denotes
shortlist membership rather than confirmed tumor specificity.

```bash
tcrsift til-select --data-dir data --rank-by marker_score_z_mean
```

| Option | Description |
| --- | --- |
| `--data-dir` | Directory containing per-timepoint files (default `./data`) |
| `--samples`, `--inputs` | Explicit `LABEL=CONSENSUS,CLONOTYPES` mappings |
| `--config` | YAML file with a `timepoints:` mapping |
| `--vdjdb`, `--iedb`, `--cedar` | Optional TCR reference databases (overrides cache) |
| `--cache-dir`, `--auto-download` | Reference cache directory and automatic fetching of missing public TCR databases |
| `--match-by` | Database matching key: `CDR3b_only` (default) or `CDR3ab` |
| `--marker-genes`, `--cytolytic-genes`, `--antigen-response-genes`, `--activation-genes`, `--exhaustion-genes` | Gene panels for branch scoring |
| `--min-cells-per-clone`, `--min-cd8-cp10k`, `--max-cd4-to-cd8-ratio` | Base CD8 selection thresholds |
| `--increase-ratio-nonzero-min`, `--immunogenic-percentile`, `--rank-by` | Branch cutoffs and top-k ranking metric |
| `--fig-dir`, `--out-table`, `--out-selected-report` | Output paths |

---

### `tcrsift annotate-gex`

Attach gene-expression columns from a standalone 10x `.h5` matrix to a VDJ CSV
table (useful when `gex_dir` was not included during `load`).

```bash
tcrsift annotate-gex \
    -i cells.csv \
    --gex-file filtered_feature_bc_matrix.h5 \
    -o clonotype_gex.csv \
    --aggregate \
    --cd4-cd8-counts
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input and output CSV paths (required) |
| `--gex-file` | 10x `filtered_feature_bc_matrix.h5` file (required) |
| `--barcode-col` | Barcode column name (default `barcode`) |
| `--genes` | Comma-separated gene symbols to extract |
| `--aggregate` | Aggregate expression per clonotype (`gex.{GENE}.sum`, `gex.{GENE}.mean`) |
| `--cd4-cd8-counts` | Add exclusive `CD4_only.count` and `CD8_only.count` columns |

---

### `tcrsift load-sct`

Load paired TCR tables from single-cell trimer (SCT) Excel workbooks.

```bash
tcrsift load-sct -i sct_data.xlsx -o sct_clonotypes.csv --aggregate
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input `.xlsx` and output `.csv` paths (required) |
| `--aggregate` | Collapse cell rows to clonotype level |
| `--min-snr` | Minimum signal-to-noise ratio (default `2.0`) |
| `--min-reads` | Minimum reads per chain (default `10`) |

---

### `tcrsift unify`

Merge clonotype tables from independent runs (`CDR3ab` or SCT `CDR3_pair`
tables) with source prefixes, `occurs_in_*` flags, and combined cell counts.

```bash
tcrsift unify \
    -i til_results/clonotypes.csv culture_results/clonotypes.csv \
    -o unified.csv
```

---

### `tcrsift mnemonic`

Generate pronounceable two-word labels from CDR3 sequences so related clonotypes
are easy to spot in tables and discussions.

```bash
tcrsift mnemonic -i clonotypes.csv -o clonotypes_with_names.csv
```

| Option | Description |
| --- | --- |
| `-i`, `--input` / `-o`, `--output` | Input and output CSV paths (required) |
| `--cdr3-col` | Source CDR3 column (auto-detected if omitted) |
| `--name-col` | Output column name (default `tcr_name`) |

---

### `tcrsift generate-config`

Write a documented YAML configuration file containing all default `tcrsift run`
parameters.

```bash
tcrsift generate-config -o config.yaml
```

Omit `-o` to print to stdout; add `--force` to overwrite an existing file.
