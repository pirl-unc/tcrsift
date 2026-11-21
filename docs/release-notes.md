# Prioritization release notes

Behavior changes in `tcrsift prioritize` and `til-prioritize` that can alter a
shortlist. Record the package version alongside every analysis; the
[prioritization guide](user-guide/til-signatures.md) documents the current
release. Full release history is on
[GitHub](https://github.com/pirl-unc/tcrsift/releases) and
[PyPI](https://pypi.org/project/tcrsift/#history).

## 3.28.0

Bug fixes affecting receptor assignment and shortlist selection. Re-run
analyses generated with earlier versions.

- **Non-productive contigs are dropped on load.** Cell Ranger's filtered contig
  table retains non-productive contigs, and the loader previously ranked a
  cell's contigs by UMI count alone. A high-UMI non-productive contig could
  therefore displace a cell's productive alpha or beta chain or leave the cell
  without a CDR3 (7% of primary chains in the public lung dataset were
  non-productive). Applies to every command that loads Cell Ranger VDJ output,
  including `run`.
- **Barcodes are matched exactly unless no exact matches exist.** Stripping
  `-1`/`-2` suffixes per barcode could previously pair GEX cell `X-2` with the
  TCR of an unrelated VDJ cell `X-1`. Core-barcode fallback now runs only when
  a sample has zero exact barcode matches, skips ambiguous core barcodes, and
  warns when nothing matches.
- **`.h5` expression matrices are indexed by Ensembl ID**, matching matrix
  directories. Duplicate gene symbols previously caused `.h5`-only samples to
  fail on load, and mixing `.h5` and matrix-directory samples duplicated the
  gene axis.
- **Mid-round budget sharing and `--max-clones 0` fixes.** A list whose sole
  remaining clone had already been selected by another list could still consume
  a turn, violating prefix consistency across `--max-clones` budgets. In
  addition, `--max-clones 0` was treated as a cap equal to the eligible clone
  count, altering `selection_rank` and `selected_signature` relative to a
  non-binding cap.
- **V and J gene filters with repeated barcodes across samples.** When a
  clone's highest-UMI cell shared its barcode string with another cell in a
  different sample, representative gene lookup returned multiple rows and
  caused `--exclude-trav12-2`, `--exclude-v-gene`, and `--include-v-gene` to
  misclassify the clone.
- **Sample-sheet parsing.** Blank `source` cells fall back to the default
  instead of raising an error, numeric YAML sample names are coerced to strings,
  and duplicate sample names are rejected on load.

## 3.27.0

Two shortlist-changing updates motivated by the
[public TIL benchmark](benchmarks.md#4-what-do-the-cell-filters-drop):

- **Default `--min-cd3` lowered from 10 to 3**, matching `run` and `phenotype`.
  Across public TIL datasets, a floor of 10 discarded 65% to 92% of CD8 cells
  without improving ranking among surviving clones. Because more cells and
  clones enter scoring, signature percentiles, background floors, and
  shortlists change. Pass `--min-cd3 10` to reproduce earlier filtering.
- **Fair mid-round budget allocation.** When `--max-clones` expires partway
  through a round, remaining turns are assigned to the signature (and within
  it, the sample) with the fewest picks so far in that round. Previously, partial
  rounds were truncated in preset order, which could starve signatures listed
  near the end of a preset when qualifying lists outnumbered the budget. Full
  rounds retain the existing order.

## 3.26.0

Loader and scoring fixes uncovered while benchmarking on public datasets
([benchmarks](benchmarks.md)). Existing runs that completed without error
produce unchanged scores and shortlists.

- Multi-sample runs reading `filtered_feature_bc_matrix/` directories preserve
  the gene symbols supplied with the matrices instead of failing with missing
  signature genes under `--context solid-tumor`.
- Signature lookup resolves common HGNC symbol renames via a built-in alias
  table (such as `PLA2G16`/`PLAAT3` and `CARS`/`CARS1`), allowing `NeoTCR8`,
  `NeoTCR4`, and `NeoTCR_PBL` to score across older and newer Cell Ranger
  references. Genuinely absent genes still raise an error.
- When only `all_contig_annotations.csv` is present, the VDJ loader filters to
  called cells and high-confidence contigs with valid chain calls.
- `prioritize` logs a warning whenever any single non-lineage cell filter
  removes more than half of the remaining cells.

## 3.25.0

Background cutoff mode defaults to `--signature-evidence either`, allowing
clones to qualify via elevated mean expression or hypergeometric enrichment of
high-scoring cells, each contributing its own ranked list within the shared
200-clone budget. Output tables add high-cell counts, fractions, exact tail
probabilities, and `selected_evidence`. Pass `--signature-evidence mean` for
3.24.0 behavior; `legacy` cutoff mode continues to use mean scores only.

## 3.24.0

Replaced the fixed `score > 0` floor with a size-matched empirical resampling
floor (`--signature-cutoff background`), keeping the top-decile rank gate and
200-clone cap. Pass `--signature-cutoff legacy` for the previous rule. NeoTCR
module scoring now always reads normalized `.X` even when `.raw` is present,
and bare-DataFrame fallback scoring applies `log1p` by default (pass
`log1p=False` if inputs are already log-transformed).

## 3.23.1

Changed `--max-clones` from 100 per patient to **200 across the entire run**
with global `selection_rank`. Required percentile ≥ 0.90 and score > 0,
skipping numerically flat lists and yielding turns from exhausted lists without
padding below cutoffs.

## 3.23.0 (changes from 3.22)

Unified `prioritize` and `til-prioritize` around context presets, CD8 selection
by default, MART-1 exclusion on, viral exclusion off, no minimum mitochondrial
floor, and single-sample support.
