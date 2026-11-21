# Prioritization release notes

Behavior changes to `tcrsift prioritize` / `til-prioritize` that can alter a
shortlist. Record the package version with every analysis; the
[prioritization guide](user-guide/til-signatures.md) always describes the
current release. The complete per-release history is on
[GitHub](https://github.com/pirl-unc/tcrsift/releases) and
[PyPI](https://pypi.org/project/tcrsift/#history).

## 3.27.0

Two changes that alter shortlists, both prompted by the
[public-data evaluation](benchmarks.md#4-what-do-the-cell-gates-cost).

- The default `--min-cd3` is 3, down from 10, matching the `run` and
  `phenotype` commands. A floor of 10 removed 65% to 92% of CD8 cells in the
  datasets evaluated and did not rank the same clones better. More cells and
  clones are scored, so scores, cutoffs and shortlists change. Pass
  `--min-cd3 10` for the previous behavior.
- When the clone budget ends partway through a round, the remaining turns are
  shared across signatures and samples: each goes to the signature, then the
  sample, with the fewest turns in that round. Previously they went to the
  first lists in signature order, so in a run with more qualifying lists than
  budget the signatures listed last in the preset could supply no candidates.
  Only the final, partial round of a selection changes; full rounds are
  visited in the same order as before.

## 3.26.0

Fixes found by running on public data ([evaluation](benchmarks.md)). None
changes scores or shortlists for inputs that already ran.

- Multi-sample runs that read `filtered_feature_bc_matrix/` directories keep
  the gene symbols supplied with the matrices. Previously
  `--context solid-tumor` could stop with "genes absent from matrix".
- Signature genes are matched across reference-annotation renames through a
  small synonym table (for example PLA2G16/PLAAT3 and CARS/CARS1), so NeoTCR8,
  NeoTCR4 and NeoTCR_PBL score on current Cell Ranger references. A gene with
  no synonym present is still an error.
- When a sample has only `all_contig_annotations.csv`, the loader keeps called
  cells and high-confidence contigs and drops contigs with no chain call.
- The command warns when any single cell filter other than lineage removes
  more than half of the remaining cells. The default `--min-cd3 10` is
  unchanged.

## 3.25.0

Background selection now defaults to either mean or high-cell subset evidence,
with separate ranked lists inside the same 200-total budget. Exports include
high-cell counts, fractions, exact enrichment tails, and the evidence route
that supplied each candidate. `--signature-evidence mean` restores the 3.24.0
selection rule; legacy mode retains its previous behavior.

## 3.24.0

The default stopping rule now uses a size-matched empirical background floor
instead of score > 0. The total cap remains 200 and the rank gate remains the
top decile. Use `--signature-cutoff legacy` for the previous score rule. This
can substantially reduce shortlists, especially with very small samples.
NeoTCR module scores now explicitly read normalized `.X` even when `.raw`
exists. Their registry metadata and descriptions correctly identify log-space
control subtraction; the bare-DataFrame proxy follows the other signatures'
default `log1p` transform (pass `log1p=False` for already logged frames).

## 3.23.1

The cap is now **200 total**, replacing 100 per patient. Selection ranks are
global. That release used percentile ≥ 0.90 and score > 0; numerically flat
lists are excluded. Exhausted lists yield their turns to the remaining lists.
There is no fallback below these cutoffs to fill the budget.

## Changes from 3.22

Both command names now use context-dependent stratification, default CD8
selection, viral exclusion off, and a disabled mitochondrial floor. MART-1
exclusion remains on. One sample is sufficient. These changes intentionally
change shortlists from the original TIL example; the old ranking algorithm
is not restored by adjusting thresholds. Use the recorded configuration and
package version when reproducing a prior analysis.
