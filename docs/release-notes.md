# Prioritization release notes

Behavior changes to `tcrsift prioritize` / `til-prioritize` that can alter a
shortlist. Record the package version with every analysis; the
[prioritization guide](user-guide/til-signatures.md) always describes the
current release. The complete per-release history is on
[GitHub](https://github.com/pirl-unc/tcrsift/releases) and
[PyPI](https://pypi.org/project/tcrsift/#history).

## 3.25.4

When a sample has only `all_contig_annotations.csv`, the loader now keeps
called cells and high-confidence contigs, as Cell Ranger's filtered table does,
and drops contigs with no chain call. Previously every row was used and the
run stopped with "Missing TRA_1_umis". Samples with
`filtered_contig_annotations.csv` are unaffected.

## 3.25.3

The command now warns when any single cell filter other than lineage removes
more than half of the remaining cells. The default `--min-cd3 10` is unchanged.
It is a raw-count floor, so it behaves as a sequencing-depth filter: in public
lung-tumor data (GSE176021, about 3,600 UMIs per cell) it removed 72% of
CD8-classified cells. Lowering it to 3 kept 86% of them but ranked validated
neoantigen-specific clones worse (TumorReactive AUROC 0.69 against 0.97, with
10 to 18 positive clones), so the stricter default stays and the loss is made
visible instead.

## 3.25.2

Signature genes are now matched across reference-annotation renames. NeoTCR8
lists PLA2G16 and CARS, which Cell Ranger 2020-A and later references call
PLAAT3 and CARS1, so `--context solid-tumor` stopped with a missing-gene error
on current references. A gene with no synonym present is still an error.
Scores do not change for inputs that already ran.

## 3.25.1

Multi-sample runs that read `filtered_feature_bc_matrix/` directories now keep
the gene symbols supplied with the matrices. Previously the symbols were lost
when samples were combined and looked up in a different Ensembl release, so
`--context solid-tumor` could stop with "genes absent from matrix" (NeoTCR8)
for genes the matrices contained. Single-sample and `.h5` inputs were not
affected. Found by the [public benchmark](benchmarks.md).

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
