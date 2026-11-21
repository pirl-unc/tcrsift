# Filtering Strategies

`tcrsift filter` assigns clonotypes from antigen-stimulated cultures to
confidence tiers (`tier1`–`tier5`) using cell counts, within-sample
frequencies, T-cell lineage, and optional viral exclusions. Tiering ranks
clonal expansion for experimental follow-up; it does not prove epitope
specificity.

## Sources of clones in stimulated cultures

Peptide- or antigen-stimulated T-cell cultures typically contain a mixture of:

- **Antigen-responsive clones** that expanded upon recognizing the target
  antigen.
- **Bystander clones** (often common anti-viral memory T cells) carried over or
  expanded non-specifically by cytokines.
- **Background singletons and low-count clones** from unexpanded input T cells
  or ambient/doublet artifacts.

## Threshold method (default)

The `threshold` method assigns each clone to the strictest tier whose minimum
cell count (`cell_count`) and peak within-sample frequency (`max_frequency`) it
meets:

| Tier | Min cells | Min frequency | Priority |
| --- | --- | --- | --- |
| `tier1` | 10 | 1.0% (`0.01`) | Highest |
| `tier2` | 5 | 0.5% (`0.005`) | High |
| `tier3` | 3 | 0.1% (`0.001`) | Medium |
| `tier4` | 2 | 0.05% (`0.0005`) | Low |
| `tier5` | 2 | 0% | Lowest |

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold --tcell-type cd8
```

Bundled tier definitions use abundance only. In Python, custom
`tier_definitions` dictionaries may also include `max_conditions` to cap how
many stimulation conditions a clone may appear in:

```python
from tcrsift import filter_clonotypes

custom_tiers = {
    "tier1": {"min_cells": 20, "min_frequency": 0.02, "max_conditions": 1},
    "tier2": {"min_cells": 10, "min_frequency": 0.01, "max_conditions": 2},
    "tier3": {"min_cells": 5, "min_frequency": 0.005, "max_conditions": 3},
}

filtered = filter_clonotypes(clonotypes, tier_definitions=custom_tiers)
```

## Logistic (FDR) method

The `logistic` method fits a logistic regression model relating clone frequency
to non-viral (or single-condition) status among abundant clones, converts model
probabilities into estimated false discovery rates, and assigns clones to FDR
tiers:

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method logistic \
    --fdr-tiers 0.0001,0.001,0.01,0.1,0.15
```

| Tier | Default FDR cutoff |
| --- | --- |
| `tier1` | 0.0001 (0.01%) |
| `tier2` | 0.001 (0.1%) |
| `tier3` | 0.01 (1%) |
| `tier4` | 0.10 (10%) |
| `tier5` | 0.15 (15%) |

!!! note "When logistic fitting is unstable"
    Logistic calibration requires annotated viral/background clones or multiple
    conditions across a sufficient range of frequencies. On small or weakly
    expanded samples where fitting cannot converge, `filter_clonotypes_logistic`
    falls back to a fixed frequency threshold; prefer `--method threshold` for
    single-sample or unannotated runs.

In multi-donor `tcrsift run` cohorts, `--fdr-scope` defaults to `auto`, which
computes separate tierings per donor (`filtered_tier{1..5}_<donor>.csv`) unless
`donors_share_antigen: true` is set on the sample sheet or `--fdr-scope global`
is passed explicitly.

## Multi-method overlap and recovery (`enrichment_method`)

When `enrichment_method` (such as `AIMpos`, `tetpos`, or `CTYneg`) is populated
on the sample sheet, `tcrsift run` writes per-method rankings and cross-method
comparison tables alongside the tier files:

- `data/filtered_by_method/<donor>__<method>.csv`: top clones per donor and
  method (`--per-method-top-n`).
- `data/method_overlap_<donor>.csv` and `plots/method_overlap_<donor>.png`:
  pairwise method overlap matrix (`--method-overlap-similarity
  {jaccard,dice,count}`) and heatmap.
- `data/method_overlap_sets.csv` and `plots/method_overlap_upset.png`: multi-way
  clone intersection counts and UpSet/bar chart.
- `data/method_recovery.csv` and `plots/method_recovery.png`: fraction of the
  strictest populated tier recovered by each method per donor.

## Lineage and viral exclusions

Restrict candidates to the T-cell lineage expected for your antigen format
(based on each clone's majority cell phenotype, `Tcell_type_consensus`):

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type cd8   # MHC-I / short peptides
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type cd4   # MHC-II
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type both  # long peptides / whole protein
```

To remove cataloged bystander viral TCRs (CMV, EBV, influenza, HIV, SARS-CoV-2),
annotate against VDJdb/IEDB first and pass `--exclude-viral`:

```bash
tcrsift annotate -i clonotypes.csv -o annotated.csv \
  --vdjdb /path/to/vdjdb --iedb /path/to/iedb --flag-only

tcrsift filter -i annotated.csv -o filtered/ --exclude-viral
```

(`tcrsift annotate ... --exclude-viral` also performs annotation and removal in
a single step.) Absence of a database hit means the sequence is uncataloged,
not proven non-viral; paired `CDR3ab` matching is much more specific than
`CDR3b_only`.

## Cross-checking against matched TIL

When tumor repertoire data are available for the same patient, clones that both
expand in antigen-stimulated culture and appear in TIL carry stronger evidence
for experimental testing:

```bash
tcrsift filter -i culture_clonotypes.csv -o filtered/ --method threshold
tcrsift match-til -i filtered/tier1.csv --til-csv til_clonotypes.csv -o matched.csv
```

(In `tcrsift run`, listing a `source: til` sample on the same sample sheet runs
this step automatically.)
