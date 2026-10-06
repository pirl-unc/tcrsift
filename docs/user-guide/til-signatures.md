# Clone prioritization by expression state

Use `tcrsift prioritize` for paired CellRanger VDJ + GEX from one or more
samples. The command selects clones from the top of several expression-state
lists. It does not infer antigen specificity or prove that a clone was induced
by vaccination. Keep patient identities in the sample sheet.

```bash
# Vaccine / peripheral-blood study; CD8 is the default.
tcrsift prioritize samples.yaml -o candidates/ --context blood

# Malignant pleural effusion, including both CD4 and CD8.
tcrsift prioritize samples.yaml -o candidates/ --context mpe --tcell-type both

# AML / other hematologic malignancy.
tcrsift prioritize samples.yaml -o candidates/ --context heme

# Solid-tumor TILs: these two commands are equivalent.
tcrsift prioritize samples.yaml -o candidates/ --context solid-tumor
tcrsift til-prioritize samples.yaml -o candidates/
```

`til-prioritize` and `examples/multi_sample_til.py` use the same workflow with
`solid-tumor` as their default context. `--context` can override it. The older
`til-select` command remains a separate workflow for its legacy ordered-input
layout; independent samples should not be interpreted as a time series.

## Simple selection rule

1. Apply cell QC, CD4/CD8 selection, paired-chain quality, and optional marker
   expression gates. Use this same cell population for counts and scores.
2. Score each enabled signature separately within each **sample and lineage**,
   using log1p(CP10K). Average scores over cells of each clone in that stratum.
3. Remove clones failing abundance or explicit exclusion filters. Rank each
   signature/sample/lineage list by score, then frequency, cell count, and CDR3
   sequence to break ties deterministically.
4. Take the next unused clone from each list in turn, skipping previously
   selected clones. Repeat until the patient has **100 unique clones** or the
   lists are exhausted. `--max-clones N` changes that budget; `0` removes it.

Every signature gets a turn; there is no weighted composite or fitted ranking
model. Lists are visited in the signature order below (or your `--signatures`
order), then lexical sample and lineage order. A small budget can end partway
through a round. Overlapping signatures are correlated views, not independent
votes. A clone shared between samples occupies one slot within a patient; the
same sequence in two patients occupies a slot for each patient.

There is **no percentile cutoff by default**: this is top-of-list sampling.
Use `--signature-quantile 0.90` to restrict lists to their top decile, optionally
with `--min-signature-support 2`. Support counts distinct signatures with at
least one qualifying sample/lineage; it does not require them in the same cell.
Tied scores use average percentile ranks. Small or constant-score strata can
therefore have no clone above a requested percentile cutoff.

## Context presets

The shared core, in order, is **Differentiated, AntigenExperienced, Cytolytic,
AcuteActivation, Proliferation, AIM**. This covers differentiation, effector
function, recent activation, proliferation, and an AIM-like state. We use one
panel per named signature, without also adding each `Broad` variant by default.

| Context | Additions to the core | Rationale |
| --- | --- | --- |
| `generic` (default) | None | Broad state coverage when the setting is unspecified. |
| `blood` | CirculatingMemory | Blood/vaccine studies need memory as well as activated and effector strata. No tumor-trained scores. |
| `blood-tumor` | CirculatingMemory, NeoTCR_PBL | Explicit search for circulating antitumor CD8 cells; not a generic vaccine preset. |
| `solid-tumor` | TumorReactive, MANAscore, NeoTCR8, NeoTCR4 | Include published tumor-reactivity programs alongside broad states. |
| `mpe` | CirculatingMemory, MANAscore, NeoTCR8, NeoTCR4 | Exploratory transfer of tumor-reactivity programs to malignant pleural effusions, plus memory; omit the epithelial-residency-containing TumorReactive panel. |
| `heme` | CirculatingMemory | Conservative state coverage for AML/hematologic malignancies. Intentionally the same gene sets as blood; no validated leukemia-specific ranking is shipped. |

These are explicit, literature-informed starting points, not validated
context-specific classifiers. Context applies to the entire invocation; it is
not inferred from filenames, `source`, or free-text metadata. Run distinct
contexts separately or supply an explicit common signature set for a mixed study.

**Lineage rule:** `--tcell-type cd8` is the default. `cd4` and `both` are also
available. Confident and likely cells of the requested type are retained;
unknown cells are excluded. NeoTCR4 is CD4-only; NeoTCR8, NeoTCR_PBL, and
MANAscore are CD8-only. Incompatible signatures are removed from presets.
With `both`, each is scored only in its applicable lineage, and CD4/CD8 have
separate lists. Explicitly requesting an incompatible signature is an error.

Override a preset or remove a panel:

```bash
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --signatures CirculatingMemory AcuteActivation Cytolytic --max-clones 60

tcrsift prioritize samples.yaml -o candidates/ --context solid-tumor \
  --exclude-signatures MANAscore --tcell-type both
```

Names are case-insensitive. Other registered panels, including `Broad` variants,
can be selected explicitly. Required genes missing from the expression matrix
cause an error; signatures are not silently shortened. The [signature registry](../api/signatures.md)
documents gene sets and scoring methods. In particular, MANAscore is TCRsift's
transparent signed-z **proxy**, not the original trained predictor.

## Input and outputs

```yaml
samples:
  - sample: patient1_pre
    patient_id: patient1
    tissue: blood
    timepoint: pre
    vdj_dir: /data/patient1_pre/vdj
    gex_dir: /data/patient1_pre/gex
  - sample: patient1_post
    patient_id: patient1
    tissue: blood
    timepoint: post
    vdj_dir: /data/patient1_post/vdj
    gex_dir: /data/patient1_post/gex
```

CSV sample sheets also work. Sample names must be unique and every sample must
have both `vdj_dir` and `gex_dir`. If `patient_id` is used, populate it for every
sample. If omitted entirely, all samples are treated as one patient/cohort.
`tissue` is descriptive; `--context blood` chooses the preset. `source` can be
omitted: its legacy pipeline enum is not a tissue or prioritization setting.

- `candidate_clones.csv`: the shortlist, with per-patient `selection_rank` and
  the `selected_signature`, `selected_sample`, and `selected_lineage` that
  supplied each selection.
- `all_scored_clones.csv`: all clones after cell filtering, their scores,
  `eligible_for_review`, `selected_for_review`, risk flags, and `excluded_reason`.
  An eligible clone can be unselected because of the clone budget.
- `clone_sample_scores.csv`: clone/sample/lineage scores, percentiles, counts,
  and frequencies. Frequencies here use retained paired cells of that lineage
  in that sample. Clone-level `max_frequency` uses all retained requested
  lineages in the sample; these coincide when only one lineage is selected.
- `prioritization.json`: version, arguments, resolved signature set, sequential
  cell-filter counts, retained counts per sample, and database availability.
  Counts start after GEX QC; loader logs report GEX QC losses.

Cells removed by QC are not present in clone tables. All scores and thresholds
are cohort-relative; the shortlist is for experimental review.

## Filters and units

| Options | Default and meaning |
| --- | --- |
| `--min-cells`, `--min-frequency` | 2 cells across a patient's samples; frequency at least 0.001 in one sample. |
| `--min-genes`, `--max-genes` | 250–15,000 detected genes per cell. |
| `--min-counts`, `--max-counts` | 500–100,000 total GEX UMI counts per cell, not sequencing reads. |
| `--min-mito-pct`, `--max-mito-pct` | 0–8% mitochondrial counts; no mitochondrial floor by default. |
| `--min-vdj-umis`, `--min-vdj-reads` | At least 2 UMIs in each primary alpha/beta chain; read gate disabled (0). |
| `--min-cd3` | At least 10 summed CD3D/E/G raw GEX counts per cell; 0 disables. |
| `--min-expression GENE=VALUE`, `--max-expression GENE=VALUE` | Optional per-cell log1p(CP10K) thresholds. Repeat for AND conditions. Missing genes are errors. |
| `--include-v-gene`, `--exclude-v-gene` | Optional exact V calls in either chain, ignoring allele suffix and case. Repeat for multiple genes. |
| `--include-j-gene`, `--exclude-j-gene` | Corresponding J-gene gates. Inclusion requires at least one listed gene; exclusion rejects any listed gene. Missing calls fail active gene filters. |
| `--min-alpha-cdr3-length`, `--min-beta-cdr3-length` | Optional amino-acid length floors, including conserved anchors; both disabled (0). |

Gene-call filters use the representative primary-chain calls in the clone table.
Only complete paired alpha/beta cells with canonical amino-acid CDR3s enter the
workflow; dual-chain cells use their primary chains. VDJ and expression gates
are applied before aggregation, so rejected cells cannot inflate abundance or
contribute signature scores.

```bash
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --min-genes 300 --max-mito-pct 15 --min-vdj-reads 10 \
  --min-expression GZMB=1 --max-expression CCR7=2 \
  --min-alpha-cdr3-length 10 --min-beta-cdr3-length 12
```

These thresholds illustrate the syntax, not universal biological cutoffs.
Expression gates can remove the memory/state diversity that stratification is
intended to preserve; use them when that restriction is intentional.

## MART-1, viral matches, and publicness

**Known MART-1 matches are excluded by default in every context.** Supply one
or more reference files with `--vdjdb`, `--iedb`, or `--cedar`. Matches to
MART-1/Melan-A/MLANA and EAAGIGILTV/AAGIGILTV/ELAGIGILTV are flagged. Without a
reference database the command warns: it cannot identify unannotated MART-1
clones. A missing database match means unknown specificity. Override explicitly
with `--no-exclude-known-mart1` for studies that target MART-1.

**Known viral matches are retained by default.** Add `--exclude-known-viral`
when they should be excluded. This is intentionally opt-in for vaccine studies
and studies of viral tumor antigens. Database match strictness is controlled by
`--database-match strict_ab|ab_with_partial|b_only` (default `ab_with_partial`).

`--exclude-trav12-2` is an optional aggressive heuristic. TRAV12-2 bias is not
proof of MART-1 specificity, and CDR3 length cannot establish it either.
`--exclude-public-quantile 0.90` optionally removes the top publicness decile
within each patient (ties can remove a different fraction). Scores and flags
are retained in the audit table even when these exclusions are disabled.

Length filters and the bundled k-mer publicness model are **not a germline
distance or an N-insertion count**. Exact junctional analysis would require
nucleotide sequences and germline V(D)J alignment; beta junctions also contain
D-segment contributions. This command does not estimate that quantity. Short
or readily generated receptors can still be antigen-specific, so length and
publicness exclusions remain opt-in.

## Evidence behind the presets

- [Lowery et al., Science 2022](https://pubmed.ncbi.nlm.nih.gov/35113651/):
  separate CD4 and CD8 neoantigen-reactive TIL signatures.
- [Yossef et al., Cancer Cell 2023](https://pubmed.ncbi.nlm.nih.gov/38039963/):
  circulating antitumor CD8 cells have a memory-like phenotype distinct from
  their TIL counterparts; basis for the explicit `blood-tumor` preset.
- [Vaccination/infection profiling, Nature Immunology 2023](https://www.nature.com/articles/s41590-023-01608-9):
  supports examining multiple antigen-specific CD8 states rather than treating
  cytotoxicity alone as vaccine specificity.
- [MPE single-cell study, Nature Communications 2021](https://www.nature.com/articles/s41467-021-27026-9):
  heterogeneous effusion immune states. Our MPE panel choice is an exploratory
  extrapolation, not a validated classifier from that study.
- [AML differentiation/dysfunction study, Blood 2024](https://pubmed.ncbi.nlm.nih.gov/38776511/):
  motivates preserving memory and differentiated strata; does not validate our
  `heme` preset or establish leukemia specificity.
- [MART-1 germline recognition study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2785656/)
  and [Stitchr](https://pmc.ncbi.nlm.nih.gov/articles/PMC9262623/): context for
  distinguishing germline V-gene bias, CDR3 length, and junctional diversity.

## Changes from 3.22

Both command names now use stratification, default CD8 selection, a 100-clone
per-patient budget, viral exclusion off, and a disabled mitochondrial floor.
MART-1 exclusion remains on. One sample is now sufficient. The signature sets
are context-dependent and the default top-decile requirement is removed.
These changes intentionally change shortlists from the original TIL example;
`--tcell-type both --max-clones 0 --signature-quantile 0.90` restores those
individual controls but does not restore the old ranking algorithm. Use the
recorded configuration and package version when reproducing a prior analysis.
