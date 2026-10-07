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
3. Keep only the informative head of each signature/sample/lineage list:
   by default **percentile ≥ 0.90 and score > 0**. A numerically flat list
   supplies no clones. Cutoffs are evaluated on all scored clones in that
   stratum, before abundance/database exclusions or selection.
4. Remove clones failing abundance or explicit exclusion filters. Rank each
   remaining list by score, then frequency, cell count, and CDR3 sequence.
5. Take the next unused clone from each list in turn. When a list runs out of
   qualifying clones, skip it and keep taking from the others. Stop at
   **200 candidate rows total across the entire run**, or when every list is
   exhausted. Never pad the shortlist with below-cutoff clones.

`--max-clones N` changes the total cap; `0` removes the cap but keeps the
signature cutoffs. This is one shared budget across patients, samples,
signatures, and lineages. A clone shared between samples occupies one slot
within a patient; the same sequence in two patients remains two separate rows,
both counted toward the total cap.

Every active list gets a turn; there are no fixed signature quotas or weighted
composite scores. For example, if one list has only 8 qualifying clones, it
stops after those 8; another can supply the remaining slots while it still has
qualifying clones. If only 137 distinct candidates qualify overall, the result
has 137 rows, even with a cap of 200.

Lists are visited in signature order (preset order below, or your
`--signatures` order), then lexical patient, sample, and lineage order. A small
budget can end partway through a round. Overlapping signatures are correlated
views, not independent votes. Duplicate selections are skipped before taking
the next qualifying clone from the same list.

The cutoff is a **transparent ranking heuristic**, not a significance test or
proof of antigen specificity. Positive scores and a varying list avoid taking
zero/flat programs just to fill the budget; the percentile floor limits each
list to its strongest relative scores. A noisy signature can still pass.

Adjust the stopping rule explicitly:

```bash
# Broaden each list to its top 20%, retaining the positive-score requirement.
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --max-clones 200 --signature-quantile 0.80 --min-signature-score 0
```

`--signature-quantile 0` disables only the percentile gate. The score cutoff
is strict (`score > --min-signature-score`), accepts finite negative values
for an intentional relaxation, and is applied within each signature's native
score units. Flat and single-clone strata are always skipped because they
provide no within-stratum ranking contrast. Tied scores use average percentile
ranks, so a tied top group can fall below a high percentile cutoff.

`--min-signature-support 2` optionally requires two distinct signatures passing
**both** cutoffs. Each supporting signature must pass both in the same
sample/lineage; different signatures may qualify in different strata. A clone
can be taken only from a list where it passes that list's cutoffs.

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

- `candidate_clones.csv`: the shortlist, with global `selection_rank` and
  the `selected_signature`, `selected_sample`, and `selected_lineage` that
  supplied each selection.
- `all_scored_clones.csv`: all clones after cell filtering, their scores,
  `eligible_for_review`, `selected_for_review`, risk flags, and `excluded_reason`.
  An eligible clone can be unselected because of the clone budget.
- `clone_sample_scores.csv`: clone/sample/lineage scores, percentiles, counts,
  and frequencies. Each `signature_NAME_passes_cutoff` column records
  qualification before clone exclusions and budgeting. Frequencies here use retained paired cells of that lineage
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

## Changes from 3.23.0

The cap is now **200 total**, replacing 100 per patient. Selection ranks are
global. The default cutoff is percentile ≥ 0.90 and score > 0; numerically flat
lists are excluded. Exhausted lists yield their turns to the remaining lists.
There is no fallback below these cutoffs to fill the budget.

## Changes from 3.22

Both command names now use context-dependent stratification, default CD8
selection, viral exclusion off, and a disabled mitochondrial floor. MART-1
exclusion remains on. One sample is sufficient. These changes intentionally
change shortlists from the original TIL example; the old ranking algorithm
is not restored by adjusting thresholds. Use the recorded configuration and
package version when reproducing a prior analysis.
