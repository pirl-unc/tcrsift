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
3. Keep clones with **percentile ≥ 0.90 and score above a size-matched
   background floor** (details below). A numerically flat list supplies no
   clones. Calibration uses all retained cells in that stratum, before clone
   abundance/database exclusions or selection.
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

## Per-signature background floors

The default `--signature-cutoff background` estimates a separate floor for
each **signature, patient, sample, lineage, and clone cell count**:

1. For a clone represented by `n` cells in that sample/lineage, draw random
   groups of `n` retained cells from the same patient/sample/lineage, uniformly
   **without replacement within a draw**. All retained cells, including the
   candidate clone, are part of this reference population.
2. Average the already computed per-cell scores in each random group, exactly
   as for the real clone. Use 2,000 draws, shared across signatures. Enumerate
   all possible groups when there are fewer; singleton references use all cells.
3. Require the clone mean to **exceed the 99th percentile** of those reference
   means, using the higher observed order statistic and treating numerical ties
   as failures. The top-decile rank gate is an additional requirement.

This preserves the observed score distribution, including skew, dropout and
correlated genes, without assuming Gaussian scores or inventing a universal
native-unit cutoff. A clone with two cells is compared with two-cell groups;
a clone with twenty cells is compared with twenty-cell groups. Larger clones
usually have narrower reference distributions. Seed 0, stable cell identifiers,
and sorted strata make repeated runs deterministic within a software environment.

**Interpretation:** this is a competitive expression-enrichment reference:
"higher than random groups of T cells in this sample." It does not estimate
pure technical noise or establish antigen experience, vaccine induction, or
antigen specificity. Cell exchangeability within a sample/lineage is an
assumption; remaining cell-state, quality and batch effects can drive enrichment.
A uniformly activated population has no contrast to itself. Including truly
activated clones in the reference can make the rule conservative. An extremely
small sample may supply no clones even when its top clone looks distinctive.

The 99th percentile is a **configurable operating choice**, not a validated
biological boundary. Noise can still pass. The exported one-sided background
tail probability includes ties and uses `(exceedances + 1) / (draws + 1)` for
Monte Carlo sampling (the exact fraction for enumeration). These probabilities
are **unadjusted across clones/signatures/samples**; the shortlist has no claimed
false-discovery-rate control. Stable identifiers and the recorded package
version, seed, draw count and floor table support reproduction.

### What the papers support

| Source | Published approach | Consequence for TCRsift |
| --- | --- | --- |
| [Zeng et al., MANAscore, 2025](https://www.nature.com/articles/s41467-024-55059-3) | Patient-specific last troughs of the imputed and non-imputed trained ensemble score distributions; clone calls required at least five MANAscore-high cells. | Supports adapting thresholds to a population. Neither the numeric thresholds nor the fitted probability scale transfer to TCRsift's three-gene signed-z **proxy**. |
| [Lowery et al., NeoTCR4/8, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8996692/) | scGSEA scores, validated against experimentally tested TCRs; figures highlight the 95th percentile. | A percentile is not a universal noise boundary. TCRsift uses the published gene sets with a different scoring implementation. |
| [Yossef et al., NeoTCR_PBL, 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10843665/) | A circulating tumor-reactive gene signature with functional validation in cancer patients. | Blood-specific biological support does not calibrate thresholds for generic blood, vaccine or heme samples. |

TCRsift's NeoTCR AnnData scores use [Scanpy `score_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html):
mean signature-gene expression minus expression-matched control-gene expression.
This is **not scGSEA or a rank-based enrichment statistic**. No paper-derived
numeric cutoff is imported into this different scale. The resampling rule above
is TCRsift's implementation choice, not a reproduction of those papers' methods.

Adjust the stopping rule explicitly:

```bash
# Broaden the rank gate while keeping the default background floor.
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --max-clones 200 --signature-quantile 0.80

# More conservative background tail, with more draws for resolution.
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --signature-background-quantile 0.995 --signature-background-draws 4000

# Explicitly reproduce the 3.23.1 cutoff behavior.
tcrsift prioritize samples.yaml -o candidates/ --context blood \
  --signature-cutoff legacy --signature-quantile 0.90 --min-signature-score 0
```

`--signature-quantile 0` disables only the percentile gate.
`--min-signature-score VALUE` adds a strict native-unit floor on top of the
background floor. There is no default universal zero cutoff in background mode;
a negative score can exceed a more-negative reference floor. In `legacy` mode,
the minimum replaces the default zero floor and accepts finite negative values.
`--signature-background-seed` changes the random seed. Draw count must be at
least 1,000 and leave at least ten expected draws above the requested quantile.
Flat and single-clone strata are always skipped because they
provide no within-stratum ranking contrast. Tied scores use average percentile
ranks, so a tied top group can fall below a high percentile cutoff.

`--min-signature-support 2` optionally requires two distinct signatures passing
**all** cutoffs. Each supporting signature must pass them in the same
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
  and frequencies. In background mode, `signature_NAME_noise_floor` and
  `signature_NAME_background_tail_probability` record the native-unit floor
  and unadjusted competitive tail probability. Each `signature_NAME_passes_cutoff` column records
  qualification before clone exclusions and budgeting. Frequencies here use retained paired cells of that lineage
  in that sample. Clone-level `max_frequency` uses all retained requested
  lineages in the sample; these coincide when only one lineage is selected.
- `signature_background.csv` (background mode): one calibration record per
  signature/patient/sample/lineage/clone-size, including background cell count,
  floor, quantile, actual draw count, exact versus Monte Carlo method, and status.
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

## Changes in 3.24.0

The default stopping rule now uses a size-matched empirical background floor
instead of score > 0. The total cap remains 200 and the rank gate remains the
top decile. Use `--signature-cutoff legacy` for the previous score rule. This
can substantially reduce shortlists, especially with very small samples.

## Changes in 3.23.1

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
