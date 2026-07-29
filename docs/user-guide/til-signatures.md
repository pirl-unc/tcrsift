# Multi-sample TIL prioritization

Use this workflow when you have paired 10x VDJ + gene-expression data from
multiple TIL samples and want a reviewable list of expanded,
antigen-experienced clones.

The result is a **prioritized hypothesis list**, not a call of antigen
specificity. Expansion and T-cell state are indirect evidence; public-database
matches are incomplete; final candidates still need peptide/tumor recognition,
functional-avidity, and off-target testing.

## Choose the path

| Design | Path | What it does |
| --- | --- | --- |
| Two or more ordered/longitudinal TIL samples in the legacy per-timepoint layout | [`tcrsift til-select`](cli.md#tcrsift-til-select) | Harmonizes clones, scores compact TIL programs, and adds frequency-change branches. |
| Multiple standard CellRanger VDJ + GEX sample directories | [`examples/multi_sample_til.py`](https://github.com/pirl-unc/tcrsift/blob/main/examples/multi_sample_til.py) | Loads a sample sheet, scores several published/curated signatures within each sample, annotates publicness, and writes an auditable shortlist. |

Do not interpret an arbitrary ordering of independent samples as an expansion
trajectory. For independent tumors, rank by within-sample frequency and
signature scores; use increase flags only for true longitudinal data.

## End-to-end example

Create a normal TCRsift sample sheet. All samples may have `source: til`; this
Python workflow does not remove them as the culture-oriented `tcrsift run`
command does. Populate `patient_id` for every sample in a multi-patient
analysis: samples from one patient are combined, while the same public CDR3
pair remains a separate candidate in different patients.

```yaml title="til_samples.yaml"
samples:
  - sample: patient1_pre
    patient_id: patient1
    timepoint: pre
    source: til
    vdj_dir: /data/patient1_pre/vdj
    gex_dir: /data/patient1_pre/gex

  - sample: patient1_on_treatment
    patient_id: patient1
    timepoint: on_treatment
    source: til
    vdj_dir: /data/patient1_on_treatment/vdj
    gex_dir: /data/patient1_on_treatment/gex

  - sample: patient2_pre
    patient_id: patient2
    timepoint: pre
    source: til
    vdj_dir: /data/patient2_pre/vdj
    gex_dir: /data/patient2_pre/gex
```

From a source checkout (or after downloading the linked script), run the
example with one or more curated TCR databases:

```bash
python examples/multi_sample_til.py til_samples.yaml \
  -o til_candidates/ \
  --vdjdb /references/vdjdb.txt \
  --iedb /references/iedb.tsv \
  --cedar /references/cedar.tsv \
  --min-cells 2 \
  --min-frequency 0.001 \
  --signature-quantile 0.90
```

The three outputs serve different purposes:

- `all_scored_clones.csv`: every clone, score, risk flag, and exclusion reason.
- `candidate_clones.csv`: clones passing abundance, signature, and configured
  exclusion rules.
- `clone_sample_scores.csv`: per-(clone, sample) frequency and expression
  scores, so a pooled maximum can always be traced to its source sample.

Scores are computed separately within each sample after log1p(CP10K)
normalization. This prevents a high-depth or high-baseline sample from winning
simply because it has larger expression values. The example stops with a clear
error if a required signature gene is absent rather than silently scoring a
different, partial signature.

## Evidence used

The example requires observed abundance plus at least one signature in the
configured top quantile. It includes:

| Score | Intended evidence | Important limit |
| --- | --- | --- |
| `TumorReactive` | CD39/CXCL13, exhaustion, and residency program in fresh TIL | Curated composite, not direct antigen binding. |
| `Cytolytic` | PRF1/GZMB effector function | Also marks cytotoxic bystanders. |
| `Differentiated` | Effector-minus-naïve state | Expansion-associated, not tumor-specific. |
| `MANAscore` | CXCL13/ENTPD1 up and IL7R down | TCRsift implements a transparent signed-z proxy; the published trained ensemble is not shipped. |
| `NeoTCR8` / `NeoTCR4` | Published CD8/CD4 neoantigen-reactive gene sets | Gene-set enrichment, not a calibrated probability in a new cohort. |

The published registries and their exact scoring methods are documented in
[Signatures](../api/signatures.md). NeoTCR_PBL is available in the API but is
deliberately absent from this TIL example because it was derived for
circulating blood cells.

## Viral, MART-1, and public-receptor filters

The example separates facts from heuristics:

- `known_viral_match` is a public-database annotation. Known viral matches are
  excluded by default. A missing match means **unknown**, not non-viral.
- `known_mart1_match` is a database match to MART-1/Melan-A/MLANA or the common
  EAAGIGILTV/AAGIGILTV or altered ELAGIGILTV epitope forms. These matches are
  excluded by default in this example.
- `uses_trav12_2` is only a germline-bias flag. TRAV12-2 is strongly enriched
  among HLA-A*02:01/MART-1 receptors, but it is also used by other
  specificities. It is not excluded unless `--exclude-trav12-2` is explicitly
  supplied.
- `publicness_percentile` ranks the more public α/β chain against the shipped
  observed-repertoire background. Optionally remove the cohort's most-public
  tail with `--exclude-public-quantile 0.90`.
- `alpha_promiscuous` flags an α CDR3 paired with at least three distinct β
  CDR3s in the observed cohort. It is a review flag, not a default exclusion.

For an aggressive HLA-A*02:01 melanoma screen:

```bash
python examples/multi_sample_til.py til_samples.yaml \
  -o til_candidates_strict/ \
  --vdjdb /references/vdjdb.txt \
  --iedb /references/iedb.tsv \
  --exclude-trav12-2 \
  --exclude-public-quantile 0.90
```

Review `risk_flags` and `excluded_reason` before accepting this stricter
shortlist. TRAV12-2 bias is well established for MART-1, but the same V gene is
also biased in other HLA-A2 responses; removing the whole family trades
sensitivity for caution.

## Biological interpretation

Several TIL studies motivate combining expression state with clonotype
expansion, but none makes the combination a specificity assay. The
[Lowery et al. NeoTCR study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8996692/)
derived separate CD4 and CD8 expression signatures. The
[MANAscore study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11791090/) used a
minimal CXCL13/ENTPD1/IL7R program. Conversely,
[Simoni et al.](https://pubmed.ncbi.nlm.nih.gov/29769722/) showed that
virus-specific bystander T cells can be abundant in tumors and overlap
phenotypically with tumor-specific T cells.

For MART-1 specifically, TRAV12-2 bias has a structural and repertoire basis,
but it is not exclusive to MART-1
([Madura et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2785656/)).
That is why the example exposes the V-gene flag separately from exact
database matches and publicness scores.
