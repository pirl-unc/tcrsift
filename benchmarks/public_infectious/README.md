# Public Infectious-Disease & Peripheral-Blood Benchmark

Reproducible evaluation of `tcrsift` across three open-access human infectious-disease and peripheral-blood single-cell RNA + TCR cohorts spanning Influenza A, SARS-CoV-2, CMV, EBV, and naive self-antigen (`MART-1`, encoded by *MLANA*) precursors:

1. **10x Genomics Multiplexed 50-Dextramer Peripheral Blood Cohort** (Healthy Donors 1 & 2; `102,009` GEX-QC blood CD8+ T cells / `63,201` paired-TCR retained cells across 44 viral/cancer pMHC dextramers, 6 negative-control multimers, and CITE-seq TotalSeq-C CD45RA/CD45RO surface antibodies).
2. **Minervina, Pogorelyy, Kirk et al., *Nature Immunology* 2022** ([Zenodo 6232103](https://zenodo.org/records/6232103); `13,614` paired-TCR CD8+ T cells / `7,317` paired αβ clonotypes including `4,212` curated dextramer-verified SARS-CoV-2 Spike and non-Spike clones across primary infection, mRNA vaccination, acute post-dose-2 recall, and breakthrough infection).
3. **Fischer, Ansari, Wagner et al., *Nature Communications* 2021** ([GSE171037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171037); 4 paired unstimulated control vs 16-hour SARS-CoV-2 Spike peptide-stimulated PBMC libraries from hospitalized severe COVID-19 patients `GT_2` and `GT_6`, with Cell Ranger 3 `.h5` GEX, `.csv.gz` contigs, and `all_contig_annotations.json.gz`).

See [`PROVENANCE.md`](PROVENANCE.md) for complete ground-truth assay descriptions, barcode/epitope label-mapping rules, VDJdb study-holdout construction, and circularity disclosures, and [`docs/benchmarks-infectious.md`](../../docs/benchmarks-infectious.md) for full results and figures.

## Reproducing

```bash
tcrsift data download --db vdjdb
python benchmarks/public_infectious/download.py --data-dir bench_infectious_data
python benchmarks/public_infectious/run.py \
    --data-dir bench_infectious_data \
    --work-dir bench_infectious_data/runs
python docs/figures/make_infectious_benchmark_figures.py
```

`run.py` writes compact summary tables and `summary.json` (including `tcrsift_version`, `vdjdb_release`, and `ground_truth_provenance`) to `benchmarks/public_infectious/results/`:

- `blood_signature_auroc.csv` — AUROC (with 2,000-draw clone bootstrap 95% intervals) for `tcrsift prioritize --context blood` signatures separating bona fide viral memory clones from naive `MART-1` precursors, separating expanded viral clones from unlabeled expanded blood clones, separating Influenza M1 vs EBV memory subtypes, and separating acute post-dose-2 SARS-CoV-2 recall from convalescent memory.
- `blood_cutoffs_and_shortlist.csv` — Default cutoff (`eligible_for_review`) and 200-clone round-robin shortlist (`selected_for_review`) counts from unmodified `tcrsift prioritize --context blood` across 10x Donors 1 & 2 and Minervina 2022, with and without `--exclude-known-viral` (`ab_with_partial` vs `strict_ab` on study-held-out VDJdb).
- `exposure_state_profiles.csv` — Mean signature trajectories across SARS-CoV-2 exposure stages (`naive`, `post_inf`, `inf_after_1`, `inf_after_2`, `break`).
- `vdjdb_matching.csv` — Exact-epitope and virus-species precision and recall across four matching rules (`strict_ab`, `ab_with_partial`, `b_only`, `levenshtein b_only`), comparing unmodified VDJdb (**circular**) against study-held-out VDJdb (**independent**).
- `multimer_calling.csv` — Per-multimer clonal concordance (`summarize_clone_multimer_specificity`), single-cell agreement with 10x's `binarized_matrix.csv`, sticky-reagent (`KLGGALQAK`) rejection, and CITE-seq `CLR(CD45RO) - CLR(CD45RA)` memory validation.
- `fischer_stimulation.csv` — Paired unstimulated vs 16 h Spike-stimulated persistent β-clonotypes (`score_stimulation_response`), contrasting Fisher exact count enrichment against within-clone cytokine/AIM transcriptional induction.
