# Public Infectious-Disease Benchmark

Reproducible evaluation of `tcrsift` across three open-access human infectious-disease single-cell RNA + TCR cohorts spanning Influenza A, SARS-CoV-2, CMV, and EBV:

1. **10x Genomics Multiplexed 50-Dextramer Peripheral Blood Cohort** (Donors 1 and 2; `102,009` GEX-QC cells / `63,201` retained paired-TCR CD8 T cells across `22,873` clonotypes, stained with 44 viral and self pMHC dextramers, 6 negative-control multimers, and TotalSeq-C CD45RA/CD45RO surface antibodies).
2. **Minervina, Pogorelyy, Kirk et al., *Nature Immunology* 2022** ([Zenodo 6232103](https://zenodo.org/records/6232103); `21,573` GEX-QC cells / `13,614` retained paired-TCR CD8 T cells across `7,317` paired α/β clonotypes, including `4,212` dextramer-verified SARS-CoV-2 Spike and non-Spike clones spanning 19 epitopes across pre-exposure, convalescent post-infection, post-dose-1, acute post-dose-2 recall, and breakthrough infection).
3. **Fischer et al., *Nature Communications* 2021** ([GSE171037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171037); 4 paired unstimulated vs 16 h SARS-CoV-2 Spike peptide-stimulated PBMC libraries from hospitalized severe COVID-19 patients `GT_2` and `GT_6`, totaling `29,184` loaded PBMCs / `27,113` CD3-gated T cells across `3,474` β-clonotypes, with Cell Ranger 3 `.h5` GEX, `.csv.gz` contigs, and `all_contig_annotations.json.gz`).

## Running the benchmark

```bash
tcrsift data download --db vdjdb
python benchmarks/public_infectious/download.py --data-dir bench_infectious_data
python benchmarks/public_infectious/run.py \
    --data-dir bench_infectious_data \
    --work-dir bench_infectious_data/runs
python docs/figures/make_infectious_benchmark_figures.py
```

`run.py` writes compact summary tables to `benchmarks/public_infectious/results/`:

- `blood_signature_auroc.csv` — AUROC (with 2,000-draw clone bootstrap 95% intervals) for `tcrsift prioritize --context blood` signatures separating dextramer-verified viral memory clones from naive self-antigen (MART-1, encoded by *MLANA*) precursors and from other expanded blood clones, plus Influenza M1 vs EBV memory subtype contrasts and acute vaccine-recall ranking.
- `exposure_state_profiles.csv` — mean per-cell signature trajectories across SARS-CoV-2 exposure states (`naive`, `post_inf`, `inf_after_1`, `inf_after_2`, and `break`) in Minervina et al. 2022.
- `fischer_stimulation.csv` — paired unstimulated vs 16 h Spike-stimulated clone response (`score_stimulation_response`) across persistent CD3-gated β-clonotypes in Fischer et al. 2021, comparing abundance fold-change against within-clone activation and cytokine-transcript deltas (`AIM`, `AIMBroad`, `AcuteActivation`, and `IFNG`/`TNF`/`IL2`/`CCL4`).
- `vdjdb_matching.csv` — precision, recall, and off-target viral cross-hit rates for `match_clonotypes` (`strict_ab`, `ab_with_partial`, `b_only`, and Levenshtein-1 `b_only`) under both full VDJdb (circular for ingested cohorts) and leave-study-out VDJdb.
- `multimer_calling.csv` — negative-control-calibrated pMHC multimer calling (`tcrsift.multimer`) vs 10x `binarized_matrix.csv`, quantifying clonal concordance and CLR-normalized CITE-seq CD45RO/CD45RA profiles across viral multimers and promiscuous `A0301_KLGGALQAK_IE-1_CMV`.
- `summary.json` — run metadata, cell-filter counts, MART-1 singleton and CITE-seq summary statistics, 5'-truncated contig fractions (showing 98.7%–99.7% 5'-partial contigs and 0 paired full-length constructs in Fischer 2021), and VDJdb off-target species counts.
