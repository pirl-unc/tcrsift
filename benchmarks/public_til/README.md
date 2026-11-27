# Public Solid-Tumor TIL Prioritization Benchmark

Reproducible evaluation of `tcrsift prioritize --context solid-tumor` across five open-access human tumor-infiltrating lymphocyte (TIL) scRNA + scTCR cohorts spanning 33 patients:

1. **Meng et al., *Science Translational Medicine* 2023** ([GSE254250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254250); 9 pancreatic ductal adenocarcinoma patients, `14,396` CD8+ TILs, `94` autologous-tumor-reactive vs `55` experimentally non-reactive clonotypes).
2. **Chiffelle et al., *Immunity* 2024 / Barras et al., *Science Immunology* 2024** ([GSE222448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222448); 11 metastatic melanoma patients, 13 baseline `T0` biopsies, `72` autologous-tumor-reactive vs `94` experimentally non-reactive clonotypes).
3. **Caushi et al., *Nature* 2021** ([GSE176021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176021); 8 NSCLC patients post-anti-PD-1, 57 tumor/metastasis/LN/normal lung libraries, MANAFEST-validated neoantigen `MANA` vs `viral` β-chains).
4. **Eberhardt et al., *Nature* 2021** ([GSE180268](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180268); 3 HPV+ HNSCC patients, 6 bulk PD-1+ CD8+ primary tumor and metastatic lymph node libraries paired with HLA-A\*01:01 HPV16 E2 pMHC tetramer sorts).
5. **Pulliam et al., *Cell Reports Medicine* 2024** ([GSE227708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227708); 2 Merkel cell carcinoma patients, hashtag-demultiplexed tumor biopsies with MCPyV pMHC multimer calls).

See [`PROVENANCE.md`](PROVENANCE.md) for complete ground-truth assay descriptions, label-mapping rules, and circularity disclosures, and [`docs/benchmarks.md`](../../docs/benchmarks.md) for full results and figures.

## Reproducing

```bash
tcrsift data download --db vdjdb
python benchmarks/public_til/download.py --data-dir bench_til_data
python benchmarks/public_til/run.py --data-dir bench_til_data --work-dir bench_til_data/runs
python scripts/figures/make_benchmark_figures.py
```

`run.py` writes compact CSV tables and `summary.json` (including `tcrsift_version` and `ground_truth_provenance`) to `benchmarks/public_til/results/`.
