# Public Paired Blood-Tumor (`blood-tumor`), Post-ACT (`T30`), & CD4+ TIL (`--tcell-type cd4`) Benchmark

Reproducible evaluation of `tcrsift` across four open-access human cancer immunology cohorts testing `--context blood-tumor` (`NeoTCR_PBL`), longitudinal post-adoptive-cell-therapy (`T30`) TIL persistence, and human CD4+ TIL prioritization (`--tcell-type cd4`, `NeoTCR4`, `pre_sorted=CD4`, and `--exclude-treg`):

1. **Pulliam et al., *Cell Reports Medicine* 2024** ([GSE227708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227708), `PMID:38340723`; 10 demultiplexed peripheral blood libraries + 2 matched Merkel cell carcinoma tumor libraries with `>=20` cells across 5 hashed 10x runs / 4 patients, with multimer-verified `MCPyV` antitumor vs `CMV`/`EBV`/`Flu` viral control TCRs).
2. **Chiffelle / Barras et al., *Immunity* & *Science Immunology* 2024** ([GSE222448](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222448), `PMID:39321802`; 8 pre-ACT `T0` and 6 Day-30 post-ACT `T30` melanoma biopsies across Patients 1, 5, 7, 8, 10, and 13 with autologous tumor-reactive `TR` vs `NTR` labels).
3. **Veatch et al., *Cancer Cell* 2022** ([GSE198265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198265), `PMID:35413271`; human melanoma intratumoral CD4+ T cells across 4 patients with experimentally verified autologous neoantigen-specific `CD4_neo` and self-antigen `CD4_self` clonotypes).
4. **Oh et al., *Cell* 2020** ([GSE149652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149652), `PMID:32497499`; whole-transcriptome 10x scRNA-seq + paired α/β TCR across 7 human bladder cancer CD4+ TIL libraries and 6 matched normal bladder CD4+ libraries).

See [`PROVENANCE.md`](PROVENANCE.md) for ground-truth assay provenance, sample-to-patient mappings, and targeted-panel / FACS `pre_sorted` disclosures, and [`docs/benchmarks-blood-cd4.md`](../../docs/benchmarks-blood-cd4.md) for full results and figures.

## Reproducing

```bash
tcrsift data download --db vdjdb
python benchmarks/public_blood_cd4/download.py --data-dir bench_blood_cd4_data
python benchmarks/public_blood_cd4/run.py \
    --data-dir bench_blood_cd4_data \
    --work-dir bench_blood_cd4_data/runs
python docs/figures/make_blood_cd4_benchmark_figures.py
```

`run.py` writes compact CSV result tables and `summary.json` to `benchmarks/public_blood_cd4/results/`:

- `auroc_blood_cd4.csv` — AUROCs (with 2,000-draw clone bootstrap 95% CIs) across Pulliam peripheral blood (`MCPyV` vs `CMV`/`EBV`/`Flu`, `MCPyV` vs unlabeled blood, and tumor-infiltrating vs blood-only clones), Chiffelle `T0` vs `T30` (`TR` vs `NTR`), Veatch melanoma CD4+ TILs (`CD4_neo` vs unlabeled/non-Treg/self), and Oh whole-transcriptome bladder CD4+ TILs (`CD4CXCL13` + cytotoxic vs `CD4CM`/`TH17` and vs `CD4IL2RAHI` Treg).
- `shortlist_pulliam_blood.csv` — Shortlist yields on Pulliam peripheral blood comparing `--context blood-tumor`, `--context blood`, and `--exclude-known-viral` (`strict_ab` vs `ab_with_partial`).
- `shortlist_chiffelle_act.csv` — Shortlist yields and tested precision (`TR / (TR + NTR)`) on matched `T0` pre-ACT vs `T30` post-ACT melanoma tumors.
- `shortlist_veatch_cd4.csv` — Shortlist yields on Veatch melanoma CD4+ TILs with and without `FOXP3+` Treg exclusion.
- `shortlist_oh_cd4.csv` — Shortlist yields and subtype composition on Oh whole-transcriptome bladder CD4+ TILs comparing `pre_sorted=CD4` default, `--exclude-treg` (`0.50` and `0.25`), and `pre_sorted=None` (`CD4 mRNA > 0` only).
