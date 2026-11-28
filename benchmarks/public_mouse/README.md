# Public Murine (`--species mouse`) Benchmark

Reproducible evaluation of `tcrsift` on murine (*Mus musculus*) single-cell RNA + paired α/β TCR sequencing across 19 open-access 10x Genomics libraries spanning syngeneic solid tumor + peripheral blood/lymphoid compartments, autoimmune target-organ vs splenic MHC-I tetramer sorts, and acute viral MHC-II tetramer+ CD4+ T cells:

1. **Multi-Compartment MC38 Colon Adenocarcinoma** ([GSE168944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168944); 8 paired 10x GEX + VDJ libraries across matched Tumor, Peripheral Blood, Draining Lymph Node, and Spleen in biological replicate C57BL/6 mice `Mouse01` & `Mouse02`).
2. **Kasmani et al., *Life Science Alliance* 2022** ([GSE200608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200608), `PMID:35667687`; 6 paired 10x GEX + VDJ libraries of NOD mouse pancreatic islet and splenic `H-2Kd IGRP206-214 / NRP-V7` [`KYNKANVFL`] tetramer+ vs islet tetramer- CD8+ T cells).
3. **Khatun et al., *J Exp Med* 2021** ([GSE158896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158896), `PMID:33201171`; 5 paired 10x GEX + VDJ libraries of Day-10 acute LCMV Armstrong `I-Ab GP66-77` [`DIYKGVYQFKSV`] tetramer+ splenic CD4+ T cells across C57BL/6 mice `M1`–`M5`).

See [`PROVENANCE.md`](PROVENANCE.md) for ground-truth assay details, accession mappings, and VDJdb circularity controls, and [`docs/benchmarks-mouse.md`](../../docs/benchmarks-mouse.md) for full results and figures.

## Reproducing

```bash
tcrsift data download --db vdjdb
tcrsift data download --group germline_mouse
python benchmarks/public_mouse/download.py --data-dir bench_mouse_data
python benchmarks/public_mouse/run.py \
    --data-dir bench_mouse_data \
    --work-dir bench_mouse_data/runs
python scripts/figures/make_mouse_benchmark_figures.py
```

`run.py` writes compact CSV result tables and `summary.json` to `benchmarks/public_mouse/results/`:

- `auroc_mouse.csv` — AUROCs (with 2,000-draw clone bootstrap 95% CIs) for murine MGI signatures on MC38 tumor-enriched vs lymphoid-shared clones (including MGI `TumorReactive` vs naive human `CXCL13` `TumorReactive`), MC38 peripheral-blood tracking of tumor-infiltrating clones, and Kasmani `IGRP` target-organ (`Islet_TetPos`) vs circulating (`Spleen_TetPos`) and within-islet (`Islet_TetNeg`) contrasts.
- `vdjdb_mouse_recovery.csv` — `MusMusculus` VDJdb recall, `IGRP`/`LCMV` precision, and out-of-strain MC38 false-positive rates across `strict_ab`, `ab_with_partial`, `b_only`, and `b_only + levenshtein-1`, comparing unmodified VDJdb (**circular** on Kasmani `PMID:35667687`) against study-held-out VDJdb (**independent**).
- `shortlist_mc38.csv` — Default shortlist yields on MC38 tumor (`--context solid-tumor`) and matched blood (`--context blood`).
- `mouse_assembly.csv` — Full-length murine single-chain TCR construct assembly (`Trac`, `Trbc1`, `Trbc2`) via `assemble_full_sequences(..., species="mouse", germline_vdj_fallback=True)`.
