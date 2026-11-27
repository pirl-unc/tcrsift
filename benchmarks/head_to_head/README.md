# Head-to-Head Comparator Benchmarks (`benchmarks/head_to_head`)

Reproducible head-to-head evaluation of `tcrsift` against open-source baselines across four stages:

1. **Ex vivo prioritization**: `tcrsift` vs `clone_size`, single-gene `CXCL13`, `scanpy.tl.score_genes` (`NeoTCR8` / `NeoTCR4` / `NeoTCR_PBL`), and `scanpy` composite top-$K$ across 5 human TIL and cancer-blood cohorts.
2. **Full-length TCR construct assembly**: `tcrsift assemble` vs `stitchr` / `thimble` (`v1.3.1`, Heather et al. *NAR* 2022) across 5'-complete human 10x contigs, 5'-truncated human deposits, and murine Cell Ranger VDJ outputs.
3. **TCR generation probability & publicness**: `tcrsift` bundled k-mer `Pgen`/`Ppost` vs `OLGA` (`v1.3.0`, Sethna et al. *Bioinformatics* 2019) across human (`n = 1,500`) and murine (`n = 1,781`) repertoires.
4. **Multiplexed pMHC multimer calling**: `tcrsift --call-multimers` clonal concordance vs 10x per-cell threshold, CLR `z >= 3`, and bimodal GMM callers on 10x Healthy Donors 1 & 2.

See [`PROVENANCE.md`](PROVENANCE.md) for methodological details and [`docs/benchmarks-overview.md`](../../docs/benchmarks-overview.md) for the summary tables and head-to-head figures.

## Reproducing

Run the four cohort benchmark suites (`public_til`, `public_infectious`, `public_blood_cd4`, `public_mouse`) first so their cached inputs and scored clone tables are available, then install the external comparators and run:

```bash
pip install stitchr==1.3.1 olga==1.3.0
stitchrdl -s HUMAN && stitchrdl -s MOUSE

python benchmarks/head_to_head/run.py
python scripts/figures/make_overview_benchmark_figures.py
```
