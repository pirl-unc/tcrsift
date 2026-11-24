# Head-to-Head Comparator Benchmarks (`benchmarks/head_to_head`)

Reproducible head-to-head evaluation of `tcrsift` against open-source baselines across four stages:

1. **Ex vivo prioritization**: `tcrsift` vs `clone_size`, single-gene `CXCL13`, `scanpy.tl.score_genes` (`NeoTCR8` / `NeoTCR4` / `NeoTCR_PBL`), and `scanpy` composite top-$K$ across 5 human TIL and cancer-blood cohorts.
2. **Full-length TCR construct assembly**: `tcrsift assemble` vs `stitchr` / `thimble` (`v1.3.1`, Heather et al. *NAR* 2022) across 5'-complete human 10x contigs, 5'-truncated human deposits, and murine Cell Ranger VDJ outputs.
3. **TCR generation probability & publicness**: `tcrsift` bundled k-mer `Pgen`/`Ppost` vs `OLGA` (`v1.3.0`, Sethna et al. *Bioinformatics* 2019) across human (`n = 1,500`) and murine (`n = 1,781`) repertoires.
4. **Multiplexed pMHC multimer calling**: `tcrsift --call-multimers` clonal concordance vs 10x per-cell threshold, CLR `z >= 3`, and bimodal GMM callers on 10x Healthy Donors 1 & 2.

See [`PROVENANCE.md`](PROVENANCE.md) for methodological details and [`docs/benchmarks-overview.md`](../../docs/benchmarks-overview.md) for the master Failure-Mode Scorecard and figures.

## Reproducing

```bash
.venv/bin/python -m venv --system-site-packages /tmp/h2h_venv
/tmp/h2h_venv/bin/pip install stitchr==1.3.1 olga==1.3.0
/tmp/h2h_venv/bin/stitchrdl -s HUMAN && /tmp/h2h_venv/bin/stitchrdl -s MOUSE

MPLCONFIGDIR=/tmp/mpl_cache /tmp/h2h_venv/bin/python benchmarks/head_to_head/run.py
MPLCONFIGDIR=/tmp/mpl_cache .venv/bin/python docs/figures/make_overview_benchmark_figures.py
```
