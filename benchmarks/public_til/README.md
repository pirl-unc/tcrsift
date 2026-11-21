# Public TIL Prioritization Benchmark

Evaluates `tcrsift prioritize --context solid-tumor` against experimentally
tested clones across five public datasets: Meng 2023 (GSE254250), Caushi 2021
(GSE176021), Eberhardt 2021 (GSE180268), Chiffelle 2024 (GSE222448), and
Pulliam 2024 (GSE227708). Methods, results, and limitations are documented in
`docs/benchmarks.md`.

```bash
pip install openpyxl                                             # Reads supplementary Excel tables
python benchmarks/public_til/download.py --data-dir bench_data   # ~3 GB download
python benchmarks/public_til/run.py --data-dir bench_data --work-dir bench_data/runs
python docs/figures/make_benchmark_figures.py
```

`results/` stores the summary CSV and JSON files plotted by
`make_benchmark_figures.py`. Downloaded raw matrices and intermediate run
directories under `bench_data/` are not committed to git.
