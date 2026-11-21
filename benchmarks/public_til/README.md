# Public TIL prioritization benchmark

Scores `tcrsift prioritize --context solid-tumor` against experimentally
labelled clones in five open datasets (Meng 2023, GSE254250; Caushi 2021,
GSE176021; Eberhardt 2021, GSE180268; Chiffelle 2024, GSE222448; Pulliam 2024,
GSE227708). Results and caveats are written up in `docs/benchmarks.md`.

```bash
pip install openpyxl                      # to read the label tables
python benchmarks/public_til/download.py --data-dir bench_data   # about 3 GB
python benchmarks/public_til/run.py --data-dir bench_data --work-dir bench_data/runs
python docs/figures/make_benchmark_figures.py
```

`results/` holds the small tables the documentation figures are drawn from.
The downloaded data and the `tcrsift` output directories are not committed.
