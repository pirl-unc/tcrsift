# Public TIL prioritization benchmark

Scores `tcrsift prioritize --context solid-tumor` against experimentally
validated clones in two open datasets (Caushi 2021, GSE176021; Meng 2023,
GSE254250). Results and caveats are written up in `docs/benchmarks.md`.

```bash
pip install openpyxl                      # to read the Caushi label table
python benchmarks/public_til/download.py --data-dir bench_data   # about 1.4 GB
python benchmarks/public_til/run.py --data-dir bench_data --work-dir bench_data/runs
python docs/figures/make_benchmark_figures.py
```

`results/` holds the small tables the documentation figures are drawn from.
The downloaded data and the `tcrsift` output directories are not committed.
