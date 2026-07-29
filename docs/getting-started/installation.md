# Installation

## Requirements

- Python 3.9 or higher
- Enough memory and temporary disk for the combined single-cell dataset

The normal install includes the scientific stack, PDF sequence reports,
constant-region references, Excel input support, and the k-mer/TCRpeg sequence
probability backends.

## Install from PyPI

```bash
pip install tcrsift
```

## Install from Source

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
pip install -e .
```

## Optional Dependencies

Install everything, including atlas and optional plotting/report helpers:

```bash
pip install "tcrsift[full]"
```

Available extras:

| Extra | Adds |
| --- | --- |
| `atlas` | Harmony, Leiden, and graph dependencies for `tcrsift cells embed` |
| `plots` | UpSet plots (a built-in bar-chart fallback is used without it) |
| `reports` | `pdfkit` for legacy HTML-to-PDF report paths |
| `docs` | MkDocs and the documentation theme/plugins |
| `dev` | pytest, coverage, Ruff, and Pylint |
| `full` | `atlas`, `plots`, and `reports` together |

There are no `assembly`, `excel`, or `tcrpeg` extras. Those capabilities and
dependencies are part of the core package. Constant-region assembly uses
packaged canonical references; it does not require a separate
`pyensembl install` download.

## Verify Installation

```bash
tcrsift --version
```

Or in Python:

```python
import tcrsift
print(tcrsift.__version__)
```
