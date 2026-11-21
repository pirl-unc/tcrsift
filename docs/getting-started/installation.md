# Installation

## Requirements

- Python 3.9–3.12
- Sufficient RAM and temporary disk space for your combined single-cell matrices

The standard `pip install tcrsift` package includes Scanpy/AnnData, PDF
sequence reports, packaged human TCR constant-region and leader references,
Excel import/export, cancer-testis-antigen panels (`oncoref`), and both the
k-mer and TCRpeg/PyTorch sequence-probability backends.

## Install from PyPI

```bash
pip install tcrsift
```

## Install from source

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
pip install -e .
```

## Optional extras

Install all optional clustering, plotting, and legacy report dependencies:

```bash
pip install "tcrsift[full]"
```

Individual extras:

| Extra | Packages added | Needed for |
| --- | --- | --- |
| `atlas` | `harmonypy`, `igraph`, `leidenalg` | `tcrsift cells embed` and `embed_cells` |
| `plots` | `upsetplot` | UpSet set-overlap plots (falls back to an intersection bar chart when absent) |
| `reports` | `pdfkit` | Legacy HTML-to-PDF conversion (`reportlab` PDF reports are already in core) |
| `docs` | `mkdocs`, `mkdocs-material`, `mkdocstrings[python]`, `pymdown-extensions` | Building the documentation site |
| `dev` | `pytest`, `pytest-cov`, `pytest-xdist`, `ruff`, `pylint` | Running tests and linters |
| `full` | `atlas` + `plots` + `reports` | All runtime extras |

Constant-region assembly uses packaged canonical FASTA references and does not
require downloading an Ensembl release via `pyensembl install`.

## Verify installation

```bash
tcrsift --version
```

```python
import tcrsift
print(tcrsift.__version__)
```
