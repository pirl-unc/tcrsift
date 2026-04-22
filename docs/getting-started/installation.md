# Installation

## Requirements

- Python 3.9 or higher
- numpy, pandas, scanpy, anndata (core dependencies)
- matplotlib, seaborn (visualization)

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

Install all optional feature bundles:

```bash
pip install "tcrsift[reports,assembly,excel]"
```

### PDF Report Generation

For generating PDF reports:

```bash
pip install "tcrsift[reports]"
```

You'll also need wkhtmltopdf:

=== "macOS"
    ```bash
    brew install wkhtmltopdf
    ```

=== "Ubuntu/Debian"
    ```bash
    sudo apt-get install wkhtmltopdf
    ```

=== "Windows"
    Download from [wkhtmltopdf.org](https://wkhtmltopdf.org/downloads.html)

### Ensembl Data for Constant Regions

For assembling full-length TCR sequences with constant regions:

```bash
pip install "tcrsift[assembly]"
pyensembl install --release 93 --species human
```

### Excel Support

For SCT Excel input files:

```bash
pip install "tcrsift[excel]"
```

## Verify Installation

```bash
tcrsift --version
```

Or in Python:

```python
import tcrsift
print(tcrsift.__version__)
```
