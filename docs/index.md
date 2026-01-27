# TCRsift

T-cell receptor selection for TCR-T studies from antigen-specific culture and scRNA/VDJ sequencing.

## Overview

TCRsift is a comprehensive pipeline for identifying antigen-specific T cell receptor clones from single-cell sequencing data. It supports:

- **Data Loading**: Parse CellRanger VDJ and GEX outputs
- **CD4/CD8 Phenotyping**: Classify cells based on gene expression
- **Clonotype Aggregation**: Group cells by CDR3 sequences
- **Tiered Filtering**: Apply biology-aware filtering to identify antigen-specific clones
- **Database Annotation**: Match against VDJdb, IEDB, and CEDAR
- **TIL Matching**: Find culture-validated TCRs in tumor samples
- **Sequence Assembly**: Build full-length TCR sequences with constant regions

## Quick Example

```python
import tcrsift

# Load sample sheet
sample_sheet = tcrsift.load_sample_sheet("samples.yaml")

# Load samples
adata = tcrsift.load_samples(sample_sheet)

# Phenotype cells
adata = tcrsift.phenotype_cells(adata)

# Aggregate clonotypes
clonotypes = tcrsift.aggregate_clonotypes(adata)

# Filter clonotypes
filtered = tcrsift.filter_clonotypes(clonotypes, method="threshold", tcell_type="cd8")

# Annotate with VDJdb
annotated = tcrsift.annotate_clonotypes(filtered, vdjdb_path="/path/to/vdjdb")
```

Or via command line:

```bash
tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/ \
    --vdjdb /path/to/vdjdb
```

## Installation

```bash
pip install tcrsift
```

See [Installation](getting-started/installation.md) for detailed instructions.

## Next Steps

- [Quick Start Guide](getting-started/quickstart.md)
- [Sample Sheet Format](getting-started/sample-sheets.md)
- [Pipeline Overview](user-guide/pipeline.md)
- [API Reference](api/sample_sheet.md)
