# Sample Sheets

TCRsift uses sample sheets to define your input data and metadata.

## Supported Formats

Both YAML and CSV formats are supported.

### YAML Format

```yaml title="samples.yaml"
samples:
  - sample: "Patient1_CMV"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_description: "CMV pp65 495-503"
    culture_days: 14
    source: "culture"
    tcell_type_expected: "CD8"

  - sample: "Patient1_KRAS"
    vdj_dir: "/data/patient1_kras/vdj"
    gex_dir: "/data/patient1_kras/gex"
    antigen_type: "long_peptide"
    antigen_description: "KRAS G12D 25-mer"
    culture_days: 21
    source: "culture"

  - sample: "Patient1_TIL"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
    tissue: "tumor"
```

### CSV Format

```csv title="samples.csv"
sample,vdj_dir,gex_dir,antigen_type,antigen_description,source
Patient1_CMV,/data/patient1/vdj,/data/patient1/gex,short_peptide,CMV pp65,culture
Patient1_KRAS,/data/patient1_kras/vdj,/data/patient1_kras/gex,long_peptide,KRAS G12D,culture
Patient1_TIL,/data/patient1_til/vdj,,,til
```

## Required Fields

| Field | Description |
|-------|-------------|
| `sample` | Unique sample identifier |
| `vdj_dir` or `gex_dir` | At least one data directory is required |

## Optional Fields

### Data Paths

| Field | Description |
|-------|-------------|
| `vdj_dir` | Path to CellRanger VDJ output directory |
| `gex_dir` | Path to CellRanger GEX output directory |

### Antigen Information

| Field | Description |
|-------|-------------|
| `antigen_type` | Type of antigen (see below) |
| `antigen_description` | Free-text description |
| `epitope_sequence` | Peptide sequence if known |
| `mhc_allele` | MHC restriction if known |

### Culture Conditions

| Field | Description |
|-------|-------------|
| `culture_days` | Duration of culture in days |
| `source` | Sample source: `culture`, `til`, `tetramer`, `sct` |

### T Cell Type Expectations

| Field | Description |
|-------|-------------|
| `tcell_type_expected` | Expected type: `CD4`, `CD8`, `mixed` |
| `pre_sorted` | If cells were pre-sorted: `CD4`, `CD8` |
| `mhc_blocking` | If MHC was blocked: `MHC-I`, `MHC-II` |

## Antigen Types

TCRsift uses antigen type to infer biology-aware defaults:

| Type | Expected T Cells | Description |
|------|-----------------|-------------|
| `short_peptide` | CD8 | 8-11aa peptides bind MHC-I directly |
| `long_peptide` | mixed (favors CD4) | 15-25+aa requires processing |
| `minigene` | mixed | Minigene expression constructs |
| `whole_protein` | mixed | Full protein antigens |
| `tetramer_mhc1` | CD8 | MHC-I tetramer-sorted cells |
| `tetramer_mhc2` | CD4 | MHC-II tetramer-sorted cells |
| `sct_mhc1` | CD8 | Single-cell tetramer, MHC-I |
| `sct_mhc2` | CD4 | Single-cell tetramer, MHC-II |

## Sample Sources

| Source | Description |
|--------|-------------|
| `culture` | Antigen-stimulated culture (default) |
| `til` | Tumor-infiltrating lymphocytes |
| `tetramer` | Tetramer-sorted cells |
| `sct` | Single-cell tetramer |

## T Cell Type Inference

TCRsift automatically infers expected T cell type based on available metadata:

1. **Direct specification** (`tcell_type_expected`) takes priority
2. **Pre-sorting** (`pre_sorted`) is definitive
3. **MHC blocking** infers the opposite type:
    - `MHC-I` blocking → expect CD4
    - `MHC-II` blocking → expect CD8
4. **Antigen type** provides biological expectations

## Validation

TCRsift validates sample sheets and warns about:

- Duplicate sample names
- Missing data directories
- Conflicting metadata (e.g., short peptide expecting CD4)

```python
from tcrsift import load_sample_sheet, validate_sample_sheet

sample_sheet = load_sample_sheet("samples.yaml")
warnings = validate_sample_sheet(sample_sheet)

for warning in warnings:
    print(f"Warning: {warning}")
```
