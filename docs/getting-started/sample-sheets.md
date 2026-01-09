# Sample Sheets

TCRsift uses sample sheets to define your input data and metadata.

## Supported Formats

Both YAML and CSV formats are supported.

### YAML Format

```yaml title="samples.yaml"
samples:
  # Minimal peptide culture (antigen == epitope)
  - sample: "Patient1_CMV"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"  # same as antigen for minimal peptides
    mhc_allele: "HLA-A*02:01"
    culture_days: 14
    source: "culture"
    tcell_type_expected: "CD8"

  # Long peptide culture (epitope unknown until processing)
  - sample: "Patient1_KRAS"
    vdj_dir: "/data/patient1_kras/vdj"
    gex_dir: "/data/patient1_kras/gex"
    antigen_type: "long_peptide"
    antigen_name: "KRAS G12D 25-mer"
    antigen_sequence: "TEYKLVVVGADGVGKSALTIQLIQ"  # full long peptide
    # epitope_sequence unknown - depends on processing
    culture_days: 21
    source: "culture"

  # Peptide pool stimulation
  - sample: "Patient1_Pool"
    vdj_dir: "/data/patient1_pool/vdj"
    gex_dir: "/data/patient1_pool/gex"
    antigen_type: "peptide_pool"
    antigen_names:  # required when >1 antigen
      - "KRAS_G12D"
      - "TP53_R175H"
      - "PIK3CA_H1047R"
    antigen_sequences:  # optional
      - "GADGVGKSAL"
      - "HMTEVVRHC"
      - "ARHGGWTTKM"
    culture_days: 14
    source: "culture"

  # SCT selection (epitope is known from the construct)
  - sample: "Patient1_SCT"
    vdj_dir: "/data/patient1_sct/vdj"
    antigen_type: "sct"
    antigen_name: "PRAME"  # source protein
    epitope_sequence: "SLLQHLIGL"  # what's in the SCT
    mhc_allele: "HLA-A*02:01"
    source: "sct"

  # TIL sample (no antigen info needed)
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
| `antigen_name` | Name of source antigen (e.g., "PRAME", "CMV pp65") |
| `antigen_sequence` | Sequence of source antigen (may be long) |
| `epitope_sequence` | Minimal peptide AA sequence that binds MHC |
| `mhc_allele` | MHC restriction (e.g., "HLA-A*02:01") |
| `antigen_names` | List of source antigen names (for pools/libraries) |
| `antigen_sequences` | List of source antigen sequences (for pools/libraries) |
| `epitope_sequences` | List of minimal epitope sequences (for pools, if known) |

**Antigen vs Epitope:**

- **Antigen** = what you gave to APCs (whole protein, long peptide, minigene, mRNA, etc.)
- **Epitope** = the minimal peptide that sits in the MHC groove (8-11aa for MHC-I, 13-25aa for MHC-II)

For minimal peptide stimulation, antigen == epitope. For whole proteins, the antigen is much larger than the processed epitope. Only the epitope can have an MHC restriction.

**When to use each field:**

- **Single antigen**: Use `antigen_name` and optionally `antigen_sequence`. If you know the minimal epitope, add `epitope_sequence` and `mhc_allele`.

- **Tetramer/SCT**: The epitope is known. Provide `epitope_sequence`, `mhc_allele`, and `antigen_name` to describe the source.

- **Pools/libraries**: Use `antigen_names` (required if >1 antigen) and optionally `antigen_sequences`.

**Naming rules**: At least a name or sequence must be provided. If only a sequence is given, it becomes the name. If multiple sequences are given without names, that's ambiguous.

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
| `peptide_pool` | mixed | Pool of peptides for stimulation |
| `minigene` | mixed | Single minigene expression construct |
| `minigene_library` | mixed | Library of multiple minigene constructs |
| `whole_protein` | mixed | Full protein antigens |
| `mrna` | mixed | mRNA encoding one or more antigens |
| `tetramer_mhc1` | CD8 | MHC-I tetramer-sorted cells (single antigen) |
| `tetramer_mhc2` | CD4 | MHC-II tetramer-sorted cells (single antigen) |
| `sct` | CD8 | Single-chain trimer (pMHC-I: alpha-B2M-peptide fusion) |

## Sample Sources

| Source | Description |
|--------|-------------|
| `culture` | Antigen-stimulated culture (default) |
| `til` | Tumor-infiltrating lymphocytes |
| `tetramer` | Tetramer-sorted cells (MHC-I or MHC-II) |
| `sct` | Single-chain trimer-sorted cells (MHC-I only) |

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
