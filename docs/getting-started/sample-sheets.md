# Sample Sheets

Sample sheets describe input files and experimental metadata for `tcrsift run`,
`tcrsift load`, and `tcrsift prioritize`. Both YAML and CSV formats are
supported.

## Example sample sheets

### YAML format

```yaml title="samples.yaml"
# Optional cohort-level flag: pool multi-donor FDR rankings only when donors
# share the same target antigen and MHC allele (default: false -> per-donor).
donors_share_antigen: false

samples:
  # Short-peptide culture (antigen == minimal epitope)
  - sample: "P1_CMV_AIM"
    patient_id: "P1"
    enrichment_method: "AIMpos"
    timepoint: "D14"
    vdj_dir: "/data/p1_cmv/vdj"
    gex_dir: "/data/p1_cmv/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"
    mhc_allele: "HLA-A*02:01"
    culture_days: 14
    source: "culture"
    tcell_type_expected: "CD8"

  # Long-peptide culture (minimal epitope depends on APC processing)
  - sample: "P1_KRAS"
    patient_id: "P1"
    vdj_dir: "/data/p1_kras/vdj"
    gex_dir: "/data/p1_kras/gex"
    antigen_type: "long_peptide"
    antigen_name: "KRAS G12D 25-mer"
    antigen_sequence: "TEYKLVVVGADGVGKSALTIQLIQ"
    culture_days: 21
    source: "culture"

  # Multi-peptide pool stimulation
  - sample: "P1_Pool"
    patient_id: "P1"
    vdj_dir: "/data/p1_pool/vdj"
    gex_dir: "/data/p1_pool/gex"
    antigen_type: "peptide_pool"
    antigen_names:
      - "KRAS_G12D"
      - "TP53_R175H"
      - "PIK3CA_H1047R"
    antigen_sequences:
      - "GADGVGKSAL"
      - "HMTEVVRHC"
      - "ARHGGWTTKM"
    culture_days: 14
    source: "culture"

  # Matched TIL sample (used automatically for TIL overlap in `tcrsift run`)
  - sample: "P1_TIL"
    patient_id: "P1"
    vdj_dir: "/data/p1_til/vdj"
    gex_dir: "/data/p1_til/gex"
    source: "til"
    tissue: "tumor"
```

### CSV format

```csv title="samples.csv"
sample,patient_id,vdj_dir,gex_dir,antigen_type,antigen_name,source
P1_CMV_AIM,P1,/data/p1_cmv/vdj,/data/p1_cmv/gex,short_peptide,CMV pp65 495-503,culture
P1_KRAS,P1,/data/p1_kras/vdj,/data/p1_kras/gex,long_peptide,KRAS G12D 25-mer,culture
P1_TIL,P1,/data/p1_til/vdj,/data/p1_til/gex,,,til
```

## Field reference

### Required fields

| Field | Requirement |
| --- | --- |
| `sample` | Unique sample identifier (numeric YAML values are coerced to text; duplicates raise an error). |
| Data path | At least one of `vdj_dir`, `gex_dir`, `sct_path`, `til_csv`, or `til_h5ad`. (`tcrsift prioritize` requires both `vdj_dir` and `gex_dir`.) |

### Input path fields

| Field | Description |
| --- | --- |
| `vdj_dir` | Cell Ranger VDJ output directory (`filtered_contig_annotations.csv` or `all_contig_annotations.csv`). |
| `gex_dir` | Cell Ranger GEX output directory (`filtered_feature_bc_matrix/` or `filtered_feature_bc_matrix.h5`). |
| `gex_path` | Standalone 10x `.h5` expression matrix when a full `gex_dir` is not used. |
| `sct_path`, `sct_sheet` | SCT platform Excel workbook path and sheet name (default sheet `"Cell"`). |
| `til_csv`, `til_h5ad` | Pre-processed TIL CSV or `.h5ad` file (alternatives to `vdj_dir` for `source: til`). |

### Study-design and cohort fields

| Field | Description |
| --- | --- |
| `patient_id` | Donor or patient identifier. Enables per-donor tiering, per-patient clone deduplication in `prioritize`, and multi-donor sharing buckets. |
| `source` | Pipeline role: `culture` (default; blank CSV cells also resolve to `culture`), `til`, `tetramer`, or `sct`. Ignored by `prioritize`, which uses `--context` instead. |
| `tissue` | Free-text anatomical site (such as `tumor`, `blood`, `LN`, `normal`). |
| `enrichment_method` | Free-text sorting or stimulation label (such as `AIMpos`, `tetpos`, `CTYneg`). Populating this in `tcrsift run` enables `filtered_by_method/` rankings, per-donor method overlap matrices, and method recovery plots. |
| `timepoint` | Free-text longitudinal label (such as `D7`, `D14`, `pre`, `post`). |
| `apc_type` | Antigen-presenting cell label for multi-APC designs (such as `mDC`, `B-LCL`). |
| `culture_days` | Integer duration of in-vitro culture in days. |
| `experiment` | Source label used when merging runs via `tcrsift unify`. |

### Antigen and epitope fields

- **Antigen** refers to the reagent delivered to APCs (short peptide, long
  peptide, minigene, mRNA, or whole protein).
- **Epitope** refers to the minimal MHC-bound peptide (typically 8–11 aa for
  MHC-I, 13–25 aa for MHC-II). Only minimal epitopes take an `mhc_allele`.

| Field | Description |
| --- | --- |
| `antigen_type` | Stimulus format (table below) |
| `antigen_name`, `antigen_sequence` | Single source antigen name and amino-acid sequence |
| `epitope_sequence`, `mhc_allele` | Known minimal epitope sequence and MHC allele (e.g. `HLA-A*02:01`) |
| `antigen_names`, `antigen_sequences`, `epitope_sequences` | Parallel lists for peptide pools or minigene libraries (`antigen_names` required when >1 sequence is given) |
| `antigen_description` | Optional free-text notes |

### Expected T-cell lineage fields

TCRsift infers the expected T-cell lineage for a culture sample using the first
available signal in this order:

1. `tcell_type_expected`: explicit `CD4`, `CD8`, or `mixed`
2. `pre_sorted`: FACS pre-sort (`CD4` or `CD8`)
3. `mhc_blocking`: `MHC-I` block implies `CD4`; `MHC-II` block implies `CD8`
4. `antigen_type`: biological default from the table below

| `antigen_type` | Default expected lineage | Typical input |
| --- | --- | --- |
| `short_peptide` | `CD8` | Minimal 8–11 aa MHC-I peptides |
| `tetramer_mhc1` | `CD8` | MHC-I multimer sort |
| `sct` | `CD8` | Single-chain pMHC-I trimer sort |
| `tetramer_mhc2` | `CD4` | MHC-II multimer sort |
| `long_peptide` | `mixed` | 15–25+ aa peptides requiring processing |
| `peptide_pool` | `mixed` | Multi-peptide stimulation pool |
| `minigene`, `minigene_library` | `mixed` | Single or pooled minigene constructs |
| `whole_protein`, `mrna` | `mixed` | Full-length protein or mRNA |

## Programmatic validation

`load_sample_sheet` enforces required fields, unique sample names, and valid
enums on load; `validate_sample_sheet` checks directory existence and flags
conflicting metadata (such as `short_peptide` paired with `CD4` expectation):

```python
from tcrsift import load_sample_sheet, validate_sample_sheet

sheet = load_sample_sheet("samples.yaml")
for warning in validate_sample_sheet(sheet):
    print(f"Warning: {warning}")
```
