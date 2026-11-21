# Quick Start

This walkthrough runs `tcrsift run` on an antigen-stimulated T-cell culture
with a matched TIL sample. If you are starting from fresh tumor, blood,
effusion, or marrow single-cell VDJ + GEX samples without in-vitro culture, go
to [Clone Prioritization](../user-guide/til-signatures.md) (`tcrsift
prioritize`) instead.

## 1. Create a sample sheet

List your Cell Ranger output directories and sample metadata in YAML (or CSV):

```yaml title="samples.yaml"
samples:
  - sample: "Patient1_CMV"
    patient_id: "Patient1"
    vdj_dir: "/data/patient1/vdj"
    gex_dir: "/data/patient1/gex"
    antigen_type: "short_peptide"
    antigen_name: "CMV pp65 495-503"
    epitope_sequence: "NLVPMVATV"
    mhc_allele: "HLA-A*02:01"
    source: "culture"

  - sample: "Patient1_TIL"
    patient_id: "Patient1"
    vdj_dir: "/data/patient1_til/vdj"
    source: "til"
```

## 2. Run the pipeline

### Command line

```bash
# Optional one-time step: download VDJdb, IEDB, CEDAR, and IMGT V-region references
tcrsift data download

tcrsift run \
    --sample-sheet samples.yaml \
    --output-dir results/
```

`tcrsift run` executes seven stages:

1. Load Cell Ranger VDJ and GEX data (`loaded.h5ad`)
2. Classify CD4+ and CD8+ T cells (`phenotyped.h5ad`)
3. Collapse cells into paired α/β clonotypes and score `Pgen`/`Ppost` (`clonotypes.csv`)
4. Assign abundance tiers (`filtered_tier1.csv`–`filtered_tier5.csv`)
5. Flag public-database matches in VDJdb/IEDB/CEDAR (`annotated.csv`)
6. Cross-reference culture clones against `source: til` samples (`til_matched.csv`)
7. Assemble full-length constructs (`full_sequences.csv`) and summary plots/PDFs

### Python API

```python
import tcrsift

sample_sheet = tcrsift.load_sample_sheet("samples.yaml")
adata = tcrsift.load_samples(sample_sheet)
adata = tcrsift.phenotype_cells(adata, cd4_cd8_ratio=3.0)
clonotypes = tcrsift.aggregate_clonotypes(adata, group_by="CDR3ab")

filtered = tcrsift.filter_clonotypes(
    clonotypes,
    method="threshold",
    tcell_type="cd8",
)

annotated = tcrsift.annotate_clonotypes(
    filtered,
    vdjdb_path="/path/to/vdjdb",
    exclude_viral=True,
)
annotated.to_csv("results/annotated_clonotypes.csv", index=False)
```

## 3. Inspect outputs

```
results/
├── data/
│   ├── loaded.h5ad
│   ├── phenotyped.h5ad
│   ├── clonotypes.csv
│   ├── filtered_tier1.csv    # Strictest abundance tier (≥10 cells, ≥1% freq)
│   ├── filtered_tier2.csv
│   ├── filtered_tier3.csv
│   ├── filtered_tier4.csv
│   ├── filtered_tier5.csv
│   ├── annotated.csv         # Written when reference databases are available
│   ├── til_matched.csv       # Written when TIL samples are supplied
│   └── full_sequences.csv    # Assembled α, β, and β-2A-α sequences
├── plots/
│   ├── *.png
│   └── all-figures.pdf
├── config.yaml               # Resolved configuration for the run
├── provenance.json           # Package and dependency versions
└── provenance.txt
```

Core clonotype columns:

| Column | Description |
| --- | --- |
| `CDR3ab` | Paired identifier (`CDR3_alpha` + `_` + `CDR3_beta`) |
| `CDR3_alpha`, `CDR3_beta` | Chain CDR3 amino-acid sequences |
| `cell_count` | Total cells carrying this clonotype |
| `max_frequency` | Highest within-sample frequency |
| `tier` | Strictest abundance tier cleared (`tier1` = highest) |
| `Tcell_type_consensus` | Majority cell phenotype (`Confident CD8+`, etc.) |
| `db_match`, `is_viral` | Public TCR database hit and viral-epitope flag |

## Next steps

- [Sample Sheets](sample-sheets.md): all donor, enrichment-method, antigen, and
  path fields
- [Culture Pipeline](../user-guide/pipeline.md): stage-by-stage options
- [Filtering Strategies](../user-guide/filtering.md): custom tiers, FDR scope,
  and viral exclusion
- [Clone Prioritization](../user-guide/til-signatures.md): expression-state
  shortlisting with `tcrsift prioritize`
