# Gene-Expression Signatures

TCRsift defines gene-expression signatures across two modules:

- `tcrsift.signatures` (`T_CELL_SIGNATURES`): compact unsigned marker panels
  for per-cell and per-clone state scoring, plus cell-type registries for
  `annotate_cells`.
- `tcrsift.signature_methods` (`SIGNATURES` and `NEOANTIGEN_SIGNATURES`):
  signed z-score panels and published neoantigen-reactivity gene sets used by
  `tcrsift prioritize` and `score_by_name`.

Signature scores quantify transcriptional state rather than epitope
specificity, and depend on tissue context: `AIM` and `AcuteActivation` reflect
recent stimulation in culture, whereas `TumorReactive` and exhaustion panels
reflect chronic antigen exposure in fresh TIL.

## Compact functional programs (`T_CELL_SIGNATURES`)

Human signatures follow HGNC nomenclature (`T_CELL_SIGNATURES`); mouse
counterparts (`MOUSE_T_CELL_SIGNATURES`, retrieved via
`t_cell_signatures_for_species("mouse")`) adapt markers without 1:1 murine
T-cell orthologs (dropping primate-restricted `GNLY` and replacing `CXCL13`
with `PDCD1, HAVCR2, ITGAE` in `tumor_reactive` and `ENTPD1` in
`expansion_core`; in `MOUSE_SIGNATURES`, `GNLY` is replaced by `GZMA` and
`CXCL13` by `CCL3, CCL4`).

| Name | Genes | Intended readout |
| --- | --- | --- |
| `effector` | `IFNG, GZMB, PRF1, GNLY, NKG7` | Cytotoxic effector differentiation |
| `activation` | Deprecated alias of `effector` | Retained for backward compatibility |
| `naive_stem` | `TCF7, LEF1, CCR7, SELL, IL7R, CD27, CD28` | Naive and stem-like memory state |
| `antigen_response` | `TNFRSF9, MKI67` | Recent activation and proliferation |
| `aim` | `TNFRSF9, TNFRSF4, IL2RA, MKI67` | Activation-induced-marker response in stimulated culture |
| `cytolytic` | `PRF1, GZMB` | Minimal perforin/granzyme B readout |
| `exhaustion` | `PDCD1, TOX, LAG3, HAVCR2, TIGIT` | Chronic stimulation / exhaustion in TIL |
| `tumor_reactive` | `CXCL13, ENTPD1` | Minimal two-gene tumor-reactive TIL readout |
| `expansion_core` | `MKI67, TNFRSF9, EGR2, IFNG, CXCL13, HAVCR2` | Cross-compartment expansion-associated markers |

`MARKER_PANEL_HGNC` is the combined gene extraction list used by `til-select`
rather than a standalone biological signature.

## Signatures scored by `tcrsift prioritize` (`SIGNATURES`)

`tcrsift prioritize` selects panels by these names (`--signatures`,
`--context`), scoring each cell as the mean of within-sample z-scored positive
genes minus negative genes:

| Name | Up genes | Down genes | Provenance |
| --- | --- | --- | --- |
| `Differentiated` | `IFNG, GZMB, PRF1, GNLY, NKG7` | `TCF7, LEF1, CCR7, SELL` | Canonical effector minus naive/stem markers, fixed on the two-donor MART-1 pilot ([reproducibility map](../reproducibility_map.md)) |
| `AntigenExperienced` | `IFNG, GZMB, PRF1, GNLY, NKG7` | — | Positive arm of `Differentiated` |
| `Cytolytic` | `PRF1, GZMB` | — | Perforin and granzyme B ([Caushi 2021](https://doi.org/10.1038/s41586-021-03752-4), [Krishna 2020](https://doi.org/10.1126/science.abb9847), [Hanada 2022](https://doi.org/10.1016/j.ccell.2022.03.012)) |
| `AcuteActivation` | `TNFRSF9, MKI67` | — | 4-1BB and Ki-67, separated from `AntigenExperienced` after the MART-1 pilot |
| `Proliferation` | `MKI67` | — | Ki-67 alone |
| `AIM` | `TNFRSF9, TNFRSF4` | — | 4-1BB and OX40 mRNAs corresponding to surface AIM sorting |
| `CirculatingMemory` | `GZMK, IL7R` | `CCR7` | Antigen-experienced circulating memory CD8/CD4 readout used by `blood`, `blood-tumor`, `mpe`, and `heme` presets |
| `TumorReactive` | `CXCL13, ENTPD1, PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4, ITGAE` | — | CD39 and CD103 ([Simoni 2018](https://doi.org/10.1038/s41586-018-0130-2), [Duhen 2018](https://doi.org/10.1038/s41467-018-05072-0)) plus `CXCL13` and six TIL exhaustion markers |

`Broad` variants are also registered in `SIGNATURES`.
Note that nine-gene `TumorReactive` above is broader than the two-gene
`tumor_reactive` entry in `T_CELL_SIGNATURES`. Empirical AUROCs across human
and murine cohorts are reported in the [Benchmark overview](../benchmarks-overview.md)
and cohort pages.

## Published neoantigen-reactivity registry (`NEOANTIGEN_SIGNATURES`)

`score_by_name` dispatches each signature using its registered scoring method:

| Name | Reference | TCRsift scoring rule |
| --- | --- | --- |
| `MANAscore` | Zeng et al., *Nat Commun* 2025 (PMID 39900903) | Unweighted three-gene signed-z proxy: `(z(CXCL13) + z(ENTPD1) - z(IL7R)) / sqrt(3)`. Does not run the paper's fitted random-forest/linear ensemble. |
| `NeoTCR8` | Lowery et al., *Science* 2022 (PMID 35113651) | Published 243-gene CD8 set via Scanpy `score_genes` |
| `NeoTCR4` | Lowery et al., *Science* 2022 (PMID 35113651) | Published 40-gene CD4 set via Scanpy `score_genes` |
| `NeoTCR_PBL` | Yossef et al., *Cancer Cell* 2023 (PMID 38039963) | Published 151-gene circulating CD8 set via Scanpy `score_genes` (`blood-tumor` preset) |

Scoring notes for `NeoTCR4`, `NeoTCR8`, and `NeoTCR_PBL`:

- **AnnData inputs:** pass an `AnnData` whose `.X` contains the full
  `log1p(CP10K)`-normalized transcriptome with `log1p=False` (`score_by_name`
  always reads `.X` rather than `.raw`, and does not apply a second logarithm).
  Scores are Scanpy `score_genes` module scores (signature mean minus
  expression-matched control-gene mean), **not** single-cell GSEA; numeric
  scGSEA cutoffs from the original papers do not transfer.
- **DataFrame inputs:** because a subset `DataFrame` lacks the full gene
  universe needed for expression-matched control bins, `score_by_name` warns
  and falls back to mean per-gene z-scores (applying `log1p` first by default;
  pass `log1p=False` if values are already log-normalized).

## Usage

Compact unsigned panels:

```python
from tcrsift import compute_signature_scores_per_cell

compute_signature_scores_per_cell(adata)
print([col for col in adata.obs if col.startswith("signature_")])
```

Signed and published signatures:

```python
import scanpy as sc
from tcrsift import score_by_name

sc.pp.normalize_total(adata, target_sum=10_000)
sc.pp.log1p(adata)

for name in ("MANAscore", "NeoTCR8", "NeoTCR4"):
    score_by_name(adata, name, log1p=False, key_added=f"signature_{name}")
```

In multi-sample cohorts, normalize and score within each sample before
aggregating clonotypes so sample-wide baseline shifts do not dominate rankings;
see [Clone Prioritization](../user-guide/til-signatures.md).

## Cell-type and state registries

Used by [`tcrsift.annotate_cells`](../user-guide/atlas.md):

| Registry constant | Contents |
| --- | --- |
| `CELL_TYPE_SIGNATURES` / `PBMC_CELL_TYPE_SIGNATURES` | Default blood/PBMC immune and basic stromal cell-type markers |
| `SOLID_TUMOR_CELL_TYPE_SIGNATURES` | PBMC lineages plus solid-tissue stromal, neural, muscle, mesothelial, and osteoclast markers |
| `T_STATE_SIGNATURES`, `B_STATE_SIGNATURES` | T-cell and B-cell functional sub-state markers |
| `PBMC_CULTURE_TYPES` | Cell-type allow-list for Ficoll-isolated PBMC cultures |

Malignant clusters are typed with `tumor_override` / `MarkerCountOverride`
rather than a fixed signature in these registries.

::: tcrsift.signatures
    options:
      members:
        - EFFECTOR_GENES_HGNC
        - ACTIVATION_GENES_HGNC
        - NAIVE_STEM_GENES_HGNC
        - ANTIGEN_RESPONSE_GENES_HGNC
        - AIM_GENES_HGNC
        - CYTOLYTIC_GENES_HGNC
        - EXHAUSTION_GENES_HGNC
        - TUMOR_REACTIVE_GENES_HGNC
        - EXPANSION_CORE_GENES_HGNC
        - MARKER_PANEL_HGNC
        - T_CELL_SIGNATURES

::: tcrsift.signature_methods
    options:
      members:
        - Signature
        - SIGNATURES
        - NEOANTIGEN_SIGNATURES
        - score_by_name
