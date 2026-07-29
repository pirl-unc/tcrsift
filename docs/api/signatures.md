# Gene-expression signatures

TCRsift has two related registries:

- `T_CELL_SIGNATURES`: compact functional programs used for exploratory
  per-cell/per-clone scoring and TIL selection.
- `NEOANTIGEN_SIGNATURES`: published neoantigen-reactivity signatures with
  their own scoring methods, input units, and citations.

Scores describe cell state. They do not establish which antigen a TCR binds.
Interpret the same program in its source compartment: AIM is most useful after
in-vitro stimulation, while tumor-reactive/exhaustion programs are most useful
in fresh TIL.

## Compact functional programs

All symbols are HGNC human gene names.

| Name | Genes / contrast | Intended use |
| --- | --- | --- |
| `effector` | IFNG, GZMB, PRF1, GNLY, NKG7 | Cytotoxic-effector state |
| `activation` | Deprecated alias of `effector` | Backward compatibility only |
| `naive_stem` | TCF7, LEF1, CCR7, SELL, IL7R, CD27, CD28 | Naïve/stem-memory state |
| `antigen_response` | TNFRSF9, MKI67 | Focal recent-response/proliferation readout |
| `aim` | TNFRSF9, TNFRSF4, IL2RA, MKI67 | Activation-induced-marker program in stimulated culture |
| `cytolytic` | PRF1, GZMB | Minimal cytolytic readout |
| `exhaustion` | PDCD1, TOX, LAG3, HAVCR2, TIGIT | Chronic-antigen/exhaustion state in TIL |
| `tumor_reactive` | CXCL13, ENTPD1 | Compact tumor-reactive TIL phenotype |
| `expansion_core` | MKI67, TNFRSF9, EGR2, IFNG, CXCL13, HAVCR2 | Cross-compartment expansion-associated panel |

`MARKER_PANEL_HGNC` is the union extracted by `til-select`; it is a display and
scoring panel, not another biological signature.

## Published neoantigen-reactivity registry

`score_by_name` dispatches each entry according to its registered method
instead of treating unlike signatures as interchangeable weighted sums.

| Name | Source | TCRsift implementation |
| --- | --- | --- |
| `MANAscore` | Zeng/Smith, Nature Communications 2025 (PMID 39900903) | Transparent signed-z proxy: `(z(CXCL13) + z(ENTPD1) - z(IL7R)) / sqrt(3)`. The paper's trained RF+linear ensemble is not shipped and this proxy must not be described as the original fitted model. |
| `NeoTCR8` | Lowery/Rosenberg, Science 2022 (PMID 35113651) | Published 243-gene CD8 set, unweighted gene-set enrichment |
| `NeoTCR4` | Lowery/Rosenberg, Science 2022 (PMID 35113651) | Published 40-gene CD4 set, unweighted gene-set enrichment |
| `NeoTCR_PBL` | Yossef/Rosenberg, Cancer Cell 2023 (PMID 38039963) | Published 151-gene circulating CD8 set, unweighted gene-set enrichment; use for blood, not as the default TIL signature |

For `NeoTCR4`, `NeoTCR8`, and `NeoTCR_PBL`, pass an AnnData object containing
the full log-normalized gene universe. A bare DataFrame cannot reproduce
rank-enrichment controls, so TCRsift warns and falls back to a mean-z proxy.

## Usage

Compact unsigned programs:

```python
from tcrsift import compute_signature_scores_per_cell

compute_signature_scores_per_cell(adata)
print([column for column in adata.obs if column.startswith("signature_")])
```

Published/signed methods:

```python
import scanpy as sc
from tcrsift import score_by_name

sc.pp.normalize_total(adata, target_sum=10_000)
sc.pp.log1p(adata)

for name in ("MANAscore", "NeoTCR8", "NeoTCR4"):
    # X is already log1p(CP10K), so do not log it again.
    score_by_name(adata, name, log1p=False, key_added=f"signature_{name}")
```

For multiple samples, score within each sample or donor before aggregating
clones; otherwise sample-level expression shifts can dominate the result. See
[Multi-sample TIL prioritization](../user-guide/til-signatures.md) for a
complete example.

## Cell-type and state registries

The per-cell annotator (`tcrsift.annotate_cells`) also exposes:

| Constant | What it contains |
| --- | --- |
| `CELL_TYPE_SIGNATURES` | Blood/PBMC-oriented default cell-type registry |
| `PBMC_CELL_TYPE_SIGNATURES` | Explicit alias for that default |
| `SOLID_TUMOR_CELL_TYPE_SIGNATURES` | Default plus solid-tissue lineages |
| `T_STATE_SIGNATURES`, `B_STATE_SIGNATURES` | T- and B-cell sub-state registries |
| `PBMC_CULTURE_TYPES` | Type allow-list for PBMC/culture annotation |

Malignant-cell calls are handled by `MarkerCountOverride`/`tumor_override`,
not by adding a “malignant” expression signature to this registry.

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
