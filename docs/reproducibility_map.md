# Signature reproducibility map

This page summarizes which gene-expression signatures replicate across donors
in stimulated T-cell cultures, and when to use TCR sequence probability instead
of RNA. These recommendations come from TCRsift's two-donor in-vitro MART-1
peptide-expansion pilot (paired scRNA-seq and VDJ sequencing across multiple
sorting and enrichment conditions).

You can query the same guidance metadata programmatically:

```python
from tcrsift.signature_methods import (
    SIGNATURE_GUIDANCE,
    RNA_REPRODUCIBILITY_NOTE,
    recommended_signatures,
)

recommended_signatures()                 # ['Differentiated', 'AcuteActivation']
SIGNATURE_GUIDANCE["Differentiated"].use_for
```

## Recommended for stimulated cultures

| Signature | Primary use | Evidence from the two-donor pilot |
| --- | --- | --- |
| **Differentiated** (`IFNG, GZMB, PRF1, GNLY, NKG7` minus `TCF7, LEF1, CCR7, SELL`) | Identifying clones that expanded in culture | Strongest and most consistent RNA axis: separates expanded from unexpanded clones at AUROC up to 0.85 in the same direction across both donors. |
| **AcuteActivation** (`TNFRSF9`, `MKI67`) | Recent antigen stimulation and proliferation | Both genes tracked culture expansion across compartments; scores were also lower among public `TRAV12-2` clones in both donors. Best used alongside expansion or sequence publicness rather than alone. |

## Situational signatures

- **Cytolytic (`PRF1, GZMB`) and AntigenExperienced (`IFNG, GZMB, PRF1, GNLY,
  NKG7`).** Measure effector differentiation without subtracting naive/stem
  markers, and correlate strongly with AIM-sorted populations by design. Prefer
  `Differentiated` when ranking clonal expansion in culture.
- **AIM and AIMBroad (`TNFRSF9, TNFRSF4` and extended panel).** Transcriptional
  counterpart of surface 4-1BB/OX40 sorting; useful for checking concordance
  between FACS gating and mRNA capture rather than as an independent selection
  axis.

## Mismatched to short-term PBMC culture

- **TumorReactive and chronic exhaustion panels (`ENTPD1/CD39`, `CXCL13`,
  `TOX`, `ITGAE/CD103`, `PDCD1`, `LAG3`, `HAVCR2`, `TIGIT`).** Reflect chronic
  in-vivo antigen exposure and tissue residency in fresh TIL
  ([public TIL benchmark](benchmarks.md)). In short-term peptide-stimulated
  PBMC cultures, these genes behave largely as noise.

## What RNA does not recover: sequence publicness and cross-reactivity

Transcriptional state does **not** reliably separate private, high-specificity
receptors from high-precursor-frequency public receptors (such as germline-biased
`TRAV12-2` MART-1 binders). Most gene sets flip sign between donors on that
task, and a cross-donor transcriptome classifier (with TCR V/J genes excluded)
reaches only **AUROC 0.62–0.69** for sequence publicness using a diffuse mix of
interferon, stress, and mucosal genes (`ZBTB16`, `ISG15`, `IL7R`, `CCL20`)
rather than a coherent T-cell state.

!!! important "Use sequence probability, not RNA, to filter public receptors"
    To deprioritize common or easily generated CDR3 sequences, use repertoire
    generation and observation probabilities (`Pgen` and `Ppost` in
    [`tcrsift.seqprob`](api/seqprob.md)), not gene expression. Gene expression
    informs **differentiation and activation state**, whereas TCR sequence
    models inform **precursor frequency and publicness**.

## Avoiding overfitting on small pilots

Do not select or trim signature gene lists to maximize separation on the cohort
being analyzed. With two donors, half of any arbitrary seven-gene panel will
agree in sign between donors by coin flip alone (~3.5 of 7 genes), so pruning a
signature on pilot labels overfits quickly. Define gene panels from external
biological literature and evaluate them on independent cohorts.
