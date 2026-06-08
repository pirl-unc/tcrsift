# Signature reproducibility map — what RNA can and can't do

The curated gene-signature set for clone selection, with an **honest map of
what each is good for** (#145) — so you don't reach for RNA where the *TCR
sequence* is the right tool. Derived from the B1-2/B1-3 MART-1 in-vitro
expansion pilot (paired scRNA + VDJ, two donors).

This map is queryable in code:

```python
from tcrsift.signature_methods import (
    SIGNATURE_GUIDANCE, recommended_signatures, RNA_REPRODUCIBILITY_NOTE,
)
recommended_signatures()        # ['Differentiated', 'AcuteActivation']
SIGNATURE_GUIDANCE["Differentiated"].use_for
```

## Recommended — reproducible across donors

| Signature | Use for | Why it's reproducible |
| --- | --- | --- |
| **Differentiated** (effector − naïve/stem, signed; #141) | **clonal expansion** — which clones expanded in culture | Best RNA axis: separates expanded vs rare clones at AUROC up to 0.85, same direction in both donors. |
| **AcuteActivation** (TNFRSF9/4-1BB, MKI67; #142) | the **one reproducible RNA correlate of publicness** | Consistently *lower* in public (TRAV12-2) clones across all three normalizations and both donors (modest but stable). Private/specific clones carry more recent cognate activation + proliferation; the public pool is more bystander-like. Runs *inverse* to expansion. |

## Situational — document, don't default

- **Cytolytic / AntigenExperienced** — effector program; correlated with an
  AIM sort by construction (#142). Prefer Differentiated for the expansion
  question.
- **AIM / AIMBroad** — the transcriptional analogue of an AIMpos sort; use
  to check sort/transcriptome concordance, not as an independent axis.

## Wrong biology for short PBMC culture

- **TumorReactive, exhaustion (CD39 / CXCL13 / TOX / CD103)** — chronic
  in-situ exhaustion biology. Reads as **noise** in short PBMC culture;
  meaningful only in fresh tissue TILs.

## Documented limitation

**RNA does not reproducibly recover the precursor-frequency /
cross-reactivity axis.** Most signatures *flip direction between donors*; a
cross-donor-validated transcriptome classifier (TCR genes removed, per
#144) reaches only **AUROC ~0.62–0.69** for publicness, on a diffuse,
non-interpretable gene set (ZBTB16, ISG15, IL7R, CCL20, stress genes — not
a naïve program).

> For **high-specificity / high-avidity / low-cross-reactivity** selection,
> use **Pgen/Ppost + sequence features** ([`tcrsift.seqprob`](api/seqprob.md),
> #143), **not RNA**. RNA is for **differentiation / expansion state**, not
> specificity.

## Methodology note

Do **not** fit/trim signatures on the selection target. With n=2 donors,
per-gene direction "consistency" is ~chance (≈3.5/7 agree by coin-flip), so
trimming a signature to the genes that separate a label *here* is
overfitting. Refine signatures from curated external references and
validate cross-donor / cross-cohort.
