# TCR annotation API + PRISM (#157/#158)

A uniform, pluggable way to add optional, **always-defined** per-clone
columns to any TCR table, plus the **PRISM** selection method.

## Pgen / Ppost — robust, never-zero (#157)

The OLGA/SONIA runtime path returned `Pgen = 0` on ~3.5–5% of α junctions
(CDR3s not ending in the conserved F/W — an anchor parse failure), and
Ppost inherited it, mis-ranking exactly the rare clones of interest. The
data-driven backends in [`seqprob`](seqprob.md) are **finite for every
valid junction**.

Two **roles** per chain:

- **Pgen** — model fit on an OLGA-generated (synthetic) repertoire =
  generation probability.
- **Ppost** *(default)* — model fit on an **observed** repertoire =
  post-selection publicness. `log_q = log Ppost − log Pgen` is the
  data-driven selection factor (`Q = Ppost/Pgen`).

```python
from tcrsift.annotate_tcrs import add_pgen_ppost
clones = add_pgen_ppost(clones)   # pgen_/ppost_/log_q_ {alpha,beta}
```

Shipped k-mer defaults: Pgen (α+β, OLGA-generated) and Ppost (β, observed).
**No observed α reference is bundled** → α-Ppost falls back to α-Pgen (finite,
not circular) with a logged note until one is supplied. TCRpeg is the
higher-accuracy backend (`backend="tcrpeg"`, needs `pip install
tcrsift[tcrpeg]`).

## GEX signature scores, z-scored per sample/donor

`antigen_response_score` (TNFRSF9, MKI67) and `naive_score` (TCF7, LEF1,
CCR7, SELL, IL7R, CD27, CD28) from per-cell expression, **z-scored within
each sample or donor group** (#144/#145) before averaging per clone:

```python
from tcrsift.annotate_tcrs import add_gex_signature_scores
clones = add_gex_signature_scores(clones, per_cell, group_col="sample")
```

`per_cell` has `gex.<SYMBOL>` columns + a clone column + the group column.

## PRISM — Percentile-Rank In-Silico Multi-criterion

Rank every clone by the **mean of percentile ranks** of: low `ppost_alpha`,
low `ppost_beta`, high `antigen_response_score`, low `naive_score`; take the
top-K. Pilot result: ~33% / 45% condition-private clones vs 6% / 17% for
frequency selection.

```python
from tcrsift.annotate_tcrs import annotate_tcrs, select_prism

annotated = annotate_tcrs(clones, per_cell=per_cell, group_col="sample",
                          backend="kmer")
picks = select_prism(annotated, k=20, group_col="enrichment_method")
```

Per-dimension `weights` and `K` are tunable; PRISM composes with the
in-silico filter layer ([`insilico_filter`](insilico_filter.md)).

::: tcrsift.annotate_tcrs
    options:
      members:
        - annotate_tcrs
        - add_pgen_ppost
        - add_gex_signature_scores
        - score_gex_signature_per_clone
        - naive_signature
        - prism_score
        - select_prism
