# TCR Annotation and PRISM Selection

`tcrsift.annotate_tcrs` attaches standardized sequence-probability and
gene-expression columns to clonotype tables and implements the **PRISM**
(Percentile-Rank In-Silico Multi-criterion) ranking score.

## Sequence probability (`Pgen` and `Ppost`)

Runtime OLGA/SONIA scoring previously returned `Pgen = 0` on 3.5–5% of α-chain
CDR3s lacking a canonical C-terminal `F`/`W` anchor, mis-ranking rare clones.
The k-mer and TCRpeg models in [`tcrsift.seqprob`](seqprob.md) assign finite
natural-log probabilities to every valid CDR3 sequence across two reference
roles:

- **`Pgen` (generation probability):** fitted offline on OLGA-simulated human
  α and β repertoires.
- **`Ppost` (observed publicness, default):** fitted on pooled healthy-donor
  10x Genomics PBMC repertoires. `log_q = log Ppost - log Pgen` estimates the
  log post-thymic selection factor (`Q = Ppost / Pgen`).

For backward compatibility, output columns are named `pgen_alpha`,
`ppost_alpha`, `pgen_beta`, and `ppost_beta`, and store **natural-log
probabilities** (`ln P`). Less negative values indicate more common (public)
sequences; more negative values indicate rarer (private) sequences.

```python
from tcrsift.annotate_tcrs import add_pgen_ppost

# Adds pgen_alpha, ppost_alpha, log_q_alpha, pgen_beta, ppost_beta, log_q_beta
clones = add_pgen_ppost(clones)
```

Bundled k-mer models run with NumPy alone; pass `backend="tcrpeg"` to score
with the autoregressive PyTorch model instead.

## Per-sample z-scored GEX signatures

`add_gex_signature_scores` computes `antigen_response_score` (`TNFRSF9`,
`MKI67`) and `naive_score` (`TCF7`, `LEF1`, `CCR7`, `SELL`, `IL7R`, `CD27`,
`CD28`) from per-cell expression, standardizing each gene within `group_col`
(sample or donor) before averaging across cells of each clone:

```python
from tcrsift.annotate_tcrs import add_gex_signature_scores

# per_cell contains gex.<SYMBOL> columns, a clone identifier, and group_col
clones = add_gex_signature_scores(clones, per_cell, group_col="sample")
```

## PRISM multi-criterion ranking

`prism_score` ranks clones by the weighted mean of within-group percentile
ranks across four criteria (`ppost_alpha` low, `ppost_beta` low,
`antigen_response_score` high, and `naive_score` low), and `select_prism`
flags the top `k` clones per group in `prism_selected`. In the two-donor
MART-1 culture pilot, top-`k` PRISM selection yielded 33% and 45%
condition-private clones versus 6% and 17% for raw frequency ranking.

```python
from tcrsift.annotate_tcrs import annotate_tcrs, select_prism

annotated = annotate_tcrs(
    clones,
    per_cell=per_cell,
    group_col="sample",
    backend="kmer",
)
scored = select_prism(annotated, k=20, group_col="enrichment_method")
picks = scored[scored["prism_selected"]]
```

Per-dimension `weights` and `k` are configurable, and PRISM composes with
[`tcrsift.insilico_filter`](insilico_filter.md).

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
