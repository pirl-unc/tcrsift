# Heuristic Generation Probability (`tcrsift.pgen`)

`tcrsift.pgen` implements a fast, dependency-free generation-probability
heuristic over CDR3 amino-acid length, V/J gene marginals, and estimated
junctional N-insertions, used by `annotate_clonotypes(add_publicness=True)` to
flag likely public CDR3s when reviewing database matches.

For new ranking workflows, prefer the fitted k-mer and TCRpeg models in
[`tcrsift.seqprob`](seqprob.md) (especially observed-repertoire `Ppost`). Do
not transfer numeric cutoffs between `tcrsift.pgen` and `tcrsift.seqprob`, as
their log-probability scales differ.

## Rationale and decomposition

Exact V(D)J dynamic programming in OLGA (Sethna et al., 2019) requires ~30 MB
of model tables and a GPL runtime dependency. When only approximate rank order
between germline-like public CDR3s and high-insertion private CDR3s is needed,
`tcrsift.pgen` decomposes base-10 log probability (`log10_pgen`) as:

```
log10 Pgen(CDR3, V, J, chain)
  ≈ log10 P(V)              # V-gene marginal
  + log10 P(J)              # J-gene marginal
  + log10 P(length)         # Chain-specific CDR3 AA length
  + log10 P(n_inserted)     # Geometric N-insertion model (1 junction in α, 2 in β)
  + log10 P(CDR3 | length)  # Uniform 20-amino-acid baseline
```

Estimated non-templated nucleotide insertions (`n_inserted`, inferred from CDR3
length minus typical germline V/D/J contribution) carry the largest weight, as
zero-insertion rearrangements arise convergently across donors. Both human
(`species="human"`) and mouse (`species="mouse"`, using murine `TRAV`/`TRAJ` and
`TRBV`/`TRBJ` usage priors and case-insensitive gene normalization) are
supported.

## Scale and calibration

Because of the uniform amino-acid baseline, scores sit roughly 10 `log10` units
(orders of magnitude) below exact OLGA `log10 Pgen`. Default `publicness_score`
bounds (`-30` / `-18`) are calibrated to this module's `log10` scale; pass
`-20` / `-8` when scoring exact OLGA `log10 Pgen` outputs, or set
`auto_quantile=True` for cohort-quantile scaling.

## Usage

```python
from tcrsift import annotate_clonotypes
from tcrsift.pgen import annotate_publicness

clones = annotate_publicness(clones)
public = clones[clones["publicness"] >= 0.5]

# Or enable directly during database annotation:
annotated = annotate_clonotypes(
    clones,
    vdjdb_path="vdjdb.tsv",
    add_publicness=True,
)
```

::: tcrsift.pgen
    options:
      members:
        - pgen_single
        - pgen_components
        - compute_pgen
        - publicness_score
        - annotate_publicness
