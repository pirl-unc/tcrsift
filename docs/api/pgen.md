# Generation Probability (Pgen) API

Lightweight generation-probability estimator for TCR β/α CDR3
sequences. Used to flag likely-public CDR3s so DB matches against
them can be appropriately discounted (#58).

This is the older, interpretable V/J/length/insertion heuristic behind
`annotate_clonotypes(add_publicness=True)`. New selection workflows normally
use the externally fitted Pgen/Ppost backgrounds in
[Sequence probability](seqprob.md), especially observed-repertoire Ppost.
Their numerical scales are different; do not transfer thresholds between
them.

## Why a proxy, not OLGA

The gold-standard tool is **OLGA** (Sethna 2019), which fits explicit
V/D/J/insertion models and computes Pgen by dynamic programming.
OLGA gives the right answer but ships ~30 MB of pre-trained model
files plus a runtime dependency chain.

Most tcrsift use cases (publicness discounting, ranking convergent
vs. private sequences) just need the **rank order**, not numerically
calibrated Pgens. This module provides that proxy in ~250 lines
with no runtime dependency beyond numpy.

## Decomposition

```
log Pgen(CDR3, V, J, chain)
  ≈ log P(V)             # V-gene usage marginal
  + log P(J)             # J-gene usage marginal
  + log P(length)        # CDR3 AA length, per chain
  + log P(n_inserted)    # total non-templated nt
  + log P(CDR3 | length) # uniform-AA composition baseline
```

The **N-insertion term** is the largest contributor to publicness:
sequences with 0 non-templated nucleotides are convergently
generatable across donors and dominate public-match noise. We model
N-insertion length as geometric per junction (β has two junctions,
α has one); the total N-insertion count is estimated from observed
CDR3 nt length minus a typical templated V+J+D contribution.

## Calibration

The estimator runs ~10 log units lower than real OLGA Pgens because
of the uniform-over-20 AA composition baseline. Default
`publicness_score` cutoffs are calibrated to **this estimator's
scale** (`-30 / -18`); pass OLGA cutoffs (`-20 / -8`) if feeding in
real OLGA values, or use `auto_quantile=True` for distribution-based
auto-calibration.

## Usage

```python
from tcrsift.pgen import annotate_publicness

clones = annotate_publicness(clones)
public = clones[clones["publicness"] >= 0.5]
print(f"{len(public)} likely-public clones; DB matches should be discounted")
```

Wired into `annotate_clonotypes`:

```python
from tcrsift import annotate_clonotypes

annotated = annotate_clonotypes(
    clones, vdjdb_path="vdjdb.tsv", add_publicness=True,
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
