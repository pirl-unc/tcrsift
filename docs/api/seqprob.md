# Sequence probability (data-driven Pgen)

Data-driven CDR3 generation/occurrence probability — the **precursor-
frequency / publicness** axis. This replaces the brittle OLGA/SONIA runtime
path ([`olga_ppost`](olga_ppost.md)) with a model **fit once on an external
reference repertoire**, so `log_pgen(seq)` is a fast, calibrated score for
"how generatable / common is this CDR3" — lower = rarer precursor / more
private.

On the B1-2 pilot the default k-mer model recovers the same alpha-chain
publicness signal as OLGA (TRAV12-2 vs other TRAV **AUROC ≈ 0.64**, matching
OLGA's 0.65–0.67) — with **no GPL and no runtime dependency beyond numpy**.

## Two backends, one interface

Both implement `SequenceProbabilityModel` (`fit` / `log_prob` / `save` /
`load`):

| Backend | Deps | Notes |
| --- | --- | --- |
| `KmerProbabilityModel` (**default**) | numpy only | Order-`k` Markov model over CDR3 AAs (length captured via an EOS symbol; add-`alpha` smoothing). Ships a default for each chain. |
| `TCRpegProbabilityModel` | `pip install tcrsift[tcrpeg]` | Wraps TCRpeg (autoregressive, PyTorch). Heavier, better-calibrated. |

## Shipped defaults & the GPL boundary

The default k-mer models (`tcrsift/refseqs/kmer_background_{alpha,beta}.npz`,
~300 KB each) are fit **offline at build time** on OLGA-generated synthetic
repertoires (`scripts/generate_kmer_background.py`). OLGA (GPL-3.0) is used
*only* to produce training sequences; tcrsift never imports it at runtime,
so the package stays Apache-2.0.

To retrain on a different reference, fit and `save` your own:

```python
from tcrsift.seqprob import KmerProbabilityModel
model = KmerProbabilityModel(order=3, chain="beta").fit(my_reference_cdr3s)
model.save("my_beta_background.npz")
```

## Usage

```python
from tcrsift.seqprob import score_log_pgen, load_background_model

# default shipped k-mer background:
clones["log_pgen_beta"] = score_log_pgen(clones, chain="beta")

# explicit model / TCRpeg backend:
clones["log_pgen_alpha"] = score_log_pgen(clones, chain="alpha", backend="tcrpeg")
```

CLI:

```bash
tcrsift log-pgen clones.csv -o clones_pgen.csv --chain both          # k-mer
tcrsift log-pgen clones.csv -o out.csv --backend tcrpeg --chain beta # TCRpeg
```

The resulting `log_pgen_<chain>` column plugs directly into the in-silico
filter ([`insilico_filter`](insilico_filter.md)) as a `Ppost^low`-style
predicate.

::: tcrsift.seqprob
    options:
      members:
        - SequenceProbabilityModel
        - KmerProbabilityModel
        - TCRpegProbabilityModel
        - load_background_model
        - score_log_pgen
