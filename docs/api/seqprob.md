# Sequence probability: Pgen and Ppost

Sequence-background scores provide a **precursor-frequency/publicness axis**:

- Pgen asks how probable a CDR3 is under a generated-repertoire background.
- Ppost asks how common it is under an observed healthy-repertoire background.

Both are ranking proxies, not specificity or avidity measurements. The shipped
k-mer scores are natural-log probabilities; lower (more negative) values mean
a rarer/private sequence.

## Backends and shipped models

Both backends implement `SequenceProbabilityModel` (`fit`, `log_prob`, `save`,
and `load`):

| Backend | Notes |
| --- | --- |
| `KmerProbabilityModel` (default) | Fast order-k Markov model over CDR3 amino acids, with length represented by an end symbol and optional V/J marginals |
| `TCRpegProbabilityModel` | Autoregressive PyTorch model; heavier and useful when trained on a sufficiently large reference |

TCRpeg and PyTorch are core dependencies. No `tcrsift[tcrpeg]` extra is
required.

Packaged files are role- and chain-specific:

- `kmer_pgen_alpha.npz`, `kmer_pgen_beta.npz`
- `kmer_ppost_alpha.npz`, `kmer_ppost_beta.npz`

Pgen references were generated offline with OLGA; tcrsift does not import OLGA
at runtime. Ppost references come from pooled observed healthy PBMC
repertoires. See `tcrsift/refseqs/PROVENANCE.md` in the source distribution for
the exact data and calibration notes.

## Usage

Add both roles per chain:

```python
from tcrsift import add_pgen_ppost

clones = add_pgen_ppost(clones, backend="kmer")
# pgen_alpha, ppost_alpha, pgen_beta, ppost_beta are ln(probability)
```

Score only one role:

```python
from tcrsift.seqprob import score_log_prob

clones["log_ppost_beta"] = score_log_prob(
    clones,
    chain="beta",
    role="ppost",
)
```

CLI:

```bash
tcrsift pgen annotate clones.csv -o clones_probability.csv --chain both
```

Use Ppost as a cohort-relative ranking signal or through
`select_specificity_candidates`; do not copy a threshold calibrated on a
different model/reference without validation. Ppost also partly reflects CDR3
length, so review length and V/J usage alongside it.

::: tcrsift.seqprob
    options:
      members:
        - SequenceProbabilityModel
        - KmerProbabilityModel
        - GeneAwareKmerModel
        - TCRpegProbabilityModel
        - load_background_model
        - score_log_prob
        - score_log_pgen
