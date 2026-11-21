# Sequence Probability (`tcrsift.seqprob`)

`tcrsift.seqprob` scores CDR3 amino-acid sequences under pre-fitted repertoire
background models to estimate precursor frequency and sequence publicness:

- **`Pgen`:** generation probability under an OLGA-simulated naive repertoire.
- **`Ppost`:** observation probability under a pooled healthy-donor PBMC
  repertoire.

Both scores are cohort-ranking signals rather than antigen-specificity or
avidity measurements. Returned values are natural-log probabilities (`ln P`);
more negative values indicate rarer, more private sequences.

## Backends and packaged models

Both backends implement `SequenceProbabilityModel` (`fit`, `log_prob`, `save`,
`load`) and are included in the standard `pip install tcrsift` installation:

| Backend | Description |
| --- | --- |
| `KmerProbabilityModel` (default) | Order-`k` amino-acid Markov model with explicit stop symbol (encoding CDR3 length) and optional V/J marginals |
| `TCRpegProbabilityModel` | Autoregressive PyTorch model for training on larger custom repertoires |

Packaged k-mer parameter files cover both roles and chains:

- `kmer_pgen_alpha.npz`, `kmer_pgen_beta.npz` (fitted offline on 500,000
  OLGA-generated sequences per chain; OLGA is not imported at runtime)
- `kmer_ppost_alpha.npz`, `kmer_ppost_beta.npz` (fitted on pooled 10x Genomics
  healthy-donor PBMC datasets; see `tcrsift/refseqs/PROVENANCE.md`)

## Usage

Score both roles on α and β chains:

```python
from tcrsift import add_pgen_ppost

clones = add_pgen_ppost(clones, backend="kmer")
# Columns added: pgen_alpha, ppost_alpha, log_q_alpha, pgen_beta, ppost_beta, log_q_beta
```

Score a single chain and role:

```python
from tcrsift.seqprob import score_log_prob

clones["log_ppost_beta"] = score_log_prob(clones, chain="beta", role="ppost")
```

Command line:

```bash
tcrsift pgen annotate clones.csv -o clones_probability.csv --chain both
```

Use `Ppost` as a relative rank within a cohort (for example via PRISM or
`--exclude-public-quantile`), rather than applying fixed numeric cutoffs across
different models. Because log sequence probability correlates strongly with CDR3
length, inspect CDR3 length and V/J usage alongside `Ppost`.

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
