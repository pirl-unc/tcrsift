# Pgen + Ppost — OLGA/SONIA gold-standard (optional)

The reproducible **precursor-frequency / publicness** axis for clone
selection (#143). To pick high-specificity rare-precursor clones and
de-prioritize broadly cross-reactive high-precursor ones, the
discriminating quantity is set by V(D)J recombination + thymic selection —
the TCR *sequence*, not the transcriptome.

On the B1-2/B1-3 MART-1 pilot, `log10 Ppost` reproducibly separated the
public (TRAV12-2) repertoire from private clones across both donors
(AUROC 0.72 / 0.72), while every RNA signature flipped direction.

- **Pgen** (OLGA, Sethna 2019) — generation probability of the CDR3.
- **Ppost** (SONIA, Sethna 2020) — post-selection probability `Pgen · Q`.
  `Ppost > Pgen` for public clones: selection *amplifies* their publicness.

## Optional GPL extra

OLGA and SONIA are **GPL-3.0**. To keep tcrsift Apache-2.0 they are a
**runtime-optional dependency, not vendored**:

```bash
pip install tcrsift[pgen]
```

Without it, the numpy-only proxy in [`tcrsift.pgen`](pgen.md) remains the
default for ranking. `tcrsift.olga_ppost.olga_sonia_available()` reports
whether the extra is installed; the compute functions raise an
`ImportError` with the install hint when it isn't.

## Gene-name handling

OLGA/SONIA anchors are allele-suffixed (`TRAV12-2*01`); CellRanger gene
calls often lack the allele. `normalize_gene_name` appends `*01` when
missing (validated).

Alleles not recognized **even after `*01`** — novel alleles from the
audit, or genes outside the model — are mapped to the **nearest supported
allele** before scoring (#150), instead of OLGA silently falling back to a
default usage mask. `compute_pgen_ppost(..., map_unsupported=True)` does
this and logs every substitution (`original → substituted`, reason); pass
`used_v_col` / `used_j_col` to record the allele actually used.

The mapping tiers (`nearest_supported_allele`):

1. `exact` — supported after appending `*01`; no substitution.
2. `nearest_allele` — same gene, different allele → lowest supported allele.
3. `nearest_gene` — gene unsupported → closest supported gene by
   locus/family/subgroup name distance, as its `*01`.
4. `unmapped` — no candidate (→ NaN Pgen, never a silent default).

`annotate_nearest_supported_allele(df, chain=, segment=, gene_col=)` adds
the mapping as columns so a detected novel allele carries its
downstream-Pgen substitution auditably.

## Caveats

- Pgen/Ppost estimate **naive precursor frequency**, a *proxy* for
  cross-reactivity risk — not a direct measure.
- They say nothing about **functional avidity** (tetramer MFI / koff) —
  that needs a wet-lab assay and is out of scope.

## Usage

```python
from tcrsift.olga_ppost import compute_pgen_ppost, flag_private_candidates

clones = compute_pgen_ppost(clones, chain="beta")
clones["private_candidate"] = flag_private_candidates(
    clones, score_col="log10_ppost", freq_col="frequency",
)
```

CLI:

```bash
tcrsift ppost clones.csv -o clones_ppost.csv --chain both
```

::: tcrsift.olga_ppost
    options:
      members:
        - olga_sonia_available
        - load_chain_model
        - supported_alleles
        - normalize_gene_name
        - nearest_supported_allele
        - annotate_nearest_supported_allele
        - compute_pgen_ppost
        - flag_private_candidates
