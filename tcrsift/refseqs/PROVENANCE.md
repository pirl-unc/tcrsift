# Reference model provenance

Shipped background models for data-driven Pgen/Ppost (`tcrsift.seqprob`,
gene-agnostic order-2 k-mer). Regenerate with
`scripts/generate_kmer_background.py`.

## Ppost — observed healthy-PBMC background (publicness)

`kmer_ppost_alpha.npz` / `kmer_ppost_beta.npz` are fit on
**`observed_pbmc_10x.csv.gz`** (bundled here): the pooled, de-identified unique
CDR3 amino-acid sequences (with V/J gene calls — no barcodes or donor IDs) from
five **10x Genomics healthy-donor VDJ-T** datasets, CC-licensed. Both α and β
are drawn from the **same** pooled source.

Filter applied: `productive == True & high_confidence == True`, then deduplicated
to unique CDR3 per chain → **15,948 α / 18,110 β** unique CDR3.

| component file | 10x dataset (`filtered_contig_annotations.csv`) |
|----------------|--------------------------------------------------|
| pbmc10k        | `cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t` |
| ref_pbmc1k     | `cell-vdj/5.0.0/sc5p_v2_hs_PBMC_1k_multi_5gex_5fb_b_t` |
| ref_pbmc3      | `cell-vdj/3.1.0/vdj_v1_hs_pbmc3` |
| ref_pbmc3ng    | `cell-vdj/3.1.0/vdj_nextgem_hs_pbmc3` |
| ref_pbmct      | `cell-vdj/2.2.0/vdj_v1_hs_pbmc_t` |

Source: 10x Genomics, Single Cell Immune Profiling public datasets
(https://www.10xgenomics.com/datasets). This is the same neutral background the
downstream `selection_analysis/build_feat.py` pilot uses, so tcrsift's native
Ppost is consistent with the validated selection recipe.

Note: order-2 CDR3 k-mer log-P is strongly anti-correlated with CDR3 length
(longer CDR3 → lower Pgen/Ppost, ~0.9 |corr|), so the `ppost_*` PRISM axes
partly encode length, not only rarity. This is inherent to generation-probability
scores; see `prism_score` / `select_freq_prism_per_condition` docs.

## Pgen — synthetic generation background

`kmer_pgen_alpha.npz` / `kmer_pgen_beta.npz` are fit on 500,000 OLGA-generated
synthetic CDR3s per chain (OLGA default human T α/β models). OLGA (GPL-3.0) is
used **offline at generation time only** — never imported at runtime — so the
shipped package stays Apache-2.0. The fitted numpy tables are committed; the
synthetic sequences are not bundled.
