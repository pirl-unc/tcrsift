# Reference Model and Sequence Provenance

Provenance for packaged background probability models (`tcrsift.seqprob`),
canonical TCR constant-region FASTAs (`tcrsift.assemble`), and IMGT V-gene
leader references (`tcrsift.assemble.leaders`).

## `Ppost`: observed healthy-PBMC background

`kmer_ppost_alpha.npz` and `kmer_ppost_beta.npz` are order-2 amino-acid Markov
models fitted on `observed_pbmc_10x.csv.gz` (bundled in this directory),
regenerated via `scripts/generate_kmer_background.py`. `tcrsift pgen fetch` and
`tcrsift pgen train` also read `observed_pbmc_10x.csv.gz` directly so human α
and β `kmer`, `kmer_gene`, and `tcrpeg` models can be trained offline.

`observed_pbmc_10x.csv.gz` pools de-identified unique CDR3 amino-acid sequences
and V/J gene calls (no cell barcodes or donor identifiers) from five
CC-licensed **10x Genomics healthy-donor VDJ-T** datasets, filtered to
`productive == True & high_confidence == True` and deduplicated per chain
(**15,948 α** and **18,110 β** unique CDR3s):

| Short name | 10x Genomics dataset (`filtered_contig_annotations.csv`) |
| --- | --- |
| `pbmc10k` | `cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t` |
| `ref_pbmc1k` | `cell-vdj/5.0.0/sc5p_v2_hs_PBMC_1k_multi_5gex_5fb_b_t` |
| `ref_pbmc3` | `cell-vdj/3.1.0/vdj_v1_hs_pbmc3` |
| `ref_pbmc3ng` | `cell-vdj/3.1.0/vdj_nextgem_hs_pbmc3` |
| `ref_pbmct` | `cell-vdj/2.2.0/vdj_v1_hs_pbmc_t` |

Source: 10x Genomics Single Cell Immune Profiling public datasets
(<https://www.10xgenomics.com/datasets>).

Note: order-2 CDR3 k-mer log-probability correlates strongly with CDR3 length
(`|r| ≈ 0.9`; longer CDR3s receive lower `Pgen` and `Ppost`), an inherent
property of sequence generation probability models.

## `Pgen`: synthetic generation background

`kmer_pgen_alpha.npz` and `kmer_pgen_beta.npz` are fitted on 500,000 synthetic
CDR3 sequences per chain generated offline using OLGA's default human T-cell α
and β models. OLGA (GPL-3.0) is invoked **offline during model generation
only** and is never imported at runtime, keeping the distributed `tcrsift`
package under Apache-2.0. Only the fitted NumPy transition tables are bundled.

## Canonical TCR constant regions (human and mouse)

`canonical_constants.fasta` and `canonical_constants_alleles.fasta` store bare
mature human `TRAC`, `TRBC1`, and `TRBC2` constant-region amino-acid sequences
(excluding the 1-aa J→C junction splice residue supplied by the rearranged J
segment or donor contig) from Ensembl release 110 (`GRCh38`:
`ENST00000611116` / `TRAC` / UniProt `P01848`, `ENST00000633705` / `TRBC1` /
UniProt `P01850`), `NCBI:AAA60662` (`TRBC2*01` population-major default), and
IMGT allele tables (`TRAC*01`–`*04`, `TRBC1*01`–`*03`, `TRBC2*01`–`*05`).

`canonical_constants_mouse.fasta` and `canonical_constants_alleles_mouse.fasta`
store the corresponding bare mature mouse (*Mus musculus*, C57BL/6J `GRCm39` /
Ensembl release 112) constant-region amino-acid sequences:

- `TRAC` (`Trac`, `ENSMUST00000198398` / UniProt `P01849`, 136 aa, C-terminus
  `LLMTLRLWSS`)
- `TRBC1` (`Trbc1`, `ENSMUST00000192856` / UniProt `P01852`, 172 aa, C-terminus
  `VKRKNS`)
- `TRBC2` (`Trbc2`, `ENSMUST00000103299` / UniProt `A0A075B5J4`, 172 aa,
  C-terminus `VKKKNS`)

Regenerate from local Ensembl caches via
`python scripts/regen_canonical_constants.py --species {human,mouse}`.

## Germline V-gene signal peptides (IMGT/GENE-DB)

`vgene_leaders.csv.gz` stores per-allele human `TRAV` and `TRBV` signal-peptide
references used to anchor `from_contig` leader extraction.
`scripts/build_vgene_leaders.py` builds this table from the vendored
IMGT/GENE-DB `L-PART1+L-PART2` FASTAs in the source repository under
`tcrsift/refseqs/imgt_vgene_leaders/` (`TRAV_L-PART.fasta`, `TRBV_L-PART.fasta`),
retaining clean `ATG`-initiated, in-frame, stop-free leaders across 270 alleles
and 123 genes.

Mouse (*Mus musculus*) IMGT `L-PART1+L-PART2` and `V-REGION` FASTAs can be
fetched on demand into the local cache via
`tcrsift data download --group germline_mouse`.

Source: IMGT/GENE-DB (<https://www.imgt.org/genedb/>, *Homo sapiens* and
*Mus musculus*, queries `8.1 TRAV` and `8.1 TRBV`, labels `L-PART1+L-PART2` and
`V-REGION`, CC BY-NC-ND 4.0; cite Lefranc et al., *Nucleic Acids Res.*). For
valid `from_contig` leaders, IMGT sequences serve as an alignment anchor to
locate start codons and leader boundaries while retaining the donor's native
contig sequence; rejected contig leaders fall back to the germline V-gene leader
by default (`--leader-fallback germline`).
