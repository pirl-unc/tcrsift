# Known Limitations and Roadmap

## Species support

TCRsift supports **human** (*Homo sapiens*) and **mouse** (*Mus musculus*) paired
α/β TCR and single-cell gene-expression datasets across loading, mitochondrial
QC, CD4/CD8 phenotyping, gene-expression signature scoring, clone
prioritization, single-cell atlas annotation, public TCR database matching, and
full-length sequence assembly (including chimeric human-VDJ / murine-constant
constructs via `--constant-species mouse`).

Current species boundaries:

1. **Other vertebrates and γδ VDJ assembly:** non-human primates, rats, swine,
   and γδ TCR contig assembly are not bundled (`--constant-source from-data` can
   supply custom constant-region sequences for α/β tables).
2. **Published human neoantigen transcriptomic gene sets (`NeoTCR4`, `NeoTCR8`,
   `NeoTCR_PBL`):** derived from human TIL and peripheral blood cohorts; mouse
   `prioritize` presets omit them by default rather than projecting multi-hundred
   gene sets across paralogous families.
3. **Mouse repertoire background models (`tcrsift.seqprob`):** packaged k-mer
   weights (`kmer_pgen_*.npz`, `kmer_ppost_*.npz`) are trained on human OLGA and
   10x PBMC repertoires. Fit a dedicated murine background model from a local
   repertoire TSV via `tcrsift pgen train --species mouse --repertoire <path>`.

## Future directions

- Bundled murine naive/splenocyte `Pgen` and `Ppost` background models
- TCR sequence-similarity clustering across clonotypes
- Depth-normalized CD3 gating and adaptive shortlist budgets for large cohorts
- Expanded validation cohorts for `blood-tumor`, `mpe`, and `heme`
  prioritization presets (for peripheral-blood, pMHC-multimer, public-database,
  and peptide-restimulation benchmarks across Influenza, SARS-CoV-2, EBV, and
  CMV, see [docs/benchmarks-infectious.md](docs/benchmarks-infectious.md))
