# Known Limitations and Roadmap

## Species support

TCRsift currently targets **human** α/β TCR and GEX datasets:

1. **Constant-region references (`tcrsift/assemble.py`):** packaged canonical
   sequences (`tcrsift/refseqs/canonical_constants.fasta`) cover human `TRAC`,
   `TRBC1`, and `TRBC2`, and C-terminus QC checks validate human sequence
   endings.
2. **Leader sequences (`tcrsift/assemble.py`, `tcrsift/leaders.py`):** bundled
   IMGT V-gene leader tables (`vgene_leaders.csv.gz`) are *Homo sapiens*
   `TRAV`/`TRBV` alleles, and predefined signal peptides (`CD8A`, `CD28`,
   `TRAC`, `TRBC`, plus murine `IgK`) use human codon optimization.
3. **Gene symbols and background models:** `loader.py`, `phenotype.py`,
   `signatures.py`, and `seqprob.py` assume human HGNC symbols (`CD3D/E/G`,
   `CD4`, `CD8A/B`, `TRAV`/`TRBV`) and human PBMC/OLGA CDR3 backgrounds.

### Multi-species design sketch

Supporting mouse (*Mus musculus*) or non-human primate data requires
parameterizing species-specific reference tables:

- Constant-region FASTAs (`Trac`, `Trbc1`, `Trbc2`) and murine codon tables
- IMGT leader tables and `Tr[ab][vdjc]` gene-name regexes
- Mouse marker symbols (`Cd3d`, `Cd3e`, `Cd3g`, `Cd4`, `Cd8a`, `Cd8b1`) and
  signature ortholog mappings
- Species filtering when matching VDJdb and IEDB

If you need support for a non-human species, open an issue at
<https://github.com/pirl-unc/tcrsift/issues> describing your assay format and
reference build.

## Future directions

- TCR sequence-similarity clustering across clonotypes
- Depth-normalized CD3 gating and adaptive shortlist budgets for large cohorts
- Expanded validation cohorts for `blood`, `mpe`, and `heme` prioritization
  presets
