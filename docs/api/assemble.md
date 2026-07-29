# Assembly API

Build full-length TCR sequences and synthesis-ready β-2A-α constructs.

## Overview

The input must contain assembled V(D)J segment sequence, as produced by
TCRsift clonotype aggregation; CDR3 and V/J names alone are not enough to
reconstruct the full variable region.

The module supports:

- **Leader sequences**: Per-chain configuration - extract from contig FASTAs or use standard signal peptides
- **Constant regions**: Packaged canonical TRAC/TRBC references or explicit sequence columns from the input
- **2A linkers**: Join alpha and beta chains with self-cleaving peptides (T2A, P2A, E2A, F2A)
- **Single-chain constructs**: Generate β-linker-α format for expression

## Leader Sequence Options

Each chain (alpha and beta) can have its own leader configuration:

| Option | Description |
|--------|-------------|
| `None` | No leader sequence |
| `"from_contig"` | Extract native leader from CellRanger FASTA |
| `"CD8A"`, `"CD28"`, etc. | Use a standard signal peptide |

**Default configuration**: CD28 on alpha chain, CD8A on beta chain (distinct leaders for identification).

## Available Leader Sequences

| Leader | Source | Species | Sequence |
|--------|--------|---------|----------|
| CD8A | CD8A signal peptide (UniProt P01732) | Human | MALPVTALLLPLALLLHAARP |
| CD28 | CD28 signal peptide (UniProt P10747) | Human | MLRLLLALNLFPSIQVTG |
| IgK | IgGκ light chain signal peptide | Mouse | METDTLLLWVLLLWVPGSTG |
| TRAC | TCR alpha constant signal peptide | Human | MAGTWLLLLLALGCPALPTG |
| TRBC | TCR beta constant signal peptide | Human | MGTSLLCWMALCLLGADHADG |

## Available Linkers

| Linker | Source | Sequence |
|--------|--------|----------|
| T2A | Thosea asigna virus | EGRGSLLTCGDVEENPGP |
| P2A | Porcine teschovirus-1 | GSGATNFSLLKQAGDVEENPGP |
| E2A | Equine rhinitis A virus | QCTNYALLKLAGDVESNPGP |
| F2A | Foot-and-mouth disease virus | VKQTLNFDLLKLAGDVESNPGP |

## Usage Examples

```python
from tcrsift import assemble_full_sequences

# Default: CD28 on alpha, CD8A on beta
assembled = assemble_full_sequences(clonotypes)

# Custom leaders
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader="CD8A",
    beta_leader="CD28",
)

# Leader only on beta chain (first in 2A construct)
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader=None,
    beta_leader="CD8A",
)

# Extract native leaders from contigs
assembled = assemble_full_sequences(
    clonotypes,
    contigs_dir="/path/to/contigs",
    alpha_leader="from_contig",
    beta_leader="from_contig",
)

# No leaders at all
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader=None,
    beta_leader=None,
)
```

## Constant-region construction

Canonical constant-region amino-acid sequences are packaged in
`tcrsift/refseqs/canonical_constants.fasta`. The default
`constant_source="ensembl"` is a backward-compatible name for this packaged
canonical path; no runtime Ensembl download occurs. Use
`constant_source="from-data"` to read explicit
`{chain}_constant_aa`/`{chain}_constant_nt` input columns.

With CellRanger contigs, nucleotide assembly is a hybrid:

1. Preserve the donor's J→C boundary and following complete codons while their
   translation agrees with the selected canonical protein.
2. Complete a 1- or 2-base terminal codon only when a compatible canonical
   amino-acid codon exists.
3. Use canonical codon-optimized nucleotides after contig coverage ends or
   translation diverges.

The protein remains canonical. `{chain}_constant_source` records how much
donor sequence was retained, and `qc_warnings` reports any divergence. Omit
`contigs_dir` to use canonical codon-optimized nucleotides for the entire
constant region.

## API Reference

::: tcrsift.assemble
    options:
      members:
        - assemble_full_sequences
        - translate_dna
        - find_longest_orf
        - parse_fasta
        - load_contigs
        - get_constant_region_sequences
        - validate_sequences
        - export_fasta
        - CODON_TABLE
        - LINKERS
        - DEFAULT_LEADERS
