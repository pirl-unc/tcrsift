# Assembly API

Module for full-length TCR sequence assembly.

## Overview

The assembly module builds full-length TCR sequences from CDR3 and V/J gene information. It supports:

- **Leader sequences**: Per-chain configuration - extract from contig FASTAs or use standard signal peptides
- **Constant regions**: Fetch from Ensembl (TRAC, TRBC1, TRBC2) or use built-in sequences
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
