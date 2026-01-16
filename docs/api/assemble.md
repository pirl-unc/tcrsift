# Assembly API

Module for full-length TCR sequence assembly.

## Overview

The assembly module builds full-length TCR sequences from CDR3 and V/J gene information. It supports:

- **Leader sequences**: Extract from contig FASTAs or use standard signal peptides (CD8A, CD28, IgK, TRAC, TRBC)
- **Constant regions**: Fetch from Ensembl (TRAC, TRBC1, TRBC2) or use built-in sequences
- **2A linkers**: Join alpha and beta chains with self-cleaving peptides (T2A, P2A, E2A, F2A)
- **Single-chain constructs**: Generate β-linker-α format for expression

## Available Linkers

| Linker | Source | Sequence |
|--------|--------|----------|
| T2A | Thosea asigna virus | EGRGSLLTCGDVEENPGP |
| P2A | Porcine teschovirus-1 | GSGATNFSLLKQAGDVEENPGP |
| E2A | Equine rhinitis A virus | QCTNYALLKLAGDVESNPGP |
| F2A | Foot-and-mouth disease virus | VKQTLNFDLLKLAGDVESNPGP |

## Default Leader Sequences

When contig FASTAs are not available, use `--default-leader` with one of:

| Leader | Source | Use Case |
|--------|--------|----------|
| CD8A | Human CD8A signal peptide | Common choice for TCR expression |
| CD28 | Human CD28 signal peptide | Alternative signal peptide |
| IgK | Human Ig kappa light chain | High secretion efficiency |
| TRAC | Native TCR alpha constant | Native-like expression |
| TRBC | Native TCR beta constant | Native-like expression |

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
