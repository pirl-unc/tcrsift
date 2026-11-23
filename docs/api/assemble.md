# Assembly API

`tcrsift.assemble` constructs full-length α, β, and synthesis-ready
single-chain (`β-2A-α`) protein and DNA sequences from clonotype tables.

## Overview

By default, input tables carry full-length V(D)J variable-region sequences
(`VDJ_alpha_aa`, `VDJ_beta_aa`) from `aggregate_clonotypes`. When working with
5'-truncated contigs (`partial_VDJ_{chain}_aa`/`_nt`, where sequencing begins
downstream of the leader start codon) or bare `(v_gene, CDR3, j_gene)` tables,
pass `germline_vdj_fallback=True` (`--germline-vdj-fallback` on the CLI, after
`tcrsift data download --group germline` or `--auto-download-germline`):

- **5'-partial contig rescue (`contig_5prime_germline_rescued`):** prepends only
  the missing N-terminal FR1/CDR1 segment from the best-matching IMGT `V-REGION`
  allele while preserving all sequenced donor framework SNPs, CDR3 codons, and
  FR4 codons verbatim.
- **Bare V+CDR3+J reconstruction (`germline_vj_reconstructed`):** splices IMGT
  germline `V-REGION` (`FR1–FR3` up to conserved `Cys 104`) + `CDR3` + IMGT
  germline `J-REGION` (`FR4` past `[FW]GXG`).
- Every germline-imputed chain records `{chain}_vdj_source`,
  `{chain}_vdj_germline_allele`, and `qc_warnings` (`{chain}_vdj_germline_reconstructed`)
  and sets `construct_contig_verified = False` so germline imputation is never
  conflated with a 5'-to-3' contig-verified sequence.

Assembly configures four components:

- **Signal peptides (leaders):** standard human or murine leader peptides per
  chain, native V-gene leaders extracted from Cell Ranger FASTA contigs, or no
  leader.
- **Constant regions:** packaged canonical human (`TRAC`, `TRBC1`, `TRBC2`) and
  mouse (`Trac`, `Trbc1`, `Trbc2`) references (including chimeric human-VDJ +
  mouse-constant pairing via `constant_species="mouse"`), or explicit
  `{chain}_constant_aa`/`_nt` columns from input data.
- **Self-cleaving 2A linkers:** `T2A`, `P2A`, `E2A`, or `F2A`.
- **Single-chain orientation:** `β-linker-α` constructs with validation for
  start/stop codons, reading frame, J→C junction seam, and repeated k-mers.

## Signal-peptide options

Each chain takes its own leader setting (defaults: `CD28` on α and `CD8A` on β
so chains remain distinguishable):

| Value | Behavior |
| --- | --- |
| `"CD8A"`, `"CD28"`, `"mCD8B"`, `"mCD28"`, `"H2Kb"`, `"IgK"`, `"TRAC"`, `"TRBC"` | Predefined signal peptide (table below) |
| `"from_contig"` | Native V-gene leader extracted from Cell Ranger contig FASTAs (germline-anchored against IMGT L-PART1+L-PART2) |
| `None` | Omit signal peptide |

| Name | Origin | Species | Amino-acid sequence |
| --- | --- | --- | --- |
| `CD8A` | `CD8A` signal peptide (UniProt P01732) | Human | `MALPVTALLLPLALLLHAARP` |
| `CD28` | `CD28` signal peptide (UniProt P10747) | Human | `MLRLLLALNLFPSIQVTG` |
| `mCD8B` | `Cd8b1` signal peptide (UniProt P10300) | Mouse | `MQPWLWLVFSMKLAALWSSSA` |
| `mCD28` | `Cd28` signal peptide (UniProt P31041) | Mouse | `MTLRLLFLALNFFSVQVT` |
| `H2Kb` | `H2-K1` (H-2Kb) signal peptide (UniProt P01901) | Mouse | `MVPCTLLLLLAAALAPTQTRA` |
| `IgK` | IgGκ light-chain signal peptide | Mouse | `METDTLLLWVLLLWVPGSTG` |
| `TRAC` | Synthetic TCR α-chain signal peptide | Human | `MAGTWLLLLLALGCPALPTG` |
| `TRBC` | Synthetic TCR β-chain signal peptide | Human | `MGTSLLCWMALCLLGADHADG` |

## 2A linkers

| Name | Virus | Amino-acid sequence |
| --- | --- | --- |
| `T2A` | *Thosea asigna* virus | `EGRGSLLTCGDVEENPGP` |
| `P2A` | Porcine teschovirus-1 | `GSGATNFSLLKQAGDVEENPGP` |
| `E2A` | Equine rhinitis A virus | `QCTNYALLKLAGDVESNPGP` |
| `F2A` | Foot-and-mouth disease virus | `VKQTLNFDLLKLAGDVESNPGP` |

## Usage examples

```python
from tcrsift import assemble_full_sequences

# Default: CD28 leader on alpha, CD8A leader on beta
assembled = assemble_full_sequences(clonotypes)

# Native V-gene leaders from Cell Ranger contig FASTAs
assembled = assemble_full_sequences(
    clonotypes,
    contigs_dir="/path/to/contigs",
    alpha_leader="from_contig",
    beta_leader="from_contig",
)

# Leader on beta chain only (upstream chain in a β-2A-α construct)
assembled = assemble_full_sequences(
    clonotypes,
    alpha_leader=None,
    beta_leader="CD8A",
)
```

## Constant-region construction

Canonical human (`tcrsift/refseqs/canonical_constants.fasta`) and mouse
(`tcrsift/refseqs/canonical_constants_mouse.fasta`) `TRAC`, `TRBC1`, and `TRBC2`
sequences are bundled in the package. Pass `species="mouse"` for murine TCRs, or
`species="human", constant_species="mouse"` for murinized chimeric TCR-T
constructs (`codon_species` defaults to `constant_species` and can be set
independently, e.g. `codon_species="human"` when expressing a murinized TCR in
human T cells). The default `constant_source="ensembl"` is retained as a
backward-compatible alias for `"canonical"` and performs no network or Ensembl
lookup. Pass `constant_source="from-data"` to read explicit
`{chain}_constant_aa` and `{chain}_constant_nt` columns instead.

When `contigs_dir` is provided, nucleotide assembly blends donor contig 3′
reads into the canonical constant region:

1. Retain donor nucleotides across the J→C junction and subsequent complete
   codons while their translation matches the canonical constant-region protein.
2. Complete a trailing 1- or 2-nt partial contig codon when a compatible codon
   exists for the canonical amino acid.
3. Fill the remainder of the constant region with codon-optimized canonical
   nucleotides once contig coverage ends or diverges.

The translated constant-region protein stays canonical; `{chain}_constant_source`
records how much donor sequence was kept, and `qc_warnings` flags any
divergence. Without `contigs_dir`, codon-optimized canonical sequence is used
across the entire constant region.

## API reference

::: tcrsift.assemble
    options:
      members:
        - assemble_full_sequences
        - reconstruct_germline_vdj
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
