# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Full-length TCR sequence assembly for TCRsift.

Builds complete TCR sequences including leader peptides and constant regions.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_directory_exists,
)

logger = logging.getLogger(__name__)


# Standard codon table
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}

# Self-cleaving 2A peptide linkers
# These are viral 2A sequences that cause ribosomal skipping during translation
LINKERS = {
    "T2A": {
        "dna": "GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGCCCG",
        "aa": "EGRGSLLTCGDVEENPGP",
        "source": "Thosea asigna virus",
    },
    "P2A": {
        "dna": "GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT",
        "aa": "GSGATNFSLLKQAGDVEENPGP",
        "source": "Porcine teschovirus-1",
    },
    "E2A": {
        "dna": "CAGTGTACTAATTATGCTCTCTTGAAATTGGCTGGAGATGTTGAGAGCAACCCAGGTCCC",
        "aa": "QCTNYALLKLAGDVESNPGP",
        "source": "Equine rhinitis A virus",
    },
    "F2A": {
        "dna": "GTGAAACAGACTTTGAATTTTGACCTTCTCAAGTTGGCGGGAGACGTGGAGTCCAACCCAGGGCCC",
        "aa": "VKQTLNFDLLKLAGDVESNPGP",
        "source": "Foot-and-mouth disease virus",
    },
}

# Backwards compatibility aliases
T2A_LINKER_DNA = LINKERS["T2A"]["dna"]
T2A_LINKER_AA = LINKERS["T2A"]["aa"]
P2A_LINKER_DNA = LINKERS["P2A"]["dna"]
P2A_LINKER_AA = LINKERS["P2A"]["aa"]

# Default leader/signal peptide sequences for TCR expression
# These can be used when contig FASTA files are not available
DEFAULT_LEADERS = {
    "CD8A": {
        "aa": "MALPVTALLLPLALLLHAARP",
        "dna": "ATGGCCCTGCCTGTGACAGCCCTGCTGCTGCCTCTGGCTCTGCTGCTGCATGCCGCTAGACCC",
        "source": "Human CD8A signal peptide (UniProt P01732)",
        "species": "human",
    },
    "CD28": {
        "aa": "MLRLLLALNLFPSIQVTG",
        "dna": "ATGCTCCGCCTGCTGCTGGCCCTGAACCTGTTCCCCAGCATCCAGGTGACCGGC",
        "source": "Human CD28 signal peptide (UniProt P10747)",
        "species": "human",
    },
    "IgK": {
        "aa": "METDTLLLWVLLLWVPGSTG",
        "dna": "ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGT",
        "source": "Murine IgGκ light chain signal peptide",
        "species": "mouse",
        "note": "Widely used for high secretion efficiency in mammalian expression",
    },
    "TRAC": {
        "aa": "MAGTWLLLLLALGCPALPTG",
        "dna": "ATGGCTGGCACCTGGCTGCTGCTGCTGCTGGCCCTGGGATGCCCAGCACTGCCCACAGGC",
        "source": "Human TRAC native signal peptide",
        "species": "human",
    },
    "TRBC": {
        "aa": "MGTSLLCWMALCLLGADHADG",
        "dna": "ATGGGCACCAGCCTGCTGTGCTGGATGGCCCTGTGCCTGCTGGGAGCAGACCACGCCGATGGC",
        "source": "Human TRBC native signal peptide",
        "species": "human",
    },
}

# Standard constant region endings for QC
CONSTANT_REGION_ENDINGS = {
    "TRAC": "LLMTLRLWSS",
    "TRBC1": "VKRKDF",
    "TRBC2": "VKRKDSRG",
}

# Minimum plausible length of the mature constant region per chain.
# TRAC mature ≈ 138 aa; TRBC1/TRBC2 mature ≈ 178 aa. The #66 bug produced
# 2–11 aa constants. Anything below these floors is a truncation.
CONSTANT_AA_FLOOR = {"alpha": 100, "beta": 150}

# In-cis J→C pairing for β chain: TRBJ1-* segments rearrange with TRBC1,
# TRBJ2-* with TRBC2. A row whose J-gene parity contradicts the called
# C-gene is either an annotation conflict or an assembly bug.
BETA_JC_PARITY = {"TRBJ1": "TRBC1", "TRBJ2": "TRBC2"}

# Allowed residue alphabet for assembled protein chains: the 20 standard
# AAs plus X (unknown — sometimes emitted by translation) and * (stop).
# Anything outside this set in a ``full_*_aa`` is corruption.
_VALID_AA_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYX*")


# Canonical human TCR constant-region amino acid sequences.
#
# These are the mature constant regions (post signal-cleavage), spanning
# the connecting peptide, transmembrane helix, and cytoplasmic tail —
# everything that goes downstream of the J-gene boundary in a cloned
# TCR construct. Hardcoding rather than fetching from Ensembl because:
#
#   * Pyensembl introduces an external dependency and a frame-handling
#     bug surface (#66) for sequences that don't actually change.
#   * The protein sequences are well-established (UniProt P01848 /
#     P01850 / A0A075B6Y0) and stable across releases.
#   * The Sarah TIL pipeline already substitutes canonical strings
#     downstream as a workaround; this brings the pipeline upstream.
#
# Each sequence ends with the canonical C-terminus stored in
# ``CONSTANT_REGION_ENDINGS`` — see ``validate_sequences``.
HUMAN_TRAC_AA = (
    "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFAC"
    "ANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS"
)
HUMAN_TRBC1_AA = (
    "EDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALND"
    "SRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGV"
    "LSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
)
HUMAN_TRBC2_AA = (
    "EDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALND"
    "SRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGV"
    "LSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG"
)

HUMAN_CONSTANT_REGIONS_AA: dict[str, str] = {
    "TRAC": HUMAN_TRAC_AA,
    "TRBC1": HUMAN_TRBC1_AA,
    "TRBC2": HUMAN_TRBC2_AA,
}


# Codon table for back-translating canonical AA constants to NT for
# downstream consumers that need DNA. Picks the most frequent human
# codon per residue (CCDS-weighted, GenScript table). This is
# representative — not authoritative for nucleotide-level work.
# Stop codon "*" maps to TAA.
HUMAN_PREFERRED_CODONS: dict[str, str] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA",
}


def back_translate(aa: str) -> str:
    """Reverse-translate a polypeptide to DNA via :data:`HUMAN_PREFERRED_CODONS`.

    The result is one of many valid back-translations; tcrsift uses it
    for NT columns that downstream pipelines expect (``*_constant_nt``)
    but doesn't claim it matches any particular Ensembl transcript.
    Unknown residues fall back to NNN.
    """
    fallback = "NNN"
    out = [HUMAN_PREFERRED_CODONS.get(r, fallback) for r in aa]
    return "".join(out)


def pick_canonical_constant(
    chain: str, c_gene: str | None = None, j_gene: str | None = None
) -> tuple[str, str]:
    """Choose the canonical constant-region AA for a chain.

    - Alpha → always TRAC.
    - Beta → J family decides when known (TRBJ1-* with TRBC1, TRBJ2-*
      with TRBC2), because the in-cis J→C pairing is locus-determined
      and reliable. ``c_gene`` is only consulted when ``j_gene`` is
      missing or its family can't be parsed; CellRanger's TRBC1/TRBC2
      call is unreliable (the two C genes share ~95% NT identity) and
      gets silently overridden when it contradicts J (#90). Aggregate
      override counts are surfaced at the assembly orchestrator level
      so a per-row override doesn't produce a log flood.

    Returns ``(gene_name, canonical_aa)``.
    """
    if chain == "alpha":
        return "TRAC", HUMAN_TRAC_AA

    j_family = _beta_j_family(j_gene)
    expected_c = BETA_JC_PARITY.get(j_family) if j_family else None
    if expected_c == "TRBC1":
        return "TRBC1", HUMAN_TRBC1_AA
    if expected_c == "TRBC2":
        return "TRBC2", HUMAN_TRBC2_AA

    # J family unknown — fall back to the (less reliable) c_gene call.
    c_base = _beta_c_base(c_gene)
    if c_base == "TRBC2":
        return "TRBC2", HUMAN_TRBC2_AA
    if c_base == "TRBC1":
        return "TRBC1", HUMAN_TRBC1_AA
    # Default to TRBC1 when nothing else is known.
    return "TRBC1", HUMAN_TRBC1_AA


def _beta_j_family(j_gene: str | None) -> str:
    """Parse the J-family prefix (``"TRBJ1"`` / ``"TRBJ2"``) from a
    raw J-gene call. Tolerates lowercase, whitespace, and allele
    suffixes (``"TRBJ1-1*02"``). Returns an empty string when the
    input doesn't look like a TRBJ call."""
    if not isinstance(j_gene, str):
        return ""
    token = j_gene.strip().upper().split("*", 1)[0]
    return token.split("-", 1)[0] if "-" in token else token


def _trbc_override_mask(df: pd.DataFrame) -> pd.Series:
    """Boolean Series identifying rows where ``beta_c_gene_canonical``
    differs from raw ``beta_c_gene`` after gene-family normalisation —
    the rows where :func:`pick_canonical_constant` overrode
    CellRanger's TRBC call.

    Returns an all-False Series (aligned to ``df.index``) when either
    column is missing. Used both by :func:`_count_trbc_overrides` and
    to emit the per-row ``beta_c_gene_overridden`` audit column.
    """
    if "beta_c_gene" not in df.columns or "beta_c_gene_canonical" not in df.columns:
        return pd.Series(False, index=df.index)
    raw_base = df["beta_c_gene"].apply(_beta_c_base)
    canon_base = df["beta_c_gene_canonical"].apply(_beta_c_base)
    return (raw_base != "") & (canon_base != "") & (raw_base != canon_base)


def _count_trbc_overrides(df: pd.DataFrame) -> int:
    """Count rows where :func:`pick_canonical_constant` overrode
    CellRanger's TRBC call. Thin wrapper over
    :func:`_trbc_override_mask` for callers that just want a tally."""
    return int(_trbc_override_mask(df).sum())


def _beta_c_base(c_gene: str | None) -> str:
    """Parse the C-gene base (``"TRBC1"`` / ``"TRBC2"``) from a raw
    C-gene call. Tolerates lowercase, whitespace, and allele suffixes
    (``"TRBC2*01"``). Returns an empty string when the input doesn't
    resolve to a known TRBC."""
    if not isinstance(c_gene, str):
        return ""
    token = c_gene.strip().upper().split("*", 1)[0]
    return token if token in {"TRBC1", "TRBC2"} else ""


def _resolve_c_gene(row, chain: str) -> str:
    """Coalesce ``{chain}_c_gene_canonical`` → ``{chain}_c_gene`` for a
    single row, treating NaN as missing.

    The naive ``row.get(canonical) or row.get(raw) or ""`` is wrong
    because ``float('nan')`` is truthy: when the canonical column is
    present but the cell is NaN, the ``or`` short-circuits to NaN and
    the raw column is never consulted. That's silent: validation
    downgrades structural failures to informational notes, and
    :func:`fix_jc_parity` silently skips the correction.
    """
    canonical = row.get(f"{chain}_c_gene_canonical")
    if isinstance(canonical, str) and canonical:
        return canonical
    raw = row.get(f"{chain}_c_gene")
    if isinstance(raw, str) and raw:
        return raw
    return ""


def verify_canonical_constant_start(
    observed: str, canonical_aa: str, min_match: int = 8
) -> bool:
    """Sanity-check the observed C-region start against the canonical.

    Returns True if at least ``min_match`` of the first ``len(observed)``
    residues match position-wise. Used to verify a CellRanger contig's
    post-J residues against the picked canonical before splicing the
    canonical full sequence into the assembled construct.
    """
    n = min(len(observed), len(canonical_aa))
    if n == 0:
        return False
    matches = sum(1 for a, b in zip(observed[:n], canonical_aa[:n]) if a == b)
    return matches >= min_match


def _describe_leader(param_value: str | None, resolved: str | dict | None) -> str:
    """Generate a description of a leader configuration for logging."""
    if param_value is None:
        return "None (no leader)"
    if resolved == "from_contig":
        return "from_contig (extract from FASTA)"
    if isinstance(resolved, dict):
        return f"{param_value.upper()} ({resolved['source']})"
    return str(param_value)


def translate_dna(dna_seq: str) -> tuple[str, str]:
    """
    Translate DNA sequence to amino acids.

    Returns
    -------
    tuple
        (amino_acid_sequence, ragged_3p_nucleotides)
    """
    seq_len = len(dna_seq)
    seq_len_trimmed = (seq_len // 3) * 3

    if seq_len != seq_len_trimmed:
        ragged_nt = dna_seq[seq_len_trimmed:]
        dna_seq = dna_seq[:seq_len_trimmed]
    else:
        ragged_nt = ""

    aa_seq = "".join([CODON_TABLE.get(dna_seq[i : i + 3], "X") for i in range(0, len(dna_seq), 3)])

    # Stop at first stop codon
    if "*" in aa_seq:
        ragged_nt = ""
        aa_seq = aa_seq[: aa_seq.index("*")]

    return aa_seq, ragged_nt


def find_longest_orf(dna_seq: str) -> tuple[str, int, str]:
    """
    Find and translate the longest open reading frame.

    Returns
    -------
    tuple
        (amino_acid_sequence, start_offset, ragged_3p_nucleotides)
    """
    start_positions = [i for i in range(len(dna_seq)) if dna_seq[i : i + 3] == "ATG"]

    longest_aa = ""
    longest_offset = 0
    longest_ragged = ""

    for start in start_positions:
        subseq = dna_seq[start:]
        aa, ragged = translate_dna(subseq)
        if len(aa) > len(longest_aa):
            longest_aa = aa
            longest_offset = start
            longest_ragged = ragged

    return longest_aa, longest_offset, longest_ragged


def parse_fasta(path: str | Path) -> dict[str, str]:
    """
    Parse a FASTA file.

    Returns
    -------
    dict
        Mapping from sequence ID to sequence
    """
    path = Path(path)
    results = {}
    curr_id = None
    lines = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id and lines:
                    results[curr_id] = "".join(lines)
                    lines = []
                curr_id = line[1:].split()[0]  # Take first word after >
            else:
                lines.append(line)

        # Don't forget last entry
        if curr_id and lines:
            results[curr_id] = "".join(lines)

    return results


def load_contigs(contig_dir: str | Path) -> dict[str, dict[str, str]]:
    """
    Load contig sequences from CellRanger output directories.

    Parameters
    ----------
    contig_dir : str or Path
        Directory containing sample subdirectories with FASTA files

    Returns
    -------
    dict
        Nested dict: sample -> contig_id -> sequence
    """
    contig_dir = Path(contig_dir)
    sample_contigs = {}

    # Look for FASTA files in subdirectories
    for fasta_path in contig_dir.rglob("*contig*.fasta"):
        sample_name = fasta_path.parent.name
        if sample_name not in sample_contigs:
            sample_contigs[sample_name] = {}
        sample_contigs[sample_name].update(parse_fasta(fasta_path))

    # Also check direct files
    for fasta_path in contig_dir.glob("*.fasta"):
        sample_name = fasta_path.stem.split("_")[0]
        if sample_name not in sample_contigs:
            sample_contigs[sample_name] = {}
        sample_contigs[sample_name].update(parse_fasta(fasta_path))

    logger.info(f"Loaded contigs from {len(sample_contigs)} samples")
    return sample_contigs


def get_constant_region_sequences() -> dict[str, str]:
    """
    Return human TCR constant-region CDS (DNA) via back-translation.

    Sources NT from :data:`HUMAN_CONSTANT_REGIONS_AA` and
    :func:`back_translate`. The earlier pyensembl-backed implementation
    of this function read the full mRNA at frame offset 2 and silently
    truncated TRAC / TRBC1 / TRBC2 to 2–11 residues for every
    assembled clonotype (#66). Hardcoding eliminates the frame bug,
    drops the pyensembl dependency, and matches the canonical
    sequences that downstream cloning constructs need (#67).

    Returns
    -------
    dict
        Gene name → CDS (DNA, ATG-prepended … stop). Sourced from the
        canonical AA via :func:`back_translate`.
    """
    out: dict[str, str] = {}
    for name, aa in HUMAN_CONSTANT_REGIONS_AA.items():
        # Construct ends at the canonical C-terminus; back-translate
        # the polypeptide and append a stop codon for completeness.
        out[name] = back_translate(aa) + HUMAN_PREFERRED_CODONS["*"]
    return out


def assemble_full_sequences(
    clonotypes: pd.DataFrame,
    contigs_dir: str | Path | None = None,
    alpha_leader: str | None = "CD28",
    beta_leader: str | None = "CD8A",
    include_constant: bool = True,
    constant_source: str = "ensembl",
    linker: str = "T2A",
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Assemble full-length TCR sequences.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame with VDJ sequences (from fwr1/cdr1/fwr2/cdr2/fwr3/cdr3/fwr4)
    contigs_dir : str or Path, optional
        Directory with CellRanger contig FASTA files. Required if alpha_leader or
        beta_leader is set to "from_contig".
    alpha_leader : str or None
        Leader sequence for alpha chain. Options:
        - None: No leader sequence
        - "from_contig": Extract native leader from contig FASTA (requires contigs_dir)
        - Key from DEFAULT_LEADERS: "CD8A", "CD28", "IgK", "TRAC", "TRBC"
        Default is "CD28" to provide distinct sequences from beta chain.
    beta_leader : str or None
        Leader sequence for beta chain. Same options as alpha_leader.
        Default is "CD8A" to provide distinct sequences from alpha chain.
    include_constant : bool
        Include constant region sequences.
    constant_source : str
        Source for constant regions:

        - ``"canonical"`` (default): splice in the hardcoded canonical
          human TRAC / TRBC1 / TRBC2 (:data:`HUMAN_CONSTANT_REGIONS_AA`).
          When CellRanger contigs are available, the canonical's first
          ~15 residues are verified against the observed contig start
          and a ``qc_warnings`` entry is emitted on mismatch.
        - ``"ensembl"``: back-compat alias for ``"canonical"``. The
          earlier pyensembl-backed path was removed in #66 — it read
          the full mRNA at the wrong frame offset and silently
          truncated constants to 2–11 residues.
        - ``"from-data"``: read ``{chain}_constant_aa`` /
          ``{chain}_constant_nt`` directly from the input frame.
    linker : str
        Linker sequence for single-chain constructs: "T2A", "P2A", "E2A", "F2A"
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with full sequences added

    Examples
    --------
    >>> # Default: CD28 on alpha, CD8A on beta (distinct leaders)
    >>> assembled = assemble_full_sequences(clonotypes)

    >>> # No leader sequences
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader=None)

    >>> # Leader only on beta chain (first in 2A construct)
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader="CD8A")

    >>> # Extract native leaders from contig FASTAs
    >>> assembled = assemble_full_sequences(
    ...     clonotypes,
    ...     contigs_dir="/path/to/contigs",
    ...     alpha_leader="from_contig",
    ...     beta_leader="from_contig",
    ... )
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_assembly=True)

    valid_constant_sources = ["canonical", "ensembl", "from-data"]
    if constant_source not in valid_constant_sources:
        raise TCRsiftValidationError(
            f"Invalid constant_source: '{constant_source}'",
            hint=f"Valid options are: {valid_constant_sources}",
        )

    # Validate and resolve leader options for each chain
    leader_config = {}
    for chain, leader_param in [("alpha", alpha_leader), ("beta", beta_leader)]:
        if leader_param is None:
            leader_config[chain] = None
        elif leader_param.lower() == "from_contig":
            if not contigs_dir:
                raise TCRsiftValidationError(
                    f"{chain}_leader='from_contig' requires contigs_dir to be specified",
                    hint="Provide contigs_dir with CellRanger FASTA files, or use a default leader like 'CD8A'",
                )
            leader_config[chain] = "from_contig"
        elif leader_param.upper() in DEFAULT_LEADERS:
            leader_config[chain] = DEFAULT_LEADERS[leader_param.upper()]
        else:
            raise TCRsiftValidationError(
                f"Unknown {chain}_leader: '{leader_param}'",
                hint=f"Valid options are: None, 'from_contig', or one of {list(DEFAULT_LEADERS.keys())}",
            )

    if verbose:
        alpha_desc = _describe_leader(alpha_leader, leader_config["alpha"])
        beta_desc = _describe_leader(beta_leader, leader_config["beta"])
        logger.info(f"Assembling full sequences for {len(clonotypes):,} clonotypes")
        logger.info(f"  Alpha leader: {alpha_desc}")
        logger.info(f"  Beta leader: {beta_desc}")
        logger.info(f"  Constant regions: {include_constant} (source: {constant_source})")
        logger.info(f"  Linker: {linker}")

    df = clonotypes.copy()

    # Load constant regions if needed. "canonical" / "ensembl" both
    # hit the same hardcoded-AA path now (#66, #67); "ensembl" is
    # retained as a back-compat alias and emits an info note.
    constant_seqs = {}
    if include_constant and constant_source in ("canonical", "ensembl"):
        if constant_source == "ensembl" and verbose:
            logger.info(
                "  constant_source='ensembl' is now an alias for 'canonical' "
                "(#66) — splicing canonical TRAC / TRBC1 / TRBC2 from "
                "HUMAN_CONSTANT_REGIONS_AA."
            )
        constant_seqs = get_constant_region_sequences()
        if verbose:
            logger.info(
                f"    Loaded {len(constant_seqs)} canonical constant region sequences"
            )

    # Warn if from-data constants requested but not present
    if include_constant and constant_source == "from-data":
        constant_cols = [
            "alpha_constant_aa",
            "alpha_constant_nt",
            "beta_constant_aa",
            "beta_constant_nt",
        ]
        if not any(col in df.columns for col in constant_cols):
            logger.warning(
                "  constant_source='from-data' but no constant region columns found in input. "
                "Constants will be omitted."
            )

    # Load contigs if needed for leader extraction
    sample_contigs = {}
    needs_contigs = (
        leader_config["alpha"] == "from_contig" or leader_config["beta"] == "from_contig"
    )
    if contigs_dir and needs_contigs:
        contigs_dir = validate_directory_exists(Path(contigs_dir), "contigs directory")
        if verbose:
            logger.info(f"  Loading contigs from {contigs_dir}...")
        sample_contigs = load_contigs(contigs_dir)
        if verbose:
            total_contigs = sum(len(c) for c in sample_contigs.values())
            logger.info(f"    Loaded {total_contigs:,} contigs from {len(sample_contigs)} samples")

    # Process each clonotype
    if verbose:
        logger.info("  Assembling sequences...")

    assembly_results = []

    # Create iterator with optional progress bar
    row_iter = df.iterrows()
    if show_progress:
        row_iter = tqdm(
            list(df.iterrows()),
            desc="Assembling sequences",
            unit="clone",
        )

    for idx, row in row_iter:
        result = _assemble_clone(
            row,
            sample_contigs,
            constant_seqs,
            leader_config,
            include_constant,
            constant_source,
        )
        assembly_results.append(result)

    # Add assembly columns to dataframe
    result_df = pd.DataFrame(assembly_results)
    for col in result_df.columns:
        df[col] = result_df[col].values

    # Add single-chain construct if requested
    if linker and "full_beta_aa" in df.columns and "full_alpha_aa" in df.columns:
        if verbose:
            logger.info(f"  Creating single-chain constructs with {linker} linker...")
        df = _add_single_chain(df, linker)

    # Summary
    if verbose:
        n_with_alpha = df["full_alpha_aa"].notna().sum() if "full_alpha_aa" in df.columns else 0
        n_with_beta = df["full_beta_aa"].notna().sum() if "full_beta_aa" in df.columns else 0
        n_single_chain = (
            df["single_chain_aa"].notna().sum() if "single_chain_aa" in df.columns else 0
        )
        logger.info("  Assembly complete:")
        logger.info(f"    With full alpha: {n_with_alpha:,}")
        logger.info(f"    With full beta: {n_with_beta:,}")
        logger.info(f"    Single-chain constructs: {n_single_chain:,}")

    # Per-row override audit column (#90). True where the canonical
    # picked by the J-family rule disagreed with CellRanger's raw
    # `beta_c_gene` call. Lets users filter / inspect post-hoc
    # without going through the log. Written even when no overrides
    # fired so the column is always present (defaults to False).
    if "beta_c_gene" in df.columns:
        df["beta_c_gene_overridden"] = _trbc_override_mask(df)

    # One aggregate line for J-family overrides of CellRanger TRBC
    # (#90). pick_canonical_constant is silent per-call to avoid a
    # log flood on cohorts where every other clone hits the case.
    # Gated by `verbose` so silent runs stay silent — the override
    # itself still happens, this is the audit trail. The per-row
    # `beta_c_gene_overridden` column above is the data-side audit
    # that survives the log.
    if verbose:
        n_overrides = (
            int(df["beta_c_gene_overridden"].sum())
            if "beta_c_gene_overridden" in df.columns
            else 0
        )
        if n_overrides:
            n_beta = df["beta_c_gene"].notna().sum() if "beta_c_gene" in df.columns else 0
            logger.warning(
                "  Overrode CellRanger TRBC call on %d / %d β clones to match "
                "J-family parity (TRBJ1→TRBC1, TRBJ2→TRBC2). CellRanger's "
                "TRBC1/2 discrimination is unreliable; the locus rule is "
                "authoritative.",
                n_overrides, n_beta,
            )

    return df


def _assemble_clone(
    row: pd.Series,
    sample_contigs: dict,
    constant_seqs: dict,
    leader_config: dict,
    include_constant: bool,
    constant_source: str,
) -> dict:
    """Assemble full sequence for a single clone."""
    result = {}

    # Try to get full sequence from VDJ columns if available
    for chain in ["alpha", "beta"]:
        vdj_col = f"VDJ_{chain}_aa"
        vdj_nt_col = f"VDJ_{chain}_nt"

        if vdj_col in row and pd.notna(row.get(vdj_col)):
            result[f"vdj_{chain}_aa"] = row[vdj_col]
        if vdj_nt_col in row and pd.notna(row.get(vdj_nt_col)):
            result[f"vdj_{chain}_nt"] = row[vdj_nt_col]

        # Get C gene and J gene for constant-region lookup. J gene is
        # needed by ``pick_canonical_constant`` to override
        # CellRanger's unreliable TRBC call via locus-determined
        # J-family parity (#90); without it the J-based override is
        # silently bypassed in the real assembly path.
        for gene in ("c_gene", "j_gene"):
            col = f"{chain}_{gene}"
            if col in row:
                result[col] = row[col]

    # Add leader sequences based on per-chain config
    for chain in ["alpha", "beta"]:
        chain_leader = leader_config.get(chain)
        if chain_leader is None:
            continue
        elif chain_leader == "from_contig":
            # Extract native leader from contigs
            _extract_leader_from_contigs_single(row, sample_contigs, result, chain)
        elif isinstance(chain_leader, dict):
            # Use specified default leader
            result[f"{chain}_leader_aa"] = chain_leader["aa"]
            result[f"{chain}_leader_nt"] = chain_leader["dna"]

    # Add constant regions
    if include_constant:
        if constant_source == "from-data":
            _add_constant_from_row(row, result)
        else:
            # ``canonical`` / ``ensembl`` (back-compat alias) flow:
            # pick the canonical TRAC / TRBC1 / TRBC2, splice the AA
            # in, and verify against the contig start when available.
            _add_constant_regions(result, constant_seqs, row=row, sample_contigs=sample_contigs)

    # Determine which chains have leaders for building full sequences
    include_alpha_leader = leader_config.get("alpha") is not None
    include_beta_leader = leader_config.get("beta") is not None

    # Build full sequences
    _build_full_sequences(result, include_alpha_leader, include_beta_leader, include_constant)

    return result


def _extract_leader_from_contigs_single(
    row: pd.Series,
    sample_contigs: dict,
    result: dict,
    chain: str,
):
    """Extract leader peptide from contig sequences for a single chain."""
    samples = str(row.get("samples", "")).split(";")

    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return

    contig_ids = str(row[contig_col]).split(";")
    vdj_aa = result.get(f"vdj_{chain}_aa", "")

    leader_counter = Counter()
    leader_dna_counter = Counter()

    for sample in samples:
        if sample not in sample_contigs:
            continue

        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue

            contig_seq = sample_contigs[sample][contig_id]
            translated, offset, ragged = find_longest_orf(contig_seq)

            if vdj_aa and vdj_aa in translated:
                parts = translated.split(vdj_aa)
                leader = parts[0]
                leader_counter[leader] += 1

                # Get leader DNA
                if offset is not None:
                    leader_dna = contig_seq[offset : offset + len(leader) * 3]
                    leader_dna_counter[leader_dna] += 1

    if leader_counter:
        result[f"{chain}_leader_aa"] = leader_counter.most_common(1)[0][0]
    if leader_dna_counter:
        result[f"{chain}_leader_nt"] = leader_dna_counter.most_common(1)[0][0]


def _add_constant_regions(
    result: dict,
    constant_seqs: dict,
    row: pd.Series | None = None,
    sample_contigs: dict | None = None,
):
    """Add constant-region sequences and verify against the observed
    CellRanger contig start where possible.

    Strategy: pick the canonical TRAC / TRBC1 / TRBC2 via
    :func:`pick_canonical_constant` using ``c_gene`` (with J-gene
    fallback for β), splice the canonical AA in as the constant, and —
    if contigs are available — pull the first ~15 residues observed
    after the VDJ in the contig and verify they match the canonical
    start. Mismatches surface as ``{chain}_constant_qc`` warnings on
    the result; the canonical is still used so downstream
    cloning-construct outputs are at least biologically valid.
    """
    for chain in ["alpha", "beta"]:
        c_gene = result.get(f"{chain}_c_gene")
        j_gene = result.get(f"{chain}_j_gene")
        canonical_name, canonical_aa = pick_canonical_constant(chain, c_gene, j_gene)

        # Splice in canonical AA + back-translated NT.
        result[f"{chain}_constant_aa"] = canonical_aa
        # `constant_seqs` already holds back-translated DNA for the
        # canonical AA; keep it pluggable in case a caller patches
        # the dict (the keys are gene names).
        result[f"{chain}_constant_nt"] = constant_seqs.get(
            canonical_name, back_translate(canonical_aa)
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name

        # Verify the observed contig C-region start against canonical.
        observed = (
            _extract_c_region_start_from_contig(
                row, sample_contigs, result.get(f"vdj_{chain}_aa", ""), chain
            )
            if (row is not None and sample_contigs)
            else None
        )
        if observed is None:
            result[f"{chain}_constant_source"] = f"canonical:{canonical_name}"
        elif verify_canonical_constant_start(observed, canonical_aa):
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (contig-verified)"
            )
        else:
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (UNVERIFIED — start mismatch)"
            )
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant start mismatch: observed "
                f"{observed!r} differs from canonical {canonical_name} "
                f"(expected start {canonical_aa[:len(observed)]!r}). "
                "Using canonical anyway; c_gene assignment may be wrong."
            )


def _extract_c_region_start_from_contig(
    row: pd.Series,
    sample_contigs: dict,
    vdj_aa: str,
    chain: str,
    n_aa: int = 15,
) -> str | None:
    """Pull the first ``n_aa`` residues observed after the VDJ in the
    CellRanger contig, if extractable.

    Mirrors the logic in :func:`_extract_leader_from_contigs_single` —
    longest-ORF translation, split on the VDJ AA, take what follows.
    Returns ``None`` when no contig contains the VDJ AA or the
    post-VDJ tail is shorter than ``n_aa``.
    """
    if not vdj_aa:
        return None
    samples = str(row.get("samples", "")).split(";")
    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return None

    contig_ids = str(row[contig_col]).split(";")
    starts = Counter()
    for sample in samples:
        if sample not in sample_contigs:
            continue
        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue
            contig_seq = sample_contigs[sample][contig_id]
            translated, _, _ = find_longest_orf(contig_seq)
            if vdj_aa in translated:
                _, after = translated.split(vdj_aa, 1)
                if len(after) >= n_aa:
                    starts[after[:n_aa]] += 1
    if starts:
        return starts.most_common(1)[0][0]
    return None


def _add_constant_from_row(row: pd.Series, result: dict):
    """Add constant region sequences directly from row columns if present.

    Also writes ``{chain}_c_gene_canonical`` so the from-data path
    surfaces the same audit column as the canonical path — without it,
    downstream code (``validate_sequences``, ``fix_jc_parity``,
    ``_count_trbc_overrides``) only sees the raw CellRanger
    ``{chain}_c_gene`` and the J-family override / parity checks
    silently degrade.
    """
    for chain in ["alpha", "beta"]:
        aa_col = f"{chain}_constant_aa"
        nt_col = f"{chain}_constant_nt"

        aa_val = row.get(aa_col) if aa_col in row else None
        nt_val = row.get(nt_col) if nt_col in row else None

        if pd.notna(aa_val):
            result[f"{chain}_constant_aa"] = aa_val
        if pd.notna(nt_val):
            result[f"{chain}_constant_nt"] = nt_val
            # If AA not provided, translate from NT
            if pd.isna(aa_val):
                const_aa, _ = translate_dna(str(nt_val))
                result[f"{chain}_constant_aa"] = const_aa

        # Record the canonical gene name so the from-data path produces
        # the same downstream audit column as the canonical path.
        canonical_name, _ = pick_canonical_constant(
            chain,
            c_gene=row.get(f"{chain}_c_gene"),
            j_gene=row.get(f"{chain}_j_gene"),
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name


def _build_full_sequences(
    result: dict,
    include_alpha_leader: bool,
    include_beta_leader: bool,
    include_constant: bool,
):
    """Build complete sequences from parts."""
    include_leader_map = {"alpha": include_alpha_leader, "beta": include_beta_leader}

    for chain in ["alpha", "beta"]:
        parts_aa = []
        parts_nt = []
        include_leader = include_leader_map[chain]

        if include_leader and f"{chain}_leader_aa" in result:
            parts_aa.append(result[f"{chain}_leader_aa"])
        if include_leader and f"{chain}_leader_nt" in result:
            parts_nt.append(result[f"{chain}_leader_nt"])

        if f"vdj_{chain}_aa" in result:
            parts_aa.append(result[f"vdj_{chain}_aa"])
        if f"vdj_{chain}_nt" in result:
            parts_nt.append(result[f"vdj_{chain}_nt"])

        if include_constant and f"{chain}_constant_aa" in result:
            parts_aa.append(result[f"{chain}_constant_aa"])
        if include_constant and f"{chain}_constant_nt" in result:
            parts_nt.append(result[f"{chain}_constant_nt"])

        if parts_aa:
            result[f"full_{chain}_aa"] = "".join(parts_aa)
        if parts_nt:
            result[f"full_{chain}_nt"] = "".join(parts_nt)


def _add_single_chain(df: pd.DataFrame, linker: str) -> pd.DataFrame:
    """Add single-chain construct (beta-linker-alpha)."""
    # Check if linker is a known 2A peptide
    if linker.upper() in LINKERS:
        linker_info = LINKERS[linker.upper()]
        linker_aa = linker_info["aa"]
        linker_nt = linker_info["dna"]
    else:
        # Custom linker sequence provided as amino acids
        linker_aa = linker
        linker_nt = ""

    # Remove stop codon from beta if present
    def strip_stop(seq):
        if not isinstance(seq, str):
            return None
        if seq.endswith("*"):
            return seq[:-1]
        return seq

    beta_aa = df["full_beta_aa"].apply(strip_stop)
    alpha_aa = df["full_alpha_aa"].where(df["full_alpha_aa"].apply(lambda x: isinstance(x, str)))

    single_chain = beta_aa.fillna("") + linker_aa + alpha_aa.fillna("")
    missing_mask = df["full_beta_aa"].isna() | df["full_alpha_aa"].isna()
    df["single_chain_aa"] = single_chain.where(~missing_mask, pd.NA)

    if "full_beta_nt" in df.columns and "full_alpha_nt" in df.columns and linker_nt:
        # Remove stop codon from beta DNA
        def strip_stop_codon_dna(seq):
            if not isinstance(seq, str):
                return None
            if len(seq) >= 3:
                last_codon = seq[-3:]
                if last_codon in {"TAA", "TAG", "TGA"}:
                    return seq[:-3]
            return seq

        beta_nt = df["full_beta_nt"].apply(strip_stop_codon_dna)
        alpha_nt = df["full_alpha_nt"].where(df["full_alpha_nt"].apply(lambda x: isinstance(x, str)))
        single_chain_nt = beta_nt.fillna("") + linker_nt + alpha_nt.fillna("")
        missing_nt_mask = df["full_beta_nt"].isna() | df["full_alpha_nt"].isna()
        df["single_chain_nt"] = single_chain_nt.where(~missing_nt_mask, pd.NA)

    df["linker"] = linker_aa

    return df


class ValidationMessage(str):
    """A validation/autocorrect message tagged with its clone index
    and severity.

    Subclasses :class:`str` so existing callers that treat the return
    value of :func:`validate_sequences` as a list of strings continue
    to work (substring tests, joins, length checks). New callers can
    filter on ``.severity`` and route on ``.idx`` instead of parsing
    the ``"Clone {idx}: ..."`` prefix back out of the string.

    Severities:

    - ``"load_bearing"`` — a structural failure (CDR3 missing, length
      out of window, premature stop, etc.) that justifies raising in
      strict mode.
    - ``"informational"`` — a "didn't have enough info to check" note;
      never raises.
    - ``"autocorrect"`` — emitted by :func:`fix_jc_parity` (and by
      :func:`validate_sequences` when ``fix=True``) to record an
      in-place correction that was applied.

    **Caveat:** any :class:`str` operation that returns a new value
    (``+``, slicing, ``.replace``, f-string interpolation,
    ``json.dumps``, etc.) returns a plain :class:`str` and drops the
    ``.idx`` / ``.severity`` attributes. Read the metadata directly
    off the original instance — don't massage the text first. ``copy``
    and ``pickle`` do preserve the attributes.
    """
    idx: object
    severity: str

    def __new__(cls, text: str, *, idx: object = None, severity: str = "load_bearing"):
        obj = super().__new__(cls, text)
        obj.idx = idx
        obj.severity = severity
        return obj


def fix_jc_parity(df: pd.DataFrame) -> list[ValidationMessage]:
    """Autocorrect β J→C parity mismatches in-place (#89).

    Locus structure guarantees TRBJ1-* rearranges with TRBC1 and TRBJ2-*
    with TRBC2. When the called C-gene disagrees, J is authoritative —
    overwrite ``beta_c_gene_canonical`` to match the J family.

    This only updates ``beta_c_gene_canonical`` (the column
    :func:`validate_sequences` and :func:`assemble_qc_report` consult
    first); the raw ``beta_c_gene`` is left alone as a record of what
    CellRanger originally called. ``full_beta_aa`` itself is not
    regenerated — use this when you have an already-assembled frame
    whose J/C parity annotation is wrong but the assembled sequence
    is correct (e.g. data from an older assembler that trusted
    CellRanger's TRBC call).

    Returns a list of :class:`ValidationMessage` with
    ``severity="autocorrect"``, one per row that was changed. The
    input DataFrame is mutated.
    """
    messages: list[ValidationMessage] = []
    if "beta_j_gene" not in df.columns:
        return messages
    # Make sure the canonical column exists with an object dtype that
    # accepts string writes (a float64 NaN column would emit a
    # FutureWarning on the first `df.at[idx, col] = "TRBC1"`).
    if "beta_c_gene_canonical" not in df.columns:
        df["beta_c_gene_canonical"] = df.get("beta_c_gene")
    df["beta_c_gene_canonical"] = df["beta_c_gene_canonical"].astype(object)

    # Vectorise the family lookup so we iterate only over the rows
    # that actually have a parity rule to enforce. On a 100k-row
    # cohort with ~0.1% override rate this skips ~99% of the loop
    # iterations.
    j_family_s = df["beta_j_gene"].apply(_beta_j_family)
    expected_s = j_family_s.map(BETA_JC_PARITY)
    candidates = expected_s.dropna()
    if candidates.empty:
        return messages

    for idx in candidates.index:
        row = df.loc[idx]
        current = _resolve_c_gene(row, "beta")
        current_base = _beta_c_base(current)
        expected_c = candidates.loc[idx]
        if current_base and current_base != expected_c:
            df.at[idx, "beta_c_gene_canonical"] = expected_c
            messages.append(
                ValidationMessage(
                    f"Clone {idx}: autocorrected β c_gene {current!r} → "
                    f"{expected_c} based on {row['beta_j_gene']} (J family rules)",
                    idx=idx,
                    severity="autocorrect",
                )
            )
    return messages


def validate_sequences(
    df: pd.DataFrame, strict: bool = False, fix: bool = False
) -> list[ValidationMessage]:
    """Validate assembled sequences end-to-end per row (#67, #68).

    Per-chain checks on each ``full_{chain}_aa`` row:

    - Length window (200–450 aa).
    - CDR3 substring present in the assembled sequence.
    - Constant region's canonical C-terminus (from
      :data:`CONSTANT_REGION_ENDINGS`).
    - Constant region's canonical start (first 8+ residues of
      :data:`HUMAN_CONSTANT_REGIONS_AA`).
    - Constant length floor (:data:`CONSTANT_AA_FLOOR`).
    - No premature stop codon (``*`` mid-chain).
    - Methionine start (when a leader was included).
    - Standard residue alphabet only (:data:`_VALID_AA_CHARS`).
    - Byte-for-byte equality of ``leader + vdj + constant`` and the
      assembled ``full_{chain}_aa``, when all three parts are present.

    Cross-row checks:

    - β-chain J→C in-cis parity: ``TRBJ1*`` pairs with TRBC1,
      ``TRBJ2*`` with TRBC2.
    - Single-chain integrity: linker present in ``single_chain_aa``.
    - Per-row ``qc_warnings`` (stashed by ``_add_constant_regions``)
      are surfaced.

    Returns a flat list of :class:`ValidationMessage` (a ``str``
    subclass with ``.idx`` and ``.severity`` attributes). When
    ``strict`` is True, raises :class:`TCRsiftValidationError` if any
    load-bearing checks fail. Informational notes ("didn't have enough
    info to check") are returned but never raise — distinguishing
    those from a silent pass was the gap that hid #66.

    When ``fix`` is True, autocorrect β J→C parity mismatches in-place
    via :func:`fix_jc_parity` before validating (#89). Autocorrection
    notes are prepended to the returned messages and the input frame
    is mutated *unconditionally* — the mutation survives even when a
    later check fails and ``strict=True`` raises.
    """
    load_bearing: list[ValidationMessage] = []
    informational: list[ValidationMessage] = []
    autocorrect_notes: list[ValidationMessage] = []

    def _lb(idx, text):
        load_bearing.append(
            ValidationMessage(f"Clone {idx}: {text}", idx=idx, severity="load_bearing")
        )

    def _info(idx, text):
        informational.append(
            ValidationMessage(f"Clone {idx}: {text}", idx=idx, severity="informational")
        )

    if fix:
        autocorrect_notes = fix_jc_parity(df)

    # 1. Per-chain, per-row shape and content checks.
    for chain in ["alpha", "beta"]:
        col = f"full_{chain}_aa"
        if col not in df.columns:
            continue
        for idx, row in df.iterrows():
            seq = row.get(col, "")
            if not seq or not isinstance(seq, str):
                continue

            if len(seq) < 200:
                _lb(idx, f"{chain} chain too short ({len(seq)} aa)")
            elif len(seq) > 450:
                _lb(idx, f"{chain} chain too long ({len(seq)} aa)")

            # Standard residue alphabet only.
            invalid = sorted(set(seq) - _VALID_AA_CHARS)
            if invalid:
                _lb(idx,
                    f"full_{chain}_aa has invalid residues {''.join(invalid)!r}")

            # Premature stop codon: ``*`` anywhere except the final position.
            body = seq.rstrip("*")
            if "*" in body:
                pos = body.index("*")
                _lb(idx,
                    f"full_{chain}_aa has premature stop "
                    f"at position {pos} of {len(seq)}")

            # Methionine start (only when leader was included).
            leader = row.get(f"{chain}_leader_aa")
            if isinstance(leader, str) and leader and not seq.startswith("M"):
                _lb(idx,
                    f"full_{chain}_aa doesn't start with M "
                    f"(starts with {seq[:5]!r}; leader present)")

            cdr3_col = f"CDR3_{chain}"
            if cdr3_col in row:
                cdr3 = row[cdr3_col]
                if isinstance(cdr3, str) and cdr3 and cdr3 not in seq:
                    _lb(idx,
                        f"CDR3_{chain}={cdr3!r} not found in full sequence")

            # Constant region length floor (separate from full-chain
            # length check — catches truncations that hid behind a
            # long leader).
            const = row.get(f"{chain}_constant_aa")
            if isinstance(const, str) and const:
                floor = CONSTANT_AA_FLOOR[chain]
                if len(const) < floor:
                    _lb(idx,
                        f"{chain}_constant_aa too short "
                        f"({len(const)} aa, floor {floor})")

            # Byte-for-byte: full == leader + vdj + constant when all
            # parts are available. Catches dropped/added residues
            # during assembly (the failure mode no other check finds).
            vdj = row.get(f"vdj_{chain}_aa")
            if (
                isinstance(leader, str) and leader
                and isinstance(vdj, str) and vdj
                and isinstance(const, str) and const
            ):
                expected = leader + vdj + const
                if expected != seq:
                    _lb(idx,
                        f"full_{chain}_aa != leader+vdj+constant "
                        f"(len {len(seq)} vs {len(expected)})")

            # Canonical-ending / canonical-start checks. Prefer
            # ``{chain}_c_gene_canonical`` (written by
            # ``_add_constant_regions``); fall back to the raw
            # CellRanger ``{chain}_c_gene`` call. If neither resolves
            # to a known canonical gene name, emit an informational
            # "did-not-check" note so the user knows validation was
            # incomplete.
            c_gene = _resolve_c_gene(row, chain)
            # Strip any allele suffix (e.g. "TRBC1*01" → "TRBC1").
            c_gene_base = c_gene.split("*")[0]
            if c_gene_base in CONSTANT_REGION_ENDINGS:
                expected_end = CONSTANT_REGION_ENDINGS[c_gene_base]
                if not seq.endswith(expected_end):
                    _lb(idx,
                        f"{chain} doesn't end with canonical "
                        f"{c_gene_base} C-terminus ({expected_end!r}); "
                        f"got {seq[-len(expected_end):]!r}")
                canonical_aa = HUMAN_CONSTANT_REGIONS_AA.get(c_gene_base)
                if canonical_aa and not verify_canonical_constant_start(
                    _expected_constant_start_from_full(seq, row, chain),
                    canonical_aa,
                    min_match=8,
                ):
                    _lb(idx,
                        f"{chain} constant start doesn't match canonical "
                        f"{c_gene_base} (expected start "
                        f"{canonical_aa[:15]!r})")
            else:
                _info(idx,
                    f"{chain} c_gene={c_gene!r} not in "
                    "CONSTANT_REGION_ENDINGS — canonical ending/start unverifiable")

    # 2. β-chain J→C in-cis parity disagreement (#67).
    if "beta_j_gene" in df.columns:
        for idx, row in df.iterrows():
            j_gene = row.get("beta_j_gene", "")
            c_gene = _resolve_c_gene(row, "beta")
            j_family = _beta_j_family(j_gene)
            c_base = _beta_c_base(c_gene)
            expected_c = BETA_JC_PARITY.get(j_family) if j_family else None
            if expected_c and c_base and c_base != expected_c:
                _lb(idx,
                    f"β J→C parity mismatch — "
                    f"{j_gene} (family {j_family}) should pair with "
                    f"{expected_c}, got {c_base}")

    # 3. Single-chain (β-linker-α) construct integrity. Three checks:
    #    a) the linker appears in ``single_chain_aa`` exactly once,
    #    b) ``single_chain_aa == β.rstrip('*') + linker + α`` byte-for-byte
    #       (catches dropped residues or wrong concatenation order),
    #    c) if the linker AA matches a known 2A peptide name, it should
    #       be the canonical sequence from :data:`LINKERS`.
    if "single_chain_aa" in df.columns and "linker" in df.columns:
        for idx, row in df.iterrows():
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            if not isinstance(sc, str) or not sc:
                continue
            if not isinstance(linker, str) or not linker:
                continue
            count = sc.count(linker)
            if count == 0:
                _lb(idx, f"single_chain_aa missing linker {linker!r}")
            elif count > 1:
                _lb(idx,
                    f"single_chain_aa contains linker "
                    f"{linker!r} {count} times (expected 1)")

            beta = row.get("full_beta_aa")
            alpha = row.get("full_alpha_aa")
            if (
                isinstance(beta, str) and beta
                and isinstance(alpha, str) and alpha
                and count == 1
            ):
                expected_sc = beta.rstrip("*") + linker + alpha
                if sc != expected_sc:
                    _lb(idx,
                        f"single_chain_aa != β+linker+α "
                        f"(len {len(sc)} vs {len(expected_sc)})")

    # 4. Surface any per-row qc_warnings the assembler stashed (e.g.
    # contig-vs-canonical start mismatch detected during assembly).
    if "qc_warnings" in df.columns:
        for idx, qcs in df["qc_warnings"].items():
            if isinstance(qcs, list):
                for msg in qcs:
                    _lb(idx, str(msg))

    all_messages = autocorrect_notes + load_bearing + informational
    if strict and load_bearing:
        raise TCRsiftValidationError(
            f"Sequence validation failed ({len(load_bearing)} load-bearing issues):\n  "
            + "\n  ".join(load_bearing[:10])
            + ("\n  ..." if len(load_bearing) > 10 else "")
        )
    return all_messages


@dataclass
class AssemblyQCCheck:
    """One QC check's tally across a frame.

    Attributes:
        name: machine-friendly identifier (e.g. ``"alpha_terminal_residue"``).
        label: human-readable label for the log/plot (e.g. ``"α terminal residue"``).
        chain: ``"alpha"``, ``"beta"``, ``"single_chain"``, or ``None`` for
            cross-chain checks.
        passed: row count that passed this check.
        total: row count this check ran against.
        median_value: optional median (used for length-style checks); unitless.
        unit: optional unit for ``median_value`` (e.g. ``"aa"``).
        examples: up to a few failure exemplars as ``"clone {idx}: ..."``.
    """

    name: str
    label: str
    chain: str | None
    passed: int
    total: int
    median_value: float | None = None
    unit: str = ""
    examples: list[str] = field(default_factory=list)

    @property
    def failed(self) -> int:
        return max(0, self.total - self.passed)

    @property
    def is_passing(self) -> bool:
        return self.failed == 0

    @property
    def pass_rate(self) -> float:
        return 1.0 if self.total == 0 else self.passed / self.total


@dataclass
class AssemblyQCReport:
    """Structured report covering every check ``assemble_qc_report``
    runs. Replaces the earlier string-only return (#67)."""

    n_rows: int
    checks: list[AssemblyQCCheck] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(c.is_passing for c in self.checks)

    def format_text(self) -> str:
        """Multi-line ``[assemble QC]`` summary suitable for logging."""
        lines = [f"[assemble QC] {self.n_rows} / {self.n_rows} rows checked"]
        for c in self.checks:
            mark = "✓" if c.is_passing else "✗"
            extras = []
            if c.median_value is not None:
                extras.append(
                    f"median {c.median_value:g}"
                    + (f" {c.unit}" if c.unit else "")
                )
            extra = f" ({', '.join(extras)})" if extras else ""
            lines.append(
                f"[assemble QC]   {mark} {c.label}: {c.passed}/{c.total} pass{extra}"
            )
        lines.append(
            f"[assemble QC] {'PASS' if self.passed else 'FAIL — see warnings above'}"
        )
        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Tabular form for programmatic consumption."""
        return pd.DataFrame(
            [
                {
                    "name": c.name,
                    "label": c.label,
                    "chain": c.chain,
                    "passed": c.passed,
                    "failed": c.failed,
                    "total": c.total,
                    "pass_rate": c.pass_rate,
                    "median": c.median_value,
                    "unit": c.unit,
                }
                for c in self.checks
            ]
        )

    def __str__(self) -> str:
        return self.format_text()


def build_assembly_qc_report(df: pd.DataFrame) -> AssemblyQCReport:
    """Build a structured :class:`AssemblyQCReport` (#67).

    Each per-check tally — terminal residue, length floors, premature
    stop, standard alphabet, β J→C parity, single-chain linker — is
    one :class:`AssemblyQCCheck` in the report's ``checks`` list.
    """
    n = len(df)
    report = AssemblyQCReport(n_rows=n)

    def _add(
        name: str,
        label: str,
        chain: str | None,
        passing: pd.Series,
        total: int,
        *,
        median_value: float | None = None,
        unit: str = "",
        df_for_examples: pd.DataFrame | None = None,
        failure_formatter=None,
    ) -> None:
        passed = int(passing.sum())
        examples: list[str] = []
        if failure_formatter is not None and df_for_examples is not None:
            failing_idx = passing.index[~passing]
            for idx in failing_idx[:3]:
                msg = failure_formatter(df_for_examples.loc[idx], idx)
                if msg:
                    examples.append(msg)
        report.checks.append(
            AssemblyQCCheck(
                name=name,
                label=label,
                chain=chain,
                passed=passed,
                total=total,
                median_value=median_value,
                unit=unit,
                examples=examples,
            )
        )

    for chain in ("alpha", "beta"):
        full_col = f"full_{chain}_aa"
        if full_col not in df.columns:
            continue
        full = df[full_col]
        present_mask = full.apply(lambda s: isinstance(s, str) and bool(s))
        total = int(present_mask.sum())
        if total == 0:
            continue
        sub = df[present_mask]
        sub_full = full[present_mask]
        chain_label = "α" if chain == "alpha" else "β"

        # Terminal residue (against canonical for the row's c_gene).
        def _ends_canonical(row, _chain=chain):
            s = row.get(f"full_{_chain}_aa", "")
            if not isinstance(s, str) or not s:
                return False
            c = _resolve_c_gene(row, _chain)
            base = c.split("*")[0]
            if base in CONSTANT_REGION_ENDINGS:
                return s.endswith(CONSTANT_REGION_ENDINGS[base])
            return any(
                s.endswith(end) for end in CONSTANT_REGION_ENDINGS.values()
            )

        end_pass = sub.apply(_ends_canonical, axis=1)
        _add(
            f"{chain}_terminal_residue",
            f"{chain_label} terminal residue",
            chain,
            end_pass,
            total,
            df_for_examples=sub,
            failure_formatter=lambda row, idx, c=chain: (
                f"clone {idx}: tail {row[f'full_{c}_aa'][-10:]!r}"
            ),
        )

        # Length floor (200 aa).
        lengths = sub_full.str.len()
        len_pass = lengths >= 200
        _add(
            f"{chain}_length_floor",
            f"{chain_label} length ≥ 200",
            chain,
            len_pass,
            total,
            median_value=float(lengths.median()),
            unit="aa",
        )

        # Constant length floor.
        const_col = f"{chain}_constant_aa"
        if const_col in df.columns:
            const_lens = sub[const_col].apply(
                lambda s: len(s) if isinstance(s, str) else 0
            )
            floor = CONSTANT_AA_FLOOR[chain]
            const_pass = const_lens >= floor
            present_const_lens = const_lens[const_lens > 0]
            cmedian = (
                float(present_const_lens.median())
                if len(present_const_lens)
                else 0.0
            )
            _add(
                f"{chain}_constant_floor",
                f"{chain_label} constant ≥ {floor}",
                chain,
                const_pass,
                total,
                median_value=cmedian,
                unit="aa",
            )

        # No premature stop.
        def _no_premature_stop(s):
            if not isinstance(s, str) or not s:
                return True
            return "*" not in s.rstrip("*")

        stop_pass = sub_full.apply(_no_premature_stop)
        _add(
            f"{chain}_no_premature_stop",
            f"{chain_label} no premature stop",
            chain,
            stop_pass,
            total,
        )

        # Standard alphabet.
        def _valid_alphabet(s):
            if not isinstance(s, str) or not s:
                return True
            return not (set(s) - _VALID_AA_CHARS)

        alpha_pass = sub_full.apply(_valid_alphabet)
        _add(
            f"{chain}_standard_alphabet",
            f"{chain_label} standard alphabet",
            chain,
            alpha_pass,
            total,
        )

    # β J→C parity.
    if "beta_j_gene" in df.columns:
        def _jc_parity_ok(row):
            j = row.get("beta_j_gene", "")
            c = _resolve_c_gene(row, "beta")
            fam = _beta_j_family(j)
            cbase = _beta_c_base(c)
            expected = BETA_JC_PARITY.get(fam) if fam else None
            if not expected or not cbase:
                return True
            return cbase == expected

        jc_pass = df.apply(_jc_parity_ok, axis=1)
        _add("beta_jc_parity", "β J→C in-cis parity", "beta", jc_pass, n)

    # Single-chain linker integrity.
    if "single_chain_aa" in df.columns and "linker" in df.columns:
        def _linker_ok(row):
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            if not isinstance(sc, str) or not sc:
                return True
            if not isinstance(linker, str) or not linker:
                return True
            return sc.count(linker) == 1

        link_pass = df.apply(_linker_ok, axis=1)
        _add(
            "single_chain_linker_present",
            "single_chain linker exactly once",
            "single_chain",
            link_pass,
            n,
        )

        # Byte-for-byte: single_chain == β.rstrip('*') + linker + α.
        def _single_chain_byte_exact(row):
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            beta = row.get("full_beta_aa")
            alpha = row.get("full_alpha_aa")
            if not all(
                isinstance(x, str) and x for x in (sc, linker, beta, alpha)
            ):
                return True
            return sc == beta.rstrip("*") + linker + alpha

        sc_byte_pass = df.apply(_single_chain_byte_exact, axis=1)
        _add(
            "single_chain_byte_exact",
            "single_chain == β + linker + α",
            "single_chain",
            sc_byte_pass,
            n,
        )

        # Canonical 2A peptide AA (when the linker name matches a known
        # 2A peptide, the AA should be the canonical sequence — guards
        # against custom-substituted linkers with the right NAME but
        # wrong AA).
        def _linker_canonical(row):
            linker = row.get("linker")
            if not isinstance(linker, str) or not linker:
                return True
            # The DataFrame's ``linker`` column stores the AA, not the
            # peptide name. Verify it matches one of the canonical 2A
            # AAs (or is a non-2A custom linker — which we accept).
            canonical_aas = {info["aa"] for info in LINKERS.values()}
            non_2a = linker not in canonical_aas and not any(
                linker.endswith(a[-6:]) for a in canonical_aas
            )
            return non_2a or linker in canonical_aas

        canon_pass = df.apply(_linker_canonical, axis=1)
        _add(
            "single_chain_2a_canonical",
            "2A peptide canonical AA",
            "single_chain",
            canon_pass,
            n,
        )

    return report


def assemble_qc_report(df: pd.DataFrame) -> str:
    """Aggregate QC summary across all rows (#67).

    Returns multi-line ``[assemble QC]`` text suitable for logging at
    the end of ``assemble_full_sequences``. Reports per-check tallies
    (pass / fail + median where relevant) plus a final PASS / FAIL
    banner. For programmatic consumption use
    :func:`build_assembly_qc_report`, which returns the same data as
    an :class:`AssemblyQCReport` object.
    """
    return build_assembly_qc_report(df).format_text()


def _expected_constant_start_from_full(
    full_aa: str, row: pd.Series, chain: str
) -> str:
    """Best-effort extract of the first ~15 aa of the constant region
    from a full assembled chain.

    Strategy: the full chain is leader + VDJ + constant. The VDJ AA
    is available as ``vdj_{chain}_aa``; split on it. If the VDJ isn't
    in the full sequence (which is itself a different failure caught
    elsewhere), fall back to slicing off the canonical-length leader.
    Returns up to 15 residues; an empty string when nothing useful
    can be extracted.
    """
    if not isinstance(full_aa, str) or not full_aa:
        return ""
    vdj = row.get(f"vdj_{chain}_aa")
    if isinstance(vdj, str) and vdj and vdj in full_aa:
        after = full_aa.split(vdj, 1)[1]
        return after[:15]
    return ""


def export_fasta(df: pd.DataFrame, output_path: str | Path, sequence_col: str = "single_chain_aa"):
    """
    Export sequences to FASTA format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequences
    output_path : str or Path
        Output file path
    sequence_col : str
        Column containing sequences to export
    """
    with open(output_path, "w") as f:
        for idx, row in df.iterrows():
            seq = row.get(sequence_col, "")
            if not seq:
                continue

            # Build header
            cdr3ab = row.get("CDR3ab", idx)
            cdr3a = row.get("CDR3_alpha", "")
            cdr3b = row.get("CDR3_beta", "")

            header = f">{cdr3ab} CDR3a={cdr3a} CDR3b={cdr3b}"
            f.write(f"{header}\n{seq}\n")

    logger.info(f"Exported {len(df)} sequences to {output_path}")
