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
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
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

    aa_seq = "".join([
        CODON_TABLE.get(dna_seq[i:i+3], 'X')
        for i in range(0, len(dna_seq), 3)
    ])

    # Stop at first stop codon
    if "*" in aa_seq:
        ragged_nt = ""
        aa_seq = aa_seq[:aa_seq.index("*")]

    return aa_seq, ragged_nt


def find_longest_orf(dna_seq: str) -> tuple[str, int, str]:
    """
    Find and translate the longest open reading frame.

    Returns
    -------
    tuple
        (amino_acid_sequence, start_offset, ragged_3p_nucleotides)
    """
    start_positions = [i for i in range(len(dna_seq)) if dna_seq[i:i+3] == "ATG"]

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
    Get human TCR constant region sequences from Ensembl.

    Returns
    -------
    dict
        Gene name to coding sequence
    """
    try:
        from pyensembl import ensembl_grch38

        def find_stop_codon(seq, offset=0):
            for i in range(offset, len(seq), 3):
                codon = seq[i:i+3]
                if codon in {"TAA", "TAG", "TGA"}:
                    return i
            return None

        constants = {}

        # TRAC
        trac = ensembl_grch38.genes_by_name("TRAC")[0]
        trac_seq = trac.transcripts[0].sequence
        stop_idx = find_stop_codon(trac_seq, offset=2)
        if stop_idx:
            constants["TRAC"] = trac_seq[:stop_idx + 3]

        # TRBC1 and TRBC2
        for name in ["TRBC1", "TRBC2"]:
            gene = ensembl_grch38.genes_by_name(name)[0]
            seq = gene.transcripts[0].sequence
            stop_idx = find_stop_codon(seq, offset=2)
            if stop_idx:
                constants[name] = seq[:stop_idx + 3]

        return constants

    except ImportError:
        logger.warning("pyensembl not available, constant regions will not be included")
        return {}
    except Exception as e:
        logger.warning(f"Could not load constant regions from Ensembl: {e}")
        return {}


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
        Include constant region sequences (fetched from Ensembl or data)
    constant_source : str
        Source for constant regions: "ensembl" or "from-data"
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

    valid_constant_sources = ["ensembl", "from-data"]
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

    # Load constant regions if needed
    constant_seqs = {}
    if include_constant and constant_source == "ensembl":
        if verbose:
            logger.info("  Loading constant regions from Ensembl...")
        constant_seqs = get_constant_region_sequences()
        if not constant_seqs:
            logger.warning("  Could not load constant regions from Ensembl, will use sequences from data")
        elif verbose:
            logger.info(f"    Loaded {len(constant_seqs)} constant region sequences")

    # Load contigs if needed for leader extraction
    sample_contigs = {}
    needs_contigs = leader_config["alpha"] == "from_contig" or leader_config["beta"] == "from_contig"
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
        n_single_chain = df["single_chain_aa"].notna().sum() if "single_chain_aa" in df.columns else 0
        logger.info("  Assembly complete:")
        logger.info(f"    With full alpha: {n_with_alpha:,}")
        logger.info(f"    With full beta: {n_with_beta:,}")
        logger.info(f"    Single-chain constructs: {n_single_chain:,}")

    return df


def _assemble_clone(
    row: pd.Series,
    sample_contigs: dict,
    constant_seqs: dict,
    leader_config: dict,
    include_constant: bool,
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

        # Get C gene for constant region lookup
        c_gene_col = f"{chain}_c_gene"
        if c_gene_col in row:
            result[f"{chain}_c_gene"] = row[c_gene_col]

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
        _add_constant_regions(result, constant_seqs)

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
                    leader_dna = contig_seq[offset:offset + len(leader) * 3]
                    leader_dna_counter[leader_dna] += 1

    if leader_counter:
        result[f"{chain}_leader_aa"] = leader_counter.most_common(1)[0][0]
    if leader_dna_counter:
        result[f"{chain}_leader_nt"] = leader_dna_counter.most_common(1)[0][0]


def _add_constant_regions(result: dict, constant_seqs: dict):
    """Add constant region sequences."""
    for chain, c_gene_default in [("alpha", "TRAC"), ("beta", "TRBC1")]:
        c_gene = result.get(f"{chain}_c_gene", c_gene_default)
        if not c_gene:
            c_gene = c_gene_default

        if c_gene in constant_seqs:
            const_nt = constant_seqs[c_gene]
            const_aa, _ = translate_dna(const_nt)
            result[f"{chain}_constant_nt"] = const_nt
            result[f"{chain}_constant_aa"] = const_aa


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
        if seq and seq.endswith("*"):
            return seq[:-1]
        return seq

    df["single_chain_aa"] = (
        df["full_beta_aa"].apply(strip_stop) +
        linker_aa +
        df["full_alpha_aa"].fillna("")
    )

    if "full_beta_nt" in df.columns and "full_alpha_nt" in df.columns and linker_nt:
        # Remove stop codon from beta DNA
        def strip_stop_codon_dna(seq):
            if seq and len(seq) >= 3:
                last_codon = seq[-3:]
                if last_codon in {"TAA", "TAG", "TGA"}:
                    return seq[:-3]
            return seq

        df["single_chain_nt"] = (
            df["full_beta_nt"].apply(strip_stop_codon_dna) +
            linker_nt +
            df["full_alpha_nt"].fillna("")
        )

    df["linker"] = linker_aa

    return df


def validate_sequences(df: pd.DataFrame) -> list[str]:
    """
    Validate assembled sequences.

    Returns
    -------
    list
        List of warning messages
    """
    warnings = []

    # Check sequence lengths
    for chain in ["alpha", "beta"]:
        col = f"full_{chain}_aa"
        if col not in df.columns:
            continue

        for idx, row in df.iterrows():
            seq = row.get(col, "")
            if not seq:
                continue

            if len(seq) < 200:
                warnings.append(f"Clone {idx}: {chain} chain too short ({len(seq)} aa)")
            if len(seq) > 450:
                warnings.append(f"Clone {idx}: {chain} chain too long ({len(seq)} aa)")

            # Check CDR3 is present
            cdr3_col = f"CDR3_{chain}"
            if cdr3_col in row:
                cdr3 = row[cdr3_col]
                if cdr3 and cdr3 not in seq:
                    warnings.append(f"Clone {idx}: CDR3_{chain} not found in full sequence")

    # Check constant region endings
    for idx, row in df.iterrows():
        for chain in ["alpha", "beta"]:
            c_gene = row.get(f"{chain}_c_gene", "")
            full_seq = row.get(f"full_{chain}_aa", "")

            if c_gene and full_seq and c_gene in CONSTANT_REGION_ENDINGS:
                expected_end = CONSTANT_REGION_ENDINGS[c_gene]
                if not full_seq.endswith(expected_end):
                    warnings.append(
                        f"Clone {idx}: {chain} constant region doesn't end with expected "
                        f"sequence for {c_gene}"
                    )

    return warnings


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
            clone_id = row.get("clone_id", idx)
            cdr3a = row.get("CDR3_alpha", "")
            cdr3b = row.get("CDR3_beta", "")

            header = f">{clone_id} CDR3a={cdr3a} CDR3b={cdr3b}"
            f.write(f"{header}\n{seq}\n")

    logger.info(f"Exported {len(df)} sequences to {output_path}")
