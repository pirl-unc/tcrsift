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

# T2A self-cleaving peptide linker
T2A_LINKER_DNA = "AGAGCCGAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGGCCG"
T2A_LINKER_AA = "RAEGRGSLLTCGDVEENPGP"

# Standard constant region endings for QC
CONSTANT_REGION_ENDINGS = {
    "TRAC": "LLMTLRLWSS",
    "TRBC1": "VKRKDF",
    "TRBC2": "VKRKDSRG",
}


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
    include_leader: bool = True,
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
        Clonotype DataFrame
    contigs_dir : str or Path, optional
        Directory with CellRanger contig FASTA files
    include_leader : bool
        Include leader peptide sequences
    include_constant : bool
        Include constant region sequences
    constant_source : str
        Source for constant regions: "ensembl" or "from-data"
    linker : str
        Linker sequence for single-chain constructs
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with full sequences added
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_assembly=True)

    valid_constant_sources = ["ensembl", "from-data"]
    if constant_source not in valid_constant_sources:
        raise TCRsiftValidationError(
            f"Invalid constant_source: '{constant_source}'",
            hint=f"Valid options are: {valid_constant_sources}",
        )

    if verbose:
        logger.info(f"Assembling full sequences for {len(clonotypes):,} clonotypes")
        logger.info(f"  Options: include_leader={include_leader}, include_constant={include_constant}, linker={linker}")

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

    # Load contigs if provided
    sample_contigs = {}
    if contigs_dir:
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
            include_leader,
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
    include_leader: bool,
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

    # If we have contigs, try to extract leader and constant
    if sample_contigs and include_leader:
        _extract_leader_from_contigs(row, sample_contigs, result)

    # Add constant regions
    if include_constant:
        _add_constant_regions(result, constant_seqs)

    # Build full sequences
    _build_full_sequences(result, include_leader, include_constant)

    return result


def _extract_leader_from_contigs(
    row: pd.Series,
    sample_contigs: dict,
    result: dict,
):
    """Extract leader peptide from contig sequences."""
    samples = str(row.get("samples", "")).split(";")

    for chain, prefix in [("alpha", "TRA"), ("beta", "TRB")]:
        contig_col = f"{chain}_contig_ids"
        if contig_col not in row or pd.isna(row[contig_col]):
            continue

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


def _build_full_sequences(result: dict, include_leader: bool, include_constant: bool):
    """Build complete sequences from parts."""
    for chain in ["alpha", "beta"]:
        parts_aa = []
        parts_nt = []

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
    if linker == "T2A":
        linker_aa = T2A_LINKER_AA
        linker_nt = T2A_LINKER_DNA
    else:
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
