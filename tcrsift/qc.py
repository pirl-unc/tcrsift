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
Quality control validation for TCRsift.

Provides validation checks for TCR sequences including:
- Start/stop codon validation
- Frame validation
- Repeated k-mer detection
- CDR3 length bounds
- Chain length bounds
- Read/UMI count thresholds
"""
from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field

import pandas as pd

logger = logging.getLogger(__name__)


# Default QC parameters
DEFAULT_MIN_CDR3_LENGTH = 5
DEFAULT_MAX_CDR3_LENGTH = 40
DEFAULT_MIN_CHAIN_LENGTH = 80  # Minimum amino acids for a chain
DEFAULT_MAX_CHAIN_LENGTH = 450
DEFAULT_MIN_READS = 10
DEFAULT_MIN_UMIS = 1
DEFAULT_KMER_SIZE = 10
DEFAULT_MIN_KMER_REPEAT = 2


@dataclass
class QCResult:
    """Result of a QC check."""
    passed: bool
    check_name: str
    message: str
    details: dict = field(default_factory=dict)


@dataclass
class QCReport:
    """Complete QC report for a dataset."""
    results: list[QCResult] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        """Check if all QC checks passed."""
        return all(r.passed for r in self.results) and len(self.errors) == 0

    @property
    def num_passed(self) -> int:
        """Number of passed checks."""
        return sum(1 for r in self.results if r.passed)

    @property
    def num_failed(self) -> int:
        """Number of failed checks."""
        return sum(1 for r in self.results if not r.passed)

    def add_result(self, result: QCResult):
        """Add a QC result."""
        self.results.append(result)
        if not result.passed:
            self.errors.append(result.message)

    def add_warning(self, message: str):
        """Add a warning message."""
        self.warnings.append(message)

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            f"QC Report: {'PASSED' if self.passed else 'FAILED'}",
            f"  Checks passed: {self.num_passed}/{len(self.results)}",
        ]
        if self.warnings:
            lines.append(f"  Warnings: {len(self.warnings)}")
        if self.errors:
            lines.append(f"  Errors: {len(self.errors)}")
            for err in self.errors[:5]:
                lines.append(f"    - {err}")
            if len(self.errors) > 5:
                lines.append(f"    ... and {len(self.errors) - 5} more")
        return "\n".join(lines)


def find_repeated_kmers(
    sequence: str,
    k: int = DEFAULT_KMER_SIZE,
    min_repeat: int = DEFAULT_MIN_KMER_REPEAT,
) -> list[tuple[str, int]]:
    """
    Find k-mers that appear multiple times in a sequence.

    Parameters
    ----------
    sequence : str
        Amino acid or nucleotide sequence
    k : int
        K-mer size
    min_repeat : int
        Minimum number of occurrences to report

    Returns
    -------
    list[tuple[str, int]]
        List of (kmer, count) tuples for repeated k-mers
    """
    if len(sequence) < k:
        return []

    kmer_counts: Counter = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        kmer_counts[kmer] += 1

    return [(kmer, count) for kmer, count in kmer_counts.most_common() if count >= min_repeat]


def check_start_codon(dna_sequence: str) -> bool:
    """Check if DNA sequence starts with ATG."""
    return dna_sequence.upper().startswith("ATG") if dna_sequence else False


def check_stop_codon(dna_sequence: str) -> bool:
    """Check if DNA sequence ends with a stop codon."""
    if not dna_sequence or len(dna_sequence) < 3:
        return False
    stop_codons = {"TAA", "TAG", "TGA"}
    return dna_sequence.upper()[-3:] in stop_codons


def check_reading_frame(dna_sequence: str) -> bool:
    """Check if DNA sequence length is a multiple of 3."""
    return len(dna_sequence) % 3 == 0 if dna_sequence else False


def check_cdr3_length(
    cdr3: str,
    min_length: int = DEFAULT_MIN_CDR3_LENGTH,
    max_length: int = DEFAULT_MAX_CDR3_LENGTH,
) -> bool:
    """Check if CDR3 length is within acceptable bounds."""
    if not cdr3:
        return False
    return min_length <= len(cdr3) <= max_length


def check_chain_length(
    chain_sequence: str,
    min_length: int = DEFAULT_MIN_CHAIN_LENGTH,
    max_length: int = DEFAULT_MAX_CHAIN_LENGTH,
) -> bool:
    """Check if chain length is within acceptable bounds."""
    if not chain_sequence:
        return False
    return min_length <= len(chain_sequence) <= max_length


def validate_sequence(
    sequence: str,
    cdr3: str | None = None,
    is_dna: bool = False,
    check_codons: bool = True,
    min_cdr3_length: int = DEFAULT_MIN_CDR3_LENGTH,
    max_cdr3_length: int = DEFAULT_MAX_CDR3_LENGTH,
    min_chain_length: int = DEFAULT_MIN_CHAIN_LENGTH,
    max_chain_length: int = DEFAULT_MAX_CHAIN_LENGTH,
    kmer_size: int = DEFAULT_KMER_SIZE,
    min_kmer_repeat: int = DEFAULT_MIN_KMER_REPEAT,
) -> QCReport:
    """
    Validate a single TCR sequence.

    Parameters
    ----------
    sequence : str
        The sequence to validate (DNA or amino acid)
    cdr3 : str, optional
        CDR3 sequence to check for presence and length
    is_dna : bool
        Whether sequence is DNA (vs amino acid)
    check_codons : bool
        Whether to check start/stop codons (only for DNA)
    min_cdr3_length, max_cdr3_length : int
        CDR3 length bounds
    min_chain_length, max_chain_length : int
        Chain length bounds
    kmer_size : int
        K-mer size for repeat detection
    min_kmer_repeat : int
        Minimum k-mer occurrences to flag

    Returns
    -------
    QCReport
        QC report with all check results
    """
    report = QCReport()

    if not sequence:
        report.add_result(QCResult(False, "sequence_present", "Sequence is empty or None"))
        return report

    # Length check (for amino acids, or DNA/3 for nucleotides)
    effective_length = len(sequence) // 3 if is_dna else len(sequence)
    length_ok = check_chain_length(sequence if not is_dna else "X" * effective_length,
                                    min_chain_length, max_chain_length)
    report.add_result(QCResult(
        length_ok,
        "chain_length",
        f"Chain length {effective_length}aa {'OK' if length_ok else f'out of bounds [{min_chain_length}, {max_chain_length}]'}",
        {"length": effective_length}
    ))

    # DNA-specific checks
    if is_dna:
        # Reading frame
        frame_ok = check_reading_frame(sequence)
        report.add_result(QCResult(
            frame_ok,
            "reading_frame",
            f"Reading frame {'OK' if frame_ok else 'NOT a multiple of 3'}",
            {"length": len(sequence), "remainder": len(sequence) % 3}
        ))

        if check_codons:
            # Start codon
            start_ok = check_start_codon(sequence)
            report.add_result(QCResult(
                start_ok,
                "start_codon",
                f"Start codon {'ATG found' if start_ok else 'NOT found'}",
                {"first_codon": sequence[:3] if len(sequence) >= 3 else ""}
            ))

            # Stop codon
            stop_ok = check_stop_codon(sequence)
            report.add_result(QCResult(
                stop_ok,
                "stop_codon",
                f"Stop codon {'found' if stop_ok else 'NOT found'}",
                {"last_codon": sequence[-3:] if len(sequence) >= 3 else ""}
            ))

    # CDR3 checks
    if cdr3:
        # CDR3 length
        cdr3_length_ok = check_cdr3_length(cdr3, min_cdr3_length, max_cdr3_length)
        report.add_result(QCResult(
            cdr3_length_ok,
            "cdr3_length",
            f"CDR3 length {len(cdr3)}aa {'OK' if cdr3_length_ok else f'out of bounds [{min_cdr3_length}, {max_cdr3_length}]'}",
            {"cdr3_length": len(cdr3)}
        ))

        # CDR3 presence in sequence
        check_seq = sequence if not is_dna else None  # Can't check AA in DNA directly
        if check_seq:
            cdr3_present = cdr3 in check_seq
            cdr3_count = check_seq.count(cdr3)
            report.add_result(QCResult(
                cdr3_present and cdr3_count == 1,
                "cdr3_present",
                f"CDR3 {'present once' if cdr3_present and cdr3_count == 1 else 'missing or duplicated'} in sequence",
                {"present": cdr3_present, "count": cdr3_count}
            ))

    # K-mer repeat check (only for amino acid sequences, or translate DNA first)
    check_seq_for_kmers = sequence if not is_dna else None
    if check_seq_for_kmers and len(check_seq_for_kmers) >= kmer_size:
        repeated_kmers = find_repeated_kmers(check_seq_for_kmers, kmer_size, min_kmer_repeat)
        kmer_ok = len(repeated_kmers) == 0
        report.add_result(QCResult(
            kmer_ok,
            "repeated_kmers",
            f"{'No' if kmer_ok else len(repeated_kmers)} repeated {kmer_size}-mers found",
            {"repeated_kmers": repeated_kmers[:5]}  # Only store first 5
        ))
        if not kmer_ok:
            for kmer, count in repeated_kmers[:3]:
                report.add_warning(f"Repeated k-mer (n={count}): {kmer}")

    return report


def validate_clonotypes(
    clonotypes: pd.DataFrame,
    min_reads: int = DEFAULT_MIN_READS,
    min_umis: int = DEFAULT_MIN_UMIS,
    min_cdr3_length: int = DEFAULT_MIN_CDR3_LENGTH,
    max_cdr3_length: int = DEFAULT_MAX_CDR3_LENGTH,
    min_chain_length: int = DEFAULT_MIN_CHAIN_LENGTH,
    max_chain_length: int = DEFAULT_MAX_CHAIN_LENGTH,
    check_full_sequences: bool = False,
) -> tuple[pd.DataFrame, QCReport]:
    """
    Validate a clonotype DataFrame and add QC columns.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        DataFrame with clonotype information
    min_reads : int
        Minimum read count threshold
    min_umis : int
        Minimum UMI count threshold
    min_cdr3_length, max_cdr3_length : int
        CDR3 length bounds
    min_chain_length, max_chain_length : int
        Chain length bounds (for full sequences)
    check_full_sequences : bool
        Whether to validate full assembled sequences

    Returns
    -------
    tuple[pd.DataFrame, QCReport]
        DataFrame with QC columns added, and overall QC report
    """
    df = clonotypes.copy()
    report = QCReport()

    n_rows = len(df)
    logger.info(f"Validating {n_rows} clonotypes")

    # Initialize QC columns
    df["qc_pass"] = True
    df["qc_flags"] = ""

    def add_flag(mask, flag_name):
        """Add a flag to rows matching mask."""
        df.loc[mask, "qc_pass"] = False
        current_flags = df.loc[mask, "qc_flags"]
        df.loc[mask, "qc_flags"] = current_flags.apply(
            lambda x: f"{x};{flag_name}" if x else flag_name
        )

    # CDR3 alpha length check
    if "CDR3_alpha" in df.columns:
        cdr3a = df["CDR3_alpha"].fillna("")
        too_short = cdr3a.str.len() < min_cdr3_length
        too_long = cdr3a.str.len() > max_cdr3_length
        invalid_cdr3a = (too_short | too_long) & (cdr3a != "")

        add_flag(invalid_cdr3a, "cdr3a_length")
        n_invalid = invalid_cdr3a.sum()
        report.add_result(QCResult(
            n_invalid == 0,
            "cdr3_alpha_length",
            f"CDR3 alpha length: {n_invalid}/{n_rows} out of bounds",
            {"n_invalid": n_invalid, "min": min_cdr3_length, "max": max_cdr3_length}
        ))

    # CDR3 beta length check
    if "CDR3_beta" in df.columns:
        cdr3b = df["CDR3_beta"].fillna("")
        too_short = cdr3b.str.len() < min_cdr3_length
        too_long = cdr3b.str.len() > max_cdr3_length
        invalid_cdr3b = (too_short | too_long) & (cdr3b != "")

        add_flag(invalid_cdr3b, "cdr3b_length")
        n_invalid = invalid_cdr3b.sum()
        report.add_result(QCResult(
            n_invalid == 0,
            "cdr3_beta_length",
            f"CDR3 beta length: {n_invalid}/{n_rows} out of bounds",
            {"n_invalid": n_invalid, "min": min_cdr3_length, "max": max_cdr3_length}
        ))

    # Read count check
    read_cols = [c for c in df.columns if c.endswith("_reads") or c == "reads"]
    if read_cols:
        for col in read_cols:
            low_reads = df[col].fillna(0) < min_reads
            add_flag(low_reads, f"low_{col}")
            n_low = low_reads.sum()
            if n_low > 0:
                report.add_warning(f"{n_low} clonotypes with {col} < {min_reads}")

    # UMI count check
    umi_cols = [c for c in df.columns if c.endswith("_umis") or c == "umis"]
    if umi_cols:
        for col in umi_cols:
            low_umis = df[col].fillna(0) < min_umis
            add_flag(low_umis, f"low_{col}")
            n_low = low_umis.sum()
            if n_low > 0:
                report.add_warning(f"{n_low} clonotypes with {col} < {min_umis}")

    # Cell count check
    if "n_cells" in df.columns:
        zero_cells = df["n_cells"] == 0
        add_flag(zero_cells, "zero_cells")
        n_zero = zero_cells.sum()
        report.add_result(QCResult(
            n_zero == 0,
            "cell_count",
            f"Cell count: {n_zero}/{n_rows} have zero cells",
            {"n_zero": n_zero}
        ))

    # Full sequence checks
    if check_full_sequences:
        for chain in ["alpha", "beta"]:
            seq_col = f"full_{chain}_aa"
            if seq_col in df.columns:
                seqs = df[seq_col].fillna("")
                too_short = (seqs.str.len() < min_chain_length) & (seqs != "")
                too_long = seqs.str.len() > max_chain_length

                add_flag(too_short, f"{chain}_too_short")
                add_flag(too_long, f"{chain}_too_long")

                n_short = too_short.sum()
                n_long = too_long.sum()
                report.add_result(QCResult(
                    n_short == 0 and n_long == 0,
                    f"{chain}_chain_length",
                    f"{chain.title()} chain: {n_short} too short, {n_long} too long",
                    {"n_short": n_short, "n_long": n_long}
                ))

            dna_col = f"full_{chain}_dna"
            if dna_col in df.columns:
                dna_seqs = df[dna_col].fillna("")
                valid_dna = dna_seqs != ""

                # Frame check
                bad_frame = valid_dna & (dna_seqs.str.len() % 3 != 0)
                add_flag(bad_frame, f"{chain}_frame")
                n_bad_frame = bad_frame.sum()
                if n_bad_frame > 0:
                    report.add_warning(f"{n_bad_frame} {chain} chains not in frame")

                # Start codon check
                no_start = valid_dna & (~dna_seqs.str.upper().str.startswith("ATG"))
                add_flag(no_start, f"{chain}_no_start")
                n_no_start = no_start.sum()
                if n_no_start > 0:
                    report.add_warning(f"{n_no_start} {chain} chains missing start codon")

                # Stop codon check
                stop_codons = dna_seqs.str.upper().str[-3:]
                has_stop = stop_codons.isin(["TAA", "TAG", "TGA"])
                no_stop = valid_dna & (~has_stop)
                add_flag(no_stop, f"{chain}_no_stop")
                n_no_stop = no_stop.sum()
                if n_no_stop > 0:
                    report.add_warning(f"{n_no_stop} {chain} chains missing stop codon")

    # Summary
    n_pass = df["qc_pass"].sum()
    report.add_result(QCResult(
        n_pass == n_rows,
        "overall",
        f"Overall: {n_pass}/{n_rows} clonotypes pass QC",
        {"n_pass": n_pass, "n_fail": n_rows - n_pass}
    ))

    return df, report


def get_qc_summary(clonotypes: pd.DataFrame) -> dict:
    """
    Get a summary of QC statistics for a clonotype DataFrame.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        DataFrame with clonotype information (may include qc_pass column)

    Returns
    -------
    dict
        Dictionary with QC statistics
    """
    summary = {
        "n_clonotypes": len(clonotypes),
    }

    # QC pass rate
    if "qc_pass" in clonotypes.columns:
        summary["n_qc_pass"] = clonotypes["qc_pass"].sum()
        summary["qc_pass_rate"] = summary["n_qc_pass"] / len(clonotypes) if len(clonotypes) > 0 else 0

    # Chain completeness
    if "CDR3_alpha" in clonotypes.columns:
        summary["n_with_alpha"] = clonotypes["CDR3_alpha"].notna().sum()
    if "CDR3_beta" in clonotypes.columns:
        summary["n_with_beta"] = clonotypes["CDR3_beta"].notna().sum()
    if "CDR3_alpha" in clonotypes.columns and "CDR3_beta" in clonotypes.columns:
        summary["n_complete"] = (
            clonotypes["CDR3_alpha"].notna() & clonotypes["CDR3_beta"].notna()
        ).sum()

    # CDR3 length statistics
    for chain, col in [("alpha", "CDR3_alpha"), ("beta", "CDR3_beta")]:
        if col in clonotypes.columns:
            lengths = clonotypes[col].dropna().str.len()
            if len(lengths) > 0:
                summary[f"cdr3_{chain}_mean_length"] = lengths.mean()
                summary[f"cdr3_{chain}_min_length"] = lengths.min()
                summary[f"cdr3_{chain}_max_length"] = lengths.max()

    # Cell count statistics
    if "n_cells" in clonotypes.columns:
        summary["total_cells"] = clonotypes["n_cells"].sum()
        summary["mean_cells_per_clonotype"] = clonotypes["n_cells"].mean()
        summary["max_cells_per_clonotype"] = clonotypes["n_cells"].max()

    return summary
