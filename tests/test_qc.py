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

"""Tests for QC validation module."""

import pandas as pd
import pytest

from tcrsift.qc import (
    QCReport,
    QCResult,
    check_cdr3_length,
    check_chain_length,
    check_reading_frame,
    check_start_codon,
    check_stop_codon,
    find_repeated_kmers,
    get_qc_summary,
    validate_clonotypes,
    validate_sequence,
)


class TestFindRepeatedKmers:
    """Tests for find_repeated_kmers function."""

    def test_no_repeats(self):
        """Sequence with no repeated k-mers."""
        seq = "ABCDEFGHIJKLMNOP"
        result = find_repeated_kmers(seq, k=5, min_repeat=2)
        assert result == []

    def test_with_repeats(self):
        """Sequence with repeated k-mers."""
        seq = "ABCDEABCDEFGH"  # "ABCDE" appears twice
        result = find_repeated_kmers(seq, k=5, min_repeat=2)
        assert len(result) == 1
        assert result[0][0] == "ABCDE"
        assert result[0][1] == 2

    def test_multiple_repeats(self):
        """Sequence with multiple repeated k-mers."""
        seq = "ABCABCABC"  # "ABC" appears 3 times at different positions
        result = find_repeated_kmers(seq, k=3, min_repeat=2)
        # ABC appears at positions 0, 3, 6
        assert any(kmer == "ABC" and count >= 2 for kmer, count in result)

    def test_sequence_shorter_than_k(self):
        """Sequence shorter than k-mer size."""
        seq = "ABC"
        result = find_repeated_kmers(seq, k=5, min_repeat=2)
        assert result == []

    def test_empty_sequence(self):
        """Empty sequence."""
        result = find_repeated_kmers("", k=5, min_repeat=2)
        assert result == []

    def test_min_repeat_threshold(self):
        """Test minimum repeat threshold."""
        seq = "ABCDEABCDE"  # "ABCDE" appears twice
        # Should find it with min_repeat=2
        result = find_repeated_kmers(seq, k=5, min_repeat=2)
        assert len(result) >= 1
        # Should not find it with min_repeat=3
        result = find_repeated_kmers(seq, k=5, min_repeat=3)
        assert len(result) == 0


class TestCodonChecks:
    """Tests for codon validation functions."""

    def test_start_codon_present(self):
        """DNA with ATG start codon."""
        assert check_start_codon("ATGCGATCG") is True

    def test_start_codon_missing(self):
        """DNA without ATG start codon."""
        assert check_start_codon("GCGATCGATG") is False

    def test_start_codon_lowercase(self):
        """Lowercase DNA with start codon."""
        assert check_start_codon("atgcgatcg") is True

    def test_start_codon_empty(self):
        """Empty sequence."""
        assert check_start_codon("") is False

    def test_stop_codon_taa(self):
        """DNA ending with TAA stop codon."""
        assert check_stop_codon("ATGCGATAA") is True

    def test_stop_codon_tag(self):
        """DNA ending with TAG stop codon."""
        assert check_stop_codon("ATGCGATAG") is True

    def test_stop_codon_tga(self):
        """DNA ending with TGA stop codon."""
        assert check_stop_codon("ATGCGATGA") is True

    def test_stop_codon_missing(self):
        """DNA without stop codon."""
        assert check_stop_codon("ATGCGACCC") is False

    def test_stop_codon_lowercase(self):
        """Lowercase DNA with stop codon."""
        assert check_stop_codon("atgcgataa") is True

    def test_stop_codon_empty(self):
        """Empty sequence."""
        assert check_stop_codon("") is False

    def test_stop_codon_too_short(self):
        """Sequence too short."""
        assert check_stop_codon("AT") is False


class TestReadingFrameCheck:
    """Tests for reading frame validation."""

    def test_valid_frame(self):
        """Sequence length divisible by 3."""
        assert check_reading_frame("ATGATG") is True  # 6 nucleotides
        assert check_reading_frame("ATGATGATG") is True  # 9 nucleotides

    def test_invalid_frame_remainder_1(self):
        """Sequence with remainder 1."""
        assert check_reading_frame("ATGATGA") is False  # 7 nucleotides

    def test_invalid_frame_remainder_2(self):
        """Sequence with remainder 2."""
        assert check_reading_frame("ATGATGAT") is False  # 8 nucleotides

    def test_empty_sequence(self):
        """Empty sequence (0 is divisible by 3)."""
        assert check_reading_frame("") is False


class TestCdr3LengthCheck:
    """Tests for CDR3 length validation."""

    def test_valid_length(self):
        """CDR3 within bounds."""
        assert check_cdr3_length("CASSLTGELF", min_length=5, max_length=40) is True

    def test_too_short(self):
        """CDR3 too short."""
        assert check_cdr3_length("CASS", min_length=5, max_length=40) is False

    def test_too_long(self):
        """CDR3 too long."""
        long_cdr3 = "C" * 50
        assert check_cdr3_length(long_cdr3, min_length=5, max_length=40) is False

    def test_empty(self):
        """Empty CDR3."""
        assert check_cdr3_length("", min_length=5, max_length=40) is False

    def test_none(self):
        """None CDR3."""
        assert check_cdr3_length(None, min_length=5, max_length=40) is False

    def test_boundary_min(self):
        """CDR3 at minimum boundary."""
        assert check_cdr3_length("CASSG", min_length=5, max_length=40) is True

    def test_boundary_max(self):
        """CDR3 at maximum boundary."""
        cdr3 = "C" * 40
        assert check_cdr3_length(cdr3, min_length=5, max_length=40) is True


class TestChainLengthCheck:
    """Tests for chain length validation."""

    def test_valid_length(self):
        """Chain within bounds."""
        chain = "M" * 150
        assert check_chain_length(chain, min_length=80, max_length=450) is True

    def test_too_short(self):
        """Chain too short."""
        chain = "M" * 50
        assert check_chain_length(chain, min_length=80, max_length=450) is False

    def test_too_long(self):
        """Chain too long."""
        chain = "M" * 500
        assert check_chain_length(chain, min_length=80, max_length=450) is False

    def test_empty(self):
        """Empty chain."""
        assert check_chain_length("", min_length=80, max_length=450) is False


class TestQCResult:
    """Tests for QCResult dataclass."""

    def test_create_passed_result(self):
        """Create a passing result."""
        result = QCResult(True, "test_check", "Test passed")
        assert result.passed is True
        assert result.check_name == "test_check"
        assert result.message == "Test passed"
        assert result.details == {}

    def test_create_failed_result(self):
        """Create a failing result with details."""
        result = QCResult(False, "test_check", "Test failed", {"value": 10})
        assert result.passed is False
        assert result.details == {"value": 10}


class TestQCReport:
    """Tests for QCReport dataclass."""

    def test_empty_report(self):
        """Empty report should pass."""
        report = QCReport()
        assert report.passed is True
        assert report.num_passed == 0
        assert report.num_failed == 0

    def test_all_passed(self):
        """Report with all passed checks."""
        report = QCReport()
        report.add_result(QCResult(True, "check1", "Passed"))
        report.add_result(QCResult(True, "check2", "Passed"))
        assert report.passed is True
        assert report.num_passed == 2
        assert report.num_failed == 0

    def test_one_failed(self):
        """Report with one failed check."""
        report = QCReport()
        report.add_result(QCResult(True, "check1", "Passed"))
        report.add_result(QCResult(False, "check2", "Failed"))
        assert report.passed is False
        assert report.num_passed == 1
        assert report.num_failed == 1

    def test_add_warning(self):
        """Add warning to report."""
        report = QCReport()
        report.add_warning("This is a warning")
        assert len(report.warnings) == 1
        assert report.passed is True  # Warnings don't fail the report

    def test_summary(self):
        """Generate summary string."""
        report = QCReport()
        report.add_result(QCResult(True, "check1", "Passed"))
        summary = report.summary()
        assert "PASSED" in summary
        assert "1/1" in summary


class TestValidateSequence:
    """Tests for validate_sequence function."""

    def test_valid_amino_acid_sequence(self):
        """Valid amino acid sequence."""
        # Use a long diverse sequence without repeats (typical TCR chain)
        seq = (
            "MGLTSAVPKDQRLSGVDRIRSKGYKGKVHVIFNYGPSGQEVNLKKVQDVGSEAELLFCQAKP"
            "HFSNYSQQMKFSASPVVVSGQEGTRHTASLTFSPPSGKGGKTVSILVGKNALRQITVNDQVF"
            "GKTLTITKQSGTQVFLNDKVVLTTGTTLLLGKQVLRQI"
        )
        report = validate_sequence(seq, is_dna=False, min_chain_length=80, kmer_size=10, min_kmer_repeat=2)
        # Check that at least the chain length check passed
        chain_result = next(r for r in report.results if r.check_name == "chain_length")
        assert chain_result.passed is True

    def test_short_sequence(self):
        """Sequence too short."""
        seq = "M" * 50  # Below minimum
        report = validate_sequence(seq, is_dna=False, min_chain_length=80)
        assert report.passed is False

    def test_with_cdr3_present(self):
        """Sequence with CDR3 present."""
        cdr3 = "CASSLTGELF"
        seq = "MMMM" + cdr3 + "MMMM" * 30  # CDR3 in middle, enough length
        report = validate_sequence(seq, cdr3=cdr3, is_dna=False, min_chain_length=50)
        assert any(r.check_name == "cdr3_present" and r.passed for r in report.results)

    def test_with_cdr3_missing(self):
        """Sequence with CDR3 missing."""
        cdr3 = "CASSLTGELF"
        seq = "M" * 150  # No CDR3
        report = validate_sequence(seq, cdr3=cdr3, is_dna=False)
        assert any(r.check_name == "cdr3_present" and not r.passed for r in report.results)

    def test_with_cdr3_duplicated(self):
        """Sequence with CDR3 appearing twice."""
        cdr3 = "CASSLTGELF"
        seq = cdr3 + "MMM" + cdr3 + "M" * 100
        report = validate_sequence(seq, cdr3=cdr3, is_dna=False, min_chain_length=50)
        cdr3_result = next(r for r in report.results if r.check_name == "cdr3_present")
        assert cdr3_result.passed is False
        assert cdr3_result.details["count"] == 2

    def test_dna_sequence_valid(self):
        """Valid DNA sequence with start and stop codons."""
        seq = "ATG" + "GCA" * 100 + "TAA"  # Start + body + stop
        report = validate_sequence(seq, is_dna=True, check_codons=True, min_chain_length=50)
        assert any(r.check_name == "start_codon" and r.passed for r in report.results)
        assert any(r.check_name == "stop_codon" and r.passed for r in report.results)
        assert any(r.check_name == "reading_frame" and r.passed for r in report.results)

    def test_dna_sequence_no_start(self):
        """DNA sequence without start codon."""
        seq = "GCA" * 101 + "TAA"
        report = validate_sequence(seq, is_dna=True, check_codons=True, min_chain_length=50)
        assert any(r.check_name == "start_codon" and not r.passed for r in report.results)

    def test_dna_sequence_bad_frame(self):
        """DNA sequence not in frame."""
        seq = "ATG" + "GCAGC"  # 8 nucleotides, not divisible by 3
        report = validate_sequence(seq, is_dna=True, min_chain_length=1)
        assert any(r.check_name == "reading_frame" and not r.passed for r in report.results)

    def test_empty_sequence(self):
        """Empty sequence should fail."""
        report = validate_sequence("")
        assert report.passed is False

    def test_none_sequence(self):
        """None sequence should fail."""
        report = validate_sequence(None)
        assert report.passed is False

    def test_kmer_detection(self):
        """Detect repeated k-mers in sequence."""
        # Sequence with repeated 10-mer
        repeated_part = "ABCDEFGHIJ"
        seq = repeated_part + "XYZ" + repeated_part + "M" * 100
        report = validate_sequence(seq, is_dna=False, kmer_size=10, min_kmer_repeat=2, min_chain_length=50)
        assert any(r.check_name == "repeated_kmers" and not r.passed for r in report.results)


class TestValidateClonotypes:
    """Tests for validate_clonotypes function."""

    def test_valid_clonotypes(self):
        """Valid clonotypes should pass."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP"],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [10, 5],
            "reads": [100, 50],
            "umis": [20, 10],
        })
        result_df, report = validate_clonotypes(df)
        assert "qc_pass" in result_df.columns
        assert result_df["qc_pass"].all()

    def test_cdr3_too_short(self):
        """CDR3 too short should be flagged."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CAS", "CASSLTGELF"],  # First is too short
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [10, 5],
        })
        result_df, report = validate_clonotypes(df, min_cdr3_length=5)
        assert not result_df.loc[0, "qc_pass"]
        assert result_df.loc[1, "qc_pass"]
        assert "cdr3a_length" in result_df.loc[0, "qc_flags"]

    def test_cdr3_too_long(self):
        """CDR3 too long should be flagged."""
        df = pd.DataFrame({
            "CDR3_alpha": ["C" * 50, "CASSLTGELF"],  # First is too long
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [10, 5],
        })
        result_df, report = validate_clonotypes(df, max_cdr3_length=40)
        assert not result_df.loc[0, "qc_pass"]
        assert result_df.loc[1, "qc_pass"]

    def test_zero_cells(self):
        """Zero cells should be flagged."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP"],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [0, 5],
        })
        result_df, report = validate_clonotypes(df)
        assert not result_df.loc[0, "qc_pass"]
        assert "zero_cells" in result_df.loc[0, "qc_flags"]

    def test_low_reads_warning(self):
        """Low reads should generate warning."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP"],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "reads": [5, 100],  # First has low reads
            "n_cells": [10, 5],
        })
        result_df, report = validate_clonotypes(df, min_reads=10)
        assert len(report.warnings) > 0
        assert "low_reads" in result_df.loc[0, "qc_flags"]

    def test_empty_dataframe(self):
        """Empty dataframe should work."""
        df = pd.DataFrame(columns=["CDR3_alpha", "CDR3_beta", "n_cells"])
        result_df, report = validate_clonotypes(df)
        assert len(result_df) == 0


class TestGetQCSummary:
    """Tests for get_qc_summary function."""

    def test_basic_summary(self):
        """Basic summary statistics."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP", None],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY", "CASSYEQY"],
            "n_cells": [10, 5, 3],
        })
        summary = get_qc_summary(df)
        assert summary["n_clonotypes"] == 3
        assert summary["n_with_alpha"] == 2
        assert summary["n_with_beta"] == 3
        assert summary["n_complete"] == 2

    def test_with_qc_pass(self):
        """Summary with qc_pass column."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP"],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [10, 5],
            "qc_pass": [True, False],
        })
        summary = get_qc_summary(df)
        assert summary["n_qc_pass"] == 1
        assert summary["qc_pass_rate"] == 0.5

    def test_cdr3_length_stats(self):
        """CDR3 length statistics."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSG", "CASSLTGELF"],  # lengths 5 and 10
            "CDR3_beta": ["CSARDR", "CASSLGQAYEQY"],  # lengths 6 and 12
            "n_cells": [10, 5],
        })
        summary = get_qc_summary(df)
        assert summary["cdr3_alpha_mean_length"] == 7.5
        assert summary["cdr3_alpha_min_length"] == 5
        assert summary["cdr3_alpha_max_length"] == 10

    def test_cell_count_stats(self):
        """Cell count statistics."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CASSLTGELF", "CAVSGGSYIP"],
            "CDR3_beta": ["CSARDRVGNTIY", "CASSLGQAYEQY"],
            "n_cells": [10, 20],
        })
        summary = get_qc_summary(df)
        assert summary["total_cells"] == 30
        assert summary["mean_cells_per_clonotype"] == 15
        assert summary["max_cells_per_clonotype"] == 20
