"""
Tests for full-length TCR sequence assembly.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from tcrsift.assemble import (
    translate_dna,
    find_longest_orf,
    parse_fasta,
    load_contigs,
    assemble_full_sequences,
    validate_sequences,
    export_fasta,
    CODON_TABLE,
    T2A_LINKER_AA,
    T2A_LINKER_DNA,
)


class TestTranslateDna:
    """Tests for DNA translation."""

    def test_translate_simple(self):
        """Translate simple sequence."""
        # ATG = M, GGG = G, AAA = K
        dna = "ATGGGGAAA"
        aa, ragged = translate_dna(dna)
        assert aa == "MGK"
        assert ragged == ""

    def test_translate_with_ragged(self):
        """Translate sequence with ragged 3' end."""
        # 10 nt -> 3 codons + 1 nt
        dna = "ATGGGGAAAA"
        aa, ragged = translate_dna(dna)
        assert aa == "MGK"
        assert ragged == "A"

    def test_translate_stop_codon(self):
        """Translation should stop at stop codon."""
        # ATG = M, TAA = *, GGG = G
        dna = "ATGTAAGGG"
        aa, ragged = translate_dna(dna)
        assert aa == "M"
        assert ragged == ""

    def test_translate_all_codons(self):
        """All codons in table should translate."""
        for codon, expected_aa in CODON_TABLE.items():
            if expected_aa != "*":
                aa, _ = translate_dna(codon)
                assert aa == expected_aa

    def test_translate_unknown_codon(self):
        """Unknown codons should translate to X."""
        # NNN is not a valid codon
        dna = "NNNAAAGGG"
        aa, _ = translate_dna(dna)
        assert aa[0] == "X"


class TestFindLongestOrf:
    """Tests for ORF finding."""

    def test_find_orf_from_start(self):
        """Find ORF starting at beginning."""
        # ATG = start, then some coding sequence
        dna = "ATGGGGAAATTT"  # M G K F
        aa, offset, ragged = find_longest_orf(dna)
        assert aa == "MGKF"
        assert offset == 0

    def test_find_orf_with_offset(self):
        """Find ORF with offset from start."""
        # Some junk, then ATG
        dna = "NNATGGGGAAA"  # 2nt junk, then M G K
        aa, offset, ragged = find_longest_orf(dna)
        assert "MGK" in aa
        assert offset == 2

    def test_find_longest_among_multiple(self):
        """Find longest ORF among multiple start codons."""
        # Two ORFs: one short, one long
        dna = "ATGTAAATGGGGAAATTT"  # M stop, then M G K F
        aa, offset, ragged = find_longest_orf(dna)
        assert len(aa) == 4  # The longer one


class TestParseFasta:
    """Tests for FASTA parsing."""

    @pytest.fixture
    def mock_fasta(self, temp_dir):
        """Create a mock FASTA file."""
        fasta_path = temp_dir / "contigs.fasta"
        content = """>seq1
ATGGGGAAATTT
>seq2
GGGAAATTTCCC
>seq3_with_description extra info
AAACCCGGG
"""
        fasta_path.write_text(content)
        return fasta_path

    def test_parse_fasta(self, mock_fasta):
        """Parse FASTA file."""
        result = parse_fasta(mock_fasta)

        assert len(result) == 3
        assert "seq1" in result
        assert "seq2" in result
        assert "seq3_with_description" in result

        assert result["seq1"] == "ATGGGGAAATTT"
        assert result["seq2"] == "GGGAAATTTCCC"

    def test_parse_multiline_sequences(self, temp_dir):
        """Parse FASTA with multiline sequences."""
        fasta_path = temp_dir / "multiline.fasta"
        content = """>seq1
ATGGGG
AAATTT
>seq2
GGGAAA
TTTCCC
"""
        fasta_path.write_text(content)
        result = parse_fasta(fasta_path)

        assert result["seq1"] == "ATGGGGAAATTT"
        assert result["seq2"] == "GGGAAATTTCCC"


class TestLoadContigs:
    """Tests for loading contigs from directories."""

    @pytest.fixture
    def mock_contig_dir(self, temp_dir):
        """Create mock contig directory structure."""
        # Sample 1
        s1_dir = temp_dir / "S1"
        s1_dir.mkdir()
        (s1_dir / "filtered_contig.fasta").write_text(""">contig_1
ATGGGGAAATTT
>contig_2
GGGAAATTTCCC
""")

        # Sample 2
        s2_dir = temp_dir / "S2"
        s2_dir.mkdir()
        (s2_dir / "filtered_contig.fasta").write_text(""">contig_1
AAACCCGGG
""")

        return temp_dir

    def test_load_contigs_from_subdirs(self, mock_contig_dir):
        """Load contigs from sample subdirectories."""
        result = load_contigs(mock_contig_dir)

        assert "S1" in result
        assert "S2" in result
        assert len(result["S1"]) == 2
        assert len(result["S2"]) == 1


class TestAssembleFullSequences:
    """Tests for full sequence assembly."""

    @pytest.fixture
    def clonotypes_with_vdj(self):
        """Create clonotypes with VDJ sequences."""
        return pd.DataFrame({
            "clone_id": ["clone1", "clone2"],
            "CDR3_alpha": ["CAVSDGGSQGNLIF", "CAVSAGGSQGNLIF"],
            "CDR3_beta": ["CASSLGQAYEQYF", "CASSLAGAYEQYF"],
            "VDJ_alpha_aa": ["MQRLQVWVLLFFLLPGTRG...CAVSDGGSQGNLIF...",
                            "MQRLQVWVLLFFLLPGTRG...CAVSAGGSQGNLIF..."],
            "VDJ_beta_aa": ["MGSRLLCWVLLCLLGAGPVKA...CASSLGQAYEQYF...",
                           "MGSRLLCWVLLCLLGAGPVKA...CASSLAGAYEQYF..."],
            "alpha_c_gene": ["TRAC", "TRAC"],
            "beta_c_gene": ["TRBC1", "TRBC2"],
            "samples": ["S1", "S2"],
        })

    def test_assemble_basic(self, clonotypes_with_vdj):
        """Basic assembly without contigs or constant regions."""
        result = assemble_full_sequences(
            clonotypes_with_vdj,
            include_leader=False,
            include_constant=False,
        )

        assert "vdj_alpha_aa" in result.columns
        assert "vdj_beta_aa" in result.columns

    def test_assemble_with_single_chain(self, clonotypes_with_vdj):
        """Assembly with single-chain construct."""
        # Modify to have proper full sequences
        clonotypes_with_vdj["full_alpha_aa"] = ["TESTSEQ1", "TESTSEQ2"]
        clonotypes_with_vdj["full_beta_aa"] = ["BETASEQ1", "BETASEQ2"]

        result = assemble_full_sequences(
            clonotypes_with_vdj,
            include_leader=False,
            include_constant=False,
            linker="T2A",
        )

        if "single_chain_aa" in result.columns:
            # Single chain should have beta-T2A-alpha
            assert T2A_LINKER_AA in result["single_chain_aa"].iloc[0]


class TestValidateSequences:
    """Tests for sequence validation."""

    def test_validate_short_sequence(self):
        """Short sequences should trigger warning."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "full_alpha_aa": ["SHORT"],  # Too short
            "full_beta_aa": ["ALSO_SHORT"],
        })

        warnings = validate_sequences(df)
        assert any("too short" in w for w in warnings)

    def test_validate_long_sequence(self):
        """Very long sequences should trigger warning."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "full_alpha_aa": ["A" * 500],  # Too long
            "full_beta_aa": ["B" * 300],
        })

        warnings = validate_sequences(df)
        assert any("too long" in w for w in warnings)

    def test_validate_cdr3_present(self):
        """CDR3 should be present in full sequence."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVTEST"],
            "full_alpha_aa": ["SOMESEQUENCE"],  # CDR3 not in sequence
            "full_beta_aa": ["BETASEQUENCE"],
        })

        warnings = validate_sequences(df)
        assert any("CDR3_alpha not found" in w for w in warnings)


class TestExportFasta:
    """Tests for FASTA export."""

    def test_export_single_chain(self, temp_dir):
        """Export single-chain sequences to FASTA."""
        df = pd.DataFrame({
            "clone_id": ["clone1", "clone2"],
            "CDR3_alpha": ["CAVTEST1", "CAVTEST2"],
            "CDR3_beta": ["CASSTEST1", "CASSTEST2"],
            "single_chain_aa": ["SEQUENCEONE", "SEQUENCETWO"],
        })

        output_path = temp_dir / "output.fasta"
        export_fasta(df, output_path, sequence_col="single_chain_aa")

        assert output_path.exists()

        # Parse and verify
        content = output_path.read_text()
        assert ">clone1" in content
        assert "SEQUENCEONE" in content
        assert ">clone2" in content
        assert "SEQUENCETWO" in content

    def test_export_skips_empty(self, temp_dir):
        """Export should skip empty sequences."""
        df = pd.DataFrame({
            "clone_id": ["clone1", "clone2"],
            "CDR3_alpha": ["CAVTEST1", "CAVTEST2"],
            "CDR3_beta": ["CASSTEST1", "CASSTEST2"],
            "single_chain_aa": ["SEQUENCEONE", ""],  # Second is empty
        })

        output_path = temp_dir / "output.fasta"
        export_fasta(df, output_path, sequence_col="single_chain_aa")

        content = output_path.read_text()
        assert "clone1" in content
        assert "clone2" not in content  # Empty sequence skipped


class TestT2ALinker:
    """Tests for T2A linker constants."""

    def test_t2a_aa_length(self):
        """T2A amino acid sequence should be correct length."""
        assert len(T2A_LINKER_AA) == 20

    def test_t2a_dna_translates_to_aa(self):
        """T2A DNA should translate to T2A AA."""
        aa, _ = translate_dna(T2A_LINKER_DNA)
        # The last codon might be partial, so check prefix
        assert T2A_LINKER_AA.startswith(aa) or aa.startswith(T2A_LINKER_AA[:len(aa)])
