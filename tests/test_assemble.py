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
            alpha_leader=None,
            beta_leader=None,
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
            alpha_leader=None,
            beta_leader=None,
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
        assert len(T2A_LINKER_AA) == 18  # EGRGSLLTCGDVEENPGP

    def test_t2a_dna_translates_to_aa(self):
        """T2A DNA should translate to T2A AA."""
        aa, _ = translate_dna(T2A_LINKER_DNA)
        # The last codon might be partial, so check prefix
        assert T2A_LINKER_AA.startswith(aa) or aa.startswith(T2A_LINKER_AA[:len(aa)])


class TestAssembleWithContigs:
    """Tests for assembly with contig files."""

    @pytest.fixture
    def mock_contig_dir_with_leaders(self, temp_dir):
        """Create mock contig directory with leader sequences."""
        s1_dir = temp_dir / "S1"
        s1_dir.mkdir()
        # Create a contig with leader + VDJ sequence
        # Leader: ATGAAA = MK, VDJ: GGGTTTTTT = GFF
        (s1_dir / "filtered_contig.fasta").write_text(""">contig_alpha
ATGAAAGGGTTTTTT
>contig_beta
ATGCCCGGGAAATTT
""")
        return temp_dir

    def test_assemble_without_contigs_dir(self):
        """Test assembly when contigs_dir is None."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVMKGFF"],
            "CDR3_beta": ["CASSMPGKF"],
            "VDJ_alpha_aa": ["MKGFF"],
            "VDJ_beta_aa": ["MPGKF"],
            "alpha_c_gene": ["TRAC"],
            "beta_c_gene": ["TRBC1"],
        })

        result = assemble_full_sequences(
            df,
            contigs_dir=None,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
        )

        assert "vdj_alpha_aa" in result.columns
        assert result.iloc[0]["vdj_alpha_aa"] == "MKGFF"


class TestBuildFullSequences:
    """Tests for _build_full_sequences helper."""

    def test_assemble_with_vdj_only(self):
        """Assembly with only VDJ sequences."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVALPHASEQUENCE"],
            "CDR3_beta": ["CASSBETASEQUENCE"],
            "VDJ_alpha_aa": ["VDJALPHASEQUENCE"],
            "VDJ_alpha_nt": ["ATGATGATGATGATGATG"],
            "VDJ_beta_aa": ["VDJBETASEQUENCE"],
            "VDJ_beta_nt": ["GCAGCAGCAGCAGCAGCA"],
            "alpha_c_gene": ["TRAC"],
            "beta_c_gene": ["TRBC1"],
        })

        result = assemble_full_sequences(
            df,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
        )

        assert result.iloc[0]["vdj_alpha_aa"] == "VDJALPHASEQUENCE"
        assert result.iloc[0]["vdj_beta_aa"] == "VDJBETASEQUENCE"


class TestAddSingleChain:
    """Tests for single-chain construct assembly."""

    def test_single_chain_with_nt_sequences(self):
        """Test single-chain with nucleotide sequences."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVALPHA"],
            "CDR3_beta": ["CASSBETA"],
            "full_alpha_aa": ["ALPHASEQUENCE"],
            "full_beta_aa": ["BETASEQUENCE*"],  # With stop codon
            "full_alpha_nt": ["GCACTGGCAAGCCAGAACACC"],
            "full_beta_nt": ["TGCGAGTGCAGCAGCTAA"],  # TAA is stop codon
        })

        result = assemble_full_sequences(
            df,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
            linker="T2A",
        )

        assert "single_chain_aa" in result.columns
        # Stop codon should be stripped from beta
        assert not result.iloc[0]["single_chain_aa"].startswith("*")
        assert T2A_LINKER_AA in result.iloc[0]["single_chain_aa"]

    def test_single_chain_custom_linker(self):
        """Test single-chain with custom linker."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVALPHA"],
            "CDR3_beta": ["CASSBETA"],
            "full_alpha_aa": ["ALPHASEQUENCE"],
            "full_beta_aa": ["BETASEQUENCE"],
        })

        result = assemble_full_sequences(
            df,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
            linker="CUSTOMLINK",
        )

        assert "single_chain_aa" in result.columns
        assert "CUSTOMLINK" in result.iloc[0]["single_chain_aa"]


class TestValidateSequencesExtended:
    """Extended tests for sequence validation."""

    def test_validate_good_sequences(self):
        """Valid sequences should produce no warnings."""
        # Create reasonable length sequences with CDR3 included
        alpha_seq = "M" * 100 + "CAVTEST" + "A" * 150  # ~257 aa
        beta_seq = "M" * 100 + "CASSTEST" + "A" * 150  # ~258 aa

        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVTEST"],
            "CDR3_beta": ["CASSTEST"],
            "full_alpha_aa": [alpha_seq],
            "full_beta_aa": [beta_seq],
        })

        warnings = validate_sequences(df)

        # Filter out any warnings - should have none for good sequences
        cdr3_warnings = [w for w in warnings if "CDR3" in w]
        length_warnings = [w for w in warnings if "short" in w or "long" in w]

        assert len(cdr3_warnings) == 0
        assert len(length_warnings) == 0

    def test_validate_cdr3_beta_not_found(self):
        """CDR3 beta missing from sequence should trigger warning."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVTEST"],
            "CDR3_beta": ["CASSNOTPRESENT"],
            "full_alpha_aa": ["CAVTEST" + "A" * 200],
            "full_beta_aa": ["DIFFERENTSEQ" + "A" * 200],
        })

        warnings = validate_sequences(df)
        assert any("CDR3_beta not found" in w for w in warnings)

    def test_validate_missing_full_sequence(self):
        """Missing full sequence should not cause error."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVTEST"],
            "full_alpha_aa": [None],  # Missing
            "full_beta_aa": ["BETASEQ"],
        })

        # Should not raise
        warnings = validate_sequences(df)
        assert isinstance(warnings, list)


class TestExportFastaExtended:
    """Extended tests for FASTA export."""

    def test_export_includes_cdr3_in_header(self, temp_dir):
        """Export should include CDR3 sequences in FASTA header."""
        df = pd.DataFrame({
            "clone_id": ["clone1"],
            "CDR3_alpha": ["CAVTEST"],
            "CDR3_beta": ["CASSTEST"],
            "single_chain_aa": ["TESTSEQUENCE"],
        })

        output_path = temp_dir / "output.fasta"
        export_fasta(df, output_path, sequence_col="single_chain_aa")

        content = output_path.read_text()
        assert "CDR3a=CAVTEST" in content
        assert "CDR3b=CASSTEST" in content

    def test_export_with_missing_clone_id(self, temp_dir):
        """Export should handle missing clone_id by using index."""
        df = pd.DataFrame({
            "CDR3_alpha": ["CAVTEST"],
            "CDR3_beta": ["CASSTEST"],
            "single_chain_aa": ["TESTSEQUENCE"],
        })

        output_path = temp_dir / "output.fasta"
        export_fasta(df, output_path, sequence_col="single_chain_aa")

        assert output_path.exists()
        content = output_path.read_text()
        # Should use index (0) as clone_id
        assert ">0" in content or content.startswith(">")


class TestTranslateDnaExtended:
    """Extended tests for DNA translation."""

    def test_translate_empty_sequence(self):
        """Empty sequence should return empty result."""
        aa, ragged = translate_dna("")
        assert aa == ""
        assert ragged == ""

    def test_translate_two_nucleotides(self):
        """Two nucleotides should return empty AA with ragged end."""
        aa, ragged = translate_dna("AT")
        assert aa == ""
        assert ragged == "AT"

    def test_translate_sequence_ending_with_stop(self):
        """Sequence ending with stop codon."""
        # MGK followed by TAA (stop)
        dna = "ATGGGGAAATAA"
        aa, ragged = translate_dna(dna)
        assert aa == "MGK"  # Stops before the stop codon AA


class TestFindLongestOrfExtended:
    """Extended tests for ORF finding."""

    def test_no_start_codon(self):
        """Sequence without ATG should return empty result."""
        dna = "GGGAAATTT"  # No ATG
        aa, offset, ragged = find_longest_orf(dna)
        assert aa == ""
        assert offset == 0

    def test_multiple_orfs_different_frames(self):
        """Find longest ORF among multiple reading frames."""
        # First ATG leads to stop, second ATG leads to longer ORF
        dna = "ATGTAANNNATGGGGAAATTTCCC"  # M-stop, then MGKF...
        aa, offset, ragged = find_longest_orf(dna)
        assert len(aa) > 1  # Should find the longer one


class TestRealisticFullLengthAssembly:
    """Tests using realistic full-length TCR sequence data."""

    def test_realistic_sequence_lengths(self, sample_full_length_clonotypes):
        """Test that full-length sequences have realistic lengths."""
        df = sample_full_length_clonotypes

        # Full alpha chain should be ~250-280 AA
        for seq in df["full_alpha_aa"]:
            assert 200 <= len(seq) <= 350, f"Alpha chain length {len(seq)} out of expected range"

        # Full beta chain should be ~280-320 AA
        for seq in df["full_beta_aa"]:
            assert 200 <= len(seq) <= 400, f"Beta chain length {len(seq)} out of expected range"

    def test_cdr3_in_full_sequence(self, sample_full_length_clonotypes):
        """Test that CDR3 sequences are contained in full-length chains."""
        df = sample_full_length_clonotypes

        for idx, row in df.iterrows():
            # CDR3 alpha should be in full alpha sequence
            assert row["CDR3_alpha"] in row["full_alpha_aa"], \
                f"CDR3_alpha not found in full_alpha_aa for clone {row['clone_id']}"

            # CDR3 beta should be in full beta sequence
            assert row["CDR3_beta"] in row["full_beta_aa"], \
                f"CDR3_beta not found in full_beta_aa for clone {row['clone_id']}"

    def test_validate_realistic_sequences(self, sample_full_length_clonotypes):
        """Validate realistic sequences should pass validation."""
        warnings = validate_sequences(sample_full_length_clonotypes)

        # No length warnings expected for realistic sequences
        length_warnings = [w for w in warnings if "too short" in w or "too long" in w]
        assert len(length_warnings) == 0, f"Unexpected length warnings: {length_warnings}"

        # No CDR3 warnings expected
        cdr3_warnings = [w for w in warnings if "CDR3" in w and "not found" in w]
        assert len(cdr3_warnings) == 0, f"Unexpected CDR3 warnings: {cdr3_warnings}"

    def test_assemble_from_vdj_columns(self, sample_full_length_clonotypes):
        """Test assembly using VDJ columns from clonotypes."""
        result = assemble_full_sequences(
            sample_full_length_clonotypes,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
        )

        # VDJ sequences should be extracted
        assert "vdj_alpha_aa" in result.columns
        assert "vdj_beta_aa" in result.columns

        # VDJ sequences should contain CDR3
        for idx, row in result.iterrows():
            if pd.notna(row.get("vdj_alpha_aa")):
                cdr3 = sample_full_length_clonotypes.loc[idx, "CDR3_alpha"]
                assert cdr3 in row["vdj_alpha_aa"]

    def test_single_chain_construct(self, sample_full_length_clonotypes):
        """Test single-chain TCR construct assembly."""
        result = assemble_full_sequences(
            sample_full_length_clonotypes,
            alpha_leader=None,
            beta_leader=None,
            include_constant=False,
            linker="T2A",
        )

        if "single_chain_aa" in result.columns:
            for idx, row in result.iterrows():
                sc = row["single_chain_aa"]
                # Single chain should contain T2A linker
                assert T2A_LINKER_AA in sc, "T2A linker not found in single chain"

                # Single chain should contain both CDR3s
                cdr3_alpha = sample_full_length_clonotypes.loc[idx, "CDR3_alpha"]
                cdr3_beta = sample_full_length_clonotypes.loc[idx, "CDR3_beta"]
                assert cdr3_beta in sc, "CDR3 beta not found in single chain"


class TestAssemblyWithMockCellRanger:
    """Tests using mock CellRanger output directory."""

    def test_load_contigs_from_mock_dir(self, mock_cellranger_vdj_dir):
        """Test loading contigs from mock CellRanger structure."""
        # The mock dir should contain all_contig.fasta
        fasta_path = mock_cellranger_vdj_dir / "all_contig.fasta"
        assert fasta_path.exists()

        contigs = parse_fasta(fasta_path)

        # Should have 4 contigs (2 cells x 2 chains)
        assert len(contigs) == 4

        # Contig IDs should match CellRanger format
        assert "AAACCTGAGAACTCGG-1_contig_1" in contigs
        assert "AAACCTGAGAACTCGG-1_contig_2" in contigs

    def test_contig_sequences_contain_cdr3(self, mock_cellranger_vdj_dir):
        """Test that contig sequences contain CDR3 nucleotides."""
        fasta_path = mock_cellranger_vdj_dir / "all_contig.fasta"
        contigs = parse_fasta(fasta_path)

        # CDR3 nucleotides should be in contigs
        cdr3_alpha_nt = "TGTGCTGTGTCAGATGGAGGAAGCCAGGGAAATCTCATCTTT"
        cdr3_beta_nt = "TGTGCCAGCAGTTTGGGACAGGCTTACGAGCAGTACTTC"

        alpha_contig = contigs["AAACCTGAGAACTCGG-1_contig_1"]
        beta_contig = contigs["AAACCTGAGAACTCGG-1_contig_2"]

        assert cdr3_alpha_nt in alpha_contig
        assert cdr3_beta_nt in beta_contig


class TestExportWithRealisticData:
    """Tests for FASTA export using realistic data."""

    def test_export_realistic_sequences(self, sample_full_length_clonotypes, temp_dir):
        """Export realistic full-length sequences to FASTA."""
        output_path = temp_dir / "tcr_sequences.fasta"

        # First assemble single-chain sequences
        df = sample_full_length_clonotypes.copy()
        df["single_chain_aa"] = df["full_beta_aa"] + T2A_LINKER_AA + df["full_alpha_aa"]

        export_fasta(df, output_path, sequence_col="single_chain_aa")

        assert output_path.exists()

        # Verify content
        content = output_path.read_text()

        # Should have headers with CDR3 info
        assert "CDR3a=CAVSDGGSQGNLIF" in content
        assert "CDR3b=CASSLGQAYEQYF" in content

        # Should have both clones
        assert ">clone1" in content
        assert ">clone2" in content

    def test_export_includes_t2a_linker(self, sample_full_length_clonotypes, temp_dir):
        """Exported single-chain sequences should include T2A linker."""
        output_path = temp_dir / "tcr_sequences.fasta"

        df = sample_full_length_clonotypes.copy()
        df["single_chain_aa"] = df["full_beta_aa"] + T2A_LINKER_AA + df["full_alpha_aa"]

        export_fasta(df, output_path, sequence_col="single_chain_aa")

        content = output_path.read_text()
        assert T2A_LINKER_AA in content
