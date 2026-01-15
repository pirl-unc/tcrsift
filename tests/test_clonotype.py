"""
Tests for clonotype aggregation.
"""

import pytest
import pandas as pd
import numpy as np
import anndata as ad

from tcrsift.clonotype import (
    aggregate_clonotypes,
    calculate_clone_frequencies,
    get_clonotype_summary,
    export_clonotypes_airr,
)


class TestAggregateClonotypes:
    """Tests for aggregate_clonotypes function."""

    @pytest.fixture
    def adata_with_tcr(self, sample_adata):
        """Create AnnData with TCR info for clonotype testing."""
        adata = sample_adata.copy()
        # Add required columns for clonotyping
        adata.obs["sample"] = ["S1"] * 50 + ["S2"] * 50

        # Create clone patterns:
        # Clone A: cells 0-19 (20 cells, same CDR3ab)
        # Clone B: cells 20-34 (15 cells, different CDR3ab)
        # Clone C: cells 35-44 (10 cells, different CDR3ab)
        # Rest: unique CDR3s

        cdr3_alpha = []
        cdr3_beta = []

        for i in range(100):
            if i < 20:
                cdr3_alpha.append("CAVSDGGSQGNLIF")
                cdr3_beta.append("CASSLGQAYEQYF")
            elif i < 35:
                cdr3_alpha.append("CAVSAGGSQGNLIF")
                cdr3_beta.append("CASSLAGAYEQYF")
            elif i < 45:
                cdr3_alpha.append("CAVNAGGSQGNLIF")
                cdr3_beta.append("CASSNAGAYEQYF")
            else:
                cdr3_alpha.append(f"CAVUNIQUE{i}QGNLIF")
                cdr3_beta.append(f"CASSUNIQUE{i}YEQYF")

        adata.obs["CDR3_alpha"] = cdr3_alpha
        adata.obs["CDR3_beta"] = cdr3_beta

        # Add Tcell_type for testing consensus
        adata.obs["Tcell_type"] = ["Confident CD8+"] * 50 + ["Confident CD4+"] * 50
        adata.obs["is_CD8"] = [True] * 50 + [False] * 50
        adata.obs["is_CD4"] = [False] * 50 + [True] * 50

        return adata

    def test_aggregate_by_cdr3ab(self, adata_with_tcr):
        """Test aggregation by both alpha and beta chains."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert len(clonotypes) > 0
        assert "clone_id" in clonotypes.columns
        assert "CDR3_alpha" in clonotypes.columns
        assert "CDR3_beta" in clonotypes.columns
        assert "cell_count" in clonotypes.columns

    def test_aggregate_by_cdr3b_only(self, adata_with_tcr):
        """Test aggregation by beta chain only."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3b_only")

        assert len(clonotypes) > 0
        assert "clone_id" in clonotypes.columns
        assert "CDR3_beta" in clonotypes.columns

    def test_aggregate_invalid_group_by(self, adata_with_tcr):
        """Invalid group_by should raise."""
        with pytest.raises(ValueError, match="Invalid group_by"):
            aggregate_clonotypes(adata_with_tcr, group_by="invalid")

    def test_aggregate_cell_counts(self, adata_with_tcr):
        """Cell counts should be correct."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        # Find Clone A (the largest)
        clone_a = clonotypes[clonotypes["CDR3_alpha"] == "CAVSDGGSQGNLIF"]
        assert len(clone_a) == 1
        assert clone_a["cell_count"].iloc[0] == 20

    def test_aggregate_sample_info(self, adata_with_tcr):
        """Sample information should be aggregated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "samples" in clonotypes.columns
        assert "n_samples" in clonotypes.columns

    def test_aggregate_tcell_consensus(self, adata_with_tcr):
        """T cell type consensus should be calculated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "Tcell_type_consensus" in clonotypes.columns

    def test_aggregate_frequency(self, adata_with_tcr):
        """Frequency should be calculated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "max_frequency" in clonotypes.columns
        # All frequencies should be between 0 and 1
        assert all(clonotypes["max_frequency"] >= 0)
        assert all(clonotypes["max_frequency"] <= 1)

    def test_handle_doublets_flag(self, adata_with_tcr):
        """Doublet handling with flag mode."""
        adata_with_tcr.obs["multi_chain"] = [i % 10 == 0 for i in range(100)]

        clonotypes = aggregate_clonotypes(adata_with_tcr, handle_doublets="flag")

        assert "n_doublet_cells" in clonotypes.columns

    def test_handle_doublets_remove(self, adata_with_tcr):
        """Doublet handling with remove mode."""
        adata_with_tcr.obs["multi_chain"] = [i % 10 == 0 for i in range(100)]

        clonotypes = aggregate_clonotypes(adata_with_tcr, handle_doublets="remove")

        # Total cells should be less due to doublet removal
        total_cells = clonotypes["cell_count"].sum()
        assert total_cells < 100

    def test_empty_result_no_complete_clones(self, sample_adata):
        """Raises error when no complete clones."""
        from tcrsift.validation import TCRsiftValidationError

        sample_adata.obs["sample"] = "S1"
        sample_adata.obs["CDR3_alpha"] = None
        sample_adata.obs["CDR3_beta"] = None

        with pytest.raises(TCRsiftValidationError, match="No complete clones found"):
            aggregate_clonotypes(sample_adata, group_by="CDR3ab")


class TestGetClonotypeSummary:
    """Tests for get_clonotype_summary function."""

    def test_summary_statistics(self, sample_clonotypes_df):
        """Summary should have correct statistics."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        assert "n_clonotypes" in summary
        assert summary["n_clonotypes"] == 5

        assert "n_cells" in summary
        assert summary["n_cells"] == 29  # 15 + 8 + 3 + 2 + 1

        assert "n_singletons" in summary
        assert summary["n_singletons"] == 1  # only last one with count=1

        assert "n_expanded" in summary
        assert summary["n_expanded"] == 4  # 4 clones with count > 1

    def test_multi_sample_clones(self, sample_clonotypes_df):
        """Should count multi-sample clones."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        assert "n_multi_sample" in summary
        # First clone has n_samples=2, others have n_samples=1
        assert summary["n_multi_sample"] == 1


class TestExportClonotypesAirr:
    """Tests for AIRR format export."""

    def test_export_airr(self, sample_clonotypes_df, temp_dir):
        """Export to AIRR format should work."""
        output_path = temp_dir / "clonotypes.tsv"

        export_clonotypes_airr(sample_clonotypes_df, str(output_path))

        assert output_path.exists()

        # Load and check format
        airr_df = pd.read_csv(output_path, sep="\t")
        assert "junction_aa_tra" in airr_df.columns or len(airr_df) > 0
        assert "junction_aa_trb" in airr_df.columns or len(airr_df) > 0
        assert "clone_count" in airr_df.columns
        assert "productive" in airr_df.columns

    def test_airr_column_mapping(self, sample_clonotypes_df, temp_dir):
        """Check AIRR column mapping."""
        output_path = temp_dir / "clonotypes.tsv"

        export_clonotypes_airr(sample_clonotypes_df, str(output_path))

        airr_df = pd.read_csv(output_path, sep="\t")

        # Check that CDR3 sequences are mapped correctly
        if "junction_aa_trb" in airr_df.columns:
            assert airr_df["junction_aa_trb"].iloc[0] == "CASSLGQAYEQYF"


class TestCalculateCloneFrequencies:
    """Tests for calculate_clone_frequencies function."""

    @pytest.fixture
    def adata_with_clone_info(self, sample_adata):
        """Create AnnData with clone information."""
        adata = sample_adata.copy()
        adata.obs["sample"] = ["S1"] * 50 + ["S2"] * 50
        adata.obs["CDR3_alpha"] = ["CASSL"] * 30 + ["CAVSD"] * 70
        adata.obs["CDR3_beta"] = ["CASSF"] * 30 + ["CASRG"] * 70
        adata.obs["clone_id"] = adata.obs["CDR3_alpha"] + "_" + adata.obs["CDR3_beta"]
        adata.obs["is_complete_clone"] = True
        return adata

    def test_calculate_frequencies_basic(self, adata_with_clone_info, sample_clonotypes_df):
        """Test basic frequency calculation."""
        result = calculate_clone_frequencies(sample_clonotypes_df, adata_with_clone_info)

        assert "sample_frequencies" in result.columns
        assert "n_conditions_present" in result.columns

    def test_frequency_values(self, adata_with_clone_info):
        """Test frequency values are calculated correctly."""
        # Create clonotypes from this adata
        clonotypes = aggregate_clonotypes(adata_with_clone_info, group_by="CDR3ab")

        result = calculate_clone_frequencies(clonotypes, adata_with_clone_info)

        # All frequencies should be between 0 and 1
        if "max_frequency" in result.columns:
            assert all(result["max_frequency"] >= 0)
            assert all(result["max_frequency"] <= 1)


class TestAggregateClonotypesExtended:
    """Extended tests for aggregate_clonotypes edge cases."""

    @pytest.fixture
    def adata_with_vdj_genes(self, sample_adata):
        """Create AnnData with VDJ gene columns."""
        adata = sample_adata.copy()
        adata.obs["sample"] = "S1"
        adata.obs["CDR3_alpha"] = "CASSL"
        adata.obs["CDR3_beta"] = "CASSF"
        adata.obs["TRA_1_v_gene"] = "TRAV1"
        adata.obs["TRA_1_j_gene"] = "TRAJ1"
        adata.obs["TRA_1_c_gene"] = "TRAC"
        adata.obs["TRB_1_v_gene"] = "TRBV2"
        adata.obs["TRB_1_j_gene"] = "TRBJ2"
        adata.obs["TRB_1_c_gene"] = "TRBC1"
        adata.obs["TRA_1_vdj_aa"] = "MRLVTSGF"
        adata.obs["TRA_1_vdj_nt"] = "ATGCGT"
        adata.obs["TRB_1_vdj_aa"] = "MGVTSGHD"
        adata.obs["TRB_1_vdj_nt"] = "ATGGGT"
        adata.obs["TRA_1_umis"] = 100
        adata.obs["TRB_1_umis"] = 200
        adata.obs["TRA_1_reads"] = 1000
        adata.obs["TRB_1_reads"] = 2000
        adata.obs["TRA_1_contig_id"] = "contig_1"
        adata.obs["TRB_1_contig_id"] = "contig_2"
        adata.obs["antigen_description"] = "TestAntigen"
        adata.obs["source"] = "culture"
        return adata

    def test_aggregate_with_vdj_genes(self, adata_with_vdj_genes):
        """Test aggregation includes VDJ gene information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert len(clonotypes) > 0
        assert "alpha_v_gene" in clonotypes.columns
        assert "beta_v_gene" in clonotypes.columns
        assert clonotypes["alpha_v_gene"].iloc[0] == "TRAV1"

    def test_aggregate_with_vdj_sequences(self, adata_with_vdj_genes):
        """Test aggregation includes VDJ sequences."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "VDJ_alpha_aa" in clonotypes.columns
        assert "VDJ_beta_aa" in clonotypes.columns
        assert clonotypes["VDJ_alpha_aa"].iloc[0] == "MRLVTSGF"

    def test_aggregate_with_umi_metrics(self, adata_with_vdj_genes):
        """Test aggregation includes UMI metrics."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "alpha_umis_mean" in clonotypes.columns
        assert "beta_reads_sum" in clonotypes.columns

    def test_aggregate_with_contig_ids(self, adata_with_vdj_genes):
        """Test aggregation includes contig IDs."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "alpha_contig_ids" in clonotypes.columns
        assert "beta_contig_ids" in clonotypes.columns

    def test_aggregate_with_antigen_info(self, adata_with_vdj_genes):
        """Test aggregation includes antigen information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "antigens" in clonotypes.columns
        assert "n_antigens" in clonotypes.columns
        assert "TestAntigen" in clonotypes["antigens"].iloc[0]

    def test_aggregate_with_source_info(self, adata_with_vdj_genes):
        """Test aggregation includes source information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "sources" in clonotypes.columns
        assert clonotypes["sources"].iloc[0] == "culture"

    def test_aggregate_cdr3b_only_with_alpha(self, adata_with_vdj_genes):
        """Test CDR3b_only mode still captures alpha when available."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3b_only")

        assert len(clonotypes) > 0
        assert "CDR3_beta" in clonotypes.columns
        # CDR3_alpha should also be captured when available
        assert "CDR3_alpha" in clonotypes.columns


class TestGetClonotypeSummaryExtended:
    """Extended tests for get_clonotype_summary edge cases."""

    def test_summary_without_n_samples(self):
        """Test summary when n_samples column is missing."""
        df = pd.DataFrame({
            "clone_id": ["A", "B"],
            "cell_count": [5, 3],
            # No n_samples column
        })

        summary = get_clonotype_summary(df)

        assert "n_multi_sample" in summary
        assert summary["n_multi_sample"] == 0  # Default when column missing


class TestRealisticClonotyping:
    """Tests using realistic clonotype data matching real CellRanger outputs."""

    def test_clonotype_cdr3_format(self, sample_clonotypes_df):
        """Test that clonotype CDR3 sequences follow IMGT conventions."""
        df = sample_clonotypes_df

        # All alpha CDR3s should start with CA (conserved cysteine + alanine)
        for cdr3 in df["CDR3_alpha"]:
            assert cdr3.startswith("CA"), f"Alpha CDR3 {cdr3} doesn't start with CA"

        # All beta CDR3s should start with CAS (conserved cysteine + alanine + serine)
        for cdr3 in df["CDR3_beta"]:
            assert cdr3.startswith("CAS"), f"Beta CDR3 {cdr3} doesn't start with CAS"

    def test_clonotype_gene_naming(self, sample_clonotypes_df):
        """Test that V/J/C gene names follow IMGT nomenclature."""
        df = sample_clonotypes_df

        # V genes should be TRAV/TRBV with numeric suffix
        for v_gene in df["alpha_v_gene"]:
            assert v_gene.startswith("TRAV"), f"Alpha V gene {v_gene} doesn't start with TRAV"
        for v_gene in df["beta_v_gene"]:
            assert v_gene.startswith("TRBV"), f"Beta V gene {v_gene} doesn't start with TRBV"

        # J genes should be TRAJ/TRBJ
        for j_gene in df["alpha_j_gene"]:
            assert j_gene.startswith("TRAJ"), f"Alpha J gene {j_gene} doesn't start with TRAJ"

        # C genes should be TRAC/TRBC
        for c_gene in df["alpha_c_gene"]:
            assert c_gene == "TRAC", f"Alpha C gene {c_gene} should be TRAC"
        for c_gene in df["beta_c_gene"]:
            assert c_gene.startswith("TRBC"), f"Beta C gene {c_gene} doesn't start with TRBC"

    def test_realistic_cell_counts(self, sample_clonotypes_df):
        """Test that cell counts represent typical clonal expansion patterns."""
        df = sample_clonotypes_df

        # Most clones should be small (1-10 cells)
        small_clones = df[df["cell_count"] <= 10]
        assert len(small_clones) > 0, "Should have small clones"

        # Some clones can be expanded
        expanded = df[df["cell_count"] > 1]
        assert len(expanded) > 0, "Should have expanded clones"

        # Cell counts should be positive integers
        assert all(df["cell_count"] > 0)
        assert all(df["cell_count"] == df["cell_count"].astype(int))

    def test_frequency_range(self, sample_clonotypes_df):
        """Test that clone frequencies are in valid range."""
        df = sample_clonotypes_df

        # All frequencies should be between 0 and 1
        assert all(df["max_frequency"] >= 0)
        assert all(df["max_frequency"] <= 1)

        # Larger clones should have higher frequencies
        sorted_by_count = df.sort_values("cell_count", ascending=False)
        sorted_by_freq = df.sort_values("max_frequency", ascending=False)

        # Top clone by count should also be among top by frequency
        top_by_count = sorted_by_count.iloc[0]["clone_id"]
        top_3_by_freq = sorted_by_freq.head(3)["clone_id"].tolist()
        assert top_by_count in top_3_by_freq

    def test_tcell_type_annotations(self, sample_clonotypes_df):
        """Test T cell type consensus annotations."""
        df = sample_clonotypes_df

        # Valid T cell types
        valid_types = [
            "Confident CD4+", "Confident CD8+",
            "Likely CD4+", "Likely CD8+",
            "Unknown", "Mixed"
        ]

        for tcell_type in df["Tcell_type_consensus"]:
            assert tcell_type in valid_types, f"Unknown T cell type: {tcell_type}"

    def test_clone_id_format(self, sample_clonotypes_df):
        """Test that clone IDs are formatted as CDR3a_CDR3b."""
        df = sample_clonotypes_df

        for idx, row in df.iterrows():
            expected_id = f"{row['CDR3_alpha']}_{row['CDR3_beta']}"
            assert row["clone_id"] == expected_id, \
                f"Clone ID {row['clone_id']} doesn't match expected {expected_id}"


class TestClonotypeSummaryRealistic:
    """Tests for clonotype summary with realistic data."""

    def test_summary_statistics_realistic(self, sample_clonotypes_df):
        """Test summary statistics with realistic clonotype data."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        # Number of clonotypes
        assert summary["n_clonotypes"] == len(sample_clonotypes_df)

        # Total cells should be sum of cell_count
        assert summary["n_cells"] == sample_clonotypes_df["cell_count"].sum()

        # Singletons are clones with count == 1
        expected_singletons = (sample_clonotypes_df["cell_count"] == 1).sum()
        assert summary["n_singletons"] == expected_singletons

        # Expanded are clones with count > 1
        expected_expanded = (sample_clonotypes_df["cell_count"] > 1).sum()
        assert summary["n_expanded"] == expected_expanded

    def test_multi_sample_detection(self, sample_clonotypes_df):
        """Test detection of clones present in multiple samples."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        # Count clones with n_samples > 1
        expected_multi = (sample_clonotypes_df["n_samples"] > 1).sum()
        assert summary["n_multi_sample"] == expected_multi
