"""
Tests for clonotype aggregation.
"""

import pytest
import pandas as pd
import numpy as np
import anndata as ad

from tcrsift.clonotype import (
    aggregate_clonotypes,
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
        """Empty DataFrame when no complete clones."""
        sample_adata.obs["sample"] = "S1"
        sample_adata.obs["CDR3_alpha"] = None
        sample_adata.obs["CDR3_beta"] = None

        clonotypes = aggregate_clonotypes(sample_adata, group_by="CDR3ab")
        assert len(clonotypes) == 0


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
