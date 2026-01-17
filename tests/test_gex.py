"""
Tests for GEX augmentation module.
"""

import numpy as np
import pandas as pd
import pytest

from tcrsift.gex import (
    DEFAULT_GENE_GROUPS,
    DEFAULT_GENE_LIST,
    aggregate_gex_by_clonotype,
    compute_cd4_cd8_counts,
)


class TestDefaultGeneLists:
    """Tests for default gene list constants."""

    def test_default_gene_list_contains_t_cell_markers(self):
        """Default gene list should contain key T cell markers."""
        assert "CD3D" in DEFAULT_GENE_LIST
        assert "CD4" in DEFAULT_GENE_LIST
        assert "CD8A" in DEFAULT_GENE_LIST
        assert "GZMA" in DEFAULT_GENE_LIST

    def test_default_gene_groups_structure(self):
        """Default gene groups should have expected structure."""
        assert "CD3" in DEFAULT_GENE_GROUPS
        assert "CD8" in DEFAULT_GENE_GROUPS
        assert "CD3D" in DEFAULT_GENE_GROUPS["CD3"]
        assert "CD8A" in DEFAULT_GENE_GROUPS["CD8"]


class TestAggregateGexByClonotype:
    """Tests for aggregate_gex_by_clonotype function."""

    @pytest.fixture
    def cells_with_gex(self):
        """Create cell DataFrame with GEX columns."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
            "barcode": ["bc1", "bc2", "bc3", "bc4", "bc5"],
            "gex.CD3": [10.0, 20.0, 30.0, 5.0, 15.0],
            "gex.CD4": [1.0, 2.0, 0.0, 10.0, 20.0],
            "gex.CD8": [15.0, 25.0, 35.0, 0.0, 0.0],
            "gex.n_reads": [1000, 2000, 1500, 800, 1200],
        })

    def test_aggregate_counts_cells(self, cells_with_gex):
        """Test that cell counts are computed correctly."""
        result = aggregate_gex_by_clonotype(cells_with_gex, verbose=False)
        assert len(result) == 2
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        cd_row = result[result["CDR3_pair"] == "C/D"].iloc[0]
        assert ab_row["total_cells.count"] == 3
        assert cd_row["total_cells.count"] == 2

    def test_aggregate_sum_operation(self, cells_with_gex):
        """Test sum aggregation operation."""
        result = aggregate_gex_by_clonotype(cells_with_gex, operations=["sum"], verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # CD3: 10 + 20 + 30 = 60
        assert ab_row["gex.CD3.sum"] == 60.0

    def test_aggregate_mean_operation(self, cells_with_gex):
        """Test mean aggregation operation."""
        result = aggregate_gex_by_clonotype(cells_with_gex, operations=["mean"], verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # CD3: (10 + 20 + 30) / 3 = 20
        assert ab_row["gex.CD3.mean"] == 20.0

    def test_aggregate_multiple_operations(self, cells_with_gex):
        """Test multiple aggregation operations."""
        result = aggregate_gex_by_clonotype(cells_with_gex, operations=["sum", "mean", "max"], verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert "gex.CD3.sum" in result.columns
        assert "gex.CD3.mean" in result.columns
        assert "gex.CD3.max" in result.columns
        assert ab_row["gex.CD3.max"] == 30.0

    def test_aggregate_custom_group_column(self, cells_with_gex):
        """Test aggregation with custom group column."""
        cells_with_gex["custom_clone"] = ["X", "X", "Y", "Y", "Y"]
        result = aggregate_gex_by_clonotype(cells_with_gex, group_col="custom_clone", verbose=False)
        assert len(result) == 2
        x_row = result[result["custom_clone"] == "X"].iloc[0]
        assert x_row["total_cells.count"] == 2

    def test_aggregate_handles_missing_values(self):
        """Test aggregation handles NaN values gracefully."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B", "A/B"],
            "gex.CD3": [10.0, np.nan, 30.0],
        })
        result = aggregate_gex_by_clonotype(df, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # Sum should be 40 (ignoring NaN)
        assert ab_row["gex.CD3.sum"] == 40.0


class TestComputeCd4Cd8Counts:
    """Tests for compute_cd4_cd8_counts function."""

    @pytest.fixture
    def cells_with_cd4_cd8(self):
        """Create cell DataFrame with CD4/CD8 expression."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B", "A/B", "A/B", "C/D", "C/D"],
            "gex.CD4": [10.0, 0.0, 5.0, 0.0, 15.0, 0.0],  # 2 CD4+
            "gex.CD8": [0.0, 20.0, 0.0, 25.0, 0.0, 10.0],  # 2 CD8+
        })

    def test_compute_cd4_only_count(self, cells_with_cd4_cd8):
        """Test computing CD4-only cell counts."""
        result = compute_cd4_cd8_counts(cells_with_cd4_cd8, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # Cells with CD4>0 and CD8=0: rows 0 and 2
        assert ab_row["CD4_only.count"] == 2

    def test_compute_cd8_only_count(self, cells_with_cd4_cd8):
        """Test computing CD8-only cell counts."""
        result = compute_cd4_cd8_counts(cells_with_cd4_cd8, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # Cells with CD8>0 and CD4=0: rows 1 and 3
        assert ab_row["CD8_only.count"] == 2

    def test_compute_total_cells(self, cells_with_cd4_cd8):
        """Test total cell count is correct."""
        result = compute_cd4_cd8_counts(cells_with_cd4_cd8, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["total_cells.count"] == 4

    def test_compute_with_custom_columns(self):
        """Test with custom CD4/CD8 column names."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B"],
            "my_cd4": [10.0, 0.0],
            "my_cd8": [0.0, 20.0],
        })
        result = compute_cd4_cd8_counts(df, cd4_col="my_cd4", cd8_col="my_cd8", verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["CD4_only.count"] == 1
        assert ab_row["CD8_only.count"] == 1

    def test_compute_handles_nan_values(self):
        """Test handling of NaN values in expression."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B", "A/B"],
            "gex.CD4": [10.0, np.nan, 0.0],
            "gex.CD8": [0.0, np.nan, 20.0],
        })
        result = compute_cd4_cd8_counts(df, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # NaN is treated as 0 via fillna
        assert ab_row["CD4_only.count"] == 1
        assert ab_row["CD8_only.count"] == 1

    def test_missing_cd4_column_raises(self):
        """Test that missing CD4 column raises error."""
        from tcrsift.validation import TCRsiftValidationError
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "gex.CD8": [10.0],
        })
        with pytest.raises(TCRsiftValidationError, match="Could not find CD4 expression column"):
            compute_cd4_cd8_counts(df, verbose=False)

    def test_missing_cd8_column_raises(self):
        """Test that missing CD8 column raises error."""
        from tcrsift.validation import TCRsiftValidationError
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "gex.CD4": [10.0],
        })
        with pytest.raises(TCRsiftValidationError, match="Could not find CD8 expression column"):
            compute_cd4_cd8_counts(df, verbose=False)

    def test_autodetect_column_names(self):
        """Test auto-detection of different column name formats."""
        # Test with 'CD4' format
        df1 = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "CD4": [10.0],
            "CD8": [0.0],
        })
        result1 = compute_cd4_cd8_counts(df1, gex_prefix="", verbose=False)
        assert result1["CD4_only.count"].iloc[0] == 1

        # Test with 'gex.CD4' format
        df2 = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "gex.CD4": [10.0],
            "gex.CD8": [0.0],
        })
        result2 = compute_cd4_cd8_counts(df2, verbose=False)
        assert result2["CD4_only.count"].iloc[0] == 1
