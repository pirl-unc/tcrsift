"""
Tests for T cell phenotyping.
"""

import pytest
import pandas as pd
import numpy as np
import anndata as ad

from tcrsift.phenotype import (
    classify_tcell_type,
    phenotype_cells,
    filter_by_tcell_type,
    get_phenotype_summary,
    validate_phenotype_vs_expected,
    TCELL_TYPE_CATEGORIES,
)


class TestClassifyTcellType:
    """Tests for classify_tcell_type function."""

    def test_confident_cd8(self):
        """High CD8/CD4 ratio should give Confident CD8+."""
        # CD8 = 20, CD4 = 1 -> ratio > 3
        result = classify_tcell_type(cd4_expr=1, cd8_expr=20)
        assert result == "Confident CD8+"

    def test_confident_cd4(self):
        """High CD4/CD8 ratio should give Confident CD4+."""
        # CD4 = 20, CD8 = 1 -> ratio > 3
        result = classify_tcell_type(cd4_expr=20, cd8_expr=1)
        assert result == "Confident CD4+"

    def test_likely_cd8(self):
        """CD8 only (no CD4) with low ratio should give Likely CD8+."""
        # With +1 smoothing: cd4=0, cd8=1 -> (1+1)/(0+1) = 2 < 3 (default)
        # But cd8 > 0 and cd4 == 0, so "Likely CD8+"
        result = classify_tcell_type(cd4_expr=0, cd8_expr=1)
        assert result == "Likely CD8+"

    def test_likely_cd4(self):
        """CD4 only (no CD8) with low ratio should give Likely CD4+."""
        # With +1 smoothing: cd4=1, cd8=0 -> (1+1)/(0+1) = 2 < 3 (default)
        # But cd4 > 0 and cd8 == 0, so "Likely CD4+"
        result = classify_tcell_type(cd4_expr=1, cd8_expr=0)
        assert result == "Likely CD4+"

    def test_unknown_similar_levels(self):
        """Similar CD4/CD8 levels should give Unknown."""
        result = classify_tcell_type(cd4_expr=5, cd8_expr=5)
        assert result == "Unknown"

    def test_unknown_both_zero(self):
        """Both zero should give Unknown (via +1 smoothing)."""
        result = classify_tcell_type(cd4_expr=0, cd8_expr=0)
        assert result == "Unknown"

    def test_custom_ratio_threshold(self):
        """Custom ratio threshold should work."""
        # With ratio=2, CD8=5, CD4=2 -> (5+1)/(2+1) = 2.0, not > 2
        result = classify_tcell_type(cd4_expr=2, cd8_expr=5, cd4_cd8_ratio=2.0)
        assert result == "Unknown"

        # With ratio=1.5, same values -> 2.0 > 1.5
        result = classify_tcell_type(cd4_expr=2, cd8_expr=5, cd4_cd8_ratio=1.5)
        assert result == "Confident CD8+"


class TestPhenotypeCells:
    """Tests for phenotype_cells function."""

    def test_phenotype_adds_tcell_type(self, sample_adata):
        """phenotype_cells should add Tcell_type column."""
        # Add required columns
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]  # CD4 is 3rd from end
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]  # CD8A + CD8B

        result = phenotype_cells(sample_adata)
        assert "Tcell_type" in result.obs.columns
        assert result.obs["Tcell_type"].dtype.name == "category"

    def test_phenotype_adds_convenience_columns(self, sample_adata):
        """phenotype_cells should add is_CD4, is_CD8 columns."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]

        result = phenotype_cells(sample_adata)
        assert "is_CD8" in result.obs.columns
        assert "is_CD4" in result.obs.columns
        assert "is_confident" in result.obs.columns

    def test_phenotype_computes_cd8_from_components(self, sample_adata):
        """Should compute CD8 from CD8A + CD8B if CD8 missing."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8A"] = sample_adata.X[:, -2]
        sample_adata.obs["CD8B"] = sample_adata.X[:, -1]
        # Don't add CD8

        result = phenotype_cells(sample_adata)
        assert "Tcell_type" in result.obs.columns
        assert "CD8" in result.obs.columns  # Should be computed

    def test_phenotype_raises_without_required_columns(self, sample_adata):
        """Should raise if missing required columns."""
        # Don't add CD4 or CD8
        with pytest.raises(ValueError, match="Missing required columns"):
            phenotype_cells(sample_adata)

    def test_phenotype_cd3_filter(self, sample_adata):
        """CD3 filter should work when CD3 column present."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        sample_adata.obs["CD3"] = sample_adata.X[:, -6]  # CD3D

        result = phenotype_cells(sample_adata, min_cd3_reads=10)
        assert "filter:min_cd3" in result.obs.columns


class TestFilterByTcellType:
    """Tests for filter_by_tcell_type function."""

    @pytest.fixture
    def phenotyped_adata(self, sample_adata):
        """Create phenotyped AnnData."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        return phenotype_cells(sample_adata)

    def test_filter_cd8(self, phenotyped_adata):
        """Filter to CD8+ cells."""
        result = filter_by_tcell_type(phenotyped_adata, tcell_type="cd8")
        assert all(result.obs["is_CD8"])

    def test_filter_cd4(self, phenotyped_adata):
        """Filter to CD4+ cells."""
        result = filter_by_tcell_type(phenotyped_adata, tcell_type="cd4")
        assert all(result.obs["is_CD4"])

    def test_filter_both(self, phenotyped_adata):
        """Filter to both CD4+ and CD8+."""
        result = filter_by_tcell_type(phenotyped_adata, tcell_type="both")
        assert all(result.obs["is_CD4"] | result.obs["is_CD8"])

    def test_filter_invalid_type_raises(self, phenotyped_adata):
        """Invalid tcell_type should raise."""
        from tcrsift.validation import TCRsiftValidationError

        with pytest.raises(TCRsiftValidationError, match="Invalid tcell_type"):
            filter_by_tcell_type(phenotyped_adata, tcell_type="invalid")

    def test_filter_without_phenotype_raises(self, sample_adata):
        """Should raise if Tcell_type column missing."""
        from tcrsift.validation import TCRsiftValidationError

        with pytest.raises(TCRsiftValidationError, match="Tcell_type"):
            filter_by_tcell_type(sample_adata, tcell_type="cd8")


class TestGetPhenotypeSummary:
    """Tests for get_phenotype_summary function."""

    @pytest.fixture
    def phenotyped_adata(self, sample_adata):
        """Create phenotyped AnnData."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        return phenotype_cells(sample_adata)

    def test_summary_by_sample(self, phenotyped_adata):
        """Summary should have one row per sample."""
        summary = get_phenotype_summary(phenotyped_adata)
        assert isinstance(summary, pd.DataFrame)
        assert len(summary) == 2  # S1 and S2
        assert "sample" in summary.columns
        assert "total_cells" in summary.columns

    def test_summary_has_type_counts(self, phenotyped_adata):
        """Summary should have counts for each T cell type."""
        summary = get_phenotype_summary(phenotyped_adata)
        for tcell_type in TCELL_TYPE_CATEGORIES:
            assert tcell_type in summary.columns
            assert f"{tcell_type}_pct" in summary.columns

    def test_summary_without_phenotype_raises(self, sample_adata):
        """Should raise if Tcell_type column missing."""
        with pytest.raises(ValueError, match="must have Tcell_type column"):
            get_phenotype_summary(sample_adata)


class TestValidatePhenotypeVsExpected:
    """Tests for validate_phenotype_vs_expected function."""

    def test_no_warnings_when_no_expected(self, sample_adata):
        """No warnings if no expected_tcell_type column."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        phenotyped = phenotype_cells(sample_adata)

        warnings = validate_phenotype_vs_expected(phenotyped)
        assert len(warnings) == 0

    def test_warning_when_mismatch_expected_cd8(self, sample_adata):
        """Warning when expected CD8 but mostly CD4."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        phenotyped = phenotype_cells(sample_adata)

        # Set expected type
        phenotyped.obs["expected_tcell_type"] = "CD8"

        # Force mostly CD4
        phenotyped.obs["is_CD4"] = True
        phenotyped.obs["is_CD8"] = False

        warnings = validate_phenotype_vs_expected(phenotyped)
        assert any("Expected CD8+ but found" in w for w in warnings)

    def test_no_warning_for_mixed_expected(self, sample_adata):
        """No warning for mixed expected type."""
        sample_adata.obs["CD4"] = sample_adata.X[:, -3]
        sample_adata.obs["CD8"] = sample_adata.X[:, -2] + sample_adata.X[:, -1]
        phenotyped = phenotype_cells(sample_adata)

        phenotyped.obs["expected_tcell_type"] = "mixed"

        warnings = validate_phenotype_vs_expected(phenotyped)
        # Should not warn for mixed
        assert not any("Expected" in w for w in warnings)
