"""
Tests for multi-experiment unification module.
"""

import numpy as np
import pandas as pd
import pytest

from tcrsift.unify import (
    add_phenotype_confidence,
    compute_condition_statistics,
    find_top_condition,
    get_unify_summary,
    merge_experiments,
)


class TestMergeExperiments:
    """Tests for merge_experiments function."""

    @pytest.fixture
    def til_clonotypes(self):
        """TIL clonotype DataFrame."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "E/F"],
            "CDR3_alpha": ["A", "C", "E"],
            "CDR3_beta": ["B", "D", "F"],
            "total_cells.count": [10, 5, 3],
            "CD4_only.count": [2, 4, 0],
            "CD8_only.count": [8, 1, 3],
        })

    @pytest.fixture
    def culture_clonotypes(self):
        """Culture clonotype DataFrame."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "G/H"],
            "CDR3_alpha": ["A", "C", "G"],
            "CDR3_beta": ["B", "D", "H"],
            "total_cells.count": [20, 15, 8],
            "CD4_only.count": [5, 10, 2],
            "CD8_only.count": [15, 5, 6],
        })

    def test_merge_basic_outer_join(self, til_clonotypes, culture_clonotypes):
        """Test that merge creates outer join of all clonotypes."""
        experiments = [(til_clonotypes, "TIL"), (culture_clonotypes, "Culture")]
        result = merge_experiments(experiments, verbose=False)

        # Should have 4 unique CDR3_pairs: A/B, C/D, E/F, G/H
        assert len(result) == 4
        assert set(result["CDR3_pair"]) == {"A/B", "C/D", "E/F", "G/H"}

    def test_merge_prefixes_columns(self, til_clonotypes, culture_clonotypes):
        """Test that columns are prefixed with experiment names."""
        experiments = [(til_clonotypes, "TIL"), (culture_clonotypes, "Culture")]
        result = merge_experiments(experiments, verbose=False)

        assert "TIL.total_cells.count" in result.columns
        assert "Culture.total_cells.count" in result.columns

    def test_merge_occurrence_flags(self, til_clonotypes, culture_clonotypes):
        """Test occurrence flag columns are added."""
        experiments = [(til_clonotypes, "TIL"), (culture_clonotypes, "Culture")]
        result = merge_experiments(experiments, add_occurrence_flags=True, verbose=False)

        assert "occurs_in_TIL" in result.columns
        assert "occurs_in_Culture" in result.columns

        # A/B is in both
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["occurs_in_TIL"] == True
        assert ab_row["occurs_in_Culture"] == True

        # E/F is only in TIL
        ef_row = result[result["CDR3_pair"] == "E/F"].iloc[0]
        assert ef_row["occurs_in_TIL"] == True
        assert ef_row["occurs_in_Culture"] == False

        # G/H is only in Culture
        gh_row = result[result["CDR3_pair"] == "G/H"].iloc[0]
        assert gh_row["occurs_in_TIL"] == False
        assert gh_row["occurs_in_Culture"] == True

    def test_merge_combined_stats(self, til_clonotypes, culture_clonotypes):
        """Test combined statistics are computed."""
        experiments = [(til_clonotypes, "TIL"), (culture_clonotypes, "Culture")]
        result = merge_experiments(experiments, add_combined_stats=True, verbose=False)

        assert "combined.total_cells.count" in result.columns

        # A/B: 10 (TIL) + 20 (Culture) = 30
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["combined.total_cells.count"] == 30

    def test_merge_adds_paired_flags(self, til_clonotypes, culture_clonotypes):
        """Test that paired/unpaired flags are added."""
        experiments = [(til_clonotypes, "TIL"), (culture_clonotypes, "Culture")]
        result = merge_experiments(experiments, verbose=False)

        assert "has_CDR3_alpha" in result.columns
        assert "has_CDR3_beta" in result.columns
        assert "is_paired" in result.columns

    def test_merge_empty_experiments_raises(self):
        """Test that empty experiment list raises error."""
        from tcrsift.validation import TCRsiftValidationError
        with pytest.raises(TCRsiftValidationError, match="No experiments provided"):
            merge_experiments([])

    def test_merge_missing_cdr3_pair_raises(self):
        """Test that missing CDR3_pair column raises error."""
        from tcrsift.validation import TCRsiftValidationError
        df = pd.DataFrame({"other": [1, 2]})
        with pytest.raises(TCRsiftValidationError, match="missing CDR3_pair column"):
            merge_experiments([(df, "Test")])


class TestAddPhenotypeConfidence:
    """Tests for add_phenotype_confidence function."""

    @pytest.fixture
    def unified_df(self):
        """Create unified DataFrame with GEX columns."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "E/F", "G/H"],
            "combined.gex.CD4.sum": [0.0, 100.0, 10.0, 50.0],
            "combined.gex.CD8.sum": [50.0, 0.0, 10.0, 0.0],
        })

    def test_confident_cd8_when_only_cd8(self, unified_df):
        """Test confident CD8 classification when only CD8 expression."""
        result = add_phenotype_confidence(unified_df, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["Confident_CD8"] == True
        assert ab_row["Confident_CD4"] == False

    def test_confident_cd4_when_only_cd4(self, unified_df):
        """Test confident CD4 classification when only CD4 expression."""
        result = add_phenotype_confidence(unified_df, verbose=False)
        cd_row = result[result["CDR3_pair"] == "C/D"].iloc[0]
        assert cd_row["Confident_CD4"] == True
        assert cd_row["Confident_CD8"] == False

    def test_not_confident_when_similar_expression(self, unified_df):
        """Test no confident call when similar CD4/CD8 expression."""
        result = add_phenotype_confidence(unified_df, verbose=False)
        ef_row = result[result["CDR3_pair"] == "E/F"].iloc[0]
        # Equal expression should not be confident either way
        assert ef_row["Confident_CD4"] == False
        assert ef_row["Confident_CD8"] == False

    def test_ratio_threshold_affects_confidence(self):
        """Test that ratio threshold parameter affects classification."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "combined.gex.CD4.sum": [10.0],
            "combined.gex.CD8.sum": [50.0],
        })
        # With ratio_threshold=10, 50 > 10*(1+10) = 110 is False
        result1 = add_phenotype_confidence(df.copy(), ratio_threshold=10.0, verbose=False)
        assert result1["Confident_CD8"].iloc[0] == False

        # With ratio_threshold=3, 50 > 3*(1+10) = 33 is True
        result2 = add_phenotype_confidence(df.copy(), ratio_threshold=3.0, verbose=False)
        assert result2["Confident_CD8"].iloc[0] == True

    def test_likely_columns_added(self, unified_df):
        """Test that Likely_CD4 and Likely_CD8 columns are added."""
        result = add_phenotype_confidence(unified_df, verbose=False)
        assert "Likely_CD4" in result.columns
        assert "Likely_CD8" in result.columns

    def test_handles_missing_gex_columns(self):
        """Test handling when GEX columns are missing."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "other_col": [1],
        })
        result = add_phenotype_confidence(df, verbose=False)
        # Should run without error, using default 0 values
        assert "Confident_CD4" in result.columns
        assert "Confident_CD8" in result.columns


class TestComputeConditionStatistics:
    """Tests for compute_condition_statistics function."""

    @pytest.fixture
    def df_with_conditions(self):
        """Create DataFrame with condition fraction columns."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
            "TIL.condition_pool1.frac": [0.1, 0.2],
            "Culture.condition_pool1.frac": [0.15, 0.25],
            "TIL.condition_pool2.frac": [0.05, 0.1],
            "Culture.condition_pool2.frac": [0.05, 0.15],
        })

    def test_computes_sum_across_sources(self, df_with_conditions):
        """Test sum computation across sources."""
        result = compute_condition_statistics(
            df_with_conditions, ["pool1", "pool2"], verbose=False
        )
        assert "combined.pool1.frac.sum" in result.columns

        # A/B: 0.1 + 0.15 = 0.25
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["combined.pool1.frac.sum"] == pytest.approx(0.25)

    def test_computes_mean_across_sources(self, df_with_conditions):
        """Test mean computation across sources."""
        result = compute_condition_statistics(
            df_with_conditions, ["pool1", "pool2"], verbose=False
        )
        assert "combined.pool1.frac.mean" in result.columns

        # A/B: (0.1 + 0.15) / 2 = 0.125
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["combined.pool1.frac.mean"] == pytest.approx(0.125)

    def test_computes_std_and_cv(self, df_with_conditions):
        """Test standard deviation and CV computation."""
        result = compute_condition_statistics(
            df_with_conditions, ["pool1", "pool2"], verbose=False
        )
        assert "combined.pool1.frac.std" in result.columns
        assert "combined.pool1.frac.cv" in result.columns

    def test_handles_missing_condition_columns(self):
        """Test handling when condition columns don't exist."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "other_col": [1],
        })
        result = compute_condition_statistics(df, ["nonexistent"], verbose=False)
        # Should run without error
        assert "combined.nonexistent.frac.sum" not in result.columns


class TestFindTopCondition:
    """Tests for find_top_condition function."""

    @pytest.fixture
    def df_with_combined_stats(self):
        """Create DataFrame with combined condition statistics."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "E/F"],
            "combined.pool1.frac.mean": [0.5, 0.1, 0.3],
            "combined.pool2.frac.mean": [0.2, 0.4, 0.3],
            "combined.pool3.frac.mean": [0.1, 0.1, 0.3],
        })

    def test_finds_top_condition(self, df_with_combined_stats):
        """Test finding the top condition for each clonotype."""
        result = find_top_condition(
            df_with_combined_stats, ["pool1", "pool2", "pool3"], verbose=False
        )
        assert "top_condition" in result.columns

        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["top_condition"] == "pool1"

        cd_row = result[result["CDR3_pair"] == "C/D"].iloc[0]
        assert cd_row["top_condition"] == "pool2"

    def test_adds_top_condition_value(self, df_with_combined_stats):
        """Test that top condition value is added."""
        result = find_top_condition(
            df_with_combined_stats, ["pool1", "pool2", "pool3"], verbose=False
        )
        assert "top_condition.value" in result.columns

        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        assert ab_row["top_condition.value"] == 0.5

    def test_require_clear_winner_with_tie(self, df_with_combined_stats):
        """Test that ties don't get assigned when require_clear_winner=True."""
        result = find_top_condition(
            df_with_combined_stats, ["pool1", "pool2", "pool3"],
            require_clear_winner=True, verbose=False
        )
        # E/F has equal values (0.3, 0.3, 0.3) - should be None
        ef_row = result[result["CDR3_pair"] == "E/F"].iloc[0]
        assert pd.isna(ef_row["top_condition"])

    def test_no_clear_winner_allows_assignment(self, df_with_combined_stats):
        """Test that ties can be assigned when require_clear_winner=False."""
        result = find_top_condition(
            df_with_combined_stats, ["pool1", "pool2", "pool3"],
            require_clear_winner=False, verbose=False
        )
        ef_row = result[result["CDR3_pair"] == "E/F"].iloc[0]
        # Should get some assignment (first one in sorted order)
        assert pd.notna(ef_row["top_condition"])

    def test_handles_missing_columns(self):
        """Test handling when condition columns don't exist."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B"],
            "other_col": [1],
        })
        result = find_top_condition(df, ["pool1", "pool2"], verbose=False)
        # Should run without error
        assert "top_condition" not in result.columns or result["top_condition"].isna().all()


class TestGetUnifySummary:
    """Tests for get_unify_summary function."""

    @pytest.fixture
    def unified_result(self):
        """Create a unified result DataFrame."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "E/F"],
            "is_paired": [True, True, False],
            "occurs_in_TIL": [True, False, True],
            "occurs_in_Culture": [True, True, False],
            "combined.total_cells.count": [30, 20, 5],
            "Confident_CD4": [False, True, False],
            "Confident_CD8": [True, False, False],
            "Likely_CD4": [False, True, False],
            "Likely_CD8": [True, False, True],
        })

    def test_summary_total_clonotypes(self, unified_result):
        """Test total clonotypes count."""
        summary = get_unify_summary(unified_result)
        assert summary["total_clonotypes"] == 3

    def test_summary_paired_clonotypes(self, unified_result):
        """Test paired clonotypes count."""
        summary = get_unify_summary(unified_result)
        assert summary["paired_clonotypes"] == 2

    def test_summary_source_counts(self, unified_result):
        """Test source occurrence counts."""
        summary = get_unify_summary(unified_result)
        assert summary["from_TIL"] == 2
        assert summary["from_Culture"] == 2

    def test_summary_total_cells(self, unified_result):
        """Test total cell count."""
        summary = get_unify_summary(unified_result)
        assert summary["total_cells"] == 55

    def test_summary_phenotype_counts(self, unified_result):
        """Test phenotype classification counts."""
        summary = get_unify_summary(unified_result)
        assert summary["Confident_CD4"] == 1
        assert summary["Confident_CD8"] == 1
        assert summary["Likely_CD4"] == 1
        assert summary["Likely_CD8"] == 2

    def test_summary_handles_missing_columns(self):
        """Test handling of missing columns."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
        })
        summary = get_unify_summary(df)
        assert summary["total_clonotypes"] == 2
        assert summary["paired_clonotypes"] == 0  # No is_paired column
