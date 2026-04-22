"""
Tests for clonotype filtering.
"""

import numpy as np
import pandas as pd
import pytest

from tcrsift.filter import (
    assign_tiers_threshold,
    filter_clonotypes,
    filter_clonotypes_logistic,
    filter_clonotypes_threshold,
    get_filter_summary,
    split_by_tier,
)
from tcrsift.model import LogisticFrequencyModel


class TestFilterClonotypesThreshold:
    """Tests for threshold-based filtering."""

    def test_filter_min_cells(self, sample_clonotypes_df):
        """Filter by minimum cell count."""
        # Cell counts: 15, 8, 3, 2, 1
        result = filter_clonotypes_threshold(sample_clonotypes_df, min_cells=5)
        assert len(result) == 2  # Only clones with 15 and 8 cells

        result = filter_clonotypes_threshold(sample_clonotypes_df, min_cells=10)
        assert len(result) == 1  # Only clone with 15 cells

    def test_filter_min_frequency(self, sample_clonotypes_df):
        """Filter by minimum frequency."""
        # Frequencies: 0.15, 0.08, 0.03, 0.02, 0.01
        result = filter_clonotypes_threshold(sample_clonotypes_df, min_frequency=0.05)
        assert len(result) == 2  # 0.15 and 0.08

    def test_filter_max_conditions(self, sample_clonotypes_df):
        """Filter by maximum conditions (samples)."""
        # n_samples: 2, 1, 1, 1, 1; cell_counts: 15, 8, 3, 2, 1
        # With min_cells=0, should exclude first clone (n_samples=2)
        result = filter_clonotypes_threshold(sample_clonotypes_df, max_conditions=1, min_cells=0)
        assert len(result) == 4  # Exclude first clone with n_samples=2

    def test_filter_max_conditions_prefers_n_conditions(self, sample_clonotypes_df):
        """True condition counts should be preferred over sample count proxies."""
        df = sample_clonotypes_df.copy()
        df["n_samples"] = [3, 1, 1, 1, 1]
        df["n_conditions"] = [1, 2, 1, 1, 1]

        result = filter_clonotypes_threshold(df, max_conditions=1, min_cells=0)

        assert "CAVSDGGSQGNLIF_CASSLGQAYEQYF" in set(result["CDR3ab"])
        assert "CAVSAGGSQGNLIF_CASSLGQAYEQYF" not in set(result["CDR3ab"])

    def test_filter_require_complete(self, sample_clonotypes_df):
        """Filter to require complete TCR."""
        # Add incomplete clone
        df = sample_clonotypes_df.copy()
        df.loc[len(df)] = {
            "CDR3ab": "incomplete",
            "CDR3_alpha": "",  # missing alpha
            "CDR3_beta": "CASSTEST",
            "cell_count": 10,
        }

        # Use min_cells=0 to not filter by cell count
        result = filter_clonotypes_threshold(df, require_complete=True, min_cells=0)
        assert len(result) == 5  # Original 5, incomplete excluded

        result = filter_clonotypes_threshold(df, require_complete=False, min_cells=0)
        assert len(result) == 6  # All 6 included

    def test_filter_tcell_type_cd8(self, sample_clonotypes_df):
        """Filter to CD8+ clones only."""
        # Tcell_type_consensus: 2x Confident CD8+, 1x Confident CD4+, 1x Likely CD8+, 1x Unknown
        result = filter_clonotypes_threshold(sample_clonotypes_df, tcell_type="cd8")
        assert len(result) == 3  # Both Confident CD8+ and Likely CD8+

    def test_filter_tcell_type_cd4(self, sample_clonotypes_df):
        """Filter to CD4+ clones only."""
        result = filter_clonotypes_threshold(sample_clonotypes_df, tcell_type="cd4")
        assert len(result) == 1  # Only Confident CD4+

    def test_filter_exclude_viral(self, sample_clonotypes_df):
        """Filter to exclude viral clones."""
        df = sample_clonotypes_df.copy()
        df["is_viral"] = [True, False, False, False, False]

        # Use min_cells=0 to not filter by cell count
        result = filter_clonotypes_threshold(df, exclude_viral=True, min_cells=0)
        assert len(result) == 4  # First clone excluded

        result = filter_clonotypes_threshold(df, exclude_viral=False, min_cells=0)
        assert len(result) == 5  # All included


class TestAssignTiersThreshold:
    """Tests for tier assignment using thresholds."""

    def test_tier_assignment_default(self, sample_clonotypes_df):
        """Test tier assignment with default thresholds."""
        result = assign_tiers_threshold(sample_clonotypes_df)

        assert "tier" in result.columns
        assert result["tier"].notna().any()

    def test_tier_assignment_custom(self, sample_clonotypes_df):
        """Test tier assignment with custom thresholds."""
        custom_tiers = {
            "tier1": {"min_cells": 10, "min_frequency": 0.1, "max_conditions": 2},
            "tier2": {"min_cells": 5, "min_frequency": 0.05, "max_conditions": 5},
        }

        result = assign_tiers_threshold(sample_clonotypes_df, tier_definitions=custom_tiers)

        # First clone: 15 cells, 0.15 freq, 2 conditions -> tier1
        tier1_clones = result[result["tier"] == "tier1"]
        assert len(tier1_clones) == 1
        assert tier1_clones["cell_count"].iloc[0] == 15

    def test_higher_tiers_override_lower(self, sample_clonotypes_df):
        """Higher quality tiers should override lower ones."""
        result = assign_tiers_threshold(sample_clonotypes_df)

        # Clone that qualifies for tier1 shouldn't be in tier2, tier3, etc.
        tier1_clones = result[result["tier"] == "tier1"]
        for CDR3ab in tier1_clones["CDR3ab"]:
            # This clone should only appear in tier1, not lower tiers
            other_tiers = result[(result["CDR3ab"] == CDR3ab) & (result["tier"] != "tier1")]
            assert len(other_tiers) == 0

    def test_tier_assignment_prefers_n_conditions(self, sample_clonotypes_df):
        """Tier assignment should use condition counts when available."""
        df = sample_clonotypes_df.copy()
        df["n_samples"] = [3, 1, 1, 1, 1]
        df["n_conditions"] = [1, 2, 1, 1, 1]
        custom_tiers = {
            "tier1": {"min_cells": 5, "min_frequency": 0.05, "max_conditions": 1},
            "tier2": {"min_cells": 1, "min_frequency": 0.0, "max_conditions": 2},
        }

        result = assign_tiers_threshold(df, tier_definitions=custom_tiers)
        by_cdr3 = result.set_index("CDR3ab")["tier"].to_dict()

        assert by_cdr3["CAVSDGGSQGNLIF_CASSLGQAYEQYF"] == "tier1"
        assert by_cdr3["CAVSAGGSQGNLIF_CASSLGQAYEQYF"] == "tier2"


class TestFilterClonotypesLogistic:
    """Tests for logistic regression filtering."""

    def test_logistic_returns_tiers(self, clonotypes_for_logistic):
        """Logistic method should return tier assignments."""
        result = filter_clonotypes_logistic(clonotypes_for_logistic)

        assert "tier" in result.columns
        # Verify model actually fit and assigned varied tiers
        tier_counts = result["tier"].value_counts()
        assert len(tier_counts) >= 2  # Should have multiple tiers

    def test_logistic_model_attributes(self, clonotypes_for_logistic):
        """Logistic method should store model info in attrs."""
        result = filter_clonotypes_logistic(clonotypes_for_logistic)

        # Model should fit successfully and store weight
        assert "logistic_model_weight" in result.attrs
        assert result.attrs["logistic_model_weight"] > 0  # Positive weight expected

    def test_logistic_custom_fdr_tiers(self, clonotypes_for_logistic):
        """Custom FDR tiers should work."""
        custom_fdr = [0.05, 0.1, 0.2]
        result = filter_clonotypes_logistic(clonotypes_for_logistic, fdr_tiers=custom_fdr)

        assert "tier" in result.columns
        # Should have at most len(custom_fdr) tiers
        unique_tiers = result["tier"].dropna().unique()
        assert len(unique_tiers) <= len(custom_fdr)

    def test_logistic_uses_default_threshold_on_fit_failure(self, clonotypes_for_logistic, monkeypatch):
        """Fit failures should fall back to the documented default frequency threshold."""

        def fail_fit(*_args, **_kwargs):
            raise ValueError("synthetic fit failure")

        monkeypatch.setattr("tcrsift.filter.fit_frequency_logistic_model", fail_fit)

        result = filter_clonotypes_logistic(
            clonotypes_for_logistic,
            fdr_tiers=[0.05, 0.1, 0.2],
            default_freq_threshold=0.4,
        )

        assert result.attrs["logistic_fallback_reason"] == "fit_failed"
        assert result.attrs["fdr_to_threshold"] == {0.05: 0.4, 0.1: 0.4, 0.2: 0.4}
        assert result.loc[result["max_frequency"] >= 0.4, "tier"].eq("tier1").all()

    def test_logistic_prefers_n_conditions_for_target(self, clonotypes_for_logistic, monkeypatch):
        """Strict logistic mode should use condition counts when available."""
        df = clonotypes_for_logistic.copy()
        df["n_samples"] = 3
        df["n_conditions"] = np.where(np.arange(len(df)) % 2 == 0, 1, 2)
        captured = {}

        def fake_fit(x, y, **_kwargs):
            captured["target"] = np.asarray(y)
            return LogisticFrequencyModel(
                params=np.asarray([10.0], dtype=float),
                converged=True,
                n_iter=2,
            )

        monkeypatch.setattr("tcrsift.filter.fit_frequency_logistic_model", fake_fit)

        filter_clonotypes_logistic(df, only_avoid_viral=False)

        expected_target = ((df["max_frequency"] > 0.09) & (df["n_conditions"] == 1)).values
        assert np.array_equal(captured["target"].astype(bool), expected_target)

    def test_logistic_fallback_on_small_data(self, sample_clonotypes_df):
        """Logistic method should gracefully handle small/degenerate data."""
        # Small dataset will cause model fitting issues - test that it falls back gracefully
        result = filter_clonotypes_logistic(sample_clonotypes_df)
        assert "tier" in result.columns


class TestFilterClonotypes:
    """Tests for main filter_clonotypes function."""

    def test_filter_threshold_method(self, sample_clonotypes_df):
        """Test threshold method through main function."""
        result = filter_clonotypes(
            sample_clonotypes_df,
            method="threshold",
            min_cells=2,
        )

        assert "tier" in result.columns

    def test_filter_logistic_method(self, clonotypes_for_logistic):
        """Test logistic method through main function."""
        result = filter_clonotypes(
            clonotypes_for_logistic,
            method="logistic",
            min_cells=1,
        )

        assert "tier" in result.columns
        # Should have meaningful tier distribution
        assert result["tier"].notna().sum() > 0

    def test_filter_combined_criteria(self, sample_clonotypes_df):
        """Test combining multiple filter criteria."""
        result = filter_clonotypes(
            sample_clonotypes_df,
            method="threshold",
            tcell_type="cd8",
            min_cells=3,
            min_frequency=0.01,
        )

        # Should have CD8+ clones with at least 3 cells and 1% frequency
        assert all(result["Tcell_type_consensus"].str.contains("CD8"))
        assert all(result["cell_count"] >= 3)
        assert all(result["max_frequency"] >= 0.01)


class TestSplitByTier:
    """Tests for split_by_tier function."""

    def test_split_returns_dict(self, sample_clonotypes_df):
        """Split should return dictionary of DataFrames."""
        tiered = assign_tiers_threshold(sample_clonotypes_df)
        split = split_by_tier(tiered)

        assert isinstance(split, dict)
        for tier_name, tier_df in split.items():
            assert isinstance(tier_df, pd.DataFrame)
            assert all(tier_df["tier"] == tier_name)

    def test_split_without_tier_raises(self, sample_clonotypes_df):
        """Should raise if no tier column."""
        with pytest.raises(ValueError, match="must have 'tier' column"):
            split_by_tier(sample_clonotypes_df)

    def test_split_excludes_none_tier(self, sample_clonotypes_df):
        """None tiers should be excluded from split."""
        df = sample_clonotypes_df.copy()
        df["tier"] = ["tier1", "tier2", None, "tier1", None]

        split = split_by_tier(df)

        assert None not in split
        assert "tier1" in split
        assert "tier2" in split


class TestGetFilterSummary:
    """Tests for get_filter_summary function."""

    def test_summary_without_tiers(self, sample_clonotypes_df):
        """Summary without tiers should have basic info."""
        summary = get_filter_summary(sample_clonotypes_df)

        assert "total_clonotypes" in summary
        assert summary["total_clonotypes"] == 5

    def test_summary_with_tiers(self, sample_clonotypes_df):
        """Summary with tiers should have tier-specific info."""
        tiered = assign_tiers_threshold(sample_clonotypes_df)
        summary = get_filter_summary(tiered)

        assert "total_clonotypes" in summary
        assert "tier_counts" in summary

        # Should have cell counts for each tier
        for tier in tiered["tier"].dropna().unique():
            assert f"{tier}_cells" in summary
