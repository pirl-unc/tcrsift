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

"""Tests for model threshold calculation module."""

import numpy as np
import pandas as pd
import pytest

from tcrsift.model import (
    calc_threshold,
    calc_thresholds_and_counts,
    count_at_threshold,
)


@pytest.fixture
def sample_clonotype_df():
    """Create sample clonotype DataFrame for testing."""
    np.random.seed(42)
    n = 100

    # Create realistic frequency distribution
    max_freq = np.concatenate([
        np.random.uniform(0.01, 0.1, 40),   # Low frequency
        np.random.uniform(0.1, 0.5, 30),    # Medium frequency
        np.random.uniform(0.5, 1.0, 30),    # High frequency
    ])

    # Create specificity descriptions
    specificity = []
    for f in max_freq:
        if f > 0.3:
            specificity.append("Single-culture")
        elif f > 0.1:
            specificity.append("Single-culture" if np.random.random() > 0.3 else "Viral")
        else:
            specificity.append("Viral" if np.random.random() > 0.5 else "Multi-culture")

    return pd.DataFrame({
        "max_freq": max_freq,
        "specificity_description": specificity,
    })


class TestCountAtThreshold:
    """Tests for count_at_threshold function."""

    def test_basic_count(self, sample_clonotype_df):
        """Test basic counting at threshold."""
        count = count_at_threshold(sample_clonotype_df, 0.5)

        # Should count Single-culture clones above threshold
        expected = (
            (sample_clonotype_df["specificity_description"] == "Single-culture") &
            (sample_clonotype_df["max_freq"] >= 0.5)
        ).sum()

        assert count == expected

    def test_zero_threshold(self, sample_clonotype_df):
        """Test with zero threshold (count all)."""
        count = count_at_threshold(sample_clonotype_df, 0.0)

        # Should count all Single-culture clones
        expected = (sample_clonotype_df["specificity_description"] == "Single-culture").sum()

        assert count == expected

    def test_high_threshold(self, sample_clonotype_df):
        """Test with very high threshold."""
        count = count_at_threshold(sample_clonotype_df, 1.1)

        # No frequencies should be >= 1.1
        assert count == 0

    def test_empty_dataframe(self):
        """Test with empty DataFrame."""
        empty_df = pd.DataFrame({
            "max_freq": [],
            "specificity_description": [],
        })

        count = count_at_threshold(empty_df, 0.5)

        assert count == 0


class TestCalcThreshold:
    """Tests for calc_threshold function."""

    def test_basic_threshold_calculation(self):
        """Test basic threshold calculation."""
        df = pd.DataFrame({
            "max_freq": np.linspace(0, 1, 100),
            "specificity_description": ["Single-culture"] * 100,
        })

        x_plot = np.linspace(0, 1, 1000)
        y_plot = x_plot  # Linear probability model

        threshold, n = calc_threshold(df, x_plot, y_plot, fdr=0.1)

        # At 10% FDR, target is 0.9 probability
        assert 0.85 <= threshold <= 0.95

    def test_min_value_constraint(self):
        """Test minimum value constraint."""
        df = pd.DataFrame({
            "max_freq": np.linspace(0, 1, 100),
            "specificity_description": ["Single-culture"] * 100,
        })

        x_plot = np.linspace(0, 1, 1000)
        y_plot = np.ones(1000) * 0.95  # Always 95% probability

        threshold, n = calc_threshold(df, x_plot, y_plot, fdr=0.1, min_value=0.5)

        # Should respect minimum value even if model suggests lower
        assert threshold >= 0.5

    def test_different_fdr_levels(self):
        """Test different FDR levels give different thresholds."""
        df = pd.DataFrame({
            "max_freq": np.linspace(0, 1, 100),
            "specificity_description": ["Single-culture"] * 100,
        })

        x_plot = np.linspace(0, 1, 1000)
        y_plot = x_plot  # Linear probability

        threshold_low_fdr, _ = calc_threshold(df, x_plot, y_plot, fdr=0.01)
        threshold_high_fdr, _ = calc_threshold(df, x_plot, y_plot, fdr=0.1)

        # Lower FDR should require higher threshold
        assert threshold_low_fdr >= threshold_high_fdr


class TestCalcThresholdsAndCounts:
    """Tests for calc_thresholds_and_counts function."""

    def test_basic_threshold_calculation(self, sample_clonotype_df):
        """Test basic threshold and count calculation."""
        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=[0.1, 0.01],
            min_freq_threshold=0.05,
        )

        assert 0.1 in fdr_to_threshold
        assert 0.01 in fdr_to_threshold
        assert len(threshold_to_count) > 0
        assert model is not None

    def test_default_threshold_on_messy_data(self):
        """Test fallback to default threshold when data is messy."""
        # Create messy data where model would fail
        messy_df = pd.DataFrame({
            "max_freq": np.random.uniform(0, 1, 50),
            "specificity_description": np.random.choice(
                ["Single-culture", "Viral", "Multi-culture"],
                50
            ),
        })

        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            messy_df,
            fdrs=[0.1],
            default_freq_threshold=0.5,
        )

        # Should return some threshold
        assert 0.1 in fdr_to_threshold

    def test_only_avoid_viral_mode(self, sample_clonotype_df):
        """Test only_avoid_viral mode."""
        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=[0.1],
            only_avoid_viral=True,
        )

        assert 0.1 in fdr_to_threshold

    def test_strict_single_culture_mode(self, sample_clonotype_df):
        """Test strict single-culture mode."""
        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=[0.1],
            only_avoid_viral=False,
        )

        assert 0.1 in fdr_to_threshold

    def test_multiple_fdrs(self, sample_clonotype_df):
        """Test multiple FDR levels."""
        fdrs = [0.15, 0.1, 0.01, 0.001]

        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=fdrs,
        )

        for fdr in fdrs:
            assert fdr in fdr_to_threshold

    def test_min_freq_threshold_respected(self, sample_clonotype_df):
        """Test minimum frequency threshold is respected."""
        min_thresh = 0.2

        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=[0.1],
            min_freq_threshold=min_thresh,
        )

        # All thresholds should be >= min_freq_threshold
        for threshold in fdr_to_threshold.values():
            assert threshold >= min_thresh

    def test_returns_model(self, sample_clonotype_df):
        """Test that logistic model is returned."""
        fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
            sample_clonotype_df,
            fdrs=[0.1],
        )

        # Model should have params
        assert hasattr(model, 'params')
