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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from tcrsift.model import (
    annotate_plot_with_thresholds_and_counts,
    calc_threshold,
    calc_thresholds_and_counts,
    count_at_threshold,
)


@pytest.fixture
def sample_clonotype_df():
    """Create sample clonotype DataFrame for testing."""
    np.random.seed(42)

    # Create realistic frequency distribution
    max_freq = np.concatenate(
        [
            np.random.uniform(0.01, 0.1, 40),  # Low frequency
            np.random.uniform(0.1, 0.5, 30),  # Medium frequency
            np.random.uniform(0.5, 1.0, 30),  # High frequency
        ]
    )

    # Create specificity descriptions
    specificity = []
    for f in max_freq:
        if f > 0.3:
            specificity.append("Single-culture")
        elif f > 0.1:
            specificity.append("Single-culture" if np.random.random() > 0.3 else "Viral")
        else:
            specificity.append("Viral" if np.random.random() > 0.5 else "Multi-culture")

    return pd.DataFrame(
        {
            "max_freq": max_freq,
            "specificity_description": specificity,
        }
    )


class TestCountAtThreshold:
    """Tests for count_at_threshold function."""

    def test_basic_count(self):
        """Counts Single-culture clones at/above the threshold.

        Uses an explicit hand-counted fixture (not the function's own masking
        logic recomputed as the expectation), so a regression in
        count_at_threshold can actually fail this test.
        """
        df = pd.DataFrame({
            "max_freq": [0.6, 0.4, 0.7, 0.9, 0.55, 0.5],
            "specificity_description": [
                "Single-culture",  # 0.6 ≥ 0.5 ✓
                "Single-culture",  # 0.4 < 0.5 ✗
                "Viral",           # 0.7 but not Single-culture ✗
                "Single-culture",  # 0.9 ≥ 0.5 ✓
                "Multi-culture",   # 0.55 but not Single-culture ✗
                "Single-culture",  # 0.5 exactly — boundary, ≥ is inclusive ✓
            ],
        })
        # 0.6, 0.9, and the 0.5 boundary row → 3. The boundary row makes this
        # sensitive to >= vs > (a > would drop it and give 2).
        assert count_at_threshold(df, 0.5) == 3

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
        empty_df = pd.DataFrame(
            {
                "max_freq": [],
                "specificity_description": [],
            }
        )

        count = count_at_threshold(empty_df, 0.5)

        assert count == 0


class TestCalcThreshold:
    """Tests for calc_threshold function."""

    def test_basic_threshold_calculation(self):
        """Test basic threshold calculation."""
        df = pd.DataFrame(
            {
                "max_freq": np.linspace(0, 1, 100),
                "specificity_description": ["Single-culture"] * 100,
            }
        )

        x_plot = np.linspace(0, 1, 1000)
        y_plot = x_plot  # Linear probability model

        threshold, n = calc_threshold(df, x_plot, y_plot, fdr=0.1)

        # At 10% FDR, target is 0.9 probability
        assert 0.85 <= threshold <= 0.95

    def test_min_value_constraint(self):
        """Test minimum value constraint."""
        df = pd.DataFrame(
            {
                "max_freq": np.linspace(0, 1, 100),
                "specificity_description": ["Single-culture"] * 100,
            }
        )

        x_plot = np.linspace(0, 1, 1000)
        y_plot = np.ones(1000) * 0.95  # Always 95% probability

        threshold, n = calc_threshold(df, x_plot, y_plot, fdr=0.1, min_value=0.5)

        # Should respect minimum value even if model suggests lower
        assert threshold >= 0.5

    def test_different_fdr_levels(self):
        """Test different FDR levels give different thresholds."""
        df = pd.DataFrame(
            {
                "max_freq": np.linspace(0, 1, 100),
                "specificity_description": ["Single-culture"] * 100,
            }
        )

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

    def test_default_threshold_when_weight_nonpositive(self):
        """Fallback branch: when frequency is ANTI-correlated with the target
        (so the logistic weight is ≤ 0), every FDR maps to default_freq_threshold.

        Previously this passed random data and only asserted ``0.1 in dict`` —
        which is true by construction and never exercised the fallback. Here the
        target deterministically decreases with max_freq, forcing the weight≤0
        branch, and we assert the actual default value is returned.
        """
        mf = np.linspace(0.05, 0.9, 40)
        # high frequency → Viral (non-target) → negative logistic weight
        spec = ["Viral" if f > 0.4 else "Single-culture" for f in mf]
        df = pd.DataFrame({"max_freq": mf, "specificity_description": spec})

        fdr_to_threshold, _counts, _model = calc_thresholds_and_counts(
            df, fdrs=[0.1, 0.01], default_freq_threshold=0.5, only_avoid_viral=False,
        )
        assert fdr_to_threshold[0.1] == 0.5
        assert fdr_to_threshold[0.01] == 0.5  # all FDRs fall back to the default

    def test_mode_flag_changes_the_threshold(self, sample_clonotype_df):
        """only_avoid_viral=True vs False must actually change the result.

        The two modes build a different positive target (non-Viral vs strictly
        Single-culture), which changes the fit and the derived threshold. The old
        pair of tests asserted only ``0.1 in dict`` for each mode, so a no-op
        flag would have passed. Here we assert the modes disagree and both
        respect the min-frequency floor.
        """
        floor = 0.09
        avoid_viral, _c1, _m1 = calc_thresholds_and_counts(
            sample_clonotype_df, fdrs=[0.1], only_avoid_viral=True,
            min_freq_threshold=floor,
        )
        strict, _c2, _m2 = calc_thresholds_and_counts(
            sample_clonotype_df, fdrs=[0.1], only_avoid_viral=False,
            min_freq_threshold=floor,
        )
        assert avoid_viral[0.1] != strict[0.1]          # the flag has an effect
        assert avoid_viral[0.1] >= floor and strict[0.1] >= floor

    def test_lower_fdr_gives_higher_or_equal_threshold(self, sample_clonotype_df):
        """Across multiple FDRs, the threshold is monotonic non-increasing as the
        FDR increases (a stricter FDR demands a higher frequency cutoff). The old
        test only asserted each requested FDR was a key — guaranteed by
        construction — and never checked this ordering property.
        """
        fdrs = [0.15, 0.1, 0.01, 0.001]
        fdr_to_threshold, _counts, _model = calc_thresholds_and_counts(
            sample_clonotype_df, fdrs=fdrs,
        )
        thresholds = [fdr_to_threshold[f] for f in fdrs]  # by increasing strictness
        # increasing strictness (decreasing fdr) → non-decreasing threshold
        assert thresholds == sorted(thresholds), thresholds

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
        assert hasattr(model, "params")


class TestAnnotatePlotWithThresholdsAndCounts:
    """Tests for threshold annotation plotting."""

    def test_plot_draws_one_curve_and_one_marker_per_threshold(self, sample_clonotype_df):
        """Plot annotation should not duplicate the curve or threshold markers."""
        fig, ax = plt.subplots()
        try:
            fdr_to_threshold, threshold_to_count, model = calc_thresholds_and_counts(
                sample_clonotype_df,
                fdrs=[0.1, 0.01],
            )

            annotate_plot_with_thresholds_and_counts(
                sample_clonotype_df,
                ax,
                model,
                fdr_to_threshold,
                threshold_to_count,
            )

            assert len(ax.lines) == 3  # 1 probability curve + 2 threshold lines
            assert len(ax.texts) == 2
        finally:
            plt.close(fig)
