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

"""Tests for the rank-family plot bundle (#43, items #2/#3/#4/#5).

File-emission tests only — visual content isn't asserted. Each plot
gets a happy-path test and a skip-on-missing-column test, plus
narrower checks on the shared kernel (``_per_sample_rank_table``,
``draw_reference_fractions``).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import (  # noqa: E402
    REFERENCE_FRACTIONS,
    _per_sample_rank_table,
    draw_reference_fractions,
    plot_cumulative_coverage,
    plot_rank_curves_per_sample,
    plot_rank_family,
    plot_top_n_cumulative_stacked,
    plot_top_n_labeled_bars_per_sample,
)


@pytest.fixture
def ranked_clonotypes():
    """50 clones across 3 samples, decreasing-frequency, with a handful
    of tier-selected entries so tier-highlight code paths fire."""
    return pd.DataFrame(
        {
            "CDR3ab": [f"CASS{i:02d}F_CAVT{i:02d}F" for i in range(50)],
            "cell_count": list(range(50, 0, -1)),
            "tier": ["tier1"] * 3 + ["tier2"] * 5 + [None] * 42,
            "sample_frequencies": [
                {
                    "S1": max(0.0, 0.3 - i * 0.005),
                    "S2": max(0.0, 0.2 - i * 0.003),
                    "S3": max(0.0, 0.15 - i * 0.002),
                }
                for i in range(50)
            ],
        }
    )


class TestDrawReferenceFractions:
    """The shared log-y reference-fraction utility."""

    def test_adds_one_axhline_per_fraction_in_range(self):
        fig, ax = plt.subplots()
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 1.0)
        before = len(ax.lines)
        draw_reference_fractions(ax)
        plt.close(fig)
        # All five defaults are in [1e-5, 1.0] → 5 new lines.
        assert len(ax.lines) - before == len(REFERENCE_FRACTIONS)

    def test_skips_fractions_outside_ylim(self):
        fig, ax = plt.subplots()
        ax.set_yscale("log")
        # Only 10% and 1% fall in [0.001, 1] — three of the defaults
        # (0.1%, 0.01%, 0.001%) are below.
        ax.set_ylim(0.001, 1.0)
        before = len(ax.lines)
        draw_reference_fractions(ax)
        plt.close(fig)
        # 10%, 1%, 0.1% (=0.001, the edge — inclusive) are inside.
        # The strict-check is <=, so 0.1% counts.
        assert len(ax.lines) - before == 3


class TestPerSampleRankTable:
    """The kernel that all four plots consume."""

    def test_returns_dict_keyed_by_sample(self, ranked_clonotypes):
        tables = _per_sample_rank_table(ranked_clonotypes)
        assert set(tables) == {"S1", "S2", "S3"}

    def test_frequencies_sorted_descending(self, ranked_clonotypes):
        tables = _per_sample_rank_table(ranked_clonotypes)
        freqs = tables["S1"]["frequency"].values
        assert all(freqs[i] >= freqs[i + 1] for i in range(len(freqs) - 1))

    def test_rank_index_is_1_based(self, ranked_clonotypes):
        tables = _per_sample_rank_table(ranked_clonotypes)
        assert tables["S1"].index[0] == 1
        assert tables["S1"].index.name == "rank"

    def test_zero_frequency_clones_excluded(self, ranked_clonotypes):
        tables = _per_sample_rank_table(ranked_clonotypes)
        # In S1, frequencies > 0 for ranks 0..59 with the formula
        # max(0, 0.3 - i*0.005); rank 60 hits 0. Confirm strict-positive.
        assert (tables["S1"]["frequency"] > 0).all()

    def test_tier_column_propagates_when_present(self, ranked_clonotypes):
        tables = _per_sample_rank_table(ranked_clonotypes)
        # The top 3 clones in S1 should be tier-selected (tier1).
        top3 = tables["S1"].head(3)
        assert all(top3["tier"] == "tier1")

    def test_returns_empty_when_no_sample_frequencies(self):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert _per_sample_rank_table(df) == {}


class TestPlotRankCurvesPerSample:
    """Plot #2 — per-sample rank curves with tier highlight."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_rank_curves_per_sample(ranked_clonotypes, tmp_path)
        assert out is not None and out.exists() and out.suffix == ".png"
        assert out.name == "rank_curves_per_sample.png"

    def test_linear_x_variant_uses_distinct_filename(
        self, ranked_clonotypes, tmp_path
    ):
        """The #43 spec calls out a `linear-x` variant; orchestrator
        emits it under a separate name so both are visible side-by-side."""
        out = plot_rank_curves_per_sample(
            ranked_clonotypes,
            tmp_path,
            filename="rank_linear.png",
            x_scale="linear",
        )
        assert out is not None and out.name == "rank_linear.png"

    def test_skips_when_no_sample_frequencies(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_rank_curves_per_sample(df, tmp_path) is None
        assert list(tmp_path.glob("*.png")) == []

    def test_works_without_tier_column(self, ranked_clonotypes, tmp_path):
        """Tier highlight is optional — rank curve must render without it."""
        df = ranked_clonotypes.drop(columns=["tier"])
        out = plot_rank_curves_per_sample(df, tmp_path)
        assert out is not None and out.exists()


class TestPlotCumulativeCoverage:
    """Plot #3 — cumulative Σ(top-k) per sample."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_cumulative_coverage(ranked_clonotypes, tmp_path)
        assert out is not None and out.exists()

    def test_skips_when_no_sample_frequencies(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_cumulative_coverage(df, tmp_path) is None


class TestPlotTopNLabeledBarsPerSample:
    """Plot #4 — top-N labeled bar chart per sample."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_top_n_labeled_bars_per_sample(ranked_clonotypes, tmp_path, n=20)
        assert out is not None and out.exists()

    def test_handles_missing_cell_count_column(self, ranked_clonotypes, tmp_path):
        """Cell-count annotations on the top bars are optional —
        plot must still render when ``cell_count`` is absent."""
        df = ranked_clonotypes.drop(columns=["cell_count"])
        out = plot_top_n_labeled_bars_per_sample(df, tmp_path, n=10)
        assert out is not None and out.exists()

    def test_skips_when_no_sample_frequencies(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_top_n_labeled_bars_per_sample(df, tmp_path) is None


class TestPlotTopNCumulativeStacked:
    """Plot #5 — per-sample stacked bar of top-N contributions."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_top_n_cumulative_stacked(ranked_clonotypes, tmp_path, n=15)
        assert out is not None and out.exists()

    def test_skips_when_no_sample_frequencies(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_top_n_cumulative_stacked(df, tmp_path) is None


class TestPlotRankFamilyOrchestrator:
    """End-to-end orchestrator."""

    def test_emits_all_five_pngs(self, ranked_clonotypes, tmp_path):
        plot_rank_family(ranked_clonotypes, tmp_path, top_n=20)
        names = {p.name for p in tmp_path.glob("*.png")}
        assert names == {
            "rank_curves_per_sample.png",
            "rank_curves_per_sample_linear.png",
            "cumulative_coverage_per_sample.png",
            "top_n_labeled_bars_per_sample.png",
            "top_n_cumulative_stacked.png",
        }

    def test_safe_on_minimal_data(self, tmp_path):
        """Pre-aggregation data without ``sample_frequencies`` must
        not error — orchestrator just emits zero files."""
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        plot_rank_family(df, tmp_path)
        assert list(tmp_path.glob("*.png")) == []
