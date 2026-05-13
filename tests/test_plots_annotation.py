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

"""Tests for the structured-annotation plot family (#48 phase B).

Each test asserts file emission and graceful-skip behavior — not visual
content. The plots are exercised end-to-end against a synthetic
clonotype frame that mirrors what `match_clonotypes` produces after
#48 phase A.
"""

from __future__ import annotations

import matplotlib

# Headless backend for CI — needs to be set before pyplot imports.
matplotlib.use("Agg")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import (  # noqa: E402
    GRANULARITY_COLUMNS,
    _expand_sample_frequencies,
    _top_n_clones,
    plot_annotations_phase_b,
    plot_category_composition_bar,
    plot_match_strength_comparison,
    plot_matched_clone_heatmap,
)


@pytest.fixture
def annotated_clonotypes():
    """Synthetic post-annotation clonotypes — five clones across three
    samples, mix of categories and match strengths."""
    return pd.DataFrame(
        {
            "CDR3ab": [f"CASS{i}_CAVT{i}" for i in range(5)],
            "cell_count": [120, 80, 60, 40, 20],
            "max_frequency_per_method": [0.4, 0.3, 0.2, 0.1, 0.05],
            "sample_frequencies": [
                {"S1": 0.4, "S2": 0.1, "S3": 0.0},
                {"S1": 0.1, "S2": 0.3, "S3": 0.05},
                {"S1": 0.0, "S2": 0.2, "S3": 0.15},
                {"S1": 0.05, "S2": 0.05, "S3": 0.1},
                {"S1": 0.0, "S2": 0.0, "S3": 0.05},
            ],
            "db_match": [True, True, True, True, False],
            "db_category": ["viral", "tumor_self", "viral", "self", None],
            "db_protein_canonical": ["CMV pp65", "MART-1", "Flu M1", "Insulin", None],
            "db_epitope": ["NLVPMVATV", "ELAGIGILTV", "GILGFVFTL", "FFRKDLEEK", None],
            "db_mhc": ["HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01", None],
            "db_species": [
                "Human cytomegalovirus",
                "Homo sapiens",
                "Influenza A",
                "Homo sapiens",
                None,
            ],
            "db_match_strength": ["ab", "b_only", "ab", "b_only", None],
            "is_viral": [True, False, True, False, False],
        }
    )


class TestExpandSampleFrequencies:
    """The pivot helper that turns the per-clone dict into a wide frame."""

    def test_pivots_dict_to_wide_table(self, annotated_clonotypes):
        wide = _expand_sample_frequencies(annotated_clonotypes)
        assert list(wide.columns) == ["S1", "S2", "S3"]
        assert wide.shape == (5, 3)
        # Absent samples zero-filled, not NaN — downstream plots use this.
        assert wide.loc["CASS0_CAVT0", "S3"] == 0.0

    def test_returns_empty_when_no_frequencies_column(self):
        df = pd.DataFrame({"CDR3ab": ["X"]})
        assert _expand_sample_frequencies(df).empty


class TestTopNClones:
    """Top-N ranker — falls back through ranking columns gracefully."""

    def test_picks_top_by_max_frequency_per_method(self, annotated_clonotypes):
        top3 = _top_n_clones(annotated_clonotypes, n=3)
        assert list(top3["CDR3ab"]) == [
            "CASS0_CAVT0",
            "CASS1_CAVT1",
            "CASS2_CAVT2",
        ]

    def test_falls_back_to_cell_count_when_no_frequency(self, annotated_clonotypes):
        df = annotated_clonotypes.drop(
            columns=["max_frequency_per_method"]
        )
        top2 = _top_n_clones(df, n=2)
        assert list(top2["CDR3ab"]) == ["CASS0_CAVT0", "CASS1_CAVT1"]


class TestPlotMatchedCloneHeatmap:
    """The headline top-N × conditions heatmap with strip annotations."""

    def test_emits_png(self, annotated_clonotypes, tmp_path):
        out = plot_matched_clone_heatmap(
            annotated_clonotypes, tmp_path, top_n=4, granularity="category"
        )
        assert out is not None and out.exists() and out.suffix == ".png"
        # Filename carries the granularity for the CLI loop.
        assert "category" in out.name

    def test_skips_when_no_matched_clones(self, annotated_clonotypes, tmp_path):
        unmatched = annotated_clonotypes.copy()
        unmatched["db_match"] = False
        out = plot_matched_clone_heatmap(unmatched, tmp_path)
        assert out is None
        assert list(tmp_path.glob("*.png")) == []

    def test_skips_when_no_sample_frequencies(self, annotated_clonotypes, tmp_path):
        df = annotated_clonotypes.drop(columns=["sample_frequencies"])
        out = plot_matched_clone_heatmap(df, tmp_path)
        assert out is None


class TestPlotCategoryCompositionBar:
    """Per-condition stacked composition bar."""

    def test_emits_png_at_category_granularity(self, annotated_clonotypes, tmp_path):
        out = plot_category_composition_bar(
            annotated_clonotypes, tmp_path, granularity="category"
        )
        assert out is not None and out.exists()
        assert "composition_category" in out.name

    def test_emits_png_at_each_granularity(self, annotated_clonotypes, tmp_path):
        """All four granularities (category/organism/protein/peptide)
        should produce a separate PNG."""
        for gran in GRANULARITY_COLUMNS:
            out = plot_category_composition_bar(
                annotated_clonotypes, tmp_path, granularity=gran
            )
            assert out is not None, f"granularity {gran} produced no plot"
            assert f"composition_{gran}" in out.name

    def test_top_n_filters_to_subset_per_sample(self, annotated_clonotypes, tmp_path):
        """Passing top_n restricts the bar to the top-N clones per
        sample — should still emit a plot without erroring."""
        out = plot_category_composition_bar(
            annotated_clonotypes, tmp_path, granularity="category", top_n=2
        )
        assert out is not None

    def test_skips_when_granularity_column_missing(
        self, annotated_clonotypes, tmp_path
    ):
        df = annotated_clonotypes.drop(columns=["db_category"])
        out = plot_category_composition_bar(df, tmp_path, granularity="category")
        assert out is None


class TestPlotMatchStrengthComparison:
    """αβ-vs-β-only side-by-side renders."""

    def test_emits_both_strength_panels(self, annotated_clonotypes, tmp_path):
        result = plot_match_strength_comparison(
            annotated_clonotypes, tmp_path, top_n=2, granularity="category"
        )
        assert result is not None
        names = sorted(p.name for p in tmp_path.glob("*.png"))
        # The fixture has both `ab` and `b_only` clones, so we expect
        # one PNG for each strength.
        assert any("_ab.png" in n for n in names)
        assert any("_b_only.png" in n for n in names)

    def test_skips_when_match_strength_column_missing(
        self, annotated_clonotypes, tmp_path
    ):
        df = annotated_clonotypes.drop(columns=["db_match_strength"])
        out = plot_match_strength_comparison(df, tmp_path)
        assert out is None


class TestPlotAnnotationsPhaseBOrchestrator:
    """End-to-end loop: granularity selector × plot family."""

    def test_emits_full_family(self, annotated_clonotypes, tmp_path):
        plot_annotations_phase_b(annotated_clonotypes, tmp_path, top_n=3)
        names = {p.name for p in tmp_path.glob("*.png")}
        # 4 granularities × (1 heatmap + 1 composition + 2 strength
        # comparisons) = 16 expected files.
        for gran in GRANULARITY_COLUMNS:
            assert f"annotation_heatmap_{gran}.png" in names
            assert f"annotation_composition_{gran}.png" in names
            assert f"annotation_heatmap_{gran}_ab.png" in names
            assert f"annotation_heatmap_{gran}_b_only.png" in names

    def test_skips_silently_on_pre_phase_a_data(self, tmp_path):
        """Old clonotype frames without the structured columns must
        not error — the orchestrator should just emit fewer files."""
        legacy = pd.DataFrame(
            {
                "CDR3ab": ["X"],
                "cell_count": [1],
                # No sample_frequencies, no db_* columns.
            }
        )
        plot_annotations_phase_b(legacy, tmp_path)
        assert list(tmp_path.glob("*.png")) == []
