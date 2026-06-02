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

"""Tests for the #72 funnel extensions: αβ-pair-denominator callout,
section dividers, label normalization, and the Selected variant."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from tcrsift.plots import (  # noqa: E402
    FUNNEL_LABEL_NICE,
    create_pipeline_funnel,
    normalize_funnel_label,
    plot_funnel,
)


class TestNormalizeFunnelLabel:
    def test_min_counts_renamed_to_umis(self):
        # The pre-#72 names ``filter:min_counts`` / ``filter:max_counts``
        # read like row-count filters but actually test per-cell UMIs.
        assert "UMI" in normalize_funnel_label("filter:min_counts")
        assert "UMI" in normalize_funnel_label("filter:max_counts")

    def test_unknown_passes_through(self):
        assert normalize_funnel_label("unknown_stage") == "unknown_stage"

    def test_all_filter_aliases_present(self):
        for key in (
            "filter:min_counts", "filter:max_counts",
            "filter:min_genes", "filter:max_genes",
            "filter:min_mito", "filter:max_mito",
            "filter:min_cd3", "filter:min_umi",
        ):
            assert key in FUNNEL_LABEL_NICE


class TestPlotFunnelCallouts:
    def test_renders_with_denominator_stage(self, tmp_path):
        stages = {
            "Loaded": 10000,
            "Confident CD8+": 5000,
            "αβ-pair denominator": 4000,
            "Unique Clonotypes": 2000,
        }
        plot_funnel(stages, tmp_path, denominator_stage="αβ-pair denominator")
        assert (tmp_path / "funnel_plot.png").exists()

    def test_renders_with_section_starts(self, tmp_path):
        stages = {
            "Loaded": 10000,
            "TCR purity": 6000,
            "scRNA QC": 5000,
        }
        plot_funnel(stages, tmp_path, section_starts=("scRNA QC",))
        assert (tmp_path / "funnel_plot.png").exists()

    def test_unknown_denominator_stage_no_error(self, tmp_path):
        """If the caller names a stage that isn't in the dict, the
        callout simply doesn't fire — no crash."""
        stages = {"Loaded": 10, "Clones": 5}
        plot_funnel(stages, tmp_path, denominator_stage="not_a_stage")
        assert (tmp_path / "funnel_plot.png").exists()

    def test_custom_filename_respected(self, tmp_path):
        stages = {"Loaded": 10, "Clones": 5}
        plot_funnel(stages, tmp_path, filename="custom_funnel.png")
        assert (tmp_path / "custom_funnel.png").exists()


class TestCreatePipelineFunnel:
    def test_ab_pair_denominator_inserted(self, tmp_path):
        """When passed, the αβ-pair-denominator stage shows up between
        Phenotyped and Unique Clonotypes."""
        create_pipeline_funnel(
            raw_cells=10000,
            with_vdj=8000,
            phenotyped=7500,
            clonotypes=2000,
            filtered=1500,
            tier_counts={"tier1": 50, "tier2": 100},
            output_dir=tmp_path,
            ab_pair_denominator=6000,
        )
        assert (tmp_path / "funnel_plot.png").exists()

    def test_emit_selected_variant_writes_sibling_funnel(self, tmp_path):
        create_pipeline_funnel(
            raw_cells=10000,
            with_vdj=8000,
            phenotyped=7500,
            clonotypes=2000,
            filtered=1500,
            tier_counts={"tier1": 50, "tier2": 100, "tier3": 300},
            output_dir=tmp_path,
            ab_pair_denominator=6000,
            selected_count=159,
            emit_selected_variant=True,
        )
        assert (tmp_path / "funnel_plot.png").exists()
        assert (tmp_path / "funnel_plot_selected.png").exists()

    def test_selected_variant_off_by_default(self, tmp_path):
        create_pipeline_funnel(
            raw_cells=10, with_vdj=8, phenotyped=7,
            clonotypes=5, filtered=3,
            output_dir=tmp_path,
            selected_count=2,
            # emit_selected_variant left False
        )
        assert (tmp_path / "funnel_plot.png").exists()
        assert not (tmp_path / "funnel_plot_selected.png").exists()

    def test_non_viral_stage_renders_in_selected_funnel(self, tmp_path):
        # exclude_viral path: a "Non-viral" stage precedes "Selected".
        create_pipeline_funnel(
            raw_cells=10000, with_vdj=8000, phenotyped=7500,
            clonotypes=2000, filtered=1500,
            tier_counts={"tier1": 50, "tier2": 100},
            output_dir=tmp_path,
            selected_count=60,
            emit_selected_variant=True,
            non_viral=1490,  # 10 viral bystanders removed
        )
        assert (tmp_path / "funnel_plot_selected.png").exists()
