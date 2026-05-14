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

"""Tests for the remaining #43 plots:

- Slope chart (#1, ``plot_clone_tracking_slopes``)
- N-methods distribution (#6, ``plot_clones_by_n_methods``)
- Signature scatter (#7, ``plot_clone_freq_vs_signature_per_sample``)
- Funnel variants (#8, ``plot_funnel_{ribbon,lollipop,terrace}``)
- Method-grouping helper (``make_method_panel_grid``)

Asserts file emission + skip-on-missing-input behavior rather than
visual content.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

from tcrsift.plots import (  # noqa: E402
    ACTIVATION_GENES_HGNC,
    ANTIGEN_RESPONSE_GENES_HGNC,
    CYTOLYTIC_GENES_HGNC,
    EXHAUSTION_GENES_HGNC,
    TUMOR_REACTIVE_GENES_HGNC,
    _derive_n_methods,
    make_method_panel_grid,
    plot_clone_freq_vs_signature_per_sample,
    plot_clone_tracking_slopes,
    plot_clones_by_n_methods,
    plot_funnel_lollipop,
    plot_funnel_ribbon,
    plot_funnel_terrace,
    plot_funnel_variants,
)


@pytest.fixture
def ranked_clonotypes():
    """10 clones × 4 samples with monotone-decreasing frequencies."""
    return pd.DataFrame(
        {
            "CDR3ab": [f"CASS{i:02d}F_CAVT{i:02d}F" for i in range(10)],
            "cell_count": list(range(60, 0, -6)),
            "tier": ["tier1"] * 2 + ["tier2"] * 3 + [None] * 5,
            "sample_frequencies": [
                {
                    "S1": max(0.0, 0.3 - i * 0.03),
                    "S2": max(0.0, 0.2 - i * 0.02),
                    "S3": max(0.0, 0.15 - i * 0.015),
                    "S4": max(0.0, 0.1 - i * 0.01),
                }
                for i in range(10)
            ],
        }
    )


@pytest.fixture
def adata_with_signature_genes(ranked_clonotypes):
    """Minimal AnnData where var_names include every default signature
    gene (activation, exhaustion, antigen-response, cytolytic,
    tumor-reactive) so any scatter call has data to work with."""
    rng = np.random.default_rng(0)
    gene_ids = sorted(
        set(ACTIVATION_GENES_HGNC)
        | set(EXHAUSTION_GENES_HGNC)
        | set(ANTIGEN_RESPONSE_GENES_HGNC)
        | set(CYTOLYTIC_GENES_HGNC)
        | set(TUMOR_REACTIVE_GENES_HGNC)
    ) + ["FOO", "BAR"]
    n_cells = 60
    X = csr_matrix(rng.poisson(0.5, size=(n_cells, len(gene_ids))).astype(float))
    a = ad.AnnData(X=X)
    a.var_names = gene_ids
    a.obs["sample"] = ["S1"] * 20 + ["S2"] * 15 + ["S3"] * 15 + ["S4"] * 10
    a.obs["CDR3ab"] = rng.choice(ranked_clonotypes["CDR3ab"], size=n_cells)
    return a


class TestFocalSignatureConstants:
    """The 2-gene focal signatures in plots.py must stay in sync with
    the til_select defaults so the per-sample scatter and TIL-selection
    scores agree on gene-set membership (#70)."""

    def test_antigen_response_matches_til_select(self):
        from tcrsift.til_select import ANTIGEN_RESPONSE_GENES_DEFAULT

        assert tuple(ANTIGEN_RESPONSE_GENES_HGNC) == tuple(
            ANTIGEN_RESPONSE_GENES_DEFAULT
        )

    def test_cytolytic_matches_til_select(self):
        from tcrsift.til_select import CYTOLYTIC_GENES_DEFAULT

        assert tuple(CYTOLYTIC_GENES_HGNC) == tuple(CYTOLYTIC_GENES_DEFAULT)

    def test_tumor_reactive_matches_til_select(self):
        from tcrsift.til_select import ENRICHMENT_GENES_DEFAULT

        assert tuple(TUMOR_REACTIVE_GENES_HGNC) == tuple(ENRICHMENT_GENES_DEFAULT)


class TestMakeMethodPanelGrid:
    """The generic multi-panel grid helper."""

    def test_default_row_major_chunks(self):
        fig, axes = make_method_panel_grid(["a", "b", "c", "d", "e"], cols=3)
        assert set(axes) == {"a", "b", "c", "d", "e"}
        # 5 panels, 3 cols → 2 rows.
        assert fig.axes  # at least one axes registered
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_custom_layout_preserved(self):
        """A caller-supplied layout is used verbatim — the helper
        doesn't reorder or merge rows."""
        custom = [["a", "b", "c"], ["d", "e"]]
        fig, axes = make_method_panel_grid(["a", "b", "c", "d", "e"], layout=custom)
        assert set(axes) == {"a", "b", "c", "d", "e"}
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestPlotCloneTrackingSlopes:
    """Slope chart — #43 item #1."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_clone_tracking_slopes(ranked_clonotypes, tmp_path, top_n=5)
        assert out is not None and out.exists() and out.suffix == ".png"

    def test_skips_when_no_sample_frequencies(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_clone_tracking_slopes(df, tmp_path) is None
        assert list(tmp_path.glob("*.png")) == []


class TestDeriveNMethods:
    """Helper that backs the n-methods plot."""

    def test_prefers_explicit_n_methods_column(self, ranked_clonotypes):
        df = ranked_clonotypes.copy()
        df["n_methods"] = [3] * len(df)
        assert (_derive_n_methods(df) == 3).all()

    def test_derives_from_sample_frequencies(self, ranked_clonotypes):
        n = _derive_n_methods(ranked_clonotypes)
        assert n is not None
        # Clone 0 has all 4 samples > 0.
        assert int(n.iloc[0]) == 4

    def test_returns_none_when_no_source_available(self):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert _derive_n_methods(df) is None


class TestPlotClonesByNMethods:
    """N-methods 1-D distribution — #43 item #6."""

    def test_emits_png(self, ranked_clonotypes, tmp_path):
        out = plot_clones_by_n_methods(ranked_clonotypes, tmp_path)
        assert out is not None and out.exists()

    def test_renders_without_tier_column(self, ranked_clonotypes, tmp_path):
        """Tier highlight is optional — plot must still render."""
        df = ranked_clonotypes.drop(columns=["tier"])
        out = plot_clones_by_n_methods(df, tmp_path)
        assert out is not None

    def test_skips_when_no_source(self, tmp_path):
        df = pd.DataFrame({"CDR3ab": ["X"], "cell_count": [1]})
        assert plot_clones_by_n_methods(df, tmp_path) is None


class TestPlotCloneFreqVsSignaturePerSample:
    """Signature scatter — #43 item #7 (GEX-dependent)."""

    def test_emits_png_with_default_activation_genes(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=ACTIVATION_GENES_HGNC,
            signature_label="activation",
            output_dir=tmp_path,
        )
        assert out is not None and out.exists()

    def test_works_with_exhaustion_genes(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=EXHAUSTION_GENES_HGNC,
            signature_label="exhaustion",
            output_dir=tmp_path,
        )
        assert out is not None

    def test_works_with_antigen_response_genes(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        """The third focal signature added in #70."""
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=ANTIGEN_RESPONSE_GENES_HGNC,
            signature_label="antigen-response",
            output_dir=tmp_path,
        )
        assert out is not None and out.exists()

    def test_works_with_cytolytic_genes(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        """Parity with the experiment's three focal signatures (#70)."""
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=CYTOLYTIC_GENES_HGNC,
            signature_label="cytolytic",
            output_dir=tmp_path,
        )
        assert out is not None

    def test_works_with_tumor_reactive_genes(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        """Parity with the experiment's three focal signatures (#70)."""
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=TUMOR_REACTIVE_GENES_HGNC,
            signature_label="tumor-reactive",
            output_dir=tmp_path,
        )
        assert out is not None

    def test_skips_when_genes_not_in_adata(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        """Zero-overlap gene list → skip with no error."""
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            ranked_clonotypes,
            gene_ids=["NONEXISTENT_GENE_1", "NONEXISTENT_GENE_2"],
            signature_label="bogus",
            output_dir=tmp_path,
        )
        assert out is None
        assert list(tmp_path.glob("*.png")) == []

    def test_skips_when_adata_obs_missing_columns(
        self, ranked_clonotypes, tmp_path
    ):
        rng = np.random.default_rng(0)
        X = csr_matrix(rng.poisson(0.5, size=(20, 5)).astype(float))
        a = ad.AnnData(X=X)
        a.var_names = list(ACTIVATION_GENES_HGNC)
        # Note: no `sample` or `CDR3ab` in obs.
        out = plot_clone_freq_vs_signature_per_sample(
            a,
            ranked_clonotypes,
            gene_ids=ACTIVATION_GENES_HGNC,
            signature_label="activation",
            output_dir=tmp_path,
        )
        assert out is None

    def test_skips_when_clonotypes_missing_sample_frequencies(
        self, adata_with_signature_genes, ranked_clonotypes, tmp_path
    ):
        df = ranked_clonotypes.drop(columns=["sample_frequencies"])
        out = plot_clone_freq_vs_signature_per_sample(
            adata_with_signature_genes,
            df,
            gene_ids=ACTIVATION_GENES_HGNC,
            signature_label="activation",
            output_dir=tmp_path,
        )
        assert out is None


class TestFunnelVariants:
    """Funnel ribbon / lollipop / terrace — #43 item #8."""

    @pytest.fixture
    def stage_counts(self):
        return {
            "Loaded": 10000,
            "Complete αβ": 7500,
            "Confident CD8+": 4500,
            "≥2 cells": 1200,
            "Tier1": 80,
        }

    def test_ribbon_emits(self, stage_counts, tmp_path):
        out = plot_funnel_ribbon(stage_counts, tmp_path)
        assert out.exists() and out.name == "funnel_ribbon.png"

    def test_lollipop_emits(self, stage_counts, tmp_path):
        out = plot_funnel_lollipop(stage_counts, tmp_path)
        assert out.exists() and out.name == "funnel_lollipop.png"

    def test_terrace_emits(self, stage_counts, tmp_path):
        out = plot_funnel_terrace(stage_counts, tmp_path)
        assert out.exists() and out.name == "funnel_terrace.png"

    def test_variants_orchestrator_emits_all_four(self, stage_counts, tmp_path):
        paths = plot_funnel_variants(stage_counts, tmp_path)
        names = {p.name for p in paths}
        assert names == {
            "funnel_plot.png",
            "funnel_ribbon.png",
            "funnel_lollipop.png",
            "funnel_terrace.png",
        }
        assert all(p.exists() for p in paths)

    def test_variants_handle_single_stage(self, tmp_path):
        """One-stage funnels are degenerate but must not error."""
        out = plot_funnel_ribbon({"only": 100}, tmp_path)
        assert out.exists()
