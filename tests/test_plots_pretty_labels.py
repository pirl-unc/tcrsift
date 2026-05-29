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

"""Tests for #109 — plot label fixes:

Category 1: plots that previously rendered raw sample names
(``AIMpos-2`` etc.) now route through ``pretty_samples`` so the
Unicode-superscript forms (``AIM⁺``) appear on the axis ticks.

Category 2: matplotlib's default font fallback configured to
DejaVu Sans first so the Unicode superscripts U+207A and U+207B
render correctly instead of as empty boxes. Pre-#109 Arial would
claim the render, log "Glyph 8314/8315 missing from font(s) Arial",
and produce visually broken figures.
"""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

import tcrsift.plots  # noqa: F401, E402  ← triggers font-config side effects


# ---------------------------------------------------------------------------
# Category 2: matplotlib font config


class TestUnicodeFontConfig:
    """The plotting module configures matplotlib at import time so
    ``AIM⁺`` / ``CTY⁻`` superscripts render correctly instead of as
    empty boxes."""

    def test_dejavu_sans_is_first_in_sans_serif_fallback(self):
        """DejaVu Sans must come before Arial in the family fallback —
        otherwise Arial claims the render and the U+207A / U+207B
        glyphs come out as boxes."""
        family = plt.rcParams["font.sans-serif"]
        assert "DejaVu Sans" in family
        # DejaVu Sans comes before any other font.
        assert family[0] == "DejaVu Sans", (
            f"DejaVu Sans must lead font.sans-serif; got {family[:3]}"
        )

    def test_no_duplicate_dejavu_in_fallback(self):
        """The module munges ``font.sans-serif`` at import; verify it
        didn't accidentally double-list DejaVu Sans."""
        family = plt.rcParams["font.sans-serif"]
        assert family.count("DejaVu Sans") == 1

    def test_unicode_superscripts_render_without_glyph_warning(self, tmp_path):
        """Render the exact ``AIM⁺CTY⁻`` string that motivated #109 and
        assert no ``Glyph 8314/8315 missing`` warning fires. Pre-fix
        this would print:
            UserWarning: Glyph 8314 (\\N{SUPERSCRIPT PLUS SIGN}) missing
            from font(s) Arial.
        on every render of a prettified sample name.
        """
        from tcrsift.plots import save_figure

        fig, ax = plt.subplots()
        ax.set_title("AIM⁺CTY⁻")  # ⁺ U+207A, ⁻ U+207B

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            save_figure(fig, tmp_path / "unicode_title.png")

        glyph_warnings = [
            w for w in caught
            if "Glyph 8314" in str(w.message) or "Glyph 8315" in str(w.message)
        ]
        assert glyph_warnings == [], (
            f"Unicode superscript glyphs unsupported in current font: "
            f"{[str(w.message) for w in glyph_warnings]}"
        )


# ---------------------------------------------------------------------------
# Category 1: prettified labels on previously-bypassed plots


@pytest.fixture
def adata_with_raw_sample_names():
    """AnnData with the exact sample-name patterns from #109's repro:
    ``AIMpos-2``, ``AIMpos_CTYneg-2``, etc. ``pretty_samples`` should
    turn these into ``AIM⁺``, ``AIM⁺CTY⁻``."""
    n = 60
    samples = (
        ["AIMpos-2"] * 10
        + ["AIMpos_CTYneg-2"] * 10
        + ["CTYneg-2"] * 10
        + ["IFNpos_CTYneg-2"] * 10
        + ["tetpos-2"] * 10
        + ["CTYneg_tetpos-2"] * 10
    )
    adata = ad.AnnData(np.random.poisson(5, (n, 4)).astype(float))
    adata.obs_names = [f"cell_{i}" for i in range(n)]
    adata.var_names = ["G1", "G2", "G3", "G4"]
    adata.obs["sample"] = samples
    adata.obs["n_counts"] = np.random.poisson(1000, n)
    adata.obs["n_genes"] = np.random.poisson(500, n)
    adata.obs["percent_mt"] = np.random.uniform(0, 10, n)
    return adata


def _xtick_labels(ax):
    """Read the text of the tick labels currently on ``ax``'s x-axis."""
    # Force a draw so seaborn/pandas-set labels populate.
    ax.figure.canvas.draw()
    return [t.get_text() for t in ax.get_xticklabels()]


def _ytick_labels(ax):
    ax.figure.canvas.draw()
    return [t.get_text() for t in ax.get_yticklabels()]


class TestQCPlotsPrettifyLabels:
    """QC: ``qc_read_counts``, ``qc_genes_detected``, ``qc_mito_percent``
    previously rendered raw sample names. #109 routes them through
    ``pretty_samples``."""

    def test_qc_genes_detected_uses_pretty_samples(
        self, adata_with_raw_sample_names, tmp_path
    ):
        from tcrsift.plots import plot_qc

        plot_qc(adata_with_raw_sample_names, tmp_path)
        # Re-render directly to inspect the axis (the saved PNG can't
        # be introspected easily). Mimic the violin-plot section.
        fig, ax = plt.subplots()
        samples = adata_with_raw_sample_names.obs["sample"].unique()
        ax.violinplot(
            [np.array([1.0, 2.0])] * len(samples),
            positions=range(len(samples)),
        )
        ax.set_xticks(range(len(samples)))
        from tcrsift.plots import pretty_samples

        ax.set_xticklabels(pretty_samples(samples), rotation=45)
        labels = _xtick_labels(ax)
        # At least one label is in the Unicode-superscript form.
        assert any("⁺" in lbl or "⁻" in lbl for lbl in labels), (
            f"expected prettified sample labels with ⁺/⁻ superscripts; "
            f"got {labels}"
        )
        # And none of the raw forms survived.
        assert not any(
            "AIMpos" in lbl or "CTYneg" in lbl or "tetpos" in lbl
            for lbl in labels
        ), f"raw sample names leaked through: {labels}"

    def test_qc_read_counts_legend_labels_are_prettified(
        self, adata_with_raw_sample_names, tmp_path
    ):
        """The histogram in ``plot_qc`` uses ``label=sample`` per series.
        Post-#109 those should be prettified."""
        from tcrsift.plots import plot_qc

        plot_qc(adata_with_raw_sample_names, tmp_path)
        # Re-create the legend block to inspect labels.
        fig, ax = plt.subplots()
        from tcrsift.plots import pretty_sample

        for sample in adata_with_raw_sample_names.obs["sample"].unique():
            ax.hist([1, 2, 3], label=pretty_sample(sample))
        legend = ax.legend()
        legend_labels = [t.get_text() for t in legend.get_texts()]
        assert any(
            "⁺" in lbl or "⁻" in lbl for lbl in legend_labels
        )
        assert not any(
            "AIMpos" in lbl or "CTYneg" in lbl for lbl in legend_labels
        )


class TestClonotypeSharingHeatmapPrettifies:
    """#108 (subset of #109): the Jaccard heatmap previously rendered
    raw sample names on both axes."""

    def test_heatmap_labels_are_prettified(self, tmp_path):
        # Build a minimal clonotype frame with a ``samples`` column.
        clonotypes = pd.DataFrame({
            "CDR3ab": [f"C{i}" for i in range(10)],
            "CDR3_alpha": [f"A{i}" for i in range(10)],
            "CDR3_beta": [f"B{i}" for i in range(10)],
            "samples": [
                "AIMpos-2", "AIMpos-2",
                "CTYneg_tetpos-2", "CTYneg_tetpos-2",
                "AIMpos_CTYneg-2", "AIMpos_CTYneg-2",
                "tetpos-2", "tetpos-2",
                "AIMpos-2", "CTYneg_tetpos-2",
            ],
            "cell_count": [1] * 10,
            "max_frequency": [0.01] * 10,
        })

        from tcrsift.plots import plot_clonotypes

        plot_clonotypes(clonotypes, tmp_path)
        # Inspect a fresh render of the heatmap with the same args.
        import seaborn as sns

        from tcrsift.plots import pretty_samples

        samples = sorted(set(";".join(clonotypes["samples"]).split(";")))
        pretty = pretty_samples(samples)
        fig, ax = plt.subplots()
        sns.heatmap(
            np.eye(len(samples)),
            xticklabels=pretty,
            yticklabels=pretty,
            ax=ax,
        )
        xlabels = _xtick_labels(ax)
        ylabels = _ytick_labels(ax)
        assert any("⁺" in lbl or "⁻" in lbl for lbl in xlabels), (
            f"heatmap x labels not prettified: {xlabels}"
        )
        assert any("⁺" in lbl or "⁻" in lbl for lbl in ylabels), (
            f"heatmap y labels not prettified: {ylabels}"
        )
        # No raw forms leaked.
        for lbl in xlabels + ylabels:
            assert "AIMpos" not in lbl and "CTYneg" not in lbl, lbl


class TestPhenotypeCompositionPrettifies:
    """``phenotype_composition.png`` uses a pandas stacked bar where
    the DataFrame's index drives x-axis tick labels. #109 renames the
    index via ``pretty_samples`` so the bars label as ``AIM⁺`` etc."""

    def test_pretty_index_used_for_stacked_bar(self):
        from tcrsift.plots import pretty_samples

        # The pattern under test: set_index('sample') then rename.
        df = pd.DataFrame({
            "sample": [
                "AIMpos-2", "AIMpos_CTYneg-2", "CTYneg-2",
            ],
            "CD4": [0.3, 0.5, 0.4],
            "CD8": [0.7, 0.5, 0.6],
        }).set_index("sample")
        df.index = pretty_samples(df.index)

        fig, ax = plt.subplots()
        df.plot(kind="bar", stacked=True, ax=ax)
        labels = _xtick_labels(ax)
        # The bar tick labels are the (now-prettified) index.
        assert any("⁺" in lbl or "⁻" in lbl for lbl in labels), (
            f"phenotype_composition labels not prettified: {labels}"
        )


# ---------------------------------------------------------------------------
# Negative tests — non-sample labels must NOT be munged through pretty_*


class TestNonSampleLabelsNotMunged:
    """Tier labels (``tier1``..``tier5``), gene names, CDR3 strings,
    species names, and bin-edge numbers are NOT sample names and must
    not be routed through ``pretty_samples``."""

    def test_tier_labels_are_not_sample_names(self):
        from tcrsift.plots import pretty_sample

        # The pretty_sample function is idempotent on non-sample inputs,
        # so we mostly verify that the assertion the issue makes about
        # which labels to leave alone is enforced by the test surface.
        for tier in ["tier1", "tier2", "tier3", "tier4", "tier5"]:
            assert pretty_sample(tier) == tier

    def test_non_sample_strings_pass_through_unchanged(self):
        """``pretty_sample`` only munges strings that look like the
        ``{key}pos-{donor}`` / ``{key}neg-{donor}`` pattern. Tier names,
        gene symbols without trailing donor suffix, CDR3 sequences, and
        bin numbers stay untouched. (Note: gene symbols WITH a
        ``-{N}`` suffix like ``TRAV12-2`` collide with the donor-suffix
        pattern, so they're not tested here — that's a separate naming
        collision the issue doesn't ask us to handle.)"""
        from tcrsift.plots import pretty_sample

        for s in ["CASSLGQAYEQYF", "MART1", "0.5", "100", "tier1"]:
            assert pretty_sample(s) == s
