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

"""AnnData atlas plots (#315): UMAP, provenance, signature-vs-background, raincloud.

The first AnnData-input members of plots.py. They consume the atlas
(embed_cells → annotate_cells) and route every antigen legend through
pretty_antigen so a reactive pool is never shown by number alone.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.format import set_antigen_labels  # noqa: E402
from tcrsift.plots import (  # noqa: E402
    plot_provenance,
    plot_raincloud,
    plot_signature_vs_background,
    plot_umap,
    plot_umap_facets,
)


@pytest.fixture
def atlas():
    """A tiny atlas: 120 cells, 3 Leiden clusters / phenotypes, a UMAP embedding,
    per-cell signature columns, an antigen-specific flag, and a peptide source."""
    rng = np.random.default_rng(0)
    n = 120
    leiden = np.array([str(i % 3) for i in range(n)])
    pheno = np.where(
        leiden == "0", "CD8 T",
        np.where(leiden == "1", "CD4 T", "(unspecified)"),
    )
    genes = ["CD8A", "CD4", "CXCL13", "GZMB"]
    X = rng.poisson(3, size=(n, len(genes))).astype(float)
    coords = np.column_stack([
        np.where(leiden == "0", 0.0, np.where(leiden == "1", 5.0, 10.0))
        + rng.normal(0, 0.3, n),
        rng.normal(0, 0.3, n),
    ])
    # A bimodal signature so the raincloud has two modes to show.
    manascore = np.where(rng.random(n) < 0.3, rng.normal(2.0, 0.3, n),
                         rng.normal(-0.5, 0.3, n))
    flag = np.zeros(n, dtype=bool)
    flag[:12] = True  # 12 antigen-specific cells spread across lineages
    peptide = np.array([""] * n, dtype=object)
    peptide[:6] = "P2"
    peptide[6:12] = "P7"
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(leiden),
            "phenotype": pd.Categorical(pheno),
            "manascore": manascore,
            "cytotoxic": rng.normal(0, 1, n),
            "antigen_specific": flag,
            "peptide": peptide,
            # Facet column + two provenance layers (CD8 vs CD4 antigen-specific).
            "timepoint": np.array(["T1", "T2", "T3"])[np.arange(n) % 3],
            "is_as_cd8": flag & (leiden == "0"),
            "is_as_cd4": flag & (leiden == "1"),
        },
        index=[f"c{i}" for i in range(n)],
    )
    a = ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=obs)
    a.obsm["X_umap"] = coords
    return a


def _nonempty(path):
    return path is not None and path.exists() and path.stat().st_size > 0


class TestPlotUmap:
    def test_categorical_phenotype(self, atlas, tmp_path):
        out = plot_umap(atlas, "phenotype", tmp_path / "pheno.png",
                        grayscale_prefixes=("(unspecified)",))
        assert _nonempty(out)

    def test_continuous_signature(self, atlas, tmp_path):
        out = plot_umap(atlas, "manascore", tmp_path / "mana.png")
        assert _nonempty(out)

    def test_gene_expression(self, atlas, tmp_path):
        # A gene name (not an obs column) resolves to log1p expression.
        out = plot_umap(atlas, "CD8A", tmp_path / "cd8a.png")
        assert _nonempty(out)

    def test_leiden_categorical(self, atlas, tmp_path):
        out = plot_umap(atlas, "leiden", tmp_path / "leiden.png")
        assert _nonempty(out)

    def test_missing_basis_raises(self, atlas):
        del atlas.obsm["X_umap"]
        with pytest.raises(ValueError, match="X_umap"):
            plot_umap(atlas, "phenotype")

    def test_missing_color_raises(self, atlas):
        with pytest.raises(ValueError, match="not found"):
            plot_umap(atlas, "does_not_exist")

    def test_returns_ax_when_no_output(self, atlas):
        ax = plot_umap(atlas, "phenotype")
        assert hasattr(ax, "scatter")

    def test_pinned_layout_and_scale(self, atlas):
        # Pinned color_map / centroids / vrange / view stamp all accepted (#326).
        cmap = {"CD8 T": "#d62728", "CD4 T": "#1f77b4", "(unspecified)": (0.8, 0.8, 0.8)}
        centroids = {"CD8 T": (0.0, 0.0), "CD4 T": (5.0, 0.0)}
        ax = plot_umap(atlas, "phenotype", color_map=cmap, centroids=centroids,
                       grayscale_prefixes=("(unspecified)",), view_label="(T1)")
        assert hasattr(ax, "scatter")
        ax2 = plot_umap(atlas, "manascore", vrange=(0.0, 1.0), show_legend=False)
        assert hasattr(ax2, "scatter")


class TestPlotUmapFacets:
    def test_categorical_facets(self, atlas, tmp_path):
        out = plot_umap_facets(atlas, "timepoint", "phenotype",
                               tmp_path / "facets.png",
                               grayscale_prefixes=("(unspecified)",))
        assert _nonempty(out)

    def test_continuous_facets_shared_vrange(self, atlas, tmp_path):
        out = plot_umap_facets(atlas, "timepoint", "manascore",
                               tmp_path / "facets_cont.png", shared_vrange=True)
        assert _nonempty(out)

    def test_returns_figure_with_integrated_plus_facets(self, atlas):
        fig = plot_umap_facets(atlas, "timepoint", "phenotype")
        # integrated + 3 timepoints = 4 drawn panels (rest turned off).
        drawn = [a for a in fig.axes if a.has_data()]
        assert len(drawn) >= 4

    def test_missing_facet_col_raises(self, atlas):
        with pytest.raises(ValueError, match="facet_col"):
            plot_umap_facets(atlas, "nope", "phenotype")


class TestPlotProvenance:
    def test_single_highlight(self, atlas, tmp_path):
        out = plot_provenance(atlas, "antigen_specific", tmp_path / "prov.png")
        assert _nonempty(out)

    def test_multi_layer_overlay(self, atlas, tmp_path):
        out = plot_provenance(
            atlas, layers=[("is_as_cd8", "#d62728", "CD8"),
                           ("is_as_cd4", "#1f77b4", "CD4")],
            output_path=tmp_path / "prov_layers.png",
        )
        assert _nonempty(out)

    def test_multi_layer_legend_labels(self, atlas):
        ax = plot_provenance(atlas, layers=[("is_as_cd8", "#d62728", "CD8"),
                                            ("is_as_cd4", "#1f77b4", "CD4")])
        texts = {t.get_text() for t in ax.get_legend().get_texts()}
        assert {"CD8", "CD4"} <= texts

    def test_flag_and_layers_mutually_exclusive(self, atlas):
        with pytest.raises(ValueError, match="exactly one"):
            plot_provenance(atlas, "antigen_specific",
                            layers=[("is_as_cd8", "#d62728", "CD8")])
        with pytest.raises(ValueError, match="exactly one"):
            plot_provenance(atlas)  # neither

    def test_by_peptide_uses_pretty_antigen(self, atlas, tmp_path):
        # The legend must show the mapped antigen name, not the bare pool token.
        set_antigen_labels({"P2": "KIF1C epitope", "P7": "MAGEA4 epitope"})
        try:
            ax = plot_provenance(atlas, "antigen_specific", by="peptide")
        finally:
            set_antigen_labels({})
        legend_texts = {t.get_text() for t in ax.get_legend().get_texts()}
        assert "KIF1C epitope" in legend_texts
        assert "P2" not in legend_texts  # never a bare pool number

    def test_missing_flag_raises(self, atlas):
        with pytest.raises(ValueError, match="flag"):
            plot_provenance(atlas, "nope")


class TestSignatureVsBackground:
    def test_subset_from_flag(self, atlas, tmp_path):
        out = plot_signature_vs_background(
            atlas, "phenotype", ["manascore", "cytotoxic"],
            tmp_path / "vsbg.png", subset="antigen_specific",
        )
        assert _nonempty(out)

    def test_subset_inferred_from_by(self, atlas, tmp_path):
        # No explicit subset: cells with a non-empty peptide are the subset,
        # and their dots are colored by antigen.
        out = plot_signature_vs_background(
            atlas, "phenotype", ["manascore"], tmp_path / "vsbg_by.png",
            by="peptide",
        )
        assert _nonempty(out)

    def test_requires_a_subset(self, atlas):
        with pytest.raises(ValueError, match="subset"):
            plot_signature_vs_background(atlas, "phenotype", ["manascore"])

    def test_missing_signature_col_raises(self, atlas):
        with pytest.raises(ValueError, match="not in adata.obs"):
            plot_signature_vs_background(
                atlas, "phenotype", ["no_such_score"], subset="antigen_specific"
            )


class TestRaincloud:
    def test_bimodal_groups(self, atlas, tmp_path):
        out = plot_raincloud(
            atlas.obs["phenotype"], atlas.obs["manascore"],
            tmp_path / "rain.png", value_label="MANAscore",
        )
        assert _nonempty(out)

    def test_handles_tiny_and_constant_group(self, tmp_path):
        # A singleton group and a constant group must not crash the KDE path.
        groups = ["a", "a", "a", "b", "c", "c"]
        values = [0.1, 0.2, 0.3, 5.0, 2.0, 2.0]  # b: singleton, c: constant
        out = plot_raincloud(groups, values, tmp_path / "rain_edge.png")
        assert _nonempty(out)

    def test_respects_order(self, tmp_path):
        ax = plot_raincloud(
            ["x", "y", "x"], [1.0, 2.0, 3.0], order=["y", "x"],
        )
        assert [t.get_text() for t in ax.get_yticklabels()] == ["y", "x"]
