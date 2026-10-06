"""String inference must not change SCT flags or per-cell plot annotations (#361)."""

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from tcrsift.plots import (
    _resolve_cell_values,
    _resolve_subset_mask,
    plot_signature_vs_background,
    plot_umap,
    plot_umap_facets,
)
from tcrsift.sct import aggregate_sct


@pytest.fixture(params=[
    pytest.param(("object", "python", False), id="object"),
    pytest.param(("string", "python", False), id="nullable-python"),
    pytest.param(("string", "pyarrow", False), id="nullable-pyarrow"),
    pytest.param((None, "python", False), id="inferred-default"),
    pytest.param((None, "python", True), id="inferred-python"),
    pytest.param((None, "pyarrow", True), id="inferred-pyarrow"),
])
def text_series(request):
    dtype, storage, infer = request.param
    if storage == "pyarrow":
        pytest.importorskip("pyarrow")
    with pd.option_context("mode.string_storage", storage, "future.infer_string", infer):
        def make(values):
            return pd.Series(values, dtype=dtype)
        yield make


def _atlas(values, key):
    obs = pd.DataFrame({key: values})
    obs.index = [f"c{i}" for i in range(len(obs))]
    atlas = ad.AnnData(obs=obs)
    atlas.obsm["X_umap"] = np.column_stack([np.arange(len(obs)), np.zeros(len(obs))])
    return atlas


@pytest.mark.parametrize("column", ["high_quality", "chosen", "complete", "custom_flag"])
def test_sct_yes_no_aggregation(text_series, column):
    flags = text_series(["No", "Yes", "Yes", "No", None])
    frame = pd.DataFrame({
        "CDR3_pair": ["A/B", "C/D", "E/F", "E/F", "E/F"],
        "CDR3_alpha": ["A", "C", "E", "E", "E"],
        "CDR3_beta": ["B", "D", "F", "F", "F"],
        column: flags,
    })
    original = frame.copy(deep=True)
    result = aggregate_sct(frame, numeric_cols=[], boolean_cols=[column], verbose=False)
    result = result.set_index("CDR3_pair")
    assert result[f"{column}.any"].to_dict() == {"A/B": False, "C/D": True, "E/F": True}
    assert result[f"{column}.all"].to_dict() == {"A/B": False, "C/D": True, "E/F": False}
    pd.testing.assert_frame_equal(frame, original)


def test_resolve_text_cell_values(text_series):
    atlas = _atlas(text_series(["CD4", "CD8"]), "cell_type")
    values, kind, label = _resolve_cell_values(atlas, "cell_type")
    assert kind == "categorical"
    assert label == "cell_type"
    assert values.tolist() == ["CD4", "CD8"]
    assert values.index.equals(atlas.obs_names)


@pytest.mark.parametrize("facets", [False, True])
def test_public_umap_preserves_text_labels(text_series, facets):
    atlas = _atlas(text_series(["CD4", "CD8", "CD4", "CD8"]), "cell_type")
    atlas.obs["sample"] = ["T1", "T1", "T2", "T2"]
    if facets:
        fig = plot_umap_facets(atlas, "sample", "cell_type")
        legend = fig.legends[0] if fig.legends else None
        axes = fig.axes
    else:
        ax = plot_umap(atlas, "cell_type")
        fig, legend, axes = ax.figure, ax.get_legend(), [ax]
    try:
        assert legend is not None
        assert {t.get_text() for t in legend.get_texts()} == {"CD4", "CD8"}
        assert sum(len(c.get_offsets()) for c in axes[0].collections) == 4
        assert {t.get_text() for t in axes[0].texts} >= {"CD4", "CD8"}
    finally:
        plt.close(fig)


def test_missing_peptides_are_not_assigned(text_series):
    atlas = _atlas(text_series(["AAAAAAAAA", None, np.nan, pd.NA, "", "nan", "None", "NaN", "<NA>"]), "peptide")
    mask = _resolve_subset_mask(atlas, None, "peptide")
    assert mask.dtype == np.dtype(bool)
    assert mask.tolist() == [True] + [False] * 8

    atlas.obs["phenotype"] = "CD8"
    atlas.obs["score"] = np.arange(atlas.n_obs, dtype=float)
    fig = plot_signature_vs_background(atlas, "phenotype", ["score"], by="peptide")
    try:
        offsets = fig.axes[0].collections[0].get_offsets()
        assert len(offsets) == 1
        assert offsets[0, 1] == 0
    finally:
        plt.close(fig)
