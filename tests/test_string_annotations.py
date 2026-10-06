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


@pytest.mark.parametrize("labels", [["CD4", "CD8"], ["1", "2"]])
def test_resolve_text_cell_values(text_series, labels):
    atlas = _atlas(text_series(labels), "cell_type")
    values, kind, label = _resolve_cell_values(atlas, "cell_type")
    assert kind == "categorical"
    assert label == "cell_type"
    assert values.tolist() == labels
    assert values.index.equals(atlas.obs_names)


@pytest.mark.parametrize("facets", [False, True])
@pytest.mark.parametrize("missing", [False, True])
def test_public_umap_preserves_text_labels(text_series, facets, missing):
    labels = ["CD4", "CD8", None if missing else "CD4", "CD8"]
    expected = {"CD4", "CD8", "nan"} if missing else {"CD4", "CD8"}
    atlas = _atlas(text_series(labels), "cell_type")
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
        assert {t.get_text() for t in legend.get_texts()} == expected
        assert sum(len(c.get_offsets()) for c in axes[0].collections) == 4
        assert {t.get_text() for t in axes[0].texts} >= expected
    finally:
        plt.close(fig)


def test_missing_facet_labels_keep_cells(text_series):
    atlas = _atlas(text_series(["T1", None, "T2"]), "sample")
    atlas.obs["cell_type"] = "CD8"
    fig = plot_umap_facets(atlas, "sample", "cell_type", include_integrated=False)
    try:
        assert len(fig.axes) == 3
        assert [sum(len(c.get_offsets()) for c in ax.collections) for ax in fig.axes] == [1, 1, 1]
        assert {t.get_text() for ax in fig.axes for t in ax.texts} >= {"(T1)", "(T2)", "(nan)"}
    finally:
        plt.close(fig)


@pytest.mark.parametrize("dtype", ["int64", "float64"])
def test_numeric_cell_values_stay_continuous(dtype):
    atlas = _atlas(pd.Series([1, 2], dtype=dtype), "score")
    values, kind, _ = _resolve_cell_values(atlas, "score")
    assert kind == "continuous"
    assert values.tolist() == [1, 2]


@pytest.mark.parametrize("dtype", ["bool", "boolean", "category"])
def test_boolean_and_categorical_cell_values_stay_categorical(dtype):
    atlas = _atlas(pd.Series([True, False], dtype=dtype), "flag")
    values, kind, _ = _resolve_cell_values(atlas, "flag")
    assert kind == "categorical"
    assert values.tolist() == ["True", "False"]


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
