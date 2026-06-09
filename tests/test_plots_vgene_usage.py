# Licensed under the Apache License, Version 2.0 (the "License").
"""Tests for the per-method V-gene usage heatmap (#251)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import plot_vgene_usage_by_method  # noqa: E402


@pytest.fixture
def usage_data():
    long_df = pd.DataFrame({
        "CDR3ab": ["c0", "c1", "c2", "c3"], "method": ["A", "A", "B", "B"],
        "frequency": [0.4, 0.3, 0.2, 0.1],
    })
    clonotypes = pd.DataFrame({
        "CDR3ab": ["c0", "c1", "c2", "c3"],
        "alpha_v_gene": ["TRAV1", "TRAV1", "TRAV2", "TRAV2"],
        "beta_v_gene": ["TRBV1", "TRBV2", "TRBV1", "TRBV2"],
    })
    return long_df, clonotypes


def test_emits_nonempty_png(usage_data, tmp_path):
    out = plot_vgene_usage_by_method(*usage_data, tmp_path / "v.png")
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_none_on_empty_long(usage_data, tmp_path):
    _, clo = usage_data
    assert plot_vgene_usage_by_method(pd.DataFrame(columns=["CDR3ab", "method", "frequency"]), clo, tmp_path / "v.png") is None


def test_none_when_no_vgene_cols(usage_data, tmp_path):
    long_df, clo = usage_data
    assert plot_vgene_usage_by_method(long_df, clo.drop(columns=["alpha_v_gene", "beta_v_gene"]), tmp_path / "v.png") is None


def test_one_vgene_col_only(usage_data, tmp_path):
    long_df, clo = usage_data
    out = plot_vgene_usage_by_method(long_df, clo.drop(columns=["beta_v_gene"]), tmp_path / "v.png")
    assert out is not None and out.exists()


def test_row_normalization(usage_data, monkeypatch):
    import tcrsift.plots as plots
    captured = []
    real = plots.sns.heatmap
    def spy(data, *a, **k):
        captured.append(data.copy())
        return real(data, *a, **k)
    monkeypatch.setattr(plots.sns, "heatmap", spy)
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        plots.plot_vgene_usage_by_method(*usage_data, Path(d) / "v.png")
    assert captured
    for frac in captured:
        assert (frac.sum(axis=1).round(6) == 1.0).all()
    assert captured[0].loc["A", "TRAV1"] == pytest.approx(1.0)
