# Licensed under the Apache License, Version 2.0 (the "License").
"""Tests for the frequency × PRISM selection scatter (#248)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import plot_freq_prism_scatter  # noqa: E402


@pytest.fixture
def selection_candidates():
    rows = []
    for cond in ("flu", "cmv"):
        rows += [
            {"CDR3ab": f"{cond}_a", "method": cond, "frequency": 0.20, "prism_score": 0.5, "selection_route": "both"},
            {"CDR3ab": f"{cond}_b", "method": cond, "frequency": 0.10, "prism_score": 1.2, "selection_route": "freq"},
            {"CDR3ab": f"{cond}_c", "method": cond, "frequency": 0.01, "prism_score": 0.8, "selection_route": "prism"},
            {"CDR3ab": f"{cond}_d", "method": cond, "frequency": 0.001, "prism_score": 3.0, "selection_route": "unselected"},
            {"CDR3ab": f"{cond}_e", "method": cond, "frequency": 0.002, "prism_score": float("nan"), "selection_route": "unselected"},
        ]
    return pd.DataFrame(rows)


def test_emits_nonempty_png(selection_candidates, tmp_path):
    out = plot_freq_prism_scatter(selection_candidates, tmp_path / "s.png", gate=0.05)
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_single_condition(selection_candidates, tmp_path):
    one = selection_candidates[selection_candidates["method"] == "flu"]
    out = plot_freq_prism_scatter(one, tmp_path / "one.png")
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_empty_input_is_noop(tmp_path):
    empty = pd.DataFrame(columns=["CDR3ab", "method", "frequency", "prism_score", "selection_route"])
    assert plot_freq_prism_scatter(empty, tmp_path / "x.png") is None
    assert list(tmp_path.glob("*.png")) == []
