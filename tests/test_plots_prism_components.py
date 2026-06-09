# Licensed under the Apache License, Version 2.0 (the "License").
"""Tests for the PRISM-components small-multiples plot (#249)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import plot_prism_components  # noqa: E402


@pytest.fixture
def scores():
    rng = np.random.default_rng(0)
    n = 30
    df = pd.DataFrame({
        "CDR3ab": [f"C{i:02d}" for i in range(n)],
        "ppost_alpha": rng.random(n), "ppost_beta": rng.random(n),
        "antigen_response_score": rng.random(n), "naive_score": rng.random(n),
    })
    df.loc[0:2, "ppost_alpha"] = np.nan
    return df


def test_emits_nonempty_png(scores, tmp_path):
    out = plot_prism_components(scores, set(scores["CDR3ab"].head(8)), tmp_path / "p.png")
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_renders_with_a_term_missing(scores, tmp_path):
    out = plot_prism_components(scores.drop(columns=["naive_score"]), set(), tmp_path / "m.png")
    assert out is not None and out.exists() and out.stat().st_size > 0


def test_no_op_on_empty(tmp_path):
    empty = pd.DataFrame(columns=["CDR3ab", "ppost_alpha"])
    assert plot_prism_components(empty, set(), tmp_path / "e.png") is None
    assert plot_prism_components(None, set(), tmp_path / "n.png") is None
