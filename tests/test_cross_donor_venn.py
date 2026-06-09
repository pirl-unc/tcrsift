# Licensed under the Apache License, Version 2.0 (the "License").
"""Pairwise cross-donor clonotype-sharing Venn diagrams (#250)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import builtins  # noqa: E402

import pandas as pd  # noqa: E402

from tcrsift.plots import plot_cross_donor_venn  # noqa: E402


def _donor_frames():
    return {
        "d1": pd.DataFrame({"CDR3ab": ["A_x", "B_y", "C_z"], "CDR3_beta": ["x", "y", "z"]}),
        "d2": pd.DataFrame({"CDR3ab": ["A_x", "D_w"], "CDR3_beta": ["x", "w"]}),
        "d3": pd.DataFrame({"CDR3ab": ["B_y", "D_w", "E_v"], "CDR3_beta": ["y", "w", "v"]}),
    }


def test_writes_png(tmp_path):
    out = tmp_path / "c.png"
    assert plot_cross_donor_venn(_donor_frames(), out) == out
    assert out.exists() and out.stat().st_size > 0


def test_two_donors_single_pair(tmp_path):
    out = tmp_path / "p.png"
    plot_cross_donor_venn({k: _donor_frames()[k] for k in ("d1", "d2")}, out)
    assert out.exists() and out.stat().st_size > 0


def test_noop_fewer_than_two(tmp_path):
    out = tmp_path / "n.png"
    assert plot_cross_donor_venn({"d1": _donor_frames()["d1"]}, out) is None
    assert not out.exists()


def test_fallback_without_matplotlib_venn(tmp_path, monkeypatch):
    real_import = builtins.__import__
    def no_venn(name, *a, **k):
        if name == "matplotlib_venn":
            raise ImportError("forced")
        return real_import(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", no_venn)
    out = tmp_path / "f.png"
    plot_cross_donor_venn(_donor_frames(), out)
    assert out.exists() and out.stat().st_size > 0
