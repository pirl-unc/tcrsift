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

"""Faithful per-condition freq∪PRISM selection — matches the downstream recipe (#193).

The reference is the validated selection_analysis/make_plots.py:
    R = lambda s, a: s.rank(pct=True, ascending=a)
    prism = (R(ppost_a, True) + R(ppost_b, True)
             + R(antigen, False) + R(naive, True)) / 4        # lowest = best
    select = top-K by freq  UNION  top-K by prism,  per condition.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from tcrsift.selection import select_freq_prism_per_condition


def _downstream_expected(feat, *, gate, top_freq, top_prism):
    """The exact make_plots.py selection, as a set of (clone, condition)."""
    R = lambda s, a: s.rank(pct=True, ascending=a)  # noqa: E731
    chosen = set()
    for cond, grp in feat.groupby("condition"):
        p = grp[grp["frequency"] > gate].copy()
        if p.empty:
            continue
        p["prism"] = (
            R(p["ppost_alpha"], True) + R(p["ppost_beta"], True)
            + R(p["antigen_response_score"], False) + R(p["naive_score"], True)
        ) / 4
        top_f = list(p.sort_values("frequency", ascending=False)["CDR3ab"].head(top_freq))
        top_p = list(p.sort_values("prism")["CDR3ab"].head(top_prism))
        for c in set(top_f) | set(top_p):
            chosen.add((c, cond))
    return chosen


def _fixture(seed=0, n=40, conds=("AIMpos", "CTYneg", "tetpos")):
    rng = np.random.default_rng(seed)
    rows = []
    for cond in conds:
        for i in range(n):
            rows.append({
                "CDR3ab": f"{cond}_C{i}",
                "condition": cond,
                "frequency": float(rng.random() * 0.05),
                "ppost_alpha": float(10 ** (-rng.random() * 10)),
                "ppost_beta": float(10 ** (-rng.random() * 10)),
                "antigen_response_score": float(rng.normal()),
                "naive_score": float(rng.normal()),
            })
    return pd.DataFrame(rows)


class TestSelectFreqPrismPerCondition:
    def test_matches_downstream_formula(self):
        feat = _fixture()
        gate, top_freq, top_prism = 0.001, 10, 5
        out = select_freq_prism_per_condition(
            feat, condition_col="condition",
            gate=gate, top_freq=top_freq, top_prism=top_prism,
        )
        got = set(zip(out["CDR3ab"], out["condition"]))
        expected = _downstream_expected(
            feat, gate=gate, top_freq=top_freq, top_prism=top_prism
        )
        assert got == expected

    def test_route_labels(self):
        feat = _fixture(seed=3)
        out = select_freq_prism_per_condition(feat, condition_col="condition")
        assert set(out["selection_route"]) <= {"freq", "prism", "both"}
        # A 'both' clone is in both the freq and prism top-K of its condition.
        for _, r in out[out["selection_route"] == "both"].iterrows():
            assert r["rank_within_route"] >= 1

    def test_gate_excludes_low_freq(self):
        feat = pd.DataFrame([
            {"CDR3ab": "A", "condition": "c", "frequency": 0.02,
             "ppost_alpha": 1e-9, "ppost_beta": 1e-9,
             "antigen_response_score": 2.0, "naive_score": 0.0},
            {"CDR3ab": "B", "condition": "c", "frequency": 0.0001,  # below 0.1% gate
             "ppost_alpha": 1e-12, "ppost_beta": 1e-12,
             "antigen_response_score": 3.0, "naive_score": -1.0},
        ])
        out = select_freq_prism_per_condition(feat, condition_col="condition")
        assert set(out["CDR3ab"]) == {"A"}  # B gated out despite better PRISM

    def test_rescue_adds_subthreshold(self):
        feat = pd.DataFrame([
            {"CDR3ab": "A", "condition": "c", "frequency": 0.02, "cov": 5,
             "ppost_alpha": 1e-9, "ppost_beta": 1e-9,
             "antigen_response_score": 2.0, "naive_score": 0.0},
            {"CDR3ab": "B", "condition": "c", "frequency": 0.0001, "cov": 100,
             "ppost_alpha": 1e-12, "ppost_beta": 1e-12,
             "antigen_response_score": 3.0, "naive_score": -1.0},
        ])
        out = select_freq_prism_per_condition(
            feat, condition_col="condition",
            rescue_target=2, rescue_rank_col="cov",
        )
        assert "B" in set(out["CDR3ab"])  # rescued (best-covered sub-threshold)

    def test_empty_input(self):
        out = select_freq_prism_per_condition(pd.DataFrame(), condition_col="condition")
        assert out.empty

    def test_incomplete_term_clone_not_prism_picked(self):
        # require_complete policy: a clone missing a PRISM term (e.g. no GEX) is
        # NOT eligible for the PRISM route, even if strong on the terms it has.
        # Complete-term clones fill the PRISM slots; the deliverable is unchanged
        # when no pick has missing terms.
        from tcrsift.selection import select_freq_prism_per_condition

        feat = pd.DataFrame({
            "CDR3ab": ["C0", "C1", "C2"],
            "method": ["m1"] * 3,
            "frequency": [0.02, 0.015, 0.012],
            "ppost_alpha": [1e-9, 1e-8, 1e-7],
            "ppost_beta": [1e-9, 1e-8, 1e-7],
            # C0 has the best Ppost but is MISSING antigen/naive GEX.
            "antigen_response_score": [float("nan"), 1.0, 0.5],
            "naive_score": [float("nan"), 0.0, 1.0],
        })
        sel = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=0, top_prism=2,
        )
        prism_clones = set(sel[sel["selection_route"] == "prism"]["CDR3ab"])
        assert "C0" not in prism_clones  # incomplete -> not PRISM-picked
        assert prism_clones == {"C1", "C2"}

    def test_dedups_clone_rows_within_condition(self):
        # Per-(clone, sample) input: clone C0 appears in 3 samples of method m1.
        # It must take ONE top-K slot, not three, so top_freq=3 yields 3 DISTINCT
        # clones with correct ranks (F2 — the duplicate-row undercount).
        feat = pd.DataFrame({
            "CDR3ab": ["C0", "C0", "C0", "C1", "C2", "C3"],
            "method": ["m1"] * 6,
            "frequency": [0.05, 0.04, 0.03, 0.02, 0.015, 0.01],
            "ppost_alpha": [1e-9] * 6,
            "ppost_beta": [1e-9] * 6,
            "antigen_response_score": [1.0] * 6,
            "naive_score": [0.0] * 6,
        })
        sel = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=3, top_prism=0,
        )
        assert sel["CDR3ab"].nunique() == 3
        assert set(sel["CDR3ab"]) == {"C0", "C1", "C2"}
        c0 = sel[sel["CDR3ab"] == "C0"]
        assert len(c0) == 1
        assert c0.iloc[0]["rank_within_route"] == 1


def test_cli_select_prism_mode(tmp_path):
    """`tcrsift select --prism` runs the recipe end-to-end."""
    from tcrsift.cli import create_parser

    long = pd.DataFrame({
        "CDR3ab": [f"C{i}" for i in range(6)] * 1,
        "sample": ["S1"] * 6,
        "method": ["AIMpos"] * 6,
        "frequency": [0.02, 0.015, 0.012, 0.01, 0.008, 0.0001],
    })
    scores = pd.DataFrame({
        "CDR3ab": [f"C{i}" for i in range(6)],
        "ppost_alpha": [1e-9, 1e-5, 1e-7, 1e-8, 1e-6, 1e-12],
        "ppost_beta": [1e-9, 1e-5, 1e-7, 1e-8, 1e-6, 1e-12],
        "antigen_response_score": [2.0, 0.1, 1.0, 1.5, 0.5, 3.0],
        "naive_score": [0.0, 2.0, 1.0, 0.5, 1.5, -1.0],
    })
    li = tmp_path / "long.csv"
    long.to_csv(li, index=False)
    sc = tmp_path / "scores.csv"
    scores.to_csv(sc, index=False)
    out = tmp_path / "selected.csv"
    args = create_parser().parse_args([
        "select", "-i", str(li), "--clone-scores", str(sc), "--prism",
        "--condition-col", "method", "--top-freq", "3", "--top-prism", "2",
        "-o", str(out),
    ])
    args.func(args)
    res = pd.read_csv(out)
    assert "selection_route" in res.columns
    assert set(res["selection_route"]) <= {"freq", "prism", "both"}
    # C5 is below the 0.1% gate -> excluded even though it has the best PRISM.
    assert "C5" not in set(res["CDR3ab"])


def test_cli_select_prism_grid(tmp_path):
    """`tcrsift select --prism --prism-grid` emits the trade-off CSV + heatmap (#207)."""
    from tcrsift.cli import create_parser

    long = pd.DataFrame({
        "CDR3ab": [f"C{i}" for i in range(6)],
        "sample": ["S1"] * 6,
        "method": ["AIMpos"] * 6,
        "frequency": [0.02, 0.015, 0.012, 0.01, 0.008, 0.005],
    })
    scores = pd.DataFrame({
        "CDR3ab": [f"C{i}" for i in range(6)],
        "ppost_alpha": [1e-9, 1e-5, 1e-7, 1e-8, 1e-6, 1e-4],
        "ppost_beta": [1e-9, 1e-5, 1e-7, 1e-8, 1e-6, 1e-4],
        "antigen_response_score": [2.0, 0.1, 1.0, 1.5, 0.5, 0.2],
        "naive_score": [0.0, 2.0, 1.0, 0.5, 1.5, 1.8],
    })
    li = tmp_path / "long.csv"
    long.to_csv(li, index=False)
    sc = tmp_path / "scores.csv"
    scores.to_csv(sc, index=False)
    out = tmp_path / "selected.csv"
    args = create_parser().parse_args([
        "select", "-i", str(li), "--clone-scores", str(sc), "--prism",
        "--condition-col", "method", "--prism-grid", "-o", str(out),
    ])
    args.func(args)
    grid_csv = tmp_path / "freq_prism_grid.csv"
    grid_png = tmp_path / "freq_prism_grid.png"
    assert grid_csv.exists()
    grid = pd.read_csv(grid_csv)
    assert set(grid.columns) == {"top_freq", "top_prism", "n_clones"}
    assert grid_png.exists()
