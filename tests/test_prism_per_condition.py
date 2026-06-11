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

    def test_rows_emitted_in_selection_rank_order(self):
        # #230: rows must be serialized in a deterministic "order of selected
        # clones" — within each condition the freq route by ascending rank
        # (incl. 'both'), then prism-only clones by ascending rank. The old
        # ``set(top_f) | set(top_p)`` iteration had no such contract and churned
        # selection_native.csv row order across processes (hash-seed dependent).
        feat = _fixture(seed=7)
        out = select_freq_prism_per_condition(feat, condition_col="condition")
        for _, grp in out.groupby("condition", sort=False):
            is_prism_only = (grp["selection_route"] == "prism").to_numpy()
            # All freq-route rows come before any prism-only row.
            split = is_prism_only.argmax() if is_prism_only.any() else len(grp)
            assert not is_prism_only[:split].any()
            # Each block is strictly ascending by rank_within_route.
            freq_ranks = list(grp["rank_within_route"].iloc[:split])
            prism_ranks = list(grp["rank_within_route"].iloc[split:])
            assert freq_ranks == sorted(freq_ranks)
            assert prism_ranks == sorted(prism_ranks)

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

    def test_freq_tiebreak_is_deterministic(self):
        # #221: equally-frequent clones at the top-K boundary must select the
        # same set regardless of input row order (was an arbitrary artifact).
        from tcrsift.selection import select_freq_prism_per_condition

        base = pd.DataFrame({
            "CDR3ab": ["Cb", "Ca", "Cd", "Cc"],
            "method": ["m"] * 4,
            "frequency": [0.02] * 4,  # all tied
        })
        picks = set()
        for seed in range(4):
            shuf = base.sample(frac=1, random_state=seed).reset_index(drop=True)
            sel = select_freq_prism_per_condition(
                shuf, condition_col="method", top_freq=2, top_prism=0,
            )
            picks.add(tuple(sorted(sel["CDR3ab"])))
        assert len(picks) == 1  # reproducible
        # CDR3ab fallback alone (no metadata) → lexical: Ca, Cb.
        assert next(iter(picks)) == ("Ca", "Cb")

    def test_freq_tiebreak_heuristic_single_alpha_then_quality(self):
        # #223: among equally-frequent clones, single-α beats dual-α, then UMIs,
        # then purity — a quality ranking, not arbitrary.
        from tcrsift.selection import select_freq_prism_per_condition

        feat = pd.DataFrame({
            "CDR3ab": ["Cb", "Ca", "Cd", "Cc"],
            "method": ["m"] * 4,
            "frequency": [0.02] * 4,
            "merged_alpha_partners": [None, "X;Y", None, None],  # Ca is dual-α
            "TRA_1_umis": [10, 99, 50, 10],
            "TRB_1_umis": [10, 99, 5, 10],
            "Tcell_type_purity": [1.0, 1.0, 1.0, 0.5],
        })
        sel = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=2, top_prism=0,
        )
        picked = set(sel["CDR3ab"])
        assert "Ca" not in picked  # dual-α de-prioritized despite equal freq
        assert picked == {"Cd", "Cb"}  # single-α, higher UMI / purity

    def test_tiebreak_configurable_cdr3ab_only(self):
        from tcrsift.selection import select_freq_prism_per_condition

        feat = pd.DataFrame({
            "CDR3ab": ["Cb", "Ca", "Cd", "Cc"],
            "method": ["m"] * 4,
            "frequency": [0.02] * 4,
            "merged_alpha_partners": [None, "X;Y", None, None],
        })
        # tie_break=[] disables the heuristic → pure lexical CDR3ab.
        sel = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=2, top_prism=0, tie_break=[],
        )
        assert set(sel["CDR3ab"]) == {"Ca", "Cb"}

    def test_all_nan_ppost_raises(self):
        # PRISM requested but a score column is entirely NaN -> loud error, not
        # a silent empty PRISM route (the "all ppost NaN = error" guard).
        import pytest

        from tcrsift.selection import select_freq_prism_per_condition
        from tcrsift.validation import TCRsiftValidationError

        feat = pd.DataFrame({
            "CDR3ab": ["C0", "C1"],
            "method": ["m1", "m1"],
            "frequency": [0.02, 0.015],
            "ppost_alpha": [float("nan"), float("nan")],  # never populated
            "ppost_beta": [1e-9, 1e-8],
            "antigen_response_score": [1.0, 0.5],
            "naive_score": [0.0, 1.0],
        })
        with pytest.raises(TCRsiftValidationError, match="missing or entirely NaN"):
            select_freq_prism_per_condition(
                feat, condition_col="method", top_freq=0, top_prism=2,
            )

    def test_all_nan_ppost_ok_when_prism_off(self):
        # top_prism=0 -> frequency route only -> the guard must NOT fire.
        from tcrsift.selection import select_freq_prism_per_condition

        feat = pd.DataFrame({
            "CDR3ab": ["C0", "C1"],
            "method": ["m1", "m1"],
            "frequency": [0.02, 0.015],
            "ppost_alpha": [float("nan"), float("nan")],
            "ppost_beta": [float("nan"), float("nan")],
            "antigen_response_score": [float("nan"), float("nan")],
            "naive_score": [float("nan"), float("nan")],
        })
        sel = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=2, top_prism=0,
        )
        assert set(sel["CDR3ab"]) == {"C0", "C1"}

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


class TestPrismCandidates:
    def test_emits_all_gated_with_routes(self):
        from tcrsift.selection import prism_candidates, select_freq_prism_per_condition

        feat = _fixture(seed=5)
        gate, top_freq, top_prism = 0.001, 10, 5
        cand = prism_candidates(
            feat, condition_col="condition", gate=gate,
            top_freq=top_freq, top_prism=top_prism,
        )
        # Columns + route vocabulary.
        assert set(cand.columns) == {
            "CDR3ab", "condition", "frequency", "prism_score", "selection_route"
        }
        assert set(cand["selection_route"]) <= {"freq", "prism", "both", "unselected"}
        # Every emitted candidate is above the gate.
        assert (cand["frequency"] > gate).all()
        # The non-'unselected' rows exactly match the actual selection.
        picked = set(zip(
            cand.loc[cand.selection_route != "unselected", "CDR3ab"],
            cand.loc[cand.selection_route != "unselected", "condition"],
        ))
        sel = select_freq_prism_per_condition(
            feat, condition_col="condition", gate=gate,
            top_freq=top_freq, top_prism=top_prism,
        )
        assert picked == set(zip(sel["CDR3ab"], sel["condition"]))

    def test_empty_input(self):
        from tcrsift.selection import prism_candidates

        out = prism_candidates(pd.DataFrame(), condition_col="condition")
        assert out.empty
        assert "selection_route" in out.columns


class TestPerRouteRanksAndSummary:
    """Independent freq/PRISM ranks + per-clone folding (#selection-cols)."""

    def test_emits_independent_freq_and_prism_ranks(self):
        # A clone in BOTH top-freq and top-prism gets both ranks; freq-only /
        # prism-only clones get just the one they placed in.
        feat = pd.DataFrame([
            {"CDR3ab": "A", "condition": "AIMpos", "frequency": 0.9, "score": 0.1},
            {"CDR3ab": "B", "condition": "AIMpos", "frequency": 0.5, "score": 0.9},
            {"CDR3ab": "C", "condition": "AIMpos", "frequency": 0.1, "score": 0.2},
        ])
        out = select_freq_prism_per_condition(
            feat, condition_col="condition", freq_col="frequency",
            score_terms=[{"col": "score", "direction": "low", "weight": 1.0}],
            gate=0.0, top_freq=2, top_prism=2, tie_break=[],
        )
        by = out.set_index("CDR3ab")
        # A: freq#1 (0.9) AND prism#1 (lowest score 0.1) → both ranks present
        assert by.loc["A", "selection_route"] == "both"
        assert by.loc["A", "freq_rank"] == 1 and by.loc["A", "prism_rank"] == 1
        # B: freq#2 only (prism rank 3 > top_prism=2)
        assert by.loc["B", "selection_route"] == "freq"
        assert by.loc["B", "freq_rank"] == 2 and pd.isna(by.loc["B", "prism_rank"])
        # C: prism#2 only (freq rank 3 > top_freq=2)
        assert by.loc["C", "selection_route"] == "prism"
        assert pd.isna(by.loc["C", "freq_rank"]) and by.loc["C", "prism_rank"] == 2

    def test_summarize_selection_per_clone_columns(self):
        from tcrsift.selection import (
            SELECTION_SUMMARY_COLS,
            summarize_selection_per_clone,
        )
        # One clone selected in two conditions: AIMpos (freq#1), CTYneg (prism#2).
        per_cond = pd.DataFrame([
            {"CDR3ab": "X", "method": "AIMpos", "selection_route": "freq",
             "rank_within_route": 1, "freq_rank": 1, "prism_rank": None,
             "prism_score": 0.42, "frequency": 0.009, "rescued": False},
            {"CDR3ab": "X", "method": "CTYneg", "selection_route": "prism",
             "rank_within_route": 2, "freq_rank": None, "prism_rank": 2,
             "prism_score": 0.18, "frequency": 0.0034, "rescued": False},
        ])
        out = summarize_selection_per_clone(per_cond)
        assert list(out.columns) == ["CDR3ab", *SELECTION_SUMMARY_COLS]
        r = out.iloc[0]
        assert r["selection_conditions"] == "AIM⁺;CTY⁻"
        assert r["selection_detail"] == "AIM⁺=freq#1(0.90%);CTY⁻=prism#2(0.34%)"
        assert r["frequency_per_condition"] == "AIM⁺=0.90%;CTY⁻=0.34%"
        assert r["prism_per_condition"] == "AIM⁺=0.420;CTY⁻=0.180"
        assert r["freq_rank_per_condition"] == "AIM⁺=1"   # only AIMpos freq-placed
        assert r["prism_rank_per_condition"] == "CTY⁻=2"  # only CTYneg prism-placed
        # tcrsift convention: plain ';' separator, no surrounding spaces
        assert " ; " not in r["selection_detail"]

    def test_summarize_empty(self):
        from tcrsift.selection import (
            SELECTION_SUMMARY_COLS,
            summarize_selection_per_clone,
        )
        out = summarize_selection_per_clone(pd.DataFrame())
        assert out.empty
        assert list(out.columns) == ["CDR3ab", *SELECTION_SUMMARY_COLS]
