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

"""Tests for partial-information clonotype attribution (#176)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift.attribution import attribute_cells
from tcrsift.config import AttributionConfig


def _cells(rows):
    """Build a cell DataFrame; fills the chain columns attribute_cells reads.

    Each row dict may set a1/a2/b1/b2 (CDR3 strings, "" = absent) and sample;
    UMIs default to 5 (passing) for any present chain.
    """
    recs = {}
    cols = ["CDR3_alpha", "TRA_2_cdr3", "CDR3_beta", "TRB_2_cdr3",
            "TRA_1_umis", "TRA_2_umis", "TRB_1_umis", "TRB_2_umis", "sample"]
    data = {c: [] for c in cols}
    idx = []
    for i, r in enumerate(rows):
        idx.append(r.get("bc", f"c{i}"))
        a1, a2 = r.get("a1", ""), r.get("a2", "")
        b1, b2 = r.get("b1", ""), r.get("b2", "")
        data["CDR3_alpha"].append(a1)
        data["TRA_2_cdr3"].append(a2)
        data["CDR3_beta"].append(b1)
        data["TRB_2_cdr3"].append(b2)
        data["TRA_1_umis"].append(r.get("au1", 5 if a1 else 0))
        data["TRA_2_umis"].append(r.get("au2", 5 if a2 else 0))
        data["TRB_1_umis"].append(r.get("bu1", 5 if b1 else 0))
        data["TRB_2_umis"].append(r.get("bu2", 5 if b2 else 0))
        data["sample"].append(r.get("sample", "S1"))
    return pd.DataFrame(data, index=idx)


def _on(**kw):
    return AttributionConfig(enabled=True, **kw)


class TestCompleteCells:
    def test_complete_single_weight_one(self):
        df = _cells([{"a1": "A", "b1": "B"}, {"a1": "A", "b1": "B"}])
        long, merge = attribute_cells(df, _on())
        assert merge == {}
        assert set(long["CDR3ab"]) == {"A_B"}
        assert long["weight"].sum() == pytest.approx(2.0)
        # Each cell contributes exactly 1.0.
        per_cell = long.groupby("cell_barcode")["weight"].sum()
        assert per_cell.to_dict() == pytest.approx({"c0": 1.0, "c1": 1.0})


class TestBetaShare:
    def test_alpha_dropout_splits_proportional_3_to_1(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior = 3
            {"a1": "A2", "b1": "B"},   # A2_B prior = 1
            {"bc": "drop", "a1": "", "b1": "B"},  # alpha-dropout, beta B
        ])
        long, _ = attribute_cells(df, _on())
        drop = long[long["cell_barcode"] == "drop"].set_index("CDR3ab")["weight"]
        assert drop["A1_B"] == pytest.approx(0.75)
        assert drop["A2_B"] == pytest.approx(0.25)
        assert drop.sum() == pytest.approx(1.0)

    def test_alpha_dropout_no_match_dropped(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "OTHER"},  # beta with no complete clone
        ])
        long, _ = attribute_cells(df, _on())
        assert "drop" not in set(long["cell_barcode"])

    def test_beta_sharing_disabled_drops(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "B"},
        ])
        long, _ = attribute_cells(df, _on(beta_sharing=False))
        assert "drop" not in set(long["cell_barcode"])

    def test_uniform_weighting_ignores_priors(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "drop", "a1": "", "b1": "B"},
        ])
        long, _ = attribute_cells(df, _on(share_weighting="uniform"))
        drop = long[long["cell_barcode"] == "drop"].set_index("CDR3ab")["weight"]
        assert drop["A1_B"] == pytest.approx(0.5)
        assert drop["A2_B"] == pytest.approx(0.5)


class TestDualAlpha:
    def test_recurrent_dual_alpha_merges(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior 3 -> canonical
            {"a1": "A2", "b1": "B"},   # A2_B prior 1
            # Two dual-alpha cells with the same (A1,A2,B) triple -> merge.
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},
            {"bc": "d2", "a1": "A2", "a2": "A1", "b1": "B"},
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_min_cells=2))
        assert merge.get("A2_B") == "A1_B"
        assert merge.get("A1_B") == "A1_B"
        # Dual cells attributed wholly to the merged clone.
        for bc in ("d1", "d2"):
            rows = long[long["cell_barcode"] == bc]
            assert list(rows["CDR3ab"]) == ["A1_B"]
            assert rows["weight"].iloc[0] == pytest.approx(1.0)
            assert rows["kind"].iloc[0] == "dual_alpha_merge"
        # A2_B no longer appears as its own clone (merged away).
        assert "A2_B" not in set(long["CDR3ab"])

    def test_singleton_dual_alpha_splits(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},
            {"a1": "A1", "b1": "B"},   # A1_B prior 3
            {"a1": "A2", "b1": "B"},   # A2_B prior 1
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},  # singleton -> split
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_min_cells=2))
        assert merge == {}
        rows = long[long["cell_barcode"] == "d1"].set_index("CDR3ab")["weight"]
        assert rows["A1_B"] == pytest.approx(0.75)
        assert rows["A2_B"] == pytest.approx(0.25)

    def test_merge_disabled_splits(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "d1", "a1": "A1", "a2": "A2", "b1": "B"},
            {"bc": "d2", "a1": "A1", "a2": "A2", "b1": "B"},
        ])
        long, merge = attribute_cells(df, _on(dual_alpha_merge=False))
        assert merge == {}
        assert (long[long["cell_barcode"] == "d1"]["kind"] == "dual_alpha_split").all()


class TestDualBeta:
    def test_dual_beta_splits(self):
        df = _cells([
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B1"},
            {"a1": "A", "b1": "B1"},   # A_B1 prior 3
            {"a1": "A", "b1": "B2"},   # A_B2 prior 1
            {"bc": "d1", "a1": "A", "b1": "B1", "b2": "B2"},
        ])
        long, _ = attribute_cells(df, _on())
        rows = long[long["cell_barcode"] == "d1"].set_index("CDR3ab")["weight"]
        assert rows["A_B1"] == pytest.approx(0.75)
        assert rows["A_B2"] == pytest.approx(0.25)


class TestWeightConservation:
    def test_every_cell_sums_to_one_or_zero(self):
        df = _cells([
            {"a1": "A1", "b1": "B"},
            {"a1": "A2", "b1": "B"},
            {"bc": "drop_ok", "a1": "", "b1": "B"},       # -> 1.0
            {"bc": "drop_no", "a1": "", "b1": "ZZZ"},     # -> 0.0 (no match)
            {"bc": "dual", "a1": "A1", "a2": "A2", "b1": "B"},  # -> 1.0
            {"bc": "nochain", "a1": "", "b1": ""},        # -> 0.0
            {"bc": "alpha_only", "a1": "A1", "b1": ""},   # -> 0.0 (no beta)
        ])
        long, _ = attribute_cells(df, _on())
        sums = long.groupby("cell_barcode")["weight"].sum()
        for bc in ("c0", "c1", "drop_ok", "dual"):
            assert sums.get(bc, 0.0) == pytest.approx(1.0), bc
        for bc in ("drop_no", "nochain", "alpha_only"):
            assert sums.get(bc, 0.0) == pytest.approx(0.0), bc

    def test_umi_gate_excludes_low_umi_chain(self):
        # Beta present but below UMI gate -> treated as no beta -> dropped.
        df = _cells([
            {"a1": "A", "b1": "B"},
            {"bc": "low", "a1": "A", "b1": "B", "bu1": 1},
        ])
        long, _ = attribute_cells(df, _on(), min_umi=2)
        assert "low" not in set(long["cell_barcode"])


class TestDisabledAndGuards:
    def test_cdr3b_only_returns_empty(self):
        df = _cells([{"a1": "A", "b1": "B"}])
        long, merge = attribute_cells(df, _on(), group_by="CDR3b_only")
        assert long.empty and merge == {}

    def test_empty_input(self):
        df = _cells([])
        long, merge = attribute_cells(df, _on())
        assert long.empty and merge == {}
