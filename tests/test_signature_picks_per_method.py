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

"""Tests for #84: ``compute_signature_picks_per_method`` and
``signature_picks_clone_to_methods``."""

from __future__ import annotations

import pandas as pd

from tcrsift.candidate import (
    compute_signature_picks_per_method,
    signature_picks_clone_to_methods,
)


def _frame(rows):
    return pd.DataFrame(rows)


class TestComputeSignaturePicksPerMethod:
    def test_top_n_per_method_per_signature(self):
        df = _frame([
            {"CDR3ab": "A", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": 5.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "B", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": 3.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "C", "method": "M2", "cells_in_method": 10,
             "signature_antigen_response": 7.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        picks = compute_signature_picks_per_method(df, top_n=1)
        # Within M1, A has higher antigen_response than B → A is picked.
        assert picks["M1"]["antigen_response"] == ["A"]
        # Within M2, only C → C picked.
        assert picks["M2"]["antigen_response"] == ["C"]

    def test_min_cells_filter(self):
        df = _frame([
            # B has higher score but fails the min_cells floor.
            {"CDR3ab": "A", "method": "M1", "cells_in_method": 5,
             "signature_antigen_response": 1.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "B", "method": "M1", "cells_in_method": 1,
             "signature_antigen_response": 10.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        picks = compute_signature_picks_per_method(df, top_n=1, min_cells=2)
        assert picks["M1"]["antigen_response"] == ["A"]

    def test_restrict_to_pool(self):
        df = _frame([
            {"CDR3ab": "TIER1_HIT", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": 99.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "TIER3_HIT", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": 5.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        # Only tier3+ pool — the high-scoring tier1 clone is excluded.
        picks = compute_signature_picks_per_method(
            df, pool_clones={"TIER3_HIT"}, top_n=1
        )
        assert picks["M1"]["antigen_response"] == ["TIER3_HIT"]

    def test_nan_scores_drop_to_bottom(self):
        df = _frame([
            {"CDR3ab": "A", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": float("nan"),
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "B", "method": "M1", "cells_in_method": 10,
             "signature_antigen_response": 1.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        picks = compute_signature_picks_per_method(df, top_n=1)
        assert picks["M1"]["antigen_response"] == ["B"]

    def test_missing_score_column_returns_empty_list(self):
        df = _frame([
            {"CDR3ab": "A", "method": "M1", "cells_in_method": 10},
        ])
        picks = compute_signature_picks_per_method(
            df, signatures=("antigen_response",), top_n=1
        )
        assert picks["M1"]["antigen_response"] == []

    def test_empty_input(self):
        picks = compute_signature_picks_per_method(
            pd.DataFrame(columns=["CDR3ab", "method", "cells_in_method"])
        )
        assert picks == {}


class TestSignaturePicksCloneToMethods:
    def test_inverts_picks(self):
        per_method = {
            "M1": {"antigen_response": ["A", "B"], "cytolytic": ["C"]},
            "M2": {"antigen_response": ["A"], "cytolytic": ["D"]},
        }
        out = signature_picks_clone_to_methods(per_method)
        assert sorted(out["A"]["antigen_response"]) == ["M1", "M2"]
        assert out["B"]["antigen_response"] == ["M1"]
        assert out["C"]["cytolytic"] == ["M1"]
        assert out["D"]["cytolytic"] == ["M2"]

    def test_empty_input(self):
        assert signature_picks_clone_to_methods({}) == {}
