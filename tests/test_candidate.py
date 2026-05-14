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

"""Tests for ``select_candidates`` (#75) — signature-based shortlist."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift.candidate import select_candidates


def _clones(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


class TestSelectCandidates:
    def test_tier1_and_tier2_always_selected(self):
        df = _clones([
            {"CDR3ab": "A", "tier": "tier1", "signature_antigen_response": 0.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "B", "tier": "tier2", "signature_antigen_response": 0.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "C", "tier": "tier5", "signature_antigen_response": 0.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        out = select_candidates(df, top_n=0)
        # top_n=0 so no signature picks; only tier1/tier2 selected.
        assert out.set_index("CDR3ab")["is_selected"].to_dict() == {
            "A": True, "B": True, "C": False,
        }

    def test_top_n_per_signature_from_tier3plus(self):
        # 5 tier3 clones with varying antigen_response scores; top 2 picked.
        df = _clones([
            {"CDR3ab": f"C{i}", "tier": "tier3",
             "signature_antigen_response": float(i),
             "signature_cytolytic": 0.0,
             "signature_tumor_reactive": 0.0}
            for i in range(5)
        ])
        out = select_candidates(df, top_n=2).set_index("CDR3ab")
        # Top 2 by antigen_response: C4 (5th) and C3 (4th).
        assert out.loc["C4", "tier3plus_top_by_antigen_response_signature"]
        assert out.loc["C3", "tier3plus_top_by_antigen_response_signature"]
        assert not out.loc["C2", "tier3plus_top_by_antigen_response_signature"]
        assert out.loc["C4", "is_selected"]
        assert out.loc["C3", "is_selected"]
        assert not out.loc["C2", "is_selected"]

    def test_tier1_excluded_from_signature_pool(self):
        """Tier1 clones with high signature score should NOT consume a
        top-N slot in the tier3+ pool — the pool is strictly tier3+."""
        df = _clones([
            # Tier1 clone with the highest signature score.
            {"CDR3ab": "Q", "tier": "tier1",
             "signature_antigen_response": 999.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            # Tier3 clones with lower scores.
            {"CDR3ab": "C0", "tier": "tier3",
             "signature_antigen_response": 1.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
            {"CDR3ab": "C1", "tier": "tier3",
             "signature_antigen_response": 2.0,
             "signature_cytolytic": 0.0, "signature_tumor_reactive": 0.0},
        ])
        out = select_candidates(df, top_n=1).set_index("CDR3ab")
        # Q is selected via tier1, not via signature.
        assert out.loc["Q", "tier1_or_tier2"]
        assert not out.loc["Q", "tier3plus_top_by_antigen_response_signature"]
        # The tier3+ slot goes to C1 (highest of the tier3 pool).
        assert out.loc["C1", "tier3plus_top_by_antigen_response_signature"]
        assert not out.loc["C0", "tier3plus_top_by_antigen_response_signature"]

    def test_different_signatures_can_pick_different_clones(self):
        df = _clones([
            {"CDR3ab": "A", "tier": "tier3",
             "signature_antigen_response": 10.0,
             "signature_cytolytic": 0.0,
             "signature_tumor_reactive": 0.0},
            {"CDR3ab": "B", "tier": "tier3",
             "signature_antigen_response": 0.0,
             "signature_cytolytic": 10.0,
             "signature_tumor_reactive": 0.0},
            {"CDR3ab": "C", "tier": "tier3",
             "signature_antigen_response": 0.0,
             "signature_cytolytic": 0.0,
             "signature_tumor_reactive": 10.0},
        ])
        out = select_candidates(df, top_n=1).set_index("CDR3ab")
        assert out.loc["A", "tier3plus_top_by_antigen_response_signature"]
        assert out.loc["B", "tier3plus_top_by_cytolytic_signature"]
        assert out.loc["C", "tier3plus_top_by_tumor_reactive_signature"]
        assert all(out["is_selected"])

    def test_missing_score_column_emits_false(self):
        df = _clones([
            {"CDR3ab": "A", "tier": "tier3"},
        ])
        out = select_candidates(df, top_n=1).set_index("CDR3ab")
        assert not out.loc["A", "tier3plus_top_by_antigen_response_signature"]
        assert not out.loc["A", "is_selected"]

    def test_nan_scores_dont_pick(self):
        df = _clones([
            {"CDR3ab": "A", "tier": "tier3",
             "signature_antigen_response": np.nan,
             "signature_cytolytic": 0.0,
             "signature_tumor_reactive": 0.0},
        ])
        out = select_candidates(df, top_n=1).set_index("CDR3ab")
        assert not out.loc["A", "tier3plus_top_by_antigen_response_signature"]

    def test_missing_tier_column_raises(self):
        df = _clones([{"CDR3ab": "A"}])
        with pytest.raises(ValueError, match="tier"):
            select_candidates(df)

    def test_does_not_mutate_input(self):
        df = _clones([
            {"CDR3ab": "A", "tier": "tier1",
             "signature_antigen_response": 0.0,
             "signature_cytolytic": 0.0,
             "signature_tumor_reactive": 0.0},
        ])
        before = df.columns.tolist()
        select_candidates(df)
        assert df.columns.tolist() == before
