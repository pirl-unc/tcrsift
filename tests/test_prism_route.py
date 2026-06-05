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

"""GEX-inclusive percentile-rank (PRISM) ranking in the select route engine (#193)."""

from __future__ import annotations

import pandas as pd

from tcrsift.selection import (
    build_selection_rules,
    percentile_rank_score,
    select_from_clone_sample_long,
)

# A=private+antigen-responding+non-naive (best PRISM); B=opposite (worst); C=middle.
_SCORES = pd.DataFrame({
    "CDR3ab": ["A", "B", "C"],
    "ppost_alpha": [1e-9, 1e-5, 1e-7],          # low is best
    "ppost_beta": [1e-9, 1e-5, 1e-7],
    "antigen_response_score": [3.0, 0.1, 1.0],  # high is best
    "naive_score": [0.1, 3.0, 1.0],             # low is best
})
_PRISM_TERMS = [
    {"col": "ppost_alpha", "direction": "low"},
    {"col": "ppost_beta", "direction": "low"},
    {"col": "antigen_response_score", "direction": "high"},
    {"col": "naive_score", "direction": "low"},
]


class TestPercentileRankScore:
    def test_prism_composite_ordering(self):
        comp = percentile_rank_score(_SCORES, _PRISM_TERMS)
        assert comp["A"] < comp["C"] < comp["B"]  # low composite = best

    def test_weight_and_direction(self):
        # Only the antigen-response term (high is best): B (0.1) worst, A best.
        comp = percentile_rank_score(
            _SCORES, [{"col": "antigen_response_score", "direction": "high"}]
        )
        assert comp["A"] < comp["B"]

    def test_empty_inputs(self):
        assert percentile_rank_score(pd.DataFrame(), _PRISM_TERMS) == {}
        assert percentile_rank_score(_SCORES, []) == {}


class TestPrismRoute:
    def _cml(self):
        return pd.DataFrame({
            "CDR3ab": ["A", "B", "C"],
            "method": ["m", "m", "m"],
            "tier": ["tier1", "tier1", "tier1"],
            "max_freq_in_method": [0.1, 0.1, 0.1],
            "cells_in_method": [10, 10, 10],
        })

    def test_shared_route_ranks_by_prism(self):
        config = {"rules": {"shared": {
            "include_tiers": ["tier1", "tier2"],
            "rank_by": {"percentile_rank": _PRISM_TERMS},
        }}}
        out = build_selection_rules(self._cml(), config, clone_scores=_SCORES)
        # Ranked by PRISM composite ascending: A (best) -> C -> B.
        assert list(out["CDR3ab"]) == ["A", "C", "B"]
        assert (out["ranking_metric"] == "percentile_rank").all()
        assert list(out["global_rank"]) == [1, 2, 3]

    def test_freq_route_unchanged_without_scores(self):
        # No clone_scores + a frequency rank_by -> existing behavior (desc freq).
        cml = self._cml()
        cml["max_freq_in_method"] = [0.05, 0.2, 0.1]
        config = {"rules": {"shared": {
            "include_tiers": ["tier1"], "rank_by": "max_frequency",
        }}}
        out = build_selection_rules(cml, config)
        assert list(out["CDR3ab"]) == ["B", "C", "A"]  # highest freq first
        assert (out["ranking_metric"] == "max_frequency").all()

    def test_end_to_end_with_clone_scores(self):
        long = pd.DataFrame({
            "CDR3ab": ["A", "B", "C"],
            "sample": ["S1", "S1", "S1"],
            "method": ["m", "m", "m"],
            "cells": [10, 10, 10],
            "frequency": [0.1, 0.1, 0.1],
        })
        config = {"rules": {"shared": {
            "include_tiers": ["tier1", "tier2"],
            "rank_by": {"percentile_rank": _PRISM_TERMS},
        }}}
        out = select_from_clone_sample_long(long, config, clone_scores=_SCORES)
        assert list(out["CDR3ab"]) == ["A", "C", "B"]
