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

"""Tests for the composable in-silico filter layer (#149)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift import insilico_filter as isf


class TestPercentileRank:
    def test_lower_is_better_zero_is_best(self):
        # Ppost^low: lowest value ranks best (0). Ppost==0 is a valid best.
        s = pd.Series([0.0, 1.0, 2.0, 3.0])
        pr = isf.percentile_rank(s, lower_is_better=True)
        assert pr.iloc[0] == 0.0
        assert pr.iloc[-1] == 1.0
        assert pr.is_monotonic_increasing

    def test_high_is_better_flips(self):
        s = pd.Series([0.0, 1.0, 2.0, 3.0])
        pr = isf.percentile_rank(s, lower_is_better=False)
        assert pr.iloc[-1] == 0.0  # highest is best
        assert pr.iloc[0] == 1.0

    def test_ties_grouped_not_dropped(self):
        # A cluster of equal zeros shares a rank (method='average'), and is
        # never dropped — the Ppost==0 case the issue calls out.
        s = pd.Series([0.0, 0.0, 0.0, 1.0, 2.0])
        pr = isf.percentile_rank(s, lower_is_better=True, method="average")
        # the three zeros share one averaged rank, all strictly best
        assert pr.iloc[0] == pr.iloc[1] == pr.iloc[2]
        assert pr.iloc[0] < pr.iloc[3] < pr.iloc[4]
        assert pr.notna().all()

    def test_nan_stays_nan(self):
        s = pd.Series([0.0, np.nan, 2.0])
        pr = isf.percentile_rank(s, lower_is_better=True)
        assert np.isnan(pr.iloc[1])
        assert pr.iloc[0] == 0.0

    def test_single_value_is_best(self):
        assert isf.percentile_rank(pd.Series([5.0])).iloc[0] == 0.0


class TestInsilicoMask:
    def _df(self):
        return pd.DataFrame({
            "CDR3ab": [f"c{i}" for i in range(6)],
            "log10_ppost": [-12, -11, -5, -4, -3, np.nan],
            "sig_response": [9, 8, 1, 0, 5, 5],
        })

    def test_single_predicate_keeps_best_fraction(self):
        preds = [isf.FilterPredicate("log10_ppost", "low", 0.5)]
        mask = isf.insilico_mask(self._df(), preds)
        # bottom-half (lowest ppost), NaN excluded
        assert list(mask) == [True, True, True, False, False, False]

    def test_stacked_predicates_all_must_pass(self):
        # Ppost^low (best half) AND response^high (best half).
        preds = [
            isf.FilterPredicate("log10_ppost", "low", 0.5),
            isf.FilterPredicate("sig_response", "high", 0.5),
        ]
        mask = isf.insilico_mask(self._df(), preds)
        # rows passing BOTH: low ppost (c0,c1,c2) and high response.
        # response ranks (high best): c0=9,c1=8 best; c2=1 worst-ish.
        assert mask.iloc[0] and mask.iloc[1]
        assert not mask.iloc[5]  # NaN ppost never passes

    def test_missing_score_never_passes(self):
        preds = [isf.FilterPredicate("log10_ppost", "low", 1.0)]
        mask = isf.insilico_mask(self._df(), preds)
        assert not mask.iloc[5]  # the NaN row

    def test_missing_column_raises(self):
        preds = [isf.FilterPredicate("nope", "low", 0.5)]
        with pytest.raises(ValueError, match="not in frame"):
            isf.insilico_mask(self._df(), preds)

    def test_empty_predicates_all_pass(self):
        mask = isf.insilico_mask(self._df(), [])
        assert mask.all()

    def test_within_group_ranking(self):
        # Two assay groups with different score scales; the filter refines
        # each against its own distribution.
        df = pd.DataFrame({
            "CDR3ab": list("abcdef"),
            "group": ["A", "A", "A", "B", "B", "B"],
            "ppost": [-10, -9, -8, -3, -2, -1],
        })
        preds = [isf.FilterPredicate("ppost", "low", 0.4)]
        mask = isf.insilico_mask(df, preds, group_col="group")
        # best ~40% within each group: the lowest-ppost member of each.
        assert mask.iloc[0] and not mask.iloc[2]   # group A
        assert mask.iloc[3] and not mask.iloc[5]   # group B


class TestPredicateValidation:
    def test_bad_direction(self):
        with pytest.raises(ValueError, match="low.*high"):
            isf.FilterPredicate("x", "sideways", 0.5)

    def test_bad_percentile(self):
        with pytest.raises(ValueError, match="percentile"):
            isf.FilterPredicate("x", "low", 1.5)

    def test_predicates_from_config(self):
        cfg = {"predicates": [
            {"score": "log10_ppost", "direction": "low", "percentile": 0.3},
            {"score": "sig_response", "direction": "high"},
        ]}
        preds = isf.predicates_from_config(cfg)
        assert len(preds) == 2
        assert preds[0].percentile == 0.3
        assert preds[1].direction == "high"
        assert preds[1].percentile == 0.5  # default

    def test_predicate_missing_score_raises(self):
        with pytest.raises(ValueError, match="missing 'score'"):
            isf.predicates_from_config({"predicates": [{"direction": "low"}]})


class TestApplyAndTwins:
    def _df(self):
        return pd.DataFrame({
            "CDR3ab": list("abcdef"),
            "method": ["CTYneg"] * 3 + ["AIMpos"] * 3,
            "cells": [10, 5, 1, 8, 4, 2],
            "ppost": [-10, -9, -1, -11, -2, -1],
        })

    def test_apply_adds_columns_and_named_composite(self):
        preds = [isf.FilterPredicate("ppost", "low", 0.5)]
        out = isf.apply_insilico_filter(
            self._df(), preds, group_col="method", label="IS",
        )
        assert {"insilico_pass", "insilico_avg_rank", "insilico_group"} <= set(out.columns)
        # passing rows are named "<method>+IS"
        passing = out[out["insilico_pass"]]
        assert set(passing["insilico_group"]) <= {"CTYneg+IS", "AIMpos+IS"}
        # non-passing → NaN composite
        assert out.loc[~out["insilico_pass"], "insilico_group"].isna().all()

    def test_expand_twins_doubles_groups(self):
        preds = [isf.FilterPredicate("ppost", "low", 0.5)]
        expanded = isf.expand_insilico_twins(
            self._df(), preds, group_col="method", label="IS",
        )
        groups = set(expanded["overlap_group"])
        assert {"CTYneg", "AIMpos"} <= groups
        assert {"CTYneg+IS", "AIMpos+IS"} <= groups
        # twin rows are a subset of the base rows
        assert (expanded["overlap_group"] == "CTYneg+IS").sum() <= 3

    def test_average_rank_for_ranking(self):
        preds = [
            isf.FilterPredicate("ppost", "low", 1.0),
        ]
        ar = isf.average_percentile_rank(self._df(), preds, group_col="method")
        # within CTYneg, 'a' (lowest ppost) is best (0.0)
        assert ar.iloc[0] == 0.0


class TestInsilicoOverlap:
    def test_overlap_includes_twin_groups(self):
        df = pd.DataFrame({
            "CDR3ab": ["x", "y", "z", "x", "w"],
            "method": ["CTYneg", "CTYneg", "CTYneg", "AIMpos", "AIMpos"],
            "cells": [10, 5, 1, 9, 3],
            "ppost": [-10, -9, -1, -11, -2],
        })
        preds = [isf.FilterPredicate("ppost", "low", 0.6)]
        mats = isf.insilico_overlap_long(df, preds, group_col="method")
        jac = mats["jaccard"]
        # both base and +IS groups appear as overlap axes
        assert "CTYneg" in jac.index and "CTYneg+IS" in jac.index
        assert "AIMpos+IS" in jac.index
        # diagonal is 1.0
        assert all(jac.loc[g, g] == 1.0 for g in jac.index)


class TestConfigIntegration:
    def test_insilico_filter_block_parsed(self):
        from tcrsift.config import TCRsiftConfig

        cfg = TCRsiftConfig._from_dict({
            "insilico_filter": {
                "label": "IS",
                "predicates": [{"score": "log10_ppost", "direction": "low"}],
            }
        })
        assert cfg.insilico_filter["label"] == "IS"
        preds = isf.predicates_from_config(cfg.insilico_filter)
        assert preds[0].score == "log10_ppost"

    def test_default_is_empty(self):
        from tcrsift.config import TCRsiftConfig

        assert TCRsiftConfig().insilico_filter == {}
