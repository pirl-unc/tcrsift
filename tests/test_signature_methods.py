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

"""Tests for gene-signature synthetic selection methods."""

from __future__ import annotations

import pandas as pd

from tcrsift import signature_methods as sm


class TestSelectionSignatures:
    def test_two_axes_defined(self):
        assert set(sm.SELECTION_SIGNATURES) == {"TumorReactive", "AntigenExperienced"}

    def test_tumor_reactive_composes_canonical_sets(self):
        tr = set(sm.TUMOR_REACTIVE_SIGNATURE)
        assert {"CXCL13", "ENTPD1"} <= tr      # tumor_reactive
        assert {"PDCD1", "LAG3", "TOX"} <= tr  # exhaustion
        assert "ITGAE" in tr                   # CD103

    def test_antigen_experienced_composes_canonical_sets(self):
        ae = set(sm.ANTIGEN_EXPERIENCED_SIGNATURE)
        assert {"TNFRSF9", "MKI67"} <= ae      # antigen_response
        assert {"IFNG", "GZMB", "NKG7"} <= ae  # activation

    def test_no_duplicate_genes(self):
        for genes in sm.SELECTION_SIGNATURES.values():
            assert len(genes) == len(set(genes))


class TestScoreSignature:
    def _expr(self):
        cells = [f"c{i}" for i in range(4)]
        return pd.DataFrame(
            {"GZMB": [0.0, 0.0, 5.0, 5.0], "NKG7": [0.0, 1.0, 4.0, 6.0]},
            index=cells,
        )

    def test_zscored_mean_orders_cells(self):
        s = sm.score_signature(self._expr(), ["GZMB", "NKG7"])
        # high-expression cells score above low-expression cells
        assert s["c3"] > s["c0"]

    def test_missing_genes_dropped_but_scores(self):
        s = sm.score_signature(self._expr(), ["GZMB", "NKG7", "ABSENT"])
        assert len(s) == 4  # scored on the 2 present genes

    def test_too_few_present_returns_zeros(self):
        s = sm.score_signature(self._expr(), ["ONLY_ONE_PRESENT", "GZMB"], min_genes_present=2)
        # only GZMB present (<2) -> zeros
        assert (s == 0).all()


class TestCallPositive:
    def test_quantile_default_top_quartile(self):
        scores = pd.Series([1.0, 2.0, 3.0, 4.0], index=list("abcd"))
        pos = sm.call_positive(scores, quantile=0.75)
        assert pos["d"] and not pos["a"]

    def test_explicit_threshold(self):
        scores = pd.Series([1.0, 2.0, 3.0], index=list("abc"))
        pos = sm.call_positive(scores, threshold=1.5)
        assert list(pos) == [False, True, True]


class TestSignatureMethodsLong:
    def test_emits_method_rows_with_frequency(self):
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "B", "B"], "sample": ["S1"] * 4},
            index=[f"c{i}" for i in range(4)],
        )
        positive = {"TumorReactive": pd.Series([True, True, True, False], index=obs.index)}
        out = sm.signature_methods_long(obs, positive)
        assert set(out["method"]) == {"TumorReactive"}
        a = out[out["CDR3ab"] == "A"].iloc[0]
        assert a["cells"] == 2
        # 3 positive cells in S1; A has 2 -> 2/3
        assert abs(a["frequency"] - 2 / 3) < 1e-9

    def test_empty_when_no_positive(self):
        obs = pd.DataFrame({"CDR3ab": ["A"], "sample": ["S1"]}, index=["c0"])
        out = sm.signature_methods_long(
            obs, {"X": pd.Series([False], index=["c0"])}
        )
        assert out.empty
        assert list(out.columns) == ["CDR3ab", "sample", "method", "cells", "frequency"]


class TestBuildSignatureMethodsIntoSelection:
    def test_signature_method_feeds_selection_rules(self):
        from tcrsift.selection import select_from_clone_sample_long

        cells = [f"c{i}" for i in range(6)]
        expr = pd.DataFrame(0.0, index=cells, columns=["CXCL13", "ENTPD1", "TOX", "PDCD1", "LAG3"])
        expr.loc[["c0", "c1", "c2"], :] = 5.0  # clone A tumor-reactive-high
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "A", "B", "B", "B"], "sample": ["S1"] * 6},
            index=cells,
        )
        sig_long = sm.build_signature_methods(
            expr, obs, signatures={"TumorReactive": sm.TUMOR_REACTIVE_SIGNATURE},
            quantile=0.4,
        )
        # The synthetic method is usable by a private rule.
        rules = select_from_clone_sample_long(
            sig_long,
            {"rules": {"private": {"include_tier": "tier3", "top_n": 5,
                                   "apply_to_methods": ["TumorReactive"]}}},
        )
        # A is tumor-reactive-positive; with enough cells it can clear a tier
        # and be selected by the TumorReactive private rule (or none if too
        # sparse) — the point is the pipeline runs end-to-end with a
        # signature acting as a method.
        assert set(rules.columns) >= {"CDR3ab", "selection_rule", "global_rank"}


class TestAdaptiveCutoffs:
    """Adaptive cutoffs select the variable, small high cluster (#sig)."""

    def _scores(self):
        # 16 smooth background in [0,1], 4 clearly clustered high.
        import numpy as np
        vals = list(np.linspace(0.0, 1.0, 16)) + [5.0, 5.2, 4.8, 5.1]
        return pd.Series(vals, index=[f"cl{i}" for i in range(20)])

    def test_gap_isolates_the_high_cluster(self):
        pos = sm.call_positive(self._scores(), method="gap")
        assert int(pos.sum()) == 4
        assert set(self._scores().index[pos]) == {"cl16", "cl17", "cl18", "cl19"}

    def test_otsu_isolates_the_high_cluster(self):
        pos = sm.call_positive(self._scores(), method="otsu")
        assert int(pos.sum()) == 4

    def test_quantile_is_fixed_count(self):
        # top quartile of 20 = 5, regardless of where the gap is
        pos = sm.call_positive(self._scores(), method="quantile", quantile=0.75)
        assert int(pos.sum()) == 5

    def test_gap_count_varies_with_cluster_size(self):
        # Two high clones instead of four -> gap picks two.
        import numpy as np
        vals = list(np.linspace(0.0, 1.0, 16)) + [5.0, 5.1]
        s = pd.Series(vals, index=[f"c{i}" for i in range(18)])
        assert int(sm.call_positive(s, method="gap").sum()) == 2

    def test_explicit_threshold_overrides_method(self):
        pos = sm.call_positive(self._scores(), method="gap", threshold=2.0)
        assert int(pos.sum()) == 4  # the four >2.0

    def test_unknown_method_raises(self):
        import pytest
        with pytest.raises(ValueError, match="unknown method"):
            sm.call_positive(self._scores(), method="bogus")
