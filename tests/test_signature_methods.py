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

"""Tests for context-aware gene-signature selection methods."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift import signature_methods as sm


class TestSignatureRegistry:
    def test_invariant_signatures_valid_everywhere(self):
        prolif = sm.SIGNATURES["Proliferation"]
        assert prolif.contexts == frozenset()
        for ctx in (sm.CULTURE, sm.CIRCULATING, sm.TISSUE, None):
            assert prolif.valid_in(ctx)

    def test_tissue_signature_invalid_in_culture(self):
        til = sm.SIGNATURES["NeoantigenExperiencedTIL"]
        assert til.contexts == frozenset({sm.TISSUE})
        assert til.valid_in(sm.TISSUE)
        assert not til.valid_in(sm.CULTURE)

    def test_signed_signature_all_genes_dedup(self):
        diff = sm.SIGNATURES["Differentiated"]
        assert diff.genes_up == ()
        assert "CCR7" in diff.genes_down
        assert len(diff.all_genes) == len(set(diff.all_genes))

    def test_small_and_large_panels_exist(self):
        assert sm.SIGNATURES["Cytolytic"].panel == "focal"
        assert sm.SIGNATURES["CytolyticBroad"].panel == "broad"
        assert len(sm.SIGNATURES["CytolyticBroad"].genes_up) > len(
            sm.SIGNATURES["Cytolytic"].genes_up
        )


class TestInferContext:
    def test_known_sources(self):
        assert sm.infer_context("culture") == sm.CULTURE
        assert sm.infer_context("TIL") == sm.TISSUE
        assert sm.infer_context("tetramer") == sm.CIRCULATING

    def test_unknown_source_is_none(self):
        assert sm.infer_context("mystery") is None
        assert sm.infer_context(None) is None


class TestScoreSignature:
    def _expr(self):
        return pd.DataFrame(
            {"GZMB": [0.0, 0.0, 5.0, 5.0], "PRF1": [0.0, 1.0, 4.0, 6.0],
             "CCR7": [5.0, 5.0, 0.0, 0.0]},
            index=[f"c{i}" for i in range(4)],
        )

    def test_signed_score_subtracts_down_genes(self):
        # Differentiated = -(CCR7): naive-high cells score low, others high.
        s = sm.score_signature(self._expr(), sm.SIGNATURES["Differentiated"],
                               min_genes_present=1)
        assert s["c3"] > s["c0"]

    def test_combine_mean_uses_raw_values(self):
        s = sm.score_signature(self._expr(), sm.SIGNATURES["Cytolytic"], combine="mean")
        # raw mean of GZMB,PRF1: c3=(5+6)/2=5.5 highest
        assert s["c3"] == 5.5

    def test_combine_zscore_default(self):
        s = sm.score_signature(self._expr(), sm.SIGNATURES["Cytolytic"])
        assert s["c3"] > s["c0"]

    def test_unknown_combine_raises(self):
        with pytest.raises(ValueError, match="combine"):
            sm.score_signature(self._expr(), sm.SIGNATURES["Cytolytic"], combine="bogus")

    def test_too_few_present_returns_zeros(self):
        expr = pd.DataFrame({"GZMB": [1.0, 2.0]}, index=["a", "b"])
        s = sm.score_signature(expr, sm.SIGNATURES["Cytolytic"], min_genes_present=2)
        assert (s == 0).all()  # only GZMB present (<2)

    def test_explicit_background_changes_zscore(self):
        # z-scoring against a low-expression background lifts all scores
        # vs z-scoring against the (higher) input itself.
        expr = self._expr()
        low_bg = pd.DataFrame(
            {"GZMB": [0.0, 0.0, 1.0, 1.0], "PRF1": [0.0, 1.0, 0.0, 1.0]},
            index=expr.index,
        )
        self_scored = sm.score_signature(expr, sm.SIGNATURES["Cytolytic"])
        bg_scored = sm.score_signature(
            expr, sm.SIGNATURES["Cytolytic"], background=low_bg,
        )
        assert bg_scored.mean() > self_scored.mean()


class TestCallPositiveAdaptive:
    def _scores(self):
        vals = list(np.linspace(0.0, 1.0, 16)) + [5.0, 5.2, 4.8, 5.1]
        return pd.Series(vals, index=[f"cl{i}" for i in range(20)])

    def test_gap_isolates_cluster(self):
        assert int(sm.call_positive(self._scores(), method="gap").sum()) == 4

    def test_otsu_isolates_cluster(self):
        assert int(sm.call_positive(self._scores(), method="otsu").sum()) == 4

    def test_quantile_fixed(self):
        assert int(sm.call_positive(self._scores(), method="quantile", quantile=0.75).sum()) == 5

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="unknown method"):
            sm.call_positive(self._scores(), method="bogus")


class TestContextGuard:
    def _data(self):
        cells = [f"c{i}" for i in range(4)]
        expr = pd.DataFrame(
            {"ENTPD1": [5.0, 5.0, 0.0, 0.0], "CXCL13": [4.0, 6.0, 0.0, 1.0]},
            index=cells,
        )
        obs = pd.DataFrame({"CDR3ab": ["A", "A", "B", "B"], "sample": ["S1"] * 4},
                          index=cells)
        return expr, obs

    def test_default_warns_but_still_scores(self):
        # Off-context is interesting, not forbidden: default warns + scores.
        expr, obs = self._data()
        out = sm.build_signature_methods(
            expr, obs, signatures=["NeoantigenExperiencedTIL"],
            context=sm.CULTURE,  # mismatch; default on_context_mismatch="warn"
        )
        assert "NeoantigenExperiencedTIL" in set(out["method"])

    def test_raise_is_opt_in(self):
        expr, obs = self._data()
        with pytest.raises(sm.SignatureContextError, match="tissue"):
            sm.build_signature_methods(
                expr, obs, signatures=["NeoantigenExperiencedTIL"],
                context=sm.CULTURE, on_context_mismatch="raise",
            )

    def test_tissue_signature_on_tissue_runs(self):
        expr, obs = self._data()
        out = sm.build_signature_methods(
            expr, obs, signatures=["NeoantigenExperiencedTIL"],
            context=sm.TISSUE, positive_method="gap",
        )
        assert set(out["method"]) <= {"NeoantigenExperiencedTIL"}


class TestSignatureMethodsLong:
    def test_emits_method_rows_with_frequency(self):
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "B", "B"], "sample": ["S1"] * 4},
            index=[f"c{i}" for i in range(4)],
        )
        positive = {"Cytolytic": pd.Series([True, True, True, False], index=obs.index)}
        out = sm.signature_methods_long(obs, positive)
        a = out[out["CDR3ab"] == "A"].iloc[0]
        assert a["cells"] == 2
        assert abs(a["frequency"] - 2 / 3) < 1e-9

    def test_feeds_selection_rules(self):
        from tcrsift.selection import select_from_clone_sample_long

        cells = [f"c{i}" for i in range(6)]
        expr = pd.DataFrame(
            {"PRF1": [5, 5, 5, 0, 0, 0], "GZMB": [5, 5, 5, 0, 0, 0]}, index=cells,
        )
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "A", "B", "B", "B"], "sample": ["S1"] * 6}, index=cells,
        )
        sig_long = sm.build_signature_methods(
            expr, obs, signatures=["Cytolytic"], positive_method="gap",
        )
        rules = select_from_clone_sample_long(
            sig_long,
            {"rules": {"private": {"include_tier": "tier3", "top_n": 5,
                                   "apply_to_methods": ["Cytolytic"]}}},
        )
        assert set(rules.columns) >= {"CDR3ab", "selection_rule", "global_rank"}
