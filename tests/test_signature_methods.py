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

import numpy as np
import pandas as pd
import pytest

from tcrsift import signature_methods as sm


class TestSignatureRegistry:
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

    def test_no_context_tagging_on_signatures(self):
        # Tagging was removed — signatures are plain gene sets now.
        sig = sm.SIGNATURES["TumorReactive"]
        assert not hasattr(sig, "contexts")
        assert not hasattr(sm, "infer_context")

    def test_selection_signatures_named_composites(self):
        # The two headline composites the pilot uses (no context tags).
        assert set(sm.SELECTION_SIGNATURES) == {"TumorReactive", "AntigenExperienced"}
        tr = sm.SELECTION_SIGNATURES["TumorReactive"]
        assert {"CXCL13", "ENTPD1", "TOX", "ITGAE", "CTLA4"} <= set(tr.genes_up)
        ae = sm.SELECTION_SIGNATURES["AntigenExperienced"]
        assert {"TNFRSF9", "MKI67", "GZMB", "NKG7"} <= set(ae.genes_up)

    def test_single_gene_focal_signature_scores(self):
        # Regression: a 1-gene focal panel (Proliferation) must score, not
        # silently return zeros under the default min_genes_present.
        expr = pd.DataFrame(
            {"MKI67": [0.0, 1.0, 5.0, 9.0]}, index=[f"c{i}" for i in range(4)],
        )
        s = sm.score_signature(expr, sm.SIGNATURES["Proliferation"])
        assert not (s == 0).all()
        assert s["c3"] > s["c0"]


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

    def test_combine_mean_raw_when_log1p_off(self):
        s = sm.score_signature(
            self._expr(), sm.SIGNATURES["Cytolytic"], combine="mean", log1p=False,
        )
        # raw mean of GZMB,PRF1: c3=(5+6)/2=5.5 highest
        assert s["c3"] == 5.5

    def test_default_applies_log1p(self):
        import numpy as np
        # combine="mean" with default log1p -> mean of log1p(TPM)
        s = sm.score_signature(self._expr(), sm.SIGNATURES["Cytolytic"], combine="mean")
        expected = (np.log1p(5.0) + np.log1p(6.0)) / 2  # c3
        assert abs(s["c3"] - expected) < 1e-9

    def test_combine_zscore_default(self):
        # default: z-score across cells of log1p(TPM), then mean across genes
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

    def test_gap_returns_none_when_no_distinct_cluster(self):
        # Smooth, evenly-spaced scores → no cluster "above the others" →
        # gap selects nothing rather than an arbitrary top slice.
        smooth = pd.Series(np.linspace(0.0, 1.0, 30), index=[f"c{i}" for i in range(30)])
        assert int(sm.call_positive(smooth, method="gap").sum()) == 0

    def test_otsu_isolates_cluster(self):
        assert int(sm.call_positive(self._scores(), method="otsu").sum()) == 4

    def test_quantile_fixed(self):
        assert int(sm.call_positive(self._scores(), method="quantile", quantile=0.75).sum()) == 5

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="unknown method"):
            sm.call_positive(self._scores(), method="bogus")


class TestBuildSignatureMethods:
    def _data(self):
        cells = [f"c{i}" for i in range(4)]
        expr = pd.DataFrame(
            {"ENTPD1": [5.0, 5.0, 0.0, 0.0], "CXCL13": [4.0, 6.0, 0.0, 1.0]},
            index=cells,
        )
        obs = pd.DataFrame({"CDR3ab": ["A", "A", "B", "B"], "sample": ["S1"] * 4},
                          index=cells)
        return expr, obs

    def test_scores_any_signature_on_any_data(self):
        # No context guard: a signature can be scored on any data.
        expr, obs = self._data()
        out = sm.build_signature_methods(
            expr, obs, signatures=["TumorReactive"], positive_method="gap",
        )
        assert set(out["method"]) <= {"TumorReactive"}

    def test_defaults_to_selection_signatures(self):
        expr, obs = self._data()
        out = sm.build_signature_methods(expr, obs, positive_method="gap")
        # default signatures = SELECTION_SIGNATURES (Tumor/AntigenExp)
        assert set(out["method"]) <= {"TumorReactive", "AntigenExperienced"}


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


class TestAimPanel:
    """AIM (activation-induced marker) GEX panel — the transcriptional
    analogue of an AIMpos sort (#pilot)."""

    def test_aim_focal_is_4_1bb_plus_ox40(self):
        aim = sm.SIGNATURES["AIM"]
        assert set(aim.genes_up) == {"TNFRSF9", "TNFRSF4"}  # 4-1BB + OX40
        assert aim.genes_down == ()
        assert aim.panel == "focal"

    def test_aim_broad_adds_surface_markers(self):
        broad = sm.SIGNATURES["AIMBroad"]
        assert {"TNFRSF9", "TNFRSF4", "CD69", "IL2RA", "CD40LG"} <= set(broad.genes_up)
        assert len(broad.genes_up) > len(sm.SIGNATURES["AIM"].genes_up)

    def test_aim_scores_activated_cells_higher(self):
        import pandas as pd
        cells = [f"c{i}" for i in range(4)]
        expr = pd.DataFrame(
            {"TNFRSF9": [0.0, 0.0, 5.0, 6.0], "TNFRSF4": [0.0, 1.0, 4.0, 7.0]},
            index=cells,
        )
        s = sm.score_signature(expr, sm.SIGNATURES["AIM"])
        assert s["c3"] > s["c0"]
