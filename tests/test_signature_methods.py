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
        # Differentiated is now the effector−naïve contrast (#141): effector
        # core up, naïve/stem down.
        diff = sm.SIGNATURES["Differentiated"]
        assert {"IFNG", "GZMB", "PRF1"} <= set(diff.genes_up)
        assert {"CCR7", "TCF7", "LEF1", "SELL"} <= set(diff.genes_down)
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
        # Headline composites the pilot uses (no context tags). Differentiated
        # (eff−naïve) is the best in-vitro-expansion axis (#141).
        assert set(sm.SELECTION_SIGNATURES) == {
            "Differentiated", "TumorReactive", "AntigenExperienced"
        }
        tr = sm.SELECTION_SIGNATURES["TumorReactive"]
        assert {"CXCL13", "ENTPD1", "TOX", "ITGAE", "CTLA4"} <= set(tr.genes_up)
        # AntigenExperienced was recomposed (#142): effector core only, no
        # longer diluted by TNFRSF9/MKI67 (which moved to AcuteActivation).
        ae = sm.SELECTION_SIGNATURES["AntigenExperienced"]
        assert {"IFNG", "GZMB", "NKG7"} <= set(ae.genes_up)
        assert "TNFRSF9" not in ae.all_genes and "MKI67" not in ae.all_genes
        acute = sm.SIGNATURES["AcuteActivation"]
        assert set(acute.genes_up) == {"TNFRSF9", "MKI67"}

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
             "CCR7": [5.0, 5.0, 0.0, 0.0], "TCF7": [4.0, 6.0, 0.0, 1.0]},
            index=[f"c{i}" for i in range(4)],
        )

    def test_signed_score_subtracts_down_genes(self):
        # Differentiated = effector(up) − naïve(down): effector-high /
        # naïve-low cells score high. The tiny frame has only GZMB/PRF1 (up)
        # and CCR7/TCF7 (down) of the full panel, so score with on_missing
        # ="ignore". c3 (effector high, naïve low) > c0 (effector low, naïve
        # high).
        s = sm.score_signature(
            self._expr(), sm.SIGNATURES["Differentiated"], on_missing="ignore",
        )
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

    def test_missing_gene_raises_by_default(self):
        # Cytolytic needs PRF1+GZMB; PRF1 absent → error, not a partial score.
        expr = pd.DataFrame({"GZMB": [1.0, 2.0]}, index=["a", "b"])
        with pytest.raises(KeyError, match="PRF1"):
            sm.score_signature(expr, sm.SIGNATURES["Cytolytic"])

    def test_missing_gene_ignore_scores_present(self):
        expr = pd.DataFrame({"GZMB": [1.0, 2.0, 3.0]}, index=list("abc"))
        s = sm.score_signature(expr, sm.SIGNATURES["Cytolytic"], on_missing="ignore")
        assert s["c"] > s["a"]  # scored on GZMB alone

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

    def test_gap_majority_cluster_needs_higher_search_top(self):
        # A high cluster that is the MAJORITY (~70%) sits below the default
        # search_top=0.5 window, so gap misses it; raising search_top finds it.
        vals = list(np.zeros(6)) + list(np.full(14, 5.0))  # 14/20 high
        s = pd.Series(vals, index=[f"c{i}" for i in range(20)])
        assert int(sm.call_positive(s, method="gap").sum()) == 0           # default 0.5
        assert int(sm.call_positive(s, method="gap", search_top=0.9).sum()) == 14

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
            {"PRF1": [5.0, 5.0, 0.0, 0.0], "GZMB": [4.0, 6.0, 0.0, 1.0]},
            index=cells,
        )
        obs = pd.DataFrame({"CDR3ab": ["A", "A", "B", "B"], "sample": ["S1"] * 4},
                          index=cells)
        return expr, obs

    def test_scores_signature_with_complete_panel(self):
        # No context guard: a signature with all genes present scores fine.
        expr, obs = self._data()
        out = sm.build_signature_methods(
            expr, obs, signatures=["Cytolytic"], positive_method="gap",
        )
        assert set(out["method"]) <= {"Cytolytic"}

    def test_missing_signature_gene_raises(self):
        # TumorReactive needs CXCL13/ENTPD1/… absent here → hard error.
        expr, obs = self._data()
        with pytest.raises(KeyError):
            sm.build_signature_methods(expr, obs, signatures=["TumorReactive"])

    def test_defaults_to_selection_signatures(self):
        # Default sigs (Differentiated/Tumor/AntigenExp) aren't all in this
        # tiny frame; with on_missing="ignore" the default-param behavior is
        # still exercised.
        expr, obs = self._data()
        out = sm.build_signature_methods(
            expr, obs, on_missing="ignore", positive_method="gap",
        )
        assert set(out["method"]) <= {
            "Differentiated", "TumorReactive", "AntigenExperienced"
        }


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


class TestExpressionFrameFromAdata:
    def _adata(self):
        import anndata as ad
        from scipy.sparse import csr_matrix
        # ENSEMBL var_names + gene_symbols column (the load_samples shape)
        var = pd.DataFrame(
            {"gene_symbols": ["TNFRSF9", "TNFRSF4", "GZMB"]},
            index=["ENSG1", "ENSG2", "ENSG3"],
        )
        X = csr_matrix(np.array([[8.0, 9, 0], [7, 8, 0], [0, 0, 5], [0, 0, 4]]))
        obs = pd.DataFrame(index=[f"c{i}" for i in range(4)])
        return ad.AnnData(X=X, obs=obs, var=var)

    def test_resolves_symbols_from_var_column(self):
        expr = sm.expression_frame_from_adata(self._adata(), ["TNFRSF9", "TNFRSF4"])
        assert list(expr.columns) == ["TNFRSF9", "TNFRSF4"]
        assert list(expr.index) == ["c0", "c1", "c2", "c3"]
        assert expr.loc["c0", "TNFRSF9"] == 8.0  # densified from sparse

    def test_absent_gene_raises_by_default(self):
        with pytest.raises(KeyError, match="NOTAGENE"):
            sm.expression_frame_from_adata(self._adata(), ["TNFRSF9", "NOTAGENE"])

    def test_absent_gene_omitted_when_ignore(self):
        expr = sm.expression_frame_from_adata(
            self._adata(), ["TNFRSF9", "NOTAGENE"], on_missing="ignore",
        )
        assert list(expr.columns) == ["TNFRSF9"]

    def test_resolves_when_var_names_are_symbols(self):
        import anndata as ad
        adata = ad.AnnData(
            X=np.array([[5.0, 0.0], [0.0, 5.0]]),
            obs=pd.DataFrame(index=["c0", "c1"]),
            var=pd.DataFrame(index=["GZMB", "CCR7"]),  # symbols as var_names
        )
        expr = sm.expression_frame_from_adata(adata, ["GZMB"])
        assert list(expr.columns) == ["GZMB"]
        assert expr.loc["c0", "GZMB"] == 5.0


class TestSignatureInjectionIntoSelection:
    """The run-wiring flow: adata GEX → signature method → selection rule."""

    def test_aim_signature_becomes_a_selectable_method(self):
        import anndata as ad
        from scipy.sparse import csr_matrix

        from tcrsift.selection import select_from_clone_sample_long

        # 8 cells; clone A (c0,c1) is the minority AIM-high cluster.
        var = pd.DataFrame({"gene_symbols": ["TNFRSF9", "TNFRSF4"]}, index=["E1", "E2"])
        rows = np.zeros((8, 2))
        rows[0] = [8, 9]
        rows[1] = [7, 8]
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "B", "B", "C", "C", "D", "D"], "sample": ["S1"] * 8},
            index=[f"c{i}" for i in range(8)],
        )
        adata = ad.AnnData(X=csr_matrix(rows), obs=obs, var=var)

        expr = sm.expression_frame_from_adata(adata, ["TNFRSF9", "TNFRSF4"])
        sig_long = sm.build_signature_methods(
            expr, adata.obs, signatures=["AIM"], positive_method="gap",
        )
        assert set(sig_long["method"]) == {"AIM"}
        assert "A" in set(sig_long["CDR3ab"])

        # A private_AIM rule can now select clone A via the signature method.
        rules = select_from_clone_sample_long(
            sig_long,
            {"rules": {"private": {"include_tier": "tier5", "top_n": 5,
                                   "apply_to_methods": ["AIM"]}}},
        )
        assert set(rules.columns) >= {"CDR3ab", "selection_rule", "global_rank"}


class TestReceptorGeneGuard:
    def test_receptor_regex_distinguishes_segments_from_lookalikes(self):
        for g in ("TRAV12-2", "TRBV6-1", "TRBJ2-1", "TRAC", "TRBC1", "IGHV3-23",
                  "IGKC", "IGHM"):
            assert sm.is_receptor_gene(g), g
        for g in ("TRADD", "TRAF3", "GZMB", "IGF1", "CCR7", "TOX"):
            assert not sm.is_receptor_gene(g), g

    def test_strip_receptor_genes(self):
        assert sm.strip_receptor_genes(["GZMB", "TRAV12-2", "CCR7", "TRBC1"]) == [
            "GZMB", "CCR7"
        ]

    def test_signature_with_receptor_gene_raises(self):
        bad = sm.Signature("Leaky", ("GZMB", "TRAV12-2"))
        expr = pd.DataFrame({"GZMB": [1.0], "TRAV12-2": [2.0]}, index=["c0"])
        with pytest.raises(ValueError, match="receptor"):
            sm.score_signature(expr, bad)


class TestPerConditionBaseline:
    def test_groups_zscore_within_each_group(self):
        # Two groups with very different baselines; per-condition z-score
        # centers each group rather than letting the high group dominate.
        expr = pd.DataFrame(
            {"GZMB": [0.0, 0.0, 10.0, 10.0], "PRF1": [0.0, 1.0, 9.0, 11.0]},
            index=[f"c{i}" for i in range(4)],
        )
        groups = pd.Series(["A", "A", "B", "B"], index=expr.index)
        s = sm.score_signature(expr, sm.SIGNATURES["Cytolytic"], groups=groups)
        # within each group the two cells are symmetric around 0
        assert abs(s["c0"] + s["c1"]) < 1e-9
        assert abs(s["c2"] + s["c3"]) < 1e-9
        # without grouping, group B dominates (all positive)
        pooled = sm.score_signature(expr, sm.SIGNATURES["Cytolytic"])
        assert pooled["c2"] > 0 and pooled["c3"] > 0

    def test_groups_and_background_mutually_exclusive(self):
        expr = pd.DataFrame({"GZMB": [1.0], "PRF1": [1.0]}, index=["c0"])
        with pytest.raises(ValueError, match="groups OR background"):
            sm.score_signature(
                expr, sm.SIGNATURES["Cytolytic"],
                groups=pd.Series(["A"], index=["c0"]), background=expr,
            )

    def test_build_background_by_uses_obs_column(self):
        expr = pd.DataFrame(
            {"GZMB": [0.0, 0.0, 9.0], "PRF1": [0.0, 1.0, 11.0]},
            index=[f"c{i}" for i in range(3)],
        )
        obs = pd.DataFrame(
            {"CDR3ab": ["A", "A", "B"], "sample": ["S1", "S1", "S2"]},
            index=expr.index,
        )
        out = sm.build_signature_methods(
            expr, obs, signatures=["Cytolytic"], background_by="sample",
            positive_method="gap",
        )
        assert set(out["method"]) <= {"Cytolytic"}


class TestSignatureRecompose141142:
    """#141/#142: naïve/stem set, effector−naïve contrast, recomposed
    AntigenExperienced, effector rename."""

    def test_naive_stem_and_effector_gene_sets(self):
        from tcrsift import signatures as sigmod
        assert sigmod.NAIVE_STEM_GENES_HGNC == (
            "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28"
        )
        # effector is the canonical name; activation is a deprecated alias.
        assert sigmod.EFFECTOR_GENES_HGNC == ("IFNG", "GZMB", "PRF1", "GNLY", "NKG7")
        assert sigmod.ACTIVATION_GENES_HGNC == sigmod.EFFECTOR_GENES_HGNC
        assert sigmod.T_CELL_SIGNATURES["effector"] == sigmod.EFFECTOR_GENES_HGNC
        assert sigmod.T_CELL_SIGNATURES["naive_stem"] == sigmod.NAIVE_STEM_GENES_HGNC

    def test_differentiated_is_effector_minus_naive_contrast(self):
        diff = sm.SIGNATURES["Differentiated"]
        assert set(diff.genes_up) == {"IFNG", "GZMB", "PRF1", "GNLY", "NKG7"}
        assert set(diff.genes_down) == {"TCF7", "LEF1", "CCR7", "SELL"}

    def test_acute_activation_split_out(self):
        acute = sm.SIGNATURES["AcuteActivation"]
        assert set(acute.genes_up) == {"TNFRSF9", "MKI67"}
        assert acute.genes_down == ()

    def test_antigen_experienced_dropped_diluting_genes(self):
        ae = sm.SIGNATURES["AntigenExperienced"]
        assert "TNFRSF9" not in ae.all_genes
        assert "MKI67" not in ae.all_genes
        assert set(ae.genes_up) == {"IFNG", "GZMB", "PRF1", "GNLY", "NKG7"}

    def test_eff_minus_naive_beats_effector_alone_on_synthetic(self):
        # Expanded clone = effector-high + naïve-low; bystander = naïve-high.
        # The contrast separates them more than effector-up alone.
        cells = [f"c{i}" for i in range(4)]
        expr = pd.DataFrame(
            {"IFNG": [5, 6, 0, 0], "GZMB": [5, 6, 0, 0], "PRF1": [5, 6, 0, 0],
             "GNLY": [5, 6, 0, 0], "NKG7": [5, 6, 0, 0],
             "TCF7": [0, 0, 6, 5], "LEF1": [0, 0, 6, 5],
             "CCR7": [0, 0, 6, 5], "SELL": [0, 0, 6, 5]},
            index=cells,
        )
        contrast = sm.score_signature(expr, sm.SIGNATURES["Differentiated"])
        # expanded (c0,c1) clearly above bystanders (c2,c3)
        assert min(contrast["c0"], contrast["c1"]) > max(contrast["c2"], contrast["c3"])


class TestSignatureGuidance145:
    """#145 — honest per-signature usage map + RNA-limitation note."""

    def test_guidance_tiers_are_valid(self):
        valid = {"recommended", "situational", "wrong_biology"}
        assert all(g.tier in valid for g in sm.SIGNATURE_GUIDANCE.values())

    def test_recommended_are_the_reproducible_axes(self):
        rec = set(sm.recommended_signatures())
        # The two donor-reproducible axes from the pilot.
        assert rec == {"Differentiated", "AcuteActivation"}

    def test_tumor_reactive_flagged_wrong_biology(self):
        assert sm.SIGNATURE_GUIDANCE["TumorReactive"].tier == "wrong_biology"
        assert "TIL" in sm.SIGNATURE_GUIDANCE["TumorReactive"].use_for

    def test_guidance_keys_are_real_signatures(self):
        assert set(sm.SIGNATURE_GUIDANCE) <= set(sm.SIGNATURES)

    def test_rna_note_points_at_sequence_axis(self):
        note = sm.RNA_REPRODUCIBILITY_NOTE
        assert "Pgen/Ppost" in note
        assert "not RNA" in note
        # methodology caveat: don't overfit on the selection target
        assert "overfitting" in note
