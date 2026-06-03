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

"""Tests for the TCR annotation API + PRISM selection (#158)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift import seqprob
from tcrsift.annotate_tcrs import (
    PRISM_DEFAULT_PREDICATES,
    add_gex_signature_scores,
    add_pgen_ppost,
    naive_signature,
    prism_score,
    score_gex_signature_per_clone,
    select_prism,
)


class TestPrismScore:
    def _df(self):
        # Four clones; clone A is the ideal PRISM pick: low ppost (private),
        # high antigen response, low naive.
        return pd.DataFrame({
            "clone": list("ABCD"),
            "ppost_alpha": [-20.0, -5.0, -6.0, -4.0],
            "ppost_beta": [-22.0, -5.0, -7.0, -3.0],
            "antigen_response_score": [2.0, 0.0, 1.0, -1.0],
            "naive_score": [-2.0, 1.0, 0.0, 2.0],
        })

    def test_ideal_clone_ranks_first(self):
        out = prism_score(self._df())
        assert {"prism_score", "prism_rank"} <= set(out.columns)
        best = out.loc[out.prism_rank == 1, "clone"].iloc[0]
        assert best == "A"
        # PRISM score is a mean of percentile ranks in [0,1]
        assert (out.prism_score.dropna().between(0, 1)).all()

    def test_missing_predicate_column_raises(self):
        df = self._df().drop(columns=["naive_score"])
        with pytest.raises(ValueError, match="missing predicate columns"):
            prism_score(df)

    def test_weights_change_ranking(self):
        df = self._df()
        # Weight antigen-response and naive to zero → rank on ppost only.
        out = prism_score(df, predicates=PRISM_DEFAULT_PREDICATES,
                          weights=[1.0, 1.0, 0.0, 0.0])
        order = out.sort_values("prism_score")["clone"].tolist()
        assert order[0] == "A"
        assert order.index("C") < order.index("D")

    def test_within_group_ranking(self):
        df = self._df()
        df["cond"] = ["X", "X", "Y", "Y"]
        out = prism_score(df, group_col="cond")
        for cond in ("X", "Y"):
            sub = out[out.cond == cond]
            assert sub.prism_score.min() == 0.0  # best-in-group hits rank 0

    def test_select_prism_top_k(self):
        out = select_prism(self._df(), k=2)
        assert out.prism_selected.sum() == 2
        assert out.loc[out.clone == "A", "prism_selected"].iloc[0]


class TestGexSignaturePerClone:
    def test_group_zscore_then_per_clone_mean(self):
        # Two samples with different baselines; within-sample z-scoring keeps
        # the high-baseline sample from dominating. Clone A is antigen-high.
        per_cell = pd.DataFrame({
            "CDR3ab": ["A", "A", "B", "B", "A", "B"],
            "sample": ["s1", "s1", "s1", "s1", "s2", "s2"],
            "gex.TNFRSF9": [5.0, 6.0, 0.0, 1.0, 10.0, 5.0],
            "gex.MKI67":   [4.0, 5.0, 0.0, 1.0, 9.0, 4.0],
        })
        per_clone = score_gex_signature_per_clone(
            per_cell, "AcuteActivation", clone_col="CDR3ab",
            group_col="sample", gex_prefix="gex",
        )
        assert per_clone["A"] > per_clone["B"]

    def test_missing_genes_raises(self):
        per_cell = pd.DataFrame({"CDR3ab": ["A"], "gex.FOO": [1.0]})
        with pytest.raises(ValueError, match="none of"):
            score_gex_signature_per_clone(per_cell, "AcuteActivation")

    def test_naive_signature_genes(self):
        assert set(naive_signature().all_genes) == {
            "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28"
        }

    def test_add_gex_scores_joins_on_clone(self):
        per_cell = pd.DataFrame({
            "CDR3ab": ["A", "A", "B"],
            "sample": ["s1", "s1", "s1"],
            "gex.TNFRSF9": [5.0, 6.0, 0.0],
            "gex.MKI67": [4.0, 5.0, 0.0],
            "gex.TCF7": [0.0, 1.0, 5.0], "gex.LEF1": [0.0, 0.0, 5.0],
            "gex.CCR7": [0.0, 1.0, 6.0], "gex.SELL": [1.0, 0.0, 5.0],
            "gex.IL7R": [0.0, 1.0, 5.0], "gex.CD27": [1.0, 0.0, 4.0],
            "gex.CD28": [0.0, 1.0, 5.0],
        })
        df = pd.DataFrame({"CDR3ab": ["A", "B"]})
        out = add_gex_signature_scores(
            df, per_cell, clone_col="CDR3ab", group_col="sample",
        )
        assert {"antigen_response_score", "naive_score"} <= set(out.columns)
        assert out.loc[out.CDR3ab == "A", "antigen_response_score"].iloc[0] > \
            out.loc[out.CDR3ab == "B", "antigen_response_score"].iloc[0]


def _models_built() -> bool:
    try:
        seqprob.load_background_model("beta", "kmer", "pgen")
        return True
    except FileNotFoundError:
        return False


class TestAddPgenPpost:
    def test_adds_pgen_ppost_columns(self):
        if not _models_built():
            pytest.skip("shipped k-mer models not built yet")
        df = pd.DataFrame({
            "CDR3_alpha": ["CAVRDSNYQLIW", "CAGHTGNQFYF"],
            "CDR3_beta": ["CASSLGTGELFF", "CASSPGQGAYEQYF"],
        })
        out = add_pgen_ppost(df, backend="kmer")
        # Pgen + Ppost (both chains) ship from neutral references (α Ppost now
        # fit on the pooled 10x healthy-PBMC reference) → all finite.
        for col in ("pgen_alpha", "pgen_beta", "ppost_alpha", "ppost_beta",
                    "log_q_alpha", "log_q_beta"):
            assert col in out.columns
            assert np.isfinite(out[col]).all()
        # Ppost is a distinct observed-repertoire model, NOT a Pgen masquerade.
        assert not np.allclose(out["ppost_alpha"], out["pgen_alpha"])

    def test_skips_absent_chain(self):
        if not _models_built():
            pytest.skip("shipped k-mer models not built yet")
        df = pd.DataFrame({"CDR3_beta": ["CASSLGTGELFF"]})
        out = add_pgen_ppost(df, chains=("alpha", "beta"))
        assert "pgen_beta" in out.columns
        assert "pgen_alpha" not in out.columns  # no CDR3_alpha column


class TestSeqprobRolesAndQ:
    def test_score_log_q_is_ppost_minus_pgen(self):
        if not _models_built():
            pytest.skip("shipped k-mer models not built yet")
        ppost = seqprob.load_background_model("beta", "kmer", "ppost")
        pgen = seqprob.load_background_model("beta", "kmer", "pgen")
        df = pd.DataFrame({"CDR3_beta": ["CASSLGTGELFF"]})
        log_q = seqprob.score_log_q(df, chain="beta")
        expected = (ppost.log_prob(["CASSLGTGELFF"])[0]
                    - pgen.log_prob(["CASSLGTGELFF"])[0])
        assert abs(log_q.iloc[0] - expected) < 1e-6

    def test_role_validation(self):
        with pytest.raises(ValueError, match="role"):
            seqprob._default_model_path("beta", "kmer", "bogus")


class TestPairedPpost:
    def test_sum_of_logs_and_either(self):
        from tcrsift.annotate_tcrs import add_paired_ppost
        df = pd.DataFrame({
            "ppost_alpha": [-10.0, -3.0, np.nan],
            "ppost_beta": [-12.0, -5.0, -4.0],
        })
        out = add_paired_ppost(df)
        # paired = sum of logs (joint, independent recombination)
        assert out.loc[0, "ppost_paired"] == -22.0
        assert np.isnan(out.loc[2, "ppost_paired"])  # missing chain → NaN
        # either = max (public if EITHER chain common) — don't hide public α
        assert out.loc[0, "ppost_either"] == -10.0
        assert out.loc[1, "ppost_either"] == -3.0
