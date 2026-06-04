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

"""Tests for the specificity/output-polish batch: #165, #169, #148, #146, #160."""

from __future__ import annotations

import logging

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# --------------------------------------------------------------------------- #
# #165 — high multi-chain / doublet-rate warning
# --------------------------------------------------------------------------- #
class TestDoubletRateWarning:
    def _adata(self, multi_flags, sample="S1"):
        n = len(multi_flags)
        obs = pd.DataFrame({
            "CDR3_alpha": ["CAAA"] * n,
            "CDR3_beta": ["CBBB"] * n,
            "TRA_1_umis": [5] * n,
            "TRB_1_umis": [5] * n,
            "multi_chain": multi_flags,
            "sample": [sample] * n,
        }, index=[f"c{i}" for i in range(n)])
        return ad.AnnData(X=csr_matrix((n, 1)), obs=obs)

    def test_warns_when_rate_high(self, caplog):
        from tcrsift.clonotype import aggregate_clonotypes

        # 3/10 multi-chain = 30% >= 10% default threshold.
        flags = [True, True, True] + [False] * 7
        with caplog.at_level(logging.WARNING, logger="tcrsift.clonotype"):
            aggregate_clonotypes(
                self._adata(flags), verbose=False, show_progress=False
            )
        warns = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("High multi-chain rate" in m for m in warns)

    def test_silent_when_rate_low(self, caplog):
        from tcrsift.clonotype import aggregate_clonotypes

        flags = [True] + [False] * 99  # 1% < 10%
        with caplog.at_level(logging.WARNING, logger="tcrsift.clonotype"):
            aggregate_clonotypes(
                self._adata(flags), verbose=False, show_progress=False
            )
        warns = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert not any("multi-chain" in m for m in warns)

    def test_per_sample_breakdown(self, caplog):
        from tcrsift.clonotype import aggregate_clonotypes

        # Sample HOT is 50% multi-chain; sample OK is 0%.
        n = 12
        obs = pd.DataFrame({
            "CDR3_alpha": ["CAAA"] * n,
            "CDR3_beta": ["CBBB"] * n,
            "TRA_1_umis": [5] * n,
            "TRB_1_umis": [5] * n,
            "multi_chain": [True, True, True] + [False] * 3 + [False] * 6,
            "sample": ["HOT"] * 6 + ["OK"] * 6,
        }, index=[f"c{i}" for i in range(n)])
        adata = ad.AnnData(X=csr_matrix((n, 1)), obs=obs)
        with caplog.at_level(logging.WARNING, logger="tcrsift.clonotype"):
            aggregate_clonotypes(adata, verbose=False, show_progress=False)
        warns = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("sample HOT" in m for m in warns)
        assert not any("sample OK" in m for m in warns)

    def test_disabled_with_zero_threshold(self, caplog):
        from tcrsift.clonotype import aggregate_clonotypes

        flags = [True] * 5 + [False] * 5  # 50%
        with caplog.at_level(logging.WARNING, logger="tcrsift.clonotype"):
            aggregate_clonotypes(
                self._adata(flags), doublet_warn_rate=0, verbose=False,
                show_progress=False,
            )
        warns = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert not any("multi-chain" in m for m in warns)


# --------------------------------------------------------------------------- #
# #169 — vector plot output
# --------------------------------------------------------------------------- #
class TestPlotFormat:
    def _fig(self):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        return fig

    def test_default_png_only(self, tmp_path):
        from tcrsift.plots import save_figure, set_plot_format

        set_plot_format("png")
        try:
            save_figure(self._fig(), tmp_path / "x.png")
        finally:
            set_plot_format("png")
        assert (tmp_path / "x.png").exists()
        assert not (tmp_path / "x.pdf").exists()

    def test_pdf_emits_vector_plus_png(self, tmp_path):
        from tcrsift.plots import save_figure, set_plot_format

        set_plot_format("pdf")
        try:
            out = save_figure(self._fig(), tmp_path / "x.png")
        finally:
            set_plot_format("png")
        # Both written; primary return is the vector file. The PNG keeps the
        # raster-embedding PDF report working.
        assert (tmp_path / "x.pdf").exists()
        assert (tmp_path / "x.png").exists()
        assert out.suffix == ".pdf"

    def test_invalid_format_rejected(self):
        from tcrsift.plots import set_plot_format

        with pytest.raises(ValueError):
            set_plot_format("tiff")
        set_plot_format("png")


# --------------------------------------------------------------------------- #
# #148 — alpha-beta pairing promiscuity
# --------------------------------------------------------------------------- #
class TestPairingPromiscuity:
    def test_counts_distinct_beta_partners(self):
        from tcrsift.annotate_tcrs import add_pairing_promiscuity

        df = pd.DataFrame({
            # Public alpha A1 pairs with 3 distinct betas; A2 with 1.
            "CDR3_alpha": ["A1", "A1", "A1", "A2"],
            "CDR3_beta": ["B1", "B2", "B3", "B9"],
        })
        out = add_pairing_promiscuity(df, promiscuous_min_partners=3)
        prom = out.set_index("CDR3_beta")["alpha_beta_promiscuity"].to_dict()
        assert prom["B1"] == 3 and prom["B2"] == 3 and prom["B3"] == 3
        assert prom["B9"] == 1
        flag = out.set_index("CDR3_beta")["alpha_promiscuous"].to_dict()
        assert flag["B1"]  # promiscuous
        assert not flag["B9"]
        assert out["alpha_cdr3_length"].tolist() == [2, 2, 2, 2]

    def test_empty_beta_not_counted(self):
        from tcrsift.annotate_tcrs import add_pairing_promiscuity

        df = pd.DataFrame({
            "CDR3_alpha": ["A1", "A1", "A1"],
            "CDR3_beta": ["B1", "B2", ""],  # the empty β must not count
        })
        out = add_pairing_promiscuity(df)
        assert out["alpha_beta_promiscuity"].tolist() == [2, 2, 2]

    def test_noop_without_cdr3_columns(self):
        from tcrsift.annotate_tcrs import add_pairing_promiscuity

        df = pd.DataFrame({"CDR3ab": ["A_B"]})
        out = add_pairing_promiscuity(df)
        assert "alpha_beta_promiscuity" not in out.columns

    def test_wired_into_annotate_tcrs_orchestrator(self):
        from tcrsift.annotate_tcrs import annotate_tcrs

        df = pd.DataFrame({
            "CDR3ab": ["A1_B1", "A1_B2"],
            "CDR3_alpha": ["A1", "A1"],
            "CDR3_beta": ["B1", "B2"],
        })
        out = annotate_tcrs(df, methods=["promiscuity"])
        assert "alpha_beta_promiscuity" in out.columns
        assert out["alpha_beta_promiscuity"].tolist() == [2, 2]


# --------------------------------------------------------------------------- #
# #146 — private / specificity-candidate selection route
# --------------------------------------------------------------------------- #
class TestSpecificityRoute:
    def _clones(self):
        # 4 complete clones above the freq gate, varying ppost_alpha; 1 below gate.
        return pd.DataFrame({
            "CDR3_alpha": ["A1", "A2", "A3", "A4", "A5"],
            "CDR3_beta": ["B1", "B2", "B3", "B4", "B5"],
            "max_frequency_per_method": [0.05, 0.05, 0.05, 0.05, 0.0001],
            "ppost_alpha": [1e-9, 1e-8, 1e-6, 1e-5, 1e-12],
        })

    def test_bottom_quartile_low_ppost_flagged(self):
        from tcrsift.selection import select_specificity_candidates

        out = select_specificity_candidates(self._clones(), percentile=25)
        sel = out.set_index("CDR3_alpha")
        # Lowest-ppost gated clone (A1) is the candidate; the high-freq-failing
        # A5 (lower ppost but below the freq gate) is excluded.
        assert sel.loc["A1", "specificity_candidate"]
        assert not sel.loc["A5", "specificity_gated"]
        assert not sel.loc["A4", "specificity_candidate"]
        # Rank 1 = lowest ppost among gated.
        assert sel.loc["A1", "specificity_ppost_rank"] == 1

    def test_missing_ppost_raises(self):
        from tcrsift.selection import select_specificity_candidates

        with pytest.raises(ValueError):
            select_specificity_candidates(pd.DataFrame({"CDR3_alpha": ["A"]}))

    def test_falls_back_to_max_frequency(self):
        from tcrsift.selection import select_specificity_candidates

        # No max_frequency_per_method column — must gate on max_frequency
        # instead of silently excluding every clone.
        df = pd.DataFrame({
            "CDR3_alpha": ["A1", "A2"],
            "CDR3_beta": ["B1", "B2"],
            "max_frequency": [0.05, 0.05],
            "ppost_alpha": [1e-9, 1e-4],
        })
        out = select_specificity_candidates(df, percentile=50)
        assert out["specificity_gated"].all()
        assert out.set_index("CDR3_alpha").loc["A1", "specificity_candidate"]


# --------------------------------------------------------------------------- #
# #160 — backend auto-select + never-zero / canonicalization acceptance
# --------------------------------------------------------------------------- #
class TestSeqprobBackendAndGuarantees:
    def test_backend_selection_by_size(self):
        from tcrsift.seqprob import select_backend_for_reference_size as sb

        assert sb(1000) == "kmer"
        assert sb(29_999) == "kmer"
        assert sb(50_000) == "tcrpeg"

    def test_kmer_never_zero_for_standard_cdr3(self):
        from tcrsift.seqprob import KmerProbabilityModel

        m = KmerProbabilityModel(order=2).fit(["CASSLGQAYEQYF", "CASSIRSSYEQYF"])
        # A never-before-seen but standard-AA CDR3 still scores finite (#160 #1).
        lp = m.log_prob(["CAWWWWWWWWF"])
        assert np.isfinite(lp[0])

    def test_gene_marginal_unseen_gene_finite(self):
        from tcrsift.seqprob import GeneAwareKmerModel

        m = GeneAwareKmerModel(order=2).fit(
            ["CASSLGQAYEQYF", "CASSIRSSYEQYF"],
            v_genes=["TRBV20-1", "TRBV20-1"],
            j_genes=["TRBJ2-7", "TRBJ2-7"],
        )
        # An unseen V gene falls back to the Laplace tail, never "gene not found".
        lp = m.log_prob(["CASSLGQAYEQYF"], v_genes=["TRBV99-9"], j_genes=["TRBJ2-7"])
        assert np.isfinite(lp[0])

    def test_canonicalize_resolves_10x_and_olga_forms(self):
        from tcrsift.genes import canonicalize_gene

        # The DV slash/no-slash variants and Adaptive forms collapse (#160 #2/#3).
        assert canonicalize_gene("TRAV14/DV4") == canonicalize_gene("TRAV14DV4")
        assert canonicalize_gene("TRBV20-1*01") == "TRBV20-1"
        assert canonicalize_gene("TCRBV20-01") == "TRBV20-1"
