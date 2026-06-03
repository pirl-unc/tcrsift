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

"""Tests for the OLGA/SONIA-backed Pgen + Ppost features (#143).

The heavy model-loading tests are gated behind ``olga``+``sonia`` being
importable (the GPL-3.0 optional extra). The pure-Python helpers
(``normalize_gene_name``, ``flag_private_candidates``, the missing-deps
error path) are always exercised.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift import olga_ppost as op

_HAS_PGEN = op.olga_sonia_available()
_needs_pgen = pytest.mark.skipif(
    not _HAS_PGEN, reason="olga+sonia not installed (optional [pgen] extra)"
)


class TestPureHelpers:
    """No OLGA/SONIA needed."""

    def test_normalize_gene_name_appends_allele(self):
        assert op.normalize_gene_name("TRBV20-1") == "TRBV20-1*01"
        assert op.normalize_gene_name("TRAV12-2") == "TRAV12-2*01"

    def test_normalize_gene_name_keeps_existing_allele(self):
        assert op.normalize_gene_name("TRBV20-1*03") == "TRBV20-1*03"

    def test_normalize_gene_name_handles_empty(self):
        assert op.normalize_gene_name("") is None
        assert op.normalize_gene_name(None) is None
        assert op.normalize_gene_name(float("nan")) is None

    def test_flag_private_candidates_rare_and_expanded(self):
        df = pd.DataFrame({
            "log10_ppost": [-5.0, -12.0, -12.5, -4.0],
            "frequency": [0.5, 0.5, 0.0001, 0.0001],
        })
        flag = op.flag_private_candidates(
            df, freq_col="frequency", score_quantile=0.5, min_freq=0.01,
        )
        # row 1: rare (bottom half) AND expanded → True.
        # row 2: rare but not expanded → False. rows 0/3: not rare → False.
        assert list(flag) == [False, True, False, False]

    def test_flag_private_candidates_without_freq_is_rare_only(self):
        df = pd.DataFrame({"log10_ppost": [-5.0, -12.0, -12.5, -4.0]})
        flag = op.flag_private_candidates(df, score_quantile=0.5, freq_col=None)
        assert list(flag) == [False, True, True, False]

    def test_flag_private_candidates_all_nan(self):
        df = pd.DataFrame({"log10_ppost": [np.nan, np.nan]})
        flag = op.flag_private_candidates(df)
        assert not flag.any()

    def test_unknown_chain_raises(self):
        with pytest.raises(ValueError, match="alpha.*beta"):
            op.load_chain_model("gamma")


class TestNearestSupportedAllele:
    """#150 — map unsupported/novel V/J alleles to nearest supported.

    Pure logic: tested against an explicit ``supported`` set, no model load.
    """

    SUPPORTED = frozenset({
        "TRBV20-1*01", "TRBV20-1*02", "TRBV28*01", "TRBV19*01",
        "TRBV6-1*01", "TRBV6-5*01", "TRBJ2-2*01",
    })

    def test_exact_match_no_substitution(self):
        m, r = op._nearest_in_supported("TRBV20-1*01", self.SUPPORTED)
        assert (m, r) == ("TRBV20-1*01", "exact")

    def test_bare_gene_appends_allele_and_is_exact(self):
        # CellRanger name without allele → *01 → already supported.
        m, r = op._nearest_in_supported("TRBV28", self.SUPPORTED)
        assert (m, r) == ("TRBV28*01", "exact")

    def test_novel_allele_maps_to_nearest_allele_same_gene(self):
        # *99 of a supported gene → lowest supported allele of that gene.
        m, r = op._nearest_in_supported("TRBV20-1*99", self.SUPPORTED)
        assert (m, r) == ("TRBV20-1*01", "nearest_allele")

    def test_unsupported_gene_maps_to_nearest_gene_by_family(self):
        # TRBV6-3 unsupported → nearest TRBV6-* by subgroup distance.
        m, r = op._nearest_in_supported("TRBV6-3", self.SUPPORTED)
        assert r == "nearest_gene"
        assert m in {"TRBV6-1*01", "TRBV6-5*01"}
        # subgroup 3 is closer to 1 (dist 2) than 5 (dist 4)
        assert m == "TRBV6-1*01"

    def test_unsupported_family_maps_to_nearest_family(self):
        # TRBV21 unsupported, no TRBV21-* → nearest family number (20 or 28).
        m, r = op._nearest_in_supported("TRBV21-1", self.SUPPORTED)
        assert r == "nearest_gene"
        assert m == "TRBV20-1*01"  # |21-20|=1 < |21-28|=7

    def test_different_locus_is_unmapped(self):
        # No TRAV in a TRBV-only supported set → unmapped.
        m, r = op._nearest_in_supported("TRAV12-2", self.SUPPORTED)
        assert (m, r) == (None, "unmapped")

    def test_empty_gene_unmapped(self):
        assert op._nearest_in_supported(None, self.SUPPORTED) == (None, "unmapped")
        assert op._nearest_in_supported("", self.SUPPORTED) == (None, "unmapped")

    def test_annotate_column(self):
        df = pd.DataFrame({"v_gene": ["TRBV20-1*99", "TRBV28", "TRAV1-1"]})
        # Patch supported_alleles so no model load is needed.
        import tcrsift.olga_ppost as mod
        orig = mod.supported_alleles
        mod.supported_alleles = lambda chain, seg: self.SUPPORTED
        try:
            out = op.annotate_nearest_supported_allele(
                df, chain="beta", segment="v", gene_col="v_gene",
            )
        finally:
            mod.supported_alleles = orig
        assert list(out["nearest_supported_reason"]) == [
            "nearest_allele", "exact", "unmapped",
        ]
        assert out.loc[0, "nearest_supported_allele"] == "TRBV20-1*01"


@pytest.mark.skipif(
    _HAS_PGEN, reason="only meaningful when the extra is NOT installed"
)
def test_missing_deps_raises_with_hint():
    with pytest.raises(ImportError, match=r"pip install tcrsift\[pgen\]"):
        op.load_chain_model("beta")


@_needs_pgen
class TestWithModels:
    """Exercised only where olga+sonia are importable (e.g. local dev)."""

    def test_supported_alleles_nonempty(self):
        v = op.supported_alleles("beta", "v")
        j = op.supported_alleles("beta", "j")
        assert any(a.startswith("TRBV") for a in v)
        assert any(a.startswith("TRBJ") for a in j)
        assert "TRBV20-1*01" in v

    def test_alpha_model_is_vj_not_beta(self):
        # Regression: load_dir alone loads alpha weights as a beta (VDJ)
        # model; chain_type must be passed explicitly. Alpha V alleles must
        # be TRAV, not TRBV.
        v = op.supported_alleles("alpha", "v")
        assert any(a.startswith("TRAV") for a in v)
        assert not any(a.startswith("TRBV") for a in v)

    def test_compute_pgen_ppost_beta(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLGTGELFF", "", "NOTAREALSEQ123"],
            "beta_v_gene": ["TRBV20-1", "TRBV20-1", "TRBV20-1"],
            "beta_j_gene": ["TRBJ2-2", "TRBJ2-2", "TRBJ2-2"],
        })
        out = op.compute_pgen_ppost(df, chain="beta")
        assert {"log10_pgen_olga", "sonia_q", "log10_ppost"} <= set(out.columns)
        # First row is a real productive CDR3 → finite Pgen/Ppost.
        assert np.isfinite(out.loc[0, "log10_pgen_olga"])
        assert np.isfinite(out.loc[0, "log10_ppost"])
        # Empty CDR3 → NaN, never a silent default.
        assert np.isnan(out.loc[1, "log10_pgen_olga"])

    def test_ppost_equals_pgen_plus_logq(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLGTGELFF"],
            "beta_v_gene": ["TRBV20-1"],
            "beta_j_gene": ["TRBJ2-2"],
        })
        out = op.compute_pgen_ppost(df, chain="beta")
        expected = out.loc[0, "log10_pgen_olga"] + np.log10(out.loc[0, "sonia_q"])
        assert abs(out.loc[0, "log10_ppost"] - expected) < 1e-9

    def test_alpha_compute_runs(self):
        df = pd.DataFrame({
            "CDR3_alpha": ["CAVRDSNYQLIW"],
            "alpha_v_gene": ["TRAV12-2"],
            "alpha_j_gene": ["TRAJ33"],
        })
        out = op.compute_pgen_ppost(
            df, chain="alpha", cdr3_col="CDR3_alpha",
            v_gene_col="alpha_v_gene", j_gene_col="alpha_j_gene",
        )
        assert np.isfinite(out.loc[0, "log10_pgen_olga"])

    def test_real_supported_alleles_resolve_against_model(self):
        # A novel allele of a real gene resolves to a supported allele.
        m, r = op.nearest_supported_allele("TRBV20-1*99", "beta", "v")
        assert r in {"nearest_allele", "exact"}
        assert m in op.supported_alleles("beta", "v")

    def test_unsupported_allele_substituted_before_pgen(self):
        # With map_unsupported, a *99 novel allele is mapped and still scores;
        # the used allele is recorded and is finite Pgen (no silent NaN).
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLGTGELFF"],
            "beta_v_gene": ["TRBV20-1*99"],
            "beta_j_gene": ["TRBJ2-2"],
        })
        out = op.compute_pgen_ppost(
            df, chain="beta", map_unsupported=True,
            used_v_col="beta_v_used",
        )
        assert out.loc[0, "beta_v_used"] in op.supported_alleles("beta", "v")
        assert np.isfinite(out.loc[0, "log10_pgen_olga"])
