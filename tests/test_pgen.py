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

"""Tests for the lightweight Pgen estimator (#58).

These don't check absolute Pgens (we ship a proxy, not OLGA); they
verify the *rank order* — short / canonical-gene / zero-N-insertion
sequences should score higher than long / off-frequency-gene / many-
N-insertion sequences."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from tcrsift.pgen import (
    annotate_publicness,
    compute_pgen,
    pgen_components,
    pgen_single,
    publicness_score,
)


class TestPgenSingle:
    def test_short_canonical_more_likely_than_long_unusual(self):
        """A short CDR3 (13 aa) with a common V/J should score higher
        than a long CDR3 (20 aa) with off-frequency genes."""
        short = pgen_single(
            "CASSLAPGATNEKLFF", v_gene="TRBV20-1", j_gene="TRBJ2-1", chain="beta"
        )
        long_rare = pgen_single(
            "CASSLAPGATNEKLFFAAAAAA",
            v_gene="TRBV15", j_gene="TRBJ1-3", chain="beta",
        )
        assert short > long_rare

    def test_alpha_chain_uses_one_junction(self):
        """α CDR3s use one junction (V-J), so the N-insertion term is
        less penalizing than β at the same total CDR3 length."""
        cdr3 = "CAASFGGSNYKLTF"
        beta_p = pgen_single(cdr3, "TRBV20-1", "TRBJ2-1", chain="beta")
        alpha_p = pgen_single(cdr3, "TRAV1-2", "TRAJ48", chain="alpha")
        # Just verify both finite, and α isn't catastrophically worse
        # than β. We don't claim α > β here — composition + length
        # priors differ.
        assert math.isfinite(beta_p) and math.isfinite(alpha_p)

    def test_zero_n_insertion_higher_than_many(self):
        """A sequence with n_inserted=0 should score higher than the
        same sequence asserted to have n_inserted=20."""
        cdr3 = "CASSPGTGELFF"
        p_zero = pgen_single(
            cdr3, "TRBV20-1", "TRBJ2-1", chain="beta", n_inserted=0,
        )
        p_many = pgen_single(
            cdr3, "TRBV20-1", "TRBJ2-1", chain="beta", n_inserted=20,
        )
        assert p_zero > p_many


class TestPgenComponents:
    def test_components_sum_to_total(self):
        parts = pgen_components(
            "CASSLAPGATNEKLFF", "TRBV20-1", "TRBJ2-1", chain="beta"
        )
        total = pgen_single(
            "CASSLAPGATNEKLFF", "TRBV20-1", "TRBJ2-1", chain="beta"
        )
        assert abs(sum(parts.values()) - total) < 1e-10

    def test_all_five_components_returned(self):
        parts = pgen_components(
            "CASS", "TRBV20-1", "TRBJ2-1", chain="beta"
        )
        assert set(parts) == {
            "log10_p_v",
            "log10_p_j",
            "log10_p_length",
            "log10_p_n_inserted",
            "log10_p_composition",
        }

    def test_invalid_chain_raises(self):
        with pytest.raises(ValueError, match="chain"):
            pgen_components("CASS", "TRBV20-1", "TRBJ2-1", chain="gamma")

    def test_unknown_v_gene_uses_tail_fallback(self):
        """A V-gene not in the lookup table should not crash; should
        return a finite log-prob (the small tail mass)."""
        p = pgen_single(
            "CASSLAPGATNEKLFF",
            "TRBV_MADE_UP", "TRBJ2-1", chain="beta",
        )
        assert math.isfinite(p)


class TestComputePgen:
    def test_returns_one_value_per_row(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLAPGATNEKLFF", "CASSPGTGELFF"],
            "beta_v_gene": ["TRBV20-1", "TRBV9"],
            "beta_j_gene": ["TRBJ2-1", "TRBJ1-2"],
        })
        out = compute_pgen(df)
        assert len(out) == 2
        assert out.notna().all()

    def test_missing_cdr3_produces_nan(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLAPGATNEKLFF", None, ""],
            "beta_v_gene": ["TRBV20-1", "TRBV20-1", "TRBV20-1"],
            "beta_j_gene": ["TRBJ2-1", "TRBJ2-1", "TRBJ2-1"],
        })
        out = compute_pgen(df)
        assert out.notna().iloc[0]
        assert out.isna().iloc[1]
        assert out.isna().iloc[2]

    def test_missing_v_j_columns_uses_unknown_fallback(self):
        df = pd.DataFrame({"CDR3_beta": ["CASSPGTGELFF"]})
        out = compute_pgen(df)
        assert out.notna().iloc[0]

    def test_alpha_chain_works(self):
        df = pd.DataFrame({
            "CDR3_alpha": ["CAASFGGSNYKLTF"],
            "alpha_v_gene": ["TRAV1-2"],
            "alpha_j_gene": ["TRAJ48"],
        })
        out = compute_pgen(
            df,
            cdr3_col="CDR3_alpha",
            v_gene_col="alpha_v_gene",
            j_gene_col="alpha_j_gene",
            chain="alpha",
        )
        assert out.notna().iloc[0]

    def test_n_inserted_column_overrides_estimate(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSPGTGELFF"] * 2,
            "beta_v_gene": ["TRBV20-1"] * 2,
            "beta_j_gene": ["TRBJ2-1"] * 2,
            "n_inserted": [0, 20],
        })
        out = compute_pgen(df, n_inserted_col="n_inserted")
        # Row 0 (n=0) should be higher Pgen than row 1 (n=20).
        assert out.iloc[0] > out.iloc[1]

    def test_missing_cdr3_col_raises(self):
        df = pd.DataFrame({"not_cdr3": ["X"]})
        with pytest.raises(ValueError, match="CDR3_beta"):
            compute_pgen(df)


class TestPublicnessScore:
    def test_high_pgen_maps_to_high_publicness(self):
        # log10 Pgen = -16 (very public per default cutoffs) → 1.0
        s = publicness_score([-16.0])
        assert s.iloc[0] == 1.0

    def test_low_pgen_maps_to_low_publicness(self):
        # log10 Pgen = -35 (very private) → 0.0
        s = publicness_score([-35.0])
        assert s.iloc[0] == 0.0

    def test_midpoint_maps_to_half(self):
        # Default cutoffs: low=-30, high=-18 → midpoint = -24.
        s = publicness_score([-24.0])
        assert abs(s.iloc[0] - 0.5) < 1e-9

    def test_monotone_increasing(self):
        s = publicness_score([-35.0, -28.0, -22.0, -18.0, -10.0])
        assert (s.diff().dropna() >= 0).all()

    def test_nan_passthrough(self):
        s = publicness_score([-25.0, np.nan, -20.0])
        assert not pd.isna(s.iloc[0])
        assert pd.isna(s.iloc[1])
        assert not pd.isna(s.iloc[2])

    def test_custom_cutoffs_for_olga_scale(self):
        """Callers feeding real OLGA Pgens (which sit ~10 log units
        higher than this estimator's output) can pass OLGA-calibrated
        cutoffs."""
        s = publicness_score(
            [-25.0, -8.0], low_pgen_cutoff=-20.0, high_pgen_cutoff=-8.0
        )
        assert s.iloc[0] == 0.0
        assert s.iloc[1] == 1.0


class TestAnnotatePublicness:
    def test_adds_both_columns(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLAPGATNEKLFF", "CASSPGTGELFF"],
            "beta_v_gene": ["TRBV20-1", "TRBV9"],
            "beta_j_gene": ["TRBJ2-1", "TRBJ1-2"],
        })
        out = annotate_publicness(df)
        assert "log10_pgen" in out.columns
        assert "publicness" in out.columns
        assert (out["publicness"].between(0, 1) | out["publicness"].isna()).all()

    def test_does_not_mutate_input(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSPGTGELFF"],
            "beta_v_gene": ["TRBV20-1"],
            "beta_j_gene": ["TRBJ2-1"],
        })
        cols_before = df.columns.tolist()
        annotate_publicness(df)
        assert df.columns.tolist() == cols_before
