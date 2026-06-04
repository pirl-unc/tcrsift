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

"""Regression tests for #80: ``build_clone_sample_long`` must use the
same UMI-gated ``is_complete_clone`` rule as ``aggregate_clonotypes``,
and must drop NaN-valued CDR3ab rows from the denominator."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from tcrsift.clonotype import build_clone_sample_long


def _adata(rows):
    """Build a minimal AnnData from a list of obs dicts."""
    obs = pd.DataFrame(rows)
    X = csr_matrix((len(obs), 1))
    return ad.AnnData(X=X, obs=obs)


class TestBuildCloneSampleLongDenominator:
    def test_nan_cdr3ab_excluded_from_denominator(self):
        """NaN CDR3ab rows used to leak through ``!= ""`` because in
        pandas ``NaN != ""`` is True. Verify they're now filtered."""
        adata = _adata([
            {"sample": "S1", "CDR3ab": "CASS_A", "CDR3_alpha": "CAS", "CDR3_beta": "CAS",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "CASS_A", "CDR3_alpha": "CAS", "CDR3_beta": "CAS",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            # Chainless cell — CDR3ab NaN. Should NOT count toward
            # the sample-S1 denominator.
            {"sample": "S1", "CDR3ab": np.nan, "CDR3_alpha": np.nan, "CDR3_beta": np.nan,
             "TRA_1_umis": 0, "TRB_1_umis": 0},
        ])
        out = build_clone_sample_long(adata)
        # Only one clone, two cells in S1. Denominator = 2 (not 3).
        row = out[(out["CDR3ab"] == "CASS_A") & (out["sample"] == "S1")].iloc[0]
        assert row["cells"] == 2
        assert abs(row["frequency"] - 1.0) < 1e-9

    def test_umi_gated_numerator_and_denominator(self):
        """Cells with TRA_1_umis < 2 or TRB_1_umis < 2 are excluded from
        BOTH the numerator and the denominator — same complete-clone rule
        as ``aggregate_clonotypes``. Pre-#175 the numerator counted the
        UMI-failing cell (cells=3, freq=1.5), inflating per-sample
        frequencies above 1.0; now both sides use the gated count."""
        adata = _adata([
            # Two passing cells.
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            # One cell failing the UMI gate (TRA < 2): excluded from both
            # the numerator and the complete-clone denominator (#175).
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 0, "TRB_1_umis": 5},
        ])
        out = build_clone_sample_long(adata)
        row = out[(out["CDR3ab"] == "A_A") & (out["sample"] == "S1")].iloc[0]
        assert row["cells"] == 2
        assert abs(row["frequency"] - 1.0) < 1e-9


class TestBuildCloneSampleLongFrequencyInvariant:
    """#175: numerator and denominator must use the same complete-cell set
    so per-sample frequencies sum to 1.0 and single-chain clones never leak
    into the table."""

    def test_per_sample_frequencies_sum_to_one(self):
        adata = _adata([
            # Two complete clones in S1 (3 + 1 cells).
            {"sample": "S1", "CDR3ab": "A_B", "CDR3_alpha": "A", "CDR3_beta": "B",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "A_B", "CDR3_alpha": "A", "CDR3_beta": "B",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "A_B", "CDR3_alpha": "A", "CDR3_beta": "B",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "C_D", "CDR3_alpha": "C", "CDR3_beta": "D",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            # A single-chain (alpha-dropout) cell sharing beta B, and a
            # both-chain cell failing the UMI gate. Pre-#175 both inflated
            # the numerator without entering the denominator.
            {"sample": "S1", "CDR3ab": "_B", "CDR3_alpha": "", "CDR3_beta": "B",
             "TRA_1_umis": 0, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "E_F", "CDR3_alpha": "E", "CDR3_beta": "F",
             "TRA_1_umis": 1, "TRB_1_umis": 5},
        ])
        out = build_clone_sample_long(adata)
        freq_sum = out[out["sample"] == "S1"]["frequency"].sum()
        assert abs(freq_sum - 1.0) < 1e-9

    def test_single_chain_clones_absent(self):
        adata = _adata([
            {"sample": "S1", "CDR3ab": "A_B", "CDR3_alpha": "A", "CDR3_beta": "B",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            # Alpha-dropout single-chain "clone" — must not appear in output.
            {"sample": "S1", "CDR3ab": "_B", "CDR3_alpha": "", "CDR3_beta": "B",
             "TRA_1_umis": 0, "TRB_1_umis": 5},
        ])
        out = build_clone_sample_long(adata)
        assert (out["CDR3ab"] == "_B").sum() == 0
        assert "A_B" in set(out["CDR3ab"])

    def test_no_umi_columns_falls_back_to_pass_all(self):
        """Datasets without TRA_1_umis / TRB_1_umis columns shouldn't
        crash; the UMI gate degrades to always-pass."""
        adata = _adata([
            {"sample": "S1", "CDR3ab": "X_Y", "CDR3_alpha": "X", "CDR3_beta": "Y"},
            {"sample": "S1", "CDR3ab": "X_Y", "CDR3_alpha": "X", "CDR3_beta": "Y"},
        ])
        out = build_clone_sample_long(adata)
        # Both cells pass; denominator = 2; freq = 1.0.
        row = out[out["sample"] == "S1"].iloc[0]
        assert abs(row["frequency"] - 1.0) < 1e-9

    def test_underscore_cdr3ab_excluded(self):
        """``CDR3ab == "_"`` (alpha and beta both empty) was already
        filtered before #80 — verify it still is."""
        adata = _adata([
            {"sample": "S1", "CDR3ab": "_", "CDR3_alpha": "", "CDR3_beta": "",
             "TRA_1_umis": 0, "TRB_1_umis": 0},
            {"sample": "S1", "CDR3ab": "X_Y", "CDR3_alpha": "X", "CDR3_beta": "Y",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
        ])
        out = build_clone_sample_long(adata)
        assert (out["CDR3ab"] != "_").all()
