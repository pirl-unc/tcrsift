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

    def test_umi_gated_denominator(self):
        """Cells with TRA_1_umis < 2 or TRB_1_umis < 2 should be
        excluded from the denominator — same rule as
        ``aggregate_clonotypes``."""
        adata = _adata([
            # Two passing cells.
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 5, "TRB_1_umis": 5},
            # One cell failing the UMI gate (TRA < 2). Should still be
            # listed as cells=1 for the clone but excluded from the
            # complete-clone denominator.
            {"sample": "S1", "CDR3ab": "A_A", "CDR3_alpha": "A", "CDR3_beta": "A",
             "TRA_1_umis": 0, "TRB_1_umis": 5},
        ])
        out = build_clone_sample_long(adata)
        # Cells count includes all rows with non-empty CDR3ab (3);
        # denominator is the UMI-gated count (2).
        row = out[(out["CDR3ab"] == "A_A") & (out["sample"] == "S1")].iloc[0]
        assert row["cells"] == 3
        assert abs(row["frequency"] - 3 / 2) < 1e-9

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
