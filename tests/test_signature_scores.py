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

"""Tests for ``compute_signature_scores_per_clonotype`` (#74).

The per-clone signature scores need to be available in ``clonotypes.csv``
so downstream candidate-selection (#75), 2D maps, and pairwise plots
don't have to re-load the h5ad and recompute. This module verifies:

- mean-log1p computation across a gene set,
- per-clonotype aggregation,
- CD8+ restriction when a CD8 column is available,
- graceful handling of missing genes / empty results.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from tcrsift.gex import compute_signature_scores_per_clonotype


def _per_cell(rows):
    return pd.DataFrame(rows)


class TestComputeSignatureScores:
    def test_mean_log1p_across_gene_set(self):
        # Two cells in one clone with TNFRSF9=expr1, MKI67=expr2 each.
        df = _per_cell([
            {"CDR3_pair": "A", "gex.TNFRSF9": 4.0, "gex.MKI67": 0.0, "gex.CD8A": 1},
            {"CDR3_pair": "A", "gex.TNFRSF9": 0.0, "gex.MKI67": 4.0, "gex.CD8A": 1},
        ])
        sigs = {"antigen_response": ("TNFRSF9", "MKI67")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, verbose=False
        )
        # Cell 1: mean(log1p(4), log1p(0)) = mean(1.609, 0) = 0.8047
        # Cell 2: mean(log1p(0), log1p(4)) = 0.8047
        # Clone mean = 0.8047
        assert result.shape == (1, 2)
        assert result["CDR3_pair"].iloc[0] == "A"
        assert abs(result["signature_antigen_response"].iloc[0] - np.log(5) / 2) < 1e-6

    def test_separates_clonotypes(self):
        df = _per_cell([
            {"CDR3_pair": "A", "gex.TNFRSF9": 10.0, "gex.MKI67": 10.0, "gex.CD8A": 1},
            {"CDR3_pair": "B", "gex.TNFRSF9": 0.0, "gex.MKI67": 0.0, "gex.CD8A": 1},
        ])
        sigs = {"antigen_response": ("TNFRSF9", "MKI67")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, verbose=False
        ).set_index("CDR3_pair")
        assert result.loc["A", "signature_antigen_response"] > 0
        assert result.loc["B", "signature_antigen_response"] == 0

    def test_three_default_signatures_emitted(self):
        """When ``signatures`` is None, the five canonical T-cell sets
        from ``tcrsift.signatures`` are scored — but the issue #74 ask
        focuses on three for clonotype output. Verify the focal three
        all appear in the result."""
        df = _per_cell([
            {
                "CDR3_pair": "A",
                "gex.TNFRSF9": 1.0, "gex.MKI67": 1.0,
                "gex.PRF1": 1.0,    "gex.GZMB": 1.0,
                "gex.CXCL13": 1.0,  "gex.ENTPD1": 1.0,
                "gex.IFNG": 1.0, "gex.GNLY": 1.0, "gex.NKG7": 1.0,
                "gex.PDCD1": 1.0, "gex.LAG3": 1.0, "gex.HAVCR2": 1.0,
                "gex.TIGIT": 1.0, "gex.TOX": 1.0, "gex.CTLA4": 1.0,
                "gex.CD8A": 1,
            },
        ])
        result = compute_signature_scores_per_clonotype(df, verbose=False)
        for col in (
            "signature_antigen_response",
            "signature_cytolytic",
            "signature_tumor_reactive",
        ):
            assert col in result.columns

    def test_cd8_only_restriction_active_by_default(self):
        """A clone with mixed CD4/CD8 cells should only count the
        CD8+ cell's expression when ``cd8_only=True`` (the default)."""
        df = _per_cell([
            # CD4 cell (CD8A=0) with high signature.
            {"CDR3_pair": "A", "gex.TNFRSF9": 100.0, "gex.MKI67": 100.0, "gex.CD8A": 0},
            # CD8 cell with zero signature.
            {"CDR3_pair": "A", "gex.TNFRSF9": 0.0,   "gex.MKI67": 0.0,   "gex.CD8A": 1},
        ])
        sigs = {"antigen_response": ("TNFRSF9", "MKI67")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, verbose=False
        )
        # CD4 cell filtered out → only the CD8 cell scores, signature = 0.
        assert result["signature_antigen_response"].iloc[0] == 0

    def test_cd8_only_off_includes_all_cells(self):
        df = _per_cell([
            {"CDR3_pair": "A", "gex.TNFRSF9": 100.0, "gex.MKI67": 100.0, "gex.CD8A": 0},
            {"CDR3_pair": "A", "gex.TNFRSF9": 0.0,   "gex.MKI67": 0.0,   "gex.CD8A": 1},
        ])
        sigs = {"antigen_response": ("TNFRSF9", "MKI67")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, cd8_only=False, verbose=False
        )
        assert result["signature_antigen_response"].iloc[0] > 0

    def test_missing_signature_genes_yields_nan(self):
        """A signature whose genes aren't in the augmented frame
        emits NaN for every clone (not zero, not an error)."""
        df = _per_cell([
            {"CDR3_pair": "A", "gex.UNRELATED": 1.0, "gex.CD8A": 1},
        ])
        sigs = {"missing": ("NOT_THERE_1", "NOT_THERE_2")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, verbose=False
        )
        assert result["signature_missing"].isna().all()

    def test_empty_after_cd8_filter_returns_nan_keyed_clones(self):
        """If CD8 filter kills every row, we should still emit one row
        per clonotype (NaN scores) so the merge into clonotypes.csv
        doesn't silently drop clones."""
        df = _per_cell([
            {"CDR3_pair": "A", "gex.TNFRSF9": 1.0, "gex.MKI67": 1.0, "gex.CD8A": 0},
            {"CDR3_pair": "B", "gex.TNFRSF9": 1.0, "gex.MKI67": 1.0, "gex.CD8A": 0},
        ])
        sigs = {"antigen_response": ("TNFRSF9", "MKI67")}
        result = compute_signature_scores_per_clonotype(
            df, signatures=sigs, verbose=False
        )
        assert set(result["CDR3_pair"]) == {"A", "B"}
        assert result["signature_antigen_response"].isna().all()

    def test_missing_group_col_raises(self):
        df = _per_cell([{"gex.TNFRSF9": 1.0, "gex.CD8A": 1}])
        try:
            compute_signature_scores_per_clonotype(df, verbose=False)
        except ValueError as e:
            assert "CDR3_pair" in str(e)
        else:
            raise AssertionError("expected ValueError")
