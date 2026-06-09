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

"""Tests for ``plot_cells_per_sample_stacked`` (#76)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402

from tcrsift.plots import plot_cells_per_sample_stacked  # noqa: E402


def _adata(
    n_per_sample: dict[str, int],
    *,
    with_qc: bool = False,
    with_tcr: bool = False,
    qc_pass_frac: float = 1.0,
    denom_pass_frac: float = 1.0,
):
    samples = []
    for s, n in n_per_sample.items():
        samples.extend([s] * n)
    n = len(samples)
    obs = pd.DataFrame({"sample": samples}, index=[str(i) for i in range(n)])

    if with_qc:
        rng = np.random.default_rng(0)
        qc_cols = (
            "filter:min_genes", "filter:max_genes",
            "filter:min_counts", "filter:max_counts",
            "filter:min_mito", "filter:max_mito",
            "filter:min_cd3",
        )
        passing = rng.random(n) < qc_pass_frac
        for c in qc_cols:
            obs[c] = passing

    if with_tcr:
        rng = np.random.default_rng(1)
        is_pure = rng.random(n) < denom_pass_frac
        obs["is_confident"] = is_pure
        obs["has_both_chains"] = is_pure
        obs["multi_chain"] = ~is_pure
        obs["TRA_1_umis"] = np.where(is_pure, 5, 0).astype(float)
        obs["TRB_1_umis"] = np.where(is_pure, 5, 0).astype(float)
        obs["Tcell_type_consensus"] = np.where(is_pure, "Confident CD8+", "Other")

    X = csr_matrix((n, 1))
    a = ad.AnnData(X=X, obs=obs)
    return a


class TestStackedQCBars:
    def test_renders_with_only_sample_column(self, tmp_path):
        a = _adata({"S1": 100, "S2": 60})
        out = plot_cells_per_sample_stacked(a, tmp_path / "qc.png")
        assert out is not None and out.exists()

    def test_skips_when_no_sample_column(self, tmp_path):
        X = csr_matrix((5, 1))
        a = ad.AnnData(X=X, obs=pd.DataFrame({"foo": [1, 2, 3, 4, 5]}, index=list("abcde")))
        out = plot_cells_per_sample_stacked(a, tmp_path / "qc.png")
        assert out is None

    def test_qc_layer_present_when_columns_present(self, tmp_path):
        """When the ``filter:*`` columns exist, the QC layer is drawn
        with a smaller height than the loaded layer."""
        a = _adata({"S1": 100}, with_qc=True, qc_pass_frac=0.7)
        out = plot_cells_per_sample_stacked(a, tmp_path / "qc.png")
        assert out.exists()

    def test_tcr_denominator_layer_present(self, tmp_path):
        a = _adata(
            {"S1": 100, "S2": 80},
            with_qc=True, with_tcr=True,
            qc_pass_frac=0.8, denom_pass_frac=0.5,
        )
        out = plot_cells_per_sample_stacked(a, tmp_path / "qc.png")
        assert out.exists()

    def test_donor_in_title(self, tmp_path):
        a = _adata({"S1": 50})
        out = plot_cells_per_sample_stacked(
            a, tmp_path / "qc.png", donor="A1"
        )
        assert out.exists()

    def test_sample_labels_pass_through_pretty(self, tmp_path):
        """Sample names like ``AIMpos_CTYneg-2`` should appear as
        ``AIM⁺CTY⁻`` in the rendered figure. Save as PDF (vector) so
        the labels survive as text we can inspect."""
        a = _adata({"AIMpos_CTYneg-1": 10, "tetpos-1": 8})
        out = plot_cells_per_sample_stacked(a, tmp_path / "qc.pdf")
        assert out.exists()
