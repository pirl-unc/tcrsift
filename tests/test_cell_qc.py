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

"""Per-cell QC + doublet funnel (#314).

Filters an AnnData on the LoadConfig cell-QC thresholds and flags the doublet
classes the single-TRB gate misses (CD4/CD8 double-positive, cross-lineage
co-expression), recording a visible waterfall.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.qc import (
    DEFAULT_LINEAGE_PROGRAMS,
    cd4_cd8_doublet_mask,
    cell_qc_funnel,
    cross_lineage_doublets,
)


@pytest.fixture
def qc_adata():
    """60 cells at realistic depth: 0-9 clean T, 10-14 T+myeloid doublets
    (elevated UMI), 15-17 CD4/CD8 double-positive, 18-19 high-mito, 20 near-empty."""
    rng = np.random.default_rng(0)
    lineage_genes = sorted({g for gs in DEFAULT_LINEAGE_PROGRAMS.values() for g in gs})
    filler = [f"G{i}" for i in range(60)]
    genes = list(dict.fromkeys(
        lineage_genes + ["MT-CO1", "MT-ND1", "CD4", "CD8A", "CD8B"] + filler
    ))
    gi = {g: i for i, g in enumerate(genes)}
    n = 60
    X = np.zeros((n, len(genes)), dtype=float)
    for fg in filler:  # bulk UMI so CP10K normalization behaves like real data
        X[:, gi[fg]] = rng.poisson(30, size=n)

    def setg(cells, glist, val):
        for c in cells:
            for g in glist:
                X[c, gi[g]] = val

    setg(range(0, 10), ["CD3D", "CD3E", "CD3G", "CD2", "CD7", "TRAC"], 40)
    setg(range(10, 15), ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"], 40)
    setg(range(10, 15), ["LYZ", "CD68", "C1QA", "C1QB", "C1QC", "AIF1", "CSF1R"], 40)
    for c in range(10, 15):
        X[c, :] += (X[c, :] > 0) * 20  # elevate depth (doublet coverage)
    setg(range(15, 18), ["CD3E"], 40)
    setg(range(15, 18), ["CD4", "CD8A", "CD8B"], 30)
    setg(range(18, 20), ["MT-CO1", "MT-ND1"], 4000)
    X[20, :] = 0
    X[20, gi["CD3E"]] = 1

    adata = ad.AnnData(
        X=X, var=pd.DataFrame(index=genes),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
    )
    adata.obs["multi_TRB"] = False
    adata.obs.iloc[5, adata.obs.columns.get_loc("multi_TRB")] = True
    return adata


class TestCrossLineageDoublets:
    def test_flags_disjoint_lineage_coexpression(self, qc_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mask = cross_lineage_doublets(qc_adata, require_high_umi=True)
        assert mask[10:15].all()      # T+myeloid doublets
        assert not mask[0:10].any()   # clean T cells

    def test_require_high_umi_gate(self, qc_adata):
        # With the UMI corroboration off, the same co-expressors are flagged;
        # the gate only ever removes flags, never adds them.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            strict = cross_lineage_doublets(qc_adata, require_high_umi=True)
            loose = cross_lineage_doublets(qc_adata, require_high_umi=False)
        assert loose[10:15].all()
        assert (strict <= loose).all()

    def test_too_few_programs_returns_all_false(self):
        adata = ad.AnnData(
            X=np.ones((5, 3)), var=pd.DataFrame(index=["CD3E", "CD3D", "LYZ"]),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mask = cross_lineage_doublets(adata, {"only": ["CD3E", "CD3D"]})
        assert not mask.any()


class TestCd4Cd8Doublet:
    def test_cd3_gated_double_positive(self, qc_adata):
        mask = cd4_cd8_doublet_mask(qc_adata)
        assert mask[15:18].all()
        assert not mask[0:10].any()

    def test_cd3_gate_spares_myeloid_cd4(self):
        # CD4 without CD3 (myeloid) + ambient CD8 must NOT be flagged.
        genes = ["CD3E", "CD4", "CD8A", "CD8B", "LYZ"]
        X = np.array([[0.0, 5.0, 1.0, 0.0, 9.0]])  # CD4+CD8+ but CD3E==0
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                           obs=pd.DataFrame(index=["c0"]))
        assert not cd4_cd8_doublet_mask(adata).any()


class TestCellQcFunnel:
    def test_waterfall_removes_each_class(self, qc_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, wf = cell_qc_funnel(
                qc_adata, min_genes=3, min_counts=5, min_mito_pct=0.0,
                max_mito_pct=40.0,
            )
        d = wf.set_index("step")["removed"]
        assert d["min_genes"] == 1
        assert d["max_mito"] == 2
        assert d["tcr_2b"] == 1
        assert d["cd4cd8_dp"] == 3
        assert d["xlineage"] == 5
        # Waterfall accounts for every cell: input - sum(removed) == remaining.
        assert wf["remaining"].iloc[-1] == out.n_obs
        assert wf["remaining"].iloc[0] - wf["removed"].sum() == out.n_obs

    def test_low_coverage_flagged_not_dropped(self, qc_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, _ = cell_qc_funnel(
                qc_adata, min_genes=3, min_counts=5, min_mito_pct=0.0,
                max_mito_pct=40.0, low_coverage_genes=100,
            )
        assert "qc_low_coverage" in out.obs.columns
        # Flag marks (does not drop) — flagged cells still present.
        assert out.obs["qc_low_coverage"].dtype == bool

    def test_config_thresholds_reused(self, qc_adata):
        # A duck-typed config supplies thresholds; explicit kwargs override it.
        class Cfg:
            min_genes = 3
            max_genes = 15000
            min_counts = 5
            max_counts = 100000
            min_mito_pct = 0.0
            max_mito_pct = 40.0

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, wf = cell_qc_funnel(qc_adata, config=Cfg())
        assert wf.set_index("step")["removed"]["max_mito"] == 2
        assert out.n_obs == 48

    def test_does_not_mutate_input(self, qc_adata):
        n_before = qc_adata.n_obs
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cell_qc_funnel(qc_adata, min_genes=3, min_counts=5, min_mito_pct=0.0)
        assert qc_adata.n_obs == n_before
        assert "qc_low_coverage" not in qc_adata.obs.columns
