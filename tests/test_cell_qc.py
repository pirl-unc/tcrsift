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
    SOLID_TUMOR_LINEAGE_PROGRAMS,
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


class TestFunnelDoubletTunables:
    """Cross-lineage stringency + CD gene symbols are reachable through the
    funnel (#327), not fixed internally."""

    def _xlineage_removed(self, adata, **kw):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, wf = cell_qc_funnel(
                adata, min_genes=3, min_counts=5, min_mito_pct=0.0, **kw
            )
        row = wf[wf["step"] == "xlineage"]
        return int(row["removed"].iloc[0]) if len(row) else 0

    def test_xlineage_min_programs_threaded(self, qc_adata):
        # The 5 T+myeloid doublets are a 2-program event: default min_programs=2
        # flags them; requiring 3 disjoint programs lets them through.
        assert self._xlineage_removed(qc_adata, xlineage_min_programs=2) == 5
        assert self._xlineage_removed(qc_adata, xlineage_min_programs=3) == 0

    def test_xlineage_score_min_threaded(self, qc_adata):
        # A punishingly high per-program score floor calls nothing "on".
        assert self._xlineage_removed(qc_adata, xlineage_score_min=99.0) == 0

    def test_cd_gene_symbols_threaded(self):
        # A double-positive T cell whose CD8 sits under a non-standard symbol is
        # missed by the defaults but flagged when cd8_genes points at it.
        genes = ["CD3E", "CD4", "CD8x"] + [f"G{i}" for i in range(40)]
        gi = {g: i for i, g in enumerate(genes)}
        X = np.zeros((4, len(genes)))
        for fg in [f"G{i}" for i in range(40)]:
            X[:, gi[fg]] = 30
        X[0, [gi["CD3E"], gi["CD4"], gi["CD8x"]]] = 40  # DP under CD8x
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                           obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]))
        assert not cd4_cd8_doublet_mask(adata).any()            # default CD8A/B
        assert cd4_cd8_doublet_mask(adata, cd8_genes=("CD8x",))[0]


class TestSolidTumorPreset:
    def test_keeps_fibroblast_folds_osteoclast(self):
        # The solid-tumor preset KEEPS the fibroblast program (#334): tumor
        # collagen is the only handle cross-lineage detection has on tumor+immune
        # doublets. It folds osteoclast markers into myeloid.
        assert "fibroblast" in SOLID_TUMOR_LINEAGE_PROGRAMS
        assert SOLID_TUMOR_LINEAGE_PROGRAMS["fibroblast"] == [
            "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"
        ]
        for g in ("ACP5", "CTSK", "MMP9"):
            assert g in SOLID_TUMOR_LINEAGE_PROGRAMS["myeloid"]
        # Same lineages as the default (plus osteoclast markers on myeloid).
        assert set(SOLID_TUMOR_LINEAGE_PROGRAMS) == set(DEFAULT_LINEAGE_PROGRAMS)

    def test_tumor_immune_doublet_needs_fibroblast_handle(self):
        # A tumor+T doublet (collagen + CD3) trips fibroblast+T = 2 programs and
        # is flagged; a PURE tumor cell (collagen only) is 1 program and is not.
        # Include the full lineage-gene universe (mostly zero) so score_genes has
        # control genes, as in the qc_adata fixture.
        lineage_genes = sorted(
            {g for gs in SOLID_TUMOR_LINEAGE_PROGRAMS.values() for g in gs}
        )
        filler = [f"G{i}" for i in range(120)]
        genes = list(dict.fromkeys(lineage_genes + filler))
        gi = {g: i for i, g in enumerate(genes)}
        collagen = ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"]
        tcell = ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "CD7"]
        rng = np.random.default_rng(0)
        n = 30  # 0-9 pure tumor, 10-19 tumor+T doublets, 20-29 clean T
        X = np.zeros((n, len(genes)))
        for fg in filler:  # spread background levels so score_genes control bins fill
            X[:, gi[fg]] = rng.poisson(rng.integers(2, 40), size=n)
        for g in collagen:
            X[0:20, gi[g]] = 60                # tumor + doublets carry collagen
        for g in tcell:
            X[10:30, gi[g]] = 60               # doublets + clean T carry CD3
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                           obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
        # require_high_umi=False isolates the program-counting logic (#334): the
        # UMI-outlier corroboration is a separate lever, not the fibroblast handle.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mask = cross_lineage_doublets(
                adata, SOLID_TUMOR_LINEAGE_PROGRAMS, require_high_umi=False
            )
        assert not mask[0:10].any()   # pure tumor: 1 program (fibroblast), not flagged
        assert mask[10:20].all()      # tumor+T doublet: caught via fibroblast handle
        assert not mask[20:30].any()  # clean T cells


class TestQcNoneDisableAndWarn:
    """#333: None means 'disable this bound', gentle defaults, loud cull warning."""

    def test_min_mito_floor_off_by_default(self, qc_adata):
        # qc_adata's non-doublet cells are ~0% mito; the OLD 2% default floor
        # would cull them. Bare-ish call must NOT apply a mito floor now.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, wf = cell_qc_funnel(qc_adata, min_genes=3, min_counts=5)
        assert "min_mito" not in set(wf["step"])   # floor off by default
        # The ~0%-mito majority survives (only doublet gates + max_mito trim);
        # the old 2% floor would have culled nearly all of them.
        assert out.n_obs > 40

    def test_none_disables_max_genes(self):
        genes = [f"G{i}" for i in range(50)]
        X = np.zeros((4, len(genes)))
        X[0, :40] = 5.0          # cell 0: 40 genes (high)
        X[1:, :8] = 5.0          # cells 1-3: 8 genes each
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                           obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # max_genes unset (None) → the high-gene cell is kept.
            keep, _ = cell_qc_funnel(adata, min_genes=1, min_counts=1)
            # max_genes=20 → it is dropped.
            drop, wf = cell_qc_funnel(adata, min_genes=1, min_counts=1, max_genes=20)
        assert keep.n_obs == 4
        assert drop.n_obs == 3
        assert "max_genes" in set(wf["step"])

    def test_loud_warning_on_large_single_step_cull(self, qc_adata, caplog):
        import logging
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with caplog.at_level(logging.WARNING, logger="tcrsift.qc"):
                out, wf = cell_qc_funnel(
                    qc_adata, min_genes=3, min_counts=5, min_mito_pct=2.0
                )
        # The 2% floor is the reported footgun: it culls the ~0%-mito majority.
        assert wf.set_index("step")["removed"]["min_mito"] > 10
        assert any(
            "min_mito" in r.getMessage() and "removed" in r.getMessage()
            for r in caplog.records
        )

    def test_config_restores_opinionated_floor(self, qc_adata):
        from tcrsift.config import LoadConfig

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            off, _ = cell_qc_funnel(qc_adata, min_genes=3, min_counts=5)
            cfg, wf = cell_qc_funnel(
                qc_adata, config=LoadConfig(min_genes=3, min_counts=5)
            )
        # LoadConfig carries the opinionated PBMC min_mito_pct=2.0 → floor applies.
        assert "min_mito" in set(wf["step"])
        assert cfg.n_obs < off.n_obs

    def test_low_coverage_flag_survives_disabled_min_genes(self):
        # min_genes=None must not crash the 2×min_genes low-coverage default.
        genes = [f"G{i}" for i in range(30)]
        X = np.zeros((3, len(genes)))
        X[:, :20] = 3.0
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                           obs=pd.DataFrame(index=[f"c{i}" for i in range(3)]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, _ = cell_qc_funnel(adata, min_genes=None, min_counts=1)
        assert out.n_obs == 3
        assert "qc_low_coverage" in out.obs.columns


class TestPctCountsMtPreserved:
    """cell_qc_funnel leaves a caller's pre-existing pct_counts_mt untouched and
    still gates the same cells (#338)."""

    def test_preexisting_pct_counts_mt_not_clobbered(self, qc_adata):
        # A caller-supplied pct_counts_mt (e.g. from scanpy calculate_qc_metrics)
        # must survive the funnel unchanged on the returned subset.
        sentinel = np.arange(qc_adata.n_obs, dtype=float) / 1000.0
        qc_adata.obs["pct_counts_mt"] = sentinel
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, _ = cell_qc_funnel(
                qc_adata, min_genes=3, min_counts=5, max_mito_pct=None
            )
        # No cells dropped for mito (max_mito disabled), so every survivor keeps
        # the exact value the caller stored — not tcrsift's own arithmetic.
        expected = pd.Series(sentinel, index=qc_adata.obs.index).loc[out.obs.index]
        assert np.array_equal(out.obs["pct_counts_mt"].to_numpy(), expected.to_numpy())

    def test_computes_pct_counts_mt_when_absent(self, qc_adata):
        # With no caller value, the funnel still populates the column.
        assert "pct_counts_mt" not in qc_adata.obs.columns
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, _ = cell_qc_funnel(qc_adata, min_genes=3, min_counts=5)
        assert "pct_counts_mt" in out.obs.columns
        assert "_pct_mt" not in out.obs.columns  # internal gate column cleaned up

    def test_mito_gate_independent_of_caller_column(self, qc_adata):
        # Cells 18-19 are high-mito. A garbage caller pct_counts_mt (all zero) must
        # NOT change which cells the mito gate drops — the gate uses tcrsift's own
        # arithmetic, so survivors match the no-column run exactly.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref, _ = cell_qc_funnel(qc_adata, min_genes=3, min_counts=5, max_mito_pct=10.0)
        poisoned = qc_adata.copy()
        poisoned.obs["pct_counts_mt"] = 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, _ = cell_qc_funnel(poisoned, min_genes=3, min_counts=5, max_mito_pct=10.0)
        assert list(out.obs.index) == list(ref.obs.index)


class TestConfigUnwrap:
    """cell_qc_funnel(config=TCRsiftConfig()) reads the nested .load section
    instead of silently ignoring it (#341)."""

    def test_top_level_tcrsift_config_unwrapped(self, qc_adata):
        from tcrsift.config import LoadConfig, TCRsiftConfig

        cfg = TCRsiftConfig(
            load=LoadConfig(
                min_genes=3, max_genes=15000, min_counts=5, max_counts=100000,
                min_mito_pct=0.0, max_mito_pct=40.0,
            )
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out, wf = cell_qc_funnel(qc_adata, config=cfg)
        # Same outcome as passing the flat LoadConfig directly (2 high-mito dropped).
        assert wf.set_index("step")["removed"]["max_mito"] == 2
        assert out.n_obs == 48

    def test_flat_load_config_still_works(self, qc_adata):
        from tcrsift.config import LoadConfig

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            flat, _ = cell_qc_funnel(
                qc_adata,
                config=LoadConfig(
                    min_genes=3, max_genes=15000, min_counts=5, max_counts=100000,
                    min_mito_pct=0.0, max_mito_pct=40.0,
                ),
            )
        assert flat.n_obs == 48
