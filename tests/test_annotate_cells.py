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

"""Generalized cell-type annotator (#312).

Per-cluster background-subtracted signature argmax with biology-aware gates
(Treg=FOXP3, γδ=TRDC-only, MAIT=SLC4A10+KLRB1, NK-vs-CD8), superseding the
CD4/CD8-only phenotype path.
"""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import pandas as pd

from tcrsift.annotate_cells import (
    AnnotationGates,
    MarkerCountOverride,
    annotate_cells,
    annotate_clusters,
    compose_phenotype_labels,
    gates_from_phenotype_config,
)
from tcrsift.config import PhenotypeConfig
from tcrsift.signatures import (
    CELL_TYPE_SIGNATURES,
    LINEAGE_GENES,
    PBMC_CULTURE_TYPES,
    T_STATE_SIGNATURES,
)

_MARKERS = sorted(
    {g for v in CELL_TYPE_SIGNATURES.values() for g in v}
    | {g for v in T_STATE_SIGNATURES.values() for g in v}
    | {g for v in LINEAGE_GENES.values() for g in v}
    | {"KLRB1", "SLC4A10", "TRDC", "IL13", "IL17A", "IFNG"}
)
_FILLER = [f"F{i}" for i in range(120)]
_GENES = list(dict.fromkeys(_MARKERS + _FILLER))
_GI = {g: i for i, g in enumerate(_GENES)}


def _build(cluster_specs):
    """AnnData at realistic depth (bulk filler ~3000 UMI so CP10K behaves), each
    cluster defined by marker bumps; genes not named are 0."""
    rng = np.random.default_rng(0)
    rows, leiden = [], []
    for cl, bumps in cluster_specs.items():
        for _ in range(40):
            x = np.zeros(len(_GENES))
            for f in _FILLER:
                x[_GI[f]] = rng.poisson(25)
            for g, v in bumps.items():
                if g in _GI:
                    x[_GI[g]] = v
            rows.append(x)
            leiden.append(cl)
    X = np.array(rows)
    return ad.AnnData(
        X=X, var=pd.DataFrame(index=_GENES),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(leiden)},
            index=[f"c{i}" for i in range(len(rows))],
        ),
    )


def _labels(adata, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = annotate_cells(adata, allowed_types=PBMC_CULTURE_TYPES, **kw)
    return out.obs.groupby("leiden", observed=True)["phenotype_base"].first().to_dict()


class TestLineageArgmax:
    def test_cd8_cd4_nk_b(self):
        labels = _labels(_build({
            "0": {"CD3D": 30, "CD3E": 30, "CD3G": 25, "CD8A": 30, "CD8B": 25,
                  "GZMB": 40, "PRF1": 30, "GZMA": 30, "GZMH": 20, "NKG7": 20},
            "1": {"CD3D": 30, "CD3E": 30, "CD3G": 25, "CD4": 30, "CCR7": 25,
                  "SELL": 20, "TCF7": 20, "IL7R": 18},
            "2": {"NKG7": 40, "KLRD1": 30, "KLRF1": 25, "NCR1": 20, "GNLY": 40,
                  "NCAM1": 20, "PRF1": 30, "TYROBP": 25, "KLRC1": 15},
            "3": {"MS4A1": 40, "CD79A": 30, "CD79B": 30, "CD19": 20, "BANK1": 18},
        }))
        assert labels["0"].startswith("CD8")
        assert labels["1"].startswith("CD4")
        assert labels["2"] == "NK cell"
        assert labels["3"] == "B cell"


class TestStateGates:
    def test_treg_requires_foxp3(self):
        # Cluster A: CD4 with FOXP3 → Treg. Cluster B: CD4 activated with
        # IL2RA/CTLA4 HIGH but FOXP3 ZERO → NOT Treg (gate fails).
        labels = _labels(_build({
            "A": {"CD3D": 30, "CD3E": 30, "CD4": 30, "FOXP3": 20, "IL2RA": 18,
                  "CTLA4": 18, "IKZF2": 15, "CCR8": 12},
            "B": {"CD3D": 30, "CD3E": 30, "CD4": 30, "IL2RA": 25, "CTLA4": 25,
                  "TNFRSF9": 20, "TNFRSF4": 18, "CD69": 15},
        }))
        assert "Treg" in labels["A"]
        assert "Treg" not in labels["B"]  # no FOXP3 → not Treg
        assert labels["B"].startswith("CD4")

    def test_gd_t_called_on_trdc_only(self):
        # A CD3+ cluster with the δ-constant TRDC high is γδ T.
        labels = _labels(_build({
            "g": {"CD3D": 30, "CD3E": 30, "CD3G": 25, "TRDC": 20, "CD8A": 10,
                  "GZMB": 20, "NKG7": 15},
        }))
        assert labels["g"] == "gd T"

    def test_mait_on_slc4a10_and_klrb1(self):
        labels = _labels(_build({
            "m": {"CD3D": 30, "CD3E": 30, "CD3G": 25, "CD8A": 20, "SLC4A10": 20,
                  "KLRB1": 20, "GZMK": 15},
        }))
        assert labels["m"] == "MAIT"


class TestAllowedTypes:
    def test_restriction_prevents_out_of_panel_win(self):
        # A myeloid cluster: with PBMC_CULTURE_TYPES it must be a myeloid type,
        # never a stromal label (fibroblast/pericyte are excluded from the panel).
        labels = _labels(_build({
            "myeloid": {"CD68": 30, "C1QA": 25, "C1QB": 25, "C1QC": 20,
                        "LYZ": 30, "CD14": 20, "FCN1": 20},
        }))
        assert labels["myeloid"] in {"Macrophage", "Monocyte", "Dendritic cell"}


class TestComposeAndConfig:
    def test_compose_merges_tiny_fragments(self):
        base = {"0": "Macrophage", "1": "Macrophage"}
        tops = {"0": "C1QA", "1": "SPP1"}
        sizes = {"0": 1000, "1": 5}  # cluster 1 is 0.5% of macrophages
        out = compose_phenotype_labels(base, tops, sizes, min_subtype_frac=0.01)
        assert out["0"] == "Macrophage · C1QA"  # big fragment keeps suffix
        assert out["1"] == "Macrophage"          # <1% fragment loses suffix

    def test_compose_protects_rare_subtypes(self):
        out = compose_phenotype_labels(
            {"0": "pDC"}, {"0": "LILRA4"}, {"0": 3}, min_subtype_frac=0.5,
            protect=frozenset({"pDC"}),
        )
        assert out["0"] == "pDC · LILRA4"

    def test_gates_from_phenotype_config(self):
        cfg = PhenotypeConfig(treg_foxp3_min=0.9, gd_tcr_min=0.7, mait_min=0.4)
        gates = gates_from_phenotype_config(cfg)
        assert isinstance(gates, AnnotationGates)
        assert gates.state_defining_genes["Treg"] == ("FOXP3", 0.9)
        assert gates.gd_tcr_min == 0.7
        assert gates.mait_min == 0.4


class TestBackCompatShim:
    def test_classify_tcell_type_still_works(self):
        # The CD4/CD8-ratio classifier stays available (thin shim, #312).
        from tcrsift.phenotype import classify_tcell_type

        assert "CD8" in classify_tcell_type(0.0, 10.0)
        assert "CD4" in classify_tcell_type(10.0, 0.0)


# --------------------------------------------------------------------------- #
# Marker-count cluster override (#325)
# --------------------------------------------------------------------------- #
_CTA = tuple(f"CTA{i}" for i in range(1, 7))
_TF = ("TF1", "TF2")
_OV_GENES = list(
    dict.fromkeys(
        ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]  # fibroblast (collagen)
        + list(_CTA) + list(_TF)
        + ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB"]
        + [f"F{i}" for i in range(120)]
    )
)
_OVI = {g: i for i, g in enumerate(_OV_GENES)}


def _build_tumor():
    """Three clusters at realistic depth: a collagen+CTA tumor (mis-argmaxes to
    Fibroblast), a low-CTA-coverage-but-high-TF tumor (needs rescue), and a clean
    CD8 T cluster with no CTA load."""
    rng = np.random.default_rng(0)
    rows, leiden = [], []

    def cell(bumps):
        x = np.zeros(len(_OV_GENES))
        for i in range(len(_OV_GENES)):
            if _OV_GENES[i].startswith("F"):
                x[i] = rng.poisson(25)
        for g, v in bumps.items():
            x[_OVI[g]] = v
        return x

    collagen = {"COL1A1": 8, "COL1A2": 8, "COL3A1": 8, "DCN": 6, "LUM": 6}
    # tumor: every cell carries 3 distinct CTAs → frac 1.0
    for _ in range(40):
        rows.append(cell({**collagen, "CTA1": 5, "CTA2": 5, "CTA3": 5,
                          "TF1": 4, "TF2": 4}))
        leiden.append("tumor")
    # tumor_lo: only 6/40 cells carry ≥2 CTAs (frac 0.15 < 0.4) but TF is high
    for j in range(40):
        b = {**collagen, "TF1": 6, "TF2": 6}
        if j < 6:
            b.update({"CTA1": 5, "CTA2": 5})
        rows.append(cell(b))
        leiden.append("tumor_lo")
    # tcell: clean CD8, no CTAs
    for _ in range(40):
        rows.append(cell({"CD3D": 6, "CD3E": 6, "CD8A": 6, "CD8B": 6, "GZMB": 6}))
        leiden.append("tcell")

    return ad.AnnData(
        X=np.array(rows), var=pd.DataFrame(index=_OV_GENES),
        obs=pd.DataFrame({"leiden": pd.Categorical(leiden)},
                         index=[f"c{i}" for i in range(len(rows))]),
    )


def _ov_labels(adata, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = annotate_cells(adata, **kw)
    return out, out.obs.groupby("leiden", observed=True)["phenotype_base"].first().to_dict()


class TestMarkerCountOverride:
    def test_relabels_collagen_cluster_as_tumor(self):
        # Without the override the collagen+CTA cluster mis-argmaxes (Fibroblast
        # or Other); with it, positive CTA evidence relabels it Tumor.
        base, base_labels = _ov_labels(_build_tumor())
        assert base_labels["tumor"] != "Tumor"

        ov = MarkerCountOverride("Tumor", _CTA, min_distinct=2, min_cluster_frac=0.4)
        out, labels = _ov_labels(_build_tumor(), overrides=[ov])
        assert labels["tumor"] == "Tumor"
        assert labels["tcell"].startswith("CD8")  # immune cluster untouched

    def test_exposes_per_cell_marker_count(self):
        ov = MarkerCountOverride("Tumor", _CTA, min_distinct=2)
        out, _ = _ov_labels(_build_tumor(), overrides=[ov])
        assert "n_markers_tumor" in out.obs.columns
        counts = out.obs.groupby("leiden", observed=True)["n_markers_tumor"].mean()
        assert counts["tumor"] >= 2.0     # 3 distinct CTAs per cell
        assert counts["tcell"] == 0.0     # no CTA load on T cells

    def test_low_coverage_needs_rescue(self):
        # tumor_lo has only 15% of cells carrying ≥2 CTAs — below min_cluster_frac.
        no_rescue = MarkerCountOverride("Tumor", _CTA, min_distinct=2,
                                        min_cluster_frac=0.4)
        _, labels = _ov_labels(_build_tumor(), overrides=[no_rescue])
        assert labels["tumor_lo"] != "Tumor"

        # With a TF-mean rescue at a lowered fraction floor, it is saved.
        rescued = MarkerCountOverride("Tumor", _CTA, min_distinct=2,
                                      min_cluster_frac=0.4,
                                      rescue=(list(_TF), 1.0, 0.1))
        _, labels2 = _ov_labels(_build_tumor(), overrides=[rescued])
        assert labels2["tumor_lo"] == "Tumor"

    def test_custom_count_col_and_first_match_wins(self):
        # Explicit count_col; two overlapping overrides → first listed wins.
        first = MarkerCountOverride("TumorA", _CTA, min_distinct=2,
                                    min_cluster_frac=0.4, count_col="cta_load")
        second = MarkerCountOverride("TumorB", _CTA, min_distinct=2,
                                     min_cluster_frac=0.4)
        out, labels = _ov_labels(_build_tumor(), overrides=[first, second])
        assert labels["tumor"] == "TumorA"
        assert "cta_load" in out.obs.columns

    def test_no_override_no_count_column(self):
        out, _ = _ov_labels(_build_tumor())
        assert not any(c.startswith("n_markers") for c in out.obs.columns)

    def test_resolved_count_col_default(self):
        assert MarkerCountOverride("Tumor", _CTA).resolved_count_col() == "n_markers_tumor"
        assert MarkerCountOverride("ISG-high", _CTA).resolved_count_col() == "n_markers_isg_high"

    def test_rescue_can_gate_on_any_marker_fraction(self):
        # Low-coverage tumor: sparse dropout leaves NO cell clearing >=2 CTAs, but
        # most cells still show >=1 CTA and the lineage-TF is high. The >=2 rescue
        # (3-tuple, rescue_min_distinct defaults to min_distinct) can't fire; a
        # rescue_min_distinct=1 (4-tuple) rescue gates on the broad any-marker
        # fraction and does (#339) — h37's real osteosarcoma rule.
        rng = np.random.default_rng(1)
        collagen = {"COL1A1": 8, "COL1A2": 8, "COL3A1": 8, "DCN": 6, "LUM": 6}

        def cell(bumps):
            x = np.zeros(len(_OV_GENES))
            for i, g in enumerate(_OV_GENES):
                if g.startswith("F"):
                    x[i] = rng.poisson(25)
            for g, v in bumps.items():
                x[_OVI[g]] = v
            return x

        rows, leiden = [], []
        for j in range(40):  # 34/40 carry exactly ONE CTA: >=1 frac 0.85, >=2 frac 0
            b = {**collagen, "TF1": 6, "TF2": 6}
            if j < 34:
                b["CTA1"] = 5
            rows.append(cell(b)); leiden.append("tumor_sparse")
        for _ in range(40):  # a clean CD8 cluster so there is >1 cluster
            rows.append(cell({"CD3D": 6, "CD3E": 6, "CD8A": 6, "CD8B": 6, "GZMB": 6}))
            leiden.append("tcell")
        adata = ad.AnnData(
            X=np.array(rows), var=pd.DataFrame(index=_OV_GENES),
            obs=pd.DataFrame({"leiden": pd.Categorical(leiden)},
                             index=[f"c{i}" for i in range(len(rows))]),
        )

        # >=2 rescue never fires (no cell clears the 2-CTA bar, so rescue_frac is 0).
        strict = MarkerCountOverride("Tumor", _CTA, min_distinct=2,
                                     min_cluster_frac=0.4, rescue=(list(_TF), 1.0, 0.1))
        _, labels = _ov_labels(adata.copy(), overrides=[strict])
        assert labels["tumor_sparse"] != "Tumor"

        # rescue_min_distinct=1 fires on the broad >=1 fraction.
        broad = MarkerCountOverride("Tumor", _CTA, min_distinct=2,
                                    min_cluster_frac=0.4, rescue=(list(_TF), 1.0, 0.1, 1))
        _, labels2 = _ov_labels(adata.copy(), overrides=[broad])
        assert labels2["tumor_sparse"] == "Tumor"

    def test_three_tuple_rescue_unchanged(self):
        # Back-compat: the existing 3-tuple rescue still saves the tumor_lo cluster.
        rescued = MarkerCountOverride("Tumor", _CTA, min_distinct=2,
                                      min_cluster_frac=0.4, rescue=(list(_TF), 1.0, 0.1))
        _, labels = _ov_labels(_build_tumor(), overrides=[rescued])
        assert labels["tumor_lo"] == "Tumor"


class TestInjectableStateRegistries:
    """Caller T-state/B-state registries are considered instead of being silently
    dropped by the module-level defaults (#337)."""

    def _t_cell_obs(self):
        # A single CD8 T cluster with a caller-scored state column already present.
        n = 12
        obs = pd.DataFrame(index=[f"c{i}" for i in range(n)])
        obs["leiden"] = pd.Categorical(["0"] * n)
        obs["ct::T cell"] = 1.0
        obs["lin::CD8"] = 1.0
        obs["lin::CD4"] = 0.0
        obs["lin::CD3"] = 2.0
        obs["lin::NKspec"] = 0.0
        obs["st::MyState"] = 1.0  # scored by a caller, not in tcrsift's registry
        return ad.AnnData(X=np.zeros((n, 1)), var=pd.DataFrame(index=["G0"]), obs=obs)

    def test_default_registry_drops_caller_state(self):
        adata = self._t_cell_obs()
        labels = annotate_clusters(adata, reference={"T cell": ["G0"]})
        assert "MyState" not in labels["0"]  # dropped by the built-in T-state keys

    def test_injected_registry_reads_caller_state(self):
        adata = self._t_cell_obs()
        labels = annotate_clusters(
            adata, reference={"T cell": ["G0"]},
            t_state_reference={"MyState": ["G0"]},
        )
        assert labels["0"] == "CD8 MyState"

    def test_warns_on_scored_state_absent_from_registry(self, caplog):
        adata = self._t_cell_obs()
        with caplog.at_level(logging.WARNING, logger="tcrsift.annotate_cells"):
            annotate_clusters(adata, reference={"T cell": ["G0"]})
        assert any("st::MyState" in r.message and "ignored" in r.message
                   for r in caplog.records)

    def test_all_defaults_do_not_warn(self, caplog):
        # The all-defaults path (scored registry == read registry) is quiet.
        labels = _labels(_build({
            "0": {"CD3D": 30, "CD3E": 30, "CD8A": 30, "CD8B": 25, "GZMB": 30},
            "1": {"MS4A1": 40, "CD79A": 30, "CD79B": 30, "CD19": 20},
        }))
        assert labels["0"].startswith("CD8")
        with caplog.at_level(logging.WARNING, logger="tcrsift.annotate_cells"):
            _labels(_build({
                "0": {"CD3D": 30, "CD3E": 30, "CD8A": 30, "CD8B": 25, "GZMB": 30},
                "1": {"MS4A1": 40, "CD79A": 30, "CD79B": 30, "CD19": 20},
            }))
        assert not any("will be ignored" in r.message for r in caplog.records)

    def test_annotate_cells_threads_custom_t_state(self):
        # End to end through the one-call path: a custom T-state registry is both
        # scored (score_reference) and read (annotate_clusters).
        custom_t = {"Cyto2": ["GZMB", "PRF1", "GZMA", "NKG7", "GZMH"]}
        adata = _build({
            "0": {"CD3D": 30, "CD3E": 30, "CD8A": 30, "CD8B": 25,
                  "GZMB": 40, "PRF1": 30, "GZMA": 30, "NKG7": 20, "GZMH": 20},
            "1": {"CD3D": 30, "CD3E": 30, "CD4": 30, "CCR7": 25, "IL7R": 18},
        })
        labels = _labels(adata, t_state_reference=custom_t)
        assert labels["0"] == "CD8 Cyto2"
