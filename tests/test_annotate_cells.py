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

import warnings

import anndata as ad
import numpy as np
import pandas as pd

from tcrsift.annotate_cells import (
    AnnotationGates,
    annotate_cells,
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
