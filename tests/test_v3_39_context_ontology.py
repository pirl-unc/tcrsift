"""Tests for v3.39.0 biological-context ontology, new contexts, and CD4 TIL Treg warning."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import tcrsift
from tcrsift.cli import create_parser
from tcrsift.prioritize import (
    CONTEXT_ONTOLOGY,
    CONTEXT_SIGNATURES,
    CONTEXT_SIGNATURES_MOUSE,
    PRIMARY_CONTEXTS,
    PRIORITY_CONTEXT_SIGNATURES,
    ContextSpec,
    resolve_signatures,
)


def test_context_ontology_covers_all_presets_and_species():
    assert set(CONTEXT_ONTOLOGY) == set(CONTEXT_SIGNATURES)
    assert set(CONTEXT_SIGNATURES_MOUSE) == set(CONTEXT_SIGNATURES)
    assert set(PRIORITY_CONTEXT_SIGNATURES) == set(CONTEXT_SIGNATURES)
    assert tcrsift.CONTEXT_ONTOLOGY is CONTEXT_ONTOLOGY
    assert tcrsift.PRIMARY_CONTEXTS == PRIMARY_CONTEXTS
    assert tcrsift.ContextSpec is ContextSpec
    assert PRIMARY_CONTEXTS == (
        "solid-tumor",
        "blood-tumor",
        "blood",
        "peptide-stim",
        "generic",
    )
    assert CONTEXT_ONTOLOGY["heme"].primary is False
    assert CONTEXT_ONTOLOGY["heme"].alias_of == "blood"
    assert CONTEXT_ONTOLOGY["mpe"].primary is False
    for name in PRIMARY_CONTEXTS:
        assert CONTEXT_ONTOLOGY[name].primary is True
    for name, spec in CONTEXT_ONTOLOGY.items():
        assert spec.name == name
        assert spec.category and spec.compartment and spec.antigen_exposure and spec.summary

    # peptide-stim prioritizes acute activation/AIM without chronic tumor exhaustion
    stim_cd8 = resolve_signatures("peptide-stim", "cd8")
    assert "TumorReactive" not in stim_cd8
    assert stim_cd8[:4] == ["AIM", "AcuteActivation", "Proliferation", "Cytolytic"]
    assert PRIORITY_CONTEXT_SIGNATURES["peptide-stim"] == (
        "AIM",
        "AcuteActivation",
        "Proliferation",
        "Cytolytic",
    )


def test_cli_accepts_new_contexts():
    parser = create_parser()
    for ctx in ("peptide-stim", "solid-tumor", "blood-tumor", "blood", "mpe", "heme", "generic"):
        args = parser.parse_args(["prioritize", "samples.yaml", "-o", "out", "--context", ctx])
        assert args.context == ctx


def test_cd4_solid_tumor_warns_when_tregs_present_and_exclude_treg_disabled(tmp_path, monkeypatch, caplog):
    from tcrsift import til_prioritize as tp

    genes = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4", "TNFRSF9", "CXCL13", "ENTPD1", "ITGAE", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "GZMB", "PRF1", "NKG7", "IFNG", "MKI67", "TOP2A", "CCL3", "CCL4", "CCL5", "TCF7", "CCR7", "SELL", "IL7R", "LEF1", "GZMK", "CD69", "NR4A1", "FOS", "JUN", "CD40LG", "ICOS", "TNFRSF4", "TNFRSF18", "BATF", "PRDM1", "TOX", "EOMES", "TBX21", "KLRG1", "CX3CR1", "S1PR1", "KLF2", "MALAT1", "ACTB"]
    # 6 cells: 3 FOXP3+ Treg clone, 3 FOXP3- effector clone
    n_cells = 6
    X = np.full((n_cells, len(genes)), 2.0, dtype=np.float32)
    gidx = {g: i for i, g in enumerate(genes)}
    X[:, gidx["CD3D"]] = 10.0
    X[:, gidx["CD4"]] = 10.0
    X[:, gidx["CD8A"]] = 0.0
    X[:, gidx["CD8B"]] = 0.0
    X[:3, gidx["FOXP3"]] = 12.0
    X[3:, gidx["FOXP3"]] = 0.0
    obs = pd.DataFrame(
        {
            "sample": ["S1"] * n_cells,
            "patient_id": ["P1"] * n_cells,
            "CDR3_alpha": ["CAVTREG"] * 3 + ["CAVEFF"] * 3,
            "CDR3_beta": ["CASSTREG"] * 3 + ["CASSEFF"] * 3,
            "CDR3ab": ["CAVTREG_CASSTREG"] * 3 + ["CAVEFF_CASSEFF"] * 3,
            "TRA_1_umis": [10] * n_cells,
            "TRB_1_umis": [10] * n_cells,
            "TRA_1_v_gene": ["TRAV1-2"] * n_cells,
            "TRB_1_v_gene": ["TRBV7-2"] * n_cells,
            "TRA_1_j_gene": ["TRAJ33"] * n_cells,
            "TRB_1_j_gene": ["TRBJ2-1"] * n_cells,
            "n_genes": [500] * n_cells,
            "total_counts": [2000] * n_cells,
            "pct_mito": [1.0] * n_cells,
            "CD3": [10.0] * n_cells,
            "CD4": [10.0] * n_cells,
            "CD8": [0.0] * n_cells,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))

    class _DummySheet:
        species = "human"
        samples = []

        def __iter__(self):
            return iter([])

    monkeypatch.setattr("tcrsift.sample_sheet.load_sample_sheet", lambda _p: _DummySheet())
    monkeypatch.setattr("tcrsift.loader.load_samples", lambda _s, **_kw: adata.copy())

    sheet_path = tmp_path / "sheet.yaml"
    sheet_path.write_text("samples: []\n")
    out_dir = tmp_path / "out"
    parser = create_parser()
    args = parser.parse_args(
        [
            "prioritize",
            str(sheet_path),
            "-o",
            str(out_dir),
            "--context",
            "solid-tumor",
            "--tcell-type",
            "cd4",
            "--signatures",
            "AIM",
            "--signature-cutoff",
            "legacy",
            "--min-signature-score",
            "-10",
            "--signature-quantile",
            "0",
            "--min-cd3",
            "1",
            "--no-exclude-treg",
        ]
    )
    with caplog.at_level(logging.WARNING, logger="tcrsift.til_prioritize"):
        tp.run_til_prioritize(args)

    assert any("FOXP3+ Treg-like" in rec.message and "--exclude-treg" in rec.message for rec in caplog.records)
    assert (Path(out_dir) / "candidate_clones.csv").exists()
