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

"""Tests for TIL timepoint selection workflow."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from tcrsift.til_select import (
    compute_marker_scores_for_timepoint,
    default_timepoint_inputs,
    load_from_consensus,
    parse_config,
    parse_sample_args,
    run_selection_pipeline,
    run_til_select,
)


def _write_timepoint_inputs(data_dir: Path, tp: str, clonotype_counts: dict[str, int]) -> None:
    """Write minimal v2-compatible input files for one timepoint."""
    consensus_rows = []
    for clonotype_id in clonotype_counts:
        clone_num = clonotype_id.replace("clonotype", "")
        consensus_rows.append(
            {
                "clonotype_id": clonotype_id,
                "consensus_id": f"{clonotype_id}_consensus1",
                "chain": "TRA",
                "productive": True,
                "full_length": True,
                "cdr3": f"CAV{clone_num}TEST",
                "fwr1": "AAA",
                "cdr1": "BBB",
                "fwr2": "CCC",
                "cdr2": "DDD",
                "fwr3": "EEE",
                "fwr4": "FFF",
                "fwr1_nt": "aaa",
                "cdr1_nt": "bbb",
                "fwr2_nt": "ccc",
                "cdr2_nt": "ddd",
                "fwr3_nt": "eee",
                "cdr3_nt": "ggg",
                "fwr4_nt": "fff",
                "v_gene": "TRAV12-2",
                "j_gene": "TRAJ33",
                "c_gene": "TRAC",
            }
        )
        consensus_rows.append(
            {
                "clonotype_id": clonotype_id,
                "consensus_id": f"{clonotype_id}_consensus2",
                "chain": "TRB",
                "productive": True,
                "full_length": True,
                "cdr3": f"CASS{clone_num}TEST",
                "fwr1": "AAA",
                "cdr1": "BBB",
                "fwr2": "CCC",
                "cdr2": "DDD",
                "fwr3": "EEE",
                "fwr4": "FFF",
                "fwr1_nt": "aaa",
                "cdr1_nt": "bbb",
                "fwr2_nt": "ccc",
                "cdr2_nt": "ddd",
                "fwr3_nt": "eee",
                "cdr3_nt": "ggg",
                "fwr4_nt": "fff",
                "v_gene": "TRBV28",
                "d_gene": "TRBD1",
                "j_gene": "TRBJ2-7",
                "c_gene": "TRBC1",
            }
        )

    consensus_df = pd.DataFrame(consensus_rows)
    consensus_df.to_csv(data_dir / f"consensus_annotations.{tp}.csv", index=False)

    clonotypes_df = pd.DataFrame(
        {
            "clonotype_id": list(clonotype_counts.keys()),
            "frequency": list(clonotype_counts.values()),
            "proportion": np.array(list(clonotype_counts.values())) / max(sum(clonotype_counts.values()), 1),
        }
    )
    clonotypes_df.to_csv(data_dir / f"clonotypes.{tp}.csv", index=False)

    contig_rows = []
    for clonotype_id, n_cells in clonotype_counts.items():
        for i in range(n_cells):
            barcode = f"{tp}_{clonotype_id}_{i:02d}"
            contig_rows.append(
                {
                    "barcode": barcode,
                    "raw_clonotype_id": clonotype_id,
                    "chain": "TRA",
                    "productive": True,
                    "high_confidence": True,
                }
            )
    pd.DataFrame(contig_rows).to_csv(data_dir / f"filtered_contig_annotations.{tp}.csv", index=False)

    # Placeholder path; tests monkeypatch scanpy.read_10x_h5
    (data_dir / f"sample_filtered_feature_bc_matrix.{tp}.h5").write_bytes(b"placeholder")


def _make_mock_adata(barcodes: list[str]) -> ad.AnnData:
    """Create small AnnData for marker scoring tests."""
    genes = ["CD4", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "MKI67", "TNFRSF9"]
    X = []
    for i, _bc in enumerate(barcodes):
        # Alternate higher CD8 and low CD4 to ensure at least some base-selected clones
        if i % 2 == 0:
            X.append([0.0, 8.0, 7.0, 2.0, 2.0, 1.0, 1.0, 1.0])
        else:
            X.append([2.0, 1.0, 1.0, 0.5, 0.5, 0.2, 0.2, 0.2])
    adata = ad.AnnData(sp.csr_matrix(np.array(X, dtype=float)))
    adata.var_names = genes
    adata.obs_names = barcodes
    return adata


class TestTilSelectParsing:
    """Tests for til-select input parsing helpers."""

    def test_parse_sample_args(self, tmp_path):
        """Labelled sample specs should parse as LABEL=consensus,clonotypes."""
        cons = tmp_path / "cons.csv"
        clono = tmp_path / "clono.csv"
        cons.write_text("clonotype_id,chain,cdr3\n")
        clono.write_text("clonotype_id,frequency\n")

        parsed = parse_sample_args([f"T1={cons},{clono}"])
        assert list(parsed.keys()) == ["T1"]
        assert parsed["T1"]["consensus"] == cons
        assert parsed["T1"]["clonotypes"] == clono

    def test_parse_config_with_timepoints(self, tmp_path):
        """YAML config with timepoints mapping should parse and resolve paths."""
        cons = tmp_path / "cons.csv"
        clono = tmp_path / "clono.csv"
        cons.write_text("clonotype_id,chain,cdr3\n")
        clono.write_text("clonotype_id,frequency\n")
        cfg = tmp_path / "config.yaml"
        cfg.write_text(
            "timepoints:\n"
            "  T1:\n"
            "    consensus: cons.csv\n"
            "    clonotypes: clono.csv\n"
        )

        parsed = parse_config(cfg)
        assert list(parsed.keys()) == ["T1"]
        assert parsed["T1"]["consensus"] == cons
        assert parsed["T1"]["clonotypes"] == clono

    def test_default_timepoint_inputs(self, tmp_path):
        """Auto-discovery should detect paired T1/T2 consensus+clonotypes files."""
        for tp in ["T1", "T2"]:
            (tmp_path / f"consensus_annotations.{tp}.csv").write_text("clonotype_id,chain,cdr3\n")
            (tmp_path / f"clonotypes.{tp}.csv").write_text("clonotype_id,frequency\n")

        parsed = default_timepoint_inputs(tmp_path)
        assert parsed is not None
        assert list(parsed.keys()) == ["T1", "T2"]


class TestTilSelectCore:
    """Tests for core harmonization/scoring/selection functions."""

    def test_load_from_consensus(self, tmp_path):
        """Consensus + clonotypes should produce paired CDR3ab with cell counts."""
        _write_timepoint_inputs(tmp_path, "T1", {"clonotype1": 3, "clonotype2": 1})
        df, _stats = load_from_consensus(
            tmp_path / "consensus_annotations.T1.csv",
            tmp_path / "clonotypes.T1.csv",
        )

        assert {"CDR3_alpha", "CDR3_beta", "CDR3ab", "cell_count"}.issubset(df.columns)
        assert len(df) == 2
        assert int(df["cell_count"].sum()) == 4

    def test_compute_marker_scores_for_timepoint(self, tmp_path, monkeypatch):
        """Marker scoring should emit per-cell and per-clonotype score tables."""
        import scanpy as sc

        _write_timepoint_inputs(tmp_path, "T1", {"clonotype1": 2, "clonotype2": 2})
        contig_df = pd.read_csv(tmp_path / "filtered_contig_annotations.T1.csv")
        adata = _make_mock_adata(contig_df["barcode"].tolist())
        monkeypatch.setattr(sc, "read_10x_h5", lambda _p: adata)

        cell_df, score_df = compute_marker_scores_for_timepoint(
            consensus_path=tmp_path / "consensus_annotations.T1.csv",
            contig_csv_path=tmp_path / "filtered_contig_annotations.T1.csv",
            gex_h5_path=tmp_path / "sample_filtered_feature_bc_matrix.T1.h5",
            marker_genes=["CD4", "CD8A", "CD8B", "GZMB"],
        )

        assert len(cell_df) > 0
        assert len(score_df) > 0
        assert "CDR3ab" in score_df.columns
        assert "score_cp10k_CD8A" in score_df.columns
        assert "marker_score_cp10k" in score_df.columns

    def test_run_selection_pipeline(self):
        """Selection pipeline should create mask columns and candidate subset."""
        df = pd.DataFrame(
            {
                "CDR3ab": ["A_B", "C_D", "E_F"],
                "frequency_T1": [0.01, 0.02, 0.01],
                "frequency_T2": [0.03, 0.005, 0.01],
                "total_cells": [10, 3, 1],
                "score_cp10k_CD8A_mean": [2.0, 0.2, 0.0],
                "score_cp10k_CD8B_mean": [2.0, 0.1, 0.0],
                "score_cp10k_CD4_mean": [0.1, 1.0, 0.0],
                "score_cp10k_GZMB_mean": [1.5, 0.1, 0.0],
                "score_cp10k_PRF1_mean": [1.5, 0.1, 0.0],
                "score_cp10k_IFNG_mean": [0.8, 0.0, 0.0],
                "score_cp10k_MKI67_mean": [0.5, 0.0, 0.0],
                "score_cp10k_TNFRSF9_mean": [0.4, 0.0, 0.0],
                "score_z_GZMB_T1": [0.5, -0.2, 0.0],
                "score_z_GZMB_T2": [1.1, -0.3, 0.0],
                "score_z_PRF1_T1": [0.6, -0.2, 0.0],
                "score_z_PRF1_T2": [1.0, -0.3, 0.0],
            }
        )

        out_df, subsets, _seq, _ind, _branch = run_selection_pipeline(
            df,
            timepoint_order=["T1", "T2"],
            marker_genes=["CD4", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "MKI67", "TNFRSF9"],
        )

        assert "is_base_selected" in out_df.columns
        assert "is_candidate_tumor_reactive" in out_df.columns
        assert "subset_candidate_tumor_reactive" in subsets


class TestTilSelectEndToEnd:
    """Small end-to-end test for file outputs and compatibility paths."""

    def test_run_til_select_writes_v2_style_outputs(self, tmp_path, monkeypatch):
        """run_til_select should write core tables/subsets with v2-style names."""
        import scanpy as sc

        data_dir = tmp_path / "data"
        data_dir.mkdir()
        _write_timepoint_inputs(data_dir, "T1", {"clonotype1": 2, "clonotype2": 1})
        _write_timepoint_inputs(data_dir, "T2", {"clonotype1": 3, "clonotype2": 1})

        # Use one shared mock AnnData per requested path
        by_path = {}
        for tp in ["T1", "T2"]:
            contig_df = pd.read_csv(data_dir / f"filtered_contig_annotations.{tp}.csv")
            by_path[str(data_dir / f"sample_filtered_feature_bc_matrix.{tp}.h5")] = _make_mock_adata(
                contig_df["barcode"].tolist()
            )

        monkeypatch.setattr(sc, "read_10x_h5", lambda p: by_path[str(p)])

        fig_dir = tmp_path / "figures"
        out_table = tmp_path / "abTCR_harmonized.csv"
        out_report = fig_dir / "selected_clones_report.pdf"

        args = SimpleNamespace(
            config=None,
            data_dir=data_dir,
            samples=None,
            count_column=None,
            verbose=False,
            top_k=10,
            min_cells_per_clone=1,
            min_cd8_cp10k=0.0,
            max_cd4_to_cd8_ratio=2.0,
            increase_ratio_nonzero_min=1.1,
            increase_ratio_all_timepoints_min=1.1,
            immunogenic_percentile=0.5,
            immunogenic_percentile_slack_frac=0.0,
            immunogenic_min_cp10k=0.0,
            immunogenic_require_above_median=False,
            cytotoxic_last_min_z=0.0,
            cytotoxic_last_min_cp10k=0.0,
            became_cytotoxic_min_delta_z=0.0,
            trend_increase_ratio_min=1.2,
            trend_decrease_ratio_max=0.8,
            marker_genes="CD4,CD8A,CD8B,GZMB,PRF1,IFNG,MKI67,TNFRSF9",
            immunogenic_genes="GZMB,PRF1,IFNG,MKI67,TNFRSF9",
            cytotoxic_genes="GZMB,PRF1,IFNG,MKI67,TNFRSF9",
            cytolytic_genes="GZMB,PRF1",
            antigen_response_genes="TNFRSF9,MKI67",
            pyensembl_release=110,
            rank_by="mean_frequency",
            fig_dir=fig_dir,
            out_table=out_table,
            out_heatmap=None,
            out_top=None,
            vdjdb=None,
            iedb=None,
            cedar=None,
            match_by="CDR3b_only",
            out_annotated=None,
            out_annotated_heatmap=None,
            out_selected_report=out_report,
        )

        run_til_select(args)

        assert out_table.exists()
        assert (fig_dir / "abTCR_master_table.csv").exists()
        assert (fig_dir / "abTCR_annotated.csv").exists()
        assert (fig_dir / "selection_masks.csv").exists()
        assert (fig_dir / "subset_candidate_tumor_reactive.csv").exists()
        assert (fig_dir / "marker_cells_T1.csv").exists()
        assert (fig_dir / "marker_clonotype_scores_T1.csv").exists()
        assert (fig_dir / "abTCR_topk.csv").exists()
        assert out_report.exists()
