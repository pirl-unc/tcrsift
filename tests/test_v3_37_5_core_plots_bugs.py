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

"""Regression tests for v3.37.5 core I/O, QC, GEX, multimer, report, CLI, and plot fixes."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from tcrsift._dtypes import rehydrate_obs
from tcrsift.attribution import attribute_cells
from tcrsift.config import AttributionConfig
from tcrsift.datacache import DownloadError, _extract_zip
from tcrsift.genes import resolve_adata_var_map
from tcrsift.gex import (
    augment_with_gex,
    compute_cd4_cd8_counts,
    compute_signature_scores_per_cell,
)
from tcrsift.multimer import call_multimers_for_adata, classify_feature_columns
from tcrsift.plots import (
    plot_category_composition_bar,
    plot_clone_tracking_slopes,
    plot_filter,
    plot_matched_clone_heatmap,
    plot_provenance,
)
from tcrsift.qc import cross_lineage_doublets, validate_clonotypes
from tcrsift.report import bundle_figure_pdf, expand_dual_alpha_variants
from tcrsift.sct import aggregate_sct
from tcrsift.til_select import (
    dominant_specificity_from_master,
    plot_annotated_heatmap,
    select_top_monotonic_increase_clones,
)


def test_aggregate_sct_boolean_normalization_does_not_invert():
    df = pd.DataFrame(
        [
            {"CDR3_pair": "A_B", "flag": "yes"},
            {"CDR3_pair": "A_B", "flag": "no"},
            {"CDR3_pair": "C_D", "flag": "True"},
            {"CDR3_pair": "C_D", "flag": "False"},
            {"CDR3_pair": "E_F", "flag": None},
        ]
    )
    out = aggregate_sct(df, group_cols=["CDR3_pair"], numeric_cols=[], boolean_cols=["flag"])
    out = out.set_index("CDR3_pair")
    assert bool(out.loc["A_B", "flag.any"]) is True
    assert bool(out.loc["A_B", "flag.all"]) is False
    assert bool(out.loc["C_D", "flag.any"]) is True
    assert bool(out.loc["C_D", "flag.all"]) is False
    assert bool(out.loc["E_F", "flag.any"]) is False
    assert bool(out.loc["E_F", "flag.all"]) is False


def test_gex_cp10k_normalizes_low_count_signatures_and_resolves_synonym_pairs_and_empty_counts(tmp_path):
    # Two raw-count cells with identical PDCD1=10 (<=25) but 10x different total library depth
    X = np.array(
        [
            [10.0, 990.0],
            [10.0, 9990.0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(index=["PDCD1", "ACTB"]),
    )
    compute_signature_scores_per_cell(adata, signatures={"Exh": ["PDCD1"]})
    scores = adata.obs["signature_Exh"].to_numpy()
    assert scores[0] > scores[1] + 1.5  # CP10K depth-normalized

    # Synonym pair in gene_list with only mouse Cd8b1 in matrix
    adata_syn = ad.AnnData(
        X=np.array([[5.0, 1.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0"]),
        var=pd.DataFrame(index=["Cd8b1", "Actb"]),
    )
    gex_h5ad = tmp_path / "syn.h5ad"
    adata_syn.write_h5ad(gex_h5ad)
    df_cells = pd.DataFrame([{"barcode": "c0", "CDR3_pair": "A_B"}])
    aug = augment_with_gex(df_cells, gex_h5ad, gene_list=["CD8B", "CD8B1"])
    assert "gex.CD8B" in aug.columns and aug.loc[0, "gex.CD8B"] == 5.0

    # Empty DataFrame in compute_cd4_cd8_counts
    empty_res = compute_cd4_cd8_counts(pd.DataFrame(columns=["CDR3_pair", "gex.CD4", "gex.CD8A", "gex.CD8B"]), verbose=True)
    assert list(empty_res.columns) == ["CDR3_pair", "total_cells.count", "CD4_only.count", "CD8_only.count"]


def test_resolve_adata_var_map_handles_zero_cell_sparse_and_nan_duplicates():
    empty_sparse = ad.AnnData(
        X=csr_matrix((0, 2), dtype=np.float32),
        obs=pd.DataFrame(index=pd.Index([], dtype=str)),
        var=pd.DataFrame({"gene_symbols": ["CD8A", "CD8A"]}, index=["ENSG1", "ENSG2"]),
    )
    _, sym_to_var = resolve_adata_var_map(empty_sparse, ["CD8A"])
    assert sym_to_var["CD8A"] == "ENSG1"

    nan_dense = ad.AnnData(
        X=np.array([[np.nan, 5.0], [1.0, 5.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame({"gene_symbols": ["CD8A", "CD8A"]}, index=["ENSG1", "ENSG2"]),
    )
    _, sym_to_var2 = resolve_adata_var_map(nan_dense, ["CD8A"])
    assert sym_to_var2["CD8A"] == "ENSG2"


def test_qc_and_attribution_handle_secondary_chain_nans_categoricals_and_small_panels():
    df = pd.DataFrame(
        {
            "CDR3_alpha": pd.Categorical(["CAVRDSNYQLIW"]),
            "CDR3_beta": pd.Categorical(["CASSLGGNTEVFF"]),
            "TRA_1_reads": [100.0],
            "TRA_2_reads": [np.nan],
            "TRA_1_umis": [10.0],
            "TRA_2_umis": [np.nan],
        }
    )
    validated, _ = validate_clonotypes(df)
    assert bool(validated.loc[0, "qc_pass"]) is True

    # Attribution with Categorical sample column
    obs_cat = pd.DataFrame(
        {
            "sample": pd.Categorical(["S1", "S1"]),
            "CDR3_alpha": ["CAVF", "CAVF"],
            "CDR3_beta": ["CASSF", "CASSF"],
        },
        index=["c0", "c1"],
    )
    long_df, _, _ = attribute_cells(obs_cat, AttributionConfig())
    assert len(long_df) == 2

    # cross_lineage_doublets on tiny 3-gene panel does not crash
    tiny = ad.AnnData(
        X=np.array([[5.0, 0.0, 1.0], [0.0, 5.0, 1.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["c0", "c1"]),
        var=pd.DataFrame(index=["CD3E", "MS4A1", "ACTB"]),
    )
    mask = cross_lineage_doublets(tiny)
    assert mask.shape == (2,)


def test_multimer_nc_regex_and_sample_prefixed_barcodes(monkeypatch):
    fb = pd.DataFrame(
        {
            "NC_1_TotalSeqC": [1.0, 1.0],
            "dCODE_NC_1_PE": [1.0, 1.0],
            "NC_ctrl": [1.0, 1.0],
            "NC_045512_Spike": [10.0, 0.0],
            "GILGFVFTL_Flu": [80.0, 2.0],
        },
        index=["ACGT-1", "TGCA-1"],
    )
    cls = classify_feature_columns(fb)
    assert set(cls["negative_controls"]) == {"NC_1_TotalSeqC", "dCODE_NC_1_PE", "NC_ctrl"}
    assert "NC_045512_Spike" in cls["multimers"]

    # Sample-prefixed barcodes in AnnData align to raw CellRanger barcodes
    adata = ad.AnnData(
        X=np.zeros((2, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"sample": ["S1", "S1"], "CDR3ab": ["A1_B1", "A1_B1"]},
            index=["S1_ACGT-1", "S1_TGCA-1"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )

    class _Sample:
        sample = "S1"
        gex_dir = "dummy"

    import tcrsift.multimer as mm

    monkeypatch.setattr(mm, "extract_feature_barcodes", lambda _p: fb.copy())
    summary = call_multimers_for_adata(
        adata, [_Sample()], min_umi=10.0, min_nc_ratio=5.0, min_specific_cells=1
    )
    assert len(summary) == 1
    assert summary.loc[0, "multimer_call"] == "GILGFVFTL_Flu"
    assert int(summary.loc[0, "multimer_call_cells"]) == 1


def test_report_datacache_and_dtypes_fixes(tmp_path):
    # bundle_figure_pdf raises FileNotFoundError on nonexistent extra_images
    with pytest.raises(FileNotFoundError):
        bundle_figure_pdf(
            [],
            tmp_path / "out.pdf",
            extra_images=[("Header", tmp_path / "missing.png")],
        )

    # expand_dual_alpha_variants includes selected_clone even when obs is None
    sel = pd.DataFrame([{"CDR3ab": "A_B", "CDR3_beta": "B"}])
    exp, _ = expand_dual_alpha_variants(sel, obs=None)
    assert "selected_clone" in exp.columns
    assert exp.loc[0, "selected_clone"] == "A_B"

    # datacache._extract_zip converts BadZipFile into DownloadError
    bad_zip = tmp_path / "corrupt.zip"
    bad_zip.write_text("<html>503 Service Unavailable</html>")
    with pytest.raises(DownloadError, match="not a valid zip"):
        _extract_zip(bad_zip, tmp_path / "dest", None)

    # rehydrate_obs makes mixed-type object columns writable to h5ad
    adata = ad.AnnData(
        X=np.zeros((3, 2), dtype=np.float32),
        obs=pd.DataFrame({"mixed": pd.Series([1, "foo", None], dtype=object)}, index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    rehydrate_obs(adata)
    adata.write_h5ad(tmp_path / "mixed.h5ad")


def test_til_select_monotonic_lineage_and_annotated_heatmap(tmp_path):
    df_mono = pd.DataFrame(
        [
            {"CDR3ab": "A1_B1", "frequency_T1": 0.0, "frequency_T2": 0.05, "frequency_T3": 0.1},
            {"CDR3ab": "A2_B2", "frequency_T1": 0.01, "frequency_T2": 0.03, "frequency_T3": 0.06},
        ]
    )
    sel_all, _ = select_top_monotonic_increase_clones(
        df_mono, ["T1", "T2", "T3"], top_k=10, require_positive_all_timepoints=True
    )
    assert list(sel_all["CDR3ab"]) == ["A2_B2"]

    # Double-negative (0 CD8, 0 CD4) clone must NOT be classified as CD8
    master = pd.DataFrame(
        [
            {
                "CDR3ab": "A1_B1",
                "CDR3_beta": "B1",
                "cell_count_T1": 3,
                "cell_count_T2": 0,
                "cd8_cp10k_max": 0.0,
                "cd4_cp10k_max": 0.0,
                "cd4_to_cd8_ratio": 0.0,
            }
        ]
    )
    dom = dominant_specificity_from_master(
        master,
        timepoint_order=["T1", "T2"],
        per_lineage_min={"CD8": 2, "CD4": 5},
        min_dominant_cells=5,
    )
    assert bool(dom.loc[0, "is_dominant_specific"]) is False

    # plot_annotated_heatmap runs and renders annotation tags
    df_hm = pd.DataFrame(
        [{"CDR3ab": "A1_B1", "cell_count_T1": 5, "cell_count_T2": 10, "is_viral": True, "db_vdjdb": True}]
    )
    hm_path = tmp_path / "annotated_hm.png"
    plot_annotated_heatmap(
        df_hm, ["cell_count_T1", "cell_count_T2"], hm_path, top_k=5,
        annotation_cols=["db_vdjdb", "is_viral"], title="Test",
    )
    assert hm_path.exists()


def test_plots_duplicate_cdr3ab_empty_filter_and_na_provenance(tmp_path):
    dup = pd.DataFrame(
        [
            {"CDR3ab": "A1_B1", "sample_frequencies": {"S1": 0.2, "S2": 0.1}, "db_match": True, "db_category": "viral", "max_frequency": 0.2},
            {"CDR3ab": "A1_B1", "sample_frequencies": {}, "db_match": True, "db_category": "viral", "max_frequency": 0.0},
            {"CDR3ab": "A2_B2", "sample_frequencies": {"S1": 0.1, "S2": 0.3}, "db_match": True, "db_category": "tumor", "max_frequency": 0.3},
        ]
    )
    assert plot_clone_tracking_slopes(dup, tmp_path) is not None
    assert plot_matched_clone_heatmap(dup, tmp_path) is not None
    assert plot_category_composition_bar(dup, tmp_path) is not None

    # Empty & missing cell_count in plot_filter
    plot_filter(pd.DataFrame(columns=["max_frequency", "n_samples"]), tmp_path / "empty_f")
    plot_filter(pd.DataFrame([{"max_frequency": 0.1, "n_samples": 2}]), tmp_path / "nocell_f")
    assert (tmp_path / "nocell_f" / "filter_frequency_specificity.png").exists()

    # plot_provenance does not leak '<NA>' on StringDtype columns
    adata = ad.AnnData(
        X=np.zeros((3, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "flag": [True, True, False],
                "antigen": pd.Series(["Flu", pd.NA, pd.NA], dtype="string"),
            },
            index=["c0", "c1", "c2"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    adata.obsm["X_umap"] = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.5]])
    out_p = plot_provenance(adata, flag="flag", by="antigen", output_path=tmp_path / "prov.png")
    assert out_p.exists()
