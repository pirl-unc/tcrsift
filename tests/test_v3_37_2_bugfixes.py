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

"""Regression tests for v3.37.2 audit bugfixes."""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.annotate import _annotate_match
from tcrsift.annotate_cells import top_markers
from tcrsift.assemble import MOUSE_TRAC_AA, assemble_full_sequences, optimize_codons
from tcrsift.cli import create_parser
from tcrsift.clonotype import aggregate_clonotypes, calculate_clone_frequencies
from tcrsift.config import AssembleConfig, TCRsiftConfig
from tcrsift.filter import filter_clonotypes_stimulation, filter_clonotypes_threshold
from tcrsift.phenotype import filter_by_tcell_type, get_phenotype_summary
from tcrsift.plots import (
    plot_pgen_distribution,
    plot_rank_curves_per_sample,
    plot_top_n_labeled_bars_per_sample,
    plot_vgene_usage_by_method,
)
from tcrsift.prioritize import (
    PrioritizeConfig,
    resolve_signatures,
    select_round_robin,
)
from tcrsift.selection import select_freq_prism_per_condition
from tcrsift.signature_methods import score_by_name
from tcrsift.til import identify_til_specific_clones
from tcrsift.til_select import match_clonotypes
from tcrsift.validation import (
    TCRsiftValidationError,
    validate_assemble_args,
    validate_run_args,
)


def test_prioritize_config_preserves_allow_canonical_fallback():
    ns = argparse.Namespace(
        sample_sheet="samples.csv",
        output_dir="out",
        allow_canonical_fallback=False,
    )
    cfg = PrioritizeConfig.from_namespace(ns)
    assert cfg.allow_canonical_fallback is False


def test_avoid_ncoi_removed_in_mouse_and_chimeric_and_allele_override_trac():
    # Direct codon optimization of MOUSE_TRAC_AA must eliminate NcoI (CCATGG) across T-M-E.
    nt = optimize_codons(MOUSE_TRAC_AA, species="mouse", forbidden_motifs=("CCATGG",))
    assert "CCATGG" not in nt

    # Mouse native + chimeric human-VDJ/mouse-C + explicit allele override must all honor avoid_enzymes.
    df_mouse = pd.DataFrame(
        [
            {
                "CDR3ab": "CAVRDSNYQLIW_CASSLGGNTEVFF",
                "CDR3_alpha": "CAVRDSNYQLIW",
                "CDR3_beta": "CASSLGGNTEVFF",
                "VDJ_alpha_aa": "AQKVTQAQTEISVVEKEDVTLDCVYETRDTTYYLFWYKQPPSGELVFLIRRNSFDEQNEISGRYSWNFQKSTSSFNFTITASQVVDSAVYFCAVRDSNYQLIWGAGTKLIIKP",
                "VDJ_beta_aa": "DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSLGGNTEVFFGKGTRLTVV",
                "alpha_v_gene": "TRAV12-2",
                "alpha_j_gene": "TRAJ33",
                "alpha_c_gene": "Trac",
                "beta_v_gene": "TRBV13-1",
                "beta_j_gene": "TRBJ1-1",
                "beta_c_gene": "Trbc1",
            }
        ]
    )
    out_mouse = assemble_full_sequences(
        df_mouse,
        species="mouse",
        include_constant=True,
        alpha_leader=None,
        beta_leader=None,
        linker="T2A",
        avoid_enzymes=("NcoI",),
        codon_optimization="full",
        show_progress=False,
    )
    assert "CCATGG" not in str(out_mouse.loc[0, "alpha_constant_nt"]).upper()

    out_chimeric = assemble_full_sequences(
        df_mouse.assign(alpha_c_gene="TRAC", beta_c_gene="TRBC1"),
        species="human",
        constant_species="mouse",
        include_constant=True,
        alpha_leader=None,
        beta_leader=None,
        linker="T2A",
        avoid_enzymes=("NcoI",),
        codon_optimization="full",
        show_progress=False,
    )
    assert "CCATGG" not in str(out_chimeric.loc[0, "alpha_constant_nt"]).upper()

    out_override = assemble_full_sequences(
        df_mouse,
        species="mouse",
        trac_allele="01",
        include_constant=True,
        alpha_leader=None,
        beta_leader=None,
        linker="T2A",
        avoid_enzymes=("NcoI",),
        codon_optimization="full",
        show_progress=False,
    )
    assert "CCATGG" not in str(out_override.loc[0, "alpha_constant_nt"]).upper()


def test_annotate_cells_top_markers_resolves_ensembl_gene_symbols():
    rng = np.random.default_rng(0)
    n_cells = 40
    # Cluster 0 strongly expresses ENSG_TRBV (should be excluded!) and ENSG_GZMK;
    # Cluster 1 strongly expresses ENSG_CCR7.
    X = rng.uniform(0.0, 0.2, size=(n_cells, 4)).astype(np.float32)
    X[:20, 0] = 6.0  # TRBV20-1 -> excluded
    X[:20, 1] = 4.0  # GZMK -> top non-TCR marker for cluster 0
    X[20:, 2] = 5.0  # CCR7 -> top marker for cluster 1
    obs = pd.DataFrame({"leiden": ["0"] * 20 + ["1"] * 20}, index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(
        {"gene_symbols": ["TRBV20-1", "GZMK", "CCR7", "ACTB"]},
        index=["ENSG00000211751", "ENSG00000113088", "ENSG00000126353", "ENSG00000075624"],
    )
    adata = ad.AnnData(X=X, obs=obs, var=var)
    markers = top_markers(adata, "leiden")
    assert markers["0"] == "GZMK"
    assert markers["1"] == "CCR7"


def test_annotate_match_ignores_empty_antigen_gene_mode():
    df = pd.DataFrame([{"CDR3ab": "CA_CB"}])
    matches = pd.DataFrame(
        [
            {"antigen_gene": "", "db_protein_canonical": "", "epitope": "", "species": ""},
            {"antigen_gene": "", "db_protein_canonical": "", "epitope": "", "species": ""},
            {"antigen_gene": "pp65", "db_protein_canonical": "CMV pp65", "epitope": "NLVPMVATV", "species": "CMV"},
        ]
    )
    _annotate_match(df, 0, matches)
    assert df.loc[0, "db_protein_canonical"] == "CMV pp65"
    assert df.loc[0, "db_epitope"] == "NLVPMVATV"


def test_select_round_robin_preserves_unselected_rows_with_duplicate_index():
    clones = pd.DataFrame(
        [
            {"donor": "D1", "CDR3ab": "A1_B1", "eligible_for_review": True},
            {"donor": "D2", "CDR3ab": "A2_B2", "eligible_for_review": False},
        ],
        index=[0, 0],
    )
    scores = pd.DataFrame(
        [
            {
                "donor": "D1",
                "sample": "S1",
                "lineage": "cd8",
                "CDR3ab": "A1_B1",
                "cells": 10,
                "frequency": 0.1,
                "signature_TumorReactive": 2.5,
                "signature_TumorReactive_noise_floor": 0.5,
                "signature_TumorReactive_percentile": 1.0,
            },
            {
                "donor": "D1",
                "sample": "S1",
                "lineage": "cd8",
                "CDR3ab": "A3_B3",
                "cells": 5,
                "frequency": 0.05,
                "signature_TumorReactive": 0.1,
                "signature_TumorReactive_noise_floor": 0.5,
                "signature_TumorReactive_percentile": 0.2,
            },
        ]
    )
    out = select_round_robin(clones, scores, ["TumorReactive"], max_clones=10)
    selected_row = out[out["CDR3ab"] == "A1_B1"].iloc[0]
    unselected_row = out[out["CDR3ab"] == "A2_B2"].iloc[0]
    assert bool(selected_row["selected_for_review"]) is True
    assert bool(unselected_row["selected_for_review"]) is False
    assert pd.isna(unselected_row["selection_rank"])


def test_resolve_signatures_is_case_insensitive_on_tcell_type():
    assert resolve_signatures("solid-tumor", "CD8") == resolve_signatures("solid-tumor", "cd8")
    assert resolve_signatures("solid-tumor", "BOTH") == resolve_signatures("solid-tumor", "both")


def test_score_by_name_does_not_double_log_log_normalized_adata_without_uns():
    X = np.array(
        [
            [0.0, 2.0, 0.5, 1.0],
            [3.0, 0.1, 2.5, 0.2],
            [1.5, 1.2, 0.8, 0.4],
            [4.0, 0.3, 3.2, 0.1],
        ],
        dtype=np.float32,
    )
    var = pd.DataFrame(index=["PRF1", "GZMB", "NKG7", "ACTB"])
    obs = pd.DataFrame(index=[f"c{i}" for i in range(4)])
    adata_no_uns = ad.AnnData(X=X.copy(), obs=obs.copy(), var=var.copy())
    adata_with_uns = ad.AnnData(X=X.copy(), obs=obs.copy(), var=var.copy())
    adata_with_uns.uns["log1p"] = {"base": None}

    s_no_uns = score_by_name(adata_no_uns, "Cytolytic")
    s_with_uns = score_by_name(adata_with_uns, "Cytolytic")
    np.testing.assert_allclose(s_no_uns.to_numpy(), s_with_uns.to_numpy(), rtol=1e-6, atol=1e-6)


def test_cli_config_contig_override_and_sample_sheet_validation(tmp_path):
    # Bug #10: CLI --cellranger-dir overrides YAML contigs_dir without dual-source conflict
    dir_a = tmp_path / "dir_a"
    dir_b = tmp_path / "dir_b"
    dir_a.mkdir()
    dir_b.mkdir()
    cfg = TCRsiftConfig(assemble=AssembleConfig(contigs_dir=str(dir_a)))
    merged = cfg.merge_with_args(argparse.Namespace(cellranger_dir=str(dir_b)))
    assert merged.assemble.cellranger_dir == str(dir_b)
    assert merged.assemble.contigs_dir is None

    # Bug #9: assemble & report selected parsers expose --sample-sheet and validate_assemble_args requires it for sheet mode
    parser = create_parser()
    asm_args = parser.parse_args(
        ["assemble", "-i", "in.csv", "-o", "out.csv", "--sample-name-from", "sheet", "--sample-sheet", "s.yaml"]
    )
    assert asm_args.sample_sheet == "s.yaml"
    rep_args = parser.parse_args(
        ["report", "selected", "--selection", "sel.csv", "--clonotypes", "c.csv", "-o", "out", "--sample-name-from", "sheet", "--sample-sheet", "s.yaml"]
    )
    assert rep_args.sample_sheet == "s.yaml"
    with pytest.raises(TCRsiftValidationError, match="--sample-sheet is required"):
        validate_assemble_args(
            argparse.Namespace(
                no_leaders=False,
                leaders_from_contigs=False,
                alpha_leader="CD8A",
                beta_leader="CD8A",
                contigs_dir=str(dir_a),
                cellranger_dir=None,
                sample_name_from="sheet",
                sample_sheet=None,
            )
        )

    # Bug #11: validate_run_args accepts --leaders-from-contigs when sample sheet provides vdj_dir
    vdj_dir = tmp_path / "S1"
    vdj_dir.mkdir()
    sheet_path = tmp_path / "samples.csv"
    pd.DataFrame(
        [{"sample": "S1", "vdj_dir": str(vdj_dir), "source": "culture"}]
    ).to_csv(sheet_path, index=False)
    validate_run_args(
        argparse.Namespace(
            sample_sheet=str(sheet_path),
            leaders_from_contigs=True,
            no_leaders=False,
            no_assemble=False,
            alpha_leader=None,
            beta_leader=None,
            contigs_dir=None,
            cellranger_dir=None,
            no_sheet_contigs=False,
            til_samples=None,
            til_sample=None,
        )
    )


def test_til_select_match_clonotypes_handles_empty_after_filter_db():
    clones = pd.DataFrame([{"CDR3_alpha": "CAVF", "CDR3_beta": "CASSF"}])
    beta_only_db = pd.DataFrame(
        [{"cdr3_alpha": "", "cdr3_beta": "CASSF", "epitope": "GIL", "species": "Flu", "database": "VDJdb", "is_viral": True}]
    )
    out = match_clonotypes(clones, beta_only_db, match_by="CDR3ab")
    assert list(out["db_match"]) == [False]
    assert list(out["is_viral"]) == [False]


def test_filter_stimulation_duplicate_index_and_threshold_beta_only():
    stim_df = pd.DataFrame(
        [
            {"CDR3ab": "A1_B1", "CDR3_alpha": "A1", "CDR3_beta": "B1", "cell_count": 5, "stim_fdr": 0.001, "stim_log2_fc": 3.0},
            {"CDR3ab": "A2_B2", "CDR3_alpha": "A2", "CDR3_beta": "B2", "cell_count": 5, "stim_fdr": 0.5, "stim_log2_fc": 0.0},
            {"CDR3ab": "A3_B3", "CDR3_alpha": "A3", "CDR3_beta": "B3", "cell_count": 5, "stim_fdr": 0.002, "stim_log2_fc": 2.0},
        ],
        index=[0, 0, 1],
    )
    res = filter_clonotypes_stimulation(stim_df)
    assert len(res) == 2
    assert set(res["CDR3ab"]) == {"A1_B1", "A3_B3"}

    beta_only = pd.DataFrame([{"CDR3ab": "CASSF", "CDR3_beta": "CASSF", "cell_count": 5}])
    with pytest.raises(TCRsiftValidationError, match="require_complete=True requires both"):
        filter_clonotypes_threshold(beta_only, require_complete=True)


def test_clonotype_all_doublets_and_empty_frequencies():
    adata = ad.AnnData(
        X=np.zeros((2, 2)),
        obs=pd.DataFrame(
            {
                "sample": ["S1", "S1"],
                "CDR3_alpha": ["CAVF", "CAVF"],
                "CDR3_beta": ["CASSF", "CASSF"],
                "multi_chain": [True, True],
            },
            index=["c0", "c1"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    with pytest.raises(TCRsiftValidationError, match="No complete clones found"):
        aggregate_clonotypes(adata, handle_doublets="remove", verbose=True)

    empty_clones = pd.DataFrame(columns=["CDR3ab", "CDR3_alpha", "CDR3_beta", "cell_count"])
    freq_out = calculate_clone_frequencies(empty_clones, adata)
    assert "sample_frequencies" in freq_out.columns
    assert len(freq_out) == 0


def test_select_freq_prism_skips_nan_clone_ids():
    feat = pd.DataFrame(
        [
            {"CDR3ab": np.nan, "condition": "C1", "frequency": 0.5},
            {"CDR3ab": "A1_B1", "condition": "C1", "frequency": 0.2},
        ]
    )
    out = select_freq_prism_per_condition(
        feat, condition_col="condition", gate=0.01, top_freq=5, top_prism=0
    )
    assert list(out["CDR3ab"]) == ["A1_B1"]


def test_phenotype_empty_adata_and_til_without_sample_column():
    empty_adata = ad.AnnData(
        X=np.zeros((0, 2)),
        obs=pd.DataFrame(
            {
                "Tcell_type": pd.Series([], dtype="object"),
                "expected_tcell_type": pd.Series([], dtype="object"),
                "is_CD8": pd.Series([], dtype="bool"),
                "is_CD4": pd.Series([], dtype="bool"),
            },
            index=pd.Index([], dtype=str),
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    summary = get_phenotype_summary(empty_adata)
    assert len(summary) == 1
    assert summary.loc[0, "total_cells"] == 0
    assert pd.isna(summary.loc[0, "expected_type"])

    with pytest.raises(TCRsiftValidationError, match="No cells remain"):
        filter_by_tcell_type(empty_adata, "cd8", verbose=True)

    til_adata = ad.AnnData(
        X=np.zeros((2, 2)),
        obs=pd.DataFrame(
            {"CDR3_alpha": ["CAVF", "CAVF"], "CDR3_beta": ["CASSF", "CASSF"]},
            index=["t0", "t1"],
        ),
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    til_clones = identify_til_specific_clones(til_adata, min_cells=2)
    assert list(til_clones["CDR3ab"]) == ["CAVF_CASSF"]
    assert list(til_clones["sample"]) == ["unknown"]


def test_plots_publicness_duplicate_cdr3ab_and_missing_cdr3ab(tmp_path):
    # Bug #20: plot_pgen_distribution colors bins by publicness_col
    df_pgen = pd.DataFrame(
        {
            "log10_pgen": [-32.0, -31.0, -25.0, -24.0, -16.0, -15.0],
            "publicness": [0.05, 0.1, 0.5, 0.55, 0.95, 0.98],
        }
    )
    p1 = plot_pgen_distribution(df_pgen, tmp_path / "pgen.png")
    assert p1 is not None and p1.exists()

    # Bug #21: duplicate CDR3ab rows in ranks plots do not crash
    dup_clones = pd.DataFrame(
        [
            {"CDR3ab": "A1_B1", "sample_frequencies": {"S1": 0.2}, "tier": "tier1", "cell_count": 10},
            {"CDR3ab": "A1_B1", "sample_frequencies": {"S1": 0.2}, "tier": "tier2", "cell_count": 5},
            {"CDR3ab": "A2_B2", "sample_frequencies": {"S1": 0.1}, "tier": None, "cell_count": 3},
        ]
    )
    assert plot_rank_curves_per_sample(dup_clones, tmp_path) is not None
    assert plot_top_n_labeled_bars_per_sample(dup_clones, tmp_path) is not None

    # Bug #22: plot_vgene_usage_by_method returns None when CDR3ab is absent
    long_df = pd.DataFrame([{"method": "M1", "CDR3_beta": "CASSF"}])
    clo_df = pd.DataFrame([{"CDR3_beta": "CASSF", "beta_v_gene": "TRBV20-1"}])
    assert plot_vgene_usage_by_method(long_df, clo_df, tmp_path / "vgene.png") is None


def test_pipeline_forwards_clone_scores_to_selection():
    """Verify run_pipeline passes clone_scores=til_matched to select_from_clone_sample_long."""
    import inspect

    from tcrsift import pipeline

    src = inspect.getsource(pipeline.run_pipeline)
    assert "clone_scores=til_matched" in src
    assert '_sample_name_from in ("parent", "sheet")' in src
