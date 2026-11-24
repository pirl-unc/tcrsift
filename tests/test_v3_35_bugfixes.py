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

"""Regression tests for the v3.35.0 library-wide bug fixes."""

from __future__ import annotations

import argparse
import inspect
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Batch 1: Loading, QC, Phenotyping, Validation, Genes, Multimer, SampleSheet
# ---------------------------------------------------------------------------


def test_cell_qc_funnel_skips_mito_gate_when_no_mt_genes():
    """Mito-free matrices must not lose 100% of cells to min_mito_pct > 0."""
    from tcrsift.qc import cell_qc_funnel

    X = np.full((10, 5), 50.0)
    adata = ad.AnnData(X)
    adata.obs_names = [f"c{i}" for i in range(10)]
    adata.var_names = ["CD3D", "CD3E", "CD8A", "TRAC", "B2M"]
    adata.obs["sample"] = "S1"

    filtered, waterfall = cell_qc_funnel(
        adata,
        min_counts=10,
        max_counts=10000,
        min_genes=2,
        max_genes=100,
        min_mito_pct=2.0,
        max_mito_pct=15.0,
    )
    assert filtered.n_obs == 10
    assert not waterfall.empty


def test_mouse_v_genes_not_flagged_as_human_only():
    """Mouse TRAV12-2 / TRBV28 must not trigger human-only species detection."""
    from tcrsift.species import _HUMAN_ONLY_V_PREFIXES

    assert "TRAV12-2" not in _HUMAN_ONLY_V_PREFIXES
    assert "TRBV28" not in _HUMAN_ONLY_V_PREFIXES


def test_phenotype_synthesizes_cd3_and_applies_numeric_cd3_gate():
    """phenotype_cells synthesizes CD3 from CD3D/E/G; apply_cd3_gate enforces min_cd3_reads."""
    from tcrsift.phenotype import apply_cd3_gate, phenotype_cells

    X = np.array(
        [
            [5.0, 5.0, 5.0, 10.0, 0.0],  # CD3=15, CD8
            [0.0, 1.0, 0.0, 10.0, 0.0],  # CD3=1, below min_cd3_reads=3
        ]
    )
    adata = ad.AnnData(X)
    adata.obs_names = ["c1", "c2"]
    adata.var_names = ["CD3D", "CD3E", "CD3G", "CD8A", "CD4"]
    adata.obs["CD3D"] = X[:, 0]
    adata.obs["CD3E"] = X[:, 1]
    adata.obs["CD3G"] = X[:, 2]
    adata.obs["CD8"] = X[:, 3]
    adata.obs["CD4"] = X[:, 4]
    adata.obs["filter:min_cd3"] = [True, True]

    out = phenotype_cells(adata)
    assert "CD3" in out.obs.columns
    assert list(out.obs["CD3"]) == [15.0, 1.0]

    gated = apply_cd3_gate(out, min_cd3_reads=3)
    assert list(gated.obs_names) == ["c1"]


def test_cd4_cd8_doublet_mask_spares_ambient_single_positive():
    """Single-positive CD8 cell with 1 ambient CD4 UMI should not be called a doublet."""
    from tcrsift.qc import cd4_cd8_doublet_mask

    X = np.array(
        [
            [10.0, 1.0, 50.0, 0.0],  # clear CD8+ with 1 ambient CD4 UMI (ratio 50x)
            [10.0, 10.0, 10.0, 0.0],  # true CD4/CD8 doublet (ratio 1x)
        ]
    )
    adata = ad.AnnData(X)
    adata.var_names = ["CD3E", "CD4", "CD8A", "CD8B"]
    mask = cd4_cd8_doublet_mask(adata, max_imbalance_ratio=5.0)
    assert list(mask) == [False, True]


def test_validation_rejects_bad_cdr3_and_bool_nan_numeric():
    """validate_cdr3_sequence rejects non-str/*; validate_numeric_param rejects bool/NaN/Inf."""
    from tcrsift.validation import (
        TCRsiftValidationError,
        validate_cdr3_sequence,
        validate_numeric_param,
    )

    assert not validate_cdr3_sequence(12345)
    assert not validate_cdr3_sequence("CASS*EQYF")
    with pytest.raises(TCRsiftValidationError):
        validate_numeric_param(True, "min_cells")
    with pytest.raises(TCRsiftValidationError):
        validate_numeric_param(float("nan"), "min_freq")
    with pytest.raises(TCRsiftValidationError):
        validate_numeric_param(float("inf"), "min_freq")


def test_ensembl_ids_to_symbols_auto_species():
    """species='auto' should not crash ensembl_ids_to_symbols."""
    from tcrsift.genes import ensembl_ids_to_symbols

    res = ensembl_ids_to_symbols(["ENSG_NONEXISTENT_1"], species="auto")
    assert isinstance(res, dict)


def test_multimer_tie_is_deterministic_and_unconfident():
    """Tied multimer counts should break deterministically with multimer_confident=False."""
    from tcrsift.multimer import summarize_clone_multimer_specificity

    obs = pd.DataFrame(
        {
            "CDR3ab": ["A_B"] * 4,
            "multimer_call": ["PepB", "PepA", "PepB", "PepA"],
            "is_multimer_specific": [True] * 4,
        }
    )
    res = summarize_clone_multimer_specificity(obs)
    assert len(res) == 1
    assert res.iloc[0]["multimer_call"] == "PepA"
    assert bool(res.iloc[0]["multimer_confident"]) is False


def test_sample_sheet_to_dataframe_includes_all_fields():
    """SampleSheet.to_dataframe should include antigen_sequence, til_csv, til_h5ad."""
    from tcrsift.sample_sheet import Sample, SampleSheet

    sheet = SampleSheet(
        samples=[
            Sample(
                sample="S1",
                antigen_sequence="ELAGIGILTV",
                til_csv="/tmp/til.csv",
                til_h5ad="/tmp/til.h5ad",
            )
        ]
    )
    df = sheet.to_dataframe()
    assert df.loc[0, "antigen_sequence"] == "ELAGIGILTV"
    assert df.loc[0, "til_csv"] == "/tmp/til.csv"
    assert df.loc[0, "til_h5ad"] == "/tmp/til.h5ad"


# ---------------------------------------------------------------------------
# Batch 2: Scoring, Filtering, Selection, TIL, Attribution
# ---------------------------------------------------------------------------


def test_analytical_logistic_threshold_uncapped_and_unreachable():
    """analytical_logistic_threshold computes exact closed-form threshold and inf."""
    from tcrsift.model import analytical_logistic_threshold

    # Within [0, 1] frequency domain
    w = 10.0
    fdr = 0.15
    expected = np.log((1.0 - fdr) / fdr) / w
    assert analytical_logistic_threshold(w, fdr) == pytest.approx(expected)
    # Uncapped when max_frequency=inf (no 100-cap artifact)
    assert analytical_logistic_threshold(0.01, fdr, max_frequency=float("inf")) == pytest.approx(
        np.log((1.0 - fdr) / fdr) / 0.01
    )
    # Unreachable when exact threshold exceeds max_frequency=1.0
    assert np.isinf(analytical_logistic_threshold(0.01, fdr, max_frequency=1.0))
    assert np.isinf(analytical_logistic_threshold(0.0, fdr))
    assert np.isinf(analytical_logistic_threshold(-1.0, fdr))
    assert np.isinf(analytical_logistic_threshold(1.0, 0.0))


def test_subset_percentile_logsf_underflow_and_zero_high_ties():
    """calibrate_signature_background uses -logsf (no tail collapse) and gives 0.0 when high_cells==0."""
    from tcrsift.signature_background import calibrate_signature_background

    clones = (
        ["C_hi"] * 15
        + ["C_med"] * 10
        + ["C_zero1"] * 10
        + ["C_zero2"] * 10
        + [f"C_bg{i}" for i in range(400)]
    )
    sig_vals = [10.0] * 15 + [10.0] * 10 + [0.0] * 10 + [0.0] * 10 + [0.0] * 400
    per_cell = pd.DataFrame(
        {
            "donor": "D1",
            "sample": "TIL",
            "lineage": "CD8",
            "CDR3ab": clones,
            "signature_Neo": sig_vals,
        }
    )
    scores = (
        per_cell.groupby(["donor", "sample", "lineage", "CDR3ab"], as_index=False)
        .agg(cells=("signature_Neo", "size"), signature_Neo=("signature_Neo", "mean"))
    )
    calibrated, _ = calibrate_signature_background(
        scores, per_cell, signatures=["Neo"], quantile=0.95, draws=1000, seed=0
    )
    by_clone = calibrated.set_index("CDR3ab")["signature_Neo_subset_percentile"]
    assert by_clone["C_zero1"] == 0.0
    assert by_clone["C_zero2"] == 0.0
    assert by_clone["C_hi"] > by_clone["C_med"]


def test_receptor_gene_guard_precision():
    """IGHMBP2 must not match RECEPTOR_GENE_RE; TRBV19 must be rejected by score_signature."""
    from tcrsift.signature_methods import Signature, is_receptor_gene, score_signature

    assert not is_receptor_gene("IGHMBP2")
    assert is_receptor_gene("TRBV19")
    assert is_receptor_gene("TRAV12-2")

    expr = pd.DataFrame(
        {"TRBV19": [100.0, 0.0], "CXCL13": [2.0, 0.0], "GZMB": [2.0, 0.0]},
        index=["c1", "c2"],
    )
    sig = Signature(name="test", genes_up=("TRBV19", "CXCL13"))
    with pytest.raises(ValueError, match="receptor transcripts"):
        score_signature(expr, sig, combine="mean", log1p=False)


def test_insilico_filter_constant_score_does_not_pass_all():
    """A constant score column must not pass 100% of clones via best_tie rescue."""
    from tcrsift.insilico_filter import FilterPredicate, apply_insilico_filter

    df = pd.DataFrame(
        {
            "CDR3ab": [f"A{i}_B{i}" for i in range(10)],
            "const_score": [5.0] * 10,
        }
    )
    out = apply_insilico_filter(
        df, [FilterPredicate(score="const_score", direction="high", percentile=0.1)]
    )
    assert not out["insilico_pass"].all()


def test_percentile_rank_score_require_complete_missing_column():
    """require_complete=True must return NaN/empty when a requested metric column is absent."""
    from tcrsift.selection import percentile_rank_score

    df = pd.DataFrame({"CDR3ab": ["A_B", "C_D"], "present_metric": [1.0, 2.0]})
    out = percentile_rank_score(
        df,
        terms=[{"col": "present_metric"}, {"col": "missing_metric"}],
        require_complete=True,
    )
    assert all(np.isnan(v) for v in out.values())


def test_attribution_dual_alpha_merge_feeds_prior():
    """dual_alpha_merge should fold merged partner counts into canon_prior."""
    from tcrsift.attribution import attribute_cells
    from tcrsift.config import AttributionConfig

    obs = pd.DataFrame(
        {
            "TRA_1_cdr3": ["CAA1"] * 5 + ["CAA2"] * 3,
            "TRA_2_cdr3": ["CAA2"] * 5 + [""] * 3,
            "TRB_1_cdr3": ["CBB1"] * 8,
            "sample": ["S1"] * 8,
            "TRA_1_umis": [10] * 8,
            "TRA_2_umis": [8] * 5 + [0] * 3,
            "TRB_1_umis": [10] * 8,
        }
    )
    cfg = AttributionConfig(enabled=True, dual_alpha_merge=True, dual_alpha_split=True)
    _weights, merge_map, stats = attribute_cells(obs, cfg)
    assert len(merge_map) >= 1
    assert isinstance(stats, dict)


# ---------------------------------------------------------------------------
# Batch 3: Annotation, Pgen, SeqProb, Cache, Unify, Cohort
# ---------------------------------------------------------------------------


def test_annotate_iedb_uniprot_tumor_self_and_bacterial_haemophilus():
    """MART-1/HER2 UniProt labels -> tumor_self; Haemophilus influenzae -> bacterial (not viral)."""
    from tcrsift.annotate import (
        CATEGORY_BACTERIAL,
        CATEGORY_TUMOR_SELF,
        _flag_viral,
        canonicalize_antigens,
        classify_category,
    )

    df = pd.DataFrame(
        {
            "species": ["Homo sapiens", "Homo sapiens", "Haemophilus influenzae"],
            "antigen_name": [
                "Melanoma antigen recognized by T-cells 1 (UniProt:Q16655)",
                "Receptor tyrosine-protein kinase erbB-2 (UniProt:P04626)",
                "Outer membrane protein P2",
            ],
            "antigen_gene": ["", "", "ompP2"],
            "epitope": ["ELAGIGILTV", "KIFGSLAFL", "VAAKLEK"],
        }
    )
    canon = canonicalize_antigens(df["antigen_name"], df["species"])
    viral = _flag_viral(df)
    cat = classify_category(df["species"], df["antigen_gene"], canonical_proteins=canon)

    assert not bool(viral.iloc[0])
    assert not bool(viral.iloc[1])
    assert not bool(viral.iloc[2]), "Haemophilus influenzae must not be flagged as viral influenza"
    assert cat.iloc[0] == CATEGORY_TUMOR_SELF
    assert cat.iloc[1] == CATEGORY_TUMOR_SELF
    assert cat.iloc[2] == CATEGORY_BACTERIAL


def test_pgen_single_handles_pd_na():
    """pgen_single must not raise TypeError when V/J gene is pd.NA."""
    from tcrsift.pgen import pgen_single

    val = pgen_single("CASSLGQAYEQYF", v_gene=pd.NA, j_gene=pd.NA, chain="beta")
    assert np.isfinite(val) and val < 0


def test_train_model_default_order_is_2():
    """pgen_models.train_model default order must match KmerProbabilityModel (order=2)."""
    from tcrsift.pgen_models import train_model

    sig = inspect.signature(train_model)
    assert sig.parameters["order"].default == 2


def test_datacache_extract_zip_rejects_zip_slip(tmp_path):
    """_extract_zip must reject archive members attempting path traversal."""
    from tcrsift.datacache import DownloadError, _extract_zip

    zpath = tmp_path / "evil.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("../escaped.txt", "pwned")

    with pytest.raises(DownloadError, match="unsafe archive member"):
        _extract_zip(zpath, tmp_path / "out", member=None)


def test_unify_dotted_conditions():
    """compute_condition_statistics handles dotted condition names."""
    from tcrsift.unify import compute_condition_statistics

    df = pd.DataFrame(
        {
            "CDR3ab": ["A_B", "C_D"],
            "src1.condition_pool.1.frac": [0.10, 0.0],
            "src1.condition_pool.2.frac": [0.0, 0.05],
        }
    )
    out = compute_condition_statistics(df, conditions=["pool.1", "pool.2"], verbose=False)
    assert out.loc[0, "combined.pool.1.frac.sum"] == pytest.approx(0.10)
    assert out.loc[1, "combined.pool.2.frac.sum"] == pytest.approx(0.05)


def test_cohort_overlap_ignores_empty_cdr3_beta():
    """Empty CDR3_beta strings must not create bogus shared clones across donors."""
    from tcrsift.cohort import cohort_overlap

    cohort_long = pd.DataFrame(
        {
            "donor": ["D1", "D1", "D2", "D2"],
            "CDR3ab": ["A1_B1", "A2_", "A3_B2", "A4_"],
            "CDR3_beta": ["B1", "", "B2", ""],
            "cells": [5, 5, 5, 5],
        }
    )
    out = cohort_overlap(cohort_long, include_beta_only=True)
    assert out["jaccard_beta_only"].loc["D1", "D2"] == 0.0


# ---------------------------------------------------------------------------
# Batch 4: Construct Assembly
# ---------------------------------------------------------------------------


def test_assemble_qc_catches_premature_stop_and_bare_mature_offset():
    """_nt_contig_in_frame_no_premature_stop catches internal stops; validate_sequences handles bare mature."""
    from tcrsift.assemble.constants import HUMAN_TRAC_AA
    from tcrsift.assemble.qc import _nt_contig_in_frame_no_premature_stop, validate_sequences

    # ATG (M) + GCT (A) + TAA (*) + GCT (A) -> internal stop at codon 3
    assert not _nt_contig_in_frame_no_premature_stop("ATGGCTTAAGCT", "MAGA")

    # Bare mature chain (no leader) with valid FGxG J-motif and human TRAC constant
    vdj = "AQSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQLLLKYYSGDPVVQGVNGFEAEFSKSNSSFHLRKASVHWSDSAVYFCAVSGGYQKVTFGTGTKLQVIP"
    df = pd.DataFrame(
        {
            "CDR3ab": ["CAVSGGYQKVTF_CASSLGQAYEQYF"],
            "CDR3_alpha": ["CAVSGGYQKVTF"],
            "alpha_v_gene": ["TRAV12-2"],
            "alpha_j_gene": ["TRAJ13"],
            "vdj_alpha_aa": [vdj],
            "full_alpha_aa": [vdj + HUMAN_TRAC_AA],
        }
    )
    msgs = validate_sequences(df, strict=False, species="human")
    junction_errs = [m for m in msgs if "junction" in str(m).lower() and m.severity == "load_bearing"]
    assert not junction_errs


def test_assemble_codons_forbids_golden_gate_and_ncoi():
    """DEFAULT_FORBIDDEN_MOTIFS includes NcoI, BsaI, and BsmBI."""
    from tcrsift.assemble.codons import DEFAULT_FORBIDDEN_MOTIFS

    for motif in ("CCATGG", "GGTCTC", "GAGACC", "CGTCTC", "GAGACG"):
        assert motif in DEFAULT_FORBIDDEN_MOTIFS


def test_assemble_skips_missing_chain_on_unpaired_clone():
    """Beta-only clone must not get a leader+constant fabricated for missing alpha chain."""
    from tcrsift.assemble import assemble_full_sequences

    df = pd.DataFrame(
        {
            "CDR3ab": ["_CASSLGQAYEQYF"],
            "CDR3_alpha": [None],
            "CDR3_beta": ["CASSLGQAYEQYF"],
            "VDJ_alpha_aa": [None],
            "VDJ_beta_aa": [
                "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSLGQAYEQYFGPGTRLTVT"
            ],
            "beta_v_gene": ["TRBV6-5"],
            "beta_j_gene": ["TRBJ2-7"],
        }
    )
    out = assemble_full_sequences(
        df,
        alpha_leader="CD28",
        beta_leader="CD8A",
        force_alpha_leader="CD28",
        include_constant=True,
        linker=None,
        verbose=False,
        show_progress=False,
    )
    assert (
        "full_alpha_aa" not in out.columns
        or pd.isna(out.loc[0, "full_alpha_aa"])
        or out.loc[0, "full_alpha_aa"] == ""
    )
    assert pd.notna(out.loc[0, "full_beta_aa"]) and len(out.loc[0, "full_beta_aa"]) > 100


def test_beta_c_base_unpads_zero_padded_trbc():
    """_beta_c_base maps TRBC01/TRBC02 to TRBC1/TRBC2."""
    from tcrsift.assemble.constants import _beta_c_base

    assert _beta_c_base("TRBC01*01") == "TRBC1"
    assert _beta_c_base("TRBC02") == "TRBC2"


# ---------------------------------------------------------------------------
# Batch 5: __init__ shadowing, Format, Plots, Report, CLI, Pipeline
# ---------------------------------------------------------------------------


def test_top_level_annotate_cells_and_tcrs_remain_callable_after_submodule_import():
    """Importing submodules must not break tcrsift.annotate_cells / tcrsift.annotate_tcrs calls."""
    import tcrsift
    from tcrsift.annotate_cells import score_reference  # noqa: F401
    from tcrsift.annotate_tcrs import prism_score  # noqa: F401

    assert callable(tcrsift.annotate_cells)
    assert callable(tcrsift.annotate_tcrs)


def test_pretty_sample_and_method_preserve_non_marker_names():
    """Non-marker hyphenated names (Patient1-TIL, MART-1) and non-marker methods are preserved."""
    from tcrsift.format import pretty_method, pretty_sample

    assert pretty_sample("Patient1-TIL") == "Patient1-TIL"
    assert pretty_sample("Patient1-PBMC") == "Patient1-PBMC"
    assert pretty_sample("MART-1") == "MART-1"
    assert pretty_sample("TRAV12-2") == "TRAV12-2"
    assert pretty_sample("AIMpos_CTYneg-2") == "AIM⁺CTY⁻"
    assert pretty_method("tumor_blood") == "tumorblood"


def test_plots_handle_stringified_sample_frequencies_and_duplicate_var_names():
    """Plots parse CSV-round-tripped sample_frequencies and duplicate var_names cleanly."""
    from tcrsift.plots import _derive_n_methods, _expand_sample_frequencies
    from tcrsift.plots.state import _per_cell_signature

    df = pd.DataFrame(
        {
            "CDR3ab": ["A1_B1", "A2_B2"],
            "sample_frequencies": ["{'S1': 0.25, 'S2': 0.10}", "{'S1': 0.05}"],
        }
    )
    wide = _expand_sample_frequencies(df)
    assert list(wide.columns) == ["S1", "S2"]
    assert wide.loc["A1_B1", "S1"] == pytest.approx(0.25)
    assert list(_derive_n_methods(df)) == [2, 1]

    # Duplicate var_names + already log1p matrix
    adata = ad.AnnData(np.array([[2.0, 4.0, 1.0], [0.0, 2.0, 3.0]]))
    adata.var_names = ["GZMB", "GZMB", "PRF1"]
    adata.uns["log1p"] = {"base": None}
    sig = _per_cell_signature(adata, ["GZMB", "PRF1"])
    assert sig[0] == pytest.approx(1.5)
    assert sig[1] == pytest.approx(1.5)


def test_generate_report_html_escapes_report_name(tmp_path):
    """HTML figure report must escape report_name."""
    from tcrsift.plots.reports import generate_report

    generate_report(tmp_path, format="html", report_name='<script>alert(1)</script>')
    html_text = (tmp_path / "all-figures.html").read_text()
    assert "<script>" not in html_text
    assert "&lt;script&gt;" in html_text


def test_expand_dual_alpha_variants_empty_preserves_columns():
    """expand_dual_alpha_variants on an empty DataFrame preserves columns."""
    from tcrsift.report import expand_dual_alpha_variants

    empty = pd.DataFrame(columns=["CDR3ab", "CDR3_beta", "merged_alpha_partners"])
    out, vmap = expand_dual_alpha_variants(empty, obs=pd.DataFrame())
    assert out.empty
    assert {"CDR3ab", "dual_alpha_variant", "selected_clone"} <= set(out.columns)
    assert vmap == {}


def test_validate_run_args_accepts_config_contigs_dir(tmp_path):
    """validate_run_args should accept --leaders-from-contigs when config sets contigs_dir."""
    from tcrsift.config import TCRsiftConfig
    from tcrsift.validation import validate_run_args

    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,vdj_dir\nS1,/tmp\n")
    cfg = TCRsiftConfig._from_dict({"assemble": {"contigs_dir": str(tmp_path)}})
    args = argparse.Namespace(
        sample_sheet=str(sheet),
        output_dir=str(tmp_path / "out"),
        leaders_from_contigs=True,
        alpha_leader=None,
        beta_leader=None,
        contigs_dir=None,
        cellranger_dir=None,
    )
    validate_run_args(args, config=cfg)


def test_pipeline_no_assemble_still_runs_selection_and_multidonor_fdr_scope(tmp_path, capsys):
    """run_pipeline detects multi-donor private cohorts for fdr_scope='auto' and runs selection with no_assemble=True."""
    import gzip
    import io

    from scipy import io as spio
    from scipy import sparse

    from tcrsift.config import TCRsiftConfig
    from tcrsift.pipeline import run_pipeline
    from tcrsift.sample_sheet import Sample, SampleSheet

    def _write_sample(base_dir, prefix, cdr3a, cdr3b, n_cells):
        vdj_dir = base_dir / "vdj"
        gex_dir = base_dir / "filtered_feature_bc_matrix"
        vdj_dir.mkdir(parents=True)
        gex_dir.mkdir(parents=True)
        barcodes = [f"{prefix}_BC{i}-1" for i in range(n_cells)]

        rows = []
        for bc in barcodes:
            rows.append(
                {
                    "barcode": bc,
                    "is_cell": True,
                    "contig_id": f"{bc}_c1",
                    "high_confidence": True,
                    "chain": "TRA",
                    "v_gene": "TRAV12-2",
                    "d_gene": "",
                    "j_gene": "TRAJ13",
                    "c_gene": "TRAC",
                    "full_length": True,
                    "productive": True,
                    "cdr3": cdr3a,
                    "cdr3_nt": "TGT" + "GCT" * (len(cdr3a) - 2) + "TTC",
                    "reads": 100,
                    "umis": 10,
                }
            )
            rows.append(
                {
                    "barcode": bc,
                    "is_cell": True,
                    "contig_id": f"{bc}_c2",
                    "high_confidence": True,
                    "chain": "TRB",
                    "v_gene": "TRBV6-5",
                    "d_gene": "TRBD1",
                    "j_gene": "TRBJ2-7",
                    "c_gene": "TRBC2",
                    "full_length": True,
                    "productive": True,
                    "cdr3": cdr3b,
                    "cdr3_nt": "TGT" + "GCT" * (len(cdr3b) - 2) + "TTC",
                    "reads": 100,
                    "umis": 10,
                }
            )
        pd.DataFrame(rows).to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        genes = ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "B2M", "ACTB"]
        # n_genes x n_cells matrix: strong CD3/CD8 expression, zero CD4
        mat = np.zeros((len(genes), n_cells), dtype=np.int32)
        mat[0:3, :] = 10  # CD3D/E/G
        mat[4:6, :] = 15  # CD8A/B
        mat[6:, :] = 500  # housekeeping genes so counts clear QC
        with gzip.open(gex_dir / "barcodes.tsv.gz", "wt") as f:
            f.write("\n".join(barcodes) + "\n")
        with gzip.open(gex_dir / "features.tsv.gz", "wt") as f:
            f.write(
                "\n".join(f"G{i}\t{g}\tGene Expression" for i, g in enumerate(genes)) + "\n"
            )
        buf = io.BytesIO()
        spio.mmwrite(buf, sparse.csc_matrix(mat))
        with gzip.open(gex_dir / "matrix.mtx.gz", "wb") as f:
            f.write(buf.getvalue())
        return vdj_dir, gex_dir

    vdj1, gex1 = _write_sample(tmp_path / "d1", "D1", "CAVSGGYQKVTF", "CASSLGQAYEQYF", 5)
    vdj2, gex2 = _write_sample(tmp_path / "d2", "D2", "CAVRDSSYKLIF", "CASSQETQYF", 5)

    sheet = SampleSheet(
        samples=[
            Sample(
                sample="D1_AIMpos",
                patient_id="Donor1",
                enrichment_method="AIMpos",
                vdj_dir=str(vdj1),
                gex_dir=str(gex1),
            ),
            Sample(
                sample="D2_AIMpos",
                patient_id="Donor2",
                enrichment_method="AIMpos",
                vdj_dir=str(vdj2),
                gex_dir=str(gex2),
            ),
        ]
    )
    cfg = TCRsiftConfig._from_dict(
        {
            "load": {"min_genes": 2, "min_counts": 10},
            "clonotype": {"add_pgen_ppost": False},
            "filter": {"fdr_scope": "auto", "min_cells": 2, "min_frequency": 0.01},
            "output": {"generate_plots": False, "generate_report": False},
            "selection": {
                "exclude_viral": False,
                "rules": {"shared": {"include_tiers": ["tier1", "tier2", "tier3"]}},
            },
        }
    )
    out_dir = tmp_path / "run_out"
    res = run_pipeline(sheet, out_dir, config=cfg, no_assemble=True)

    # Multi-donor private cohort must trigger per-donor FDR scope.
    captured = capsys.readouterr().out
    assert "FDR scope: per-donor" in captured
    # Selection must still run and emit selected_clones.csv even with no_assemble=True.
    sel_csv = res.data_dir / "selected_clones.csv"
    assert sel_csv.is_file()
    sel_df = pd.read_csv(sel_csv)
    assert len(sel_df) == 2
