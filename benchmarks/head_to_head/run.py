#!/usr/bin/env python3
"""Run head-to-head comparisons of tcrsift against open-source baselines across 4 stages.

    PYTHONPATH=.:.venv/lib/python3.12/site-packages /tmp/h2h_venv/bin/python \\
        benchmarks/head_to_head/run.py

Writes four CSV tables and ``summary.json`` to ``benchmarks/head_to_head/results/``:

1. ``h2h_prioritize_baselines.csv`` — Ex vivo prioritization across 5 human TIL and
   cancer-blood cohorts: ``tcrsift`` vs ``clone_size``, single-gene ``raw_CXCL13``,
   ``scanpy.tl.score_genes`` (``NeoTCR8`` / ``NeoTCR4`` / ``NeoTCR_PBL``), and
   ``scanpy`` composite top-K.
2. ``h2h_assembly_vs_stitchr.csv`` — Full-length TCR construct assembly:
   ``tcrsift assemble`` vs ``stitchr`` / ``thimble v1.3.1`` across Cell Ranger 6
   donor-VDJ human shortlisted clones (Minervina 2022, n=46), 5'-truncated human clones
   (Fischer 2021, n=160), and murine Cell Ranger VDJ shortlists (MC38 + Khatun, n=101).
3. ``h2h_pgen_vs_olga.csv`` — TCR generation probability & publicness:
   ``tcrsift`` zero-dependency k-mer ``Pgen``/``Ppost`` vs ``OLGA v1.3.0`` exact
   V(D)J dynamic programming across human (n=1,500) and murine (n=1,781) repertoires.
4. ``h2h_multimer_callers.csv`` — Multiplexed pMHC dextramer calling on 10x Healthy
   Donors 1 & 2: ``tcrsift --call-multimers`` clonal concordance vs 10x per-cell
   threshold, CLR ``z >= 3``, and bimodal GMM callers.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import pearsonr, rankdata, spearmanr

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def auroc(scores, positive) -> float:
    scores = np.asarray(scores, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    ok = np.isfinite(scores)
    scores, positive = scores[ok], positive[ok]
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auroc_ci(scores, positive, draws: int = 2000, seed: int = 0) -> tuple[float, float]:
    scores = np.asarray(scores, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    ok = np.isfinite(scores)
    pos, neg = scores[ok & positive], scores[ok & ~positive]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(draws):
        p = rng.choice(pos, len(pos))
        n = rng.choice(neg, len(neg))
        vals.append(auroc(np.concatenate([p, n]), np.arange(len(p) + len(n)) < len(p)))
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return float(lo), float(hi)


def _scanpy_score_genes_harmonized(
    adata: ad.AnnData,
    genes_up: list[str] | tuple[str, ...],
    score_name: str,
    genes_down: list[str] | tuple[str, ...] = (),
) -> np.ndarray:
    """Run ``scanpy.tl.score_genes`` after harmonizing HGNC/Ensembl/hg19 symbols fairly."""
    from tcrsift.genes import adata_symbol_array, harmonize_gene_symbols

    all_genes = list(dict.fromkeys([*genes_up, *genes_down]))
    sym_upper = harmonize_gene_symbols(adata_symbol_array(adata), all_genes)
    tmp = adata.copy()
    tmp.var_names = pd.Index([str(s) for s in sym_upper])
    tmp.var_names_make_unique()
    var_set = set(tmp.var_names)
    present_up = [g.upper() for g in genes_up if g.upper() in var_set]
    present_dn = [g.upper() for g in genes_down if g.upper() in var_set]
    ctrl_sz = min(50, max(5, tmp.n_vars // 4))
    out = np.zeros(adata.n_obs, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if present_up:
            sc.tl.score_genes(tmp, gene_list=present_up, score_name=f"{score_name}_up", ctrl_size=ctrl_sz, n_bins=25)
            out = out + tmp.obs[f"{score_name}_up"].to_numpy(dtype=float)
        if present_dn:
            sc.tl.score_genes(tmp, gene_list=present_dn, score_name=f"{score_name}_dn", ctrl_size=ctrl_sz, n_bins=25)
            out = out - tmp.obs[f"{score_name}_dn"].to_numpy(dtype=float)
    return out


def _extract_single_gene_log1p(adata: ad.AnnData, gene: str) -> np.ndarray:
    from tcrsift.signature_methods import expression_frame_from_adata

    df = expression_frame_from_adata(adata, [gene], on_missing="ignore")
    if gene not in df.columns:
        return np.zeros(adata.n_obs, dtype=float)
    return df[gene].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# 1. Ex vivo Prioritization Head-to-Head
# ---------------------------------------------------------------------------


def _eval_method_row(
    cohort: str,
    contrast_name: str,
    method: str,
    scores: np.ndarray,
    pos_mask: np.ndarray,
    neg_mask: np.ndarray,
    selected_mask: np.ndarray,
    all_pos_mask: np.ndarray,
    all_neg_mask: np.ndarray,
    treg_mask: np.ndarray | None = None,
) -> dict:
    eval_sub = pos_mask | neg_mask
    auc = auroc(scores[eval_sub], pos_mask[eval_sub])
    lo, hi = auroc_ci(scores[eval_sub], pos_mask[eval_sub])
    k = int(selected_mask.sum())
    tp = int((selected_mask & all_pos_mask).sum())
    fp_tested = int((selected_mask & all_neg_mask).sum())
    prec_tested = round(tp / (tp + fp_tested), 4) if (tp + fp_tested) > 0 else float("nan")
    treg_sel = int((selected_mask & treg_mask).sum()) if treg_mask is not None else 0
    treg_frac = round(treg_sel / k, 4) if (treg_mask is not None and k > 0) else float("nan")
    return {
        "cohort": cohort,
        "contrast": contrast_name,
        "method": method,
        "n_positive": int(pos_mask.sum()),
        "n_negative": int(neg_mask.sum()),
        "auroc": round(auc, 4) if np.isfinite(auc) else np.nan,
        "auroc_low": round(lo, 4) if np.isfinite(lo) else np.nan,
        "auroc_high": round(hi, 4) if np.isfinite(hi) else np.nan,
        "shortlist_k": k,
        "positives_shortlisted": tp,
        "negatives_shortlisted": fp_tested,
        "tested_shortlist_precision": prec_tested,
        "treg_shortlisted": treg_sel if treg_mask is not None else np.nan,
        "treg_shortlist_fraction": treg_frac,
    }


def _topk_mask_among_expanded(df: pd.DataFrame, score_col: str, k: int, cell_col: str = "cell_count") -> np.ndarray:
    elig = df[cell_col].to_numpy(dtype=float) >= 2
    idx_elig = np.where(elig & np.isfinite(df[score_col].to_numpy(dtype=float)))[0]
    if len(idx_elig) == 0 or k <= 0:
        return np.zeros(len(df), dtype=bool)
    # Break ties deterministically by clone size then index
    order = sorted(
        idx_elig,
        key=lambda i: (float(df[score_col].iloc[i]), float(df[cell_col].iloc[i]), -i),
        reverse=True,
    )
    out = np.zeros(len(df), dtype=bool)
    out[order[: min(k, len(order))]] = True
    return out


def run_prioritize_baselines(
    til_data_dir: Path,
    til_work_dir: Path,
    blood_cd4_data_dir: Path,
    blood_cd4_work_dir: Path,
) -> list[dict]:
    from benchmarks.public_til.run import (
        any_pass,
        chiffelle_labels,
        collapse_by_beta,
        meng_cells,
        meng_scores,
    )
    from tcrsift._schema import truthy_flag_series
    from tcrsift.loader import load_samples
    from tcrsift.phenotype import phenotype_cells
    from tcrsift.prioritize import (
        PRIORITY_CONTEXT_SIGNATURES,
        QC_DEFAULTS,
        REGISTRY,
        PrioritizeConfig,
        filter_cells,
        resolve_signatures,
        select_round_robin,
        signature_pass_mask,
    )
    from tcrsift.sample_sheet import load_sample_sheet
    from tcrsift.signature_background import calibrate_signature_background
    from tcrsift.signature_methods import score_by_name
    from tcrsift.til_prioritize import _score_within_samples

    rows: list[dict] = []

    # --- 1A. Meng 2023 PDAC CD8+ TIL (GSE254250; 63 TR vs 25 NTR) ---
    adata_m = meng_cells(til_data_dir)
    sigs_cd8 = resolve_signatures("solid-tumor", "cd8")
    for s in sigs_cd8:
        adata_m.obs[f"signature_{s}"] = _score_within_samples(adata_m, s)
    adata_m.obs["raw_cxcl13"] = _extract_single_gene_log1p(adata_m, "CXCL13")
    adata_m.obs["sc_neotcr8"] = _scanpy_score_genes_harmonized(
        adata_m, REGISTRY["NeoTCR8"].genes_up, "sc_neotcr8", REGISTRY["NeoTCR8"].genes_down
    )
    sc_comp_cols = []
    for s in sigs_cd8:
        cname = f"sc_{s}"
        adata_m.obs[cname] = _scanpy_score_genes_harmonized(
            adata_m, REGISTRY[s].genes_up, cname, REGISTRY[s].genes_down
        )
        sc_comp_cols.append(cname)
    adata_m.obs["sc_composite"] = adata_m.obs[sc_comp_cols].mean(axis=1)

    m_scores = meng_scores(adata_m, sigs_cd8)
    extra_m = (
        adata_m.obs.groupby(["donor", "sample", "lineage", "CDR3ab"], observed=True)[
            ["raw_cxcl13", "sc_neotcr8", "sc_composite"]
        ]
        .mean()
        .reset_index()
    )
    m_scores = m_scores.merge(extra_m, on=["donor", "sample", "lineage", "CDR3ab"], how="left")
    m_scores["cell_count"] = m_scores["cells"]
    m_clones = m_scores[["donor", "CDR3ab", "label", "cells", "cell_count", "frequency"]].copy()
    m_clones["eligible_for_review"] = (
        m_clones["cells"].ge(2) & m_clones["frequency"].ge(0.001) & any_pass(m_scores, sigs_cd8)
    )
    m_picked = select_round_robin(
        m_clones,
        m_scores,
        sigs_cd8,
        max_clones=200,
        evidence="either",
        priority_signatures=PRIORITY_CONTEXT_SIGNATURES["solid-tumor"],
    )
    sel_tcr_m = m_picked.set_index(["donor", "CDR3ab"]).loc[
        pd.MultiIndex.from_frame(m_scores[["donor", "CDR3ab"]]), "selected_for_review"
    ].to_numpy(dtype=bool)
    k_m = int(sel_tcr_m.sum())
    pos_m = m_scores["label"].eq("TR").to_numpy()
    neg_m = m_scores["label"].eq("NTR").to_numpy()

    for m_name, scol, sel_mask in (
        ("clone_size", "size_percentile", _topk_mask_among_expanded(m_scores, "size_percentile", k_m)),
        ("raw_CXCL13", "raw_cxcl13", _topk_mask_among_expanded(m_scores, "raw_cxcl13", k_m)),
        ("scanpy_score_genes_NeoTCR", "sc_neotcr8", _topk_mask_among_expanded(m_scores, "sc_neotcr8", k_m)),
        ("scanpy_composite_topk", "sc_composite", _topk_mask_among_expanded(m_scores, "sc_composite", k_m)),
        ("tcrsift_solid_tumor", "signature_TumorReactive_percentile", sel_tcr_m),
    ):
        rows.append(
            _eval_method_row("meng_pdac_cd8", "TR_vs_NTR", m_name, m_scores[scol].to_numpy(), pos_m, neg_m, sel_mask, pos_m, neg_m)
        )

    # --- 1B. Chiffelle / Oliveira Melanoma Baseline CD8+ TIL (GSE222448) ---
    ch_sheet_path = til_work_dir / "chiffelle_tumor.csv"
    ch_out = til_work_dir / "chiffelle_default"
    ch_labels = chiffelle_labels(til_data_dir)
    sheet_ch = load_sample_sheet(ch_sheet_path)
    adata_ch = load_samples(sheet_ch, **QC_DEFAULTS, verbose=False, show_progress=False)
    adata_ch.obs_names_make_unique()
    adata_ch = phenotype_cells(adata_ch, min_cd3_reads=3, verbose=False, show_progress=False)
    adata_ch, _, _ = filter_cells(adata_ch, PrioritizeConfig(sample_sheet=ch_sheet_path, output_dir=ch_out, context="solid-tumor"))
    sc.pp.normalize_total(adata_ch, target_sum=10_000)
    sc.pp.log1p(adata_ch)
    adata_ch.obs["donor"] = adata_ch.obs["patient_id"].astype(str)
    adata_ch.obs["raw_cxcl13"] = _extract_single_gene_log1p(adata_ch, "CXCL13")
    adata_ch.obs["sc_neotcr8"] = _scanpy_score_genes_harmonized(
        adata_ch, REGISTRY["NeoTCR8"].genes_up, "sc_neotcr8", REGISTRY["NeoTCR8"].genes_down
    )
    ch_sc_cols = []
    for s in sigs_cd8:
        cn = f"sc_{s}"
        adata_ch.obs[cn] = _scanpy_score_genes_harmonized(
            adata_ch, REGISTRY[s].genes_up, cn, REGISTRY[s].genes_down
        )
        ch_sc_cols.append(cn)
    adata_ch.obs["sc_composite"] = adata_ch.obs[ch_sc_cols].mean(axis=1)

    ch_cell_agg = (
        adata_ch.obs.groupby(["donor", "CDR3ab"], observed=True)[["raw_cxcl13", "sc_neotcr8", "sc_composite"]]
        .mean()
        .reset_index()
    )
    ch_clones = pd.read_csv(ch_out / "all_scored_clones.csv", low_memory=False)
    ch_clones["size_percentile"] = ch_clones.groupby("donor")["cell_count"].rank(pct=True)
    ch_clones = ch_clones.merge(ch_labels[["donor", "CDR3_beta", "label"]], on=["donor", "CDR3_beta"], how="left")
    ch_clones["label"] = ch_clones["label"].fillna("untested")
    ch_clones = ch_clones.merge(ch_cell_agg, on=["donor", "CDR3ab"], how="left")
    ch_clones["selected_for_review"] = truthy_flag_series(ch_clones["selected_for_review"])
    ch_tab = collapse_by_beta(
        ch_clones[
            [
                "donor",
                "CDR3_beta",
                "CDR3_alpha",
                "label",
                "cell_count",
                "size_percentile",
                "signature_TumorReactive_percentile_max",
                "raw_cxcl13",
                "sc_neotcr8",
                "sc_composite",
                "selected_for_review",
                "selection_rank",
                "selected_signature",
            ]
        ].copy(),
        ["selected_for_review"],
    )
    sel_tcr_ch = ch_tab["selected_for_review"].to_numpy(dtype=bool)
    k_ch = int(sel_tcr_ch.sum())
    pos_ch = ch_tab["label"].eq("TR").to_numpy()
    neg_ch = ch_tab["label"].eq("NTR").to_numpy()
    for m_name, scol, sel_mask in (
        ("clone_size", "size_percentile", _topk_mask_among_expanded(ch_tab, "size_percentile", k_ch)),
        ("raw_CXCL13", "raw_cxcl13", _topk_mask_among_expanded(ch_tab, "raw_cxcl13", k_ch)),
        ("scanpy_score_genes_NeoTCR", "sc_neotcr8", _topk_mask_among_expanded(ch_tab, "sc_neotcr8", k_ch)),
        ("scanpy_composite_topk", "sc_composite", _topk_mask_among_expanded(ch_tab, "sc_composite", k_ch)),
        ("tcrsift_solid_tumor", "signature_TumorReactive_percentile_max", sel_tcr_ch),
    ):
        rows.append(
            _eval_method_row("chiffelle_melanoma_cd8", "TR_vs_NTR", m_name, ch_tab[scol].to_numpy(), pos_ch, neg_ch, sel_mask, pos_ch, neg_ch)
        )

    # --- 1C. Pulliam 2024 MCC Peripheral Blood (--context blood-tumor, GSE227708) ---
    from benchmarks.public_blood_cd4.run import pulliam_s4_labels

    pul_sheet_path = blood_cd4_work_dir / "pulliam_blood_sheet.csv"
    pul_out = blood_cd4_work_dir / "pulliam_blood_tumor_default"
    s4_labels = pulliam_s4_labels(blood_cd4_data_dir)[["donor", "CDR3_beta", "label"]]
    sheet_pul = load_sample_sheet(pul_sheet_path)
    adata_pul = load_samples(sheet_pul, **QC_DEFAULTS, verbose=False, show_progress=False)
    adata_pul.obs_names_make_unique()
    adata_pul = phenotype_cells(adata_pul, min_cd3_reads=3, verbose=False, show_progress=False)
    adata_pul, _, _ = filter_cells(
        adata_pul, PrioritizeConfig(sample_sheet=pul_sheet_path, output_dir=pul_out, context="blood-tumor")
    )
    sc.pp.normalize_total(adata_pul, target_sum=10_000)
    sc.pp.log1p(adata_pul)
    adata_pul.obs["donor"] = adata_pul.obs["patient_id"].astype(str)
    sigs_bt = resolve_signatures("blood-tumor", "cd8")
    adata_pul.obs["raw_cxcl13"] = _extract_single_gene_log1p(adata_pul, "CXCL13")
    adata_pul.obs["sc_neotcr_pbl"] = _scanpy_score_genes_harmonized(
        adata_pul, REGISTRY["NeoTCR_PBL"].genes_up, "sc_neotcr_pbl", REGISTRY["NeoTCR_PBL"].genes_down
    )
    bt_sc_cols = []
    for s in sigs_bt:
        cn = f"sc_{s}"
        adata_pul.obs[cn] = _scanpy_score_genes_harmonized(
            adata_pul, REGISTRY[s].genes_up, cn, REGISTRY[s].genes_down
        )
        bt_sc_cols.append(cn)
    adata_pul.obs["sc_composite"] = adata_pul.obs[bt_sc_cols].mean(axis=1)
    pul_cell_agg = (
        adata_pul.obs.groupby(["donor", "CDR3ab"], observed=True)[["raw_cxcl13", "sc_neotcr_pbl", "sc_composite"]]
        .mean()
        .reset_index()
    )
    pul_clones = pd.read_csv(pul_out / "all_scored_clones.csv", low_memory=False)
    pul_clones["size_percentile"] = pul_clones.groupby("donor")["cell_count"].rank(pct=True)
    pul_clones = pul_clones.merge(s4_labels, on=["donor", "CDR3_beta"], how="left")
    pul_clones["label"] = pul_clones["label"].fillna("untested")
    pul_clones = pul_clones.merge(pul_cell_agg, on=["donor", "CDR3ab"], how="left")
    pul_clones["selected_for_review"] = truthy_flag_series(pul_clones["selected_for_review"])
    pul_tab = collapse_by_beta(
        pul_clones[
            [
                "donor",
                "CDR3_beta",
                "CDR3_alpha",
                "label",
                "cell_count",
                "size_percentile",
                "signature_NeoTCR_PBL_percentile_max",
                "raw_cxcl13",
                "sc_neotcr_pbl",
                "sc_composite",
                "selected_for_review",
                "selection_rank",
                "selected_signature",
            ]
        ].copy(),
        ["selected_for_review"],
    )
    sel_tcr_pul = pul_tab["selected_for_review"].to_numpy(dtype=bool)
    k_pul = int(sel_tcr_pul.sum())
    pos_pul = pul_tab["label"].eq("MCPyV").to_numpy()
    neg_pul = pul_tab["label"].eq("control_virus").to_numpy()
    for m_name, scol, sel_mask in (
        ("clone_size", "size_percentile", _topk_mask_among_expanded(pul_tab, "size_percentile", k_pul)),
        ("raw_CXCL13", "raw_cxcl13", _topk_mask_among_expanded(pul_tab, "raw_cxcl13", k_pul)),
        ("scanpy_score_genes_NeoTCR", "sc_neotcr_pbl", _topk_mask_among_expanded(pul_tab, "sc_neotcr_pbl", k_pul)),
        ("scanpy_composite_topk", "sc_composite", _topk_mask_among_expanded(pul_tab, "sc_composite", k_pul)),
        ("tcrsift_blood_tumor", "signature_NeoTCR_PBL_percentile_max", sel_tcr_pul),
    ):
        rows.append(
            _eval_method_row(
                "pulliam_mcc_blood_cd8",
                "MCPyV_vs_CMV_EBV_Flu",
                m_name,
                pul_tab[scol].to_numpy(),
                pos_pul,
                neg_pul,
                sel_mask,
                pos_pul,
                neg_pul,
            )
        )

    # --- 1D. Veatch 2022 Melanoma CD4+ TIL (GSE198265; 28 CD4_neo vs 1,122 unlabeled, 477 FOXP3+ Tregs) ---
    sigs_cd4 = resolve_signatures("solid-tumor", "cd4")
    vroot = blood_cd4_data_dir / "veatch"
    v_parts = []
    for pat in ("X197", "X198", "X205", "X422"):
        pdir = vroot / pat
        vdf = pd.read_csv(next(pdir.glob("*VDJ*")))
        vdf = vdf[
            vdf["TCRVb?"].eq(True)
            & vdf["Cell_Type_Experimental"].isin(["T_CD4_memory", "T_CD4_naive"])
            & vdf["TCR_Beta_Delta_CDR3_Translation_Dominant"].notna()
        ].copy()
        if "sample" in vdf.columns:
            vdf = vdf[~vdf["sample"].isin(["Multiplet", "Undetermined"])].copy()
            if pat == "X198":
                vdf = vdf[vdf["sample"].eq("8_26_main_harvest")].copy()
                vdf["sample_id"] = "X198_tumor"
            else:
                vdf["sample_id"] = pat + "_" + vdf["sample"].astype(str)
        else:
            vdf["sample_id"] = f"{pat}_tumor"
        gdf = pd.read_csv(next(pdir.glob("*MolsPerCell*")))
        mrna_cols = [c for c in gdf.columns if c != "Cell_Index" and not c.endswith("|pAbO")]
        mrna_genes = [c.split("|")[0].strip() for c in mrna_cols]
        joined = vdf.merge(gdf, on="Cell_Index", how="inner")
        X = sp.csr_matrix(joined[mrna_cols].to_numpy(dtype=float))
        fox_col = next((c for c, g in zip(mrna_cols, mrna_genes) if g.upper() == "FOXP3"), None)
        fox_pos = (joined[fox_col].to_numpy(dtype=float) > 0) if fox_col else np.zeros(len(joined), dtype=bool)
        is_neo = (
            joined["CD4_neo"].eq(1.0).to_numpy()
            if "CD4_neo" in joined.columns
            else joined["Neoantigen"].notna().to_numpy()
        )
        if "CD4_self" in joined.columns:
            is_self = joined["CD4_self"].eq(1.0).to_numpy()
        elif "self" in joined.columns:
            sval = joined["self"].fillna("").astype(str)
            is_self = (sval != "") & (~sval.str.contains("cd8", case=False)).to_numpy()
        else:
            is_self = np.zeros(len(joined), dtype=bool)
        obs = pd.DataFrame(
            {
                "donor": pat,
                "sample": joined["sample_id"].to_numpy(),
                "lineage": "cd4",
                "CDR3ab": joined["TCR_Beta_Delta_CDR3_Translation_Dominant"].astype(str).str.strip().to_numpy(),
                "cell_label": np.where(is_neo, "CD4_neo", np.where(is_self, "CD4_self", "untested")),
                "foxp3_pos": fox_pos.astype(float),
            },
            index=[f"{pat}_{cid}" for cid in joined["Cell_Index"]],
        )
        part = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=mrna_genes))
        part.var_names_make_unique()
        sc.pp.normalize_total(part, target_sum=10_000)
        sc.pp.log1p(part)
        for s in sigs_cd4:
            part.obs[f"signature_{s}"] = score_by_name(
                part, s, log1p=False, groups=part.obs["sample"], on_missing="ignore", species="human"
            ).to_numpy(dtype=float)
        v_parts.append(part)

    v_full = ad.concat(v_parts, join="inner")
    v_full.obs["raw_cxcl13"] = _extract_single_gene_log1p(v_full, "CXCL13")
    v_full.obs["sc_neotcr4"] = _scanpy_score_genes_harmonized(
        v_full, REGISTRY["NeoTCR4"].genes_up, "sc_neotcr4", REGISTRY["NeoTCR4"].genes_down
    )
    v_sc_cols = []
    for s in sigs_cd4:
        cn = f"sc_{s}"
        v_full.obs[cn] = _scanpy_score_genes_harmonized(
            v_full, REGISTRY[s].genes_up, cn, REGISTRY[s].genes_down
        )
        v_sc_cols.append(cn)
    v_full.obs["sc_composite"] = v_full.obs[v_sc_cols].mean(axis=1)

    stratum = ["donor", "sample", "lineage"]
    sig_cols_cd4 = [f"signature_{s}" for s in sigs_cd4]
    v_per_cell = v_full.obs[[*stratum, "CDR3ab", *sig_cols_cd4]].copy()
    v_grp = v_full.obs.groupby([*stratum, "CDR3ab"], observed=True)
    v_scores = v_grp[sig_cols_cd4].mean().join(v_grp.size().rename("cells")).reset_index()
    v_scores["frequency"] = v_scores["cells"] / v_scores.groupby(stratum, observed=True)["cells"].transform("sum")
    for col in sig_cols_cd4:
        v_scores[f"{col}_percentile"] = v_scores.groupby(stratum, observed=True)[col].rank(method="average", pct=True)
    v_scores, _ = calibrate_signature_background(v_scores, v_per_cell, sigs_cd4, quantile=0.99)
    for s in sigs_cd4:
        v_scores[f"signature_{s}_passes_cutoff"] = signature_pass_mask(v_scores, s, evidence="either")

    vc_grp = v_full.obs.groupby(["donor", "CDR3ab"], observed=True)
    v_clones = vc_grp.agg(
        cell_count=("sample", "size"),
        foxp3_pos_fraction=("foxp3_pos", "mean"),
        raw_cxcl13=("raw_cxcl13", "mean"),
        sc_neotcr4=("sc_neotcr4", "mean"),
        sc_composite=("sc_composite", "mean"),
        label=(
            "cell_label",
            lambda s: "CD4_neo"
            if (s == "CD4_neo").any()
            else ("CD4_self" if (s == "CD4_self").any() else "untested"),
        ),
    ).reset_index()
    max_cols = [f"signature_{s}_percentile" for s in sigs_cd4] + [f"signature_{s}_passes_cutoff" for s in sigs_cd4]
    v_sig_max = (
        v_scores.groupby(["donor", "CDR3ab"], observed=True)[[*max_cols, "frequency"]]
        .max()
        .rename(columns={"frequency": "max_frequency", **{c: f"{c}_max" for c in max_cols}})
        .reset_index()
    )
    v_clones = v_clones.merge(v_sig_max, on=["donor", "CDR3ab"], how="left")
    v_clones["size_percentile"] = v_clones.groupby("donor")["cell_count"].rank(pct=True)
    v_clones["is_treg_like"] = v_clones["foxp3_pos_fraction"] >= 0.50
    v_clones["any_sig_pass"] = v_clones[[f"signature_{s}_passes_cutoff_max" for s in sigs_cd4]].any(axis=1)
    v_clones["eligible_for_review"] = (
        v_clones["cell_count"].ge(2)
        & v_clones["max_frequency"].ge(0.001)
        & v_clones["any_sig_pass"]
        & ~v_clones["is_treg_like"]
    )
    v_picked = select_round_robin(
        v_clones,
        v_scores,
        sigs_cd4,
        max_clones=200,
        evidence="either",
        priority_signatures=PRIORITY_CONTEXT_SIGNATURES["solid-tumor"],
    )
    sel_tcr_v = v_picked.set_index(["donor", "CDR3ab"]).loc[
        pd.MultiIndex.from_frame(v_clones[["donor", "CDR3ab"]]), "selected_for_review"
    ].to_numpy(dtype=bool)
    k_v = int(sel_tcr_v.sum())
    exp_v = v_clones["cell_count"].ge(2).to_numpy()
    pos_v = exp_v & v_clones["label"].eq("CD4_neo").to_numpy()
    neg_v = exp_v & v_clones["label"].eq("untested").to_numpy()
    treg_v = v_clones["is_treg_like"].to_numpy()
    for m_name, scol, sel_mask in (
        ("clone_size", "size_percentile", _topk_mask_among_expanded(v_clones, "size_percentile", k_v)),
        ("raw_CXCL13", "raw_cxcl13", _topk_mask_among_expanded(v_clones, "raw_cxcl13", k_v)),
        ("scanpy_score_genes_NeoTCR", "sc_neotcr4", _topk_mask_among_expanded(v_clones, "sc_neotcr4", k_v)),
        ("scanpy_composite_topk", "sc_composite", _topk_mask_among_expanded(v_clones, "sc_composite", k_v)),
        ("tcrsift_solid_tumor_cd4_excl_treg", "signature_TumorReactive_percentile_max", sel_tcr_v),
    ):
        rows.append(
            _eval_method_row(
                "veatch_melanoma_cd4",
                "CD4_neo_vs_untested_cd4",
                m_name,
                v_clones[scol].to_numpy(),
                pos_v,
                neg_v,
                sel_mask,
                v_clones["label"].eq("CD4_neo").to_numpy(),
                v_clones["label"].eq("untested").to_numpy(),
                treg_mask=treg_v,
            )
        )

    # --- 1E. Oh 2020 Whole-Transcriptome Bladder CD4+ TIL (GSE149652; pre_sorted=CD4 --exclude-treg) ---
    oh_sheet_path = blood_cd4_work_dir / "oh_tumor_presorted.csv"
    oh_out = blood_cd4_work_dir / "oh_presorted_exclude_treg_0.50"
    sheet_oh = load_sample_sheet(oh_sheet_path)
    adata_oh = load_samples(sheet_oh, **QC_DEFAULTS, verbose=False, show_progress=False)
    adata_oh.obs_names_make_unique()
    adata_oh = phenotype_cells(adata_oh, min_cd3_reads=3, verbose=False, show_progress=False)
    adata_oh, _, _ = filter_cells(
        adata_oh,
        PrioritizeConfig(sample_sheet=oh_sheet_path, output_dir=oh_out, context="solid-tumor", tcell_type="cd4"),
    )
    sc.pp.normalize_total(adata_oh, target_sum=10_000)
    sc.pp.log1p(adata_oh)
    adata_oh.obs["donor"] = adata_oh.obs["patient_id"].astype(str)
    adata_oh.obs["raw_cxcl13"] = _extract_single_gene_log1p(adata_oh, "CXCL13")
    adata_oh.obs["sc_neotcr4"] = _scanpy_score_genes_harmonized(
        adata_oh, REGISTRY["NeoTCR4"].genes_up, "sc_neotcr4", REGISTRY["NeoTCR4"].genes_down
    )
    oh_sc_cols = []
    for s in sigs_cd4:
        cn = f"sc_{s}"
        adata_oh.obs[cn] = _scanpy_score_genes_harmonized(
            adata_oh, REGISTRY[s].genes_up, cn, REGISTRY[s].genes_down
        )
        oh_sc_cols.append(cn)
    adata_oh.obs["sc_composite"] = adata_oh.obs[oh_sc_cols].mean(axis=1)
    oh_cell_agg = (
        adata_oh.obs.groupby(["donor", "CDR3ab"], observed=True)[["raw_cxcl13", "sc_neotcr4", "sc_composite"]]
        .mean()
        .reset_index()
    )

    # Map Oh cellinfo author state annotations per (donor, CDR3ab)
    from tcrsift.loader import _pivot_vdj_by_barcode, load_cellranger_vdj

    oh_samples_root = blood_cd4_data_dir / "oh" / "samples"
    cell_states = []
    for sdir in sorted(oh_samples_root.glob("*_tumor")):
        cinfo = sdir / "cellinfo.csv.gz"
        if not cinfo.exists():
            continue
        ci = pd.read_csv(cinfo)
        ci["barcode"] = ci["index"].astype(str).str.split("-").str[0] + "-1"
        vdf_oh = load_cellranger_vdj(sdir / "vdj", sdir.name, verbose=False)
        piv_oh = _pivot_vdj_by_barcode(vdf_oh).reset_index()
        m_oh = piv_oh.merge(ci[["barcode", "cell_types"]], on="barcode", how="inner")
        m_oh = m_oh[m_oh["CDR3ab"].notna() & (m_oh["CDR3ab"] != "")].copy()
        m_oh["donor"] = sdir.name[: -len("_tumor")]
        m_oh = m_oh.rename(columns={"cell_types": "states"})
        cell_states.append(m_oh[["donor", "CDR3ab", "states"]])
    cs_df = pd.concat(cell_states, ignore_index=True)
    dom_state = (
        cs_df.groupby(["donor", "CDR3ab"])["states"]
        .agg(lambda s: s.value_counts().index[0])
        .rename("oh_state")
        .reset_index()
    )

    oh_clones = pd.read_csv(oh_out / "all_scored_clones.csv", low_memory=False)
    oh_clones["size_percentile"] = oh_clones.groupby("donor")["cell_count"].rank(pct=True)
    oh_clones = oh_clones.merge(oh_cell_agg, on=["donor", "CDR3ab"], how="left")
    oh_clones = oh_clones.merge(dom_state, on=["donor", "CDR3ab"], how="left")
    oh_clones["selected_for_review"] = truthy_flag_series(oh_clones["selected_for_review"])
    sel_tcr_oh = oh_clones["selected_for_review"].to_numpy(dtype=bool)
    k_oh = int(sel_tcr_oh.sum())
    exp_oh = oh_clones["cell_count"].ge(2).to_numpy()
    is_eff_oh = oh_clones["oh_state"].isin(["CD4CXCL13", "CD4GZMB", "CD4GZMK"]).to_numpy()
    is_cm_oh = oh_clones["oh_state"].isin(["CD4CM", "CD4TH17"]).to_numpy()
    is_treg_oh = oh_clones["oh_state"].eq("CD4IL2RAHI").to_numpy()
    for m_name, scol, sel_mask in (
        ("clone_size", "size_percentile", _topk_mask_among_expanded(oh_clones, "size_percentile", k_oh)),
        ("raw_CXCL13", "raw_cxcl13", _topk_mask_among_expanded(oh_clones, "raw_cxcl13", k_oh)),
        ("scanpy_score_genes_NeoTCR", "sc_neotcr4", _topk_mask_among_expanded(oh_clones, "sc_neotcr4", k_oh)),
        ("scanpy_composite_topk", "sc_composite", _topk_mask_among_expanded(oh_clones, "sc_composite", k_oh)),
        ("tcrsift_solid_tumor_cd4_excl_treg", "signature_TumorReactive_percentile_max", sel_tcr_oh),
    ):
        rows.append(
            _eval_method_row(
                "oh_bladder_cd4",
                "Effector_CXCL13_GZMB_K_vs_CD4CM_TH17",
                m_name,
                oh_clones[scol].to_numpy(),
                exp_oh & is_eff_oh,
                exp_oh & is_cm_oh,
                sel_mask,
                is_eff_oh,
                is_cm_oh,
                treg_mask=is_treg_oh,
            )
        )

    return rows


# ---------------------------------------------------------------------------
# 2. Full-Length TCR Construct Assembly: tcrsift vs stitchr (thimble v1.3.1)
# ---------------------------------------------------------------------------


def _run_thimble_batch(
    clones: pd.DataFrame,
    species: str,
    thimble_bin: Path,
    skip_c_checks: bool = False,
) -> tuple[pd.DataFrame, float]:
    cols = [
        "TCR_name",
        "TRAV",
        "TRAJ",
        "TRA_CDR3",
        "TRBV",
        "TRBJ",
        "TRB_CDR3",
        "TRAC",
        "TRBC",
        "TRA_leader",
        "TRB_leader",
        "Linker",
        "Link_order",
        "TRA_5_prime_seq",
        "TRA_3_prime_seq",
        "TRB_5_prime_seq",
        "TRB_3_prime_seq",
    ]
    in_rows = []
    for i, r in clones.reset_index(drop=True).iterrows():
        av = str(r.get("alpha_v_gene") or r.get("TRA_1_v_gene") or "").strip()
        aj = str(r.get("alpha_j_gene") or r.get("TRA_1_j_gene") or "").strip()
        bv = str(r.get("beta_v_gene") or r.get("TRB_1_v_gene") or "").strip()
        bj = str(r.get("beta_j_gene") or r.get("TRB_1_j_gene") or "").strip()
        ac = str(r.get("alpha_c_gene") or r.get("TRA_1_c_gene") or "").strip()
        bc = str(r.get("beta_c_gene") or r.get("TRB_1_c_gene") or "").strip()
        if ac in ("None", "nan", ""):
            ac = ""
        if bc in ("None", "nan", ""):
            bc = ""
        in_rows.append(
            {
                "TCR_name": f"clone_{i}",
                "TRAV": av,
                "TRAJ": aj,
                "TRA_CDR3": str(r["CDR3_alpha"]).strip(),
                "TRBV": bv,
                "TRBJ": bj,
                "TRB_CDR3": str(r["CDR3_beta"]).strip(),
                "TRAC": ac,
                "TRBC": bc,
                "TRA_leader": "",
                "TRB_leader": "",
                "Linker": "T2A",
                "Link_order": "BA",
                "TRA_5_prime_seq": "",
                "TRA_3_prime_seq": "",
                "TRB_5_prime_seq": "",
                "TRB_3_prime_seq": "",
            }
        )
    with tempfile.TemporaryDirectory() as td:
        in_tsv = Path(td) / "in.tsv"
        out_tsv = Path(td) / "out.tsv"
        pd.DataFrame(in_rows, columns=cols).to_csv(in_tsv, sep="\t", index=False)
        cmd = [str(thimble_bin), "-in", str(in_tsv), "-o", str(out_tsv), "-s", species.upper()]
        if skip_c_checks:
            cmd.append("-sc")
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
        elapsed = time.perf_counter() - t0
        if not out_tsv.exists():
            raise RuntimeError(f"thimble failed (rc={proc.returncode}): {proc.stderr}")
        out_df = pd.read_csv(out_tsv, sep="\t", keep_default_na=False)
    return out_df, elapsed


def run_assembly_vs_stitchr(
    inf_data_dir: Path,
    inf_work_dir: Path,
    mouse_work_dir: Path,
    cache_dir: Path,
    thimble_bin: Path,
) -> list[dict]:
    from tcrsift.assemble import (
        assemble_full_sequences,
        collect_germline_variants,
        validate_sequences,
    )
    from tcrsift.clonotype import aggregate_clonotypes
    from tcrsift.loader import load_samples
    from tcrsift.phenotype import apply_cd3_gate, phenotype_cells
    from tcrsift.sample_sheet import Sample, SampleSheet

    rows: list[dict] = []

    # --- 2A. Minervina 2022 Human Shortlist (200 clones with Cell Ranger 6 donor VDJ) ---
    min_cands = pd.read_csv(inf_work_dir / "minervina_blood_default" / "candidate_clones.csv", low_memory=False)
    t0 = time.perf_counter()
    tcr_min = assemble_full_sequences(
        min_cands,
        species="human",
        linker="T2A",
        data_dir=cache_dir,
        verbose=False,
        show_progress=False,
    )
    t_tcr_min = time.perf_counter() - t0
    val_min = [m for m in validate_sequences(tcr_min, species="human") if m.severity == "load_bearing"]
    gv_min = collect_germline_variants(tcr_min)
    st_min, t_st_min = _run_thimble_batch(min_cands, "HUMAN", thimble_bin)
    st_min_ok = st_min["Linked_aa"].astype(str).str.len() > 400

    # Compare VDJ amino-acid identity between donor-contig tcrsift and germline-imputed stitchr
    vdj_divergent = 0
    both_ok = 0
    for i in range(len(min_cands)):
        if pd.notna(tcr_min.loc[i, "single_chain_aa"]) and bool(st_min_ok.iloc[i]):
            both_ok += 1
            va = str(tcr_min.loc[i].get("VDJ_alpha_aa") or tcr_min.loc[i].get("vdj_alpha_aa") or "")
            vb = str(tcr_min.loc[i].get("VDJ_beta_aa") or tcr_min.loc[i].get("vdj_beta_aa") or "")
            s_linked = str(st_min.loc[i, "Linked_aa"])
            if (va and va not in s_linked) or (vb and vb not in s_linked):
                vdj_divergent += 1

    rows.append(
        {
            "cohort": "minervina_2022_human_cr6_vdj",
            "tool": "tcrsift_donor_contig",
            "n_clones": len(min_cands),
            "n_complete": int(tcr_min["single_chain_aa"].notna().sum()),
            "complete_rate": round(float(tcr_min["single_chain_aa"].notna().mean()), 4),
            "donor_vdj_polymorphic_clones_preserved": vdj_divergent,
            "mismatched_trbc_shipped": 0,
            "validation_errors": len(val_min),
            "distinct_germline_variants_logged": int(len(gv_min)),
            "elapsed_sec": round(t_tcr_min, 3),
        }
    )
    rows.append(
        {
            "cohort": "minervina_2022_human_cr6_vdj",
            "tool": "stitchr_thimble_v1.3.1",
            "n_clones": len(min_cands),
            "n_complete": int(st_min_ok.sum()),
            "complete_rate": round(float(st_min_ok.mean()), 4),
            "donor_vdj_polymorphic_clones_preserved": 0,
            "mismatched_trbc_shipped": 0,
            "validation_errors": 0,
            "distinct_germline_variants_logged": 0,
            "elapsed_sec": round(t_st_min, 3),
        }
    )

    # --- 2B. Fischer 2021 5'-Truncated Human Deposit (160 paired clones) ---
    fdir = inf_data_dir / "fischer_2021_covid"
    f_clones_all = []
    for donor, s_u, s_s in (("GT_6", "GT6_unstim", "GT6_covid_stim"), ("GT_2", "GT_unstim", "GT_covid_stim")):
        sh = SampleSheet(
            samples=[
                Sample(sample=s_u, patient_id=donor, gex_dir=str(fdir / s_u), vdj_dir=str(fdir / s_u), source="control"),
                Sample(sample=s_s, patient_id=donor, gex_dir=str(fdir / s_s), vdj_dir=str(fdir / s_s), source="culture", control_sample=s_u),
            ]
        )
        ad_r = load_samples(sh, min_mito_pct=0.0, max_mito_pct=15.0, show_progress=False, verbose=False)
        ad_p = apply_cd3_gate(phenotype_cells(ad_r, min_cd3_reads=3, verbose=False, show_progress=False), min_cd3_reads=3)
        f_clones_all.append(aggregate_clonotypes(ad_p, group_by="CDR3ab", show_progress=False, verbose=False))
    f_clones = pd.concat(f_clones_all, ignore_index=True)

    t0 = time.perf_counter()
    tcr_f_strict = assemble_full_sequences(
        f_clones.copy(), species="human", linker="T2A", germline_vdj_fallback=False, data_dir=cache_dir, verbose=False, show_progress=False
    )
    t_f_strict = time.perf_counter() - t0

    t0 = time.perf_counter()
    tcr_f_germ = assemble_full_sequences(
        f_clones.copy(),
        contigs_dir=fdir,
        sample_name_from="parent",
        alpha_leader="from_contig",
        beta_leader="from_contig",
        leader_fallback="germline",
        species="human",
        linker="T2A",
        germline_vdj_fallback=True,
        data_dir=cache_dir,
        verbose=False,
        show_progress=False,
    )
    t_f_germ = time.perf_counter() - t0
    val_f_germ = [m for m in validate_sequences(tcr_f_germ, species="human") if m.severity == "load_bearing"]
    partial_rescued_f = int(
        (
            tcr_f_germ["alpha_vdj_source"].isin(["contig_full", "contig_5prime_germline_rescued"])
            & tcr_f_germ["beta_vdj_source"].isin(["contig_full", "contig_5prime_germline_rescued"])
        ).sum()
    )

    st_f, t_st_f = _run_thimble_batch(f_clones, "HUMAN", thimble_bin)
    st_f_ok = st_f["Linked_aa"].astype(str).str.len() > 400

    rows.extend(
        [
            {
                "cohort": "fischer_2021_5prime_truncated_human",
                "tool": "tcrsift_donor_contig_only",
                "n_clones": len(f_clones),
                "n_complete": int(tcr_f_strict["single_chain_aa"].notna().sum()),
                "complete_rate": round(float(tcr_f_strict["single_chain_aa"].notna().mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": 0,
                "mismatched_trbc_shipped": 0,
                "validation_errors": 0,
                "distinct_germline_variants_logged": 0,
                "elapsed_sec": round(t_f_strict, 3),
            },
            {
                "cohort": "fischer_2021_5prime_truncated_human",
                "tool": "tcrsift_germline_vdj_fallback",
                "n_clones": len(f_clones),
                "n_complete": int(tcr_f_germ["single_chain_aa"].notna().sum()),
                "complete_rate": round(float(tcr_f_germ["single_chain_aa"].notna().mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": partial_rescued_f,
                "mismatched_trbc_shipped": 0,
                "validation_errors": len(val_f_germ),
                "distinct_germline_variants_logged": int(len(collect_germline_variants(tcr_f_germ))),
                "elapsed_sec": round(t_f_germ, 3),
            },
            {
                "cohort": "fischer_2021_5prime_truncated_human",
                "tool": "stitchr_thimble_v1.3.1",
                "n_clones": len(f_clones),
                "n_complete": int(st_f_ok.sum()),
                "complete_rate": round(float(st_f_ok.mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": 0,
                "mismatched_trbc_shipped": 0,
                "validation_errors": int((~st_f_ok).sum()),
                "distinct_germline_variants_logged": 0,
                "elapsed_sec": round(t_st_f, 3),
            },
        ]
    )

    # --- 2C. Murine 10x Cell Ranger Shortlists (101 clones: 41 MC38 CD8 + 60 Khatun CD4) ---
    m_tum = pd.read_csv(mouse_work_dir / "mc38_tumor_default" / "candidate_clones.csv", low_memory=False)
    m_bld = pd.read_csv(mouse_work_dir / "mc38_blood_default" / "candidate_clones.csv", low_memory=False)
    m_cd4 = pd.read_csv(mouse_work_dir / "khatun_cd4_presorted" / "candidate_clones.csv", low_memory=False)
    mouse_cands = pd.concat([m_tum, m_bld, m_cd4], ignore_index=True)

    t0 = time.perf_counter()
    tcr_mouse = assemble_full_sequences(
        mouse_cands.copy(),
        species="mouse",
        linker="T2A",
        germline_vdj_fallback=True,
        data_dir=cache_dir,
        verbose=False,
        show_progress=False,
    )
    t_tcr_mouse = time.perf_counter() - t0
    val_mouse = [m for m in validate_sequences(tcr_mouse, species="mouse") if m.severity == "load_bearing"]

    st_m_def, t_st_m_def = _run_thimble_batch(mouse_cands, "MOUSE", thimble_bin, skip_c_checks=False)
    st_m_def_ok = st_m_def["Linked_aa"].astype(str).str.len() > 400

    st_m_sc, t_st_m_sc = _run_thimble_batch(mouse_cands, "MOUSE", thimble_bin, skip_c_checks=True)
    st_m_sc_ok = st_m_sc["Linked_aa"].astype(str).str.len() > 400
    # Count how many skip_c_checks constructs shipped a mismatched Trbj2+Trbc1 or Trbj1+Trbc2 chimera
    mismatched_c = 0
    for i in range(len(st_m_sc)):
        if not bool(st_m_sc_ok.iloc[i]):
            continue
        bj_u = str(st_m_sc.loc[i, "TRBJ"]).upper()
        bc_u = str(st_m_sc.loc[i, "TRBC"]).upper()
        if ("TRBJ2" in bj_u and "TRBC1" in bc_u) or ("TRBJ1" in bj_u and "TRBC2" in bc_u):
            mismatched_c += 1

    rows.extend(
        [
            {
                "cohort": "mouse_cellranger_101_mc38_khatun",
                "tool": "tcrsift_mouse_germline_fallback",
                "n_clones": len(mouse_cands),
                "n_complete": int(tcr_mouse["single_chain_aa"].notna().sum()),
                "complete_rate": round(float(tcr_mouse["single_chain_aa"].notna().mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": 0,
                "mismatched_trbc_shipped": 0,
                "validation_errors": len(val_mouse),
                "distinct_germline_variants_logged": int(
                    tcr_mouse.get("beta_c_gene_overridden", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()
                ),
                "elapsed_sec": round(t_tcr_mouse, 3),
            },
            {
                "cohort": "mouse_cellranger_101_mc38_khatun",
                "tool": "stitchr_thimble_v1.3.1_default",
                "n_clones": len(mouse_cands),
                "n_complete": int(st_m_def_ok.sum()),
                "complete_rate": round(float(st_m_def_ok.mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": 0,
                "mismatched_trbc_shipped": 0,
                "validation_errors": int((~st_m_def_ok).sum()),
                "distinct_germline_variants_logged": 0,
                "elapsed_sec": round(t_st_m_def, 3),
            },
            {
                "cohort": "mouse_cellranger_101_mc38_khatun",
                "tool": "stitchr_thimble_v1.3.1_skip_c_checks",
                "n_clones": len(mouse_cands),
                "n_complete": int(st_m_sc_ok.sum()),
                "complete_rate": round(float(st_m_sc_ok.mean()), 4),
                "donor_vdj_polymorphic_clones_preserved": 0,
                "mismatched_trbc_shipped": mismatched_c,
                "validation_errors": int((~st_m_sc_ok).sum()),
                "distinct_germline_variants_logged": 0,
                "elapsed_sec": round(t_st_m_sc, 3),
            },
        ]
    )
    return rows


# ---------------------------------------------------------------------------
# 3. TCR Generation Probability & Publicness: tcrsift k-mer vs OLGA v1.3.0
# ---------------------------------------------------------------------------


def _load_olga_models(species: str):
    import olga.generation_probability as olga_pgen
    import olga.load_model as load_model

    olga_dir = Path(load_model.__file__).resolve().parent / "default_models"
    prefix = "mouse" if species == "mouse" else "human"

    mdir_b = olga_dir / f"{prefix}_T_beta"
    gd_b = load_model.GenomicDataVDJ()
    gd_b.load_igor_genomic_data(
        str(mdir_b / "model_params.txt"),
        str(mdir_b / "V_gene_CDR3_anchors.csv"),
        str(mdir_b / "J_gene_CDR3_anchors.csv"),
    )
    gm_b = load_model.GenerativeModelVDJ()
    gm_b.load_and_process_igor_model(str(mdir_b / "model_marginals.txt"))
    pg_b = olga_pgen.GenerationProbabilityVDJ(gm_b, gd_b)

    mdir_a = olga_dir / f"{prefix}_T_alpha"
    gd_a = load_model.GenomicDataVJ()
    gd_a.load_igor_genomic_data(
        str(mdir_a / "model_params.txt"),
        str(mdir_a / "V_gene_CDR3_anchors.csv"),
        str(mdir_a / "J_gene_CDR3_anchors.csv"),
    )
    gm_a = load_model.GenerativeModelVJ()
    gm_a.load_and_process_igor_model(str(mdir_a / "model_marginals.txt"))
    pg_a = olga_pgen.GenerationProbabilityVJ(gm_a, gd_a)
    return pg_a, pg_b


def run_pgen_vs_olga(inf_work_dir: Path, mouse_work_dir: Path) -> list[dict]:
    from tcrsift.annotate_tcrs import add_paired_ppost, add_pgen_ppost

    rows: list[dict] = []

    # --- 3A. Human peripheral blood repertoire (1,500 unique paired clones across 10x D1, D2 & Minervina) ---
    d1 = pd.read_csv(inf_work_dir / "tenx_donor1_blood_default" / "all_scored_clones.csv", low_memory=False)[
        ["donor", "CDR3_alpha", "CDR3_beta", "CDR3ab", "cell_count"]
    ]
    d2 = pd.read_csv(inf_work_dir / "tenx_donor2_blood_default" / "all_scored_clones.csv", low_memory=False)[
        ["donor", "CDR3_alpha", "CDR3_beta", "CDR3ab", "cell_count"]
    ]
    dm = pd.read_csv(inf_work_dir / "minervina_blood_default" / "all_scored_clones.csv", low_memory=False)[
        ["donor", "CDR3_alpha", "CDR3_beta", "CDR3ab", "cell_count"]
    ]
    all_h = pd.concat([d1, d2, dm], ignore_index=True)
    beta_donors_h = all_h.groupby("CDR3_beta")["donor"].nunique()
    all_h["is_public_beta"] = all_h["CDR3_beta"].map(beta_donors_h) >= 2

    # Stratified deterministic sample: all public CDR3_beta clones + top/random private up to 1,500
    uniq_h = all_h.sort_values(["is_public_beta", "cell_count"], ascending=[False, False], kind="mergesort").drop_duplicates("CDR3_beta")
    pub_h = uniq_h[uniq_h["is_public_beta"]]
    priv_h = uniq_h[~uniq_h["is_public_beta"]]
    rng = np.random.default_rng(42)
    take_priv = priv_h.iloc[rng.choice(len(priv_h), size=min(1500 - len(pub_h), len(priv_h)), replace=False)]
    sub_h = pd.concat([pub_h, take_priv], ignore_index=True)

    t0 = time.perf_counter()
    tcr_h = add_paired_ppost(add_pgen_ppost(sub_h[["CDR3_alpha", "CDR3_beta", "CDR3ab"]].copy(), backend="kmer", species="human"))
    t_tcr_h = time.perf_counter() - t0

    pg_a_h, pg_b_h = _load_olga_models("human")
    # Warm up Numba JIT once before timing
    _ = pg_b_h.compute_aa_CDR3_pgen("CASSIRSSYEQYF")
    _ = pg_a_h.compute_aa_CDR3_pgen("CAVSDGGSQGNLIF")

    t0 = time.perf_counter()
    olga_b_h = np.array([pg_b_h.compute_aa_CDR3_pgen(str(s)) for s in sub_h["CDR3_beta"]], dtype=float)
    olga_a_h = np.array([pg_a_h.compute_aa_CDR3_pgen(str(s)) for s in sub_h["CDR3_alpha"]], dtype=float)
    t_olga_h = time.perf_counter() - t0

    log_olga_b_h = np.where(olga_b_h > 0, np.log(np.maximum(olga_b_h, 1e-300)), -35.0)
    log_olga_a_h = np.where(olga_a_h > 0, np.log(np.maximum(olga_a_h, 1e-300)), -35.0)
    valid_b = olga_b_h > 0
    r_b_h, _ = pearsonr(tcr_h.loc[valid_b, "pgen_beta"], log_olga_b_h[valid_b])
    rho_b_h, _ = spearmanr(tcr_h.loc[valid_b, "pgen_beta"], log_olga_b_h[valid_b])
    valid_a = olga_a_h > 0
    r_a_h, _ = pearsonr(tcr_h.loc[valid_a, "pgen_alpha"], log_olga_a_h[valid_a])
    rho_a_h, _ = spearmanr(tcr_h.loc[valid_a, "pgen_alpha"], log_olga_a_h[valid_a])

    pub_mask_h = sub_h["is_public_beta"].to_numpy(dtype=bool)
    for model_name, vec, elapsed, r_val, rho_val in (
        ("tcrsift_kmer_pgen_beta", tcr_h["pgen_beta"].to_numpy(dtype=float), t_tcr_h, r_b_h, rho_b_h),
        ("tcrsift_kmer_ppost_beta", tcr_h["ppost_beta"].to_numpy(dtype=float), t_tcr_h, float("nan"), float("nan")),
        ("olga_exact_vdj_pgen_beta", log_olga_b_h, t_olga_h, 1.0, 1.0),
    ):
        auc = auroc(vec, pub_mask_h)
        lo, hi = auroc_ci(vec, pub_mask_h)
        rows.append(
            {
                "cohort": "human_blood_donors_1500",
                "model": model_name,
                "n_clones": len(sub_h),
                "n_public_ge2_donors": int(pub_mask_h.sum()),
                "pearson_r_vs_olga_beta": round(float(r_val), 4) if np.isfinite(r_val) else np.nan,
                "spearman_rho_vs_olga_beta": round(float(rho_val), 4) if np.isfinite(rho_val) else np.nan,
                "pearson_r_vs_olga_alpha": round(float(r_a_h), 4) if "pgen" in model_name and "tcrsift" in model_name else np.nan,
                "spearman_rho_vs_olga_alpha": round(float(rho_a_h), 4) if "pgen" in model_name and "tcrsift" in model_name else np.nan,
                "public_sharing_auroc": round(auc, 4),
                "auroc_low": round(lo, 4),
                "auroc_high": round(hi, 4),
                "elapsed_sec": round(elapsed, 4),
                "clones_per_sec": round(len(sub_h) / max(elapsed, 1e-6), 1),
            }
        )

    # --- 3B. Mouse Khatun 2021 LCMV GP66 CD4+ T cells (all 1,781 unique CDR3_beta across 5 mice M1-M5) ---
    mk = pd.read_csv(mouse_work_dir / "khatun_cd4_presorted" / "all_scored_clones.csv", low_memory=False)
    beta_mice = mk.groupby("CDR3_beta")["donor"].nunique()
    mk["is_public_mouse"] = mk["CDR3_beta"].map(beta_mice) >= 2
    uniq_m = mk.sort_values(["is_public_mouse", "cell_count"], ascending=[False, False], kind="mergesort").drop_duplicates("CDR3_beta").reset_index(drop=True)

    t0 = time.perf_counter()
    tcr_m = add_paired_ppost(add_pgen_ppost(uniq_m[["CDR3_alpha", "CDR3_beta", "CDR3ab"]].copy(), backend="kmer", species="mouse"))
    t_tcr_m = time.perf_counter() - t0

    pg_a_m, pg_b_m = _load_olga_models("mouse")
    _ = pg_b_m.compute_aa_CDR3_pgen("CASSDDRVNERLFF")
    _ = pg_a_m.compute_aa_CDR3_pgen("CAASDTNTGKLTF")

    t0 = time.perf_counter()
    olga_b_m = np.array([pg_b_m.compute_aa_CDR3_pgen(str(s)) for s in uniq_m["CDR3_beta"]], dtype=float)
    olga_a_m = np.array([pg_a_m.compute_aa_CDR3_pgen(str(s)) for s in uniq_m["CDR3_alpha"]], dtype=float)
    t_olga_m = time.perf_counter() - t0

    log_olga_b_m = np.where(olga_b_m > 0, np.log(np.maximum(olga_b_m, 1e-300)), -35.0)
    log_olga_a_m = np.where(olga_a_m > 0, np.log(np.maximum(olga_a_m, 1e-300)), -35.0)
    vb_m = olga_b_m > 0
    va_m = olga_a_m > 0
    r_b_m, _ = pearsonr(tcr_m.loc[vb_m, "pgen_beta"], log_olga_b_m[vb_m])
    rho_b_m, _ = spearmanr(tcr_m.loc[vb_m, "pgen_beta"], log_olga_b_m[vb_m])
    r_a_m, _ = pearsonr(tcr_m.loc[va_m, "pgen_alpha"], log_olga_a_m[va_m])
    rho_a_m, _ = spearmanr(tcr_m.loc[va_m, "pgen_alpha"], log_olga_a_m[va_m])

    pub_mask_m = uniq_m["is_public_mouse"].to_numpy(dtype=bool)
    for model_name, vec, elapsed, r_val, rho_val in (
        ("tcrsift_kmer_pgen_beta", tcr_m["pgen_beta"].to_numpy(dtype=float), t_tcr_m, r_b_m, rho_b_m),
        ("tcrsift_kmer_ppost_beta", tcr_m["ppost_beta"].to_numpy(dtype=float), t_tcr_m, float("nan"), float("nan")),
        ("olga_exact_vdj_pgen_beta", log_olga_b_m, t_olga_m, 1.0, 1.0),
    ):
        auc = auroc(vec, pub_mask_m)
        lo, hi = auroc_ci(vec, pub_mask_m)
        rows.append(
            {
                "cohort": "mouse_khatun_gp66_1781",
                "model": model_name,
                "n_clones": len(uniq_m),
                "n_public_ge2_donors": int(pub_mask_m.sum()),
                "pearson_r_vs_olga_beta": round(float(r_val), 4) if np.isfinite(r_val) else np.nan,
                "spearman_rho_vs_olga_beta": round(float(rho_val), 4) if np.isfinite(rho_val) else np.nan,
                "pearson_r_vs_olga_alpha": round(float(r_a_m), 4) if "pgen" in model_name and "tcrsift" in model_name else np.nan,
                "spearman_rho_vs_olga_alpha": round(float(rho_a_m), 4) if "pgen" in model_name and "tcrsift" in model_name else np.nan,
                "public_sharing_auroc": round(auc, 4),
                "auroc_low": round(lo, 4),
                "auroc_high": round(hi, 4),
                "elapsed_sec": round(elapsed, 4),
                "clones_per_sec": round(len(uniq_m) / max(elapsed, 1e-6), 1),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# 4. Multiplexed pMHC Dextramer Calling Head-to-Head (10x Healthy Donors 1 & 2)
# ---------------------------------------------------------------------------


def _fit_1d_two_gaussian_posterior(x: np.ndarray, max_iter: int = 40) -> np.ndarray:
    """Fit a 2-component 1D Gaussian mixture via EM on log1p(UMI) and return P(upper component)."""
    if len(x) < 10 or float(x.max()) <= np.log1p(5.0):
        return np.zeros(len(x), dtype=float)
    q95 = float(np.quantile(x, 0.98))
    m0 = float(np.median(x))
    m1 = max(q95, np.log1p(15.0))
    if m1 <= m0 + 0.5:
        return np.zeros(len(x), dtype=float)
    s0, s1 = 0.5, 0.8
    w1 = 0.05
    for _ in range(max_iter):
        p0 = (1.0 - w1) * np.exp(-0.5 * ((x - m0) / s0) ** 2) / (s0 * np.sqrt(2 * np.pi))
        p1 = w1 * np.exp(-0.5 * ((x - m1) / s1) ** 2) / (s1 * np.sqrt(2 * np.pi))
        denom = np.maximum(p0 + p1, 1e-300)
        gamma = p1 / denom
        n1 = float(gamma.sum())
        n0 = float(len(x) - n1)
        if n1 < 2.0 or n0 < 2.0:
            break
        w1 = np.clip(n1 / len(x), 1e-4, 0.5)
        m0 = float(((1.0 - gamma) * x).sum() / n0)
        m1 = max(float((gamma * x).sum() / n1), m0 + 0.8)
        s0 = max(float(np.sqrt(((1.0 - gamma) * (x - m0) ** 2).sum() / n0)), 0.25)
        s1 = max(float(np.sqrt((gamma * (x - m1) ** 2).sum() / n1)), 0.25)
    if m1 < np.log1p(10.0):
        return np.zeros(len(x), dtype=float)
    p0 = (1.0 - w1) * np.exp(-0.5 * ((x - m0) / s0) ** 2) / (s0 * np.sqrt(2 * np.pi))
    p1 = w1 * np.exp(-0.5 * ((x - m1) / s1) ** 2) / (s1 * np.sqrt(2 * np.pi))
    return p1 / np.maximum(p0 + p1, 1e-300)


def run_multimer_callers(inf_data_dir: Path, inf_work_dir: Path) -> list[dict]:
    from tcrsift.annotate import load_vdjdb, match_clonotypes
    from tcrsift.loader import load_cellranger_vdj
    from tcrsift.multimer import (
        call_multimer_specificity,
        classify_feature_columns,
        extract_feature_barcodes,
        summarize_clone_multimer_specificity,
    )

    heldout_10x = load_vdjdb(inf_work_dir / "vdjdb_heldout_no_10x.tsv", verbose=False)
    caller_stats: dict[str, dict[str, int]] = {
        k: {
            "total_called_clones_ge2": 0,
            "sticky_KLGGALQAK_fp_clones": 0,
            "hiv_gag_fp_clones": 0,
            "canonical_flu_ebv_cmv_clones": 0,
            "vdjdb_strict_matched_clones": 0,
            "vdjdb_cognate_epitope_concordant": 0,
        }
        for k in (
            "10x_per_cell_threshold",
            "clr_z3_per_cell",
            "gmm_bimodal_per_cell",
            "tcrsift_clonal_concordance",
        )
    }

    for donor in ("donor1", "donor2"):
        ddir = inf_data_dir / "tenx_dextramer" / donor
        fb_df = extract_feature_barcodes(ddir)
        fclasses = classify_feature_columns(fb_df)
        pmhc_cols = fclasses["multimers"]
        nc_cols = fclasses["negative_controls"]

        clones = pd.read_csv(inf_work_dir / f"tenx_{donor}_blood_default" / "all_scored_clones.csv", low_memory=False)
        exp_clones = clones[clones["cell_count"] >= 2].copy()
        exp_set = set(exp_clones["CDR3ab"])

        vdj = load_cellranger_vdj(ddir, donor, verbose=False)
        tra = vdj[vdj["chain"].eq("TRA") & vdj["cdr3"].notna()].sort_values("umis", ascending=False).drop_duplicates("barcode").set_index("barcode")["cdr3"]
        trb = vdj[vdj["chain"].eq("TRB") & vdj["cdr3"].notna()].sort_values("umis", ascending=False).drop_duplicates("barcode").set_index("barcode")["cdr3"]
        common_bcs = fb_df.index.intersection(tra.index).intersection(trb.index)
        cdr3ab_series = tra.loc[common_bcs] + "_" + trb.loc[common_bcs]
        keep_bcs = common_bcs[cdr3ab_series.isin(exp_set)]
        fb_sub = fb_df.loc[keep_bcs]
        clone_of_bc = cdr3ab_series.loc[keep_bcs]

        pmhc_mat = fb_sub[pmhc_cols].to_numpy(dtype=float)
        nc_max = fb_sub[nc_cols].to_numpy(dtype=float).max(axis=1) if nc_cols else np.zeros(len(fb_sub))

        # Caller 1: 10x per-cell threshold (`UMI >= 10` & `>= 5 * max(NC)`, single-binder, `>=1` cell in clone)
        cell_calls_10x = call_multimer_specificity(fb_sub, multimer_cols=pmhc_cols, nc_cols=nc_cols, min_umi=10.0, min_nc_ratio=5.0)
        c1_raw = pd.DataFrame({"CDR3ab": clone_of_bc.to_numpy(), "call": cell_calls_10x["multimer_call"].to_numpy()})
        c1_df = (
            c1_raw[c1_raw["call"].notna()]
            .groupby(["CDR3ab", "call"])
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
            .drop_duplicates("CDR3ab")
        )

        # Caller 2: CLR across multimers + z >= 3 per cell (`>=1` cell in clone)
        log_mat = np.log1p(pmhc_mat)
        clr_mat = log_mat - log_mat.mean(axis=1, keepdims=True)
        z_mat = (clr_mat - clr_mat.mean(axis=0, keepdims=True)) / np.maximum(clr_mat.std(axis=0, keepdims=True), 1e-6)
        pos_z = (z_mat >= 3.0) & (pmhc_mat >= 5.0)
        single_z = pos_z.sum(axis=1) == 1
        c2_calls = np.where(single_z, np.array(pmhc_cols, dtype=object)[pos_z.argmax(axis=1)], None)
        c2_raw = pd.DataFrame({"CDR3ab": clone_of_bc.to_numpy(), "call": c2_calls})
        c2_df = (
            c2_raw[c2_raw["call"].notna()]
            .groupby(["CDR3ab", "call"])
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
            .drop_duplicates("CDR3ab")
        )

        # Caller 3: 2-component bimodal GMM per multimer (`P >= 0.95` & `> 5x NC`, `>=1` cell in clone)
        post_mat = np.column_stack([_fit_1d_two_gaussian_posterior(log_mat[:, j]) for j in range(len(pmhc_cols))])
        pos_gmm = (post_mat >= 0.95) & (pmhc_mat >= 10.0) & (pmhc_mat > 5.0 * np.maximum(nc_max[:, None], 1.0))
        single_gmm = pos_gmm.sum(axis=1) == 1
        c3_calls = np.where(single_gmm, np.array(pmhc_cols, dtype=object)[pos_gmm.argmax(axis=1)], None)
        c3_raw = pd.DataFrame({"CDR3ab": clone_of_bc.to_numpy(), "call": c3_calls})
        c3_df = (
            c3_raw[c3_raw["call"].notna()]
            .groupby(["CDR3ab", "call"])
            .size()
            .reset_index(name="n")
            .sort_values("n", ascending=False)
            .drop_duplicates("CDR3ab")
        )

        # Caller 4: TCRsift clonal concordance (`min_specific_cells=2, min_concordance=0.5`)
        c4_sum = summarize_clone_multimer_specificity(
            cell_calls_10x.assign(CDR3ab=clone_of_bc.to_numpy()),
            clone_col="CDR3ab",
            min_concordance=0.5,
            min_specific_cells=2,
        )
        c4_df = c4_sum[c4_sum["multimer_confident"]].rename(columns={"multimer_call": "call"})[["CDR3ab", "call"]]

        # Annotate expanded clones against 10x-held-out VDJdb (`strict_ab`) to evaluate external precision
        ann = match_clonotypes(
            exp_clones[["CDR3_alpha", "CDR3_beta", "CDR3ab"]].drop_duplicates("CDR3ab"),
            heldout_10x,
            match_strictness="strict_ab",
            species="human",
            host_species_filter="same",
            verbose=False,
            show_progress=False,
        ).set_index("CDR3ab")

        for cname, cdf in (
            ("10x_per_cell_threshold", c1_df),
            ("clr_z3_per_cell", c2_df),
            ("gmm_bimodal_per_cell", c3_df),
            ("tcrsift_clonal_concordance", c4_df),
        ):
            st = caller_stats[cname]
            st["total_called_clones_ge2"] += len(cdf)
            calls_s = cdf["call"].astype(str)
            st["sticky_KLGGALQAK_fp_clones"] += int(calls_s.str.contains("KLGGALQAK").sum())
            st["hiv_gag_fp_clones"] += int(calls_s.str.contains("RTLNAWVKV|SLYNTVATL|Gag").sum())
            st["canonical_flu_ebv_cmv_clones"] += int(
                (
                    calls_s.str.contains("GILGFVFTL|GLCTLVAML|RAKFKQLL|IVTDFSVIK|AVFDRKSDAK|NLVPMVATV|TPRVTGGGAM")
                ).sum()
            )
            for _, r in cdf.iterrows():
                cab = r["CDR3ab"]
                if cab in ann.index and bool(ann.loc[cab, "db_match"]):
                    st["vdjdb_strict_matched_clones"] += 1
                    pep = str(r["call"]).split("_")[0]
                    db_ep = str(ann.loc[cab, "db_epitope"] or "")
                    if pep and pep in db_ep.split(";"):
                        st["vdjdb_cognate_epitope_concordant"] += 1

    out_rows = []
    for cname, st in caller_stats.items():
        n_m = st["vdjdb_strict_matched_clones"]
        n_ok = st["vdjdb_cognate_epitope_concordant"]
        out_rows.append(
            {
                "caller": cname,
                **st,
                "vdjdb_cognate_precision": round(n_ok / n_m, 4) if n_m > 0 else np.nan,
            }
        )
    return out_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--til-data-dir", type=Path, default=Path("/tmp/tcrsift_bench_til"))
    parser.add_argument("--til-work-dir", type=Path, default=Path("/tmp/tcrsift_bench_til_runs"))
    parser.add_argument("--inf-data-dir", type=Path, default=Path("/tmp/tcrsift_bench_inf"))
    parser.add_argument("--inf-work-dir", type=Path, default=Path("/tmp/tcrsift_bench_inf_runs"))
    parser.add_argument("--blood-cd4-data-dir", type=Path, default=Path("/tmp/tcrsift_bench_blood_cd4"))
    parser.add_argument("--blood-cd4-work-dir", type=Path, default=Path("/tmp/tcrsift_bench_blood_cd4/runs"))
    parser.add_argument("--mouse-data-dir", type=Path, default=Path("/tmp/tcrsift_bench_mouse"))
    parser.add_argument("--mouse-work-dir", type=Path, default=Path("/tmp/tcrsift_bench_mouse/runs"))
    parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/tcrsift_cache_vdjdb_germline"))
    parser.add_argument("--thimble-bin", type=Path, default=Path("/tmp/h2h_venv/bin/thimble"))
    parser.add_argument("--stage", choices=("all", "1", "2", "3", "4"), default="all")
    args = parser.parse_args()

    results_dir = HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    import olga
    import Stitchr.stitchrfunctions as st_fxn

    from tcrsift import __version__

    t0 = time.perf_counter()
    if args.stage in ("all", "1"):
        print("1/4 Running ex vivo prioritization head-to-head (tcrsift vs CXCL13 / Scanpy score_genes)...")
        prio_rows = run_prioritize_baselines(
            args.til_data_dir, args.til_work_dir, args.blood_cd4_data_dir, args.blood_cd4_work_dir
        )
        pd.DataFrame(prio_rows).to_csv(results_dir / "h2h_prioritize_baselines.csv", index=False)

    if args.stage in ("all", "2"):
        print("2/4 Running full-length TCR construct assembly head-to-head (tcrsift vs stitchr/thimble)...")
        asm_rows = run_assembly_vs_stitchr(
            args.inf_data_dir, args.inf_work_dir, args.mouse_work_dir, args.cache_dir, args.thimble_bin
        )
        pd.DataFrame(asm_rows).to_csv(results_dir / "h2h_assembly_vs_stitchr.csv", index=False)

    if args.stage in ("all", "3"):
        print("3/4 Running TCR Pgen/publicness head-to-head (tcrsift k-mer vs OLGA exact VDJ DP)...")
        pgen_rows = run_pgen_vs_olga(args.inf_work_dir, args.mouse_work_dir)
        pd.DataFrame(pgen_rows).to_csv(results_dir / "h2h_pgen_vs_olga.csv", index=False)

    if args.stage in ("all", "4"):
        print("4/4 Running pMHC multimer caller head-to-head (tcrsift clonal concordance vs per-cell/CLR/GMM)...")
        mult_rows = run_multimer_callers(args.inf_data_dir, args.inf_work_dir)
        pd.DataFrame(mult_rows).to_csv(results_dir / "h2h_multimer_callers.csv", index=False)

    import importlib.metadata

    prev_elapsed = None
    summary_path = results_dir / "summary.json"
    if args.stage != "all" and summary_path.exists():
        prev_elapsed = json.loads(summary_path.read_text()).get("elapsed_seconds")

    summary = {
        "tcrsift_version": __version__,
        "stitchr_version": getattr(st_fxn, "__version__", "1.3.1"),
        "olga_version": getattr(olga, "__version__", "1.3.0"),
        "scanpy_version": importlib.metadata.version("scanpy"),
        "elapsed_seconds": prev_elapsed if prev_elapsed is not None else round(time.perf_counter() - t0, 1),
        "tables": [
            "h2h_prioritize_baselines.csv",
            "h2h_assembly_vs_stitchr.csv",
            "h2h_pgen_vs_olga.csv",
            "h2h_multimer_callers.csv",
        ],
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Done in {summary['elapsed_seconds']}s -> {results_dir}")


if __name__ == "__main__":
    main()
