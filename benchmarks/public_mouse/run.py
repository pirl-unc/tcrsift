#!/usr/bin/env python3
"""Run the public murine (`--species mouse`) benchmark and write summary tables to ``results/``.

    python benchmarks/public_mouse/download.py --data-dir CACHE
    python benchmarks/public_mouse/run.py --data-dir CACHE --work-dir CACHE/runs

Exercises ``tcrsift`` CLI and library workflows across three open-access murine
single-cell VDJ + GEX cohorts (19 paired 10x libraries):

1. Multi-compartment MC38 colon adenocarcinoma (GSE168944; Mouse01 & Mouse02 across
   matched Tumor, Peripheral Blood, Draining Lymph Node, and Spleen):
   - Unmodified ``tcrsift prioritize --species mouse --context solid-tumor`` on
     MC38 tumors and ``--context blood`` on matched peripheral blood.
   - Multi-compartment clonal tracking (tumor-enriched vs lymphoid-shared clones
     in tumor; tumor-infiltrating vs blood-restricted clones in peripheral blood).
   - Direct ablation of MGI species-adapted ``TumorReactive`` (``Ccl3``/``Ccl4``)
     against naive human ``CXCL13``-containing ``TumorReactive`` on mouse CD8 TILs.
   - Full-length murine TCR construct assembly (``Trac``/``Trbc1``/``Trbc2``) via
     ``assemble_full_sequences(..., species="mouse", germline_vdj_fallback=True)``.
2. Kasmani et al. 2022 (PMID:35667687; GSE200608 — NOD mouse pancreatic islet and
   spleen H-2Kd IGRP206-214 / NRP-V7 ``KYNKANVFL`` tetramer+ vs islet tetramer-
   CD8+ T cells across Rounds 1 & 2):
   - Same cognate epitope across compartments (target-organ Islet Tet+ vs
     circulating Spleen Tet+, including exact sister ``CDR3ab`` clones) and
     within-islet Tet+ vs Tet- CD8+ T cells.
   - VDJdb ``MusMusculus`` specificity lookup comparing unmodified VDJdb
     (**circular**: contains 780 entries from PMID:35667687) against study-held-out
     VDJdb across ``strict_ab``, ``ab_with_partial``, ``b_only``, and
     ``b_only + levenshtein`` (distance <= 1), with C57BL/6 MC38 clones as
     out-of-strain specificity controls.
3. Khatun et al. 2021 (PMID:33201171; GSE158896 — C57BL/6 Day-10 acute LCMV
   Armstrong I-Ab GP66-77 ``DIYKGVYQFKSV`` tetramer+ splenic CD4+ T cells across
   5 biological replicate mice M1-M5):
   - Unmodified ``tcrsift prioritize --species mouse --tcell-type cd4 --context blood``.
   - Murine CD4 phenotyping (with and without ``pre_sorted=CD4``), ``Foxp3`` Treg
     audit, cross-mouse public ``GP66`` receptor sharing, VDJdb ``LCMV`` lookup,
     and full-length murine CD4 TCR construct assembly.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent


def is_current(out: Path) -> bool:
    from tcrsift import __version__

    cfg = out / "prioritization.json"
    return cfg.exists() and json.loads(cfg.read_text()).get("version") == __version__


def auroc(scores, positive) -> float:
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
    ok = np.isfinite(scores)
    scores, positive = scores[ok], positive[ok]
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if not n_pos or not n_neg:
        return float("nan")
    ranks = rankdata(scores)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auroc_ci(scores, positive, draws: int = 2000, seed: int = 0) -> tuple[float, float]:
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
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
    low, high = np.quantile(vals, [0.025, 0.975])
    return float(low), float(high)


def eval_auroc_block(
    df: pd.DataFrame,
    pos_mask: np.ndarray,
    signatures: list[str],
    *,
    col_fmt: str = "signature_{}_percentile_max",
    size_col: str = "size_percentile",
    **tags,
) -> list[dict]:
    rows = []
    pos_mask = np.asarray(pos_mask, dtype=bool)
    for sig in [*signatures, "clone_size"]:
        col = size_col if sig == "clone_size" else col_fmt.format(sig)
        if col not in df.columns:
            continue
        vals = df[col].to_numpy(dtype=float)
        ok = np.isfinite(vals)
        est = auroc(vals[ok], pos_mask[ok])
        low, high = auroc_ci(vals[ok], pos_mask[ok])
        rows.append(
            {
                **tags,
                "signature": sig,
                "n_positive": int(pos_mask[ok].sum()),
                "n_negative": int((~pos_mask[ok]).sum()),
                "auroc": round(est, 4) if np.isfinite(est) else np.nan,
                "auroc_low": round(low, 4) if np.isfinite(low) else np.nan,
                "auroc_high": round(high, 4) if np.isfinite(high) else np.nan,
            }
        )
    return rows


def _paired_vdj_counts_by_sample(vdj_dir: Path, sample_name: str) -> pd.DataFrame:
    """Load productive paired alpha/beta CDR3ab cell counts from a 10x VDJ folder."""
    from tcrsift.loader import _pivot_vdj_by_barcode, load_cellranger_vdj

    vdj_df = load_cellranger_vdj(vdj_dir, sample_name, verbose=False)
    piv = _pivot_vdj_by_barcode(vdj_df)
    paired = piv[
        piv["CDR3_alpha"].notna()
        & (piv["CDR3_alpha"] != "")
        & piv["CDR3_beta"].notna()
        & (piv["CDR3_beta"] != "")
    ].copy()
    counts = paired.groupby("CDR3ab").size().rename("cells").reset_index()
    total = max(int(counts["cells"].sum()), 1)
    counts["freq"] = counts["cells"] / total
    return counts


# ---------------------------------------------------- 1. MC38 (GSE168944)


def run_mc38(data_dir: Path, work_dir: Path, cache_dir: Path | None):
    from tcrsift.assemble import assemble_full_sequences, validate_sequences
    from tcrsift.cli import main as tcrsift_main
    from tcrsift.genes import adata_symbol_array
    from tcrsift.loader import load_samples
    from tcrsift.phenotype import phenotype_cells
    from tcrsift.signature_methods import SIGNATURES, score_signature

    samples_root = data_dir / "mc38" / "samples"
    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []

    # Build multi-tissue paired CDR3ab repertoire map across Tumor, Blood, Spleen, dLN
    comp_frames = []
    for mouse in ("Mouse01", "Mouse02"):
        for tissue in ("Tumor", "Blood", "Spleen", "dLN"):
            sname = f"{mouse}_{tissue}"
            c_df = _paired_vdj_counts_by_sample(samples_root / sname / "vdj", sname)
            c_df["donor"] = mouse
            c_df["tissue"] = tissue.lower()
            comp_frames.append(c_df)
    comp_long = pd.concat(comp_frames, ignore_index=True)
    comp_wide = (
        comp_long.pivot_table(
            index=["donor", "CDR3ab"],
            columns="tissue",
            values=["cells", "freq"],
            fill_value=0.0,
        )
    )
    comp_wide.columns = [f"{metric}_{tis}" for metric, tis in comp_wide.columns]
    comp_wide = comp_wide.reset_index()
    comp_wide["cells_lymphoid"] = comp_wide["cells_spleen"] + comp_wide["cells_dln"]
    lymph_tot = comp_wide.groupby("donor")["cells_lymphoid"].transform("sum").clip(lower=1.0)
    comp_wide["freq_lymphoid"] = comp_wide["cells_lymphoid"] / lymph_tot

    # 1a. Default tumor prioritization (`--species mouse --context solid-tumor`)
    tumor_sheet = work_dir / "mc38_tumor_sheet.csv"
    tumor_sheet.write_text(
        "sample,patient_id,tissue,vdj_dir,gex_dir\n"
        f"Mouse01_Tumor,Mouse01,tumor,{samples_root / 'Mouse01_Tumor/vdj'},{samples_root / 'Mouse01_Tumor/gex'}\n"
        f"Mouse02_Tumor,Mouse02,tumor,{samples_root / 'Mouse02_Tumor/vdj'},{samples_root / 'Mouse02_Tumor/gex'}\n"
    )
    out_tum = work_dir / "mc38_tumor_default"
    if not is_current(out_tum):
        tcrsift_main(
            [
                "prioritize",
                str(tumor_sheet),
                "-o",
                str(out_tum),
                "--species",
                "mouse",
                "--context",
                "solid-tumor",
                *cache_args,
            ]
        )
    cfg_tum = json.loads((out_tum / "prioritization.json").read_text())
    sigs_tum = cfg_tum["resolved_signatures"]
    clones_tum = pd.read_csv(out_tum / "all_scored_clones.csv", low_memory=False)
    clones_tum["size_percentile"] = clones_tum.groupby("donor")["cell_count"].rank(pct=True)
    clones_tum = clones_tum.merge(comp_wide, on=["donor", "CDR3ab"], how="left")
    for col in ("cells_blood", "cells_spleen", "cells_dln", "cells_lymphoid", "freq_lymphoid"):
        clones_tum[col] = clones_tum[col].fillna(0.0)

    # MGI vs naive Human CXCL13-containing TumorReactive ablation on retained MC38 tumor CD8 cells
    from tcrsift.prioritize import PrioritizeConfig, filter_cells

    adata_tum = load_samples(tumor_sheet, min_mito_pct=0.0, verbose=False, show_progress=False)
    adata_tum = phenotype_cells(adata_tum, min_cd3_reads=3, verbose=False, show_progress=False)
    adata_cd8, _, _ = filter_cells(
        adata_tum,
        PrioritizeConfig(
            sample_sheet=tumor_sheet,
            output_dir=out_tum,
            species="mouse",
            context="solid-tumor",
        ),
    )
    syms_upper = adata_symbol_array(adata_cd8)
    cxcl13_idx = np.where(syms_upper == "CXCL13")[0]
    if cxcl13_idx.size > 0:
        cx_raw = adata_cd8[:, cxcl13_idx[0]].X
        cx_arr = np.asarray(cx_raw.todense()).ravel() if hasattr(cx_raw, "todense") else np.asarray(cx_raw).ravel()
        cxcl13_pos_pct = float((cx_arr > 0).mean())
    else:
        cxcl13_pos_pct = 0.0
    ccl3_idx = np.where(syms_upper == "CCL3")[0]
    ccl4_idx = np.where(syms_upper == "CCL4")[0]
    ccl3_pos_pct = (
        float((np.asarray(adata_cd8[:, ccl3_idx[0]].X.todense()).ravel() > 0).mean())
        if ccl3_idx.size > 0
        else 0.0
    )
    ccl4_pos_pct = (
        float((np.asarray(adata_cd8[:, ccl4_idx[0]].X.todense()).ravel() > 0).mean())
        if ccl4_idx.size > 0
        else 0.0
    )

    sc.pp.normalize_total(adata_cd8, target_sum=10_000)
    sc.pp.log1p(adata_cd8)
    from tcrsift.signature_methods import expression_frame_from_adata

    human_tr_sig = SIGNATURES["TumorReactive"]  # includes CXCL13
    expr_htr = expression_frame_from_adata(adata_cd8, human_tr_sig.all_genes, on_missing="ignore")
    adata_cd8.obs["sig_human_tr"] = score_signature(
        expr_htr,
        human_tr_sig,
        log1p=False,
        groups=adata_cd8.obs["sample"],
        on_missing="ignore",
    ).to_numpy(dtype=float)
    htr_clone = (
        adata_cd8.obs.assign(donor=adata_cd8.obs["patient_id"].astype(str))
        .groupby(["donor", "sample", "CDR3ab"], observed=True)["sig_human_tr"]
        .mean()
        .reset_index()
    )
    htr_clone["signature_TumorReactive_HumanCXCL13_percentile"] = htr_clone.groupby(
        ["donor", "sample"], observed=True
    )["sig_human_tr"].rank(method="average", pct=True)
    htr_max = (
        htr_clone.groupby(["donor", "CDR3ab"], observed=True)[
            "signature_TumorReactive_HumanCXCL13_percentile"
        ]
        .max()
        .rename("signature_TumorReactive_HumanCXCL13_percentile_max")
        .reset_index()
    )
    clones_tum = clones_tum.merge(htr_max, on=["donor", "CDR3ab"], how="left")

    # Define multi-compartment contrasts:
    # (A) Between-organ MC38 Tumor (>=2 cells) vs matched dLN + Spleen non-tumor clones (>=2 cells, 0 tumor cells)
    #     scored on a shared per-mouse gene baseline (testing MGI TumorReactive vs naive human CXCL13 TumorReactive!).
    organ_sheet = work_dir / "mc38_organs_sheet.csv"
    organ_lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    for mouse in ("Mouse01", "Mouse02"):
        for tis in ("Tumor", "dLN", "Spleen"):
            sname = f"{mouse}_{tis}"
            organ_lines.append(
                f"{sname},{mouse},{tis.lower()},{samples_root / sname / 'vdj'},{samples_root / sname / 'gex'}"
            )
    organ_sheet.write_text("\n".join(organ_lines) + "\n")
    adata_org = load_samples(organ_sheet, min_mito_pct=0.0, verbose=False, show_progress=False)
    adata_org = phenotype_cells(adata_org, min_cd3_reads=3, verbose=False, show_progress=False)
    keep_org = (
        adata_org.obs["is_CD8"].to_numpy(dtype=bool)
        & adata_org.obs["filter:min_cd3"].to_numpy(dtype=bool)
        & adata_org.obs["CDR3ab"].notna().to_numpy()
        & pd.to_numeric(adata_org.obs["TRA_1_umis"], errors="coerce").ge(2).to_numpy()
        & pd.to_numeric(adata_org.obs["TRB_1_umis"], errors="coerce").ge(2).to_numpy()
    )
    adata_org_cd8 = adata_org[keep_org].copy()
    sc.pp.normalize_total(adata_org_cd8, target_sum=10_000)
    sc.pp.log1p(adata_org_cd8)
    adata_org_cd8.uns["species"] = "mouse"
    from tcrsift.signature_methods import score_by_name

    for sig in sigs_tum:
        adata_org_cd8.obs[f"org_{sig}"] = score_by_name(
            adata_org_cd8,
            sig,
            log1p=False,
            groups=adata_org_cd8.obs["patient_id"],
            species="mouse",
        ).to_numpy(dtype=float)
    expr_org_htr = expression_frame_from_adata(
        adata_org_cd8, human_tr_sig.all_genes, on_missing="ignore"
    )
    adata_org_cd8.obs["org_TumorReactive_HumanCXCL13"] = score_signature(
        expr_org_htr,
        human_tr_sig,
        log1p=False,
        groups=adata_org_cd8.obs["patient_id"],
        on_missing="ignore",
    ).to_numpy(dtype=float)

    org_sigs = [*sigs_tum, "TumorReactive_HumanCXCL13"]
    org_obs = adata_org_cd8.obs[
        ["patient_id", "sample", "tissue", "CDR3ab", *[f"org_{s}" for s in org_sigs]]
    ].copy()
    org_obs["donor"] = org_obs["patient_id"].astype(str)
    org_obs["compartment_group"] = np.where(org_obs["tissue"].eq("tumor"), "tumor", "lymphoid")
    org_grp = org_obs.groupby(["donor", "compartment_group", "CDR3ab"], observed=True)
    org_clones = (
        org_grp[[f"org_{s}" for s in org_sigs]].mean().join(org_grp.size().rename("cells")).reset_index()
    )
    org_clones = org_clones.merge(comp_wide[["donor", "CDR3ab", "cells_tumor"]], on=["donor", "CDR3ab"], how="left")
    # Disjoint expanded (>=2 cells) clones: Tumor clones vs non-tumor lymphoid (dLN/Spleen, 0 tumor cells) clones
    org_disjoint = org_clones[
        (org_clones["cells"] >= 2)
        & (
            org_clones["compartment_group"].eq("tumor")
            | (org_clones["cells_tumor"].fillna(0) == 0)
        )
    ].copy()
    org_disjoint["size_percentile"] = org_disjoint.groupby("donor", observed=True)["cells"].rank(pct=True)
    for s in org_sigs:
        org_disjoint[f"signature_{s}_org_pct"] = org_disjoint.groupby("donor", observed=True)[
            f"org_{s}"
        ].rank(method="average", pct=True)
    auroc_rows = eval_auroc_block(
        org_disjoint,
        org_disjoint["compartment_group"].eq("tumor").to_numpy(),
        org_sigs,
        col_fmt="signature_{}_org_pct",
        dataset="mc38",
        compartment="tumor_vs_lymphoid",
        contrast="tumor_til_vs_dln_spleen_nontumor",
    )

    # (B) Within-tumor exclusivity contrast among >=2-cell tumor clones:
    # Positives: tumor-exclusive expanded (`cell_count >= 5` and `cells_lymphoid == 0`)
    # Negatives: lymphoid-shared non-enriched (`cells_lymphoid >= 1` and `max_frequency <= 2 * freq_lymphoid`)
    exp_tum = clones_tum[clones_tum["cell_count"] >= 2].copy()
    is_tum_enriched = (exp_tum["cell_count"] >= 5) & (exp_tum["cells_lymphoid"] == 0)
    is_lymph_bystander = (exp_tum["cells_lymphoid"] >= 1) & (
        exp_tum["max_frequency"] <= 2.0 * exp_tum["freq_lymphoid"]
    )
    sub_contrast = exp_tum[is_tum_enriched | is_lymph_bystander].copy()
    pos_tum = is_tum_enriched.loc[sub_contrast.index].to_numpy()
    auroc_rows += eval_auroc_block(
        sub_contrast,
        pos_tum,
        org_sigs,
        dataset="mc38",
        compartment="within_tumor",
        contrast="tumor_exclusive_ge5_vs_lymphoid_shared",
    )

    # 1b. Peripheral blood prioritization (`--species mouse --context blood`)
    blood_sheet = work_dir / "mc38_blood_sheet.csv"
    blood_sheet.write_text(
        "sample,patient_id,tissue,vdj_dir,gex_dir\n"
        f"Mouse01_Blood,Mouse01,blood,{samples_root / 'Mouse01_Blood/vdj'},{samples_root / 'Mouse01_Blood/gex'}\n"
        f"Mouse02_Blood,Mouse02,blood,{samples_root / 'Mouse02_Blood/vdj'},{samples_root / 'Mouse02_Blood/gex'}\n"
    )
    out_bld = work_dir / "mc38_blood_default"
    if not is_current(out_bld):
        tcrsift_main(
            [
                "prioritize",
                str(blood_sheet),
                "-o",
                str(out_bld),
                "--species",
                "mouse",
                "--context",
                "blood",
                *cache_args,
            ]
        )
    cfg_bld = json.loads((out_bld / "prioritization.json").read_text())
    sigs_bld = cfg_bld["resolved_signatures"]
    clones_bld = pd.read_csv(out_bld / "all_scored_clones.csv", low_memory=False)
    clones_bld["size_percentile"] = clones_bld.groupby("donor")["cell_count"].rank(pct=True)
    clones_bld = clones_bld.merge(comp_wide, on=["donor", "CDR3ab"], how="left")
    clones_bld["cells_tumor"] = clones_bld["cells_tumor"].fillna(0.0)
    exp_bld = clones_bld[clones_bld["cell_count"] >= 2].copy()
    in_tum_exp = exp_bld["cells_tumor"] >= 5
    not_in_tum = exp_bld["cells_tumor"] == 0
    sub_bld = exp_bld[in_tum_exp | not_in_tum].copy()
    auroc_rows += eval_auroc_block(
        sub_bld,
        in_tum_exp.loc[sub_bld.index].to_numpy(),
        sigs_bld,
        dataset="mc38",
        compartment="blood",
        contrast="tumor_infiltrating_ge5_vs_blood_only",
    )

    # 1c. Murine full-length construct assembly on shortlisted MC38 tumor + blood candidates
    cands_tum = pd.read_csv(out_tum / "candidate_clones.csv", low_memory=False)
    cands_bld = pd.read_csv(out_bld / "candidate_clones.csv", low_memory=False)
    asm_input = pd.concat(
        [cands_tum.assign(cohort_source="mc38_tumor"), cands_bld.assign(cohort_source="mc38_blood")],
        ignore_index=True,
    )
    asm_out = assemble_full_sequences(
        asm_input,
        species="mouse",
        germline_vdj_fallback=True,
        data_dir=cache_dir,
        verbose=False,
        show_progress=False,
    )
    val_msgs = validate_sequences(asm_out, species="mouse")
    val_errors = [m for m in val_msgs if m.severity == "load_bearing"]

    shortlist_rows = [
        {
            "dataset": "mc38_tumor",
            "setting": "solid-tumor default",
            "retained_cells": int(cfg_tum["cell_filters"][-1]["after"]),
            "total_clones": len(clones_tum),
            "expanded_ge2": int((clones_tum["cell_count"] >= 2).sum()),
            "shortlisted": int(clones_tum["selected_for_review"].sum()),
            "tumor_enriched_total": int(is_tum_enriched.sum()),
            "tumor_enriched_shortlisted": int(
                exp_tum.loc[is_tum_enriched, "selected_for_review"].sum()
            ),
            "lymphoid_shared_total": int(is_lymph_bystander.sum()),
            "lymphoid_shared_shortlisted": int(
                exp_tum.loc[is_lymph_bystander, "selected_for_review"].sum()
            ),
        },
        {
            "dataset": "mc38_blood",
            "setting": "blood default",
            "retained_cells": int(cfg_bld["cell_filters"][-1]["after"]),
            "total_clones": len(clones_bld),
            "expanded_ge2": int((clones_bld["cell_count"] >= 2).sum()),
            "shortlisted": int(clones_bld["selected_for_review"].sum()),
            "tumor_enriched_total": int((exp_bld["cells_tumor"] >= 5).sum()),
            "tumor_enriched_shortlisted": int(
                exp_bld.loc[exp_bld["cells_tumor"] >= 5, "selected_for_review"].sum()
            ),
            "lymphoid_shared_total": int((exp_bld["cells_tumor"] == 0).sum()),
            "lymphoid_shared_shortlisted": int(
                exp_bld.loc[exp_bld["cells_tumor"] == 0, "selected_for_review"].sum()
            ),
        },
    ]

    facts = {
        "tumor_retained_cd8_cells": int(cfg_tum["cell_filters"][-1]["after"]),
        "tumor_clones": len(clones_tum),
        "tumor_expanded_clones": int((clones_tum["cell_count"] >= 2).sum()),
        "tumor_shortlisted": int(clones_tum["selected_for_review"].sum()),
        "blood_retained_cd8_cells": int(cfg_bld["cell_filters"][-1]["after"]),
        "blood_clones": len(clones_bld),
        "blood_expanded_clones": int((clones_bld["cell_count"] >= 2).sum()),
        "blood_shortlisted": int(clones_bld["selected_for_review"].sum()),
        "mouse_cd8_cxcl13_pos_cell_fraction": round(cxcl13_pos_pct, 5),
        "mouse_cd8_ccl3_pos_cell_fraction": round(ccl3_pos_pct, 4),
        "mouse_cd8_ccl4_pos_cell_fraction": round(ccl4_pos_pct, 4),
        "assembled_candidates": len(asm_out),
        "assembled_complete": int(asm_out["single_chain_aa"].notna().sum()),
        "assembly_validation_errors": len(val_errors),
    }
    return auroc_rows, shortlist_rows, asm_out, clones_tum, facts


# ---------------------------------------------------- 2. Kasmani 2022 (GSE200608)


def run_kasmani(
    data_dir: Path,
    work_dir: Path,
    cache_dir: Path | None,
    vdjdb_full: pd.DataFrame,
    mc38_clones: pd.DataFrame,
):
    from tcrsift.annotate import match_clonotypes
    from tcrsift.cli import main as tcrsift_main

    sroot = data_dir / "kasmani" / "samples"
    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []

    samples_meta = [
        ("R1_Islet_TetNeg", "Round1", "islet_tet_neg"),
        ("R1_Islet_TetPos", "Round1", "islet_tet_pos"),
        ("R1_Spleen_TetPos", "Round1", "spleen_tet_pos"),
        ("R2_Islet_TetNeg", "Round2", "islet_tet_neg"),
        ("R2_Islet_TetPos", "Round2", "islet_tet_pos"),
        ("R2_Spleen_TetPos", "Round2", "spleen_tet_pos"),
    ]
    # Run each sample as its own stratum while grouping donor by biological round;
    # ALSO score raw per-cell expression across pooled Round1+Round2 samples so
    # between-compartment Islet Tet+ vs Spleen Tet+ differences are directly measured
    # on a shared normalization/z-score baseline (within-sample percentiles have mean
    # 0.50 in every sample by construction!).
    sheet = work_dir / "kasmani_sheet.csv"
    lines = ["sample,patient_id,tissue,pre_sorted,vdj_dir,gex_dir"]
    for sname, rnd, tis in samples_meta:
        lines.append(f"{sname},{rnd},{tis},CD8,{sroot / sname / 'vdj'},{sroot / sname / 'gex'}")
    sheet.write_text("\n".join(lines) + "\n")

    out_kas = work_dir / "kasmani_default"
    if not is_current(out_kas):
        tcrsift_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_kas),
                "--species",
                "mouse",
                "--context",
                "solid-tumor",
                "--signatures",
                "TumorReactive",
                "Differentiated",
                "AntigenExperienced",
                "Cytolytic",
                "CirculatingMemory",
                "AcuteActivation",
                "Proliferation",
                "AIM",
                *cache_args,
            ]
        )
    cfg = json.loads((out_kas / "prioritization.json").read_text())
    sigs = cfg["resolved_signatures"]

    # Pooled within-round scoring across the 3 compartments per round so Islet_TetPos,
    # Spleen_TetPos, and Islet_TetNeg share a single per-round gene baseline!
    from tcrsift.loader import load_samples
    from tcrsift.phenotype import phenotype_cells
    from tcrsift.signature_methods import score_by_name

    adata = load_samples(sheet, min_mito_pct=0.0, verbose=False, show_progress=False)
    adata = phenotype_cells(adata, min_cd3_reads=3, verbose=False, show_progress=False)
    keep = (
        adata.obs["is_CD8"].to_numpy(dtype=bool)
        & adata.obs["filter:min_cd3"].to_numpy(dtype=bool)
        & adata.obs["CDR3ab"].notna().to_numpy()
        & (pd.to_numeric(adata.obs["TRA_1_umis"], errors="coerce").ge(2).to_numpy())
        & (pd.to_numeric(adata.obs["TRB_1_umis"], errors="coerce").ge(2).to_numpy())
    )
    adata_cd8 = adata[keep].copy()
    sc.pp.normalize_total(adata_cd8, target_sum=10_000)
    sc.pp.log1p(adata_cd8)
    adata_cd8.uns["species"] = "mouse"

    for sig in sigs:
        adata_cd8.obs[f"pooled_{sig}"] = score_by_name(
            adata_cd8,
            sig,
            log1p=False,
            groups=adata_cd8.obs["patient_id"],
            species="mouse",
        ).to_numpy(dtype=float)

    obs_df = adata_cd8.obs[
        ["patient_id", "sample", "tissue", "CDR3ab", "CDR3_alpha", "CDR3_beta", *[f"pooled_{s}" for s in sigs]]
    ].copy()
    obs_df["donor"] = obs_df["patient_id"].astype(str)
    grp = obs_df.groupby(["donor", "sample", "tissue", "CDR3ab", "CDR3_alpha", "CDR3_beta"], observed=True)
    comp_clones = grp[[f"pooled_{s}" for s in sigs]].mean().join(grp.size().rename("cells")).reset_index()
    comp_clones["size_percentile"] = comp_clones.groupby(["donor", "sample"], observed=True)["cells"].rank(pct=True)
    for sig in sigs:
        comp_clones[f"signature_{sig}_pooled_pct"] = comp_clones.groupby("donor", observed=True)[
            f"pooled_{sig}"
        ].rank(method="average", pct=True)

    exp_comp = comp_clones[comp_clones["cells"] >= 2].copy()
    auroc_rows: list[dict] = []

    # Contrast A: Target-organ Islet IGRP Tet+ vs Splenic circulating IGRP Tet+ (same KYNKANVFL epitope!)
    sub_isl_spl = exp_comp[exp_comp["tissue"].isin(["islet_tet_pos", "spleen_tet_pos"])].copy()
    auroc_rows += eval_auroc_block(
        sub_isl_spl,
        sub_isl_spl["tissue"].eq("islet_tet_pos").to_numpy(),
        sigs,
        col_fmt="signature_{}_pooled_pct",
        dataset="kasmani",
        compartment="islet_vs_spleen_igrp",
        contrast="islet_tet_pos_vs_spleen_tet_pos",
    )

    # Contrast B: Within-islet IGRP Tet+ vs Islet Tet- CD8 T cells
    sub_isl_pn = exp_comp[exp_comp["tissue"].isin(["islet_tet_pos", "islet_tet_neg"])].copy()
    auroc_rows += eval_auroc_block(
        sub_isl_pn,
        sub_isl_pn["tissue"].eq("islet_tet_pos").to_numpy(),
        sigs,
        col_fmt="signature_{}_pooled_pct",
        dataset="kasmani",
        compartment="within_islet",
        contrast="islet_tet_pos_vs_islet_tet_neg",
    )

    # Sister-clone analysis: identical CDR3ab present in BOTH Spleen_TetPos and Islet_TetPos in the same round
    piv_sis = (
        comp_clones[comp_clones["tissue"].isin(["islet_tet_pos", "spleen_tet_pos"])]
        .pivot_table(
            index=["donor", "CDR3ab"],
            columns="tissue",
            values=["cells", *[f"pooled_{s}" for s in sigs]],
        )
        .dropna()
    )
    sister_deltas = {}
    for sig in sigs:
        diff = piv_sis[("pooled_" + sig, "islet_tet_pos")] - piv_sis[("pooled_" + sig, "spleen_tet_pos")]
        sister_deltas[sig] = {
            "n_sister_clones": int(len(diff)),
            "mean_delta_islet_minus_spleen": round(float(diff.mean()), 4),
            "frac_higher_in_islet": round(float((diff > 0).mean()), 4),
        }

    # VDJdb MusMusculus evaluation: Circular (full VDJdb) vs Held-out (excluding PMID:35667687)
    igrp_clones = (
        comp_clones[comp_clones["tissue"].isin(["islet_tet_pos", "spleen_tet_pos"])]
        .sort_values("cells", ascending=False)
        .drop_duplicates(["donor", "CDR3ab"])[
            ["donor", "CDR3ab", "CDR3_alpha", "CDR3_beta", "cells"]
        ]
        .copy()
    )
    mc38_ctrl = mc38_clones[["donor", "CDR3ab", "CDR3_alpha", "CDR3_beta", "cell_count"]].rename(
        columns={"cell_count": "cells"}
    )
    db_heldout = vdjdb_full[
        ~vdjdb_full["reference"].fillna("").astype(str).str.contains("35667687", case=False)
    ].copy()

    vdjdb_rows = []
    for db_label, db_df in (("vdjdb_full_circular", vdjdb_full), ("vdjdb_heldout_no_kasmani", db_heldout)):
        for strictness, mode in (
            ("strict_ab", "exact"),
            ("ab_with_partial", "exact"),
            ("b_only", "exact"),
            ("b_only", "levenshtein"),
        ):
            m_igrp = match_clonotypes(
                igrp_clones,
                db_df,
                match_strictness=strictness,
                match_mode=mode,
                species="mouse",
                host_species_filter="same",
                verbose=False,
                show_progress=False,
            )
            m_ctrl = match_clonotypes(
                mc38_ctrl,
                db_df,
                match_strictness=strictness,
                match_mode=mode,
                species="mouse",
                host_species_filter="same",
                verbose=False,
                show_progress=False,
            )
            hit = m_igrp["db_match"].to_numpy(dtype=bool)
            igrp_ep = (
                m_igrp["db_epitope"].fillna("").astype(str).str.contains("KYNKANVFL", case=False)
                | m_igrp["db_protein"].fillna("").astype(str).str.contains("IGRP|G6PC2", case=False)
            ).to_numpy(dtype=bool)
            ctrl_igrp_fp = (
                m_ctrl["db_match"].to_numpy(dtype=bool)
                & (
                    m_ctrl["db_epitope"].fillna("").astype(str).str.contains("KYNKANVFL", case=False)
                    | m_ctrl["db_protein"].fillna("").astype(str).str.contains("IGRP|G6PC2", case=False)
                ).to_numpy(dtype=bool)
            )
            n_hit = int(hit.sum())
            n_igrp = int((hit & igrp_ep).sum())
            vdjdb_rows.append(
                {
                    "cohort": "kasmani_igrp_nod",
                    "database": db_label,
                    "strictness": strictness if mode == "exact" else f"{strictness}_lev1",
                    "n_query_clones": len(igrp_clones),
                    "n_matched": n_hit,
                    "any_match_rate": round(n_hit / max(len(igrp_clones), 1), 4),
                    "n_igrp_matched": n_igrp,
                    "igrp_recall": round(n_igrp / max(len(igrp_clones), 1), 4),
                    "igrp_precision_among_matched": round(n_igrp / n_hit, 4) if n_hit else np.nan,
                    "mc38_control_clones": len(mc38_ctrl),
                    "mc38_any_matched": int(m_ctrl["db_match"].sum()),
                    "mc38_false_igrp_matches": int(ctrl_igrp_fp.sum()),
                }
            )

    facts = {
        "retained_cd8_cells": int(len(adata_cd8)),
        "islet_tet_pos_expanded_clones": int((exp_comp["tissue"] == "islet_tet_pos").sum()),
        "spleen_tet_pos_expanded_clones": int((exp_comp["tissue"] == "spleen_tet_pos").sum()),
        "islet_tet_neg_expanded_clones": int((exp_comp["tissue"] == "islet_tet_neg").sum()),
        "unique_igrp_tet_pos_clones": len(igrp_clones),
        "vdjdb_pmid_35667687_entries": int(
            vdjdb_full["reference"].fillna("").astype(str).str.contains("35667687").sum()
        ),
        "sister_clones_islet_and_spleen": sister_deltas,
    }
    return auroc_rows, vdjdb_rows, facts


# ---------------------------------------------------- 3. Khatun 2021 (GSE158896)


def run_khatun(
    data_dir: Path,
    work_dir: Path,
    cache_dir: Path | None,
    vdjdb_full: pd.DataFrame,
):
    from tcrsift.annotate import match_clonotypes
    from tcrsift.assemble import assemble_full_sequences, validate_sequences
    from tcrsift.cli import main as tcrsift_main

    sroot = data_dir / "khatun" / "samples"
    cache_args = ["--cache-dir", str(cache_dir)] if cache_dir else []

    # Compare with `pre_sorted=CD4` vs without `pre_sorted` across M1-M5
    sheet_sorted = work_dir / "khatun_presorted.csv"
    sheet_unsorted = work_dir / "khatun_unsorted.csv"
    rows_s = ["sample,patient_id,tissue,pre_sorted,vdj_dir,gex_dir"]
    rows_u = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    for i in range(1, 6):
        sname = f"M{i}_GP66_CD4"
        rows_s.append(f"{sname},M{i},spleen,CD4,{sroot / sname / 'vdj'},{sroot / sname / 'gex'}")
        rows_u.append(f"{sname},M{i},spleen,{sroot / sname / 'vdj'},{sroot / sname / 'gex'}")
    sheet_sorted.write_text("\n".join(rows_s) + "\n")
    sheet_unsorted.write_text("\n".join(rows_u) + "\n")

    out_s = work_dir / "khatun_cd4_presorted"
    out_u = work_dir / "khatun_cd4_unsorted"
    for sh, out in ((sheet_sorted, out_s), (sheet_unsorted, out_u)):
        if not is_current(out):
            tcrsift_main(
                [
                    "prioritize",
                    str(sh),
                    "-o",
                    str(out),
                    "--species",
                    "mouse",
                    "--tcell-type",
                    "cd4",
                    "--context",
                    "blood",
                    *cache_args,
                ]
            )

    cfg_s = json.loads((out_s / "prioritization.json").read_text())
    cfg_u = json.loads((out_u / "prioritization.json").read_text())
    clones_s = pd.read_csv(out_s / "all_scored_clones.csv", low_memory=False)
    clones_u = pd.read_csv(out_u / "all_scored_clones.csv", low_memory=False)
    cands_s = pd.read_csv(out_s / "candidate_clones.csv", low_memory=False)

    # Cross-mouse public receptor sharing across genetically identical C57BL/6 mice (M1-M5)
    beta_mice = clones_s.groupby("CDR3_beta")["donor"].nunique()
    ab_mice = clones_s.groupby("CDR3ab")["donor"].nunique()

    # VDJdb LCMV matching across M1-M5 unique CDR3ab clones
    uniq_gp66 = clones_s.sort_values("cell_count", ascending=False).drop_duplicates("CDR3ab").copy()
    db_no_khatun = vdjdb_full[
        ~vdjdb_full["reference"].fillna("").astype(str).str.contains("33201171", case=False)
    ].copy()
    vdjdb_rows = []
    for strictness, mode in (
        ("strict_ab", "exact"),
        ("ab_with_partial", "exact"),
        ("b_only", "exact"),
        ("b_only", "levenshtein"),
    ):
        m = match_clonotypes(
            uniq_gp66,
            db_no_khatun,
            match_strictness=strictness,
            match_mode=mode,
            species="mouse",
            host_species_filter="same",
            verbose=False,
            show_progress=False,
        )
        hit = m["db_match"].to_numpy(dtype=bool)
        lcmv_hit = hit & m["db_species"].fillna("").astype(str).str.contains("LCMV", case=False).to_numpy(dtype=bool)
        n_hit = int(hit.sum())
        n_lcmv = int(lcmv_hit.sum())
        vdjdb_rows.append(
            {
                "cohort": "khatun_gp66_cd4",
                "database": "vdjdb_independent",
                "strictness": strictness if mode == "exact" else f"{strictness}_lev1",
                "n_query_clones": len(uniq_gp66),
                "n_matched": n_hit,
                "any_match_rate": round(n_hit / max(len(uniq_gp66), 1), 4),
                "n_igrp_matched": n_lcmv,  # target pathogen = LCMV
                "igrp_recall": round(n_lcmv / max(len(uniq_gp66), 1), 4),
                "igrp_precision_among_matched": round(n_lcmv / n_hit, 4) if n_hit else np.nan,
                "mc38_control_clones": 0,
                "mc38_any_matched": 0,
                "mc38_false_igrp_matches": 0,
            }
        )

    # Full-length murine CD4 construct assembly across all shortlisted candidates
    asm_cd4 = assemble_full_sequences(
        cands_s,
        species="mouse",
        germline_vdj_fallback=True,
        data_dir=cache_dir,
        verbose=False,
        show_progress=False,
    )
    val_cd4 = [
        m for m in validate_sequences(asm_cd4, species="mouse") if m.severity == "load_bearing"
    ]

    sig_breakdown = {
        str(k): int(v) for k, v in cands_s["selected_signature"].value_counts().items()
    }
    facts = {
        "vdjdb_pmid_33201171_entries": int(
            vdjdb_full["reference"].fillna("").astype(str).str.contains("33201171").sum()
        ),
        "unsorted_retained_cells": int(cfg_u["cell_filters"][-1]["after"]),
        "presorted_retained_cells": int(cfg_s["cell_filters"][-1]["after"]),
        "unsorted_clones": len(clones_u),
        "presorted_clones": len(clones_s),
        "presorted_expanded_clones": int((clones_s["cell_count"] >= 2).sum()),
        "shortlisted_candidates": len(cands_s),
        "shortlisted_by_signature": sig_breakdown,
        "foxp3_treg_like_clones": int(clones_s["is_treg_like"].sum()),
        "foxp3_treg_like_shortlisted": int(cands_s["is_treg_like"].sum()),
        "unique_cdr3_beta_across_5_mice": int(len(beta_mice)),
        "public_cdr3_beta_ge2_mice": int((beta_mice >= 2).sum()),
        "public_cdr3_ab_ge2_mice": int((ab_mice >= 2).sum()),
        "assembled_cd4_candidates": len(asm_cd4),
        "assembled_cd4_complete": int(asm_cd4["single_chain_aa"].notna().sum()),
        "assembled_cd4_validation_errors": len(val_cd4),
    }
    return vdjdb_rows, asm_cd4, facts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional tcrsift cache directory containing vdjdb and mouse IMGT V/J references",
    )
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    from tcrsift import __version__
    from tcrsift.annotate import load_vdjdb
    from tcrsift.datacache import ensure_cached

    vdjdb_path = ensure_cached("vdjdb", data_dir=args.cache_dir, auto_download=False)
    if vdjdb_path is None:
        raise RuntimeError(
            "VDJdb not found in cache; run `tcrsift data download --db vdjdb --group germline_mouse`"
        )
    vdjdb_subs = sorted(p.name for p in Path(vdjdb_path).glob("vdjdb-*") if p.is_dir())
    vdjdb_release = vdjdb_subs[-1] if vdjdb_subs else Path(vdjdb_path).name
    vdjdb_full = load_vdjdb(vdjdb_path, verbose=False)
    mus_rows = vdjdb_full[
        vdjdb_full["host_species"].fillna("").astype(str).str.contains("Mus|mouse", case=False)
    ]

    t0 = time.perf_counter()
    print("1/3 Running MC38 multi-compartment murine benchmark (GSE168944)...")
    mc38_aurocs, mc38_shortlists, asm_mc38, mc38_tum_clones, mc38_facts = run_mc38(
        args.data_dir, args.work_dir, args.cache_dir
    )

    print("2/3 Running Kasmani 2022 NOD islet/spleen IGRP CD8 benchmark (GSE200608)...")
    kas_aurocs, kas_vdjdb, kas_facts = run_kasmani(
        args.data_dir, args.work_dir, args.cache_dir, vdjdb_full, mc38_tum_clones
    )

    print("3/3 Running Khatun 2021 LCMV GP66 murine CD4 benchmark (GSE158896)...")
    kha_vdjdb, asm_cd4, kha_facts = run_khatun(
        args.data_dir, args.work_dir, args.cache_dir, vdjdb_full
    )

    auroc_df = pd.DataFrame(mc38_aurocs + kas_aurocs)
    auroc_df.to_csv(results_dir / "auroc_mouse.csv", index=False)

    vdjdb_df = pd.DataFrame(kas_vdjdb + kha_vdjdb)
    vdjdb_df.to_csv(results_dir / "vdjdb_mouse_recovery.csv", index=False)

    shortlist_df = pd.DataFrame(mc38_shortlists)
    shortlist_df.to_csv(results_dir / "shortlist_mc38.csv", index=False)

    asm_summary_rows = []
    for label, adf in (("mc38_cd8", asm_mc38), ("khatun_gp66_cd4", asm_cd4)):
        c_col = "beta_c_gene_canonical" if "beta_c_gene_canonical" in adf.columns else "beta_c_gene"
        trbc1 = int(adf[c_col].fillna("").astype(str).str.upper().str.startswith("TRBC1").sum())
        trbc2 = int(adf[c_col].fillna("").astype(str).str.upper().str.startswith("TRBC2").sum())
        overridden = int(adf["beta_c_gene_overridden"].sum()) if "beta_c_gene_overridden" in adf.columns else 0
        asm_summary_rows.append(
            {
                "cohort": label,
                "n_candidates": len(adf),
                "n_complete_single_chain": int(adf["single_chain_aa"].notna().sum()),
                "complete_rate": round(float(adf["single_chain_aa"].notna().mean()), 4),
                "alpha_germline_reconstructed": int(
                    (adf["alpha_vdj_source"] == "germline_vj_reconstructed").sum()
                ),
                "beta_germline_reconstructed": int(
                    (adf["beta_vdj_source"] == "germline_vj_reconstructed").sum()
                ),
                "trbc1_constructs": trbc1,
                "trbc2_constructs": trbc2,
                "cellranger_trbc_miscalls_corrected": overridden,
                "median_single_chain_aa_len": int(adf["single_chain_aa"].dropna().str.len().median()),
            }
        )
    pd.DataFrame(asm_summary_rows).to_csv(results_dir / "mouse_assembly.csv", index=False)

    summary = {
        "tcrsift_version": __version__,
        "vdjdb_release": vdjdb_release,
        "vdjdb_mus_musculus_trb_entries": int(len(mus_rows)),
        "vdjdb_mus_musculus_paired_ab_entries": int(
            (mus_rows["cdr3_alpha"].notna() & (mus_rows["cdr3_alpha"] != "")).sum()
        ),
        "elapsed_seconds": round(time.perf_counter() - t0, 1),
        "ground_truth_provenance": {
            "mc38": {
                "study": "Multi-tissue murine MC38 colon adenocarcinoma (GEO GSE168944)",
                "accession": "GSE168944 (GSM5172668-GSM5172675 GEX; GSM5172684-GSM5172691 VDJ)",
                "ground_truth": "Matched 4-organ paired alpha/beta TCR tracking across Tumor, Blood, dLN, and Spleen in biological replicate C57BL/6 mice (Mouse01 & Mouse02)",
            },
            "kasmani": {
                "study": "Kasmani et al., Life Science Alliance 2022 (PMID:35667687)",
                "accession": "GSE200608 (GSM6038943-GSM6038954)",
                "ground_truth": "FACS-sorted H-2Kd IGRP206-214 / NRP-V7 (KYNKANVFL) MHC-I tetramer+ CD8+ T cells in NOD mouse pancreatic islets vs spleen, and islet tetramer- CD8+ T cells",
                "circularity_note": "PMID:35667687 contributes 780 MusMusculus TRB rows in VDJdb 2026-11-06; both circular full-VDJdb and study-held-out VDJdb evaluations are reported.",
            },
            "khatun": {
                "study": "Khatun et al., J Exp Med 2021 (PMID:33201171)",
                "accession": "GSE158896 (GSM4815363-GSM4815372)",
                "ground_truth": "FACS-sorted I-Ab LCMV GP66-77 (DIYKGVYQFKSV) MHC-II tetramer+ splenic CD4+ T cells at Day 10 post-acute LCMV Armstrong infection across 5 C57BL/6 mice (M1-M5)",
            },
        },
        "mc38": mc38_facts,
        "kasmani": kas_facts,
        "khatun": kha_facts,
    }
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Done in {summary['elapsed_seconds']}s -> {results_dir}")


if __name__ == "__main__":
    main()
