#!/usr/bin/env python3
"""Run the public infectious-disease benchmark and write summary tables to ``results/``.

    python benchmarks/public_infectious/download.py --data-dir CACHE
    python benchmarks/public_infectious/run.py --data-dir CACHE --work-dir CACHE/runs

Exercises ``tcrsift`` CLI and library workflows across three public human
infectious-disease cohorts:

1. 10x Genomics Multiplexed 50-Dextramer Donors 1 & 2 (Influenza M1, EBV, CMV,
   MART-1 + 6 negative-control dextramers): ``tcrsift prioritize --context blood``,
   ``tcrsift.multimer`` negative-control-calibrated pMHC calling, and VDJdb
   specificity recovery (both full VDJdb and holding out the 10x study).
2. Minervina et al. 2022 (Zenodo 6232103 — SARS-CoV-2 primary infection,
   vaccination, booster recall, and breakthrough CD8+ T cells): ``tcrsift
   prioritize --context blood`` longitudinal state trajectories and VDJdb
   specificity recovery (full VDJdb and holding out PMID:35383307).
3. Fischer et al. 2021 (GSE171037 — paired unstimulated vs 16 h SARS-CoV-2 Spike
   peptide-stimulated PBMC from hospitalized severe COVID-19 patients):
   ``load_samples``, ``phenotype_cells`` + ``apply_cd3_gate``, and
   ``score_stimulation_response`` (abundance fold-change vs paired
   activation-signature and cytokine transcript delta across persistent clones).
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
VIRAL_BONA_FIDE_RE = r"Influenza|EBV|pp65_CMV"
FLU_EPITOPE = "A0201_GILGFVFTL_Flu-MP_Influenza"
MART1_EPITOPE = "A0201_ELAGIGILTV_MART-1_Cancer"
STICKY_KLG_EPITOPE = "A0301_KLGGALQAK_IE-1_CMV"


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


# ------------------------------------------------ 1. 10x Dextramer Donors 1 & 2

def run_tenx_dextramer(data_dir: Path, work_dir: Path, vdjdb_full: pd.DataFrame):
    from tcrsift._schema import truthy_flag_series
    from tcrsift.annotate import match_clonotypes
    from tcrsift.cli import main as tcrsift_main
    from tcrsift.loader import _pivot_vdj_by_barcode, load_cellranger_vdj
    from tcrsift.multimer import (
        call_multimer_specificity,
        classify_feature_columns,
        extract_feature_barcodes,
        summarize_clone_multimer_specificity,
    )

    auroc_rows: list[dict] = []
    vdjdb_rows: list[dict] = []
    multimer_rows: list[dict] = []
    shortlist_rows: list[dict] = []
    facts: dict[str, dict] = {}

    db_no_10x = vdjdb_full[
        ~vdjdb_full["reference"].fillna("").astype(str).str.contains("10xgenomics|tenx", case=False)
    ].copy()
    heldout_10x_tsv = work_dir / "vdjdb_heldout_no_10x.tsv"
    if not heldout_10x_tsv.exists():
        # Emit a paired-format TSV so unmodified CLI `prioritize --vdjdb ...` can consume study-held-out VDJdb
        export_df = pd.DataFrame(
            {
                "cdr3.alpha": db_no_10x["cdr3_alpha"],
                "cdr3.beta": db_no_10x["cdr3_beta"],
                "v.alpha": db_no_10x.get("v_alpha", pd.NA),
                "v.beta": db_no_10x.get("v_beta", pd.NA),
                "species": db_no_10x.get("host_species", "HomoSapiens"),
                "antigen.epitope": db_no_10x["epitope"],
                "antigen.gene": db_no_10x["antigen_gene"],
                "antigen.species": db_no_10x["species"],
                "mhc.a": db_no_10x.get("mhc_allele", pd.NA),
                "mhc.class": db_no_10x.get("mhc_class", pd.NA),
                "reference.id": db_no_10x.get("reference", pd.NA),
            }
        )
        export_df.to_csv(heldout_10x_tsv, sep="\t", index=False)

    for donor in ("donor1", "donor2"):
        ddir = data_dir / "tenx_dextramer" / donor
        sheet = work_dir / f"tenx_{donor}_sheet.csv"
        sheet.write_text(
            f"sample,patient_id,tissue,vdj_dir,gex_dir\n{donor},{donor},blood,{ddir},{ddir}\n"
        )
        # Score the 7 blood-preset signatures plus NeoTCR_PBL and TumorReactive
        # across all paired-TCR CD8 clones (--min-cells 1 so singleton MART-1
        # naive precursors are retained in all_scored_clones.csv; expanded >=2
        # cell clones are evaluated separately below).
        out = work_dir / f"tenx_{donor}_blood"
        sec_file = out / "benchmark_seconds.txt"
        if not is_current(out):
            t0 = time.perf_counter()
            tcrsift_main(
                [
                    "prioritize",
                    str(sheet),
                    "-o",
                    str(out),
                    "--context",
                    "blood",
                    "--signatures",
                    "Differentiated",
                    "AntigenExperienced",
                    "Cytolytic",
                    "CirculatingMemory",
                    "AcuteActivation",
                    "Proliferation",
                    "AIM",
                    "NeoTCR_PBL",
                    "TumorReactive",
                    "--min-cells",
                    "1",
                ]
            )
            elapsed = round(time.perf_counter() - t0, 1)
            sec_file.write_text(f"{elapsed}\n")
        else:
            elapsed = float(sec_file.read_text().strip()) if sec_file.exists() else None
        cfg = json.loads((out / "prioritization.json").read_text())
        sigs = cfg["resolved_signatures"]
        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones["cell_count"].rank(pct=True)

        # Load 10x ground-truth binarized_matrix.csv + raw feature-barcode counts
        bin_f = next(ddir.glob("*binarized_matrix.csv"))
        bdf = pd.read_csv(bin_f).set_index("barcode")
        binder_cols = [c for c in bdf.columns if c.endswith("_binder")]
        bmat = bdf[binder_cols].infer_objects(copy=False).fillna(False).astype(bool)
        nb = bmat.sum(axis=1)
        single_bcs = bmat.index[nb == 1]
        tenx_cell_call = bmat.loc[single_bcs].idxmax(axis=1).str.replace("_binder$", "", regex=True)

        # Map barcodes -> paired CDR3ab from VDJ contigs with >=50% clonal concordance
        vdj_df = load_cellranger_vdj(ddir, donor, verbose=False)
        piv = _pivot_vdj_by_barcode(vdj_df)
        piv_paired = piv[
            piv["CDR3_alpha"].notna()
            & (piv["CDR3_alpha"] != "")
            & piv["CDR3_beta"].notna()
            & (piv["CDR3_beta"] != "")
        ].copy()
        piv_paired["tenx_call"] = tenx_cell_call.reindex(piv_paired.index)
        total_by_clone = piv_paired.groupby("CDR3ab").size().rename("vdj_clone_cells")
        cc = (
            piv_paired.dropna(subset=["tenx_call"])
            .groupby(["CDR3ab", "tenx_call"])
            .size()
            .reset_index(name="tenx_cells")
            .sort_values(["tenx_cells", "tenx_call"], ascending=[False, True], kind="mergesort")
            .drop_duplicates("CDR3ab")
            .merge(total_by_clone, on="CDR3ab", how="left")
        )
        cc["tenx_concordance"] = cc["tenx_cells"] / cc["vdj_clone_cells"].clip(lower=1)
        cc = cc[cc["tenx_concordance"] >= 0.5].copy()
        clones = clones.merge(
            cc[["CDR3ab", "tenx_call", "tenx_cells", "tenx_concordance"]],
            on="CDR3ab",
            how="left",
        )

        # Extract feature barcodes + compute centered log-ratio (CLR) CITE-seq CD45RO/CD45RA
        fb_df = extract_feature_barcodes(ddir)
        fclasses = classify_feature_columns(fb_df)
        ro_log = np.log1p(fb_df["CD45RO_TotalSeqC"])
        ra_log = np.log1p(fb_df["CD45RA_TotalSeqC"])
        ro_clr = ro_log - float(ro_log.mean())
        ra_clr = ra_log - float(ra_log.mean())
        piv_paired["cd45ro_clr"] = ro_clr.reindex(piv_paired.index)
        piv_paired["cd45ra_clr"] = ra_clr.reindex(piv_paired.index)
        piv_paired["cd45ro_log1p"] = ro_log.reindex(piv_paired.index)
        piv_paired["cd45ra_log1p"] = ra_log.reindex(piv_paired.index)
        prot_by_clone = piv_paired.groupby("CDR3ab")[
            ["cd45ro_clr", "cd45ra_clr", "cd45ro_log1p", "cd45ra_log1p"]
        ].mean()
        clones = clones.merge(prot_by_clone, left_on="CDR3ab", right_index=True, how="left")

        # Comparison A: Viral memory (Flu/EBV/CMV pp65) vs MART-1 naive precursors (all clones >=1 cell)
        is_viral = clones["tenx_call"].fillna("").str.contains(VIRAL_BONA_FIDE_RE)
        is_mart1 = clones["tenx_call"].eq(MART1_EPITOPE)
        sub_vm = clones[is_viral | is_mart1].copy()
        auroc_rows += eval_auroc_block(
            sub_vm,
            is_viral.loc[sub_vm.index].to_numpy(),
            sigs,
            dataset=f"tenx_{donor}",
            comparison="viral_memory_vs_naive_mart1",
            min_cells=1,
        )

        mart1_clones = clones[is_mart1]
        viral_clones = clones[is_viral]
        mart1_singleton_frac = (
            float((mart1_clones["cell_count"] == 1).mean()) if len(mart1_clones) else np.nan
        )
        mart1_memory_clr = float(
            (mart1_clones["cd45ro_clr"] - mart1_clones["cd45ra_clr"]).mean()
        )
        viral_memory_clr = float(
            (viral_clones["cd45ro_clr"] - viral_clones["cd45ra_clr"]).mean()
        )

        # Comparison B: Expanded (>=2 cells) viral clones vs unlabeled expanded blood clones
        exp = clones[clones["cell_count"] >= 2].copy()
        exp["size_percentile"] = exp["cell_count"].rank(pct=True)
        exp_viral = exp["tenx_call"].fillna("").str.contains(VIRAL_BONA_FIDE_RE) & (
            exp["tenx_cells"].fillna(0) >= 2
        )
        exp_unlab = exp["tenx_call"].isna()
        sub_vu = exp[exp_viral | exp_unlab].copy()
        auroc_rows += eval_auroc_block(
            sub_vu,
            exp_viral.loc[sub_vu.index].to_numpy(),
            sigs,
            dataset=f"tenx_{donor}",
            comparison="expanded_viral_vs_unlabeled_expanded",
            min_cells=2,
        )

        # Comparison C: Cleared respiratory virus (Flu M1 GILGFVFTL) vs latent herpesvirus (EBV) among expanded clones
        exp_flu = exp["tenx_call"].eq(FLU_EPITOPE) & (exp["tenx_cells"].fillna(0) >= 2)
        exp_ebv = exp["tenx_call"].fillna("").str.contains("_EBV$") & (
            exp["tenx_cells"].fillna(0) >= 2
        )
        sub_fe = exp[exp_flu | exp_ebv].copy()
        auroc_rows += eval_auroc_block(
            sub_fe,
            exp_flu.loc[sub_fe.index].to_numpy(),
            sigs,
            dataset=f"tenx_{donor}",
            comparison="influenza_m1_vs_ebv_expanded",
            min_cells=2,
        )

        # Question 3: VDJdb specificity recovery on expanded viral + MART-1 clones
        is_eval_spec = exp["tenx_call"].fillna("").str.contains(VIRAL_BONA_FIDE_RE) | exp[
            "tenx_call"
        ].eq(MART1_EPITOPE)
        eval_clones = exp[(exp["tenx_cells"].fillna(0) >= 2) & is_eval_spec].copy()
        eval_clones["true_peptide"] = eval_clones["tenx_call"].str.split("_").str[1]
        for db_scope, db_frame, circular in (
            ("full_vdjdb", vdjdb_full, True),
            ("leave_study_out", db_no_10x, False),
        ):
            for mode, strictness in (
                ("exact", "strict_ab"),
                ("exact", "ab_with_partial"),
                ("exact", "b_only"),
                ("levenshtein", "b_only"),
            ):
                ann = match_clonotypes(
                    eval_clones,
                    db_frame,
                    match_mode=mode,
                    match_strictness=strictness,
                    species="human",
                    host_species_filter="same",
                    verbose=False,
                    show_progress=False,
                )
                has_hit = ann["db_epitope"].notna()
                n_hits = int(has_hit.sum())
                exact_pep = sum(
                    tp in str(de).split(";") if pd.notna(de) else False
                    for tp, de in zip(ann["true_peptide"], ann["db_epitope"])
                )
                vdjdb_rows.append(
                    {
                        "dataset": f"tenx_{donor}",
                        "db_scope": db_scope,
                        "circular": circular,
                        "match_mode": mode,
                        "match_strictness": strictness,
                        "n_query_clones": len(eval_clones),
                        "n_matched": n_hits,
                        "n_exact_epitope": int(exact_pep),
                        "epitope_precision": round(exact_pep / n_hits, 4) if n_hits else np.nan,
                        "epitope_recall": round(exact_pep / len(eval_clones), 4)
                        if len(eval_clones)
                        else np.nan,
                    }
                )

        # Question 4: Feature-barcode multimer calling + clonal concordance + CITE-seq CLR memory QC
        mcalls = call_multimer_specificity(
            fb_df,
            multimer_cols=fclasses["multimers"],
            nc_cols=fclasses["negative_controls"],
            min_umi=10.0,
            min_nc_ratio=5.0,
            min_dominance_ratio=2.0,
        )
        clone_mm = summarize_clone_multimer_specificity(
            piv_paired, mcalls, min_specific_cells=2, min_concordance=0.5
        )
        clone_mm = clone_mm.merge(prot_by_clone, left_on="CDR3ab", right_index=True, how="left")

        # Compare cell-level agreement with 10x binarized_matrix on single-binder cells
        shared_bcs = fb_df.index.intersection(bdf.index)
        tenx_pep = tenx_cell_call.reindex(shared_bcs).str.split("_").str[1]
        our_pep = mcalls.reindex(shared_bcs)["multimer_call"].dropna().astype(str).str.split("_").str[0]
        both_called = tenx_pep.notna() & mcalls.reindex(shared_bcs)["multimer_call"].notna()
        cell_agreement = float((tenx_pep[both_called] == our_pep.reindex(shared_bcs)[both_called]).mean())

        # Summarize per major multimer specificity across expanded clones (>=2 cells)
        exp_mm = clone_mm[(clone_mm["clone_cells"] >= 2) & clone_mm["multimer_call"].notna()].copy()
        for mname, grp in exp_mm.groupby("multimer_call"):
            if len(grp) < 5:
                continue
            multimer_rows.append(
                {
                    "dataset": f"tenx_{donor}",
                    "multimer": str(mname),
                    "is_known_sticky": "KLGGALQAK" in str(mname),
                    "expanded_clones_with_any_call": len(grp),
                    "confident_clones_concordance_ge_50pct": int(grp["multimer_confident"].sum()),
                    "median_clone_concordance": round(float(grp["multimer_concordance"].median()), 4),
                    "mean_clone_cells": round(float(grp["clone_cells"].mean()), 2),
                    "mean_cd45ro_clr": round(float(grp["cd45ro_clr"].mean()), 3),
                    "mean_cd45ra_clr": round(float(grp["cd45ra_clr"].mean()), 3),
                    "mean_memory_clr_diff": round(
                        float((grp["cd45ro_clr"] - grp["cd45ra_clr"]).mean()), 3
                    ),
                    "mean_cd45ro_log1p": round(float(grp["cd45ro_log1p"].mean()), 3),
                    "mean_cd45ra_log1p": round(float(grp["cd45ra_log1p"].mean()), 3),
                    "cell_peptide_agreement_vs_10x": round(cell_agreement, 4),
                }
            )

        # Question 2b: Unmodified default `tcrsift prioritize --context blood` (--min-cells 2, 7-sig blood preset)
        # evaluated with held-out VDJdb under default, --exclude-known-viral (strict_ab default), and --exclude-known-viral (ab_with_partial).
        empty_cache = work_dir / "_empty_db_cache"
        empty_cache.mkdir(parents=True, exist_ok=True)
        label_lookup = cc.set_index("CDR3ab")[["tenx_call", "tenx_cells"]]
        for run_tag, extra_cli in (
            ("blood_default", ["--vdjdb", str(heldout_10x_tsv), "--cache-dir", str(empty_cache)]),
            (
                "blood_exclude_viral_strict_ab",
                ["--vdjdb", str(heldout_10x_tsv), "--cache-dir", str(empty_cache), "--exclude-known-viral"],
            ),
            (
                "blood_exclude_viral_ab_with_partial",
                [
                    "--vdjdb",
                    str(heldout_10x_tsv),
                    "--cache-dir",
                    str(empty_cache),
                    "--exclude-known-viral",
                    "--database-match",
                    "ab_with_partial",
                ],
            ),
        ):
            out_def = work_dir / f"tenx_{donor}_{run_tag}"
            if not is_current(out_def):
                tcrsift_main(["prioritize", str(sheet), "-o", str(out_def), "--context", "blood", *extra_cli])
            def_clones = pd.read_csv(out_def / "all_scored_clones.csv", low_memory=False)
            def_clones = def_clones.join(label_lookup, on="CDR3ab", how="left")
            def_exp = def_clones[def_clones["cell_count"] >= 2].copy()
            d_viral = def_exp["tenx_call"].fillna("").str.contains(VIRAL_BONA_FIDE_RE) & (
                def_exp["tenx_cells"].fillna(0) >= 2
            )
            d_flu = def_exp["tenx_call"].eq(FLU_EPITOPE) & (def_exp["tenx_cells"].fillna(0) >= 2)
            d_ebv = def_exp["tenx_call"].fillna("").str.contains("_EBV$") & (
                def_exp["tenx_cells"].fillna(0) >= 2
            )
            d_mart1 = def_exp["tenx_call"].eq(MART1_EPITOPE)
            d_unlab = def_exp["tenx_call"].isna()
            elig = truthy_flag_series(def_exp["eligible_for_review"])
            sel = truthy_flag_series(def_exp["selected_for_review"])
            shortlist_rows.append(
                {
                    "dataset": f"tenx_{donor}",
                    "run": run_tag,
                    "expanded_clones_ge2": len(def_exp),
                    "eligible_total": int(elig.sum()),
                    "selected_total": int(sel.sum()),
                    "viral_total": int(d_viral.sum()),
                    "viral_eligible": int((d_viral & elig).sum()),
                    "viral_selected": int((d_viral & sel).sum()),
                    "flu_total": int(d_flu.sum()),
                    "flu_eligible": int((d_flu & elig).sum()),
                    "flu_selected": int((d_flu & sel).sum()),
                    "ebv_total": int(d_ebv.sum()),
                    "ebv_eligible": int((d_ebv & elig).sum()),
                    "ebv_selected": int((d_ebv & sel).sum()),
                    "mart1_expanded_total": int(d_mart1.sum()),
                    "mart1_expanded_eligible": int((d_mart1 & elig).sum()),
                    "mart1_expanded_selected": int((d_mart1 & sel).sum()),
                    "unlabeled_expanded_total": int(d_unlab.sum()),
                    "unlabeled_expanded_eligible": int((d_unlab & elig).sum()),
                    "unlabeled_expanded_selected": int((d_unlab & sel).sum()),
                }
            )

        facts[f"tenx_{donor}"] = {
            "seconds": elapsed,
            "cells_after_gex_qc": cfg["cell_filters"][0]["before"],
            "cells_retained": cfg["cell_filters"][-1]["after"],
            "clones_scored": len(clones),
            "expanded_clones_ge2": len(exp),
            "mart1_clones_total": len(mart1_clones),
            "mart1_singleton_clones": int((mart1_clones["cell_count"] == 1).sum()),
            "mart1_singleton_fraction": round(mart1_singleton_frac, 4),
            "mart1_mean_memory_clr_diff": round(mart1_memory_clr, 3),
            "viral_mean_memory_clr_diff": round(viral_memory_clr, 3),
            "multimer_features": len(fclasses["multimers"]),
            "negative_control_features": len(fclasses["negative_controls"]),
            "cell_peptide_agreement_vs_10x": round(cell_agreement, 4),
        }

    return auroc_rows, vdjdb_rows, multimer_rows, shortlist_rows, facts


# ------------------------------------------------ 2. Minervina et al. 2022

def run_minervina(data_dir: Path, work_dir: Path, vdjdb_full: pd.DataFrame):
    from tcrsift._schema import truthy_flag_series
    from tcrsift.annotate import match_clonotypes
    from tcrsift.cli import main as tcrsift_main
    from tcrsift.loader import combine_gex_and_vdj, load_cellranger_gex, load_cellranger_vdj
    from tcrsift.phenotype import apply_cd3_gate, phenotype_cells
    from tcrsift.signature_methods import score_by_name

    mdir = data_dir / "minervina_2022_covid_flu"
    agg_dir = mdir / "extracted" / "libs_aggregate_rev"
    tsv = pd.read_csv(mdir / "cd8_only_dextr_rev_clean.tsv", sep="\t", low_memory=False)

    sheet = work_dir / "minervina_sheet.csv"
    sheet.write_text(
        f"sample,patient_id,tissue,vdj_dir,gex_dir\nminervina,cohort,blood,{agg_dir},{agg_dir}\n"
    )
    out = work_dir / "minervina_blood"
    if not is_current(out):
        tcrsift_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out),
                "--context",
                "blood",
                "--min-cells",
                "1",
            ]
        )
    cfg = json.loads((out / "prioritization.json").read_text())
    sigs = cfg["resolved_signatures"]
    clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)

    adata_m = load_cellranger_gex(agg_dir, "minervina", min_mito_pct=0.0, verbose=False)
    vdj_m = load_cellranger_vdj(agg_dir, "minervina", verbose=False)
    comb_m = apply_cd3_gate(phenotype_cells(combine_gex_and_vdj(adata_m, vdj_m), min_cd3_reads=3), 3)
    bc_meta = tsv.set_index("barcode")[["donor", "epitope", "spike", "donor_category", "sample_category"]]
    comb_m.obs = comb_m.obs.join(bc_meta, how="left")

    obs_paired = comb_m.obs[
        comb_m.obs["CDR3ab"].notna()
        & ~comb_m.obs["CDR3ab"].astype(str).str.startswith("_")
        & ~comb_m.obs["CDR3ab"].astype(str).str.endswith("_")
        & (comb_m.obs["CDR3ab"].astype(str) != "")
    ].copy()
    clone_ep = (
        obs_paired.dropna(subset=["epitope"])
        .groupby(["CDR3ab", "epitope"])
        .size()
        .reset_index(name="ep_cells")
        .sort_values(["ep_cells", "epitope"], ascending=[False, True], kind="mergesort")
        .drop_duplicates("CDR3ab")
    )
    clones = clones.merge(clone_ep, on="CDR3ab", how="left")

    # Longitudinal exposure state profiles across paired-TCR SARS-CoV-2 CD8+ T cells
    paired_mask = (
        comb_m.obs["donor_category"].notna()
        & comb_m.obs["is_CD8"].astype(bool)
        & comb_m.obs["CDR3ab"].notna()
        & ~comb_m.obs["CDR3ab"].astype(str).str.startswith("_")
        & ~comb_m.obs["CDR3ab"].astype(str).str.endswith("_")
    )
    sub_cells = comb_m[paired_mask].copy()
    sc.pp.normalize_total(sub_cells, target_sum=10_000)
    sc.pp.log1p(sub_cells)
    eval_sigs = (*sigs, "TumorReactive")
    for sname in eval_sigs:
        sub_cells.obs[f"sig_{sname}"] = score_by_name(sub_cells, sname, log1p=False, species="human").to_numpy()

    exposure_rows = []
    for cat, grp in sub_cells.obs.groupby("donor_category", observed=True):
        row = {
            "dataset": "minervina_2022",
            "donor_category": str(cat),
            "n_cells": len(grp),
            "n_clones": int(grp["CDR3ab"].dropna().nunique()),
        }
        for sname in eval_sigs:
            row[f"mean_{sname}"] = round(float(grp[f"sig_{sname}"].mean()), 4)
        exposure_rows.append(row)

    # AUROC: acute booster recall (inf_after_2) vs resting post-infection convalescent (post_inf)
    # Exclude longitudinal clones sampled in both categories so the two clone sets are disjoint.
    clone_cat_scores = (
        sub_cells.obs.groupby(["donor_category", "CDR3ab"], observed=True)[[f"sig_{s}" for s in eval_sigs]]
        .mean()
        .reset_index()
    )
    clone_cat_counts = (
        sub_cells.obs.groupby(["donor_category", "CDR3ab"], observed=True)
        .size()
        .reset_index(name="cat_cells")
    )
    clone_cat = clone_cat_scores.merge(clone_cat_counts, on=["donor_category", "CDR3ab"])
    clone_cat = clone_cat[
        clone_cat["donor_category"].isin(["inf_after_2", "post_inf"]) & (clone_cat["cat_cells"] >= 2)
    ].copy()
    shared_longitudinal = set(
        clone_cat.loc[clone_cat["donor_category"].eq("inf_after_2"), "CDR3ab"]
    ) & set(clone_cat.loc[clone_cat["donor_category"].eq("post_inf"), "CDR3ab"])
    clone_cat = clone_cat[~clone_cat["CDR3ab"].isin(shared_longitudinal)].copy()
    clone_cat["size_percentile"] = clone_cat["cat_cells"].rank(pct=True)
    auroc_rows = eval_auroc_block(
        clone_cat,
        clone_cat["donor_category"].eq("inf_after_2").to_numpy(),
        list(eval_sigs),
        col_fmt="sig_{}",
        size_col="size_percentile",
        dataset="minervina_2022",
        comparison="acute_post_vax_dose2_vs_resting_convalescent",
        min_cells=2,
    )

    # VDJdb matching: full VDJdb vs holding out Minervina 2022 (PMID:35383307)
    sars_clones = clones[clones["epitope"].notna()].copy()
    db_no_minervina = vdjdb_full[
        ~vdjdb_full["reference"].fillna("").astype(str).str.contains("35383307", case=False)
    ].copy()

    sars_clones["pep_prefix"] = sars_clones["epitope"].astype(str).str.split("_").str[1]
    vdjdb_rows = []
    off_target_species: dict[str, int] = {}
    for db_scope, db_frame, circular in (
        ("full_vdjdb", vdjdb_full, True),
        ("leave_study_out", db_no_minervina, False),
    ):
        for mode, strictness in (
            ("exact", "strict_ab"),
            ("exact", "ab_with_partial"),
            ("exact", "b_only"),
            ("levenshtein", "b_only"),
        ):
            ann = match_clonotypes(
                sars_clones,
                db_frame,
                match_mode=mode,
                match_strictness=strictness,
                species="human",
                host_species_filter="same",
                verbose=False,
                show_progress=False,
            )
            has_hit = ann["db_species"].notna()
            n_hits = int(has_hit.sum())
            sars_sp = int(ann["db_species"].fillna("").str.contains("SARS", case=False).sum())
            pep_ok = sum(
                any(ep.startswith(pref) for ep in str(de).split(";"))
                if (pd.notna(de) and isinstance(pref, str) and pref)
                else False
                for pref, de in zip(ann["pep_prefix"], ann["db_epitope"])
            )
            if db_scope == "leave_study_out" and mode == "exact" and strictness == "ab_with_partial":
                wrong_sp = ann.loc[
                    has_hit & ~ann["db_species"].fillna("").str.contains("SARS", case=False),
                    "db_species",
                ]
                off_target_species = {
                    str(k): int(v) for k, v in wrong_sp.value_counts().head(5).items()
                }
            vdjdb_rows.append(
                {
                    "dataset": "minervina_2022_sarscov2",
                    "db_scope": db_scope,
                    "circular": circular,
                    "match_mode": mode,
                    "match_strictness": strictness,
                    "n_query_clones": len(sars_clones),
                    "n_matched": n_hits,
                    "n_exact_epitope": int(pep_ok),
                    "n_correct_virus": sars_sp,
                    "epitope_precision": round(pep_ok / n_hits, 4) if n_hits else np.nan,
                    "virus_precision": round(sars_sp / n_hits, 4) if n_hits else np.nan,
                    "epitope_recall": round(pep_ok / len(sars_clones), 4),
                    "virus_recall": round(sars_sp / len(sars_clones), 4),
                }
            )

    # Also evaluate unmodified default `prioritize --context blood` viral-exclusion false-drop rate on Minervina SARS-CoV-2 clones
    heldout_min_tsv = work_dir / "vdjdb_heldout_no_minervina.tsv"
    if not heldout_min_tsv.exists():
        pd.DataFrame(
            {
                "cdr3.alpha": db_no_minervina["cdr3_alpha"],
                "cdr3.beta": db_no_minervina["cdr3_beta"],
                "v.alpha": db_no_minervina.get("v_alpha", pd.NA),
                "v.beta": db_no_minervina.get("v_beta", pd.NA),
                "species": db_no_minervina.get("host_species", "HomoSapiens"),
                "antigen.epitope": db_no_minervina["epitope"],
                "antigen.gene": db_no_minervina["antigen_gene"],
                "antigen.species": db_no_minervina["species"],
                "mhc.a": db_no_minervina.get("mhc_allele", pd.NA),
                "mhc.class": db_no_minervina.get("mhc_class", pd.NA),
                "reference.id": db_no_minervina.get("reference", pd.NA),
            }
        ).to_csv(heldout_min_tsv, sep="\t", index=False)

    empty_cache = work_dir / "_empty_db_cache"
    empty_cache.mkdir(parents=True, exist_ok=True)
    min_shortlist_rows: list[dict] = []
    min_ep_lookup = clone_ep.set_index("CDR3ab")["epitope"]
    for run_tag, extra_cli in (
        ("blood_default", ["--vdjdb", str(heldout_min_tsv), "--cache-dir", str(empty_cache)]),
        (
            "blood_exclude_viral_strict_ab",
            ["--vdjdb", str(heldout_min_tsv), "--cache-dir", str(empty_cache), "--exclude-known-viral"],
        ),
        (
            "blood_exclude_viral_ab_with_partial",
            [
                "--vdjdb",
                str(heldout_min_tsv),
                "--cache-dir",
                str(empty_cache),
                "--exclude-known-viral",
                "--database-match",
                "ab_with_partial",
            ],
        ),
    ):
        out_mdef = work_dir / f"minervina_{run_tag}"
        if not is_current(out_mdef):
            tcrsift_main(["prioritize", str(sheet), "-o", str(out_mdef), "--context", "blood", *extra_cli])
        mdef_clones = pd.read_csv(out_mdef / "all_scored_clones.csv", low_memory=False)
        mdef_clones["sars_ep"] = mdef_clones["CDR3ab"].map(min_ep_lookup)
        m_exp = mdef_clones[mdef_clones["cell_count"] >= 2].copy()
        is_sars = m_exp["sars_ep"].notna()
        elig = truthy_flag_series(m_exp["eligible_for_review"])
        sel = truthy_flag_series(m_exp["selected_for_review"])
        min_shortlist_rows.append(
            {
                "dataset": "minervina_2022",
                "run": run_tag,
                "expanded_clones_ge2": len(m_exp),
                "eligible_total": int(elig.sum()),
                "selected_total": int(sel.sum()),
                "viral_total": int(is_sars.sum()),
                "viral_eligible": int((is_sars & elig).sum()),
                "viral_selected": int((is_sars & sel).sum()),
                "flu_total": 0,
                "flu_eligible": 0,
                "flu_selected": 0,
                "ebv_total": 0,
                "ebv_eligible": 0,
                "ebv_selected": 0,
                "mart1_expanded_total": 0,
                "mart1_expanded_eligible": 0,
                "mart1_expanded_selected": 0,
                "unlabeled_expanded_total": int((~is_sars).sum()),
                "unlabeled_expanded_eligible": int(((~is_sars) & elig).sum()),
                "unlabeled_expanded_selected": int(((~is_sars) & sel).sum()),
            }
        )

    ep_counts = sars_clones["epitope"].value_counts()
    facts = {
        "cells_after_gex_qc": cfg["cell_filters"][0]["before"],
        "cells_retained": cfg["cell_filters"][-1]["after"],
        "clones_scored": len(clones),
        "sars_cov_2_labelled_clones": len(sars_clones),
        "n_distinct_epitopes": int(sars_clones["epitope"].nunique()),
        "top_epitope": str(ep_counts.index[0]) if len(ep_counts) else "",
        "top_epitope_clones": int(ep_counts.iloc[0]) if len(ep_counts) else 0,
        "shared_longitudinal_clones_excluded_from_auroc": len(shared_longitudinal),
        "leave_study_out_ab_with_partial_off_target_species": off_target_species,
    }
    return auroc_rows, exposure_rows, vdjdb_rows, min_shortlist_rows, facts


# ------------------------------------------------ 3. Fischer et al. 2021

def run_fischer(data_dir: Path, work_dir: Path):
    from tcrsift.assemble import assemble_full_sequences
    from tcrsift.clonotype import aggregate_clonotypes
    from tcrsift.filter import score_stimulation_response
    from tcrsift.loader import load_samples
    from tcrsift.phenotype import apply_cd3_gate, phenotype_cells
    from tcrsift.sample_sheet import Sample, SampleSheet
    from tcrsift.signature_methods import score_by_name

    fdir = data_dir / "fischer_2021_covid"
    donor_pairs = [
        ("GT_6", "GT6_unstim", "GT6_covid_stim"),
        ("GT_2", "GT_unstim", "GT_covid_stim"),
    ]
    stim_rows = []
    facts = {}

    for donor, s_unstim, s_stim in donor_pairs:
        # Measure raw 5'-partial (full_length == False) rate across productive/rescued TCR contigs
        raw_partial_counts = []
        for s_lbl in (s_unstim, s_stim):
            raw_csv = pd.read_csv(next((fdir / s_lbl).glob("*filtered_contig_annotations.csv*")))
            tcr_rows = raw_csv[
                raw_csv["is_cell"].astype(str).str.lower().eq("true")
                & raw_csv["high_confidence"].astype(str).str.lower().eq("true")
                & raw_csv["chain"].isin(["TRA", "TRB"])
                & ~raw_csv["productive"].astype(str).str.lower().eq("false")
            ]
            if len(tcr_rows) and "full_length" in tcr_rows.columns:
                raw_partial_counts.append(
                    float((~tcr_rows["full_length"].astype(str).str.lower().eq("true")).mean())
                )

        sheet = SampleSheet(
            samples=[
                Sample(
                    sample=s_unstim,
                    patient_id=donor,
                    gex_dir=str(fdir / s_unstim),
                    vdj_dir=str(fdir / s_unstim),
                    source="control",
                    tissue="blood",
                ),
                Sample(
                    sample=s_stim,
                    patient_id=donor,
                    gex_dir=str(fdir / s_stim),
                    vdj_dir=str(fdir / s_stim),
                    source="culture",
                    control_sample=s_unstim,
                    antigen_type="peptide_pool",
                    antigen_name="SARS-CoV-2 Spike",
                    tissue="blood",
                ),
            ]
        )
        adata_raw = load_samples(sheet, min_mito_pct=0.0, max_mito_pct=15.0, show_progress=False, verbose=False)
        adata_pheno = phenotype_cells(adata_raw, min_cd3_reads=3)
        # Gate to CD3-passing T cells before log-normalization and signature/cytokine scoring
        # so non-T PBMC populations (monocytes, NK, B cells) do not bias reference bins or TNF/CCL4.
        adata = apply_cd3_gate(adata_pheno, min_cd3_reads=3)

        norm = adata.copy()
        sc.pp.normalize_total(norm, target_sum=10_000)
        sc.pp.log1p(norm)
        for sname in ("AIM", "AIMBroad", "AcuteActivation", "Differentiated", "Cytolytic"):
            norm.obs[f"sig_{sname}"] = score_by_name(norm, sname, log1p=False, species="human").to_numpy()

        sym_map = dict(zip(norm.var["gene_symbols"], norm.var_names))
        cyto_genes = [g for g in ("IFNG", "TNF", "IL2", "CCL4") if g in sym_map]
        cyto_mat = norm[:, [sym_map[g] for g in cyto_genes]].X
        cyto_arr = cyto_mat.toarray() if hasattr(cyto_mat, "toarray") else np.asarray(cyto_mat)
        norm.obs["effector_cytokines"] = cyto_arr.mean(axis=1)

        resp = score_stimulation_response(
            norm,
            stim_samples=[s_stim],
            control_samples=[s_unstim],
            clone_col="CDR3_beta",
            score_cols=["sig_AIM", "sig_AIMBroad", "sig_AcuteActivation", "sig_Differentiated", "effector_cytokines"],
        )
        # Persistent clones observed with >=2 CD3-gated T cells in both unstimulated and stimulated wells
        pers = resp[(resp["stim_cells"] >= 2) & (resp["control_cells"] >= 2)].copy()
        for _, r in pers.iterrows():
            stim_rows.append(
                {
                    "donor": donor,
                    "CDR3_beta": r["CDR3_beta"],
                    "stim_cells": int(r["stim_cells"]),
                    "control_cells": int(r["control_cells"]),
                    "log2_fold_change": round(float(r["log2_fold_change"]), 4),
                    "fisher_pvalue": round(float(r["fisher_pvalue"]), 4),
                    "fisher_fdr": round(float(r["fisher_fdr"]), 4),
                    "stim_effector_cytokines": round(float(r["stim_mean_effector_cytokines"]), 4),
                    "control_effector_cytokines": round(float(r["control_mean_effector_cytokines"]), 4),
                    "delta_effector_cytokines": round(float(r["delta_effector_cytokines"]), 4),
                    "delta_AIMBroad": round(float(r["delta_sig_AIMBroad"]), 4),
                    "delta_AIM": round(float(r["delta_sig_AIM"]), 4),
                    "delta_AcuteActivation": round(float(r["delta_sig_AcuteActivation"]), 4),
                }
            )

        clones_ab = aggregate_clonotypes(adata, group_by="CDR3ab", show_progress=False, verbose=False)
        has_vdj = (
            clones_ab["VDJ_alpha_aa"].notna()
            & (clones_ab["VDJ_alpha_aa"].astype(str) != "")
            & clones_ab["VDJ_beta_aa"].notna()
            & (clones_ab["VDJ_beta_aa"].astype(str) != "")
            if {"VDJ_alpha_aa", "VDJ_beta_aa"} <= set(clones_ab.columns)
            else pd.Series(False, index=clones_ab.index)
        )
        n_assembled_ok = 0
        n_contig_verified = 0
        if has_vdj.any():
            asm = assemble_full_sequences(
                clones_ab[has_vdj].copy(),
                contigs_dir=fdir,
                sample_name_from="parent",
                alpha_leader="from_contig",
                beta_leader="from_contig",
                include_constant=True,
                linker="T2A",
                show_progress=False,
                verbose=False,
            )
            n_assembled_ok = int(
                (asm["single_chain_aa"].notna() & asm["qc_warnings"].map(len).eq(0)).sum()
            )
            n_contig_verified = int(asm["construct_contig_verified"].sum())

        has_partial_or_full = pd.Series(True, index=clones_ab.index)
        for ch in ("alpha", "beta"):
            c_full = (
                clones_ab[f"VDJ_{ch}_aa"].notna() & (clones_ab[f"VDJ_{ch}_aa"].astype(str) != "")
                if f"VDJ_{ch}_aa" in clones_ab.columns
                else pd.Series(False, index=clones_ab.index)
            )
            c_part = (
                clones_ab[f"partial_VDJ_{ch}_aa"].notna()
                & (clones_ab[f"partial_VDJ_{ch}_aa"].astype(str) != "")
                if f"partial_VDJ_{ch}_aa" in clones_ab.columns
                else pd.Series(False, index=clones_ab.index)
            )
            has_partial_or_full &= c_full | c_part

        n_germline_rescued = 0
        n_germline_5prime_both = 0
        n_germline_junction_contig_both = 0
        asm_germline = assemble_full_sequences(
            clones_ab.copy(),
            contigs_dir=fdir,
            sample_name_from="parent",
            alpha_leader="from_contig",
            beta_leader="from_contig",
            leader_fallback="germline",
            include_constant=True,
            linker="T2A",
            germline_vdj_fallback=True,
            show_progress=False,
            verbose=False,
        )
        valid_sc = asm_germline["single_chain_aa"].notna() & ~asm_germline["single_chain_aa"].astype(
            str
        ).str.contains(r"\*", regex=True, na=True)
        n_germline_rescued = int(valid_sc.sum())
        if {"alpha_vdj_source", "beta_vdj_source"} <= set(asm_germline.columns):
            n_germline_5prime_both = int(
                (
                    valid_sc
                    & asm_germline["alpha_vdj_source"].isin(
                        ["contig_full", "contig_5prime_germline_rescued"]
                    )
                    & asm_germline["beta_vdj_source"].isin(
                        ["contig_full", "contig_5prime_germline_rescued"]
                    )
                ).sum()
            )
        if {"alpha_junction_residue_source", "beta_junction_residue_source"} <= set(
            asm_germline.columns
        ):
            n_germline_junction_contig_both = int(
                (
                    valid_sc
                    & asm_germline["alpha_junction_residue_source"].eq("contig")
                    & asm_germline["beta_junction_residue_source"].eq("contig")
                ).sum()
            )

        facts[donor] = {
            "total_pbmc_loaded": int(adata_raw.n_obs),
            "cd3_gated_t_cells": int(adata.n_obs),
            "cd4_cells": int(adata.obs["is_CD4"].sum()),
            "cd8_cells": int(adata.obs["is_CD8"].sum()),
            "five_prime_truncated_contig_fraction_min": round(min(raw_partial_counts), 4),
            "five_prime_truncated_contig_fraction_max": round(max(raw_partial_counts), 4),
            "paired_ab_clones": len(clones_ab),
            "beta_clones_evaluated": len(resp),
            "persistent_clones_ge2_both": len(pers),
            "cytokine_activated_persistent_clones_delta_gt_0_5": int(
                (pers["delta_effector_cytokines"] > 0.5).sum()
            ),
            "fisher_expanded_persistent_clones_fdr_lt_0_05": int(
                (pers["fisher_fdr"] < 0.05).sum()
            ),
            "json_full_vdj_paired_clones": int(has_vdj.sum()),
            "json_assembled_valid_constructs": n_assembled_ok,
            "json_contig_verified_constructs": n_contig_verified,
            "partial_or_full_vdj_paired_clones": int(has_partial_or_full.sum()),
            "germline_fallback_assembled_constructs": n_germline_rescued,
            "germline_fallback_5prime_contig_rescued_constructs": n_germline_5prime_both,
            "germline_fallback_contig_junction_verified_constructs": n_germline_junction_contig_both,
        }

    return stim_rows, facts


GROUND_TRUTH_PROVENANCE = {
    "tenx_dextramer": {
        "citation": "10x Genomics Application Note (2019): 'A new way of exploring immunity: linking highly multiplexed antigen recognition to immune repertoire and phenotype' (Cell Ranger 3.0.2 Healthy Donors 1 & 2)",
        "setting": "Resting steady-state peripheral blood CD8+ T cells from healthy HLA-typed donors",
        "gex_vdj_accession": "https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor{1,2}/ (*_sc5p_v2_hs_B_T_all_feature_bc_matrix_filtered.h5, saved locally as filtered_feature_bc_matrix.h5, + *_all_contig_annotations.csv)",
        "label_source": "10x Genomics published binarized_matrix.csv (44 viral/cancer pMHC dextramer '*_binder' columns) + raw Antibody Capture UMI counts for 6 negative-control ('*_NC') dextramers and TotalSeq-C CD45RO/CD45RA surface antibodies.",
        "experimental_assay": "Simultaneous staining with 50 DNA-barcoded dCODE MHC class I dextramers (44 viral/self epitopes + 6 negative-control peptides) + CITE-seq TotalSeq-C antibodies + 5' scRNA-seq + scTCR-seq.",
        "positive_definition": "Paired CDR3ab clonotypes where >=50% of cells in the clone carry the same single-specificity binder call in binarized_matrix.csv (plus >=2 single-binder cells for all expanded >=2-cell evaluations): bona fide viral memory = Influenza A M1 (GILGFVFTL), EBV (BMLF1/BZLF1/EBNA-3A/3B/LMP2A), and CMV pp65 (NLVPMVATV). Sticky reagent A0301_KLGGALQAK_IE-1_CMV is excluded from ground-truth viral positives and evaluated separately as a multimer QC stress test.",
        "negative_definition": "(A) Naive antigen-specific precursors: HLA-A*02:01 MART-1 (ELAGIGILTV) dextramer+ clones (verified phenotypically naive via TotalSeq-C CLR(CD45RO) - CLR(CD45RA) < -2.9 and >99% singletons). (B) Unlabeled expanded blood clones: paired CD8 clones (>=2 cells) negative for all 44 panel dextramers (representing memory T cells against unpanelled antigens). (C) Memory subtype contrast: Influenza M1 GILGFVFTL vs EBV expanded clones.",
        "circularity_controls": "10x Genomics dextramer data were ingested into VDJdb (reference containing '10xgenomics'); VDJdb matching reports both unmodified VDJdb (circular) and 10x-study-held-out VDJdb (independent).",
    },
    "minervina_2022": {
        "citation": "Minervina, Pogorelyy, Kirk et al., Nature Immunology 2022 (23:781-790; DOI:10.1038/s41590-022-01184-4; PMID:35383307)",
        "setting": "Longitudinal human peripheral blood SARS-CoV-2 CD8+ T cells across pre-exposure, convalescence, mRNA vaccination, acute post-dose-2 recall, and breakthrough infection",
        "gex_vdj_accession": "Zenodo 6232103 (libs_aggregate_rev.zip: aggregated Cell Ranger filtered_feature_bc_matrix.h5 + filtered_contig_annotations.csv)",
        "label_source": "Curated single-cell dextramer & cohort metadata table pogorely/COVID_vax_CD8 (cd8_only_dextr_rev_clean.tsv; columns 'barcode', 'donor', 'epitope', 'spike', 'donor_category', 'sample_category').",
        "experimental_assay": "FACS and DNA-barcoded HLA class I dextramer profiling across 18 SARS-CoV-2 Spike and non-Spike (ORF1ab, N, M, ORF3a) epitopes.",
        "positive_definition": "(1) State recall AUROC: paired-TCR CD8 clones (>=2 cells) sampled during acute post-vaccine dose-2 recall (donor_category == 'inf_after_2'). (2) VDJdb specificity recovery: paired-TCR CD8 clones carrying a curated SARS-CoV-2 epitope call in cd8_only_dextr_rev_clean.tsv.",
        "negative_definition": "State recall AUROC: disjoint paired-TCR SARS-CoV-2 CD8 clones (>=2 cells) sampled during resting post-infection convalescence (donor_category == 'post_inf'), excluding longitudinal clones observed in both categories.",
        "circularity_controls": "PMID:35383307 is indexed in VDJdb; VDJdb matching reports both unmodified VDJdb (circular) and PMID:35383307-held-out VDJdb (independent).",
    },
    "fischer_2021": {
        "citation": "Fischer, Ansari, Wagner et al., Nature Communications 2021 (12:4515; DOI:10.1038/s41467-021-25548-w; PMID:34301934)",
        "setting": "Paired unstimulated control vs 16-hour SARS-CoV-2 Spike peptide-pool-stimulated PBMCs from hospitalized severe COVID-19 patients (Gauting cohort)",
        "gex_vdj_accession": "GEO GSE171037 (GSM5216975 GT6_unstim, GSM5216976 GT6_covid_stim, GSM5216977 GT_unstim, GSM5216978 GT_covid_stim: Cell Ranger 3 .h5 + .csv.gz + all_contig_annotations.json.gz)",
        "label_source": "Experimental well assignment from GEO sample metadata (unstimulated control well vs autologous 16 h SARS-CoV-2 Spike overlapping peptide pool well) across persistent CDR3_beta clonotypes (>=2 CD3-gated T cells in both wells).",
        "experimental_assay": "Ex vivo 16 h overlapping SARS-CoV-2 Spike peptide pool restimulation ('reverse phenotyping') vs parallel unstimulated PBMC control from the same blood draw.",
        "positive_definition": "Persistent T-cell clones exhibiting paired within-clone effector cytokine transcript induction (delta mean log1p CP10K of IFNG, TNF, IL2, CCL4 > 0.5 and/or AIM/AIMBroad induction) upon Spike peptide stimulation.",
        "negative_definition": "Persistent bystander T-cell clones present at >=2 cells in both unstimulated and Spike-stimulated wells showing no transcriptional cytokine/AIM induction.",
        "circularity_controls": "None of TCRsift's signatures (AIM, AIMBroad, AcuteActivation, Differentiated, Cytolytic) were derived from GSE171037.",
    },
}


def main() -> None:
    warnings.filterwarnings("ignore", message="Variable names are not unique")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--vdjdb", type=Path, default=None, help="Optional path to VDJdb dir/file")
    args = parser.parse_args()

    from tcrsift import __version__
    from tcrsift.annotate import _pick_vdjdb_source, load_vdjdb
    from tcrsift.datacache import ensure_cached

    args.work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    vdjdb_path = args.vdjdb or ensure_cached("vdjdb", auto_download=False)
    if vdjdb_path is None:
        raise RuntimeError(
            "VDJdb is required in the TCRsift cache (`tcrsift data download --db vdjdb` or --vdjdb)."
        )
    vdjdb_src = _pick_vdjdb_source(Path(vdjdb_path)) if Path(vdjdb_path).is_dir() else Path(vdjdb_path)
    vdjdb_src_label = (
        vdjdb_src[0].parent.name if isinstance(vdjdb_src, tuple) else vdjdb_src.name
    )
    vdjdb_full = load_vdjdb(vdjdb_path, verbose=False)

    print("Running Cohort 1: 10x Genomics Multiplexed Dextramer Donors 1 & 2 ...")
    tenx_auroc, tenx_vdjdb, multimer_rows, tenx_shortlist, tenx_facts = run_tenx_dextramer(
        args.data_dir, args.work_dir, vdjdb_full
    )

    print("Running Cohort 2: Minervina et al. 2022 SARS-CoV-2 CD8 Dextramers ...")
    min_auroc, exposure_rows, min_vdjdb, min_shortlist, min_facts = run_minervina(
        args.data_dir, args.work_dir, vdjdb_full
    )

    print("Running Cohort 3: Fischer et al. 2021 COVID-19 PBMC Stimulation ...")
    fischer_rows, fischer_facts = run_fischer(args.data_dir, args.work_dir)

    pd.DataFrame(tenx_auroc + min_auroc).to_csv(
        results_dir / "blood_signature_auroc.csv", index=False
    )
    pd.DataFrame(tenx_shortlist + min_shortlist).to_csv(
        results_dir / "blood_cutoffs_and_shortlist.csv", index=False
    )
    pd.DataFrame(exposure_rows).to_csv(
        results_dir / "exposure_state_profiles.csv", index=False
    )
    pd.DataFrame(tenx_vdjdb + min_vdjdb).to_csv(
        results_dir / "vdjdb_matching.csv", index=False
    )
    pd.DataFrame(multimer_rows).to_csv(
        results_dir / "multimer_calling.csv", index=False
    )
    pd.DataFrame(fischer_rows).to_csv(
        results_dir / "fischer_stimulation.csv", index=False
    )

    summary = {
        "tcrsift_version": __version__,
        "vdjdb_release": vdjdb_src_label,
        "vdjdb_entries": len(vdjdb_full),
        "vdjdb_paired_ab_entries": int(vdjdb_full["cdr3_alpha"].notna().sum()),
        "ground_truth_provenance": GROUND_TRUTH_PROVENANCE,
        "tenx_dextramer": tenx_facts,
        "minervina_2022": min_facts,
        "fischer_2021": fischer_facts,
    }
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote infectious-disease benchmark results to {results_dir}")


if __name__ == "__main__":
    main()
