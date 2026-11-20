#!/usr/bin/env python3
"""Run the public TIL prioritization benchmark and write small result tables.

    python benchmarks/public_til/download.py --data-dir CACHE
    python benchmarks/public_til/run.py --data-dir CACHE --work-dir CACHE/runs

Two datasets, two questions:

* Caushi 2021 (NSCLC): ``tcrsift prioritize --context solid-tumor`` is run
  unchanged from a sample sheet. Labels are CDR3-beta sequences the authors
  showed to be neoantigen (MANA) specific or virus specific.
* Meng 2023 (PDAC): GEO has no VDJ contigs, only a per-cell clonotype id, so
  the same scoring, background calibration and cutoff functions are called
  directly on the authors' CD8 cells. Labels are clonotypes the authors
  tested as tumor reactive (TR) or not (NTR).

For both we ask (1) how well each signature's within-sample clone percentile
separates positives from negatives (AUROC) and (2) how many positives and
negatives pass TCRsift's default cutoffs. Results go to ``results/``.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
CAUSHI_RUNS = {
    "default": [],
    # The default floor of ten CD3 counts removes most CD8-classified cells in
    # this dataset. Same command with that gate off.
    "no_cd3_gate": ["--min-cd3", "0"],
}


def auroc(scores, positive):
    """Probability a random positive outranks a random negative (ties = 0.5)."""
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
    n_pos, n_neg = positive.sum(), (~positive).sum()
    if not n_pos or not n_neg:
        return np.nan
    ranks = rankdata(scores)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auroc_ci(scores, positive, draws=2000, seed=0):
    """Percentile bootstrap over clones, resampling each class separately."""
    scores, positive = np.asarray(scores, float), np.asarray(positive, bool)
    pos, neg = scores[positive], scores[~positive]
    if len(pos) < 2 or len(neg) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(draws):
        p, n = rng.choice(pos, len(pos)), rng.choice(neg, len(neg))
        values.append(auroc(np.concatenate([p, n]), np.arange(len(p) + len(n)) < len(p)))
    return tuple(float(v) for v in np.quantile(values, [0.025, 0.975]))


def summarize(table, signatures, dataset, run, positive_label, negative_label):
    """Per-signature AUROC and cutoff pass rates for one labelled clone table."""
    labelled = table[table.label.isin([positive_label, negative_label])].copy()
    positive = labelled.label.eq(positive_label).to_numpy()
    rows = []
    for name in [*signatures, "clone_size"]:
        col = "size_percentile" if name == "clone_size" else f"signature_{name}_percentile"
        ok = labelled[col].notna().to_numpy()
        low, high = auroc_ci(labelled.loc[ok, col], positive[ok])
        row = {"dataset": dataset, "run": run, "signature": name,
               "n_positive": int(positive[ok].sum()), "n_negative": int((~positive[ok]).sum()),
               "auroc": auroc(labelled.loc[ok, col], positive[ok]), "auroc_low": low, "auroc_high": high}
        if name != "clone_size":
            passed = labelled[f"signature_{name}_passes_cutoff"].fillna(False).astype(bool).to_numpy()
            row |= {"positive_pass": int(passed[positive].sum()), "negative_pass": int(passed[~positive].sum())}
        rows.append(row)
    return rows


# ---------------------------------------------------------------- Caushi 2021

def caushi_labels(data_dir):
    raw = pd.read_excel(data_dir / "caushi" / "supplementary_tables.xlsx",
                        sheet_name="Supplementary Table 8", header=None)
    table = raw.iloc[3:, :3].dropna()
    table.columns = ["donor", "CDR3_beta", "antigen_type"]
    table["CDR3_beta"] = table.CDR3_beta.str.replace("*", "", regex=False).str.strip()
    table["label"] = np.where(table.antigen_type.str.contains("MANA", case=False), "MANA", "viral")
    # A handful of sequences are listed for two epitopes; one row per clone.
    return table.drop_duplicates(["donor", "CDR3_beta"])[["donor", "CDR3_beta", "label"]]


def run_caushi(data_dir, work_dir):
    from tcrsift.cli import main as tcrsift_main

    samples = sorted((data_dir / "caushi" / "samples").iterdir())
    sheet = work_dir / "caushi_samples.csv"
    lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    for sample in samples:
        patient, tissue, _ = sample.name.split("_")
        lines.append(f"{sample.name},{patient},{tissue},{sample / 'vdj'},{sample / 'gex'}")
    sheet.write_text("\n".join(lines) + "\n")
    labels = caushi_labels(data_dir)

    summary, tables, facts = [], [], {}
    for run, extra in CAUSHI_RUNS.items():
        out = work_dir / f"caushi_{run}"
        if not (out / "prioritization.json").exists():
            tcrsift_main(["prioritize", str(sheet), "-o", str(out), "--context", "solid-tumor", *extra])
        config = json.loads((out / "prioritization.json").read_text())
        signatures = config["resolved_signatures"]
        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones.groupby("donor").cell_count.rank(pct=True)
        clones = clones.merge(labels, on=["donor", "CDR3_beta"], how="left")
        clones["label"] = clones.label.fillna("untested")
        clones = clones.rename(columns={
            f"signature_{n}{suffix}_max": f"signature_{n}{suffix}"
            for n in signatures for suffix in ("_percentile", "_passes_cutoff")})
        keep = ["donor", "CDR3_beta", "CDR3_alpha", "label", "cell_count", "size_percentile",
                "eligible_for_review", "selected_for_review", "selection_rank", "selected_signature",
                *[f"signature_{n}{s}" for n in signatures for s in ("_percentile", "_passes_cutoff")]]
        table = clones[keep].assign(run=run)
        # Labels are beta chains. One labelled beta can appear with several
        # alpha chains (dual alphas, dropout); count it once, represented by
        # its largest alpha-beta clone, and as selected if any of them was.
        labelled = table[table.label.ne("untested")].sort_values("cell_count", ascending=False)
        any_cols = ["eligible_for_review", "selected_for_review", *[c for c in keep if c.endswith("_passes_cutoff")]]
        grouped = labelled.groupby(["donor", "CDR3_beta"], sort=False)
        collapsed = grouped.head(1).drop(columns=any_cols).merge(
            grouped[any_cols].any().reset_index(), on=["donor", "CDR3_beta"])
        collapsed["alpha_beta_clones"] = collapsed.set_index(["donor", "CDR3_beta"]).index.map(grouped.size())
        table = pd.concat([table[table.label.eq("untested")], collapsed], ignore_index=True)
        tables.append(table)
        summary += summarize(table, signatures, "caushi", run, "MANA", "viral")
        filters = config["cell_filters"]
        facts[run] = {
            "samples": len(samples), "patients": int(clones.donor.nunique()),
            "cells_after_gex_qc": filters[0]["before"], "cells_retained": filters[-1]["after"],
            "cell_filters": filters, "clones_scored": len(clones),
            "candidates": int(clones.selected_for_review.sum()),
            "labelled_clones_scored": table.label.value_counts().to_dict(),
            "labelled_clones_published": labels.label.value_counts().to_dict(),
            "selected": table[table.selected_for_review].label.value_counts().to_dict(),
            "eligible": table[table.eligible_for_review].label.value_counts().to_dict(),
        }
    return summary, pd.concat(tables, ignore_index=True), facts


# ------------------------------------------------------------------ Meng 2023

def run_meng(data_dir):
    import anndata as ad
    import scanpy as sc

    from tcrsift.prioritize import resolve_signatures, signature_pass_mask
    from tcrsift.signature_background import calibrate_signature_background
    from tcrsift.til_prioritize import _score_within_samples

    root = data_dir / "meng"
    meta = pd.read_csv(root / "GSE254250_PDAC_9sample_CD8_RNA_metadata_matrix.csv.gz").iloc[:, [0, 1, 5, 8]]
    meta.columns = ["cell", "patient", "clonotype", "label"]
    meta = meta.dropna(subset=["clonotype"])
    meta["barcode"] = meta.cell.str.split("_").str[1]
    meta["label"] = meta.label.replace("-", "untested")
    files = {patient: "Pool" if patient in ("TIPC301", "TIPC413", "TIPC416") else patient[4:]
             for patient in sorted(meta.patient.unique())}
    parts = []
    for patient, tag in files.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            matrix = sc.read_10x_h5(root / f"GSE254250_{tag}_GEX_filtered_feature_bc_matrix.h5")
        matrix.var_names_make_unique()
        # These matrices use the GRCh38-2020-A annotation, where PLA2G16 is
        # named PLAAT3. NeoTCR8 lists the older symbol and TCRsift refuses to
        # score a signature with an absent gene, so restore the published name.
        matrix.var_names = matrix.var_names.where(matrix.var_names != "PLAAT3", "PLA2G16")
        cells = meta[meta.patient.eq(patient)].set_index("barcode")
        part = matrix[cells.index].copy()
        part.obs = pd.DataFrame({"donor": patient, "sample": patient, "lineage": "cd8",
                                 "CDR3ab": cells.clonotype.to_numpy(), "label": cells.label.to_numpy()},
                                index=cells.cell.to_numpy())
        parts.append(part)
    adata = ad.concat(parts)
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    signatures = resolve_signatures("solid-tumor", "cd8")
    columns = [f"signature_{name}" for name in signatures]
    for name, col in zip(signatures, columns):
        adata.obs[col] = _score_within_samples(adata, name)
    stratum = ["donor", "sample", "lineage"]
    per_cell = adata.obs[[*stratum, "CDR3ab", *columns]]
    groups = per_cell.groupby([*stratum, "CDR3ab"], observed=True)
    scores = groups[columns].mean().join(groups.size().rename("cells")).reset_index()
    scores["frequency"] = scores.cells / scores.groupby(stratum).cells.transform("sum")
    for col in columns:
        scores[f"{col}_percentile"] = scores.groupby(stratum)[col].rank(method="average", pct=True)
    scores, _ = calibrate_signature_background(scores, per_cell, signatures)
    for name in signatures:
        for route in ("mean", "subset", "either"):
            suffix = "passes_cutoff" if route == "either" else f"{route}_passes_cutoff"
            scores[f"signature_{name}_{suffix}"] = signature_pass_mask(scores, name, evidence=route)
    scores["size_percentile"] = scores.groupby(stratum).cells.rank(pct=True)
    labels = adata.obs.groupby("CDR3ab").label.first()
    scores["label"] = scores.CDR3ab.map(labels)
    table = scores.rename(columns={"CDR3ab": "clonotype", "cells": "cell_count"})
    table["passes_any"] = table[[f"signature_{n}_passes_cutoff" for n in signatures]].any(axis=1)
    summary = summarize(table, signatures, "meng", "default", "TR", "NTR")
    facts = {
        "patients": int(table.donor.nunique()), "cells": int(adata.n_obs), "clones_scored": len(table),
        "labelled_clones": table.label.value_counts().to_dict(),
        "passes_any_default_cutoff": table[table.passes_any].label.value_counts().to_dict(),
        "passes_any_default_cutoff_min2_cells":
            table[table.passes_any & table.cell_count.ge(2)].label.value_counts().to_dict(),
        "min2_cells": table[table.cell_count.ge(2)].label.value_counts().to_dict(),
    }
    keep = ["donor", "clonotype", "label", "cell_count", "size_percentile", "passes_any",
            *[f"signature_{n}{s}" for n in signatures
              for s in ("_percentile", "_passes_cutoff", "_mean_passes_cutoff", "_subset_passes_cutoff")]]
    per_patient = []
    for patient, group in table.groupby("donor"):
        for name in signatures:
            g = group[group.label.isin(["TR", "NTR"])]
            per_patient.append({"donor": patient, "signature": name, "n_positive": int(g.label.eq("TR").sum()),
                                "n_negative": int(g.label.eq("NTR").sum()),
                                "auroc": auroc(g[f"signature_{name}_percentile"], g.label.eq("TR"))})
    return summary, table[keep], facts, pd.DataFrame(per_patient)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path, help="tcrsift output directories")
    parser.add_argument("--out-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    from tcrsift import __version__

    caushi_summary, caushi_table, caushi_facts = run_caushi(args.data_dir, args.work_dir)
    meng_summary, meng_table, meng_facts, meng_patients = run_meng(args.data_dir)
    pd.DataFrame(caushi_summary + meng_summary).to_csv(args.out_dir / "signature_auroc.csv", index=False)
    labelled = caushi_table[caushi_table.label.ne("untested")]
    labelled.to_csv(args.out_dir / "caushi_labelled_clones.csv", index=False)
    meng_table[meng_table.label.ne("untested")].to_csv(args.out_dir / "meng_labelled_clones.csv", index=False)
    meng_patients.to_csv(args.out_dir / "meng_per_patient_auroc.csv", index=False)
    (args.out_dir / "summary.json").write_text(json.dumps(
        {"tcrsift_version": __version__, "caushi": caushi_facts, "meng": meng_facts}, indent=2, default=int) + "\n")
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
