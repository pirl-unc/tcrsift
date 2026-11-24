#!/usr/bin/env python3
"""Run the public TIL prioritization benchmark and write small result tables.

    python benchmarks/public_til/download.py --data-dir CACHE
    python benchmarks/public_til/run.py --data-dir CACHE --work-dir CACHE/runs

Datasets and how TCRsift is exercised on each:

* Caushi 2021 (NSCLC, GSE176021): the unmodified ``tcrsift prioritize
  --context solid-tumor`` command on sample sheets. Labels are CDR3-beta
  sequences shown to be neoantigen (MANA) specific or virus specific.
* Eberhardt 2021 (HPV+ head and neck, GSE180268): the same command on bulk
  PD-1+ CD8 TIL. Positives are beta chains found in HPV-tetramer-sorted cells
  of the same patient; there are no tested negatives, so they are compared
  with every other (unlabelled) clone.
* Chiffelle 2024 (melanoma, GSE222448): the same command on baseline tumors.
  Labels are clonotypes tested as tumor reactive (TR) or not (NTR).
* Pulliam 2024 (Merkel cell carcinoma, GSE227708): the same command on the
  tumor cells of two hashed runs. Positives are beta chains of clones that
  bound Merkel cell polyomavirus multimers; there are no tested negatives in
  the same tumor.
* Meng 2023 (PDAC, GSE254250): GEO has no VDJ contigs, only a per-cell
  clonotype id, so the same scoring, background calibration, cutoff and
  round-robin functions are called on the authors' CD8 cells. Labels are
  clonotypes tested as tumor reactive (TR) or not (NTR).

Questions answered, each written to its own table in ``results/``:

1. Ranking: per-signature AUROC of the within-sample clone percentile.
2. Cutoffs: who passes the default cutoffs, by evidence route and signature
   set, and how that moves with the two main thresholds.
3. Selection: what the round-robin shortlist contains at several budgets,
   against choosing the largest clones or one signature's top clones.
4. Cell gates: what the CD3 floor does to cells, clones and the shortlist.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
STRATUM = ["donor", "sample", "lineage"]
TUMOR_TISSUES = ("tumor", "mettumor")
# (run name, tissues, extra arguments). The first is the default command.
CAUSHI_RUNS = [
    ("tumor_default", TUMOR_TISSUES, []),
    ("tumor_cd3_10", TUMOR_TISSUES, ["--min-cd3", "10"]),
    ("tumor_cd3_0", TUMOR_TISSUES, ["--min-cd3", "0"]),
    # Is the CD3 floor special, or is any depth filter as good? No CD3 floor,
    # but require 4,000 UMIs per cell instead of 500.
    ("tumor_depth_gate", TUMOR_TISSUES, ["--min-cd3", "0", "--min-counts", "4000"]),
    ("other_tissue_default", ("LN", "normal"), []),
    ("other_tissue_cd3_10", ("LN", "normal"), ["--min-cd3", "10"]),
]
TUMOR_SIGNATURES = ("TumorReactive", "MANAscore", "NeoTCR8")
# The default command (CD3 floor of 3), the former floor of 10, no floor, and
# the tumor-reactivity signatures alone, for the datasets small enough to rerun cheaply.
SENSITIVITY_RUNS = [
    ("default", []),
    ("cd3_10", ["--min-cd3", "10"]),
    ("cd3_0", ["--min-cd3", "0"]),
    ("tumor_signatures", ["--signatures", *TUMOR_SIGNATURES]),
]
BUDGETS = (25, 50, 100, 200)


def is_current(out):
    """True when ``out`` holds a finished run of the installed tcrsift version.
    Runs left by another version are redone, so every table is from one version."""
    from tcrsift import __version__

    config = out / "prioritization.json"
    return config.exists() and json.loads(config.read_text()).get("version") == __version__


# ------------------------------------------------------------------- metrics

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


def auroc_rows(table, signatures, positive_label, negative_label, **tags):
    """One row per signature (plus clone size) for a labelled clone table."""
    labelled = table[table.label.isin([positive_label, negative_label])]
    positive = labelled.label.eq(positive_label).to_numpy()
    rows = []
    for name in [*signatures, "clone_size"]:
        col = "size_percentile" if name == "clone_size" else f"signature_{name}_percentile"
        ok = labelled[col].notna().to_numpy()
        low, high = auroc_ci(labelled.loc[ok, col], positive[ok])
        rows.append({**tags, "signature": name,
                     "n_positive": int(positive[ok].sum()), "n_negative": int((~positive[ok]).sum()),
                     "auroc": auroc(labelled.loc[ok, col], positive[ok]), "auroc_low": low, "auroc_high": high})
    return rows


def same_clone_rows(labelled, signatures, positive_label, negative_label, default_run, strict_run, **tags):
    """Does a stricter CD3 floor rank clones better, or only score fewer of them?

    Splits the default run's labelled clones into those the stricter run also
    scored and those it dropped, and reports AUROC for each part."""
    key = ["donor", "CDR3_beta"]
    strict = labelled.loc[labelled.run.eq(strict_run), key].drop_duplicates().assign(kept=True)
    default = labelled[labelled.run.eq(default_run)].merge(strict, on=key, how="left")
    kept = default.kept.eq(True)
    rows = []
    for name, part in ((f"{default_run}_also_scored_at_{strict_run}", default[kept]),
                       (f"{default_run}_dropped_at_{strict_run}", default[~kept])):
        rows += [r for r in auroc_rows(part, signatures, positive_label, negative_label, **tags, run=name)
                 if r["signature"] != "clone_size"]
    return rows


def counts(table, mask=None):
    """Label counts as a plain dict, optionally within a boolean mask."""
    labels = table.label if mask is None else table.label[np.asarray(mask, bool)]
    return {str(k): int(v) for k, v in labels.value_counts().items()}


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


def collapse_by_beta(table, any_cols):
    """Labels are beta chains. One labelled beta can appear with several alpha
    chains (dual alphas, dropout): count it once, represented by its largest
    alpha-beta clone, and as selected or passing if any of them was."""
    labelled = table[table.label.ne("untested")].sort_values(
        ["cell_count", "CDR3_alpha"], ascending=[False, True], kind="mergesort"
    )
    grouped = labelled.groupby(["donor", "CDR3_beta"], sort=False)
    collapsed = grouped.head(1).drop(columns=any_cols).merge(
        grouped[any_cols].any().reset_index(), on=["donor", "CDR3_beta"]
    )
    if "selection_rank" in labelled.columns and "selected_signature" in labelled.columns:
        sel = labelled[labelled["selection_rank"].notna()].sort_values(
            ["selection_rank", "CDR3_alpha"], ascending=[True, True], kind="mergesort"
        )
        if not sel.empty:
            best_sel = sel.drop_duplicates(["donor", "CDR3_beta"])[
                ["donor", "CDR3_beta", "selection_rank", "selected_signature"]
            ].rename(
                columns={
                    "selection_rank": "_best_sel_rank",
                    "selected_signature": "_best_sel_sig",
                }
            )
            collapsed = collapsed.merge(best_sel, on=["donor", "CDR3_beta"], how="left")
            collapsed["selection_rank"] = collapsed["selection_rank"].combine_first(
                collapsed.pop("_best_sel_rank")
            )
            collapsed["selected_signature"] = collapsed["selected_signature"].combine_first(
                collapsed.pop("_best_sel_sig")
            )
    return pd.concat([table[table.label.eq("untested")], collapsed], ignore_index=True)


def run_caushi(data_dir, work_dir):
    from tcrsift.cli import main as tcrsift_main

    samples = sorted((data_dir / "caushi" / "samples").iterdir())
    labels = caushi_labels(data_dir)
    aurocs, clone_tables, facts = [], [], {}
    for run, tissues, extra in CAUSHI_RUNS:
        chosen = [s for s in samples if s.name.split("_")[1] in tissues]
        sheet = work_dir / f"caushi_{run}.csv"
        lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
        lines += [f"{s.name},{s.name.split('_')[0]},{s.name.split('_')[1]},{s / 'vdj'},{s / 'gex'}" for s in chosen]
        sheet.write_text("\n".join(lines) + "\n")
        out = work_dir / f"caushi_{run}"
        seconds_file = out / "benchmark_seconds.txt"
        if not is_current(out):
            start = time.perf_counter()
            tcrsift_main(["prioritize", str(sheet), "-o", str(out), "--context", "solid-tumor", *extra])
            seconds_file.write_text(f"{time.perf_counter() - start:.0f}\n")
        config = json.loads((out / "prioritization.json").read_text())
        signatures = config["resolved_signatures"]
        suffixes = ("_percentile", "_passes_cutoff")
        sig_cols = [f"signature_{n}{s}" for n in signatures for s in suffixes]

        clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
        clones["size_percentile"] = clones.groupby("donor").cell_count.rank(pct=True)
        clones = clones.merge(labels, on=["donor", "CDR3_beta"], how="left")
        clones["label"] = clones.label.fillna("untested")
        clones = clones.rename(columns={f"{c}_max": c for c in sig_cols})
        any_cols = ["eligible_for_review", "selected_for_review", *[c for c in sig_cols if c.endswith("cutoff")]]
        keep = ["donor", "CDR3_beta", "CDR3_alpha", "label", "cell_count", "size_percentile",
                "eligible_for_review", "selected_for_review", "selection_rank", "selected_signature", *sig_cols]
        table = collapse_by_beta(clones[keep], any_cols)
        clone_tables.append(table[table.label.ne("untested")].assign(run=run))
        aurocs += auroc_rows(table, signatures, "MANA", "viral", dataset="caushi", run=run, tissue="all")

        # The same question within each tissue, from the per-sample scores.
        per_sample = pd.read_csv(out / "clone_sample_scores.csv", low_memory=False)
        per_sample["tissue"] = per_sample["sample"].str.split("_").str[1]
        per_sample["CDR3_beta"] = per_sample.CDR3ab.str.split("_").str[-1]
        per_sample = per_sample.merge(labels, on=["donor", "CDR3_beta"], how="inner")
        per_sample["size_percentile"] = np.nan
        for tissue, group in per_sample.groupby("tissue"):
            best = group.groupby(["donor", "CDR3_beta", "label"], as_index=False)[
                [f"signature_{n}_percentile" for n in signatures] + ["size_percentile"]].max()
            rows = auroc_rows(best, signatures, "MANA", "viral", dataset="caushi", run=run, tissue=tissue)
            aurocs += [r for r in rows if r["signature"] != "clone_size"]

        filters = config["cell_filters"]
        scored = table.selected_for_review
        facts[run] = {
            "arguments": extra, "samples": len(chosen), "patients": int(clones.donor.nunique()),
            "tissues": list(tissues), "seconds": int(seconds_file.read_text()) if seconds_file.exists() else None,
            "cells_after_gex_qc": filters[0]["before"], "cells_retained": filters[-1]["after"],
            "cell_filters": filters, "clones_scored": len(clones), "candidates": int(clones.selected_for_review.sum()),
            "eligible_clones": int(clones.eligible_for_review.sum()),
            "labelled_scored": counts(table), "labelled_min2_cells": counts(table, table.cell_count.ge(2)),
            "labelled_eligible": counts(table, table.eligible_for_review), "labelled_selected": counts(table, scored),
        }
    published = {str(k): int(v) for k, v in labels.label.value_counts().items()}
    labelled = pd.concat(clone_tables, ignore_index=True)
    aurocs += same_clone_rows(labelled, signatures, "MANA", "viral", "tumor_default", "tumor_cd3_10",
                              dataset="caushi", tissue="all")
    return aurocs, labelled, {"labels_published": published, "runs": facts}


# ------------------------------------------------------------- Eberhardt 2021

def top_beta_per_cell(contigs_csv):
    """Each cell's productive beta CDR3 with the most UMIs."""
    contigs = pd.read_csv(contigs_csv, keep_default_na=False)
    beta = contigs[contigs.chain.eq("TRB") & contigs.productive.eq("True") & contigs.cdr3.ne("None")]
    return (
        beta.sort_values(["umis", "cdr3"], ascending=[False, True], kind="mergesort")
        .drop_duplicates("barcode")
        .set_index("barcode")
        .cdr3
    )


def eberhardt_labels(samples):
    """Beta chains with >= 5 cells in a patient's tetramer sorts, at least as
    frequent there as in that patient's bulk sample (drops sort carry-over)."""
    rows = []
    for patient in sorted({s.name.split("_")[0] for s in samples}):
        mine = [s for s in samples if s.name.split("_")[0] == patient]
        betas = {s.name: top_beta_per_cell(s / "vdj" / "filtered_contig_annotations.csv") for s in mine}
        sorts = [n for n in betas if "PD1pos" not in n]
        bulks = [n for n in betas if "PD1pos" in n]
        cells = pd.concat([betas[n] for n in sorts]).value_counts()
        sort_frac = pd.concat([betas[n].value_counts(normalize=True) for n in sorts], axis=1).max(axis=1)
        bulk_frac = pd.concat([betas[n].value_counts(normalize=True) for n in bulks], axis=1).max(axis=1)
        keep = cells[cells.ge(5)].index
        keep = [b for b in keep if sort_frac[b] >= bulk_frac.get(b, 0)]
        rows += [{"donor": patient, "CDR3_beta": b, "label": "HPV", "tetramer_cells": int(cells[b])} for b in keep]
    return pd.DataFrame(rows)


def run_eberhardt(data_dir, work_dir):
    from tcrsift.cli import main as tcrsift_main

    samples = sorted((data_dir / "eberhardt" / "samples").iterdir())
    labels = eberhardt_labels(samples)
    bulk = [s for s in samples if "PD1pos" in s.name]
    sheet = work_dir / "eberhardt_bulk.csv"
    lines = ["sample,patient_id,tissue,vdj_dir,gex_dir"]
    lines += [f"{s.name},{s.name.split('_')[0]},{s.name.split('_')[2]},{s / 'vdj'},{s / 'gex'}" for s in bulk]
    sheet.write_text("\n".join(lines) + "\n")
    out = work_dir / "eberhardt_default"
    if not is_current(out):
        tcrsift_main(["prioritize", str(sheet), "-o", str(out), "--context", "solid-tumor"])
    config = json.loads((out / "prioritization.json").read_text())
    signatures = config["resolved_signatures"]
    sig_cols = [f"signature_{n}{s}" for n in signatures for s in ("_percentile", "_passes_cutoff")]
    clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
    clones["size_percentile"] = clones.groupby("donor").cell_count.rank(pct=True)
    clones = clones.merge(labels[["donor", "CDR3_beta", "label"]], on=["donor", "CDR3_beta"], how="left")
    clones["label"] = clones.label.fillna("untested")
    clones = clones.rename(columns={f"{c}_max": c for c in sig_cols})
    any_cols = ["eligible_for_review", "selected_for_review", *[c for c in sig_cols if c.endswith("cutoff")]]
    keep = ["donor", "CDR3_beta", "CDR3_alpha", "label", "cell_count", "size_percentile",
            "eligible_for_review", "selected_for_review", "selection_rank", "selected_signature", *sig_cols]
    table = collapse_by_beta(clones[keep], any_cols)
    # Positives against all unlabelled clones, among clones the command can
    # select (two or more cells), so singletons do not dominate the negatives.
    expanded = table[table.cell_count.ge(2)]
    aurocs = auroc_rows(expanded, signatures, "HPV", "untested", dataset="eberhardt", run="default", tissue="all")
    by_patient = {}
    for patient, group in table.groupby("donor"):
        by_patient[patient] = {"scored": counts(group), "min2_cells": counts(group, group.cell_count.ge(2)),
                               "selected": counts(group, group.selected_for_review)}
    filters = config["cell_filters"]
    facts = {
        "samples": len(bulk), "patients": int(clones.donor.nunique()),
        "tetramer_positive_betas": {str(k): int(v) for k, v in labels.donor.value_counts().items()},
        "cells_after_gex_qc": filters[0]["before"], "cells_retained": filters[-1]["after"], "cell_filters": filters,
        "clones_scored": len(clones), "candidates": int(clones.selected_for_review.sum()),
        "labelled_scored": counts(table), "labelled_min2_cells": counts(table, table.cell_count.ge(2)),
        "labelled_selected": counts(table, table.selected_for_review), "by_patient": by_patient,
    }
    return aurocs, table[table.label.ne("untested")], facts


# ------------------------------------------- shared by the full-command runs

def one_label_per_beta(table):
    """One row per patient and beta chain. A beta listed with two different
    labels (for example reactive with one alpha, not with another) is dropped."""
    distinct = table.groupby(["donor", "CDR3_beta"]).label.transform("nunique")
    return table[distinct.eq(1)].drop_duplicates(["donor", "CDR3_beta"])


def run_prioritize(sheet_rows, sheet, out, extra=()):
    """Write a sample sheet, run the unmodified command once, return its config."""
    from tcrsift.cli import main as tcrsift_main

    sheet.write_text("\n".join(["sample,patient_id,tissue,vdj_dir,gex_dir", *sheet_rows]) + "\n")
    if not is_current(out):
        tcrsift_main(["prioritize", str(sheet), "-o", str(out), "--context", "solid-tumor", *extra])
    return json.loads((out / "prioritization.json").read_text())


def labelled_clone_table(out, labels, signatures):
    """Scored clones joined to beta-chain labels, one row per labelled beta."""
    sig_cols = [f"signature_{n}{s}" for n in signatures for s in ("_percentile", "_passes_cutoff")]
    clones = pd.read_csv(out / "all_scored_clones.csv", low_memory=False)
    clones = clones.merge(labels[["donor", "CDR3_beta", "label"]], on=["donor", "CDR3_beta"], how="left")
    clones = clones.rename(columns={f"{c}_max": c for c in sig_cols}).copy()
    clones["label"] = clones.label.fillna("untested")
    clones["size_percentile"] = clones.groupby("donor").cell_count.rank(pct=True)
    from tcrsift._schema import truthy_flag_series

    cutoff_cols = [c for c in sig_cols if c.endswith("cutoff")]
    clones["passes_any"] = pd.DataFrame(
        {c: truthy_flag_series(clones[c]) for c in cutoff_cols}, index=clones.index
    ).any(axis=1)
    clones["passes_tumor_signature"] = pd.DataFrame(
        {f"signature_{n}_passes_cutoff": truthy_flag_series(clones[f"signature_{n}_passes_cutoff"]) for n in TUMOR_SIGNATURES},
        index=clones.index,
    ).any(axis=1)
    # The simple alternative to the shortlist: a patient's top clones by one score.
    rank = clones[clones.cell_count.ge(2)].groupby("donor")[
        "signature_TumorReactive_percentile"].rank(ascending=False, method="first")
    clones["top10_tumor_reactive"] = rank.le(10).reindex(clones.index, fill_value=False)
    clones["top20_tumor_reactive"] = rank.le(20).reindex(clones.index, fill_value=False)
    any_cols = ["eligible_for_review", "selected_for_review", "passes_any", "passes_tumor_signature",
                "top10_tumor_reactive", "top20_tumor_reactive", *cutoff_cols]
    keep = ["donor", "CDR3_beta", "CDR3_alpha", "label", "cell_count", "size_percentile",
            "selection_rank", "selected_signature", *any_cols,
            *[c for c in sig_cols if c.endswith("percentile")]]
    return clones, collapse_by_beta(clones[keep], any_cols)


def command_facts(config, clones, table):
    filters = config["cell_filters"]
    expanded = table.cell_count.ge(2)
    return {
        "patients": int(clones.donor.nunique()), "cells_after_gex_qc": filters[0]["before"],
        "cells_retained": filters[-1]["after"], "cell_filters": filters, "clones_scored": len(clones),
        "candidates": int(clones.selected_for_review.sum()), "labelled_scored": counts(table),
        "labelled_min2_cells": counts(table, expanded),
        "labelled_min2_pass_any": counts(table, expanded & table.passes_any),
        "labelled_min2_pass_tumor_signature": counts(table, expanded & table.passes_tumor_signature),
        "labelled_eligible": counts(table, table.eligible_for_review),
        "labelled_selected": counts(table, table.selected_for_review),
        "labelled_top10_tumor_reactive_per_patient": counts(table, table.top10_tumor_reactive),
        "labelled_top20_tumor_reactive_per_patient": counts(table, table.top20_tumor_reactive),
    }


# ------------------------------------------------------------- Chiffelle 2024

def chiffelle_labels(data_dir):
    """Table S3: one clonotype per row, tested for tumor reactivity. Chains are
    written V_CDR3_J with one extra trailing residue on the CDR3."""
    raw = pd.read_excel(data_dir / "chiffelle" / "table_s3.xlsx", header=2)
    raw["Patient"] = raw["Patient"].ffill()
    raw = raw[raw["Tumor-reactivity"].isin(["yes", "no"]) & raw.iloc[:, 1].notna()]
    table = pd.DataFrame({
        "donor": "P" + raw["Patient"].astype(int).astype(str),
        "CDR3_beta": raw.iloc[:, 1].str.split("_").str[1].str[:-1],
        "label": np.where(raw["Tumor-reactivity"].eq("yes"), "TR", "NTR"),
    })
    return one_label_per_beta(table)


def run_chiffelle(data_dir, work_dir):
    samples = sorted((data_dir / "chiffelle" / "samples").iterdir())
    labels = chiffelle_labels(data_dir)
    rows = [f"{s.name},{s.name.split('_')[0]},tumor,{s / 'vdj'},{s / 'gex'}" for s in samples]
    aurocs, tables, per_patient, facts = [], [], [], {}
    for run, extra in SENSITIVITY_RUNS:
        out = work_dir / f"chiffelle_{run}"
        config = run_prioritize(rows, work_dir / "chiffelle_tumor.csv", out, extra)
        signatures = config["resolved_signatures"]
        clones, table = labelled_clone_table(out, labels, signatures)
        aurocs += auroc_rows(table, signatures, "TR", "NTR", dataset="chiffelle", run=run, tissue="tumor")
        tables.append(table[table.label.ne("untested")].assign(run=run))
        facts[run] = {"arguments": extra, **command_facts(config, clones, table)}
        if run == "default":
            preset_signatures = signatures
            for patient, group in table.groupby("donor"):
                for row in auroc_rows(group, signatures, "TR", "NTR", donor=patient):
                    per_patient.append({k: row[k] for k in ("donor", "signature", "n_positive", "n_negative", "auroc")})
    summary = {"samples": len(samples), "labels_published": counts(labels), "runs": facts}
    labelled = pd.concat(tables, ignore_index=True)
    aurocs += same_clone_rows(labelled, preset_signatures, "TR", "NTR",
                              "default", "cd3_10", dataset="chiffelle", tissue="tumor")
    return aurocs, labelled, summary, pd.DataFrame(per_patient)


# --------------------------------------------------------------- Pulliam 2024

def pulliam_tumor_samples(data_dir, work_dir):
    """Split each hashed run by hashtag and write its tumor cells as a sample.
    A cell is assigned to its top hashtag when that has at least 10 counts and
    at least three times the runner-up."""
    import gzip

    import scipy.io

    root = data_dir / "pulliam"
    meta = pd.read_csv(root / "meta_data.csv.gz", encoding="utf-8-sig")
    made, demux = [], {}
    for run_dir in sorted((root / "runs").iterdir()):
        run = run_dir.name
        tubes = meta[meta.Run_id.eq(run)]
        features = pd.read_csv(run_dir / "features.tsv.gz", sep="\t", header=None)
        barcodes = pd.read_csv(run_dir / "barcodes.tsv.gz", header=None)[0].to_numpy()
        matrix = scipy.io.mmread(run_dir / "matrix.mtx.gz").tocsr()
        tags = features.index[features[0].str.startswith("Hashtag")].to_numpy()
        hto = matrix[tags].toarray().T
        ordered = np.sort(hto, axis=1)
        clear = (ordered[:, -1] >= 10) & (ordered[:, -1] >= 3 * np.maximum(ordered[:, -2], 1))
        number = features.loc[tags[hto.argmax(axis=1)], 0].str[len("Hashtag"):].astype(int).to_numpy()
        tissue = pd.Series(np.where(clear, number, 0)).map(dict(zip(tubes.Hashtag, tubes.tissue))).fillna("unassigned")
        demux[run] = {str(k): int(v) for k, v in tissue.value_counts().items()}
        tumor = tissue.eq("Tumor").to_numpy()
        genes = features[2].eq("Gene Expression").to_numpy()
        sample = work_dir / "pulliam_samples" / f"{run}_tumor"
        gex = sample / "gex" / "filtered_feature_bc_matrix"
        # Rebuilt on every run (a few seconds), so a changed rule or an interrupted
        # run cannot leave a stale sample behind.
        gex.mkdir(parents=True, exist_ok=True)
        (sample / "vdj").mkdir(exist_ok=True)
        with gzip.open(gex / "matrix.mtx.gz", "wb") as handle:
            scipy.io.mmwrite(handle, matrix[genes][:, tumor].tocoo(), field="integer")
        # The reference holds human and virus genomes; drop the genome prefix.
        features[genes].replace(r"^GRCh38_", "", regex=True).to_csv(
            gex / "features.tsv.gz", sep="\t", header=False, index=False)
        pd.Series(barcodes[tumor]).to_csv(gex / "barcodes.tsv.gz", header=False, index=False)
        contigs = pd.read_csv(run_dir / "all_contig_annotations.csv.gz", low_memory=False)
        called = contigs.is_cell.astype(str).str.lower().eq("true")
        confident = contigs.high_confidence.astype(str).str.lower().eq("true")
        contigs[called & confident & contigs.barcode.isin(barcodes[tumor])].to_csv(
            sample / "vdj" / "filtered_contig_annotations.csv", index=False)
        made.append((sample, tubes.Patient_id.iloc[0]))
    return made, demux


def pulliam_labels(data_dir):
    """Table S4 clones with their virus, assigned to a patient through the
    tube (specimen) they were found in. Clones from the other cohort are dropped."""
    root = data_dir / "pulliam"
    meta = pd.read_csv(root / "meta_data.csv.gz", encoding="utf-8-sig")
    plain = lambda names: names.astype(str).str.lower().str.replace(r"[^a-z0-9]", "", regex=True)  # noqa: E731
    patient = dict(zip(plain(meta.Tube), meta.Patient_id))
    raw = pd.read_excel(root / "table_s4.xlsx", header=1)
    table = pd.DataFrame({"donor": plain(raw.subject).map(patient), "CDR3_beta": raw.cdr3_b_aa,
                          "label": np.where(raw.virus.eq("MCPyV"), "MCPyV", "other_virus")})
    return one_label_per_beta(table.dropna(subset=["donor"]))


def run_pulliam(data_dir, work_dir):
    samples, demux = pulliam_tumor_samples(data_dir, work_dir)
    labels = pulliam_labels(data_dir)
    rows = [f"{s.name},{patient},tumor,{s / 'vdj'},{s / 'gex'}" for s, patient in samples]
    aurocs, tables, facts = [], [], {}
    for run, extra in SENSITIVITY_RUNS[:3]:
        out = work_dir / f"pulliam_{run}"
        config = run_prioritize(rows, work_dir / "pulliam_tumor.csv", out, extra)
        signatures = config["resolved_signatures"]
        clones, table = labelled_clone_table(out, labels, signatures)
        # As for Eberhardt: positives against unlabelled clones of two or more
        # cells, within the one tumor that has positives.
        positive_donors = table.loc[table.label.eq("MCPyV"), "donor"].unique()
        expanded = table[table.cell_count.ge(2) & table.donor.isin(positive_donors)]
        aurocs += auroc_rows(expanded, signatures, "MCPyV", "untested", dataset="pulliam", run=run, tissue="tumor")
        tables.append(table[table.label.ne("untested")].assign(run=run))
        facts[run] = {"arguments": extra, **command_facts(config, clones, table)}
    in_run = set(clones.donor)
    summary = {"samples": len(samples), "cells_by_hashtag_tissue": demux, "runs": facts,
               "labels_in_these_patients": {d: counts(g) for d, g in labels.groupby("donor") if d in in_run}}
    return aurocs, pd.concat(tables, ignore_index=True), summary


# ------------------------------------------------------------------ Meng 2023

def meng_cells(data_dir):
    """Log-normalized AnnData of the authors' CD8 cells with clone ids and labels."""
    import anndata as ad
    import scanpy as sc

    root = data_dir / "meng"
    meta = pd.read_csv(root / "GSE254250_PDAC_9sample_CD8_RNA_metadata_matrix.csv.gz").iloc[:, [0, 1, 5, 8]]
    meta.columns = ["cell", "patient", "clonotype", "label"]
    meta = meta.dropna(subset=["clonotype"])
    meta["barcode"] = meta.cell.str.split("_").str[1]
    meta["label"] = meta.label.replace("-", "untested")
    parts = []
    for patient in sorted(meta.patient.unique()):
        tag = "Pool" if patient in ("TIPC301", "TIPC413", "TIPC416") else patient[4:]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            matrix = sc.read_10x_h5(root / f"GSE254250_{tag}_GEX_filtered_feature_bc_matrix.h5")
        matrix.var_names_make_unique()
        cells = meta[meta.patient.eq(patient)].set_index("barcode")
        part = matrix[cells.index].copy()
        part.obs = pd.DataFrame({"donor": patient, "sample": patient, "lineage": "cd8",
                                 "CDR3ab": cells.clonotype.to_numpy(), "label": cells.label.to_numpy()},
                                index=cells.cell.to_numpy())
        parts.append(part)
    adata = ad.concat(parts)
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    return adata


def meng_scores(adata, signatures, background_quantile=0.99):
    """Clone-by-sample score table exactly as ``prioritize`` builds it."""
    from tcrsift.signature_background import calibrate_signature_background

    columns = [f"signature_{name}" for name in signatures]
    per_cell = adata.obs[[*STRATUM, "CDR3ab", *columns]]
    groups = per_cell.groupby([*STRATUM, "CDR3ab"], observed=True)
    scores = groups[columns].mean().join(groups.size().rename("cells")).reset_index()
    scores["frequency"] = scores.cells / scores.groupby(STRATUM).cells.transform("sum")
    for col in columns:
        scores[f"{col}_percentile"] = scores.groupby(STRATUM)[col].rank(method="average", pct=True)
    scores, _ = calibrate_signature_background(scores, per_cell, signatures, quantile=background_quantile)
    scores["size_percentile"] = scores.groupby(STRATUM).cells.rank(pct=True)
    scores["label"] = scores.CDR3ab.map(adata.obs.groupby("CDR3ab").label.first())
    return scores


def any_pass(scores, signatures, evidence="either", quantile=0.9):
    from tcrsift.prioritize import signature_pass_mask

    passed = np.zeros(len(scores), bool)
    for name in signatures:
        passed |= signature_pass_mask(scores, name, quantile=quantile, evidence=evidence).to_numpy()
    return passed


def rate_row(scores, passed, **tags):
    row = dict(tags)
    for label in ("TR", "NTR", "untested"):
        of = scores.label.eq(label).to_numpy()
        row |= {f"{label}_pass": int((passed & of).sum()), f"{label}_total": int(of.sum())}
    row["precision_tested"] = row["TR_pass"] / max(row["TR_pass"] + row["NTR_pass"], 1)
    return row


def run_meng(data_dir):
    from tcrsift.prioritize import (
        PRIORITY_CONTEXT_SIGNATURES,
        resolve_signatures,
        select_round_robin,
        signature_pass_mask,
    )
    from tcrsift.til_prioritize import _score_within_samples

    adata = meng_cells(data_dir)
    signatures = resolve_signatures("solid-tumor", "cd8")
    for name in signatures:
        adata.obs[f"signature_{name}"] = _score_within_samples(adata, name)
    scores = meng_scores(adata, signatures)

    aurocs = auroc_rows(scores, signatures, "TR", "NTR", dataset="meng", run="default", tissue="tumor")
    per_patient = []
    for patient, group in scores.groupby("donor"):
        for row in auroc_rows(group, signatures, "TR", "NTR", donor=patient):
            per_patient.append({k: row[k] for k in ("donor", "signature", "n_positive", "n_negative", "auroc")})

    # 2. Cutoffs: per signature and route, then signature sets, then thresholds.
    cutoffs = []
    for name in signatures:
        for route in ("mean", "subset", "either"):
            passed = signature_pass_mask(scores, name, evidence=route).to_numpy()
            cutoffs.append(rate_row(scores, passed, comparison="signature", setting=name, evidence=route))
    core = [n for n in signatures if n not in TUMOR_SIGNATURES]
    sets = {"solid-tumor preset (9)": signatures, "general state only (6)": core,
            "tumor-reactivity only (3)": list(TUMOR_SIGNATURES)}
    for label, names in sets.items():
        for route in ("mean", "subset", "either"):
            cutoffs.append(rate_row(scores, any_pass(scores, names, route),
                                    comparison="signature_set", setting=label, evidence=route))
    for quantile in (0.0, 0.5, 0.7, 0.8, 0.9, 0.95):
        cutoffs.append(rate_row(scores, any_pass(scores, signatures, quantile=quantile),
                                comparison="rank_gate", setting=quantile, evidence="either"))
    for background in (0.9, 0.95, 0.99, 0.995):
        swept = scores if background == 0.99 else meng_scores(adata, signatures, background)
        cutoffs.append(rate_row(swept, any_pass(swept, signatures),
                                comparison="background_quantile", setting=background, evidence="either"))

    # 3. Selection: the real round-robin against two simple alternatives.
    clones = scores[["donor", "CDR3ab", "label", "cells", "frequency"]].copy()
    clones["eligible_for_review"] = (clones.cells.ge(2) & clones.frequency.ge(0.001)
                                     & any_pass(scores, signatures))
    abundant = clones[clones.cells.ge(2) & clones.frequency.ge(0.001)]
    selection = []
    for budget in (*BUDGETS, 0):
        picked = select_round_robin(
            clones,
            scores,
            signatures,
            max_clones=budget,
            evidence="either",
            priority_signatures=PRIORITY_CONTEXT_SIGNATURES["solid-tumor"],
        )
        picked = picked[picked.selected_for_review]
        size = len(picked)
        picked_flat = select_round_robin(
            clones, scores, signatures, max_clones=size, evidence="either"
        )
        narrow = clones.assign(eligible_for_review=clones.cells.ge(2) & clones.frequency.ge(0.001)
                               & any_pass(scores, TUMOR_SIGNATURES))
        narrow = select_round_robin(narrow, scores, list(TUMOR_SIGNATURES), max_clones=size, evidence="either")
        methods = {
            "TCRsift round-robin": picked.label,
            "TCRsift round-robin (flat)": picked_flat[picked_flat.selected_for_review].label,
            "TCRsift round-robin, tumor-reactivity signatures only": narrow[narrow.selected_for_review].label,
            "largest clones": abundant.sort_values(
                ["frequency", "cells", "CDR3ab"], ascending=[False, False, True], kind="mergesort"
            ).head(size).label,
            "TumorReactive score only": abundant.assign(
                s=scores.loc[abundant.index, "signature_TumorReactive_percentile"]
            ).sort_values(["s", "CDR3ab"], ascending=[False, True], kind="mergesort").head(size).label,
        }
        rng = np.random.default_rng(0)
        draws = [abundant.label.iloc[rng.choice(len(abundant), min(size, len(abundant)), replace=False)]
                 for _ in range(500)]
        for method, chosen in methods.items():
            tr, ntr = int(chosen.eq("TR").sum()), int(chosen.eq("NTR").sum())
            selection.append({"budget": budget or "all qualifying", "method": method, "selected": len(chosen),
                              "TR": tr, "NTR": ntr, "untested": len(chosen) - tr - ntr})
        selection.append({"budget": budget or "all qualifying", "method": "random abundant clones (mean of 500)",
                          "selected": size, "TR": float(np.mean([d.eq("TR").sum() for d in draws])),
                          "NTR": float(np.mean([d.eq("NTR").sum() for d in draws])),
                          "untested": float(np.mean([d.eq("untested").sum() for d in draws]))})

    default_pass = any_pass(scores, signatures)
    facts = {
        "patients": int(scores.donor.nunique()), "cells": int(adata.n_obs), "clones_scored": len(scores),
        "labelled_clones": counts(scores), "passes_any_default_cutoff": counts(scores, default_pass),
        "abundant_clones": counts(scores, scores.index.isin(abundant.index)),
        "eligible_clones": counts(clones, clones.eligible_for_review),
    }
    keep = ["donor", "CDR3ab", "label", "cells", "frequency", "size_percentile",
            *[f"signature_{n}_percentile" for n in signatures]]
    labelled = scores.loc[scores.label.ne("untested"), keep].rename(columns={"CDR3ab": "clonotype"})
    labelled["passes_any"] = default_pass[scores.label.ne("untested").to_numpy()]
    return aurocs, labelled, facts, pd.DataFrame(per_patient), pd.DataFrame(cutoffs), pd.DataFrame(selection)


GROUND_TRUTH_PROVENANCE = {
    "meng": {
        "citation": "Meng et al., Science Translational Medicine 2023 (15:eadh9562; DOI:10.1126/scitranslmed.adh9562; PMID:37992152)",
        "cancer_type": "Pancreatic ductal adenocarcinoma (PDAC)",
        "gex_vdj_accession": "GEO GSE254250 (7 Cell Ranger filtered_feature_bc_matrix.h5 files across 9 patients)",
        "label_source": "GEO GSE254250_PDAC_9sample_CD8_RNA_metadata_matrix.csv.gz (columns 0, 1, 5, 8: 'Unnamed: 0', 'orig.ident', 'clonotype_id', 'Reactivity (n=149)')",
        "experimental_assay": "Cloned TCRs retrovirally/lentivirally expressed in reporter/autologous T cells and co-cultured directly against each patient's autologous PDAC organoids/tumor lines.",
        "positive_definition": "TR: clonotypes experimentally reactive against autologous tumor cells (94 clonotypes in deposited metadata).",
        "negative_definition": "NTR: clonotypes from the same tumors experimentally tested against autologous tumor cells with no response (55 clonotypes in deposited metadata).",
        "chain_matching": "Exact per-cell 'clonotype_id' provided by authors (raw VDJ contig tables not deposited on GEO).",
        "circular_signatures": [],
    },
    "chiffelle": {
        "citation": "Chiffelle et al., Immunity 2024 (57:2466-2482; DOI:10.1016/j.immuni.2024.08.014; PMID:39321802); scRNA/scTCR deposited with Barras et al., Science Immunology 2024 (DOI:10.1126/sciimmunol.adg7995)",
        "cancer_type": "Metastatic melanoma (baseline pre-ACT tumor biopsies T0/T01/T02)",
        "gex_vdj_accession": "GEO GSE222448 (13 baseline tumor CD45+ count matrices and unfiltered contig CSVs across 11 tested patients: Patients 1-5, 7-11, 13 from GSM6885352-GSM6885368 and GSM6930579-GSM6930593)",
        "label_source": "Immunity Supplementary Table S3 (1-s2.0-S1074761324004138-mmc3.xlsx, sheet 'clonotypes', header row 3, columns 'Patient', 'TCR β chain' in V_CDR3_J format, 'Tumor-reactivity' in {'yes', 'no'}).",
        "experimental_assay": "TCRs cloned from baseline tumors or TIL infusion products and screened against autologous melanoma cell lines (CD137/4-1BB upregulation and/or IFN-gamma secretion).",
        "positive_definition": "TR: clonotypes with 'Tumor-reactivity == yes' in Table S3 (123 published; 72 scored in baseline tumor after default cell filters).",
        "negative_definition": "NTR: clonotypes experimentally tested against autologous tumor with 'Tumor-reactivity == no' (124 published; 94 scored in baseline tumor). Ambiguous beta chains carrying conflicting labels across alpha partners are dropped.",
        "chain_matching": "Within-patient exact CDR3_beta amino-acid match (stripping the extra trailing J-anchor residue encoded in Table S3's V_CDR3_J notation); collapsed to one row per (donor, CDR3_beta).",
        "circular_signatures": [],
    },
    "caushi": {
        "citation": "Caushi et al., Nature 2021 (596:126-132; DOI:10.1038/s41586-021-03752-4; PMID:34290408)",
        "cancer_type": "Resectable non-small-cell lung cancer (NSCLC) post-neoadjuvant anti-PD-1 (nivolumab)",
        "gex_vdj_accession": "GEO GSE176021 (57 paired Cell Ranger GEX + filtered VDJ samples across 8 antigen-annotated patients: 33 tumor/metastasis, 24 lymph node/adjacent normal lung)",
        "label_source": "Nature Supplementary Table 8 (41586_2021_3752_MOESM8_ESM.xlsx, sheet 'Supplementary Table 8', header row 3, columns 'Patient', 'TCR clonotype (AA)', 'Type of antigen recognized (viral or MANA)').",
        "experimental_assay": "MANAFEST (Mutation-Associated Neoantigen Functional Expansion of Specific T cells): 10-day peptide-stimulated PBMC/TIL TCRB sequencing identifying patient-specific neoantigen (MANA) vs CEF viral (Influenza, EBV, CMV) peptide expansion.",
        "positive_definition": "MANA: CDR3_beta clonotypes expanding specifically to predicted patient neoantigen peptides (74 published; 18 scored in tumor).",
        "negative_definition": "viral: CDR3_beta clonotypes from the same patients expanding to Influenza, EBV, or CMV epitopes (82 published; 13 scored in tumor). Note: negatives are verified viral bystanders, not tumor-non-reactive random TILs.",
        "chain_matching": "Within-patient exact CDR3_beta match collapsed across alpha partners.",
        "circular_signatures": ["MANAscore"],
    },
    "eberhardt": {
        "citation": "Eberhardt et al., Nature 2021 (597:279-284; DOI:10.1038/s41586-021-03862-z; PMID:34471285)",
        "cancer_type": "HPV-positive head and neck squamous cell carcinoma (HNSCC)",
        "gex_vdj_accession": "GEO GSE180268 (6 bulk PD-1+ CD8+ libraries [3 primary tumor, 3 metastatic lymph node] paired with HLA-A*01:01 HPV16 E2 tetramer-sorted libraries across 3 patients: HPV7, HPV34, HPV51)",
        "label_source": "Derived directly from paired GSE180268 HLA-A*01:01 HPV16 E2 (QVDYYGLYY and KSAIVTLTY) tetramer-sorted vs bulk PD-1+ VDJ contig tables.",
        "experimental_assay": "FACS isolation with HLA-A*01:01 HPV16 E2 peptide-MHC class I tetramers vs paired PD-1+ CD8+ bulk tumor/metLN libraries.",
        "positive_definition": "HPV: productive CDR3_beta chains observed in >=5 cells in a patient's HPV16 E2 tetramer-sorted libraries and at least as frequent in tetramer sort as in paired bulk PD-1+ libraries (guards against sort carryover).",
        "negative_definition": "untested: all other expanded (>=2 cells) paired-TCR CD8 clonotypes in bulk PD-1+ libraries (no tested tumor-non-reactive negatives exist; bulk samples were pre-gated on PD-1 and depleted of tetramer+ cells).",
        "chain_matching": "Within-patient exact CDR3_beta match collapsed across alpha partners.",
        "circular_signatures": [],
    },
    "pulliam": {
        "citation": "Pulliam et al., Cell Reports Medicine 2024 (5:101412; DOI:10.1016/j.xcrm.2024.101412; PMID:38340723)",
        "cancer_type": "Merkel cell carcinoma (MCPyV-positive virus-driven skin neuroendocrine carcinoma)",
        "gex_vdj_accession": "GEO GSE227708 (hashed runs Z1249 and Z1308 demultiplexed to Tumor cells via GSE227708_meta_data.csv.gz: top Hashtag >=10 UMI and >=3x runner-up)",
        "label_source": "Cell Reports Medicine Table S4 (1-s2.0-S2666379124000351-mmc5.xlsx, sheet 'Table 4', header row 2, columns 'subject', 'cdr3_b_aa', 'virus').",
        "experimental_assay": "DNA-barcoded and fluorescently labeled MCPyV T-antigen and control viral (CMV/EBV/Flu) peptide-MHC class I multimers.",
        "positive_definition": "MCPyV: CDR3_beta clonotypes annotated with virus == 'MCPyV' in Table S4 (22 in Z1249 patient; 5 with >=2 cells in demultiplexed tumor).",
        "negative_definition": "untested: remaining expanded (>=2 cells) paired-TCR CD8 clonotypes in the same demultiplexed tumor (53 clones; Z1249 has no non-MCPyV viral controls in tumor).",
        "chain_matching": "Within-patient exact CDR3_beta match mapped via specimen Tube -> Patient_id and collapsed across alpha partners.",
        "circular_signatures": [],
    },
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path, help="tcrsift output directories")
    parser.add_argument("--out-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    from tcrsift import __version__
    from tcrsift.datacache import cached_path

    vdjdb_release = "uncached"
    vdjdb_entries = 0
    vdjdb_paired_ab_entries = 0
    vdj_p = cached_path("vdjdb")
    if vdj_p is not None:
        meta_candidates = [
            (vdj_p if vdj_p.is_dir() else vdj_p.parent) / "_meta.json",
            vdj_p.parent / "_meta.json",
        ]
        for meta_json in meta_candidates:
            if meta_json.exists():
                meta_obj = json.loads(meta_json.read_text())
                url = str(meta_obj.get("download_url") or "")
                if "/releases/download/" in url:
                    vdjdb_release = "vdjdb-" + url.split("/releases/download/")[1].split("/")[0]
                break
        if vdj_p.is_dir():
            txt_hits = sorted(vdj_p.rglob("vdjdb.txt"))
            txt_p = txt_hits[0] if txt_hits else vdj_p / "vdjdb.txt"
            if vdjdb_release == "uncached" and txt_hits:
                vdjdb_release = txt_hits[0].parent.name
        else:
            txt_p = vdj_p
        if txt_p.exists():
            vdf = pd.read_csv(txt_p, sep="\t", usecols=["gene", "complex.id"], low_memory=False)
            trb = vdf[vdf["gene"] == "TRB"]
            vdjdb_entries = int(len(trb))
            tra_cids = set(vdf.loc[(vdf["gene"] == "TRA") & (vdf["complex.id"] > 0), "complex.id"])
            vdjdb_paired_ab_entries = int(trb["complex.id"].isin(tra_cids).sum())

    meng_auroc, meng_table, meng_facts, meng_patients, meng_cutoffs, meng_selection = run_meng(args.data_dir)
    caushi_auroc, caushi_table, caushi_facts = run_caushi(args.data_dir, args.work_dir)
    hpv_auroc, hpv_table, hpv_facts = run_eberhardt(args.data_dir, args.work_dir)
    mel_auroc, mel_table, mel_facts, mel_patients = run_chiffelle(args.data_dir, args.work_dir)
    mcc_auroc, mcc_table, mcc_facts = run_pulliam(args.data_dir, args.work_dir)
    out = args.out_dir
    mel_table.to_csv(out / "chiffelle_labelled_clones.csv", index=False)
    mel_patients.to_csv(out / "chiffelle_per_patient_auroc.csv", index=False)
    mcc_table.to_csv(out / "pulliam_labelled_clones.csv", index=False)
    pd.DataFrame(caushi_auroc + meng_auroc + hpv_auroc + mel_auroc + mcc_auroc).to_csv(out / "signature_auroc.csv", index=False)
    caushi_table.to_csv(out / "caushi_labelled_clones.csv", index=False)
    hpv_table.to_csv(out / "eberhardt_labelled_clones.csv", index=False)
    meng_table.to_csv(out / "meng_labelled_clones.csv", index=False)
    meng_patients.to_csv(out / "meng_per_patient_auroc.csv", index=False)
    meng_cutoffs.to_csv(out / "meng_cutoffs.csv", index=False)
    meng_selection.to_csv(out / "meng_selection.csv", index=False)
    (out / "summary.json").write_text(json.dumps(
        {
            "tcrsift_version": __version__,
            "vdjdb_release": vdjdb_release,
            "vdjdb_entries": vdjdb_entries,
            "vdjdb_paired_ab_entries": vdjdb_paired_ab_entries,
            "ground_truth_provenance": GROUND_TRUTH_PROVENANCE,
            "caushi": caushi_facts,
            "meng": meng_facts,
            "eberhardt": hpv_facts,
            "chiffelle": mel_facts,
            "pulliam": mcc_facts,
        },
        indent=2,
    ) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
