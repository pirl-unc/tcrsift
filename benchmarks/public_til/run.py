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
    ("tumor_cd3_3", TUMOR_TISSUES, ["--min-cd3", "3"]),
    ("tumor_cd3_0", TUMOR_TISSUES, ["--min-cd3", "0"]),
    # Is the CD3 floor special, or is any depth filter as good? No CD3 floor,
    # but require 4,000 UMIs per cell instead of 500.
    ("tumor_depth_gate", TUMOR_TISSUES, ["--min-cd3", "0", "--min-counts", "4000"]),
    ("other_tissue_default", ("LN", "normal"), []),
]
TUMOR_SIGNATURES = ("TumorReactive", "MANAscore", "NeoTCR8")
# The default command, two lower CD3 floors, and the tumor-reactivity
# signatures alone, for the datasets small enough to rerun cheaply.
SENSITIVITY_RUNS = [
    ("default", []),
    ("cd3_3", ["--min-cd3", "3"]),
    ("cd3_0", ["--min-cd3", "0"]),
    ("tumor_signatures", ["--signatures", *TUMOR_SIGNATURES]),
]
BUDGETS = (25, 50, 100, 200)


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
    labelled = table[table.label.ne("untested")].sort_values("cell_count", ascending=False)
    grouped = labelled.groupby(["donor", "CDR3_beta"], sort=False)
    collapsed = grouped.head(1).drop(columns=any_cols).merge(
        grouped[any_cols].any().reset_index(), on=["donor", "CDR3_beta"])
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
        if not (out / "prioritization.json").exists():
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
    return aurocs, pd.concat(clone_tables, ignore_index=True), {"labels_published": published, "runs": facts}


# ------------------------------------------------------------- Eberhardt 2021

def top_beta_per_cell(contigs_csv):
    """Each cell's productive beta CDR3 with the most UMIs."""
    contigs = pd.read_csv(contigs_csv, keep_default_na=False)
    beta = contigs[contigs.chain.eq("TRB") & contigs.productive.eq("True") & contigs.cdr3.ne("None")]
    return beta.sort_values("umis", ascending=False).drop_duplicates("barcode").set_index("barcode").cdr3


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
    if not (out / "prioritization.json").exists():
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

def run_prioritize(sheet_rows, sheet, out, extra=()):
    """Write a sample sheet, run the unmodified command once, return its config."""
    from tcrsift.cli import main as tcrsift_main

    sheet.write_text("\n".join(["sample,patient_id,tissue,vdj_dir,gex_dir", *sheet_rows]) + "\n")
    if not (out / "prioritization.json").exists():
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
    cutoff_cols = [c for c in sig_cols if c.endswith("cutoff")]
    clones["passes_any"] = clones[cutoff_cols].fillna(False).astype(bool).any(axis=1)
    clones["passes_tumor_signature"] = clones[
        [f"signature_{n}_passes_cutoff" for n in TUMOR_SIGNATURES]].fillna(False).astype(bool).any(axis=1)
    clones["tumor_reactive_rank"] = clones[clones.cell_count.ge(2)].groupby("donor")[
        "signature_TumorReactive_percentile"].rank(ascending=False, method="first")
    any_cols = ["eligible_for_review", "selected_for_review", "passes_any", "passes_tumor_signature", *cutoff_cols]
    keep = ["donor", "CDR3_beta", "CDR3_alpha", "label", "cell_count", "size_percentile",
            "selection_rank", "selected_signature", "tumor_reactive_rank", *any_cols,
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
        # The simple alternative: each patient's top clones by one score.
        "labelled_top10_tumor_reactive_per_patient": counts(table, table.tumor_reactive_rank.le(10)),
        "labelled_top20_tumor_reactive_per_patient": counts(table, table.tumor_reactive_rank.le(20)),
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
    return table.drop_duplicates(["donor", "CDR3_beta"])


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
            for patient, group in table.groupby("donor"):
                for row in auroc_rows(group, signatures, "TR", "NTR", donor=patient):
                    per_patient.append({k: row[k] for k in ("donor", "signature", "n_positive", "n_negative", "auroc")})
    summary = {"samples": len(samples), "labels_published": counts(labels), "runs": facts}
    return aurocs, pd.concat(tables, ignore_index=True), summary, pd.DataFrame(per_patient)


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
        if not (gex / "matrix.mtx.gz").exists():
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
    return table.dropna(subset=["donor"]).drop_duplicates(["donor", "CDR3_beta"])


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
    from tcrsift.prioritize import resolve_signatures, select_round_robin, signature_pass_mask
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
        picked = select_round_robin(clones, scores, signatures, max_clones=budget, evidence="either")
        picked = picked[picked.selected_for_review]
        size = len(picked)
        narrow = clones.assign(eligible_for_review=clones.cells.ge(2) & clones.frequency.ge(0.001)
                               & any_pass(scores, TUMOR_SIGNATURES))
        narrow = select_round_robin(narrow, scores, list(TUMOR_SIGNATURES), max_clones=size, evidence="either")
        methods = {
            "TCRsift round-robin": picked.label,
            "TCRsift round-robin, tumor-reactivity signatures only": narrow[narrow.selected_for_review].label,
            "largest clones": abundant.sort_values(["frequency", "cells"], ascending=False).head(size).label,
            "TumorReactive score only": abundant.assign(s=scores.loc[abundant.index, "signature_TumorReactive_percentile"])
                                                 .sort_values("s", ascending=False).head(size).label,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path, help="tcrsift output directories")
    parser.add_argument("--out-dir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    from tcrsift import __version__

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
        {"tcrsift_version": __version__, "caushi": caushi_facts, "meng": meng_facts,
         "eberhardt": hpv_facts, "chiffelle": mel_facts, "pulliam": mcc_facts}, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
