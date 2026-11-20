#!/usr/bin/env python3
"""Draw the benchmark figures from ``benchmarks/public_til/results``.

    python docs/figures/make_benchmark_figures.py

Unlike ``make_figures.py`` these are measurements: every number comes from the
result tables written by ``benchmarks/public_til/run.py`` on public data.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from make_figures import BLUE, INK, INK2, NEUTRAL, ORANGE, SURFACE, save

RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "public_til" / "results"
LABELS = {"clone_size": "Clone size (baseline)", "MANAscore": "MANAscore proxy"}


def auroc_forest(ax, rows, title, note=None):
    rows = rows.sort_values("auroc")
    y = range(len(rows))
    for i, row in zip(y, rows.itertuples()):
        baseline = row.signature == "clone_size"
        color = NEUTRAL if baseline else BLUE
        ax.plot([row.auroc_low, row.auroc_high], [i, i], color=color, lw=2, solid_capstyle="round")
        ax.scatter([row.auroc], [i], s=55, color=color, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text(1.03, i, f"{row.auroc:.2f}", va="center", fontsize=9, color=INK)
    ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_yticks(list(y), [LABELS.get(s, s) for s in rows.signature])
    ax.set_xlim(0.1, 1.0)
    ax.set_xlabel("AUROC (95% bootstrap interval)")
    ax.set_title(title)
    ax.grid(axis="x", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    if note:
        ax.text(0.49, len(rows) - 0.55, "chance", ha="right", va="center", fontsize=8.5, color=INK2)
        ax.set_ylim(-0.6, len(rows) - 0.3)


def signature_auroc():
    table = pd.read_csv(RESULTS / "signature_auroc.csv")
    meng = table[table.dataset.eq("meng")]
    caushi = table[table.dataset.eq("caushi") & table.run.eq("default")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    m, c = meng.iloc[0], caushi.iloc[0]
    auroc_forest(ax1, meng, f"PDAC (Meng 2023): {m.n_positive} reactive vs "
                            f"{m.n_negative} non-reactive clones", note=True)
    auroc_forest(ax2, caushi, f"NSCLC (Caushi 2021): {c.n_positive} neoantigen vs "
                              f"{c.n_negative} viral clones", note=True)
    save(fig, "benchmark_auroc")


def cutoff_pass_rates():
    summary = json.loads((RESULTS / "summary.json").read_text())["meng"]
    passed, total = summary["passes_any_default_cutoff"], summary["labelled_clones"]
    groups = [("TR", "Tumor-reactive", BLUE), ("NTR", "Tested, non-reactive", ORANGE),
              ("untested", "Untested", NEUTRAL)]
    table = pd.read_csv(RESULTS / "signature_auroc.csv")
    meng = table[table.dataset.eq("meng") & table.signature.ne("clone_size")].copy()
    meng["tr"] = meng.positive_pass / meng.n_positive
    meng["ntr"] = meng.negative_pass / meng.n_negative
    meng = meng.sort_values("tr")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1, 1.5]})
    for i, (key, label, color) in enumerate(groups):
        rate = passed.get(key, 0) / total[key]
        ax1.bar(i, rate, width=0.6, color=color)
        ax1.text(i, rate + 0.02, f"{rate:.0%}\n{passed.get(key, 0):,} of {total[key]:,}",
                 ha="center", va="bottom", fontsize=9, color=INK)
    ax1.set_xticks(range(3), [g[1] for g in groups])
    ax1.set_ylim(0, 1)
    ax1.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax1.set_ylabel("clones passing any signature cutoff")
    ax1.set_title("Any of nine signatures")
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)

    height = 0.36
    for i, row in enumerate(meng.itertuples()):
        ax2.barh(i + height / 2 + 0.02, row.tr, height=height, color=BLUE)
        ax2.barh(i - height / 2 - 0.02, row.ntr, height=height, color=ORANGE)
        ax2.text(row.tr + 0.01, i + height / 2 + 0.02, f"{row.tr:.0%}", va="center", fontsize=8.5, color=INK)
        ax2.text(row.ntr + 0.01, i - height / 2 - 0.02, f"{row.ntr:.0%}", va="center", fontsize=8.5, color=INK)
    ax2.set_yticks(range(len(meng)), [LABELS.get(s, s) for s in meng.signature])
    ax2.set_xlim(0, 0.75)
    ax2.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax2.set_xlabel("clones passing that signature's cutoff")
    ax2.set_title("By signature")
    ax2.tick_params(axis="y", length=0)
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="x", lw=0.6)
    ax2.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    ax2.legend(handles, ["Tumor-reactive", "Tested, non-reactive"], frameon=False, fontsize=9,
               loc="lower right")
    save(fig, "benchmark_cutoffs")


def cell_attrition():
    summary = json.loads((RESULTS / "summary.json").read_text())["caushi"]["default"]
    steps = [("Cells after expression QC", summary["cells_after_gex_qc"])]
    names = {"lineage": "Classified CD8", "cd3": "≥10 CD3 counts", "alpha_cdr3": "Alpha CDR3 found",
             "beta_cdr3": "Beta CDR3 found", "TRA_1_umis": "≥2 alpha UMIs", "TRB_1_umis": "≥2 beta UMIs"}
    steps += [(names[f["filter"]], f["after"]) for f in summary["cell_filters"]]
    fig, ax = plt.subplots(figsize=(8.2, 3.6), constrained_layout=True)
    total = steps[0][1]
    for i, (label, count) in enumerate(steps):
        y = len(steps) - 1 - i
        ax.barh(y, count, height=0.62, color=BLUE)
        ax.text(count + total * 0.012, y, f"{count:,}  ({count / total:.0%})", va="center", fontsize=9, color=INK)
    ax.set_yticks(range(len(steps)), [s[0] for s in reversed(steps)])
    ax.set_xlim(0, total * 1.22)
    ax.set_xticks([])
    ax.tick_params(axis="y", length=0)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title("NSCLC tumor samples: cells remaining after each default gate")
    save(fig, "benchmark_attrition")


if __name__ == "__main__":
    signature_auroc()
    cutoff_pass_rates()
    cell_attrition()
    print("wrote benchmark figures")
