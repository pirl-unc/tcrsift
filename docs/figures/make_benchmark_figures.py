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
import matplotlib.patches
import matplotlib.pyplot as plt
import pandas as pd
from make_figures import BLUE, INK, INK2, NEUTRAL, ORANGE, SURFACE, save

RESULTS = Path(__file__).resolve().parents[2] / "benchmarks" / "public_til" / "results"
LABELS = {"clone_size": "Clone size (baseline)", "MANAscore": "MANAscore proxy"}


def auroc_forest(ax, rows, title, note=None, circular=()):
    """Dot and interval per signature. ``circular`` names signatures that were
    derived from the dataset being plotted: drawn hollow, in orange, and
    labeled, because their score there is not an independent test."""
    rows = rows.sort_values("auroc")
    y = range(len(rows))
    for i, row in zip(y, rows.itertuples()):
        baseline, tainted = row.signature == "clone_size", row.signature in circular
        color = NEUTRAL if baseline else ORANGE if tainted else BLUE
        ax.plot([row.auroc_low, row.auroc_high], [i, i], color=color, lw=2, solid_capstyle="round",
                ls=(0, (2, 2)) if tainted else "-")
        ax.scatter([row.auroc], [i], s=55, color=SURFACE if tainted else color,
                   edgecolor=color if tainted else SURFACE, linewidth=1.8 if tainted else 1.5, zorder=3)
        ax.text(1.03, i, f"{row.auroc:.2f}" + ("  circular" if tainted else ""), va="center", fontsize=9,
                color=INK)
    ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_yticks(list(y), [LABELS.get(s, s) + (" *" if s in circular else "") for s in rows.signature])
    if circular:
        ax.text(0.0, -0.2, "* derived from this dataset, so not an independent test", transform=ax.transAxes,
                fontsize=8.5, color=INK2, va="top")
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
    caushi = table[table.dataset.eq("caushi") & table.run.eq("tumor_default") & table.tissue.eq("all")]
    melanoma = table[table.dataset.eq("chiffelle") & table.run.eq("default")]
    fig1, ax1 = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    fig2, ax2 = plt.subplots(figsize=(7.2, 4.7), constrained_layout=True)
    fig3, ax3 = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    m, c, k = meng.iloc[0], caushi.iloc[0], melanoma.iloc[0]
    auroc_forest(ax1, meng, f"PDAC (Meng 2023): {m.n_positive} reactive vs "
                            f"{m.n_negative} non-reactive clones", note=True)
    auroc_forest(ax2, caushi, f"NSCLC (Caushi 2021): {c.n_positive} neoantigen vs "
                              f"{c.n_negative} viral clones", note=True, circular=("MANAscore",))
    auroc_forest(ax3, melanoma, f"Melanoma (Chiffelle 2024): {k.n_positive} reactive vs "
                                f"{k.n_negative} non-reactive clones", note=True)
    save(fig1, "benchmark_auroc_pancreatic")
    save(fig2, "benchmark_auroc_lung")
    save(fig3, "benchmark_auroc_melanoma")


def auroc_by_dataset():
    """One independent signature across all five datasets. Datasets without
    tested negatives (positives against unlabeled clones) are drawn hollow."""
    table = pd.read_csv(RESULTS / "signature_auroc.csv")
    table = table[table.signature.eq("TumorReactive")]
    # (dataset, run, tissue, label, has tested negatives)
    wanted = [("meng", "default", "tumor", "Pancreatic (Meng 2023)", True),
              ("chiffelle", "default", "tumor", "Melanoma (Chiffelle 2024)", True),
              ("caushi", "tumor_default", "all", "Lung (Caushi 2021)", True),
              ("eberhardt", "default", "all", "Head and neck, HPV+ (Eberhardt 2021)", False),
              ("pulliam", "default", "tumor", "Merkel cell (Pulliam 2024)", False)]
    fig, ax = plt.subplots(figsize=(8.2, 3.6), constrained_layout=True)
    for i, (dataset, run, tissue, label, tested) in enumerate(reversed(wanted)):
        row = table[table.dataset.eq(dataset) & table.run.eq(run) & table.tissue.eq(tissue)].iloc[0]
        color = BLUE if tested else INK2
        ax.plot([row.auroc_low, row.auroc_high], [i, i], color=color, lw=2, solid_capstyle="round")
        ax.scatter([row.auroc], [i], s=55, color=color if tested else SURFACE, edgecolor=SURFACE if tested else color,
                   linewidth=1.5 if tested else 1.8, zorder=3)
        versus = "tested negatives" if tested else "unlabeled"
        ax.text(1.03, i, f"{row.auroc:.2f}  ({row.n_positive} vs {row.n_negative} {versus})", va="center",
                fontsize=9, color=INK)
    ax.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.set_yticks(range(len(wanted)), [w[3] for w in reversed(wanted)])
    ax.set_xlim(0.1, 1.0)
    ax.set_ylim(-0.6, len(wanted) - 0.4)
    ax.set_xlabel("TumorReactive AUROC (95% bootstrap interval)")
    ax.set_title("TumorReactive across the five datasets")
    ax.grid(axis="x", lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    save(fig, "benchmark_auroc_by_dataset")


def percent_axis(ax, axis="x"):
    getattr(ax, f"{axis}axis").set_major_formatter(matplotlib.ticker.PercentFormatter(1))


def cutoff_pass_rates():
    cutoffs = pd.read_csv(RESULTS / "meng_cutoffs.csv")
    either = cutoffs[cutoffs.evidence.eq("either")]
    preset = either[either.setting.eq("solid-tumor preset (9)")].iloc[0]
    groups = [("TR", "Tumor-reactive", BLUE), ("NTR", "Tested, non-reactive", ORANGE),
              ("untested", "Untested", NEUTRAL)]
    by_signature = either[either.comparison.eq("signature")].copy()
    by_signature["tr"] = by_signature.TR_pass / by_signature.TR_total
    by_signature["ntr"] = by_signature.NTR_pass / by_signature.NTR_total
    by_signature = by_signature.sort_values("tr")

    fig1, ax1 = plt.subplots(figsize=(5.4, 4.2), constrained_layout=True)
    fig2, ax2 = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    for i, (key, label, color) in enumerate(groups):
        passed, total = preset[f"{key}_pass"], preset[f"{key}_total"]
        ax1.bar(i, passed / total, width=0.6, color=color)
        ax1.text(i, passed / total + 0.02, f"{passed / total:.0%}\n{passed:,} of {total:,}",
                 ha="center", va="bottom", fontsize=9, color=INK)
    ax1.set_xticks(range(3), [g[1] for g in groups])
    ax1.set_ylim(0, 1)
    percent_axis(ax1, "y")
    ax1.set_ylabel("clones passing any signature cutoff")
    ax1.set_title("Any of nine signatures")
    ax1.grid(axis="y", lw=0.6)
    ax1.set_axisbelow(True)

    height = 0.36
    for i, row in enumerate(by_signature.itertuples()):
        ax2.barh(i + height / 2 + 0.02, row.tr, height=height, color=BLUE)
        ax2.barh(i - height / 2 - 0.02, row.ntr, height=height, color=ORANGE)
        ax2.text(row.tr + 0.01, i + height / 2 + 0.02, f"{row.tr:.0%}", va="center", fontsize=8.5, color=INK)
        ax2.text(row.ntr + 0.01, i - height / 2 - 0.02, f"{row.ntr:.0%}", va="center", fontsize=8.5, color=INK)
    ax2.set_yticks(range(len(by_signature)), [LABELS.get(s, s) for s in by_signature.setting])
    ax2.set_xlim(0, 0.75)
    percent_axis(ax2)
    ax2.set_xlabel("clones passing that signature's cutoff")
    ax2.set_title("By signature")
    ax2.tick_params(axis="y", length=0)
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="x", lw=0.6)
    ax2.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    ax2.legend(handles, ["Tumor-reactive", "Tested, non-reactive"], frameon=False, fontsize=9,
               loc="lower right")
    save(fig1, "benchmark_cutoffs_any")
    save(fig2, "benchmark_cutoffs_by_signature")


def selection_by_budget():
    """Tested clones on the list as the budget grows, for each way of choosing."""
    table = pd.read_csv(RESULTS / "meng_selection.csv")
    table = table[table.budget.astype(str).isin(["25", "50", "100", "200"])].assign(
        budget=lambda d: d.budget.astype(int))
    methods = [("TCRsift round-robin", "Round-robin, tiered solid-tumor preset", BLUE, "-"),
               ("TCRsift round-robin (flat)", "Round-robin, flat 9-signature preset", BLUE, ":"),
               ("TCRsift round-robin, tumor-reactivity signatures only", "Round-robin, 3 tumor signatures", BLUE, "--"),
               ("TumorReactive score only", "Top TumorReactive score", ORANGE, "-"),
               ("largest clones", "Largest clones", INK2, "-"),
               ("random abundant clones (mean of 500)", "Random clones", NEUTRAL, "-")]
    panels = [("TR", "Tumor-reactive clones chosen", "benchmark_selection_reactive"),
              ("NTR", "Tested non-reactive clones chosen", "benchmark_selection_nonreactive")]
    for column, title, name in panels:
        fig, ax = plt.subplots(figsize=(7.4, 4.2), constrained_layout=True)
        for method, label, color, style in methods:
            rows = table[table.method.eq(method)].sort_values("budget")
            ax.plot(rows.budget, rows[column], color=color, ls=style, lw=2, marker="o", ms=5, label=label)
        ax.set_xticks([25, 50, 100, 200])
        ax.set_ylim(0, table[["TR", "NTR"]].to_numpy().max() * 1.08)
        ax.set_xlabel("clones chosen (budget)")
        ax.set_title(title)
        ax.grid(axis="y", lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        save(fig, name)


def threshold_sweeps():
    """Reactive against non-reactive pass rate as each threshold moves."""
    cutoffs = pd.read_csv(RESULTS / "meng_cutoffs.csv")
    panels = [("rank_gate", "Rank gate (--signature-quantile)", 0.9),
              ("background_quantile", "Background floor (--signature-background-quantile)", 0.99)]
    names = ["benchmark_threshold_rank", "benchmark_threshold_background"]
    for name, (comparison, title, default) in zip(names, panels):
        fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
        rows = cutoffs[cutoffs.comparison.eq(comparison)].copy()
        rows["setting"] = rows.setting.astype(float)
        rows["tr"], rows["ntr"] = rows.TR_pass / rows.TR_total, rows.NTR_pass / rows.NTR_total
        rows = rows.drop_duplicates(["tr", "ntr"], keep="last").sort_values("setting")
        ax.plot(rows.ntr, rows.tr, color=NEUTRAL, lw=1.5, zorder=1)
        for row in rows.itertuples():
            chosen = abs(row.setting - default) < 1e-9
            ax.scatter([row.ntr], [row.tr], s=70 if chosen else 45, color=BLUE if chosen else INK2,
                       edgecolor=SURFACE, linewidth=1.5, zorder=3)
            ax.annotate(f"{row.setting:g}" + (" (default)" if chosen else ""), (row.ntr, row.tr),
                        textcoords="offset points", xytext=(8, -12), fontsize=8.5, color=INK)
        ax.plot([0, 1], [0, 1], color=INK2, lw=1, ls=(0, (4, 3)))
        ax.set_xlim(0, 0.7)
        ax.set_ylim(0.5, 1)
        percent_axis(ax)
        percent_axis(ax, "y")
        ax.set_xlabel("tested non-reactive clones passing")
        ax.set_title(title, fontsize=10)
        ax.grid(lw=0.6)
        ax.set_axisbelow(True)
        ax.set_ylabel("tumor-reactive clones passing")
        save(fig, name)


def cd3_gate(dataset, tissue, order, versus, suffix=""):
    """Cells kept and ranking of validated clones under three CD3 floors."""
    runs = json.loads((RESULTS / "summary.json").read_text())[dataset]["runs"]
    aurocs = pd.read_csv(RESULTS / "signature_auroc.csv")
    aurocs = aurocs[aurocs.dataset.eq(dataset) & aurocs.tissue.eq(tissue) & aurocs.signature.eq("TumorReactive")]
    labels = [label for _, label in reversed(order)]
    fig1, ax1 = plt.subplots(figsize=(6.6, 3.0), constrained_layout=True)
    fig2, ax2 = plt.subplots(figsize=(7.4, 3.2), constrained_layout=True)
    for i, (key, _) in enumerate(order):
        run = runs[key]
        lineage = next(f["after"] for f in run["cell_filters"] if f["filter"] == "lineage")
        y = len(order) - 1 - i
        ax1.barh(y, run["cells_retained"], height=0.6, color=BLUE)
        ax1.text(run["cells_retained"] + lineage * 0.015, y,
                 f"{run['cells_retained']:,}  ({run['cells_retained'] / lineage:.0%} of CD8)",
                 va="center", fontsize=9, color=INK)
        row = aurocs[aurocs.run.eq(key)].iloc[0]
        ax2.plot([row.auroc_low, row.auroc_high], [y, y], color=BLUE, lw=2, solid_capstyle="round")
        ax2.scatter([row.auroc], [y], s=55, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax2.text(1.04, y, f"{row.auroc:.2f}  ({row.n_positive} vs {row.n_negative} clones)", va="center",
                 fontsize=9, color=INK)
    ax1.set_yticks(range(len(order)), labels)
    ax1.set_xticks([])
    ax1.set_xlim(0, ax1.get_xlim()[1] * 1.45)
    ax1.tick_params(axis="y", length=0)
    for side in ("left", "bottom"):
        ax1.spines[side].set_visible(False)
    ax1.set_title("CD8 cells retained with a paired TCR", fontsize=10)
    ax2.axvline(0.5, color=INK2, lw=1, ls=(0, (4, 3)))
    ax2.set_xlim(0.3, 1.0)
    ax2.set_ylim(-0.6, len(order) - 0.4)
    ax2.set_yticks(range(len(order)), labels)
    ax2.tick_params(axis="y", length=0)
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="x", lw=0.6)
    ax2.set_axisbelow(True)
    ax2.set_xlabel(f"TumorReactive AUROC, {versus} clones")
    ax2.set_title("Ranking of validated clones", fontsize=10)
    save(fig1, f"benchmark_cd3_cells{suffix}")
    save(fig2, f"benchmark_cd3_ranking{suffix}")


def dataset_designs():
    """Schematic of the five source studies: who, what was sequenced, how
    specificity was established, and which labels the benchmark uses.
    Drawn here from each paper's abstract and deposit; not a paper figure."""
    columns = ["Patients", "Sequenced", "How specificity was established", "Labels used here"]
    rows = [
        ("Meng 2023\npancreatic", [
            "9 patients with\npancreatic ductal\nadenocarcinoma",
            "Single-cell RNA + TCR\nof tumor-infiltrating\nCD8 T cells",
            "TCRs tested for\nreactivity against the\npatient's own tumor cells",
            "94 reactive and\n55 non-reactive\nclonotypes"]),
        ("Caushi 2021\nlung", [
            "8 of 16 patients with\nlung cancer treated with\nanti-PD-1 before surgery",
            "Single-cell RNA + TCR\nof T cells from tumor,\nlymph node, normal lung",
            "MANAFEST: T cells\nexpanded with neoantigen\nor viral peptides",
            "Neoantigen-specific vs\nvirus-specific\nbeta chains"]),
        ("Eberhardt 2021\nhead and neck", [
            "3 of 6 patients with\nHPV-positive head\nand neck cancer",
            "Single-cell RNA + TCR\nof tetramer-sorted and\nbulk PD-1+ CD8 TIL",
            "MHC class I tetramers\nloaded with HPV E2\npeptides",
            "Tetramer-positive\nbeta chains;\nno negatives"]),
        ("Chiffelle 2024\nmelanoma", [
            "11 of 13 patients with\nmelanoma given TIL\ntherapy",
            "Single-cell RNA + TCR\nof CD45+ cells from the\ntumor before treatment",
            "TCRs cloned and tested\nagainst the patient's\nown tumor cells",
            "72 reactive and\n94 non-reactive\nclonotypes"]),
        ("Pulliam 2024\nMerkel cell", [
            "2 of 7 patients with\nMerkel cell carcinoma",
            "Single-cell RNA + TCR\nof hashed tumor and\nblood cells",
            "MHC class I multimers\nloaded with polyomavirus\nor other viral peptides",
            "Polyomavirus-specific\nbeta chains;\nno negatives"]),
    ]
    fig, ax = plt.subplots(figsize=(11, 7.4), constrained_layout=True)
    width, height, x0, dx, dy = 2.25, 1.0, 1.75, 2.6, 1.35
    for c, name in enumerate(columns):
        ax.text(x0 + c * dx + width / 2, 0.42, name, ha="center", va="center", fontsize=9.5,
                fontweight="bold", color=INK)
    for r, (study, cells) in enumerate(rows):
        y = -r * dy - 0.25
        ax.text(1.5, y - height / 2, study, ha="right", va="center", fontsize=9.5, fontweight="bold", color=INK)
        for c, words in enumerate(cells):
            x = x0 + c * dx
            last = c == len(cells) - 1
            ax.add_patch(matplotlib.patches.FancyBboxPatch(
                (x, y - height), width, height, boxstyle="round,pad=0,rounding_size=0.08",
                facecolor=BLUE if last else SURFACE, edgecolor=BLUE if last else NEUTRAL, linewidth=1.3,
                alpha=1.0))
            ax.text(x + width / 2, y - height / 2, words, ha="center", va="center", fontsize=8.6,
                    color=SURFACE if last else INK, linespacing=1.25)
            if not last:
                ax.annotate("", xy=(x + dx - 0.04, y - height / 2), xytext=(x + width + 0.04, y - height / 2),
                            arrowprops={"arrowstyle": "->", "color": INK2, "lw": 1.1})
    ax.set_xlim(0, x0 + 3 * dx + width + 0.1)
    ax.set_ylim(-(len(rows) - 1) * dy - height - 0.4, 0.7)
    ax.axis("off")
    save(fig, "benchmark_datasets")


if __name__ == "__main__":
    dataset_designs()
    signature_auroc()
    auroc_by_dataset()
    cutoff_pass_rates()
    selection_by_budget()
    threshold_sweeps()
    floors = ["No floor", "3 counts (default)", "10 counts"]
    cd3_gate("caushi", "all", list(zip(["tumor_cd3_0", "tumor_default", "tumor_cd3_10"], floors)),
             "neoantigen vs viral")
    cd3_gate("chiffelle", "tumor", list(zip(["cd3_0", "default", "cd3_10"], floors)), "reactive vs non-reactive",
             suffix="_melanoma")
    print("wrote benchmark figures")
