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
    caushi = table[table.dataset.eq("caushi") & table.run.eq("tumor_default") & table.tissue.eq("all")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    m, c = meng.iloc[0], caushi.iloc[0]
    auroc_forest(ax1, meng, f"PDAC (Meng 2023): {m.n_positive} reactive vs "
                            f"{m.n_negative} non-reactive clones", note=True)
    auroc_forest(ax2, caushi, f"NSCLC (Caushi 2021): {c.n_positive} neoantigen vs "
                              f"{c.n_negative} viral clones", note=True)
    save(fig, "benchmark_auroc")


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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1, 1.5]})
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
    save(fig, "benchmark_cutoffs")


def selection_by_budget():
    """What each way of picking N clones puts on the list (tested clones only)."""
    table = pd.read_csv(RESULTS / "meng_selection.csv")
    names = {"TCRsift round-robin": "Round-robin, solid-tumor preset",
             "TCRsift round-robin, tumor-reactivity signatures only": "Round-robin, 3 tumor signatures",
             "TumorReactive score only": "Top TumorReactive score",
             "largest clones": "Largest clones",
             "random abundant clones (mean of 500)": "Random clones"}
    budgets = ["25", "50", "100", "200"]
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.6), constrained_layout=True, sharey=True)
    limit = table[table.budget.astype(str).isin(budgets)][["TR", "NTR"]].to_numpy().max() * 1.18
    height = 0.36
    for ax, budget in zip(axes, budgets):
        rows = table[table.budget.astype(str).eq(budget)].set_index("method").loc[list(names)][::-1]
        for i, row in enumerate(rows.itertuples()):
            ax.barh(i + height / 2 + 0.02, row.TR, height=height, color=BLUE)
            ax.barh(i - height / 2 - 0.02, row.NTR, height=height, color=ORANGE)
            ax.text(row.TR + limit * 0.02, i + height / 2 + 0.02, f"{row.TR:.0f}", va="center", fontsize=8.5, color=INK)
            ax.text(row.NTR + limit * 0.02, i - height / 2 - 0.02, f"{row.NTR:.0f}", va="center", fontsize=8.5, color=INK)
        ax.set_yticks(range(len(rows)), [names[m] for m in rows.index])
        ax.set_xlim(0, limit)
        ax.set_title(f"{budget} clones chosen", fontsize=10)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x", lw=0.6)
        ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=BLUE), plt.Rectangle((0, 0), 1, 1, color=ORANGE)]
    fig.legend(handles, ["Tumor-reactive clones chosen", "Tested non-reactive clones chosen"], frameon=False,
               fontsize=9, loc="outside lower center", ncol=2)
    save(fig, "benchmark_selection")


def threshold_sweeps():
    """Reactive against non-reactive pass rate as each threshold moves."""
    cutoffs = pd.read_csv(RESULTS / "meng_cutoffs.csv")
    panels = [("rank_gate", "Rank gate (--signature-quantile)", 0.9),
              ("background_quantile", "Background floor (--signature-background-quantile)", 0.99)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True, sharex=True, sharey=True)
    for ax, (comparison, title, default) in zip(axes, panels):
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
    axes[0].set_ylabel("tumor-reactive clones passing")
    save(fig, "benchmark_thresholds")


def cd3_gate():
    """Cells kept and ranking of validated clones under three CD3 floors (lung data)."""
    runs = json.loads((RESULTS / "summary.json").read_text())["caushi"]["runs"]
    aurocs = pd.read_csv(RESULTS / "signature_auroc.csv")
    aurocs = aurocs[aurocs.dataset.eq("caushi") & aurocs.tissue.eq("all") & aurocs.signature.eq("TumorReactive")]
    order = [("tumor_cd3_0", "No floor"), ("tumor_cd3_3", "3 counts"), ("tumor_default", "10 counts (default)")]
    labels = [label for _, label in reversed(order)]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.2), constrained_layout=True)
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
    ax2.set_xlabel("TumorReactive AUROC, neoantigen vs viral clones")
    ax2.set_title("Ranking of validated clones", fontsize=10)
    save(fig, "benchmark_cd3_gate")


def dataset_designs():
    """Schematic of the three source studies: who, what was sequenced, how
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
    ]
    fig, ax = plt.subplots(figsize=(11, 4.6), constrained_layout=True)
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
    ax.set_ylim(-2 * dy - height - 0.4, 0.7)
    ax.axis("off")
    save(fig, "benchmark_datasets")


if __name__ == "__main__":
    dataset_designs()
    signature_auroc()
    cutoff_pass_rates()
    selection_by_budget()
    threshold_sweeps()
    cd3_gate()
    print("wrote benchmark figures")
