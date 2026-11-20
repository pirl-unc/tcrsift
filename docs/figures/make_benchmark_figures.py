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
    """Cells kept and labelled clones scored under three CD3 floors (lung data)."""
    runs = json.loads((RESULTS / "summary.json").read_text())["caushi"]["runs"]
    order = [("tumor_cd3_0", "No floor"), ("tumor_default", "3 counts (default)"), ("tumor_cd3_10", "10 counts")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True)
    for i, (key, label) in enumerate(order):
        run = runs[key]
        lineage = next(f["after"] for f in run["cell_filters"] if f["filter"] == "lineage")
        y = len(order) - 1 - i
        ax1.barh(y, run["cells_retained"], height=0.6, color=BLUE)
        ax1.text(run["cells_retained"] + lineage * 0.015, y,
                 f"{run['cells_retained']:,}  ({run['cells_retained'] / lineage:.0%} of CD8)",
                 va="center", fontsize=9, color=INK)
        scored = run["labelled_scored"]
        selected = run["labelled_selected"]
        ax2.barh(y + 0.19, scored.get("MANA", 0), height=0.34, color=NEUTRAL)
        ax2.barh(y - 0.19, selected.get("MANA", 0), height=0.34, color=BLUE)
        ax2.text(scored.get("MANA", 0) + 0.4, y + 0.19, f"{scored.get('MANA', 0)} scored", va="center", fontsize=8.5)
        ax2.text(selected.get("MANA", 0) + 0.4, y - 0.19, f"{selected.get('MANA', 0)} shortlisted",
                 va="center", fontsize=8.5)
    for ax, title in ((ax1, "CD8 cells retained with a paired TCR"), (ax2, "Neoantigen-specific clonotypes")):
        ax.set_yticks(range(len(order)), [label for _, label in reversed(order)])
        ax.set_xticks([])
        ax.tick_params(axis="y", length=0)
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, ax.get_xlim()[1] * 1.45)
    save(fig, "benchmark_cd3_gate")


if __name__ == "__main__":
    signature_auroc()
    cutoff_pass_rates()
    selection_by_budget()
    threshold_sweeps()
    cd3_gate()
    print("wrote benchmark figures")
