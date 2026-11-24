#!/usr/bin/env python3
"""Regenerate the explanatory figures used in the documentation.

    python docs/figures/make_figures.py

Every panel here is drawn from SIMULATED per-cell scores (fixed seed). The
figures illustrate how the selection rule works; they are not measurements of
performance on real samples. The background floors are computed with the same
helper the ``prioritize`` command uses (``signature_background._background_means``)
and the subset tail with the same hypergeometric survival function, so the
pictures track the implementation rather than a re-description of it.
"""

from __future__ import annotations

from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom

from tcrsift.signature_background import _background_means

OUT = Path(__file__).resolve().parent.parent / "assets" / "figures"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
NEUTRAL = "#b9b8b2"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 10, "text.color": INK,
    "axes.edgecolor": NEUTRAL, "axes.labelcolor": INK2, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "grid.color": GRID,
    "svg.fonttype": "none", "svg.hashsalt": "tcrsift-docs",
})


def simulated_cell_scores(n=3000, seed=0):
    """A right-skewed score population: most cells low, a minority high."""
    rng = np.random.default_rng(seed)
    high = rng.random(n) < 0.12
    return np.where(high, rng.normal(1.1, 0.45, n), rng.normal(-0.15, 0.3, n))


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", pad_inches=0.08, metadata={"Date": None})
    plt.close(fig)


def background_floor():
    """Mean route: the floor a clone must exceed depends on its cell count."""
    cells = simulated_cell_scores()
    matrix = cells[:, None]
    sizes = [2, 3, 5, 8, 12, 20, 35, 60]
    nulls = {n: _background_means(matrix, n, 2000, seed=n)[0][:, 0] for n in sizes}
    floors = {n: float(np.quantile(nulls[n], 0.99, method="higher")) for n in sizes}

    fig1, ax1 = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    fig2, ax2 = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    bins = np.linspace(-1.0, 2.4, 70)
    ax1.hist(cells, bins=bins, color=NEUTRAL, density=True, label="single cells")
    for n, color in ((2, BLUE), (20, ORANGE)):
        ax1.hist(nulls[n], bins=bins, density=True, histtype="step", lw=2, color=color,
                 label=f"means of random {n}-cell groups")
        ax1.axvline(floors[n], color=color, lw=1.2, ls=(0, (4, 3)))
        ax1.text(floors[n] + 0.04, ax1.get_ylim()[1] * (0.42 if n == 2 else 0.6),
                 f"99th pct.\n(n={n})", color=INK2, fontsize=8.5, va="top")
    ax1.set(title="Random groups of n cells set the reference",
            xlabel="signature score (simulated)", ylabel="density")
    ax1.set_ylim(0, ax1.get_ylim()[1] * 1.05)
    ax1.legend(frameon=False, fontsize=8.5, loc="upper right")

    ax2.plot(sizes, [floors[n] for n in sizes], color=INK2, lw=2, marker="o", ms=4,
             label="background floor (99th pct.)")
    ax2.fill_between(sizes, [floors[n] for n in sizes], -0.2, color=GRID, lw=0)
    examples = [(2, 1.05, True, "A"), (5, 0.42, False, "B"), (20, 0.42, True, "C"), (60, 0.08, False, "D")]
    for n, mean, passes, tag in examples:
        ax2.scatter([n], [mean], s=70, color=AQUA if passes else SURFACE,
                    edgecolor=AQUA if passes else INK2, linewidth=1.8, zorder=3)
        ax2.annotate(f"{tag}: {'above' if passes else 'below'} floor", (n, mean),
                     textcoords="offset points", xytext=(8, 5) if passes else (8, -14),
                     fontsize=8.5, color=INK)
    ax2.set_xscale("log")
    ax2.set_xticks(sizes, [str(n) for n in sizes])
    ax2.minorticks_off()
    ax2.set_ylim(-0.2, 1.45)
    ax2.set(title="Larger clones face a lower floor",
            xlabel="cells in clone (n)", ylabel="clone mean score")
    ax2.grid(axis="y", lw=0.6)
    ax2.legend(frameon=False, fontsize=8.5, loc="upper right")
    save(fig1, "background_reference")
    save(fig2, "background_floor")


def subset_evidence():
    """Subset route: 10 of 25 cells high even though the clone mean is ordinary."""
    cells = simulated_cell_scores()
    threshold = float(np.quantile(cells, 0.90, method="higher"))
    total, total_high = len(cells), int((cells > threshold).sum())
    rng = np.random.default_rng(3)
    clone = np.concatenate([rng.normal(1.25, 0.2, 10).clip(threshold + 0.05),
                            rng.normal(-0.45, 0.15, 15)])
    n, k = len(clone), int((clone > threshold).sum())
    floor = float(np.quantile(_background_means(cells[:, None], n, 2000, seed=1)[0][:, 0],
                              0.99, method="higher"))
    tail = hypergeom.sf(k - 1, total, total_high, n)

    fig1, ax1 = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    fig2, ax2 = plt.subplots(figsize=(6.4, 3.9), constrained_layout=True)
    jitter = rng.uniform(-0.28, 0.28, n)
    high = clone > threshold
    ax1.scatter(clone[~high], jitter[~high], s=42, color=NEUTRAL, edgecolor=SURFACE,
                label=f"other cells ({n - k})")
    ax1.scatter(clone[high], jitter[high], s=42, color=BLUE, edgecolor=SURFACE,
                label=f"high-scoring cells ({k})")
    ymax = 0.62
    ax1.axvline(threshold, color=INK2, lw=1.2, ls=(0, (4, 3)))
    ax1.text(threshold + 0.03, ymax, "sample's 90th-pct.\ncell score", fontsize=8.5,
             color=INK2, va="top")
    ax1.axvline(clone.mean(), color=ORANGE, lw=2)
    ax1.text(clone.mean() - 0.03, ymax, f"clone mean\n{clone.mean():.2f}", fontsize=8.5,
             color=INK, va="top", ha="right")
    ax1.axvline(floor, color=INK2, lw=1.2)
    ax1.text(floor + 0.03, -ymax, f"mean-route floor\n(n={n}): {floor:.2f}", fontsize=8.5,
             color=INK2, va="bottom")
    ax1.set_ylim(-0.68, 0.68)
    ax1.set_yticks([])
    ax1.spines["left"].set_visible(False)
    ax1.set(title=f"One clone, {n} cells: mean fails, subset qualifies",
            xlabel="signature score (simulated)")
    ax1.legend(frameon=False, fontsize=8.5, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.3))

    ks = np.arange(0, 13)
    pmf = hypergeom.pmf(ks, total, total_high, n)
    ax2.bar(ks, pmf, width=0.72, color=[BLUE if x >= k else NEUTRAL for x in ks])
    ax2.annotate(f"observed k = {k}\nP(X ≥ {k}) = {tail:.1e}", (k, pmf[k]),
                 xytext=(k - 0.4, 0.14), fontsize=9, color=INK, ha="center",
                 arrowprops={"arrowstyle": "-", "color": INK2, "lw": 1})
    ax2.set_xticks(ks)
    ax2.set(title=f"High cells expected in a random {n}-cell group",
            xlabel=f"high-scoring cells among {n}", ylabel="probability")
    ax2.grid(axis="y", lw=0.6)
    ax2.set_axisbelow(True)
    save(fig1, "subset_cells")
    save(fig2, "subset_probability")


def round_robin():
    """Priority-tiered budget: Tier-1 reactivity lists exhaust before Tier-2 backfill."""
    lists = {
        "NeoTCR8\nmean": ["c1", "c2", "c3", "c4", "c5"],
        "NeoTCR8\nsubset": ["c6", "c1"],
        "MANAscore\nmean": ["c2", "c7", "c8"],
        "Proliferation\nmean": ["c9", "c3", "c10", "c11"],
    }
    tier1_cols = {0, 1, 2}
    budget = 9
    queues = [deque(v) for v in lists.values()]
    picked, order = {}, {}
    for pool in ([0, 1, 2], [3]):
        while len(picked) < budget:
            before = len(picked)
            for j in pool:
                queue = queues[j]
                while queue and queue[0] in picked:
                    queue.popleft()
                if not queue:
                    continue
                clone = queue.popleft()
                picked[clone] = j
                order[(j, clone)] = len(picked)
                if len(picked) == budget:
                    break
            if len(picked) == before:
                break

    fig, ax = plt.subplots(figsize=(8.6, 4.5), constrained_layout=True)
    ax.text(1.0, 1.38, "Tier 1: context reactivity lists (exhausted first)",
            ha="center", va="bottom", fontsize=8.5, color=BLUE, fontweight="bold")
    ax.plot([-0.38, 2.38], [1.32, 1.32], color=BLUE, lw=1.2)
    ax.text(3.0, 1.38, "Tier 2: backfill",
            ha="center", va="bottom", fontsize=8.5, color=INK2, fontweight="bold")
    ax.plot([2.62, 3.38], [1.32, 1.32], color=NEUTRAL, lw=1.2)
    for j, (name, clones) in enumerate(lists.items()):
        ax.text(j, 0.72, name, ha="center", va="bottom", fontsize=9.5, fontweight="bold")
        is_t1 = j in tier1_cols
        for i, clone in enumerate(clones):
            rank = order.get((j, clone))
            taken_elsewhere = rank is None and clone in picked
            face = (BLUE if is_t1 else "#4B6A88") if rank else SURFACE
            edge = (BLUE if is_t1 else "#4B6A88") if rank else NEUTRAL
            ax.add_patch(plt.Rectangle((j - 0.4, -i - 0.36), 0.8, 0.72, facecolor=face,
                                       edgecolor=edge, lw=1.5,
                                       ls="-" if rank or not taken_elsewhere else (0, (3, 2))))
            label = f"{clone}\npick {rank}" if rank else (
                f"{clone}\nalready picked" if taken_elsewhere else f"{clone}\nover budget")
            ax.text(j, -i, label, ha="center", va="center", fontsize=8.5,
                    color=SURFACE if rank else INK2)
    ax.text(-0.75, 0, "best", ha="right", va="center", fontsize=9, color=INK2)
    ax.annotate("", xy=(-0.9, -4.2), xytext=(-0.9, -0.5),
                arrowprops={"arrowstyle": "->", "color": INK2, "lw": 1})
    ax.text(-1.0, -2.35, "rank within list", rotation=90, ha="right", va="center",
            fontsize=9, color=INK2)
    ax.set_xlim(-1.3, len(lists) - 0.4)
    ax.set_ylim(-4.6, 1.8)
    ax.axis("off")
    ax.set_title(f"Priority-tiered round-robin over qualifying lists (budget = {budget} clones)", pad=4)
    save(fig, "round_robin")


if __name__ == "__main__":
    background_floor()
    subset_evidence()
    round_robin()
    print(f"wrote figures to {OUT}")
