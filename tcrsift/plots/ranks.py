# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Rank-frequency, cumulative coverage, slope-tracking, and per-sample scatter plots."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..format import pretty_sample, pretty_samples
from .annotations import _expand_sample_frequencies
from .state import (
    _grid_shape,
    _parse_sample_frequencies,
    _per_cell_signature,
    save_figure,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Rank-family Plots — #43 (#2/#3/#4/#5)
# =============================================================================
#
# Four plots built on a single "rank clones within a sample" kernel:
#
#   1. Per-sample log-log rank curves with optional tier highlight.
#   2. Cumulative coverage by clone rank across samples.
#   3. Top-N labeled bar chart per sample.
#   4. Top-N cumulative stacked bar per sample (mono/oligo-clonality).
#
# Plus a `draw_reference_fractions` log-y utility that all four reuse.
# Reference fractions: 10/1/0.1/0.01/0.001%.


# Reference-fraction defaults for log-y frequency axes. Powers of ten
# from 10% down to 0.001%; helpers above can override.
REFERENCE_FRACTIONS = (0.1, 0.01, 0.001, 1e-4, 1e-5)


def draw_reference_fractions(
    ax: plt.Axes,
    fractions: tuple[float, ...] = REFERENCE_FRACTIONS,
    label: bool = True,
    color: str = "#999999",
    linestyle: str = ":",
    linewidth: float = 0.8,
) -> None:
    """Draw labeled horizontal dashed lines at the given fractions.

    Designed for log-y frequency axes — labels render as percent
    strings on the right edge so they don't fight with data lines.
    Skips fractions outside the current y-limits.
    """
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    for f in fractions:
        if not (ymin <= f <= ymax):
            continue
        ax.axhline(f, color=color, linestyle=linestyle, linewidth=linewidth, zorder=1)
        if label:
            pct = f"{f * 100:g}%"
            ax.text(
                xmax,
                f,
                f" {pct}",
                ha="left",
                va="center",
                fontsize=8,
                color=color,
                clip_on=False,
            )


def _per_sample_rank_table(clonotypes: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build a per-sample sorted-by-frequency view.

    Returns a dict keyed by sample name; each value is a DataFrame
    indexed by 1-based rank with columns ``CDR3ab``, ``frequency``,
    and ``tier`` (if the column exists; else None). Clones with zero
    frequency in that sample are excluded.
    """
    wide = _expand_sample_frequencies(clonotypes)
    if wide.empty or wide.shape[1] == 0:
        return {}

    tier_lookup = (
        clonotypes.dropna(subset=["tier"])
        .drop_duplicates(subset=["CDR3ab"], keep="first")
        .set_index("CDR3ab")["tier"]
        if "tier" in clonotypes.columns
        else None
    )

    out: dict[str, pd.DataFrame] = {}
    for sample in wide.columns:
        col = wide[sample]
        present = col[col > 0].sort_values(ascending=False)
        if len(present) == 0:
            continue
        sample_df = pd.DataFrame(
            {
                "CDR3ab": present.index.tolist(),
                "frequency": present.values,
            }
        )
        if tier_lookup is not None:
            sample_df["tier"] = sample_df["CDR3ab"].map(tier_lookup)
        sample_df.index = range(1, len(sample_df) + 1)
        sample_df.index.name = "rank"
        out[sample] = sample_df
    return out


def plot_rank_curves_per_sample(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    filename: str = "rank_curves_per_sample.png",
    x_scale: str = "log",
) -> Path | None:
    """Log-log rank vs within-sample frequency, one panel per sample.

    Tier-selected clones (``tier`` column non-null) get highlighted as
    red dots overlaid on the rank curve so it's obvious which clones
    are picked from where on the distribution. Reference-fraction
    lines (10/1/0.1%/...) make the y-scale readable.

    ``x_scale`` is ``"log"`` (default — best for whole-distribution
    views) or ``"linear"`` (matches the "top-X linear-x" variant
    mentioned in #43).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = _per_sample_rank_table(clonotypes)
    if not tables:
        logger.info("plot_rank_curves_per_sample: no sample_frequencies; skipping.")
        return None

    samples = list(tables)
    rows, cols = _grid_shape(len(samples))
    fig, axes = plt.subplots(
        rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False, sharey=True
    )

    for idx, sample in enumerate(samples):
        ax = axes[idx // cols][idx % cols]
        df = tables[sample]
        ax.plot(df.index, df["frequency"], color="#1f77b4", linewidth=1.2)

        if "tier" in df.columns:
            tier_hits = df[df["tier"].notna()]
            if len(tier_hits) > 0:
                ax.scatter(
                    tier_hits.index,
                    tier_hits["frequency"],
                    color="#d62728",
                    s=20,
                    zorder=3,
                    label=f"tier-selected (n={len(tier_hits)})",
                )
                ax.legend(loc="upper right", fontsize=9)

        ax.set_xscale(x_scale)
        ax.set_yscale("log")
        ax.set_xlabel("clone rank")
        ax.set_ylabel("within-sample frequency")

        ax.set_title(pretty_sample(sample), fontsize=11)
        draw_reference_fractions(ax)

    # Blank any unused panels.
    for unused in range(len(samples), rows * cols):
        axes[unused // cols][unused % cols].axis("off")

    fig.suptitle(
        f"Clone-frequency rank curves per sample ({x_scale}-x)", y=1.02, fontsize=13
    )
    out_path = output_dir / filename
    return save_figure(fig, out_path)


def plot_cumulative_coverage(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    filename: str = "cumulative_coverage_per_sample.png",
    reference_levels: tuple[float, ...] = (0.5, 0.9),
) -> Path | None:
    """Σ(top-k frequency) vs k, one line per sample on a shared axes.

    Reference horizontal lines at 50% / 90% (configurable) make
    clonality differences across samples visually unambiguous: a
    sample that hits 90% coverage by clone 10 is much more clonal
    than one that needs 1000.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = _per_sample_rank_table(clonotypes)
    if not tables:
        logger.info("plot_cumulative_coverage: no sample_frequencies; skipping.")
        return None

    fig, ax = plt.subplots(figsize=(10, 6))
    palette = plt.get_cmap("tab10")
    for i, (sample, df) in enumerate(tables.items()):
        cum = df["frequency"].cumsum().values
        ax.plot(
            range(1, len(cum) + 1),
            cum,
            label=sample,
            color=palette(i % 10),
            linewidth=1.5,
        )

    for level in reference_levels:
        ax.axhline(level, color="#999999", linestyle=":", linewidth=0.8, zorder=1)
        ax.text(
            ax.get_xlim()[1],
            level,
            f" {level * 100:g}%",
            ha="left",
            va="center",
            fontsize=8,
            color="#666666",
            clip_on=False,
        )

    ax.set_xscale("log")
    ax.set_xlabel("top-k clones")
    ax.set_ylabel("cumulative frequency")
    ax.set_ylim(0, 1.05)
    ax.set_title("Cumulative coverage by clone rank")
    ax.legend(loc="lower right", fontsize=9, title="sample")

    out_path = output_dir / filename
    return save_figure(fig, out_path)


def plot_top_n_labeled_bars_per_sample(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    n: int = 40,
    annotate_top: int = 5,
    filename: str = "top_n_labeled_bars_per_sample.png",
) -> Path | None:
    """Top-N clones per sample as a labeled bar chart, multi-panel.

    Each panel: a bar per clone in the sample's top-N. X labels are
    CDR3αβ (rotated and truncated for readability). Bar color
    indicates tier-selected vs. not when the ``tier`` column is
    available. The top ``annotate_top`` bars carry a small cell-count
    annotation above them when ``cell_count`` is present.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = _per_sample_rank_table(clonotypes)
    if not tables:
        logger.info(
            "plot_top_n_labeled_bars_per_sample: no sample_frequencies; skipping."
        )
        return None

    cell_count_lookup = (
        clonotypes.drop_duplicates(subset=["CDR3ab"], keep="first")
        .set_index("CDR3ab")["cell_count"]
        if "cell_count" in clonotypes.columns
        else None
    )

    samples = list(tables)
    rows, cols = _grid_shape(len(samples), max_cols=2)
    fig, axes = plt.subplots(
        rows, cols, figsize=(8 * cols, 4 * rows), squeeze=False
    )

    for idx, sample in enumerate(samples):
        ax = axes[idx // cols][idx % cols]
        df = tables[sample].head(n)
        positions = range(len(df))

        if "tier" in df.columns:
            colors = [
                "#d62728" if pd.notna(t) else "#7f7f7f" for t in df["tier"]
            ]
        else:
            colors = "#1f77b4"

        ax.bar(positions, df["frequency"].values, color=colors, edgecolor="white")
        ax.set_yscale("log")
        ax.set_xticks(positions)
        # Show CDR3αβ but truncate to avoid an unreadable axis on long
        # paired strings (e.g. ``CASS...F_CAVR...F``).
        labels = [
            (cdr3 if len(cdr3) <= 20 else cdr3[:17] + "…")
            for cdr3 in df["CDR3ab"]
        ]
        ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
        ax.set_ylabel("within-sample freq")

        ax.set_title(f"{pretty_sample(sample)} — top {len(df)} clones", fontsize=11)
        draw_reference_fractions(ax)

        # Cell-count annotation on the leading bars.
        if cell_count_lookup is not None:
            top_to_label = df.head(annotate_top)
            for pos, (_, row) in zip(positions, top_to_label.iterrows()):
                count = cell_count_lookup.get(row["CDR3ab"])
                if pd.notna(count):
                    ax.annotate(
                        f"{int(count)}",
                        xy=(pos, row["frequency"]),
                        xytext=(0, 4),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=7,
                        color="#333333",
                    )

    for unused in range(len(samples), rows * cols):
        axes[unused // cols][unused % cols].axis("off")

    fig.suptitle(f"Top-{n} clones per sample", y=1.02, fontsize=13)
    out_path = output_dir / filename
    return save_figure(fig, out_path)


def plot_top_n_cumulative_stacked(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    n: int = 20,
    filename: str = "top_n_cumulative_stacked.png",
) -> Path | None:
    """Per-sample stacked bar showing the contribution of clones 1..N.

    One bar per sample, each segment is a single clone's within-sample
    frequency. Visually conveys mono/oligo-clonality at a glance: a
    sample dominated by one segment is monoclonal, one with many even
    segments is polyclonal.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = _per_sample_rank_table(clonotypes)
    if not tables:
        logger.info("plot_top_n_cumulative_stacked: no sample_frequencies; skipping.")
        return None

    # Build a tidy long table of (sample, rank, frequency) for the top N.
    rows_data = []
    for sample, df in tables.items():
        head = df.head(n)
        for rank, freq in zip(head.index, head["frequency"]):
            rows_data.append({"sample": sample, "rank": rank, "frequency": freq})
    if not rows_data:
        return None
    long = pd.DataFrame(rows_data)
    pivot = long.pivot(index="sample", columns="rank", values="frequency").fillna(0.0)
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)

    fig, ax = plt.subplots(figsize=(max(8, 1.0 * len(pivot)), 6))
    bottom = pd.Series(0.0, index=pivot.index)
    # Rank-1 darkest → rank-N lightest, so the eye reads the leading
    # contributor first.
    cmap = plt.get_cmap("viridis_r")
    for rank in pivot.columns:
        ax.bar(
            pivot.index,
            pivot[rank].values,
            bottom=bottom.values,
            color=cmap((rank - 1) / max(1, n - 1)),
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += pivot[rank]

    ax.set_ylabel("cumulative within-sample frequency")
    ax.set_xlabel("sample")
    ax.set_title(f"Top-{n} cumulative frequency per sample")
    ax.set_ylim(0, max(1.0, bottom.max() * 1.05))
    plt.xticks(rotation=45, ha="right")

    # Colorbar for clone rank.
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=n))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("clone rank within sample")

    out_path = output_dir / filename
    return save_figure(fig, out_path)


def plot_rank_family(
    clonotypes: pd.DataFrame, output_dir: str | Path, top_n: int = 40
) -> None:
    """One-shot orchestrator for the rank-family plots.

    Emits the four plots described in #43 (#2/#3/#4/#5) under
    ``output_dir``. Each helper skips itself when prerequisites are
    missing, so this is safe on partial data.
    """
    plot_rank_curves_per_sample(clonotypes, output_dir, x_scale="log")
    plot_rank_curves_per_sample(
        clonotypes,
        output_dir,
        filename="rank_curves_per_sample_linear.png",
        x_scale="linear",
    )
    plot_cumulative_coverage(clonotypes, output_dir)
    plot_top_n_labeled_bars_per_sample(clonotypes, output_dir, n=top_n)
    plot_top_n_cumulative_stacked(clonotypes, output_dir, n=min(20, top_n))


# =============================================================================
# Multi-Panel Grid Helper — #43 method-grouping helper
# =============================================================================
#
# A generic version of the wrapper's `make_method_grid` / `method_panel_layout`
# that doesn't hardcode any cohort's specific method set. Callers who want
# a cohort-specific layout (e.g. the MART-1 cohort's 3-2-2 layout) can pass
# it as `layout`; otherwise the helper falls back to row-major chunks.


def make_method_panel_grid(
    methods: list[str],
    panel_size: tuple[float, float] = (4.5, 3.5),
    cols: int = 3,
    layout: list[list[str]] | None = None,
    sharex: bool = False,
    sharey: bool = False,
) -> tuple[plt.Figure, dict[str, plt.Axes]]:
    """Build a multi-panel grid keyed by method/sample name.

    Returns ``(fig, {method: ax})``. When ``layout`` is given (a
    list-of-lists), it's used verbatim — useful for cohort-specific
    layouts where panel groupings carry meaning. Otherwise methods are
    chunked row-major into ``cols`` columns.

    Designed so panel-keyed plots (slope chart, signature scatter,
    overlap matrices) all share the same multi-panel skeleton.
    """
    if layout is None:
        layout = [methods[i : i + cols] for i in range(0, len(methods), cols)]
    n_rows = len(layout)
    n_cols = max((len(row) for row in layout), default=1)
    w, h = panel_size
    fig = plt.figure(figsize=(w * n_cols, h * n_rows))
    gs = fig.add_gridspec(n_rows, n_cols)

    axes_map: dict[str, plt.Axes] = {}
    first_ax: plt.Axes | None = None
    for r, row in enumerate(layout):
        for c, method in enumerate(row):
            kw: dict = {}
            if sharex and first_ax is not None:
                kw["sharex"] = first_ax
            if sharey and first_ax is not None:
                kw["sharey"] = first_ax
            ax = fig.add_subplot(gs[r, c], **kw)
            if first_ax is None:
                first_ax = ax
            axes_map[method] = ax
    return fig, axes_map


# =============================================================================
# Slope chart — #43 item #1
# =============================================================================


def plot_clone_tracking_slopes(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    top_n: int = 20,
    floor: float = 1e-5,
    filename: str = "clone_tracking_slopes.png",
    layout: list[list[str]] | None = None,
) -> Path | None:
    """Slope chart: per source-sample, top-N clones' frequency across all samples.

    For each source sample, picks the top-N clones by within-source
    frequency and draws their frequencies across every sample as
    connected lines. Clones unique to the source plummet to the
    log-y floor; broadly shared clones stay high. The visual delta
    makes it obvious which top clones are source-private vs.
    cohort-public (#43 item #1).

    Skips silently when ``sample_frequencies`` is missing.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    wide = _expand_sample_frequencies(clonotypes)
    if wide.empty or wide.shape[1] == 0:
        logger.info("plot_clone_tracking_slopes: no sample_frequencies; skipping.")
        return None

    samples = list(wide.columns)
    fig, axes_map = make_method_panel_grid(
        samples, panel_size=(5.5, 3.4), layout=layout, sharey=True
    )

    palette = sns.color_palette("viridis", n_colors=top_n)
    y_max = max(float(wide.values.max()), floor * 2) * 1.5

    for source, ax in axes_map.items():
        # Pick top-N clones by frequency in this source.
        col = wide[source]
        present = col[col > 0].sort_values(ascending=False).head(top_n)
        if len(present) == 0:
            ax.axis("off")
            continue

        # Order x-axis: source first, then others by mean freq across
        # the source's top-N (decreasing) — same ordering as the wrapper.
        other = [s for s in samples if s != source]
        if other:
            other_means = wide.loc[present.index, other].mean(axis=0)
            other = list(other_means.sort_values(ascending=False).index)
        x_order = [source] + other
        x_pos = np.arange(len(x_order))

        mat = wide.loc[present.index, x_order].clip(lower=floor)
        for i in range(len(mat)):
            ax.plot(
                x_pos,
                mat.iloc[i].to_numpy(),
                color=palette[i],
                linewidth=1.1,
                alpha=0.85,
                marker="o",
                markersize=3.5,
            )

        ax.set_yscale("log")
        ax.set_ylim(floor / 2, y_max)
        ax.axvline(0, color="#dc2626", linewidth=0.7, alpha=0.5, zorder=1)
        ax.set_xticks(x_pos)
        # x_order is a list of sample names; pretty-print so AIM⁺ etc.
        # render instead of AIMpos-2 (#109).
        ax.set_xticklabels(pretty_samples(x_order), rotation=30, ha="right", fontsize=8)
        ax.set_title(pretty_sample(source), loc="left", fontsize=10)
        ax.grid(True, which="both", linewidth=0.3, alpha=0.4)
        draw_reference_fractions(ax)

    fig.suptitle(
        f"Top-{top_n} clones traced across samples (source marked at left of each panel)",
        fontsize=12,
        y=1.01,
    )
    fig.tight_layout()
    out_path = output_dir / filename
    return save_figure(fig, out_path)


# =============================================================================
# N-methods distribution — #43 item #6
# =============================================================================


def _derive_n_methods(clonotypes: pd.DataFrame) -> pd.Series | None:
    """Return per-clone count of distinct samples present.

    Prefers an explicit ``n_methods`` / ``n_conditions_present`` column
    if present; otherwise derives by counting non-zero entries in the
    per-clone ``sample_frequencies`` dict. Returns ``None`` when
    neither source is available.
    """
    if "n_methods" in clonotypes.columns:
        return clonotypes["n_methods"]
    if "n_conditions_present" in clonotypes.columns:
        return clonotypes["n_conditions_present"]
    if "sample_frequencies" in clonotypes.columns:
        return clonotypes["sample_frequencies"].apply(
            lambda d: sum(1 for v in _parse_sample_frequencies(d).values() if v > 0)
        )
    return None


def plot_clones_by_n_methods(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    filename: str = "clones_by_n_methods.png",
    highlight_tiers: tuple[str, ...] = ("tier1", "tier2"),
) -> Path | None:
    """1-D distribution of clones × number of distinct methods present (#43 #6).

    Two side-by-side bars per bin:
    - "all surviving clones" — every clone in the input frame
    - tier-selected subset (default tier1 + tier2) when a ``tier``
      column is present

    Uses a symlog y-axis so a single highly-public clone and the tail
    of singletons both stay visible. Returns ``None`` when neither
    ``n_methods`` nor ``sample_frequencies`` is available.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_methods = _derive_n_methods(clonotypes)
    if n_methods is None:
        logger.info(
            "plot_clones_by_n_methods: no n_methods / sample_frequencies; skipping."
        )
        return None
    n_methods = n_methods.dropna().astype(int)
    if len(n_methods) == 0:
        return None

    max_m = int(n_methods.max())
    bins = list(range(1, max_m + 1))
    surviving = n_methods.value_counts().reindex(bins, fill_value=0)

    if "tier" in clonotypes.columns:
        tier_mask = clonotypes["tier"].isin(highlight_tiers)
        tier_subset = _derive_n_methods(clonotypes[tier_mask])
        tier_subset = (
            tier_subset.dropna().astype(int).value_counts().reindex(bins, fill_value=0)
            if tier_subset is not None
            else pd.Series(0, index=bins)
        )
    else:
        tier_subset = None

    fig, ax = plt.subplots(figsize=(max(6, 0.7 * max_m + 3), 5))
    x = np.arange(len(bins))
    if tier_subset is not None:
        ax.bar(x - 0.2, surviving.values, width=0.38, label="all surviving clones",
               color="#1f77b4")
        ax.bar(x + 0.2, tier_subset.values, width=0.38,
               label=f"{'+'.join(highlight_tiers)}", color="#d62728")
        for xi, v in zip(x - 0.2, surviving.values):
            if v > 0:
                ax.text(xi, v, f"{int(v)}", ha="center", va="bottom", fontsize=8)
        for xi, v in zip(x + 0.2, tier_subset.values):
            if v > 0:
                ax.text(xi, v, f"{int(v)}", ha="center", va="bottom",
                        fontsize=8, color="#7f1d1d")
    else:
        ax.bar(x, surviving.values, width=0.6, color="#1f77b4",
               label="all surviving clones")
        for xi, v in zip(x, surviving.values):
            if v > 0:
                ax.text(xi, v, f"{int(v)}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_xlabel("Distinct samples / enrichment methods clone appears in")
    ax.set_ylabel("Number of clones")
    ax.set_yscale("symlog", linthresh=10)
    ax.set_title("Clones × number of distinct samples")
    ax.legend()

    out_path = output_dir / filename
    return save_figure(fig, out_path)


# =============================================================================
# Signature scatter — #43 item #7
# =============================================================================
#
# Gene-set defaults live in :mod:`tcrsift.signatures` (top-of-file
# import above). The per-sample scatter takes any ``gene_ids`` tuple,
# so callers can drop in custom signatures via that argument.


def plot_clone_freq_vs_signature_per_sample(
    adata: ad.AnnData,
    clonotypes: pd.DataFrame,
    gene_ids: list[str] | tuple[str, ...],
    signature_label: str,
    output_dir: str | Path,
    filename: str | None = None,
    min_cells_per_clone: int = 2,
    layout: list[list[str]] | None = None,
) -> Path | None:
    """Per-sample scatter of clone freq vs. mean log1p signature expression.

    For each sample, plots one point per (sample, clone) pair where
    ``x = within-sample frequency`` and ``y = mean log1p expression``
    of the genes in ``gene_ids`` averaged across that clone's cells
    in the sample. Spearman correlation between log-freq and signature
    rendered as a corner annotation per panel (#43 item #7).

    Requires:
    - ``adata`` with ``obs["sample"]``, ``obs["CDR3ab"]``
    - ``clonotypes`` with ``sample_frequencies`` (so we can look up
      per-(sample, clone) frequency) and ``cell_count`` for point
      sizing.

    Skips silently when any prerequisite is missing or when the gene
    list doesn't overlap ``adata.var_names`` at all.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    needed = {"sample", "CDR3ab"}
    if not needed.issubset(set(adata.obs.columns)):
        logger.info(
            "plot_clone_freq_vs_signature_per_sample: adata.obs missing "
            f"required columns ({needed - set(adata.obs.columns)}); skipping."
        )
        return None
    if "sample_frequencies" not in clonotypes.columns:
        logger.info(
            "plot_clone_freq_vs_signature_per_sample: no sample_frequencies; skipping."
        )
        return None

    gene_ids = list(gene_ids)
    valid_genes = [g for g in gene_ids if g in adata.var_names]
    if not valid_genes:
        logger.info(
            "plot_clone_freq_vs_signature_per_sample: none of the requested "
            f"genes ({gene_ids[:5]}...) match adata.var_names; skipping."
        )
        return None

    sig = _per_cell_signature(adata, gene_ids)
    cell_df = pd.DataFrame(
        {
            "sample": adata.obs["sample"].astype(str).values,
            "CDR3ab": adata.obs["CDR3ab"].astype(str).values,
            "signature": sig,
        }
    )
    cell_df = cell_df[cell_df["CDR3ab"] != "nan"]
    clone_sig = (
        cell_df.groupby(["sample", "CDR3ab"], observed=True)["signature"]
        .agg(["mean", "size"])
        .rename(columns={"size": "cells_in_sample"})
        .reset_index()
    )
    if clone_sig.empty:
        return None

    # Long-format frequency lookup from sample_frequencies.
    freq_rows = []
    for _, row in clonotypes.iterrows():
        freqs = _parse_sample_frequencies(row.get("sample_frequencies"))
        for sample, freq in freqs.items():
            if freq > 0:
                freq_rows.append(
                    {
                        "sample": str(sample),
                        "CDR3ab": str(row["CDR3ab"]),
                        "frequency": float(freq),
                    }
                )
    if not freq_rows:
        return None
    freq_long = pd.DataFrame(freq_rows)
    merged = clone_sig.merge(freq_long, on=["sample", "CDR3ab"], how="inner")
    merged = merged[merged["cells_in_sample"] >= min_cells_per_clone]
    if merged.empty:
        return None

    samples = sorted(merged["sample"].unique())
    fig, axes_map = make_method_panel_grid(
        samples, panel_size=(5.0, 3.8), layout=layout
    )

    for sample, ax in axes_map.items():
        m = merged[merged["sample"] == sample]
        if len(m) == 0:
            ax.axis("off")
            continue
        ax.scatter(
            m["frequency"],
            m["mean"],
            s=np.clip(m["cells_in_sample"], 4, 80),
            alpha=0.55,
            color="#1f77b4",
            edgecolor="white",
            linewidth=0.4,
        )
        if len(m) >= 5 and m["frequency"].std() > 0:
            # Spearman = Pearson of ranks (rank-based, so no log transform needed
            # and monotonic-robust). Previously labeled "Spearman" but computed
            # Pearson on log-frequency — now the label matches the statistic.
            r = np.corrcoef(
                m["frequency"].rank().to_numpy(), m["mean"].rank().to_numpy()
            )[0, 1]
            ax.text(
                0.02,
                0.97,
                f"Spearman freq vs sig: r = {r:.2f}\nn clones = {len(m)}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8,
                bbox=dict(
                    boxstyle="round,pad=0.3",
                    facecolor="white",
                    alpha=0.7,
                    edgecolor="none",
                ),
            )
        ax.set_xscale("log")
        ax.set_xlabel("clone frequency within sample")
        ax.set_ylabel(f"mean log1p {signature_label}")

        ax.set_title(pretty_sample(sample))
        ax.grid(True, which="both", linewidth=0.3, alpha=0.4)

    fig.suptitle(
        f"Clone frequency vs {signature_label} signature per sample",
        fontsize=12,
        y=1.01,
    )
    fig.tight_layout()
    filename = filename or f"clone_freq_vs_{signature_label.split()[0]}.png"
    out_path = output_dir / filename
    return save_figure(fig, out_path)

