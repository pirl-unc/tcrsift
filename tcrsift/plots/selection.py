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

"""TIL, method-overlap, PRISM, cross-donor, V-gene, leader, and assembly plots."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..format import pretty_method, pretty_methods
from .state import _category_color_strip, save_figure

logger = logging.getLogger(__name__)


# =============================================================================
# TIL Plots
# =============================================================================


def plot_til(matched_clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate TIL matching plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if "til_match" not in matched_clonotypes.columns:
        logger.warning("No TIL match data found")
        return

    from .._schema import truthy_flag_series

    til_mask = truthy_flag_series(matched_clonotypes["til_match"])

    # TIL recovery by tier
    if "tier" in matched_clonotypes.columns:
        tiers = sorted(matched_clonotypes["tier"].dropna().unique())
        if tiers:
            fig, ax = plt.subplots(figsize=(10, 6))
            recovery_rates = []

            for tier in tiers:
                tier_sel = matched_clonotypes["tier"] == tier
                n_tier = int(tier_sel.sum())
                recovery = (
                    float(til_mask[tier_sel].sum()) / n_tier * 100 if n_tier > 0 else 0.0
                )
                recovery_rates.append(recovery)

            ax.bar(
                range(len(tiers)), recovery_rates, color=plt.cm.viridis(np.linspace(0, 1, len(tiers)))
            )
            ax.set_xticks(range(len(tiers)))
            ax.set_xticklabels(tiers)
            ax.set_ylabel("TIL Recovery Rate (%)")
            ax.set_title("Culture→TIL Recovery by Confidence Tier")
            for i, v in enumerate(recovery_rates):
                ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=10)
            save_figure(fig, output_dir / "til_recovery_by_tier.png")

    # Culture vs TIL frequency scatter
    matched = matched_clonotypes.loc[til_mask]
    if (
        len(matched) > 0
        and "max_frequency" in matched.columns
        and "til_frequency" in matched.columns
    ):
        x_vals = pd.to_numeric(matched["max_frequency"], errors="coerce") * 100
        y_vals = pd.to_numeric(matched["til_frequency"], errors="coerce") * 100
        pos_mask = (x_vals > 0) & (y_vals > 0) & np.isfinite(x_vals) & np.isfinite(y_vals)
        if pos_mask.any():
            fig, ax = plt.subplots(figsize=(10, 10))
            sizes = (
                pd.to_numeric(matched.loc[pos_mask, "cell_count"], errors="coerce").fillna(4.0) * 10
                if "cell_count" in matched.columns
                else 40
            )
            ax.scatter(
                x_vals[pos_mask],
                y_vals[pos_mask],
                alpha=0.6,
                s=sizes,
            )
            ax.set_xlabel("Culture Frequency (%)")
            ax.set_ylabel("TIL Frequency (%)")
            ax.set_title("Culture vs TIL Frequency\n(size = culture cell count)")

            # Add positive-domain 1:1 diagonal for log-log axes
            min_pos = float(min(x_vals[pos_mask].min(), y_vals[pos_mask].min()))
            max_val = float(max(x_vals[pos_mask].max(), y_vals[pos_mask].max()))
            ax.plot([min_pos, max_val], [min_pos, max_val], "k--", alpha=0.3, label="1:1 line")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend()
            save_figure(fig, output_dir / "til_frequency_scatter.png")


# =============================================================================
# Method-stratified plots (#27)
# =============================================================================


def plot_set_overlap(
    sets: dict[str, set],
    output_path: str | Path,
    *,
    title: str = "Selection overlap",
    min_subset_size: int = 1,
    max_bars: int = 30,
):
    """N-way set-overlap (UpSet) plot over selected-clone sets (#208).

    ``sets`` is ``{set_name: {clones}}`` (e.g. from
    :func:`tcrsift.clonotype.build_selection_sets`, keyed by method / condition
    / donor). Renders an UpSet plot when the optional ``upsetplot`` dependency
    is installed; otherwise falls back to a bar chart of the largest
    intersection patterns so the plot still renders without the extra dep.

    No-op (writes nothing) when fewer than two non-empty sets are given.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sets = {k: set(v) for k, v in (sets or {}).items() if v}
    if len(sets) < 2:
        return

    # Prettify set names (AIMpos -> AIM⁺) and order them consistently so a
    # baseline marker (CTY⁻) always reads last — matching every other figure.
    from ..format import order_conditions

    ordered_names = order_conditions(sets)
    pretty = {name: pretty_method(name) for name in ordered_names}
    # Guard against two raw names collapsing to one pretty label.
    if len(set(pretty.values())) == len(pretty):
        sets = {pretty[name]: sets[name] for name in ordered_names}

    upset_fig = None
    try:
        import warnings

        import upsetplot
        from pandas.errors import PerformanceWarning

        # upsetplot ≤0.9 uses deprecated pandas idioms internally (chained
        # ``fillna(..., inplace=True)``, object-dtype downcasting, fragmenting
        # inserts → FutureWarning/PerformanceWarning) AND builds the figure with
        # array-valued text positions that warn at DRAW time (NumPy-1.25
        # "ndim>0 to scalar" DeprecationWarning, fired in save_figure). None are
        # ours to fix. Scope-suppress around the library calls AND the savefig
        # (the draw is where the DeprecationWarning fires); our own code doesn't
        # run inside this block.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            warnings.simplefilter("ignore", PerformanceWarning)
            warnings.simplefilter("ignore", DeprecationWarning)
            data = upsetplot.from_contents(sets)
            upset_fig = plt.figure(figsize=(max(7, 0.6 * len(sets) + 5), 5))
            upsetplot.plot(
                data, fig=upset_fig, sort_by="cardinality", show_counts=True,
                min_subset_size=min_subset_size,
            )
            upset_fig.suptitle(title)
            save_figure(upset_fig, output_path, dpi=150)
        return
    except ImportError:
        fallback_note = " (install 'upsetplot' for UpSet view)"
    except Exception as error:  # noqa: BLE001 - any failure inside the optional library
        # upsetplot 0.9.0 fills missing dot styles with chained in-place
        # ``fillna``, which pandas 3 ignores, so its matrix reaches matplotlib
        # with NaN colors. Draw the bar chart instead of failing the report.
        if upset_fig is not None:
            plt.close(upset_fig)
        logger.warning("UpSet plot failed (%s: %s); drawing a bar chart of overlaps instead",
                       type(error).__name__, error)
        fallback_note = ""

    # Fallback: bar chart of the top intersection patterns.
    from ..clonotype import set_overlap_table

    table = set_overlap_table(sets)
    table = table[table["n_clones"] >= min_subset_size].head(max_bars)
    if table.empty:
        return
    # Labels are already prettified+ordered above; just space the ∩ for reading.
    labels = [s.replace("+", " ∩ ") for s in table["sets"].values]
    fig, ax = plt.subplots(figsize=(max(7, 0.4 * len(table) + 3), 5))
    ax.bar(range(len(table)), table["n_clones"].values, color="#4c72b0")
    ax.set_xticks(range(len(table)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("clones")
    ax.set_title(f"{title}{fallback_note}")
    for i, n in enumerate(table["n_clones"].values):
        ax.text(i, n, str(int(n)), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    save_figure(fig, output_path, dpi=150)


def plot_method_overlap(
    matrix: pd.DataFrame,
    output_path: str | Path,
    similarity: str = "jaccard",
    donor: str | None = None,
):
    """Heatmap of method × method overlap for a single donor.

    ``matrix`` is a square DataFrame from ``build_method_overlap_matrices``.
    For ``similarity == 'count'`` cell labels show integer intersection
    counts; for ``jaccard`` / ``dice`` they show two-decimal floats.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(matrix)
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * n + 2), max(5, 0.7 * n + 1)))
    vals = np.asarray(matrix.values, dtype=float)
    finite_vals = vals[np.isfinite(vals)]
    vmax = float(finite_vals.max()) if finite_vals.size else 1.0
    is_int_matrix = (
        similarity == "count"
        and not np.isnan(vals).any()
        and all(np.issubdtype(dt, np.integer) for dt in matrix.dtypes)
    )
    fmt = "d" if is_int_matrix else (".0f" if similarity == "count" else ".2f")
    sns.heatmap(
        matrix,
        annot=True,
        fmt=fmt,
        cmap="viridis",
        ax=ax,
        cbar_kws={"label": similarity},
        square=True,
        vmin=0,
        vmax=vmax if similarity == "count" else 1.0,
    )
    title = f"Method × method overlap ({similarity})"
    if donor:
        title += f" — donor {donor}"
    ax.set_title(title)
    ax.set_xlabel("method")
    ax.set_ylabel("method")
    ax.set_xticklabels(pretty_methods(matrix.columns))
    ax.set_yticklabels(pretty_methods(matrix.index))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=150)


def plot_freq_prism_grid(
    grid: pd.DataFrame,
    output_path: str | Path,
    *,
    chosen: tuple[int, int] | None = None,
    title: str = "Distinct clones selected: frequency × PRISM",
):
    """Heatmap of total distinct clones over a ``top_freq`` × ``top_prism`` grid (#207).

    ``grid`` is the tidy DataFrame from
    :func:`tcrsift.selection.freq_prism_grid` (columns ``top_freq``,
    ``top_prism``, ``n_clones``). Rows are ``top_prism`` (PRISM picks per
    condition), columns are ``top_freq`` (frequency picks per condition); each
    cell is annotated with the distinct-clone count. Pass ``chosen=(top_freq,
    top_prism)`` to box the operating point (e.g. the ``(10, 5)`` default).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if grid is None or grid.empty:
        return

    pivot = grid.pivot(index="top_prism", columns="top_freq", values="n_clones")
    # Ascending so the origin (0, 0) sits at the lower-left like the downstream.
    pivot = pivot.sort_index(ascending=True).sort_index(axis=1, ascending=True)

    n_rows, n_cols = pivot.shape
    fig, ax = plt.subplots(
        figsize=(max(6, 0.8 * n_cols + 2), max(5, 0.8 * n_rows + 1))
    )
    pivot_int = (
        not pivot.isna().any().any()
        and all(np.issubdtype(dt, np.integer) for dt in pivot.dtypes)
    )
    sns.heatmap(
        pivot,
        annot=True,
        fmt="d" if pivot_int else ".0f",
        cmap="viridis",
        ax=ax,
        cbar_kws={"label": "distinct clones"},
        square=True,
    )
    ax.invert_yaxis()  # PRISM increases upward
    ax.set_title(title)
    ax.set_xlabel("top clones by frequency (per condition)")
    ax.set_ylabel("top clones by PRISM (per condition)")

    if chosen is not None:
        cf, cp = int(chosen[0]), int(chosen[1])
        if cf in list(pivot.columns) and cp in list(pivot.index):
            col_i = list(pivot.columns).index(cf)
            row_i = list(pivot.index).index(cp)
            ax.add_patch(
                plt.Rectangle(
                    (col_i, row_i), 1, 1, fill=False,
                    edgecolor="red", lw=3, clip_on=False,
                )
            )

    fig.tight_layout()
    save_figure(fig, output_path, dpi=150)


def plot_freq_prism_scatter(
    candidates: pd.DataFrame,
    output_path: str | Path,
    *,
    condition_col: str = "method",
    gate: float | None = None,
    title: str = "Selection: frequency × PRISM",
) -> Path | None:
    """Faceted scatter of frequency × PRISM, colored by selection route (#248).

    ``candidates`` is the per-clone selection table (see
    :func:`tcrsift.selection.prism_candidates`) with columns ``CDR3ab``,
    ``<condition_col>``, ``frequency``, ``prism_score`` (lower = better; may be
    NaN), and ``selection_route`` (``freq``/``prism``/``both``/``unselected``).
    One small-multiple subplot per condition: ``unselected`` clones in light
    grey, selected routes overlaid in distinct colors, frequency on a log
    x-axis, with the ``gate`` drawn as a vertical dashed red line. Returns the
    saved path, or None when there's nothing to plot.
    """
    output_path = Path(output_path)
    if candidates is None or candidates.empty:
        return None
    route_colors = {"freq": "tab:blue", "prism": "tab:orange", "both": "tab:green"}
    conditions = sorted(candidates[condition_col].dropna().unique(), key=str)
    if not conditions:
        return None
    ncols = min(3, len(conditions))
    nrows = int(np.ceil(len(conditions) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5.0 * ncols, 4.0 * nrows), squeeze=False,
    )
    flat_axes = axes.flatten()
    for ax, condition in zip(flat_axes, conditions):
        sub = candidates[candidates[condition_col] == condition]
        sub = sub[sub["prism_score"].notna()]
        bg = sub[sub["selection_route"] == "unselected"]
        ax.scatter(
            bg["frequency"], bg["prism_score"], s=12, color="lightgrey",
            alpha=0.6, linewidths=0, zorder=1,
        )
        for route, color in route_colors.items():
            sel = sub[sub["selection_route"] == route]
            ax.scatter(
                sel["frequency"], sel["prism_score"], s=45, color=color,
                edgecolors="black", linewidths=0.4, label=route, zorder=2,
            )
        if gate is not None:
            ax.axvline(gate, color="red", linestyle="--", linewidth=1)
        ax.set_xscale("log")
        ax.set_xlabel("clone frequency (within condition)")
        ax.set_ylabel("PRISM score (lower = better)")
        # Invert so better-PRISM clones rise to the TOP (#259): the tick labels
        # still read "lower = better", but now both axes agree that up/right =
        # stronger candidate, so the selected prism/both points cluster in one
        # corner. invert_yaxis() is idempotent per-axis (one call each).
        ax.invert_yaxis()
        ax.set_title(pretty_method(condition))
    for ax in flat_axes[len(conditions):]:
        ax.set_visible(False)
    handles = [
        plt.Line2D(
            [], [], marker="o", linestyle="", color=color,
            markeredgecolor="black", markeredgewidth=0.4, label=route,
        )
        for route, color in route_colors.items()
    ]
    # Legend on the first panel (not fig-level) so it can't overlap a subplot
    # in a dense grid.
    flat_axes[0].legend(handles=handles, fontsize="small", frameon=False)
    fig.suptitle(title)
    fig.tight_layout()
    return save_figure(fig, output_path)


def plot_prism_components(
    scores: pd.DataFrame,
    selected_ids,
    output_path: str | Path,
    *,
    clone_col: str = "CDR3ab",
    terms=("ppost_alpha", "ppost_beta", "antigen_response_score", "naive_score"),
    title: str = "PRISM components: selected vs background",
) -> Path | None:
    """Small-multiples of PRISM component distributions: selected vs background (#249).

    One histogram per present ``terms`` column — all gated candidates' values in
    grey, with the PRISM-selected subset (``selected_ids``) overlaid in orange —
    so it reads at a glance which components separate the picks from the
    background (and surfaces the length-confound / weighting concerns of #238).
    Absent term columns are skipped; NaNs dropped per term. Returns the saved
    path, or None when there's nothing to plot.
    """
    output_path = Path(output_path)
    if scores is None or scores.empty:
        return None
    present = [t for t in terms if t in scores.columns]
    if not present:
        return None
    selected_ids = set(selected_ids)
    is_selected = scores[clone_col].isin(selected_ids)
    n = len(present)
    ncols = min(n, 4)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(max(4, 3.5 * ncols), 3.5 * nrows), squeeze=False,
    )
    flat = axes.flatten()
    for i, term in enumerate(present):
        ax = flat[i]
        all_vals = scores[term].dropna()
        sel_vals = scores.loc[is_selected, term].dropna()
        if len(all_vals):
            # Shared bin edges (from the full candidate range) so the selected
            # overlay is bin-aligned with the background — the comparison this
            # plot exists to make.
            edges = np.histogram_bin_edges(all_vals, bins=40)
            ax.hist(all_vals, bins=edges, color="grey", label="candidates")
            if len(sel_vals):
                ax.hist(sel_vals, bins=edges, color="tab:orange", alpha=0.8,
                        label="PRISM-selected")
        ax.set_title(term)
        ax.set_xlabel("score")
        ax.set_ylabel("clones")
        if i == 0:
            ax.legend(fontsize="small", frameon=False)
    for ax in flat[n:]:
        ax.set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    return save_figure(fig, output_path)


def plot_cross_donor_venn(
    clonotypes_by_donor: dict,
    output_path: str | Path,
    *,
    keys: Sequence = (
        ("CDR3ab", "CDR3αβ (paired)"),
        ("CDR3_beta", "CDR3β only"),
    ),
    title: str = "Cross-donor clonotype sharing",
    max_pairs: int = 15,
) -> Path | None:
    """Pairwise cross-donor clonotype-sharing Venn diagrams (#250).

    For every unordered pair of donors and each granularity in ``keys`` (paired
    αβ and β-only), draws a 2-way Venn of shared clonotypes — rows = keys,
    columns = donor-pairs. β-only sharing (the permissive "public TCR" view) is
    typically much larger than αβ. Uses the optional ``matplotlib_venn`` when
    installed, else falls back to a grouped |A only|/|shared|/|B only| bar so it
    still renders. No-op (returns None) for fewer than two donors.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    donors = list(clonotypes_by_donor or {})
    if len(donors) < 2:
        logger.info("plot_cross_donor_venn: need >= 2 donors; skipping.")
        return None
    sets_by_key: dict = {}
    for col, _label in keys:
        sets_by_key[col] = {
            d: (set(df[col].dropna().astype(str)) if col in df.columns else set())
            for d, df in clonotypes_by_donor.items()
        }
    pairs = list(combinations(donors, 2))
    if len(pairs) > max_pairs:
        # Avoid an unreadably wide figure for many donors; an UpSet over donors
        # is the better view at that scale (#250). Cap + log what was dropped.
        logger.info(
            "plot_cross_donor_venn: %d donor-pairs > max_pairs=%d; plotting the "
            "first %d (consider an UpSet over donors instead).",
            len(pairs), max_pairs, max_pairs,
        )
        pairs = pairs[:max_pairs]
    n_rows, n_cols = len(keys), len(pairs)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(max(4.0, 4.0 * n_cols), max(4.0, 4.0 * n_rows)),
        squeeze=False,
    )
    try:
        import warnings

        import matplotlib_venn
        from pandas.errors import PerformanceWarning

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            warnings.simplefilter("ignore", PerformanceWarning)
            for r, (col, label) in enumerate(keys):
                for c, (a, b) in enumerate(pairs):
                    matplotlib_venn.venn2(
                        [sets_by_key[col][a], sets_by_key[col][b]],
                        set_labels=(a, b), ax=axes[r][c],
                    )
                    axes[r][c].set_title(f"{label}: {a} vs {b}", fontsize=11)
    except ImportError:
        for r, (col, label) in enumerate(keys):
            for c, (a, b) in enumerate(pairs):
                ax = axes[r][c]
                sa, sb = sets_by_key[col][a], sets_by_key[col][b]
                heights = [len(sa - sb), len(sa & sb), len(sb - sa)]
                ax.bar(range(3), heights, color=["#4c72b0", "#55a868", "#c44e52"])
                ax.set_xticks(range(3))
                ax.set_xticklabels([f"{a} only", "shared", f"{b} only"],
                                   rotation=30, ha="right", fontsize=8)
                ax.set_ylabel("clones")
                for i, h in enumerate(heights):
                    ax.text(i, h, str(h), ha="center", va="bottom", fontsize=8)
                ax.set_title(
                    f"{label}: {a} vs {b}\n(install 'tcrsift[plots]' / 'matplotlib-venn' for Venn view)",
                    fontsize=9,
                )
    fig.suptitle(title)
    fig.tight_layout()
    save_figure(fig, output_path, dpi=150)
    return output_path


def plot_vgene_usage_by_method(
    long_df: pd.DataFrame,
    clonotypes: pd.DataFrame,
    output_path: str | Path,
    *,
    vgene_cols=(("alpha_v_gene", "TRAV"), ("beta_v_gene", "TRBV")),
    method_col: str = "method",
    clone_col: str = "CDR3ab",
    min_prevalence: float = 0.01,
    title: str = "V-gene usage by method",
) -> Path | None:
    """Per-method V-gene usage heatmap(s), one subplot per V-gene column (#251).

    Merges ``long_df`` (``CDR3ab``, ``method``, ``frequency``) with
    ``clonotypes`` V-gene columns, builds a method × V-gene matrix
    **row-normalized** per method (fraction of that method's clones using each
    gene), applies a ``min_prevalence`` floor, and orders columns by descending
    usage. Skips absent/empty V-gene columns; returns None when nothing renders.
    """
    output_path = Path(output_path)
    if long_df is None or clonotypes is None or long_df.empty or clonotypes.empty:
        return None
    if (
        method_col not in long_df.columns
        or clone_col not in long_df.columns
        or clone_col not in clonotypes.columns
    ):
        return None
    present = [(c, p) for c, p in vgene_cols if c in clonotypes.columns]
    if not present:
        return None
    matrices = []
    long_dedup = long_df[[method_col, clone_col]].drop_duplicates()
    for vgcol, pretty in present:
        clo_vg = clonotypes[[clone_col, vgcol]].dropna(subset=[vgcol]).drop_duplicates(subset=[clone_col])
        merged = long_dedup.merge(clo_vg, on=clone_col, how="inner")
        if merged.empty:
            continue
        counts = merged.groupby([method_col, vgcol]).size().unstack(fill_value=0)
        if counts.empty or counts.shape[1] == 0:
            continue
        row_sums = counts.sum(axis=1).replace(0, np.nan)
        frac = counts.div(row_sums, axis=0).fillna(0.0)
        keep = frac.columns[frac.max(axis=0) >= min_prevalence]
        frac = frac[keep]
        if frac.empty or frac.shape[1] == 0:
            continue
        col_order = frac.sum(axis=0).sort_values(ascending=False).index
        matrices.append((frac[col_order], pretty))
    if not matrices:
        return None
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(matrices)
    widths = [max(6, 0.45 * m.shape[1] + 3) for m, _ in matrices]
    height = max(4, 0.5 * max(m.shape[0] for m, _ in matrices) + 2)
    fig, axes = plt.subplots(1, n, figsize=(sum(widths), height), squeeze=False)
    for ax, (frac, pretty) in zip(axes[0], matrices):
        sns.heatmap(
            frac, cmap="viridis", linewidths=0.3, ax=ax,
            cbar_kws={"label": "fraction of method's clones"},
        )
        ax.set_yticklabels([pretty_method(m) for m in frac.index], rotation=0)
        ax.set_title(f"{pretty} usage by method")
        ax.set_xlabel("V gene")
        ax.set_ylabel("method")
    fig.suptitle(title)
    fig.tight_layout()
    return save_figure(fig, output_path)


_FREQ_HEATMAP_PAGE_ROWS = 60  # clones per page — keeps every y-label legible


# Heatmap geometry (inches). Cells are deliberately large and the CDR3αβ label
# font small — the prior clustermap squeezed the cells to a sliver while the long
# labels and dendrogram gutters dominated. Laying axes out in absolute inches
# (not clustermap's auto-sizing) is what keeps the cells big and the condition
# labels from overlapping regardless of clone count.
_FREQ_CELL_W = 0.62   # per condition column
_FREQ_CELL_H = 0.24   # per clone row
_FREQ_LABEL_W = 3.0   # right-hand gutter for CDR3αβ labels
_FREQ_STRIP_W = 0.16  # public-DB colour strip
_FREQ_CBAR_W = 0.18
_FREQ_YFONT = 5.5     # CDR3αβ labels — small, they are long
_FREQ_XFONT = 9


def _render_frequency_heatmap_page(
    data, output_path, *, title, annotations, page, n_pages, vmin, vmax
):
    """Render one page of the selected-frequency heatmap (#selected-freq).

    Explicit inch-based axes layout: a left colour-bar (log scale), an optional
    public-DB colour strip, the heatmap with big cells, and CDR3αβ labels in a
    right-hand gutter. Frequencies are log-scaled (LogNorm) so the low-frequency
    bulk is distinguishable instead of washing out under a few dominant clones.
    """
    import matplotlib.patches as mpatches
    from matplotlib.colors import LogNorm, to_rgb

    n_rows, n_cols = data.shape

    # Public-DB strip: only when a real match is present on this page (an all-grey
    # "none" strip is just noise).
    strip_colors = None
    legend_handles = None
    if annotations is not None:
        cats = annotations.reindex(data.index)
        colors = _category_color_strip(cats)
        seen: dict[str, str] = {}
        for v, c in zip(cats.fillna("none").astype(str), colors):
            seen.setdefault(v, c)
        if set(seen) - {"none"}:
            strip_colors = colors
            legend_handles = [
                mpatches.Patch(color=c, label=k)
                for k, c in seen.items() if k != "none"
            ]
    has_strip = strip_colors is not None

    # Absolute inch layout → fractions, so cells stay a fixed size no matter how
    # many clones/conditions a page has. A minimum figure height keeps sparse
    # pages from collapsing into a tiny figure where the (absolute-point) title /
    # colourbar fonts look comically large; the heatmap block is centred in the
    # slack rather than stretching the cells.
    heat_w = n_cols * _FREQ_CELL_W
    heat_h = n_rows * _FREQ_CELL_H
    left = 0.55          # left margin → colourbar
    gap = 0.7            # colourbar → strip/heatmap (room for the rotated label)
    strip_w = _FREQ_STRIP_W if has_strip else 0.0
    strip_gap = 0.06 if has_strip else 0.0
    bottom = 1.25        # x-tick labels + xlabel
    top = 0.7            # title
    fig_w = left + _FREQ_CBAR_W + gap + strip_w + strip_gap + heat_w + _FREQ_LABEL_W
    fig_h = max(4.5, bottom + heat_h + top)
    y0 = bottom + (fig_h - bottom - top - heat_h) / 2  # centre the block

    fig = plt.figure(figsize=(fig_w, fig_h))

    def _ax(x, y, w, h):
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])

    cbar_ax = _ax(left, y0, _FREQ_CBAR_W, heat_h)
    hx = left + _FREQ_CBAR_W + gap
    if has_strip:
        strip_ax = _ax(hx, y0, strip_w, heat_h)
        rgb = np.array([to_rgb(c) for c in strip_colors]).reshape((n_rows, 1, 3))
        strip_ax.imshow(rgb, aspect="auto", origin="upper")
        strip_ax.set_xticks([])
        strip_ax.set_yticks([])
        strip_ax.set_title("DB", fontsize=7, pad=3)
        hx += strip_w + strip_gap
    heat_ax = _ax(hx, y0, heat_w, heat_h)

    # Log colour scale: mask zeros (absent in a condition) so they read as blank
    # rather than being forced to the floor of the log scale.
    masked = data.mask(data <= 0)
    norm = LogNorm(vmin=vmin, vmax=vmax) if vmin and vmax and vmax > vmin else None
    sns.heatmap(
        masked, ax=heat_ax, cbar_ax=cbar_ax, cmap="rocket_r", norm=norm,
        vmin=None if norm else vmin, vmax=None if norm else vmax,
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "within-condition frequency (%, log)"},
    )
    cbar_ax.yaxis.label.set_size(8)
    cbar_ax.tick_params(labelsize=7)
    # CDR3αβ labels in the right gutter, small font.
    heat_ax.yaxis.set_label_position("right")
    heat_ax.yaxis.tick_right()
    heat_ax.set_yticks(np.arange(n_rows) + 0.5)
    heat_ax.set_yticklabels(
        list(data.index), rotation=0, fontsize=_FREQ_YFONT, fontfamily="monospace",
    )
    heat_ax.set_ylabel("selected clone (CDR3αβ)", fontsize=9)
    heat_ax.set_xticks(np.arange(n_cols) + 0.5)
    heat_ax.set_xticklabels(
        list(data.columns), rotation=45, ha="right", fontsize=_FREQ_XFONT,
    )
    heat_ax.set_xlabel("condition", fontsize=10)
    heat_ax.tick_params(length=2)

    # Left-align the title at the heatmap's left edge (right of the colourbar) so
    # it never collides with the public-DB legend in the top-left corner.
    suffix = f" (page {page}/{n_pages})" if n_pages > 1 else ""
    fig.text(
        hx / fig_w, 1 - 0.35 / fig_h, ha="left", va="top", fontsize=11,
        s=f"{title} — selected-clone frequency by condition{suffix}",
    )
    if legend_handles:
        fig.legend(
            handles=legend_handles, title="public-DB match",
            loc="upper left", bbox_to_anchor=(0.01, 1 - 0.04 / fig_h),
            fontsize=7, title_fontsize=8, frameon=True, framealpha=0.9,
        )
    save_figure(fig, output_path, dpi=200)
    return output_path


def plot_selected_frequency_heatmap(
    pivot, output_path, *, title: str = "selected", annotations=None
):
    """Heatmap of the selected clones × conditions frequency table (#selected-freq).

    ``pivot`` is rows = CDR3ab, columns = conditions, cells = within-condition
    frequency (%). When ``annotations`` (a Series keyed by CDR3ab giving each
    clone's public-DB ``db_category`` — ``viral`` / ``tumor_self`` / … / NaN) is
    supplied, a left-side colour strip + legend marks the public-DB match so
    viral/public TCR-T candidates are readable straight off the heatmap.

    Rows are sorted by total frequency (descending) and **paginated** at
    ~60 clones/page so every clone's y-label stays legible — no row is ever
    label-less. The first page is written to ``output_path``; further pages get a
    ``_p2`` / ``_p3`` … stem suffix. Returns the page-1 path (a list of all page
    paths is available via the side effect of the written files), or None on
    empty input.
    """
    output_path = Path(output_path)
    if pivot is None or getattr(pivot, "empty", True):
        return None
    data = pivot.fillna(0.0)
    # Sort by total frequency so the strongest clones lead page 1.
    order = data.sum(axis=1).sort_values(ascending=False, kind="stable").index
    data = data.loc[order]

    n_rows = data.shape[0]
    # Balance pages evenly (62 clones → 31 + 31, never 60 + 2) so no page is a
    # tiny stub with comically large relative fonts.
    n_pages = max(1, (n_rows + _FREQ_HEATMAP_PAGE_ROWS - 1) // _FREQ_HEATMAP_PAGE_ROWS)
    page_rows = (n_rows + n_pages - 1) // n_pages
    # Shared LOG color scale across all pages so a cell's shade means the same
    # frequency everywhere (per-page autoscaling would make page 2's max look as
    # hot as page 1's). vmin = smallest positive frequency (log can't take 0).
    flat = data.to_numpy().ravel()
    positive = flat[flat > 0]
    vmin = float(positive.min()) if positive.size else None
    vmax = float(positive.max()) if positive.size else None

    paths = []
    for page in range(1, n_pages + 1):
        chunk = data.iloc[(page - 1) * page_rows : page * page_rows]
        if page == 1:
            page_path = output_path
        else:
            page_path = output_path.with_name(
                f"{output_path.stem}_p{page}{output_path.suffix}"
            )
        paths.append(
            _render_frequency_heatmap_page(
                chunk, page_path, title=title, annotations=annotations,
                page=page, n_pages=n_pages, vmin=vmin, vmax=vmax,
            )
        )
    return paths[0]


def plot_leader_summary(
    df: pd.DataFrame,
    output_path: str | Path,
    *,
    chains=("alpha", "beta"),
    top_n: int = 20,
    typical_range: tuple[int, int] = (15, 25),
    title: str = "Signal-peptide (leader) summary",
) -> Path | None:
    """Leader-sequence frequency + α/β leader length distribution (#262).

    Two views per chain that make ``from_contig`` signal peptides legible: the
    top-``top_n`` distinct leader peptides by construct count (dominant SPs vs the
    one-off tail that's often a mis-extraction), and the leader length histogram
    with the typical signal-peptide window (~15-25 aa) shaded so over-long /
    too-short outliers pop. ``df`` is ``full_sequences`` / ``selected_clones``
    (needs ``{chain}_leader_aa``). When a ``{chain}_leader_qc`` column is present
    the length bars are split SP-sound (``ok``-like) vs flagged, so QC problems
    read off the length axis. Returns the saved path, or None if no leader data.
    """
    output_path = Path(output_path)
    present = [c for c in chains if f"{c}_leader_aa" in getattr(df, "columns", [])]
    present = [c for c in present if df[f"{c}_leader_aa"].notna().any()]
    if df is None or df.empty or not present:
        return None

    # ok-like QC verdicts (sound SPs); anything else is a flagged leader.
    _ok = {"ok", "weak_kozak_start", "long_leader", "hregion_trimmed",
           "germline_reference", "substituted"}
    fig, axes = plt.subplots(2, len(present), figsize=(5.5 * len(present), 7), squeeze=False)
    for j, ch in enumerate(present):
        s = df[f"{ch}_leader_aa"].dropna().astype(str)
        # Row 0: top-N distinct leaders by count.
        vc = s.value_counts().head(top_n)
        ax = axes[0][j]
        ax.barh(range(len(vc)), vc.values, color="#4c78a8")
        ax.set_yticks(range(len(vc)))
        ax.set_yticklabels(
            [f"{x[:18]}{'…' if len(x) > 18 else ''} (len {len(x)})" for x in vc.index],
            fontsize=7,
        )
        ax.invert_yaxis()
        ax.set_title(f"{ch} leader frequency (top {min(top_n, len(vc))})")
        ax.set_xlabel("constructs")

        # Row 1: length histogram, split by QC when available.
        ax = axes[1][j]
        lengths = s.str.len()
        bins = range(0, max(40, int(lengths.max()) + 4), 2)
        qc_col = f"{ch}_leader_qc"
        if qc_col in df.columns:
            ok_mask = df[qc_col].isin(_ok)
            ok_len = df.loc[ok_mask, f"{ch}_leader_aa"].dropna().astype(str).str.len()
            bad_len = df.loc[~ok_mask, f"{ch}_leader_aa"].dropna().astype(str).str.len()
            ax.hist([ok_len, bad_len], bins=bins, stacked=True,
                    color=["#54a24b", "#e45756"], label=["SP-sound", "flagged"])
            ax.legend(fontsize=7, frameon=False)
        else:
            ax.hist(lengths, bins=bins, color="#4c78a8")
        ax.axvspan(*typical_range, color="green", alpha=0.10)
        ax.set_title(f"{ch} leader length")
        ax.set_xlabel("aa")
        ax.set_ylabel("constructs")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(title)
    fig.tight_layout()
    return save_figure(fig, output_path)


def plot_method_recovery(
    recovery: pd.DataFrame,
    output_path: str | Path,
    tier_label: str = "tier1",
):
    """Bar plot of per-method recovery of a target tier, paired by donor.

    ``recovery`` is the long-form DataFrame from
    ``build_method_recovery_table`` with columns
    ``[donor, method, recovered, total, fraction]``.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if recovery.empty:
        return

    # Sort methods by overall (across-donor) recovery fraction so the worst
    # / best methods read off the chart at a glance.
    method_order = (
        recovery.groupby("method", observed=True)["fraction"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )
    donor_order = sorted(recovery["donor"].astype(str).unique())

    n_methods = len(method_order)
    n_donors = len(donor_order)
    width = 0.8 / max(n_donors, 1)
    x = np.arange(n_methods)

    fig, ax = plt.subplots(figsize=(max(6, 0.6 * n_methods + 2), 4.5))
    palette = sns.color_palette("Set2", n_colors=max(n_donors, 1))
    for i, donor in enumerate(donor_order):
        sub = recovery[recovery["donor"].astype(str) == donor].set_index("method")
        sub = sub.reindex(method_order)
        offsets = (i - (n_donors - 1) / 2) * width
        ax.bar(
            x + offsets,
            sub["fraction"].fillna(0).values,
            width,
            label=str(donor),
            color=palette[i % len(palette)],
        )
        # Annotate "recovered/total" above each bar.
        for j, (_, row) in enumerate(sub.iterrows()):
            if pd.isna(row.get("total")):
                continue
            ax.text(
                x[j] + offsets,
                (row.get("fraction") or 0) + 0.02,
                f"{int(row['recovered'])}/{int(row['total'])}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(pretty_methods(method_order), rotation=30, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(f"fraction of {tier_label} clones recovered")
    ax.set_xlabel("method")
    ax.set_title(f"Method recovery of {tier_label} clones (per donor)")
    if n_donors > 1:
        ax.legend(title="donor", loc="best")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=150)


# =============================================================================
# Assembly Plots
# =============================================================================


def plot_assembly(clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate sequence assembly plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Chain length distributions
    for chain in ["alpha", "beta"]:
        col = f"full_{chain}_aa"
        if col in clonotypes.columns:
            lengths = clonotypes[col].dropna().str.len()
            if len(lengths) > 0:
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.hist(lengths, bins=30, edgecolor="black")
                ax.set_xlabel("Sequence Length (aa)")
                ax.set_ylabel("Number of Clonotypes")
                ax.set_title(f"{chain.upper()} Chain Length Distribution")
                ax.axvline(
                    lengths.median(),
                    color="red",
                    linestyle="--",
                    label=f"Median: {lengths.median():.0f}",
                )
                ax.legend()
                save_figure(fig, output_dir / f"assembly_{chain}_length.png")

    # CDR3 length distributions
    for chain in ["alpha", "beta"]:
        col = f"CDR3_{chain}"
        if col in clonotypes.columns:
            lengths = clonotypes[col].dropna().str.len()
            if len(lengths) > 0:
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.hist(lengths, bins=range(5, 30), edgecolor="black")
                ax.set_xlabel("CDR3 Length (aa)")
                ax.set_ylabel("Number of Clonotypes")
                ax.set_title(f"CDR3 {chain.upper()} Length Distribution")
                save_figure(fig, output_dir / f"assembly_cdr3_{chain}_length.png")

