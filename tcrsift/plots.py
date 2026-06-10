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
"""
Visualization functions for TCRsift.

Generates plots for QC, phenotyping, clonotype analysis, and filtering.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .format import pretty_method, pretty_methods, pretty_sample, pretty_samples  # noqa: F401, E402

# Re-export canonical T-cell signatures so callers of the per-sample
# scatter can keep importing them from ``tcrsift.plots``. Source of
# truth is :mod:`tcrsift.signatures`.
from .signatures import (  # noqa: F401, E402
    ACTIVATION_GENES_HGNC,
    ANTIGEN_RESPONSE_GENES_HGNC,
    CYTOLYTIC_GENES_HGNC,
    EXHAUSTION_GENES_HGNC,
    TUMOR_REACTIVE_GENES_HGNC,
)

logger = logging.getLogger(__name__)

# Set default style
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.facecolor"] = "#f8f9fa"

# Configure matplotlib to use a Unicode-complete sans-serif so the
# ``pretty_sample`` / ``pretty_method`` strings (which use
# U+207A SUPERSCRIPT PLUS SIGN and U+207B SUPERSCRIPT MINUS for the
# ``AIM⁺`` / ``CTY⁻`` style) render correctly (#109 category 2).
# Pre-fix the fallback was Arial, which lacks both glyphs, producing
# empty boxes plus "Glyph 8314/8315 missing from font(s) Arial"
# warnings on every render. DejaVu Sans is bundled with matplotlib,
# carries the full Unicode range we need, and is the canonical
# matplotlib default — so the only thing to do is make it explicit
# at the head of the family-fallback list so Arial doesn't get a
# chance to claim the render first.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "DejaVu Sans",
    *[
        f for f in plt.rcParams.get("font.sans-serif", [])
        if f != "DejaVu Sans"
    ],
]


VALID_PLOT_FORMATS = ("png", "pdf", "svg", "both")

# Module-level output format for save_figure (#169). Default "png" keeps the
# historical behavior byte-for-byte. cmd_run sets it from output.plot_format.
_PLOT_FORMAT = "png"


def set_plot_format(fmt: str | None) -> None:
    """Set the vector/raster format save_figure emits (#169): png/pdf/svg/both.

    A PNG is ALWAYS written (it's what the raster-embedding PDF report consumes);
    the format selects which vector copy is *also* emitted alongside it:
    ``pdf``→+PDF, ``svg``→+SVG, ``both``→+PDF +SVG (#258), ``png``→nothing extra.

    Process-global (not thread-safe): it stays in effect until set again, so a
    library caller invoking ``plot_*`` after a pdf ``run`` inherits ``pdf``.
    ``cmd_run`` sets it from ``output.plot_format`` at the start of every run.
    """
    global _PLOT_FORMAT
    fmt = (fmt or "png").lower()
    if fmt not in VALID_PLOT_FORMATS:
        raise ValueError(
            f"Invalid plot_format: {fmt!r}. Valid: {VALID_PLOT_FORMATS}"
        )
    _PLOT_FORMAT = fmt


def set_polished_style(style: str = "clean-white") -> None:
    """Apply a publication style profile to matplotlib (#123 report polished).

    ``clean-white`` (default): DejaVu Sans (renders the Unicode superscript
    method labels, #109) on a white background. ``paper`` additionally bumps the
    base font size. Process-global; pair with ``set_plot_format("pdf")`` for
    selectable-text vector output.
    """
    import matplotlib as mpl

    mpl.rcParams["font.family"] = "DejaVu Sans"
    for k in ("figure.facecolor", "axes.facecolor", "savefig.facecolor"):
        mpl.rcParams[k] = "white"
    if style == "paper":
        mpl.rcParams["font.size"] = 11


def save_figure(fig: plt.Figure, output_path: str | Path, dpi: int = 300):
    """Save figure with consistent settings.

    Honors the configured plot format (#169, #258). A PNG is always written —
    it's the default and what the (raster-embedding) PDF report consumes — and
    a vector copy is written alongside it: ``pdf``→+PDF, ``svg``→+SVG,
    ``both``→+PDF +SVG. With the default ``png`` format only the PNG is written,
    unchanged. Returns the primary output path (the first vector file when one
    was requested, else the PNG).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if _PLOT_FORMAT == "png":
        # Default / back-compat: honor the caller's exact path and suffix
        # (byte-identical to pre-#169 behavior, including explicit non-png
        # suffixes a caller passes deliberately).
        targets = [output_path]
    else:
        # Vector requested: emit the vector file(s) FIRST (so the returned
        # primary path is the vector), then always the PNG (keeps the
        # raster-embedding PDF report working). "both" → PDF + SVG.
        vector_exts = ["pdf", "svg"] if _PLOT_FORMAT == "both" else [_PLOT_FORMAT]
        targets = [output_path.with_suffix(f".{ext}") for ext in vector_exts]
        targets.append(output_path.with_suffix(".png"))
    for target in targets:
        fig.savefig(target, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    logger.info(f"Saved plot to {', '.join(str(t) for t in targets)}")
    return targets[0]


def plot_pgen_distribution(
    df: pd.DataFrame,
    output_path: str | Path,
    *,
    pgen_col: str = "log10_pgen",
    publicness_col: str | None = "publicness",
    group_col: str | None = None,
    bins: int = 40,
    title: str | None = None,
) -> Path | None:
    """Histogram of log10 Pgen with the publicness cutoff overlaid (#58).

    Parameters
    ----------
    df : pd.DataFrame
        Clonotypes with a ``log10_pgen`` (or named via ``pgen_col``)
        column produced by :func:`tcrsift.pgen.compute_pgen` or
        :func:`tcrsift.pgen.annotate_publicness`.
    output_path : str | Path
        File to write.
    pgen_col : str
        Column name with log10 Pgen.
    publicness_col : str | None
        If a publicness column is present, color bars by mean
        publicness in each bin (red = public, blue = private).
    group_col : str | None
        Optional categorical column to facet on (e.g. ``tier``). When
        set, overlays one translucent histogram per group.
    bins : int
        Number of histogram bins.
    title : str | None
        Optional title override.

    Returns
    -------
    Path | None
        Output path, or None if there's no data to plot.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if pgen_col not in df.columns:
        logger.info(
            f"plot_pgen_distribution: {pgen_col!r} not in frame; skipping."
        )
        return None
    values = df[pgen_col].dropna()
    if values.empty:
        logger.info("plot_pgen_distribution: no non-NaN Pgen values; skipping.")
        return None

    fig, ax = plt.subplots(figsize=(8, 4.5))

    if group_col and group_col in df.columns:
        palette = sns.color_palette("Set2", n_colors=df[group_col].nunique())
        for color, (g, sub) in zip(palette, df.groupby(group_col, observed=True)):
            ax.hist(
                sub[pgen_col].dropna(), bins=bins, alpha=0.55,
                label=str(g), color=color, edgecolor="white",
            )
        ax.legend(title=group_col, frameon=False, fontsize=9)
    else:
        ax.hist(values, bins=bins, color="#3b82f6", edgecolor="white", alpha=0.85)

    # Annotate the publicness cutoffs (matches the
    # ``publicness_score`` defaults — pass override lines via title if
    # the caller is using a different scale).
    for cutoff, label, color in (
        (-30.0, "low (private)", "#1d4ed8"),
        (-18.0, "high (public)", "#dc2626"),
    ):
        ax.axvline(cutoff, color=color, linestyle="--", linewidth=1.0, alpha=0.7)
        ax.text(
            cutoff, ax.get_ylim()[1] * 0.95, f" {label}",
            color=color, fontsize=8.5, va="top", ha="left",
        )

    ax.set_xlabel(r"log$_{10}$ Pgen (proxy)")
    ax.set_ylabel("clones")
    ax.set_title(title or "Generation-probability distribution")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    save_figure(fig, output_path)
    return output_path


def plot_publicness_vs_match_score(
    df: pd.DataFrame,
    output_path: str | Path,
    *,
    publicness_col: str = "publicness",
    match_score_col: str = "n_db_matches",
    tier_col: str | None = "tier",
    title: str | None = None,
) -> Path | None:
    """Scatter of publicness vs. raw DB match count (#58).

    Helps spot DB matches that are inflated by public sequences:
    points in the high-publicness, high-match-score corner are clones
    that look "well-annotated" only because they're easy to generate
    by chance.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if publicness_col not in df.columns or match_score_col not in df.columns:
        logger.info(
            f"plot_publicness_vs_match_score: missing "
            f"{publicness_col!r} or {match_score_col!r}; skipping."
        )
        return None
    sub = df.dropna(subset=[publicness_col, match_score_col])
    if sub.empty:
        return None

    fig, ax = plt.subplots(figsize=(7, 5))
    if tier_col and tier_col in sub.columns:
        palette = sns.color_palette("rocket_r", n_colors=sub[tier_col].nunique())
        for color, (t, g) in zip(palette, sub.groupby(tier_col, observed=True)):
            ax.scatter(
                g[publicness_col], g[match_score_col],
                color=color, alpha=0.7, s=22, edgecolor="white",
                linewidth=0.4, label=str(t),
            )
        ax.legend(title=tier_col, frameon=False, fontsize=9, loc="upper left")
    else:
        ax.scatter(
            sub[publicness_col], sub[match_score_col],
            color="#3b82f6", alpha=0.7, s=22, edgecolor="white",
            linewidth=0.4,
        )
    ax.set_xlabel("publicness (1 = generatable by chance)")
    ax.set_ylabel(match_score_col)
    ax.set_title(title or "DB match count vs. publicness")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_figure(fig, output_path)
    return output_path


def plot_assembly_qc(
    report,
    output_path: str | Path | None = None,
    *,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Stacked horizontal bar of pass/fail per assembly QC check (#67).

    Parameters
    ----------
    report : tcrsift.assemble.AssemblyQCReport
        Structured report from :func:`tcrsift.assemble.build_assembly_qc_report`.
    output_path : str | Path | None
        If given, save the figure here and return it after closing.
    figsize : tuple[float, float] | None
        Override default sizing. Default scales height to check count.
    """
    checks = list(report.checks)
    if not checks:
        fig, ax = plt.subplots(figsize=(6, 1.5))
        ax.text(
            0.5, 0.5, "No QC checks ran (empty input)",
            ha="center", va="center", transform=ax.transAxes,
        )
        ax.set_axis_off()
        if output_path:
            save_figure(fig, output_path)
        return fig

    if figsize is None:
        figsize = (9, max(2.5, 0.4 * len(checks) + 1.5))
    fig, ax = plt.subplots(figsize=figsize)

    labels = [c.label for c in checks]
    passed = np.array([c.passed for c in checks])
    failed = np.array([c.failed for c in checks])
    y = np.arange(len(checks))

    ax.barh(y, passed, color="#2a9d8f", label="pass")
    ax.barh(y, failed, left=passed, color="#e63946", label="fail")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("clones")
    status = "PASS" if report.passed else "FAIL"
    ax.set_title(f"Assembly QC — {report.n_rows} clones — {status}")
    ax.legend(loc="lower right", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate each row with the pass/total fraction (and median where set).
    for i, c in enumerate(checks):
        right_edge = c.total
        text_parts = [f"{c.passed}/{c.total}"]
        if c.median_value is not None:
            text_parts.append(
                f"med {c.median_value:g}"
                + (f" {c.unit}" if c.unit else "")
            )
        ax.text(
            right_edge,
            i,
            "  " + " · ".join(text_parts),
            va="center",
            fontsize=9,
            color="#444",
        )

    fig.tight_layout()
    if output_path:
        save_figure(fig, output_path)
    return fig


# =============================================================================
# QC Plots (load command)
# =============================================================================


def plot_cells_per_sample_stacked(
    adata: ad.AnnData,
    output_path: str | Path,
    *,
    qc_cols: tuple[str, ...] = (
        "filter:min_genes", "filter:max_genes",
        "filter:min_counts", "filter:max_counts",
        "filter:min_mito", "filter:max_mito",
        "filter:min_cd3",
    ),
    cd8_label: str = "Confident CD8+",
    donor: str | None = None,
) -> Path | None:
    """Per-sample cells bar with three layered cohorts (#76).

    Stacks three bars at the same baseline so the smaller bars sit
    visually *inside* the larger one (a "filled subset" view):

    1. **Loaded barcodes** (pale gray) — all rows in ``adata.obs``.
    2. **Passing scRNA QC** (medium blue) — when ``qc_cols`` are all
       present, the AND of those mask columns; otherwise this layer
       equals the loaded count.
    3. **αβ-pair denominator** (deep blue) — confident CD8+ ∧ complete
       αβ ∧ no doublet ∧ TRA_1_umis ≥ 2 ∧ TRB_1_umis ≥ 2. This is the
       denominator we report per-condition cell fractions against
       (#72). When the TCR-purity columns aren't present, this layer
       is omitted.

    Parameters
    ----------
    adata : AnnData
        Must have ``obs["sample"]``. Optionally ``filter:*`` mask
        columns (boolean) and the TCR-purity flags
        ``is_confident``, ``has_both_chains``, ``multi_chain``,
        ``Tcell_type`` / ``Tcell_type_consensus``, ``TRA_1_umis``,
        ``TRB_1_umis``.
    output_path : str | Path
        File to write the figure to.
    qc_cols : tuple[str, ...]
        Per-cell boolean QC columns to AND together for the QC layer.
    cd8_label : str
        Value of ``Tcell_type`` / ``Tcell_type_consensus`` indicating a
        confident-CD8+ cell.
    donor : str | None
        Optional donor name added to the title.

    Returns
    -------
    Path | None
        Output path on success.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    obs = adata.obs
    if "sample" not in obs.columns:
        logger.info("plot_cells_per_sample_stacked: no 'sample' column; skipping.")
        return None

    samples = sorted(obs["sample"].astype(str).unique())
    samp = obs["sample"].astype(str).values
    loaded_n = (
        pd.Series(np.ones(len(obs), dtype=int), index=samp)
        .groupby(level=0).sum().reindex(samples, fill_value=0)
    )

    has_qc = all(c in obs.columns for c in qc_cols)
    if has_qc:
        qc_mask = np.logical_and.reduce(
            [obs[c].fillna(False).astype(bool).values for c in qc_cols]
        )
    else:
        qc_mask = np.ones(len(obs), dtype=bool)
    qc_n = (
        pd.Series(qc_mask.astype(int), index=samp)
        .groupby(level=0).sum().reindex(samples, fill_value=0)
    )

    tcell_col = None
    for cand in ("Tcell_type_consensus", "Tcell_type"):
        if cand in obs.columns:
            tcell_col = cand
            break
    has_tcr = (
        "has_both_chains" in obs.columns
        and "multi_chain" in obs.columns
        and "is_confident" in obs.columns
        and "TRA_1_umis" in obs.columns
        and "TRB_1_umis" in obs.columns
        and tcell_col is not None
    )
    if has_tcr:
        is_cd8 = (obs[tcell_col].astype(str) == cd8_label).values
        tcr_mask = (
            obs["is_confident"].fillna(False).astype(bool).values
            & is_cd8
            & obs["has_both_chains"].fillna(False).astype(bool).values
            & (~obs["multi_chain"].fillna(False).astype(bool).values)
            & (obs["TRA_1_umis"].fillna(0).astype(float).values >= 2)
            & (obs["TRB_1_umis"].fillna(0).astype(float).values >= 2)
        )
        denom_n = (
            pd.Series((qc_mask & tcr_mask).astype(int), index=samp)
            .groupby(level=0).sum().reindex(samples, fill_value=0)
        )
    else:
        denom_n = None

    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(samples) + 2), 4.8))
    x = np.arange(len(samples))
    ax.bar(
        x, loaded_n.values, color="#cbd5e1", edgecolor="white",
        label="Loaded cells (all barcodes)",
    )
    if has_qc:
        ax.bar(
            x, qc_n.values, color="#60a5fa", edgecolor="white",
            label="Passing scRNA QC",
        )
    if denom_n is not None:
        ax.bar(
            x, denom_n.values, color="#1d4ed8", edgecolor="white",
            label="αβ-pair denominator (confident CD8+, complete αβ, "
                  "no doublet, min_umi ≥ 2)",
        )
    for xi, lv in zip(x, loaded_n.values):
        ax.text(
            xi, lv, f"{int(lv):,}", ha="center", va="bottom",
            fontsize=8.5, color="#475569",
        )
    if denom_n is not None:
        for xi, dv in zip(x, denom_n.values):
            if dv > 0:
                ax.text(
                    xi, dv, f"{int(dv):,}", ha="center", va="bottom",
                    fontsize=8, color="white", fontweight="bold",
                )
    ax.set_xticks(x)
    ax.set_xticklabels(pretty_samples(samples), rotation=30, ha="right")
    ax.set_ylabel("Cells")
    title = "Cells per sample"
    if donor:
        title += f" — donor {donor}"
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    save_figure(fig, output_path)
    return output_path


def plot_qc(adata: ad.AnnData, output_dir: str | Path):
    """Generate QC plots for loaded data."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Reads per cell distribution
    if "n_counts" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        for sample in adata.obs["sample"].unique():
            sample_data = adata.obs[adata.obs["sample"] == sample]["n_counts"]
            ax.hist(sample_data, bins=50, alpha=0.5, label=pretty_sample(sample))
        ax.set_xlabel("Total Counts per Cell")
        ax.set_ylabel("Number of Cells")
        ax.set_title("Read Count Distribution by Sample")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xscale("log")
        save_figure(fig, output_dir / "qc_read_counts.png")

    # Genes detected distribution
    if "n_genes" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        samples = adata.obs["sample"].unique()
        data = [adata.obs[adata.obs["sample"] == s]["n_genes"].values for s in samples]
        ax.violinplot(data, positions=range(len(samples)))
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(pretty_samples(samples), rotation=45, ha="right")
        ax.set_ylabel("Genes Detected")
        ax.set_title("Gene Detection by Sample")
        save_figure(fig, output_dir / "qc_genes_detected.png")

    # Mitochondrial percentage
    if "percent_mt" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        samples = adata.obs["sample"].unique()
        data = [adata.obs[adata.obs["sample"] == s]["percent_mt"].values for s in samples]
        ax.violinplot(data, positions=range(len(samples)))
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(pretty_samples(samples), rotation=45, ha="right")
        ax.set_ylabel("Mitochondrial %")
        ax.set_title("Mitochondrial Content by Sample")
        ax.axhline(y=8, color="red", linestyle="--", alpha=0.5, label="Max threshold")
        ax.axhline(y=2, color="orange", linestyle="--", alpha=0.5, label="Min threshold")
        ax.legend()
        save_figure(fig, output_dir / "qc_mito_percent.png")

    # Cells per sample — stacked layout when QC/TCR masks are present
    # (#76); falls back to a single-bar count when only ``sample`` is
    # available on ``adata.obs``.
    plot_cells_per_sample_stacked(adata, output_dir / "qc_cells_per_sample.png")


# =============================================================================
# Phenotype Plots
# =============================================================================


def plot_phenotype(adata: ad.AnnData, output_dir: str | Path):
    """Generate phenotype plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # CD4 vs CD8 scatter
    if "CD4" in adata.obs.columns and "CD8" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 10))

        if "Tcell_type" in adata.obs.columns:
            for tcell_type in adata.obs["Tcell_type"].unique():
                mask = adata.obs["Tcell_type"] == tcell_type
                ax.scatter(
                    adata.obs.loc[mask, "CD4"] + 0.1,
                    adata.obs.loc[mask, "CD8"] + 0.1,
                    alpha=0.3,
                    s=10,
                    label=tcell_type,
                )
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax.scatter(adata.obs["CD4"] + 0.1, adata.obs["CD8"] + 0.1, alpha=0.3, s=10)

        ax.set_xlabel("CD4 Expression")
        ax.set_ylabel("CD8 Expression (CD8A + CD8B)")
        ax.set_title("CD4 vs CD8 Expression")
        ax.set_xscale("log")
        ax.set_yscale("log")
        save_figure(fig, output_dir / "phenotype_cd4_cd8_scatter.png")

    # T cell type composition by sample
    if "Tcell_type" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(14, 8))

        # Create stacked bar chart
        samples = adata.obs["sample"].unique()
        tcell_types = adata.obs["Tcell_type"].cat.categories

        # Calculate percentages
        data = []
        for sample in samples:
            sample_data = adata.obs[adata.obs["sample"] == sample]
            total = len(sample_data)
            row = {"sample": sample}
            for tt in tcell_types:
                row[tt] = (sample_data["Tcell_type"] == tt).sum() / total * 100
            data.append(row)

        df_plot = pd.DataFrame(data)
        df_plot = df_plot.set_index("sample")
        # pandas bar plot uses the index as tick labels — pretty-rename
        # so the x-axis renders ``AIM⁺`` etc. instead of ``AIMpos-2`` (#109).
        df_plot.index = pretty_samples(df_plot.index)

        df_plot.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
        ax.set_ylabel("Percentage")
        ax.set_title("T Cell Type Composition by Sample")
        ax.legend(title="T Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.xticks(rotation=45, ha="right")
        save_figure(fig, output_dir / "phenotype_composition.png")


# =============================================================================
# Clonotype Plots
# =============================================================================


def plot_clonotypes(clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate clonotype analysis plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Clone size distribution (rank plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    sorted_counts = clonotypes["cell_count"].sort_values(ascending=False).values
    ax.plot(range(1, len(sorted_counts) + 1), sorted_counts, linewidth=2)
    ax.set_xlabel("Clone Rank")
    ax.set_ylabel("Cell Count")
    ax.set_title("Clone Size Distribution (Rank Plot)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.fill_between(range(1, len(sorted_counts) + 1), sorted_counts, alpha=0.3)
    save_figure(fig, output_dir / "clonotype_rank_plot.png")

    # Clone size histogram. cell_count is float under attribution (#176), so
    # ceil to an int range for integer histogram bins.
    fig, ax = plt.subplots(figsize=(10, 6))
    max_size = int(min(50, np.ceil(clonotypes["cell_count"].max())))
    ax.hist(
        clonotypes["cell_count"].clip(upper=max_size),
        bins=range(1, max_size + 2),
        edgecolor="black",
    )
    ax.set_xlabel("Clone Size (cells)")
    ax.set_ylabel("Number of Clones")
    ax.set_title("Clone Size Distribution")
    save_figure(fig, output_dir / "clonotype_size_histogram.png")

    # V gene usage heatmap (if available)
    for chain, gene_col in [("alpha", "alpha_v_gene"), ("beta", "beta_v_gene")]:
        if gene_col in clonotypes.columns:
            v_genes = clonotypes[gene_col].dropna()
            if len(v_genes) > 0:
                fig, ax = plt.subplots(figsize=(14, 8))
                gene_counts = v_genes.value_counts().head(20)
                ax.barh(range(len(gene_counts)), gene_counts.values)
                ax.set_yticks(range(len(gene_counts)))
                ax.set_yticklabels(gene_counts.index)
                ax.set_xlabel("Number of Clonotypes")
                ax.set_title(f"{chain.upper()} V Gene Usage (Top 20)")
                ax.invert_yaxis()
                save_figure(fig, output_dir / f"clonotype_{chain}_v_gene_usage.png")

    # Sample sharing matrix
    if "samples" in clonotypes.columns:
        samples = set()
        for s in clonotypes["samples"].dropna():
            samples.update(s.split(";"))
        samples = sorted(samples)

        if len(samples) > 1:
            # Calculate Jaccard similarity
            sample_clones = {}
            for sample in samples:
                mask = clonotypes["samples"].fillna("").str.contains(sample)
                sample_clones[sample] = set(clonotypes.loc[mask, "CDR3ab"])

            jaccard_matrix = np.zeros((len(samples), len(samples)))
            for i, s1 in enumerate(samples):
                for j, s2 in enumerate(samples):
                    intersection = len(sample_clones[s1] & sample_clones[s2])
                    union = len(sample_clones[s1] | sample_clones[s2])
                    jaccard_matrix[i, j] = intersection / union if union > 0 else 0

            fig, ax = plt.subplots(figsize=(10, 8))
            pretty = pretty_samples(samples)
            sns.heatmap(
                jaccard_matrix,
                xticklabels=pretty,
                yticklabels=pretty,
                annot=True,
                fmt=".2f",
                cmap="viridis",
                ax=ax,
            )
            ax.set_title("Sample Sharing (Jaccard Similarity)")
            plt.xticks(rotation=45, ha="right")
            save_figure(fig, output_dir / "clonotype_sharing_heatmap.png")

    # Rank-family bundle (#43 — #2/#3/#4/#5). Each helper skips itself
    # when prerequisites are missing, so callers without
    # ``sample_frequencies`` get only the legacy clonotype plots above.
    plot_rank_family(clonotypes, output_dir)


# =============================================================================
# Filter Plots
# =============================================================================


def plot_filter(clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate filtering analysis plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Frequency vs condition scatter (the key plot)
    if "max_frequency" in clonotypes.columns and "n_samples" in clonotypes.columns:
        fig, ax = plt.subplots(figsize=(12, 8))

        if "tier" in clonotypes.columns:
            tiers = sorted(clonotypes["tier"].dropna().unique())
            colors = plt.cm.viridis(np.linspace(0, 1, len(tiers)))

            for tier, color in zip(tiers, colors):
                mask = clonotypes["tier"] == tier
                ax.scatter(
                    clonotypes.loc[mask, "max_frequency"] * 100,
                    clonotypes.loc[mask, "n_samples"],
                    c=[color],
                    alpha=0.6,
                    s=clonotypes.loc[mask, "cell_count"] * 5,
                    label=tier,
                )
            ax.legend(title="Tier", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax.scatter(
                clonotypes["max_frequency"] * 100,
                clonotypes["n_samples"],
                alpha=0.5,
                s=clonotypes["cell_count"] * 5,
            )

        ax.set_xlabel("Max Frequency (%)")
        ax.set_ylabel("Number of Conditions")
        ax.set_title("Clone Frequency vs Condition Specificity\n(size = cell count)")
        ax.set_xscale("log")
        save_figure(fig, output_dir / "filter_frequency_specificity.png")

    # Tier distribution
    if "tier" in clonotypes.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        tier_counts = clonotypes["tier"].value_counts().sort_index()
        ax.bar(
            range(len(tier_counts)),
            tier_counts.values,
            color=plt.cm.viridis(np.linspace(0, 1, len(tier_counts))),
        )
        ax.set_xticks(range(len(tier_counts)))
        ax.set_xticklabels(tier_counts.index)
        ax.set_ylabel("Number of Clonotypes")
        ax.set_title("Clonotype Distribution by Tier")
        for i, v in enumerate(tier_counts.values):
            ax.text(i, v + 1, str(v), ha="center", fontsize=12, fontweight="bold")
        save_figure(fig, output_dir / "filter_tier_distribution.png")

    # CD4/CD8 by tier
    if "tier" in clonotypes.columns and "Tcell_type_consensus" in clonotypes.columns:
        fig, ax = plt.subplots(figsize=(12, 6))

        tiers = sorted(clonotypes["tier"].dropna().unique())
        cd8_counts = []
        cd4_counts = []

        for tier in tiers:
            tier_data = clonotypes[clonotypes["tier"] == tier]
            cd8_counts.append(tier_data["Tcell_type_consensus"].str.contains("CD8", na=False).sum())
            cd4_counts.append(tier_data["Tcell_type_consensus"].str.contains("CD4", na=False).sum())

        x = np.arange(len(tiers))
        width = 0.35

        ax.bar(x - width / 2, cd8_counts, width, label="CD8+")
        ax.bar(x + width / 2, cd4_counts, width, label="CD4+")
        ax.set_xticks(x)
        ax.set_xticklabels(tiers)
        ax.set_ylabel("Number of Clonotypes")
        ax.set_title("T Cell Type Distribution by Tier")
        ax.legend()
        save_figure(fig, output_dir / "filter_tcell_type_by_tier.png")


# =============================================================================
# Annotation Plots
# =============================================================================


def plot_annotations(clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate annotation plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if "db_match" not in clonotypes.columns:
        logger.warning("No annotation data found in clonotypes")
        return

    # Viral vs non-viral pie chart
    if "is_viral" in clonotypes.columns:
        fig, ax = plt.subplots(figsize=(8, 8))

        matched = clonotypes["db_match"].sum()
        viral = clonotypes["is_viral"].sum()
        non_viral_matched = matched - viral
        unmatched = len(clonotypes) - matched

        sizes = [viral, non_viral_matched, unmatched]
        labels = [
            f"Viral ({viral})",
            f"Non-viral matched ({non_viral_matched})",
            f"Unmatched ({unmatched})",
        ]
        colors = ["red", "orange", "gray"]

        ax.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90)
        ax.set_title("Database Match Summary")
        save_figure(fig, output_dir / "annotation_summary_pie.png")

    # Top species/epitopes
    if "db_species" in clonotypes.columns:
        species_counts = clonotypes["db_species"].dropna().value_counts().head(10)
        if len(species_counts) > 0:
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.barh(range(len(species_counts)), species_counts.values)
            ax.set_yticks(range(len(species_counts)))
            ax.set_yticklabels(species_counts.index)
            ax.set_xlabel("Number of Clonotypes")
            ax.set_title("Top 10 Matched Species/Antigens")
            ax.invert_yaxis()
            save_figure(fig, output_dir / "annotation_top_species.png")

    # Phase-B granularity-aware plot family. Each helper skips itself
    # when prerequisites are missing, so it's safe to call on
    # pre-#48-phase-A annotations too.
    plot_annotations_phase_b(clonotypes, output_dir)


# =============================================================================
# Structured Annotation Plots — #48 phase B
# =============================================================================
#
# Four plots that consume the structured-annotation columns produced by
# match_clonotypes (db_category / db_protein_canonical / db_species /
# db_epitope / db_mhc / db_match_strength), built on a small set of
# shared helpers. The granularity selector is wired as a single
# ``granularity`` parameter that maps a human-friendly label to the
# underlying column.


# Categorical color palette used across all phase-B plots. Mapping is
# stable across runs so the same category always gets the same color
# in publication figures. Anything not in this dict falls back to a
# muted grey ("other"/"unknown").
DB_CATEGORY_PALETTE: dict[str, str] = {
    "tumor_self": "#d62728",   # red — the cohort's target axis
    "viral": "#1f77b4",        # blue — bystander viral
    "bacterial": "#2ca02c",    # green
    "self": "#ff7f0e",         # orange — non-tumor self
    "other": "#7f7f7f",        # grey
    "unknown": "#bcbcbc",      # light grey
    None: "#e5e5e5",           # very light grey (unmatched)
}


# Human-friendly granularity labels → underlying column name. Order is
# coarsest → finest so the CLI loop produces a sensible filename order.
GRANULARITY_COLUMNS: dict[str, str] = {
    "category": "db_category",
    "organism": "db_species",
    "protein": "db_protein_canonical",
    "peptide": "db_epitope",
}


def _expand_sample_frequencies(clonotypes: pd.DataFrame) -> pd.DataFrame:
    """Pivot the per-clone ``sample_frequencies`` dict into a wide table.

    Returns an N_clones × N_samples DataFrame indexed by ``CDR3ab`` with
    one column per sample (zero-filled where a clone is absent).
    """
    if "sample_frequencies" not in clonotypes.columns:
        return pd.DataFrame(index=clonotypes.get("CDR3ab", pd.Index([], name="CDR3ab")))
    freqs = clonotypes["sample_frequencies"].dropna()
    if len(freqs) == 0:
        return pd.DataFrame(index=clonotypes.get("CDR3ab", pd.Index([], name="CDR3ab")))
    expanded = pd.DataFrame(list(freqs.values), index=clonotypes.loc[freqs.index, "CDR3ab"])
    return expanded.fillna(0.0).sort_index(axis=1)


def _top_n_clones(
    clonotypes: pd.DataFrame, n: int = 20, by: str = "max_frequency_per_method"
) -> pd.DataFrame:
    """Return the top-N clones by the chosen ranking column.

    Falls back through ranking columns when the preferred one is
    missing (different upstream paths produce different column sets).
    """
    df = clonotypes
    for candidate in (by, "max_frequency_per_method", "max_frequency", "cell_count"):
        if candidate in df.columns:
            return df.nlargest(n, candidate)
    # No usable ranking column — return the head as a last resort so
    # callers still get a non-empty plot rather than an exception.
    return df.head(n)


def _category_color_strip(
    series: pd.Series, palette: dict[str | None, str] | None = None
) -> list[str]:
    """Map each value in ``series`` to a color, honouring the shared
    category palette where possible, otherwise hashing the value to a
    matplotlib tab20 slot."""
    palette = palette or DB_CATEGORY_PALETTE
    fallback_cmap = plt.get_cmap("tab20")
    unique = [v for v in series.dropna().unique() if v not in palette]
    unique_to_color = {
        val: fallback_cmap(i % 20) for i, val in enumerate(sorted(map(str, unique)))
    }
    return [
        palette[v] if v in palette else unique_to_color.get(str(v), "#cccccc")
        for v in series
    ]


def _resolve_granularity(granularity: str) -> str:
    """Map a human-friendly granularity label to its underlying column."""
    if granularity in GRANULARITY_COLUMNS:
        return GRANULARITY_COLUMNS[granularity]
    # Allow passing the raw column name directly for advanced callers.
    return granularity


def plot_matched_clone_heatmap(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    top_n: int = 20,
    granularity: str = "category",
    filename: str | None = None,
    title_suffix: str = "",
) -> Path | None:
    """Top-N matched clones × conditions heatmap with annotation strip.

    Rows are the top-N matched clones by within-sample frequency.
    Columns are the samples present in ``sample_frequencies``. A
    left-side colored strip carries category / protein / peptide / MHC
    annotations so the eye can read specificity directly off the heatmap.

    Skips silently (and returns None) when the input lacks either the
    structured-annotation columns or any matched clones.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # ``granularity`` shows up in the title + filename. The heatmap's
    # annotation strip always carries every available label (category,
    # protein, peptide, MHC); the parameter is informational so the
    # CLI loop produces one heatmap PNG per granularity.
    _ = _resolve_granularity(granularity)

    if "db_match" not in clonotypes.columns:
        logger.info("plot_matched_clone_heatmap: no db_match column; skipping.")
        return None
    matched = clonotypes[clonotypes["db_match"].fillna(False).astype(bool)]
    if len(matched) == 0:
        logger.info("plot_matched_clone_heatmap: no matched clones; skipping.")
        return None

    top = _top_n_clones(matched, n=top_n)
    if len(top) == 0:
        return None

    freq_wide = _expand_sample_frequencies(top)
    if freq_wide.empty or freq_wide.shape[1] == 0:
        logger.info("plot_matched_clone_heatmap: no sample_frequencies; skipping.")
        return None

    # Annotation strip — only include columns that exist + carry data.
    strip_candidates = ("db_category", "db_protein_canonical", "db_epitope", "db_mhc")
    strip_cols = [c for c in strip_candidates if c in top.columns and top[c].notna().any()]
    row_colors: pd.DataFrame | None = None
    if strip_cols:
        # Build a DataFrame the same height as freq_wide, indexed the
        # same way, with a colour per (clone, annotation field).
        strip_indexed = top.set_index("CDR3ab")[strip_cols].reindex(freq_wide.index)
        row_colors = pd.DataFrame(
            {col: _category_color_strip(strip_indexed[col]) for col in strip_cols},
            index=freq_wide.index,
        )

    grid = sns.clustermap(
        freq_wide,
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        cmap="rocket_r",
        cbar_kws={"label": "within-sample frequency"},
        linewidths=0,
        figsize=(max(10, 1 + 0.6 * freq_wide.shape[1]), max(8, 0.35 * top_n)),
        xticklabels=True,
        yticklabels=True,
    )
    grid.ax_heatmap.set_xlabel("sample / condition")
    grid.ax_heatmap.set_ylabel("clone (CDR3αβ)")
    title = f"Top-{top_n} matched clones × conditions"
    if title_suffix:
        title += f" — {title_suffix}"
    grid.figure.suptitle(title, y=1.02, fontsize=14)

    filename = filename or f"annotation_heatmap_{granularity}.png"
    out_path = output_dir / filename
    save_figure(grid.figure, out_path, dpi=300)
    logger.info(f"Saved plot to {out_path}")
    return out_path


def plot_category_composition_bar(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    granularity: str = "category",
    top_n: int | None = None,
    filename: str | None = None,
) -> Path | None:
    """Stacked composition bar per condition.

    For each condition (sample), shows the fraction of clones falling
    into each value of the chosen granularity (category / protein /
    organism / peptide). Unmatched clones become a dedicated ``unmatched``
    bucket so the y-axis sums to 1.0.

    When ``top_n`` is set, only the top-N clones (by max frequency) per
    condition are counted — matches the heatmap's "top-N" framing.
    Otherwise all matched clones are counted.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gran_col = _resolve_granularity(granularity)

    if gran_col not in clonotypes.columns:
        logger.info(
            f"plot_category_composition_bar: column {gran_col} not present; skipping."
        )
        return None
    freq_wide = _expand_sample_frequencies(clonotypes)
    if freq_wide.empty or freq_wide.shape[1] == 0:
        logger.info(
            "plot_category_composition_bar: no sample_frequencies; skipping."
        )
        return None

    # For each sample, pick the clones present (frequency > 0); take
    # top-N within that sample if requested.
    samples = list(freq_wide.columns)
    annotation_by_cdr3ab = clonotypes.set_index("CDR3ab")[gran_col]
    annotation_by_cdr3ab = annotation_by_cdr3ab.fillna("unmatched").astype(str)

    rows: list[dict] = []
    for sample in samples:
        present = freq_wide[freq_wide[sample] > 0][sample]
        if top_n is not None:
            present = present.nlargest(top_n)
        if len(present) == 0:
            continue
        labels = annotation_by_cdr3ab.reindex(present.index).fillna("unmatched")
        for label, count in labels.value_counts().items():
            rows.append({"sample": sample, "group": label, "count": int(count)})
    if not rows:
        return None

    long = pd.DataFrame(rows)
    pivot = long.pivot_table(
        index="sample", columns="group", values="count", aggfunc="sum", fill_value=0
    )
    fractions = pivot.div(pivot.sum(axis=1), axis=0)

    # Stable column order — categories first in palette order, then
    # everything else alphabetically.
    palette_order = [k for k in DB_CATEGORY_PALETTE if isinstance(k, str)]
    ordered = [c for c in palette_order if c in fractions.columns]
    ordered += sorted(c for c in fractions.columns if c not in ordered)
    fractions = fractions[ordered]

    # Anything missing from the curated palette falls back to a muted
    # grey rather than None (pandas's plot backend rejects None).
    colors = [DB_CATEGORY_PALETTE.get(c) or "#cccccc" for c in fractions.columns]
    fig, ax = plt.subplots(figsize=(max(8, 1.0 * len(fractions)), 6))
    fractions.plot(kind="bar", stacked=True, ax=ax, color=colors, edgecolor="white")
    ax.set_ylabel("fraction of clones")
    ax.set_xlabel("condition / sample")
    suffix = f" (top {top_n} per sample)" if top_n else ""
    ax.set_title(f"Composition by {granularity}{suffix}")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), title=granularity)
    ax.set_ylim(0, 1.0)

    filename = filename or f"annotation_composition_{granularity}.png"
    out_path = output_dir / filename
    save_figure(fig, out_path)
    return out_path


def plot_match_strength_comparison(
    clonotypes: pd.DataFrame,
    output_dir: str | Path,
    top_n: int = 20,
    granularity: str = "category",
) -> Path | None:
    """αβ matches vs β-only fallback, side-by-side.

    Renders the heatmap twice: once filtered to ``db_match_strength == "ab"``
    (full αβ pairing), once to ``"b_only"`` (β-only fallback). The
    visual delta makes it obvious how much of the headline match
    signal disappears when full pairing is required. Skips when
    ``db_match_strength`` is absent.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if "db_match_strength" not in clonotypes.columns:
        logger.info(
            "plot_match_strength_comparison: db_match_strength not present; skipping."
        )
        return None

    written: list[Path] = []
    for strength in ("ab", "b_only"):
        subset = clonotypes[clonotypes["db_match_strength"] == strength]
        if len(subset) == 0:
            continue
        out = plot_matched_clone_heatmap(
            subset,
            output_dir=output_dir,
            top_n=top_n,
            granularity=granularity,
            filename=f"annotation_heatmap_{granularity}_{strength}.png",
            title_suffix=f"match strength = {strength}",
        )
        if out is not None:
            written.append(out)
    return written[0] if written else None


def plot_annotations_phase_b(
    clonotypes: pd.DataFrame, output_dir: str | Path, top_n: int = 20
) -> None:
    """One-shot orchestrator: loop the granularity selector across every
    granularity and emit the heatmap, composition bar, and αβ/β-only
    comparison for each. Wired into the CLI ``plot-annotations`` flag
    so users get the full phase-B family with one switch.

    Skips silently when prerequisites are missing — works fine on
    pre-#48-phase-A data, just produces fewer files.
    """
    for gran in GRANULARITY_COLUMNS:
        plot_matched_clone_heatmap(clonotypes, output_dir, top_n=top_n, granularity=gran)
        plot_category_composition_bar(clonotypes, output_dir, granularity=gran, top_n=top_n)
        plot_match_strength_comparison(
            clonotypes, output_dir, top_n=top_n, granularity=gran
        )


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
        clonotypes.set_index("CDR3ab")["tier"]
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


def _grid_shape(n_panels: int, max_cols: int = 3) -> tuple[int, int]:
    """Pick (rows, cols) for a multi-panel layout."""
    import math

    cols = min(max(1, n_panels), max_cols)
    rows = math.ceil(n_panels / cols)
    return rows, cols


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
    save_figure(fig, out_path)
    return out_path


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
    save_figure(fig, out_path)
    return out_path


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
        clonotypes.set_index("CDR3ab")["cell_count"]
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
    save_figure(fig, out_path)
    return out_path


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
    save_figure(fig, out_path)
    return out_path


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
        for i, cdr3 in enumerate(mat.index):
            ax.plot(
                x_pos,
                mat.loc[cdr3].values,
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
    save_figure(fig, out_path)
    return out_path


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
            lambda d: sum(1 for v in (d or {}).values() if v > 0)
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
    save_figure(fig, out_path)
    return out_path


# =============================================================================
# Signature scatter — #43 item #7
# =============================================================================
#
# Gene-set defaults live in :mod:`tcrsift.signatures` (top-of-file
# import above). The per-sample scatter takes any ``gene_ids`` tuple,
# so callers can drop in custom signatures via that argument.


def _per_cell_signature(adata: ad.AnnData, gene_ids: list[str]) -> np.ndarray:
    """Mean log1p expression across ``gene_ids`` per cell.

    Drops any IDs not present in ``adata.var_names``; returns a
    zero-vector if none of the requested genes are in the matrix.
    """
    valid = [g for g in gene_ids if g in adata.var_names]
    if not valid:
        return np.zeros(adata.n_obs)
    idx = [adata.var_names.get_loc(g) for g in valid]
    X = adata.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.log1p(X).mean(axis=1)


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
        freqs = row.get("sample_frequencies") or {}
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
    save_figure(fig, out_path)
    return out_path


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

    # TIL recovery by tier
    if "tier" in matched_clonotypes.columns:
        fig, ax = plt.subplots(figsize=(10, 6))

        tiers = sorted(matched_clonotypes["tier"].dropna().unique())
        recovery_rates = []

        for tier in tiers:
            tier_data = matched_clonotypes[matched_clonotypes["tier"] == tier]
            recovery = (
                tier_data["til_match"].sum() / len(tier_data) * 100 if len(tier_data) > 0 else 0
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
    matched = matched_clonotypes[matched_clonotypes["til_match"]]
    if (
        len(matched) > 0
        and "max_frequency" in matched.columns
        and "til_frequency" in matched.columns
    ):
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(
            matched["max_frequency"] * 100,
            matched["til_frequency"] * 100,
            alpha=0.6,
            s=matched["cell_count"] * 10,
        )
        ax.set_xlabel("Culture Frequency (%)")
        ax.set_ylabel("TIL Frequency (%)")
        ax.set_title("Culture vs TIL Frequency\n(size = culture cell count)")

        # Add diagonal line
        max_val = max(matched["max_frequency"].max(), matched["til_frequency"].max()) * 100
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="1:1 line")
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
    from .format import order_conditions

    ordered_names = order_conditions(sets)
    pretty = {name: pretty_method(name) for name in ordered_names}
    # Guard against two raw names collapsing to one pretty label.
    if len(set(pretty.values())) == len(pretty):
        sets = {pretty[name]: sets[name] for name in ordered_names}

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
            fig = plt.figure(figsize=(max(7, 0.6 * len(sets) + 5), 5))
            upsetplot.plot(
                data, fig=fig, sort_by="cardinality", show_counts=True,
                min_subset_size=min_subset_size,
            )
            fig.suptitle(title)
            save_figure(fig, output_path, dpi=150)
        return
    except ImportError:
        pass

    # Fallback: bar chart of the top intersection patterns.
    from .clonotype import set_overlap_table

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
    ax.set_title(f"{title} (install 'upsetplot' for UpSet view)")
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
    vmax = float(matrix.values.max()) if matrix.size else 1.0
    fmt = "d" if similarity == "count" else ".2f"
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
    sns.heatmap(
        pivot,
        annot=True,
        fmt="d",
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
                ax.set_title(f"{label}: {a} vs {b}\n(install 'matplotlib_venn' for Venn view)",
                             fontsize=9)
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
    present = [(c, p) for c, p in vgene_cols if c in clonotypes.columns]
    if not present:
        return None
    matrices = []
    for vgcol, pretty in present:
        merged = long_df.merge(clonotypes[[clone_col, vgcol]], on=clone_col, how="inner")
        merged = merged[merged[vgcol].notna()]
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


# =============================================================================
# Funnel Plot
# =============================================================================


FUNNEL_LABEL_NICE: dict[str, str] = {
    # Raw mask-column names → readable funnel labels (#72). The
    # ``filter:min_counts`` / ``filter:max_counts`` names are
    # particularly misleading — they read like row-count filters but
    # actually test per-cell *UMI counts*.
    "filter:min_counts": "min UMIs per cell",
    "filter:max_counts": "max UMIs per cell",
    "filter:min_genes": "min genes per cell",
    "filter:max_genes": "max genes per cell",
    "filter:min_mito": "min mito %",
    "filter:max_mito": "max mito %",
    "filter:min_cd3": "min CD3 reads",
    "filter:min_umi": "min TCR UMIs per chain",
}


def normalize_funnel_label(name: str) -> str:
    """Map a raw mask/stage key to a reader-friendly funnel label
    (#72). Pass-through for names not in :data:`FUNNEL_LABEL_NICE`."""
    return FUNNEL_LABEL_NICE.get(name, name)


def plot_funnel(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
    *,
    denominator_stage: str | None = None,
    section_starts: tuple[str, ...] | None = None,
    filename: str = "funnel_plot.png",
):
    """
    Generate a funnel plot showing TCR counts at each filtering stage.

    Parameters
    ----------
    stage_counts : dict
        Ordered dictionary mapping stage names to counts.
        Example: ``{"Raw Cells": 10000, "With VDJ": 8000, ...}``.
        Caller controls ordering — use this to put TCR-purity gates
        ahead of scRNA-QC gates (or vice versa) per #72.
    output_dir : str or Path
        Directory to save the plot.
    title : str
        Plot title.
    denominator_stage : str | None
        When set, annotate this stage as the per-condition cell-fraction
        denominator (the αβ-pair denominator, per #72). Draws an arrow
        + label to the right of that bar.
    section_starts : tuple[str, ...] | None
        Stage names that *start* a new logical section (e.g. the first
        scRNA-QC gate after the TCR-purity block). A thin horizontal
        divider line is drawn above each named stage so a reader can
        see where one gate group ends and the next begins.
    filename : str
        Output filename (default ``funnel_plot.png``).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = list(stage_counts.keys())
    counts = list(stage_counts.values())
    # `or 1` guards both the empty-dict and all-zero-counts cases (a sample
    # where no cells clear any gate) — otherwise width = count / max_count
    # divides by zero.
    max_count = (max(counts) if counts else 0) or 1

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create funnel bars
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(stages)))

    for i, (stage, count, color) in enumerate(zip(stages, counts, colors)):
        # Calculate bar width proportional to count
        width = count / max_count * 0.8
        left = (1 - width) / 2

        # Draw bar
        ax.barh(
            len(stages) - i - 1,
            width,
            left=left,
            height=0.7,
            color=color,
            edgecolor="white",
            linewidth=2,
        )

        # Add count label
        ax.text(
            0.5,
            len(stages) - i - 1,
            f"{stage}\n{count:,}",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            color="white" if count / max_count > 0.3 else "black",
        )

        # Add percentage from previous stage
        if i > 0 and counts[i - 1] > 0:
            pct = count / counts[i - 1] * 100
            ax.text(
                0.92,
                len(stages) - i - 1,
                f"{pct:.0f}%",
                ha="left",
                va="center",
                fontsize=10,
                color="gray",
            )

        # αβ-pair denominator callout (#72).
        if denominator_stage and stage == denominator_stage:
            y = len(stages) - i - 1
            ax.annotate(
                "← αβ-pair denominator",
                xy=(0.9, y),
                xytext=(1.15, y),
                ha="left", va="center",
                fontsize=10, color="#1d4ed8", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#1d4ed8", lw=1.4),
            )

    # Section divider lines (#72). Drawn ABOVE the named stage row.
    if section_starts:
        for name in section_starts:
            if name not in stages:
                continue
            i = stages.index(name)
            y = len(stages) - i - 1 + 0.45
            ax.axhline(y=y, xmin=0.02, xmax=0.98, color="#94a3b8",
                       linestyle="--", linewidth=0.9, alpha=0.7)

    ax.set_xlim(0, 1.3 if denominator_stage else 1.1)
    ax.set_ylim(-0.5, len(stages) - 0.5)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    ax.axis("off")

    # Add overall retention
    if len(counts) >= 2 and counts[0] > 0:
        overall_pct = counts[-1] / counts[0] * 100
        ax.text(
            0.5,
            -0.3,
            f"Overall retention: {overall_pct:.1f}%",
            ha="center",
            va="top",
            fontsize=12,
            style="italic",
        )

    save_figure(fig, output_dir / filename)


# =============================================================================
# Funnel variants — #43 item #8
# =============================================================================
#
# Ribbon, lollipop, and terrace siblings to the existing bars funnel.
# All four siblings consume the same ``stage_counts: dict[str, int]``
# input as ``plot_funnel`` for API consistency. The shared
# `_log_stage_widths` helper makes the visual narrowing rate match
# between bars / ribbon / terrace.


def _log_stage_widths(counts: list[int], floor: float = 0.10) -> np.ndarray:
    """Map stage counts to relative widths in [floor, 1.0].

    Log-spaced so a 10× narrowing reads as a 10× width drop visually
    (rather than 10% as it would on a linear scale). The smallest
    stage gets ``floor`` width to stay visible; the largest gets 1.0.
    """
    arr = np.array([max(c, 1) for c in counts], dtype=float)
    log_counts = np.log10(arr)
    lo, hi = log_counts.min(), log_counts.max()
    if hi == lo:
        return np.ones(len(counts))
    return floor + (1.0 - floor) * (log_counts - lo) / (hi - lo)


def _draw_funnel_side_labels(
    ax: plt.Axes, stages: list[str], counts: list[int]
) -> None:
    """Draw stage labels on the left, count + retention% on the right."""
    initial = counts[0] if counts else 1
    for i, (stage, count) in enumerate(zip(stages, counts)):
        ax.text(
            -0.04, i, stage, ha="right", va="center",
            fontsize=10.5, fontweight="semibold", color="#1a1a1a",
        )
        frac = count / initial if initial else 0
        ax.text(
            1.04, i, f"n = {count:,}  ·  {frac * 100:.1f}% of initial",
            ha="left", va="center", fontsize=9.5, color="#444444",
        )


def plot_funnel_ribbon(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
    filename: str = "funnel_ribbon.png",
) -> Path:
    """Funnel as a smooth narrowing polygon (#43 #8).

    Interpolates stage widths along the vertical axis to render a
    continuous taper rather than discrete bars. Inner tick marks at
    each stage's exact width preserve the per-stage information.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = list(stage_counts.keys())
    counts = list(stage_counts.values())
    n = len(stages)
    widths = _log_stage_widths(counts)

    fig, ax = plt.subplots(figsize=(12, max(5, 0.7 * n + 2)))
    ys = np.arange(n)
    interp_y = np.linspace(0, n - 1, (n - 1) * 30 + 1) if n > 1 else np.array([0.0])
    interp_w = np.interp(interp_y, ys, widths)
    left = (1 - interp_w) / 2
    right = (1 + interp_w) / 2
    ax.fill_betweenx(interp_y, left, right, color="#fda4a4", alpha=0.85, linewidth=0)
    ax.plot(left, interp_y, color="#b91c1c", linewidth=2.0)
    ax.plot(right, interp_y, color="#b91c1c", linewidth=2.0)
    for i, w in enumerate(widths):
        ax.plot([(1 - w) / 2, (1 + w) / 2], [i, i],
                color="#b91c1c", linewidth=0.8, alpha=0.4)

    _draw_funnel_side_labels(ax, stages, counts)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(n - 0.5, -0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, pad=14)

    out_path = output_dir / filename
    save_figure(fig, out_path)
    return out_path


def plot_funnel_lollipop(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
    filename: str = "funnel_lollipop.png",
) -> Path:
    """Funnel as a horizontal lollipop on a log-x axis (#43 #8).

    One dot per stage at its actual count; lollipop "stick" runs from
    the dot out to a common right edge so the visual fall-off is
    obvious. Best for funnels with wide dynamic range where the
    other variants compress the small stages.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = list(stage_counts.keys())
    counts = np.array([max(c, 1) for c in stage_counts.values()], dtype=float)
    n = len(stages)

    fig, ax = plt.subplots(figsize=(10.5, max(4, 0.6 * n + 2)))
    y = np.arange(n)
    ax.plot(counts, y, color="#b91c1c", linewidth=2.0, zorder=2)
    ax.scatter(
        counts, y,
        s=np.linspace(95, 42, n),
        color="#b91c1c", edgecolor="white", linewidth=1.0, zorder=3,
    )
    for yi, stage, count in zip(y, stages, counts):
        ax.hlines(yi, xmin=count, xmax=counts.max() * 2.4,
                  color="#fda4a4", linewidth=2.3)
        ax.text(count / 1.10, yi, f"n = {int(count):,}",
                va="center", ha="left", fontsize=9.5)

    ax.set_xscale("log")
    ax.set_xlim(counts.max() * 2.4, counts.min() * 0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.set_xlabel("count (log scale)")
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_title(title, pad=14)

    out_path = output_dir / filename
    save_figure(fig, out_path)
    return out_path


def plot_funnel_terrace(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
    filename: str = "funnel_terrace.png",
) -> Path:
    """Funnel as stepped trapezoids between adjacent stages (#43 #8).

    Each transition is a filled trapezoid that visually emphasizes
    the *drop* between consecutive stages — useful for funnels where
    the per-step retention matters more than absolute counts.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = list(stage_counts.keys())
    counts = list(stage_counts.values())
    n = len(stages)
    widths = _log_stage_widths(counts)
    cmap = plt.get_cmap("Reds")

    fig, ax = plt.subplots(figsize=(12, max(5, 0.7 * n + 2)))
    for i in range(n - 1):
        w0, w1 = widths[i], widths[i + 1]
        y0, y1 = i + 0.22, i + 0.78
        poly = np.array(
            [
                [(1 - w0) / 2, y0], [(1 + w0) / 2, y0],
                [(1 + w1) / 2, y1], [(1 - w1) / 2, y1],
            ]
        )
        ax.fill(
            poly[:, 0], poly[:, 1],
            color=cmap(0.34 + 0.55 * i / max(1, n - 2)),
            alpha=0.9, edgecolor="white", linewidth=1.1,
        )
    for i, w in enumerate(widths):
        ax.plot([(1 - w) / 2, (1 + w) / 2], [i, i],
                color="#7f1d1d", linewidth=1.0, alpha=0.6)

    _draw_funnel_side_labels(ax, stages, counts)
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(n - 0.5, -0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, pad=14)

    out_path = output_dir / filename
    save_figure(fig, out_path)
    return out_path


def plot_funnel_variants(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
) -> list[Path]:
    """Emit all four funnel renderings (bars / ribbon / lollipop / terrace).

    Useful when prepping figures for an audience-uncertain context —
    each variant exposes a different aspect (absolute counts vs.
    per-step retention vs. dynamic range).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    # The legacy bars renderer writes to ``funnel_plot.png``; use the
    # same path so existing pipelines pick it up unchanged.
    plot_funnel(stage_counts, output_dir, title=title)
    paths.append(output_dir / "funnel_plot.png")
    paths.append(plot_funnel_ribbon(stage_counts, output_dir, title=title))
    paths.append(plot_funnel_lollipop(stage_counts, output_dir, title=title))
    paths.append(plot_funnel_terrace(stage_counts, output_dir, title=title))
    return paths


def create_pipeline_funnel(
    raw_cells: int,
    with_vdj: int,
    phenotyped: int,
    clonotypes: int,
    filtered: int,
    tier_counts: dict[str, int] | None = None,
    output_dir: str | Path = ".",
    *,
    ab_pair_denominator: int | None = None,
    selected_count: int | None = None,
    emit_selected_variant: bool = False,
    non_viral: int | None = None,
):
    """
    Create a funnel plot for the TCRsift pipeline stages.

    Parameters
    ----------
    raw_cells : int
        Number of cells after loading.
    with_vdj : int
        Number of cells with VDJ data.
    phenotyped : int
        Number of cells after phenotyping.
    clonotypes : int
        Number of unique clonotypes.
    filtered : int
        Number of clonotypes passing filters.
    tier_counts : dict | None
        Counts per confidence tier (tier1 / tier2 / …).
    output_dir : str | Path
        Output directory for the plot.
    ab_pair_denominator : int | None
        Per-cell count after the TCR-purity gates (confident CD8+ ∧
        complete αβ ∧ no doublet ∧ min UMI ≥ 2 per chain). When given,
        inserted between ``Phenotyped`` and ``Unique Clonotypes`` and
        annotated as the αβ-pair denominator (#72).
    selected_count : int | None
        Number of clones in the shortlist from
        :func:`tcrsift.candidate.select_candidates`. Only used when
        ``emit_selected_variant=True``.
    emit_selected_variant : bool
        When True, also emit a sibling funnel ``funnel_plot_selected.png``
        that collapses the tier cascade into a single ``Selected``
        stage (tier1 ∪ tier2 ∪ top-N-by-signature from tier3+).
    """
    stage_counts = {
        "Raw Cells": raw_cells,
        "With VDJ": with_vdj,
        "Phenotyped (CD4/CD8)": phenotyped,
    }
    if ab_pair_denominator is not None:
        stage_counts["αβ-pair denominator"] = ab_pair_denominator
    stage_counts["Unique Clonotypes"] = clonotypes
    stage_counts["Passing Filters"] = filtered

    tier_stage_counts = dict(stage_counts)
    if tier_counts:
        for tier, count in tier_counts.items():
            tier_stage_counts[f"Tier: {tier}"] = count

    section_starts = (
        ("Unique Clonotypes",)
        if ab_pair_denominator is not None
        else None
    )
    # Emit the full cascade in all four complementary styles (#255): bars
    # carries the αβ-pair denominator callout + section dividers; ribbon /
    # lollipop / terrace are alternate views of the SAME stage counts (each
    # exposes a different aspect — absolute counts vs. per-step retention vs.
    # dynamic range), so a reader / figure-assembler can pick what reads best.
    plot_funnel(
        tier_stage_counts,
        output_dir,
        denominator_stage=(
            "αβ-pair denominator" if ab_pair_denominator is not None else None
        ),
        section_starts=section_starts,
    )
    plot_funnel_ribbon(tier_stage_counts, output_dir)
    plot_funnel_lollipop(tier_stage_counts, output_dir)
    plot_funnel_terrace(tier_stage_counts, output_dir)

    if emit_selected_variant and selected_count is not None:
        selected_stage_counts = dict(stage_counts)
        # When viral bystanders were excluded from selection, show the
        # post-viral-filter survivor count as its own stage so the drop
        # is visible in the shortlist funnel (#122 exclude_viral).
        if non_viral is not None:
            selected_stage_counts["Non-viral"] = non_viral
        selected_stage_counts["Selected"] = selected_count
        # Selected-shortlist overlay variant in all four styles too (#255):
        # collapses the tier cascade into a single "Selected" stage so the
        # reader sees how many clones survived each gate AND where the final
        # picks land.
        sel_title = "TCR Selection Funnel (Selected shortlist)"
        plot_funnel(
            selected_stage_counts,
            output_dir,
            denominator_stage=(
                "αβ-pair denominator" if ab_pair_denominator is not None else None
            ),
            section_starts=section_starts,
            filename="funnel_plot_selected.png",
            title=sel_title,
        )
        plot_funnel_ribbon(
            selected_stage_counts, output_dir,
            title=sel_title, filename="funnel_ribbon_selected.png",
        )
        plot_funnel_lollipop(
            selected_stage_counts, output_dir,
            title=sel_title, filename="funnel_lollipop_selected.png",
        )
        plot_funnel_terrace(
            selected_stage_counts, output_dir,
            title=sel_title, filename="funnel_terrace_selected.png",
        )


# =============================================================================
# Color-Coded TCR Sequence PDF
# =============================================================================


# Canonical assembled-sequence columns rendered by the TCR sequence
# PDF, listed in biological assembly order. The linker direction in
# tcrsift is β-T2A-α, so the β block comes first. Earlier versions
# had (1) the α leader before the β block, (2) lowercase
# ``vdj_{chain}_aa`` and uppercase ``VDJ_{chain}_aa`` as aliases both
# matched, so both columns rendered when both existed, and (3) the
# ``full_{chain}_aa`` summary columns listed alongside their
# constituent parts, duplicating the chain content. All three are
# fixed here (#65).
#
# The canonical column names ``assemble.py`` writes to
# ``full_sequences.csv`` are all lowercase ``vdj_*_aa``. The
# ``full_*_aa`` summary columns are intentionally not rendered: their
# content equals leader + VDJ + constant, which are already shown.
_TCR_SEQUENCE_COLUMN_CANDIDATES: tuple[tuple[str, str], ...] = (
    ("beta_leader_aa", "Beta Leader"),
    ("vdj_beta_aa", "Beta VDJ"),
    ("beta_constant_aa", "Beta Constant"),
    ("linker", "Linker"),
    ("alpha_leader_aa", "Alpha Leader"),
    ("vdj_alpha_aa", "Alpha VDJ"),
    ("alpha_constant_aa", "Alpha Constant"),
)


def _default_tcr_sequence_columns(df: pd.DataFrame) -> dict[str, str]:
    """Auto-detect sequence columns for the TCR PDF in assembly order.

    Returns a dict mapping column-name → display label, biological
    assembly order preserved (β leader → β VDJ → β constant → linker
    → α leader → α VDJ → α constant). When none of the canonical
    assembled columns are present but ``single_chain_aa`` is, returns
    just that as a fallback so callers with a pre-concatenated single
    chain still get a useful page.
    """
    detected = {
        col: label for col, label in _TCR_SEQUENCE_COLUMN_CANDIDATES if col in df.columns
    }
    if not detected and "single_chain_aa" in df.columns:
        return {"single_chain_aa": "Single Chain"}
    return detected


def _expand_annotation_lines(lines, *, max_lines: int = 18) -> list[str]:
    """Lay out per-clone annotation lines for the sequence PDF (#202).

    Makes the block readable: each ``;``-separated item gets its own indented
    line (so one condition per line, not a long run-on), and all text is run
    through :func:`tcrsift.format.pdf_safe` so superscript ``⁺``/``⁻`` and other
    non-WinAnsi glyphs render as ASCII instead of missing-glyph boxes. A line
    of the form ``"label: a ; b ; c"`` becomes a ``label:`` header followed by
    indented ``a`` / ``b`` / ``c``. Capped at ``max_lines`` with an explicit
    ``... (+N more)`` marker rather than a silent truncation.
    """
    from .format import pdf_safe

    out: list[str] = []
    for raw in lines:
        line = pdf_safe(raw).rstrip()
        head, sep, rest = line.partition(":")
        if sep and ";" in rest:
            out.append(f"{head}:")
            out.extend(
                f"    {p.strip()}" for p in rest.split(";") if p.strip()
            )
        elif ";" in line and not sep:
            out.extend(p.strip() for p in line.split(";") if p.strip())
        else:
            out.append(line)
    if len(out) > max_lines:
        kept = out[: max_lines - 1]
        kept.append(f"... (+{len(out) - (max_lines - 1)} more)")
        return kept
    return out


def create_tcr_sequence_pdf(
    df: pd.DataFrame,
    output_path: str | Path,
    sequence_columns: dict[str, str] | None = None,
    title_column: str | None = None,
    sequence_font_size: int = 14,
    label_font_size: int = 11,
    title_font_size: int = 12,
    chars_per_line: int = 60,
    strict: bool | str = True,
    annotations: dict | None = None,
    annotation_key_column: str = "CDR3ab",
):
    """
    Create a PDF with color-coded TCR sequences.

    Each TCR is displayed on a separate page with:
    - TCR identifier and metadata
    - Color-coded sequence showing different regions
    - Color legend

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with TCR sequence data
    output_path : str or Path
        Path for output PDF
    sequence_columns : dict, optional
        Mapping of column names to display labels for sequence parts.
        Default: beta_leader, beta_VDJ, beta_constant, linker, alpha_leader, alpha_VDJ, alpha_constant
    title_column : str, optional
        Column to use for TCR title (default: auto-detect)
    sequence_font_size : int
        Font size for sequences
    label_font_size : int
        Font size for labels
    title_font_size : int
        Font size for titles
    chars_per_line : int
        Characters per line before wrapping
    strict : bool or str
        How to handle :func:`tcrsift.assemble.validate_sequences`
        load-bearing failures (CDR3 missing, canonical C-terminus
        mismatch, length window, etc.):

        - ``True`` (default) — autocorrect β J→C parity mismatches
          in-place (locus rule overrides CellRanger TRBC ambiguity;
          #89), then raise on any remaining load-bearing failures.
        - ``"skip"`` — autocorrect J/C parity, drop rows that still
          have load-bearing failures from the PDF, and render the rest.
        - ``False`` — render every row with a warning banner; no
          autocorrect or filtering.

        Defaults to ``True`` so users don't accidentally export PDFs
        of structurally-invalid TCRs (#68) while still tolerating the
        common upstream-CellRanger TRBC ambiguity case (#89).
    """
    # Sanity-gate the input. The PDF rendered nonsense in earlier
    # versions because nothing checked that the sequences it was
    # iterating over were biologically valid (#68). With #66 and #67
    # fixed, this is the last line of defense against future
    # construction regressions silently producing broken figures.
    from .assemble import validate_sequences

    apply_fix = strict is True or strict == "skip"
    qc_messages = validate_sequences(df, strict=False, fix=apply_fix)
    autocorrect_msgs = [m for m in qc_messages if m.severity == "autocorrect"]
    for m in autocorrect_msgs:
        logger.warning("create_tcr_sequence_pdf: %s", m)

    load_bearing = [m for m in qc_messages if m.severity == "load_bearing"]

    if load_bearing:
        preview = "\n  ".join(load_bearing[:10])
        more = f"\n  ... ({len(load_bearing) - 10} more)" if len(load_bearing) > 10 else ""
        msg = (
            f"create_tcr_sequence_pdf: {len(load_bearing)} load-bearing "
            f"validation failures:\n  {preview}{more}"
        )
        if strict == "skip":
            # Each ValidationMessage carries its source row index, so
            # filtering is a simple set membership check (no string
            # parsing of the message text).
            failing_labels = {m.idx for m in load_bearing if m.idx is not None}
            n_before = len(df)
            df = df.loc[~df.index.isin(failing_labels)]
            logger.warning(
                "create_tcr_sequence_pdf: dropping %d / %d clones with "
                "load-bearing validation failures (strict='skip'); "
                "rendering remaining %d.\n%s",
                n_before - len(df), n_before, len(df), msg,
            )
            if df.empty:
                from .validation import TCRsiftValidationError

                raise TCRsiftValidationError(
                    "create_tcr_sequence_pdf: every clone failed "
                    "validation (strict='skip' filtered all rows); "
                    "nothing to render.\n" + msg,
                    hint="Inspect the validation failures above. "
                    "Pass strict=False to render anyway with a warning "
                    "banner, or fix the upstream assembly.",
                )
        elif strict:
            from .validation import TCRsiftValidationError

            raise TCRsiftValidationError(
                msg,
                hint="Pass strict='skip' to drop failing clones and "
                "render the rest, or strict=False to render every row "
                "with a warning banner.",
            )
        else:
            logger.warning(msg)
    try:
        from itertools import cycle

        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfgen import canvas
    except ImportError:
        logger.warning("reportlab not installed, cannot generate sequence PDF")
        return

    if sequence_columns is None:
        sequence_columns = _default_tcr_sequence_columns(df)

    if not sequence_columns:
        logger.warning("No sequence columns found in DataFrame")
        return

    # Color palette for different regions
    color_list = [
        colors.HexColor("#1f77b4"),  # blue - beta leader
        colors.HexColor("#2ca02c"),  # green - beta VDJ
        colors.HexColor("#17becf"),  # cyan - beta constant
        colors.HexColor("#ff7f0e"),  # orange - linker
        colors.HexColor("#9467bd"),  # purple - alpha leader
        colors.HexColor("#d62728"),  # red - alpha VDJ
        colors.HexColor("#e377c2"),  # pink - alpha constant
        colors.HexColor("#7f7f7f"),  # gray
        colors.HexColor("#bcbd22"),  # olive
        colors.HexColor("#8c564b"),  # brown
    ]

    color_cycle = cycle(color_list)
    color_map = {}
    for col in sequence_columns.keys():
        color_map[col] = next(color_cycle)

    # Create PDF
    output_path = Path(output_path)
    c = canvas.Canvas(str(output_path), pagesize=letter)
    width, height = letter

    # Ensure Courier is available
    pdfmetrics.getFont("Courier")

    char_width = sequence_font_size * 0.6  # Monospace character width

    for idx, row in df.iterrows():
        y_position = height - 50

        def write_text(text, x=30, y_offset=18, font="Helvetica", size=None, color="black"):
            nonlocal y_position
            c.setFont(font, size or title_font_size)
            c.setFillColor(color)
            c.drawString(x, y_position, str(text))
            y_position -= y_offset

        def blank(space=25):
            nonlocal y_position
            y_position -= space

        # Title
        c.setFont("Helvetica-Bold", title_font_size + 2)
        c.setFillColor(colors.HexColor("#333333"))

        # Find a good title
        title = None
        if title_column and title_column in row:
            title = f"TCR: {row[title_column]}"
        elif "tcr_name" in row:
            title = f"TCR: {row['tcr_name']}"
        elif "CDR3ab" in row:
            title = f"Clone: {row['CDR3ab']}"
        else:
            title = f"TCR #{idx}"

        c.drawString(30, y_position, title)
        y_position -= 25

        # Metadata
        c.setFont("Helvetica", label_font_size)
        c.setFillColor("black")

        # CDR3 sequences
        if "CDR3_alpha" in row and pd.notna(row.get("CDR3_alpha")):
            write_text(f"CDR3 Alpha: {row['CDR3_alpha']}")
        if "CDR3_beta" in row and pd.notna(row.get("CDR3_beta")):
            write_text(f"CDR3 Beta: {row['CDR3_beta']}")

        # Gene information
        gene_cols = [
            c
            for c in row.index
            if c.endswith("_gene") or c.endswith("_v_gene") or c.endswith("_j_gene")
        ]
        for col in gene_cols[:6]:  # Limit to 6 genes
            if pd.notna(row.get(col)):
                write_text(f"{col}: {row[col]}", x=40)

        # Sequence length
        total_len = sum(
            len(str(row.get(col, ""))) for col in sequence_columns.keys() if pd.notna(row.get(col))
        )
        write_text(f"Total Length: {total_len} aa", x=40)

        blank(20)

        # Color legend
        legend_x = width - 180
        legend_y = height - 50
        c.setFont("Helvetica-Bold", label_font_size)
        c.setFillColor("black")
        c.drawString(legend_x, legend_y, "Legend:")
        legend_y -= 18

        c.setFont("Helvetica", label_font_size - 1)
        for col, label in sequence_columns.items():
            if col in color_map:
                c.setFillColor(color_map[col])
                c.drawString(legend_x, legend_y, f"  {label}")
                legend_y -= 15

        # Write sequence with color coding
        x_position = 30
        current_line_width = 0

        for col in sequence_columns.keys():
            sequence = row.get(col, "")
            if pd.isna(sequence) or not sequence:
                continue

            sequence = str(sequence)
            color = color_map.get(col, colors.black)

            for char in sequence:
                # Check for page break
                if y_position < 80:
                    c.showPage()
                    y_position = height - 50
                    x_position = 30
                    current_line_width = 0

                # Check for line wrap
                if current_line_width >= chars_per_line:
                    y_position -= sequence_font_size * 1.4
                    x_position = 30
                    current_line_width = 0

                c.setFont("Courier", sequence_font_size)
                c.setFillColor(color)
                c.drawString(x_position, y_position, char)
                x_position += char_width
                current_line_width += 1

        # Per-method evidence / selection-rule overlay (#123/#125). Drawn
        # at the bottom of the page when an annotation block is supplied
        # for this clone (keyed by ``annotation_key_column``).
        if annotations:
            key = row.get(annotation_key_column)
            ann_lines = annotations.get(key)
            if ann_lines is None and key is not None:
                ann_lines = annotations.get(str(key))
            if ann_lines:
                disp = _expand_annotation_lines(ann_lines)
                ann_y = 60 + 13 * min(len(disp), 18)
                c.setFillColor(colors.black)
                c.setFont("Helvetica-Bold", label_font_size)
                c.drawString(30, ann_y, "Selection:")
                ann_y -= 13
                c.setFont("Helvetica", label_font_size - 1)
                for line in disp:
                    c.drawString(36, ann_y, line)
                    ann_y -= 13

        c.showPage()

    c.save()
    logger.info(f"Generated TCR sequence PDF: {output_path}")


# =============================================================================
# Combined Report Generation
# =============================================================================


def generate_report(
    output_dir: str | Path,
    format: str = "pdf",
):
    """
    Generate combined report from all plots in output directory.

    Parameters
    ----------
    output_dir : str or Path
        Directory containing plot PNG files
    format : str
        Output format: "pdf" or "html"
    """
    output_dir = Path(output_dir)

    if format == "pdf":
        try:
            from reportlab.lib.pagesizes import letter
            from reportlab.lib.utils import ImageReader
            from reportlab.pdfgen import canvas

            pdf_path = output_dir / "tcrsift_report.pdf"
            c = canvas.Canvas(str(pdf_path), pagesize=letter)
            width, height = letter

            # Title page
            c.setFont("Helvetica-Bold", 24)
            c.drawCentredString(width / 2, height - 100, "TCRsift Analysis Report")
            c.showPage()

            # Add each plot
            for png_file in sorted(output_dir.glob("*.png")):
                img = ImageReader(str(png_file))
                img_width, img_height = img.getSize()

                # Scale to fit page
                scale = min((width - 100) / img_width, (height - 100) / img_height)
                scaled_width = img_width * scale
                scaled_height = img_height * scale

                x = (width - scaled_width) / 2
                y = (height - scaled_height) / 2

                c.drawImage(
                    img,
                    x,
                    y,
                    width=scaled_width,
                    height=scaled_height,
                )
                c.showPage()

            c.save()
            logger.info(f"Generated PDF report: {pdf_path}")

        except ImportError:
            logger.warning("reportlab not installed, cannot generate PDF report")

    elif format == "html":
        html_path = output_dir / "tcrsift_report.html"

        html_content = ["<html><head><title>TCRsift Report</title></head><body>"]
        html_content.append("<h1>TCRsift Analysis Report</h1>")

        for png_file in sorted(output_dir.glob("*.png")):
            title = png_file.stem.replace("_", " ").title()
            html_content.append(f"<h2>{title}</h2>")
            html_content.append(f'<img src="{png_file.name}" style="max-width:100%;">')
            html_content.append("<hr>")

        html_content.append("</body></html>")

        with open(html_path, "w") as f:
            f.write("\n".join(html_content))

        logger.info(f"Generated HTML report: {html_path}")
