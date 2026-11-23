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
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ._plots_atlas import (  # noqa: F401
    _FREQ_CBAR_W,
    _FREQ_CELL_H,
    _FREQ_CELL_W,
    _FREQ_HEATMAP_PAGE_ROWS,
    _FREQ_LABEL_W,
    _FREQ_STRIP_W,
    _FREQ_XFONT,
    _FREQ_YFONT,
    _categorical_colors,
    _coerce_row_mask,
    _render_frequency_heatmap_page,
    _resolve_cell_values,
    _resolve_subset_mask,
    _shared_umap_scale,
    _umap_coords,
    _violin_horizontal,
    plot_assembly,
    plot_cross_donor_venn,
    plot_freq_prism_grid,
    plot_freq_prism_scatter,
    plot_leader_summary,
    plot_method_overlap,
    plot_method_recovery,
    plot_prism_components,
    plot_provenance,
    plot_raincloud,
    plot_selected_frequency_heatmap,
    plot_set_overlap,
    plot_signature_map,
    plot_signature_vs_background,
    plot_til,
    plot_umap,
    plot_umap_facets,
    plot_vgene_usage_by_method,
)
from ._plots_reports import (  # noqa: F401
    _TCR_SEQUENCE_COLUMN_CANDIDATES,
    FUNNEL_LABEL_NICE,
    _default_tcr_sequence_columns,
    _draw_funnel_side_labels,
    _expand_annotation_lines,
    _log_stage_widths,
    create_pipeline_funnel,
    create_tcr_sequence_pdf,
    generate_report,
    normalize_funnel_label,
    plot_funnel,
    plot_funnel_lollipop,
    plot_funnel_ribbon,
    plot_funnel_terrace,
    plot_funnel_variants,
)
from ._plots_state import (  # noqa: F401
    DB_CATEGORY_PALETTE,
    VALID_PLOT_FORMATS,
    _category_color_strip,
    _grid_shape,
    _per_cell_signature,
    get_plot_format,
    save_figure,
    set_plot_format,
    set_polished_style,
)
from .format import pretty_method, pretty_methods, pretty_sample, pretty_samples  # noqa: F401, E402

# Re-export canonical T-cell signatures so callers of the per-sample
# scatter can keep importing them from ``tcrsift.plots``. Source of
# truth is :mod:`tcrsift.signatures`.
from .signatures import (  # noqa: F401, E402
    ACTIVATION_GENES_HGNC,
    AIM_GENES_HGNC,
    ANTIGEN_RESPONSE_GENES_HGNC,
    CYTOLYTIC_GENES_HGNC,
    EXHAUSTION_GENES_HGNC,
    EXPANSION_CORE_GENES_HGNC,
    TUMOR_REACTIVE_GENES_HGNC,
)

logger = logging.getLogger(__name__)

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
    if adata.n_obs == 0:
        logger.info("plot_qc: empty AnnData; skipping.")
        return
    obs = adata.obs
    samples = obs["sample"].dropna().unique() if "sample" in obs.columns else np.array(["all"])

    # Reads per cell distribution
    if "n_counts" in obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        for sample in samples:
            sample_data = (
                obs[obs["sample"] == sample]["n_counts"]
                if "sample" in obs.columns
                else obs["n_counts"]
            ).dropna()
            if len(sample_data) > 0:
                ax.hist(sample_data, bins=50, alpha=0.5, label=pretty_sample(sample))
        ax.set_xlabel("Total Counts per Cell")
        ax.set_ylabel("Number of Cells")
        ax.set_title("Read Count Distribution by Sample")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        if (obs["n_counts"].dropna() > 0).any():
            ax.set_xscale("log")
        save_figure(fig, output_dir / "qc_read_counts.png")

    # Genes detected distribution
    if "n_genes" in obs.columns:
        valid_pairs = [
            (s, (obs[obs["sample"] == s]["n_genes"] if "sample" in obs.columns else obs["n_genes"]).dropna().values)
            for s in samples
        ]
        valid_pairs = [(s, vals) for s, vals in valid_pairs if len(vals) > 0]
        if valid_pairs:
            fig, ax = plt.subplots(figsize=(12, 6))
            v_samples, data = zip(*valid_pairs)
            ax.violinplot(data, positions=range(len(v_samples)))
            ax.set_xticks(range(len(v_samples)))
            ax.set_xticklabels(pretty_samples(v_samples), rotation=45, ha="right")
            ax.set_ylabel("Genes Detected")
            ax.set_title("Gene Detection by Sample")
            save_figure(fig, output_dir / "qc_genes_detected.png")

    # Mitochondrial percentage
    mt_col = "percent_mt" if "percent_mt" in obs.columns else ("pct_counts_mt" if "pct_counts_mt" in obs.columns else None)
    if mt_col is not None:
        valid_pairs = [
            (s, (obs[obs["sample"] == s][mt_col] if "sample" in obs.columns else obs[mt_col]).dropna().values)
            for s in samples
        ]
        valid_pairs = [(s, vals) for s, vals in valid_pairs if len(vals) > 0]
        if valid_pairs:
            fig, ax = plt.subplots(figsize=(12, 6))
            v_samples, data = zip(*valid_pairs)
            ax.violinplot(data, positions=range(len(v_samples)))
            ax.set_xticks(range(len(v_samples)))
            ax.set_xticklabels(pretty_samples(v_samples), rotation=45, ha="right")
            ax.set_ylabel("Mitochondrial %")
            ax.set_title("Mitochondrial Content by Sample")
            ax.axhline(y=8, color="red", linestyle="--", alpha=0.5, label="Max threshold")
            ax.axhline(y=2, color="orange", linestyle="--", alpha=0.5, label="Min threshold")
            ax.legend()
            save_figure(fig, output_dir / "qc_mito_percent.png")

    plot_cells_per_sample_stacked(adata, output_dir / "qc_cells_per_sample.png")


# =============================================================================
# Phenotype Plots
# =============================================================================


def plot_phenotype(adata: ad.AnnData, output_dir: str | Path):
    """Generate phenotype plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if adata.n_obs == 0:
        logger.info("plot_phenotype: empty AnnData; skipping.")
        return

    # CD4 vs CD8 scatter
    if "CD4" in adata.obs.columns and "CD8" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 10))

        if "Tcell_type" in adata.obs.columns:
            for tcell_type in adata.obs["Tcell_type"].dropna().unique():
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
        tt_series = adata.obs["Tcell_type"]
        tcell_types = (
            list(tt_series.cat.categories)
            if isinstance(tt_series.dtype, pd.CategoricalDtype)
            else sorted(tt_series.dropna().astype(str).unique())
        )
        samples = (
            adata.obs["sample"].dropna().unique()
            if "sample" in adata.obs.columns
            else ["all"]
        )
        data = []
        for sample in samples:
            sample_data = (
                adata.obs[adata.obs["sample"] == sample]
                if "sample" in adata.obs.columns
                else adata.obs
            )
            total = len(sample_data)
            if total == 0:
                continue
            row = {"sample": sample}
            for tt in tcell_types:
                row[tt] = (sample_data["Tcell_type"] == tt).sum() / total * 100
            data.append(row)

        if data and tcell_types:
            fig, ax = plt.subplots(figsize=(14, 8))
            df_plot = pd.DataFrame(data).set_index("sample")
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
    if clonotypes is None or len(clonotypes) == 0 or "cell_count" not in clonotypes.columns:
        logger.info("plot_clonotypes: empty or missing cell_count; skipping.")
        return

    # Clone size distribution (rank plot)
    fig, ax = plt.subplots(figsize=(12, 6))
    sorted_counts = pd.to_numeric(clonotypes["cell_count"], errors="coerce").dropna().sort_values(ascending=False).values
    if len(sorted_counts) == 0:
        plt.close(fig)
        return
    ax.plot(range(1, len(sorted_counts) + 1), sorted_counts, linewidth=2)
    ax.set_xlabel("Clone Rank")
    ax.set_ylabel("Cell Count")
    ax.set_title("Clone Size Distribution (Rank Plot)")
    ax.set_xscale("log")
    if (sorted_counts > 0).any():
        ax.set_yscale("log")
    ax.fill_between(range(1, len(sorted_counts) + 1), sorted_counts, alpha=0.3)
    save_figure(fig, output_dir / "clonotype_rank_plot.png")

    # Clone size histogram. cell_count is float under attribution (#176), so
    # ceil to an int range for integer histogram bins.
    fig, ax = plt.subplots(figsize=(10, 6))
    raw_max = float(pd.to_numeric(clonotypes["cell_count"], errors="coerce").max() or 1.0)
    max_size = max(1, int(min(50, np.ceil(raw_max))))
    ax.hist(
        pd.to_numeric(clonotypes["cell_count"], errors="coerce").fillna(0).clip(upper=max_size),
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

    # Sample sharing matrix (exact token match on `;`-delimited sample names, not regex substring)
    if "samples" in clonotypes.columns and "CDR3ab" in clonotypes.columns:
        parsed_sample_sets = [
            {tok.strip() for tok in str(val).split(";") if tok.strip()}
            if pd.notna(val)
            else set()
            for val in clonotypes["samples"]
        ]
        samples = sorted(set().union(*parsed_sample_sets)) if parsed_sample_sets else []

        if len(samples) > 1:
            cdr3_vals = clonotypes["CDR3ab"].tolist()
            sample_clones = {
                sample: {
                    cdr3
                    for cdr3, sset in zip(cdr3_vals, parsed_sample_sets)
                    if sample in sset and pd.notna(cdr3)
                }
                for sample in samples
            }

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
    if "is_viral" in clonotypes.columns and len(clonotypes) > 0:
        matched_mask = clonotypes["db_match"].fillna(False).astype(bool)
        viral_mask = matched_mask & clonotypes["is_viral"].fillna(False).astype(bool)
        matched = int(matched_mask.sum())
        viral = int(viral_mask.sum())
        non_viral_matched = max(0, matched - viral)
        unmatched = max(0, len(clonotypes) - matched)

        sizes = [viral, non_viral_matched, unmatched]
        if sum(sizes) > 0:
            fig, ax = plt.subplots(figsize=(8, 8))
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
    return save_figure(fig, out_path)


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
    return save_figure(fig, out_path)

