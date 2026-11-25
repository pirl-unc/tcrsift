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

"""Core repertoire QC, phenotyping, clonotype, filter, Pgen, and assembly QC plots."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..format import pretty_sample, pretty_samples
from .ranks import plot_rank_family
from .state import save_figure

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
    elif publicness_col and publicness_col in df.columns:
        pub_series = pd.to_numeric(df.loc[values.index, publicness_col], errors="coerce")
        counts, bin_edges, patches = ax.hist(
            values, bins=bins, color="#3b82f6", edgecolor="white", alpha=0.85
        )
        valid_pub = pub_series.dropna()
        if not valid_pub.empty:
            vmin = float(min(0.0, valid_pub.min()))
            vmax = float(max(1.0, valid_pub.max()))
            if vmax <= vmin:
                vmax = vmin + 1.0
            norm = plt.Normalize(vmin=vmin, vmax=vmax)
            cmap = plt.cm.coolwarm
            vals_arr = values.to_numpy(dtype=float)
            pub_arr = pub_series.to_numpy(dtype=float)
            bin_idx = np.clip(
                np.digitize(vals_arr, bin_edges[:-1]) - 1, 0, len(patches) - 1
            )
            for idx_b, patch in enumerate(patches):
                in_bin = pub_arr[bin_idx == idx_b]
                in_bin = in_bin[np.isfinite(in_bin)]
                if len(in_bin) > 0:
                    patch.set_facecolor(cmap(norm(float(np.mean(in_bin)))))
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.02)
            cbar.set_label(publicness_col, fontsize=9)
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
    if clonotypes is None or clonotypes.empty:
        return

    # Frequency vs condition scatter (the key plot)
    if "max_frequency" in clonotypes.columns and "n_samples" in clonotypes.columns:
        fig, ax = plt.subplots(figsize=(12, 8))
        has_cell_count = "cell_count" in clonotypes.columns

        if "tier" in clonotypes.columns:
            tiers = sorted(clonotypes["tier"].dropna().unique())
            colors = plt.cm.viridis(np.linspace(0, 1, max(len(tiers), 1)))

            for tier, color in zip(tiers, colors):
                mask = clonotypes["tier"] == tier
                sizes = clonotypes.loc[mask, "cell_count"] * 5 if has_cell_count else 30
                ax.scatter(
                    clonotypes.loc[mask, "max_frequency"] * 100,
                    clonotypes.loc[mask, "n_samples"],
                    c=[color],
                    alpha=0.6,
                    s=sizes,
                    label=tier,
                )
            ax.legend(title="Tier", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            sizes = clonotypes["cell_count"] * 5 if has_cell_count else 30
            ax.scatter(
                clonotypes["max_frequency"] * 100,
                clonotypes["n_samples"],
                alpha=0.5,
                s=sizes,
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

