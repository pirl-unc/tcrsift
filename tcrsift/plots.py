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

logger = logging.getLogger(__name__)

# Set default style
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams["figure.facecolor"] = "#f8f9fa"


def save_figure(fig: plt.Figure, output_path: str | Path, dpi: int = 300):
    """Save figure with consistent settings."""
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    logger.info(f"Saved plot to {output_path}")


# =============================================================================
# QC Plots (load command)
# =============================================================================

def plot_qc(adata: ad.AnnData, output_dir: str | Path):
    """Generate QC plots for loaded data."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Reads per cell distribution
    if "n_counts" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        for sample in adata.obs["sample"].unique():
            sample_data = adata.obs[adata.obs["sample"] == sample]["n_counts"]
            ax.hist(sample_data, bins=50, alpha=0.5, label=sample)
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
        ax.set_xticklabels(samples, rotation=45, ha="right")
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
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_ylabel("Mitochondrial %")
        ax.set_title("Mitochondrial Content by Sample")
        ax.axhline(y=8, color="red", linestyle="--", alpha=0.5, label="Max threshold")
        ax.axhline(y=2, color="orange", linestyle="--", alpha=0.5, label="Min threshold")
        ax.legend()
        save_figure(fig, output_dir / "qc_mito_percent.png")

    # Cells per sample
    fig, ax = plt.subplots(figsize=(12, 6))
    sample_counts = adata.obs["sample"].value_counts()
    ax.bar(range(len(sample_counts)), sample_counts.values)
    ax.set_xticks(range(len(sample_counts)))
    ax.set_xticklabels(sample_counts.index, rotation=45, ha="right")
    ax.set_ylabel("Number of Cells")
    ax.set_title("Cells per Sample")
    for i, v in enumerate(sample_counts.values):
        ax.text(i, v + 0.01 * max(sample_counts.values), str(v), ha="center", fontsize=10)
    save_figure(fig, output_dir / "qc_cells_per_sample.png")


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

    # Clone size histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    max_size = min(50, clonotypes["cell_count"].max())
    ax.hist(clonotypes["cell_count"].clip(upper=max_size), bins=range(1, max_size + 2), edgecolor="black")
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
                sample_clones[sample] = set(clonotypes.loc[mask, "clone_id"])

            jaccard_matrix = np.zeros((len(samples), len(samples)))
            for i, s1 in enumerate(samples):
                for j, s2 in enumerate(samples):
                    intersection = len(sample_clones[s1] & sample_clones[s2])
                    union = len(sample_clones[s1] | sample_clones[s2])
                    jaccard_matrix[i, j] = intersection / union if union > 0 else 0

            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(
                jaccard_matrix,
                xticklabels=samples,
                yticklabels=samples,
                annot=True,
                fmt=".2f",
                cmap="viridis",
                ax=ax,
            )
            ax.set_title("Clone Sharing Between Samples (Jaccard Similarity)")
            plt.xticks(rotation=45, ha="right")
            save_figure(fig, output_dir / "clonotype_sharing_heatmap.png")


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
        ax.bar(range(len(tier_counts)), tier_counts.values, color=plt.cm.viridis(np.linspace(0, 1, len(tier_counts))))
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

        ax.bar(x - width/2, cd8_counts, width, label="CD8+")
        ax.bar(x + width/2, cd4_counts, width, label="CD4+")
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
        labels = [f"Viral ({viral})", f"Non-viral matched ({non_viral_matched})", f"Unmatched ({unmatched})"]
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
            recovery = tier_data["til_match"].sum() / len(tier_data) * 100 if len(tier_data) > 0 else 0
            recovery_rates.append(recovery)

        ax.bar(range(len(tiers)), recovery_rates, color=plt.cm.viridis(np.linspace(0, 1, len(tiers))))
        ax.set_xticks(range(len(tiers)))
        ax.set_xticklabels(tiers)
        ax.set_ylabel("TIL Recovery Rate (%)")
        ax.set_title("Culture→TIL Recovery by Confidence Tier")
        for i, v in enumerate(recovery_rates):
            ax.text(i, v + 1, f"{v:.1f}%", ha="center", fontsize=10)
        save_figure(fig, output_dir / "til_recovery_by_tier.png")

    # Culture vs TIL frequency scatter
    matched = matched_clonotypes[matched_clonotypes["til_match"]]
    if len(matched) > 0 and "max_frequency" in matched.columns and "til_frequency" in matched.columns:
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
                ax.axvline(lengths.median(), color="red", linestyle="--", label=f"Median: {lengths.median():.0f}")
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

def plot_funnel(
    stage_counts: dict[str, int],
    output_dir: str | Path,
    title: str = "TCR Selection Funnel",
):
    """
    Generate a funnel plot showing TCR counts at each filtering stage.

    Parameters
    ----------
    stage_counts : dict
        Ordered dictionary mapping stage names to counts.
        Example: {"Raw Cells": 10000, "With VDJ": 8000, "Phenotyped": 7500, ...}
    output_dir : str or Path
        Directory to save the plot
    title : str
        Plot title
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stages = list(stage_counts.keys())
    counts = list(stage_counts.values())
    max_count = max(counts) if counts else 1

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

    ax.set_xlim(0, 1.1)
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

    save_figure(fig, output_dir / "funnel_plot.png")


def create_pipeline_funnel(
    raw_cells: int,
    with_vdj: int,
    phenotyped: int,
    clonotypes: int,
    filtered: int,
    tier_counts: dict[str, int] | None = None,
    output_dir: str | Path = ".",
):
    """
    Create a funnel plot for the TCRsift pipeline stages.

    Parameters
    ----------
    raw_cells : int
        Number of cells after loading
    with_vdj : int
        Number of cells with VDJ data
    phenotyped : int
        Number of cells after phenotyping
    clonotypes : int
        Number of unique clonotypes
    filtered : int
        Number of clonotypes passing filters
    tier_counts : dict, optional
        Counts per confidence tier
    output_dir : str or Path
        Output directory for the plot
    """
    stage_counts = {
        "Raw Cells": raw_cells,
        "With VDJ": with_vdj,
        "Phenotyped (CD4/CD8)": phenotyped,
        "Unique Clonotypes": clonotypes,
        "Passing Filters": filtered,
    }

    if tier_counts:
        for tier, count in tier_counts.items():
            stage_counts[f"Tier: {tier}"] = count

    plot_funnel(stage_counts, output_dir)


# =============================================================================
# Color-Coded TCR Sequence PDF
# =============================================================================

def create_tcr_sequence_pdf(
    df: pd.DataFrame,
    output_path: str | Path,
    sequence_columns: dict[str, str] | None = None,
    title_column: str | None = None,
    sequence_font_size: int = 14,
    label_font_size: int = 11,
    title_font_size: int = 12,
    chars_per_line: int = 60,
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
    """
    try:
        from itertools import cycle

        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfgen import canvas
    except ImportError:
        logger.warning("reportlab not installed, cannot generate sequence PDF")
        return

    # Default sequence columns
    if sequence_columns is None:
        sequence_columns = {}
        # Try to find sequence columns in order
        column_candidates = [
            ("beta_leader_aa", "Beta Leader"),
            ("alpha_leader_aa", "Alpha Leader"),
            ("vdj_beta_aa", "Beta VDJ"),
            ("VDJ_beta_aa", "Beta VDJ"),
            ("full_beta_aa", "Beta Full"),
            ("beta_constant_aa", "Beta Constant"),
            ("linker", "Linker"),
            ("vdj_alpha_aa", "Alpha VDJ"),
            ("VDJ_alpha_aa", "Alpha VDJ"),
            ("full_alpha_aa", "Alpha Full"),
            ("alpha_constant_aa", "Alpha Constant"),
        ]
        for col, label in column_candidates:
            if col in df.columns:
                sequence_columns[col] = label

        # If we have single_chain_aa, use that instead
        if "single_chain_aa" in df.columns and not sequence_columns:
            sequence_columns = {"single_chain_aa": "Single Chain"}

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
        elif "clone_id" in row:
            title = f"Clone: {row['clone_id']}"
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
        gene_cols = [c for c in row.index if c.endswith("_gene") or c.endswith("_v_gene") or c.endswith("_j_gene")]
        for col in gene_cols[:6]:  # Limit to 6 genes
            if pd.notna(row.get(col)):
                write_text(f"{col}: {row[col]}", x=40)

        # Sequence length
        total_len = sum(len(str(row.get(col, ""))) for col in sequence_columns.keys() if pd.notna(row.get(col)))
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
                    x, y,
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
