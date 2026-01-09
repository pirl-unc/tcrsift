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

from pathlib import Path
from typing import Optional
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import anndata as ad

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
            from reportlab.pdfgen import canvas
            from reportlab.lib.utils import ImageReader

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
