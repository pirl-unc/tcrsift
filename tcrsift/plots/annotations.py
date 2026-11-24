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

"""Database annotation summary and granularity-aware structured annotation plots."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .state import (
    DB_CATEGORY_PALETTE,
    _category_color_strip,
    _parse_sample_frequencies,
    save_figure,
)

logger = logging.getLogger(__name__)


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
    one column per sample (zero-filled where a clone is absent). Accepts
    either Python dicts or stringified dicts produced by CSV round-trips.
    """
    if "sample_frequencies" not in clonotypes.columns or "CDR3ab" not in clonotypes.columns:
        return pd.DataFrame(index=clonotypes.get("CDR3ab", pd.Index([], name="CDR3ab")))
    parsed_rows = []
    idx_labels = []
    for cdr3, raw in zip(clonotypes["CDR3ab"].tolist(), clonotypes["sample_frequencies"].tolist()):
        d = _parse_sample_frequencies(raw)
        if d:
            parsed_rows.append(d)
            idx_labels.append(cdr3)
    if not parsed_rows:
        return pd.DataFrame(index=pd.Index([], name="CDR3ab"))
    expanded = pd.DataFrame(parsed_rows, index=pd.Index(idx_labels, name="CDR3ab"))
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


def plot_annotations(clonotypes: pd.DataFrame, output_dir: str | Path):
    """Generate annotation summary and granularity-aware phase-B plots."""
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

