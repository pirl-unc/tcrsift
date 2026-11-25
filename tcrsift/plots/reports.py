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

"""Pipeline funnel plots, color-coded TCR sequence PDFs, and figure report generation."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .state import save_figure

logger = logging.getLogger(__name__)

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
    if not stage_counts:
        return output_dir / filename

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

    return save_figure(fig, output_dir / filename)


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
    """Map stage counts to relative widths in [floor, 1.0]."""
    if not counts:
        return np.array([], dtype=float)
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
    if not stage_counts:
        return output_dir / filename

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
    return save_figure(fig, out_path)


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
    if not stage_counts:
        return output_dir / filename

    stages = list(stage_counts.keys())
    raw_counts = [int(c) for c in stage_counts.values()]
    counts = np.array([max(c, 1) for c in raw_counts], dtype=float)
    n = len(stages)

    fig, ax = plt.subplots(figsize=(10.5, max(4, 0.6 * n + 2)))
    y = np.arange(n)
    ax.plot(counts, y, color="#b91c1c", linewidth=2.0, zorder=2)
    ax.scatter(
        counts, y,
        s=np.linspace(95, 42, n),
        color="#b91c1c", edgecolor="white", linewidth=1.0, zorder=3,
    )
    for yi, stage, count, raw_count in zip(y, stages, counts, raw_counts):
        ax.hlines(yi, xmin=count, xmax=counts.max() * 2.4,
                  color="#fda4a4", linewidth=2.3)
        ax.text(count / 1.10, yi, f"n = {raw_count:,}",
                va="center", ha="left", fontsize=9.5)

    ax.set_xscale("log")
    ax.set_xlim(counts.max() * 2.4, counts.min() * 0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(stages)
    ax.set_xlabel("count (log scale)")
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_title(title, pad=14)

    out_path = output_dir / filename
    return save_figure(fig, out_path)


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
    if not stage_counts:
        return output_dir / filename

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
    return save_figure(fig, out_path)


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
    if not stage_counts:
        return []
    paths: list[Path] = []
    # The legacy bars renderer writes to ``funnel_plot.<ext>``; capture its
    # returned Path so plot_format="pdf"/"svg" returns the vector path too.
    paths.append(plot_funnel(stage_counts, output_dir, title=title))
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
# The canonical column names ``tcrsift.assemble`` writes to
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
    from ..format import pdf_safe

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
    species: str | None = None,
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
    from ..assemble import validate_sequences

    apply_fix = strict is True or strict == "skip"
    qc_messages = validate_sequences(df, strict=False, fix=apply_fix, species=species)
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
                from ..validation import TCRsiftValidationError

                raise TCRsiftValidationError(
                    "create_tcr_sequence_pdf: every clone failed "
                    "validation (strict='skip' filtered all rows); "
                    "nothing to render.\n" + msg,
                    hint="Inspect the validation failures above. "
                    "Pass strict=False to render anyway with a warning "
                    "banner, or fix the upstream assembly.",
                )
        elif strict:
            from ..validation import TCRsiftValidationError

            raise TCRsiftValidationError(
                msg,
                hint="Pass strict='skip' to drop failing clones and "
                "render the rest, or strict=False to render every row "
                "with a warning banner.",
            )
        else:
            logger.warning(msg)
    from itertools import cycle

    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfgen import canvas

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

    from ..format import pdf_safe

    # Ensure Courier is available
    pdfmetrics.getFont("Courier")

    char_width = sequence_font_size * 0.6  # Monospace character width

    for idx, row in df.iterrows():
        y_position = height - 50

        def write_text(text, x=30, y_offset=18, font="Helvetica", size=None, color="black"):
            nonlocal y_position
            c.setFont(font, size or title_font_size)
            c.setFillColor(color)
            c.drawString(x, y_position, pdf_safe(str(text)))
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

        c.drawString(30, y_position, pdf_safe(title))
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
        # below the sequence (or on a fresh page if space is tight) so it
        # never overlaps long amino-acid/nucleotide sequences.
        if annotations:
            key = row.get(annotation_key_column)
            ann_lines = annotations.get(key)
            if ann_lines is None and key is not None:
                ann_lines = annotations.get(str(key))
            if ann_lines:
                disp = _expand_annotation_lines(ann_lines)
                ann_needed = 13 * (min(len(disp), 18) + 1)
                ann_y = min(y_position - 25, 60 + 13 * min(len(disp), 18))
                if y_position - 25 < 60 + ann_needed:
                    c.showPage()
                    ann_y = height - 50
                c.setFillColor(colors.black)
                c.setFont("Helvetica-Bold", label_font_size)
                c.drawString(30, ann_y, "Selection:")
                ann_y -= 13
                c.setFont("Helvetica", label_font_size - 1)
                for line in disp:
                    if ann_y < 40:
                        c.showPage()
                        ann_y = height - 50
                        c.setFillColor(colors.black)
                        c.setFont("Helvetica", label_font_size - 1)
                    c.drawString(36, ann_y, pdf_safe(line))
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
    *,
    report_name: str | None = None,
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
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.utils import ImageReader
        from reportlab.pdfgen import canvas

        # all-figures.pdf carries the donor/run name on the cover so a
        # per-donor bundle is identifiable on its own (#262 follow-up).
        pdf_path = output_dir / "all-figures.pdf"
        c = canvas.Canvas(str(pdf_path), pagesize=letter)
        width, height = letter

        # Title page
        c.setFont("Helvetica-Bold", 24)
        c.drawCentredString(width / 2, height - 100, "TCRsift — all figures")
        if report_name:
            c.setFont("Helvetica", 16)
            c.drawCentredString(width / 2, height - 140, str(report_name))
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

    elif format == "html":
        import html as _html

        html_path = output_dir / "all-figures.html"

        safe_name = f" — {_html.escape(str(report_name))}" if report_name else ""
        html_content = [f"<html><head><title>TCRsift figures{safe_name}</title></head><body>"]
        html_content.append(f"<h1>TCRsift — all figures{safe_name}</h1>")

        for png_file in sorted(output_dir.glob("*.png")):
            title = _html.escape(png_file.stem.replace("_", " ").title())
            src = _html.escape(png_file.name, quote=True)
            html_content.append(f"<h2>{title}</h2>")
            html_content.append(f'<img src="{src}" style="max-width:100%;">')
            html_content.append("<hr>")

        html_content.append("</body></html>")

        with open(html_path, "w") as f:
            f.write("\n".join(html_content))

        logger.info(f"Generated HTML report: {html_path}")
