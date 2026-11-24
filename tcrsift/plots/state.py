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

"""Shared Matplotlib/Seaborn style configuration and figure-saving state."""

from __future__ import annotations

import logging
import math
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


def get_plot_format() -> str:
    """Return the active module-level plot format (``'png'``, ``'pdf'``, ``'svg'``, or ``'both'``)."""
    return _PLOT_FORMAT


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


DB_CATEGORY_PALETTE: dict[str | None, str] = {
    "tumor_self": "#d62728",   # red — the cohort's target axis
    "viral": "#1f77b4",        # blue — bystander viral
    "bacterial": "#2ca02c",    # green
    "self": "#ff7f0e",         # orange — non-tumor self
    "other": "#7f7f7f",        # grey
    "unknown": "#bcbcbc",      # light grey
    None: "#e5e5e5",           # very light grey (unmatched)
}


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


def _grid_shape(n_panels: int, max_cols: int = 3) -> tuple[int, int]:
    """Pick (rows, cols) for a multi-panel layout."""
    cols = min(max(1, n_panels), max_cols)
    rows = math.ceil(n_panels / cols)
    return rows, cols


def _parse_sample_frequencies(val) -> dict[str, float]:
    """Parse a per-clone ``sample_frequencies`` cell (dict or stringified dict from CSV)."""
    import ast
    import json

    if isinstance(val, dict):
        parsed = val
    elif isinstance(val, str):
        s = val.strip()
        if not s or s.lower() in {"nan", "none", "{}"}:
            return {}
        try:
            parsed = json.loads(s)
        except Exception:
            try:
                parsed = ast.literal_eval(s)
            except Exception:
                return {}
        if not isinstance(parsed, dict):
            return {}
    else:
        return {}

    out: dict[str, float] = {}
    for k, v in parsed.items():
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if np.isfinite(fv):
            out[str(k)] = fv
    return out


def _per_cell_signature(adata: ad.AnnData, gene_ids: list[str]) -> np.ndarray:
    """Mean log1p expression across ``gene_ids`` per cell.

    Drops any IDs not present in ``adata.var_names``; returns a
    zero-vector if none of the requested genes are in the matrix.
    Safe against duplicate ``var_names`` and already-log-transformed
    or residual matrices.
    """
    var_names = adata.var_names
    idx = []
    for g in gene_ids:
        matches = np.flatnonzero(var_names == g)
        if len(matches) > 0:
            idx.append(int(matches[0]))
    if not idx:
        return np.zeros(adata.n_obs)
    X = adata.X[:, idx]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    already_logged = "log1p" in getattr(adata, "uns", {}) or (
        X.size > 0 and float(np.nanmin(X)) < 0
    )
    if already_logged:
        return np.nanmean(X, axis=1)
    return np.nanmean(np.log1p(np.clip(X, 0.0, None)), axis=1)


