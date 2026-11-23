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

"""Single-cell UMAP, signature-map, PRISM, method-overlap, and assembly plots."""

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

from ._plots_state import (
    _category_color_strip,
    _grid_shape,
    _per_cell_signature,
    save_figure,
)
from .format import pretty_method, pretty_methods

logger = logging.getLogger(__name__)

def _coerce_row_mask(
    mask: pd.Series | np.ndarray | Sequence[bool] | None,
    index: pd.Index,
    default: bool,
) -> np.ndarray:
    """Coerce a row mask to a boolean array aligned to ``index``.

    A :class:`pandas.Series` is reindexed (label-aligned, missing→False); an
    array/sequence is taken positionally; ``None`` fills with ``default``.
    """
    if mask is None:
        return np.full(len(index), default, dtype=bool)
    if isinstance(mask, pd.Series):
        return mask.reindex(index).fillna(False).to_numpy(dtype=bool)
    arr = np.asarray(mask, dtype=bool)
    if arr.shape[0] != len(index):
        raise ValueError(
            f"mask length {arr.shape[0]} != table length {len(index)}"
        )
    return arr


def plot_signature_map(
    table: pd.DataFrame,
    x_score: str,
    y_score: str,
    *,
    highlight_mask: pd.Series | np.ndarray | None = None,
    label_by: str | None = None,
    background_mask: pd.Series | np.ndarray | None = None,
    ax: plt.Axes | None = None,
    x_label: str | None = None,
    y_label: str | None = None,
    title: str | None = None,
    highlight_color: str = "#d62728",
    output_path: str | Path | None = None,
) -> plt.Axes | Path:
    """Two-score signature scatter with a highlighted subset (#302).

    The standard antigen-response (e.g. ``TNFRSF9``/``MKI67``) × cytolytic
    (``PRF1``/``GZMB``) signature map, generalized: plot any two per-clone
    score columns against each other, emphasize a chosen subset, and
    optionally color points by a category. Reused across the peptide-culture
    and TIL figures so the layout doesn't drift per paper.

    Three layers, drawn back-to-front:

    1. **background** — rows in ``background_mask``, faint gray context.
    2. **main, de-emphasized** — non-background rows not in
       ``highlight_mask``: light gray, or their ``label_by`` color.
    3. **main, highlighted** — non-background rows in ``highlight_mask``:
       larger, edged, on top, in ``highlight_color`` (or their ``label_by``
       color). ``highlight_mask=None`` emphasizes every main row, giving a
       plain colored scatter.

    This single primitive covers the issue's three views — *chosen vs
    background* (pass ``background_mask``), *chosen vs all* (just
    ``highlight_mask``), and *chosen + overlap pool* (``highlight_mask`` for
    the picks plus ``label_by`` for specificity/viral category) — by varying
    the masks, not the code. Pair it with
    :func:`tcrsift.gex.annotate_chosen`, which produces ``highlight_mask``.

    Parameters
    ----------
    table : pd.DataFrame
        Per-clone table containing ``x_score``, ``y_score`` and, if used,
        ``label_by``.
    x_score, y_score : str
        Numeric score columns for the x / y axes.
    highlight_mask : pd.Series | np.ndarray | None
        Rows to emphasize (e.g. from :func:`tcrsift.gex.annotate_chosen`). A
        Series is label-aligned to ``table.index``; an array is positional.
        ``None`` emphasizes all non-background rows.
    label_by : str | None
        Categorical column to color points by (legend rendered). When None,
        highlighted points use ``highlight_color`` and the rest light gray.
    background_mask : pd.Series | np.ndarray | None
        Rows to draw as faint gray context beneath everything.
    ax : matplotlib.axes.Axes | None
        Axis to draw into (for composing into a multi-panel figure). When
        None a new figure/axis is created.
    x_label, y_label, title : str | None
        Axis labels / title. Axis labels default to the column name with
        underscores turned to spaces.
    highlight_color : str
        Accent color for highlighted points when ``label_by`` is None.
    output_path : str | Path | None
        When given *and* ``ax`` is None, save via :func:`save_figure` and
        return the written path. Otherwise the Axes is returned.

    Returns
    -------
    matplotlib.axes.Axes | Path
        The Axes drawn into, or the saved path when ``output_path`` is given
        and a figure was created here.
    """
    for col in (x_score, y_score):
        if col not in table.columns:
            raise ValueError(f"plot_signature_map: missing score column {col!r}")
    if label_by is not None and label_by not in table.columns:
        raise ValueError(f"plot_signature_map: missing label_by column {label_by!r}")

    bg = _coerce_row_mask(background_mask, table.index, default=False)
    hl = _coerce_row_mask(highlight_mask, table.index, default=True)
    main = ~bg

    created_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.figure

    x = pd.to_numeric(table[x_score], errors="coerce").to_numpy()
    y = pd.to_numeric(table[y_score], errors="coerce").to_numpy()

    # Layer 1: background context.
    if bg.any():
        ax.scatter(
            x[bg], y[bg], s=14, color="#dddddd", alpha=0.6,
            linewidth=0, zorder=1, label="_nolegend_",
        )

    cat_color: dict[str, tuple] = {}
    if label_by is not None:
        cats = sorted(
            table.loc[main, label_by].dropna().astype(str).unique()
        )
        palette = sns.color_palette("tab10", n_colors=max(len(cats), 1))
        cat_color = {c: palette[i % len(palette)] for i, c in enumerate(cats)}
        labels = table[label_by].astype(str)

        def _colors_for(rows: np.ndarray) -> list:
            return [cat_color.get(labels.iloc[i], "#bdbdbd")
                    for i in np.flatnonzero(rows)]

        main_dim = main & ~hl
        if main_dim.any():
            ax.scatter(
                x[main_dim], y[main_dim], s=22, c=_colors_for(main_dim),
                alpha=0.55, linewidth=0, zorder=2,
            )
        main_hi = main & hl
        if main_hi.any():
            ax.scatter(
                x[main_hi], y[main_hi], s=90, c=_colors_for(main_hi),
                alpha=0.95, edgecolor="black", linewidth=0.6, zorder=3,
            )
        # Category legend (proxy handles) + optional emphasis note.
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", markersize=8,
                       markerfacecolor=cat_color[c], markeredgecolor="none",
                       label=c)
            for c in cats
        ]
        if highlight_mask is not None and (main & hl).any():
            handles.append(
                plt.Line2D([0], [0], marker="o", linestyle="", markersize=9,
                           markerfacecolor="white", markeredgecolor="black",
                           markeredgewidth=0.9, label="highlighted")
            )
        if handles:
            ax.legend(handles=handles, title=label_by, fontsize=8,
                      title_fontsize=9, loc="best", framealpha=0.85)
    else:
        main_dim = main & ~hl
        if main_dim.any():
            ax.scatter(
                x[main_dim], y[main_dim], s=20, color="#c7c7c7", alpha=0.6,
                linewidth=0, zorder=2,
                label="other" if (main & hl).any() else "_nolegend_",
            )
        main_hi = main & hl
        if main_hi.any():
            ax.scatter(
                x[main_hi], y[main_hi], s=90, color=highlight_color, alpha=0.95,
                edgecolor="white", linewidth=0.6, zorder=3,
                label="highlighted" if (main & ~hl).any() else "_nolegend_",
            )
        if (main & hl).any() and (main & ~hl).any():
            ax.legend(fontsize=8, loc="best", framealpha=0.85)

    ax.set_xlabel(x_label if x_label is not None else x_score.replace("_", " "))
    ax.set_ylabel(y_label if y_label is not None else y_score.replace("_", " "))
    if title:
        ax.set_title(title)
    ax.grid(True, which="both", linewidth=0.3, alpha=0.4)

    if output_path is not None and created_fig:
        return save_figure(fig, output_path)
    return ax


# =============================================================================
# AnnData atlas plots (#315) — the first AnnData-input members of plots.py.
#
# These consume the single-cell atlas (embed_cells → annotate_cells, #311/#312):
# an AnnData with obsm["X_umap"], obs["leiden"]/obs["phenotype"], per-cell
# signature-score columns, and TCR joins (CDR3ab / sample). They reuse
# save_figure / set_plot_format / set_polished_style — not a parallel stack —
# and route every antigen legend through pretty_antigen so a reactive pool is
# never shown by number alone.
# =============================================================================


def _violin_horizontal(ax, dataset, positions, widths):
    """``ax.violinplot`` drawn horizontally, across the matplotlib 3.7–4.0 pin.

    matplotlib 3.10 renamed the ``vert=False`` flag to
    ``orientation="horizontal"`` and deprecated ``vert``; older pinned versions
    only know ``vert``. Pick the kwarg the installed version accepts.
    """
    import matplotlib

    kw = dict(positions=positions, showextrema=False, widths=widths)
    if tuple(int(p) for p in matplotlib.__version__.split(".")[:2]) >= (3, 10):
        return ax.violinplot(dataset, orientation="horizontal", **kw)
    return ax.violinplot(dataset, vert=False, **kw)


def _umap_coords(adata: ad.AnnData, basis: str = "X_umap") -> np.ndarray:
    """The 2-D embedding from ``adata.obsm[basis]`` as an (n, 2) float array."""
    if basis not in adata.obsm:
        raise ValueError(
            f"{basis!r} not in adata.obsm — run embed_cells first "
            f"(have: {list(adata.obsm)})"
        )
    coords = np.asarray(adata.obsm[basis])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(f"adata.obsm[{basis!r}] must be (n, >=2); got {coords.shape}")
    return coords[:, :2].astype(float)


def _resolve_cell_values(
    adata: ad.AnnData, key: str
) -> tuple[pd.Series, str, str]:
    """Resolve a per-cell coloring ``key`` to (values, kind, label).

    ``key`` is looked up first in ``adata.obs`` (a phenotype/leiden category or
    a continuous signature-score column), then in ``adata.var_names`` (mean
    log1p expression of that one gene). ``kind`` is ``"categorical"`` for
    object/string/category/bool obs columns and ``"continuous"`` for numeric columns
    or a gene. Missing categorical values use the display label ``"nan"``.
    ``label`` is a display string.
    """
    if key in adata.obs.columns:
        s = adata.obs[key]
        if (
            isinstance(s.dtype, pd.CategoricalDtype)
            or pd.api.types.is_string_dtype(s.dtype)
            or pd.api.types.is_bool_dtype(s.dtype)
        ):
            # astype(str) preserves missing values under pandas 3, which would
            # silently drop those cells when matching the string color keys.
            return s.astype("string").fillna("nan"), "categorical", key
        return pd.to_numeric(s, errors="coerce"), "continuous", key
    if key in adata.var_names:
        vec = _per_cell_signature(adata, [key])
        return pd.Series(vec, index=adata.obs_names), "continuous", f"{key} (log1p)"
    raise ValueError(
        f"color {key!r} not found in adata.obs columns or var_names"
    )


def _categorical_colors(
    categories: Sequence[str],
    *,
    order: Sequence[str] | None = None,
    grayscale_prefixes: Sequence[str] | None = None,
) -> tuple[list[str], dict[str, tuple]]:
    """Order categories and build a color map, greying out de-emphasized ones.

    Categories whose label starts with any of ``grayscale_prefixes`` (e.g.
    ``"(unspecified)"`` / ``"(activated)"`` from the annotator's catch-all
    labels) are mapped to a neutral grey and dropped from the color legend so
    the specific phenotypes carry the palette. Returns ``(ordered, color_map)``
    where ``color_map`` covers every category (grey for the de-emphasized).
    """
    uniq = list(dict.fromkeys(str(c) for c in categories))
    if order is not None:
        ordered = [c for c in order if c in uniq]
        ordered += [c for c in uniq if c not in set(ordered)]
    else:
        ordered = sorted(uniq)
    prefixes = tuple(grayscale_prefixes or ())
    grey = [c for c in ordered if prefixes and c.startswith(prefixes)]
    colored = [c for c in ordered if c not in set(grey)]
    if len(colored) <= 20:
        palette = sns.color_palette("tab20", n_colors=max(len(colored), 1))
    else:
        palette = sns.color_palette("husl", n_colors=len(colored))
    color_map: dict[str, tuple] = {
        c: palette[i % len(palette)] for i, c in enumerate(colored)
    }
    for c in grey:
        color_map[c] = (0.78, 0.78, 0.78)
    return ordered, color_map


def _resolve_subset_mask(
    adata: ad.AnnData, subset: str | None, by: str | None
) -> np.ndarray:
    """Boolean per-cell mask for the highlighted subset.

    ``subset`` (a boolean obs column) wins when given; otherwise the subset is
    every cell with a non-empty ``by`` value (so ``by="peptide"`` implies the
    antigen-assigned cells). Empty tokens (``""``/``"nan"``/``"None"``) count as
    unassigned. Raises when neither is provided.
    """
    if subset is not None:
        if subset not in adata.obs.columns:
            raise ValueError(f"subset column {subset!r} not in adata.obs")
        return adata.obs[subset].fillna(False).to_numpy(dtype=bool)
    if by is not None:
        if by not in adata.obs.columns:
            raise ValueError(f"by column {by!r} not in adata.obs")
        vals = adata.obs[by]
        assigned = vals.notna() & ~vals.astype(str).isin(["", "nan", "None", "NaN", "<NA>"])
        return assigned.to_numpy(dtype=bool)
    raise ValueError(
        "plot needs a subset: pass subset=<bool obs column> or by=<obs column>"
    )


def plot_umap(
    adata: ad.AnnData,
    color: str,
    output_path: str | Path | None = None,
    *,
    basis: str = "X_umap",
    categorical_order: Sequence[str] | None = None,
    grayscale_prefixes: Sequence[str] | None = None,
    label_centroids: bool = True,
    cmap: str = "viridis",
    point_size: float = 6.0,
    vrange: tuple[float, float] | None = None,
    color_map: dict | None = None,
    centroids: dict | None = None,
    view_label: str | None = None,
    show_legend: bool = True,
    title: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes | Path:
    """UMAP scatter colored by a phenotype/cluster category or a continuous score.

    The standard atlas coloring (#308/#315). ``color`` names an ``obs`` column
    (a categorical phenotype/``leiden`` label, or a continuous per-cell signature
    score) or a gene in ``var_names`` (mean log1p expression). Categorical
    colorings draw a per-category legend and, with ``label_centroids``, place
    each category's name at its median position; ``grayscale_prefixes`` greys out
    de-emphasized labels (e.g. ``("(unspecified)", "(activated)")``) so the
    specific phenotypes carry the palette. Continuous colorings render a
    percentile-clipped colorbar.

    To keep a facet grid comparable across panels (#326), pin the layout/scale so
    every call renders identically: ``color_map`` (category→color) and
    ``centroids`` (category→(x, y)) for categorical, ``vrange`` (lo, hi) for
    continuous. ``view_label`` stamps a small view tag (e.g. ``"(T1)"``) in the
    corner; ``show_legend`` suppresses the per-panel legend/colorbar so a facet
    grid can carry one shared legend. :func:`plot_umap_facets` wires these up.

    Returns the saved path when ``output_path`` is given and a figure was created
    here, otherwise the Axes (so panels compose into a multi-panel figure).
    """
    coords = _umap_coords(adata, basis)
    values, kind, label = _resolve_cell_values(adata, color)

    created_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.figure

    if kind == "categorical":
        cats = values.to_numpy()
        prefixes = tuple(grayscale_prefixes or ())
        if color_map is None:
            ordered, color_map = _categorical_colors(
                cats, order=categorical_order, grayscale_prefixes=prefixes
            )
        else:
            # Pinned palette (facet grid): honor a requested order, then append
            # any categories present here but not in the order.
            present = list(dict.fromkeys(str(c) for c in cats))
            ordered = [c for c in (categorical_order or color_map) if c in color_map]
            ordered += [c for c in present if c not in set(ordered)]
        grey_cats = {c for c in ordered if prefixes and c.startswith(prefixes)}
        # Draw greyed (de-emphasized) categories underneath the colored ones.
        for c in list(grey_cats) + [c for c in ordered if c not in grey_cats]:
            m = cats == c
            if not m.any():
                continue
            is_grey = c in grey_cats
            ax.scatter(
                coords[m, 0], coords[m, 1], s=point_size,
                color=color_map.get(c, (0.78, 0.78, 0.78)),
                alpha=0.35 if is_grey else 0.8,
                linewidth=0, zorder=1 if is_grey else 2,
                label="_nolegend_" if is_grey else c,
            )
        if label_centroids:
            for c in ordered:
                if c in grey_cats:
                    continue
                if centroids is not None and c in centroids:
                    cx, cy = centroids[c]
                else:
                    m = cats == c
                    if not m.any():
                        continue
                    cx, cy = np.median(coords[m, 0]), np.median(coords[m, 1])
                ax.text(
                    cx, cy, c, fontsize=7, fontweight="bold", ha="center",
                    va="center", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              alpha=0.7, edgecolor="none"),
                )
        n_leg = sum(1 for c in ordered if c not in grey_cats)
        if show_legend and 0 < n_leg <= 25:
            ax.legend(
                markerscale=2.0, fontsize=7, loc="center left",
                bbox_to_anchor=(1.01, 0.5), framealpha=0.9,
                title=label, title_fontsize=8,
            )
    else:
        v = pd.to_numeric(values, errors="coerce").to_numpy()
        if vrange is not None:
            vmin, vmax = vrange
        else:
            finite = v[np.isfinite(v)]
            vmin, vmax = (
                (np.percentile(finite, 2), np.percentile(finite, 98))
                if finite.size else (0.0, 1.0)
            )
        order = np.argsort(np.nan_to_num(v, nan=-np.inf))  # high values on top
        sc = ax.scatter(
            coords[order, 0], coords[order, 1], c=v[order], s=point_size,
            cmap=cmap, vmin=vmin, vmax=max(vmax, vmin + 1e-9), linewidth=0,
            alpha=0.85,
        )
        if show_legend:
            cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
            cbar.set_label(label, fontsize=9)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title if title is not None else label)
    if view_label:
        ax.text(
            0.02, 0.98, view_label, transform=ax.transAxes, ha="left", va="top",
            fontsize=8, fontweight="bold", zorder=5,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7,
                      edgecolor="none"),
        )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)

    if output_path is not None and created_fig:
        return save_figure(fig, output_path)
    return ax


def _shared_umap_scale(adata, color, basis, *, categorical_order, grayscale_prefixes):
    """Compute the pinned color_map + centroids (categorical) or vrange
    (continuous) from the FULL atlas, so facet panels stay comparable (#326)."""
    coords = _umap_coords(adata, basis)
    values, kind, _ = _resolve_cell_values(adata, color)
    if kind == "categorical":
        cats = values.to_numpy()
        ordered, color_map = _categorical_colors(
            cats, order=categorical_order, grayscale_prefixes=grayscale_prefixes
        )
        prefixes = tuple(grayscale_prefixes or ())
        centroids = {}
        for c in ordered:
            if prefixes and c.startswith(prefixes):
                continue
            m = cats == c
            if m.any():
                centroids[c] = (float(np.median(coords[m, 0])),
                                float(np.median(coords[m, 1])))
        return kind, ordered, color_map, centroids, None
    v = pd.to_numeric(values, errors="coerce").to_numpy()
    finite = v[np.isfinite(v)]
    vrange = (
        (float(np.percentile(finite, 2)), float(np.percentile(finite, 98)))
        if finite.size else (0.0, 1.0)
    )
    return kind, None, None, None, vrange


def plot_umap_facets(
    adata: ad.AnnData,
    facet_col: str,
    color: str,
    output_path: str | Path | None = None,
    *,
    basis: str = "X_umap",
    shared_layout: bool = True,
    shared_vrange: bool = True,
    include_integrated: bool = True,
    view_labels: dict | None = None,
    integrated_label: str | None = None,
    categorical_order: Sequence[str] | None = None,
    grayscale_prefixes: Sequence[str] | None = None,
    cmap: str = "viridis",
    point_size: float = 6.0,
    max_cols: int = 3,
    title: str | None = None,
) -> Path | plt.Figure:
    """A grid of per-``facet_col`` UMAPs that share layout + color scale (#326).

    Loop-``plot_umap``-per-subset recomputes the palette / centroid layout / color
    range independently, so panels don't line up. This computes them ONCE from the
    full atlas and pins them across every panel: the same category→color mapping
    and centroid positions (``shared_layout``) and the same continuous vrange
    (``shared_vrange``). With ``include_integrated`` the first panel is the whole
    atlas; each panel gets a corner view stamp (``view_labels`` maps facet value →
    label, default ``"(<value>)"``; the integrated panel uses ``integrated_label``
    or ``"(<all>/<values>)"``). One shared legend/colorbar is drawn on the first
    panel only.
    """
    if facet_col not in adata.obs.columns:
        raise ValueError(f"facet_col {facet_col!r} not in adata.obs")
    # Keep missing annotations as a visible facet, just like missing colors.
    fvals = adata.obs[facet_col].astype("string").fillna("nan")
    groups = sorted(fvals.unique())

    kind = color_map = centroids = vrange = ordered = None
    if shared_layout or shared_vrange:
        kind, ordered, color_map, centroids, vrange = _shared_umap_scale(
            adata, color, basis,
            categorical_order=categorical_order,
            grayscale_prefixes=grayscale_prefixes,
        )
    if not shared_layout:
        color_map = centroids = None
    if not shared_vrange:
        vrange = None

    view_labels = view_labels or {}
    panels = []
    if include_integrated:
        integ = integrated_label or f"({'/'.join(groups)})"
        panels.append((adata, integ))
    for g in groups:
        panels.append((adata[(fvals == g).to_numpy()], view_labels.get(g, f"({g})")))

    # One shared legend/colorbar for the whole grid (from the pinned scale) beats
    # a per-panel legend that would jut into its neighbor.
    shared_key = (kind == "categorical" and color_map is not None) or (
        kind == "continuous" and vrange is not None
    )
    n_rows, n_cols = _grid_shape(len(panels), max_cols=max_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6.0 * n_cols, 5.0 * n_cols / max_cols + 0.5 * n_rows),
        squeeze=False,
    )
    flat = axes.flatten()
    for i, (panel_adata, vlabel) in enumerate(panels):
        plot_umap(
            panel_adata, color, ax=flat[i], basis=basis,
            categorical_order=ordered or categorical_order,
            grayscale_prefixes=grayscale_prefixes, cmap=cmap, point_size=point_size,
            color_map=color_map, centroids=centroids, vrange=vrange,
            view_label=vlabel, show_legend=(not shared_key and i == 0), title=" ",
        )
    for j in range(len(panels), len(flat)):
        flat[j].axis("off")

    if shared_key and kind == "categorical":
        prefixes = tuple(grayscale_prefixes or ())
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", markersize=8,
                       markerfacecolor=color_map[c], markeredgecolor="none", label=c)
            for c in ordered if not (prefixes and c.startswith(prefixes))
        ]
        if handles:
            fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
                       fontsize=8, framealpha=0.9, title=color, title_fontsize=9)
    elif shared_key and kind == "continuous":
        sm = plt.cm.ScalarMappable(
            norm=plt.Normalize(vrange[0], max(vrange[1], vrange[0] + 1e-9)), cmap=cmap
        )
        sm.set_array([])
        fig.colorbar(sm, ax=flat[:len(panels)].tolist(), fraction=0.02, pad=0.02,
                     label=color)

    fig.suptitle(title if title is not None else f"{color} by {facet_col}",
                 fontsize=13, y=1.01)
    # An external shared legend/colorbar defeats tight_layout; save_figure's
    # bbox_inches="tight" already crops the final image cleanly.
    if not shared_key:
        fig.tight_layout()
    if output_path is not None:
        return save_figure(fig, output_path)
    return fig


def plot_provenance(
    adata: ad.AnnData,
    flag: str | None = None,
    output_path: str | Path | None = None,
    *,
    layers: Sequence[tuple] | None = None,
    basis: str = "X_umap",
    by: str | None = None,
    group_colors: dict[str, str] | None = None,
    antigen_labels: dict | None = None,
    highlight_color: str = "#d62728",
    background_color: str = "#dddddd",
    point_size: float = 7.0,
    title: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes | Path:
    """Highlight a TCR subset over a grey UMAP background (the provenance panel).

    ``flag`` is a boolean ``obs`` column marking the subset to emphasize (e.g.
    tetramer/culture-validated or antigen-specific cells); every other cell is
    drawn as faint grey context. With ``by`` (an ``obs`` column such as the
    antigen/peptide source) the subset is recolored by that column, its legend
    routed through :func:`pretty_antigen` (never a bare pool number) and honoring
    a ``group_colors`` override keyed by the *pretty* label. ``antigen_labels``
    is an optional per-call ``condition → antigen`` map.

    Pass ``layers=[(flag_col, color, label), ...]`` to overlay MULTIPLE subsets on
    one panel with distinct colors (#326) — e.g. CD4-antigen-specific (blue) and
    CD8-antigen-specific (red) stars on the same UMAP. Background is every cell in
    no layer; later layers draw on top. ``flag`` and ``layers`` are mutually
    exclusive.
    """
    from .format import pretty_antigen

    coords = _umap_coords(adata, basis)
    if (flag is None) == (layers is None):
        raise ValueError("plot_provenance needs exactly one of flag= or layers=")

    created_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.figure

    # Multi-layer overlay: grey background under N distinctly-colored subsets.
    if layers is not None:
        masks = []
        for spec in layers:
            fcol = spec[0]
            if fcol not in adata.obs.columns:
                raise ValueError(f"layer flag column {fcol!r} not in adata.obs")
            masks.append(adata.obs[fcol].fillna(False).to_numpy(dtype=bool))
        any_hi = np.logical_or.reduce(masks) if masks else np.zeros(adata.n_obs, bool)
        ax.scatter(
            coords[~any_hi, 0], coords[~any_hi, 1], s=point_size * 0.6,
            color=background_color, alpha=0.5, linewidth=0, zorder=1,
            label="_nolegend_",
        )
        palette = sns.color_palette("tab10", n_colors=max(len(layers), 1))
        for i, (spec, m) in enumerate(zip(layers, masks)):
            col = spec[1] if len(spec) > 1 and spec[1] else palette[i % len(palette)]
            lbl = spec[2] if len(spec) > 2 and spec[2] else spec[0]
            ax.scatter(
                coords[m, 0], coords[m, 1], s=point_size * 3, color=col,
                alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3 + i,
                label=lbl,
            )
        ax.legend(markerscale=1.6, fontsize=7, loc="center left",
                  bbox_to_anchor=(1.01, 0.5), framealpha=0.9)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(title if title is not None else "provenance")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        if output_path is not None and created_fig:
            return save_figure(fig, output_path)
        return ax

    if flag not in adata.obs.columns:
        raise ValueError(f"flag column {flag!r} not in adata.obs")
    hi = adata.obs[flag].fillna(False).to_numpy(dtype=bool)

    ax.scatter(
        coords[~hi, 0], coords[~hi, 1], s=point_size * 0.6,
        color=background_color, alpha=0.5, linewidth=0, zorder=1,
        label="_nolegend_",
    )
    if by is None:
        ax.scatter(
            coords[hi, 0], coords[hi, 1], s=point_size * 3, color=highlight_color,
            alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3,
            label=flag,
        )
    else:
        if by not in adata.obs.columns:
            raise ValueError(f"by column {by!r} not in adata.obs")
        raw = adata.obs[by].astype(str).to_numpy()
        sub_raw = raw[hi]
        sub_coords = coords[hi]
        tokens = sorted(t for t in set(sub_raw) if t not in ("", "nan", "None"))
        pretty = {t: pretty_antigen(t, labels=antigen_labels) for t in tokens}
        palette = sns.color_palette("tab10", n_colors=max(len(tokens), 1))
        default_color = {
            t: palette[i % len(palette)] for i, t in enumerate(tokens)
        }
        for t in tokens:
            m = sub_raw == t
            plabel = pretty[t]
            col = (group_colors or {}).get(plabel, default_color[t])
            ax.scatter(
                sub_coords[m, 0], sub_coords[m, 1], s=point_size * 3, color=col,
                alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3,
                label=plabel,
            )
        if tokens:
            ax.legend(
                markerscale=1.6, fontsize=7, loc="center left",
                bbox_to_anchor=(1.01, 0.5), framealpha=0.9,
                title=by, title_fontsize=8,
            )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title if title is not None else f"provenance: {flag}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)

    if output_path is not None and created_fig:
        return save_figure(fig, output_path)
    return ax


def plot_signature_vs_background(
    adata: ad.AnnData,
    group_col: str,
    signature_cols: Sequence[str],
    output_path: str | Path | None = None,
    *,
    subset: str | None = None,
    by: str | None = None,
    antigen_labels: dict | None = None,
    group_colors: dict[str, str] | None = None,
    subset_color: str = "#d62728",
    max_cols: int = 3,
    title: str | None = None,
) -> Path | plt.Figure:
    """Selected/antigen-specific subset vs. per-lineage background, one signature per panel.

    For each column in ``signature_cols`` (per-cell score columns in ``obs``) a
    panel shows the background distribution of every cell grouped by ``group_col``
    (the lineage/phenotype background) as boxes, with the small highlighted subset
    overlaid as jittered dots at its own lineage's position — so you can read where
    the antigen-specific cells sit relative to the bulk of their lineage. The
    subset is ``subset`` (a boolean obs column) or, when omitted, every cell with a
    non-empty ``by`` value. With ``by`` the subset dots are recolored by that
    column via :func:`pretty_antigen`.
    """
    from .format import pretty_antigen

    signature_cols = list(signature_cols)
    if group_col not in adata.obs.columns:
        raise ValueError(f"group_col {group_col!r} not in adata.obs")
    missing = [c for c in signature_cols if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"signature_cols not in adata.obs: {missing}")
    if not signature_cols:
        raise ValueError("signature_cols is empty")

    sub_mask = _resolve_subset_mask(adata, subset, by)
    groups = adata.obs[group_col].astype(str)
    order = sorted(groups.unique())
    gpos = {g: i for i, g in enumerate(order)}

    # Antigen coloring of the subset dots.
    if by is not None:
        by_raw = adata.obs[by].astype(str).to_numpy()
        tokens = sorted(
            {by_raw[i] for i in np.flatnonzero(sub_mask)}
            - {"", "nan", "None", "NaN", "<NA>"}
        )
        pretty = {t: pretty_antigen(t, labels=antigen_labels) for t in tokens}
        palette = sns.color_palette("tab10", n_colors=max(len(tokens), 1))
        tok_color = {
            t: (group_colors or {}).get(pretty[t], palette[i % len(palette)])
            for i, t in enumerate(tokens)
        }

    n_rows, n_cols = _grid_shape(len(signature_cols), max_cols=max_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5.0 * n_cols, 4.0 * n_rows), squeeze=False
    )
    flat_axes = axes.flatten()
    rng = np.random.default_rng(0)

    for panel_i, sig in enumerate(signature_cols):
        ax = flat_axes[panel_i]
        vals = pd.to_numeric(adata.obs[sig], errors="coerce").to_numpy()
        # Background boxes: one per lineage over all its cells (empty groups
        # draw nothing at their position, which boxplot tolerates).
        box_data = [vals[(groups == g).to_numpy() & np.isfinite(vals)] for g in order]
        ax.boxplot(
            box_data, positions=range(len(order)), widths=0.6,
            showfliers=False, patch_artist=True,
            boxprops=dict(facecolor="#eeeeee", edgecolor="#999999"),
            medianprops=dict(color="#666666"),
            whiskerprops=dict(color="#999999"), capprops=dict(color="#999999"),
        )
        # Subset dots at each lineage's x, jittered.
        sub_idx = np.flatnonzero(sub_mask & np.isfinite(vals))
        gx = np.array([gpos.get(groups.iloc[i], np.nan) for i in sub_idx], dtype=float)
        keep = np.isfinite(gx)
        sub_idx, gx = sub_idx[keep], gx[keep]
        jitter = (rng.random(len(sub_idx)) - 0.5) * 0.35
        if by is not None:
            colors = [tok_color.get(by_raw[i], subset_color) for i in sub_idx]
        else:
            colors = subset_color
        ax.scatter(
            gx + jitter, vals[sub_idx], s=26, c=colors, alpha=0.9,
            edgecolor="white", linewidth=0.4, zorder=3,
        )
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(sig.replace("_", " "), fontsize=9)
        ax.set_title(sig.replace("_", " "), fontsize=10)
        ax.grid(True, axis="y", linewidth=0.3, alpha=0.4)

    for j in range(len(signature_cols), len(flat_axes)):
        flat_axes[j].axis("off")

    if by is not None and tokens:
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", markersize=7,
                       markerfacecolor=tok_color[t], markeredgecolor="white",
                       label=pretty[t])
            for t in tokens
        ]
        fig.legend(
            handles=handles, title=by, fontsize=8, title_fontsize=9,
            loc="center left", bbox_to_anchor=(1.0, 0.5), framealpha=0.9,
        )
    fig.suptitle(
        title if title is not None
        else f"Subset vs {group_col} background",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()

    if output_path is not None:
        return save_figure(fig, output_path)
    return fig


def plot_raincloud(
    groups: Sequence,
    values: Sequence[float],
    output_path: str | Path | None = None,
    *,
    order: Sequence[str] | None = None,
    palette: dict[str, str] | str | None = None,
    title: str | None = None,
    value_label: str | None = None,
    ax: plt.Axes | None = None,
) -> plt.Axes | Path:
    """Raincloud (half-violin + jittered raw dots + median/IQR) per group.

    Built for the bimodal / zero-inflated per-cell signatures (e.g. MANAscore),
    where a box or bar hides the two modes: each group gets a half-violin density
    above its row, the raw per-cell dots jittered below, and a median dot with an
    IQR bar. ``groups`` and ``values`` are tidy per-observation parallel arrays.
    """
    g = pd.Series(np.asarray(groups)).astype(str).to_numpy()
    v = pd.to_numeric(pd.Series(np.asarray(values)), errors="coerce").to_numpy()
    finite = np.isfinite(v)
    g, v = g[finite], v[finite]

    if order is not None:
        cats = [c for c in order if c in set(g)]
    else:
        cats = sorted(set(g))
    if isinstance(palette, dict):
        color_of = lambda c: palette.get(c, "#4C72B0")  # noqa: E731
    else:
        pal = sns.color_palette(palette or "Set2", n_colors=max(len(cats), 1))
        _pmap = {c: pal[i % len(pal)] for i, c in enumerate(cats)}
        color_of = lambda c: _pmap[c]  # noqa: E731

    created_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 0.9 * len(cats) + 2))
    else:
        fig = ax.figure

    rng = np.random.default_rng(0)
    for i, c in enumerate(cats):
        cv = v[g == c]
        if cv.size == 0:
            continue
        col = color_of(c)
        # Half-violin above the row (needs variance; skip for tiny/constant groups).
        if cv.size >= 2 and np.ptp(cv) > 0:
            parts = _violin_horizontal(ax, [cv], positions=[i], widths=0.9)
            for body in parts["bodies"]:
                verts = body.get_paths()[0].vertices
                verts[:, 1] = np.clip(verts[:, 1], i, np.inf)  # upper half only
                body.set_facecolor(col)
                body.set_edgecolor("none")
                body.set_alpha(0.5)
        # Raw dots jittered below the row.
        jitter = i - 0.18 - rng.random(cv.size) * 0.18
        ax.scatter(cv, jitter, s=10, color=col, alpha=0.5, linewidth=0, zorder=2)
        # Median dot + IQR bar.
        med = np.median(cv)
        q1, q3 = np.percentile(cv, [25, 75])
        ax.plot([q1, q3], [i - 0.05, i - 0.05], color="black", linewidth=1.5, zorder=3)
        ax.scatter([med], [i - 0.05], s=30, color="black", zorder=4)

    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats)
    ax.set_ylim(-0.7, len(cats) - 0.3)
    ax.set_xlabel(value_label or "value")
    if title:
        ax.set_title(title)
    ax.grid(True, axis="x", linewidth=0.3, alpha=0.4)

    if output_path is not None and created_fig:
        return save_figure(fig, output_path)
    return ax


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
        fallback_note = " (install 'upsetplot' for UpSet view)"
    except Exception as error:  # noqa: BLE001 - any failure inside the optional library
        # upsetplot 0.9.0 fills missing dot styles with chained in-place
        # ``fillna``, which pandas 3 ignores, so its matrix reaches matplotlib
        # with NaN colors. Draw the bar chart instead of failing the report.
        plt.close("all")
        logger.warning("UpSet plot failed (%s: %s); drawing a bar chart of overlaps instead",
                       type(error).__name__, error)
        fallback_note = ""

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

