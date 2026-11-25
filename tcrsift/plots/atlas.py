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
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .state import (
    _grid_shape,
    _per_cell_signature,
    save_figure,
)

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
# AnnData atlas plots (#315) — the first AnnData-input members of tcrsift.plots.
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
    if adata.n_obs == 0:
        logger.warning("plot_umap_facets: empty AnnData; skipping")
        return None
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
    if not panels:
        return None

    # One shared legend/colorbar for the whole grid (from the pinned scale) beats
    # a per-panel legend that would jut into its neighbor.
    shared_key = (kind == "categorical" and color_map is not None) or (
        kind == "continuous" and vrange is not None
    )
    n_rows, n_cols = _grid_shape(len(panels), max_cols=max_cols)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6.0 * n_cols, 5.0 * n_rows),
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
    from ..format import pretty_antigen

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
        tokens = sorted(t for t in set(sub_raw) if t not in ("", "nan", "NaN", "None", "<NA>"))
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
    from ..format import pretty_antigen

    signature_cols = list(signature_cols)
    if group_col not in adata.obs.columns:
        raise ValueError(f"group_col {group_col!r} not in adata.obs")
    missing = [c for c in signature_cols if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"signature_cols not in adata.obs: {missing}")
    if not signature_cols:
        raise ValueError("signature_cols is empty")
    if adata.n_obs == 0:
        logger.warning("plot_signature_vs_background: empty AnnData; skipping")
        return None

    sub_mask = _resolve_subset_mask(adata, subset, by)
    groups = adata.obs[group_col].astype(str)
    order = sorted(groups.unique())
    if not order:
        return None
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

