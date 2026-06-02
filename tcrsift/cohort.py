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

"""Multi-donor cohort analysis (#123 / #125 increment 7).

Reads N per-donor ``tcrsift run`` output directories and computes
cross-donor clone-overlap matrices (paired CDR3αβ and β-only), reusing
the tested overlap math in :func:`tcrsift.clonotype.compute_sample_overlap_matrices`
with donor as the grouping axis. This is the one genuinely multi-input
analysis that can't fold into a per-donor command.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .clonotype import compute_sample_overlap_matrices

# Candidate locations within a donor output dir, in priority order.
_CLONOTYPE_FILES = ("data/clonotypes.csv", "clonotypes.csv")
_SELECTED_FILES = ("data/selected_clones.csv", "selected_clones.csv")
_CELL_COUNT_COLS = ("cell_count", "cells", "n_cells")


def _first_existing(donor_dir: Path, names: tuple[str, ...]) -> Path | None:
    for name in names:
        p = donor_dir / name
        if p.exists():
            return p
    return None


def read_donor_clones(
    donor_dir: str | Path,
    *,
    selected_only: bool = False,
) -> pd.DataFrame:
    """Read one donor's clone table as ``(CDR3ab, cells, CDR3_beta)``.

    Looks for ``selected_clones.csv`` (when ``selected_only``) or
    ``clonotypes.csv`` under the donor dir (``data/`` or top level). The
    cell count is taken from the first of ``cell_count`` / ``cells`` /
    ``n_cells`` present, defaulting to 1 (presence) when none exist.
    """
    donor_dir = Path(donor_dir)
    names = _SELECTED_FILES if selected_only else _CLONOTYPE_FILES
    path = _first_existing(donor_dir, names)
    if path is None:
        raise FileNotFoundError(
            f"No {'selected_clones' if selected_only else 'clonotypes'}.csv "
            f"found under {donor_dir} (looked in {list(names)})."
        )
    df = pd.read_csv(path)
    if "CDR3ab" not in df.columns:
        raise ValueError(f"{path} has no CDR3ab column")

    out = pd.DataFrame({"CDR3ab": df["CDR3ab"].astype(str)})
    cells_col = next((c for c in _CELL_COUNT_COLS if c in df.columns), None)
    out["cells"] = (
        pd.to_numeric(df[cells_col], errors="coerce").fillna(0.0)
        if cells_col else 1.0
    )
    # β-only key: explicit CDR3_beta column, else parse from CDR3ab
    # (``alpha_beta`` convention from build_clone_sample_long).
    if "CDR3_beta" in df.columns:
        out["CDR3_beta"] = df["CDR3_beta"].astype(str)
    else:
        out["CDR3_beta"] = out["CDR3ab"].str.split("_").str[-1]
    return out


def build_cohort_long(
    donor_dirs: dict[str, str | Path],
    *,
    selected_only: bool = False,
) -> pd.DataFrame:
    """Concatenate donor clone tables into a ``(donor, …)`` long frame."""
    frames = []
    for donor, ddir in donor_dirs.items():
        d = read_donor_clones(ddir, selected_only=selected_only)
        d.insert(0, "donor", str(donor))
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["donor", "CDR3ab", "cells", "CDR3_beta"])
    return pd.concat(frames, ignore_index=True)


def cohort_overlap(
    cohort_long: pd.DataFrame,
    *,
    include_beta_only: bool = True,
) -> dict[str, pd.DataFrame]:
    """Donor × donor overlap matrices.

    Reuses :func:`compute_sample_overlap_matrices` with ``donor`` as the
    grouping axis. Returns ``jaccard`` / ``cell_weighted_jaccard`` on the
    paired CDR3αβ clone key, and (when ``include_beta_only``) the same on
    the β-only key — public β-chains shared across donors are a common
    cross-donor signal even when the paired clone differs.
    """
    result: dict[str, pd.DataFrame] = {}
    paired = compute_sample_overlap_matrices(
        cohort_long, clone_col="CDR3ab", sample_col="donor", cells_col="cells",
    )
    result["jaccard"] = paired["jaccard"]
    result["cell_weighted_jaccard"] = paired["cell_weighted_jaccard"]

    if include_beta_only and "CDR3_beta" in cohort_long.columns:
        beta = compute_sample_overlap_matrices(
            cohort_long, clone_col="CDR3_beta", sample_col="donor", cells_col="cells",
        )
        result["jaccard_beta_only"] = beta["jaccard"]
        result["cell_weighted_jaccard_beta_only"] = beta["cell_weighted_jaccard"]
    return result


def run_cohort_analysis(
    donor_dirs: dict[str, str | Path],
    output_dir: str | Path,
    *,
    selected_only: bool = False,
    include_beta_only: bool = True,
    emit_tables: bool = True,
    emit_plots: bool = False,
) -> dict[str, pd.DataFrame]:
    """Orchestrate the cohort overlap analysis and write outputs.

    Returns the matrices dict; writes a CSV per matrix when
    ``emit_tables`` and a heatmap PDF per matrix when ``emit_plots``
    (heatmaps are skipped gracefully if matplotlib is unavailable).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_long = build_cohort_long(donor_dirs, selected_only=selected_only)
    matrices = cohort_overlap(cohort_long, include_beta_only=include_beta_only)

    if emit_tables:
        for name, mat in matrices.items():
            mat.to_csv(output_dir / f"cohort_{name}.csv")

    if emit_plots:
        _emit_heatmaps(matrices, output_dir)

    return matrices


def _emit_heatmaps(matrices: dict[str, pd.DataFrame], output_dir: Path) -> None:
    """Render one heatmap PDF per matrix; no-op if matplotlib missing."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - matplotlib optional
        return

    for name, mat in matrices.items():
        if mat.empty:
            continue
        fig, ax = plt.subplots(figsize=(1.2 + 0.6 * len(mat), 1.2 + 0.6 * len(mat)))
        im = ax.imshow(mat.to_numpy(dtype=float), vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(mat.columns)))
        ax.set_yticks(range(len(mat.index)))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right")
        ax.set_yticklabels(mat.index)
        ax.set_title(name)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(output_dir / f"cohort_{name}.pdf")
        plt.close(fig)
