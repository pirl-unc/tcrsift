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

"""pMHC multimer / feature-barcode specificity extraction and calling.

Supports 10x Genomics Feature Barcode matrices (``Antibody Capture`` and
``Antigen Capture`` modalities, including dCODE Dextramer and BEAM-T libraries)
stored alongside gene expression in ``filtered_feature_bc_matrix.h5`` or
Matrix Market directories.
"""

from __future__ import annotations

import logging
import re
import warnings
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from .validation import TCRsiftValidationError, resolve_cellranger_gex_source

logger = logging.getLogger(__name__)

_DEFAULT_NC_RE = re.compile(
    r"(?:"
    r"(?:^|[_.-])(?:neg(?:ative)?[-_]?ctrl|neg(?:ative)?[-_]?control|blank)(?:$|[_.-])"
    r"|(?:^|[_.-])nc(?:$|[_.-](?!\d{5,})(?:total\w*|\d{0,3})?$)"
    r")",
    re.IGNORECASE,
)
_DEFAULT_HASHTAG_RE = re.compile(
    r"(?:^|[_.-])(?:fb_hash\d*|hashtag[-_]?\d*|hash\d*|hto[-_]?\d*|cmo[-_]?\d*|cellplex\d*)(?:$|[_.-])",
    re.IGNORECASE,
)
_DEFAULT_SURFACE_RE = re.compile(
    r"(?:totalseq|cite[-_]?seq|isotype|(?:mouse|rat|hamster)[-_]?igg\d*[ab]?|"
    r"(?:^|[_.-])fb_(?:cd\d|ccr\d|cxcr\d|hla|pd-?1|tigit|lag3|tim-?3|klrg1)|"
    r"^(?:cd\d{1,3}[a-z]?|ccr\d+|cxcr\d+|hla[-_]?dr|pd-?1|pd-?l1|tigit|lag-?3|tim-?3|"
    r"klrg1|tcrab|tcrgd|igg\d*[ab]?)(?:$|[_.-]))",
    re.IGNORECASE,
)


def extract_feature_barcodes(
    path: str | Path,
    *,
    feature_types: Sequence[str] = ("Antibody Capture", "Antigen Capture"),
) -> pd.DataFrame:
    """Load non-GEX feature-barcode UMI counts from a 10x ``.h5`` or MTX source.

    Parameters
    ----------
    path : str or Path
        CellRanger ``outs`` directory, sample-prefixed ``.h5`` directory,
        direct ``.h5`` file, or 10x Matrix Market directory.
    feature_types : sequence of str
        ``feature_types`` values to extract (default:
        ``("Antibody Capture", "Antigen Capture")``).

    Returns
    -------
    pd.DataFrame
        Cells × features float64 DataFrame of raw UMI counts, with
        ``df.attrs["feature_types"]`` mapping each column name to its 10x
        modality string. Returns an empty 0-column DataFrame indexed by cell
        barcode when no matching modality is present.
    """
    import scanpy as sc

    source_kind, resolved = resolve_cellranger_gex_source(path)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Variable names are not unique",
            category=UserWarning,
        )
        if source_kind == "h5":
            adata = sc.read_10x_h5(str(resolved), gex_only=False)
        else:
            adata = sc.read_10x_mtx(str(resolved), gex_only=False)

    if "feature_types" not in adata.var.columns:
        return pd.DataFrame(index=adata.obs_names.astype(str))

    wanted = {str(ft) for ft in feature_types}
    mask = adata.var["feature_types"].astype(str).isin(wanted).to_numpy()
    if not mask.any():
        return pd.DataFrame(index=adata.obs_names.astype(str))

    sub = adata[:, mask].copy()
    sub.var_names_make_unique()
    mat = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    cols = sub.var_names.astype(str)
    out = pd.DataFrame(
        np.asarray(mat, dtype=float),
        index=sub.obs_names.astype(str),
        columns=cols,
    )
    out.attrs["feature_types"] = dict(zip(cols, sub.var["feature_types"].astype(str)))
    return out


def classify_feature_columns(
    columns: Iterable[str] | pd.DataFrame,
    *,
    feature_types: Mapping[str, str] | None = None,
    nc_regex: str | re.Pattern[str] = _DEFAULT_NC_RE,
    hashtag_regex: str | re.Pattern[str] = _DEFAULT_HASHTAG_RE,
    surface_regex: str | re.Pattern[str] = _DEFAULT_SURFACE_RE,
) -> dict[str, list[str]]:
    """Partition non-GEX feature names into multimers, negative controls, surface antibodies, and hashtags.

    Parameters
    ----------
    columns : iterable of str or pd.DataFrame
        Feature names (or a DataFrame returned by :func:`extract_feature_barcodes`,
        whose ``attrs["feature_types"]`` will be consulted automatically).
    feature_types : mapping of str to str, optional
        Optional ``{feature_name: modality}`` map; features explicitly typed as
        ``"Antigen Capture"`` (10x BEAM-T / dCODE Antigen Capture) are routed to
        ``negative_controls`` (when matching ``nc_regex``) or ``multimers`` even
        if their label includes ``TotalSeq``.
    nc_regex : str or Pattern
        Pattern matching negative-control multimers.
    hashtag_regex : str or Pattern
        Pattern matching sample-multiplexing hashtags / CMOs.
    surface_regex : str or Pattern
        Pattern matching CITE-seq / TotalSeq surface protein antibodies.

    Returns
    -------
    dict of str to list of str
        Dictionary with keys ``"multimers"``, ``"negative_controls"``,
        ``"surface_proteins"``, and ``"hashtags"``.
    """
    if isinstance(columns, pd.DataFrame):
        if feature_types is None:
            feature_types = columns.attrs.get("feature_types")
        col_iter: Iterable[str] = columns.columns
    else:
        col_iter = columns

    nc_rx = re.compile(nc_regex) if isinstance(nc_regex, str) else nc_regex
    ht_rx = re.compile(hashtag_regex) if isinstance(hashtag_regex, str) else hashtag_regex
    sf_rx = re.compile(surface_regex) if isinstance(surface_regex, str) else surface_regex

    out: dict[str, list[str]] = {
        "multimers": [],
        "negative_controls": [],
        "surface_proteins": [],
        "hashtags": [],
    }
    for col in col_iter:
        name = str(col)
        modality = (feature_types or {}).get(name, "")
        if nc_rx.search(name):
            out["negative_controls"].append(name)
        elif ht_rx.search(name):
            out["hashtags"].append(name)
        elif modality == "Antigen Capture":
            out["multimers"].append(name)
        elif sf_rx.search(name):
            out["surface_proteins"].append(name)
        else:
            out["multimers"].append(name)
    return out


def call_multimer_specificity(
    counts_df: pd.DataFrame,
    *,
    multimer_cols: Sequence[str] | None = None,
    nc_cols: Sequence[str] | None = None,
    min_umi: float = 10.0,
    min_nc_ratio: float = 5.0,
    nc_pseudocount: float = 1.0,
    min_dominance_ratio: float = 2.0,
) -> pd.DataFrame:
    """Call per-cell pMHC multimer specificity calibrated against negative-control multimers.

    A cell is called specific for its top pMHC multimer when:

    1. ``multimer_umi >= min_umi``,
    2. ``multimer_umi / max(max_nc_umi, nc_pseudocount) >= min_nc_ratio``, and
    3. ``multimer_umi / max(second_multimer_umi, nc_pseudocount) >= min_dominance_ratio``
       (rejecting promiscuous / sticky poly-binders, which are flagged in
       ``is_multi_reactive``).

    Parameters
    ----------
    counts_df : pd.DataFrame
        Cells × features UMI count table from :func:`extract_feature_barcodes`.
    multimer_cols : sequence of str, optional
        Explicit pMHC multimer columns. Inferred via :func:`classify_feature_columns`
        when ``None``.
    nc_cols : sequence of str, optional
        Explicit negative-control multimer columns. Inferred via
        :func:`classify_feature_columns` when ``None``.
    min_umi : float
        Minimum UMI count on the top multimer (default ``10.0``).
    min_nc_ratio : float
        Minimum ratio over ``max(max_nc_umi, nc_pseudocount)`` (default ``5.0``).
    nc_pseudocount : float
        Pseudocount floor on negative-control and runner-up UMI denominators
        (default ``1.0``).
    min_dominance_ratio : float
        Minimum ratio of top multimer UMI over ``max(second_multimer_umi, nc_pseudocount)``
        (default ``2.0``).

    Returns
    -------
    pd.DataFrame
        Per-cell table indexed identically to ``counts_df`` with ``top_multimer``,
        ``multimer_call``, ``multimer_umi``, ``second_multimer_umi``, ``max_nc_umi``,
        ``nc_fold_change``, ``dominance_ratio``, ``is_multimer_specific``, and
        ``is_multi_reactive``.
    """
    if min_umi < 0 or min_nc_ratio <= 0 or nc_pseudocount <= 0 or min_dominance_ratio < 1.0:
        raise TCRsiftValidationError(
            "Invalid multimer calling threshold(s)",
            hint="Require min_umi >= 0, min_nc_ratio > 0, nc_pseudocount > 0, and min_dominance_ratio >= 1.0.",
        )

    if multimer_cols is None or nc_cols is None:
        classified = classify_feature_columns(counts_df)
        if multimer_cols is None:
            multimer_cols = classified["multimers"]
        if nc_cols is None:
            nc_cols = classified["negative_controls"]

    multimer_cols = [c for c in multimer_cols if c in counts_df.columns]
    nc_cols = [c for c in nc_cols if c in counts_df.columns]
    if not multimer_cols:
        raise TCRsiftValidationError(
            "No pMHC multimer columns found in feature-barcode DataFrame",
            hint=f"Available columns: {list(counts_df.columns)[:15]}",
        )
    if not nc_cols:
        logger.warning(
            "No negative-control multimer columns provided or detected; "
            "nc_fold_change will divide by nc_pseudocount (%g) only.",
            nc_pseudocount,
        )

    n_cells = len(counts_df)
    if n_cells == 0:
        return pd.DataFrame(
            columns=[
                "top_multimer",
                "multimer_call",
                "multimer_umi",
                "second_multimer_umi",
                "max_nc_umi",
                "nc_fold_change",
                "dominance_ratio",
                "is_multimer_specific",
                "is_multi_reactive",
            ],
            index=counts_df.index,
        )

    mat = counts_df[multimer_cols].fillna(0.0).to_numpy(dtype=float)
    nc_max = (
        counts_df[nc_cols].fillna(0.0).to_numpy(dtype=float).max(axis=1)
        if nc_cols
        else np.zeros(n_cells, dtype=float)
    )

    top_idx = np.argmax(mat, axis=1)
    top_umi = mat[np.arange(n_cells), top_idx]
    raw_top_names = np.asarray(multimer_cols, dtype=object)[top_idx]
    top_names = np.where(top_umi > 0, raw_top_names, None)

    if mat.shape[1] >= 2:
        part = np.partition(mat, -2, axis=1)
        second_umi = part[:, -2]
    else:
        second_umi = np.zeros(n_cells, dtype=float)

    nc_fc = top_umi / np.maximum(nc_max, nc_pseudocount)
    dom_ratio = top_umi / np.maximum(second_umi, nc_pseudocount)

    above_bg = (top_umi >= min_umi) & (nc_fc >= min_nc_ratio)
    is_specific = above_bg & (dom_ratio >= min_dominance_ratio)
    is_multi = above_bg & (second_umi >= min_umi) & (dom_ratio < min_dominance_ratio)

    calls = np.where(is_specific, top_names, None)
    return pd.DataFrame(
        {
            "top_multimer": top_names,
            "multimer_call": calls,
            "multimer_umi": top_umi,
            "second_multimer_umi": second_umi,
            "max_nc_umi": nc_max,
            "nc_fold_change": nc_fc,
            "dominance_ratio": dom_ratio,
            "is_multimer_specific": is_specific,
            "is_multi_reactive": is_multi,
        },
        index=counts_df.index,
    )


def summarize_clone_multimer_specificity(
    obs_df: pd.DataFrame,
    calls_df: pd.DataFrame | None = None,
    *,
    clone_col: str = "CDR3ab",
    min_specific_cells: int = 2,
    min_concordance: float = 0.5,
) -> pd.DataFrame:
    """Aggregate per-cell multimer calls onto clonotypes with clonal concordance QC.

    True pMHC-specific clonotypes bind the same multimer reproducibly across
    sibling cells of the clone (high ``multimer_concordance``), whereas sticky
    pMHC reagents bind isolated cells across unrelated clones at low per-clone
    concordance.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Per-cell metadata table containing ``clone_col`` (and optionally the
        output columns of :func:`call_multimer_specificity` when ``calls_df``
        is omitted).
    calls_df : pd.DataFrame, optional
        Per-cell output of :func:`call_multimer_specificity`. Must share unique
        cell barcode indices with ``obs_df`` (or match ``obs_df.index``
        positionally 1-to-1).
    clone_col : str
        Clonotype identifier column (default ``"CDR3ab"``; partial single-chain
        ``_beta`` and ``alpha_`` entries are excluded automatically when
        ``clone_col == "CDR3ab"``).
    min_specific_cells : int
        Minimum cells calling the modal multimer for ``multimer_confident``
        (default ``2``).
    min_concordance : float
        Minimum fraction of cells in the clone calling the modal multimer for
        ``multimer_confident`` (default ``0.5``).

    Returns
    -------
    pd.DataFrame
        One row per clonotype with ``clone_cells``, ``multimer_specific_cells``,
        ``multimer_call``, ``multimer_call_cells``, ``multimer_concordance``,
        ``multi_reactive_fraction``, ``mean_multimer_umi``,
        ``mean_nc_fold_change``, and ``multimer_confident``.
    """
    if min_specific_cells < 1 or not (0.0 <= min_concordance <= 1.0):
        raise TCRsiftValidationError(
            "Invalid clone multimer QC threshold(s)",
            hint="Require min_specific_cells >= 1 and 0.0 <= min_concordance <= 1.0.",
        )
    if clone_col not in obs_df.columns:
        raise TCRsiftValidationError(
            f"Missing clone identifier column {clone_col!r} in obs_df",
            hint=f"Available columns: {list(obs_df.columns)[:15]}",
        )

    call_fields = (
        "multimer_call",
        "multimer_umi",
        "max_nc_umi",
        "nc_fold_change",
        "is_multimer_specific",
        "is_multi_reactive",
    )
    work = obs_df[[clone_col]].copy()
    if calls_df is not None:
        if len(calls_df) == len(work) and calls_df.index.equals(work.index):
            for col in call_fields:
                if col in calls_df.columns:
                    work[col] = calls_df[col].to_numpy()
        else:
            if not calls_df.index.is_unique or not work.index.is_unique:
                raise TCRsiftValidationError(
                    "Duplicate cell barcodes prevent unambiguous index alignment between obs_df and calls_df",
                    hint="Prefix cell barcodes by sample (e.g. sample + '_' + barcode) or pass pre-aligned inputs.",
                )
            aligned = calls_df.reindex(work.index)
            for col in call_fields:
                if col in aligned.columns:
                    work[col] = aligned[col]
    else:
        for col in call_fields:
            if col in obs_df.columns:
                work[col] = obs_df[col]

    if "multimer_call" not in work.columns:
        raise TCRsiftValidationError(
            "Missing 'multimer_call' column; pass calls_df from call_multimer_specificity()."
        )

    c_str = work[clone_col].astype(str).str.strip()
    valid_clone = work[clone_col].notna() & ~c_str.isin({"", "_", "nan", "None", "<NA>"})
    if clone_col == "CDR3ab":
        valid_clone &= ~c_str.str.startswith("_") & ~c_str.str.endswith("_")
    work = work[valid_clone].copy()
    if work.empty:
        return pd.DataFrame(
            columns=[
                clone_col,
                "clone_cells",
                "multimer_specific_cells",
                "multimer_call",
                "multimer_call_cells",
                "multimer_concordance",
                "multi_reactive_fraction",
                "mean_multimer_umi",
                "mean_nc_fold_change",
                "multimer_confident",
            ]
        )

    work["_is_spec"] = work.get(
        "is_multimer_specific", work["multimer_call"].notna()
    ).fillna(False).astype(bool)
    work["_is_multi"] = work.get(
        "is_multi_reactive", pd.Series(False, index=work.index)
    ).fillna(False).astype(bool)

    rows = []
    for clone_id, grp in work.groupby(clone_col, sort=False, observed=True):
        n_total = len(grp)
        spec_grp = grp[grp["_is_spec"] & grp["multimer_call"].notna()]
        n_spec = len(spec_grp)
        if n_spec > 0:
            vc = spec_grp["multimer_call"].astype(str).value_counts()
            top_m = str(vc.index[0])
            top_m_cells = int(vc.iloc[0])
            concordance = top_m_cells / n_total
            top_sub = spec_grp[spec_grp["multimer_call"].astype(str) == top_m]
            mean_umi = (
                float(top_sub["multimer_umi"].mean())
                if "multimer_umi" in top_sub.columns
                else np.nan
            )
            mean_nc_fc = (
                float(top_sub["nc_fold_change"].mean())
                if "nc_fold_change" in top_sub.columns
                else np.nan
            )
        else:
            top_m = None
            top_m_cells = 0
            concordance = 0.0
            mean_umi = np.nan
            mean_nc_fc = np.nan

        rows.append(
            {
                clone_col: clone_id,
                "clone_cells": n_total,
                "multimer_specific_cells": n_spec,
                "multimer_call": top_m,
                "multimer_call_cells": top_m_cells,
                "multimer_concordance": concordance,
                "multi_reactive_fraction": float(grp["_is_multi"].mean()),
                "mean_multimer_umi": mean_umi,
                "mean_nc_fold_change": mean_nc_fc,
                "multimer_confident": bool(
                    top_m is not None
                    and top_m_cells >= min_specific_cells
                    and concordance >= min_concordance
                ),
            }
        )

    out = pd.DataFrame(rows)
    return out.sort_values(
        ["multimer_confident", "multimer_call_cells", "multimer_concordance"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
