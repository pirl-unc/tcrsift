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
TIL timepoint harmonization and clone prioritization.

This module powers the ``tcrsift til-select`` command. It is designed to
support the same input model as legacy ad-hoc scripts:

- per-timepoint ``consensus_annotations.<TP>.csv``
- per-timepoint ``clonotypes.<TP>.csv``
- per-timepoint ``filtered_contig_annotations.<TP>.csv``
- per-timepoint ``sample_filtered_feature_bc_matrix.<TP>.h5``

The workflow produces a harmonized table across timepoints, marker-expression
scores from 10x GEX, optional public-DB annotation, and v2-style selection
subset files.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy.sparse import issparse

from .annotate import annotate_clonotypes
from .validation import TCRsiftValidationError, validate_file_exists

logger = logging.getLogger(__name__)

MARKER_GENES_DEFAULT = ["CD4", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "MKI67", "TNFRSF9"]
CYTOLYTIC_GENES_DEFAULT = ("GZMB", "PRF1")
ANTIGEN_RESPONSE_GENES_DEFAULT = ("TNFRSF9", "MKI67")
RANK_METRICS = (
    "mean_frequency",
    "max_frequency",
    "total_cells",
    "marker_score_cp10k_mean",
    "marker_score_z_mean",
)

VDJ_SEGMENT_COLS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
VDJ_SEGMENT_NT_COLS = [f"{c}_nt" for c in VDJ_SEGMENT_COLS]


def _timepoint_sort_key(label: str) -> tuple[int, str]:
    """Sort key that keeps T1/T2/T10 in numeric order when possible."""
    m = re.search(r"(\d+)", str(label))
    if m:
        return (int(m.group(1)), str(label))
    return (10**9, str(label))


def _truthy_mask(series: pd.Series) -> pd.Series:
    """Convert a mixed-typed boolean-like series to bool mask."""
    if series.dtype == bool:
        return series.fillna(False)
    lowered = series.astype(str).str.lower()
    return lowered.isin({"true", "1", "yes", "y", "t"})


def _strip_barcode_suffix(barcode: str) -> str:
    """Strip gem-group style suffix (e.g. ``-1``) when present."""
    if "-" in barcode:
        return barcode.rsplit("-", 1)[0]
    return barcode


def parse_sample_args(items: list[str]) -> dict[str, dict[str, Path]]:
    """
    Parse ``--samples`` / ``--inputs`` entries.

    Expected item format:
    ``LABEL=CONSENSUS_PATH,CLONOTYPES_PATH`` (comma or semicolon separator).
    """
    parsed: dict[str, dict[str, Path]] = {}
    for item in items:
        if "=" not in item:
            raise TCRsiftValidationError(
                f"Invalid sample mapping: '{item}'",
                hint="Use LABEL=CONSENSUS_PATH,CLONOTYPES_PATH",
            )
        label, rest = item.split("=", 1)
        label = label.strip()
        if not label:
            raise TCRsiftValidationError(
                f"Invalid sample mapping: '{item}'",
                hint="Label before '=' cannot be empty.",
            )
        if label in parsed:
            raise TCRsiftValidationError(
                f"Duplicate timepoint label: '{label}'",
                hint="Use unique labels (e.g. T1, T2, T3).",
            )

        if "," in rest:
            consensus_raw, clonotypes_raw = rest.split(",", 1)
        elif ";" in rest:
            consensus_raw, clonotypes_raw = rest.split(";", 1)
        else:
            raise TCRsiftValidationError(
                f"Invalid sample mapping: '{item}'",
                hint="Expected CONSENSUS_PATH,CLONOTYPES_PATH after label.",
            )

        consensus_path = validate_file_exists(Path(consensus_raw.strip()), f"{label} consensus")
        clonotypes_path = validate_file_exists(Path(clonotypes_raw.strip()), f"{label} clonotypes")
        parsed[label] = {"consensus": consensus_path, "clonotypes": clonotypes_path}

    return dict(sorted(parsed.items(), key=lambda kv: _timepoint_sort_key(kv[0])))


def parse_config(path: str | Path) -> dict[str, dict[str, Path]]:
    """Parse YAML config mapping timepoints to consensus/clonotypes files."""
    path = validate_file_exists(path, "til-select config")
    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    if "timepoints" in raw:
        raw = raw["timepoints"]

    if not isinstance(raw, dict):
        raise TCRsiftValidationError(
            "Invalid til-select config format",
            hint="Config must be a mapping of label -> {consensus, clonotypes}.",
        )

    base = path.parent
    parsed: dict[str, dict[str, Path]] = {}
    for label, entry in raw.items():
        if not isinstance(entry, dict):
            raise TCRsiftValidationError(
                f"Invalid config entry for {label}",
                hint="Each label must map to a dict with consensus/clonotypes paths.",
            )
        consensus = entry.get("consensus")
        clonotypes = entry.get("clonotypes")
        if not consensus or not clonotypes:
            raise TCRsiftValidationError(
                f"Config entry '{label}' missing paths",
                hint="Each label needs both 'consensus' and 'clonotypes'.",
            )
        consensus_path = Path(consensus)
        clonotypes_path = Path(clonotypes)
        if not consensus_path.is_absolute():
            consensus_path = base / consensus_path
        if not clonotypes_path.is_absolute():
            clonotypes_path = base / clonotypes_path
        parsed[str(label)] = {
            "consensus": validate_file_exists(consensus_path, f"{label} consensus"),
            "clonotypes": validate_file_exists(clonotypes_path, f"{label} clonotypes"),
        }

    return dict(sorted(parsed.items(), key=lambda kv: _timepoint_sort_key(kv[0])))


def default_timepoint_inputs(data_dir: str | Path) -> dict[str, dict[str, Path]] | None:
    """Discover paired ``consensus_annotations.<TP>.csv`` + ``clonotypes.<TP>.csv``."""
    data_dir = Path(data_dir)
    if not data_dir.exists():
        return None

    discovered: dict[str, dict[str, Path]] = {}
    for consensus_path in sorted(data_dir.glob("consensus_annotations.*.csv")):
        m = re.match(r"^consensus_annotations\.(.+)\.csv$", consensus_path.name)
        if not m:
            continue
        label = m.group(1)
        clonotypes_path = data_dir / f"clonotypes.{label}.csv"
        if clonotypes_path.exists():
            discovered[label] = {
                "consensus": consensus_path,
                "clonotypes": clonotypes_path,
            }

    if not discovered:
        return None

    return dict(sorted(discovered.items(), key=lambda kv: _timepoint_sort_key(kv[0])))


def discover_gex_inputs(data_dir: str | Path, timepoint_order: list[str]) -> dict[str, dict[str, Path]]:
    """Discover per-timepoint marker scoring inputs from standard file names."""
    data_dir = Path(data_dir)
    missing: list[str] = []
    out: dict[str, dict[str, Path]] = {}
    for tp in timepoint_order:
        contig_csv = data_dir / f"filtered_contig_annotations.{tp}.csv"
        gex_h5 = data_dir / f"sample_filtered_feature_bc_matrix.{tp}.h5"
        if not contig_csv.exists():
            missing.append(str(contig_csv))
        if not gex_h5.exists():
            missing.append(str(gex_h5))
        out[tp] = {"contigs": contig_csv, "gex_h5": gex_h5}

    if missing:
        preview = "\n".join(f"- {p}" for p in missing[:10])
        raise TCRsiftValidationError(
            "Missing files required for marker scoring",
            hint=f"Expected per-timepoint contigs + 10x H5 files:\n{preview}",
        )
    return out


def _detect_count_column(clonotypes_df: pd.DataFrame, count_column: str | None = None) -> str:
    """Pick the clonotype count column used for harmonization."""
    if count_column is not None:
        if count_column not in clonotypes_df.columns:
            raise TCRsiftValidationError(
                f"Count column not found: '{count_column}'",
                hint=f"Available columns: {list(clonotypes_df.columns)}",
            )
        return count_column

    candidates = [
        "frequency",
        "cell_count",
        "count",
        "cells",
        "n_cells",
        "size",
        "total_cells",
    ]
    for cand in candidates:
        if cand in clonotypes_df.columns:
            return cand

    numeric_cols = clonotypes_df.select_dtypes(include=[np.number]).columns.tolist()
    if numeric_cols:
        return numeric_cols[0]

    raise TCRsiftValidationError(
        "Could not auto-detect clonotype count column",
        hint="Use --count-column to specify which numeric column contains clonotype counts.",
    )


def _build_consensus_pairs(consensus_df: pd.DataFrame, require_paired: bool = True) -> pd.DataFrame:
    """Build one row per clonotype with paired CDR3 alpha/beta sequences."""
    required = {"clonotype_id", "chain", "cdr3"}
    missing = required - set(consensus_df.columns)
    if missing:
        raise TCRsiftValidationError(
            f"Consensus file missing required columns: {sorted(missing)}",
            hint="Expected at least clonotype_id, chain, cdr3 columns.",
        )

    df = consensus_df.copy()
    if "productive" in df.columns:
        df = df[_truthy_mask(df["productive"])]
    if "full_length" in df.columns:
        df = df[_truthy_mask(df["full_length"])]

    df = df[df["chain"].isin(["TRA", "TRB"])]
    if len(df) == 0:
        raise TCRsiftValidationError(
            "Consensus file has no productive full-length TRA/TRB rows",
            hint="Check filtered consensus annotations and chain labels.",
        )

    sort_cols = [c for c in ["clonotype_id", "chain", "umis", "reads"] if c in df.columns]
    ascending = [True, True] + [False] * max(len(sort_cols) - 2, 0)
    if sort_cols:
        df = df.sort_values(sort_cols, ascending=ascending)

    keep_cols = {"clonotype_id", "chain", "cdr3"}
    keep_cols.update([c for c in ["v_gene", "d_gene", "j_gene", "c_gene"] if c in df.columns])
    keep_cols.update([c for c in VDJ_SEGMENT_COLS if c in df.columns])
    keep_cols.update([c for c in VDJ_SEGMENT_NT_COLS if c in df.columns])
    trimmed = df[list(keep_cols)].copy()

    rows = []
    for clonotype_id, grp in trimmed.groupby("clonotype_id", sort=False):
        tra = grp[grp["chain"] == "TRA"].head(1)
        trb = grp[grp["chain"] == "TRB"].head(1)
        if require_paired and (len(tra) == 0 or len(trb) == 0):
            continue
        if len(trb) == 0:
            continue

        tra_row = tra.iloc[0] if len(tra) > 0 else None
        trb_row = trb.iloc[0]
        cdr3_alpha = str(tra_row["cdr3"]) if tra_row is not None else ""
        cdr3_beta = str(trb_row["cdr3"])
        if not cdr3_beta:
            continue
        cdr3ab = f"{cdr3_alpha}_{cdr3_beta}"

        row = {
            "clonotype_id": clonotype_id,
            "CDR3_alpha": cdr3_alpha,
            "CDR3_beta": cdr3_beta,
            "CDR3ab": cdr3ab,
        }

        for base_col in ["v_gene", "d_gene", "j_gene", "c_gene", *VDJ_SEGMENT_COLS, *VDJ_SEGMENT_NT_COLS]:
            alpha_key = f"alpha_{base_col}"
            beta_key = f"beta_{base_col}"
            row[alpha_key] = tra_row[base_col] if (tra_row is not None and base_col in tra_row) else np.nan
            row[beta_key] = trb_row[base_col] if base_col in trb_row else np.nan

        rows.append(row)

    out = pd.DataFrame(rows)
    if len(out) == 0:
        raise TCRsiftValidationError(
            "No paired clonotypes were extracted from consensus annotations",
            hint="Ensure consensus file includes productive full-length TRA and TRB rows per clonotype.",
        )
    return out


def load_from_consensus(
    consensus_path: str | Path,
    clonotypes_path: str | Path,
    count_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load one timepoint from consensus/clonotypes CSV files.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        First DataFrame is clone-level rows with CDR3ab and cell_count.
        Second DataFrame is a tiny timepoint summary table.
    """
    consensus_path = validate_file_exists(consensus_path, "consensus annotations")
    clonotypes_path = validate_file_exists(clonotypes_path, "clonotypes")

    consensus_df = pd.read_csv(consensus_path)
    clonotypes_df = pd.read_csv(clonotypes_path)

    pairs = _build_consensus_pairs(consensus_df, require_paired=True)
    detected_count_col = _detect_count_column(clonotypes_df, count_column=count_column)

    counts = clonotypes_df[["clonotype_id", detected_count_col]].copy()
    counts = counts.rename(columns={detected_count_col: "cell_count"})
    counts["cell_count"] = pd.to_numeric(counts["cell_count"], errors="coerce").fillna(0.0)

    merged = pairs.merge(counts, on="clonotype_id", how="left")
    merged["cell_count"] = pd.to_numeric(merged["cell_count"], errors="coerce").fillna(0.0)
    merged = merged.sort_values(["cell_count", "CDR3ab"], ascending=[False, True]).reset_index(drop=True)

    stats = pd.DataFrame(
        {
            "n_clonotypes": [len(merged)],
            "total_cells": [float(merged["cell_count"].sum())],
            "count_column": [detected_count_col],
        }
    )
    return merged, stats


def _safe_zscore(x: pd.Series) -> pd.Series:
    """Stable z-score (returns zeros when variance is zero)."""
    vals = pd.to_numeric(x, errors="coerce")
    mean = float(vals.mean()) if len(vals) else 0.0
    std = float(vals.std(ddof=0)) if len(vals) else 0.0
    if std <= 0 or not np.isfinite(std):
        return pd.Series(np.zeros(len(vals), dtype=float), index=vals.index)
    return (vals - mean) / std


def _extract_gene_vector(adata, gene: str) -> np.ndarray:
    """Extract one gene vector from AnnData by exact symbol match, else zeros."""
    if gene in adata.var_names:
        vec = adata[:, gene].X
        if issparse(vec):
            return np.asarray(vec.toarray()).ravel()
        return np.asarray(vec).ravel()
    return np.zeros(adata.n_obs, dtype=float)


def _build_barcode_map(contig_barcodes: pd.Series, adata_barcodes: list[str]) -> pd.DataFrame:
    """Map contig barcodes to available AnnData barcodes with suffix-tolerant fallback."""
    adata_set = set(adata_barcodes)
    adata_by_base: dict[str, str] = {}
    for bc in adata_barcodes:
        base = _strip_barcode_suffix(bc)
        if base not in adata_by_base:
            adata_by_base[base] = bc

    rows = []
    for bc in contig_barcodes.astype(str):
        if bc in adata_set:
            rows.append({"barcode": bc, "adata_barcode": bc})
            continue
        base = _strip_barcode_suffix(bc)
        mapped = adata_by_base.get(base)
        if mapped is not None:
            rows.append({"barcode": bc, "adata_barcode": mapped})
    return pd.DataFrame(rows)


def compute_marker_scores_for_timepoint(
    consensus_path: str | Path,
    contig_csv_path: str | Path,
    gex_h5_path: str | Path,
    marker_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute per-cell and per-clonotype marker scores for one timepoint.

    Marker scores are computed from CP10K and log1p(CP10K)-z values per marker.
    """
    consensus_path = validate_file_exists(consensus_path, "consensus annotations")
    contig_csv_path = validate_file_exists(contig_csv_path, "filtered contig annotations")
    gex_h5_path = validate_file_exists(gex_h5_path, "10x filtered_feature_bc_matrix.h5")

    consensus_df = pd.read_csv(consensus_path)
    clone_map = _build_consensus_pairs(consensus_df, require_paired=False)
    clone_map = clone_map[["clonotype_id", "CDR3_alpha", "CDR3_beta", "CDR3ab"]].drop_duplicates()

    contigs = pd.read_csv(contig_csv_path)
    clonotype_col = "raw_clonotype_id" if "raw_clonotype_id" in contigs.columns else "clonotype_id"
    if clonotype_col not in contigs.columns:
        raise TCRsiftValidationError(
            f"Contig file missing clonotype id column: {contig_csv_path}",
            hint="Expected raw_clonotype_id or clonotype_id.",
        )
    if "barcode" not in contigs.columns:
        raise TCRsiftValidationError(
            f"Contig file missing barcode column: {contig_csv_path}",
            hint="Expected CellRanger filtered_contig_annotations with barcode column.",
        )

    keep = contigs[[clonotype_col, "barcode"]].copy()
    if "productive" in contigs.columns:
        keep = keep[_truthy_mask(contigs["productive"])]
    if "high_confidence" in contigs.columns:
        keep = keep[_truthy_mask(contigs["high_confidence"])]
    keep = keep.rename(columns={clonotype_col: "clonotype_id"}).drop_duplicates()

    keep = keep.merge(clone_map, on="clonotype_id", how="inner")
    if len(keep) == 0:
        raise TCRsiftValidationError(
            f"No contig barcodes mapped to paired clonotypes for {contig_csv_path}",
            hint="Check barcode/clonotype columns and consensus file consistency.",
        )

    adata = sc.read_10x_h5(str(gex_h5_path))
    barcode_map = _build_barcode_map(keep["barcode"], list(adata.obs_names))
    if len(barcode_map) == 0:
        raise TCRsiftValidationError(
            f"No barcodes from {contig_csv_path.name} were found in {gex_h5_path.name}",
            hint="Check if barcode suffixes or sample IDs differ across VDJ and GEX files.",
        )

    keep = keep.merge(barcode_map, on="barcode", how="inner")
    keep = keep.drop_duplicates(subset=["barcode", "clonotype_id", "adata_barcode"])

    adata_subset = adata[keep["adata_barcode"].tolist(), :].copy()
    if issparse(adata_subset.X):
        total_reads = np.asarray(adata_subset.X.sum(axis=1)).ravel()
    else:
        total_reads = np.asarray(adata_subset.X.sum(axis=1)).ravel()
    total_reads = np.maximum(total_reads.astype(float), 1.0)

    cell_df = keep.copy().reset_index(drop=True)
    cell_df["n_reads"] = total_reads

    marker_genes = [g.strip().upper() for g in marker_genes if g.strip()]
    marker_genes = list(dict.fromkeys(marker_genes))
    cp10k_cols = []
    z_cols = []
    for gene in marker_genes:
        counts = _extract_gene_vector(adata_subset, gene)
        cp10k = counts / total_reads * 1e4
        z = _safe_zscore(pd.Series(np.log1p(cp10k))).to_numpy()
        cell_df[f"count_{gene}"] = counts
        cell_df[f"cp10k_{gene}"] = cp10k
        cell_df[f"z_{gene}"] = z
        cp10k_cols.append(f"cp10k_{gene}")
        z_cols.append(f"z_{gene}")

    if cp10k_cols:
        cell_df["marker_mean_cp10k"] = cell_df[cp10k_cols].mean(axis=1)
    else:
        cell_df["marker_mean_cp10k"] = 0.0
    if z_cols:
        cell_df["marker_mean_z"] = cell_df[z_cols].mean(axis=1)
    else:
        cell_df["marker_mean_z"] = 0.0

    agg_spec: dict[str, tuple[str, str]] = {
        "n_cells_with_gex": ("barcode", "count"),
        "marker_score_cp10k": ("marker_mean_cp10k", "mean"),
        "marker_score_z": ("marker_mean_z", "mean"),
    }
    for gene in marker_genes:
        agg_spec[f"score_cp10k_{gene}"] = (f"cp10k_{gene}", "mean")
        agg_spec[f"score_z_{gene}"] = (f"z_{gene}", "mean")

    score_df = (
        cell_df.groupby(["CDR3_alpha", "CDR3_beta", "CDR3ab"], as_index=False)
        .agg(**agg_spec)
        .sort_values("n_cells_with_gex", ascending=False)
        .reset_index(drop=True)
    )

    return cell_df, score_df


def infer_timepoint_metric_bases(df: pd.DataFrame, timepoint_order: list[str]) -> list[str]:
    """Infer metric bases from ``<base>_<TP>`` column patterns."""
    suffixes = [f"_{tp}" for tp in timepoint_order]
    bases: set[str] = set()
    for col in df.columns:
        for suffix in suffixes:
            if col.endswith(suffix):
                bases.add(col[: -len(suffix)])
    return sorted(bases)


def add_timepoint_summaries(
    df: pd.DataFrame,
    timepoint_order: list[str],
    metric_bases: list[str],
) -> pd.DataFrame:
    """Add mean/max/min/sum summaries for all metric bases across timepoints."""
    out = df.copy()
    for base in metric_bases:
        cols = [f"{base}_{tp}" for tp in timepoint_order if f"{base}_{tp}" in out.columns]
        if not cols:
            continue
        numeric = out[cols].apply(pd.to_numeric, errors="coerce")
        out[f"{base}_mean"] = numeric.mean(axis=1)
        out[f"{base}_max"] = numeric.max(axis=1)
        out[f"{base}_min"] = numeric.min(axis=1)
        out[f"{base}_sum"] = numeric.sum(axis=1)
    return out


def add_marker_scores_to_harmonized(
    harmonized: pd.DataFrame,
    marker_scores_by_timepoint: dict[str, pd.DataFrame],
    timepoint_order: list[str],
) -> pd.DataFrame:
    """Merge per-timepoint marker score tables into harmonized clone table."""
    out = harmonized.copy()
    for tp in timepoint_order:
        score_df = marker_scores_by_timepoint.get(tp)
        if score_df is None or len(score_df) == 0:
            continue
        keys = ["CDR3ab", "CDR3_alpha", "CDR3_beta"]
        metric_cols = [c for c in score_df.columns if c not in keys]
        renamed = score_df[keys + metric_cols].copy()
        rename_map = {c: f"{c}_{tp}" for c in metric_cols}
        renamed = renamed.rename(columns=rename_map)
        out = out.merge(renamed, on=["CDR3ab", "CDR3_alpha", "CDR3_beta"], how="left")

    bases = infer_timepoint_metric_bases(out, timepoint_order)
    out = add_timepoint_summaries(out, timepoint_order, bases)
    return out


def sort_harmonized_by_rank(df: pd.DataFrame, rank_by: str = "mean_frequency") -> pd.DataFrame:
    """Sort clone table by ranking metric with sensible fallbacks."""
    out = df.copy()
    metric = rank_by if rank_by in out.columns else "mean_frequency"
    sort_cols = [metric]
    for backup in ["max_frequency", "total_cells", "n_timepoints"]:
        if backup in out.columns and backup not in sort_cols:
            sort_cols.append(backup)
    if not sort_cols:
        return out
    return out.sort_values(sort_cols, ascending=[False] * len(sort_cols)).reset_index(drop=True)


def build_harmonized_table(
    samples: dict[str, pd.DataFrame],
    timepoint_order: list[str],
    rank_by: str = "mean_frequency",
) -> pd.DataFrame:
    """Build one harmonized clone table across timepoints."""
    if not timepoint_order:
        raise TCRsiftValidationError("No timepoints loaded", hint="Check your inputs.")

    all_ids: set[str] = set()
    for tp in timepoint_order:
        if tp not in samples:
            raise TCRsiftValidationError(
                f"Missing loaded sample for timepoint '{tp}'",
                hint="Ensure every discovered timepoint was loaded successfully.",
            )
        all_ids.update(samples[tp]["CDR3ab"].dropna().astype(str).tolist())

    out = pd.DataFrame({"CDR3ab": sorted(all_ids)})
    split = out["CDR3ab"].str.split("_", n=1, expand=True)
    out["CDR3_alpha"] = split[0].fillna("")
    out["CDR3_beta"] = split[1].fillna("")

    for tp in timepoint_order:
        sample_df = samples[tp].copy()
        by_id = sample_df.set_index("CDR3ab")
        out[f"cell_count_{tp}"] = (
            out["CDR3ab"].map(by_id["cell_count"]).astype(float).fillna(0.0)
        )

    for tp in timepoint_order:
        count_col = f"cell_count_{tp}"
        denom = float(out[count_col].sum())
        if denom > 0:
            out[f"frequency_{tp}"] = out[count_col] / denom
        else:
            out[f"frequency_{tp}"] = 0.0

    count_cols = [f"cell_count_{tp}" for tp in timepoint_order]
    freq_cols = [f"frequency_{tp}" for tp in timepoint_order]
    out["total_cells"] = out[count_cols].sum(axis=1)
    out["n_timepoints"] = (out[count_cols] > 0).sum(axis=1)
    out["mean_frequency"] = out[freq_cols].mean(axis=1)
    out["max_frequency"] = out[freq_cols].max(axis=1)
    out["min_frequency"] = out[freq_cols].min(axis=1)

    # Carry over static metadata when available from any timepoint.
    metadata_cols: set[str] = set()
    for tp in timepoint_order:
        metadata_cols.update(samples[tp].columns)
    metadata_cols -= {
        "CDR3ab",
        "CDR3_alpha",
        "CDR3_beta",
        "cell_count",
    }
    for col in sorted(metadata_cols):
        out[col] = np.nan
        for tp in timepoint_order:
            sample_df = samples[tp]
            if col not in sample_df.columns:
                continue
            mapped = out["CDR3ab"].map(sample_df.set_index("CDR3ab")[col])
            out[col] = out[col].where(out[col].notna(), mapped)

    return sort_harmonized_by_rank(out, rank_by=rank_by)


def compute_increasing_masks_from_frequencies(
    freq_df: pd.DataFrame,
    ratio_nonzero_min: float = 1.5,
    ratio_all_timepoints_min: float = 1.25,
) -> tuple[pd.Series, pd.Series]:
    """Compute nonzero and all-timepoint increasing masks from frequency table."""
    if ratio_nonzero_min <= 0:
        raise ValueError(f"ratio_nonzero_min must be > 0, got {ratio_nonzero_min}")
    if ratio_all_timepoints_min <= 0:
        raise ValueError(f"ratio_all_timepoints_min must be > 0, got {ratio_all_timepoints_min}")

    arr = freq_df.apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    nonzero_mask = np.zeros(arr.shape[0], dtype=bool)
    all_positive_mask = np.zeros(arr.shape[0], dtype=bool)

    if arr.shape[1] < 2:
        return pd.Series(nonzero_mask, index=freq_df.index), pd.Series(all_positive_mask, index=freq_df.index)

    for i in range(arr.shape[0]):
        row = arr[i]
        prev = row[:-1]
        last = row[-1]
        nz_prev = prev[prev > 0]

        nonzero_mask[i] = bool(
            last > 0
            and len(nz_prev) > 0
            and last >= (ratio_nonzero_min * float(np.min(nz_prev)))
        )
        all_positive_mask[i] = bool(
            np.all(prev > 0)
            and last > 0
            and last >= (ratio_all_timepoints_min * float(np.max(prev)))
        )

    return pd.Series(nonzero_mask, index=freq_df.index), pd.Series(all_positive_mask, index=freq_df.index)


def _resolve_best_metric_series(
    df: pd.DataFrame,
    candidates: list[str],
) -> tuple[str | None, pd.Series]:
    """Resolve first available metric column and return a numeric series."""
    for col in candidates:
        if col in df.columns:
            return col, pd.to_numeric(df[col], errors="coerce")
    return None, pd.Series(np.nan, index=df.index, dtype=float)


def _resolve_gene_cp10k_series(
    df: pd.DataFrame,
    gene: str,
    timepoint_order: list[str],
) -> tuple[str | None, pd.Series]:
    """Resolve best CP10K series for a gene using summary-then-timepoint fallback."""
    candidates = [
        f"score_cp10k_{gene}_mean",
        f"score_cp10k_{gene}_max",
        f"score_cp10k_{gene}_sum",
        *[f"score_cp10k_{gene}_{tp}" for tp in timepoint_order],
    ]
    return _resolve_best_metric_series(df, candidates)


def _select_top_panel_score(
    df: pd.DataFrame,
    score_col: str,
    flag_col: str,
    immunogenic_percentile: float,
    immunogenic_percentile_slack_frac: float,
    immunogenic_min_cp10k: float,
    immunogenic_require_above_median: bool,
) -> int:
    """Apply percentile-based top-score selection for one aggregate panel score."""
    aux_prefix = flag_col.replace("is_", "score_sel_", 1)
    df[flag_col] = False
    if score_col not in df.columns:
        df[f"{aux_prefix}_median"] = np.nan
        df[f"{aux_prefix}_percentile"] = np.nan
        df[f"{aux_prefix}_cutoff"] = np.nan
        return 0

    score_values = pd.to_numeric(df[score_col], errors="coerce")
    base_mask = df["is_base_selected"] & score_values.notna()
    if int(base_mask.sum()) == 0:
        df[f"{aux_prefix}_median"] = np.nan
        df[f"{aux_prefix}_percentile"] = np.nan
        df[f"{aux_prefix}_cutoff"] = np.nan
        return 0

    median_value = float(score_values[base_mask].median())
    df[f"{aux_prefix}_median"] = median_value
    floor = max(float(immunogenic_min_cp10k), 0.0)
    eligible = df["is_base_selected"] & score_values.gt(floor)
    if immunogenic_require_above_median:
        eligible &= score_values.gt(median_value)
    if int(eligible.sum()) == 0:
        df[f"{aux_prefix}_percentile"] = np.nan
        df[f"{aux_prefix}_cutoff"] = np.nan
        return 0

    eligible_values = score_values[eligible]
    percentile_value = float(np.nanquantile(eligible_values.to_numpy(dtype=float), immunogenic_percentile))
    cutoff_value = percentile_value * (1.0 - immunogenic_percentile_slack_frac)
    df[f"{aux_prefix}_percentile"] = percentile_value
    df[f"{aux_prefix}_cutoff"] = cutoff_value
    selected = eligible & score_values.ge(cutoff_value)
    df.loc[selected, flag_col] = True
    return int((df["is_base_selected"] & df[flag_col]).sum())


def run_selection_pipeline(
    master_df: pd.DataFrame,
    timepoint_order: list[str],
    marker_genes: list[str],
    min_cells_per_clone: int = 2,
    min_cd8_cp10k: float = 0.0,
    max_cd4_to_cd8_ratio: float = 1.0,
    increase_ratio_nonzero_min: float = 1.5,
    increase_ratio_all_timepoints_min: float = 1.25,
    immunogenic_percentile: float = 0.90,
    immunogenic_percentile_slack_frac: float = 0.01,
    immunogenic_min_cp10k: float = 0.0,
    immunogenic_require_above_median: bool = True,
    cytotoxic_last_min_z: float = 0.25,
    cytotoxic_last_min_cp10k: float = 0.05,
    became_cytotoxic_min_delta_z: float = 0.5,
    immunogenic_genes_preferred: list[str] | None = None,
    cytotoxic_genes_preferred: list[str] | None = None,
    cytolytic_genes_preferred: list[str] | None = None,
    antigen_response_genes_preferred: list[str] | None = None,
) -> tuple[
    pd.DataFrame,
    dict[str, pd.DataFrame],
    list[tuple[str, int]],
    list[tuple[str, int]],
    list[tuple[str, int]],
]:
    """Compute selection masks/subsets for promising TIL clonotypes."""
    df = master_df.copy()
    if not timepoint_order:
        raise ValueError("No timepoints provided for selection pipeline")
    if min_cells_per_clone < 1:
        raise ValueError(f"min_cells_per_clone must be >= 1, got {min_cells_per_clone}")
    if increase_ratio_nonzero_min <= 0:
        raise ValueError(f"increase_ratio_nonzero_min must be > 0, got {increase_ratio_nonzero_min}")
    if increase_ratio_all_timepoints_min <= 0:
        raise ValueError(
            "increase_ratio_all_timepoints_min must be > 0, "
            f"got {increase_ratio_all_timepoints_min}"
        )
    if not (0 < immunogenic_percentile <= 1):
        raise ValueError(f"immunogenic_percentile must be in (0,1], got {immunogenic_percentile}")
    if not (0 <= immunogenic_percentile_slack_frac < 1):
        raise ValueError(
            "immunogenic_percentile_slack_frac must be in [0,1), "
            f"got {immunogenic_percentile_slack_frac}"
        )

    freq_cols = [f"frequency_{tp}" for tp in timepoint_order if f"frequency_{tp}" in df.columns]
    if len(freq_cols) < 2:
        raise ValueError("Need at least two frequency_<TP> columns for selection.")
    nonzero_increasing, all_positive_increasing = compute_increasing_masks_from_frequencies(
        df[freq_cols],
        ratio_nonzero_min=increase_ratio_nonzero_min,
        ratio_all_timepoints_min=increase_ratio_all_timepoints_min,
    )
    df["is_increasing"] = nonzero_increasing
    df["is_increasing_nonzero"] = nonzero_increasing
    df["is_increasing_positive"] = all_positive_increasing
    df["is_increasing_all_timepoints"] = all_positive_increasing
    df["is_enriched"] = df["is_increasing_nonzero"]

    if "total_cells" not in df.columns:
        count_cols = [f"cell_count_{tp}" for tp in timepoint_order if f"cell_count_{tp}" in df.columns]
        if not count_cols:
            raise ValueError("Missing total_cells and per-timepoint cell_count columns.")
        df["total_cells"] = df[count_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=1)
    else:
        df["total_cells"] = pd.to_numeric(df["total_cells"], errors="coerce").fillna(0.0)
    df["is_min_cells"] = df["total_cells"] >= float(min_cells_per_clone)

    first_tp = timepoint_order[0]
    last_tp = timepoint_order[-1]
    first_freq = pd.to_numeric(df.get(f"frequency_{first_tp}", 0.0), errors="coerce").fillna(0.0)
    last_freq = pd.to_numeric(df.get(f"frequency_{last_tp}", 0.0), errors="coerce").fillna(0.0)
    df["delta_frequency_last_vs_first"] = last_freq - first_freq
    df["fold_change_last_vs_first"] = np.where(
        first_freq > 0,
        last_freq / first_freq,
        np.where(last_freq > 0, np.inf, 1.0),
    )

    if "is_viral" in df.columns:
        df["is_non_viral"] = ~df["is_viral"].fillna(False).astype(bool)
    else:
        df["is_non_viral"] = True

    cd8_candidate_cols = [
        c for c in df.columns if c.startswith("score_cp10k_CD8A") or c.startswith("score_cp10k_CD8B")
    ]
    if cd8_candidate_cols:
        df["cd8_cp10k_max"] = (
            df[cd8_candidate_cols].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)
        )
    else:
        df["cd8_cp10k_max"] = 0.0

    cd4_candidate_cols = [c for c in df.columns if c.startswith("score_cp10k_CD4")]
    if cd4_candidate_cols:
        df["cd4_cp10k_max"] = (
            df[cd4_candidate_cols].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)
        )
    else:
        df["cd4_cp10k_max"] = 0.0

    df["is_cd8_positive"] = df["cd8_cp10k_max"] > float(min_cd8_cp10k)
    df["cd4_to_cd8_ratio"] = np.where(
        df["cd8_cp10k_max"] > 0,
        df["cd4_cp10k_max"] / df["cd8_cp10k_max"],
        np.where(df["cd4_cp10k_max"] > 0, np.inf, 0.0),
    )
    df["is_cd8_not_cd4"] = df["cd4_to_cd8_ratio"] < float(max_cd4_to_cd8_ratio)
    df["is_base_cd8_non_viral"] = df["is_non_viral"] & df["is_cd8_positive"] & df["is_cd8_not_cd4"]
    df["is_base_selected"] = df["is_base_cd8_non_viral"] & df["is_min_cells"]

    requested_genes = [g.strip().upper() for g in marker_genes if g.strip()]
    requested_genes = list(dict.fromkeys(requested_genes))

    if immunogenic_genes_preferred:
        immunogenic_default = [g.strip().upper() for g in immunogenic_genes_preferred if g.strip()]
    else:
        immunogenic_default = ["GZMB", "PRF1", "IFNG", "MKI67", "TNFRSF9"]
    immunogenic_genes = [g for g in immunogenic_default if g in requested_genes]
    if not immunogenic_genes:
        immunogenic_genes = [g for g in requested_genes if g not in {"CD4", "CD8A", "CD8B"}]
    if not immunogenic_genes:
        immunogenic_genes = immunogenic_default

    def _dedupe_valid(candidate: list[str], fallback: tuple[str, ...]) -> list[str]:
        deduped = list(dict.fromkeys([g.strip().upper() for g in candidate if g and g.strip()]))
        valid = [g for g in deduped if g in requested_genes]
        if valid:
            return valid
        return [g for g in fallback if g in requested_genes]

    cytotoxic_seed = (
        [g.strip().upper() for g in cytotoxic_genes_preferred if g.strip()]
        if cytotoxic_genes_preferred
        else ["GZMB", "PRF1", "IFNG", "MKI67", "TNFRSF9"]
    )
    cytotoxic_genes = _dedupe_valid(cytotoxic_seed, tuple(cytotoxic_seed))
    if not cytotoxic_genes:
        cytotoxic_genes = immunogenic_genes

    cytolytic_seed = (
        [g.strip().upper() for g in cytolytic_genes_preferred if g.strip()]
        if cytolytic_genes_preferred
        else list(CYTOLYTIC_GENES_DEFAULT)
    )
    cytolytic_genes = _dedupe_valid(cytolytic_seed, CYTOLYTIC_GENES_DEFAULT)

    antigen_seed = (
        [g.strip().upper() for g in antigen_response_genes_preferred if g.strip()]
        if antigen_response_genes_preferred
        else list(ANTIGEN_RESPONSE_GENES_DEFAULT)
    )
    antigen_response_genes = _dedupe_valid(antigen_seed, ANTIGEN_RESPONSE_GENES_DEFAULT)

    for tp in timepoint_order:
        z_cols = [f"score_z_{g}_{tp}" for g in cytotoxic_genes if f"score_z_{g}_{tp}" in df.columns]
        cp10k_cols = [
            f"score_cp10k_{g}_{tp}" for g in cytotoxic_genes if f"score_cp10k_{g}_{tp}" in df.columns
        ]
        if z_cols:
            df[f"cytotoxic_score_z_{tp}"] = (
                df[z_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).fillna(0.0)
            )
        else:
            df[f"cytotoxic_score_z_{tp}"] = 0.0
        if cp10k_cols:
            df[f"cytotoxic_score_cp10k_{tp}"] = (
                df[cp10k_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).fillna(0.0)
            )
        else:
            df[f"cytotoxic_score_cp10k_{tp}"] = 0.0

    def _add_panel_scores(panel_name: str, genes: list[str]) -> None:
        per_tp_cols = []
        for tp in timepoint_order:
            src_cols = [f"score_cp10k_{g}_{tp}" for g in genes if f"score_cp10k_{g}_{tp}" in df.columns]
            out_col = f"{panel_name}_score_cp10k_{tp}"
            if src_cols:
                df[out_col] = (
                    df[src_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).fillna(0.0)
                )
            else:
                # Fallback to resolved summary columns when per-timepoint cols are absent.
                series_list = []
                for g in genes:
                    _source, series = _resolve_gene_cp10k_series(df, g, timepoint_order)
                    if series.notna().any():
                        series_list.append(series.fillna(0.0))
                if series_list:
                    df[out_col] = pd.concat(series_list, axis=1).mean(axis=1)
                else:
                    df[out_col] = 0.0
            per_tp_cols.append(out_col)
        df[f"{panel_name}_score_cp10k_mean"] = df[per_tp_cols].mean(axis=1)
        df[f"{panel_name}_score_cp10k_max"] = df[per_tp_cols].max(axis=1)
        df[f"{panel_name}_score_cp10k_sum"] = df[per_tp_cols].sum(axis=1)

    _add_panel_scores("cytolytic", cytolytic_genes)
    _add_panel_scores("antigen_response", antigen_response_genes)

    first_z = pd.to_numeric(df.get(f"cytotoxic_score_z_{first_tp}", 0.0), errors="coerce").fillna(0.0)
    last_z = pd.to_numeric(df.get(f"cytotoxic_score_z_{last_tp}", 0.0), errors="coerce").fillna(0.0)
    first_cp10k = pd.to_numeric(
        df.get(f"cytotoxic_score_cp10k_{first_tp}", 0.0),
        errors="coerce",
    ).fillna(0.0)
    last_cp10k = pd.to_numeric(
        df.get(f"cytotoxic_score_cp10k_{last_tp}", 0.0),
        errors="coerce",
    ).fillna(0.0)

    df["cytotoxic_score_z_delta_last_vs_first"] = last_z - first_z
    df["cytotoxic_score_cp10k_delta_last_vs_first"] = last_cp10k - first_cp10k
    df["is_cytotoxic_at_last"] = (last_z >= cytotoxic_last_min_z) | (
        last_cp10k >= cytotoxic_last_min_cp10k
    )
    df["is_became_cytotoxic"] = (
        df["cytotoxic_score_z_delta_last_vs_first"] >= became_cytotoxic_min_delta_z
    )
    df["is_increasing_and_became_cytotoxic"] = (
        df["is_increasing_nonzero"] & df["is_became_cytotoxic"] & df["is_cytotoxic_at_last"]
    )

    top_flag_cols = []
    immunogenic_branch_counts_within_base: list[tuple[str, int]] = []
    for gene in immunogenic_genes:
        ranking_source, ranking_values = _resolve_gene_cp10k_series(df, gene, timepoint_order)
        flag_col = f"is_top_immunogenic_{gene}"
        df[flag_col] = False
        df[f"immunogenic_cp10k_source_{gene}"] = ranking_source if ranking_source else ""

        base_mask = df["is_base_selected"] & ranking_values.notna()
        if int(base_mask.sum()) == 0:
            df[f"immunogenic_cp10k_median_{gene}"] = np.nan
            immunogenic_branch_counts_within_base.append((gene, 0))
            continue

        median_value = float(ranking_values[base_mask].median())
        df[f"immunogenic_cp10k_median_{gene}"] = median_value
        floor = max(float(immunogenic_min_cp10k), 0.0)
        eligible = df["is_base_selected"] & ranking_values.gt(floor)
        if immunogenic_require_above_median:
            eligible &= ranking_values.gt(median_value)
        if int(eligible.sum()) == 0:
            df[f"immunogenic_cp10k_percentile_{gene}"] = np.nan
            df[f"immunogenic_cp10k_cutoff_{gene}"] = np.nan
            immunogenic_branch_counts_within_base.append((gene, 0))
            continue

        values = ranking_values[eligible]
        percentile_value = float(np.nanquantile(values.to_numpy(dtype=float), immunogenic_percentile))
        cutoff_value = percentile_value * (1.0 - immunogenic_percentile_slack_frac)
        df[f"immunogenic_cp10k_percentile_{gene}"] = percentile_value
        df[f"immunogenic_cp10k_cutoff_{gene}"] = cutoff_value
        selected = eligible & ranking_values.ge(cutoff_value)
        df.loc[selected, flag_col] = True
        top_flag_cols.append(flag_col)
        n_selected = int((df["is_base_selected"] & df[flag_col]).sum())
        immunogenic_branch_counts_within_base.append((gene, n_selected))

    if top_flag_cols:
        df["is_top_immunogenic_any"] = df[top_flag_cols].any(axis=1)
    else:
        df["is_top_immunogenic_any"] = False

    n_top_cytolytic = _select_top_panel_score(
        df,
        "cytolytic_score_cp10k_mean",
        "is_top_cytolytic_score",
        immunogenic_percentile,
        immunogenic_percentile_slack_frac,
        immunogenic_min_cp10k,
        immunogenic_require_above_median,
    )
    n_top_antigen_response = _select_top_panel_score(
        df,
        "antigen_response_score_cp10k_mean",
        "is_top_antigen_response_score",
        immunogenic_percentile,
        immunogenic_percentile_slack_frac,
        immunogenic_min_cp10k,
        immunogenic_require_above_median,
    )

    df["is_branch_top_immunogenic"] = df["is_top_immunogenic_any"]
    df["is_branch_cytolytic"] = df["is_top_cytolytic_score"]
    df["is_branch_antigen_response"] = df["is_top_antigen_response_score"]
    df["is_branch_enriched"] = df["is_increasing_nonzero"]
    df["is_branch_increasing"] = df["is_increasing_all_timepoints"]
    df["is_branch_any"] = (
        df["is_branch_top_immunogenic"]
        | df["is_branch_cytolytic"]
        | df["is_branch_antigen_response"]
        | df["is_branch_enriched"]
        | df["is_branch_increasing"]
    )
    df["is_immunogenic"] = (
        df["is_top_immunogenic_any"]
        | df["is_top_cytolytic_score"]
        | df["is_top_antigen_response_score"]
    )
    df["is_candidate_tumor_reactive"] = df["is_base_selected"] & df["is_branch_any"]
    df["is_branch_union_within_base"] = df["is_candidate_tumor_reactive"]
    df["is_priority_cytotoxic_expander"] = (
        df["is_base_selected"] & df["is_increasing_and_became_cytotoxic"]
    )

    subsets = {
        "subset_non_viral": df[df["is_non_viral"]].copy(),
        "subset_base_cd8_nonviral": df[df["is_base_cd8_non_viral"]].copy(),
        "subset_min_cells": df[df["is_min_cells"]].copy(),
        "subset_base_selected": df[df["is_base_selected"]].copy(),
        "subset_immunogenic": df[df["is_immunogenic"]].copy(),
        "subset_top_immunogenic_any": df[df["is_top_immunogenic_any"]].copy(),
        "subset_top_cytolytic_score": df[df["is_top_cytolytic_score"]].copy(),
        "subset_top_antigen_response_score": df[df["is_top_antigen_response_score"]].copy(),
        "subset_enriched": df[df["is_enriched"]].copy(),
        "subset_increasing": df[df["is_increasing_nonzero"]].copy(),
        "subset_increasing_positive": df[df["is_increasing_all_timepoints"]].copy(),
        "subset_increasing_nonzero": df[df["is_increasing_nonzero"]].copy(),
        "subset_increasing_all_timepoints": df[df["is_increasing_all_timepoints"]].copy(),
        "subset_increasing_became_cytotoxic": df[df["is_increasing_and_became_cytotoxic"]].copy(),
        "subset_priority_cytotoxic_expander": df[df["is_priority_cytotoxic_expander"]].copy(),
        "subset_candidate_tumor_reactive": df[df["is_candidate_tumor_reactive"]].copy(),
    }
    for gene in immunogenic_genes:
        flag_col = f"is_top_immunogenic_{gene}"
        subsets[f"subset_top_immunogenic_{gene}"] = df[df[flag_col]].copy()

    terminal_freq_col = f"frequency_{last_tp}"
    sort_cols = [c for c in [terminal_freq_col, "mean_frequency", "total_cells"] if c in df.columns]
    if sort_cols:
        for key, subset_df in list(subsets.items()):
            if len(subset_df) == 0:
                continue
            subsets[key] = subset_df.sort_values(sort_cols, ascending=[False] * len(sort_cols))

    immunogenic_branch_rows = [(f"Imm.{gene}", count) for gene, count in immunogenic_branch_counts_within_base]
    immunogenic_independent_rows = []
    for gene in immunogenic_genes:
        flag_col = f"is_top_immunogenic_{gene}"
        immunogenic_independent_rows.append((f"Imm.{gene}", int(df[flag_col].sum())))

    cd4_stage_label = (
        "CD4<CD8" if np.isclose(max_cd4_to_cd8_ratio, 1.0) else f"CD4/CD8<{max_cd4_to_cd8_ratio:g}"
    )
    base_non_viral = df["is_non_viral"]
    base_cd8 = base_non_viral & df["is_cd8_positive"]
    base_cd4_ratio = base_cd8 & df["is_cd8_not_cd4"]
    base_min_cells = base_cd4_ratio & df["is_min_cells"]

    sequential_stages: list[tuple[str, int]] = []
    branch_union_stages = [
        ("All", len(df)),
        ("NonViral", int(base_non_viral.sum())),
        ("CD8", int(base_cd8.sum())),
        (cd4_stage_label, int(base_cd4_ratio.sum())),
        (f"Cells>={min_cells_per_clone}", int(base_min_cells.sum())),
        *immunogenic_branch_rows,
        ("Cytolytic", n_top_cytolytic),
        ("AgResp", n_top_antigen_response),
        (f"Inc>={increase_ratio_nonzero_min:g}x", int((df["is_base_selected"] & df["is_branch_enriched"]).sum())),
        ("Union", int(df["is_branch_union_within_base"].sum())),
        ("Final", int(df["is_candidate_tumor_reactive"].sum())),
    ]

    independent_stages = [
        ("All", len(df)),
        ("NonViral", int(df["is_non_viral"].sum())),
        ("CD8", int(df["is_cd8_positive"].sum())),
        (cd4_stage_label, int(df["is_cd8_not_cd4"].sum())),
        (f"Cells>={min_cells_per_clone}", int(df["is_min_cells"].sum())),
        ("Base", int(df["is_base_selected"].sum())),
        ("ImmAny", int(df["is_top_immunogenic_any"].sum())),
        ("Cytolytic", int(df["is_top_cytolytic_score"].sum())),
        ("AgResp", int(df["is_top_antigen_response_score"].sum())),
        *immunogenic_independent_rows,
        (f"Inc>={increase_ratio_nonzero_min:g}x", int(df["is_increasing_nonzero"].sum())),
        ("Inc+Cyto", int(df["is_increasing_and_became_cytotoxic"].sum())),
        ("Priority", int(df["is_priority_cytotoxic_expander"].sum())),
        ("Final", int(df["is_candidate_tumor_reactive"].sum())),
    ]

    return df, subsets, sequential_stages, independent_stages, branch_union_stages


def plot_selection_funnel(stages: list[tuple[str, int]], output_path: Path, title: str) -> None:
    """Plot a simple selection funnel bar chart."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels = [s[0] for s in stages]
    values = [s[1] for s in stages]
    plt.figure(figsize=(10, max(4, 0.35 * len(labels))))
    ax = sns.barplot(x=values, y=labels, color="#4c72b0")
    ax.set_title(title)
    ax.set_xlabel("n clonotypes")
    ax.set_ylabel("")
    for i, v in enumerate(values):
        ax.text(v + max(values) * 0.01 if max(values) > 0 else 0.1, i, f"{v:,}", va="center")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_heatmap(df: pd.DataFrame, count_cols: list[str], output_path: Path, top_k: int, title: str) -> None:
    """Plot a compact top-k heatmap for per-timepoint counts."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top = df.head(top_k).copy()
    if len(top) == 0 or not count_cols:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "No data for heatmap", ha="center", va="center")
        plt.axis("off")
        plt.savefig(output_path, dpi=200)
        plt.close()
        return
    mat = top.set_index("CDR3ab")[count_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    plt.figure(figsize=(max(6, 1.2 * len(count_cols)), max(4, 0.25 * len(top) + 2)))
    sns.heatmap(mat, cmap="YlOrRd")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def create_selected_clone_pdf_report(
    selected_df: pd.DataFrame,
    output_path: Path,
    timepoint_order: list[str],
    marker_genes: list[str],
    pyensembl_release: int = 110,
) -> None:
    """
    Write a compact PDF report for selected clones.

    The report intentionally remains lightweight and dependency-free
    (matplotlib backend only) while keeping the expected output artifact name.
    """
    del pyensembl_release  # kept for CLI compatibility
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output_path) as pdf:
        page = 0
        if len(selected_df) == 0:
            fig = plt.figure(figsize=(8.5, 11))
            plt.axis("off")
            plt.title("Selected Clones Report")
            plt.text(
                0.5,
                0.5,
                "No clones met current selection criteria.",
                ha="center",
                va="center",
            )
            pdf.savefig(fig)
            plt.close(fig)
            return

        cols = [
            c
            for c in [
                "CDR3ab",
                "total_cells",
                "mean_frequency",
                f"frequency_{timepoint_order[-1]}" if timepoint_order else None,
                "is_candidate_tumor_reactive",
                "is_priority_cytotoxic_expander",
            ]
            if c and c in selected_df.columns
        ]
        for gene in marker_genes:
            col = f"score_cp10k_{gene}_mean"
            if col in selected_df.columns:
                cols.append(col)

        show = selected_df[cols].head(200).copy()
        chunk_size = 35
        for start in range(0, len(show), chunk_size):
            chunk = show.iloc[start : start + chunk_size]
            fig, ax = plt.subplots(figsize=(11, 8.5))
            ax.axis("off")
            ax.set_title(f"Selected Clones ({start + 1}-{start + len(chunk)} of {len(show)})")
            table = ax.table(
                cellText=chunk.astype(str).values,
                colLabels=chunk.columns.tolist(),
                cellLoc="left",
                loc="center",
            )
            table.auto_set_font_size(False)
            table.set_fontsize(6)
            table.scale(1, 1.2)
            pdf.savefig(fig)
            plt.close(fig)
            page += 1
        if page == 0:
            fig = plt.figure(figsize=(8.5, 11))
            plt.axis("off")
            plt.text(0.5, 0.5, "No report rows.", ha="center", va="center")
            pdf.savefig(fig)
            plt.close(fig)


def _write_subset_tables(subset_dfs: dict[str, pd.DataFrame], fig_dir: Path) -> None:
    """Write subset tables to ``subset_*.csv`` files."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    for name, subset_df in subset_dfs.items():
        subset_df.to_csv(fig_dir / f"{name}.csv", index=False)


def run_til_select(args: argparse.Namespace) -> pd.DataFrame:
    """
    Execute TIL-only clone selection workflow.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments from ``til-select`` command.

    Returns
    -------
    pd.DataFrame
        Final selected master table with masks.
    """
    if getattr(args, "config", None):
        discovered = parse_config(args.config)
    elif getattr(args, "samples", None):
        discovered = parse_sample_args(args.samples)
    else:
        discovered = default_timepoint_inputs(args.data_dir)
        if not discovered:
            raise TCRsiftValidationError(
                f"No timepoint inputs discovered under {args.data_dir}",
                hint="Provide --samples LABEL=CONSENSUS,CLONOTYPES or --config YAML.",
            )

    timepoint_order = list(discovered.keys())
    samples: dict[str, pd.DataFrame] = {}
    for tp in timepoint_order:
        counts_df, _stats_df = load_from_consensus(
            discovered[tp]["consensus"],
            discovered[tp]["clonotypes"],
            count_column=getattr(args, "count_column", None),
        )
        samples[tp] = counts_df

    marker_genes = [g.strip().upper() for g in str(args.marker_genes).split(",") if g.strip()]
    marker_genes = list(dict.fromkeys(marker_genes))
    if not marker_genes:
        marker_genes = MARKER_GENES_DEFAULT.copy()

    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    gex_inputs = discover_gex_inputs(args.data_dir, timepoint_order)
    marker_scores_by_tp: dict[str, pd.DataFrame] = {}
    for tp in timepoint_order:
        cell_df, score_df = compute_marker_scores_for_timepoint(
            consensus_path=discovered[tp]["consensus"],
            contig_csv_path=gex_inputs[tp]["contigs"],
            gex_h5_path=gex_inputs[tp]["gex_h5"],
            marker_genes=marker_genes,
        )
        cell_df.to_csv(fig_dir / f"marker_cells_{tp}.csv", index=False)
        score_df.to_csv(fig_dir / f"marker_clonotype_scores_{tp}.csv", index=False)
        marker_scores_by_tp[tp] = score_df

    harmonized = build_harmonized_table(samples, timepoint_order, rank_by="mean_frequency")
    harmonized = add_marker_scores_to_harmonized(harmonized, marker_scores_by_tp, timepoint_order)

    annotated = annotate_clonotypes(
        harmonized,
        vdjdb_path=getattr(args, "vdjdb", None),
        iedb_path=getattr(args, "iedb", None),
        cedar_path=getattr(args, "cedar", None),
        match_by=getattr(args, "match_by", "CDR3b_only"),
        exclude_viral=False,
        flag_only=True,
    )
    annotated["db_vdjdb"] = annotated["db_database"].fillna("").str.contains("VDJdb")
    annotated["db_iedb"] = annotated["db_database"].fillna("").str.contains("IEDB")
    annotated["db_cedar"] = annotated["db_database"].fillna("").str.contains("CEDAR")

    metric_bases = infer_timepoint_metric_bases(annotated, timepoint_order)
    master_df = add_timepoint_summaries(annotated, timepoint_order, metric_bases)

    immunogenic_genes = [g.strip().upper() for g in str(args.immunogenic_genes).split(",") if g.strip()]
    cytotoxic_genes = [g.strip().upper() for g in str(args.cytotoxic_genes).split(",") if g.strip()]
    cytolytic_genes = [g.strip().upper() for g in str(args.cytolytic_genes).split(",") if g.strip()]
    antigen_response_genes = [
        g.strip().upper() for g in str(args.antigen_response_genes).split(",") if g.strip()
    ]

    master_df, subset_dfs, _seq, _ind, branch_stages = run_selection_pipeline(
        master_df,
        timepoint_order=timepoint_order,
        marker_genes=marker_genes,
        min_cells_per_clone=int(args.min_cells_per_clone),
        min_cd8_cp10k=float(args.min_cd8_cp10k),
        max_cd4_to_cd8_ratio=float(args.max_cd4_to_cd8_ratio),
        increase_ratio_nonzero_min=float(args.increase_ratio_nonzero_min),
        increase_ratio_all_timepoints_min=float(args.increase_ratio_all_timepoints_min),
        immunogenic_percentile=float(args.immunogenic_percentile),
        immunogenic_percentile_slack_frac=float(args.immunogenic_percentile_slack_frac),
        immunogenic_min_cp10k=float(args.immunogenic_min_cp10k),
        immunogenic_require_above_median=bool(args.immunogenic_require_above_median),
        cytotoxic_last_min_z=float(args.cytotoxic_last_min_z),
        cytotoxic_last_min_cp10k=float(args.cytotoxic_last_min_cp10k),
        became_cytotoxic_min_delta_z=float(args.became_cytotoxic_min_delta_z),
        immunogenic_genes_preferred=immunogenic_genes,
        cytotoxic_genes_preferred=cytotoxic_genes,
        cytolytic_genes_preferred=cytolytic_genes,
        antigen_response_genes_preferred=antigen_response_genes,
    )

    _write_subset_tables(subset_dfs, fig_dir)

    out_table = Path(args.out_table)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    master_df.to_csv(out_table, index=False)
    master_df.to_csv(fig_dir / "abTCR_master_table.csv", index=False)

    out_annotated = Path(args.out_annotated) if args.out_annotated else fig_dir / "abTCR_annotated.csv"
    out_annotated.parent.mkdir(parents=True, exist_ok=True)
    master_df.to_csv(out_annotated, index=False)

    mask_cols = sorted([c for c in master_df.columns if c.startswith("is_")])
    id_cols = [c for c in ["CDR3ab", "CDR3_alpha", "CDR3_beta"] if c in master_df.columns]
    metric_cols = [
        c
        for c in [
            *[f"cell_count_{tp}" for tp in timepoint_order],
            *[f"frequency_{tp}" for tp in timepoint_order],
            "cd8_cp10k_max",
            "cd4_cp10k_max",
            "cd4_to_cd8_ratio",
            "delta_frequency_last_vs_first",
            "fold_change_last_vs_first",
            "cytotoxic_score_z_delta_last_vs_first",
            "cytotoxic_score_cp10k_delta_last_vs_first",
        ]
        if c in master_df.columns
    ]
    master_df[id_cols + metric_cols + mask_cols].to_csv(fig_dir / "selection_masks.csv", index=False)

    out_selected_report = (
        Path(args.out_selected_report)
        if args.out_selected_report
        else fig_dir / "selected_clones_report.pdf"
    )
    selected = master_df[master_df["is_candidate_tumor_reactive"]].copy()
    create_selected_clone_pdf_report(
        selected,
        out_selected_report,
        timepoint_order=timepoint_order,
        marker_genes=marker_genes,
        pyensembl_release=int(args.pyensembl_release),
    )

    rank_source = sort_harmonized_by_rank(master_df, rank_by=args.rank_by)
    out_top = Path(args.out_top) if args.out_top else fig_dir / "abTCR_topk.csv"
    out_heatmap = Path(args.out_heatmap) if args.out_heatmap else fig_dir / "abTCR_topk_heatmap.png"
    out_top.parent.mkdir(parents=True, exist_ok=True)
    out_heatmap.parent.mkdir(parents=True, exist_ok=True)
    rank_source.head(int(args.top_k)).to_csv(out_top, index=False)
    count_cols = [f"cell_count_{tp}" for tp in timepoint_order if f"cell_count_{tp}" in rank_source.columns]
    plot_heatmap(
        rank_source,
        count_cols=count_cols,
        output_path=out_heatmap,
        top_k=int(args.top_k),
        title=f"Top {int(args.top_k)} abTCRs (ranked by {args.rank_by})",
    )

    for rank_metric in RANK_METRICS:
        ranked = sort_harmonized_by_rank(master_df, rank_by=rank_metric)
        ranked.head(int(args.top_k)).to_csv(fig_dir / f"abTCR_topk_{rank_metric}.csv", index=False)
        plot_heatmap(
            ranked,
            count_cols=count_cols,
            output_path=fig_dir / f"abTCR_topk_heatmap_{rank_metric}.png",
            top_k=int(args.top_k),
            title=f"Top {int(args.top_k)} abTCRs (ranked by {rank_metric})",
        )

    plot_selection_funnel(
        branch_stages,
        output_path=fig_dir / "selection_funnel.png",
        title="Selection Funnel",
    )

    logger.info(
        "til-select complete: %s clones total, %s final candidates",
        len(master_df),
        int(master_df["is_candidate_tumor_reactive"].sum()),
    )
    return master_df
