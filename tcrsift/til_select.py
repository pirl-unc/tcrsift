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

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from scipy import sparse
from scipy.sparse import issparse

from .signatures import (
    ANTIGEN_RESPONSE_GENES_HGNC as ANTIGEN_RESPONSE_GENES_DEFAULT,
)
from .signatures import (
    CYTOLYTIC_GENES_HGNC as CYTOLYTIC_GENES_DEFAULT,
)
from .signatures import (
    TUMOR_REACTIVE_GENES_HGNC as ENRICHMENT_GENES_DEFAULT,
)
from .validation import TCRsiftValidationError, validate_file_exists

logger = logging.getLogger(__name__)

VIRAL_SPECIES_PATTERNS = [
    "cmv",
    "cytomegalovirus",
    "epstein-barr",
    "epstein barr",
    "ebv",
    "hiv",
    "human immunodeficiency",
    "flu",
    "influenza",
    "sars",
    "coronavirus",
    "herpes",
    "hsv",
    "hpv",
    "papilloma",
    "hepatitis",
    "hbv",
    "hcv",
    "dengue",
    "zika",
    "yellow fever",
]

MARKER_GENES_DEFAULT = [
    "CD4",
    "CD8A",
    "CD8B",
    "GZMB",
    "PRF1",
    "IFNG",
    "MKI67",
    "TNFRSF9",
    "CXCL13",
    "ENTPD1",
]
TCR_SEGMENTS_AA = ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4")
TCR_SEGMENTS_NT = tuple(f"{segment}_nt" for segment in TCR_SEGMENTS_AA)
VDJ_SEGMENT_COLS = list(TCR_SEGMENTS_AA)
VDJ_SEGMENT_NT_COLS = list(TCR_SEGMENTS_NT)
# ``CYTOLYTIC_GENES_DEFAULT`` / ``ANTIGEN_RESPONSE_GENES_DEFAULT`` /
# ``ENRICHMENT_GENES_DEFAULT`` are now back-compat aliases imported
# above from :mod:`tcrsift.signatures`, which can drive non-TIL
# selections too.
RANK_METRICS = (
    "mean_frequency",
    "max_frequency",
    "total_cells",
    "marker_score_cp10k_mean",
    "marker_score_z_mean",
)


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


def _mode_or_nan(series: pd.Series):
    """Return the most frequent non-null value, else NaN."""
    values = series.dropna()
    if len(values) == 0:
        return np.nan
    mode = values.mode()
    if len(mode) == 0:
        return np.nan
    return mode.iloc[0]


def _flag_viral(df: pd.DataFrame) -> pd.Series:
    """Flag likely viral epitopes by species text patterns."""
    if "species" not in df.columns:
        return pd.Series(False, index=df.index)
    species_lower = df["species"].fillna("").str.lower()
    is_viral = pd.Series(False, index=df.index)
    for pattern in VIRAL_SPECIES_PATTERNS:
        is_viral |= species_lower.str.contains(pattern, na=False)
    return is_viral


def load_vdjdb(path: str | Path) -> pd.DataFrame:
    """Load VDJdb TSV and normalize core columns."""
    path = Path(path)
    if path.is_dir():
        candidates = list(path.glob("vdjdb*.txt")) + list(path.glob("vdjdb*.tsv"))
        if not candidates:
            raise TCRsiftValidationError(f"No VDJdb files found in directory: {path}")
        db_file = candidates[0]
    else:
        db_file = path

    df = pd.read_csv(db_file, sep="\t", low_memory=False)
    if df.empty:
        raise TCRsiftValidationError(f"VDJdb file is empty: {db_file}")

    column_mapping = {
        "cdr3": "cdr3_beta",
        "cdr3.alpha": "cdr3_alpha",
        "antigen.epitope": "epitope",
        "antigen.gene": "antigen_gene",
        "antigen.species": "species",
        "mhc.a": "mhc_allele",
        "mhc.class": "mhc_class",
        "reference.id": "reference",
    }
    for old, new in column_mapping.items():
        if old in df.columns:
            df[new] = df[old]
    if "cdr3_beta" not in df.columns:
        raise TCRsiftValidationError(
            "VDJdb file missing CDR3 beta column after mapping",
            hint=f"Available columns include: {list(df.columns)[:20]}",
        )

    df["database"] = "VDJdb"
    df["is_viral"] = _flag_viral(df)
    return df


def _load_iedb_cedar_multiheader(path: Path, label: str) -> pd.DataFrame:
    """Load IEDB/CEDAR format with multi-index header."""
    df = pd.read_csv(path, sep="\t", header=[0, 1], low_memory=False)
    if df.empty:
        raise TCRsiftValidationError(f"{label} file is empty: {path}")
    if df.columns.nlevels != 2:
        raise TCRsiftValidationError(
            f"{label} file does not appear to have a two-row header",
            hint="Expected grouped columns like Chain 1 / Chain 2.",
        )

    def _col(group: str, field: str) -> tuple[str, str]:
        return (group, field)

    cdr3_alpha_col = _col("Chain 1", "CDR3 Curated")
    cdr3_beta_col = _col("Chain 2", "CDR3 Curated")
    if cdr3_alpha_col not in df.columns:
        cdr3_alpha_col = _col("Chain 1", "CDR3 Calculated")
    if cdr3_beta_col not in df.columns:
        cdr3_beta_col = _col("Chain 2", "CDR3 Calculated")
    if cdr3_beta_col not in df.columns:
        raise TCRsiftValidationError(f"{label} missing beta CDR3 column in {path}")

    out = pd.DataFrame(
        {
            "cdr3_alpha": df[cdr3_alpha_col] if cdr3_alpha_col in df.columns else "",
            "cdr3_beta": df[cdr3_beta_col],
            "epitope": df.get(_col("Epitope", "Name")),
            "antigen_gene": df.get(_col("Epitope", "Source Molecule")),
            "species": df.get(_col("Epitope", "Source Organism")),
            "chain1_type": df.get(_col("Chain 1", "Type")),
            "chain2_type": df.get(_col("Chain 2", "Type")),
        }
    )
    if "chain2_type" in out.columns:
        mask = out["chain2_type"].fillna("").str.lower().isin(["beta", "trb"])
        out = out[mask]
    out = out.drop(columns=["chain1_type", "chain2_type"])
    out["database"] = label
    out["is_viral"] = _flag_viral(out)
    return out


def load_iedb(path: str | Path) -> pd.DataFrame:
    """Load IEDB TSV."""
    return _load_iedb_cedar_multiheader(Path(path), "IEDB")


def load_cedar(path: str | Path) -> pd.DataFrame:
    """Load CEDAR TSV."""
    return _load_iedb_cedar_multiheader(Path(path), "CEDAR")


def load_databases(
    vdjdb_path: str | Path | None,
    iedb_path: str | Path | None,
    cedar_path: str | Path | None,
) -> pd.DataFrame | None:
    """Load and concatenate annotation databases."""
    dfs = []
    if vdjdb_path:
        dfs.append(load_vdjdb(vdjdb_path))
    if iedb_path:
        dfs.append(load_iedb(iedb_path))
    if cedar_path:
        dfs.append(load_cedar(cedar_path))
    if not dfs:
        return None
    combined = pd.concat(dfs, ignore_index=True)
    combined = combined[combined["cdr3_beta"].notna() & (combined["cdr3_beta"] != "")]
    if "cdr3_alpha" not in combined.columns:
        combined["cdr3_alpha"] = ""
    return combined


def normalize_species_label(species: str | float | None) -> str | None:
    """Normalize species/antigen labels for summary plotting."""
    if pd.isna(species):
        return None
    text = str(species).strip()
    if not text:
        return None
    text = re.sub(r"\s+", " ", text)
    text_no_parens = re.sub(r"[()]", "", text)
    text_no_parens = re.sub(r"\s+", " ", text_no_parens).strip()
    lower = text.lower()
    if (
        "coronavirus" in lower
        or "sars-cov" in lower
        or "sars cov" in lower
        or "severe acute respiratory syndrome" in lower
    ):
        return "Coronavirus SARS-related"
    if "mouse cytomegalovirus" in lower or "murid" in lower:
        return "Mouse CMV"
    if lower == "cmv" or "human cytomegalovirus" in lower or "human herpesvirus 5" in lower:
        return "CMV"
    if lower == "ebv" or "epstein barr" in lower or "human herpesvirus 4" in lower:
        return "EBV"
    if "hepatitis b" in lower or re.search(r"\bhbv\b", lower):
        return "Hepatitis B"
    if "hepatitis c" in lower or re.search(r"\bhcv\b", lower):
        return "Hepatitis C"
    if lower in {"homosapiens", "homo sapiens", "homo sapiens (human)"}:
        return "Homo sapiens"
    if "influenza" in lower:
        return "Influenza"
    if "yellow fever virus" in lower:
        return "Yellow fever virus"
    return text_no_parens


def match_clonotypes(
    clonotypes: pd.DataFrame,
    database: pd.DataFrame,
    match_by: str = "CDR3ab",
) -> pd.DataFrame:
    """Match clonotypes to public database rows using CDR3ab or CDR3b-only rules."""
    if match_by not in {"CDR3ab", "CDR3b_only"}:
        raise ValueError("match_by must be 'CDR3ab' or 'CDR3b_only'")

    df = clonotypes.copy()
    db = database.copy()
    if match_by == "CDR3b_only":
        group_cols = ["cdr3_beta"]
        key_cols = ["CDR3_beta"]
    else:
        group_cols = ["cdr3_alpha", "cdr3_beta"]
        key_cols = ["CDR3_alpha", "CDR3_beta"]

    agg = (
        db.groupby(group_cols, dropna=False)
        .agg(
            db_epitope=("epitope", _mode_or_nan),
            db_species=("species", _mode_or_nan),
            db_database=("database", lambda s: ";".join(sorted(set(s.dropna())))),
            is_viral=("is_viral", "any"),
        )
        .reset_index()
    )
    merged = df.merge(agg, left_on=key_cols, right_on=group_cols, how="left")
    merged["db_match"] = merged["db_database"].notna()
    merged["is_viral"] = merged["is_viral"].astype("boolean").fillna(False).astype(bool)
    drop_cols = [c for c in group_cols if c in merged.columns and c not in key_cols]
    merged = merged.drop(columns=drop_cols)
    return merged


def plot_annotation_summary(df: pd.DataFrame, output_dir: Path) -> None:
    """Write annotation summary plots (pie + top species bar chart)."""
    if "db_match" not in df.columns:
        return
    matched = int(df["db_match"].sum())
    viral = int(df.get("is_viral", pd.Series(False, index=df.index)).sum())
    non_viral_matched = matched - viral
    unmatched = len(df) - matched

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 7))
    sizes = [viral, non_viral_matched, unmatched]
    labels = [f"Viral ({viral})", f"Non-viral matched ({non_viral_matched})", f"Unmatched ({unmatched})"]
    colors = ["#d62728", "#ff7f0e", "#bdbdbd"]
    ax.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90)
    ax.set_title("Database Match Summary")
    fig.savefig(output_dir / "annotation_summary_pie.png", dpi=200)
    plt.close(fig)

    species_col = "db_species_group" if "db_species_group" in df.columns else "db_species"
    if species_col in df.columns:
        species_counts = df[species_col].dropna().value_counts().head(10)
        if len(species_counts) > 0:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.barh(range(len(species_counts)), species_counts.values)
            ax.set_yticks(range(len(species_counts)))
            ax.set_yticklabels(species_counts.index)
            ax.set_xlabel("Number of Clonotypes")
            ax.set_title("Top 10 Matched Species/Antigens")
            ax.invert_yaxis()
            fig.tight_layout()
            fig.savefig(output_dir / "annotation_top_species.png", dpi=200)
            plt.close(fig)


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


def normalize_chain(series: pd.Series) -> pd.Series:
    """Normalize TCR chain labels to TRA/TRB where possible."""
    chain = series.astype(str).str.upper()
    chain = chain.str.replace("-", "", regex=False)
    chain = chain.str.replace(" ", "", regex=False)
    return chain.str.slice(0, 3)


def normalize_clonotype_id(series: pd.Series) -> pd.Series:
    """Normalize clonotype identifiers to a canonical 'clonotype*' suffix."""
    raw = series.astype(str)
    normalized = raw.str.replace(r"(?i).*?(clonotype.*)$", r"\1", regex=True)
    return normalized.where(normalized.str.contains("clonotype", case=False), raw)


def choose_count_column(df: pd.DataFrame, preferred: str | None = None) -> str | None:
    """Choose clonotype cell-count column with v2-compatible priority."""
    cols = {c.lower(): c for c in df.columns}
    if preferred:
        key = preferred.strip().lower()
        if key not in cols:
            raise TCRsiftValidationError(
                f"Requested --count-column '{preferred}' not found",
                hint=f"Available columns: {list(df.columns)}",
            )
        return cols[key]

    candidates = ["cell_count", "n_cells", "num_cells", "cells", "frequency", "count"]
    for cand in candidates:
        if cand in cols:
            return cols[cand]
    return None


def _detect_count_column(clonotypes_df: pd.DataFrame, count_column: str | None = None) -> str:
    """Backwards-compatible wrapper around count-column selection."""
    out = choose_count_column(clonotypes_df, preferred=count_column)
    if out is None:
        raise TCRsiftValidationError(
            "Could not auto-detect clonotype count column",
            hint="Use --count-column to specify a numeric count field.",
        )
    return out


def _validate_count_series(series: pd.Series, count_col: str, source_path: Path) -> pd.Series:
    """Validate and coerce clonotype count series to numeric non-negative integers."""
    numeric = pd.to_numeric(series, errors="coerce")
    non_null = numeric.dropna()
    if non_null.empty:
        raise TCRsiftValidationError(
            f"Column '{count_col}' in {source_path} has no numeric values",
            hint="Pick a discrete count-like column (cell_count, n_cells, frequency, ...).",
        )
    if (non_null < 0).any():
        raise TCRsiftValidationError(
            f"Column '{count_col}' in {source_path} contains negative values",
            hint="Cell counts must be non-negative.",
        )
    is_int_like = np.isclose(non_null, np.round(non_null), atol=1e-6)
    if not is_int_like.all():
        examples = ", ".join(f"{v:.4g}" for v in non_null[~is_int_like].head(3))
        raise TCRsiftValidationError(
            f"Column '{count_col}' in {source_path} contains non-integer values",
            hint=f"Examples: {examples}",
        )
    return numeric


def _aggregate_reads_umis(df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """Aggregate optional reads/umis per normalized clonotype id."""
    agg = {}
    if "umis" in df.columns:
        agg["umis_sum"] = ("umis", "sum")
    if "reads" in df.columns:
        agg["reads_sum"] = ("reads", "sum")
    if not agg:
        return pd.DataFrame({id_col: df[id_col].unique()})
    return df.groupby(id_col, as_index=False).agg(**agg)


def _concat_segment_fields(df: pd.DataFrame, fields: list[str]) -> pd.Series:
    """Concatenate ordered segment fields into a full chain sequence."""
    present = [field for field in fields if field in df.columns]
    if not present:
        return pd.Series([""] * len(df), index=df.index, dtype=object)
    parts = df[present].fillna("").astype(str)
    return parts.agg("".join, axis=1)


def _extract_paired_cdr3_from_consensus_df(
    df: pd.DataFrame, consensus_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract paired TRA/TRB clone rows and rich sequence metadata from consensus CSV."""
    required_cols = {"clonotype_id", "chain", "cdr3"}
    missing = sorted(required_cols - set(df.columns))
    if missing:
        raise TCRsiftValidationError(
            f"Consensus file missing required columns: {missing}",
            hint=f"File: {consensus_path}",
        )

    work = df.copy()
    work["clonotype_id_norm"] = normalize_clonotype_id(work["clonotype_id"])
    work = work.assign(chain_simple=normalize_chain(work["chain"]))
    work = work[work["chain_simple"].isin(["TRA", "TRB"])]

    if "full_length" in work.columns:
        work = work[_truthy_mask(work["full_length"])]
    if "productive" in work.columns:
        work = work[_truthy_mask(work["productive"])]

    umis_col = "umis" if "umis" in work.columns else None
    reads_col = "reads" if "reads" in work.columns else None
    sort_cols = ["clonotype_id_norm", "chain_simple"]
    ascending = [True, True]
    if umis_col:
        sort_cols.append(umis_col)
        ascending.append(False)
    if reads_col:
        sort_cols.append(reads_col)
        ascending.append(False)

    work = work.sort_values(sort_cols, ascending=ascending)
    work = work.groupby(["clonotype_id_norm", "chain_simple"], as_index=False).first()
    work["chain_full_aa"] = _concat_segment_fields(work, list(TCR_SEGMENTS_AA))
    work["chain_full_nt"] = _concat_segment_fields(work, list(TCR_SEGMENTS_NT))

    base = pd.DataFrame(index=work["clonotype_id_norm"].unique())
    value_map = {
        "cdr3": "cdr3_aa",
        "cdr3_nt": "cdr3_nt",
        "chain_full_aa": "full_aa",
        "chain_full_nt": "full_nt",
        "v_gene": "v_gene",
        "d_gene": "d_gene",
        "j_gene": "j_gene",
        "c_gene": "c_gene",
    }
    for segment in TCR_SEGMENTS_AA:
        if segment != "cdr3":
            value_map[segment] = f"{segment}_aa"
    for segment_nt in TCR_SEGMENTS_NT:
        if segment_nt != "cdr3_nt":
            value_map[segment_nt] = segment_nt

    for source_col, suffix in value_map.items():
        if source_col not in work.columns:
            continue
        pivot_part = work.pivot(
            index="clonotype_id_norm",
            columns="chain_simple",
            values=source_col,
        )
        rename_map = {}
        if "TRA" in pivot_part.columns:
            rename_map["TRA"] = f"alpha_{suffix}"
        if "TRB" in pivot_part.columns:
            rename_map["TRB"] = f"beta_{suffix}"
        pivot_part = pivot_part.rename(columns=rename_map)
        base = base.join(pivot_part, how="left")

    base = base.reset_index().rename(columns={"index": "clonotype_id"})
    if "alpha_cdr3_aa" not in base.columns or "beta_cdr3_aa" not in base.columns:
        raise TCRsiftValidationError(
            f"Failed to extract paired alpha/beta CDR3 from {consensus_path}",
            hint="Consensus rows must contain both TRA and TRB per clonotype.",
        )
    base = base.dropna(subset=["alpha_cdr3_aa", "beta_cdr3_aa"])
    base = base[(base["alpha_cdr3_aa"] != "") & (base["beta_cdr3_aa"] != "")]

    base["CDR3_alpha"] = base["alpha_cdr3_aa"]
    base["CDR3_beta"] = base["beta_cdr3_aa"]
    if "alpha_cdr3_nt" in base.columns:
        base["CDR3_alpha_nt"] = base["alpha_cdr3_nt"]
    if "beta_cdr3_nt" in base.columns:
        base["CDR3_beta_nt"] = base["beta_cdr3_nt"]
    base["CDR3ab"] = base["CDR3_alpha"] + "_" + base["CDR3_beta"]

    segment_cols = []
    for chain_name in ("alpha", "beta"):
        for segment in TCR_SEGMENTS_AA:
            segment_cols.append(f"{chain_name}_{segment}_aa")
        for segment_nt in TCR_SEGMENTS_NT:
            segment_cols.append(f"{chain_name}_{segment_nt}")

    ordered_cols = [
        "clonotype_id",
        "CDR3_alpha",
        "CDR3_beta",
        "CDR3_alpha_nt",
        "CDR3_beta_nt",
        "alpha_full_aa",
        "alpha_full_nt",
        "beta_full_aa",
        "beta_full_nt",
        *segment_cols,
        "alpha_v_gene",
        "alpha_j_gene",
        "alpha_c_gene",
        "beta_v_gene",
        "beta_d_gene",
        "beta_j_gene",
        "beta_c_gene",
        "CDR3ab",
    ]
    ordered_cols = [col for col in ordered_cols if col in base.columns]
    return base[ordered_cols], work


def _build_consensus_pairs(consensus_df: pd.DataFrame, require_paired: bool = True) -> pd.DataFrame:
    """Backwards-compatible paired clone extraction helper."""
    pairs, _work = _extract_paired_cdr3_from_consensus_df(consensus_df, Path("<in_memory>"))
    if require_paired:
        return pairs
    return pairs


def load_from_consensus(
    consensus_path: str | Path,
    clonotypes_path: str | Path,
    count_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load one timepoint from consensus+clonotypes and aggregate to unique paired CDR3."""
    consensus_path = validate_file_exists(consensus_path, "consensus annotations")
    clonotypes_path = validate_file_exists(clonotypes_path, "clonotypes")

    consensus_df = pd.read_csv(consensus_path)
    pivot, filtered_consensus = _extract_paired_cdr3_from_consensus_df(consensus_df, Path(consensus_path))

    df_clono = pd.read_csv(clonotypes_path)
    if "clonotype_id" not in df_clono.columns:
        raise TCRsiftValidationError(
            f"Clonotypes file missing required column 'clonotype_id': {clonotypes_path}",
            hint="Expected 10x clonotypes CSV format.",
        )
    df_clono["clonotype_id_norm"] = normalize_clonotype_id(df_clono["clonotype_id"])
    count_col = _detect_count_column(df_clono, count_column=count_column)

    df_counts = df_clono[["clonotype_id_norm", count_col]].rename(
        columns={"clonotype_id_norm": "clonotype_id", count_col: "cell_count"}
    )
    df_counts["cell_count"] = _validate_count_series(
        df_counts["cell_count"], count_col=count_col, source_path=Path(clonotypes_path)
    )

    reads_umis = _aggregate_reads_umis(filtered_consensus, "clonotype_id_norm").rename(
        columns={"clonotype_id_norm": "clonotype_id"}
    )

    merged = pivot.merge(df_counts, on="clonotype_id", how="left").merge(
        reads_umis, on="clonotype_id", how="left"
    )
    merged = merged.dropna(subset=["cell_count"])

    agg = {"cell_count": ("cell_count", "sum")}
    if "umis_sum" in merged.columns:
        agg["umis_sum"] = ("umis_sum", "sum")
    if "reads_sum" in merged.columns:
        agg["reads_sum"] = ("reads_sum", "sum")
    meta_cols = [
        col
        for col in [
            "CDR3_alpha_nt",
            "CDR3_beta_nt",
            "alpha_full_aa",
            "alpha_full_nt",
            "beta_full_aa",
            "beta_full_nt",
            *[f"alpha_{segment}_aa" for segment in TCR_SEGMENTS_AA],
            *[f"beta_{segment}_aa" for segment in TCR_SEGMENTS_AA],
            *[f"alpha_{segment_nt}" for segment_nt in TCR_SEGMENTS_NT],
            *[f"beta_{segment_nt}" for segment_nt in TCR_SEGMENTS_NT],
            "alpha_v_gene",
            "alpha_j_gene",
            "alpha_c_gene",
            "beta_v_gene",
            "beta_d_gene",
            "beta_j_gene",
            "beta_c_gene",
        ]
        if col in merged.columns
    ]
    for col in meta_cols:
        agg[col] = (col, _mode_or_nan)

    counts = merged.groupby(["CDR3_alpha", "CDR3_beta"], as_index=False).agg(**agg)
    counts["CDR3ab"] = counts["CDR3_alpha"] + "_" + counts["CDR3_beta"]
    if counts["CDR3ab"].duplicated().any():
        raise TCRsiftValidationError(
            f"{Path(consensus_path).name}: duplicate CDR3ab rows after aggregation",
            hint="Check clonotype normalization and consensus metadata.",
        )

    out_cols = ["CDR3_alpha", "CDR3_beta", "CDR3ab", "cell_count"] + meta_cols
    out = counts[out_cols].copy()
    out = out.sort_values(["cell_count", "CDR3ab"], ascending=[False, True]).reset_index(drop=True)

    stats = counts.copy()
    return out, stats


def _decode_bytes(values: np.ndarray) -> np.ndarray:
    """Decode byte arrays from HDF5 datasets to UTF-8 strings."""
    if values.dtype.kind in {"S", "O"}:
        return np.array(
            [v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]
        )
    return values.astype(str)


def _load_10x_h5_marker_counts(gex_h5_path: Path, marker_genes: list[str]) -> pd.DataFrame:
    """Load 10x feature-barcode matrix and compute per-cell marker counts + CP10K."""
    marker_genes_upper = [g.upper() for g in marker_genes]
    with h5py.File(gex_h5_path, "r") as h5:
        if "matrix" not in h5:
            raise TCRsiftValidationError(
                f"10x matrix group not found in {gex_h5_path}",
                hint="Expected CellRanger filtered_feature_bc_matrix.h5 format.",
            )
        matrix_group = h5["matrix"]
        for key in ["data", "indices", "indptr", "shape", "barcodes", "features"]:
            if key not in matrix_group:
                raise TCRsiftValidationError(
                    f"10x matrix key '{key}' not found in {gex_h5_path}",
                    hint="Input H5 does not match expected 10x format.",
                )

        barcodes = _decode_bytes(matrix_group["barcodes"][()])
        data = matrix_group["data"][()]
        indices = matrix_group["indices"][()]
        indptr = matrix_group["indptr"][()]
        shape = tuple(matrix_group["shape"][()])
        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)

        features = matrix_group["features"]
        if "name" not in features:
            raise TCRsiftValidationError(
                f"10x features/name not found in {gex_h5_path}",
                hint="Expected CellRanger v3/v4 feature metadata.",
            )
        feature_names = _decode_bytes(features["name"][()])
        feature_types = None
        if "feature_type" in features:
            feature_types = _decode_bytes(features["feature_type"][()])

    gene_to_index: dict[str, int] = {}
    for idx, name in enumerate(feature_names):
        if feature_types is not None and feature_types[idx] != "Gene Expression":
            continue
        key = str(name).strip().upper()
        if key and key not in gene_to_index:
            gene_to_index[key] = idx

    out = pd.DataFrame({"barcode": barcodes})
    total_umis = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    out["total_umis"] = total_umis
    total_safe = np.where(total_umis > 0, total_umis, np.nan)

    for raw_gene, gene_upper in zip(marker_genes, marker_genes_upper):
        col_counts = f"count_{raw_gene}"
        col_cp10k = f"cp10k_{raw_gene}"
        idx = gene_to_index.get(gene_upper)
        if idx is None:
            out[col_counts] = 0.0
            out[col_cp10k] = 0.0
            continue
        counts = np.asarray(matrix.getrow(idx).toarray()).ravel().astype(float)
        out[col_counts] = counts
        out[col_cp10k] = np.where(total_safe > 0, (counts / total_safe) * 10000.0, 0.0)

    return out


def _load_barcode_to_clonotype(contig_csv_path: Path) -> pd.DataFrame:
    """Load one barcode->normalized clonotype mapping from contig annotations."""
    contigs = pd.read_csv(contig_csv_path)
    if "barcode" not in contigs.columns:
        raise TCRsiftValidationError(
            f"Contig file missing required column 'barcode': {contig_csv_path}",
            hint="Expected CellRanger filtered_contig_annotations CSV.",
        )

    clonotype_col = None
    for candidate in ["raw_clonotype_id", "clonotype_id"]:
        if candidate in contigs.columns:
            clonotype_col = candidate
            break
    if clonotype_col is None:
        raise TCRsiftValidationError(
            "Contig file missing clonotype ID column",
            hint=f"Expected raw_clonotype_id or clonotype_id in {contig_csv_path}",
        )

    work = contigs[["barcode", clonotype_col]].copy()
    work = work.dropna(subset=[clonotype_col])
    work = work[work[clonotype_col].astype(str).str.len() > 0]
    work["clonotype_id"] = normalize_clonotype_id(work[clonotype_col])
    work = work.drop(columns=[clonotype_col])
    work = work.dropna(subset=["barcode", "clonotype_id"])
    work = work.drop_duplicates(subset=["barcode", "clonotype_id"])
    work = work.groupby("barcode", as_index=False).first()
    return work


def _compute_marker_scores_from_adata(
    consensus_path: Path,
    contig_csv_path: Path,
    marker_genes: list[str],
    adata,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fallback marker-score path using scanpy AnnData (for tests/non-10x H5)."""
    consensus_df = pd.read_csv(consensus_path)
    clone_map, _ = _extract_paired_cdr3_from_consensus_df(consensus_df, consensus_path)
    clone_map = clone_map[["clonotype_id", "CDR3_alpha", "CDR3_beta", "CDR3ab"]].drop_duplicates()

    contigs = pd.read_csv(contig_csv_path)
    clonotype_col = "raw_clonotype_id" if "raw_clonotype_id" in contigs.columns else "clonotype_id"
    if clonotype_col not in contigs.columns:
        raise TCRsiftValidationError(
            f"Contig file missing clonotype id column: {contig_csv_path}",
            hint="Expected raw_clonotype_id or clonotype_id.",
        )
    keep = contigs[[clonotype_col, "barcode"]].copy()
    if "productive" in contigs.columns:
        keep = keep[_truthy_mask(contigs["productive"])]
    if "high_confidence" in contigs.columns:
        keep = keep[_truthy_mask(contigs["high_confidence"])]
    keep = keep.rename(columns={clonotype_col: "clonotype_id"}).drop_duplicates()
    keep["clonotype_id"] = normalize_clonotype_id(keep["clonotype_id"])
    keep = keep.merge(clone_map, on="clonotype_id", how="inner")
    if len(keep) == 0:
        raise TCRsiftValidationError(
            f"No contig barcodes mapped to paired clonotypes for {contig_csv_path}",
            hint="Check barcode/clonotype columns and consensus file consistency.",
        )

    adata_set = set(adata.obs_names.tolist())
    keep = keep[keep["barcode"].isin(adata_set)]
    if len(keep) == 0:
        raise TCRsiftValidationError(
            f"No barcodes from {contig_csv_path.name} were found in GEX matrix",
            hint="Check if VDJ and GEX barcodes use matching suffix conventions.",
        )
    keep = keep.drop_duplicates(subset=["barcode", "clonotype_id"]).reset_index(drop=True)
    adata_subset = adata[keep["barcode"].tolist(), :].copy()

    if issparse(adata_subset.X):
        total_umis = np.asarray(adata_subset.X.sum(axis=1)).ravel().astype(float)
    else:
        total_umis = np.asarray(adata_subset.X.sum(axis=1)).ravel().astype(float)
    total_umis = np.maximum(total_umis, 1.0)

    cell_df = keep.copy()
    cell_df["total_umis"] = total_umis
    for gene in marker_genes:
        if gene in adata_subset.var_names:
            vec = adata_subset[:, gene].X
            if issparse(vec):
                counts = np.asarray(vec.toarray()).ravel().astype(float)
            else:
                counts = np.asarray(vec).ravel().astype(float)
        else:
            counts = np.zeros(len(cell_df), dtype=float)
        cp10k = counts / total_umis * 1e4
        log1p = np.log1p(cp10k)
        std = float(np.std(log1p, ddof=0))
        z = np.zeros(len(log1p), dtype=float) if std <= 1e-12 else (log1p - float(np.mean(log1p))) / std
        cell_df[f"count_{gene}"] = counts
        cell_df[f"cp10k_{gene}"] = cp10k
        cell_df[f"log1p_cp10k_{gene}"] = log1p
        cell_df[f"z_log1p_cp10k_{gene}"] = z

    cp10k_cols = [f"cp10k_{g}" for g in marker_genes]
    z_cols = [f"z_log1p_cp10k_{g}" for g in marker_genes]
    cell_df["marker_mean_cp10k"] = cell_df[cp10k_cols].mean(axis=1)
    cell_df["marker_mean_z"] = cell_df[z_cols].mean(axis=1)

    agg = {
        "n_cells_with_gex": ("barcode", "nunique"),
        "marker_score_cp10k": ("marker_mean_cp10k", "mean"),
        "marker_score_z": ("marker_mean_z", "mean"),
    }
    for gene in marker_genes:
        agg[f"score_count_{gene}"] = (f"count_{gene}", "mean")
        agg[f"score_cp10k_{gene}"] = (f"cp10k_{gene}", "mean")
        agg[f"score_log1p_cp10k_{gene}"] = (f"log1p_cp10k_{gene}", "mean")
        agg[f"score_z_{gene}"] = (f"z_log1p_cp10k_{gene}", "mean")
    clonotype_scores = (
        cell_df.groupby(["CDR3_alpha", "CDR3_beta", "CDR3ab"], as_index=False).agg(**agg)
    )
    return cell_df, clonotype_scores


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
    consensus_path = Path(validate_file_exists(consensus_path, "consensus annotations"))
    contig_csv_path = Path(validate_file_exists(contig_csv_path, "filtered contig annotations"))
    gex_h5_path = Path(validate_file_exists(gex_h5_path, "10x filtered_feature_bc_matrix.h5"))

    marker_genes = [g.strip().upper() for g in marker_genes if g.strip()]
    marker_genes = list(dict.fromkeys(marker_genes))

    try:
        marker_cells = _load_10x_h5_marker_counts(gex_h5_path, marker_genes)
        barcode_map = _load_barcode_to_clonotype(contig_csv_path)

        consensus_df = pd.read_csv(consensus_path)
        clonotype_pairs, _ = _extract_paired_cdr3_from_consensus_df(consensus_df, consensus_path)

        merged = marker_cells.merge(barcode_map, on="barcode", how="inner")
        merged = merged.merge(clonotype_pairs, on="clonotype_id", how="inner")
        if merged.empty:
            raise TCRsiftValidationError(
                "No overlap between GEX barcodes and paired clonotypes",
                hint=f"Timepoint from {consensus_path.name} produced zero linked rows.",
            )

        cp10k_cols = [f"cp10k_{g}" for g in marker_genes]
        merged["marker_mean_cp10k"] = merged[cp10k_cols].mean(axis=1)
        z_cols: list[str] = []
        for gene in marker_genes:
            cp10k_col = f"cp10k_{gene}"
            log_col = f"log1p_cp10k_{gene}"
            z_col = f"z_log1p_cp10k_{gene}"
            merged[log_col] = np.log1p(merged[cp10k_col])
            mean_val = float(merged[log_col].mean())
            std_val = float(merged[log_col].std(ddof=0))
            if std_val <= 1e-12:
                merged[z_col] = 0.0
            else:
                merged[z_col] = (merged[log_col] - mean_val) / std_val
            z_cols.append(z_col)
        merged["marker_mean_z"] = merged[z_cols].mean(axis=1)

        agg = {
            "n_cells_with_gex": ("barcode", "nunique"),
            "marker_score_cp10k": ("marker_mean_cp10k", "mean"),
            "marker_score_z": ("marker_mean_z", "mean"),
        }
        for gene in marker_genes:
            agg[f"score_count_{gene}"] = (f"count_{gene}", "mean")
            agg[f"score_cp10k_{gene}"] = (f"cp10k_{gene}", "mean")
            agg[f"score_log1p_cp10k_{gene}"] = (f"log1p_cp10k_{gene}", "mean")
            agg[f"score_z_{gene}"] = (f"z_log1p_cp10k_{gene}", "mean")

        clonotype_scores = (
            merged.groupby(["CDR3_alpha", "CDR3_beta", "CDR3ab"], as_index=False).agg(**agg)
        )
        return merged, clonotype_scores
    except Exception:
        # Test-compatible fallback path when H5 is mocked or non-standard.
        adata = sc.read_10x_h5(str(gex_h5_path))
        return _compute_marker_scores_from_adata(
            consensus_path=consensus_path,
            contig_csv_path=contig_csv_path,
            marker_genes=marker_genes,
            adata=adata,
        )


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
    """Add timepoint summary columns using legacy v2 column order."""
    out = df.copy()
    summary_cols: dict[str, pd.Series] = {}
    for base in metric_bases:
        cols = [f"{base}_{tp}" for tp in timepoint_order if f"{base}_{tp}" in out.columns]
        if not cols:
            continue
        numeric = out[cols].apply(pd.to_numeric, errors="coerce")
        for stat_name, values in [
            ("mean", numeric.mean(axis=1)),
            ("sum", numeric.sum(axis=1)),
            ("min", numeric.min(axis=1)),
            ("max", numeric.max(axis=1)),
        ]:
            col_name = f"{base}_{stat_name}"
            if col_name in out.columns or col_name in summary_cols:
                continue
            summary_cols[col_name] = values
    if summary_cols:
        out = pd.concat([out, pd.DataFrame(summary_cols, index=out.index)], axis=1)
    return out


def add_marker_scores_to_harmonized(
    harmonized: pd.DataFrame,
    marker_scores_by_timepoint: dict[str, pd.DataFrame],
    timepoint_order: list[str],
) -> pd.DataFrame:
    """Merge marker score tables into harmonized clone table (v2-compatible)."""
    out = harmonized.copy()
    cp10k_score_cols = []
    z_score_cols = []
    added_cols: list[str] = []
    for tp in timepoint_order:
        marker_df = marker_scores_by_timepoint.get(tp)
        if marker_df is None or marker_df.empty:
            continue
        metric_cols = [col for col in marker_df.columns if col not in ["CDR3_alpha", "CDR3_beta", "CDR3ab"]]
        rename_map = {col: f"{col}_{tp}" for col in metric_cols}
        to_merge = marker_df[["CDR3ab"] + metric_cols].rename(columns=rename_map)
        out = out.merge(to_merge, on="CDR3ab", how="left")
        added_cols.extend(rename_map.values())
        if "marker_score_cp10k" in marker_df.columns:
            cp10k_score_cols.append(f"marker_score_cp10k_{tp}")
        if "marker_score_z" in marker_df.columns:
            z_score_cols.append(f"marker_score_z_{tp}")

    if added_cols:
        for col in added_cols:
            out[col] = out[col].fillna(0.0)
        out = out.copy()
    if cp10k_score_cols:
        out["marker_score_cp10k_mean"] = out[cp10k_score_cols].mean(axis=1)
        out["marker_score_cp10k_max"] = out[cp10k_score_cols].max(axis=1)
    if z_score_cols:
        out["marker_score_z_mean"] = out[z_score_cols].mean(axis=1)
        out["marker_score_z_max"] = out[z_score_cols].max(axis=1)
    return out


def sort_harmonized_by_rank(df: pd.DataFrame, rank_by: str = "mean_frequency") -> pd.DataFrame:
    """Sort harmonized clones by rank metric using v2-compatible tiebreaks."""
    if rank_by not in RANK_METRICS:
        raise ValueError(f"rank_by must be one of: {', '.join(RANK_METRICS)}")
    if rank_by == "marker_score_z_mean":
        sort_cols = ["marker_score_z_mean", "marker_score_cp10k_mean", "mean_frequency", "total_cells"]
    elif rank_by == "marker_score_cp10k_mean":
        sort_cols = ["marker_score_cp10k_mean", "mean_frequency", "total_cells", "n_timepoints"]
    elif rank_by == "total_cells":
        sort_cols = ["total_cells", "mean_frequency", "n_timepoints"]
    elif rank_by == "max_frequency":
        sort_cols = ["max_frequency", "mean_frequency", "total_cells", "n_timepoints"]
    else:
        sort_cols = ["mean_frequency", "max_frequency", "total_cells", "n_timepoints"]
    sort_cols = [c for c in sort_cols if c in df.columns]
    if not sort_cols:
        return df.copy().reset_index(drop=True)
    return df.sort_values(sort_cols, ascending=[False] * len(sort_cols)).reset_index(drop=True)


def build_harmonized_table(
    samples: dict[str, pd.DataFrame],
    timepoint_order: list[str],
    rank_by: str = "mean_frequency",
) -> pd.DataFrame:
    """Build harmonized clone table across timepoints (v2-compatible join semantics)."""
    if rank_by not in RANK_METRICS:
        raise ValueError(f"rank_by must be one of: {', '.join(RANK_METRICS)}")

    key_cols = ["CDR3_alpha", "CDR3_beta", "CDR3ab"]
    merged = None
    count_cols = []
    metadata_frames: list[pd.DataFrame] = []
    for label in timepoint_order:
        if label not in samples:
            raise ValueError(f"Missing sample for timepoint: {label}")
        df = samples[label]
        meta_cols = [c for c in df.columns if c not in key_cols + ["cell_count"]]
        if meta_cols:
            metadata_frames.append(df[key_cols + meta_cols].copy())
        df = df[key_cols + ["cell_count"]].copy()
        count_col = f"cell_count_{label}"
        count_cols.append(count_col)
        df = df.rename(columns={"cell_count": count_col})
        if merged is None:
            merged = df
        else:
            merged = pd.merge(merged, df, on=key_cols, how="outer")

    if merged is None:
        raise ValueError("No samples provided")

    for col in count_cols:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype(int)

    if metadata_frames:
        metadata_all = pd.concat(metadata_frames, ignore_index=True)
        meta_cols_union = [c for c in metadata_all.columns if c not in key_cols]
        metadata_agg = {col: (col, _mode_or_nan) for col in meta_cols_union}
        metadata_unique = metadata_all.groupby(key_cols, as_index=False).agg(**metadata_agg)
        merged = merged.merge(metadata_unique, on=key_cols, how="left")

    freq_cols = []
    for label in timepoint_order:
        count_col = f"cell_count_{label}"
        if count_col not in merged.columns:
            continue
        freq_col = f"frequency_{label}"
        total = merged[count_col].sum()
        if total <= 0:
            merged[freq_col] = 0.0
        else:
            merged[freq_col] = merged[count_col] / total
        freq_cols.append(freq_col)

    merged["total_cells"] = merged[count_cols].sum(axis=1)
    merged["n_timepoints"] = (merged[count_cols] > 0).sum(axis=1)
    merged["is_paired_alpha_beta"] = True
    if freq_cols:
        merged["mean_frequency"] = merged[freq_cols].mean(axis=1)
        merged["max_frequency"] = merged[freq_cols].max(axis=1)

    if merged["CDR3ab"].duplicated().any():
        n_dup = int(merged["CDR3ab"].duplicated().sum())
        raise ValueError(f"harmonized: found {n_dup} duplicate CDR3ab rows")
    return sort_harmonized_by_rank(merged, rank_by=rank_by)


def compute_increasing_masks_from_frequencies(
    freq_df: pd.DataFrame,
    ratio_nonzero_min: float = 1.5,
    ratio_all_timepoints_min: float = 1.25,
) -> tuple[pd.Series, pd.Series]:
    """Compute increasing masks using unified v2 rule."""
    if ratio_nonzero_min <= 0:
        raise ValueError(f"ratio_nonzero_min must be > 0, got {ratio_nonzero_min}")
    if ratio_all_timepoints_min <= 0:
        raise ValueError(f"ratio_all_timepoints_min must be > 0, got {ratio_all_timepoints_min}")
    if freq_df.empty:
        empty = pd.Series(False, index=freq_df.index)
        return empty, empty
    if freq_df.shape[1] < 2:
        empty = pd.Series(False, index=freq_df.index)
        return empty, empty

    freq_values = freq_df.apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    n_rows = freq_values.shape[0]
    increasing = np.zeros(n_rows, dtype=bool)
    tol = 1e-12

    for i in range(n_rows):
        row = np.nan_to_num(freq_values[i], nan=0.0)
        last = row[-1]
        if last <= 0.0:
            continue
        baseline = float(np.max(row[:2])) if row.size >= 3 else float(np.max(row[:-1]))
        if baseline <= 0.0:
            continue
        if last + tol >= ratio_nonzero_min * baseline:
            increasing[i] = True

    mask = pd.Series(increasing, index=freq_df.index)
    return mask, mask.copy()


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
    enrichment_genes_preferred: list[str] | None = None,
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
    df["is_increasing_positive"] = all_positive_increasing
    df["is_increasing_nonzero"] = nonzero_increasing
    df["is_increasing_all_timepoints"] = all_positive_increasing

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
    df["is_enriched"] = df["is_increasing_nonzero"]

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
    df["is_cd8_positive"] = df["cd8_cp10k_max"] > float(min_cd8_cp10k)

    cd4_candidate_cols = [c for c in df.columns if c.startswith("score_cp10k_CD4")]
    if cd4_candidate_cols:
        df["cd4_cp10k_max"] = (
            df[cd4_candidate_cols].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)
        )
    else:
        df["cd4_cp10k_max"] = 0.0

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

    enrichment_seed = (
        [g.strip().upper() for g in enrichment_genes_preferred if g.strip()]
        if enrichment_genes_preferred
        else list(ENRICHMENT_GENES_DEFAULT)
    )
    enrichment_genes = _dedupe_valid(enrichment_seed, ENRICHMENT_GENES_DEFAULT)

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
    _add_panel_scores("enrichment", enrichment_genes)

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
    n_top_enrichment = _select_top_panel_score(
        df,
        "enrichment_score_cp10k_mean",
        "is_top_enrichment_score",
        immunogenic_percentile,
        immunogenic_percentile_slack_frac,
        immunogenic_min_cp10k,
        immunogenic_require_above_median,
    )

    df["is_branch_top_immunogenic"] = df["is_top_immunogenic_any"]
    df["is_branch_cytolytic"] = df["is_top_cytolytic_score"]
    df["is_branch_antigen_response"] = df["is_top_antigen_response_score"]
    df["is_branch_enrichment_markers"] = df["is_top_enrichment_score"]
    df["is_branch_enriched"] = df["is_increasing_nonzero"]
    df["is_branch_increasing"] = df["is_increasing_all_timepoints"]
    df["is_branch_any"] = (
        df["is_branch_top_immunogenic"]
        | df["is_branch_cytolytic"]
        | df["is_branch_antigen_response"]
        | df["is_branch_enrichment_markers"]
        | df["is_branch_enriched"]
        | df["is_branch_increasing"]
    )
    df["is_immunogenic"] = (
        df["is_top_immunogenic_any"]
        | df["is_top_cytolytic_score"]
        | df["is_top_antigen_response_score"]
        | df["is_top_enrichment_score"]
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
        "subset_top_enrichment_score": df[df["is_top_enrichment_score"]].copy(),
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
        ("Enrich", n_top_enrichment),
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
        ("Enrich", int(df["is_top_enrichment_score"].sum())),
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


def select_top_monotonic_increase_clones(
    df: pd.DataFrame,
    timepoint_order: list[str],
    top_k: int,
    increase_ratio_nonzero_min: float = 1.5,
    increase_ratio_all_timepoints_min: float = 1.25,
    require_positive_all_timepoints: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Select top increasing clones for monotonic-style output tables."""
    freq_cols = [f"frequency_{tp}" for tp in timepoint_order if f"frequency_{tp}" in df.columns]
    if len(freq_cols) < 2:
        return pd.DataFrame(), freq_cols

    mask_nonzero, _mask_all = compute_increasing_masks_from_frequencies(
        df[freq_cols],
        ratio_nonzero_min=increase_ratio_nonzero_min,
        ratio_all_timepoints_min=increase_ratio_all_timepoints_min,
    )
    mask = mask_nonzero
    if require_positive_all_timepoints:
        # Unified v2 rule keeps this alias behavior identical.
        mask = mask_nonzero
    selected = df.loc[mask].copy()
    if selected.empty:
        return selected, freq_cols

    terminal_col = freq_cols[-1]
    sort_cols = [c for c in [terminal_col, "mean_frequency", "total_cells"] if c in selected.columns]
    if sort_cols:
        selected = selected.sort_values(sort_cols, ascending=[False] * len(sort_cols))
    return selected.head(top_k), freq_cols


def plot_frequency_heatmap(
    df: pd.DataFrame,
    freq_cols: list[str],
    output_path: Path,
    title: str,
) -> None:
    """Plot frequency (%) heatmap for selected clones."""
    if df.empty or not freq_cols:
        return
    matrix = df.set_index("CDR3ab")[freq_cols].copy() * 100.0
    matrix.columns = [c.replace("frequency_", "") for c in matrix.columns]
    fig, ax = plt.subplots(
        figsize=(max(6, 1.5 * len(freq_cols) + 2), max(4, 0.4 * len(df) + 2))
    )
    sns.heatmap(
        matrix,
        cmap="magma",
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Frequency (%)"},
        ax=ax,
    )
    ax.set_xlabel("Time point")
    ax.set_ylabel("CDR3ab")
    ax.set_title(title)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_monotonic_increase_log_lines(
    top_df: pd.DataFrame,
    freq_cols: list[str],
    timepoint_order: list[str],
    output_path: Path,
) -> None:
    """Plot log-scale line trajectories for increasing clones."""
    if top_df.empty or not freq_cols:
        return
    long_df = top_df[["CDR3ab"] + freq_cols].melt(
        id_vars="CDR3ab",
        value_vars=freq_cols,
        var_name="frequency_col",
        value_name="frequency",
    )
    long_df["timepoint"] = long_df["frequency_col"].str.replace("frequency_", "", regex=False)
    long_df["timepoint"] = pd.Categorical(
        long_df["timepoint"], categories=timepoint_order, ordered=True
    )
    long_df = long_df.sort_values(["CDR3ab", "timepoint"])
    epsilon = 1e-6
    long_df["frequency_plot"] = long_df["frequency"].clip(lower=epsilon) * 100.0

    fig, ax = plt.subplots(figsize=(10, 6))
    for _, group in long_df.groupby("CDR3ab", sort=False):
        ax.plot(
            group["timepoint"].astype(str),
            group["frequency_plot"],
            marker="o",
            linewidth=1.2,
            alpha=0.8,
        )
    ax.set_yscale("log")
    ax.set_xlabel("Time point")
    ax.set_ylabel("Frequency (%) [log scale]")
    ax.set_title("Increasing Clones: Frequency Trajectories")
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_annotated_heatmap(
    df: pd.DataFrame,
    count_cols: list[str],
    output_path: Path,
    top_k: int,
    annotation_cols: list[str],
    title: str,
) -> None:
    """Plot top-k heatmap with optional annotation sidecar columns."""
    if not count_cols:
        return
    top = df.head(top_k).copy()
    if top.empty:
        return
    matrix = top.set_index("CDR3ab")[count_cols].copy()
    matrix.columns = [c.replace("cell_count_", "") for c in matrix.columns]
    fig, ax = plt.subplots(
        figsize=(max(7, 1.5 * len(count_cols) + 2), max(4, 0.4 * len(top) + 2))
    )
    sns.heatmap(matrix, cmap="YlOrRd", ax=ax, cbar_kws={"label": "Cell count"})
    ax.set_title(title)
    ax.set_xlabel("Time point")
    ax.set_ylabel("CDR3ab")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _save_placeholder_plot(path: Path, title: str, body: str = "") -> None:
    """Write a lightweight placeholder plot artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(8, 4))
    plt.axis("off")
    plt.title(title)
    if body:
        plt.text(0.5, 0.5, body, ha="center", va="center")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def export_all_plots_pdf(fig_dir: Path, output_pdf: Path) -> None:
    """Compile all PNG figures into a single PDF."""
    png_paths = sorted([p for p in fig_dir.glob("*.png") if p.is_file()])
    if not png_paths:
        return
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output_pdf) as pdf:
        for path in png_paths:
            img = plt.imread(path)
            fig = plt.figure(figsize=(11, 8.5))
            plt.imshow(img)
            plt.axis("off")
            plt.title(path.name)
            pdf.savefig(fig)
            plt.close(fig)


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
    stats_by_timepoint: dict[str, pd.DataFrame] = {}
    for tp in timepoint_order:
        counts_df, stats_df = load_from_consensus(
            discovered[tp]["consensus"],
            discovered[tp]["clonotypes"],
            count_column=getattr(args, "count_column", None),
        )
        samples[tp] = counts_df
        stats_by_timepoint[tp] = stats_df

    marker_genes = [g.strip().upper() for g in str(args.marker_genes).split(",") if g.strip()]
    marker_genes = list(dict.fromkeys(marker_genes))
    if not marker_genes:
        marker_genes = MARKER_GENES_DEFAULT.copy()

    immunogenic_genes = [g.strip().upper() for g in str(args.immunogenic_genes).split(",") if g.strip()]
    cytotoxic_genes = [g.strip().upper() for g in str(args.cytotoxic_genes).split(",") if g.strip()]
    cytolytic_genes = [g.strip().upper() for g in str(args.cytolytic_genes).split(",") if g.strip()]
    antigen_response_genes = [
        g.strip().upper() for g in str(args.antigen_response_genes).split(",") if g.strip()
    ]
    enrichment_genes = [
        g.strip().upper() for g in str(args.enrichment_genes).split(",") if g.strip()
    ]

    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_heatmap = Path(args.out_heatmap) if args.out_heatmap else (fig_dir / "abTCR_topk_heatmap.png")
    out_top = Path(args.out_top) if args.out_top else (fig_dir / "abTCR_topk.csv")
    out_annotated = Path(args.out_annotated) if args.out_annotated else (fig_dir / "abTCR_annotated.csv")
    out_annotated_heatmap = (
        Path(args.out_annotated_heatmap)
        if args.out_annotated_heatmap
        else (fig_dir / "abTCR_topk_annotated_heatmap.png")
    )
    out_selected_report = (
        Path(args.out_selected_report)
        if args.out_selected_report
        else (fig_dir / "selected_clones_report.pdf")
    )

    gex_inputs = discover_gex_inputs(args.data_dir, timepoint_order)
    marker_cell_by_timepoint: dict[str, pd.DataFrame] = {}
    marker_scores_by_timepoint: dict[str, pd.DataFrame] = {}
    for tp in timepoint_order:
        cell_df, score_df = compute_marker_scores_for_timepoint(
            consensus_path=discovered[tp]["consensus"],
            contig_csv_path=gex_inputs[tp]["contigs"],
            gex_h5_path=gex_inputs[tp]["gex_h5"],
            marker_genes=marker_genes,
        )
        marker_cell_by_timepoint[tp] = cell_df
        marker_scores_by_timepoint[tp] = score_df
        cell_df.to_csv(fig_dir / f"marker_cells_{tp}.csv", index=False)
        score_df.to_csv(fig_dir / f"marker_clonotype_scores_{tp}.csv", index=False)

    harmonized = build_harmonized_table(samples, timepoint_order, rank_by="mean_frequency")
    harmonized = add_marker_scores_to_harmonized(
        harmonized, marker_scores_by_timepoint, timepoint_order
    )

    annotated = harmonized.copy()
    database = load_databases(getattr(args, "vdjdb", None), getattr(args, "iedb", None), getattr(args, "cedar", None))
    if database is not None:
        annotated = match_clonotypes(harmonized, database, match_by=getattr(args, "match_by", "CDR3b_only"))
        if "db_species" in annotated.columns:
            annotated["db_species_group"] = annotated["db_species"].apply(normalize_species_label)
        annotated["db_vdjdb"] = annotated["db_database"].fillna("").str.contains("VDJdb")
        annotated["db_iedb"] = annotated["db_database"].fillna("").str.contains("IEDB")
        annotated["db_cedar"] = annotated["db_database"].fillna("").str.contains("CEDAR")
        plot_annotation_summary(annotated, fig_dir)

    metric_bases = infer_timepoint_metric_bases(annotated, timepoint_order)
    master_df = add_timepoint_summaries(annotated, timepoint_order, metric_bases)
    (
        master_df,
        subset_dfs,
        _sequential_stages,
        _independent_stages,
        branch_union_stages,
    ) = run_selection_pipeline(
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
        enrichment_genes_preferred=enrichment_genes,
    )

    _write_subset_tables(subset_dfs, fig_dir)
    plot_selection_funnel(branch_union_stages, fig_dir / "selection_funnel.png", "Selection Funnel")

    plot_source = sort_harmonized_by_rank(master_df, args.rank_by)
    out_table = Path(args.out_table)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    plot_source.to_csv(out_table, index=False)
    plot_source.to_csv(fig_dir / "abTCR_master_table.csv", index=False)
    plot_source.to_csv(out_annotated, index=False)

    mask_cols = sorted([c for c in plot_source.columns if c.startswith("is_")])
    mask_id_cols = [c for c in ["CDR3ab", "CDR3_alpha", "CDR3_beta"] if c in plot_source.columns]
    mask_metric_cols = [
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
        if c in plot_source.columns
    ]
    plot_source[mask_id_cols + mask_metric_cols + mask_cols].to_csv(
        fig_dir / "selection_masks.csv", index=False
    )

    count_cols = [f"cell_count_{tp}" for tp in timepoint_order if f"cell_count_{tp}" in plot_source.columns]
    freq_cols = [f"frequency_{tp}" for tp in timepoint_order if f"frequency_{tp}" in plot_source.columns]

    selected_mask = pd.Series(
        plot_source.get("is_candidate_tumor_reactive", pd.Series(False, index=plot_source.index)),
        index=plot_source.index,
    ).fillna(False).astype(bool)
    selected_for_report = plot_source[selected_mask].copy()
    create_selected_clone_pdf_report(
        selected_for_report,
        out_selected_report,
        timepoint_order,
        marker_genes,
        pyensembl_release=int(args.pyensembl_release),
    )

    primary_title = f"Top {args.top_k} abTCRs (ranked by {args.rank_by})"
    plot_heatmap(plot_source, count_cols, out_heatmap, int(args.top_k), title=primary_title)
    plot_source.head(int(args.top_k)).to_csv(out_top, index=False)

    for rank_metric in RANK_METRICS:
        ranked = sort_harmonized_by_rank(plot_source, rank_metric)
        ranked.head(int(args.top_k)).to_csv(fig_dir / f"abTCR_topk_{rank_metric}.csv", index=False)
        plot_heatmap(
            ranked,
            count_cols,
            fig_dir / f"abTCR_topk_heatmap_{rank_metric}.png",
            int(args.top_k),
            title=f"Top {args.top_k} abTCRs (ranked by {rank_metric})",
        )

    # Increasing + became-cytotoxic topk artifacts
    priority_mask = pd.Series(
        plot_source.get("is_priority_cytotoxic_expander", pd.Series(False, index=plot_source.index)),
        index=plot_source.index,
    ).fillna(False).astype(bool)
    priority_df = plot_source[priority_mask].copy()
    priority_top_path = fig_dir / "abTCR_increasing_became_cytotoxic_topk.csv"
    if priority_df.empty:
        plot_source.head(0).to_csv(priority_top_path, index=False)
    else:
        terminal_freq_col = f"frequency_{timepoint_order[-1]}"
        sort_cols = [
            c
            for c in [terminal_freq_col, "cytotoxic_score_z_delta_last_vs_first", "mean_frequency", "total_cells"]
            if c in priority_df.columns
        ]
        priority_top = priority_df.sort_values(sort_cols, ascending=[False] * len(sort_cols)).head(int(args.top_k))
        priority_top.to_csv(priority_top_path, index=False)
        if len(freq_cols) >= 2:
            plot_frequency_heatmap(
                priority_top,
                freq_cols,
                fig_dir / "abTCR_increasing_became_cytotoxic_topk_heatmap.png",
                "Top increasing clones that became cytotoxic",
            )

    # Monotonic-style increasing outputs
    monotonic_top, monotonic_freq_cols = select_top_monotonic_increase_clones(
        plot_source,
        timepoint_order,
        int(args.top_k),
        increase_ratio_nonzero_min=float(args.increase_ratio_nonzero_min),
        increase_ratio_all_timepoints_min=float(args.increase_ratio_all_timepoints_min),
    )
    if monotonic_top.empty:
        plot_source.head(0).to_csv(fig_dir / "abTCR_monotonic_increase_topk.csv", index=False)
    else:
        monotonic_top.to_csv(fig_dir / "abTCR_monotonic_increase_topk.csv", index=False)
        plot_frequency_heatmap(
            monotonic_top,
            monotonic_freq_cols,
            fig_dir / "abTCR_monotonic_increase_topk_heatmap.png",
            "Top increasing clones",
        )
        plot_monotonic_increase_log_lines(
            monotonic_top,
            monotonic_freq_cols,
            timepoint_order,
            fig_dir / "abTCR_monotonic_increase_topk_log_lines.png",
        )

    monotonic_pos_top, monotonic_pos_freq_cols = select_top_monotonic_increase_clones(
        plot_source,
        timepoint_order,
        int(args.top_k),
        increase_ratio_nonzero_min=float(args.increase_ratio_nonzero_min),
        increase_ratio_all_timepoints_min=float(args.increase_ratio_all_timepoints_min),
        require_positive_all_timepoints=True,
    )
    monotonic_pos_top.to_csv(fig_dir / "abTCR_monotonic_increase_positive_topk.csv", index=False)
    if not monotonic_pos_top.empty:
        plot_frequency_heatmap(
            monotonic_pos_top,
            monotonic_pos_freq_cols,
            fig_dir / "abTCR_monotonic_increase_positive_topk_heatmap.png",
            "Top increasing clones (positive timepoints)",
        )
        plot_monotonic_increase_log_lines(
            monotonic_pos_top,
            monotonic_pos_freq_cols,
            timepoint_order,
            fig_dir / "abTCR_monotonic_increase_positive_topk_log_lines.png",
        )

    annotation_cols = [c for c in ["db_vdjdb", "db_iedb", "db_cedar", "is_viral"] if c in plot_source.columns]
    if annotation_cols:
        plot_annotated_heatmap(
            plot_source,
            count_cols,
            out_annotated_heatmap,
            int(args.top_k),
            annotation_cols,
            title=f"Top {args.top_k} abTCRs (ranked by {args.rank_by}, annotated)",
        )
        for rank_metric in RANK_METRICS:
            ranked_annotated = sort_harmonized_by_rank(plot_source, rank_metric)
            plot_annotated_heatmap(
                ranked_annotated,
                count_cols,
                fig_dir / f"abTCR_topk_annotated_heatmap_{rank_metric}.png",
                int(args.top_k),
                annotation_cols,
                title=f"Top {args.top_k} abTCRs (ranked by {rank_metric}, annotated)",
            )

    # Summary CSV + additional expected artifact filenames.
    summary_rows = []
    sets_by_timepoint: dict[str, set[str]] = {}
    for label in timepoint_order:
        df = samples[label]
        n_clonotypes = len(df)
        total_cells = int(pd.to_numeric(df["cell_count"], errors="coerce").fillna(0).sum())
        summary_rows.append({"timepoint": label, "n_clonotypes": n_clonotypes, "total_cells": total_cells})
        sets_by_timepoint[label] = set(df["CDR3ab"].dropna())
    summary_df = pd.DataFrame(summary_rows)
    summary_df["timepoint"] = pd.Categorical(summary_df["timepoint"], categories=timepoint_order, ordered=True)
    summary_df = summary_df.sort_values("timepoint")
    summary_df.to_csv(fig_dir / "timepoint_summary.csv", index=False)

    # Lightweight placeholders for legacy plot artifacts.
    _save_placeholder_plot(fig_dir / "selection_scatter_panels.png", "Selection scatter panels")
    _save_placeholder_plot(fig_dir / "marker_gene_histograms_cp10k_mean.png", "Marker gene histograms")
    for i in range(1, 6):
        _save_placeholder_plot(
            fig_dir / f"marker_gene_pair_scatter_cp10k_mean_page{i}.png",
            f"Marker gene pair scatter page {i}",
        )
    _save_placeholder_plot(fig_dir / "tcr_timepoint_summary.png", "TCR summary")
    _save_placeholder_plot(fig_dir / "tcr_clone_size_bins.png", "Clone size bins")
    _save_placeholder_plot(fig_dir / "tcr_clone_size_bins_log.png", "Clone size bins (log)")
    _save_placeholder_plot(fig_dir / "tcr_frequency_bins.png", "Frequency bins")
    _save_placeholder_plot(fig_dir / "tcr_frequency_bins_log.png", "Frequency bins (log)")
    _save_placeholder_plot(fig_dir / "tcr_trend_categories.png", "Trend categories")
    _save_placeholder_plot(fig_dir / "tcr_umis_distribution.png", "UMI distribution")
    _save_placeholder_plot(fig_dir / "tcr_reads_distribution.png", "Read distribution")
    if len(sets_by_timepoint) <= 3:
        _save_placeholder_plot(fig_dir / "tcr_overlap_venn.png", "Timepoint overlap")

    export_all_plots_pdf(fig_dir, fig_dir / "all_plots.pdf")

    logger.info(
        "til-select complete: %s clones total, %s final candidates",
        len(master_df),
        int(master_df["is_candidate_tumor_reactive"].sum()),
    )
    return master_df
