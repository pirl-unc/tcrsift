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
SCT (single-cell TCR) platform data loading and processing.

Functions for loading TCR data from PACT-style single-cell platforms,
which provide paired TCR sequencing with antigen specificity information
via pMHC tetramer staining.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .validation import TCRsiftValidationError, validate_file_exists

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def load_sct(
    path: str | Path,
    *,
    sheet_name: str = "Cell",
    min_snr: float = 2.0,
    min_reads_per_chain: int = 10,
    require_mutation_match: bool = True,
    require_compact_match: bool = False,
    standardize_columns: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load TCR data from SCT platform Excel file.

    The SCT (single-cell TCR) platform provides paired TCR sequencing with
    antigen specificity information from pMHC tetramer staining.

    Parameters
    ----------
    path : str or Path
        Path to SCT Excel file
    sheet_name : str
        Sheet name to read (default: "Cell")
    min_snr : float
        Minimum signal-to-noise ratio (default: 2.0)
    min_reads_per_chain : int
        Minimum reads per chain (default: 10)
    require_mutation_match : bool
        Require PE and APC mutation calls to match (default: True)
    require_compact_match : bool
        Require PE and APC comPACT ID calls to match (default: False)
    standardize_columns : bool
        Add standardized CDR3_alpha, CDR3_beta, CDR3_pair columns
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        SCT data with quality flags and standardized columns
    """
    path = validate_file_exists(Path(path), "SCT file")

    if verbose:
        logger.info(f"Loading SCT data from {path}")

    # Load Excel file
    try:
        df = pd.read_excel(path, sheet_name=sheet_name)
    except Exception as e:
        raise TCRsiftValidationError(
            f"Failed to read SCT Excel file: {e}",
            hint="Ensure the file is a valid Excel file with the expected sheet name",
        )

    if verbose:
        logger.info(f"  Loaded {len(df):,} rows from sheet '{sheet_name}'")

    # Detect column naming convention
    tra_cdr3_col = _find_column(df, ["tra.CDR3", "TRA_CDR3", "CDR3_alpha", "cdr3_alpha"])
    trb_cdr3_col = _find_column(df, ["trb.CDR3", "TRB_CDR3", "CDR3_beta", "cdr3_beta"])
    tra_reads_col = _find_column(df, ["tra.readcount", "TRA_reads", "alpha_reads"])
    trb_reads_col = _find_column(df, ["trb.readcount", "TRB_reads", "beta_reads"])
    snr_col = _find_column(df, ["SNR", "snr", "signal_noise_ratio"])

    if not tra_cdr3_col or not trb_cdr3_col:
        raise TCRsiftValidationError(
            "Could not find CDR3 columns in SCT file",
            hint="Expected columns like 'tra.CDR3'/'trb.CDR3' or 'CDR3_alpha'/'CDR3_beta'",
        )

    # Add quality flags
    quality_mask = pd.Series(True, index=df.index)

    # SNR filter
    if snr_col and min_snr > 0:
        snr_mask = df[snr_col].fillna(0) >= min_snr
        quality_mask &= snr_mask
        if verbose:
            logger.info(f"  SNR >= {min_snr}: {snr_mask.sum():,} pass")

    # Read count filters
    if tra_reads_col and min_reads_per_chain > 0:
        tra_mask = df[tra_reads_col].fillna(0) >= min_reads_per_chain
        quality_mask &= tra_mask
        if verbose:
            logger.info(f"  Alpha reads >= {min_reads_per_chain}: {tra_mask.sum():,} pass")

    if trb_reads_col and min_reads_per_chain > 0:
        trb_mask = df[trb_reads_col].fillna(0) >= min_reads_per_chain
        quality_mask &= trb_mask
        if verbose:
            logger.info(f"  Beta reads >= {min_reads_per_chain}: {trb_mask.sum():,} pass")

    # CDR3 length filters
    tra_len_mask = df[tra_cdr3_col].fillna("").str.len() > 1
    trb_len_mask = df[trb_cdr3_col].fillna("").str.len() > 1
    quality_mask &= tra_len_mask & trb_len_mask

    # Mutation match filter
    mutation_pe_col = _find_column(df, ["Top1.mutation.PE.eq.Top1.mutation.APC", "mutation_match"])
    if require_mutation_match and mutation_pe_col:
        mutation_mask = df[mutation_pe_col] == "Yes"
        quality_mask &= mutation_mask
        if verbose:
            logger.info(f"  Mutation match: {mutation_mask.sum():,} pass")

    # comPACT match filter
    compact_col = _find_column(df, ["Top1.comPACT.ID.PE.eq.Top1.comPACT.ID.APC", "compact_match"])
    if require_compact_match and compact_col:
        compact_mask = df[compact_col] == "Yes"
        quality_mask &= compact_mask
        if verbose:
            logger.info(f"  comPACT match: {compact_mask.sum():,} pass")

    df["high_quality"] = quality_mask

    # Stricter "chosen" criteria (SNR >= 3.4, reads >= 50)
    chosen_mask = quality_mask.copy()
    if snr_col:
        chosen_mask &= df[snr_col].fillna(0) >= 3.4
    if tra_reads_col:
        chosen_mask &= df[tra_reads_col].fillna(0) >= 50
    if trb_reads_col:
        chosen_mask &= df[trb_reads_col].fillna(0) >= 50
    if compact_col:
        chosen_mask &= df[compact_col] == "Yes"
    df["chosen"] = chosen_mask

    if verbose:
        logger.info(f"  High quality: {df['high_quality'].sum():,}")
        logger.info(f"  Chosen (strict): {df['chosen'].sum():,}")

    # Standardize columns
    if standardize_columns:
        df["CDR3_alpha"] = df[tra_cdr3_col].fillna("")
        df["CDR3_beta"] = df[trb_cdr3_col].fillna("")
        df["CDR3_pair"] = df["CDR3_alpha"] + "/" + df["CDR3_beta"]
        df["complete"] = (df["CDR3_alpha"].str.len() > 0) & (df["CDR3_beta"].str.len() > 0)

        # Move standardized columns to front
        cols = ["CDR3_pair", "CDR3_alpha", "CDR3_beta", "complete", "high_quality", "chosen"]
        other_cols = [c for c in df.columns if c not in cols]
        df = df[cols + other_cols]

    if verbose:
        logger.info(f"  Complete pairs: {df['complete'].sum():,}")

    return df


def aggregate_sct(
    df: pd.DataFrame,
    *,
    group_cols: list[str] | None = None,
    numeric_cols: list[str] | None = None,
    boolean_cols: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Aggregate SCT data by CDR3 pair, computing statistics.

    Parameters
    ----------
    df : pd.DataFrame
        SCT data (from load_sct)
    group_cols : list of str, optional
        Columns to group by (default: CDR3_pair, CDR3_alpha, CDR3_beta)
    numeric_cols : list of str, optional
        Numeric columns to aggregate with min/median/max
    boolean_cols : list of str, optional
        Boolean columns to aggregate with any/all
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        Aggregated clonotype data with statistics
    """
    if group_cols is None:
        group_cols = ["CDR3_pair", "CDR3_alpha", "CDR3_beta"]

    if numeric_cols is None:
        numeric_cols = [
            c for c in df.columns
            if c not in group_cols
            and df[c].dtype in ["int64", "float64"]
            and not c.endswith("_bool")
        ]
        # Filter to likely SCT columns
        numeric_cols = [
            c for c in numeric_cols
            if any(x in c.lower() for x in ["snr", "readcount", "count", "ratio"])
        ]

    if boolean_cols is None:
        boolean_cols = ["high_quality", "chosen", "complete"]
        boolean_cols = [c for c in boolean_cols if c in df.columns]

    if verbose:
        logger.info(f"Aggregating {len(df):,} rows by CDR3 pair")
        logger.info(f"  Numeric columns: {numeric_cols}")
        logger.info(f"  Boolean columns: {boolean_cols}")

    # Build aggregation dictionary
    agg_dict = {"CDR3_pair": "count"}  # Count entries

    for col in numeric_cols:
        if col in df.columns:
            agg_dict[col] = ["min", "median", "max"]

    for col in boolean_cols:
        if col in df.columns:
            # Convert Yes/No to boolean if needed
            if df[col].dtype == "object":
                df[f"{col}_bool"] = df[col].map({"Yes": True, "No": False, True: True, False: False})
            else:
                df[f"{col}_bool"] = df[col].astype(bool)
            agg_dict[f"{col}_bool"] = ["any", "all"]

    # Aggregate
    result = df.groupby(group_cols).agg(agg_dict)

    # Flatten column names
    result.columns = [
        f"{col[0]}_{col[1]}" if col[1] else col[0]
        for col in result.columns
    ]

    # Rename columns
    result = result.rename(columns={"CDR3_pair_count": "num_original_entries"})
    for col in numeric_cols:
        if col in df.columns:
            result = result.rename(columns={
                f"{col}_min": f"{col}.min",
                f"{col}_median": f"{col}.median",
                f"{col}_max": f"{col}.max",
            })
    for col in boolean_cols:
        if col in df.columns:
            result = result.rename(columns={
                f"{col}_bool_any": f"{col}.any",
                f"{col}_bool_all": f"{col}.all",
            })

    result = result.reset_index()

    # Extract mutation/specificity if consistent
    mutation_col = _find_column(df, ["Top1.mutation.PE", "mutation"])
    if mutation_col:
        mutation_map = {}
        for pair, sub_df in df.groupby("CDR3_pair"):
            mutations = sub_df[mutation_col].dropna().unique()
            if len(mutations) == 1:
                mutation_map[pair] = mutations[0]
        result["mutation"] = result["CDR3_pair"].map(mutation_map)

    if verbose:
        logger.info(f"  Aggregated to {len(result):,} unique CDR3 pairs")
        if "mutation" in result.columns:
            n_with_mutation = result["mutation"].notna().sum()
            logger.info(f"  With consistent mutation: {n_with_mutation:,}")

    # Clean up temp columns
    for col in boolean_cols:
        if f"{col}_bool" in df.columns:
            del df[f"{col}_bool"]

    return result


def _find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Find first matching column from candidates."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def get_sct_specificities(df: pd.DataFrame) -> dict[str, str]:
    """
    Extract CDR3 pair to specificity mapping from SCT data.

    Parameters
    ----------
    df : pd.DataFrame
        SCT data (from load_sct or aggregate_sct)

    Returns
    -------
    dict
        Mapping of CDR3_pair to mutation/specificity
    """
    mutation_col = _find_column(df, ["mutation", "Top1.mutation.PE", "specificity"])
    if not mutation_col:
        return {}

    result = {}
    for _, row in df.iterrows():
        pair = row.get("CDR3_pair")
        mutation = row.get(mutation_col)
        if pair and pd.notna(mutation):
            result[pair] = mutation

    return result
