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
TCR annotation using public databases for TCRsift.

Matches TCRs against VDJdb, IEDB, and CEDAR to identify known specificities.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_dataframe,
    validate_file_exists,
)

logger = logging.getLogger(__name__)


# Known viral species patterns for flagging
VIRAL_SPECIES_PATTERNS = [
    "cmv", "cytomegalovirus",
    "ebv", "epstein-barr",
    "hiv", "human immunodeficiency",
    "flu", "influenza",
    "sars", "coronavirus",
    "herpes", "hsv",
    "hpv", "papilloma",
    "hepatitis", "hbv", "hcv",
    "dengue", "zika",
    "yellow fever",
]


def load_vdjdb(path: str | Path, verbose: bool = True) -> pd.DataFrame:
    """
    Load VDJdb database.

    Parameters
    ----------
    path : str or Path
        Path to VDJdb directory or file
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        VDJdb entries with standardized columns
    """
    path = Path(path)

    if path.is_dir():
        # Look for the main database file
        candidates = list(path.glob("vdjdb*.txt")) + list(path.glob("vdjdb*.tsv"))
        if not candidates:
            available = [f.name for f in path.iterdir()][:15]
            raise TCRsiftValidationError(
                f"No VDJdb files found in directory: {path}",
                hint=f"Expected files matching 'vdjdb*.txt' or 'vdjdb*.tsv'. "
                f"Available files: {available}",
            )
        db_file = candidates[0]
    else:
        db_file = validate_file_exists(path, "VDJdb database file")

    if verbose:
        logger.info(f"Loading VDJdb from {db_file}")

    try:
        df = pd.read_csv(db_file, sep="\t", low_memory=False)
    except Exception as e:
        raise TCRsiftValidationError(
            f"Failed to read VDJdb file: {db_file}",
            hint=f"Error: {e}. Make sure the file is a valid TSV file.",
        )

    if len(df) == 0:
        raise TCRsiftValidationError(
            f"VDJdb file is empty: {db_file}",
            hint="Download a fresh copy from https://vdjdb.cdr3.net/",
        )

    # Standardize columns
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

    df["database"] = "VDJdb"

    # Flag viral entries
    df["is_viral"] = _flag_viral(df)

    if verbose:
        logger.info(f"  Loaded {len(df):,} VDJdb entries ({df['is_viral'].sum():,} viral)")
    return df


def load_iedb(path: str | Path) -> pd.DataFrame:
    """
    Load IEDB TCR database.

    Parameters
    ----------
    path : str or Path
        Path to IEDB file

    Returns
    -------
    pd.DataFrame
        IEDB entries with standardized columns
    """
    path = Path(path)
    logger.info(f"Loading IEDB from {path}")

    df = pd.read_csv(path, sep="\t", low_memory=False)

    # Standardize columns (IEDB format varies)
    # Common IEDB column names
    column_mapping = {
        "Chain 2 CDR3 Curated": "cdr3_beta",
        "Chain 1 CDR3 Curated": "cdr3_alpha",
        "Epitope - Name": "epitope",
        "Epitope - Source Molecule Name": "antigen_gene",
        "Epitope - Source Organism Name": "species",
        "MHC Allele Names": "mhc_allele",
    }

    for old, new in column_mapping.items():
        if old in df.columns:
            df[new] = df[old]

    df["database"] = "IEDB"
    df["is_viral"] = _flag_viral(df)

    logger.info(f"Loaded {len(df)} IEDB entries ({df['is_viral'].sum()} viral)")
    return df


def load_cedar(path: str | Path) -> pd.DataFrame:
    """
    Load CEDAR TCR database.

    Parameters
    ----------
    path : str or Path
        Path to CEDAR file

    Returns
    -------
    pd.DataFrame
        CEDAR entries with standardized columns
    """
    path = Path(path)
    logger.info(f"Loading CEDAR from {path}")

    df = pd.read_csv(path, sep="\t", low_memory=False)

    # Standardize columns
    column_mapping = {
        "cdr3_b_aa": "cdr3_beta",
        "cdr3_a_aa": "cdr3_alpha",
        "epitope_sequence": "epitope",
        "antigen_name": "antigen_gene",
        "organism": "species",
    }

    for old, new in column_mapping.items():
        if old in df.columns:
            df[new] = df[old]

    df["database"] = "CEDAR"
    df["is_viral"] = _flag_viral(df)

    logger.info(f"Loaded {len(df)} CEDAR entries ({df['is_viral'].sum()} viral)")
    return df


def _flag_viral(df: pd.DataFrame) -> pd.Series:
    """Flag entries as viral based on species column."""
    if "species" not in df.columns:
        return pd.Series(False, index=df.index)

    species_lower = df["species"].fillna("").str.lower()

    is_viral = pd.Series(False, index=df.index)
    for pattern in VIRAL_SPECIES_PATTERNS:
        is_viral |= species_lower.str.contains(pattern, na=False)

    return is_viral


def load_databases(
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    Load and combine multiple TCR databases.

    Parameters
    ----------
    vdjdb_path : str or Path, optional
        Path to VDJdb
    iedb_path : str or Path, optional
        Path to IEDB
    cedar_path : str or Path, optional
        Path to CEDAR

    Returns
    -------
    pd.DataFrame
        Combined database with standardized columns
    """
    dfs = []

    if vdjdb_path:
        dfs.append(load_vdjdb(vdjdb_path))
    if iedb_path:
        dfs.append(load_iedb(iedb_path))
    if cedar_path:
        dfs.append(load_cedar(cedar_path))

    if not dfs:
        raise ValueError("At least one database path must be provided")

    # Combine and deduplicate
    combined = pd.concat(dfs, ignore_index=True)

    # Keep only rows with at least a beta CDR3
    combined = combined[combined["cdr3_beta"].notna() & (combined["cdr3_beta"] != "")]

    logger.info(f"Combined database has {len(combined)} entries")
    return combined


def match_clonotypes(
    clonotypes: pd.DataFrame,
    database: pd.DataFrame,
    match_by: str = "CDR3ab",
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Match clonotypes against public database.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    database : pd.DataFrame
        Combined database from load_databases
    match_by : str
        Matching strategy: "CDR3ab" (both chains) or "CDR3b_only" (beta only)
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with match annotations added
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_annotation=True)
    database = validate_dataframe(database, "database", min_rows=1)

    valid_match_by = ["CDR3ab", "CDR3b_only"]
    if match_by not in valid_match_by:
        raise TCRsiftValidationError(
            f"Invalid match_by: '{match_by}'",
            hint=f"Valid options are: {valid_match_by}",
        )

    if verbose:
        logger.info(f"Matching {len(clonotypes):,} clonotypes against {len(database):,} database entries by {match_by}")

    df = clonotypes.copy()

    # Initialize annotation columns
    df["db_match"] = False
    df["db_epitope"] = None
    df["db_species"] = None
    df["db_database"] = None
    df["is_viral"] = False

    # Build lookup sets for fast matching
    if match_by == "CDR3ab":
        # Match on both alpha and beta
        db_alpha_beta = set(
            zip(
                database["cdr3_alpha"].fillna(""),
                database["cdr3_beta"].fillna("")
            )
        )

        # Create iterator with optional progress bar
        row_iter = df.iterrows()
        if show_progress:
            row_iter = tqdm(
                list(df.iterrows()),
                desc="Matching clonotypes",
                unit="clone",
            )

        for idx, row in row_iter:
            alpha = row.get("CDR3_alpha", "") or ""
            beta = row.get("CDR3_beta", "") or ""

            if (alpha, beta) in db_alpha_beta:
                matches = database[
                    (database["cdr3_alpha"] == alpha) &
                    (database["cdr3_beta"] == beta)
                ]
                _annotate_match(df, idx, matches)

            # Also try beta-only match as fallback
            elif beta and beta in database["cdr3_beta"].values:
                matches = database[database["cdr3_beta"] == beta]
                _annotate_match(df, idx, matches, partial=True)

    else:  # CDR3b_only
        db_beta_set = set(database["cdr3_beta"].dropna())

        # Create iterator with optional progress bar
        row_iter = df.iterrows()
        if show_progress:
            row_iter = tqdm(
                list(df.iterrows()),
                desc="Matching clonotypes",
                unit="clone",
            )

        for idx, row in row_iter:
            beta = row.get("CDR3_beta", "") or ""
            if beta in db_beta_set:
                matches = database[database["cdr3_beta"] == beta]
                _annotate_match(df, idx, matches)

    n_matches = df["db_match"].sum()
    n_viral = df["is_viral"].sum()
    if verbose:
        logger.info(f"  Found {n_matches:,} matches ({n_viral:,} viral)")

    return df


def _annotate_match(
    df: pd.DataFrame,
    idx: int,
    matches: pd.DataFrame,
    partial: bool = False,
):
    """Annotate a single clonotype with match information."""
    if len(matches) == 0:
        return

    df.loc[idx, "db_match"] = True

    # Take most common epitope
    epitopes = matches["epitope"].dropna()
    if len(epitopes) > 0:
        df.loc[idx, "db_epitope"] = epitopes.mode().iloc[0]

    # Take most common species
    species = matches["species"].dropna()
    if len(species) > 0:
        df.loc[idx, "db_species"] = species.mode().iloc[0]

    # Record database sources
    df.loc[idx, "db_database"] = ";".join(matches["database"].unique())

    # Viral flag
    df.loc[idx, "is_viral"] = matches["is_viral"].any()

    # Partial match flag
    if partial:
        df.loc[idx, "db_match_partial"] = True


def annotate_clonotypes(
    clonotypes: pd.DataFrame,
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    match_by: str = "CDR3ab",
    exclude_viral: bool = False,
    flag_only: bool = False,
) -> pd.DataFrame:
    """
    Main annotation function.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    vdjdb_path, iedb_path, cedar_path : str or Path, optional
        Paths to databases
    match_by : str
        Matching strategy
    exclude_viral : bool
        Remove clones matching viral epitopes
    flag_only : bool
        Just flag viral, don't remove

    Returns
    -------
    pd.DataFrame
        Annotated clonotypes
    """
    # Load databases
    database = load_databases(
        vdjdb_path=vdjdb_path,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
    )

    # Match clonotypes
    df = match_clonotypes(clonotypes, database, match_by=match_by)

    # Handle viral exclusion
    if exclude_viral and not flag_only:
        initial = len(df)
        df = df[~df["is_viral"]]
        logger.info(f"Excluded {initial - len(df)} viral clones")

    return df


def get_annotation_summary(clonotypes: pd.DataFrame) -> dict:
    """
    Get summary of annotation results.

    Returns
    -------
    dict
        Summary statistics
    """
    summary = {
        "total": len(clonotypes),
        "matched": clonotypes["db_match"].sum() if "db_match" in clonotypes.columns else 0,
        "viral": clonotypes["is_viral"].sum() if "is_viral" in clonotypes.columns else 0,
    }

    if "db_database" in clonotypes.columns:
        db_counts = {}
        for db in ["VDJdb", "IEDB", "CEDAR"]:
            db_counts[db] = clonotypes["db_database"].fillna("").str.contains(db).sum()
        summary["database_breakdown"] = db_counts

    if "db_species" in clonotypes.columns:
        species_counts = clonotypes["db_species"].value_counts().head(10).to_dict()
        summary["top_species"] = species_counts

    return summary
