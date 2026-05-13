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
    "cmv",
    "cytomegalovirus",
    "ebv",
    "epstein-barr",
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

# Tumor-associated self antigens. Matched as case-insensitive substrings
# against the antigen_gene/Source Molecule field. Independent of species
# (these are all Homo sapiens) — the override exists so that e.g. MART-1
# doesn't end up bucketed as plain "self".
TUMOR_SELF_ANTIGEN_PATTERNS = [
    "mart-1",
    "mart1",
    "mlana",
    "melan-a",
    "ny-eso-1",
    "ny eso",
    "ctag1b",
    "cancer/testis",
    "mage-a",
    "magea",
    "gp100",
    "pmel",
    "tyrosinase",
    "ceacam5",
    "carcinoembryonic",
    "wt1",
    "wilms",
    "telomerase",
    "tert",
    "ny-br-1",
    "nybr1",
    "her2",
    "erbb2",
    "survivin",
    "birc5",
    "psa ",
    "prostate-specific antigen",
]

# Bacterial species patterns — used to derive db_category="bacterial".
# Kept separate from VIRAL_SPECIES_PATTERNS so the existing is_viral
# heuristic stays unchanged.
BACTERIAL_SPECIES_PATTERNS = [
    "mycobacterium",
    "tuberculosis",
    "staphylococcus",
    "streptococcus",
    "listeria",
    "salmonella",
    "escherichia coli",
    "borrelia",
    "bacillus",
    "clostridium",
    "neisseria",
    "helicobacter",
    "chlamydia",
    "treponema",
]

CATEGORY_VIRAL = "viral"
CATEGORY_BACTERIAL = "bacterial"
CATEGORY_SELF = "self"
CATEGORY_TUMOR_SELF = "tumor_self"
CATEGORY_OTHER = "other"
CATEGORY_UNKNOWN = "unknown"


# Valid strictness modes for match_clonotypes. Translated to the internal
# (do_ab_match, do_beta_fallback) tuple inside the function.
MATCH_STRICTNESS_MODES = ("strict_ab", "ab_with_partial", "b_only")


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

    Handles two on-disk shapes:

    1. The current `receptor_full_v3.zip` export (cached as
       ``tcr_full_v3.csv``): comma-separated with a two-row hierarchical
       header — top row is the section (``Receptor``/``Epitope``/
       ``Assay``/``Chain 1``/``Chain 2``…), second row is the field
       (e.g. ``CDR3 Curated``, ``Source Organism``).
    2. Older flat-TSV exports kept around for compatibility.

    The format is sniffed from the first byte of the file rather than
    the extension, so a user-supplied path is treated correctly
    regardless of suffix.

    Parameters
    ----------
    path : str or Path
        Path to IEDB file (v3 CSV or legacy flat TSV).

    Returns
    -------
    pd.DataFrame
        IEDB entries with standardized columns (``cdr3_alpha``,
        ``cdr3_beta``, ``epitope``, ``antigen_gene``, ``species``,
        ``mhc_allele``, ``is_viral``, ``database``).
    """
    path = Path(path)
    logger.info(f"Loading IEDB from {path}")

    if _looks_like_iedb_v3(path):
        df = _load_iedb_v3(path)
    else:
        df = _load_iedb_legacy_tsv(path)

    df["database"] = "IEDB"
    df["is_viral"] = _flag_viral(df)

    logger.info(f"Loaded {len(df)} IEDB entries ({df['is_viral'].sum()} viral)")
    return df


def _looks_like_iedb_v3(path: Path) -> bool:
    """Detect the v3 hierarchical-header CSV by sniffing the first line.

    The v3 export begins ``Receptor,Receptor,Receptor,...`` (each
    section name repeated once per column under it). Anything else
    (notably the older flat TSV) falls through to the legacy loader.
    """
    try:
        with open(path, encoding="utf-8") as f:
            first = f.readline()
    except OSError:
        return False
    return first.startswith("Receptor,Receptor")


def _load_iedb_v3(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, header=[0, 1], low_memory=False)

    # IEDB v3 puts the alpha chain in "Chain 1" and the beta chain in
    # "Chain 2" for alphabeta receptors. Drop gammadelta/construct
    # rows — they aren't matchable against αβ clonotypes.
    receptor_type = raw[("Receptor", "Type")] if ("Receptor", "Type") in raw.columns else None
    if receptor_type is not None:
        raw = raw[receptor_type == "alphabeta"].copy()

    def _pick(col_section: str, *field_options: str) -> pd.Series:
        """Pull the first present (section, field) combo, else NaN."""
        for field in field_options:
            key = (col_section, field)
            if key in raw.columns:
                return raw[key]
        return pd.Series([pd.NA] * len(raw), index=raw.index)

    # Prefer "CDR3 Curated" (manually validated) over "CDR3 Calculated"
    # (algorithmic), matching how VDJdb prioritizes curated entries.
    cdr3_alpha = _pick("Chain 1", "CDR3 Curated").fillna(
        _pick("Chain 1", "CDR3 Calculated")
    )
    cdr3_beta = _pick("Chain 2", "CDR3 Curated").fillna(
        _pick("Chain 2", "CDR3 Calculated")
    )

    out = pd.DataFrame(
        {
            "cdr3_alpha": cdr3_alpha,
            "cdr3_beta": cdr3_beta,
            "epitope": _pick("Epitope", "Name"),
            "antigen_gene": _pick("Epitope", "Source Molecule"),
            "species": _pick("Epitope", "Source Organism"),
            "mhc_allele": _pick("Assay", "MHC Allele Names"),
        }
    )
    # Drop rows with no CDR3 at all — they can't participate in matching.
    out = out[out["cdr3_alpha"].notna() | out["cdr3_beta"].notna()].reset_index(drop=True)
    return out


def _load_iedb_legacy_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)

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


def classify_category(species: pd.Series, antigen_gene: pd.Series) -> pd.Series:
    """Classify each entry into a coarse specificity category.

    Categories: ``viral``, ``bacterial``, ``tumor_self``, ``self``,
    ``other``, ``unknown``. Tumor-associated self antigens override the
    species-derived category so MART-1 / NY-ESO-1 / MAGE etc. don't
    appear as plain ``self``.

    Parameters
    ----------
    species : pd.Series
        Source organism per entry (VDJdb's ``antigen.species`` /
        IEDB's ``Epitope - Source Organism``).
    antigen_gene : pd.Series
        Source protein per entry. Used only for the tumor-self override.

    Returns
    -------
    pd.Series
        Category string per entry, same index as inputs.
    """
    species_lower = species.fillna("").astype(str).str.lower()
    antigen_lower = antigen_gene.fillna("").astype(str).str.lower()

    category = pd.Series(
        [CATEGORY_UNKNOWN] * len(species),
        index=species.index,
        dtype=object,
    )

    is_known_species = species_lower.str.len() > 0
    category[is_known_species] = CATEGORY_OTHER

    is_viral = pd.Series(False, index=species.index)
    for pattern in VIRAL_SPECIES_PATTERNS:
        is_viral |= species_lower.str.contains(pattern, na=False)
    category[is_viral] = CATEGORY_VIRAL

    is_bacterial = pd.Series(False, index=species.index)
    for pattern in BACTERIAL_SPECIES_PATTERNS:
        is_bacterial |= species_lower.str.contains(pattern, na=False)
    category[is_bacterial & ~is_viral] = CATEGORY_BACTERIAL

    # "Homo sapiens", "Homo sapiens (human)", "human" all bucket as self.
    is_self = species_lower.str.contains("homo sapiens", na=False) | (
        species_lower == "human"
    )
    category[is_self & ~is_viral & ~is_bacterial] = CATEGORY_SELF

    # Tumor-self override: applies independent of species. A handful of
    # IEDB rows store tumor peptides under a viral alias; here we trust
    # the antigen name over the species field.
    is_tumor_self = pd.Series(False, index=species.index)
    for pattern in TUMOR_SELF_ANTIGEN_PATTERNS:
        is_tumor_self |= antigen_lower.str.contains(pattern, na=False)
    category[is_tumor_self] = CATEGORY_TUMOR_SELF

    return category


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
    match_strictness: str | None = None,
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
        Legacy matching strategy. ``"CDR3ab"`` (both chains, β-only
        fallback) or ``"CDR3b_only"`` (beta only). Kept for back-compat;
        prefer ``match_strictness`` for new code.
    match_strictness : str, optional
        Explicit matching strictness. When set, takes precedence over
        ``match_by``. One of:

        - ``"strict_ab"`` — αβ-only, no β-only fallback. Best when you
          want headline match counts you can trust.
        - ``"ab_with_partial"`` — αβ first, β-only fallback per clone
          (equivalent to legacy ``match_by="CDR3ab"``).
        - ``"b_only"`` — β-only across the board (equivalent to legacy
          ``match_by="CDR3b_only"``).
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with match annotations: ``db_match``,
        ``db_match_partial`` (β-only fallback flag — back-compat),
        ``db_match_strength`` (``"ab"`` / ``"b_only"`` / None),
        ``db_epitope``, ``db_protein``, ``db_species``, ``db_mhc``,
        ``db_category``, ``db_database``, ``is_viral``.
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_annotation=True)
    database = validate_dataframe(database, "database", min_rows=1)

    if match_strictness is not None:
        if match_strictness not in MATCH_STRICTNESS_MODES:
            raise TCRsiftValidationError(
                f"Invalid match_strictness: '{match_strictness}'",
                hint=f"Valid options are: {list(MATCH_STRICTNESS_MODES)}",
            )
        strictness = match_strictness
    else:
        valid_match_by = ["CDR3ab", "CDR3b_only"]
        if match_by not in valid_match_by:
            raise TCRsiftValidationError(
                f"Invalid match_by: '{match_by}'",
                hint=f"Valid options are: {valid_match_by}",
            )
        strictness = "ab_with_partial" if match_by == "CDR3ab" else "b_only"

    if verbose:
        logger.info(
            f"Matching {len(clonotypes):,} clonotypes against "
            f"{len(database):,} database entries (strictness={strictness})"
        )

    df = clonotypes.copy()

    # Initialize annotation columns. New columns added in #48; existing
    # ones (db_match, db_match_partial, db_epitope, db_species,
    # db_database, is_viral) kept for back-compat.
    df["db_match"] = False
    df["db_match_partial"] = False
    df["db_match_strength"] = None
    df["db_epitope"] = None
    df["db_protein"] = None
    df["db_species"] = None
    df["db_mhc"] = None
    df["db_category"] = None
    df["db_database"] = None
    df["is_viral"] = False

    # Defensive check: the upstream `load_vdjdb` / `load_iedb` / `load_cedar`
    # loaders standardize column names via mappings keyed on header strings
    # they expect. When a database file's format drifts (e.g. IEDB's v3 CSV
    # export with hierarchical headers vs the older flat-header TSV the
    # loader was written for, #46), the standardization silently does
    # nothing and the resulting `database` lacks `cdr3_alpha`/`cdr3_beta`.
    # Catch that here with a clear warning rather than `KeyError` later.
    required = {"cdr3_beta"} if strictness == "b_only" else {"cdr3_alpha", "cdr3_beta"}
    missing = required - set(database.columns)
    if missing:
        dbs = sorted(set(database.get("database", pd.Series(["unknown"]))))
        logger.warning(
            f"Skipping annotation: database(s) {dbs} are missing required "
            f"columns {sorted(missing)} after standardization. The source "
            "file format may have changed since the loader was written. "
            "Inspect the column list before reporting this as a tcrsift bug."
        )
        return df

    # Build lookup sets for fast matching
    if strictness in ("strict_ab", "ab_with_partial"):
        allow_b_fallback = strictness == "ab_with_partial"
        db_alpha_beta = set(
            zip(database["cdr3_alpha"].fillna(""), database["cdr3_beta"].fillna(""))
        )
        db_beta_values = (
            set(database["cdr3_beta"].dropna()) if allow_b_fallback else set()
        )

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
                    (database["cdr3_alpha"] == alpha) & (database["cdr3_beta"] == beta)
                ]
                _annotate_match(df, idx, matches, strength="ab")
            elif allow_b_fallback and beta and beta in db_beta_values:
                matches = database[database["cdr3_beta"] == beta]
                _annotate_match(df, idx, matches, strength="b_only", partial=True)

    else:  # strictness == "b_only"
        db_beta_set = set(database["cdr3_beta"].dropna())

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
                _annotate_match(df, idx, matches, strength="b_only")

    n_matches = df["db_match"].sum()
    n_viral = df["is_viral"].sum()
    if verbose:
        logger.info(f"  Found {n_matches:,} matches ({n_viral:,} viral)")

    return df


def _annotate_match(
    df: pd.DataFrame,
    idx: int,
    matches: pd.DataFrame,
    strength: str = "ab",
    partial: bool = False,
):
    """Annotate a single clonotype with match information.

    For multi-row matches, picks the most common value per field
    (``.mode().iloc[0]``). Category is derived from the picked
    species/protein so it stays internally consistent with what's
    surfaced in ``db_species`` / ``db_protein``.
    """
    if len(matches) == 0:
        return

    df.loc[idx, "db_match"] = True
    df.loc[idx, "db_match_strength"] = strength

    epitopes = matches["epitope"].dropna()
    if len(epitopes) > 0:
        df.loc[idx, "db_epitope"] = epitopes.mode().iloc[0]

    species_val = None
    species = matches["species"].dropna() if "species" in matches.columns else pd.Series(dtype=object)
    if len(species) > 0:
        species_val = species.mode().iloc[0]
        df.loc[idx, "db_species"] = species_val

    protein_val = None
    if "antigen_gene" in matches.columns:
        proteins = matches["antigen_gene"].dropna()
        if len(proteins) > 0:
            protein_val = proteins.mode().iloc[0]
            df.loc[idx, "db_protein"] = protein_val

    if "mhc_allele" in matches.columns:
        mhcs = matches["mhc_allele"].dropna()
        if len(mhcs) > 0:
            df.loc[idx, "db_mhc"] = mhcs.mode().iloc[0]

    if species_val is not None or protein_val is not None:
        category = classify_category(
            pd.Series([species_val or ""]),
            pd.Series([protein_val or ""]),
        ).iloc[0]
        df.loc[idx, "db_category"] = category

    df.loc[idx, "db_database"] = ";".join(matches["database"].unique())
    df.loc[idx, "is_viral"] = matches["is_viral"].any()

    if partial:
        df.loc[idx, "db_match_partial"] = True


def annotate_clonotypes(
    clonotypes: pd.DataFrame,
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    match_by: str = "CDR3ab",
    match_strictness: str | None = None,
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
        Legacy matching strategy. Prefer ``match_strictness`` for new code.
    match_strictness : str, optional
        Explicit matching strictness — ``"strict_ab"`` /
        ``"ab_with_partial"`` / ``"b_only"``. Overrides ``match_by``
        when set. See :func:`match_clonotypes`.
    exclude_viral : bool
        Remove clones matching viral epitopes
    flag_only : bool
        Just flag viral, don't remove

    Returns
    -------
    pd.DataFrame
        Annotated clonotypes
    """
    # Columns that get initialized on the no-database short-circuit.
    # Mirrors the column set produced by match_clonotypes so downstream
    # code can rely on these always being present.
    _DEFAULT_ANNOTATION_COLUMNS = {
        "db_match": False,
        "db_match_partial": False,
        "db_match_strength": None,
        "db_epitope": None,
        "db_protein": None,
        "db_species": None,
        "db_mhc": None,
        "db_category": None,
        "db_database": None,
        "is_viral": False,
    }

    # Annotation is optional: if no databases are provided, return input with default annotation columns.
    if not any([vdjdb_path, iedb_path, cedar_path]):
        logger.info("No annotation database paths provided; returning input with empty annotations")
        df = clonotypes.copy()
        for col, default in _DEFAULT_ANNOTATION_COLUMNS.items():
            if col not in df.columns:
                df[col] = default
        return df

    # Load databases
    database = load_databases(
        vdjdb_path=vdjdb_path,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
    )

    # Match clonotypes
    df = match_clonotypes(
        clonotypes,
        database,
        match_by=match_by,
        match_strictness=match_strictness,
    )

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
