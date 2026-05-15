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

from .pgen import annotate_publicness as _annotate_publicness
from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_dataframe,
    validate_file_exists,
)

logger = logging.getLogger(__name__)


def _apply_publicness(
    df: pd.DataFrame,
    *,
    cdr3_col: str,
    v_gene_col: str,
    j_gene_col: str,
    log_summary: bool,
) -> pd.DataFrame:
    """Helper used by ``annotate_clonotypes`` to add publicness on
    both the no-DB short-circuit path and the post-match path."""
    if cdr3_col not in df.columns:
        logger.warning(
            f"add_publicness=True but {cdr3_col!r} not in DataFrame; "
            "skipping publicness annotation"
        )
        return df
    df = _annotate_publicness(
        df, cdr3_col=cdr3_col, v_gene_col=v_gene_col, j_gene_col=j_gene_col
    )
    if log_summary:
        n_public = int((df["publicness"] >= 0.5).sum())
        logger.info(
            f"  Publicness: {n_public:,}/{len(df):,} clones above 0.5 "
            "(consider discounting their DB matches)"
        )
    return df


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

# Tumor-associated self antigens. Patterns are case-insensitive regex
# fragments matched against antigen_gene / Source Molecule, anchored at
# word boundaries (``\b``) so short tokens like ``wt1`` / ``her2`` /
# ``tert`` / ``psa`` only fire on standalone occurrences. Without that,
# 4-letter substrings would false-positive on unrelated gene names that
# happen to contain the fragment (``interterritorial`` → ``tert``, etc.).
# Multi-word patterns embed an explicit ``[- _]`` so they match
# ``ny-eso-1`` / ``ny eso 1`` / ``ny_eso_1`` uniformly.
TUMOR_SELF_ANTIGEN_PATTERNS = [
    r"\bmart-?1\b",
    r"\bmlana\b",
    r"\bmelan-?a\b",
    r"\bny[- _]eso[- _]?1?\b",
    r"\bctag1b\b",
    r"\bcancer[/-]?testis\b",
    r"\bmage-?a\d*\b",
    r"\bgp100\b",
    r"\bpmel\b",
    r"\btyrosinase\b",
    r"\bceacam5\b",
    r"\bcarcinoembryonic\b",
    r"\bwt1\b",
    r"\bwilms\b",
    r"\btelomerase\b",
    r"\bh?tert\b",
    r"\bny-?br-?1\b",
    r"\bher2\b",
    r"\berbb2\b",
    r"\bsurvivin\b",
    r"\bbirc5\b",
    r"\bpsa\b",
    r"\bprostate[- ]specific antigen\b",
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
# Set on a single clone when its DB matches disagree on category
# (e.g. one match says ``viral``, another says ``tumor_self``). Rather
# than pick one mode arbitrarily we expose the disagreement (#83).
CATEGORY_CONTRADICTORY = "contradictory"


# Valid strictness modes for match_clonotypes. The function dispatches
# on this string directly — there's no further internal translation.
MATCH_STRICTNESS_MODES = ("strict_ab", "ab_with_partial", "b_only")


# Ordered list of (pattern, canonical) pairs for collapsing the many
# free-text Source-Molecule strings IEDB/VDJdb emit into the short
# symbols people read papers in. Order matters: first matching pattern
# wins (#55). Patterns are case-insensitive substrings, matched against
# the raw antigen_gene field — this is more forgiving than word
# boundaries given how messy these strings are (``Spike glycoprotein
# [Severe acute respiratory syndrome coronavirus 2]`` etc).
#
# Specificity rule of thumb: list LONGER and more specific patterns
# BEFORE shorter, more ambiguous ones. ``Spike glycoprotein`` should be
# checked before ``Spike`` alone, ``Cancer/testis antigen 1`` before
# ``Cancer/testis antigen``. Where short tokens collide across species
# (``pp65`` is CMV; ``M1`` is Flu), the organism prefix is baked into
# the canonical form (``CMV pp65``, ``Flu M1``) so downstream plots
# stay unambiguous.
#
# This is intentionally a static seed list. The right long-term answer
# is an HGNC/UniProt join, but until that lands a curated table covers
# the antigens that actually show up in TIL / vaccine cohorts.
CANONICAL_ANTIGEN_ALIASES: tuple[tuple[str, str], ...] = (
    # ─── Tumor-associated self antigens (TIL cohorts) ───
    ("melanoma antigen recognized by t-cells", "MART-1"),
    ("melan-a", "MART-1"),
    ("mlana", "MART-1"),
    ("mart-1", "MART-1"),
    ("mart1", "MART-1"),
    ("cancer/testis antigen 1", "NY-ESO-1"),
    ("ctag1b", "NY-ESO-1"),
    ("ny-eso-1", "NY-ESO-1"),
    ("ny eso 1", "NY-ESO-1"),
    ("ny_eso_1", "NY-ESO-1"),
    ("mage-a3", "MAGE-A3"),
    ("magea3", "MAGE-A3"),
    ("mage-a4", "MAGE-A4"),
    ("magea4", "MAGE-A4"),
    ("mage-a1", "MAGE-A1"),
    ("magea1", "MAGE-A1"),
    ("mage-a", "MAGE-A"),
    ("melanocyte protein pmel", "gp100"),
    ("melanocyte-specific secreted glycoprotein", "gp100"),
    ("gp100", "gp100"),
    ("pmel", "gp100"),
    ("tyrosinase", "Tyrosinase"),
    ("ceacam5", "CEACAM5"),
    ("carcinoembryonic", "CEA"),
    ("wilms tumor", "WT1"),
    ("wt1", "WT1"),
    ("telomerase", "TERT"),
    ("htert", "TERT"),
    ("erbb2", "HER2"),
    ("her2", "HER2"),
    ("survivin", "Survivin"),
    ("birc5", "Survivin"),
    ("prostate-specific antigen", "PSA"),
    ("ny-br-1", "NY-BR-1"),
    ("bone marrow stromal antigen 2", "BST2"),
    ("tetherin", "BST2"),
    ("bst2", "BST2"),
    # ─── HLA class I / MHC ───
    ("hla class i histocompatibility antigen", "MHC class I"),
    ("mhc class i protein", "MHC class I"),
    ("mhc class i", "MHC class I"),
    ("hla class ii", "MHC class II"),
    # ─── Viral — CMV (HHV-5) ───
    ("65 kda phosphoprotein", "CMV pp65"),
    ("phosphoprotein 65", "CMV pp65"),
    ("cmv pp65", "CMV pp65"),
    ("pp65", "CMV pp65"),
    ("55 kda immediate-early", "CMV IE1"),
    ("regulatory protein ie1", "CMV IE1"),
    ("immediate-early protein 1", "CMV IE1"),
    ("immediate early 1", "CMV IE1"),
    ("cmv ie1", "CMV IE1"),
    # ─── Viral — EBV (HHV-4) ───
    ("nuclear antigen ebna-3", "EBV EBNA-3"),
    ("ebna-3", "EBV EBNA-3"),
    ("ebna3", "EBV EBNA-3"),
    ("ebna-1", "EBV EBNA-1"),
    ("nuclear antigen ebna-1", "EBV EBNA-1"),
    ("bzlf1", "EBV BZLF1"),
    ("bmlf1", "EBV BMLF1"),
    ("brlf1", "EBV BRLF1"),
    ("lmp1", "EBV LMP1"),
    ("lmp2", "EBV LMP2"),
    # ─── Viral — Influenza ───
    ("matrix protein 1", "Flu M1"),
    ("matrix protein m1", "Flu M1"),
    ("m1 protein", "Flu M1"),
    ("influenza a virus matrix", "Flu M1"),
    ("nucleoprotein", "Flu NP"),
    ("hemagglutinin", "Flu HA"),
    # ─── Viral — SARS-CoV-2 ───
    ("spike glycoprotein", "SARS-CoV-2 Spike"),
    ("surface glycoprotein", "SARS-CoV-2 Spike"),
    ("spike protein", "SARS-CoV-2 Spike"),
    ("sars-cov-2 spike", "SARS-CoV-2 Spike"),
    ("replicase polyprotein 1ab", "SARS-CoV-2 ORF1ab"),
    ("replicase polyprotein 1a", "SARS-CoV-2 ORF1ab"),
    ("orf1ab polyprotein", "SARS-CoV-2 ORF1ab"),
    ("orf1ab", "SARS-CoV-2 ORF1ab"),
    ("nucleocapsid", "SARS-CoV-2 N"),
    # ─── Viral — HIV ───
    ("gag polyprotein", "HIV Gag"),
    ("nef protein", "HIV Nef"),
    # ─── Viral — HTLV ───
    ("protein tax-1", "HTLV-1 Tax"),
    ("transcriptional activator tax", "HTLV-1 Tax"),
    ("tax-1", "HTLV-1 Tax"),
    # ─── Viral — HCV ───
    ("ns3", "HCV NS3"),
    ("ns5b", "HCV NS5B"),
    # ─── Bacterial ───
    ("esat-6", "ESAT-6"),
    ("listeriolysin", "LLO"),
    # ─── Common self ───
    ("insulin", "Insulin"),
    ("beta-2-microglobulin", "B2M"),
    ("b2m", "B2M"),
)


def canonicalize_antigen(antigen: str | None) -> str | None:
    """Map a raw antigen string to its canonical short symbol, if known.

    Returns the canonical symbol when the input matches any pattern in
    :data:`CANONICAL_ANTIGEN_ALIASES` (case-insensitive substring,
    first match wins); returns the input unchanged otherwise. ``None``
    and empty strings pass through.

    Designed for messy free-text inputs from VDJdb's ``antigen.gene``
    and IEDB's ``Source Molecule`` — both routinely include
    ``[organism]`` suffixes, "protein" descriptors, capitalization
    variants, and parenthetical synonyms that should be ignored.
    """
    if antigen is None:
        return None
    raw = str(antigen).strip()
    if not raw:
        return raw
    lowered = raw.lower()
    for pattern, canonical in CANONICAL_ANTIGEN_ALIASES:
        if pattern in lowered:
            return canonical
    return raw


def canonicalize_antigens(antigens: pd.Series) -> pd.Series:
    """Vectorized :func:`canonicalize_antigen` over a Series."""
    return antigens.map(canonicalize_antigen)


def load_vdjdb(path: str | Path, verbose: bool = True) -> pd.DataFrame:
    """
    Load VDJdb database.

    VDJdb ships two on-disk shapes with **different row semantics**:

    - ``vdjdb_full.txt`` — *paired-chain* format: one row per
      clonotype, with ``cdr3.alpha`` and ``cdr3.beta`` as separate
      columns. This is what αβ matching needs.
    - ``vdjdb.txt`` — *long/slim* format: one row per **chain**, keyed
      by ``complex.id`` and ``gene`` (``TRA``/``TRB``). The ``cdr3``
      column holds an α-CDR3 when ``gene == TRA`` and a β-CDR3 when
      ``gene == TRB``. A naive ``cdr3 → cdr3_beta`` rename mixes the
      two and matches alpha sequences as if they were beta. To avoid
      that, the long-format path is filtered to ``gene == TRB`` rows
      so ``cdr3_beta`` only carries actual β-CDR3s; α rows are
      dropped (would need a ``complex.id`` pivot to recover αβ pairs,
      which is best deferred to the full file).

    When given a directory, looks for one of three canonical filenames
    in priority order: ``vdjdb_full.txt``, ``vdjdb.txt``,
    ``vdjdb.slim.txt``. Other files in the directory (metadata
    sidecars, processed views, etc.) are ignored — a directory that
    has none of these names raises with a hint to re-run
    ``tcrsift data download`` or pass the file path explicitly.

    Parameters
    ----------
    path : str or Path
        Path to VDJdb directory or single file.
    verbose : bool
        Print progress information.

    Returns
    -------
    pd.DataFrame
        VDJdb entries with standardized columns (``cdr3_alpha``,
        ``cdr3_beta``, ``epitope``, ``antigen_gene``, ``species``,
        ``mhc_allele``, ``database``, ``is_viral``). The
        ``cdr3_alpha`` column is empty when loaded from the
        long/slim format.
    """
    path = Path(path)

    if path.is_dir():
        db_file = _pick_vdjdb_file(path)
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

    # Format-specific normalization, dispatched by column presence:
    # the paired file has ``cdr3.beta``, the long file has
    # ``cdr3`` + ``gene``. ``cdr3.alpha`` alone (no β) doesn't occur
    # in real VDJdb files and isn't useful for matching — falls
    # through to the unrecognized-schema error.
    if "cdr3.beta" in df.columns:
        df = _normalize_vdjdb_paired(df)
    elif "cdr3" in df.columns and "gene" in df.columns:
        n_total = len(df)
        df = _normalize_vdjdb_long(df)
        if verbose:
            logger.info(
                f"  Long/slim VDJdb format: kept {len(df):,} TRB rows "
                f"out of {n_total:,} total (α rows discarded — load "
                f"vdjdb_full.txt for αβ matching)."
            )
    else:
        raise TCRsiftValidationError(
            f"Unrecognized VDJdb schema in {db_file}",
            hint=(
                "Expected either paired-chain columns "
                "(`cdr3.alpha` + `cdr3.beta`) or long-format columns "
                "(`cdr3` + `gene`). Got: "
                f"{list(df.columns)[:10]}{'...' if len(df.columns) > 10 else ''}"
            ),
        )

    df["database"] = "VDJdb"
    df["is_viral"] = _flag_viral(df)

    if verbose:
        logger.info(f"  Loaded {len(df):,} VDJdb entries ({df['is_viral'].sum():,} viral)")
    return df


def _normalize_vdjdb_paired(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize the paired-chain ``vdjdb_full.txt`` format.

    VDJdb's full export carries two species columns: the TCR donor's
    species and the epitope's source organism. We propagate
    ``antigen.species`` → ``species`` (the convention used throughout
    the codebase) and preserve the donor species as ``host_species``
    (#83 — needed for cross-species match-strength detection).
    """
    if "species" in df.columns:
        df["host_species"] = df["species"]
        if "antigen.species" in df.columns:
            df = df.drop(columns=["species"])

    column_mapping = {
        "cdr3.alpha": "cdr3_alpha",
        "cdr3.beta": "cdr3_beta",
        "antigen.epitope": "epitope",
        "antigen.gene": "antigen_gene",
        "mhc.a": "mhc_allele",
        "mhc.class": "mhc_class",
        "reference.id": "reference",
    }
    for old, new in column_mapping.items():
        if old in df.columns:
            df[new] = df[old]

    if "antigen.species" in df.columns:
        df["species"] = df["antigen.species"]
    return df


def _normalize_vdjdb_long(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize the long/slim ``vdjdb.txt`` format.

    The long format has one row per chain. We filter to ``gene == TRB``
    so the resulting ``cdr3_beta`` column only contains actual β-CDR3s
    (α-CDR3s share the same ``cdr3`` column under ``gene == TRA`` rows
    and would otherwise leak into ``cdr3_beta`` and produce wrong-chain
    matches). Recovering αβ pairs would require a pivot on
    ``complex.id`` — deferred; use ``vdjdb_full.txt`` for that.
    """
    # Exact-string match on the canonical VDJdb labels (``TRA``/``TRB``).
    # Any other gene (``TRG``/``TRD`` etc.) is silently dropped — αβ
    # matching can't consume them anyway.
    df = df[df["gene"] == "TRB"].copy()

    if "species" in df.columns:
        df["host_species"] = df["species"]
        if "antigen.species" in df.columns:
            df = df.drop(columns=["species"])

    column_mapping = {
        "cdr3": "cdr3_beta",
        "antigen.epitope": "epitope",
        "antigen.gene": "antigen_gene",
        "mhc.a": "mhc_allele",
        "mhc.class": "mhc_class",
        "reference.id": "reference",
    }
    for old, new in column_mapping.items():
        if old in df.columns:
            df[new] = df[old]

    if "antigen.species" in df.columns:
        df["species"] = df["antigen.species"]
    # No α data in this format; explicit empty column keeps the
    # downstream schema stable.
    df["cdr3_alpha"] = pd.NA
    return df


# Canonical VDJdb data filenames in resolution priority order. The
# release zip drops a swarm of processed sidecars (``vdjdb.meta.txt``,
# ``vdjdb.scored.txt``, ``vdjdb_full_*_broken.txt`` etc.) alongside
# these, but they have different schemas / are incomplete. We only
# look for the names below; everything else stays invisible to the
# loader. Users with a custom file should pass its path directly.
_VDJDB_CANONICAL_FILENAMES = ("vdjdb_full.txt", "vdjdb.txt", "vdjdb.slim.txt")


def _pick_vdjdb_file(dir_path: Path) -> Path:
    """Choose the canonical VDJdb data file from an extracted release dir.

    Priority:

    1. ``vdjdb_full.txt`` — paired αβ-per-row, the format αβ matching
       needs.
    2. ``vdjdb.txt`` — long/slim, one chain per row (β-only after the
       gene filter applied in ``_normalize_vdjdb_long``).
    3. ``vdjdb.slim.txt`` — smaller variant of the long format.

    A directory without any of these (e.g. the cache was corrupted or
    the release format changed) raises with a clear hint pointing back
    at ``tcrsift data download`` (#45).
    """
    for name in _VDJDB_CANONICAL_FILENAMES:
        candidate = dir_path / name
        if candidate.is_file():
            return candidate

    available = [f.name for f in dir_path.iterdir()][:15]
    raise TCRsiftValidationError(
        f"No canonical VDJdb data file found in directory: {dir_path}",
        hint=(
            f"Expected one of {list(_VDJDB_CANONICAL_FILENAMES)}. "
            f"Available files: {available}. "
            "Re-run `tcrsift data download --db vdjdb` to refresh, "
            "or pass the path to your VDJdb file directly with --vdjdb."
        ),
    )


def load_iedb(
    path: str | Path,
    epitope_path: str | Path | None = None,
) -> pd.DataFrame:
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
        Path to IEDB receptor file (v3 CSV or legacy flat TSV).
    epitope_path : str or Path, optional
        Path to the IEDB epitope-level CSV
        (``epitope_full_v3.csv`` — companion to the receptor file,
        cached by ``tcrsift data download --db iedb_epitope``).
        When provided, the receptor file's ``antigen_gene`` / ``species``
        fields are overridden with the epitope table's values
        wherever the epitope table has them. Empirically, the
        epitope-table names are shorter and more publication-canonical
        (e.g. ``Protein Tax-1`` vs the receptor file's
        ``transcriptional activator Tax``), which reduces downstream
        synonym sprawl (#54).

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

    if epitope_path is not None:
        df = _apply_iedb_epitope_overrides(df, Path(epitope_path))

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
            "epitope_iri": _pick("Epitope", "IEDB IRI"),
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


def _normalize_iedb_iri(s: pd.Series) -> pd.Series:
    """Normalize IEDB epitope IRIs to a comparable form.

    The receptor file (``tcr_full_v3.csv``) emits ``https://`` IRIs;
    the epitope file (``epitope_full_v3.csv``) emits ``http://``.
    Same epitope, different scheme — the join needs them aligned.
    """
    return s.fillna("").astype(str).str.replace("https://", "http://", regex=False)


def load_iedb_epitope_lookup(path: str | Path) -> pd.DataFrame:
    """Load the IEDB epitope-level table as a (canonical) name lookup.

    The epitope file has hierarchical CSV headers like the receptor
    file. Returns a deduplicated frame indexed by normalized IRI with
    columns ``antigen_gene`` and ``species`` — the values the
    receptor-level loader should defer to when present.

    Parameters
    ----------
    path : str or Path
        Path to ``epitope_full_v3.csv``.

    Returns
    -------
    pd.DataFrame
        Indexed by ``epitope_iri`` (string, ``http://...`` form),
        columns ``antigen_gene`` and ``species``.
    """
    raw = pd.read_csv(path, header=[0, 1], low_memory=False)

    iri_col = ("Epitope ID", "IEDB IRI")
    sm_col = ("Epitope", "Source Molecule")
    so_col = ("Epitope", "Source Organism")
    for required in (iri_col, sm_col, so_col):
        if required not in raw.columns:
            raise TCRsiftValidationError(
                f"IEDB epitope file at {path} is missing required column {required}",
                hint="Expected the v3 epitope CSV with hierarchical headers.",
            )

    lookup = pd.DataFrame(
        {
            "iri": _normalize_iedb_iri(raw[iri_col]),
            "antigen_gene": raw[sm_col],
            "species": raw[so_col],
        }
    )
    # Multiple rows per IRI exist (one per Related Object); keep the
    # row that actually carries a Source Molecule value (sorting NaNs
    # last and dropping later duplicates).
    lookup = (
        lookup.sort_values("antigen_gene", na_position="last")
        .drop_duplicates("iri", keep="first")
        .set_index("iri")
    )
    return lookup


def _apply_iedb_epitope_overrides(
    receptor_df: pd.DataFrame, epitope_path: Path
) -> pd.DataFrame:
    """Override receptor antigen_gene/species with epitope-table values.

    Strategy: for each receptor row that has an ``epitope_iri``, look
    it up in the epitope table. If the epitope table has a non-null
    value for that field, use it; otherwise keep the receptor value.
    This catches the ~80K receptor rows where the epitope table's
    name is shorter/more canonical, and the ~200 receptor rows where
    the receptor was blank but the epitope table had data (#54).
    """
    if "epitope_iri" not in receptor_df.columns:
        logger.warning(
            "IEDB receptor frame has no epitope_iri column; skipping "
            "epitope-table override (epitope_path was provided but "
            "the receptor loader didn't capture IRIs)."
        )
        return receptor_df

    logger.info(f"Loading IEDB epitope table from {epitope_path}")
    lookup = load_iedb_epitope_lookup(epitope_path)

    df = receptor_df.copy()
    rec_iri = _normalize_iedb_iri(df["epitope_iri"])

    for field in ("antigen_gene", "species"):
        overrides = rec_iri.map(lookup[field])
        # Keep receptor value where the epitope table is NaN.
        df[field] = overrides.where(overrides.notna(), df[field])

    n_changed = (rec_iri.isin(lookup.index)).sum()
    logger.info(
        f"  Applied epitope-table overrides to {n_changed:,} receptor rows."
    )
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

    Categories (precedence, low → high): ``unknown`` < ``other`` <
    ``viral`` / ``bacterial`` < ``self`` < ``tumor_self``. The tumor-self
    bucket is the last assignment and overrides *every* other category —
    a curated antigen-name match (MART-1, NY-ESO-1, MAGE, etc.) wins
    over whatever the species column says. In practice tumor antigens
    are Homo sapiens so this only overrides ``self``, but the rule is
    "trust the antigen name when it's on the curated list".

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
        is_viral |= species_lower.str.contains(pattern, na=False, regex=False)
    category[is_viral] = CATEGORY_VIRAL

    is_bacterial = pd.Series(False, index=species.index)
    for pattern in BACTERIAL_SPECIES_PATTERNS:
        is_bacterial |= species_lower.str.contains(pattern, na=False, regex=False)
    category[is_bacterial & ~is_viral] = CATEGORY_BACTERIAL

    # "Homo sapiens", "Homo sapiens (human)", "human" all bucket as self.
    is_self = species_lower.str.contains("homo sapiens", na=False, regex=False) | (
        species_lower == "human"
    )
    category[is_self & ~is_viral & ~is_bacterial] = CATEGORY_SELF

    # Tumor-self override: regex patterns anchored at word boundaries so
    # short tokens (her2, tert, psa, wt1) don't false-positive on
    # unrelated gene names containing them as substrings.
    is_tumor_self = pd.Series(False, index=species.index)
    for pattern in TUMOR_SELF_ANTIGEN_PATTERNS:
        is_tumor_self |= antigen_lower.str.contains(pattern, na=False, regex=True)
    category[is_tumor_self] = CATEGORY_TUMOR_SELF

    return category


def load_databases(
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    iedb_epitope_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    Load and combine multiple TCR databases.

    Parameters
    ----------
    vdjdb_path : str or Path, optional
        Path to VDJdb
    iedb_path : str or Path, optional
        Path to IEDB receptor table
    cedar_path : str or Path, optional
        Path to CEDAR
    iedb_epitope_path : str or Path, optional
        Path to IEDB epitope-level table (companion to ``iedb_path``).
        When provided, ``load_iedb`` defers to its shorter/more
        canonical antigen and organism names — see :func:`load_iedb`.

    Returns
    -------
    pd.DataFrame
        Combined database with standardized columns
    """
    dfs = []

    if vdjdb_path:
        dfs.append(load_vdjdb(vdjdb_path))
    if iedb_path:
        dfs.append(load_iedb(iedb_path, epitope_path=iedb_epitope_path))
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
        ``db_match_strength`` (``"ab"`` / ``"b_only"`` / ``"ab_cross"``
        / ``"b_only_cross"`` / None; ``_cross`` suffix indicates a
        non-human host-species match per #83),
        ``db_epitope``, ``db_protein``, ``db_protein_canonical``,
        ``db_species`` (antigen source organism),
        ``db_host_species`` (TCR donor organism, new in #83),
        ``db_mhc``, ``db_category`` (``viral`` / ``bacterial`` /
        ``self`` / ``tumor_self`` / ``other`` / ``unknown`` /
        ``contradictory`` — last when multiple matches disagree on
        an informative label per #83),
        ``db_database``, ``is_viral``.
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
    df["db_protein_canonical"] = None
    df["db_species"] = None
    df["db_host_species"] = None
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

    # Pre-classify the entire database once so the per-clone match path
    # picks db_category / db_protein_canonical as modes like any other
    # field — avoids invoking the classifiers with a 1-row Series per
    # match (#48 follow-up). `database.assign` returns a copy so the
    # caller's df isn't mutated.
    species_col = (
        database["species"]
        if "species" in database.columns
        else pd.Series([""] * len(database), index=database.index)
    )
    antigen_col = (
        database["antigen_gene"]
        if "antigen_gene" in database.columns
        else pd.Series([""] * len(database), index=database.index)
    )
    new_columns: dict[str, pd.Series] = {}
    if "db_category" not in database.columns:
        new_columns["db_category"] = classify_category(species_col, antigen_col)
    if "db_protein_canonical" not in database.columns:
        new_columns["db_protein_canonical"] = canonicalize_antigens(antigen_col)
    if new_columns:
        database = database.assign(**new_columns)

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
    (``.mode().iloc[0]``). Two exceptions (#83):

    - ``db_category`` is set to :data:`CATEGORY_CONTRADICTORY` when
      matches disagree on a non-null category (rather than silently
      picking a mode); the agreed value is used otherwise.
    - ``db_match_strength`` gets an ``_cross`` suffix when any matched
      row carries a non-human host species (``host_species``). Lets
      downstream filters distinguish a confident human-vs-human match
      from a cross-species curiosity.
    """
    if len(matches) == 0:
        return

    # Cross-species detection — applied before strength is recorded.
    host_species_match: str | None = None
    if "host_species" in matches.columns:
        hs = matches["host_species"].dropna()
        if len(hs) > 0:
            host_species_match = hs.mode().iloc[0]
            df.loc[idx, "db_host_species"] = host_species_match
            hs_lower = hs.astype(str).str.lower()
            any_non_human = (
                ~hs_lower.str.contains("homo sapiens", na=False, regex=False)
                & (hs_lower != "human")
            ).any()
            if any_non_human:
                strength = f"{strength}_cross"

    df.loc[idx, "db_match"] = True
    df.loc[idx, "db_match_strength"] = strength

    epitopes = matches["epitope"].dropna()
    if len(epitopes) > 0:
        df.loc[idx, "db_epitope"] = epitopes.mode().iloc[0]

    if "species" in matches.columns:
        species = matches["species"].dropna()
        if len(species) > 0:
            df.loc[idx, "db_species"] = species.mode().iloc[0]

    if "antigen_gene" in matches.columns:
        proteins = matches["antigen_gene"].dropna()
        if len(proteins) > 0:
            df.loc[idx, "db_protein"] = proteins.mode().iloc[0]

    if "db_protein_canonical" in matches.columns:
        canonicals = matches["db_protein_canonical"].dropna()
        if len(canonicals) > 0:
            df.loc[idx, "db_protein_canonical"] = canonicals.mode().iloc[0]

    if "mhc_allele" in matches.columns:
        mhcs = matches["mhc_allele"].dropna()
        if len(mhcs) > 0:
            df.loc[idx, "db_mhc"] = mhcs.mode().iloc[0]

    if "db_category" in matches.columns:
        cats = matches["db_category"].dropna()
        if len(cats) > 0:
            unique_cats = set(cats.unique())
            # ``unknown`` is the "we don't know" sentinel — don't let a
            # mixed bag of "unknown" + a real label trigger a
            # contradiction. Only flag when ≥ 2 *informative* categories
            # disagree.
            informative = unique_cats - {CATEGORY_UNKNOWN}
            if len(informative) > 1:
                df.loc[idx, "db_category"] = CATEGORY_CONTRADICTORY
            else:
                df.loc[idx, "db_category"] = cats.mode().iloc[0]

    df.loc[idx, "db_database"] = ";".join(matches["database"].unique())
    df.loc[idx, "is_viral"] = matches["is_viral"].any()

    if partial:
        df.loc[idx, "db_match_partial"] = True


def annotate_clonotypes(
    clonotypes: pd.DataFrame,
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    iedb_epitope_path: str | Path | None = None,
    match_by: str = "CDR3ab",
    match_strictness: str | None = None,
    exclude_viral: bool = False,
    flag_only: bool = False,
    *,
    add_publicness: bool = False,
    publicness_cdr3_col: str = "CDR3_beta",
    publicness_v_gene_col: str = "beta_v_gene",
    publicness_j_gene_col: str = "beta_j_gene",
) -> pd.DataFrame:
    """
    Main annotation function.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    vdjdb_path, iedb_path, cedar_path : str or Path, optional
        Paths to databases
    iedb_epitope_path : str or Path, optional
        Path to the IEDB epitope-level table (companion to
        ``iedb_path``). When provided, IEDB receptor-row antigen/species
        strings are overridden with the epitope-table's typically
        shorter and more canonical equivalents (#54).
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
    add_publicness : bool
        When True, add ``log10_pgen`` and ``publicness`` columns from
        :func:`tcrsift.pgen.annotate_publicness` (#58). High-publicness
        sequences are likely to be public; downstream callers can
        discount DB matches against them by multiplying any per-match
        score by ``(1 - publicness)``.
    publicness_cdr3_col, publicness_v_gene_col, publicness_j_gene_col : str
        Columns Pgen estimator should read.

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
        "db_protein_canonical": None,
        "db_species": None,
        "db_host_species": None,
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
        if add_publicness:
            df = _apply_publicness(
                df,
                cdr3_col=publicness_cdr3_col,
                v_gene_col=publicness_v_gene_col,
                j_gene_col=publicness_j_gene_col,
                log_summary=False,
            )
        return df

    # Load databases
    database = load_databases(
        vdjdb_path=vdjdb_path,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        iedb_epitope_path=iedb_epitope_path,
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

    # Publicness annotation (#58). Computed after match filtering so
    # the resulting Pgen/publicness columns align with the final row
    # set.
    if add_publicness:
        df = _apply_publicness(
            df,
            cdr3_col=publicness_cdr3_col,
            v_gene_col=publicness_v_gene_col,
            j_gene_col=publicness_j_gene_col,
            log_summary=True,
        )

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
