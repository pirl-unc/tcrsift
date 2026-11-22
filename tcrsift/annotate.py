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
import re
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from .pgen import annotate_publicness as _annotate_publicness
from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_dataframe,
    validate_file_exists,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


def _apply_publicness(
    df: pd.DataFrame,
    *,
    cdr3_col: str,
    v_gene_col: str,
    j_gene_col: str,
    log_summary: bool,
    species: str = "human",
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
        df,
        cdr3_col=cdr3_col,
        v_gene_col=v_gene_col,
        j_gene_col=j_gene_col,
        species=species,
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
    "mcmv",
    "ebv",
    "epstein-barr",
    "hiv",
    "human immunodeficiency",
    "htlv",
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
    "west nile",
    "lcmv",
    "lymphocytic choriomeningitis",
    "mammarenavirus choriomeningitidis",
    "vesicular stomatitis",
    "vsv",
    "sendai",
    "respirovirus",
    "vaccinia",
    "ectromelia",
    "murine leukemia",
    "mulv",
    "polyomavirus",
    "adenovirus",
    "norovirus",
    "rotavirus",
    "parvovirus",
]

# Tumor-associated self antigens. Patterns are case-insensitive regex
# fragments matched against antigen_gene / Source Molecule, anchored at
# word boundaries (``\b``) so short tokens like ``wt1`` / ``her2`` /
# ``tert`` / ``psa`` / ``dct`` only fire on standalone occurrences.
TUMOR_SELF_ANTIGEN_PATTERNS = [
    r"\bmart-?1\b",
    r"\bmlana\b",
    r"\bmelan-?a\b",
    r"\bny[- _]eso[- _]?1?\b",
    r"\bctag1[ab]?\b",
    r"\bctag2\b",
    r"\blage[- _]?[12]\b",
    r"\bcancer[/-]?testis\b",
    r"\bmage[- _]?[abc]\d*\b",
    r"\bgp100\b",
    r"\bpmel\b",
    r"\btyrosinase\b",
    r"\bdct\b",
    r"\bdopachrome tautomerase\b",
    r"\btrp[- _]?[12]\b",
    r"\btyrp[12]\b",
    r"\bceacam5\b",
    r"\bcarcinoembryonic antigen(?:\b|[- ]related cell adhesion molecule 5\b)",
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
    r"\bklk3\b",
    r"\bprostate[- ]specific antigen\b",
    r"\bssx2\b",
    r"\bprame\b",
    r"\bfolh1\b",
    r"\bpsma\b",
    r"\bmsln\b",
    r"\bmesothelin\b",
    r"\bafp\b",
    r"\balpha[- ]fetoprotein\b",
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

UNKNOWN_SPECIES_TOKENS = {
    "unknown",
    "unk",
    "na",
    "nan",
    "none",
    "notapplicable",
    "notavailable",
    "unclassified",
    "unspecified",
}

# Exact aliases keyed by lowercased alphanumeric text. This catches spelling,
# separator, and source-specific compact forms (e.g. VDJdb's HomoSapiens).
CANONICAL_SPECIES_ALIASES: dict[str, str] = {
    "human": "Homo sapiens",
    "homosapiens": "Homo sapiens",
    "homosapienshuman": "Homo sapiens",
    "mouse": "Mus musculus",
    "mice": "Mus musculus",
    "murine": "Mus musculus",
    "musmusculus": "Mus musculus",
    "rat": "Rattus norvegicus",
    "rattusnorvegicus": "Rattus norvegicus",
    "rhesusmacaque": "Macaca mulatta",
    "macacamulatta": "Macaca mulatta",
    "cynomolgusmacaque": "Macaca fascicularis",
    "macacafascicularis": "Macaca fascicularis",
    "pig": "Sus scrofa",
    "swine": "Sus scrofa",
    "susscrofa": "Sus scrofa",
    "cow": "Bos taurus",
    "cattle": "Bos taurus",
    "bostaurus": "Bos taurus",
    "cmv": "CMV",
    "hcmv": "CMV",
    "humancmv": "CMV",
    "humancytomegalovirus": "CMV",
    "humanherpesvirus5": "CMV",
    "hhv5": "CMV",
    "mcmv": "Mouse CMV",
    "mousecmv": "Mouse CMV",
    "murinecmv": "Mouse CMV",
    "mousecytomegalovirus": "Mouse CMV",
    "muridbetaherpesvirus1": "Mouse CMV",
    "ebv": "EBV",
    "epsteinbarr": "EBV",
    "epsteinbarrvirus": "EBV",
    "humanherpesvirus4": "EBV",
    "hhv4": "EBV",
    "hiv": "HIV",
    "humanimmunodeficiencyvirus": "HIV",
    "hiv1": "HIV-1",
    "humanimmunodeficiencyvirus1": "HIV-1",
    "hiv2": "HIV-2",
    "humanimmunodeficiencyvirus2": "HIV-2",
    "htlv1": "HTLV-1",
    "humantlymphotropicvirus1": "HTLV-1",
    "humantcellleukemiavirus1": "HTLV-1",
    "sarscov2": "SARS-CoV-2",
    "severeacuterespiratorysyndromecoronavirus2": "SARS-CoV-2",
    "2019ncov": "SARS-CoV-2",
    "sarscov": "SARS-CoV-1",
    "severeacuterespiratorysyndromecoronavirus": "SARS-CoV-1",
    "influenza": "Influenza",
    "flu": "Influenza",
    "influenzaa": "Influenza A",
    "influenzaavirus": "Influenza A",
    "influenzaavirusa": "Influenza A",
    "influenzab": "Influenza B",
    "influenzabvirus": "Influenza B",
    "hsv1": "HSV-1",
    "herpessimplexvirus1": "HSV-1",
    "humanherpesvirus1": "HSV-1",
    "hsv2": "HSV-2",
    "herpessimplexvirus2": "HSV-2",
    "humanherpesvirus2": "HSV-2",
    "hpv": "HPV",
    "humanpapillomavirus": "HPV",
    "hbv": "HBV",
    "hepatitisb": "HBV",
    "hepatitisbvirus": "HBV",
    "hcv": "HCV",
    "hepatitisc": "HCV",
    "hepatitiscvirus": "HCV",
    "dengue": "Dengue virus",
    "denguevirus": "Dengue virus",
    "zika": "Zika virus",
    "zikavirus": "Zika virus",
    "yellowfever": "Yellow fever virus",
    "yellowfevervirus": "Yellow fever virus",
    "gallusgallus": "Gallus gallus",
    "chicken": "Gallus gallus",
    "lcmv": "LCMV",
    "lymphocyticchoriomeningitisvirus": "LCMV",
    "lymphocyticchoriomeningitismammarenavirus": "LCMV",
    "mammarenaviruschoriomeningitidis": "LCMV",
    "sendai": "Sendai virus",
    "sendaivirus": "Sendai virus",
    "murinerespirovirus": "Sendai virus",
    "vsv": "VSV",
    "vesicularstomatitisvirus": "VSV",
    "vesicularstomatitisindianavirus": "VSV",
    "vaccinia": "Vaccinia virus",
    "vacciniavirus": "Vaccinia virus",
    "ectromelia": "Ectromelia virus",
    "ectromeliavirus": "Ectromelia virus",
    "mulv": "MuLV",
    "murineleukemiavirus": "MuLV",
    "friendmurineleukemiavirus": "MuLV",
    "moloneymurineleukemiavirus": "MuLV",
    "mousepolyomavirus": "Mouse polyomavirus",
}

# Substring aliases are ordered most-specific first so strains and longer
# source names collapse before their broader family names.
CANONICAL_SPECIES_CONTAINS: tuple[tuple[str, str], ...] = (
    ("severeacuterespiratorysyndromecoronavirus2", "SARS-CoV-2"),
    ("sarscov2", "SARS-CoV-2"),
    ("2019ncov", "SARS-CoV-2"),
    ("mousecytomegalovirus", "Mouse CMV"),
    ("murinecytomegalovirus", "Mouse CMV"),
    ("muridbetaherpesvirus1", "Mouse CMV"),
    ("humancytomegalovirus", "CMV"),
    ("humanherpesvirus5", "CMV"),
    ("cytomegalovirus", "CMV"),
    ("epsteinbarrvirus", "EBV"),
    ("epsteinbarr", "EBV"),
    ("humanherpesvirus4", "EBV"),
    ("humanimmunodeficiencyvirus1", "HIV-1"),
    ("humanimmunodeficiencyvirus2", "HIV-2"),
    ("humanimmunodeficiencyvirus", "HIV"),
    ("sarscov", "SARS-CoV-1"),
    ("influenzaavirusb", "Influenza B"),
    ("influenzab", "Influenza B"),
    ("influenzaavirusa", "Influenza A"),
    ("influenzaa", "Influenza A"),
    ("influenza", "Influenza"),
    ("herpessimplexvirus1", "HSV-1"),
    ("herpessimplexvirus2", "HSV-2"),
    ("humanpapillomavirus", "HPV"),
    ("hepatitisbvirus", "HBV"),
    ("hepatitiscvirus", "HCV"),
    ("denguevirus", "Dengue virus"),
    ("zikavirus", "Zika virus"),
    ("yellowfevervirus", "Yellow fever virus"),
    ("lymphocyticchoriomeningitis", "LCMV"),
    ("mammarenaviruschoriomeningitidis", "LCMV"),
    ("vesicularstomatitis", "VSV"),
    ("sendaivirus", "Sendai virus"),
    ("murinerespirovirus", "Sendai virus"),
    ("vacciniavirus", "Vaccinia virus"),
    ("ectromeliavirus", "Ectromelia virus"),
    ("murineleukemiavirus", "MuLV"),
    ("gallusgallus", "Gallus gallus"),
    ("musmusculus", "Mus musculus"),
    ("homosapiens", "Homo sapiens"),
    ("c57bl", "Mus musculus"),
    ("balbc", "Mus musculus"),
    ("dba2", "Mus musculus"),
    ("fvbn", "Mus musculus"),
    ("129sv", "Mus musculus"),
    ("sjlj", "Mus musculus"),
    ("nodshi", "Mus musculus"),
)

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
DB_CATEGORY_DETAIL_CONTRADICTORY_DOMINANT = "contradictory_dominant"


# Valid strictness modes for match_clonotypes. The function dispatches
# on this string directly — there's no further internal translation.
MATCH_STRICTNESS_MODES = ("strict_ab", "ab_with_partial", "b_only")

# Valid match modes for match_clonotypes (#57). ``"exact"`` is the
# pre-2.1 behaviour; ``"levenshtein"`` adds β-only neighbor search at
# the given ``max_distance`` for clones that don't hit exactly. TCRdist
# is intentionally out-of-scope here — it'd pull an optional dep tree.
MATCH_MODES = ("exact", "levenshtein")


def _canonical_v_gene(value: object) -> str:
    """Normalize a V-gene call so VDJdb / IEDB / CellRanger spellings
    compare cleanly (#56).

    Strips allele suffix (``TRBV6-1*01`` → ``TRBV6-1``), uppercases,
    and trims whitespace. Returns ``""`` when the input isn't a string
    or is empty / NaN so the downstream comparison can be a simple
    ``a == b and a``.
    """
    if not isinstance(value, str):
        return ""
    base = value.strip().upper().split("*", 1)[0]
    return base


def _canonical_variants(s: str):
    """Yield ``{s} ∪ {s with one character deleted at each position}``.

    Two strings within Levenshtein distance 1 of each other always
    share at least one element of this set:
      - identical strings → ``s`` itself is in both sets;
      - single substitution → the deletion at the differing position
        matches in both sets;
      - single indel → the longer string's deletion-at-the-extra-
        position equals the shorter string, which is in its own
        canonical set.

    Used to build a small candidate set per query CDR3 in O(L) lookup
    instead of brute-forcing pairwise distance against every DB entry.
    """
    if not s:
        return
    yield s
    for i in range(len(s)):
        yield s[:i] + s[i + 1 :]


def _levenshtein_distance_at_most_1(s: str, t: str) -> int | None:
    """Return the Levenshtein distance when it's ≤ 1, else ``None``.

    Fast path for the post-canonical-set candidate filtering: we only
    care about distance ≤ 1 (max_distance for #57's default fuzzy
    matching), and the candidate set is already filtered to near-
    neighbors. Returns 0 for identical strings, 1 for single
    substitution / single indel, ``None`` for everything else.
    """
    if s == t:
        return 0
    ls, lt = len(s), len(t)
    if abs(ls - lt) > 1:
        return None
    if ls == lt:
        # Substitution path: count mismatches; bail if >1.
        diffs = 0
        for a, b in zip(s, t):
            if a != b:
                diffs += 1
                if diffs > 1:
                    return None
        return diffs if diffs <= 1 else None
    # Indel path: short is one shorter than long; find the insertion
    # position and verify the surrounding chars match.
    short, long_ = (s, t) if ls < lt else (t, s)
    i = j = 0
    while i < len(short) and short[i] == long_[j]:
        i += 1
        j += 1
    # Skip one character in long_ (the inserted one).
    j += 1
    while i < len(short):
        if short[i] != long_[j]:
            return None
        i += 1
        j += 1
    return 1


def _build_cdr3_neighbor_index(cdr3_series: pd.Series) -> dict[str, set[str]]:
    """Build a deletion-canonical index over a CDR3 column for O(L)
    Levenshtein-1 candidate lookup (#57).

    Map: canonical variant → set of DB CDR3 source strings that
    produced it. To find Lev-1 neighbors of a query CDR3, compute its
    canonical variants and union the indexed sets — the result is
    every DB CDR3 within Lev distance ≤ 1 (with some false positives
    at Lev distance 2 that the caller filters out via
    :func:`_levenshtein_distance_at_most_1`).
    """
    index: dict[str, set[str]] = {}
    for cdr3 in cdr3_series.dropna().unique():
        if not isinstance(cdr3, str) or not cdr3:
            continue
        for variant in _canonical_variants(cdr3):
            index.setdefault(variant, set()).add(cdr3)
    return index


def _find_lev_neighbors(
    query: str,
    index: dict[str, set[str]],
    max_distance: int,
) -> list[tuple[str, int]]:
    """Return ``[(neighbor_cdr3, distance), ...]`` sorted by distance,
    excluding the exact match (caller handles exact separately).

    Uses the canonical-variant index for candidate generation, then
    confirms distance with :func:`_levenshtein_distance_at_most_1`.
    """
    if not query:
        return []
    candidates: set[str] = set()
    for variant in _canonical_variants(query):
        candidates.update(index.get(variant, set()))
    candidates.discard(query)  # exact match handled separately
    results: list[tuple[str, int]] = []
    for cand in candidates:
        d = _levenshtein_distance_at_most_1(query, cand)
        if d is not None and 0 < d <= max_distance:
            results.append((cand, d))
    results.sort(key=lambda x: x[1])
    return results


def _species_alias_key(value: str) -> str:
    """Collapse punctuation/case so source-specific spellings compare cleanly."""
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _fallback_species_label(value: str) -> str | None:
    """Return a stable display label for species not in the alias table."""
    text = re.sub(r"[/_]+", " ", value).strip()
    text = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", text)
    text = re.sub(r"\s*\([^)]*\)\s*", " ", text)
    text = re.sub(r"\s*\[[^\]]*\]\s*", " ", text)
    text = re.sub(r"\s+", " ", text).strip(" ;,")
    if not text:
        return None

    tokens = text.split()
    if len(tokens) >= 2 and tokens[0].isalpha() and tokens[1].isalpha():
        return " ".join([tokens[0].capitalize(), tokens[1].lower(), *tokens[2:]])
    return text


def canonicalize_species(value: object) -> str | None:
    """Canonicalize a species/source-organism label.

    Public TCR databases use a mixture of compact taxonomy names
    (``HomoSapiens``), informal host labels (``human`` / ``mouse``),
    pathogen abbreviations (``CMV`` / ``EBV``), and long NCBI organism
    names. This helper maps common aliases to one stable display symbol
    and normalizes punctuation/case for less common binomials.
    """
    if pd.isna(value):
        return None

    text = str(value).strip()
    if not text:
        return None

    text = re.sub(r"\s+", " ", text)
    key = _species_alias_key(text)
    if key in UNKNOWN_SPECIES_TOKENS:
        return None
    if key in CANONICAL_SPECIES_ALIASES:
        return CANONICAL_SPECIES_ALIASES[key]

    for needle, canonical in CANONICAL_SPECIES_CONTAINS:
        if needle in key:
            return canonical

    fallback = _fallback_species_label(text)
    if fallback is None:
        return None

    fallback_key = _species_alias_key(fallback)
    if fallback_key in UNKNOWN_SPECIES_TOKENS:
        return None
    return CANONICAL_SPECIES_ALIASES.get(fallback_key, fallback)


def canonicalize_species_labels(species: pd.Series) -> pd.Series:
    """Vectorized wrapper for :func:`canonicalize_species`."""
    return species.map(canonicalize_species)


def _canonicalize_database_species(database: pd.DataFrame) -> pd.DataFrame:
    """Normalize species columns on a match-time copy of the DB frame."""
    species_columns = [
        col for col in ("species", "host_species") if col in database.columns
    ]
    if not species_columns:
        return database
    return database.assign(
        **{col: canonicalize_species_labels(database[col]) for col in species_columns}
    )


def _format_category_counts(counts: pd.Series) -> str:
    """Return deterministic ``category:count`` audit text."""
    return ";".join(f"{category}:{int(count)}" for category, count in counts.items())


def _summarize_match_categories(
    categories: pd.Series,
) -> tuple[str | None, str | None, str | None, str | None]:
    """Summarize category calls across all matched DB rows.

    Returns ``(db_category, db_category_dominant, db_category_detail,
    db_category_counts)``. ``unknown`` rows are ignored whenever at
    least one informative category is present; they only surface when
    every matched row is unknown.
    """
    cats = categories.dropna().astype(str)
    if len(cats) == 0:
        return None, None, None, None

    counts = cats.value_counts(sort=False)
    counts = counts.sort_index().sort_values(ascending=False, kind="stable")
    informative_counts = counts.drop(labels=[CATEGORY_UNKNOWN], errors="ignore")
    if len(informative_counts) == 0:
        counts_text = _format_category_counts(counts)
        return CATEGORY_UNKNOWN, None, CATEGORY_UNKNOWN, counts_text

    counts = informative_counts
    counts_text = _format_category_counts(counts)

    top_count = informative_counts.iloc[0]
    top_categories = list(informative_counts[informative_counts == top_count].index)
    dominant = top_categories[0] if len(top_categories) == 1 else None

    if len(informative_counts) == 1:
        category = informative_counts.index[0]
        return category, category, category, counts_text

    if dominant is not None:
        detail = f"{DB_CATEGORY_DETAIL_CONTRADICTORY_DOMINANT}_{dominant}"
        return CATEGORY_CONTRADICTORY, dominant, detail, counts_text

    return CATEGORY_CONTRADICTORY, None, CATEGORY_CONTRADICTORY, counts_text


# Word-bounded regex rules for canonicalizing protein / proteoform / antigen
# labels (db_protein_canonical). Distinct from encoding gene symbols
# (db_antigen_gene, resolved by resolve_antigen_gene): e.g. CTAG1A and CTAG1B
# are two paralogous X-linked genes that encode the identical NY-ESO-1 protein,
# whereas CTAG2 encodes LAGE-1.
_HOST_SPECIES_SET = frozenset(
    {"Homo sapiens", "Mus musculus", "Macaca mulatta", "Macaca fascicularis", "Rattus norvegicus"}
)

_SELF_TUMOR_PROTEIN_RULES: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"melanoma antigen recognized by t-cells|\bmelan-?a\b|\bmlana\b|\bmart-?1\b", re.I), "MART-1"),
    (re.compile(r"cancer/testis antigen 2\b|\bctag2\b|\blage[- _]?1[ab]?\b", re.I), "LAGE-1"),
    (re.compile(r"cancer/testis antigen 1\b|\bctag1[ab]?\b|\bny[- _]eso[- _]?1\b|\blage[- _]?2[ab]?\b", re.I), "NY-ESO-1"),
    (re.compile(r"melanocyte protein pmel|melanocyte-specific secreted glycoprotein|\bgp100\b|\bpmel\b|\bsilv\b", re.I), "gp100"),
    (re.compile(r"dopachrome tautomerase|tyrosinase[- _]?related protein[- _]?2\b|\bdct\b|\btrp[- _]?2\b|\btyrp2\b", re.I), "TRP-2"),
    (re.compile(r"5,6-dihydroxyindole-2-carboxylic acid oxidase|tyrosinase[- _]?related protein[- _]?1\b|\bgp75\b|\btyrp1\b|\btrp[- _]?1\b", re.I), "TRP-1"),
    (re.compile(r"\btyrosinase\b|\btyr\b", re.I), "Tyrosinase"),
    (re.compile(r"\bwilms[' ]*tumor\b|\bwt1\b", re.I), "WT1"),
    (re.compile(r"\btelomerase\b|\bh?tert\b", re.I), "TERT"),
    (re.compile(r"receptor tyrosine-protein kinase erbb-2|\berbb2\b|\bher2(?:/neu)?\b", re.I), "HER2"),
    (re.compile(r"\bsurvivin\b|\bbirc5\b", re.I), "Survivin"),
    (re.compile(r"prostate[- ]specific antigen|\bklk3\b|\bpsa\b", re.I), "PSA"),
    (re.compile(r"\bny-?br-?1\b|\bankrd30a\b", re.I), "NY-BR-1"),
    (re.compile(r"bone marrow stromal antigen 2|\btetherin\b|\bbst2\b", re.I), "BST2"),
    (re.compile(r"preferentially expressed antigen in melanoma|\bprame\b", re.I), "PRAME"),
    (re.compile(r"\bssx2\b", re.I), "SSX2"),
    (re.compile(r"\bmesothelin\b|\bmsln\b", re.I), "Mesothelin"),
    (re.compile(r"alpha[- ]fetoprotein|\bafp\b", re.I), "AFP"),
    (re.compile(r"\bovalbumin\b|\bserpinb14\b|\bova\b", re.I), "OVA"),
    (re.compile(r"hen egg[- ]white lysozyme|\bhel\b", re.I), "HEL"),
    (re.compile(r"keyhole limpet hemocyanin|\bklh\b", re.I), "KLH"),
    (re.compile(r"\binsulin\b", re.I), "Insulin"),
    (re.compile(r"beta-2-microglobulin|\bb2m\b", re.I), "B2M"),
    (re.compile(r"hla class i histocompatibility antigen|mhc class i(?: protein)?\b", re.I), "MHC class I"),
    (re.compile(r"hla class ii|mhc class ii\b", re.I), "MHC class II"),
)

_MAGE_NUMBERED_RE = re.compile(r"\bmage[- _]?([abc])(\d+)\b", re.I)
_MAGE_A_BARE_RE = re.compile(r"\bmage[- _]?a\b", re.I)
_CEA_MATCH_RE = re.compile(r"\bceacam5\b|\bcarcinoembryonic\b", re.I)
_CEA_OTHER_PARALOG_RE = re.compile(
    r"\bceacam(?!5\b)\d+\b|adhesion molecule (?!5\b)\d+\b", re.I
)
_INLINE_BRACKET_ORG_RE = re.compile(r"\[([^\]]+)\]\s*$")

# Back-compat tuple of (substring, canonical) retained for external callers that
# inspect CANONICAL_ANTIGEN_ALIASES directly.
CANONICAL_ANTIGEN_ALIASES: tuple[tuple[str, str], ...] = (
    ("melan-a", "MART-1"),
    ("mlana", "MART-1"),
    ("mart-1", "MART-1"),
    ("ctag1a", "NY-ESO-1"),
    ("ctag1b", "NY-ESO-1"),
    ("ny-eso-1", "NY-ESO-1"),
    ("ctag2", "LAGE-1"),
    ("lage-1", "LAGE-1"),
    ("gp100", "gp100"),
    ("pmel", "gp100"),
    ("dct", "TRP-2"),
    ("trp-2", "TRP-2"),
    ("ceacam5", "CEA"),
    ("wt1", "WT1"),
    ("tert", "TERT"),
    ("erbb2", "HER2"),
    ("her2", "HER2"),
    ("birc5", "Survivin"),
    ("klk3", "PSA"),
    ("bst2", "BST2"),
    ("ovalbumin", "OVA"),
    ("pp65", "CMV pp65"),
    ("cmv ie1", "CMV IE1"),
    ("ebna-3", "EBV EBNA-3"),
    ("ebna-1", "EBV EBNA-1"),
    ("bzlf1", "EBV BZLF1"),
    ("bmlf1", "EBV BMLF1"),
    ("matrix protein 1", "Flu M1"),
    ("nucleoprotein", "Flu NP"),
    ("spike glycoprotein", "SARS-CoV-2 Spike"),
    ("orf1ab", "SARS-CoV-2 ORF1ab"),
    ("nucleocapsid", "SARS-CoV-2 N"),
    ("protein tax-1", "HTLV-1 Tax"),
    ("ns3", "HCV NS3"),
    ("b2m", "B2M"),
)


def _extract_org_hint(raw: str, species: object = None) -> str | None:
    """Resolve canonical organism context from explicit ``species`` or inline ``[...]``."""
    org = canonicalize_species(species) if species is not None and not pd.isna(species) else None
    if org:
        return org
    m = _INLINE_BRACKET_ORG_RE.search(raw)
    if m:
        return canonicalize_species(m.group(1))
    return None


def canonicalize_antigen(
    antigen: str | None, species: str | None = None
) -> str | None:
    """Map a raw antigen/protein string to its canonical protein/antigen label.

    Strictly word-bounded and organism-aware:

    - Distinguishes proteins/proteoforms from encoding genes: both ``CTAG1A`` and
      ``CTAG1B`` canonicalize to ``"NY-ESO-1"``, whereas ``CTAG2`` canonicalizes to
      ``"LAGE-1"``; ``CEACAM5`` and ``Carcinoembryonic antigen`` canonicalize to
      ``"CEA"`` while ``CEACAM1``/``CEACAM6`` stay untouched; ``ERBB2`` → ``"HER2"``
      while ``ERBB2IP`` stays untouched; ``MAGE-A10`` / ``MAGEA12`` keep their
      exact paralog numbers rather than collapsing into ``MAGE-A1``.
    - Disambiguates shared viral protein names (``Nucleoprotein``,
      ``Nucleocapsid``, ``NS3``, ``Spike``, ``IE1``) using ``species`` (or an
      inline ``[organism]`` suffix), and never maps host self proteins such as
      human ``Small nuclear ribonucleoprotein G`` (``SNRPG``) onto viral labels.
    """
    if pd.isna(antigen):
        return None
    raw = str(antigen).strip()
    if not raw:
        return ""
    org = _extract_org_hint(raw, species)
    is_host_org = org in _HOST_SPECIES_SET

    # Strip trailing [organism] before matching so organism names inside brackets
    # don't accidentally match protein patterns.
    clean = _INLINE_BRACKET_ORG_RE.sub("", raw).strip() or raw

    # 1. Dynamic numbered MAGE family (MAGE-A1..A12, MAGE-B*, MAGE-C*).
    m_mage = _MAGE_NUMBERED_RE.search(clean)
    if m_mage:
        return f"MAGE-{m_mage.group(1).upper()}{int(m_mage.group(2))}"
    if _MAGE_A_BARE_RE.search(clean):
        return "MAGE-A"

    # 2. CEA / CEACAM5 (excluding other CEACAM paralogs like CEACAM1/CEACAM6).
    if _CEA_MATCH_RE.search(clean) and not _CEA_OTHER_PARALOG_RE.search(clean):
        return "CEA"

    # 3. Curated tumor/self/model proteins.
    for rx, canonical in _SELF_TUMOR_PROTEIN_RULES:
        if rx.search(clean):
            return canonical

    # Host proteins that didn't hit a curated self/tumor alias stay verbatim
    # (never fall through to viral protein heuristics).
    if is_host_org:
        return raw

    lowered = clean.lower()

    # 4. Murine viral antigens (LCMV, MCMV, MuLV, Sendai, VSV).
    if org == "Mouse CMV" or re.search(r"\bmcmv\b", lowered):
        if re.search(r"\bm45\b", lowered):
            return "MCMV M45"
        if re.search(r"\bm38\b", lowered):
            return "MCMV M38"
        if re.search(r"\bm139\b", lowered):
            return "MCMV m139"
        if re.search(r"\bm57\b", lowered):
            return "MCMV M57"
        if re.search(r"\bie1\b|immediate[- ]early", lowered):
            return "MCMV IE1"

    if org == "LCMV" or re.search(r"\blcmv\b", lowered):
        if re.search(r"glycoprotein|\bgpc?\b|\bgp(?:33|61|276)\b", lowered):
            return "LCMV GP"
        if re.search(r"\bnucleoprotein\b|\bnucleocapsid\b|\bnp\b", lowered):
            return "LCMV NP"

    if org == "MuLV" or re.search(r"\bmulv\b|murine leukemia", lowered):
        if re.search(r"\bp15e\b", lowered):
            return "MuLV p15E"
        if re.search(r"\bgp70\b|\bah1\b", lowered):
            return "MuLV gp70"
        if re.search(r"\benv(?:elope)?\b", lowered):
            return "MuLV Env"
    elif re.search(r"\bp15e\b", lowered):
        return "MuLV p15E"

    # 5. Human & shared viral / bacterial antigens (word-bounded + species-aware).
    if re.search(r"65 kda phosphoprotein|phosphoprotein 65|\bcmv pp65\b|\bpp65\b|\bul83\b", lowered):
        return "CMV pp65"
    if re.search(
        r"55 kda immediate-early|regulatory protein ie1|immediate[- ]early (?:protein )?1\b|\bcmv ie1\b|\bul123\b",
        lowered,
    ):
        return "MCMV IE1" if org == "Mouse CMV" else "CMV IE1"

    if re.search(r"nuclear antigen ebna-?3|\bebna-?3[abc]?\b", lowered):
        return "EBV EBNA-3"
    if re.search(r"nuclear antigen ebna-?1|\bebna-?1\b", lowered):
        return "EBV EBNA-1"
    if re.search(r"\bbzlf1\b", lowered):
        return "EBV BZLF1"
    if re.search(r"\bbmlf1\b", lowered):
        return "EBV BMLF1"
    if re.search(r"\bbrlf1\b", lowered):
        return "EBV BRLF1"
    if re.search(r"\blmp1\b", lowered):
        return "EBV LMP1"
    if re.search(r"\blmp2a?\b", lowered):
        return "EBV LMP2"

    if re.search(r"matrix protein 1\b|matrix protein m1\b|\bm1 protein\b|influenza a virus matrix", lowered):
        return "Flu M1"
    if re.search(r"\bhemagglutinin\b", lowered):
        return "Flu HA"

    if re.search(r"\bnucleoprotein\b|\bnp protein\b", lowered):
        if org in (None, "", "Influenza", "Influenza A", "Influenza B"):
            return "Flu NP"
        if org == "Sendai virus":
            return "Sendai NP"
        if org == "VSV":
            return "VSV N"
        if org == "SARS-CoV-2":
            return "SARS-CoV-2 N"
        if org == "SARS-CoV-1":
            return "SARS-CoV-1 N"
        return f"{org} NP"

    if re.search(r"replicase polyprotein 1ab?\b|orf1ab(?: polyprotein)?\b", lowered):
        return "SARS-CoV-1 ORF1ab" if org == "SARS-CoV-1" else "SARS-CoV-2 ORF1ab"

    if re.search(r"spike glycoprotein|surface glycoprotein|\bspike protein\b|sars-cov-2 spike", lowered):
        if org == "SARS-CoV-1":
            return "SARS-CoV-1 Spike"
        if org == "VSV":
            return "VSV G"
        if org in (None, "", "SARS-CoV-2") or (
            species is None and org in ("Pangolin coronavirus", "Bat coronavirus")
        ):
            return "SARS-CoV-2 Spike"
        return f"{org} Spike"

    if re.search(r"\bnucleocapsid\b", lowered):
        if org == "SARS-CoV-1":
            return "SARS-CoV-1 N"
        if org == "HBV":
            return "HBV Core"
        if org == "VSV":
            return "VSV N"
        if org in (None, "", "SARS-CoV-2"):
            return "SARS-CoV-2 N"
        return f"{org} N"

    if re.search(r"\bgag polyprotein\b|\bhiv[- ]?1? gag\b", lowered):
        return "HIV Gag"
    if re.search(r"\bnef protein\b|\bhiv[- ]?1? nef\b", lowered):
        return "HIV Nef"
    if re.search(r"protein tax-1\b|transcriptional activator tax\b|\btax-1\b", lowered):
        return "HTLV-1 Tax"

    if re.search(r"\bns3\b|non[- ]?structural protein (?:ns)?3\b", lowered):
        if org == "Dengue virus":
            return "Dengue NS3"
        if org == "Zika virus":
            return "Zika NS3"
        if org == "Yellow fever virus":
            return "YFV NS3"
        if org in (None, "", "HCV"):
            return "HCV NS3"
        return f"{org} NS3"
    if re.search(r"\bns5b\b", lowered):
        return "HCV NS5B"

    if re.search(r"\besat-?6\b|\besxa\b", lowered):
        return "ESAT-6"
    if re.search(r"\blisteriolysin\b|\bhly\b", lowered):
        return "LLO"

    return raw


def canonicalize_antigens(
    antigens: pd.Series, species: pd.Series | str | None = None
) -> pd.Series:
    """Vectorized :func:`canonicalize_antigen` over a Series (with optional species)."""
    if species is None or isinstance(species, str):
        return antigens.map(lambda a: canonicalize_antigen(a, species=species))
    return pd.Series(
        [canonicalize_antigen(a, s) for a, s in zip(antigens, species)],
        index=antigens.index,
        dtype=object,
    )


# Official HGNC and MGI gene symbols for curated endogenous self/tumor proteins.
# Note NY-ESO-1 is encoded identically by two X-linked paralogs (CTAG1A, CTAG1B);
# when the raw DB string names a specific paralog we preserve it, otherwise we
# emit "CTAG1A;CTAG1B" so gene-level ambiguity is explicit rather than erased.
_CANONICAL_PROTEIN_TO_GENE_HGNC: dict[str, str] = {
    "MART-1": "MLANA",
    "LAGE-1": "CTAG2",
    "NY-ESO-1": "CTAG1A;CTAG1B",
    "gp100": "PMEL",
    "TRP-2": "DCT",
    "TRP-1": "TYRP1",
    "Tyrosinase": "TYR",
    "CEA": "CEACAM5",
    "WT1": "WT1",
    "TERT": "TERT",
    "HER2": "ERBB2",
    "Survivin": "BIRC5",
    "PSA": "KLK3",
    "NY-BR-1": "ANKRD30A",
    "BST2": "BST2",
    "PRAME": "PRAME",
    "SSX2": "SSX2",
    "Mesothelin": "MSLN",
    "AFP": "AFP",
    "Insulin": "INS",
    "B2M": "B2M",
}

_CANONICAL_PROTEIN_TO_GENE_MGI: dict[str, str | None] = {
    "MART-1": "Mlana",
    "LAGE-1": "Ctag2",
    "NY-ESO-1": None,  # Human-specific X-linked duplication (no Mus musculus ortholog)
    "gp100": "Pmel",
    "TRP-2": "Dct",
    "TRP-1": "Tyrp1",
    "Tyrosinase": "Tyr",
    "CEA": None,       # CEACAM5 is human/primate-specific
    "WT1": "Wt1",
    "TERT": "Tert",
    "HER2": "Erbb2",
    "Survivin": "Birc5",
    "PSA": None,       # KLK3 is human/primate-specific
    "NY-BR-1": None,   # ANKRD30A has no 1:1 mouse ortholog
    "BST2": "Bst2",
    "PRAME": None,     # Expanded Pramel family in mouse; no single 1:1 ortholog
    "SSX2": None,      # Human SSX2 has no 1:1 mouse ortholog
    "Mesothelin": "Msln",
    "AFP": "Afp",
    "Insulin": "Ins2",
    "B2M": "B2m",
}

_BARE_HGNC_SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9]{1,14}$")
_BARE_MGI_SYMBOL_RE = re.compile(r"^[A-Z][a-z0-9]{1,14}$")
_NON_GENE_BARE_TOKENS = frozenset(
    {"MHC", "HLA", "OVA", "HEL", "KLH", "UNKNOWN", "NONE", "NA", "NAN"}
)
_NON_GENE_CANONICAL_ANTIGENS = frozenset(
    {"OVA", "HEL", "KLH", "MAGE-A", "MHC class I", "MHC class II"}
)


def resolve_antigen_gene(
    antigen: str | None,
    canonical_protein: str | None = None,
    species: str | None = None,
) -> str | None:
    """Resolve the official encoding gene symbol (HGNC/MGI) for an endogenous antigen.

    Returns ``None`` for viral, bacterial, and foreign exogenous antigens (such as
    chicken ``OVA`` or ``CMV pp65``) so ``db_antigen_gene`` strictly holds host
    or tumor gene symbols while ``db_protein_canonical`` holds the protein or
    antigen name.
    """
    if pd.isna(antigen) and pd.isna(canonical_protein):
        return None
    raw = "" if pd.isna(antigen) else str(antigen).strip()
    clean = _INLINE_BRACKET_ORG_RE.sub("", raw).strip() if raw else ""
    org = _extract_org_hint(raw, species) if raw else canonicalize_species(species)

    # Foreign / pathogen species do not map to endogenous host gene symbols.
    if org is not None and org not in {"Homo sapiens", "Mus musculus"}:
        return None

    is_mouse = org == "Mus musculus"
    canon = (
        canonical_protein
        if (isinstance(canonical_protein, str) and canonical_protein)
        else canonicalize_antigen(raw, species=species)
    )
    if canon in _NON_GENE_CANONICAL_ANTIGENS:
        return None

    # Paralog-exact CTAG1A vs CTAG1B preservation before generic NY-ESO-1 lookup.
    if re.search(r"\bctag1a\b", clean, re.I):
        return None if is_mouse else "CTAG1A"
    if re.search(r"\bctag1b\b", clean, re.I):
        return None if is_mouse else "CTAG1B"

    # Dynamic MAGE paralogs (MAGE-A1 -> MAGEA1, MAGE-A10 -> MAGEA10, etc.).
    if isinstance(canon, str):
        m_mage = re.match(r"^MAGE-([ABC])(\d+)$", canon)
        if m_mage:
            sub, num = m_mage.group(1), m_mage.group(2)
            return f"Mage{sub.lower()}{num}" if is_mouse else f"MAGE{sub}{num}"

    if isinstance(canon, str):
        table = _CANONICAL_PROTEIN_TO_GENE_MGI if is_mouse else _CANONICAL_PROTEIN_TO_GENE_HGNC
        if canon in table:
            return table[canon]

    # When the source organism is explicitly human or mouse and the raw field is
    # already a clean bare HGNC/MGI gene symbol (e.g. SNRPG, ERBB2IP, TP53, Adpgk,
    # Trp53), preserve it verbatim.
    if org == "Homo sapiens" and clean and clean.upper() not in _NON_GENE_BARE_TOKENS:
        if _BARE_HGNC_SYMBOL_RE.match(clean):
            return clean
    if org == "Mus musculus" and clean and clean.upper() not in _NON_GENE_BARE_TOKENS:
        if _BARE_MGI_SYMBOL_RE.match(clean):
            return clean
        if _BARE_HGNC_SYMBOL_RE.match(clean):
            return clean[0].upper() + clean[1:].lower()
    return None


def resolve_antigen_genes(
    antigens: pd.Series,
    canonical_proteins: pd.Series | None = None,
    species: pd.Series | str | None = None,
) -> pd.Series:
    """Vectorized :func:`resolve_antigen_gene` over Series inputs."""
    c_iter = (
        canonical_proteins
        if canonical_proteins is not None
        else [None] * len(antigens)
    )
    if species is None or isinstance(species, str):
        return pd.Series(
            [resolve_antigen_gene(a, c, species=species) for a, c in zip(antigens, c_iter)],
            index=antigens.index,
            dtype=object,
        )
    return pd.Series(
        [resolve_antigen_gene(a, c, s) for a, c, s in zip(antigens, c_iter, species)],
        index=antigens.index,
        dtype=object,
    )


def load_vdjdb(path: str | Path, verbose: bool = True) -> pd.DataFrame:
    """
    Load VDJdb database.

    VDJdb ships two on-disk shapes with **different row semantics**:

    - ``vdjdb_full.txt`` — *paired-chain* format: one row per
      clonotype, with ``cdr3.alpha`` and ``cdr3.beta`` as separate
      columns. This is what αβ matching needs.
    - ``vdjdb.txt`` — *long/slim* format: one row per **chain**, keyed
      by ``complex.id`` and ``gene`` (``TRA``/``TRB``). Output rows are
      restricted to ``gene == TRB`` so ``cdr3_beta`` only carries actual
      β-CDR3s; when ``complex.id > 0`` links a ``TRB`` row to a paired
      ``TRA`` row, ``cdr3_alpha`` and ``v_alpha`` are reconstructed onto
      that row (unpaired ``complex.id == 0`` rows keep ``cdr3_alpha = pd.NA``).

    When given a directory, looks for canonical filenames in priority
    order (``vdjdb_full.txt``, ``vdjdb.txt``, ``vdjdb.slim.txt``, or a
    date-stamped ``vdjdb-YYYY-MM-DD.txt``). Other sidecar files in the
    directory (``.meta.txt``, ``.scored.txt``, etc.) are ignored.

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
        ``mhc_allele``, ``database``, ``is_viral``).
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
        n_paired = int(df["cdr3_alpha"].notna().sum()) if "cdr3_alpha" in df.columns else 0
        if verbose:
            logger.info(
                f"  Long/slim VDJdb format: kept {len(df):,} TRB rows "
                f"({n_paired:,} paired αβ via complex.id) out of {n_total:,} total rows."
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
        "v.alpha": "v_alpha",
        "v.beta": "v_beta",
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

    The long format has one row per chain. We filter output rows to
    ``gene == TRB`` so ``cdr3_beta`` only contains actual β-CDR3s
    (α-CDR3s share the same ``cdr3`` column under ``gene == TRA`` rows
    and must never leak into ``cdr3_beta``). When ``complex.id > 0``
    links a ``TRB`` row to a paired ``TRA`` row in the same file, we
    join the ``TRA`` row's ``cdr3`` and ``v.segm`` onto ``cdr3_alpha``
    and ``v_alpha``; unpaired rows (``complex.id == 0``) keep
    ``cdr3_alpha = pd.NA``.
    """
    tra_cdr3_by_complex: dict[int, object] = {}
    tra_v_by_complex: dict[int, object] = {}
    if "complex.id" in df.columns:
        cid_num = pd.to_numeric(df["complex.id"], errors="coerce")
        has_cdr3 = (
            df["cdr3"].notna() & (df["cdr3"].astype(str).str.strip() != "")
            if "cdr3" in df.columns
            else pd.Series(False, index=df.index)
        )
        tra_rows = df[(df["gene"] == "TRA") & cid_num.gt(0) & has_cdr3].copy()
        if not tra_rows.empty:
            tra_rows["_cid_int"] = cid_num.loc[tra_rows.index].astype(int)
            tra_rows = tra_rows.drop_duplicates("_cid_int", keep="first")
            tra_cdr3_by_complex = dict(zip(tra_rows["_cid_int"], tra_rows["cdr3"]))
            if "v.segm" in tra_rows.columns:
                tra_v_by_complex = dict(zip(tra_rows["_cid_int"], tra_rows["v.segm"]))

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
        "v.segm": "v_beta",  # long format calls the V-gene column 'v.segm'
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

    if tra_cdr3_by_complex and "complex.id" in df.columns:
        trb_cid = pd.to_numeric(df["complex.id"], errors="coerce")
        valid_cid = trb_cid.where(trb_cid.gt(0))
        df["cdr3_alpha"] = valid_cid.map(
            lambda x: tra_cdr3_by_complex.get(int(x), pd.NA) if pd.notna(x) else pd.NA
        )
        df["v_alpha"] = valid_cid.map(
            lambda x: tra_v_by_complex.get(int(x), pd.NA) if pd.notna(x) else pd.NA
        )
    else:
        df["cdr3_alpha"] = pd.NA
        df["v_alpha"] = pd.NA
    return df


# Canonical VDJdb data filenames in resolution priority order. The
# release zip drops a swarm of processed sidecars (``vdjdb.meta.txt``,
# ``vdjdb.scored.txt``, ``vdjdb_full_*_broken.txt`` etc.) alongside
# these, but they have different schemas / are incomplete. We only
# look for the names below (plus date-stamped ``vdjdb-YYYY-MM-DD.txt``
# archives when none of the canonical names exist).
_VDJDB_CANONICAL_FILENAMES = ("vdjdb_full.txt", "vdjdb.txt", "vdjdb.slim.txt")
_VDJDB_SIDECAR_SUFFIXES = (".meta.txt", ".scored.txt", "_broken.txt", ".meta.tsv", ".scored.tsv")


def _pick_vdjdb_file(dir_path: Path) -> Path:
    """Choose the canonical VDJdb data file from an extracted release dir."""
    for name in _VDJDB_CANONICAL_FILENAMES:
        candidate = dir_path / name
        if candidate.is_file():
            return candidate

    dated_candidates = sorted(
        p
        for p in (*dir_path.glob("vdjdb-*.txt"), *dir_path.glob("vdjdb-*.tsv"))
        if p.is_file() and not p.name.lower().endswith(_VDJDB_SIDECAR_SUFFIXES)
    )
    if dated_candidates:
        return dated_candidates[-1]

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


def _infer_iedb_host_species(
    host_name: pd.Series, host_iri: pd.Series, mhc_allele: pd.Series
) -> pd.Series:
    """Infer TCR host species on IEDB/CEDAR v3 exports from Host Name, Host IRI, or MHC."""
    known_hosts = {"Homo sapiens", "Mus musculus", "Macaca mulatta", "Rattus norvegicus", "Pan troglodytes"}
    out = pd.Series([pd.NA] * len(host_name), index=host_name.index, dtype=object)
    for idx in out.index:
        hn = host_name.at[idx]
        if pd.notna(hn) and str(hn).strip():
            canon = canonicalize_species(hn)
            if canon in known_hosts:
                out.at[idx] = canon
                continue
        iri = str(host_iri.at[idx]).strip() if pd.notna(host_iri.at[idx]) else ""
        if iri:
            if "NCBITaxon_9606" in iri:
                out.at[idx] = "Homo sapiens"
                continue
            if "NCBITaxon_10090" in iri:
                out.at[idx] = "Mus musculus"
                continue
            if "NCBITaxon_9544" in iri:
                out.at[idx] = "Macaca mulatta"
                continue
            if "NCBITaxon_10116" in iri:
                out.at[idx] = "Rattus norvegicus"
                continue
        mhc = str(mhc_allele.at[idx]).strip().upper() if pd.notna(mhc_allele.at[idx]) else ""
        if mhc.startswith("HLA-"):
            out.at[idx] = "Homo sapiens"
            continue
        if mhc.startswith(("H2-", "H-2")):
            out.at[idx] = "Mus musculus"
            continue
        if mhc.startswith("MAMU-"):
            out.at[idx] = "Macaca mulatta"
            continue
        # IEDB inbred mouse strain IRIs (e.g. ONTIE_0000462 for C57BL/6) when MHC is absent.
        if "/ONTIE_" in iri:
            out.at[idx] = "Mus musculus"
            continue
        if pd.notna(hn) and str(hn).strip():
            out.at[idx] = canonicalize_species(hn)
    return out


_VALID_CDR3_ANCHOR_ENDS = frozenset({"F", "W", "C", "L", "Y"})


def _upgrade_iedb_cdr3_with_junction(
    cdr3: pd.Series, junction: pd.Series
) -> pd.Series:
    """Restore C...F/W anchors on IEDB v3 CDR3s when Curated/Calculated stripped them.

    Over half of IEDB v3 entries store anchorless core CDR3s (e.g. ``ASSIRSSYEQY``)
    in ``CDR3 Curated`` / ``CDR3 Calculated`` while storing the IMGT-anchored
    sequence (``CASSIRSSYEQYF``) in ``Junction Calculated``; CellRanger and VDJdb
    always emit IMGT-anchored ``C...F/W`` CDR3s. Upgrade to ``Junction Calculated``
    only when it is a valid IMGT junction (length 7..30, starts with ``C``, ends
    with a valid J-anchor residue ``F/W/C/L/Y``, and adds at most the two N/C
    anchor residues over a non-empty core ``cdr3`` so trailing J-framework
    ``FGXG`` tails are never grafted).
    """
    out = cdr3.astype(object).copy()
    for idx in out.index:
        c_val = out.at[idx]
        j_val = junction.at[idx] if idx in junction.index else pd.NA
        c_str = str(c_val).strip().upper() if pd.notna(c_val) else ""
        j_str = str(j_val).strip().upper() if pd.notna(j_val) else ""
        if c_str:
            out.at[idx] = c_str
        if (
            not j_str
            or not (7 <= len(j_str) <= 30)
            or not j_str.startswith("C")
            or j_str[-1] not in _VALID_CDR3_ANCHOR_ENDS
        ):
            continue
        if not c_str:
            out.at[idx] = j_str
        elif not c_str.startswith("C") and c_str in j_str and len(j_str) <= len(c_str) + 2:
            out.at[idx] = j_str
    return out


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

    cdr3_alpha_raw = _pick("Chain 1", "CDR3 Curated").combine_first(
        _pick("Chain 1", "CDR3 Calculated")
    )
    cdr3_beta_raw = _pick("Chain 2", "CDR3 Curated").combine_first(
        _pick("Chain 2", "CDR3 Calculated")
    )
    cdr3_alpha = _upgrade_iedb_cdr3_with_junction(
        cdr3_alpha_raw, _pick("Chain 1", "Junction Calculated", "Junction Curated")
    )
    cdr3_beta = _upgrade_iedb_cdr3_with_junction(
        cdr3_beta_raw, _pick("Chain 2", "Junction Calculated", "Junction Curated")
    )
    v_alpha = _pick("Chain 1", "Curated V Gene").combine_first(
        _pick("Chain 1", "Calculated V Gene")
    )
    v_beta = _pick("Chain 2", "Curated V Gene").combine_first(
        _pick("Chain 2", "Calculated V Gene")
    )
    mhc_allele = _pick("Assay", "MHC Allele Names")
    host_name = _pick("Host", "Name").combine_first(
        _pick("Assay", "Host Name", "Host")
    )
    host_iri = _pick("Host", "IRI", "IEDB IRI").combine_first(
        _pick("Assay", "Host IRI")
    )
    host_species = _infer_iedb_host_species(
        host_name,
        host_iri,
        mhc_allele,
    )

    out = pd.DataFrame(
        {
            "cdr3_alpha": cdr3_alpha,
            "cdr3_beta": cdr3_beta,
            "v_alpha": v_alpha,
            "v_beta": v_beta,
            "epitope": _pick("Epitope", "Name"),
            "epitope_iri": _pick("Epitope", "IEDB IRI"),
            "antigen_gene": _pick("Epitope", "Source Molecule"),
            "species": _pick("Epitope", "Source Organism"),
            "host_species": host_species,
            "mhc_allele": mhc_allele,
            "reference": _pick("Reference", "IEDB IRI", "PMID"),
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
    Load the CEDAR TCR database.

    Handles two on-disk shapes, sniffed from the file's first line (not the
    extension):

    1. The **bulk** receptor export (``doc/receptor_full_v3.zip`` →
       ``tcr_full_v3.csv``, cached by ``tcrsift data download --db cedar``):
       comma-separated with the same v3 two-row hierarchical header as IEDB, so
       it is parsed by the shared :func:`_load_iedb_v3`.
    2. A **manual web-UI export** ``cedar.tsv`` — tab-separated, flat header with
       fields like ``cdr3_b_aa`` / ``epitope_sequence`` (the original supported
       shape; kept for back-compat).

    Parameters
    ----------
    path : str or Path
        Path to a CEDAR receptor file (bulk v3 CSV or manual flat TSV).

    Returns
    -------
    pd.DataFrame
        CEDAR entries with standardized columns (``cdr3_alpha``, ``cdr3_beta``,
        ``epitope``, ``antigen_gene``, ``species``, ``is_viral``, ``database``).
    """
    path = Path(path)
    logger.info(f"Loading CEDAR from {path}")

    if _looks_like_iedb_v3(path):
        # CEDAR's bulk export IS the IEDB v3 receptor schema (cancer-scoped).
        df = _load_iedb_v3(path)
    else:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        # Standardize the manual web-export columns.
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

    species_lower = canonicalize_species_labels(df["species"]).fillna("").str.lower()

    is_viral = pd.Series(False, index=df.index)
    for pattern in VIRAL_SPECIES_PATTERNS:
        is_viral |= species_lower.str.contains(pattern, na=False)

    return is_viral


def classify_category(
    species: pd.Series,
    antigen_gene: pd.Series,
    *,
    query_species: str = "human",
) -> pd.Series:
    """Classify each entry into a coarse specificity category.

    Categories (precedence, low → high): ``unknown`` < ``other`` <
    ``viral`` / ``bacterial`` < ``self`` < ``tumor_self``. ``self`` matches
    the query organism's binomial (``Homo sapiens`` for ``query_species='human'``,
    ``Mus musculus`` for ``query_species='mouse'``).
    """
    from .species import species_to_binomial

    self_binomial = species_to_binomial(query_species)
    species_canonical = canonicalize_species_labels(species)
    species_lower = species_canonical.fillna("").astype(str).str.lower()
    antigen_lower = antigen_gene.fillna("").astype(str).str.lower()

    category = pd.Series(
        [CATEGORY_UNKNOWN] * len(species),
        index=species.index,
        dtype=object,
    )

    is_known_species = species_canonical.notna() & (species_lower.str.len() > 0)
    category[is_known_species] = CATEGORY_OTHER

    is_viral = pd.Series(False, index=species.index)
    for pattern in VIRAL_SPECIES_PATTERNS:
        is_viral |= species_lower.str.contains(pattern, na=False, regex=False)
    category[is_viral] = CATEGORY_VIRAL

    is_bacterial = pd.Series(False, index=species.index)
    for pattern in BACTERIAL_SPECIES_PATTERNS:
        is_bacterial |= species_lower.str.contains(pattern, na=False, regex=False)
    category[is_bacterial & ~is_viral] = CATEGORY_BACTERIAL

    is_self = species_canonical == self_binomial
    category[is_self & ~is_viral & ~is_bacterial] = CATEGORY_SELF

    # Tumor-self override: regex patterns anchored at word boundaries so
    # short tokens (her2, tert, psa, wt1, dct) don't false-positive on
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
    verbose: bool = True,
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
    verbose : bool
        Print progress information.

    Returns
    -------
    pd.DataFrame
        Combined database with standardized columns
    """
    dfs = []

    if vdjdb_path:
        dfs.append(load_vdjdb(vdjdb_path, verbose=verbose))
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

    if verbose:
        logger.info(f"Combined database has {len(combined)} entries")
    return combined


def match_clonotypes(
    clonotypes: pd.DataFrame,
    database: pd.DataFrame,
    match_by: str = "CDR3ab",
    match_strictness: str | None = None,
    match_mode: str = "exact",
    max_distance: int = 1,
    verbose: bool = True,
    show_progress: bool = True,
    *,
    species: str = "human",
    host_species_filter: str = "all",
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
    match_mode : str, default ``"exact"``
        Distance regime for the β-only match (#57). When ``"exact"``,
        only exact-string β matches are recorded — pre-2.1 behaviour.
        When ``"levenshtein"``, clones that don't hit exactly fall
        through to a Levenshtein-≤``max_distance`` neighbor search
        against the database's β CDR3s. αβ matching stays strict-exact
        in both modes — fuzzy αβ is too noisy biologically.
    max_distance : int, default 1
        Maximum edit distance for fuzzy β matching (only used when
        ``match_mode == "levenshtein"``). Currently capped at 1 in
        the implementation; higher values fall back to 1 with a
        warning.
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
        / ``"b_only_cross"`` / ``"b_only_near"`` /
        ``"b_only_near_cross"`` / None; ``_cross`` suffix indicates a
        non-human host-species match per #83; ``_near`` suffix
        indicates a Levenshtein-≤1 β-only fuzzy match per #57),
        ``db_match_distance`` (int — edit distance to the matched DB
        β CDR3; 0 for exact, 1 for Levenshtein-1, ``None`` when
        unmatched),
        ``db_epitope``, ``db_protein``, ``db_protein_canonical``,
        ``db_species`` (canonical antigen source organism),
        ``db_host_species`` (canonical TCR donor organism, new in #83),
        ``db_mhc``, ``db_category`` (``viral`` / ``bacterial`` /
        ``self`` / ``tumor_self`` / ``other`` / ``unknown`` /
        ``contradictory`` — last when multiple matches disagree on
        an informative label per #83), ``db_category_dominant`` (most
        common informative label when unique), ``db_category_detail``
        (e.g. ``viral`` / ``contradictory`` /
        ``contradictory_dominant_viral``), and
        ``db_category_counts`` (``label:count`` audit string),
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

    # #57: validate fuzzy-match parameters before doing any DB work.
    if match_mode not in MATCH_MODES:
        raise TCRsiftValidationError(
            f"Invalid match_mode: {match_mode!r}",
            hint=f"Valid options are: {list(MATCH_MODES)}",
        )
    validate_numeric_param(max_distance, "max_distance", min_value=1)
    if match_mode == "levenshtein" and max_distance > 1:
        logger.warning(
            "match_mode='levenshtein' currently supports max_distance=1 "
            "only; clamping max_distance=%d → 1.", max_distance,
        )
        max_distance = 1

    if verbose:
        logger.info(
            f"Matching {len(clonotypes):,} clonotypes against "
            f"{len(database):,} database entries "
            f"(strictness={strictness}, mode={match_mode})"
        )

    from .species import resolve_species, species_to_binomial

    query_species = resolve_species(species, df=clonotypes, warn_mismatch=False)
    if host_species_filter not in ("all", "same"):
        raise TCRsiftValidationError(
            f"Invalid host_species_filter: {host_species_filter!r}",
            hint="Valid options are: ['all', 'same']",
        )

    df = clonotypes.copy()

    # Initialize annotation columns.
    df["db_match"] = False
    df["db_match_partial"] = False
    df["db_match_strength"] = None
    df["db_match_distance"] = None  # #57: edit distance to matched DB β
    df["db_match_score"] = None     # #56: composite confidence ∈ [0, 1]
    df["db_v_gene_agreement"] = None
    df["db_n_studies"] = None
    df["db_category_votes"] = None
    df["db_epitope"] = None
    df["db_protein"] = None
    df["db_protein_canonical"] = None
    df["db_antigen_gene"] = None
    df["db_species"] = None
    df["db_host_species"] = None
    df["db_mhc"] = None
    df["db_category"] = None
    df["db_category_dominant"] = None
    df["db_category_detail"] = None
    df["db_category_counts"] = None
    df["db_database"] = None
    df["is_viral"] = False

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

    database = _canonicalize_database_species(database)

    if host_species_filter == "same" and "host_species" in database.columns:
        target_host = species_to_binomial(query_species)
        database = database[
            database["host_species"].isna() | (database["host_species"] == target_host)
        ]
        if database.empty:
            if verbose:
                logger.info(
                    "  No database entries remain after host_species_filter='same' (%s)",
                    target_host,
                )
            return df

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
        new_columns["db_category"] = classify_category(
            species_col, antigen_col, query_species=query_species
        )
    if "db_protein_canonical" not in database.columns:
        new_columns["db_protein_canonical"] = canonicalize_antigens(
            antigen_col, species=species_col
        )
    canon_col = new_columns.get("db_protein_canonical", database.get("db_protein_canonical"))
    if "db_antigen_gene" not in database.columns:
        new_columns["db_antigen_gene"] = resolve_antigen_genes(
            antigen_col, canonical_proteins=canon_col, species=species_col
        )
    if new_columns:
        database = database.assign(**new_columns)

    beta_neighbor_index: dict[str, set[str]] | None = (
        _build_cdr3_neighbor_index(database["cdr3_beta"])
        if match_mode == "levenshtein"
        else None
    )

    # Pre-index positional row slices by CDR3β and paired (CDR3α, CDR3β) once
    # so per-clone lookups run in O(1) instead of scanning a 140k-row DataFrame.
    import numpy as np

    b_to_pos: dict[str, np.ndarray] = {}
    if strictness in ("ab_with_partial", "b_only"):
        valid_b = database["cdr3_beta"].notna() & (database["cdr3_beta"].astype(str) != "")
        if valid_b.any():
            b_sub = database.loc[valid_b, "cdr3_beta"].astype(str)
            b_positions = np.flatnonzero(valid_b.to_numpy())
            for beta_key, rel_idx in b_sub.groupby(b_sub, sort=False).indices.items():
                if len(rel_idx):
                    b_to_pos[str(beta_key)] = b_positions[rel_idx]

    if strictness in ("strict_ab", "ab_with_partial"):
        allow_b_fallback = strictness == "ab_with_partial"
        valid_ab = (
            database["cdr3_alpha"].notna()
            & (database["cdr3_alpha"].astype(str) != "")
            & database["cdr3_beta"].notna()
            & (database["cdr3_beta"].astype(str) != "")
        )
        ab_to_pos: dict[tuple[str, str], np.ndarray] = {}
        if valid_ab.any():
            ab_positions = np.flatnonzero(valid_ab.to_numpy())
            ab_sub = database.loc[valid_ab, ["cdr3_alpha", "cdr3_beta"]].astype(str)
            for (a_key, b_key), rel_idx in ab_sub.groupby(
                ["cdr3_alpha", "cdr3_beta"], sort=False
            ).indices.items():
                if len(rel_idx):
                    ab_to_pos[(str(a_key), str(b_key))] = ab_positions[rel_idx]

        row_iter = df.iterrows()
        if show_progress:
            row_iter = tqdm(
                list(df.iterrows()),
                desc="Matching clonotypes",
                unit="clone",
            )

        for idx, row in row_iter:
            raw_a = row.get("CDR3_alpha", "")
            raw_b = row.get("CDR3_beta", "")
            alpha = raw_a.strip() if isinstance(raw_a, str) else ""
            beta = raw_b.strip() if isinstance(raw_b, str) else ""

            if alpha and beta and (alpha, beta) in ab_to_pos:
                matches = database.iloc[ab_to_pos[(alpha, beta)]]
                _annotate_match(
                    df, idx, matches, strength="ab", distance=0, query_species=query_species
                )
            elif allow_b_fallback and beta and beta in b_to_pos:
                matches = database.iloc[b_to_pos[beta]]
                _annotate_match(
                    df, idx, matches, strength="b_only", partial=True, distance=0,
                    query_species=query_species,
                )
            elif allow_b_fallback and beta_neighbor_index is not None and beta:
                neighbors = _find_lev_neighbors(beta, beta_neighbor_index, max_distance)
                if neighbors:
                    best_d = neighbors[0][1]
                    pos_list = [b_to_pos[c] for c, d in neighbors if d == best_d and c in b_to_pos]
                    if pos_list:
                        pos = np.sort(np.concatenate(pos_list))
                        matches = database.iloc[pos]
                        _annotate_match(
                            df, idx, matches,
                            strength="b_only_near", partial=True, distance=best_d,
                            query_species=query_species,
                        )

    else:  # strictness == "b_only"
        row_iter = df.iterrows()
        if show_progress:
            row_iter = tqdm(
                list(df.iterrows()),
                desc="Matching clonotypes",
                unit="clone",
            )

        for idx, row in row_iter:
            raw_b = row.get("CDR3_beta", "")
            beta = raw_b.strip() if isinstance(raw_b, str) else ""
            if beta and beta in b_to_pos:
                matches = database.iloc[b_to_pos[beta]]
                _annotate_match(
                    df, idx, matches, strength="b_only", distance=0,
                    query_species=query_species,
                )
            elif beta_neighbor_index is not None and beta:
                neighbors = _find_lev_neighbors(beta, beta_neighbor_index, max_distance)
                if neighbors:
                    best_d = neighbors[0][1]
                    pos_list = [b_to_pos[c] for c, d in neighbors if d == best_d and c in b_to_pos]
                    if pos_list:
                        pos = np.sort(np.concatenate(pos_list))
                        matches = database.iloc[pos]
                        _annotate_match(
                            df, idx, matches,
                            strength="b_only_near", distance=best_d,
                            query_species=query_species,
                        )

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
    distance: int | None = None,
    *,
    query_species: str = "human",
):
    """Annotate a single clonotype with match information."""
    from .species import species_to_binomial

    if len(matches) == 0:
        return

    query_binomial = species_to_binomial(query_species)

    antigen_matches = matches
    for group_col in ("db_protein_canonical", "antigen_gene"):
        if group_col not in matches.columns:
            continue
        group_values = matches[group_col].dropna()
        if len(group_values) == 0:
            continue
        picked = group_values.mode().iloc[0]
        antigen_matches = matches[matches[group_col] == picked]
        break

    if (
        "db_protein_canonical" in antigen_matches.columns
        and antigen_matches["db_protein_canonical"].nunique(dropna=True) > 1
    ):
        canonicals = antigen_matches["db_protein_canonical"].dropna()
        if len(canonicals) > 0:
            picked_canonical = canonicals.mode().iloc[0]
            antigen_matches = antigen_matches[
                antigen_matches["db_protein_canonical"] == picked_canonical
            ]

    if "host_species" in antigen_matches.columns:
        hs = canonicalize_species_labels(antigen_matches["host_species"]).dropna()
        if len(hs) > 0:
            picked_hs = hs.mode().iloc[0]
            df.loc[idx, "db_host_species"] = picked_hs
            if picked_hs != query_binomial:
                strength = f"{strength}_cross"

    df.loc[idx, "db_match"] = True
    df.loc[idx, "db_match_strength"] = strength
    if distance is not None:
        df.loc[idx, "db_match_distance"] = distance

    epitopes = antigen_matches["epitope"].dropna()
    if len(epitopes) > 0:
        df.loc[idx, "db_epitope"] = epitopes.mode().iloc[0]

    if "species" in antigen_matches.columns:
        species = canonicalize_species_labels(antigen_matches["species"]).dropna()
        if len(species) > 0:
            df.loc[idx, "db_species"] = species.mode().iloc[0]

    if "antigen_gene" in antigen_matches.columns:
        proteins = antigen_matches["antigen_gene"].dropna()
        if len(proteins) > 0:
            df.loc[idx, "db_protein"] = proteins.mode().iloc[0]

    if "db_protein_canonical" in antigen_matches.columns:
        canonicals = antigen_matches["db_protein_canonical"].dropna()
        if len(canonicals) > 0:
            df.loc[idx, "db_protein_canonical"] = canonicals.mode().iloc[0]

    if "db_antigen_gene" in antigen_matches.columns:
        ag_genes = antigen_matches["db_antigen_gene"].dropna()
        if len(ag_genes) > 0:
            df.loc[idx, "db_antigen_gene"] = ag_genes.mode().iloc[0]
    else:
        df.loc[idx, "db_antigen_gene"] = resolve_antigen_gene(
            df.loc[idx, "db_protein"],
            canonical_protein=df.loc[idx, "db_protein_canonical"],
            species=df.loc[idx, "db_species"],
        )

    if "mhc_allele" in antigen_matches.columns:
        mhcs = antigen_matches["mhc_allele"].dropna()
        if len(mhcs) > 0:
            df.loc[idx, "db_mhc"] = mhcs.mode().iloc[0]

    if "db_category" in matches.columns:
        category, dominant, detail, counts = _summarize_match_categories(
            matches["db_category"]
        )
        df.loc[idx, "db_category"] = category
        df.loc[idx, "db_category_dominant"] = dominant
        df.loc[idx, "db_category_detail"] = detail
        df.loc[idx, "db_category_counts"] = counts

    df.loc[idx, "db_database"] = ";".join(matches["database"].unique())
    df.loc[idx, "is_viral"] = matches["is_viral"].any()

    if partial:
        df.loc[idx, "db_match_partial"] = True

    # #56: structured match-confidence signals + composite score. Each
    # signal is None when the underlying data isn't available so
    # downstream code can recognise "couldn't be checked" vs "checked
    # and found no agreement".
    if "v_beta" in matches.columns:
        clone_v_beta = _canonical_v_gene(df.loc[idx].get("beta_v_gene"))
        if clone_v_beta:
            db_v_betas = matches["v_beta"].dropna().map(_canonical_v_gene)
            db_v_betas = db_v_betas[db_v_betas != ""]
            if len(db_v_betas) > 0:
                df.loc[idx, "db_v_gene_agreement"] = float(
                    (db_v_betas == clone_v_beta).mean()
                )

    if "reference" in matches.columns:
        refs = matches["reference"].dropna().astype(str)
        refs = refs[refs.str.strip() != ""]
        df.loc[idx, "db_n_studies"] = int(refs.nunique())

    if "db_category" in matches.columns:
        cat_counts = matches["db_category"].dropna().value_counts()
        if len(cat_counts) > 0:
            total = int(cat_counts.sum())
            df.loc[idx, "db_category_votes"] = ";".join(
                f"{cat}:{cnt / total:.2f}"
                for cat, cnt in cat_counts.items()
            )

    df.loc[idx, "db_match_score"] = _compute_match_score(
        strength=strength,
        distance=distance,
        n_studies=df.loc[idx].get("db_n_studies"),
        v_gene_agreement=df.loc[idx].get("db_v_gene_agreement"),
    )


# Base scores per match strength label. ``_cross`` and ``_near`` suffixes
# get multiplicative penalties, applied in :func:`_compute_match_score`.
_BASE_STRENGTH_SCORE = {
    "ab": 1.0,
    "b_only": 0.7,
}
_NEAR_PENALTY = 0.7   # b_only_near = b_only * 0.7
_CROSS_PENALTY = 0.5  # *_cross = base * 0.5


def _compute_match_score(
    strength: str,
    distance: int | None,
    n_studies: object,
    v_gene_agreement: object,
) -> float:
    """Combine the structured-confidence signals into a single
    ``db_match_score`` ∈ [0, 1] (#56).

    Formula:

    1. Start from ``base = _BASE_STRENGTH_SCORE[strength_stem]`` —
       1.0 for αβ exact, 0.7 for β-only.
    2. Apply ``_NEAR_PENALTY`` (0.7×) if the strength carries the
       ``_near`` suffix (Lev-1 fuzzy β match, #57).
    3. Apply ``_CROSS_PENALTY`` (0.5×) if the strength carries the
       ``_cross`` suffix (non-human host species, #83 / #54).
    4. Reduce by up to 30% for V-gene disagreement when measurable:
       ``factor = 1 - 0.3 * (1 - v_gene_agreement)``. ``None`` (unknown
       — DB or clone missing V-gene info) is treated as 1.0 (neutral).
    5. Reduce for single-study evidence (Britanova / Emerson style
       repeats-from-one-paper concern from #56):
       ``factor = 1 - 0.3 / max(n_studies, 1)``. So 1 study → 0.7,
       2 → 0.85, 3+ → 0.9 and approaching 1.0. ``None`` and 0 are
       treated as 1 study.

    Returns a float rounded to 3 decimals so the score is stable
    across pandas/numpy float representations and reads cleanly in
    CSV output.
    """
    stem = strength
    # Suffixes can stack as ``..._near_cross``; peel ``_cross`` first
    # because it's always outermost (added in ``_annotate_match`` after
    # ``_near`` is already on the strength label).
    has_cross = stem.endswith("_cross")
    if has_cross:
        stem = stem[: -len("_cross")]
    has_near = stem.endswith("_near")
    if has_near:
        stem = stem[: -len("_near")]
    base = _BASE_STRENGTH_SCORE.get(stem, 0.5)
    if has_near:
        base *= _NEAR_PENALTY
    if has_cross:
        base *= _CROSS_PENALTY

    if isinstance(v_gene_agreement, (int, float)) and not pd.isna(v_gene_agreement):
        v_factor = 1.0 - 0.3 * (1.0 - float(v_gene_agreement))
    else:
        v_factor = 1.0

    try:
        n = int(n_studies) if n_studies is not None and not pd.isna(n_studies) else 0
    except (TypeError, ValueError):
        n = 0
    n = max(n, 1)
    studies_factor = 1.0 - 0.3 / n

    score = base * v_factor * studies_factor
    return round(min(1.0, max(0.0, score)), 3)


def annotate_clonotypes(
    clonotypes: pd.DataFrame,
    vdjdb_path: str | Path | None = None,
    iedb_path: str | Path | None = None,
    cedar_path: str | Path | None = None,
    iedb_epitope_path: str | Path | None = None,
    match_by: str = "CDR3ab",
    match_strictness: str | None = None,
    match_mode: str = "exact",
    max_distance: int = 1,
    exclude_viral: bool = False,
    flag_only: bool = False,
    *,
    verbose: bool = True,
    show_progress: bool = True,
    add_publicness: bool = False,
    publicness_cdr3_col: str = "CDR3_beta",
    publicness_v_gene_col: str = "beta_v_gene",
    publicness_j_gene_col: str = "beta_j_gene",
    species: str = "human",
    host_species_filter: str = "all",
) -> pd.DataFrame:
    """Main annotation function."""
    from .species import resolve_species

    resolved_species = resolve_species(species, df=clonotypes, warn_mismatch=False)
    _DEFAULT_ANNOTATION_COLUMNS = {
        "db_match": False,
        "db_match_partial": False,
        "db_match_strength": None,
        "db_match_distance": None,
        "db_match_score": None,
        "db_v_gene_agreement": None,
        "db_n_studies": None,
        "db_category_votes": None,
        "db_epitope": None,
        "db_protein": None,
        "db_protein_canonical": None,
        "db_antigen_gene": None,
        "db_species": None,
        "db_host_species": None,
        "db_mhc": None,
        "db_category": None,
        "db_category_dominant": None,
        "db_category_detail": None,
        "db_category_counts": None,
        "db_database": None,
        "is_viral": False,
    }

    if not any([vdjdb_path, iedb_path, cedar_path]):
        if verbose:
            logger.info("No annotation database paths provided; returning input with empty annotations")
        df = clonotypes.copy()
        if "CDR3ab" not in df.columns and {"CDR3_alpha", "CDR3_beta"} <= set(df.columns):
            df["CDR3ab"] = df["CDR3_alpha"].fillna("") + "_" + df["CDR3_beta"].fillna("")
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
                species=resolved_species,
            )
        return df

    database = load_databases(
        vdjdb_path=vdjdb_path,
        iedb_path=iedb_path,
        cedar_path=cedar_path,
        iedb_epitope_path=iedb_epitope_path,
        verbose=verbose,
    )

    df = match_clonotypes(
        clonotypes,
        database,
        match_by=match_by,
        match_strictness=match_strictness,
        match_mode=match_mode,
        max_distance=max_distance,
        verbose=verbose,
        show_progress=show_progress,
        species=resolved_species,
        host_species_filter=host_species_filter,
    )

    if exclude_viral and not flag_only:
        initial = len(df)
        df = df[~df["is_viral"]]
        if verbose:
            logger.info(f"Excluded {initial - len(df)} viral clones")

    if add_publicness:
        df = _apply_publicness(
            df,
            cdr3_col=publicness_cdr3_col,
            v_gene_col=publicness_v_gene_col,
            j_gene_col=publicness_j_gene_col,
            log_summary=verbose,
            species=resolved_species,
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
