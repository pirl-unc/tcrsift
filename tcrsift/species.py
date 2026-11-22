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

"""Centralized species normalization and auto-detection for TCRsift.

Supports ``"human"`` (*Homo sapiens*) and ``"mouse"`` (*Mus musculus*), plus
``"auto"`` detection from single-cell AnnData gene axes or V(D)J gene calls.
"""

from __future__ import annotations

import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)

VALID_SPECIES: tuple[str, ...] = ("human", "mouse")
VALID_SPECIES_WITH_AUTO: tuple[str, ...] = ("human", "mouse", "auto")

_SPECIES_ALIASES: dict[str, str] = {
    "human": "human",
    "humans": "human",
    "hs": "human",
    "hsapiens": "human",
    "homosapiens": "human",
    "homosapienshuman": "human",
    "grch38": "human",
    "grch37": "human",
    "hg38": "human",
    "hg19": "human",
    "mouse": "mouse",
    "mice": "mouse",
    "murine": "mouse",
    "mm": "mouse",
    "mmusculus": "mouse",
    "musmusculus": "mouse",
    "musmusculusmouse": "mouse",
    "grcm39": "mouse",
    "grcm38": "mouse",
    "mm10": "mouse",
    "mm39": "mouse",
}

_BINOMIAL_BY_SPECIES: dict[str, str] = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
}

_PYENSEMBL_BY_SPECIES: dict[str, str] = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
}

# Unambiguous V-gene subfamilies that appear in only one of human vs mouse IMGT.
_MOUSE_ONLY_V_PREFIXES = (
    "TRAV14D",
    "TRAV14N",
    "TRAV6D",
    "TRAV6N",
    "TRAV7D",
    "TRAV7N",
    "TRAV9D",
    "TRAV9N",
    "TRAV12D",
    "TRAV12N",
    "TRAV13D",
    "TRAV13N",
    "TRAV15",
    "TRAV16",
    "TRBV12-1",
    "TRBV12-2",
    "TRBV13-1",
    "TRBV13-2",
    "TRBV13-3",
    "TRBV14",
    "TRBV15",
    "TRBV16",
    "TRBV17",
    "TRBV26",
    "TRBV31",
)
_HUMAN_ONLY_V_PREFIXES = (
    "TRAV12-2",
    "TRAV38-2",
    "TRAV26-1",
    "TRAV26-2",
    "TRBV20-1",
    "TRBV28",
    "TRBV6-5",
    "TRBV7-9",
)

_SYMBOL_COLS = ("gene_symbols", "feature_name", "symbol", "gene_symbol", "gene_name")


def normalize_species(
    value: object,
    *,
    allow_auto: bool = False,
    default: str = "human",
) -> str:
    """Normalize a user/config species token to ``'human'`` or ``'mouse'``
    (or ``'auto'`` when ``allow_auto=True``).
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return default
    raw = str(value).strip()
    if not raw:
        return default
    key = re.sub(r"[^a-z0-9]+", "", raw.lower())
    if allow_auto and key == "auto":
        return "auto"
    if key in _SPECIES_ALIASES:
        return _SPECIES_ALIASES[key]
    valid = VALID_SPECIES_WITH_AUTO if allow_auto else VALID_SPECIES
    raise ValueError(
        f"Unsupported species {value!r}; expected one of {list(valid)} "
        "(or 'Homo sapiens' / 'Mus musculus')."
    )


def species_to_binomial(species: str | None) -> str:
    """Return ``'Homo sapiens'`` or ``'Mus musculus'`` for a normalized species."""
    norm = normalize_species(species)
    return _BINOMIAL_BY_SPECIES[norm]


def species_to_pyensembl(species: str | None) -> str:
    """Return pyensembl species identifier (``'homo_sapiens'`` or ``'mus_musculus'``)."""
    norm = normalize_species(species)
    return _PYENSEMBL_BY_SPECIES[norm]


def detect_vdj_species(df: pd.DataFrame | None) -> str | None:
    """Sniff ``'human'`` or ``'mouse'`` from V/J/C gene columns when unambiguous."""
    if df is None or len(df) == 0:
        return None
    cols = [
        c
        for c in (
            "alpha_v_gene",
            "beta_v_gene",
            "alpha_c_gene",
            "beta_c_gene",
            "TRA_1_v_gene",
            "TRB_1_v_gene",
            "TRA_1_c_gene",
            "TRB_1_c_gene",
            "v_gene",
            "c_gene",
        )
        if c in df.columns
    ]
    if not cols:
        return None

    tokens: list[str] = []
    for c in cols:
        vals = df[c].dropna().astype(str).str.strip()
        tokens.extend(v for v in vals if v and v.lower() != "nan")
    if not tokens:
        return None

    # 1. Title-case Cell Ranger GRCm39 calls (Trav*, Trbv*, Trac, Trbc1/2).
    n_title = sum(1 for t in tokens if t.startswith(("Trav", "Trbv", "Trac", "Trbc", "Traj", "Trbj")))
    n_upper = sum(1 for t in tokens if t.startswith(("TRAV", "TRBV", "TRAC", "TRBC", "TRAJ", "TRBJ")))
    if n_title > 0 and n_title >= n_upper:
        return "mouse"

    # 2. Species-unique IMGT V-gene subfamilies (works even when uppercased).
    upper_bases = [t.upper().split("*", 1)[0] for t in tokens]
    n_mouse_v = sum(1 for u in upper_bases if u.startswith(_MOUSE_ONLY_V_PREFIXES) or "D-" in u)
    n_human_v = sum(1 for u in upper_bases if u.startswith(_HUMAN_ONLY_V_PREFIXES))
    if n_mouse_v > 0 and n_human_v == 0:
        return "mouse"
    if n_human_v > 0 and n_mouse_v == 0:
        return "human"
    if n_upper > 0 and n_title == 0 and n_mouse_v == 0:
        return "human"
    return None


def detect_adata_species(adata) -> str | None:
    """Sniff ``'human'`` or ``'mouse'`` from an AnnData's gene axis and obs."""
    if adata is None or getattr(adata, "n_vars", 0) == 0:
        return None

    var_names = pd.Series(adata.var_names.astype(str))
    if (var_names.str.startswith("ENSMUSG")).mean() > 0.2:
        return "mouse"
    if (var_names.str.startswith("ENSG")).mean() > 0.2:
        return "human"

    symbols = None
    var = getattr(adata, "var", None)
    if var is not None:
        for col in _SYMBOL_COLS:
            if col in getattr(var, "columns", []):
                symbols = pd.Series(var[col].astype(str))
                break
    if symbols is None:
        symbols = var_names

    sym_set = set(symbols)
    sym_upper = {s.upper() for s in sym_set}

    # Unmistakable species-specific symbols / mitochondrial case conventions.
    has_lower_mt = any(s.startswith("mt-") for s in sym_set)
    has_upper_mt = any(s.startswith("MT-") for s in sym_set)
    if has_lower_mt and not has_upper_mt:
        return "mouse"
    if "CD8B1" in sym_upper or "LY6C2" in sym_upper or "KLRB1C" in sym_upper or "LY6G" in sym_upper:
        return "mouse"
    if "GNLY" in sym_upper or "GZMH" in sym_upper or "KLRF1" in sym_upper:
        return "human"
    if has_upper_mt and not has_lower_mt and "CD8B1" not in sym_upper:
        # Check MGI TitleCase vs HGNC UPPERCASE among alphabetic symbols.
        alpha_syms = [s for s in sym_set if len(s) >= 3 and s[0].isalpha() and not s.upper().startswith("MT-")]
        if alpha_syms:
            title_frac = sum(1 for s in alpha_syms if s[0].isupper() and any(ch.islower() for ch in s[1:])) / len(alpha_syms)
            if title_frac > 0.5:
                return "mouse"
        return "human"

    alpha_syms = [s for s in sym_set if len(s) >= 3 and s[0].isalpha()]
    if len(alpha_syms) >= 5:
        title_frac = sum(
            1 for s in alpha_syms
            if s[0].isupper() and any(ch.islower() for ch in s[1:]) and not s.startswith(("HLA-", "LINC"))
        ) / len(alpha_syms)
        upper_frac = sum(1 for s in alpha_syms if s == s.upper()) / len(alpha_syms)
        if title_frac > 0.5 and upper_frac < 0.2:
            return "mouse"
        if upper_frac > 0.7 and title_frac < 0.1:
            return "human"

    obs = getattr(adata, "obs", None)
    if obs is not None:
        return detect_vdj_species(obs)
    return None


def resolve_species(
    species: object = "human",
    *,
    adata=None,
    df: pd.DataFrame | None = None,
    warn_mismatch: bool = True,
    context: str = "",
) -> str:
    """Resolve ``species`` (handling ``'auto'``) and warn on obvious data mismatches."""
    norm = normalize_species(species, allow_auto=True)
    detected = detect_adata_species(adata) if adata is not None else None
    if detected is None and df is not None:
        detected = detect_vdj_species(df)

    if norm == "auto":
        resolved = detected or "human"
        if detected:
            prefix = f"{context}: " if context else ""
            logger.info("%sAuto-detected species=%r", prefix, resolved)
        return resolved

    if warn_mismatch and detected is not None and detected != norm:
        prefix = f"{context}: " if context else ""
        logger.warning(
            "%sInput data carries %s identifiers, while species=%r is configured. "
            "Pass --species %s (or set species: %s in config/sample sheet) if this is %s data.",
            prefix,
            detected,
            norm,
            detected,
            detected,
            detected,
        )
    return norm
