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
Centralized gene definitions and lookup utilities.

This module provides a single source of truth for gene identifiers used throughout
TCRsift, with version-robust ENSEMBL ID matching.

ENSEMBL IDs can have version suffixes (e.g., ENSG00000167286.10) that change between
genome builds. The utilities here handle both versioned and unversioned IDs.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Gene:
    """A gene with its identifiers."""

    symbol: str
    ensembl_id: str
    description: str = ""

    def matches_ensembl(self, query: str) -> bool:
        """Check if query matches this gene's ENSEMBL ID (version-robust).

        Handles both exact matches and versioned IDs:
        - "ENSG00000167286" matches "ENSG00000167286"
        - "ENSG00000167286" matches "ENSG00000167286.10"
        - "ENSG00000167286.10" matches "ENSG00000167286"
        - "ENSG00000167286.10" matches "ENSG00000167286.10"
        """
        query_base = query.split(".")[0]
        self_base = self.ensembl_id.split(".")[0]
        return query_base == self_base


# =============================================================================
# T-cell marker genes
# =============================================================================

CD3D = Gene("CD3D", "ENSG00000167286", "CD3 delta subunit")
CD3E = Gene("CD3E", "ENSG00000198851", "CD3 epsilon subunit")
CD3G = Gene("CD3G", "ENSG00000160654", "CD3 gamma subunit")
CD4 = Gene("CD4", "ENSG00000010610", "CD4 co-receptor")
CD8A = Gene("CD8A", "ENSG00000153563", "CD8 alpha subunit")
CD8B = Gene("CD8B", "ENSG00000172116", "CD8 beta subunit")

# Grouped markers
TCELL_MARKERS = [CD3D, CD3E, CD3G, CD4, CD8A, CD8B]
CD3_GENES = [CD3D, CD3E, CD3G]
CD8_GENES = [CD8A, CD8B]


# =============================================================================
# Lookup utilities
# =============================================================================


def canonicalize_gene(name: str | None, *, keep_allele: bool = False) -> str | None:
    """Canonicalize a TCR V/J gene name to a single IMGT-style token.

    The shared canonicalizer for *all* V/J gene matching (Pgen/Ppost
    marginals, allele mapping, repertoire joins) so the same string never
    silently fails to match itself across formats. Handles:

    - **Allele suffix** — dropped by default (gene-level), e.g.
      ``TRBV20-1*01`` → ``TRBV20-1`` (``keep_allele=True`` keeps ``*01``).
    - **Adaptive prefix** — ``TCRBV20-01`` → ``TRBV20-1`` (``TCR[ABDG]`` → ``TR[ABDG]``).
    - **Leading zeros** — ``TRBV09-01`` → ``TRBV9-1``, ``TRAJ08`` → ``TRAJ8``.
    - **DV/OR slash fusion** — ``TRAV14/DV4`` → ``TRAV14DV4`` (slash removed so
      the two interchangeable spellings collapse to one token).

    Returns ``None`` for empty / non-string input.
    """
    if not isinstance(name, str):
        return None
    g = name.strip().upper()
    if not g:
        return None
    allele = ""
    if "*" in g:
        g, _, allele = g.partition("*")
    g = re.sub(r"^TCR([ABDG])", r"TR\1", g)        # Adaptive TCRBV → TRBV
    g = g.replace("/", "")                          # TRAV14/DV4 → TRAV14DV4
    g = re.sub(r"(?<=[A-Z-])0+(\d)", r"\1", g)      # strip leading zeros
    if keep_allele and allele:
        # normalize allele to two digits (Adaptive '1' → '01')
        return f"{g}*{int(allele):02d}" if allele.isdigit() else f"{g}*{allele}"
    return g


def strip_ensembl_version(ensembl_id: str) -> str:
    """Remove version suffix from ENSEMBL ID.

    Examples:
        >>> strip_ensembl_version("ENSG00000167286.10")
        'ENSG00000167286'
        >>> strip_ensembl_version("ENSG00000167286")
        'ENSG00000167286'
    """
    return ensembl_id.split(".")[0]


def find_gene_by_symbol(symbol: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by its symbol.

    Args:
        symbol: Gene symbol (e.g., "CD3D")
        genes: List of genes to search (default: TCELL_MARKERS)

    Returns:
        Gene if found, None otherwise
    """
    if genes is None:
        genes = TCELL_MARKERS
    symbol_upper = symbol.upper()
    for gene in genes:
        if gene.symbol.upper() == symbol_upper:
            return gene
    return None


def find_gene_by_ensembl(ensembl_id: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by its ENSEMBL ID (version-robust).

    Args:
        ensembl_id: ENSEMBL ID with or without version (e.g., "ENSG00000167286.10")
        genes: List of genes to search (default: TCELL_MARKERS)

    Returns:
        Gene if found, None otherwise
    """
    if genes is None:
        genes = TCELL_MARKERS
    for gene in genes:
        if gene.matches_ensembl(ensembl_id):
            return gene
    return None


def find_gene(identifier: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by symbol or ENSEMBL ID.

    Args:
        identifier: Gene symbol or ENSEMBL ID
        genes: List of genes to search (default: TCELL_MARKERS)

    Returns:
        Gene if found, None otherwise
    """
    if genes is None:
        genes = TCELL_MARKERS
    # Try symbol first
    gene = find_gene_by_symbol(identifier, genes)
    if gene:
        return gene
    # Try ENSEMBL ID
    return find_gene_by_ensembl(identifier, genes)


def find_column_for_gene(
    gene: Gene,
    column_names: Sequence[str],
) -> str | None:
    """Find a column name that matches a gene (by symbol or ENSEMBL ID).

    Handles various naming conventions:
    - Direct symbol match: "CD3D"
    - ENSEMBL ID (with or without version): "ENSG00000167286", "ENSG00000167286.10"

    Args:
        gene: Gene to find
        column_names: List of column names to search

    Returns:
        Matching column name if found, None otherwise
    """
    # Try exact symbol match first (most common case)
    if gene.symbol in column_names:
        return gene.symbol

    # Try case-insensitive symbol match
    symbol_lower = gene.symbol.lower()
    for col in column_names:
        if col.lower() == symbol_lower:
            return col

    # Try ENSEMBL ID match (version-robust)
    ensembl_base = strip_ensembl_version(gene.ensembl_id)
    for col in column_names:
        col_base = strip_ensembl_version(col)
        if col_base == ensembl_base:
            return col

    return None


def build_versioned_rename_map(
    column_names: Sequence[str],
    genes: Sequence[Gene] | None = None,
) -> dict[str, str]:
    """Build a rename map from versioned ENSEMBL IDs to gene symbols.

    Useful for renaming DataFrame columns from ENSEMBL IDs to symbols.

    Args:
        column_names: List of column names (may include versioned ENSEMBL IDs)
        genes: List of genes to match against (default: TCELL_MARKERS)

    Returns:
        Dict mapping original column names to gene symbols

    Example:
        >>> cols = ["ENSG00000167286.10", "ENSG00000010610.5", "other_col"]
        >>> build_versioned_rename_map(cols)
        {'ENSG00000167286.10': 'CD3D', 'ENSG00000010610.5': 'CD4'}
    """
    if genes is None:
        genes = TCELL_MARKERS
    rename_map = {}
    for col in column_names:
        # Skip if already a symbol
        if not col.startswith("ENSG"):
            continue
        gene = find_gene_by_ensembl(col, genes)
        if gene:
            rename_map[col] = gene.symbol
    return rename_map


# =============================================================================
# AnnData symbol resolution — THE shared path for expression lookups
# =============================================================================
#
# Every part of tcrsift that needs to find genes in an expression matrix
# (signature scoring, GEX aggregation, marker lookups) must resolve symbols
# through here, so the logic — including Ensembl↔HGNC mapping — lives in one
# place rather than scattered, bespoke copies.

_SYMBOL_COLS = ("gene_symbols", "feature_name", "symbol", "gene_symbol",
                "gene_name")

_ENSEMBL_DATA_CACHE: dict[int, object] = {}


def _ensembl_data(release: int = 110):
    """Cached :class:`pyensembl.EnsemblRelease`. pyensembl is a core dep."""
    if release not in _ENSEMBL_DATA_CACHE:
        import pyensembl  # pylint: disable=import-error

        _ENSEMBL_DATA_CACHE[release] = pyensembl.EnsemblRelease(release)
    return _ENSEMBL_DATA_CACHE[release]


def ensembl_ids_to_symbols(
    ensembl_ids, *, release: int = 110,
) -> dict[str, str]:
    """Map Ensembl gene IDs → HGNC symbols via pyensembl (version-robust).

    Returns ``{ensembl_id: symbol}`` for the IDs pyensembl recognizes. IDs it
    can't resolve are omitted. Requires the Ensembl release data to be
    installed (``pyensembl install --release 110``); if it isn't, logs once
    and returns ``{}`` so callers fall back to raw names.
    """
    import logging

    logger = logging.getLogger(__name__)
    try:
        data = _ensembl_data(release)
    except Exception as exc:  # pyensembl data not downloaded, etc.
        logger.warning(
            "pyensembl release %d unavailable (%s); Ensembl→symbol mapping "
            "skipped. Run `pyensembl install --release %d` to enable it.",
            release, exc, release,
        )
        return {}
    out: dict[str, str] = {}
    for eid in ensembl_ids:
        base = strip_ensembl_version(str(eid))
        try:
            gene = data.gene_by_id(base)
        except Exception:
            continue
        if gene is not None and gene.gene_name:
            out[str(eid)] = gene.gene_name
    return out


def adata_symbol_array(adata, *, release: int = 110, use_pyensembl: bool = True):
    """Per-``var`` HGNC symbol (UPPERCASE) for an AnnData — the shared resolver.

    Resolution order (single source of truth for the whole library):

    1. a symbol column in ``adata.var`` (``gene_symbols``/``feature_name``/…);
    2. ``var_names`` themselves when they look like symbols;
    3. Ensembl ``var_names`` → HGNC via :func:`ensembl_ids_to_symbols`
       (pyensembl), falling back to the raw Ensembl id for anything unmapped.

    Returns a numpy array of uppercase symbols aligned to ``adata.var``.
    """
    import pandas as pd

    var = adata.var
    for col in _SYMBOL_COLS:
        if col in getattr(var, "columns", []):
            return var[col].astype(str).str.upper().to_numpy()

    names = pd.Series(adata.var_names.astype(str))
    looks_ensembl = names.str.startswith("ENSG")
    if use_pyensembl and looks_ensembl.mean() > 0.5:
        mapping = ensembl_ids_to_symbols(names[looks_ensembl].tolist(),
                                         release=release)
        resolved = names.map(lambda n: mapping.get(n, n))
        return resolved.astype(str).str.upper().to_numpy()
    return names.str.upper().to_numpy()
