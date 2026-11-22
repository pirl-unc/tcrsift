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
Centralized gene definitions and lookup utilities (human HGNC & mouse MGI).

Provides a single source of truth for gene identifiers across TCRsift, with:
- Version-robust Ensembl ID matching for both human (``ENSG*``) and mouse (``ENSMUSG*``).
- Offline-safe ``pyensembl`` local-release probing backed by a static Ensembl→symbol
  table for core immune/signature markers so offline machines never hang on FTP.
- 1:1 human↔mouse official gene-symbol harmonization (e.g. ``CD8B`` ↔ ``Cd8b1``).
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from .species import normalize_species, species_to_pyensembl


@dataclass(frozen=True)
class Gene:
    """A gene with its canonical output symbol, Ensembl ID, and optional aliases."""

    symbol: str
    ensembl_id: str
    description: str = ""
    aliases: tuple[str, ...] = ()

    def matches_ensembl(self, query: str) -> bool:
        """Check if query matches this gene's Ensembl ID (version-robust)."""
        query_base = query.split(".")[0]
        self_base = self.ensembl_id.split(".")[0]
        if query_base == self_base:
            return True
        for alias in self.aliases:
            if alias.startswith(("ENSG", "ENSMUSG")) and query_base == alias.split(".")[0]:
                return True
        return False


# =============================================================================
# T-cell marker genes (human GRCh38 & mouse GRCm39)
# =============================================================================
# Output column symbols stay standardized (`CD3D`, `CD3E`, `CD3G`, `CD4`, `CD8A`,
# `CD8B`) across species so downstream `CD3` and `CD8` composites work uniformly;
# `aliases` cover mouse `ENSMUSG*` IDs and mouse `Cd8b1`.

CD3D = Gene("CD3D", "ENSG00000167286", "CD3 delta subunit", aliases=("ENSMUSG00000032094", "Cd3d"))
CD3E = Gene("CD3E", "ENSG00000198851", "CD3 epsilon subunit", aliases=("ENSMUSG00000032093", "Cd3e"))
CD3G = Gene("CD3G", "ENSG00000160654", "CD3 gamma subunit", aliases=("ENSMUSG00000002033", "Cd3g"))
CD4 = Gene("CD4", "ENSG00000010610", "CD4 co-receptor", aliases=("ENSMUSG00000023274", "Cd4"))
CD8A = Gene("CD8A", "ENSG00000153563", "CD8 alpha subunit", aliases=("ENSMUSG00000053977", "Cd8a"))
CD8B = Gene(
    "CD8B",
    "ENSG00000172116",
    "CD8 beta subunit (mouse Cd8b1)",
    aliases=("ENSMUSG00000053044", "Cd8b1", "CD8B1"),
)

MOUSE_CD3D = Gene("CD3D", "ENSMUSG00000032094", "Mouse Cd3d subunit", aliases=("Cd3d",))
MOUSE_CD3E = Gene("CD3E", "ENSMUSG00000032093", "Mouse Cd3e subunit", aliases=("Cd3e",))
MOUSE_CD3G = Gene("CD3G", "ENSMUSG00000002033", "Mouse Cd3g subunit", aliases=("Cd3g",))
MOUSE_CD4 = Gene("CD4", "ENSMUSG00000023274", "Mouse Cd4 co-receptor", aliases=("Cd4",))
MOUSE_CD8A = Gene("CD8A", "ENSMUSG00000053977", "Mouse Cd8a subunit", aliases=("Cd8a",))
MOUSE_CD8B = Gene(
    "CD8B",
    "ENSMUSG00000053044",
    "Mouse Cd8b1 subunit",
    aliases=("Cd8b1", "CD8B1", "Cd8b"),
)

# Grouped markers
TCELL_MARKERS = [CD3D, CD3E, CD3G, CD4, CD8A, CD8B]
TCELL_MARKERS_MOUSE = [MOUSE_CD3D, MOUSE_CD3E, MOUSE_CD3G, MOUSE_CD4, MOUSE_CD8A, MOUSE_CD8B]
CD3_GENES = [CD3D, CD3E, CD3G]
CD8_GENES = [CD8A, CD8B]


def tcell_markers_for_species(species: str = "human") -> list[Gene]:
    """Return T-cell marker definitions for ``'human'`` or ``'mouse'``."""
    return TCELL_MARKERS_MOUSE if normalize_species(species) == "mouse" else TCELL_MARKERS


# =============================================================================
# Lookup utilities
# =============================================================================


def canonicalize_gene(name: str | None, *, keep_allele: bool = False) -> str | None:
    """Canonicalize a TCR V/J/C gene name to a single IMGT-style uppercase token.

    Handles:
    - Case normalization (mouse ``Trbv13-2*01`` / ``Trac`` → ``TRBV13-2`` / ``TRAC``)
    - Allele suffix (dropped unless ``keep_allele=True``)
    - Adaptive prefix (``TCRBV20-01`` → ``TRBV20-1``)
    - Leading zeros (``TRBV09-01`` → ``TRBV9-1``, ``TRAJ08`` → ``TRAJ8``)
    - DV/OR slash fusion (``TRAV14/DV4`` → ``TRAV14DV4``, ``Trav14d-3/dv8`` → ``TRAV14D-3DV8``)

    Returns ``None`` for empty / non-string input.
    """
    if not isinstance(name, str):
        return None
    g = name.strip().upper()
    if not g or g == "NAN":
        return None
    allele = ""
    if "*" in g:
        g, _, allele = g.partition("*")
    g = re.sub(r"^TCR([ABDG])", r"TR\1", g)
    g = g.replace("/", "")
    g = re.sub(r"(?<=[A-Z-])0+(\d)", r"\1", g)
    if keep_allele and allele:
        return f"{g}*{int(allele):02d}" if allele.isdigit() else f"{g}*{allele}"
    return g


def strip_ensembl_version(ensembl_id: str) -> str:
    """Remove version suffix from an Ensembl ID (``ENSG...`` or ``ENSMUSG...``)."""
    return ensembl_id.split(".")[0]


def find_gene_by_symbol(symbol: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by its symbol or symbol alias (case-insensitive)."""
    if genes is None:
        genes = TCELL_MARKERS
    symbol_upper = symbol.upper()
    for gene in genes:
        if gene.symbol.upper() == symbol_upper:
            return gene
        for alias in gene.aliases:
            if not alias.startswith(("ENSG", "ENSMUSG")) and alias.upper() == symbol_upper:
                return gene
    return None


def find_gene_by_ensembl(ensembl_id: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by its Ensembl ID (version-robust across human & mouse)."""
    if genes is None:
        genes = TCELL_MARKERS
    for gene in genes:
        if gene.matches_ensembl(ensembl_id):
            return gene
    return None


def find_gene(identifier: str, genes: Sequence[Gene] | None = None) -> Gene | None:
    """Find a gene by symbol, alias, or Ensembl ID."""
    if genes is None:
        genes = TCELL_MARKERS
    gene = find_gene_by_symbol(identifier, genes)
    if gene:
        return gene
    return find_gene_by_ensembl(identifier, genes)


def find_column_for_gene(
    gene: Gene,
    column_names: Sequence[str],
) -> str | None:
    """Find a column name matching a gene by symbol, alias, or Ensembl ID."""
    if gene.symbol in column_names:
        return gene.symbol

    candidate_symbols = [gene.symbol.lower()] + [
        a.lower() for a in gene.aliases if not a.startswith(("ENSG", "ENSMUSG"))
    ]
    col_lower_map: dict[str, str] = {}
    for col in column_names:
        cl = str(col).lower()
        if cl not in col_lower_map:
            col_lower_map[cl] = col
    for sym_l in candidate_symbols:
        if sym_l in col_lower_map:
            return col_lower_map[sym_l]

    for col in column_names:
        if gene.matches_ensembl(str(col)):
            return col

    return None


def build_versioned_rename_map(
    column_names: Sequence[str],
    genes: Sequence[Gene] | None = None,
) -> dict[str, str]:
    """Build a rename map from versioned ``ENSG*``/``ENSMUSG*`` IDs (or ``Cd8b1``)
    to canonical gene symbols.
    """
    if genes is None:
        genes = TCELL_MARKERS
    rename_map: dict[str, str] = {}
    for col in column_names:
        col_s = str(col)
        if col_s.startswith(("ENSG", "ENSMUSG")):
            gene = find_gene_by_ensembl(col_s, genes)
            if gene:
                rename_map[col] = gene.symbol
        else:
            gene = find_gene_by_symbol(col_s, genes)
            if gene and col_s != gene.symbol:
                rename_map[col] = gene.symbol
    return rename_map


# =============================================================================
# AnnData symbol resolution — THE shared path for expression lookups
# =============================================================================

_SYMBOL_COLS = ("gene_symbols", "feature_name", "symbol", "gene_symbol", "gene_name")

# Candidate releases probed locally (newest-first) when no exact release is cached.
_CANDIDATE_ENSEMBL_RELEASES = (115, 114, 113, 112, 111, 110, 109, 108, 104, 102)

# Static offline fallback for core human (GRCh38) and mouse (GRCm39) immune/signature
# genes so Ensembl-indexed AnnData works offline even without `pyensembl install`.
_STATIC_ENSEMBL_TO_SYMBOL: dict[str, str] = {
    # Human GRCh38
    "ENSG00000167286": "CD3D", "ENSG00000198851": "CD3E", "ENSG00000160654": "CD3G",
    "ENSG00000198821": "CD247", "ENSG00000010610": "CD4", "ENSG00000153563": "CD8A",
    "ENSG00000172116": "CD8B", "ENSG00000277734": "TRAC", "ENSG00000211751": "TRBC1",
    "ENSG00000211772": "TRBC2", "ENSG00000211829": "TRDC", "ENSG00000148773": "MKI67",
    "ENSG00000131747": "TOP2A", "ENSG00000132646": "PCNA", "ENSG00000176890": "TYMS",
    "ENSG00000134057": "CCNB1", "ENSG00000170312": "CDK1", "ENSG00000089685": "BIRC5",
    "ENSG00000180644": "PRF1", "ENSG00000100453": "GZMB", "ENSG00000145649": "GZMA",
    "ENSG00000113088": "GZMK", "ENSG00000115523": "GNLY", "ENSG00000100450": "GZMH",
    "ENSG00000105374": "NKG7", "ENSG00000134539": "KLRD1", "ENSG00000139187": "KLRG1",
    "ENSG00000111537": "IFNG", "ENSG00000168329": "CX3CR1", "ENSG00000073861": "TBX21",
    "ENSG00000057657": "PRDM1", "ENSG00000081059": "TCF7", "ENSG00000138795": "LEF1",
    "ENSG00000126353": "CCR7", "ENSG00000188404": "SELL", "ENSG00000168685": "IL7R",
    "ENSG00000139193": "CD27", "ENSG00000178562": "CD28", "ENSG00000049249": "TNFRSF9",
    "ENSG00000186827": "TNFRSF4", "ENSG00000134460": "IL2RA", "ENSG00000110848": "CD69",
    "ENSG00000102245": "CD40LG", "ENSG00000186891": "TNFRSF18", "ENSG00000123358": "NR4A1",
    "ENSG00000153234": "NR4A2", "ENSG00000122877": "EGR2", "ENSG00000156234": "CXCL13",
    "ENSG00000138185": "ENTPD1", "ENSG00000188389": "PDCD1", "ENSG00000089692": "LAG3",
    "ENSG00000135077": "HAVCR2", "ENSG00000181847": "TIGIT", "ENSG00000198846": "TOX",
    "ENSG00000163599": "CTLA4", "ENSG00000083457": "ITGAE", "ENSG00000277632": "CCL3",
    "ENSG00000275302": "CCL4", "ENSG00000049768": "FOXP3", "ENSG00000143365": "RORC",
    "ENSG00000107485": "GATA3", "ENSG00000169194": "IL13", "ENSG00000112115": "IL17A",
    "ENSG00000144290": "SLC4A10", "ENSG00000111796": "KLRB1", "ENSG00000109906": "ZBTB16",
    "ENSG00000189430": "NCR1", "ENSG00000164171": "ITGA2", "ENSG00000090382": "LYZ",
    "ENSG00000140678": "ITGAX", "ENSG00000156738": "MS4A1", "ENSG00000105369": "CD79A",
    "ENSG00000177455": "CD19", "ENSG00000129226": "CD68", "ENSG00000170458": "CD14",
    "ENSG00000203747": "FCGR3A", "ENSG00000261371": "PECAM1", "ENSG00000108821": "COL1A1",
    "ENSG00000107796": "ACTA2",
    # Mouse GRCm39 (official MGI casing; adata_symbol_array uppercases for uniform lookup)
    "ENSMUSG00000032094": "Cd3d", "ENSMUSG00000032093": "Cd3e", "ENSMUSG00000002033": "Cd3g",
    "ENSMUSG00000005763": "Cd247", "ENSMUSG00000023274": "Cd4", "ENSMUSG00000053977": "Cd8a",
    "ENSMUSG00000053044": "Cd8b1", "ENSMUSG00000076928": "Trac", "ENSMUSG00000076490": "Trbc1",
    "ENSMUSG00000076498": "Trbc2", "ENSMUSG00000104876": "Trdc", "ENSMUSG00000031004": "Mki67",
    "ENSMUSG00000020914": "Top2a", "ENSMUSG00000027342": "Pcna", "ENSMUSG00000025747": "Tyms",
    "ENSMUSG00000041431": "Ccnb1", "ENSMUSG00000019942": "Cdk1", "ENSMUSG00000017716": "Birc5",
    "ENSMUSG00000037202": "Prf1", "ENSMUSG00000015437": "Gzmb", "ENSMUSG00000023132": "Gzma",
    "ENSMUSG00000042385": "Gzmk", "ENSMUSG00000079186": "Gzmc", "ENSMUSG00000004612": "Nkg7",
    "ENSMUSG00000030165": "Klrd1", "ENSMUSG00000030114": "Klrg1", "ENSMUSG00000000817": "Fasl",
    "ENSMUSG00000055170": "Ifng", "ENSMUSG00000052336": "Cx3cr1", "ENSMUSG00000001444": "Tbx21",
    "ENSMUSG00000038151": "Prdm1", "ENSMUSG00000000782": "Tcf7", "ENSMUSG00000027985": "Lef1",
    "ENSMUSG00000037944": "Ccr7", "ENSMUSG00000026581": "Sell", "ENSMUSG00000003882": "Il7r",
    "ENSMUSG00000030336": "Cd27", "ENSMUSG00000026012": "Cd28", "ENSMUSG00000028965": "Tnfrsf9",
    "ENSMUSG00000029075": "Tnfrsf4", "ENSMUSG00000026770": "Il2ra", "ENSMUSG00000030156": "Cd69",
    "ENSMUSG00000031132": "Cd40lg", "ENSMUSG00000041954": "Tnfrsf18", "ENSMUSG00000023034": "Nr4a1",
    "ENSMUSG00000026826": "Nr4a2", "ENSMUSG00000037868": "Egr2", "ENSMUSG00000023078": "Cxcl13",
    "ENSMUSG00000048120": "Entpd1", "ENSMUSG00000026285": "Pdcd1", "ENSMUSG00000030124": "Lag3",
    "ENSMUSG00000020399": "Havcr2", "ENSMUSG00000071552": "Tigit", "ENSMUSG00000041272": "Tox",
    "ENSMUSG00000026011": "Ctla4", "ENSMUSG00000005947": "Itgae", "ENSMUSG00000000982": "Ccl3",
    "ENSMUSG00000018930": "Ccl4", "ENSMUSG00000039521": "Foxp3", "ENSMUSG00000028150": "Rorc",
    "ENSMUSG00000015619": "Gata3", "ENSMUSG00000020383": "Il13", "ENSMUSG00000025929": "Il17a",
    "ENSMUSG00000026904": "Slc4a10", "ENSMUSG00000079299": "Klrb1", "ENSMUSG00000030325": "Klrb1c",
    "ENSMUSG00000079298": "Klrb1b", "ENSMUSG00000066687": "Zbtb16", "ENSMUSG00000062524": "Ncr1",
    "ENSMUSG00000015533": "Itga2", "ENSMUSG00000022584": "Ly6c2", "ENSMUSG00000022582": "Ly6g",
    "ENSMUSG00000069516": "Lyz2", "ENSMUSG00000030789": "Itgax", "ENSMUSG00000024673": "Ms4a1",
    "ENSMUSG00000003379": "Cd79a", "ENSMUSG00000030724": "Cd19", "ENSMUSG00000018774": "Cd68",
    "ENSMUSG00000051439": "Cd14", "ENSMUSG00000059498": "Fcgr3", "ENSMUSG00000059089": "Fcgr4",
    "ENSMUSG00000020717": "Pecam1", "ENSMUSG00000001506": "Col1a1", "ENSMUSG00000035783": "Acta2",
    # Mitochondrial genes (GRCh38 & GRCm39) so bare-Ensembl mito QC works offline
    "ENSG00000198888": "MT-ND1", "ENSG00000198763": "MT-ND2", "ENSG00000198804": "MT-CO1",
    "ENSG00000198712": "MT-CO2", "ENSG00000228253": "MT-ATP8", "ENSG00000198899": "MT-ATP6",
    "ENSG00000198938": "MT-CO3", "ENSG00000198840": "MT-ND3", "ENSG00000212907": "MT-ND4L",
    "ENSG00000198886": "MT-ND4", "ENSG00000198786": "MT-ND5", "ENSG00000198695": "MT-ND6",
    "ENSG00000198727": "MT-CYB",
    "ENSMUSG00000064341": "mt-Nd1", "ENSMUSG00000064345": "mt-Nd2", "ENSMUSG00000064351": "mt-Co1",
    "ENSMUSG00000064354": "mt-Co2", "ENSMUSG00000064356": "mt-Atp8", "ENSMUSG00000064357": "mt-Atp6",
    "ENSMUSG00000064358": "mt-Co3", "ENSMUSG00000064360": "mt-Nd3", "ENSMUSG00000065947": "mt-Nd4l",
    "ENSMUSG00000064363": "mt-Nd4", "ENSMUSG00000064367": "mt-Nd5", "ENSMUSG00000064368": "mt-Nd6",
    "ENSMUSG00000064370": "mt-Cytb",
}

_ENSEMBL_DATA_CACHE: dict[tuple[str, int | None], object | None] = {}


def local_pyensembl_releases(species: str = "human") -> list[int]:
    """Return locally indexed pyensembl release numbers for ``species`` without network IO."""
    try:
        import pyensembl  # pylint: disable=import-error
    except ImportError:
        return []
    py_sp = species_to_pyensembl(species)
    found: list[int] = []
    for rel in _CANDIDATE_ENSEMBL_RELEASES:
        try:
            ens = pyensembl.EnsemblRelease(rel, species=py_sp)
            # Check `gtf.local_path_or_url` / SQLite file path via datacache without triggering download!
            db_dir = ens.gtf.cache_directory_path
            if os.path.isdir(db_dir) and any(f.endswith(".db") for f in os.listdir(db_dir)):
                found.append(rel)
        except Exception:
            continue
    return found


def _ensembl_data(release: int | None = 110, *, species: str = "human"):
    """Offline-safe cached :class:`pyensembl.EnsemblRelease`.

    Probes ``release`` first; if that exact release's SQLite DB is not installed
    locally, falls back to any locally installed release for ``species`` rather
    than triggering an implicit multi-100 MB FTP download from Ensembl.
    """
    sp = normalize_species(species)
    key = (sp, release)
    if key in _ENSEMBL_DATA_CACHE:
        return _ENSEMBL_DATA_CACHE[key]

    import pyensembl  # pylint: disable=import-error

    py_sp = species_to_pyensembl(sp)
    candidates: list[int] = []
    if release is not None:
        candidates.append(release)
    for r in _CANDIDATE_ENSEMBL_RELEASES:
        if r not in candidates:
            candidates.append(r)

    chosen = None
    for rel in candidates:
        try:
            ens = pyensembl.EnsemblRelease(rel, species=py_sp)
            db_dir = ens.gtf.cache_directory_path
            if os.path.isdir(db_dir) and any(f.endswith(".db") for f in os.listdir(db_dir)):
                chosen = ens
                break
        except Exception:
            continue

    _ENSEMBL_DATA_CACHE[key] = chosen
    return chosen


def ensembl_ids_to_symbols(
    ensembl_ids,
    *,
    release: int | None = 110,
    species: str | None = None,
) -> dict[str, str]:
    """Map Ensembl gene IDs (``ENSG*`` or ``ENSMUSG*``) → gene symbols (version-robust).

    Uses locally installed ``pyensembl`` when available (auto-detecting mouse from
    ``ENSMUSG*`` IDs when ``species`` is omitted), supplemented by the bundled
    static fallback map for core immune/signature/mitochondrial genes. Never
    initiates network downloads during analysis.
    """
    logger = logging.getLogger(__name__)
    ids_list = [str(e) for e in ensembl_ids]
    if not ids_list:
        return {}

    if species is None:
        n_mouse = sum(1 for e in ids_list if e.startswith("ENSMUSG"))
        n_human = sum(1 for e in ids_list if e.startswith("ENSG"))
        sp = "mouse" if n_mouse > n_human else "human"
    else:
        sp = normalize_species(species)

    out: dict[str, str] = {}
    data = None
    try:
        data = _ensembl_data(release, species=sp)
    except Exception as exc:
        logger.debug("pyensembl lookup skipped (%s)", exc)

    if data is not None:
        for eid in ids_list:
            base = strip_ensembl_version(eid)
            try:
                gene = data.gene_by_id(base)
            except Exception:
                gene = None
            if gene is not None and gene.gene_name:
                out[eid] = gene.gene_name
                continue
            if base in _STATIC_ENSEMBL_TO_SYMBOL:
                out[eid] = _STATIC_ENSEMBL_TO_SYMBOL[base]
    else:
        for eid in ids_list:
            base = strip_ensembl_version(eid)
            if base in _STATIC_ENSEMBL_TO_SYMBOL:
                out[eid] = _STATIC_ENSEMBL_TO_SYMBOL[base]
        if not out:
            rel_hint = release or (112 if sp == "mouse" else 110)
            logger.warning(
                "No local pyensembl %s database installed; unmapped Ensembl IDs kept raw. "
                "Install offline index with: pyensembl install --species %s --release %d",
                sp,
                species_to_pyensembl(sp),
                rel_hint,
            )
    return out


def ensembl_to_gene_name(
    ensembl_id: str,
    *,
    release: int | None = 110,
    species: str | None = None,
) -> str | None:
    """Resolve a single ``ENSG*`` or ``ENSMUSG*`` ID (with optional version suffix)
    to its official HGNC or MGI symbol offline.
    """
    eid = str(ensembl_id).strip()
    if not eid:
        return None
    mapping = ensembl_ids_to_symbols([eid], release=release, species=species)
    return mapping.get(eid)


# Signature genes whose identifier differs between GENCODE/Ensembl releases OR
# between 1:1 human HGNC and mouse MGI official symbols (uppercased).
# Each tuple is ONE orthologous gene group: Ensembl IDs + human/mouse symbols.
GENE_SYNONYMS: tuple[tuple[str, ...], ...] = (
    # Human GENCODE / HGNC version renames
    ("ENSG00000176485", "ENSMUSG00000060675", "PLA2G16", "PLAAT3"),
    ("ENSG00000110619", "ENSMUSG00000010755", "CARS", "CARS1"),
    ("ENSG00000184584", "ENSMUSG00000024349", "TMEM173", "STING1"),
    ("ENSG00000119979", "FAM45A", "DENND10"),
    ("ENSG00000104964", "AES", "TLE5"),
    ("ENSG00000135842", "FAM129A", "NIBAN1"),
    ("ENSG00000265206", "AC004687.1", "MIR142HG"),
    ("ENSG00000267364", "AC022706.1"),
    ("ENSG00000263264", "AC119396.1"),
    ("ENSG00000268027", "AC243960.1"),
    # 1:1 Human HGNC <-> Mouse MGI official symbol differences (uppercased)
    ("ENSG00000172116", "ENSMUSG00000053044", "CD8B", "CD8B1"),
    ("ENSG00000077238", "ENSMUSG00000030748", "IL4R", "IL4RA"),
    ("ENSG00000122223", "ENSMUSG00000004709", "CD244", "CD244A"),
    ("ENSG00000176083", "ENSMUSG00000049410", "ZNF683", "ZFP683"),
    ("ENSG00000100721", "ENSMUSG00000041359", "TCL1A", "TCL1"),
    ("ENSG00000104921", "ENSMUSG00000005540", "FCER2", "FCER2A"),
    ("ENSG00000118308", "ENSMUSG00000030263", "LRMP", "IRAG2"),
    ("ENSG00000133742", "ENSMUSG00000027556", "CA1", "CAR1"),
    ("ENSG00000105220", "ENSMUSG00000036427", "GPI", "GPI1"),
    ("ENSG00000101294", "ENSMUSG00000019188", "HM13", "H13"),
    ("ENSG00000178982", "C17ORF49", "BAP18", "2610507B11RIK"),
    ("ENSG00000235162", "C12ORF75", "OCC1", "1700025G04RIK"),
)
_SYNONYMS_OF = {name.upper(): tuple(x.upper() for x in group) for group in GENE_SYNONYMS for name in group}


def harmonize_gene_symbols(symbols, wanted):
    """Rename matrix symbols to the synonym a signature/caller asks for.

    ``symbols`` is the uppercase per-gene array from :func:`adata_symbol_array`;
    ``wanted`` the requested genes. Returns a new array; input is not mutated.
    """
    symbols = np.asarray(symbols, dtype=object).copy()
    present = set(symbols)
    renamed = {}
    for gene in {str(g).upper() for g in wanted} - present:
        for other in _SYNONYMS_OF.get(gene, ()):
            if other != gene and other in present:
                symbols[symbols == other] = gene
                renamed[other] = gene
                break
    if renamed:
        logging.getLogger(__name__).info(
            "Matched signature genes by synonym: %s",
            ", ".join(f"{new} (as {old})" for old, new in sorted(renamed.items())),
        )
    return symbols


def adata_symbol_array(
    adata,
    *,
    release: int | None = 110,
    use_pyensembl: bool = True,
    species: str | None = None,
):
    """Per-``var`` gene symbol (UPPERCASE) for an AnnData — the shared resolver.

    Resolution order (single source of truth across human and mouse):
    1. a symbol column in ``adata.var`` (``gene_symbols``/``feature_name``/…);
    2. ``var_names`` themselves when they look like symbols (MGI ``Gzmb`` → ``GZMB``);
    3. Ensembl ``var_names`` (``ENSG*`` or ``ENSMUSG*``) → symbols via
       :func:`ensembl_ids_to_symbols`, falling back to raw IDs when unmapped.
    """
    import pandas as pd

    var = adata.var
    for col in _SYMBOL_COLS:
        if col in getattr(var, "columns", []):
            return var[col].astype(str).str.upper().to_numpy()

    names = pd.Series(adata.var_names.astype(str))
    base_names = names.map(strip_ensembl_version)
    looks_human = base_names.str.startswith("ENSG")
    looks_mouse = base_names.str.startswith("ENSMUSG")
    looks_ensembl = looks_human | looks_mouse
    if looks_ensembl.mean() > 0.5:
        detected_sp = species or ("mouse" if looks_mouse.sum() > looks_human.sum() else "human")
        if use_pyensembl:
            mapping = ensembl_ids_to_symbols(
                names[looks_ensembl].tolist(),
                release=release,
                species=detected_sp,
            )
        else:
            mapping = {
                eid: _STATIC_ENSEMBL_TO_SYMBOL[strip_ensembl_version(eid)]
                for eid in names[looks_ensembl].tolist()
                if strip_ensembl_version(eid) in _STATIC_ENSEMBL_TO_SYMBOL
            }
        resolved = names.map(lambda n: mapping.get(n, n))
        return resolved.astype(str).str.upper().to_numpy()
    return names.str.upper().to_numpy()


def adata_mito_mask(adata) -> np.ndarray:
    """Boolean mask over ``adata.var`` for mitochondrial genes (``MT-*`` / ``mt-*``),
    working uniformly whether ``adata.var_names`` are HGNC, MGI, or ``ENSG``/``ENSMUSG``.
    """
    syms = adata_symbol_array(adata)
    return np.array([str(s).startswith("MT-") for s in syms], dtype=bool)


def resolve_adata_var_map(
    adata,
    wanted_genes: Sequence[str] | None = None,
) -> tuple[np.ndarray, dict[str, str]]:
    """Map uppercase gene symbols (after synonym/ortholog harmonization) to the
    actual ``adata.var_names`` entries in ``adata``.

    When duplicate ``var`` entries share a symbol, picks the highest-mean entry.
    Returns ``(sym_upper_array, {SYMBOL_UPPER: var_name})``.
    """
    sym_upper = adata_symbol_array(adata)
    if wanted_genes is not None:
        sym_upper = harmonize_gene_symbols(sym_upper, wanted_genes)

    var_names = np.asarray(adata.var_names.astype(str), dtype=object)
    by_sym: dict[str, list[int]] = {}
    for idx, s in enumerate(sym_upper):
        key = str(s).upper()
        if key:
            by_sym.setdefault(key, []).append(idx)

    X = getattr(adata, "X", None)
    sym_to_var: dict[str, str] = {}
    for sym, indices in by_sym.items():
        if len(indices) == 1 or X is None:
            chosen = indices[0]
        else:
            means = np.asarray(X[:, indices].mean(axis=0)).ravel()
            chosen = indices[int(np.argmax(means))]
        sym_to_var[sym] = str(var_names[chosen])
    return sym_upper, sym_to_var
