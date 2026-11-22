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

"""Germline V-gene signal-peptide (leader) reference + anchoring (#267).

The bundled reference (``tcrsift/refseqs/vgene_leaders.csv.gz``, built from IMGT
L-PART1+L-PART2 by ``scripts/build_vgene_leaders.py``) gives the canonical
leader per TRAV/TRBV allele. We use it as an **alignment anchor**, not a
template: the donor contig's native leader is kept (preserving allelic
variants); germline only picks WHERE the leader starts and HOW long it is.

Because the longest-ORF over-capture prepends 5'-UTR-encoded residues, the
donor leader is the **suffix** of the ORF-translated leader region whose length
equals the germline leader and whose first residue is the start Met. The anchor
is confidence-gated: gene absent / no allele aligns / identity below threshold
→ ``germline_anchor_leader`` returns ``None`` and the caller falls back to the
Kozak + h-region heuristic.
"""

from __future__ import annotations

import functools
import gzip
import html
import logging
import re
from importlib.resources import files as _pkg_files

logger = logging.getLogger(__name__)

# Minimum donor-vs-germline identity (over the germline leader length) for a
# confident anchor. Allelic variation is a few substitutions (TRBV13 is 27/29 =
# 0.93); a weird allele / wrong gene / truncation falls below this → heuristic.
GERMLINE_MIN_IDENTITY = 0.80

LEADER_DBS_MOUSE = {
    "TRAV": "imgt_trav_leader_mouse",
    "TRBV": "imgt_trbv_leader_mouse",
}


def _is_productive(functionality: str) -> bool:
    """True for a functional/ORF allele, tolerating IMGT bracket forms.

    IMGT marks inferred functionality with brackets/parens (``[F]``, ``(ORF)``);
    strip them before the F/ORF test so those alleles still sort as productive.
    """
    return functionality.strip("[]()") in ("F", "ORF")


def _alias_slash_genes(
    by_gene: dict[str, list[tuple[str, str, str]]],
) -> dict[str, list[tuple[str, str, str]]]:
    """Sort productive alleles first and alias unambiguous slash-named genes."""
    for recs in by_gene.values():
        recs.sort(key=lambda r: (0 if _is_productive(r[1]) else 1, r[0]))
    prefix_targets: dict[str, set[str]] = {}
    slashless_targets: dict[str, str] = {}
    for gene in list(by_gene):
        if "/" in gene:
            prefix_targets.setdefault(gene.split("/")[0], set()).add(gene)
            slashless = gene.replace("/", "")
            if slashless not in by_gene:
                slashless_targets[slashless] = gene
    for slashless, full in slashless_targets.items():
        by_gene[slashless] = by_gene[full]
    for prefix, fulls in prefix_targets.items():
        if len(fulls) == 1 and prefix not in by_gene:
            by_gene[prefix] = by_gene[next(iter(fulls))]
    return by_gene


def _parse_leader_csv_bytes(
    raw: bytes,
) -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    """Parse bundled ``vgene_leaders*.csv.gz`` bytes into AA + NT indices."""
    by_gene: dict[str, list[tuple[str, str, str]]] = {}
    nt_index: dict[tuple[str, str], str] = {}
    text = gzip.decompress(raw).decode()
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        gene, allele, func, leader_aa, leader_nt = line.split(",")
        gene_u = gene.strip().upper()
        by_gene.setdefault(gene_u, []).append((allele, func, leader_aa))
        nt_index[(gene_u, allele)] = leader_nt
    return _alias_slash_genes(by_gene), nt_index


def _parse_imgt_leader_html(
    text: str,
) -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    """Parse an IMGT GENElect ``L-PART1+L-PART2`` HTML/FASTA dump."""
    from .assemble import translate_dna

    pres = re.findall(r"<pre>(.*?)</pre>", text, re.DOTALL | re.IGNORECASE)
    fasta = (
        html.unescape(re.sub(r"<[^>]+>", "", max(pres, key=lambda p: p.count(">"))))
        if pres
        else text
    )
    by_gene: dict[str, list[tuple[str, str, str]]] = {}
    nt_index: dict[tuple[str, str], str] = {}
    header: str | None = None
    chunks: list[str] = []

    def _flush():
        if header is None:
            return
        fields = header.split("|")
        if len(fields) < 4:
            return
        gene_raw, _, allele = fields[1].partition("*")
        gene = gene_raw.strip().upper()
        allele = allele.strip() or "01"
        func = fields[3].strip()
        nt = re.sub(r"[^ACGTUacgtu]", "", "".join(chunks)).upper().replace("U", "T")
        if len(nt) < 36:
            return
        nt_trimmed = nt[: (len(nt) // 3) * 3]
        aa, _ = translate_dna(nt_trimmed)
        if not aa or not aa.startswith("M") or "*" in aa or "X" in aa:
            return
        by_gene.setdefault(gene, []).append((allele, func, aa))
        nt_index[(gene, allele)] = nt_trimmed

    for line in fasta.splitlines():
        if line.startswith(">"):
            _flush()
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    _flush()
    return _alias_slash_genes(by_gene), nt_index


def _load_mouse_leaders(
    data_dir=None,
) -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    """Load mouse leader references from optional bundled CSV or datacache."""
    try:
        path = _pkg_files("tcrsift.refseqs").joinpath("vgene_leaders_mouse.csv.gz")
        raw = path.read_bytes()
        return _parse_leader_csv_bytes(raw)
    except (FileNotFoundError, ModuleNotFoundError, OSError):
        pass

    from .datacache import cached_path

    by_gene: dict[str, list[tuple[str, str, str]]] = {}
    nt_index: dict[tuple[str, str], str] = {}
    for _locus, db in LEADER_DBS_MOUSE.items():
        cached = cached_path(db, data_dir)
        if cached is None:
            continue
        try:
            g_part, nt_part = _parse_imgt_leader_html(
                cached.read_text(errors="replace")
            )
            for gene, recs in g_part.items():
                by_gene.setdefault(gene, []).extend(recs)
            nt_index.update(nt_part)
        except Exception as exc:
            logger.warning("Mouse leader reference %s unreadable: %s", cached, exc)
    if by_gene:
        _alias_slash_genes(by_gene)
    return by_gene, nt_index


@functools.lru_cache(maxsize=4)
def _load_human_reference() -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    try:
        path = _pkg_files("tcrsift.refseqs").joinpath("vgene_leaders.csv.gz")
        raw = path.read_bytes()
    except (FileNotFoundError, ModuleNotFoundError, OSError):
        return {}, {}
    return _parse_leader_csv_bytes(raw)


_MOUSE_LEADER_CACHE: tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]] | None = None


def clear_reference_cache() -> None:
    """Clear memoized human and mouse leader reference tables."""
    global _MOUSE_LEADER_CACHE
    _MOUSE_LEADER_CACHE = None
    _load_human_reference.cache_clear()


def _load_reference(species: str = "human") -> dict[str, list[tuple[str, str, str]]]:
    """Load the leader reference for ``species`` → ``{gene: [(allele, functionality, aa), …]}``.

    Functional/ORF alleles are listed before pseudogenes so the anchor prefers a
    productive germline. Dual α/δ and orphon genes carry an IMGT slash name
    (``TRAV14/DV4``, ``TRAV14D-3/DV8``); each is ALSO aliased under its bare
    prefix (``TRAV14``, ``TRAV14D-3``) when that prefix is unambiguous and not
    already a distinct gene. Returns ``{}`` when no reference is available for
    ``species`` so mouse TCRs never cross-anchor against human leaders.
    """
    global _MOUSE_LEADER_CACHE
    from .species import normalize_species

    sp = normalize_species(species)
    if sp == "mouse":
        if _MOUSE_LEADER_CACHE is not None:
            return _MOUSE_LEADER_CACHE[0]
        by_gene, nt_idx = _load_mouse_leaders()
        if by_gene:
            _MOUSE_LEADER_CACHE = (by_gene, nt_idx)
        return by_gene
    by_gene, _ = _load_human_reference()
    return by_gene


def _leader_nt_index(species: str = "human") -> dict[tuple[str, str], str]:
    """``{(gene, allele): leader_nt}`` from the reference for ``species``."""
    global _MOUSE_LEADER_CACHE
    from .species import normalize_species

    sp = normalize_species(species)
    if sp == "mouse":
        if _MOUSE_LEADER_CACHE is not None:
            return _MOUSE_LEADER_CACHE[1]
        by_gene, nt_idx = _load_mouse_leaders()
        if by_gene:
            _MOUSE_LEADER_CACHE = (by_gene, nt_idx)
        return nt_idx
    _, nt_idx = _load_human_reference()
    return nt_idx


_load_reference.cache_clear = clear_reference_cache  # type: ignore[attr-defined]
_leader_nt_index.cache_clear = clear_reference_cache  # type: ignore[attr-defined]


def germline_leader_nt(
    v_gene: str | None,
    allele: str,
    *,
    species: str = "human",
    data_dir=None,
) -> str | None:
    """Germline leader nucleotides for a (gene, allele), or None if absent."""
    from .species import normalize_species

    gene = normalize_vgene(v_gene)
    if gene is None:
        return None
    sp = normalize_species(species)
    if data_dir is not None and sp == "mouse":
        _, nt_idx = _load_mouse_leaders(data_dir)
        return nt_idx.get((gene, allele))
    return _leader_nt_index(sp).get((gene, allele))


def normalize_vgene(v_gene: str | None) -> str | None:
    """Strip allele/suffix noise from a V-gene call → uppercase bare IMGT gene name.

    ``TRBV13*01`` / ``Trbv13-2*01`` / ``TRBV13 (gene)`` → ``TRBV13`` / ``TRBV13-2``;
    ``TRAV14/DV4*01`` → ``TRAV14/DV4``. Returns None for empty/non-string input.
    """
    if not isinstance(v_gene, str) or not v_gene.strip():
        return None
    return v_gene.strip().upper().split("*")[0].split()[0]


def germline_vgene_leaders(
    v_gene: str | None,
    *,
    species: str = "human",
    data_dir=None,
) -> list[tuple[str, str, str]]:
    """Reference leaders for a V gene → ``[(allele, functionality, aa), …]``.

    Resolves both the IMGT slash name (``TRAV14/DV4``) and a bare-prefix call
    (``TRAV14``). Empty when the gene isn't in the species reference (caller
    falls back to the heuristic).
    """
    from .species import normalize_species

    gene = normalize_vgene(v_gene)
    if gene is None:
        return []
    sp = normalize_species(species)
    if data_dir is not None and sp == "mouse":
        by_gene, _ = _load_mouse_leaders(data_dir)
        return by_gene.get(gene, [])
    return _load_reference(sp).get(gene, [])


def characterize_divergence(donor: str, germline: str) -> str:
    """Describe how a (length-mismatched) leader diverges from germline (#270).

    Detects the common single internal deletion (donor = germline's prefix +
    suffix, missing a contiguous middle — the CellRanger repeat-collapse / rare-
    germline-indel signature) as ``internal_deletion:Δ{n}@{pos}({chunk})``; a
    pure C-terminal shortfall as ``5p_truncation`` only when the donor is a clean
    suffix of germline; an insertion when longer; else a generic length note.
    """
    if len(donor) == len(germline):
        return _diff_vs_germline(donor, germline) or "identical"
    if len(donor) > len(germline):
        return f"insertion ({len(donor)} vs {len(germline)} aa)"
    # donor shorter: maximal shared prefix + suffix
    p = 0
    while p < len(donor) and donor[p] == germline[p]:
        p += 1
    s = 0
    while s < len(donor) - p and donor[-1 - s] == germline[-1 - s]:
        s += 1
    if p + s == len(donor) and p > 0 and s > 0:
        missing = germline[p:len(germline) - s]
        return f"internal_deletion:Δ{len(missing)}@{p + 1}({missing})"
    if donor == germline[-len(donor):]:
        return f"5p_truncation:Δ{len(germline) - len(donor)}(N-term)"
    return f"length_mismatch ({len(donor)} vs {len(germline)} aa)"


def _diff_vs_germline(donor: str, germline: str) -> str:
    """Per-position substitutions donor-vs-germline as ``G17R;V26G`` (1-based,
    germline residue + position + donor residue). Empty string if identical."""
    diffs = [
        f"{g}{i}{d}"
        for i, (d, g) in enumerate(zip(donor, germline), start=1)
        if d != g
    ]
    return ";".join(diffs)


def germline_anchor_leader(
    orig_leader: str,
    v_gene: str | None,
    *,
    min_identity: float = GERMLINE_MIN_IDENTITY,
    species: str = "human",
    data_dir=None,
):
    """Anchor the leader start/length to the germline reference (#267).

    ``orig_leader`` is the ORF-translated leader region (everything before the
    VDJ, possibly an over-capture). For each germline allele leader we test the
    suffix of ``orig_leader`` of the germline length: it must begin with Met and
    match the germline above ``min_identity`` (tolerating allelic variants). The
    best-matching allele wins.

    Returns ``(leader_aa, germline_allele, germline_aa, identity, diff)`` — where
    ``leader_aa`` is the DONOR suffix (native residues kept) — or ``None`` when
    no allele anchors confidently (gene absent / weird allele / too divergent /
    contig too short), so the caller uses the Kozak + h-region heuristic.

    Assumes substitution-level allelic variation (the donor leader is the same
    LENGTH as some germline allele). A genuine indel shifts the suffix frame,
    drops identity below ``min_identity``, and falls back to the heuristic — so
    an indel leader is never silently mis-anchored, just not germline-resolved.
    """
    candidates = germline_vgene_leaders(v_gene, species=species, data_dir=data_dir)
    if not candidates or not orig_leader:
        return None

    best = None  # (identity, allele, germline_aa, donor_suffix)
    for allele, _func, g_aa in candidates:
        n = len(g_aa)
        if n > len(orig_leader):
            continue  # contig leader shorter than germline — can't suffix-anchor
        cand = orig_leader[-n:]
        if not cand.startswith("M"):
            continue  # start Met must align with the germline's
        identity = sum(a == b for a, b in zip(cand, g_aa)) / n
        if best is None or identity > best[0]:
            best = (identity, allele, g_aa, cand)

    if best is None or best[0] < min_identity:
        return None
    identity, allele, g_aa, cand = best
    return cand, allele, g_aa, identity, _diff_vs_germline(cand, g_aa)


def germline_compare_leader(
    leader_aa: str | None,
    v_gene: str | None,
    *,
    species: str = "human",
    data_dir=None,
):
    """Compare a FINAL (shipped) leader to the closest germline allele (#267).

    Independent of anchoring: works for *any* leader — germline-anchored,
    Kozak/h-region heuristic, or curated default — so we always report which
    germline allele the V gene maps to and how divergent the leader we're
    shipping is from it. (A curated CD28/CD8A leader vs the native germline reads
    as a large divergence, which is exactly the useful signal.)

    Leaders are right-aligned (they share the C-terminal cleavage region; any
    over-capture / N-terminal variation differs). Identity is over the LONGER of
    the two, so a length mismatch lowers it. Returns
    ``(allele, germline_aa, identity, diff)`` — ``diff`` is the per-position
    substitution string when lengths match, else a ``length_mismatch`` marker —
    or ``None`` when the gene isn't in the reference.
    """
    candidates = germline_vgene_leaders(v_gene, species=species, data_dir=data_dir)
    if not candidates or not isinstance(leader_aa, str) or not leader_aa:
        return None
    best = None  # (identity, allele, germline_aa)
    for allele, _func, g_aa in candidates:
        n = min(len(leader_aa), len(g_aa))
        matches = sum(a == b for a, b in zip(leader_aa[-n:], g_aa[-n:]))
        identity = matches / max(len(leader_aa), len(g_aa))
        if best is None or identity > best[0]:
            best = (identity, allele, g_aa)
    identity, allele, g_aa = best
    diff = characterize_divergence(leader_aa, g_aa)
    return allele, g_aa, round(identity, 3), diff
