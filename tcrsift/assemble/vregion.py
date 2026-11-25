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

"""Germline V-REGION & J-REGION reference + donor framework comparison (#286 follow-up).

The leader reference (``leaders.py``) covers only the signal peptide, which is
cleaved off the mature TCR. The functionally-retained part of the V domain —
FR1, CDR1, FR2, CDR2, FR3 — comes from the germline V gene and survives into the
receptor, so a donor germline SNP there is biologically real (unlike a leader
variant). This module compares each clone's observed framework against the IMGT
V-REGION germline and reports the divergence, feeding the ``framework`` rows of
:func:`tcrsift.assemble.collect_germline_variants`, and also supplies clean
IMGT V-REGION and J-REGION sequences (AA + NT) for opt-in germline rescue of
5'-truncated contigs and bare ``(v_gene, CDR3, j_gene)`` tables.

The references are the IMGT/GENE-DB **V-REGION** (TRAV + TRBV) and **J-REGION**
(TRAJ + TRBJ) FASTAs, managed through the tcrsift data cache
(``tcrsift data download --group germline`` / ``--group germline_mouse``) rather
than vendored under Apache-2.0. When the cache is empty every entry point
degrades gracefully (returns ``{}`` / ``[]`` / ``None``).
"""

from __future__ import annotations

import html
import logging
import re

from .leaders import (
    GERMLINE_MIN_IDENTITY,
    _diff_vs_germline,
    normalize_vgene,
)

logger = logging.getLogger(__name__)

# Data-cache DB keys (registered in datacache.DATABASES); one HTML per locus.
VREGION_DBS = {"TRAV": "imgt_trav_vregion", "TRBV": "imgt_trbv_vregion"}
VREGION_DBS_MOUSE = {
    "TRAV": "imgt_trav_vregion_mouse",
    "TRBV": "imgt_trbv_vregion_mouse",
}
VREGION_DBS_BY_SPECIES = {
    "human": VREGION_DBS,
    "mouse": VREGION_DBS_MOUSE,
}

JREGION_DBS = {"TRAJ": "imgt_traj_jregion", "TRBJ": "imgt_trbj_jregion"}
JREGION_DBS_MOUSE = {
    "TRAJ": "imgt_traj_jregion_mouse",
    "TRBJ": "imgt_trbj_jregion_mouse",
}
JREGION_DBS_BY_SPECIES = {
    "human": JREGION_DBS,
    "mouse": JREGION_DBS_MOUSE,
}

# Memoized parsed V-REGION reference: gene -> [(allele, functionality, vregion_aa), ...].
# ``_REFERENCE`` holds human (kept under its legacy name & 3-tuple shape for test
# monkeypatching); ``_MOUSE_REFERENCE`` holds mouse so the two species never collide.
_REFERENCE: dict[str, list[tuple[str, str, str]]] | None = None
_MOUSE_REFERENCE: dict[str, list[tuple[str, str, str]]] | None = None

# Parallel nucleotide lookup maps: (gene, allele) -> clean frame-0 uppercase NT.
_VREGION_NT: dict[tuple[str, str], str] | None = None
_MOUSE_VREGION_NT: dict[tuple[str, str], str] | None = None

# Memoized J-REGION reference: gene -> [(allele, functionality, jregion_aa, jregion_nt), ...].
_JREGION_REFERENCE: dict[str, list[tuple[str, str, str, str]]] | None = None
_MOUSE_JREGION_REFERENCE: dict[str, list[tuple[str, str, str, str]]] | None = None

_FR4_MOTIF_RE = re.compile(r"[FWC][GA][ACDEFGHIKLMNPQRSTVWY]G[ACDEFGHIKLMNPQRSTVWY]{5,8}$")


def normalize_jgene(j_gene: str | None) -> str | None:
    """Strip allele/suffix noise and ImmunoSEQ zero-padding → uppercase bare IMGT J-gene name."""
    if not isinstance(j_gene, str) or not j_gene.strip():
        return None
    bare = j_gene.strip().upper().split("*")[0].split()[0]
    return re.sub(r"(?<=[A-Z/-])0+(\d)", r"\1", bare)


def _gene_lookup_candidates(gene: str | None) -> list[str]:
    """Return ordered lookup keys including zero-unpadded ImmunoSEQ/VDJdb synonyms
    (e.g. ``TRBV05-01`` → ``TRBV5-1``, ``TRBJ02-07`` → ``TRBJ2-7``)."""
    if not gene:
        return []
    out = [gene]
    unpadded = re.sub(r"(?<=[A-Z/-])0+(\d)", r"\1", gene)
    if unpadded != gene:
        out.append(unpadded)
    return out


def _func_rank_key(func: str) -> int:
    raw = str(func).strip()
    if raw == "F":
        return 0
    fc = raw.strip("[]() ")
    if fc.startswith("F"):
        return 1
    if fc.startswith("ORF"):
        return 2
    return 3


def _alias_gene_dict(
    ref: dict[str, list],
    nt_map: dict[tuple[str, str], str] | None = None,
) -> None:
    """Sort alleles (F > ORF > P, then allele id) and alias slash-named genes."""
    for gene_key in list(ref):
        ref[gene_key].sort(key=lambda item: (_func_rank_key(item[1]), item[0]))

    prefix_targets: dict[str, set[str]] = {}
    alias_pairs: list[tuple[str, str]] = []
    for gene in list(ref):
        if "/" in gene:
            prefix_targets.setdefault(gene.split("/")[0], set()).add(gene)
            slashless = gene.replace("/", "")
            if slashless not in ref:
                ref[slashless] = ref[gene]
                alias_pairs.append((slashless, gene))
    for prefix, fulls in prefix_targets.items():
        if len(fulls) == 1 and prefix not in ref:
            full = next(iter(fulls))
            ref[prefix] = ref[full]
            alias_pairs.append((prefix, full))
    if nt_map is not None:
        for alias, full in alias_pairs:
            for item in ref[full]:
                allele = item[0]
                if (full, allele) in nt_map and (alias, allele) not in nt_map:
                    nt_map[(alias, allele)] = nt_map[(full, allele)]


def _extract_imgt_fasta(text: str) -> str:
    pres = re.findall(r"<pre>(.*?)</pre>", text, re.DOTALL | re.IGNORECASE)
    if pres:
        return html.unescape(re.sub(r"<[^>]+>", "", max(pres, key=lambda p: p.count(">"))))
    return text


def _parse_vregion_fasta_with_nt(
    text: str,
) -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    """Parse an IMGT GENElect V-REGION HTML/FASTA dump into AA 3-tuples + NT map."""
    from .codons import translate_dna

    fasta = _extract_imgt_fasta(text)
    ref: dict[str, list[tuple[str, str, str]]] = {}
    nt_map: dict[tuple[str, str], str] = {}
    header: str | None = None
    chunks: list[str] = []

    def _flush() -> None:
        if header is None:
            return
        fields = header.split("|")
        if len(fields) < 4:
            return
        gene_raw, _, allele = fields[1].partition("*")
        gene = gene_raw.strip().upper()
        functionality = fields[3].strip()
        # Strip ONLY whitespace and IMGT alignment gaps; preserve IUPAC codes
        # so ambiguous bases trigger 'X' rejection instead of frameshifting.
        nt_raw = re.sub(r"[\s.\-]", "", "".join(chunks)).upper().replace("U", "T")
        if not nt_raw:
            return
        n_codons = len(nt_raw) // 3
        if n_codons == 0:
            return
        aa, _ragged = translate_dna(nt_raw[: 3 * n_codons])
        if not aa or len(aa) != n_codons or "*" in aa or "X" in aa:
            return
        clean_nt = nt_raw[: 3 * len(aa)]
        ref.setdefault(gene, []).append((allele, functionality, aa))
        nt_map[(gene, allele)] = clean_nt

    for line in fasta.splitlines():
        if line.startswith(">"):
            _flush()
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    _flush()

    _alias_gene_dict(ref, nt_map)
    return ref, nt_map


def _parse_vregion_fasta(text: str) -> dict[str, list[tuple[str, str, str]]]:
    """Parse an IMGT GENElect V-REGION HTML/FASTA dump → translated reference."""
    ref, _nt = _parse_vregion_fasta_with_nt(text)
    return ref


def _parse_jregion_fasta_with_nt(
    text: str,
) -> dict[str, list[tuple[str, str, str, str]]]:
    """Parse an IMGT GENElect J-REGION HTML/FASTA dump → ``[(allele, func, aa, nt), ...]``.

    IMGT J-REGION records begin at the germline heptamer cleavage boundary (which
    need not fall on codon frame 0) and carry a 1-nt 3' splice-donor overhang into
    C-REGION. Testing reading frames ``f in (0, 1, 2)`` against the conserved
    ``[FW]GXG`` C-terminal FR4 anchor selects the exact in-frame codons up to the
    end of FR4 (excluding the 1-nt J→C splice overhang).
    """
    from .codons import translate_dna

    fasta = _extract_imgt_fasta(text)
    ref: dict[str, list[tuple[str, str, str, str]]] = {}
    header: str | None = None
    chunks: list[str] = []

    def _flush() -> None:
        if header is None:
            return
        fields = header.split("|")
        if len(fields) < 4:
            return
        gene_raw, _, allele = fields[1].partition("*")
        gene = gene_raw.strip().upper()
        functionality = fields[3].strip()
        nt_raw = re.sub(r"[\s.\-]", "", "".join(chunks)).upper().replace("U", "T")
        if len(nt_raw) < 27:
            return

        best_candidate: tuple[tuple[int, int], str, str] | None = None
        for frame in (0, 1, 2):
            rem = len(nt_raw) - frame
            n_codons = rem // 3
            if n_codons < 9:
                continue
            sub_nt = nt_raw[frame : frame + 3 * n_codons]
            aa, _ = translate_dna(sub_nt)
            if not aa or len(aa) != n_codons or "*" in aa or "X" in aa:
                continue
            if not _FR4_MOTIF_RE.search(aa):
                continue
            # Canonical spliced J-REGIONs leave a 1-base phase-1 overhang at 3'.
            overhang = rem - 3 * n_codons
            rank = (1 if overhang == 1 else 0, len(aa))
            if best_candidate is None or rank > best_candidate[0]:
                best_candidate = (rank, aa, sub_nt)

        if best_candidate is not None:
            _rank, aa_clean, nt_clean = best_candidate
            ref.setdefault(gene, []).append((allele, functionality, aa_clean, nt_clean))

    for line in fasta.splitlines():
        if line.startswith(">"):
            _flush()
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    _flush()

    _alias_gene_dict(ref)
    return ref


def clear_reference_cache() -> None:
    """Clear memoized human and mouse V-REGION and J-REGION reference tables."""
    global _REFERENCE, _MOUSE_REFERENCE, _VREGION_NT, _MOUSE_VREGION_NT
    global _JREGION_REFERENCE, _MOUSE_JREGION_REFERENCE
    _REFERENCE = None
    _MOUSE_REFERENCE = None
    _VREGION_NT = None
    _MOUSE_VREGION_NT = None
    _JREGION_REFERENCE = None
    _MOUSE_JREGION_REFERENCE = None


def _load_vregion_bundle(
    data_dir=None, *, force: bool = False, species: str = "human"
) -> tuple[dict[str, list[tuple[str, str, str]]], dict[tuple[str, str], str]]:
    global _REFERENCE, _MOUSE_REFERENCE, _VREGION_NT, _MOUSE_VREGION_NT
    from ..datacache import cached_path
    from ..species import normalize_species

    sp = normalize_species(species)
    cached_ref = _MOUSE_REFERENCE if sp == "mouse" else _REFERENCE
    cached_nt = _MOUSE_VREGION_NT if sp == "mouse" else _VREGION_NT
    if cached_ref is not None and not force and data_dir is None:
        return cached_ref, (cached_nt or {})

    ref: dict[str, list[tuple[str, str, str]]] = {}
    nt_map: dict[tuple[str, str], str] = {}
    for _locus, db in VREGION_DBS_BY_SPECIES[sp].items():
        path = cached_path(db, data_dir)
        if path is None:
            continue
        try:
            r_part, nt_part = _parse_vregion_fasta_with_nt(path.read_text(errors="replace"))
            ref.update(r_part)
            nt_map.update(nt_part)
        except Exception as e:  # a corrupt cache file must not break assembly
            logger.warning("V-REGION reference %s unreadable: %s", path, e)
    if data_dir is None and ref:
        if sp == "mouse":
            _MOUSE_REFERENCE = ref
            _MOUSE_VREGION_NT = nt_map
        else:
            _REFERENCE = ref
            _VREGION_NT = nt_map
    return ref, nt_map


def _load_reference(
    data_dir=None, *, force: bool = False, species: str = "human"
) -> dict[str, list[tuple[str, str, str]]]:
    """Load + memoize the V-REGION reference from the data cache."""
    ref, _nt = _load_vregion_bundle(data_dir=data_dir, force=force, species=species)
    return ref


def _load_jregion_reference(
    data_dir=None, *, force: bool = False, species: str = "human"
) -> dict[str, list[tuple[str, str, str, str]]]:
    """Load + memoize the J-REGION reference from the data cache."""
    global _JREGION_REFERENCE, _MOUSE_JREGION_REFERENCE
    from ..datacache import cached_path
    from ..species import normalize_species

    sp = normalize_species(species)
    cached_memo = _MOUSE_JREGION_REFERENCE if sp == "mouse" else _JREGION_REFERENCE
    if cached_memo is not None and not force and data_dir is None:
        return cached_memo

    ref: dict[str, list[tuple[str, str, str, str]]] = {}
    for _locus, db in JREGION_DBS_BY_SPECIES[sp].items():
        path = cached_path(db, data_dir)
        if path is None:
            continue
        try:
            ref.update(_parse_jregion_fasta_with_nt(path.read_text(errors="replace")))
        except Exception as e:
            logger.warning("J-REGION reference %s unreadable: %s", path, e)
    if data_dir is None and ref:
        if sp == "mouse":
            _MOUSE_JREGION_REFERENCE = ref
        else:
            _JREGION_REFERENCE = ref
    return ref


def reference_is_available(data_dir=None, *, species: str = "human") -> bool:
    """True when at least one locus of the V-REGION reference is cached."""
    return bool(_load_reference(data_dir, species=species))


def jregion_reference_is_available(data_dir=None, *, species: str = "human") -> bool:
    """True when at least one locus of the J-REGION reference is cached."""
    return bool(_load_jregion_reference(data_dir, species=species))


def germline_vregion(
    v_gene: str | None, *, species: str = "human", data_dir=None
) -> list[tuple[str, str, str]]:
    """Reference V-REGIONs for a V gene → ``[(allele, functionality, aa), …]``."""
    gene = normalize_vgene(v_gene)
    if gene is None:
        return []
    ref = _load_reference(data_dir, species=species)
    for cand in _gene_lookup_candidates(gene):
        if cand in ref:
            return ref[cand]
    return []


def germline_vregion_sequences(
    v_gene: str | None, *, species: str = "human", data_dir=None
) -> list[tuple[str, str, str, str]]:
    """Reference V-REGIONs with nucleotides → ``[(allele, functionality, aa, nt), …]``."""
    gene = normalize_vgene(v_gene)
    if gene is None:
        return []
    ref, nt_map = _load_vregion_bundle(data_dir=data_dir, species=species)
    for cand in _gene_lookup_candidates(gene):
        if cand in ref:
            return [
                (allele, func, aa, nt_map.get((cand, allele), ""))
                for allele, func, aa in ref[cand]
            ]
    return []


def germline_jregion_sequences(
    j_gene: str | None, *, species: str = "human", data_dir=None
) -> list[tuple[str, str, str, str]]:
    """Reference J-REGIONs with nucleotides → ``[(allele, functionality, aa, nt), …]``."""
    gene = normalize_jgene(j_gene)
    if gene is None:
        return []
    ref = _load_jregion_reference(data_dir=data_dir, species=species)
    for cand in _gene_lookup_candidates(gene):
        if cand in ref:
            return ref[cand]
    return []


def _framework(vdj_aa: str | None, cdr3: str | None) -> str | None:
    """The germline-templated V framework: the V(D)J amino acids up to the CDR3
    (exclusive). ``None`` when the CDR3 can't be located, so the junction is
    never mistaken for framework divergence."""
    if not isinstance(vdj_aa, str) or not vdj_aa:
        return None
    if isinstance(cdr3, str) and cdr3 and cdr3 in vdj_aa:
        return vdj_aa[: vdj_aa.index(cdr3)]
    return None


def germline_compare_vregion(
    vdj_aa: str | None,
    cdr3: str | None,
    v_gene: str | None,
    *,
    species: str = "human",
    data_dir=None,
):
    """Compare a clone's framework (FR1–FR3) to the closest germline V-REGION.

    Returns ``(allele, germline_framework_aa, identity, diff)`` for the best
    matching allele, or ``None`` when not comparable.
    """
    fr = _framework(vdj_aa, cdr3)
    if not fr:
        return None
    candidates = germline_vregion(v_gene, species=species, data_dir=data_dir)
    if not candidates:
        return None
    best = None  # (rank_score, identity, allele, germline_overlap)
    for allele, func, g in candidates:
        n = min(len(fr), len(g))
        if n == 0:
            continue
        matches = sum(a == b for a, b in zip(fr[:n], g[:n]))
        identity = matches / len(fr)
        raw_func = str(func).strip()
        fc = raw_func.strip("[]() ")
        func_rank = (
            3
            if raw_func == "F"
            else (2 if fc.startswith("F") else (1 if fc.startswith("ORF") else 0))
        )
        score = (identity, func_rank)
        if best is None or score > best[0]:
            best = (score, identity, allele, g[: len(fr)])
    if best is None or best[1] < GERMLINE_MIN_IDENTITY:
        return None
    _score, identity, allele, g_overlap = best
    diff = _diff_vs_germline(fr, g_overlap)
    if len(g_overlap) < len(fr):
        len_note = f"len:{len(fr)}vs{len(g_overlap)}"
        diff = f"{diff},{len_note}" if diff else len_note
    elif not diff:
        diff = "identical"
    return allele, g_overlap, identity, diff
