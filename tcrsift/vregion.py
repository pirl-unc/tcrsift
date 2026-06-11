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

"""Germline V-REGION framework reference + donor comparison (#286 follow-up).

The leader reference (``leaders.py``) covers only the signal peptide, which is
cleaved off the mature TCR. The functionally-retained part of the V domain —
FR1, CDR1, FR2, CDR2, FR3 — comes from the germline V gene and survives into the
receptor, so a donor germline SNP there is biologically real (unlike a leader
variant). This module compares each clone's observed framework against the IMGT
V-REGION germline and reports the divergence, feeding the ``framework`` rows of
:func:`tcrsift.assemble.collect_germline_variants`.

The reference is the IMGT/GENE-DB **V-REGION** FASTA (TRAV + TRBV), managed
through the tcrsift data cache (``tcrsift data download --db imgt_trav_vregion
imgt_trbv_vregion``) rather than vendored — it's large and rarely changes. When
the cache is empty every entry point degrades gracefully (returns ``{}`` / ``None``),
so the framework comparison is simply inert until the user downloads it.

CDR3 is deliberately excluded from the comparison: it's the V(D)J junction, not
germline-templated, so its hypervariability is not "divergence". Only the V-gene
framework up to the CDR3 boundary is compared.
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

# Memoized parsed reference: gene -> [(allele, functionality, vregion_aa), ...].
_REFERENCE: dict[str, list[tuple[str, str, str]]] | None = None


def _parse_vregion_fasta(text: str) -> dict[str, list[tuple[str, str, str]]]:
    """Parse an IMGT GENElect V-REGION HTML/FASTA dump → translated reference.

    Pulls the FASTA out of the largest ``<pre>`` block (GENElect wraps the
    records there), translates each nucleotide record in frame 0, and keeps only
    clean, stop-free V-REGION translations. IMGT header:
    ``ACCESSION|GENE*ALLELE|species|functionality|label|...``.
    """
    from .assemble import translate_dna

    pres = re.findall(r"<pre>(.*?)</pre>", text, re.DOTALL | re.IGNORECASE)
    fasta = html.unescape(re.sub(r"<[^>]+>", "", max(pres, key=lambda p: p.count(">")))) if pres else text
    ref: dict[str, list[tuple[str, str, str]]] = {}
    header, chunks = None, []

    def _flush():
        if header is None:
            return
        fields = header.split("|")
        if len(fields) < 4:
            return
        gene, _, allele = fields[1].partition("*")
        functionality = fields[3].strip()
        nt = "".join(chunks)
        if not nt:
            return
        aa, _ragged = translate_dna(nt.upper())
        # frame-0 V-REGION must translate clean (no internal stop / unknown)
        if not aa or "*" in aa or "X" in aa:
            return
        ref.setdefault(gene, []).append((allele, functionality, aa))

    for line in fasta.splitlines():
        if line.startswith(">"):
            _flush()
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    _flush()
    return ref


def _load_reference(data_dir=None, *, force: bool = False) -> dict[str, list[tuple[str, str, str]]]:
    """Load + memoize the V-REGION reference from the data cache.

    Returns ``{}`` when the reference HTML isn't cached (download via
    ``tcrsift data download --db imgt_trav_vregion imgt_trbv_vregion``), so every
    caller degrades to a no-op rather than raising.
    """
    global _REFERENCE
    if _REFERENCE is not None and not force:
        return _REFERENCE
    from .datacache import cached_path

    ref: dict[str, list[tuple[str, str, str]]] = {}
    for _locus, db in VREGION_DBS.items():
        path = cached_path(db, data_dir)
        if path is None:
            continue
        try:
            ref.update(_parse_vregion_fasta(path.read_text(errors="replace")))
        except Exception as e:  # a corrupt cache file must not break assembly
            logger.warning("V-REGION reference %s unreadable: %s", path, e)
    _REFERENCE = ref
    return ref


def reference_is_available(data_dir=None) -> bool:
    """True when at least one locus of the V-REGION reference is cached."""
    return bool(_load_reference(data_dir))


def germline_vregion(v_gene: str | None) -> list[tuple[str, str, str]]:
    """Reference V-REGIONs for a V gene → ``[(allele, functionality, aa), …]``.

    Resolves the IMGT slash name (``TRAV14/DV4``) and a bare-prefix call. Empty
    when the gene isn't in the (possibly absent) reference.
    """
    gene = normalize_vgene(v_gene)
    if gene is None:
        return []
    return _load_reference().get(gene, [])


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
    vdj_aa: str | None, cdr3: str | None, v_gene: str | None
):
    """Compare a clone's framework (FR1–FR3) to the closest germline V-REGION.

    Returns ``(allele, germline_framework_aa, identity, diff)`` for the best
    matching allele, or ``None`` when not comparable (no reference for the gene,
    CDR3 not locatable, or best identity below :data:`GERMLINE_MIN_IDENTITY` —
    i.e. the gene call or alignment is off). ``diff`` is ``identical`` or the
    per-position substitution string (``S3N;G29V``); CDR3 is excluded.
    """
    fr = _framework(vdj_aa, cdr3)
    if not fr:
        return None
    candidates = germline_vregion(v_gene)
    if not candidates:
        return None
    best = None  # (identity, allele, germline_overlap)
    for allele, _func, g in candidates:
        n = min(len(fr), len(g))
        if n == 0:
            continue
        matches = sum(a == b for a, b in zip(fr[:n], g[:n]))
        # Normalize by the donor framework length, NOT the overlap: the IMGT
        # V-REGION dump carries many partial alleles that are short prefixes of
        # the full one. Over the overlap a truncated allele scores a spurious
        # 1.0 and hijacks the tie-break — reporting the wrong allele AND silently
        # dropping any donor SNP past its truncated end. Penalizing the residues
        # it doesn't cover lets the full-length allele win and the real
        # divergence surface (#vregion-allele).
        identity = matches / len(fr)
        key = (identity, allele)
        if best is None or key > (best[0], best[1]):
            best = (identity, allele, g[: len(fr)])
    if best is None or best[0] < GERMLINE_MIN_IDENTITY:
        return None
    identity, allele, g_overlap = best
    diff = _diff_vs_germline(fr, g_overlap) or "identical"
    return allele, g_overlap, identity, diff
