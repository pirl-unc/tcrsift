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

"""IMGT germline V+J reconstruction and 5'-truncated contig rescue.

Supports two opt-in modes (enabled via ``germline_vdj_fallback=True`` on
:func:`tcrsift.assemble.assemble_full_sequences` or directly via
:func:`reconstruct_germline_vdj`):

1. **5'-partial contig rescue** (``{chain}_vdj_source = "contig_5prime_germline_rescued"``):
   Targeted TCR / Cell Ranger 2–3 libraries (such as Fischer et al. 2021,
   ``GSE171037``) frequently start inside ``V-REGION`` downstream of the leader
   start codon while preserving 100% intact 3' ``FR2/FR3 + CDR3 + FR4 + C-REGION``
   donor sequence. When ``partial_VDJ_{chain}_aa`` / ``partial_VDJ_{chain}_nt``
   is available, ONLY the missing N-terminal germline FR1/CDR1 prefix is
   prepended from the best-matching IMGT ``V-REGION`` allele; every sequenced
   donor framework SNP, CDR3 codon, and FR4 codon is preserved verbatim.
2. **Bare ``(v_gene, CDR3, j_gene)`` reconstruction**
   (``{chain}_vdj_source = "germline_vj_reconstructed"``):
   When no contig VDJ sequence is available at all, splices germline ``FR1–FR3``
   (up to conserved ``Cys 104``) from IMGT ``V-REGION`` + ``CDR3`` + germline
   ``FR4`` (past conserved ``[FW]GXG``) from IMGT ``J-REGION``, preserving
   germline codons on templated V/J overlaps and back-translating only untemplated
   junctional residues when ``cdr3_nt`` is absent.

Every germline-imputed chain is explicitly stamped in ``{chain}_vdj_source``,
``{chain}_vdj_germline_allele``, and ``qc_warnings`` (``{chain}_vdj_germline_reconstructed``),
and marks ``construct_contig_verified = False`` so germline imputation is never
silently conflated with a 5'-to-3' contig-verified sequence.
"""

from __future__ import annotations

import logging
import re

import pandas as pd

from ._assemble_codons import back_translate, translate_dna
from .leaders import GERMLINE_MIN_IDENTITY, normalize_vgene
from .vregion import (
    germline_jregion_sequences,
    germline_vregion_sequences,
    normalize_jgene,
)

logger = logging.getLogger(__name__)

VDJ_SOURCE_CONTIG_FULL = "contig_full"
VDJ_SOURCE_5PRIME_RESCUED = "contig_5prime_germline_rescued"
VDJ_SOURCE_VJ_RECONSTRUCTED = "germline_vj_reconstructed"
GERMLINE_IMPUTED_VDJ_SOURCES = frozenset(
    {VDJ_SOURCE_5PRIME_RESCUED, VDJ_SOURCE_VJ_RECONSTRUCTED}
)

_FR4_TAIL_RE = re.compile(r"([FWC])([GA][ACDEFGHIKLMNPQRSTVWY]G[ACDEFGHIKLMNPQRSTVWY]{5,8})$")


def _v_func_rank(func: object) -> int:
    """Rank unbracketed functional ``F`` above partial ``[F]``/``(F)`` above ``ORF`` above ``P``."""
    raw = str(func).strip()
    if raw == "F":
        return 3
    fc = raw.strip("[]() ")
    if fc.startswith("F"):
        return 2
    if fc.startswith("ORF"):
        return 1
    return 0


def _explicit_allele(gene_raw: object) -> str | None:
    if not isinstance(gene_raw, str) or "*" not in gene_raw:
        return None
    allele = gene_raw.strip().split("*", 1)[1].split()[0].strip()
    return allele or None


def _find_c104_index(v_aa: str) -> int:
    """Locate conserved C-terminal Cys-104 in an IMGT V-REGION AA sequence.

    In IMGT V-REGION records, Cys-104 sits within the final 8 residues (followed
    by 0–6 germline-encoded CDR3 N-terminal residues such as ``CASS`` or ``CAV``).
    """
    tail_start = max(0, len(v_aa) - 8)
    rel = v_aa[tail_start:].rfind("C")
    if rel >= 0:
        return tail_start + rel
    # Rare non-Cys anchor fallback (e.g. F/W/Y at position 104 in atypical alleles):
    return max(0, len(v_aa) - 3)


def _ensure_nt_for_aa(aa: str, nt: str, *, codon_species: str) -> str:
    """Ensure ``nt`` is uppercase and translates cleanly to ``aa``; else back-translate.

    Allows 0–2 leading ragged bases on 5'-truncated contig segments (and ignores
    any 3' J→C splice-donor overhang base) so sequenced donor codons survive intact.
    """
    if nt and len(nt) >= 3 * len(aa):
        nt_up = nt.upper()
        for lead in (0, 1, 2):
            if len(nt_up) < lead + 3 * len(aa):
                continue
            cand = nt_up[lead : lead + 3 * len(aa)]
            tr, _ = translate_dna(cand)
            if tr == aa:
                return cand
    return back_translate(aa, species=codon_species).upper()


def _min_complete_v_len(candidates: list[tuple[str, str, str, str]]) -> int:
    """Minimum acceptable V-REGION length when full-length ``F`` alleles exist for a gene."""
    full_f_lens = [
        len(aa) for _al, fn, aa, _nt in candidates if str(fn).strip() == "F" and len(aa) >= 85
    ]
    if full_f_lens:
        return max(85, max(full_f_lens) - 3)
    return 0


def _pick_v_for_partial(
    v_gene_raw: object,
    donor_fr_suffix: str,
    *,
    species: str,
    codon_species: str,
    data_dir=None,
) -> tuple[str, str, str, int] | None:
    """Pick ``(allele_label, fr_prefix_aa, fr_prefix_nt, offset)`` for 5'-partial rescue."""
    candidates = germline_vregion_sequences(
        str(v_gene_raw) if v_gene_raw is not None else None,
        species=species,
        data_dir=data_dir,
    )
    if not candidates:
        return None

    requested_allele = _explicit_allele(v_gene_raw)
    gene_norm = normalize_vgene(str(v_gene_raw)) or str(v_gene_raw)
    min_v_len = _min_complete_v_len(candidates)
    best: tuple[tuple[float, int, int], str, str, str, int] | None = None

    for allele, func, v_aa, v_nt_raw in candidates:
        if min_v_len and str(func).strip() != "F" and len(v_aa) < min_v_len:
            continue
        c104 = _find_c104_index(v_aa)
        g_fr = v_aa[:c104]
        if not g_fr:
            continue
        v_nt = _ensure_nt_for_aa(v_aa, v_nt_raw, codon_species=codon_species)
        g_fr_nt = v_nt[: 3 * c104]

        func_rank = _v_func_rank(func)
        allele_pref = 1 if (requested_allele and allele == requested_allele) else 0

        if not donor_fr_suffix:
            score = (1.0, allele_pref, func_rank)
            if best is None or score > best[0]:
                best = (score, f"{gene_norm}*{allele}", g_fr, g_fr_nt, c104)
            continue

        # Slide ±3 residues around expected C-terminal alignment in case of
        # exonucleolytic nibbling before Cys-104 or 1-residue V-length variants.
        nominal = len(g_fr) - len(donor_fr_suffix)
        low = max(0, nominal - 3)
        high = min(len(g_fr), max(0, nominal + 3))
        for offset in range(low, high + 1):
            g_tail = g_fr[offset:]
            n_cmp = min(len(donor_fr_suffix), len(g_tail))
            if n_cmp == 0:
                continue
            matches = sum(a == b for a, b in zip(donor_fr_suffix[:n_cmp], g_tail[:n_cmp]))
            # Penalize length discrepancy slightly so exact C104 register wins ties.
            len_pen = abs(len(g_tail) - len(donor_fr_suffix)) * 0.02
            ident = (matches / max(len(donor_fr_suffix), len(g_tail))) - len_pen
            if n_cmp >= 3 and (matches / n_cmp) < GERMLINE_MIN_IDENTITY:
                continue
            score = (ident, allele_pref, func_rank)
            if best is None or score > best[0]:
                best = (
                    score,
                    f"{gene_norm}*{allele}",
                    g_fr[:offset],
                    g_fr_nt[: 3 * offset],
                    offset,
                )

    if best is None:
        return None
    _sc, label, prefix_aa, prefix_nt, offset = best
    return label, prefix_aa, prefix_nt, offset


def _reconstruct_from_bare_vj(
    v_gene_raw: object,
    cdr3_aa: str,
    j_gene_raw: object,
    *,
    cdr3_nt_raw: object = None,
    species: str = "human",
    codon_species: str = "human",
    data_dir=None,
) -> tuple[str, str, str] | None:
    """Reconstruct ``(vdj_aa, vdj_nt_lower, allele_summary)`` from ``(v_gene, CDR3, j_gene)``."""
    v_cands = germline_vregion_sequences(
        str(v_gene_raw) if v_gene_raw is not None else None,
        species=species,
        data_dir=data_dir,
    )
    j_cands = germline_jregion_sequences(
        str(j_gene_raw) if j_gene_raw is not None else None,
        species=species,
        data_dir=data_dir,
    )
    if not v_cands or not j_cands or not cdr3_aa:
        return None

    req_v = _explicit_allele(v_gene_raw)
    req_j = _explicit_allele(j_gene_raw)
    v_norm = normalize_vgene(str(v_gene_raw)) or str(v_gene_raw)
    j_norm = normalize_jgene(str(j_gene_raw)) or str(j_gene_raw)

    # Pick V allele: prefer requested full-length allele, complete F status, and CDR3 N-terminal match.
    min_v_len = _min_complete_v_len(v_cands)
    best_v = None
    for allele, func, v_aa, v_nt_raw in v_cands:
        if min_v_len and str(func).strip() != "F" and len(v_aa) < min_v_len:
            continue
        c104 = _find_c104_index(v_aa)
        if c104 == 0:
            continue
        v_cdr3_germline = v_aa[c104:]
        n_pref = 0
        for a, b in zip(cdr3_aa, v_cdr3_germline):
            if a == b:
                n_pref += 1
            else:
                break
        func_rank = _v_func_rank(func)
        score = (1 if (req_v and allele == req_v) else 0, func_rank, n_pref)
        if best_v is None or score > best_v[0]:
            best_v = (score, allele, v_aa, v_nt_raw, c104, n_pref)
    if best_v is None:
        return None
    _vs, v_allele, v_aa, v_nt_raw, c104, n_v_cdr3 = best_v
    v_nt = _ensure_nt_for_aa(v_aa, v_nt_raw, codon_species=codon_species)
    # Preserve conserved Cys-104 if the input CDR3 is IMGT anchorless (missing leading 'C').
    include_c104 = (not cdr3_aa.startswith("C")) and v_aa[c104] == "C"
    fr_end = c104 + 1 if include_c104 else c104
    fr1_fr3_aa = v_aa[:fr_end]
    fr1_fr3_nt = v_nt[: 3 * fr_end]

    # Pick J allele: locate conserved [FWC][GA]XG anchor and score C-terminal CDR3 overlap.
    best_j = None
    for allele, func, j_aa, j_nt_raw in j_cands:
        m = _FR4_TAIL_RE.search(j_aa)
        if m is None:
            continue
        anchor_pos = m.start(1)  # index of conserved F/W (or C in TRAJ35*01) in j_aa
        # Standard Cell Ranger / VDJdb CDR3 includes position 118 ([FW], or C on TRAJ35).
        cdr3_includes_anchor = cdr3_aa[-1] in ("F", "W") or (
            cdr3_aa[-1] == "C" and j_aa[anchor_pos] == "C"
        )
        j_cdr3_end = anchor_pos + (1 if cdr3_includes_anchor else 0)
        j_cdr3_germline = j_aa[:j_cdr3_end]
        post_cdr3_aa = j_aa[j_cdr3_end:]
        n_suf = 0
        for a, b in zip(reversed(cdr3_aa), reversed(j_cdr3_germline)):
            if a == b:
                n_suf += 1
            else:
                break
        func_rank = _v_func_rank(func)
        anchor_match = 1 if (cdr3_includes_anchor and cdr3_aa[-1] == j_aa[anchor_pos]) else 0
        score = (
            1 if (req_j and allele == req_j) else 0,
            anchor_match,
            func_rank,
            n_suf,
        )
        if best_j is None or score > best_j[0]:
            best_j = (score, allele, j_aa, j_nt_raw, j_cdr3_end, post_cdr3_aa, n_suf)
    if best_j is None:
        return None
    _js, j_allele, j_aa, j_nt_raw, j_cdr3_end, post_cdr3_aa, n_j_cdr3 = best_j
    j_nt = _ensure_nt_for_aa(j_aa, j_nt_raw, codon_species=codon_species)
    post_cdr3_nt = j_nt[3 * j_cdr3_end :]

    # Build CDR3 nucleotides: donor cdr3_nt if valid, else germline V-prefix + back-translated N/D + germline J-suffix.
    cdr3_nt_clean = ""
    if isinstance(cdr3_nt_raw, str) and len(cdr3_nt_raw.strip()) == 3 * len(cdr3_aa):
        cand_cdr3_nt = cdr3_nt_raw.strip().upper()
        tr_c, _ = translate_dna(cand_cdr3_nt)
        if tr_c == cdr3_aa:
            cdr3_nt_clean = cand_cdr3_nt
    if not cdr3_nt_clean:
        v_take = min(n_v_cdr3, len(cdr3_aa))
        j_take = min(n_j_cdr3, max(0, len(cdr3_aa) - v_take))
        v_cdr3_nt = v_nt[3 * c104 : 3 * (c104 + v_take)]
        mid_aa = cdr3_aa[v_take : len(cdr3_aa) - j_take]
        mid_nt = back_translate(mid_aa, species=codon_species).upper() if mid_aa else ""
        j_cdr3_nt = (
            j_nt[3 * (j_cdr3_end - j_take) : 3 * j_cdr3_end] if j_take > 0 else ""
        )
        cdr3_nt_clean = v_cdr3_nt + mid_nt + j_cdr3_nt

    vdj_aa = fr1_fr3_aa + cdr3_aa + post_cdr3_aa
    vdj_nt = (fr1_fr3_nt + cdr3_nt_clean + post_cdr3_nt).lower()
    tr_check, _ = translate_dna(vdj_nt.upper())
    if tr_check != vdj_aa or "*" in vdj_aa or "X" in vdj_aa:
        return None
    return vdj_aa, vdj_nt, f"{v_norm}*{v_allele};{j_norm}*{j_allele}"


def reconstruct_germline_vdj_chain(
    row: pd.Series | dict,
    chain: str,
    *,
    species: str = "human",
    codon_species: str = "human",
    data_dir=None,
) -> tuple[str | None, str | None, str | None, str | None]:
    """Resolve ``(vdj_aa, vdj_nt, vdj_source, germline_allele)`` for one chain."""
    existing_aa = row.get(f"VDJ_{chain}_aa")
    existing_nt = row.get(f"VDJ_{chain}_nt")
    if isinstance(existing_aa, str) and existing_aa.strip() and "*" not in existing_aa:
        aa_clean = existing_aa.strip()
        if isinstance(existing_nt, str) and existing_nt.strip():
            nt_clean = existing_nt.strip()
        else:
            nt_clean = back_translate(aa_clean, species=codon_species).lower()
        return aa_clean, nt_clean, VDJ_SOURCE_CONTIG_FULL, None

    cdr3_val = row.get(f"CDR3_{chain}")
    cdr3_aa = str(cdr3_val).strip() if isinstance(cdr3_val, str) and cdr3_val.strip() else ""
    if not cdr3_aa or "*" in cdr3_aa:
        return None, None, None, None

    v_gene = row.get(f"{chain}_v_gene")
    j_gene = row.get(f"{chain}_j_gene")

    # Path 1: 5'-partial contig rescue (preserves sequenced 3' FR2/FR3 + CDR3 + FR4).
    partial_aa_raw = row.get(f"partial_VDJ_{chain}_aa")
    partial_nt_raw = row.get(f"partial_VDJ_{chain}_nt")
    if (
        isinstance(partial_aa_raw, str)
        and partial_aa_raw.strip()
        and "*" not in partial_aa_raw
        and cdr3_aa in partial_aa_raw
        and pd.notna(v_gene)
    ):
        p_aa = partial_aa_raw.strip()
        cdr_idx = p_aa.rindex(cdr3_aa)
        # Ensure the partial contig preserves an intact 3' FR4 tail past CDR3
        # before trusting Path 1; otherwise fall through to Path 2 bare V+J.
        if _FR4_TAIL_RE.search(p_aa[cdr_idx:]) is None:
            picked = None
        else:
            donor_fr_suffix = p_aa[:cdr_idx]
            picked = _pick_v_for_partial(
                v_gene,
                donor_fr_suffix,
                species=species,
                codon_species=codon_species,
                data_dir=data_dir,
            )
        if picked is not None:
            allele_label, prefix_aa, prefix_nt, _offset = picked
            p_nt = _ensure_nt_for_aa(
                p_aa,
                str(partial_nt_raw).strip() if isinstance(partial_nt_raw, str) else "",
                codon_species=codon_species,
            )
            rescued_aa = prefix_aa + p_aa
            rescued_nt = (prefix_nt + p_nt).lower()
            tr_check, _ = translate_dna(rescued_nt.upper())
            if tr_check == rescued_aa and "*" not in rescued_aa:
                return rescued_aa, rescued_nt, VDJ_SOURCE_5PRIME_RESCUED, allele_label

    # Path 2: Bare (v_gene, CDR3, j_gene) full germline V+J reconstruction.
    if pd.notna(v_gene) and pd.notna(j_gene):
        cdr3_nt_col = row.get(f"CDR3_{chain}_nt") or row.get(f"{chain}_cdr3_nt")
        bare = _reconstruct_from_bare_vj(
            v_gene,
            cdr3_aa,
            j_gene,
            cdr3_nt_raw=cdr3_nt_col,
            species=species,
            codon_species=codon_species,
            data_dir=data_dir,
        )
        if bare is not None:
            rec_aa, rec_nt, allele_summary = bare
            return rec_aa, rec_nt, VDJ_SOURCE_VJ_RECONSTRUCTED, allele_summary

    return None, None, None, None


def reconstruct_germline_vdj(
    clonotypes: pd.DataFrame,
    *,
    species: str = "human",
    codon_species: str | None = None,
    data_dir=None,
    auto_download: bool = False,
) -> pd.DataFrame:
    """Populate missing ``VDJ_{alpha,beta}_{aa,nt}`` via IMGT germline rescue/reconstruction.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype table with ``CDR3_alpha``/``CDR3_beta``, ``{chain}_v_gene``,
        ``{chain}_j_gene``, and optional ``partial_VDJ_{chain}_{aa,nt}``.
    species : str, default ``"human"``
        VDJ organism species (``"human"``, ``"mouse"``, or ``"auto"``).
    codon_species : str or None
        Species codon table for any untemplated junctional back-translation.
    data_dir : str or Path, optional
        Override tcrsift reference cache directory.
    auto_download : bool, default False
        When True, automatically download missing IMGT V-REGION and J-REGION
        references for ``species`` into the cache before reconstructing.

    Returns
    -------
    pd.DataFrame
        Copy of ``clonotypes`` with ``VDJ_{chain}_aa``, ``VDJ_{chain}_nt``,
        ``{chain}_vdj_source``, ``{chain}_vdj_germline_allele``, and updated
        ``qc_warnings``.
    """
    from .datacache import GROUPS, ensure_cached
    from .species import normalize_species, resolve_species

    sp = resolve_species(species, df=clonotypes, warn_mismatch=False)
    cod_sp = normalize_species(codon_species) if codon_species else sp

    if auto_download:
        grp = "germline_mouse" if sp == "mouse" else "germline"
        for db in GROUPS[grp]:
            if "leader" in db:
                continue
            ensure_cached(db, data_dir=data_dir, auto_download=True)

    df = clonotypes.copy()
    if len(df) == 0:
        for chain in ("alpha", "beta"):
            df[f"{chain}_vdj_source"] = pd.Series(dtype=object)
            df[f"{chain}_vdj_germline_allele"] = pd.Series(dtype=object)
        return df

    existing_qcs: list[list[str]] = []
    if "qc_warnings" in df.columns:
        for val in df["qc_warnings"]:
            if isinstance(val, list):
                existing_qcs.append(list(val))
            elif isinstance(val, str) and val.strip():
                existing_qcs.append([val.strip()])
            else:
                existing_qcs.append([])
    else:
        existing_qcs = [[] for _ in range(len(df))]

    for chain in ("alpha", "beta"):
        aa_col = f"VDJ_{chain}_aa"
        nt_col = f"VDJ_{chain}_nt"
        aa_out: list[str | None] = []
        nt_out: list[str | None] = []
        src_out: list[str | None] = []
        allele_out: list[str | None] = []

        for row_pos, (_, row) in enumerate(df.iterrows()):
            vdj_aa, vdj_nt, src, g_allele = reconstruct_germline_vdj_chain(
                row,
                chain,
                species=sp,
                codon_species=cod_sp,
                data_dir=data_dir,
            )
            aa_out.append(vdj_aa if vdj_aa is not None else row.get(aa_col))
            nt_out.append(vdj_nt if vdj_nt is not None else row.get(nt_col))
            src_out.append(src)
            allele_out.append(g_allele)
            if src in GERMLINE_IMPUTED_VDJ_SOURCES:
                tag = f"{chain}_vdj_germline_reconstructed"
                if tag not in existing_qcs[row_pos]:
                    existing_qcs[row_pos].append(tag)

        df[aa_col] = aa_out
        df[nt_col] = nt_out
        df[f"{chain}_vdj_source"] = src_out
        df[f"{chain}_vdj_germline_allele"] = allele_out

    df["qc_warnings"] = existing_qcs
    return df
