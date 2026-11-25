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

"""Contig leader extraction, Kozak QC, and leader fallback/forcing policy."""

from __future__ import annotations

import logging
from collections import Counter

import pandas as pd

from ..validation import TCRsiftValidationError
from .codons import (
    DEFAULT_FORBIDDEN_MOTIFS,
    CodonOptimization,
    _add_single_chain,
    find_longest_orf,
    optimize_codons,
    recode_codons,
)
from .constants import (
    DEFAULT_LEADERS,
    KOZAK_ADEQUATE,
    SP_LENGTH_MAX,
    SP_LENGTH_MIN,
    _has_h_region,
    _kozak_score,
    _lookup_default_leader,
    _record_leader,
)
from .leaders import (
    characterize_divergence,
    germline_anchor_leader,
    germline_compare_leader,
    germline_leader_nt,
    germline_vgene_leaders,
    normalize_vgene,
)
from .qc import add_synthesis_qc
from .sigpep import signal_peptide_features, sp_features_summary

logger = logging.getLogger(__name__)


def _kozak_correct_leader(
    contig_seq: str, offset: int, orig_leader: str
) -> tuple[str, str, int, str]:
    """Choose the leader's start codon, auto-shortening ONLY an over-capture (#263).

    ``orig_leader`` is the ORF prefix before the VDJ (translated from the ATG at
    ``offset``). Policy:

    - **In-range leaders (len ≤ SP_LENGTH_MAX) are never altered.** The ORF
      start is kept. This avoids silently shortening a correct leader whose true
      start happens to have weak Kozak context (a downstream adequate Met must
      NOT win here). Weakness of the kept start is still surfaced as QC by the
      caller (``weak_kozak_start``) — leaky-scan informs, it doesn't act.
    - **Over-captures (len > SP_LENGTH_MAX)** trigger a leaky scan: the in-frame
      Mets are walked 5'→3' and initiation moves to the FIRST one with adequate
      Kozak context (``>= KOZAK_ADEQUATE``) whose resulting leader lands in the
      ``SP_LENGTH_MIN..MAX`` window — what the ribosome does, leaky-scanning past
      weak upstream AUGs. If no such start exists the over-capture is kept (the
      caller flags it ``too_long``); we never fabricate a worse leader.

    - **No adequate-Kozak start, but a structurally-valid nested SP exists**
      (over-capture only): fall back to a *best-available STRUCTURAL trim* —
      among in-range nested Mets whose leader carries a hydrophobic h-region
      (the SP's defining feature), take the best-available Kozak (even if weak),
      tie-broken toward the most DOWNSTREAM (shortest, "SP at the end"). This
      catches a real SP whose start is genuinely weak-Kozak — e.g. the TRAV1
      family leader ``MWG.FLLYVSMKMGGT.`` nested at the tail of a 28-34 aa
      over-capture. Flagged ``contig_hregion_trimmed`` so the weaker basis is
      visible. If no such nested SP exists (e.g. a genuinely long ~29 aa leader
      with no shorter Met, like TRBV13), the ORF start is kept and the caller
      bands it by length — ``long_leader`` (26-30 aa, accepted) or ``too_long``
      (> SP_LENGTH_HARD); we never fabricate one.

    Returns ``(leader_aa, leader_nt, kozak_score_of_chosen_start, source)`` where
    ``source`` is ``"contig"`` (ORF start kept), ``"contig_kozak_reselected"``
    (over-capture → in-range adequate-Kozak start), or
    ``"contig_hregion_trimmed"`` (over-capture → in-range h-region SP, weak Kozak).
    """
    orig_len = len(orig_leader)
    # In-frame Met starts (aa offsets into the ORF); nt position is offset+3*k.
    atgs = [k for k in range(orig_len) if orig_leader[k] == "M"]

    chosen_k, source = 0, "contig"
    if orig_len > SP_LENGTH_MAX:
        # 1. Leaky-scan for the first (most upstream) adequate Met whose leader
        #    lands in-range. atgs ascend in k, so the resulting length descends;
        #    once it drops below the floor, no later Met can qualify → stop.
        for k in atgs:
            if k == 0:
                continue  # can't shorten the over-capture to itself
            new_len = orig_len - k
            if new_len < SP_LENGTH_MIN:
                break
            if new_len <= SP_LENGTH_MAX and (
                _kozak_score(contig_seq, offset + 3 * k) >= KOZAK_ADEQUATE
            ):
                chosen_k, source = k, "contig_kozak_reselected"
                break
        # 2. No adequate-Kozak start → best-available STRUCTURAL trim: an
        #    in-range nested Met WITH a hydrophobic h-region (a real SP whose
        #    start is weak-Kozak). Best-available Kozak, tie-break most
        #    downstream (the SP sits at the tail of the over-capture).
        if source == "contig":
            h_cands = [
                k
                for k in atgs
                if k != 0
                and SP_LENGTH_MIN <= orig_len - k <= SP_LENGTH_MAX
                and _has_h_region(orig_leader[k:])
            ]
            if h_cands:
                chosen_k = max(
                    h_cands,
                    key=lambda k: (_kozak_score(contig_seq, offset + 3 * k), k),
                )
                source = "contig_hregion_trimmed"

    leader_aa = orig_leader[chosen_k:]
    leader_nt = contig_seq[offset + 3 * chosen_k : offset + orig_len * 3]
    score = _kozak_score(contig_seq, offset + 3 * chosen_k)
    return leader_aa, leader_nt, score, source


def _extract_leader_from_contigs_single(
    row: pd.Series,
    sample_contigs: dict,
    result: dict,
    chain: str,
    *,
    species: str = "human",
    data_dir=None,
):
    """Extract + QC the signal peptide from contigs for a single chain (#263).

    Each supporting contig votes its (over-capture-corrected) leader; the
    consensus (most-common) leader is recorded with QC provenance:

    - ``{chain}_leader_qc`` — ``ok`` / ``weak_kozak_start`` (in-range, kept, but
      its start has weak Kozak context) / ``too_long`` / ``too_short`` /
      ``no_met`` / ``no_h_region`` / ``internal_stop`` / ``missing``.
    - ``{chain}_leader_source`` — ``contig`` (ORF start kept) or
      ``contig_kozak_reselected`` (over-capture auto-shortened to an in-range
      adequate-Kozak start).
    - ``{chain}_leader_len``, ``_leader_kozak_score`` (of the chosen start),
      ``_leader_support`` (n agreeing / n total).

    When no contig covers the VDJ the leader is simply left unset (no QC fields).
    """
    samples = str(row.get("samples", "")).split(";")

    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return

    contig_ids = str(row[contig_col]).split(";")
    vdj_aa = result.get(f"vdj_{chain}_aa", "")
    v_gene = row.get(f"{chain}_v_gene")

    leader_counter = Counter()
    # leader_aa -> (leader_nt, kozak_score, source). Metadata is consistent across
    # contigs voting the same leader_aa, so first sighting is fine. The germline
    # columns are added later by _record_leader's universal comparison.
    leader_meta: dict = {}
    n_total = 0

    for sample in samples:
        if sample not in sample_contigs:
            continue

        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue

            contig_seq = sample_contigs[sample][contig_id]
            translated, offset, ragged = find_longest_orf(contig_seq)

            if vdj_aa and vdj_aa in translated and offset is not None:
                orig_leader = translated.split(vdj_aa)[0]
                orig_len = len(orig_leader)
                n_total += 1
                # Layer 1: germline anchor (#267) decides the START — keep the
                # donor's native suffix whose length/start match the germline
                # V-leader. Falls back (None) when the gene/allele is absent or
                # the contig is too divergent → Kozak + h-region heuristic (#263).
                anchored = germline_anchor_leader(
                    orig_leader, v_gene, species=species, data_dir=data_dir
                )
                if anchored is not None:
                    leader_aa = anchored[0]
                    start_k = orig_len - len(leader_aa)
                    leader_nt = contig_seq[offset + 3 * start_k : offset + orig_len * 3]
                    score = _kozak_score(contig_seq, offset + 3 * start_k)
                    source = "contig_germline_anchored"
                else:
                    leader_aa, leader_nt, score, source = _kozak_correct_leader(
                        contig_seq, offset, orig_leader
                    )
                leader_counter[leader_aa] += 1
                leader_meta.setdefault(leader_aa, (leader_nt, score, source))

    if not leader_counter:
        return  # no contig covered the VDJ — leader simply absent

    leader_aa, support = leader_counter.most_common(1)[0]
    leader_nt, score, source = leader_meta[leader_aa]
    _record_leader(
        result, chain, leader_aa, leader_nt, source,
        v_gene=v_gene, kozak_score=score, support=f"{support}/{n_total}",
        species=species, data_dir=data_dir,
    )


# from_contig leader sources (the policy keeps/switches these; curated leaders
# the user chose explicitly are left alone except by an explicit force/secondary).
_FROM_CONTIG_SOURCES = frozenset(
    {"contig", "contig_kozak_reselected", "contig_hregion_trimmed", "contig_germline_anchored"}
)


def _resolve_named_leader(
    name: str, row, chain: str, *, species: str = "human", data_dir=None
) -> tuple[str, str | None, str] | None:
    """Resolve a leader directive → ``(aa, nt, source_tag)`` for ``chain``."""
    if name == "germline":
        g_aa = row.get(f"{chain}_germline_leader_aa")
        g_allele = row.get(f"{chain}_germline_allele")
        if isinstance(g_aa, str) and g_aa and isinstance(g_allele, str) and g_allele:
            gene, _, allele = g_allele.partition("*")
            return (
                g_aa,
                germline_leader_nt(gene, allele, species=species, data_dir=data_dir),
                "germline_reference_leader",
            )
        v_gene = row.get(f"{chain}_v_gene")
        cands = germline_vgene_leaders(v_gene, species=species, data_dir=data_dir)
        if cands:
            allele, _func, cand_aa = cands[0]
            gene_norm = normalize_vgene(v_gene) or str(v_gene)
            return (
                cand_aa,
                germline_leader_nt(gene_norm, allele, species=species, data_dir=data_dir),
                "germline_reference_leader",
            )
        return None
    hit = _lookup_default_leader(name)
    if hit is not None:
        canon_key, info = hit
        return info["aa"], info["dna"], f"curated:{canon_key}"
    return None


def _stash_preswitch(work: dict, chain: str) -> None:
    """Record a clone's donor leader variant BEFORE it's replaced (#276).

    When the policy switches (or a force overrides) a divergent contig leader,
    :func:`_substitute_chain_leader` refreshes ``{chain}_leader_vs_germline`` to
    describe the SHIPPED (germline/curated) leader — so the discarded donor
    variant would otherwise vanish from :func:`collect_germline_variants`. Stash
    it under ``{chain}_preswitch_*`` so the variant is still reported, tagged
    ``shipped="switched"``. No-op when the pre-switch leader already matched
    germline (nothing divergent to preserve).
    """
    diff = work.get(f"{chain}_leader_vs_germline")
    donor_aa = work.get(f"{chain}_leader_aa")
    germ_aa = work.get(f"{chain}_germline_leader_aa")
    # Prefer the recorded comparison; if it's absent (some paths don't stamp it)
    # but we have both sequences, derive it so the variant isn't lost on switch.
    if (not isinstance(diff, str) or not diff) and isinstance(donor_aa, str) and isinstance(germ_aa, str):
        diff = "identical" if donor_aa == germ_aa else characterize_divergence(donor_aa, germ_aa)
    if not isinstance(diff, str) or diff == "identical":
        return
    work[f"{chain}_preswitch_variant"] = diff
    work[f"{chain}_preswitch_donor_aa"] = work.get(f"{chain}_leader_aa")
    work[f"{chain}_preswitch_germline_aa"] = work.get(f"{chain}_germline_leader_aa")
    work[f"{chain}_preswitch_germline_allele"] = work.get(f"{chain}_germline_allele")
    work[f"{chain}_preswitch_germline_identity"] = work.get(f"{chain}_germline_identity")
    work[f"{chain}_preswitch_leader_source"] = work.get(f"{chain}_leader_source")


def _substitute_chain_leader(
    target: dict,
    orig,
    chain: str,
    new_aa: str,
    new_nt: str | None,
    source: str,
    *,
    species: str = "human",
    recode: CodonOptimization | None = None,
    data_dir=None,
) -> None:
    """In ``target`` (a row dict), replace ``chain``'s leader with ``new_aa`` and
    re-stitch the affected full-sequence columns (including recoded ``_optimized``)."""
    eff_recode = recode or CodonOptimization(codon_species=species)
    donor_aa = orig.get(f"{chain}_leader_aa")
    donor_aa = donor_aa if isinstance(donor_aa, str) else ""
    donor_nt = orig.get(f"{chain}_leader_nt")
    donor_nt = donor_nt if isinstance(donor_nt, str) else ""
    eff_new_nt = (
        new_nt.upper()
        if isinstance(new_nt, str) and new_nt
        else optimize_codons(new_aa, species=eff_recode.codon_species)
    )
    target[f"{chain}_leader_aa"] = new_aa
    target[f"{chain}_leader_nt"] = eff_new_nt
    target[f"{chain}_leader_len"] = len(new_aa)
    target[f"{chain}_leader_source"] = source
    target[f"{chain}_leader_qc"] = "substituted"
    feats = signal_peptide_features(new_aa)
    target[f"{chain}_sp_features_ok"] = feats["features_ok"]
    target[f"{chain}_sp_score"] = feats["score"]
    target[f"{chain}_sp_features"] = sp_features_summary(feats)
    cmp = germline_compare_leader(
        new_aa, orig.get(f"{chain}_v_gene"), species=species, data_dir=data_dir
    )
    if cmp is not None:
        allele, g_aa, identity, diff = cmp
        target[f"{chain}_germline_allele"] = f"{normalize_vgene(orig.get(f'{chain}_v_gene'))}*{allele}"
        target[f"{chain}_germline_leader_aa"] = g_aa
        target[f"{chain}_germline_identity"] = identity
        target[f"{chain}_leader_vs_germline"] = diff

    fa = orig.get(f"full_{chain}_aa")
    if isinstance(fa, str) and donor_aa and fa.startswith(donor_aa):
        target[f"full_{chain}_aa"] = new_aa + fa[len(donor_aa) :]
    elif isinstance(fa, str) and not donor_aa:
        target[f"full_{chain}_aa"] = new_aa + fa

    fn_blend = orig.get(f"full_{chain}_nt")
    if isinstance(fn_blend, str) and donor_nt and fn_blend.upper().startswith(donor_nt.upper()):
        target[f"full_{chain}_nt"] = (eff_new_nt + fn_blend[len(donor_nt) :]) if eff_new_nt else None
    elif isinstance(fn_blend, str) and not donor_nt:
        target[f"full_{chain}_nt"] = (eff_new_nt + fn_blend) if eff_new_nt else None

    fn_contig = orig.get(f"full_{chain}_nt_contig")
    if isinstance(fn_contig, str) and donor_nt and fn_contig.upper().startswith(donor_nt.upper()):
        target[f"full_{chain}_nt_contig"] = (eff_new_nt + fn_contig[len(donor_nt) :]) if eff_new_nt else None
    elif f"full_{chain}_nt_contig" in orig:
        target[f"full_{chain}_nt_contig"] = None

    def _as_str(v) -> str:
        return v if isinstance(v, str) else ""

    # Rebuild full_{chain}_nt_optimized (and recoded position tags) from the new
    # leader + VDJ + constant_nt_optimized so translation always matches full_{chain}_aa.
    if f"full_{chain}_nt_optimized" in orig or _as_str(orig.get(f"{chain}_constant_nt_optimized")):
        const_opt = _as_str(target.get(f"{chain}_constant_nt_optimized")) or _as_str(
            orig.get(f"{chain}_constant_nt_optimized")
        )
        const_aa = _as_str(target.get(f"{chain}_constant_aa")) or _as_str(
            orig.get(f"{chain}_constant_aa")
        )
        full_aa_now = _as_str(target.get(f"full_{chain}_aa"))
        if not full_aa_now:
            return
        donor_aa_now = (
            full_aa_now[: len(full_aa_now) - len(const_aa)]
            if (const_aa and full_aa_now.endswith(const_aa))
            else full_aa_now
        )
        sites = eff_recode.restriction_sites()
        vdj_nt_raw = (
            _as_str(target.get(f"vdj_{chain}_nt"))
            or _as_str(orig.get(f"vdj_{chain}_nt"))
            or _as_str(orig.get(f"VDJ_{chain}_nt"))
        ).strip()
        if vdj_nt_raw:
            vdj_aa_raw = (
                _as_str(target.get(f"vdj_{chain}_aa"))
                or _as_str(orig.get(f"vdj_{chain}_aa"))
                or _as_str(orig.get(f"VDJ_{chain}_aa"))
            ).strip()
            if vdj_aa_raw:
                vdj_nt_raw = vdj_nt_raw[: 3 * len(vdj_aa_raw)]
        if eff_recode.mode == "full" or not vdj_nt_raw:
            full_forbidden = (
                tuple(dict.fromkeys((*DEFAULT_FORBIDDEN_MOTIFS, *sites)))
                if sites
                else DEFAULT_FORBIDDEN_MOTIFS
            )
            donor_opt = optimize_codons(
                donor_aa_now,
                species=eff_recode.codon_species,
                forbidden_motifs=full_forbidden,
            )
            if sites or (eff_recode.max_codon_repeats and eff_recode.max_codon_repeats >= 1):
                donor_opt, _ = recode_codons(
                    donor_opt,
                    max_codon_repeats=eff_recode.max_codon_repeats,
                    restriction_sites=sites,
                    species=eff_recode.codon_species,
                    right_context=const_opt,
                )
            target[f"full_{chain}_nt_optimized"] = (donor_opt + const_opt) if eff_new_nt else None
            target[f"{chain}_leader_recoded_positions"] = ""
            target[f"{chain}_vdj_recoded_positions"] = ""
        elif eff_new_nt:
            native_prefix = (eff_new_nt + vdj_nt_raw).upper()
            recoded_prefix, changed = recode_codons(
                native_prefix,
                max_codon_repeats=eff_recode.max_codon_repeats,
                restriction_sites=sites,
                species=eff_recode.codon_species,
                right_context=const_opt,
            )
            target[f"full_{chain}_nt_optimized"] = recoded_prefix + const_opt
            vdj_start = len(eff_new_nt)

            def _fmt(o, _np=native_prefix, _rp=recoded_prefix):
                return f"{o}:{_np[o : o + 3]}>{_rp[o : o + 3]}"

            target[f"{chain}_leader_recoded_positions"] = ";".join(
                _fmt(o) for o in changed if o < vdj_start
            )
            target[f"{chain}_vdj_recoded_positions"] = ";".join(
                _fmt(o) for o in changed if o >= vdj_start
            )
        else:
            target[f"full_{chain}_nt_optimized"] = None


def _unanimous_divergent_genes(df: pd.DataFrame) -> set[tuple[str, str]]:
    """``(chain, gene)`` whose from_contig leader is unanimous across the gene's
    clones AND diverges from germline — the consistent-divergent case.

    Divergence is ANY mismatch (substitution OR indel), not length only: a
    unanimous single-residue variant (e.g. a TRBV13*03-style donor SNP that the
    270-allele reference subset lacks) is the donor's likely real allele just as
    much as a unanimous indel is, so it must be protected from the keep/switch
    policy too. The SP-soundness gate in :func:`apply_leader_policy` still decides
    whether such a leader is actually kept."""
    out: set[tuple[str, str]] = set()
    for chain in ("alpha", "beta"):
        v_col, aa_col, g_col = (
            f"{chain}_v_gene", f"{chain}_leader_aa", f"{chain}_germline_leader_aa"
        )
        if not {v_col, aa_col, g_col} <= set(df.columns):
            continue
        sub = df[df[aa_col].notna() & df[v_col].notna()]
        for gene, grp in sub.groupby(v_col):
            leaders = set(grp[aa_col].astype(str))
            germs = set(grp[g_col].dropna().astype(str))
            if len(leaders) == 1 and len(germs) == 1:
                leader, germ = next(iter(leaders)), next(iter(germs))
                if leader != germ:  # any unanimous divergence (substitution or indel)
                    out.add((chain, str(gene)))
    return out


def _protected_divergent_leaders(
    df: pd.DataFrame, *, min_support: int = 2
) -> set[tuple[str, str, str]]:
    """``(chain, gene, leader_aa)`` for divergent leaders that look like a real
    donor allele rather than a one-off artifact — eligible to KEEP (the SP gate
    in :func:`apply_leader_policy` still decides the actual keep).

    A divergent leader qualifies when EITHER:

    * it is the gene's **unanimous** leader (whole-gene consistent — the original
      #270 case; covers single-clone genes), OR
    * it is carried by ``>= min_support`` clones of that gene (**het-aware**, #274):
      a heterozygous gene's minority allele has independent multi-clone support,
      so a TRBV6-2 G9A-style second allele is the donor's real allele even though
      the gene as a whole isn't unanimous. Only a singleton divergent leader in a
      non-unanimous gene (likely a sequencing artifact) is left unprotected.

    Divergence is judged per row (leader vs *that clone's* closest germline allele),
    so a het gene whose alleles resolve to different germline leaders is handled
    correctly. Returns ``leader_aa``-keyed tuples (not just gene) so only the
    supported allele is protected, not every leader the gene carries.
    """
    out: set[tuple[str, str, str]] = set()
    # Count support by distinct CLONE, not construct row: dual-α expansion
    # duplicates the β row (same β leader, two α variants), so a row-count would
    # promote a single dual-α clone's β leader to "multi-clone support" and keep
    # an artifact. `selected_clone` is shared across a clone's dual-α twins; fall
    # back to CDR3ab, then the row index, when it's absent.
    clone_col = next(
        (c for c in ("selected_clone", "CDR3ab") if c in df.columns), None
    )
    for chain in ("alpha", "beta"):
        v_col, aa_col, g_col = (
            f"{chain}_v_gene", f"{chain}_leader_aa", f"{chain}_germline_leader_aa"
        )
        if not {v_col, aa_col, g_col} <= set(df.columns):
            continue
        sub = df[df[aa_col].notna() & df[v_col].notna()]
        for gene, grp in sub.groupby(v_col):
            leaders = grp[aa_col].astype(str)
            keys = grp[clone_col].astype(str) if clone_col else grp.index.to_series()
            # distinct clones carrying each leader (de-duped over dual-α rows)
            support = keys.groupby(leaders).nunique()
            whole_gene_unanimous = leaders.nunique() == 1
            for _, r in grp.iterrows():
                leader = str(r[aa_col])
                germ = r[g_col]
                if not isinstance(germ, str) or leader == germ:
                    continue  # not germline-comparable, or matches germline (not divergent)
                if whole_gene_unanimous or support.get(leader, 0) >= min_support:
                    out.add((chain, str(gene), leader))
    return out


def _validate_leader_directive(val: str | None, param_name: str, *, allow_keep: bool = False) -> None:
    if val is None:
        return
    s = str(val).strip()
    if s == "germline" or (allow_keep and s == "keep"):
        return
    if _lookup_default_leader(s) is not None:
        return
    valid = ["germline", *(["keep"] if allow_keep else []), *DEFAULT_LEADERS.keys()]
    raise TCRsiftValidationError(
        f"Unknown {param_name}={val!r}; expected one of {valid}"
    )


def apply_leader_policy(
    df: pd.DataFrame,
    linker: str | None,
    *,
    leader_fallback: str = "germline",
    force_alpha: str | None = None,
    force_beta: str | None = None,
    secondary_alpha: str | None = None,
    secondary_beta: str | None = None,
    verbose: bool = False,
    species: str = "human",
    recode: CodonOptimization | None = None,
    data_dir=None,
) -> pd.DataFrame:
    """Decide the shipped signal peptide for every clone (von Heijne policy, #270)."""
    _validate_leader_directive(leader_fallback, "leader_fallback", allow_keep=True)
    _validate_leader_directive(force_alpha, "force_alpha")
    _validate_leader_directive(force_beta, "force_beta")
    _validate_leader_directive(secondary_alpha, "secondary_alpha")
    _validate_leader_directive(secondary_beta, "secondary_beta")

    protected = _protected_divergent_leaders(df)
    force = {"alpha": force_alpha, "beta": force_beta}
    secondary = {"alpha": secondary_alpha, "beta": secondary_beta}
    orig_index = df.index
    df = df.reset_index(drop=True)
    n_forced = n_switched = 0

    for idx, row in df.iterrows():
        work = row.to_dict()
        changed = False
        for chain in ("alpha", "beta"):
            aa = work.get(f"{chain}_leader_aa")
            has_leader = isinstance(aa, str) and bool(aa)
            if force[chain]:
                full_chain_val = work.get(f"full_{chain}_aa")
                if not isinstance(full_chain_val, str) or not full_chain_val:
                    continue
                res = _resolve_named_leader(
                    force[chain], work, chain, species=species, data_dir=data_dir
                )
                if res is not None:
                    new_aa, new_nt, src = res
                    if has_leader:
                        _stash_preswitch(work, chain)
                    _substitute_chain_leader(
                        work, work, chain, new_aa, new_nt, f"forced_{src}",
                        species=species, recode=recode, data_dir=data_dir,
                    )
                    changed = True
                    n_forced += 1
                continue
            if not has_leader:
                if (
                    work.get(f"{chain}_leader_source") == "missing_from_contig"
                    and isinstance(work.get(f"full_{chain}_aa"), str)
                    and work.get(f"full_{chain}_aa")
                    and leader_fallback != "keep"
                ):
                    res = _resolve_named_leader(
                        leader_fallback, work, chain, species=species, data_dir=data_dir
                    )
                    if res is None:
                        curated_fb = (
                            ("mCD28" if chain == "alpha" else "mCD8B")
                            if species == "mouse"
                            else ("CD28" if chain == "alpha" else "CD8A")
                        )
                        res = _resolve_named_leader(
                            curated_fb, work, chain, species=species, data_dir=data_dir
                        )
                    if res is not None:
                        new_aa, new_nt, src = res
                        _substitute_chain_leader(
                            work, work, chain, new_aa, new_nt, src,
                            species=species, recode=recode, data_dir=data_dir,
                        )
                        changed = True
                        n_switched += 1
                continue
            if work.get(f"{chain}_leader_source") not in _FROM_CONTIG_SOURCES:
                continue
            if leader_fallback == "keep":
                continue
            germ_aa = work.get(f"{chain}_germline_leader_aa")
            ident = work.get(f"{chain}_germline_identity")
            sp_ok = bool(work.get(f"{chain}_sp_features_ok"))
            matches_germline = (
                isinstance(germ_aa, str) and isinstance(ident, (int, float))
                and ident >= 0.999
            )
            divergent = isinstance(germ_aa, str) and not matches_germline
            gene = work.get(f"{chain}_v_gene")
            leader_now = str(work.get(f"{chain}_leader_aa"))
            is_consistent = (chain, str(gene), leader_now) in protected
            keep = matches_germline or (sp_ok and (not divergent or is_consistent))
            if keep:
                continue
            res = _resolve_named_leader(
                leader_fallback, work, chain, species=species, data_dir=data_dir
            )
            if res is None:
                continue
            new_aa, new_nt, src = res
            _stash_preswitch(work, chain)
            _substitute_chain_leader(
                work, work, chain, new_aa, new_nt, src,
                species=species, recode=recode, data_dir=data_dir,
            )
            changed = True
            n_switched += 1
            if verbose:
                logger.warning(
                    "  %s leader rejected (sp_ok=%s, divergent=%s, consistent=%s) → "
                    "switched to %s", chain, sp_ok, divergent, is_consistent, src,
                )
        if changed:
            for k, v in work.items():
                if k not in df.columns:
                    df[k] = pd.Series([None] * len(df), index=df.index, dtype=object)
                elif isinstance(v, str) and pd.api.types.is_numeric_dtype(df[k].dtype):
                    df[k] = df[k].astype(object)
                df.at[idx, k] = v

    twins, orig_idx = [], []
    if any(secondary.values()):
        for idx, row in df.iterrows():
            sec = []
            for chain in ("alpha", "beta"):
                if not secondary[chain]:
                    continue
                cur_aa = row.get(f"{chain}_leader_aa")
                full_chain_val = row.get(f"full_{chain}_aa")
                if not isinstance(cur_aa, str) or not isinstance(full_chain_val, str) or not full_chain_val:
                    continue
                res = _resolve_named_leader(
                    secondary[chain], row, chain, species=species, data_dir=data_dir
                )
                if res is not None and res[0] != cur_aa:
                    sec.append((chain, secondary[chain], res))
            if not sec:
                continue
            twin = row.to_dict()
            for chain, name, (new_aa, new_nt, src) in sec:
                _substitute_chain_leader(
                    twin, row, chain, new_aa, new_nt, f"secondary_{src}",
                    species=species, recode=recode, data_dir=data_dir,
                )
            twin["leader_variant"] = "secondary:" + "+".join(s[1] for s in sec)
            twins.append(twin)
            orig_idx.append(idx)

    if twins:
        if "leader_variant" not in df.columns:
            df["leader_variant"] = pd.NA
        df.loc[orig_idx, "leader_variant"] = "primary"
        df = pd.concat([df, pd.DataFrame(twins)], ignore_index=True)
    else:
        df.index = orig_index

    # Rebuild single_chain (and refresh synthesis QC if already present) wherever leaders changed.
    if n_forced or n_switched or twins:
        if linker and {"full_beta_aa", "full_alpha_aa"} <= set(df.columns):
            df = _add_single_chain(df, linker, recode=recode)
        if "synth_gc_ok" in df.columns:
            df = add_synthesis_qc(df)
    if verbose and (n_forced or n_switched or twins):
        logger.info(
            "  Leader policy: %d forced, %d switched to fallback, %d secondary constructs",
            n_forced, n_switched, len(twins),
        )
    return df

