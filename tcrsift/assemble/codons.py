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

"""Codon tables, synonymous recoding, DNA translation, and 2A single-chain linking."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)

# Standard codon table
CODON_TABLE = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGC": "C",
    "TGT": "C",
    "TGA": "*",
    "TGG": "W",
}

# Self-cleaving 2A peptide linkers
# These are viral 2A sequences that cause ribosomal skipping during translation
LINKERS = {
    "T2A": {
        "dna": "GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGCCCG",
        "aa": "EGRGSLLTCGDVEENPGP",
        "source": "Thosea asigna virus",
    },
    "P2A": {
        "dna": "GGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT",
        "aa": "GSGATNFSLLKQAGDVEENPGP",
        "source": "Porcine teschovirus-1",
    },
    "E2A": {
        "dna": "CAGTGTACTAATTATGCTCTCTTGAAATTGGCTGGAGATGTTGAGAGCAACCCAGGTCCC",
        "aa": "QCTNYALLKLAGDVESNPGP",
        "source": "Equine rhinitis A virus",
    },
    "F2A": {
        "dna": "GTGAAACAGACTTTGAATTTTGACCTTCTCAAGTTGGCGGGAGACGTGGAGTCCAACCCAGGGCCC",
        "aa": "VKQTLNFDLLKLAGDVESNPGP",
        "source": "Foot-and-mouth disease virus",
    },
}

# Backwards compatibility aliases
T2A_LINKER_DNA = LINKERS["T2A"]["dna"]
T2A_LINKER_AA = LINKERS["T2A"]["aa"]
P2A_LINKER_DNA = LINKERS["P2A"]["dna"]
P2A_LINKER_AA = LINKERS["P2A"]["aa"]

# Default leader/signal peptide sequences for TCR expression
# These can be used when contig FASTA files are not available


# Named Type-II restriction sites scanned by the synthesis-hazard QC (#206).
# (The codon optimizer already AVOIDS these; this QC flags any that survive —
# e.g. in the *unoptimized* CDS, or across a leader/constant boundary.)
SYNTHESIS_RESTRICTION_SITES: dict[str, str] = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "XhoI": "CTCGAG", "NotI": "GCGGCCGC",
    "NheI": "GCTAGC", "HindIII": "AAGCTT", "SalI": "GTCGAC", "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC", "NcoI": "CCATGG", "NdeI": "CATATG", "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
}



# Codon tables for human and mouse (CCDS / Kazusa CAI-ordered).
HUMAN_PREFERRED_CODONS: dict[str, str] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA",
}

MOUSE_PREFERRED_CODONS: dict[str, str] = {
    "A": "GCC", "R": "AGA", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA",
}

PREFERRED_CODONS_BY_SPECIES: dict[str, dict[str, str]] = {
    "human": HUMAN_PREFERRED_CODONS,
    "mouse": MOUSE_PREFERRED_CODONS,
}

HUMAN_CODON_ALTERNATIVES: dict[str, list[str]] = {
    "A": ["GCC", "GCT", "GCA", "GCG"],
    "R": ["CGG", "AGG", "CGC", "AGA", "CGA", "CGT"],
    "N": ["AAC", "AAT"],
    "D": ["GAC", "GAT"],
    "C": ["TGC", "TGT"],
    "Q": ["CAG", "CAA"],
    "E": ["GAG", "GAA"],
    "G": ["GGC", "GGA", "GGG", "GGT"],
    "H": ["CAC", "CAT"],
    "I": ["ATC", "ATT", "ATA"],
    "L": ["CTG", "CTC", "TTG", "CTT", "CTA", "TTA"],
    "K": ["AAG", "AAA"],
    "M": ["ATG"],
    "F": ["TTC", "TTT"],
    "P": ["CCC", "CCT", "CCA", "CCG"],
    "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
    "T": ["ACC", "ACA", "ACT", "ACG"],
    "W": ["TGG"],
    "Y": ["TAC", "TAT"],
    "V": ["GTG", "GTC", "GTT", "GTA"],
    "*": ["TAA", "TGA", "TAG"],
}

MOUSE_CODON_ALTERNATIVES: dict[str, list[str]] = {
    "A": ["GCC", "GCT", "GCA", "GCG"],
    "R": ["AGA", "AGG", "CGG", "CGC", "CGA", "CGT"],
    "N": ["AAC", "AAT"],
    "D": ["GAC", "GAT"],
    "C": ["TGC", "TGT"],
    "Q": ["CAG", "CAA"],
    "E": ["GAG", "GAA"],
    "G": ["GGC", "GGA", "GGG", "GGT"],
    "H": ["CAC", "CAT"],
    "I": ["ATC", "ATT", "ATA"],
    "L": ["CTG", "CTC", "CTT", "TTG", "CTA", "TTA"],
    "K": ["AAG", "AAA"],
    "M": ["ATG"],
    "F": ["TTC", "TTT"],
    "P": ["CCC", "CCT", "CCA", "CCG"],
    "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
    "T": ["ACC", "ACA", "ACT", "ACG"],
    "W": ["TGG"],
    "Y": ["TAC", "TAT"],
    "V": ["GTG", "GTC", "GTT", "GTA"],
    "*": ["TAA", "TGA", "TAG"],
}

CODON_ALTERNATIVES_BY_SPECIES: dict[str, dict[str, list[str]]] = {
    "human": HUMAN_CODON_ALTERNATIVES,
    "mouse": MOUSE_CODON_ALTERNATIVES,
}


# Default forbidden motif set for ``optimize_codons``. Four categories:
# (1) ≥5-mer mononucleotide runs that hurt expression / stall ribosomes
#     and cause synthesis QC failures.
# (2) Type-II restriction enzyme recognition sites that are widely used
#     in cloning workflows and routinely scrubbed out of synthesized
#     constructs.
# (3) Cryptic polyadenylation signals. Premature polyA causes 3' UTR
#     truncation and unstable transcripts (Tian & Manley 2013, Gruber
#     & Zavolan 2019). AATAAA is the canonical signal (~70% of human
#     polyA sites); ATTAAA is the most common variant (~10–15%).
# (4) Canonical 5' splice donor (intron donor) consensus. Adding these
#     prevents the synthesized transcript from being mis-spliced when
#     transcribed in mammalian cells. ``GTAAGT`` and ``GTGAGT`` match
#     the strongest donor consensus (Sheth et al. 2006). Note: the
#     intron acceptor (``CAG``) is too short to forbid without
#     blocking every Q codon (CAG = CAI-best for Q), so we accept the
#     residual risk on the acceptor side — donors are the higher-
#     leverage filter empirically.
DEFAULT_FORBIDDEN_MOTIFS: tuple[str, ...] = (
    # Mononucleotide runs (5+).
    "AAAAA",
    "TTTTT",
    "GGGGG",
    "CCCCC",
    # Common Type-II / Type-IIS restriction sites (both strands for asymmetric Type-IIS).
    "GAATTC",    # EcoRI
    "GGATCC",    # BamHI
    "AAGCTT",    # HindIII
    "CATATG",    # NdeI
    "GCTAGC",    # NheI
    "GCGGCCGC",  # NotI
    "GTCGAC",    # SalI
    "TCTAGA",    # XbaI
    "CTCGAG",    # XhoI
    "GGTACC",    # KpnI
    "CCATGG",    # NcoI
    "GGTCTC",    # BsaI (+)
    "GAGACC",    # BsaI (-)
    "CGTCTC",    # BsmBI (+)
    "GAGACG",    # BsmBI (-)
    # Cryptic polyadenylation signals.
    "AATAAA",
    "ATTAAA",
    # Canonical 5' splice donors (cryptic splicing prevention).
    "GTAAGT",
    "GTGAGT",
)


# GC-content target for sliding-window optimization. Mammalian coding
# regions are typically 45-55% GC; expression suffers at the extremes
# (Kudla et al. 2006, Gustafsson et al. 2004): too-high GC drives
# strong mRNA secondary structure, too-low GC drives instability
# (AU-rich element ARE-like behavior). The range below is permissive
# enough to leave most codon choices unrestricted but kicks in to
# break ties when alternatives are otherwise equivalent.
DEFAULT_GC_TARGET_WINDOW: int = 60
DEFAULT_GC_TARGET_RANGE: tuple[float, float] = (0.40, 0.65)


def back_translate(aa: str, *, species: str = "human") -> str:
    """Reverse-translate a polypeptide to DNA via preferred codons for ``species``."""
    preferred = PREFERRED_CODONS_BY_SPECIES.get(species, HUMAN_PREFERRED_CODONS)
    fallback = "NNN"
    out = [preferred.get(r, fallback) for r in aa]
    return "".join(out)


def optimize_codons(
    aa_sequence: str,
    *,
    prefix_nt: str = "",
    forbidden_motifs: tuple[str, ...] = DEFAULT_FORBIDDEN_MOTIFS,
    context_window: int = 8,
    gc_target_window: int = DEFAULT_GC_TARGET_WINDOW,
    gc_target_range: tuple[float, float] = DEFAULT_GC_TARGET_RANGE,
    species: str = "human",
) -> str:
    """Back-translate ``aa_sequence`` to DNA with motif avoidance and
    GC-content balancing for ``species`` ('human' or 'mouse')."""
    if gc_target_range[0] > gc_target_range[1]:
        raise ValueError(
            f"gc_target_range must be (low, high); got {gc_target_range}"
        )
    codon_alts = CODON_ALTERNATIVES_BY_SPECIES.get(species, HUMAN_CODON_ALTERNATIVES)
    preferred_codons = PREFERRED_CODONS_BY_SPECIES.get(species, HUMAN_PREFERRED_CODONS)
    gc_balancing_active = gc_target_range != (0.0, 1.0)
    prefix_up = str(prefix_nt).upper() if prefix_nt else ""
    motifs_up = tuple(str(m).upper() for m in forbidden_motifs if m)
    eff_window = max(context_window, max((len(m) for m in motifs_up), default=0))

    def _creates_motif(buf: str, codon: str) -> bool:
        for motif in motifs_up:
            k = len(motif) - 1
            tail = buf[-k:] if k > 0 else ""
            if motif in (tail + codon):
                return True
        return False

    def _has_valid_continuation(pos: int, new_buffer: str, depth: int) -> bool:
        if depth <= 0 or pos >= len(aa_sequence):
            return True
        next_alts = codon_alts.get(aa_sequence[pos])
        if not next_alts:
            return True
        for next_codon in next_alts:
            if not _creates_motif(new_buffer, next_codon):
                if depth == 1 or _has_valid_continuation(
                    pos + 1, (new_buffer + next_codon)[-eff_window:], depth - 1
                ):
                    return True
        return False

    def _gc_in_range(codon: str, gc_buffer: str) -> bool:
        if not gc_balancing_active:
            return True
        window = (gc_buffer + codon)[-gc_target_window:]
        if not window:
            return True
        gc = window.count("G") + window.count("C")
        frac = gc / len(window)
        return gc_target_range[0] <= frac <= gc_target_range[1]

    nt_parts: list[str] = []
    buffer = prefix_up[-eff_window:] if prefix_up else ""
    gc_buffer = prefix_up[-gc_target_window:] if prefix_up else ""
    max_motif_len = max((len(m) for m in motifs_up), default=0)
    lookahead_depth = max(1, (max_motif_len - 1 + 2) // 3) if motifs_up else 0
    for i, aa in enumerate(aa_sequence):
        alternatives = codon_alts.get(aa)
        if not alternatives:
            codon = preferred_codons.get(aa, "NNN")
            nt_parts.append(codon)
            buffer = (buffer + codon)[-eff_window:]
            gc_buffer = (gc_buffer + codon)[-gc_target_window:]
            continue
        # Filter 1 (hard): motif-free across the prefix+codon boundary
        motif_free = [
            c for c in alternatives
            if not _creates_motif(buffer, c)
        ]
        if not motif_free:
            # Every alternative triggers a motif (e.g. G after GGGGG); pick the
            # candidate that minimizes the resulting homopolymer/motif overlap,
            # breaking ties by highest CAI (first in `alternatives`).
            def _badness(cand: str) -> int:
                joined = buffer + cand
                score = sum(1 for m in motifs_up if m in joined)
                # Penalize extending a trailing homopolymer run further.
                tail_ch = cand[-1]
                run = 0
                for ch in reversed(joined):
                    if ch == tail_ch:
                        run += 1
                    else:
                        break
                return score * 100 + run

            chosen = min(alternatives, key=_badness)
        else:
            # Filter 2 (soft): bounded multi-codon lookahead so 6–8 bp sites
            # spanning 3 codons (e.g. mouse TRAC T-M-E -> ACC-ATG-GAG = NcoI)
            # are avoided proactively. Fall back to motif_free if empty.
            if lookahead_depth > 0 and i + 1 < len(aa_sequence):
                lookahead_valid = [
                    c for c in motif_free
                    if _has_valid_continuation(
                        i + 1, (buffer + c)[-eff_window:], min(lookahead_depth, 2)
                    )
                ] or motif_free
            else:
                lookahead_valid = motif_free
            # Filter 3 (soft): GC-in-range. Refine; fall back to
            # lookahead_valid if every option is out of range.
            gc_balanced = [
                c for c in lookahead_valid if _gc_in_range(c, gc_buffer)
            ] or lookahead_valid
            # CAI tiebreak: first remaining candidate is highest-CAI.
            chosen = gc_balanced[0]
        nt_parts.append(chosen)
        buffer = (buffer + chosen)[-eff_window:]
        gc_buffer = (gc_buffer + chosen)[-gc_target_window:]
    return "".join(nt_parts)


def _revcomp(s: str) -> str:
    return s.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def recode_codons(
    nt: str,
    *,
    max_codon_repeats: int = 1,
    restriction_sites: tuple[str, ...] = (),
    species: str = "human",
    right_context: str = "",
) -> tuple[str, list[int]]:
    """Minimal synonymous recoding of an in-frame coding NT (#293)."""
    if not isinstance(nt, str) or len(nt) < 3:
        return (nt.upper() if isinstance(nt, str) else ""), []
    nt = nt.upper()
    rc_upper = str(right_context).upper() if right_context else ""
    codon_alts = CODON_ALTERNATIVES_BY_SPECIES.get(species, HUMAN_CODON_ALTERNATIVES)
    forbidden: set[str] = set()
    for site in restriction_sites:
        s = str(site).upper()
        if s:
            forbidden.add(s)
            forbidden.add(_revcomp(s))
    # Avoid newly injecting common restriction sites during repeat-breaking swaps
    # even when the caller did not pass explicit `restriction_sites`.
    soft_forbidden: set[str] = set()
    for site in SYNTHESIS_RESTRICTION_SITES.values():
        s = site.upper()
        soft_forbidden.add(s)
        soft_forbidden.add(_revcomp(s))
    max_site = max(
        (len(s) for s in (forbidden | soft_forbidden)), default=0
    )
    repeat_active = bool(max_codon_repeats and max_codon_repeats >= 1)

    def _repeats_ok(prefix_codons: list[str], cand: str) -> bool:
        if not repeat_active:
            return True
        run = 1
        for pc in reversed(prefix_codons):
            if pc == cand:
                run += 1
            else:
                break
        return run <= max_codon_repeats

    def _sites_ok(trail_buf: str, cand: str, lookahead: str = "", *, check_set: set[str] | None = None) -> bool:
        target_set = forbidden if check_set is None else check_set
        if not target_set:
            return True
        for f in target_set:
            k = len(f) - 1
            left = trail_buf[-k:] if k > 0 else ""
            right = lookahead[:k] if k > 0 else ""
            if f in (left + cand + right):
                return False
        return True

    out: list[str] = []
    n_codons = len(nt) // 3
    orig_codons = [nt[3 * i : 3 * i + 3] for i in range(n_codons)]
    lookahead_tail = rc_upper[: max(0, max_site - 1)] if max_site > 1 else ""

    for i, orig in enumerate(orig_codons):
        aa = CODON_TABLE.get(orig)
        alts = codon_alts.get(aa, [])
        candidates = [orig, *[c for c in alts if c != orig]]
        trail = "".join(out)[-(max_site - 1) :] if (max_site > 1 and out) else ""
        # Look ahead into upcoming native codons (plus right_context) so a swap
        # does not create a restriction site spanning codon i -> i+1 (e.g. GCC->GCT + AGC = NheI).
        upcoming = "".join(orig_codons[i + 1 : i + 4]) + lookahead_tail
        la = upcoming[: max(0, max_site - 1)] if max_site > 1 else ""
        chosen = None
        for cand in candidates:
            if not (_repeats_ok(out, cand) and _sites_ok(trail, cand, la)):
                continue
            # When swapping away from the native codon, avoid introducing a new
            # restriction site from SYNTHESIS_RESTRICTION_SITES.
            if cand != orig and not _sites_ok(trail, cand, la, check_set=soft_forbidden):
                continue
            chosen = cand
            break
        if chosen is None:
            for cand in candidates:
                if _repeats_ok(out, cand) and _sites_ok(trail, cand, la):
                    chosen = cand
                    break
        if chosen is None and i > 0:
            # 1-codon backtracking: needed when codon i has no alternative (e.g.
            # ATG/Met or TGG/Trp at the 3' end of NcoI CCATGG or NdeI CATATG).
            prev_orig = orig_codons[i - 1]
            prev_aa = CODON_TABLE.get(prev_orig)
            prev_alts = [prev_orig, *[c for c in codon_alts.get(prev_aa, []) if c != prev_orig]]
            prefix_before_prev = out[:-1]
            trail_before_prev = (
                "".join(prefix_before_prev)[-(max_site - 1) :]
                if (max_site > 1 and prefix_before_prev)
                else ""
            )
            for prev_cand in prev_alts:
                if prev_cand == out[-1]:
                    continue
                if not (_repeats_ok(prefix_before_prev, prev_cand) and _sites_ok(trail_before_prev, prev_cand)):
                    continue
                new_prefix = [*prefix_before_prev, prev_cand]
                new_trail = "".join(new_prefix)[-(max_site - 1) :] if max_site > 1 else ""
                for cand in candidates:
                    if _repeats_ok(new_prefix, cand) and _sites_ok(new_trail, cand, la):
                        out[-1] = prev_cand
                        chosen = cand
                        break
                if chosen is not None:
                    break
        if chosen is None:
            chosen = orig
        out.append(chosen)

    changed = [3 * i for i, (o, c) in enumerate(zip(orig_codons, out)) if o != c]
    return "".join(out) + nt[3 * n_codons :], changed


@dataclass(frozen=True)
class CodonOptimization:
    """How ``full_{chain}_nt_optimized`` is generated (#293).

    * ``mode="focal"`` (default) — keep the donor's native leader+VDJ codons and
      make only the minimal synonymous swaps that clear a constraint
      (:func:`recode_codons`); the canonical constant stays codon-optimized.
      Conservative: the synthesis CDS stays byte-faithful to the contig wherever
      it can.
    * ``mode="full"`` — recode the entire ORF with :func:`optimize_codons`
      (preferred codons + motif/GC avoidance), like the constant region.

    Constraints (focal mode): ``max_codon_repeats`` (default 1 = no immediate
    identical-codon repeat) and ``avoid_enzymes`` (restriction enzyme names from
    :data:`SYNTHESIS_RESTRICTION_SITES`, or raw ACGT recognition sites; default
    none).
    """

    mode: str = "focal"
    max_codon_repeats: int = 1
    avoid_enzymes: tuple[str, ...] = ()
    codon_species: str = "human"

    def restriction_sites(self) -> tuple[str, ...]:
        sites: list[str] = []
        upper_map = {k.upper(): v for k, v in SYNTHESIS_RESTRICTION_SITES.items()}
        for e in self.avoid_enzymes:
            s = str(e).strip()
            if s in SYNTHESIS_RESTRICTION_SITES:
                sites.append(SYNTHESIS_RESTRICTION_SITES[s])
            elif s.upper() in upper_map:
                sites.append(upper_map[s.upper()])
            elif s and set(s.upper()) <= set("ACGT"):
                sites.append(s.upper())  # raw recognition site
            elif s:
                logger.warning("avoid_enzymes: %r is not a known enzyme or ACGT site; skipping", s)
        return tuple(sites)


def stop_codons_nt(stop_codons: tuple[str, ...]) -> str:
    """Validate and concatenate stop codons for appending to a CDS.

    Each entry must be a 3-character DNA codon coding for a stop
    (``TAA`` / ``TAG`` / ``TGA``). Duplicates are allowed but the
    canonical "two non-redundant stops" pattern from cloning practice
    uses two DIFFERENT stops (different release-factor recognition);
    we don't enforce non-redundancy here but the default in
    :func:`assemble_full_sequences` is ``("TAA", "TGA")``.
    """
    if not stop_codons:
        return ""
    valid = {"TAA", "TAG", "TGA"}
    bad = [c for c in stop_codons if c not in valid]
    if bad:
        raise ValueError(
            f"stop_codons entries must be one of {sorted(valid)}; "
            f"got: {bad}"
        )
    return "".join(stop_codons)



def translate_dna(dna_seq: str) -> tuple[str, str]:
    """
    Translate DNA sequence to amino acids.

    Returns
    -------
    tuple
        (amino_acid_sequence, ragged_3p_nucleotides)
    """
    dna_seq = str(dna_seq).upper()
    seq_len = len(dna_seq)
    seq_len_trimmed = (seq_len // 3) * 3

    if seq_len != seq_len_trimmed:
        ragged_nt = dna_seq[seq_len_trimmed:]
        dna_seq = dna_seq[:seq_len_trimmed]
    else:
        ragged_nt = ""

    aa_seq = "".join([CODON_TABLE.get(dna_seq[i : i + 3], "X") for i in range(0, len(dna_seq), 3)])

    # Stop at first stop codon
    if "*" in aa_seq:
        ragged_nt = ""
        aa_seq = aa_seq[: aa_seq.index("*")]

    return aa_seq, ragged_nt


def find_longest_orf(dna_seq: str) -> tuple[str, int, str]:
    """
    Find and translate the longest open reading frame.

    Returns
    -------
    tuple
        (amino_acid_sequence, start_offset, ragged_3p_nucleotides)
    """
    dna_seq = str(dna_seq).upper()
    start_positions = [i for i in range(len(dna_seq)) if dna_seq[i : i + 3] == "ATG"]

    longest_aa = ""
    longest_offset = 0
    longest_ragged = ""

    for start in start_positions:
        subseq = dna_seq[start:]
        aa, ragged = translate_dna(subseq)
        if len(aa) > len(longest_aa):
            longest_aa = aa
            longest_offset = start
            longest_ragged = ragged

    return longest_aa, longest_offset, longest_ragged


def _add_single_chain(df: pd.DataFrame, linker: str) -> pd.DataFrame:
    """Add single-chain construct (β-linker-α) including the NT triad.

    Emits four columns:

    * ``single_chain_aa``
    * ``single_chain_nt`` — built from the blend ``full_*_nt``. The
      column most callers reach for; donor-faithful at the J→C
      boundary, codon-optimized canonical otherwise.
    * ``single_chain_nt_optimized`` — built from
      ``full_*_nt_optimized``. The ready-to-order synthesis construct;
      same bytes for every donor with the same picked allele.

    No ``single_chain_nt_contig`` is emitted: the stitched contig construct is
    truncated at contig coverage (not a CDS, not synthesis-ready) and carries no
    QC signal beyond the per-chain ``{chain}_constant_nt_contig`` /
    ``full_{chain}_nt_contig`` provenance columns it was built from. Inspect
    those for the donor's raw bytes / contig backing.

    The β CDS in each variant has ALL trailing stops stripped before
    the linker — the 2A construct requires a continuous ORF across
    β → linker → α. The dual-stop default (``TAATGA``) made the
    one-codon-strip in 2.3 incorrect; the loop here is the fix.
    """
    if linker.upper() in LINKERS:
        linker_info = LINKERS[linker.upper()]
        linker_aa = linker_info["aa"]
        linker_nt = linker_info["dna"]
    else:
        # Custom linker sequence provided as amino acids; no NT.
        linker_aa = linker
        linker_nt = ""

    def strip_stop(seq):
        if not isinstance(seq, str):
            return None
        return seq.rstrip("*")

    beta_aa = df["full_beta_aa"].apply(strip_stop)
    alpha_aa = df["full_alpha_aa"].where(df["full_alpha_aa"].apply(lambda x: isinstance(x, str)))

    single_chain = beta_aa.fillna("") + linker_aa + alpha_aa.fillna("")
    missing_mask = df["full_beta_aa"].isna() | df["full_alpha_aa"].isna()
    df["single_chain_aa"] = single_chain.where(~missing_mask, pd.NA)

    def strip_stop_codon_dna(seq):
        if not isinstance(seq, str):
            return None
        # Strip ALL trailing stops. Single-stop strip pre-2.4 broke
        # frame across the linker once dual-stop became the default.
        while len(seq) >= 3 and seq[-3:] in {"TAA", "TAG", "TGA"}:
            seq = seq[:-3]
        return seq

    if linker_nt:
        # Build the triad. Each variant uses its matching beta/alpha
        # input column; if either input is missing or contains None
        # (contig view when contig didn't cover the C region) the
        # output is None — splicing a None into a 2A construct would
        # be silently wrong.
        # No ``_contig`` single-chain variant: the stitched contig construct is
        # truncated at contig coverage (not a CDS, not synthesis-ready) and adds
        # nothing over the per-chain ``{chain}_constant_nt_contig`` /
        # ``full_{chain}_nt_contig`` provenance columns it's built from.
        for suffix, out_col in (
            ("", "single_chain_nt"),
            ("_optimized", "single_chain_nt_optimized"),
        ):
            beta_col = f"full_beta_nt{suffix}"
            alpha_col = f"full_alpha_nt{suffix}"
            if beta_col not in df.columns or alpha_col not in df.columns:
                continue

            def _build_row(row, b=beta_col, a=alpha_col):
                bv = row.get(b)
                av = row.get(a)
                if not isinstance(bv, str) or not isinstance(av, str):
                    return None
                return strip_stop_codon_dna(bv) + linker_nt + av

            df[out_col] = df.apply(_build_row, axis=1)

    df["linker"] = linker_aa

    return df

