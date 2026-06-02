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
Full-length TCR sequence assembly for TCRsift.

Builds complete TCR sequences including leader peptides and constant regions.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from importlib.resources import files as _pkg_files
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_directory_exists,
)

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
DEFAULT_LEADERS = {
    "CD8A": {
        "aa": "MALPVTALLLPLALLLHAARP",
        "dna": "ATGGCCCTGCCTGTGACAGCCCTGCTGCTGCCTCTGGCTCTGCTGCTGCATGCCGCTAGACCC",
        "source": "Human CD8A signal peptide (UniProt P01732)",
        "species": "human",
    },
    "CD28": {
        "aa": "MLRLLLALNLFPSIQVTG",
        "dna": "ATGCTCCGCCTGCTGCTGGCCCTGAACCTGTTCCCCAGCATCCAGGTGACCGGC",
        "source": "Human CD28 signal peptide (UniProt P10747)",
        "species": "human",
    },
    "IgK": {
        "aa": "METDTLLLWVLLLWVPGSTG",
        "dna": "ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGT",
        "source": "Murine IgGκ light chain signal peptide",
        "species": "mouse",
        "note": "Widely used for high secretion efficiency in mammalian expression",
    },
    "TRAC": {
        "aa": "MAGTWLLLLLALGCPALPTG",
        "dna": "ATGGCTGGCACCTGGCTGCTGCTGCTGCTGGCCCTGGGATGCCCAGCACTGCCCACAGGC",
        "source": "Human TRAC native signal peptide",
        "species": "human",
    },
    "TRBC": {
        "aa": "MGTSLLCWMALCLLGADHADG",
        "dna": "ATGGGCACCAGCCTGCTGTGCTGGATGGCCCTGTGCCTGCTGGGAGCAGACCACGCCGATGGC",
        "source": "Human TRBC native signal peptide",
        "species": "human",
    },
}

# Standard constant region endings for QC
CONSTANT_REGION_ENDINGS = {
    "TRAC": "LLMTLRLWSS",
    "TRBC1": "VKRKDF",
    "TRBC2": "VKRKDSRG",
}

# Minimum plausible length of the mature constant region per chain.
# TRAC mature ≈ 138 aa; TRBC1/TRBC2 mature ≈ 178 aa. The #66 bug produced
# 2–11 aa constants. Anything below these floors is a truncation.
CONSTANT_AA_FLOOR = {"alpha": 100, "beta": 150}

# In-cis J→C pairing for β chain: TRBJ1-* segments rearrange with TRBC1,
# TRBJ2-* with TRBC2. A row whose J-gene parity contradicts the called
# C-gene is either an annotation conflict or an assembly bug.
BETA_JC_PARITY = {"TRBJ1": "TRBC1", "TRBJ2": "TRBC2"}

# Allowed residue alphabet for assembled protein chains: the 20 standard
# AAs plus X (unknown — sometimes emitted by translation) and * (stop).
# Anything outside this set in a ``full_*_aa`` is corruption.
_VALID_AA_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYX*")


# Canonical human TCR constant-region amino acid sequences.
#
# These are the mature constant regions (post signal-cleavage), spanning
# the connecting peptide, transmembrane helix, and cytoplasmic tail —
# everything that goes downstream of the J-gene boundary in a cloned
# TCR construct.
#
# Sequences live in the packaged FASTA ``tcrsift/data/canonical_constants.fasta``
# (sourced from pyensembl GRCh38 release 110 ``transcript.protein_sequence``
# for transcripts ENST00000611116 / ENST00000633705 / ENST00000466254, then
# cross-checked against UniProt P01848 / P01850 / A0A5B9). The FASTA is the
# single source of truth at runtime; CI re-runs the pyensembl fetch and asserts
# the FASTA hasn't drifted (``test_canonical_matches_pyensembl``).
#
# History: an earlier hardcoded-string implementation (PR #69, fixing #66's
# pyensembl frame bug) drifted from the cited UniProt entries — TRAC had a
# C↔T error at the conserved Thr46 position, TRBC1 had E↔V at position 135
# of mature, TRBC2 had 5 errors including the JOVI.1-distinguishing K/N swap
# at positions 3-4. #100 traced the drift to a wrong UniProt accession in
# the source-of-truth comment (the cited A0A075B6Y0 is actually TRAJ49,
# not TRBC2). This loader replaces those hand-typed strings with FASTA
# regenerated from pyensembl so the same class of drift can't recur silently.
#
# As of 2.0 (#105), all three sequences are stored as BARE MATURE protein
# matching the UniProt chain annotations exactly — no synthetic residue is
# prepended. The J→C junction residue (universally E for β chains because
# the splice combines J's terminal nt with the C-exon's first 2 nt to spell
# GAG; J-dependent N/Y/D/H for α chains) is read per-clone from the
# CellRanger contig in :func:`_add_constant_regions`. Pre-2.0 versions
# baked a synthetic E into the stored TRBC1/TRBC2; that worked but
# couldn't represent α's per-clone variation and meant the canonical
# constants didn't match UniProt directly.
#
# Each sequence ends with the canonical C-terminus stored in
# ``CONSTANT_REGION_ENDINGS`` — see ``validate_sequences``.

_CANONICAL_CONSTANTS_FASTA = "canonical_constants.fasta"
_CANONICAL_ALLELES_FASTA = "canonical_constants_alleles.fasta"


def _load_canonical_constants_fasta() -> dict[str, str]:
    """Parse the packaged canonical-constants FASTA into a {gene: AA} dict.

    The FASTA is shipped in ``tcrsift/data/`` and is the runtime source of
    truth for the human canonical TCR constant-region AA sequences. See
    the module-level note above for the provenance chain.
    """
    fasta_text = (
        _pkg_files("tcrsift.refseqs").joinpath(_CANONICAL_CONSTANTS_FASTA).read_text()
    )
    out: dict[str, str] = {}
    name: str | None = None
    seq_parts: list[str] = []
    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(seq_parts)
            # Header format: >NAME|...|... — first token after ">" is the gene name
            name = line[1:].split("|", 1)[0].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if name is not None:
        out[name] = "".join(seq_parts)
    expected = {"TRAC", "TRBC1", "TRBC2"}
    missing = expected - out.keys()
    if missing:
        raise RuntimeError(
            f"{_CANONICAL_CONSTANTS_FASTA} is missing canonical entries: {sorted(missing)}"
        )
    return out


_canonical = _load_canonical_constants_fasta()
HUMAN_TRAC_AA: str = _canonical["TRAC"]
HUMAN_TRBC1_AA: str = _canonical["TRBC1"]
HUMAN_TRBC2_AA: str = _canonical["TRBC2"]
del _canonical


def _load_canonical_alleles_fasta() -> dict[str, dict[str, str]]:
    """Parse the packaged multi-allele FASTA into a {gene: {allele: AA}}
    nested dict. Used by the per-clone allele picker (#113).

    Header convention: ``>{GENE}*{NN}|provenance|notes`` — e.g.
    ``>TRBC2*01|NCBI:AAA60662|note:major-allele``. The loader keys on
    everything before the first ``|``, splits ``GENE`` from ``allele``
    around ``*``, and ignores allele *numbers* — only the protein
    sequence matters. Multiple IMGT allele numbers can collapse to the
    same protein form (TRBC2 *01/*02/*04 all encode the same mature
    protein); we ship them as a single entry under the lowest-numbered
    allele label.
    """
    fasta_text = (
        _pkg_files("tcrsift.refseqs").joinpath(_CANONICAL_ALLELES_FASTA).read_text()
    )
    out: dict[str, dict[str, str]] = {}
    name: str | None = None
    seq_parts: list[str] = []

    def _commit():
        if name is None:
            return
        gene, _, allele = name.partition("*")
        if not allele:
            allele = "01"
        out.setdefault(gene, {})[allele] = "".join(seq_parts)

    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            _commit()
            name = line[1:].split("|", 1)[0].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    _commit()
    expected = {"TRAC", "TRBC1", "TRBC2"}
    missing = expected - out.keys()
    if missing:
        raise RuntimeError(
            f"{_CANONICAL_ALLELES_FASTA} is missing canonical gene entries: "
            f"{sorted(missing)}"
        )
    return out


# Per-gene allele dict: ``HUMAN_CONSTANT_ALLELES[gene][allele] = AA``.
# Populated at import. Users can extend in-place to register their own
# alleles before calling :func:`assemble_full_sequences`.
HUMAN_CONSTANT_ALLELES: dict[str, dict[str, str]] = _load_canonical_alleles_fasta()


# Allele-called-reason enum (#118). Surfaces the *cause* of every
# allele call (or no-call) as a single machine-readable string so
# downstream consumers don't have to parse the human-readable
# ``{chain}_constant_source`` text or the ``qc_warnings`` list. Each
# clone receives exactly one value in ``{chain}_allele_called_reason``.
ALLELE_REASON_AUTO_DETECTED = "auto_detected"
ALLELE_REASON_OVERRIDDEN = "overridden"
ALLELE_REASON_INVALID_OVERRIDE = "invalid_override"
ALLELE_REASON_DIVERGENT_CONTIG = "divergent_contig"
ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION = "divergent_at_polymorphic_position"
ALLELE_REASON_SPARSE_CONTIG = "sparse_contig"
ALLELE_REASON_NO_CONTIG = "no_contig"
ALLELE_REASON_NO_ALLELE_POOL = "no_allele_pool"
ALLELE_REASON_FRAME_ERROR = "frame_error"

# All valid values, for validation + cohort audit grouping.
ALLELE_REASON_VALUES: tuple[str, ...] = (
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_INVALID_OVERRIDE,
    ALLELE_REASON_DIVERGENT_CONTIG,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_SPARSE_CONTIG,
    ALLELE_REASON_NO_CONTIG,
    ALLELE_REASON_NO_ALLELE_POOL,
    ALLELE_REASON_FRAME_ERROR,
)


def _polymorphic_positions(allele_pool: dict[str, str]) -> set[int]:
    """Return the set of mature-residue positions (0-indexed) where
    alleles in ``allele_pool`` disagree (#120).

    These are the residues that DISTINGUISH alleles from one another —
    e.g., mature TRBC2 position 8 (0-indexed) where *01 has E and *03
    has K. A contig that disagrees with the picked allele at one of
    these positions almost certainly indicates a novel polymorphism
    (or sequencing noise); the picker should refuse to commit rather
    than silently flip to a sibling on overall score.

    Positions past the shortest allele are excluded — we can't
    distinguish "polymorphic" from "absent" there.
    """
    if len(allele_pool) < 2:
        return set()
    alleles = list(allele_pool.values())
    n_compare = min(len(a) for a in alleles)
    polymorphic: set[int] = set()
    for i in range(n_compare):
        residues = {a[i] for a in alleles}
        if len(residues) > 1:
            polymorphic.add(i)
    return polymorphic


# Number of mature C-region residues scanned for allele divergences.
# Past this depth the contig is rarely covered and the residues are
# dominated by invariant regions (#120 ask 2). Shared by the stored
# divergence computation (:func:`_divergence_positions`), the recompute
# loop, and the observed-position denominator so all three scan the same
# window — otherwise a divergence could land outside the counted range.
_MAX_DIVERGENCE_SCAN_POSITIONS = 15


def _divergence_positions(
    called_allele_aa: str,
    observed_aa: str,
    *,
    max_positions: int = _MAX_DIVERGENCE_SCAN_POSITIONS,
) -> list[tuple[int, str, str]]:
    """Return a list of ``(position_1_indexed, expected_aa, observed_aa)``
    tuples where the observed contig translation disagrees with the
    called allele.

    Caps at the first ``max_positions`` mature residues — past that
    the contig is rarely covered and the residues are dominated by
    invariant regions (#120 ask 2).
    """
    if not called_allele_aa or not observed_aa:
        return []
    n = min(len(called_allele_aa), len(observed_aa), max_positions)
    return [
        (i + 1, called_allele_aa[i], observed_aa[i])
        for i in range(n)
        if called_allele_aa[i] != observed_aa[i]
    ]


def _format_divergence_positions(positions: list[tuple[int, str, str]]) -> str | None:
    """Serialize divergence positions for the output column.

    Format: ``"3:N->K;7:V->I"``. Empty list returns ``None`` so callers
    can use ``df[col].notna()`` to filter divergent clones.
    """
    if not positions:
        return None
    return ";".join(f"{p}:{e}->{o}" for p, e, o in positions)


def _score_allele_against_contig(
    contig_aa: str,
    allele_aa: str,
) -> tuple[int, int]:
    """Score one allele's mature protein against a per-clone contig
    translation past the J→C junction codon.

    Returns ``(n_agreeing, n_compared)`` where ``n_agreeing`` is the
    number of canonical positions where contig matches the allele
    and ``n_compared`` is the number of positions actually scored
    (the lesser of ``len(contig_aa)`` and ``len(allele_aa)``).

    Callers should pass ``contig_aa`` already trimmed past the J→C
    junction codon — the picker compares against the bare-mature
    canonical AA, which by convention doesn't include the junction
    residue (that's a per-clone concern handled at NT level in
    :func:`_add_constant_regions`, #105). For both chains the slice
    happens uniformly at NT level (``contig_nt_past_vdj[3:]``).
    """
    n = min(len(contig_aa), len(allele_aa))
    if n == 0:
        return 0, 0
    n_agree = sum(1 for a, b in zip(contig_aa[:n], allele_aa[:n]) if a == b)
    return n_agree, n


# Minimum number of residues the picker needs to score before it'll
# commit to an allele call. Set just past the deepest known
# distinguishing position so the comparison window always includes the
# residue that actually disambiguates the alleles we ship.
#
# Current floor: 10. TRBC2*01 vs *03 differ only at mature pos 9
# (idx 8). With a 5-codon contig the window is [0:5] which is too
# short to see idx 8 — the picker would tie and silently fall back
# to FASTA order. 10 codons covers idx 0..9 inclusive and gives a
# small safety margin.
_DEFAULT_MIN_PICKER_POSITIONS = 10


def _verify_first_codon_frame(
    contig_nt_past_junction: str,
    allele_aa: str,
) -> bool:
    """NT-level sanity check that the contig is codon-aligned to the
    allele (#113 follow-up).

    AA-level scoring is robust to synonymous donor codon variation but
    blind to a frame-off-by-N bug upstream (e.g. if the #91 vdj_nt
    trim regressed and ``contig_nt_past_vdj`` no longer starts on a
    clean codon boundary). This check translates the contig's first
    codon past the junction and asserts it matches the allele's first
    expected residue at the protein level — failure means our offset
    is wrong, regardless of how good the downstream AA score looks.
    """
    if len(contig_nt_past_junction) < 3 or not allele_aa:
        return False
    first_codon = contig_nt_past_junction[:3]
    translated = CODON_TABLE.get(first_codon, "X")
    return translated == allele_aa[0]


def _pick_best_allele(
    contig_aa: str,
    gene: str,
    candidate_alleles: dict[str, str] | None = None,
    *,
    min_compared: int = _DEFAULT_MIN_PICKER_POSITIONS,
) -> tuple[str | None, float, dict[str, float]]:
    """Pick the best-matching allele for ``gene`` (#113).

    Scans every candidate allele in ``candidate_alleles`` (defaulting to
    :data:`HUMAN_CONSTANT_ALLELES[gene]`), scores each against the
    contig's mature C-region translation, and returns ``(allele, score,
    all_scores)`` for the highest-scoring one.

    ``score`` is ``n_agreeing / n_compared`` (a fraction in [0, 1]).
    Ties are broken by the dict-iteration order, which is FASTA file
    order — the user's preferred allele can be listed first in the
    FASTA to bias tied calls. The audit-column format
    ``"GENE*ALLELE:score;..."`` (see
    :func:`_add_constant_regions`) is sorted by score descending and
    relies on the same tie-break rule.

    When coverage is too thin to commit to a call (``n_compared <
    min_compared`` for every allele), the function returns
    ``(None, 0.0, {})`` — the caller should fall back to the
    user-specified default. Pre-2.3 the first allele in the pool was
    returned as a "best-effort" signal; this caused silent
    fall-throughs when the comparison window didn't include any
    distinguishing position. ``None`` makes the no-decision explicit.

    Parameters
    ----------
    contig_aa
        Translation of the contig PAST the J→C junction codon. The
        picker compares position-wise against the bare-mature
        canonical (no junction residue prepended).
    gene
        ``"TRAC"`` / ``"TRBC1"`` / ``"TRBC2"``.
    candidate_alleles
        Override the FASTA-derived allele set (mostly for testing).
    min_compared
        Skip the call when the comparison window for every candidate
        allele is shorter than this. Default
        :data:`_DEFAULT_MIN_PICKER_POSITIONS` (10) — past the
        deepest known distinguishing position so a tie isn't
        silently broken by FASTA order.
    """
    pool = candidate_alleles if candidate_alleles is not None else HUMAN_CONSTANT_ALLELES.get(gene, {})
    if not pool:
        return None, 0.0, {}
    all_scores: dict[str, float] = {}
    best_allele: str | None = None
    best_n_agree = -1
    for allele, allele_aa in pool.items():
        n_agree, n_compared = _score_allele_against_contig(contig_aa, allele_aa)
        if n_compared >= min_compared:
            score = n_agree / n_compared
            all_scores[allele] = round(score, 3)
            # Prefer higher absolute n_agree (more positions matched);
            # tie-break by FASTA order (dict iteration is insertion
            # order in Python 3.7+).
            if n_agree > best_n_agree:
                best_n_agree = n_agree
                best_allele = allele
    if best_allele is None:
        return None, 0.0, {}
    return best_allele, all_scores[best_allele], all_scores

# Load-time invariant: each canonical AA must end with the C-terminus
# advertised in ``CONSTANT_REGION_ENDINGS``. Catching FASTA drift here is
# cheaper than discovering it in the assembly stage.
for _gene, _aa in (
    ("TRAC", HUMAN_TRAC_AA),
    ("TRBC1", HUMAN_TRBC1_AA),
    ("TRBC2", HUMAN_TRBC2_AA),
):
    _expected_tail = CONSTANT_REGION_ENDINGS[_gene]
    if not _aa.endswith(_expected_tail):
        raise RuntimeError(
            f"{_CANONICAL_CONSTANTS_FASTA}: {_gene} canonical AA tail {_aa[-15:]!r} "
            f"doesn't match CONSTANT_REGION_ENDINGS[{_gene!r}]={_expected_tail!r}"
        )
del _gene, _aa, _expected_tail

HUMAN_CONSTANT_REGIONS_AA: dict[str, str] = {
    "TRAC": HUMAN_TRAC_AA,
    "TRBC1": HUMAN_TRBC1_AA,
    "TRBC2": HUMAN_TRBC2_AA,
}


# Codon table for back-translating canonical AA constants to NT for
# downstream consumers that need DNA. Picks the most frequent human
# codon per residue (CCDS-weighted, GenScript table). This is
# representative — not authoritative for nucleotide-level work.
# Stop codon "*" maps to TAA.
HUMAN_PREFERRED_CODONS: dict[str, str] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TAA",
}


# Synonymous codon alternatives per AA, ordered by human CAI descending.
# Used by :func:`optimize_codons` for motif-aware back-translation: if the
# CAI-best codon would create a forbidden motif in the surrounding
# context, the picker walks down the list until it finds one that
# doesn't. Falls back to the lowest-CAI option when every alternative
# fails (rare; only possible when all synonyms share a problematic
# common prefix — practically never for the motifs we forbid).
#
# Sources: CCDS-weighted human codon usage (GenScript table) cross-
# checked against the Kazusa Codon Usage Database.
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
    # Common Type-II restriction sites (sense strand only).
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


def back_translate(aa: str) -> str:
    """Reverse-translate a polypeptide to DNA via :data:`HUMAN_PREFERRED_CODONS`.

    The result is one of many valid back-translations; this is the naive
    "single best codon per AA" version. For motif-aware back-translation
    that avoids creating restriction sites and homopolymer runs at codon
    boundaries, use :func:`optimize_codons` instead — which is what the
    assembly path calls.

    This function is kept for backward compatibility with external
    callers (and a few internal sites where motif checks would be
    counter-productive — e.g. when the caller is back-translating a
    single forced codon for a junction residue). Unknown residues fall
    back to NNN.
    """
    fallback = "NNN"
    out = [HUMAN_PREFERRED_CODONS.get(r, fallback) for r in aa]
    return "".join(out)


def optimize_codons(
    aa_sequence: str,
    *,
    prefix_nt: str = "",
    forbidden_motifs: tuple[str, ...] = DEFAULT_FORBIDDEN_MOTIFS,
    context_window: int = 8,
    gc_target_window: int = DEFAULT_GC_TARGET_WINDOW,
    gc_target_range: tuple[float, float] = DEFAULT_GC_TARGET_RANGE,
) -> str:
    """Back-translate ``aa_sequence`` to DNA with motif avoidance and
    GC-content balancing.

    For each residue, score synonymous codons (see
    :data:`HUMAN_CODON_ALTERNATIVES`) by three filters applied in
    order:

    1. **Motif avoidance (hard)** — reject codons that introduce a
       new forbidden substring in the trailing ``context_window`` nt.
       Default forbidden set covers homopolymer runs, common Type-II
       restriction sites, cryptic polyA signals, and canonical 5′
       splice donors (see :data:`DEFAULT_FORBIDDEN_MOTIFS`).
    2. **1-codon lookahead (soft preference)** — prefer codons that
       leave at least one valid motif-free choice for the next
       residue. Without lookahead, F+P+P in TRBC1 collapses into a
       6-mer C run.
    3. **GC-content windowing (soft preference)** — prefer codons
       that keep the trailing ``gc_target_window`` nt within
       ``gc_target_range``. Mammalian expression suffers at GC
       extremes (Kudla et al. 2006); the default ``(0.40, 0.65)``
       window is permissive but breaks ties.

    Final tiebreaker: CAI rank (first matching candidate wins).
    Deterministic: same input always produces the same output.

    Parameters
    ----------
    aa_sequence : str
        AA sequence. Standard 20 + ``*`` (stop). Other characters
        map to ``NNN``.
    prefix_nt : str
        Existing NT to keep in the buffer. Useful when optimizing a
        constant region appended to a VDJ — pass the trailing nt of
        VDJ so boundary motifs are caught.
    forbidden_motifs : tuple of str
        Substrings to avoid. Pass ``()`` to skip motif avoidance
        (equivalent to :func:`back_translate`).
    context_window : int
        Trailing nt retained for motif checks. Default 8 catches
        common 6-mer sites + 5-mer runs.
    gc_target_window : int
        Sliding-window length (nt) for GC balancing. Default 60.
    gc_target_range : tuple of float
        ``(low, high)`` GC fraction. Default ``(0.40, 0.65)``.
        Pass ``(0.0, 1.0)`` to disable GC balancing entirely.

    Notes
    -----
    For deeper expression optimization (mRNA 5' secondary-structure
    minimization, codon-pair bias, tAI — Tuller et al. 2010; Goodman
    et al. 2013), use a constraint solver like ``dnachisel``. The
    in-house optimizer here covers what synthesis vendors filter on
    by default (motifs, GC) plus the two highest-leverage expression
    pathologies (polyA + splice donors). For short C-region
    sequences (≤180 aa), this is empirically competitive with vendor
    tools on every restriction-site / polyA benchmark we've
    checked.

    Motifs already present in the prefix buffer are treated as
    pre-existing and not blamed on subsequent codons (which can't
    rewrite them anyway).
    """
    if gc_target_range[0] > gc_target_range[1]:
        raise ValueError(
            f"gc_target_range must be (low, high); got {gc_target_range}"
        )
    gc_balancing_active = gc_target_range != (0.0, 1.0)

    def _new_motifs(buffer_with_codon: str, baseline: set[str]) -> bool:
        return any(
            motif not in baseline and motif in buffer_with_codon
            for motif in forbidden_motifs
        )

    def _has_valid_next(next_aa: str, new_buffer: str) -> bool:
        next_alts = HUMAN_CODON_ALTERNATIVES.get(next_aa)
        if not next_alts:
            return True
        baseline = {m for m in forbidden_motifs if m in new_buffer}
        for next_codon in next_alts:
            if not _new_motifs(new_buffer + next_codon, baseline):
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
    buffer = prefix_nt[-context_window:] if prefix_nt else ""
    # GC tracking uses a longer trailing window than motif checks.
    gc_buffer = prefix_nt[-gc_target_window:] if prefix_nt else ""
    preexisting = {m for m in forbidden_motifs if m in buffer}
    for i, aa in enumerate(aa_sequence):
        alternatives = HUMAN_CODON_ALTERNATIVES.get(aa)
        if not alternatives:
            codon = HUMAN_PREFERRED_CODONS.get(aa, "NNN")
            nt_parts.append(codon)
            buffer = (buffer + codon)[-context_window:]
            gc_buffer = (gc_buffer + codon)[-gc_target_window:]
            preexisting = {m for m in forbidden_motifs if m in buffer}
            continue
        next_aa = aa_sequence[i + 1] if i + 1 < len(aa_sequence) else None
        # Filter 1 (hard): motif-free
        motif_free = [
            c for c in alternatives
            if not _new_motifs(buffer + c, preexisting)
        ]
        if not motif_free:
            # Every alternative introduces a motif; pick the lowest-CAI
            # option to minimize cascade pressure on the next codon.
            chosen = alternatives[-1]
        else:
            # Filter 2 (soft): lookahead-valid. Refine motif_free; if
            # the refined set is empty, fall back to motif_free.
            lookahead_valid = [
                c for c in motif_free
                if next_aa is None
                or _has_valid_next(next_aa, (buffer + c)[-context_window:])
            ] or motif_free
            # Filter 3 (soft): GC-in-range. Refine; fall back to
            # lookahead_valid if every option is out of range.
            gc_balanced = [
                c for c in lookahead_valid if _gc_in_range(c, gc_buffer)
            ] or lookahead_valid
            # CAI tiebreak: first remaining candidate is highest-CAI.
            chosen = gc_balanced[0]
        nt_parts.append(chosen)
        buffer = (buffer + chosen)[-context_window:]
        gc_buffer = (gc_buffer + chosen)[-gc_target_window:]
        preexisting = {m for m in forbidden_motifs if m in buffer}
    return "".join(nt_parts)


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


def pick_canonical_constant(
    chain: str, c_gene: str | None = None, j_gene: str | None = None
) -> tuple[str, str]:
    """Choose the canonical constant-region AA for a chain.

    - Alpha → always TRAC.
    - Beta → J family decides when known (TRBJ1-* with TRBC1, TRBJ2-*
      with TRBC2), because the in-cis J→C pairing is locus-determined
      and reliable. ``c_gene`` is only consulted when ``j_gene`` is
      missing or its family can't be parsed; CellRanger's TRBC1/TRBC2
      call is unreliable (the two C genes share ~95% NT identity) and
      gets silently overridden when it contradicts J (#90). Aggregate
      override counts are surfaced at the assembly orchestrator level
      so a per-row override doesn't produce a log flood.

    Returns ``(gene_name, canonical_aa)``.
    """
    if chain == "alpha":
        return "TRAC", HUMAN_TRAC_AA

    j_family = _beta_j_family(j_gene)
    expected_c = BETA_JC_PARITY.get(j_family) if j_family else None
    if expected_c == "TRBC1":
        return "TRBC1", HUMAN_TRBC1_AA
    if expected_c == "TRBC2":
        return "TRBC2", HUMAN_TRBC2_AA

    # J family unknown — fall back to the (less reliable) c_gene call.
    c_base = _beta_c_base(c_gene)
    if c_base == "TRBC2":
        return "TRBC2", HUMAN_TRBC2_AA
    if c_base == "TRBC1":
        return "TRBC1", HUMAN_TRBC1_AA
    # Default to TRBC1 when nothing else is known.
    return "TRBC1", HUMAN_TRBC1_AA


def _beta_j_family(j_gene: str | None) -> str:
    """Parse the J-family prefix (``"TRBJ1"`` / ``"TRBJ2"``) from a
    raw J-gene call. Tolerates lowercase, whitespace, and allele
    suffixes (``"TRBJ1-1*02"``). Returns an empty string when the
    input doesn't look like a TRBJ call."""
    if not isinstance(j_gene, str):
        return ""
    token = j_gene.strip().upper().split("*", 1)[0]
    return token.split("-", 1)[0] if "-" in token else token


def _trbc_override_mask(df: pd.DataFrame) -> pd.Series:
    """Boolean Series identifying rows where ``beta_c_gene_canonical``
    differs from raw ``beta_c_gene`` after gene-family normalisation —
    the rows where :func:`pick_canonical_constant` overrode
    CellRanger's TRBC call.

    Returns an all-False Series (aligned to ``df.index``) when either
    column is missing. Used both by :func:`_count_trbc_overrides` and
    to emit the per-row ``beta_c_gene_overridden`` audit column.
    """
    if "beta_c_gene" not in df.columns or "beta_c_gene_canonical" not in df.columns:
        return pd.Series(False, index=df.index)
    raw_base = df["beta_c_gene"].apply(_beta_c_base)
    canon_base = df["beta_c_gene_canonical"].apply(_beta_c_base)
    return (raw_base != "") & (canon_base != "") & (raw_base != canon_base)


def _count_trbc_overrides(df: pd.DataFrame) -> int:
    """Count rows where :func:`pick_canonical_constant` overrode
    CellRanger's TRBC call. Thin wrapper over
    :func:`_trbc_override_mask` for callers that just want a tally."""
    return int(_trbc_override_mask(df).sum())


def _beta_c_base(c_gene: str | None) -> str:
    """Parse the C-gene base (``"TRBC1"`` / ``"TRBC2"``) from a raw
    C-gene call. Tolerates lowercase, whitespace, and allele suffixes
    (``"TRBC2*01"``). Returns an empty string when the input doesn't
    resolve to a known TRBC."""
    if not isinstance(c_gene, str):
        return ""
    token = c_gene.strip().upper().split("*", 1)[0]
    return token if token in {"TRBC1", "TRBC2"} else ""


def _resolve_c_gene(row, chain: str) -> str:
    """Coalesce ``{chain}_c_gene_canonical`` → ``{chain}_c_gene`` for a
    single row, treating NaN as missing.

    The naive ``row.get(canonical) or row.get(raw) or ""`` is wrong
    because ``float('nan')`` is truthy: when the canonical column is
    present but the cell is NaN, the ``or`` short-circuits to NaN and
    the raw column is never consulted. That's silent: validation
    downgrades structural failures to informational notes, and
    :func:`fix_jc_parity` silently skips the correction.
    """
    canonical = row.get(f"{chain}_c_gene_canonical")
    if isinstance(canonical, str) and canonical:
        return canonical
    raw = row.get(f"{chain}_c_gene")
    if isinstance(raw, str) and raw:
        return raw
    return ""


def verify_canonical_constant_start(
    observed: str, canonical_aa: str, min_match: int = 8
) -> bool:
    """Sanity-check the observed C-region start against the canonical.

    Returns True if at least ``min_match`` of the first ``len(observed)``
    residues match position-wise. Used to verify a CellRanger contig's
    post-J residues against the picked canonical before splicing the
    canonical full sequence into the assembled construct.
    """
    n = min(len(observed), len(canonical_aa))
    if n == 0:
        return False
    matches = sum(1 for a, b in zip(observed[:n], canonical_aa[:n]) if a == b)
    return matches >= min_match


def _describe_leader(param_value: str | None, resolved: str | dict | None) -> str:
    """Generate a description of a leader configuration for logging."""
    if param_value is None:
        return "None (no leader)"
    if resolved == "from_contig":
        return "from_contig (extract from FASTA)"
    if isinstance(resolved, dict):
        return f"{param_value.upper()} ({resolved['source']})"
    return str(param_value)


def translate_dna(dna_seq: str) -> tuple[str, str]:
    """
    Translate DNA sequence to amino acids.

    Returns
    -------
    tuple
        (amino_acid_sequence, ragged_3p_nucleotides)
    """
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


def parse_fasta(path: str | Path) -> dict[str, str]:
    """
    Parse a FASTA file.

    Returns
    -------
    dict
        Mapping from sequence ID to sequence
    """
    path = Path(path)
    results = {}
    curr_id = None
    lines = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id and lines:
                    results[curr_id] = "".join(lines)
                    lines = []
                curr_id = line[1:].split()[0]  # Take first word after >
            else:
                lines.append(line)

        # Don't forget last entry
        if curr_id and lines:
            results[curr_id] = "".join(lines)

    return results


SAMPLE_NAME_FROM_CHOICES = ("parent", "grandparent", "sheet")


def load_contigs(
    contig_dir: str | Path | None = None,
    *,
    sample_name_from: str = "parent",
    cellranger_dir: str | Path | None = None,
    sample_sheet: pd.DataFrame | None = None,
) -> dict[str, dict[str, str]]:
    """
    Load contig sequences from a CellRanger-style output tree (#124).

    Parameters
    ----------
    contig_dir : str or Path, optional
        Root directory to scan for ``*contig*.fasta``. Mutually exclusive
        with ``cellranger_dir``.
    sample_name_from : {'parent', 'grandparent', 'sheet'}, default 'parent'
        How to derive the sample name from each discovered FASTA's path:

        - ``'parent'`` (default, backward-compat) — sample = FASTA's
          immediate parent directory name. Matches the symlinked layout
          ``contig_dir/{sample}/filtered_contig.fasta``.
        - ``'grandparent'`` — sample = FASTA's grandparent directory name.
          Matches CellRanger ``multi`` output:
          ``per_sample_outs/{sample}/vdj_t/filtered_contig.fasta``.
        - ``'sheet'`` — match each FASTA against ``sample_sheet['vdj_dir']``
          and take the corresponding ``sample`` column. Most explicit;
          useful for non-standard layouts.
    cellranger_dir : str or Path, optional
        Shorthand for ``contig_dir=cellranger_dir, sample_name_from='grandparent'``.
        Pass a raw CellRanger ``per_sample_outs/`` directory and the sample
        names will resolve correctly without a symlink tree. Mutually
        exclusive with ``contig_dir``.
    sample_sheet : pandas.DataFrame, optional
        Required when ``sample_name_from='sheet'``. Must have columns
        ``sample`` and ``vdj_dir``; each FASTA's path must be inside one of
        the listed ``vdj_dir`` paths.

    Returns
    -------
    dict
        Nested dict: ``sample -> contig_id -> sequence``.
    """
    if cellranger_dir is not None:
        if contig_dir is not None:
            raise ValueError(
                "load_contigs: pass exactly one of contig_dir or cellranger_dir, not both."
            )
        contig_dir = cellranger_dir
        sample_name_from = "grandparent"
    if contig_dir is None:
        raise ValueError(
            "load_contigs: must pass either contig_dir or cellranger_dir."
        )
    if sample_name_from not in SAMPLE_NAME_FROM_CHOICES:
        raise ValueError(
            f"load_contigs: sample_name_from={sample_name_from!r} not in "
            f"{SAMPLE_NAME_FROM_CHOICES}"
        )
    if sample_name_from == "sheet" and sample_sheet is None:
        raise ValueError(
            "load_contigs: sample_name_from='sheet' requires sample_sheet=..."
        )

    contig_dir = Path(contig_dir)
    sample_contigs: dict[str, dict[str, str]] = {}

    sheet_dir_to_sample: dict[Path, str] = {}
    if sample_name_from == "sheet":
        missing = {"sample", "vdj_dir"} - set(sample_sheet.columns)
        if missing:
            raise ValueError(
                f"load_contigs: sample_sheet is missing required columns {sorted(missing)}"
            )
        sheet_dir_to_sample = {
            Path(row.vdj_dir).resolve(): row.sample
            for row in sample_sheet.itertuples(index=False)
        }

    def _sample_for(fasta_path: Path) -> str | None:
        if sample_name_from == "parent":
            return fasta_path.parent.name
        if sample_name_from == "grandparent":
            return fasta_path.parent.parent.name
        # sheet
        resolved = fasta_path.resolve()
        for sheet_path, sample in sheet_dir_to_sample.items():
            try:
                resolved.relative_to(sheet_path)
            except ValueError:
                continue
            return sample
        return None

    for fasta_path in contig_dir.rglob("*contig*.fasta"):
        sample_name = _sample_for(fasta_path)
        if sample_name is None:
            logger.warning(
                "load_contigs: %s did not match any sample_sheet vdj_dir; skipping",
                fasta_path,
            )
            continue
        sample_contigs.setdefault(sample_name, {}).update(parse_fasta(fasta_path))

    # Flat-layout fallback: bare ``*.fasta`` files directly in contig_dir.
    # Only meaningful for the 'parent' default — the grandparent/sheet
    # modes assume the CellRanger-style nested layout where this wouldn't
    # match anything.
    if sample_name_from == "parent":
        for fasta_path in contig_dir.glob("*.fasta"):
            sample_name = fasta_path.stem.split("_")[0]
            sample_contigs.setdefault(sample_name, {}).update(parse_fasta(fasta_path))

    logger.info(f"Loaded contigs from {len(sample_contigs)} samples")
    return sample_contigs


def get_constant_region_sequences(
    *,
    stop_codons: tuple[str, ...] = ("TAA", "TGA"),
) -> dict[str, str]:
    """
    Return human TCR constant-region CDS (DNA) via codon-aware back-
    translation.

    Sources NT from :data:`HUMAN_CONSTANT_REGIONS_AA` and
    :func:`optimize_codons` (motif-aware — avoids 5+ mononucleotide
    runs and common Type-II restriction sites; see #116). The earlier
    pyensembl-backed implementation read the full mRNA at frame
    offset 2 and silently truncated TRAC / TRBC1 / TRBC2 to 2–11
    residues for every assembled clonotype (#66). Hardcoding
    eliminates the frame bug, drops the pyensembl dependency, and
    matches the canonical sequences that downstream cloning
    constructs need (#67).

    Parameters
    ----------
    stop_codons : tuple of str, default ``("TAA", "TGA")``
        Stop codons to append to each CDS. Default is the two
        non-redundant stops (different release factors → reduces
        read-through in synthesized constructs).

    Returns
    -------
    dict
        Gene name → CDS (DNA, post-junction-codon-aligned … stop).
        Each CDS is the codon-optimized translation of the canonical
        AA plus the requested stop codons.
    """
    stops_nt = stop_codons_nt(stop_codons)
    out: dict[str, str] = {}
    for name, aa in HUMAN_CONSTANT_REGIONS_AA.items():
        out[name] = optimize_codons(aa) + stops_nt
    return out


def assemble_full_sequences(
    clonotypes: pd.DataFrame,
    contigs_dir: str | Path | None = None,
    alpha_leader: str | None = "CD28",
    beta_leader: str | None = "CD8A",
    include_constant: bool = True,
    constant_source: str = "ensembl",
    linker: str = "T2A",
    trac_allele: str = "auto",
    trbc1_allele: str = "auto",
    trbc2_allele: str = "auto",
    stop_codons: tuple[str, ...] = ("TAA", "TGA"),
    verbose: bool = True,
    show_progress: bool = True,
    sample_name_from: str = "parent",
    cellranger_dir: str | Path | None = None,
    sample_sheet: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Assemble full-length TCR sequences.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame with VDJ sequences (from fwr1/cdr1/fwr2/cdr2/fwr3/cdr3/fwr4)
    contigs_dir : str or Path, optional
        Directory with CellRanger contig FASTA files. Required if alpha_leader or
        beta_leader is set to "from_contig".
    alpha_leader : str or None
        Leader sequence for alpha chain. Options:
        - None: No leader sequence
        - "from_contig": Extract native leader from contig FASTA (requires contigs_dir)
        - Key from DEFAULT_LEADERS: "CD8A", "CD28", "IgK", "TRAC", "TRBC"
        Default is "CD28" to provide distinct sequences from beta chain.
    beta_leader : str or None
        Leader sequence for beta chain. Same options as alpha_leader.
        Default is "CD8A" to provide distinct sequences from alpha chain.
    include_constant : bool
        Include constant region sequences.
    constant_source : str
        Source for constant regions:

        - ``"canonical"`` (default): splice in the hardcoded canonical
          human TRAC / TRBC1 / TRBC2 (:data:`HUMAN_CONSTANT_REGIONS_AA`).
          When CellRanger contigs are available, the canonical's first
          ~15 residues are verified against the observed contig start
          and a ``qc_warnings`` entry is emitted on mismatch.
        - ``"ensembl"``: back-compat alias for ``"canonical"``. The
          earlier pyensembl-backed path was removed in #66 — it read
          the full mRNA at the wrong frame offset and silently
          truncated constants to 2–11 residues.
        - ``"from-data"``: read ``{chain}_constant_aa`` /
          ``{chain}_constant_nt`` directly from the input frame.
    linker : str
        Linker sequence for single-chain constructs: "T2A", "P2A", "E2A", "F2A"
    trac_allele, trbc1_allele, trbc2_allele : str
        Per-gene allele selection (#113). Default ``"auto"`` runs the
        allele picker over packaged alleles in
        :data:`HUMAN_CONSTANT_ALLELES` scoring each against the
        donor's contig translation; pass an explicit allele label
        (e.g. ``"01"`` or ``"03"``) to force a specific canonical.
        Invalid labels emit a QC warning and fall back to the
        default. ``constant_source="from-data"`` ignores these knobs
        (the constant comes from the input frame, not the FASTA);
        a warning is logged in that combination.
    stop_codons : tuple of str, default ``("TAA", "TGA")``
        Stop codons appended to the codon-optimized constant CDS.
        Default uses two non-redundant stops (recognized by
        different release factors → reduced read-through in
        synthesized constructs). Pass ``("TAA",)`` for single-stop
        (pre-2.4 behavior) or ``()`` to omit stops entirely. Each
        entry must be one of ``"TAA"``/``"TAG"``/``"TGA"``;
        invalid entries raise ``ValueError``. Stops are only
        appended to the ``_optimized`` columns and the final
        ``_constant_nt`` (alias) — not to ``_contig`` columns,
        which preserve the donor's CellRanger contig bytes verbatim
        (and are typically truncated past the contig coverage edge).
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with full sequences added. Allele audit columns
        emitted per clone (#113):

        - ``{chain}_allele_called`` — e.g. ``"TRBC2*01"`` or ``None``
          when no contig + non-explicit override.
        - ``{chain}_allele_score`` — float ∈ [0, 1] (1.0 = perfect AA
          agreement; 1.0 also when the user forced an allele).
        - ``{chain}_allele_alternatives`` — semicolon-joined string
          of runner-up alleles + scores, sorted by score desc with
          FASTA order as tie-break. Example: ``"TRBC2*01:1.000;TRBC2*03:0.933"``.
          ``None`` when the picker no-decided or only one allele was
          packaged.

        Constant NT triad emitted per chain (#116):

        - ``{chain}_constant_nt_contig`` — pure CellRanger contig
          NT past the J→C junction (no canonical splicing, no
          codon optimization). Truncated where contig coverage
          ends. ``None`` when no contig is available for the clone.
        - ``{chain}_constant_nt_optimized`` — codon-optimized
          canonical CDS for the picked allele (motif-aware via
          :func:`optimize_codons`; avoids common restriction sites
          and homopolymer runs) plus the requested stop codons.
        - ``{chain}_constant_nt`` — the assembly-aware blend: uses
          donor-real bytes from the contig where they agree with
          the canonical AA, falls back to ``_optimized`` for the
          rest. This is the column most callers want; kept as the
          default for back-compat.

        The same triad is also exposed for ``full_{chain}_nt``
        (leader + VDJ + constant + stop).

    Examples
    --------
    >>> # Default: CD28 on alpha, CD8A on beta (distinct leaders)
    >>> assembled = assemble_full_sequences(clonotypes)

    >>> # No leader sequences
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader=None)

    >>> # Leader only on beta chain (first in 2A construct)
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader="CD8A")

    >>> # Extract native leaders from contig FASTAs
    >>> assembled = assemble_full_sequences(
    ...     clonotypes,
    ...     contigs_dir="/path/to/contigs",
    ...     alpha_leader="from_contig",
    ...     beta_leader="from_contig",
    ... )
    """
    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_assembly=True)

    # Validate stop_codons early so users get the error before any
    # heavy lifting. Empty tuple is allowed (caller wants no stops).
    stops_nt = stop_codons_nt(stop_codons)

    valid_constant_sources = ["canonical", "ensembl", "from-data"]
    if constant_source not in valid_constant_sources:
        raise TCRsiftValidationError(
            f"Invalid constant_source: '{constant_source}'",
            hint=f"Valid options are: {valid_constant_sources}",
        )

    # #113: warn if the user combined ``constant_source="from-data"``
    # with an explicit allele override — the from-data path reads the
    # constant directly from the input frame and doesn't consult the
    # packaged FASTA, so allele kwargs are silently ignored. Surfacing
    # this loudly prevents a "I set --trbc2-allele 03 but my output is
    # still *01" confusion.
    if constant_source == "from-data" and any(
        a != "auto" for a in (trac_allele, trbc1_allele, trbc2_allele)
    ):
        logger.warning(
            "constant_source='from-data' ignores trac_allele / "
            "trbc1_allele / trbc2_allele — the constant comes from "
            "the input frame, not the packaged canonical FASTA. "
            "Reset the allele kwargs to 'auto' or switch to "
            "constant_source='canonical' to use them."
        )

    # Validate and resolve leader options for each chain
    leader_config = {}
    for chain, leader_param in [("alpha", alpha_leader), ("beta", beta_leader)]:
        if leader_param is None:
            leader_config[chain] = None
        elif leader_param.lower() == "from_contig":
            if not contigs_dir:
                raise TCRsiftValidationError(
                    f"{chain}_leader='from_contig' requires contigs_dir to be specified",
                    hint="Provide contigs_dir with CellRanger FASTA files, or use a default leader like 'CD8A'",
                )
            leader_config[chain] = "from_contig"
        elif leader_param.upper() in DEFAULT_LEADERS:
            leader_config[chain] = DEFAULT_LEADERS[leader_param.upper()]
        else:
            raise TCRsiftValidationError(
                f"Unknown {chain}_leader: '{leader_param}'",
                hint=f"Valid options are: None, 'from_contig', or one of {list(DEFAULT_LEADERS.keys())}",
            )

    if verbose:
        alpha_desc = _describe_leader(alpha_leader, leader_config["alpha"])
        beta_desc = _describe_leader(beta_leader, leader_config["beta"])
        logger.info(f"Assembling full sequences for {len(clonotypes):,} clonotypes")
        logger.info(f"  Alpha leader: {alpha_desc}")
        logger.info(f"  Beta leader: {beta_desc}")
        logger.info(f"  Constant regions: {include_constant} (source: {constant_source})")
        logger.info(f"  Linker: {linker}")

    df = clonotypes.copy()

    # Load constant regions if needed. "canonical" / "ensembl" both
    # hit the same hardcoded-AA path now (#66, #67); "ensembl" is
    # retained as a back-compat alias and emits an info note.
    constant_seqs = {}
    if include_constant and constant_source in ("canonical", "ensembl"):
        if constant_source == "ensembl" and verbose:
            logger.info(
                "  constant_source='ensembl' is now an alias for 'canonical' "
                "(#66) — splicing canonical TRAC / TRBC1 / TRBC2 from "
                "HUMAN_CONSTANT_REGIONS_AA."
            )
        constant_seqs = get_constant_region_sequences(stop_codons=stop_codons)
        if verbose:
            logger.info(
                f"    Loaded {len(constant_seqs)} canonical constant region sequences"
            )
            if stop_codons:
                logger.info(
                    f"    Appending stop codons {list(stop_codons)} to "
                    "codon-optimized constants."
                )
            else:
                logger.info(
                    "    No stop codons configured (stop_codons=()); "
                    "constants will not be terminated."
                )

    # Warn if from-data constants requested but not present
    if include_constant and constant_source == "from-data":
        constant_cols = [
            "alpha_constant_aa",
            "alpha_constant_nt",
            "beta_constant_aa",
            "beta_constant_nt",
        ]
        if not any(col in df.columns for col in constant_cols):
            logger.warning(
                "  constant_source='from-data' but no constant region columns found in input. "
                "Constants will be omitted."
            )

    # Load contigs whenever a contigs_dir is given. Pre-1.3 we only
    # loaded when a leader was "from_contig"; that gate is too narrow
    # now that the C-region NT blend (#103) also consumes contigs to
    # splice donor-real bytes into the J→C boundary. Loading is cheap
    # relative to assembly and avoids a silent "contigs given but
    # blend silently noop'd" failure mode.
    sample_contigs: dict = {}
    # cellranger_dir is shorthand for contigs_dir + sample_name_from='grandparent'
    # (#124). Resolved here so the rest of the function only sees contigs_dir.
    if cellranger_dir is not None:
        if contigs_dir is not None:
            raise ValueError(
                "assemble_full_sequences: pass exactly one of contigs_dir or "
                "cellranger_dir, not both."
            )
        contigs_dir = cellranger_dir
        sample_name_from = "grandparent"
    if contigs_dir:
        contigs_dir = validate_directory_exists(Path(contigs_dir), "contigs directory")
        if verbose:
            logger.info(
                f"  Loading contigs from {contigs_dir} (sample_name_from={sample_name_from!r})..."
            )
        sample_contigs = load_contigs(
            contigs_dir,
            sample_name_from=sample_name_from,
            sample_sheet=sample_sheet,
        )
        if verbose:
            total_contigs = sum(len(c) for c in sample_contigs.values())
            logger.info(f"    Loaded {total_contigs:,} contigs from {len(sample_contigs)} samples")

    # Process each clonotype
    if verbose:
        logger.info("  Assembling sequences...")

    assembly_results = []

    # Build per-gene allele override dict for the picker (#113).
    allele_overrides = {
        "TRAC": trac_allele,
        "TRBC1": trbc1_allele,
        "TRBC2": trbc2_allele,
    }

    # Create iterator with optional progress bar
    row_iter = df.iterrows()
    if show_progress:
        row_iter = tqdm(
            list(df.iterrows()),
            desc="Assembling sequences",
            unit="clone",
        )

    for idx, row in row_iter:
        result = _assemble_clone(
            row,
            sample_contigs,
            constant_seqs,
            leader_config,
            include_constant,
            constant_source,
            allele_overrides=allele_overrides,
            stops_nt=stops_nt,
        )
        assembly_results.append(result)

    # Add assembly columns to dataframe
    result_df = pd.DataFrame(assembly_results)
    for col in result_df.columns:
        df[col] = result_df[col].values

    # Add single-chain construct if requested
    if linker and "full_beta_aa" in df.columns and "full_alpha_aa" in df.columns:
        if verbose:
            logger.info(f"  Creating single-chain constructs with {linker} linker...")
        df = _add_single_chain(df, linker)

    # Summary
    if verbose:
        n_with_alpha = df["full_alpha_aa"].notna().sum() if "full_alpha_aa" in df.columns else 0
        n_with_beta = df["full_beta_aa"].notna().sum() if "full_beta_aa" in df.columns else 0
        n_single_chain = (
            df["single_chain_aa"].notna().sum() if "single_chain_aa" in df.columns else 0
        )
        logger.info("  Assembly complete:")
        logger.info(f"    With full alpha: {n_with_alpha:,}")
        logger.info(f"    With full beta: {n_with_beta:,}")
        logger.info(f"    Single-chain constructs: {n_single_chain:,}")

    # Per-row override audit column (#90). True where the canonical
    # picked by the J-family rule disagreed with CellRanger's raw
    # `beta_c_gene` call. Lets users filter / inspect post-hoc
    # without going through the log. Written even when no overrides
    # fired so the column is always present (defaults to False).
    if "beta_c_gene" in df.columns:
        df["beta_c_gene_overridden"] = _trbc_override_mask(df)

    # One aggregate line for J-family overrides of CellRanger TRBC
    # (#90). pick_canonical_constant is silent per-call to avoid a
    # log flood on cohorts where every other clone hits the case.
    # Gated by `verbose` so silent runs stay silent — the override
    # itself still happens, this is the audit trail. The per-row
    # `beta_c_gene_overridden` column above is the data-side audit
    # that survives the log.
    if verbose:
        n_overrides = (
            int(df["beta_c_gene_overridden"].sum())
            if "beta_c_gene_overridden" in df.columns
            else 0
        )
        if n_overrides:
            n_beta = df["beta_c_gene"].notna().sum() if "beta_c_gene" in df.columns else 0
            logger.warning(
                "  Overrode CellRanger TRBC call on %d / %d β clones to match "
                "J-family parity (TRBJ1→TRBC1, TRBJ2→TRBC2). CellRanger's "
                "TRBC1/2 discrimination is unreliable; the locus rule is "
                "authoritative.",
                n_overrides, n_beta,
            )

        # Cohort allele audit (#119, #120). Emits the per-chain
        # called/no-called breakdown + any novel-allele candidates
        # detected across the cohort. Silent on small cohorts where
        # the heuristic thresholds don't fire.
        if any(
            f"{c}_allele_called_reason" in df.columns for c in ("alpha", "beta")
        ):
            for line in allele_audit_report(df).splitlines():
                logger.info(line)

    return df


def _assemble_clone(
    row: pd.Series,
    sample_contigs: dict,
    constant_seqs: dict,
    leader_config: dict,
    include_constant: bool,
    constant_source: str,
    allele_overrides: dict[str, str] | None = None,
    stops_nt: str = "TAATGA",
) -> dict:
    """Assemble full sequence for a single clone."""
    result = {}

    # Try to get full sequence from VDJ columns if available.
    # Note (#91): the CellRanger-produced ``VDJ_*_nt`` typically has a
    # +1 trailing nucleotide past the J segment — a contig artifact
    # that does not encode any AA. If we leave it in, splicing
    # leader + VDJ + constant frame-shifts everything past the VDJ→C
    # boundary. Trim to ``3 * len(vdj_aa)`` here so the canonical
    # stored ``vdj_*_nt`` is always a clean ORF and every downstream
    # NT concatenation translates back to the assembled AA.
    for chain in ["alpha", "beta"]:
        vdj_col = f"VDJ_{chain}_aa"
        vdj_nt_col = f"VDJ_{chain}_nt"

        if vdj_col in row and pd.notna(row.get(vdj_col)):
            result[f"vdj_{chain}_aa"] = row[vdj_col]
        if vdj_nt_col in row and pd.notna(row.get(vdj_nt_col)):
            vdj_nt = row[vdj_nt_col]
            vdj_aa = result.get(f"vdj_{chain}_aa")
            if isinstance(vdj_aa, str) and vdj_aa:
                vdj_nt = vdj_nt[: 3 * len(vdj_aa)]
            result[f"vdj_{chain}_nt"] = vdj_nt

        # Get C gene and J gene for constant-region lookup. J gene is
        # needed by ``pick_canonical_constant`` to override
        # CellRanger's unreliable TRBC call via locus-determined
        # J-family parity (#90); without it the J-based override is
        # silently bypassed in the real assembly path.
        for gene in ("c_gene", "j_gene"):
            col = f"{chain}_{gene}"
            if col in row:
                result[col] = row[col]

    # Add leader sequences based on per-chain config
    for chain in ["alpha", "beta"]:
        chain_leader = leader_config.get(chain)
        if chain_leader is None:
            continue
        elif chain_leader == "from_contig":
            # Extract native leader from contigs
            _extract_leader_from_contigs_single(row, sample_contigs, result, chain)
        elif isinstance(chain_leader, dict):
            # Use specified default leader
            result[f"{chain}_leader_aa"] = chain_leader["aa"]
            result[f"{chain}_leader_nt"] = chain_leader["dna"]

    # Add constant regions
    if include_constant:
        if constant_source == "from-data":
            _add_constant_from_row(row, result)
        else:
            # ``canonical`` / ``ensembl`` (back-compat alias) flow:
            # pick the canonical TRAC / TRBC1 / TRBC2, splice the AA
            # in, and verify against the contig start when available.
            _add_constant_regions(
                result, constant_seqs,
                row=row, sample_contigs=sample_contigs,
                allele_overrides=allele_overrides,
                stops_nt=stops_nt,
            )

    # Determine which chains have leaders for building full sequences
    include_alpha_leader = leader_config.get("alpha") is not None
    include_beta_leader = leader_config.get("beta") is not None

    # Build full sequences
    _build_full_sequences(result, include_alpha_leader, include_beta_leader, include_constant)

    return result


def _extract_leader_from_contigs_single(
    row: pd.Series,
    sample_contigs: dict,
    result: dict,
    chain: str,
):
    """Extract leader peptide from contig sequences for a single chain."""
    samples = str(row.get("samples", "")).split(";")

    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return

    contig_ids = str(row[contig_col]).split(";")
    vdj_aa = result.get(f"vdj_{chain}_aa", "")

    leader_counter = Counter()
    leader_dna_counter = Counter()

    for sample in samples:
        if sample not in sample_contigs:
            continue

        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue

            contig_seq = sample_contigs[sample][contig_id]
            translated, offset, ragged = find_longest_orf(contig_seq)

            if vdj_aa and vdj_aa in translated:
                parts = translated.split(vdj_aa)
                leader = parts[0]
                leader_counter[leader] += 1

                # Get leader DNA
                if offset is not None:
                    leader_dna = contig_seq[offset : offset + len(leader) * 3]
                    leader_dna_counter[leader_dna] += 1

    if leader_counter:
        result[f"{chain}_leader_aa"] = leader_counter.most_common(1)[0][0]
    if leader_dna_counter:
        result[f"{chain}_leader_nt"] = leader_dna_counter.most_common(1)[0][0]


def _add_constant_regions(
    result: dict,
    constant_seqs: dict,
    row: pd.Series | None = None,
    sample_contigs: dict | None = None,
    allele_overrides: dict[str, str] | None = None,
    stops_nt: str = "TAATGA",
):
    """Add constant-region sequences and verify against the observed
    CellRanger contig start where possible.

    Strategy: pick the canonical TRAC / TRBC1 / TRBC2 via
    :func:`pick_canonical_constant` using ``c_gene`` (with J-gene
    fallback for β), splice the canonical AA in as the constant, and —
    if contigs are available — pull the first ~15 residues observed
    after the VDJ in the contig and verify they match the canonical
    start. Mismatches surface as ``{chain}_constant_qc`` warnings on
    the result; the canonical is still used so downstream
    cloning-construct outputs are at least biologically valid.

    ``allele_overrides`` controls per-gene allele picking (#113). Each
    key is a gene name (``"TRAC"`` / ``"TRBC1"`` / ``"TRBC2"``); each
    value is either ``"auto"`` (default — score every packaged allele
    against the contig translation and pick the best), or an explicit
    allele label (e.g. ``"01"`` or ``"03"``) to force a specific
    canonical. Missing keys default to ``"auto"``.

    ``stops_nt`` is the pre-validated, pre-concatenated NT string for
    stop codons (e.g. ``"TAATGA"`` for the default
    ``("TAA", "TGA")``). It's appended to every codon-optimized
    canonical NT — including the picker-override and J-junction
    branches that previously dropped stops (#116). Pass ``""`` for no
    stops.
    """
    allele_overrides = allele_overrides or {}
    for chain in ["alpha", "beta"]:
        c_gene = result.get(f"{chain}_c_gene")
        j_gene = result.get(f"{chain}_j_gene")
        canonical_name, canonical_aa = pick_canonical_constant(chain, c_gene, j_gene)

        # `constant_seqs` already holds codon-optimized back-translated
        # DNA for the canonical AA (built via `optimize_codons` with
        # the same stop_codons config — see #116); fall back to a
        # fresh `optimize_codons + stops_nt` build in case a caller
        # patched the dict with a missing key.
        canonical_nt_codon_opt = constant_seqs.get(
            canonical_name, optimize_codons(canonical_aa) + stops_nt
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name

        # If we have CellRanger contigs for this clone, splice the
        # donor's actual NT at the J→C junction into the front of the C
        # region NT (see :func:`_blend_constant_nt_with_contig` for the
        # full strategy). This keeps the assembled NT faithful to the
        # donor at the boundary, while the deep C region — invariant
        # across donors — stays codon-optimized for synthesis.
        contig_nt_past_vdj = (
            _extract_c_region_nt_from_contig(
                row, sample_contigs, result.get(f"vdj_{chain}_nt", ""), chain
            )
            if (row is not None and sample_contigs)
            else None
        )

        # Allele auto-detect (#113). After 2.3, ``HUMAN_TRBC2_AA``
        # defaults to the TRBC2*01-protein form (E at mature pos 9 —
        # major-allele in humans), but the auto-detect can pick the
        # TRBC2*03-protein form (K at pos 9) for the ~1% of donors who
        # carry it. ``allele_overrides[gene]`` controls the call:
        # ``"auto"`` scores all packaged alleles against contig
        # translation; an explicit allele label forces that allele.
        allele_choice = allele_overrides.get(canonical_name, "auto")
        allele_pool = HUMAN_CONSTANT_ALLELES.get(canonical_name, {})
        allele_called: str | None = None
        allele_score: float | None = None
        allele_alternatives: dict[str, float] = {}
        allele_called_reason: str | None = None

        # Translate the contig past the J→C junction codon up-front so
        # every code path can use it for #120's per-clone audit
        # (observed_constant_aa_start, allele_divergence_positions).
        contig_post_junction_nt = (
            contig_nt_past_vdj[3:] if contig_nt_past_vdj else ""
        )
        contig_aa_past_junction, _ = (
            translate_dna(contig_post_junction_nt)
            if contig_post_junction_nt else ("", False)
        )

        if not allele_pool:
            allele_called_reason = ALLELE_REASON_NO_ALLELE_POOL
        elif allele_choice != "auto":
            # Explicit user override. Validate against the pool
            # so a typo like ``trbc2_allele='3'`` (vs the right
            # ``'03'``) doesn't silently fall through to the default.
            if allele_choice in allele_pool:
                canonical_aa = allele_pool[allele_choice]
                canonical_nt_codon_opt = (
                    optimize_codons(canonical_aa) + stops_nt
                )
                allele_called = allele_choice
                allele_score = 1.0  # by definition; user chose it
                allele_called_reason = ALLELE_REASON_OVERRIDDEN
            else:
                allele_called_reason = ALLELE_REASON_INVALID_OVERRIDE
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: invalid allele override "
                    f"{canonical_name}*{allele_choice!s} (available: "
                    f"{sorted(allele_pool.keys())}); falling back to "
                    f"default allele."
                )
        elif not contig_nt_past_vdj:
            allele_called_reason = ALLELE_REASON_NO_CONTIG
            if sample_contigs:
                # contigs_dir was given but this clone has no contig
                # post-VDJ — surface it so #119 cohort aggregation can
                # split "no-contig" clones from "divergent" ones.
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — no contig coverage "
                    f"past VDJ for this clone."
                )
            # else: no contigs_dir at all — silent (user opted out).
        else:
            # Auto-detect. Pool exists, contig exists → run the picker.
            best_allele, best_score, all_scores = _pick_best_allele(
                contig_aa_past_junction, canonical_name, allele_pool,
            )
            if best_allele is None:
                # Picker punted. Either the contig translation was
                # empty (rare; handled below as "sparse") or
                # ``n_compared`` was below the floor for every allele.
                n_codons = len(contig_aa_past_junction) if contig_aa_past_junction else 0
                allele_called_reason = ALLELE_REASON_SPARSE_CONTIG
                # #118 ask 1: uniformly emit a qc_warning for sparse
                # contigs (Mode 2 was previously silent).
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — contig covers only "
                    f"{n_codons} codons past VDJ; need "
                    f"≥{_DEFAULT_MIN_PICKER_POSITIONS} for confident "
                    f"detection."
                )
            elif best_score <= 0:
                # Contig translates but disagrees at every comparable
                # position — likely a novel divergent allele or noise.
                allele_called_reason = ALLELE_REASON_DIVERGENT_CONTIG
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — contig disagrees with "
                    f"every packaged allele at every compared position."
                )
            else:
                candidate_aa = allele_pool[best_allele]
                # NT-level frame sanity check (#115).
                if not _verify_first_codon_frame(
                    contig_post_junction_nt, candidate_aa,
                ):
                    allele_called_reason = ALLELE_REASON_FRAME_ERROR
                    result.setdefault("qc_warnings", []).append(
                        f"{chain}: allele picker selected "
                        f"{canonical_name}*{best_allele} (AA score "
                        f"{best_score:.3f}) but the contig's first "
                        "codon past the junction does not translate "
                        "to that allele's expected first residue — "
                        "possible frame error upstream. Falling "
                        "back to default."
                    )
                else:
                    # #120 ask 4: tighten tolerance at allele-
                    # distinguishing positions. If the picked allele
                    # disagrees with the contig at one of the positions
                    # where alleles DIFFER from each other, that's a
                    # strong signal the donor has a novel variant —
                    # punt the *call* rather than silently committing.
                    #
                    # Note: "punt to NaN" means the allele audit
                    # columns go NaN. The constant region BYTES
                    # (canonical_aa / canonical_nt_codon_opt / the
                    # final ``{chain}_constant_*`` outputs) still
                    # use the default allele's canonical — so a
                    # downstream consumer that synthesizes from
                    # ``_constant_nt_optimized`` gets the default
                    # allele's residue at the polymorphic position
                    # even though the contig had a different one.
                    # Surfacing this is the point of
                    # ``_observed_constant_aa_start`` and the
                    # ``divergent_at_polymorphic_position`` reason;
                    # the user is expected to inspect those before
                    # ordering a construct.
                    polymorphic = _polymorphic_positions(allele_pool)
                    disagreed_at_polymorphic = sorted(
                        i for i in polymorphic
                        if i < len(contig_aa_past_junction)
                        and i < len(candidate_aa)
                        and contig_aa_past_junction[i] != candidate_aa[i]
                    )
                    if disagreed_at_polymorphic:
                        allele_called_reason = (
                            ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION
                        )
                        positions_str = ", ".join(
                            f"{i + 1}({candidate_aa[i]}/"
                            f"{contig_aa_past_junction[i]})"
                            for i in disagreed_at_polymorphic
                        )
                        result.setdefault("qc_warnings", []).append(
                            f"{chain}: allele not called — contig "
                            f"disagrees with best-fit "
                            f"{canonical_name}*{best_allele} at "
                            f"allele-distinguishing position(s) "
                            f"{positions_str} (canonical AA / observed AA). "
                            f"Possible novel allele; falling back to default."
                        )
                        # Keep allele_alternatives populated so the
                        # downstream cohort audit (#119) can see the
                        # picker's runner-up scores even when no call
                        # was committed.
                        allele_alternatives = all_scores
                    else:
                        canonical_aa = candidate_aa
                        canonical_nt_codon_opt = (
                            optimize_codons(canonical_aa) + stops_nt
                        )
                        allele_called = best_allele
                        allele_score = best_score
                        allele_alternatives = all_scores
                        allele_called_reason = ALLELE_REASON_AUTO_DETECTED

        # J→C junction residue, peeled from the contig — UNIFORMLY for
        # both chains (#105 in 2.0). The mature TCR mRNA in the donor
        # encodes a residue at the J→C splice that is NOT part of the
        # C-exon's canonical AA: the splice combines the J segment's
        # terminal partial codon (1-2 nt) with the C-exon's first 1-2
        # nt to spell the junction codon. For β chains the result is
        # universally E (GAG, across 118/118 audited B1 clones). For α
        # the junction is J-dependent — N (most common), Y, D, or H.
        #
        # Pre-2.0 versions baked an E into the stored HUMAN_TRBC1/2_AA;
        # α had no analogous handling and every α chain was 1 aa short.
        # Reading the junction codon from the contig handles both
        # chains the same way and lets β's per-clone junction NT be
        # the donor's actual bytes instead of codon-optimized AAC.
        junction_residue: str | None = None
        if contig_nt_past_vdj and len(contig_nt_past_vdj) >= 3:
            junction_codon = contig_nt_past_vdj[:3]
            translated = CODON_TABLE.get(junction_codon, "X")
            if translated and translated not in {"*", "X"}:
                junction_residue = translated
                canonical_aa = junction_residue + canonical_aa
                # Prepend the codon for the junction residue. Using
                # `back_translate` (single-best codon, no motif
                # logic) is intentional here: it's one codon at the
                # start of the construct, and `canonical_nt_codon_opt`
                # already has its own motif-optimized body. Stops are
                # already on canonical_nt_codon_opt's tail — don't
                # re-append.
                canonical_nt_codon_opt = (
                    back_translate(junction_residue) + canonical_nt_codon_opt
                )
            else:
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: junction codon {junction_codon!r} translates "
                    f"to {translated!r}; not prepended to "
                    f"{chain}_constant_aa. The assembled {chain} may be 1 "
                    "aa shorter than the donor's expressed protein at the "
                    "J→C junction."
                )
        elif sample_contigs:
            # Contigs were loaded but this clone's contig didn't
            # extend past VDJ — assembled chain is 1 aa short.
            result.setdefault("qc_warnings", []).append(
                f"{chain}: no contig coverage past VDJ; "
                f"{chain}_constant_aa will be missing the J→C junction "
                f"residue (1 aa shorter than the donor's expressed "
                f"{chain} protein)."
            )
        # else: no contigs_dir provided at all — silent. The user opted
        # out of contig-derived NT entirely; assembly uses bare mature
        # canonical, same as if the clone had no contig.
        result[f"{chain}_junction_residue"] = junction_residue

        # Record the allele-picker audit (#113, expanded #118/#120).
        result[f"{chain}_allele_called"] = (
            f"{canonical_name}*{allele_called}" if allele_called else None
        )
        result[f"{chain}_allele_score"] = allele_score
        # #118: machine-readable reason for the call (or no-call).
        result[f"{chain}_allele_called_reason"] = allele_called_reason
        if allele_alternatives:
            result[f"{chain}_allele_alternatives"] = ";".join(
                f"{canonical_name}*{a}:{s:.3f}"
                for a, s in sorted(
                    allele_alternatives.items(), key=lambda kv: -kv[1]
                )
            )
        else:
            result[f"{chain}_allele_alternatives"] = None

        # #120 ask 1: raw observed AA past the J→C junction codon.
        # Lets downstream consumers see the donor's actual residues
        # regardless of whether the picker matched a packaged allele.
        # Capped at 15 residues — past that contig coverage is rare
        # and the residues are dominated by invariant regions.
        result[f"{chain}_observed_constant_aa_start"] = (
            contig_aa_past_junction[:15] if contig_aa_past_junction else None
        )

        # #120 ask 2: per-position divergences between the called
        # allele and the observation. Format: "3:N->K;7:V->I". None
        # when no allele was called OR when there's no divergence.
        if allele_called is not None and contig_aa_past_junction:
            called_aa = allele_pool[allele_called]
            divs = _divergence_positions(called_aa, contig_aa_past_junction)
            result[f"{chain}_allele_divergence_positions"] = (
                _format_divergence_positions(divs)
            )
        else:
            result[f"{chain}_allele_divergence_positions"] = None

        result[f"{chain}_constant_aa"] = canonical_aa
        blended_nt, blend_debug = _blend_constant_nt_with_contig(
            contig_nt_past_vdj or "",
            canonical_aa,
            canonical_nt_codon_opt,
        )

        # New (#116) NT triad. Three views of the constant-region NT:
        # (1) ``_contig`` — pure CellRanger contig bytes past the J→C
        #     junction. **Truncated at contig coverage** — typically
        #     covers 30-150 nt of the ~430-530 nt C region (i.e.,
        #     does NOT reach the stop). Useful for QC, allele
        #     detection, and any downstream that wants the donor's
        #     actual NT with no canonical splice. ``None`` when no
        #     contig is available for the clone.
        # (2) ``_optimized`` — codon-optimized canonical CDS (+stops).
        #     Same byte string for every donor carrying the same
        #     picked allele; what synthesis pipelines order.
        # (3) ``_constant_nt`` — the legacy column, kept for the
        #     assembly-aware blend (donor bytes where they agree
        #     with canonical AA, ``_optimized`` for the rest).
        #     Unchanged behavior pre/post #116 for callers that
        #     didn't ask for the new views.
        result[f"{chain}_constant_nt_contig"] = (
            contig_nt_past_vdj if contig_nt_past_vdj else None
        )
        result[f"{chain}_constant_nt_optimized"] = canonical_nt_codon_opt
        result[f"{chain}_constant_nt"] = blended_nt

        # Verify the observed contig C-region start against canonical.
        observed = (
            _extract_c_region_start_from_contig(
                row, sample_contigs, result.get(f"vdj_{chain}_aa", ""), chain
            )
            if (row is not None and sample_contigs)
            else None
        )
        if observed is None:
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} [{blend_debug['source']}]"
            )
        elif verify_canonical_constant_start(observed, canonical_aa):
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (contig-verified) "
                f"[{blend_debug['source']}]"
            )
        else:
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (UNVERIFIED — start mismatch) "
                f"[{blend_debug['source']}]"
            )
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant start mismatch: observed "
                f"{observed!r} differs from canonical {canonical_name} "
                f"(expected start {canonical_aa[:len(observed)]!r}). "
                "Using canonical anyway; c_gene assignment may be wrong."
            )
        if blend_debug["aa_mismatch_at"] is not None:
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant: contig diverged from canonical at AA "
                f"position {blend_debug['aa_mismatch_at']}; switched to "
                f"codon-optimized canonical from that point."
            )
        if blend_debug.get("partial_codon_dropped"):
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant: contig provided a partial codon at the "
                "C-region boundary, but its bytes were incompatible with "
                "the canonical residue; dropped donor fidelity at that "
                "codon and used canonical."
            )


def _extract_c_region_nt_from_contig(
    row: pd.Series,
    sample_contigs: dict,
    vdj_nt: str,
    chain: str,
) -> str | None:
    """Pull NT bytes observed immediately after the VDJ end in the CellRanger
    contig — the donor's real, unedited sequence at the J→C boundary.

    Returns the contig NT starting one nt after the last nt of ``vdj_nt`` in
    the contig. The returned NT is **frame-aligned** to the next codon: in
    the assembled construct, this is the first nt of the J→C junction codon
    (whose first 1-2 nt are biologically encoded by the J's terminal partial
    codon and whose remaining nt come from the C exon's start). The donor's
    actual junction codon and the first several downstream C-region codons
    live here, so :func:`_blend_constant_nt_with_contig` can splice them in
    before falling back to a codon-optimized canonical for the rest.

    Notes
    -----
    * ``vdj_nt`` is the frame-aligned (no +1 overshoot) lowercase
      ``vdj_{chain}_nt`` written by the assembler after #91. The original
      CellRanger contig (``filtered_contig.fasta`` parsed into ``sample_contigs``)
      retains the J→C boundary nt that were trimmed during clonotype-aggregation.
    * Returns ``None`` when no contig contains ``vdj_nt`` exactly — usually a
      sign that the canonical-from-contigs path isn't available for this
      clonotype (e.g., contigs_dir wasn't loaded, or VDJ_nt was re-derived
      from a different cell than the one whose contig we have).
    * When multiple contigs contain ``vdj_nt`` (multi-sample clonotype), the
      most-common post-VDJ NT is returned.
    * **Cross-source hazard mitigation (#94 archetype).** ``contig_ids`` is
      the semicolon-joined list of EVERY cell's contig — not just the
      representative cell whose ``vdj_nt`` was kept by the clonotype
      aggregator. We rely on ``contig_seq.find(vdj_nt)`` to filter: contigs
      whose own VDJ doesn't contain the representative ``vdj_nt`` byte-for-
      byte are skipped. Cells with cell-to-cell polymorphism inside the VDJ
      therefore contribute nothing to the post-VDJ consensus, and the
      blended NT stays sourced from a coherent biological signal. Don't
      remove the ``find`` filter without re-deriving the safety argument.
    """
    if not vdj_nt:
        return None
    samples = str(row.get("samples", "")).split(";")
    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return None

    contig_ids = str(row[contig_col]).split(";")
    candidates: Counter = Counter()
    for sample in samples:
        if sample not in sample_contigs:
            continue
        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue
            contig_seq = sample_contigs[sample][contig_id]
            # Locate vdj_nt as an exact substring. Try both the contig as
            # given and its reverse complement (CellRanger contigs are
            # already oriented, but defensive).
            idx = contig_seq.find(vdj_nt)
            if idx >= 0:
                after = contig_seq[idx + len(vdj_nt):]
                if after:
                    candidates[after] += 1
    if candidates:
        return candidates.most_common(1)[0][0]
    return None


def _blend_constant_nt_with_contig(
    contig_nt_past_vdj: str,
    canonical_aa: str,
    canonical_nt_codon_opt: str,
) -> tuple[str, dict]:
    """Splice contig-derived NT for the J→C junction into a codon-optimized
    canonical C region NT.

    The result respects three principles, listed in order of precedence:

    1. **Junction-residue codon comes from the donor where possible.**
       The first nt(s) past the (trimmed) VDJ are biologically encoded by
       the J segment's terminal partial codon; the next nt(s) come from
       the C exon. Together they spell the J→C junction codon. Using the
       contig's actual bytes here keeps the assembled NT faithful to the
       donor's real sequence at this site.
    2. **Incomplete codons at the contig edge are completed using the
       canonical reference.** If the contig provides 1-2 nt of an
       otherwise-incomplete codon at its 3' end, we extend with the
       canonical codon's bytes (preferring the codon that matches both
       the contig prefix AND the canonical residue, so the protein output
       is preserved).
    3. **Everything past the contig coverage uses codon-optimized canonical
       NT.** Synthesis vendors prefer codon-optimized constructs, and the
       deep C region is invariant across donors — there's no donor-specific
       NT to preserve there.

    Stops at the first AA mismatch between the contig-translated codon and
    the canonical AA at that position — at that point the contig and
    canonical disagree about the protein, and we trust the canonical (the
    contig may be carrying a sequencing error or a non-canonical variant).
    The mismatch surfaces in the returned ``debug`` dict so callers can
    emit a QC warning.

    Parameters
    ----------
    contig_nt_past_vdj : str
        NT extracted from the contig starting one nt past the VDJ. Frame
        offset 0 — i.e., position 0 is the first nt of the J→C junction
        codon. May be empty (no contig coverage); in that case the entire
        ``canonical_nt_codon_opt`` is returned.
    canonical_aa : str
        Canonical AA sequence of the constant region (from the packaged
        FASTA, post-#100). Position 0 corresponds to the J→C junction
        residue (for β chains) or the first mature residue (for α).
    canonical_nt_codon_opt : str
        Codon-optimized back-translated NT for ``canonical_aa`` with
        stop codons already appended at the tail. Length
        ``3 * len(canonical_aa) + 3 * n_stops``. (Pre-#116 this was
        single-stop; from 2.4 onward it defaults to dual-stop.)

    Returns
    -------
    tuple[str, dict]
        ``(blended_nt, debug)`` where ``blended_nt`` is the assembled C
        region NT with stops preserved at the tail (caller does NOT
        re-append), and ``debug`` carries:

        - ``n_contig_codons``: number of full codons spliced from contig.
        - ``partial_codon_completed``: True when a partial codon at the
          contig 3' edge was extended with canonical bytes to spell the
          canonical residue (donor fidelity preserved in the partial nt).
        - ``partial_codon_dropped``: True when a partial codon WAS present
          but its bytes were incompatible with every canonical codon
          for the expected residue — partial bytes discarded and the
          full canonical codon used. Distinct from ``partial_codon_completed=False``
          when there simply was no partial codon to begin with.
        - ``aa_mismatch_at``: 1-indexed canonical AA position where the
          contig translation first disagreed with canonical (or ``None``
          if no mismatch was hit before the contig ran out).
        - ``source``: human-readable summary for the
          ``{chain}_constant_source`` column.

    Length invariant: ``len(blended_nt) == len(canonical_nt_codon_opt)``
    always, regardless of which path was taken. Equivalent to
    ``3 * len(canonical_aa) + 3 * n_stops`` since the input includes
    stops. Tested in ``TestBlendConstantNtWithContig``.
    """
    debug = {
        "n_contig_codons": 0,
        "partial_codon_completed": False,
        "partial_codon_dropped": False,
        "aa_mismatch_at": None,
        "source": "canonical-codon-opt",
    }
    if not contig_nt_past_vdj or not canonical_aa:
        return canonical_nt_codon_opt, debug

    n_full_codons_in_contig = len(contig_nt_past_vdj) // 3
    partial_nt = contig_nt_past_vdj[n_full_codons_in_contig * 3 :]

    # Find longest contig prefix whose per-codon translation matches the
    # canonical AA byte-for-byte. First mismatch ends the contig-derived
    # run; everything past that point uses canonical codon-optimized NT.
    n_matching = 0
    for i in range(min(n_full_codons_in_contig, len(canonical_aa))):
        codon = contig_nt_past_vdj[i * 3 : (i + 1) * 3]
        translated = CODON_TABLE.get(codon, "X")
        if translated == canonical_aa[i]:
            n_matching += 1
        else:
            debug["aa_mismatch_at"] = i + 1
            break

    contig_part = contig_nt_past_vdj[: n_matching * 3]
    debug["n_contig_codons"] = n_matching

    # Try to complete a partial codon at the contig 3' edge — only if
    # every full codon matched (otherwise we've already switched to
    # canonical for this position).
    completed_codon = ""
    canonical_offset = n_matching
    if (
        n_matching == n_full_codons_in_contig
        and partial_nt
        and canonical_offset < len(canonical_aa)
    ):
        expected_aa = canonical_aa[canonical_offset]
        # All codons that translate to expected_aa AND start with partial_nt.
        compatible = [
            cdn
            for cdn, aa in CODON_TABLE.items()
            if aa == expected_aa and cdn.startswith(partial_nt)
        ]
        if compatible:
            # If the codon-optimized choice is among the compatible ones,
            # prefer it (keeps downstream consistency). Otherwise pick any.
            preferred = HUMAN_PREFERRED_CODONS.get(expected_aa)
            if preferred in compatible:
                completed_codon = preferred
            else:
                completed_codon = compatible[0]
            canonical_offset += 1
            debug["partial_codon_completed"] = True
        else:
            # Contig partial bytes incompatible with canonical residue —
            # discard the partial bytes and let canonical drive this codon.
            # Recorded separately from "no partial at all" so the audit
            # trail can flag clones where donor fidelity at the boundary
            # WAS available but had to be dropped.
            debug["partial_codon_dropped"] = True

    canonical_rest = canonical_nt_codon_opt[canonical_offset * 3 :]
    blended = contig_part + completed_codon + canonical_rest

    if n_matching > 0 and debug["aa_mismatch_at"] is None:
        suffix_parts = []
        if debug["partial_codon_completed"]:
            suffix_parts.append("partial completed")
        elif debug["partial_codon_dropped"]:
            suffix_parts.append("partial dropped — incompatible with canonical AA")
        suffix = (", " + "; ".join(suffix_parts)) if suffix_parts else ""
        debug["source"] = f"contig({n_matching} codons{suffix}) + canonical-codon-opt"
    elif n_matching > 0:
        debug["source"] = (
            f"contig({n_matching} codons) + canonical-codon-opt "
            f"(switched at AA mismatch pos {debug['aa_mismatch_at']})"
        )
    elif debug["partial_codon_completed"]:
        # No full codons, but a partial was completed — still donor-derived
        # at the J→C boundary even though nothing matched at the codon level.
        debug["source"] = "contig(partial completed) + canonical-codon-opt"
    elif debug["partial_codon_dropped"]:
        # No full codons AND partial dropped — purely canonical output, but
        # surface the dropped partial so the audit trail records that donor
        # fidelity was available and discarded.
        debug["source"] = (
            "canonical-codon-opt (partial dropped — incompatible with canonical AA)"
        )
    # else: no contig codons spliced and no partial — keep default
    #       "canonical-codon-opt"

    return blended, debug


def _extract_c_region_start_from_contig(
    row: pd.Series,
    sample_contigs: dict,
    vdj_aa: str,
    chain: str,
    n_aa: int = 15,
) -> str | None:
    """Pull the first ``n_aa`` residues observed after the VDJ in the
    CellRanger contig, if extractable.

    Mirrors the logic in :func:`_extract_leader_from_contigs_single` —
    longest-ORF translation, split on the VDJ AA, take what follows.
    Returns ``None`` when no contig contains the VDJ AA or the
    post-VDJ tail is shorter than ``n_aa``.
    """
    if not vdj_aa:
        return None
    samples = str(row.get("samples", "")).split(";")
    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return None

    contig_ids = str(row[contig_col]).split(";")
    starts = Counter()
    for sample in samples:
        if sample not in sample_contigs:
            continue
        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue
            contig_seq = sample_contigs[sample][contig_id]
            translated, _, _ = find_longest_orf(contig_seq)
            if vdj_aa in translated:
                _, after = translated.split(vdj_aa, 1)
                if len(after) >= n_aa:
                    starts[after[:n_aa]] += 1
    if starts:
        return starts.most_common(1)[0][0]
    return None


def _add_constant_from_row(row: pd.Series, result: dict):
    """Add constant region sequences directly from row columns if present.

    Also writes ``{chain}_c_gene_canonical`` so the from-data path
    surfaces the same audit column as the canonical path — without it,
    downstream code (``validate_sequences``, ``fix_jc_parity``,
    ``_count_trbc_overrides``) only sees the raw CellRanger
    ``{chain}_c_gene`` and the J-family override / parity checks
    silently degrade.
    """
    for chain in ["alpha", "beta"]:
        aa_col = f"{chain}_constant_aa"
        nt_col = f"{chain}_constant_nt"

        aa_val = row.get(aa_col) if aa_col in row else None
        nt_val = row.get(nt_col) if nt_col in row else None

        if pd.notna(aa_val):
            result[f"{chain}_constant_aa"] = aa_val
        if pd.notna(nt_val):
            result[f"{chain}_constant_nt"] = nt_val
            # If AA not provided, translate from NT
            if pd.isna(aa_val):
                const_aa, _ = translate_dna(str(nt_val))
                result[f"{chain}_constant_aa"] = const_aa

        # Record the canonical gene name so the from-data path produces
        # the same downstream audit column as the canonical path.
        canonical_name, _ = pick_canonical_constant(
            chain,
            c_gene=row.get(f"{chain}_c_gene"),
            j_gene=row.get(f"{chain}_j_gene"),
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name


def _build_full_sequences(
    result: dict,
    include_alpha_leader: bool,
    include_beta_leader: bool,
    include_constant: bool,
):
    """Build complete sequences from parts.

    Emits the standard ``full_{chain}_aa`` / ``full_{chain}_nt`` plus
    (when constants are included) the NT triad introduced in #116:

    * ``full_{chain}_nt_contig`` — leader + VDJ + pure CellRanger
      contig bytes for the constant. **Truncated at contig coverage**
      (donor's contig typically reads ~30-150 nt into the C region,
      not the full ~430-530 nt). ``None`` when no contig is
      available — there is no canonical fallback for this view by
      design.
    * ``full_{chain}_nt_optimized`` — leader + VDJ + codon-optimized
      canonical constant (+stops). Always present when constants
      are included. Synthesis-ready.
    * ``full_{chain}_nt`` (legacy / blend) — leader + VDJ +
      donor bytes where they agree with canonical AA, ``_optimized``
      for the rest.
    """
    include_leader_map = {"alpha": include_alpha_leader, "beta": include_beta_leader}

    for chain in ["alpha", "beta"]:
        parts_aa: list[str] = []
        parts_nt_blend: list[str] = []
        parts_nt_contig: list[str | None] = []
        parts_nt_optimized: list[str] = []
        include_leader = include_leader_map[chain]

        if include_leader and f"{chain}_leader_aa" in result:
            parts_aa.append(result[f"{chain}_leader_aa"])
        if include_leader and f"{chain}_leader_nt" in result:
            leader_nt = result[f"{chain}_leader_nt"]
            parts_nt_blend.append(leader_nt)
            parts_nt_contig.append(leader_nt)
            parts_nt_optimized.append(leader_nt)

        if f"vdj_{chain}_aa" in result:
            parts_aa.append(result[f"vdj_{chain}_aa"])
        if f"vdj_{chain}_nt" in result:
            vdj_nt = result[f"vdj_{chain}_nt"]
            parts_nt_blend.append(vdj_nt)
            parts_nt_contig.append(vdj_nt)
            parts_nt_optimized.append(vdj_nt)

        if include_constant and f"{chain}_constant_aa" in result:
            parts_aa.append(result[f"{chain}_constant_aa"])
        if include_constant and f"{chain}_constant_nt" in result:
            parts_nt_blend.append(result[f"{chain}_constant_nt"])
            # _contig is None when contig had no post-VDJ coverage —
            # propagate the None to full_{chain}_nt_contig so callers
            # can filter on it without parsing.
            parts_nt_contig.append(result.get(f"{chain}_constant_nt_contig"))
            parts_nt_optimized.append(
                result.get(f"{chain}_constant_nt_optimized", "")
            )

        if parts_aa:
            result[f"full_{chain}_aa"] = "".join(parts_aa)
        if parts_nt_blend:
            result[f"full_{chain}_nt"] = "".join(parts_nt_blend)
        if include_constant:
            # Contig view: any None constant means the donor's bytes
            # don't reach the C region; surface None instead of a
            # silently-truncated string.
            if all(p is not None for p in parts_nt_contig) and parts_nt_contig:
                result[f"full_{chain}_nt_contig"] = "".join(parts_nt_contig)
            else:
                result[f"full_{chain}_nt_contig"] = None
            if parts_nt_optimized:
                result[f"full_{chain}_nt_optimized"] = "".join(parts_nt_optimized)


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
    * ``single_chain_nt_contig`` — built from ``full_*_nt_contig``.
      ``None`` for any clone where either chain's contig didn't cover
      its constant region (``_contig`` columns are ``None`` there).
      **Truncated at contig coverage** — useful for QC and inspecting
      the donor's actual bytes, but not a complete CDS and not
      synthesis-ready. Reach for ``_optimized`` if you want a
      complete cassette.

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
        if seq.endswith("*"):
            return seq[:-1]
        return seq

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
        for suffix, out_col in (
            ("", "single_chain_nt"),
            ("_optimized", "single_chain_nt_optimized"),
            ("_contig", "single_chain_nt_contig"),
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


class ValidationMessage(str):
    """A validation/autocorrect message tagged with its clone index
    and severity.

    Subclasses :class:`str` so existing callers that treat the return
    value of :func:`validate_sequences` as a list of strings continue
    to work (substring tests, joins, length checks). New callers can
    filter on ``.severity`` and route on ``.idx`` instead of parsing
    the ``"Clone {idx}: ..."`` prefix back out of the string.

    Severities:

    - ``"load_bearing"`` — a structural failure (CDR3 missing, length
      out of window, premature stop, etc.) that justifies raising in
      strict mode.
    - ``"informational"`` — a "didn't have enough info to check" note;
      never raises.
    - ``"autocorrect"`` — emitted by :func:`fix_jc_parity` (and by
      :func:`validate_sequences` when ``fix=True``) to record an
      in-place correction that was applied.

    **Caveat:** any :class:`str` operation that returns a new value
    (``+``, slicing, ``.replace``, f-string interpolation,
    ``json.dumps``, etc.) returns a plain :class:`str` and drops the
    ``.idx`` / ``.severity`` attributes. Read the metadata directly
    off the original instance — don't massage the text first. ``copy``
    and ``pickle`` do preserve the attributes.
    """
    idx: object
    severity: str

    def __new__(cls, text: str, *, idx: object = None, severity: str = "load_bearing"):
        obj = super().__new__(cls, text)
        obj.idx = idx
        obj.severity = severity
        return obj


def fix_jc_parity(df: pd.DataFrame) -> list[ValidationMessage]:
    """Autocorrect β J→C parity mismatches in-place (#89).

    Locus structure guarantees TRBJ1-* rearranges with TRBC1 and TRBJ2-*
    with TRBC2. When the called C-gene disagrees, J is authoritative —
    overwrite ``beta_c_gene_canonical`` to match the J family.

    This only updates ``beta_c_gene_canonical`` (the column
    :func:`validate_sequences` and :func:`assemble_qc_report` consult
    first); the raw ``beta_c_gene`` is left alone as a record of what
    CellRanger originally called. ``full_beta_aa`` itself is not
    regenerated — use this when you have an already-assembled frame
    whose J/C parity annotation is wrong but the assembled sequence
    is correct (e.g. data from an older assembler that trusted
    CellRanger's TRBC call).

    Returns a list of :class:`ValidationMessage` with
    ``severity="autocorrect"``, one per row that was changed. The
    input DataFrame is mutated.
    """
    messages: list[ValidationMessage] = []
    if "beta_j_gene" not in df.columns:
        return messages
    # Make sure the canonical column exists with an object dtype that
    # accepts string writes (a float64 NaN column would emit a
    # FutureWarning on the first `df.at[idx, col] = "TRBC1"`).
    if "beta_c_gene_canonical" not in df.columns:
        df["beta_c_gene_canonical"] = df.get("beta_c_gene")
    df["beta_c_gene_canonical"] = df["beta_c_gene_canonical"].astype(object)

    # Vectorise the family lookup so we iterate only over the rows
    # that actually have a parity rule to enforce. On a 100k-row
    # cohort with ~0.1% override rate this skips ~99% of the loop
    # iterations.
    j_family_s = df["beta_j_gene"].apply(_beta_j_family)
    expected_s = j_family_s.map(BETA_JC_PARITY)
    candidates = expected_s.dropna()
    if candidates.empty:
        return messages

    for idx in candidates.index:
        row = df.loc[idx]
        current = _resolve_c_gene(row, "beta")
        current_base = _beta_c_base(current)
        expected_c = candidates.loc[idx]
        if current_base and current_base != expected_c:
            df.at[idx, "beta_c_gene_canonical"] = expected_c
            messages.append(
                ValidationMessage(
                    f"Clone {idx}: autocorrected β c_gene {current!r} → "
                    f"{expected_c} based on {row['beta_j_gene']} (J family rules)",
                    idx=idx,
                    severity="autocorrect",
                )
            )
    return messages


def validate_sequences(
    df: pd.DataFrame, strict: bool = False, fix: bool = False
) -> list[ValidationMessage]:
    """Validate assembled sequences end-to-end per row (#67, #68).

    Per-chain checks on each ``full_{chain}_aa`` row:

    - Length window (200–450 aa).
    - CDR3 substring present in the assembled sequence.
    - Constant region's canonical C-terminus (from
      :data:`CONSTANT_REGION_ENDINGS`).
    - Constant region's canonical start (first 8+ residues of
      :data:`HUMAN_CONSTANT_REGIONS_AA`).
    - Constant length floor (:data:`CONSTANT_AA_FLOOR`).
    - No premature stop codon (``*`` mid-chain).
    - Methionine start (when a leader was included).
    - Standard residue alphabet only (:data:`_VALID_AA_CHARS`).
    - Byte-for-byte equality of ``leader + vdj + constant`` and the
      assembled ``full_{chain}_aa``, when all three parts are present.

    Cross-row checks:

    - β-chain J→C in-cis parity: ``TRBJ1*`` pairs with TRBC1,
      ``TRBJ2*`` with TRBC2.
    - Single-chain integrity: linker present in ``single_chain_aa``.
    - Per-row ``qc_warnings`` (stashed by ``_add_constant_regions``)
      are surfaced.

    Returns a flat list of :class:`ValidationMessage` (a ``str``
    subclass with ``.idx`` and ``.severity`` attributes). When
    ``strict`` is True, raises :class:`TCRsiftValidationError` if any
    load-bearing checks fail. Informational notes ("didn't have enough
    info to check") are returned but never raise — distinguishing
    those from a silent pass was the gap that hid #66.

    When ``fix`` is True, autocorrect β J→C parity mismatches in-place
    via :func:`fix_jc_parity` before validating (#89). Autocorrection
    notes are prepended to the returned messages and the input frame
    is mutated *unconditionally* — the mutation survives even when a
    later check fails and ``strict=True`` raises.
    """
    load_bearing: list[ValidationMessage] = []
    informational: list[ValidationMessage] = []
    autocorrect_notes: list[ValidationMessage] = []

    def _lb(idx, text):
        load_bearing.append(
            ValidationMessage(f"Clone {idx}: {text}", idx=idx, severity="load_bearing")
        )

    def _info(idx, text):
        informational.append(
            ValidationMessage(f"Clone {idx}: {text}", idx=idx, severity="informational")
        )

    if fix:
        autocorrect_notes = fix_jc_parity(df)

    # 1. Per-chain, per-row shape and content checks.
    for chain in ["alpha", "beta"]:
        col = f"full_{chain}_aa"
        if col not in df.columns:
            continue
        for idx, row in df.iterrows():
            seq = row.get(col, "")
            if not seq or not isinstance(seq, str):
                continue

            if len(seq) < 200:
                _lb(idx, f"{chain} chain too short ({len(seq)} aa)")
            elif len(seq) > 450:
                _lb(idx, f"{chain} chain too long ({len(seq)} aa)")

            # Standard residue alphabet only.
            invalid = sorted(set(seq) - _VALID_AA_CHARS)
            if invalid:
                _lb(idx,
                    f"full_{chain}_aa has invalid residues {''.join(invalid)!r}")

            # Premature stop codon: ``*`` anywhere except the final position.
            body = seq.rstrip("*")
            if "*" in body:
                pos = body.index("*")
                _lb(idx,
                    f"full_{chain}_aa has premature stop "
                    f"at position {pos} of {len(seq)}")

            # Methionine start (only when leader was included).
            leader = row.get(f"{chain}_leader_aa")
            if isinstance(leader, str) and leader and not seq.startswith("M"):
                _lb(idx,
                    f"full_{chain}_aa doesn't start with M "
                    f"(starts with {seq[:5]!r}; leader present)")

            cdr3_col = f"CDR3_{chain}"
            if cdr3_col in row:
                cdr3 = row[cdr3_col]
                if isinstance(cdr3, str) and cdr3 and cdr3 not in seq:
                    _lb(idx,
                        f"CDR3_{chain}={cdr3!r} not found in full sequence")

            # Constant region length floor (separate from full-chain
            # length check — catches truncations that hid behind a
            # long leader).
            const = row.get(f"{chain}_constant_aa")
            if isinstance(const, str) and const:
                floor = CONSTANT_AA_FLOOR[chain]
                if len(const) < floor:
                    _lb(idx,
                        f"{chain}_constant_aa too short "
                        f"({len(const)} aa, floor {floor})")

            # Byte-for-byte: full == leader + vdj + constant when all
            # parts are available. Catches dropped/added residues
            # during assembly (the failure mode no other check finds).
            vdj = row.get(f"vdj_{chain}_aa")
            if (
                isinstance(leader, str) and leader
                and isinstance(vdj, str) and vdj
                and isinstance(const, str) and const
            ):
                expected = leader + vdj + const
                if expected != seq:
                    _lb(idx,
                        f"full_{chain}_aa != leader+vdj+constant "
                        f"(len {len(seq)} vs {len(expected)})")

            # Canonical-ending / canonical-start checks. Prefer
            # ``{chain}_c_gene_canonical`` (written by
            # ``_add_constant_regions``); fall back to the raw
            # CellRanger ``{chain}_c_gene`` call. If neither resolves
            # to a known canonical gene name, emit an informational
            # "did-not-check" note so the user knows validation was
            # incomplete.
            c_gene = _resolve_c_gene(row, chain)
            # Strip any allele suffix (e.g. "TRBC1*01" → "TRBC1").
            c_gene_base = c_gene.split("*")[0]
            if c_gene_base in CONSTANT_REGION_ENDINGS:
                expected_end = CONSTANT_REGION_ENDINGS[c_gene_base]
                if not seq.endswith(expected_end):
                    _lb(idx,
                        f"{chain} doesn't end with canonical "
                        f"{c_gene_base} C-terminus ({expected_end!r}); "
                        f"got {seq[-len(expected_end):]!r}")
                # Prefer the row's actual ``{chain}_constant_aa`` as the
                # expected start. After 2.0 (#105) this carries the
                # per-clone junction residue prepend (e.g.
                # ``NIQNPDPAVY...`` for α with N junction;
                # ``EDLNKVFP...`` for β with universal-E junction).
                # Falling back to ``HUMAN_CONSTANT_REGIONS_AA[base]``
                # — the bare-mature canonical — was the #107 regression:
                # observed includes the junction; expected didn't.
                expected_const = row.get(f"{chain}_constant_aa")
                if not isinstance(expected_const, str) or not expected_const:
                    expected_const = HUMAN_CONSTANT_REGIONS_AA.get(c_gene_base)
                if expected_const and not verify_canonical_constant_start(
                    _expected_constant_start_from_full(seq, row, chain),
                    expected_const,
                    min_match=8,
                ):
                    _lb(idx,
                        f"{chain} constant start doesn't match canonical "
                        f"{c_gene_base} (expected start "
                        f"{expected_const[:15]!r})")
            else:
                _info(idx,
                    f"{chain} c_gene={c_gene!r} not in "
                    "CONSTANT_REGION_ENDINGS — canonical ending/start unverifiable")

    # 2. β-chain J→C in-cis parity disagreement (#67).
    if "beta_j_gene" in df.columns:
        for idx, row in df.iterrows():
            j_gene = row.get("beta_j_gene", "")
            c_gene = _resolve_c_gene(row, "beta")
            j_family = _beta_j_family(j_gene)
            c_base = _beta_c_base(c_gene)
            expected_c = BETA_JC_PARITY.get(j_family) if j_family else None
            if expected_c and c_base and c_base != expected_c:
                _lb(idx,
                    f"β J→C parity mismatch — "
                    f"{j_gene} (family {j_family}) should pair with "
                    f"{expected_c}, got {c_base}")

    # 3. Single-chain (β-linker-α) construct integrity. Three checks:
    #    a) the linker appears in ``single_chain_aa`` exactly once,
    #    b) ``single_chain_aa == β.rstrip('*') + linker + α`` byte-for-byte
    #       (catches dropped residues or wrong concatenation order),
    #    c) if the linker AA matches a known 2A peptide name, it should
    #       be the canonical sequence from :data:`LINKERS`.
    if "single_chain_aa" in df.columns and "linker" in df.columns:
        for idx, row in df.iterrows():
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            if not isinstance(sc, str) or not sc:
                continue
            if not isinstance(linker, str) or not linker:
                continue
            count = sc.count(linker)
            if count == 0:
                _lb(idx, f"single_chain_aa missing linker {linker!r}")
            elif count > 1:
                _lb(idx,
                    f"single_chain_aa contains linker "
                    f"{linker!r} {count} times (expected 1)")

            beta = row.get("full_beta_aa")
            alpha = row.get("full_alpha_aa")
            if (
                isinstance(beta, str) and beta
                and isinstance(alpha, str) and alpha
                and count == 1
            ):
                expected_sc = beta.rstrip("*") + linker + alpha
                if sc != expected_sc:
                    _lb(idx,
                        f"single_chain_aa != β+linker+α "
                        f"(len {len(sc)} vs {len(expected_sc)})")

    # 4. NT → AA round-trip invariants (#91). Every NT column the
    # assembler emits must back-translate to its corresponding AA
    # column. This is the integration check that catches splice-time
    # frame bugs — the original #91 (+1 trailing nt in VDJ_*_nt was
    # not trimmed before concatenation, frame-shifting everything
    # past the VDJ→C boundary) would have been caught here on day one.
    _validate_nt_aa_roundtrip(df, _lb)

    # 5. Surface any per-row qc_warnings the assembler stashed (e.g.
    # contig-vs-canonical start mismatch detected during assembly).
    if "qc_warnings" in df.columns:
        for idx, qcs in df["qc_warnings"].items():
            if isinstance(qcs, list):
                for msg in qcs:
                    _lb(idx, str(msg))

    all_messages = autocorrect_notes + load_bearing + informational
    if strict and load_bearing:
        raise TCRsiftValidationError(
            f"Sequence validation failed ({len(load_bearing)} load-bearing issues):\n  "
            + "\n  ".join(load_bearing[:10])
            + ("\n  ..." if len(load_bearing) > 10 else "")
        )
    return all_messages


@dataclass
class AssemblyQCCheck:
    """One QC check's tally across a frame.

    Attributes:
        name: machine-friendly identifier (e.g. ``"alpha_terminal_residue"``).
        label: human-readable label for the log/plot (e.g. ``"α terminal residue"``).
        chain: ``"alpha"``, ``"beta"``, ``"single_chain"``, or ``None`` for
            cross-chain checks.
        passed: row count that passed this check.
        total: row count this check ran against.
        median_value: optional median (used for length-style checks); unitless.
        unit: optional unit for ``median_value`` (e.g. ``"aa"``).
        examples: up to a few failure exemplars as ``"clone {idx}: ..."``.
    """

    name: str
    label: str
    chain: str | None
    passed: int
    total: int
    median_value: float | None = None
    unit: str = ""
    examples: list[str] = field(default_factory=list)

    @property
    def failed(self) -> int:
        return max(0, self.total - self.passed)

    @property
    def is_passing(self) -> bool:
        return self.failed == 0

    @property
    def pass_rate(self) -> float:
        return 1.0 if self.total == 0 else self.passed / self.total


@dataclass
class AssemblyQCReport:
    """Structured report covering every check ``assemble_qc_report``
    runs. Replaces the earlier string-only return (#67)."""

    n_rows: int
    checks: list[AssemblyQCCheck] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(c.is_passing for c in self.checks)

    def format_text(self) -> str:
        """Multi-line ``[assemble QC]`` summary suitable for logging."""
        lines = [f"[assemble QC] {self.n_rows} / {self.n_rows} rows checked"]
        for c in self.checks:
            mark = "✓" if c.is_passing else "✗"
            extras = []
            if c.median_value is not None:
                extras.append(
                    f"median {c.median_value:g}"
                    + (f" {c.unit}" if c.unit else "")
                )
            extra = f" ({', '.join(extras)})" if extras else ""
            lines.append(
                f"[assemble QC]   {mark} {c.label}: {c.passed}/{c.total} pass{extra}"
            )
        lines.append(
            f"[assemble QC] {'PASS' if self.passed else 'FAIL — see warnings above'}"
        )
        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Tabular form for programmatic consumption."""
        return pd.DataFrame(
            [
                {
                    "name": c.name,
                    "label": c.label,
                    "chain": c.chain,
                    "passed": c.passed,
                    "failed": c.failed,
                    "total": c.total,
                    "pass_rate": c.pass_rate,
                    "median": c.median_value,
                    "unit": c.unit,
                }
                for c in self.checks
            ]
        )

    def __str__(self) -> str:
        return self.format_text()


def build_assembly_qc_report(df: pd.DataFrame) -> AssemblyQCReport:
    """Build a structured :class:`AssemblyQCReport` (#67).

    Each per-check tally — terminal residue, length floors, premature
    stop, standard alphabet, β J→C parity, single-chain linker — is
    one :class:`AssemblyQCCheck` in the report's ``checks`` list.
    """
    n = len(df)
    report = AssemblyQCReport(n_rows=n)

    def _add(
        name: str,
        label: str,
        chain: str | None,
        passing: pd.Series,
        total: int,
        *,
        median_value: float | None = None,
        unit: str = "",
        df_for_examples: pd.DataFrame | None = None,
        failure_formatter=None,
    ) -> None:
        passed = int(passing.sum())
        examples: list[str] = []
        if failure_formatter is not None and df_for_examples is not None:
            failing_idx = passing.index[~passing]
            for idx in failing_idx[:3]:
                msg = failure_formatter(df_for_examples.loc[idx], idx)
                if msg:
                    examples.append(msg)
        report.checks.append(
            AssemblyQCCheck(
                name=name,
                label=label,
                chain=chain,
                passed=passed,
                total=total,
                median_value=median_value,
                unit=unit,
                examples=examples,
            )
        )

    for chain in ("alpha", "beta"):
        full_col = f"full_{chain}_aa"
        if full_col not in df.columns:
            continue
        full = df[full_col]
        present_mask = full.apply(lambda s: isinstance(s, str) and bool(s))
        total = int(present_mask.sum())
        if total == 0:
            continue
        sub = df[present_mask]
        sub_full = full[present_mask]
        chain_label = "α" if chain == "alpha" else "β"

        # Terminal residue (against canonical for the row's c_gene).
        def _ends_canonical(row, _chain=chain):
            s = row.get(f"full_{_chain}_aa", "")
            if not isinstance(s, str) or not s:
                return False
            c = _resolve_c_gene(row, _chain)
            base = c.split("*")[0]
            if base in CONSTANT_REGION_ENDINGS:
                return s.endswith(CONSTANT_REGION_ENDINGS[base])
            return any(
                s.endswith(end) for end in CONSTANT_REGION_ENDINGS.values()
            )

        end_pass = sub.apply(_ends_canonical, axis=1)
        _add(
            f"{chain}_terminal_residue",
            f"{chain_label} terminal residue",
            chain,
            end_pass,
            total,
            df_for_examples=sub,
            failure_formatter=lambda row, idx, c=chain: (
                f"clone {idx}: tail {row[f'full_{c}_aa'][-10:]!r}"
            ),
        )

        # Length floor (200 aa).
        lengths = sub_full.str.len()
        len_pass = lengths >= 200
        _add(
            f"{chain}_length_floor",
            f"{chain_label} length ≥ 200",
            chain,
            len_pass,
            total,
            median_value=float(lengths.median()),
            unit="aa",
        )

        # Constant length floor.
        const_col = f"{chain}_constant_aa"
        if const_col in df.columns:
            const_lens = sub[const_col].apply(
                lambda s: len(s) if isinstance(s, str) else 0
            )
            floor = CONSTANT_AA_FLOOR[chain]
            const_pass = const_lens >= floor
            present_const_lens = const_lens[const_lens > 0]
            cmedian = (
                float(present_const_lens.median())
                if len(present_const_lens)
                else 0.0
            )
            _add(
                f"{chain}_constant_floor",
                f"{chain_label} constant ≥ {floor}",
                chain,
                const_pass,
                total,
                median_value=cmedian,
                unit="aa",
            )

        # No premature stop.
        def _no_premature_stop(s):
            if not isinstance(s, str) or not s:
                return True
            return "*" not in s.rstrip("*")

        stop_pass = sub_full.apply(_no_premature_stop)
        _add(
            f"{chain}_no_premature_stop",
            f"{chain_label} no premature stop",
            chain,
            stop_pass,
            total,
        )

        # Standard alphabet.
        def _valid_alphabet(s):
            if not isinstance(s, str) or not s:
                return True
            return not (set(s) - _VALID_AA_CHARS)

        alpha_pass = sub_full.apply(_valid_alphabet)
        _add(
            f"{chain}_standard_alphabet",
            f"{chain_label} standard alphabet",
            chain,
            alpha_pass,
            total,
        )

    # β J→C parity.
    if "beta_j_gene" in df.columns:
        def _jc_parity_ok(row):
            j = row.get("beta_j_gene", "")
            c = _resolve_c_gene(row, "beta")
            fam = _beta_j_family(j)
            cbase = _beta_c_base(c)
            expected = BETA_JC_PARITY.get(fam) if fam else None
            if not expected or not cbase:
                return True
            return cbase == expected

        jc_pass = df.apply(_jc_parity_ok, axis=1)
        _add("beta_jc_parity", "β J→C in-cis parity", "beta", jc_pass, n)

    # Single-chain linker integrity.
    if "single_chain_aa" in df.columns and "linker" in df.columns:
        def _linker_ok(row):
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            if not isinstance(sc, str) or not sc:
                return True
            if not isinstance(linker, str) or not linker:
                return True
            return sc.count(linker) == 1

        link_pass = df.apply(_linker_ok, axis=1)
        _add(
            "single_chain_linker_present",
            "single_chain linker exactly once",
            "single_chain",
            link_pass,
            n,
        )

        # Byte-for-byte: single_chain == β.rstrip('*') + linker + α.
        def _single_chain_byte_exact(row):
            sc = row.get("single_chain_aa")
            linker = row.get("linker")
            beta = row.get("full_beta_aa")
            alpha = row.get("full_alpha_aa")
            if not all(
                isinstance(x, str) and x for x in (sc, linker, beta, alpha)
            ):
                return True
            return sc == beta.rstrip("*") + linker + alpha

        sc_byte_pass = df.apply(_single_chain_byte_exact, axis=1)
        _add(
            "single_chain_byte_exact",
            "single_chain == β + linker + α",
            "single_chain",
            sc_byte_pass,
            n,
        )

        # Canonical 2A peptide AA (when the linker name matches a known
        # 2A peptide, the AA should be the canonical sequence — guards
        # against custom-substituted linkers with the right NAME but
        # wrong AA).
        def _linker_canonical(row):
            linker = row.get("linker")
            if not isinstance(linker, str) or not linker:
                return True
            # The DataFrame's ``linker`` column stores the AA, not the
            # peptide name. Verify it matches one of the canonical 2A
            # AAs (or is a non-2A custom linker — which we accept).
            canonical_aas = {info["aa"] for info in LINKERS.values()}
            non_2a = linker not in canonical_aas and not any(
                linker.endswith(a[-6:]) for a in canonical_aas
            )
            return non_2a or linker in canonical_aas

        canon_pass = df.apply(_linker_canonical, axis=1)
        _add(
            "single_chain_2a_canonical",
            "2A peptide canonical AA",
            "single_chain",
            canon_pass,
            n,
        )

    return report


def assemble_qc_report(df: pd.DataFrame) -> str:
    """Aggregate QC summary across all rows (#67).

    Returns multi-line ``[assemble QC]`` text suitable for logging at
    the end of ``assemble_full_sequences``. Reports per-check tallies
    (pass / fail + median where relevant) plus a final PASS / FAIL
    banner. For programmatic consumption use
    :func:`build_assembly_qc_report`, which returns the same data as
    an :class:`AssemblyQCReport` object.
    """
    return build_assembly_qc_report(df).format_text()


def _parse_divergence_positions_string(s: str) -> list[tuple[int, str, str]]:
    """Inverse of :func:`_format_divergence_positions`."""
    if not isinstance(s, str) or not s:
        return []
    out: list[tuple[int, str, str]] = []
    for chunk in s.split(";"):
        if ":" not in chunk or "->" not in chunk:
            continue
        pos_str, change = chunk.split(":", 1)
        exp, obs = change.split("->", 1)
        try:
            out.append((int(pos_str), exp, obs))
        except ValueError:
            continue
    return out


_NOVEL_ALLELE_EMPTY_COLUMNS: tuple[str, ...] = (
    "chain", "gene", "reference_allele", "candidate_allele_label",
    "variant_description", "position", "expected_aa", "observed_aa",
    "n_clones", "n_observed_at_position", "pct_observed_at_position",
    "pct_chain_clones", "n_v_genes", "n_samples", "verdict",
)


def _expand_samples_unique(samples_series: pd.Series) -> int:
    """Count distinct semicolon-joined sample names in a Series."""
    sample_set: set[str] = set()
    for s in samples_series.dropna().astype(str):
        for piece in s.split(";"):
            piece = piece.strip()
            if piece:
                sample_set.add(piece)
    return len(sample_set)


def _candidate_constant_alleles_for_chain(chain: str) -> list[tuple[str, str, str]]:
    """Return ``(gene, allele_label, aa)`` references compatible with chain."""
    if chain == "alpha":
        genes = ("TRAC",)
    elif chain == "beta":
        genes = ("TRBC1", "TRBC2")
    else:
        genes = tuple(HUMAN_CONSTANT_ALLELES)

    refs: list[tuple[str, str, str]] = []
    for gene in genes:
        for allele, aa in HUMAN_CONSTANT_ALLELES.get(gene, {}).items():
            refs.append((gene, f"{gene}*{allele}", aa))
    return refs


def _default_reference_for_gene(
    canonical_name: object,
) -> tuple[str | None, str | None]:
    """Return the default packaged allele label + AA for a known gene.

    The "default" is the first allele in FASTA order, which is ``*01``
    for every packaged gene and whose AA prefix matches the legacy
    ``HUMAN_CONSTANT_REGIONS_AA`` baseline — so the recompute path's
    ``expected_aa`` is unchanged from before this helper existed.
    """
    gene = str(canonical_name)
    pool = HUMAN_CONSTANT_ALLELES.get(gene, {})
    if pool:
        allele, aa = next(iter(pool.items()))
        return f"{gene}*{allele}", aa
    aa = HUMAN_CONSTANT_REGIONS_AA.get(gene)
    if aa:
        return gene, aa
    return None, None


def _closest_reference_for_observations(
    chain: str,
    observations: pd.Series,
) -> tuple[str | None, str | None]:
    """Pick the nearest packaged constant allele for a novel gene label.

    The novel-allele audit still needs a reference sequence to define
    ``expected_aa``. For labels outside the packaged constants, choose
    the compatible chain reference with the highest aggregate prefix
    agreement across observed contig AA strings.

    O(n_obs × n_alleles × prefix_len), but only reached on the
    novel-gene path — known genes resolve via
    :func:`_default_reference_for_gene` before this is called.
    """
    obs_values = [
        str(obs)
        for obs in observations.dropna()
        if isinstance(obs, str) and obs
    ]
    if not obs_values:
        return None, None

    best_label: str | None = None
    best_aa: str | None = None
    best_key: tuple[float, int, int] = (-1.0, -1, -1)
    for _, label, allele_aa in _candidate_constant_alleles_for_chain(chain):
        n_agree_total = 0
        n_compared_total = 0
        for obs in obs_values:
            n_agree, n_compared = _score_allele_against_contig(obs, allele_aa)
            n_agree_total += n_agree
            n_compared_total += n_compared
        if n_compared_total == 0:
            continue
        score = n_agree_total / n_compared_total
        key = (score, n_agree_total, n_compared_total)
        if key > best_key:
            best_key = key
            best_label = label
            best_aa = allele_aa
    return best_label, best_aa


def _reference_for_gene_or_observations(
    chain: str,
    canonical_name: object,
    observations: pd.Series,
) -> tuple[str | None, str | None]:
    """Resolve the reference allele used to describe divergences."""
    reference_label, reference_aa = _default_reference_for_gene(canonical_name)
    if reference_aa:
        return reference_label, reference_aa
    return _closest_reference_for_observations(chain, observations)


def _variant_description(position: int, expected_aa: str, observed_aa: str) -> str:
    """Compact protein-diff descriptor, e.g. ``p3N>K``."""
    return f"p{int(position)}{expected_aa}>{observed_aa}"


def _candidate_allele_label(
    gene: object,
    reference_allele: object,
    variant_description: object,
) -> str:
    """Human-readable candidate label anchored to a packaged reference."""
    gene_s = str(gene)
    ref_s = str(reference_allele)
    desc_s = str(variant_description)
    ref_gene = ref_s.split("*", 1)[0]
    if gene_s == ref_gene:
        return f"{ref_s}:{desc_s}"
    return f"{gene_s}~{ref_s}:{desc_s}"


def _count_observed_positions_for_chain(
    df: pd.DataFrame,
    chain: str,
) -> pd.DataFrame:
    """Count clones with contig coverage at each observed C-region position.

    #119's denominator is "all contigs observed at this position", not
    only no-call clones and not every clone in the cohort. Coverage can
    drop with depth into the constant region, so each position gets its
    own denominator.
    """
    obs_col = f"{chain}_observed_constant_aa_start"
    gene_col = f"{chain}_c_gene_canonical"
    columns = ["chain", "gene", "position", "n_observed_at_position"]
    if obs_col not in df.columns or gene_col not in df.columns:
        return pd.DataFrame(columns=columns)

    observed = df.loc[df[obs_col].notna(), [gene_col, obs_col]].copy()
    if observed.empty:
        return pd.DataFrame(columns=columns)

    observed[obs_col] = observed[obs_col].astype(str)
    observed = observed[observed[obs_col].str.len() > 0]
    if observed.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    # ``observed=True`` skips unused categorical levels (AnnData/H5AD
    # round-trips leave them behind), so no empty groups reach the body.
    for canonical_name, gene_slice in observed.groupby(
        gene_col, dropna=True, observed=True,
    ):
        obs_series = gene_slice[obs_col]
        max_pos = min(int(obs_series.str.len().max()), _MAX_DIVERGENCE_SCAN_POSITIONS)
        for i in range(max_pos):
            n_observed = int(obs_series.str.get(i).notna().sum())
            if n_observed:
                rows.append({
                    "chain": chain,
                    "gene": str(canonical_name),
                    "position": i + 1,
                    "n_observed_at_position": n_observed,
                })
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def _collect_divergences_for_chain(
    df: pd.DataFrame,
    chain: str,
) -> pd.DataFrame:
    """Return a long-form DataFrame of divergence rows for ``chain``.

    Two sources are unioned:

    * **Stored** — clones with ``_allele_called`` set AND non-empty
      ``_allele_divergence_positions``. Parsed via
      :func:`_parse_divergence_positions_string`. Captures
      ``auto_detected`` clones whose contig diverges from the called
      allele at tolerated positions — the heterozygous-donor case
      from #120.
    * **Recomputed** — clones with no-call reasons (no allele was
      committed, so the stored divergence column is ``None``).
      Compare ``_observed_constant_aa_start`` vs the gene's default
      canonical AA position-by-position. If the gene label is novel,
      compare against the nearest compatible packaged allele.
    """
    reason_col = f"{chain}_allele_called_reason"
    if reason_col not in df.columns:
        return pd.DataFrame()

    called_col = f"{chain}_allele_called"
    divs_col = f"{chain}_allele_divergence_positions"
    obs_col = f"{chain}_observed_constant_aa_start"
    gene_col = f"{chain}_c_gene_canonical"
    v_col = f"{chain}_v_gene"

    if gene_col not in df.columns:
        return pd.DataFrame()

    collected_frames: list[pd.DataFrame] = []

    # ---- Source 1: stored divergences from auto_detected clones ----
    if divs_col in df.columns and called_col in df.columns:
        stored = df.loc[
            df[divs_col].notna() & df[called_col].notna(),
            [called_col, divs_col, gene_col, v_col, "samples"] if v_col in df.columns
            else [called_col, divs_col, gene_col, "samples"],
        ].copy()
        if not stored.empty:
            # Parse each clone's divergence string into a list, then
            # explode so each (clone, position) becomes its own row.
            stored["__parsed"] = stored[divs_col].apply(
                _parse_divergence_positions_string
            )
            stored = stored[stored["__parsed"].str.len() > 0]
            if not stored.empty:
                stored = stored.explode("__parsed", ignore_index=False)
                stored = stored[stored["__parsed"].notna()]
                stored["position"] = stored["__parsed"].str[0].astype(int)
                stored["expected_aa"] = stored["__parsed"].str[1]
                stored["observed_aa"] = stored["__parsed"].str[2]
                stored = stored.drop(columns=["__parsed", divs_col])
                stored = stored.rename(columns={gene_col: "gene"})
                stored["gene"] = stored["gene"].astype(str)
                stored["reference_allele"] = stored[called_col].astype(str)
                if v_col in stored.columns:
                    stored = stored.rename(columns={v_col: "v_gene"})
                else:
                    stored["v_gene"] = None
                stored = stored.rename(columns={"samples": "sample"})
                stored["clone_idx"] = stored.index
                collected_frames.append(stored[[
                    "clone_idx", "gene", "reference_allele", "position",
                    "expected_aa", "observed_aa", "v_gene", "sample",
                ]])

    # ---- Source 2: recomputed for no-call clones ----
    # Compares against ``HUMAN_CONSTANT_REGIONS_AA`` (the *default*
    # canonical), not the user-overridden allele. That's safe here
    # because user-override clones get ``reason=overridden``, which
    # is not in ``no_call_reasons`` — so we never recompute for
    # them. The stored-source path above handles overridden clones
    # via their own ``_allele_divergence_positions`` column.
    no_call_reasons = {
        ALLELE_REASON_DIVERGENT_CONTIG,
        ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
        ALLELE_REASON_SPARSE_CONTIG,
    }
    if obs_col in df.columns:
        recompute_mask = (
            df[reason_col].isin(no_call_reasons)
            & df[obs_col].notna()
        )
        # Process per-gene so we have one canonical AA string to
        # compare against; vectorized within each gene group.
        # ``observed=True`` skips unused categorical levels, so no empty
        # groups reach the body (AnnData/H5AD round-trip robustness).
        for canonical_name, gene_slice in df.loc[
            recompute_mask
        ].groupby(gene_col, dropna=True, observed=True):
            reference_label, reference_aa = _reference_for_gene_or_observations(
                chain, canonical_name, gene_slice[obs_col],
            )
            if not reference_aa or not reference_label:
                continue
            obs_series = gene_slice[obs_col]
            v_series = gene_slice[v_col] if v_col in gene_slice.columns else pd.Series(
                [None] * len(gene_slice), index=gene_slice.index,
            )
            sample_series = gene_slice["samples"] if "samples" in gene_slice.columns else pd.Series(
                [None] * len(gene_slice), index=gene_slice.index,
            )
            max_pos = min(len(reference_aa), _MAX_DIVERGENCE_SCAN_POSITIONS)
            per_position: list[pd.DataFrame] = []
            for i in range(max_pos):
                canonical_residue = reference_aa[i]
                # str.get returns None when index is past the string length.
                obs_residue = obs_series.str.get(i)
                diverged_mask = obs_residue.notna() & (obs_residue != canonical_residue)
                if not diverged_mask.any():
                    continue
                per_position.append(pd.DataFrame({
                    "clone_idx": gene_slice.index[diverged_mask],
                    "gene": str(canonical_name),
                    "reference_allele": reference_label,
                    "position": i + 1,
                    "expected_aa": canonical_residue,
                    "observed_aa": obs_residue[diverged_mask].values,
                    "v_gene": v_series[diverged_mask].values,
                    "sample": sample_series[diverged_mask].values,
                }))
            if per_position:
                collected_frames.append(pd.concat(per_position, ignore_index=True))

    if not collected_frames:
        return pd.DataFrame()
    out = pd.concat(collected_frames, ignore_index=True)
    out.insert(0, "chain", chain)
    return out


def detect_novel_alleles(
    df: pd.DataFrame,
    *,
    min_pct: float = 0.05,
    min_v_spread: int = 3,
    min_samples: int = 2,
    min_cohort_size: int = 20,
) -> pd.DataFrame:
    """Aggregate per-clone divergence positions into cohort-level
    novel-allele candidates (#119).

    For each ``(chain, position, expected_aa, observed_aa)`` tuple
    across all clones whose contig diverges from canonical (either
    stored ``_allele_divergence_positions`` on auto-called clones or
    recomputed against canonical for no-call clones), count:

    * ``n_clones`` — number of clones carrying the variant.
    * ``pct_chain_clones`` — fraction of in-cohort clones for that
      chain showing the variant.
    * ``n_observed_at_position`` — number of clones with contig AA
      coverage at that chain/gene/position.
    * ``pct_observed_at_position`` — fraction of observed contigs at
      that position showing the variant. This is the candidate-calling
      denominator.
    * ``n_v_genes`` — number of distinct V genes the variant
      appears in.
    * ``n_samples`` — number of distinct CellRanger samples the
      variant appears in.
    * ``verdict`` — ``"novel_allele_candidate"`` when
      ``pct_observed_at_position >= min_pct`` AND
      ``n_v_genes >= min_v_spread`` AND ``n_samples >= min_samples``;
      ``"likely_artifact"`` otherwise.

    Returns a long-form DataFrame sorted by ``n_clones`` descending.
    Empty DataFrame when ``df`` has fewer than ``min_cohort_size``
    clones with the chain populated (the heuristic is unreliable
    below that floor) or when no divergences are present.

    Vectorized: ~100x faster than the original ``iterrows`` loop on
    cohorts of 10k+ clones (no per-row Python). For 100k clones the
    runtime is dominated by pandas groupby — still single-threaded
    but fully numpy under the hood.

    The heuristic encodes simple population genetics: a real
    polymorphism distributes across V-genes and samples (since V-gene
    usage is independent of constant-region allele); an artifact
    concentrates in one V or one sample (PCR / assembly noise tied to
    a specific V sequence). See #119 for the full rationale and pilot
    validation.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`assemble_full_sequences`.
    min_pct : float, default 0.05
        Minimum fraction of contigs observed at that position for the
        variant to be called a candidate. 5% is conservative for
        heterozygous donors (~50% expected) and dominates over
        single-V artifacts (~1%) but still surfaces minor population
        alleles.
    min_v_spread : int, default 3
        Minimum distinct V genes the variant must appear in. V-gene
        usage is independent of constant allele, so true polymorphisms
        spread across V's; V-concentrated variants are almost always
        contig-assembly artifacts.
    min_samples : int, default 2
        Minimum distinct samples. Single-sample variants may reflect
        a batch effect rather than a real allele.
    min_cohort_size : int, default 20
        Skip the aggregator entirely when neither chain has at least
        this many clones with ``_allele_called_reason`` populated —
        ``min_v_spread`` can't exceed the cohort size, so smaller
        cohorts produce noise verdicts. Set to 0 to force evaluation
        on any cohort.
    """
    empty = pd.DataFrame(columns=list(_NOVEL_ALLELE_EMPTY_COLUMNS))

    # Compute chain totals up-front and apply the min_cohort_size gate.
    chain_totals: dict[str, int] = {}
    for chain in ("alpha", "beta"):
        reason_col = f"{chain}_allele_called_reason"
        if reason_col in df.columns:
            chain_totals[chain] = int(df[reason_col].notna().sum())
    if (
        min_cohort_size > 0
        and max(chain_totals.values(), default=0) < min_cohort_size
    ):
        return empty

    per_chain: list[pd.DataFrame] = []
    observed_position_counts: list[pd.DataFrame] = []
    for chain in ("alpha", "beta"):
        counts = _count_observed_positions_for_chain(df, chain)
        if not counts.empty:
            observed_position_counts.append(counts)
        chain_df = _collect_divergences_for_chain(df, chain)
        if not chain_df.empty:
            per_chain.append(chain_df)
    if not per_chain:
        return empty

    long_df = pd.concat(per_chain, ignore_index=True)

    # Vectorized aggregation. ``unique`` -> nunique on clone_idx and
    # v_gene; samples need to be expanded via the helper.
    grouped = long_df.groupby(
        [
            "chain", "gene", "reference_allele", "position",
            "expected_aa", "observed_aa",
        ],
        dropna=False, sort=False, observed=True,
    )
    n_clones = grouped["clone_idx"].nunique().rename("n_clones")
    n_v_genes = grouped["v_gene"].nunique(dropna=True).rename("n_v_genes")
    n_samples = grouped["sample"].agg(_expand_samples_unique).rename("n_samples")
    summary = pd.concat([n_clones, n_v_genes, n_samples], axis=1).reset_index()

    if observed_position_counts:
        position_counts = pd.concat(observed_position_counts, ignore_index=True)
        summary = summary.merge(
            position_counts,
            on=["chain", "gene", "position"],
            how="left",
        )
    else:
        summary["n_observed_at_position"] = pd.NA

    # Vectorized per-row pct: map each row's chain to its total, then
    # clip the denominator to ≥1 so zero-total chains (impossible
    # given the gate above but cheap to guard) don't produce inf.
    chain_total_series = summary["chain"].map(chain_totals).fillna(1).clip(lower=1)
    summary["pct_chain_clones"] = (
        summary["n_clones"] / chain_total_series
    ).round(4)
    # Prefer the position-specific observed-contig denominator. For
    # legacy CSVs that have divergence columns but not observed AA
    # columns, fall back to the old chain-level denominator. The
    # observed-position count and the variant count come from different
    # columns (``_observed_constant_aa_start`` vs the parsed divergence
    # string), so clamp the denominator to at least ``n_clones`` — a
    # variant clone is by definition observed at its own position, and
    # this guarantees ``pct_observed_at_position`` stays in [0, 1] even
    # if the two columns disagree.
    summary["n_observed_at_position"] = (
        summary["n_observed_at_position"]
        .fillna(chain_total_series)
        .astype(int)
    )
    summary["n_observed_at_position"] = summary[
        ["n_observed_at_position", "n_clones"]
    ].max(axis=1).clip(lower=1)
    summary["pct_observed_at_position"] = (
        summary["n_clones"] / summary["n_observed_at_position"]
    ).round(4)
    summary["variant_description"] = [
        _variant_description(pos, exp, obs)
        for pos, exp, obs in zip(
            summary["position"], summary["expected_aa"], summary["observed_aa"],
        )
    ]
    summary["candidate_allele_label"] = [
        _candidate_allele_label(gene, ref, desc)
        for gene, ref, desc in zip(
            summary["gene"],
            summary["reference_allele"],
            summary["variant_description"],
        )
    ]
    summary["verdict"] = (
        (summary["pct_observed_at_position"] >= min_pct)
        & (summary["n_v_genes"] >= min_v_spread)
        & (summary["n_samples"] >= min_samples)
    ).map(
        {True: "novel_allele_candidate", False: "likely_artifact"}
    )
    summary["position"] = summary["position"].astype(int)
    summary = summary[list(_NOVEL_ALLELE_EMPTY_COLUMNS)].sort_values(
        ["n_clones", "pct_observed_at_position"], ascending=[False, False],
    ).reset_index(drop=True)
    return summary


def allele_audit_report(
    df: pd.DataFrame,
    *,
    min_pct: float = 0.05,
    min_v_spread: int = 3,
    min_samples: int = 2,
    min_cohort_size: int = 20,
) -> str:
    """Human-readable cohort-level allele audit (#120 ask 3).

    Combines:
    * per-chain confident allele tally (calls vs. no-calls broken out
      by ``allele_called_reason``)
    * cohort-level novel-allele candidates from
      :func:`detect_novel_alleles`

    ``min_cohort_size`` (default 20) is forwarded to
    :func:`detect_novel_alleles` — smaller cohorts skip the
    aggregator entirely (it adds noticeable latency on every
    verbose assembly and produces noise verdicts below the
    ``min_v_spread`` floor).

    Output is suitable for printing at the end of a ``tcrsift run``
    pipeline or as the body of ``tcrsift audit-alleles``.
    """
    lines: list[str] = ["[allele audit]"]
    for chain in ("alpha", "beta"):
        called_col = f"{chain}_allele_called"
        reason_col = f"{chain}_allele_called_reason"
        if called_col not in df.columns or reason_col not in df.columns:
            continue
        chain_mask = df[reason_col].notna()
        n_total = int(chain_mask.sum())
        if n_total == 0:
            continue
        n_called = int(df[called_col].notna().sum())
        called_counts = df.loc[
            df[called_col].notna(), called_col
        ].value_counts()
        reason_counts = df[reason_col].value_counts(dropna=False)
        lines.append(
            f"  {chain} chain: {n_called}/{n_total} called "
            f"({100 * n_called / max(n_total, 1):.1f}%)"
        )
        for allele, n in called_counts.items():
            lines.append(f"    {allele}: {n}")
        no_calls = n_total - n_called
        if no_calls > 0:
            lines.append(f"    No-call breakdown ({no_calls} clones):")
            for reason, n in reason_counts.items():
                if reason in {ALLELE_REASON_AUTO_DETECTED, ALLELE_REASON_OVERRIDDEN}:
                    continue
                if pd.isna(reason):
                    continue
                lines.append(f"      {reason}: {n}")

    novel = detect_novel_alleles(
        df,
        min_pct=min_pct,
        min_v_spread=min_v_spread,
        min_samples=min_samples,
        min_cohort_size=min_cohort_size,
    )
    if not novel.empty:
        lines.append("")
        lines.append("  Observed polymorphisms (contig vs. canonical):")
        for _, row in novel.iterrows():
            verdict_tag = (
                "← novel_allele_candidate"
                if row["verdict"] == "novel_allele_candidate"
                else "→ likely_artifact"
            )
            lines.append(
                f"    {row['chain']} {row['candidate_allele_label']}  "
                f"{row['expected_aa']}→{row['observed_aa']}  "
                f"n={row['n_clones']}/{row['n_observed_at_position']} "
                f"observed ({100 * row['pct_observed_at_position']:.1f}%)  "
                f"V-genes={row['n_v_genes']}  samples={row['n_samples']}  "
                f"{verdict_tag}"
            )
        candidates = novel[novel["verdict"] == "novel_allele_candidate"]
        if candidates.empty:
            lines.append("  No novel-allele candidates detected.")
    else:
        lines.append("")
        lines.append("  No contig-vs-canonical divergences observed.")
    return "\n".join(lines)


def _nt_translates_to(nt: str, aa: str, strip_trailing_stop: bool = True) -> bool:
    """Check that ``translate(nt)`` matches ``aa``, optionally dropping
    *all* trailing stop codons from ``nt`` first.

    Used by the NT→AA round-trip invariants (#91). Returns False on
    any of: NT length not a multiple of 3 (after optional stop trim),
    translation contains a stop codon mid-chain, translated AA differs
    from the assembled AA. The caller composes the failure message
    so it can include which column failed.

    The ``strip_trailing_stop`` flag strips *every* trailing stop, not
    just one. This was a one-stop strip pre-#116, which broke when the
    constant CDS picked up its second non-redundant stop (e.g. the
    default ``"TAATGA"`` suffix).
    """
    if not isinstance(nt, str) or not nt:
        return False
    if not isinstance(aa, str) or not aa:
        return False
    work = nt
    if strip_trailing_stop:
        while len(work) >= 3 and work[-3:] in {"TAA", "TAG", "TGA"}:
            work = work[:-3]
    if len(work) != 3 * len(aa):
        return False
    translated, ragged = translate_dna(work)
    if ragged:
        return False
    # ``translate_dna`` stops at the first internal stop, returning a
    # shorter aa string than expected — equality with ``aa`` (which
    # has no stop) already catches that case, but assert explicitly
    # to make intent clear.
    return translated == aa


def _nt_contig_in_frame_no_premature_stop(nt: str, aa: str) -> bool:
    """Check that a ``_contig`` NT column (donor's CellRanger contig
    bytes past the J→C junction) is in-frame and free of premature
    stops.

    We DON'T require AA agreement with ``aa`` — the contig view is
    raw donor bytes; polymorphisms vs. the canonical reference are
    expected (that's what the blender handles by switching at the
    mismatch). We DO require the largest 3·N prefix to:

    * translate cleanly (no ragged byte that ``translate_dna`` flags),
    * contain no stop codons in the first ``len(aa)`` positions
      (only the 3' UTR past the coding region may carry stops, but
      capping at ``len(aa)`` ignores anything past the coding region).

    Returns False if the contig NT has zero in-frame codons before
    the first stop, or if it appears non-coding.
    """
    if not isinstance(nt, str) or not nt:
        return False
    if not isinstance(aa, str) or not aa:
        return False
    n_codons = min(len(nt) // 3, len(aa))
    if n_codons == 0:
        return False
    work = nt[: 3 * n_codons]
    translated, ragged = translate_dna(work)
    if ragged:
        return False
    if "*" in translated:
        return False
    return True


def _validate_nt_aa_roundtrip(df: pd.DataFrame, _lb) -> None:
    """Run NT→AA back-translation checks on every NT column the
    assembler emits, appending load-bearing messages via ``_lb`` for
    mismatches.

    Checks (per chain α/β where columns are present):

    - ``vdj_{chain}_nt`` translates to ``vdj_{chain}_aa`` (no stop;
      length must be exactly 3×AA).
    - ``{chain}_constant_nt`` translates to ``{chain}_constant_aa``
      after dropping any trailing stop codons.
    - ``{chain}_constant_nt_optimized`` translates to
      ``{chain}_constant_aa`` after dropping trailing stops (#116).
    - ``{chain}_constant_nt_contig`` (when not None) is in-frame
      and free of premature stops — donor polymorphisms vs. canonical
      AA are NOT flagged (the blender handles them); only frame
      errors or non-coding contigs surface here (#116).
    - ``full_{chain}_nt`` translates to ``full_{chain}_aa`` after
      dropping the trailing stop codon. **This is the integration
      check that would have caught #91 (+1 nt overshoot at the
      VDJ→C boundary) on day one.**
    - ``full_{chain}_nt_optimized`` translates to ``full_{chain}_aa``
      after dropping trailing stops (#116).
    - ``single_chain_nt`` translates to ``single_chain_aa`` (trailing
      stop dropped).
    """
    for chain in ("alpha", "beta"):
        nt_col = f"vdj_{chain}_nt"
        aa_col = f"vdj_{chain}_aa"
        if nt_col in df.columns and aa_col in df.columns:
            for idx, row in df.iterrows():
                nt = row.get(nt_col)
                aa = row.get(aa_col)
                if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                    if not _nt_translates_to(nt, aa, strip_trailing_stop=False):
                        _lb(idx,
                            f"{nt_col} does not translate to {aa_col} "
                            f"(NT length {len(nt)}, expected {3 * len(aa)})")

        # Both `_constant_nt` (the legacy blend column) and
        # `_constant_nt_optimized` should translate exactly to
        # `{chain}_constant_aa` after dropping all trailing stops.
        aa_col = f"{chain}_constant_aa"
        for nt_col in (f"{chain}_constant_nt", f"{chain}_constant_nt_optimized"):
            if nt_col in df.columns and aa_col in df.columns:
                for idx, row in df.iterrows():
                    nt = row.get(nt_col)
                    aa = row.get(aa_col)
                    if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                        if not _nt_translates_to(nt, aa, strip_trailing_stop=True):
                            _lb(idx,
                                f"{nt_col} does not translate to {aa_col} "
                                f"(after stripping trailing stops)")

        # `_constant_nt_contig` (#116): pure donor NT. May end mid-
        # codon or extend past the AA — translate the largest 3·N
        # prefix and check it matches the AA prefix of equal length.
        nt_col = f"{chain}_constant_nt_contig"
        if nt_col in df.columns and aa_col in df.columns:
            for idx, row in df.iterrows():
                nt = row.get(nt_col)
                aa = row.get(aa_col)
                if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                    if not _nt_contig_in_frame_no_premature_stop(nt, aa):
                        _lb(idx,
                            f"{nt_col} is out-of-frame or hits a premature "
                            f"stop within the first {3 * len(aa)} nt — the "
                            f"donor's contig past the J→C boundary appears "
                            f"non-coding")

        for nt_col in (f"full_{chain}_nt", f"full_{chain}_nt_optimized"):
            aa_col = f"full_{chain}_aa"
            if nt_col in df.columns and aa_col in df.columns:
                for idx, row in df.iterrows():
                    nt = row.get(nt_col)
                    aa = row.get(aa_col)
                    if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                        if not _nt_translates_to(nt, aa, strip_trailing_stop=True):
                            _lb(idx,
                                f"{nt_col} does not translate to {aa_col} "
                                f"— frame likely broken at a splice boundary (#91)")

    # ``single_chain_nt`` and ``_optimized`` must translate exactly to
    # ``single_chain_aa`` (canonical AA throughout — these views use
    # canonical bytes for the C region). ``_contig`` uses donor
    # bytes for the C region and may legitimately differ from
    # canonical at polymorphisms — only require it to be in-frame
    # and free of mid-chain stops (#116).
    for nt_col in ("single_chain_nt", "single_chain_nt_optimized"):
        if nt_col in df.columns and "single_chain_aa" in df.columns:
            for idx, row in df.iterrows():
                nt = row.get(nt_col)
                aa = row.get("single_chain_aa")
                if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                    if not _nt_translates_to(nt, aa, strip_trailing_stop=True):
                        _lb(idx,
                            f"{nt_col} does not translate to single_chain_aa "
                            f"— frame likely broken at a splice boundary (#91)")
    nt_col = "single_chain_nt_contig"
    if nt_col in df.columns and "single_chain_aa" in df.columns:
        for idx, row in df.iterrows():
            nt = row.get(nt_col)
            aa = row.get("single_chain_aa")
            if isinstance(nt, str) and nt and isinstance(aa, str) and aa:
                if not _nt_contig_in_frame_no_premature_stop(nt, aa):
                    _lb(idx,
                        f"{nt_col} is out-of-frame or hits a premature "
                        f"stop within the coding region — the donor's "
                        f"contig appears non-coding past the J→C boundary")


def _expected_constant_start_from_full(
    full_aa: str, row: pd.Series, chain: str
) -> str:
    """Best-effort extract of the first ~15 aa of the constant region
    from a full assembled chain.

    Strategy: the full chain is leader + VDJ + constant. The VDJ AA
    is available as ``vdj_{chain}_aa``; split on it. If the VDJ isn't
    in the full sequence (which is itself a different failure caught
    elsewhere), fall back to slicing off the canonical-length leader.
    Returns up to 15 residues; an empty string when nothing useful
    can be extracted.
    """
    if not isinstance(full_aa, str) or not full_aa:
        return ""
    vdj = row.get(f"vdj_{chain}_aa")
    if isinstance(vdj, str) and vdj and vdj in full_aa:
        after = full_aa.split(vdj, 1)[1]
        return after[:15]
    return ""


def export_fasta(df: pd.DataFrame, output_path: str | Path, sequence_col: str = "single_chain_aa"):
    """
    Export sequences to FASTA format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with sequences
    output_path : str or Path
        Output file path
    sequence_col : str
        Column containing sequences to export
    """
    with open(output_path, "w") as f:
        for idx, row in df.iterrows():
            seq = row.get(sequence_col, "")
            if not seq:
                continue

            # Build header
            cdr3ab = row.get("CDR3ab", idx)
            cdr3a = row.get("CDR3_alpha", "")
            cdr3b = row.get("CDR3_beta", "")

            header = f">{cdr3ab} CDR3a={cdr3a} CDR3b={cdr3b}"
            f.write(f"{header}\n{seq}\n")

    logger.info(f"Exported {len(df)} sequences to {output_path}")
