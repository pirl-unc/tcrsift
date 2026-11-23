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

"""Canonical TCR constant-region alleles, J->C parity rules, and leader QC helpers."""

from __future__ import annotations

import logging
import re
from importlib.resources import files as _pkg_files

import pandas as pd

from ._assemble_codons import (
    CODON_TABLE,
    back_translate,
)
from .leaders import (
    germline_compare_leader,
    normalize_vgene,
)
from .sigpep import signal_peptide_features, sp_features_summary

logger = logging.getLogger(__name__)

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
    "mCD8B": {
        "aa": "MQPWLWLVFSMKLAALWSSSA",
        "dna": "ATGCAGCCTTGGCTGTGGCTCGTGTTCAGCATGAAGCTGGCCGCTCTGTGGAGCTCCAGCGCC",
        "source": "Murine Cd8b1 signal peptide (UniProt P10300)",
        "species": "mouse",
    },
    "mCD28": {
        "aa": "MTLRLLFLALNFFSVQVT",
        "dna": "ATGACCCTGAGACTGCTCTTCCTGGCCCTCAACTTCTTTAGCGTGCAGGTGACC",
        "source": "Murine Cd28 signal peptide (UniProt P31041)",
        "species": "mouse",
    },
    "H2Kb": {
        "aa": "MVPCTLLLLLAAALAPTQTRA",
        "dna": "ATGGTGCCCTGCACCCTGCTCCTGCTTCTGGCCGCTGCCCTCGCTCCTACCCAGACAAGAGCC",
        "source": "Murine H2-K1 (H-2Kb) signal peptide (UniProt P01901)",
        "species": "mouse",
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
        "source": "Synthetic TCR alpha-chain leader peptide (legacy key TRAC)",
        "species": "human",
    },
    "TRBC": {
        "aa": "MGTSLLCWMALCLLGADHADG",
        "dna": "ATGGGCACCAGCCTGCTGTGCTGGATGGCCCTGTGCCTGCTGGGAGCAGACCACGCCGATGGC",
        "source": "Synthetic TCR beta-chain leader peptide (legacy key TRBC)",
        "species": "human",
    },
}

_DEFAULT_LEADERS_CASEFOLD: dict[str, tuple[str, dict]] = {
    k.upper(): (k, v) for k, v in DEFAULT_LEADERS.items()
}


def _lookup_default_leader(name: str) -> tuple[str, dict] | None:
    """Case-insensitive lookup into :data:`DEFAULT_LEADERS` → ``(canonical_key, info)``."""
    return _DEFAULT_LEADERS_CASEFOLD.get(str(name).strip().upper())

# --- Signal-peptide (leader) QC + Kozak start selection (#263) -------------
#
# `from_contig` leader extraction takes the ORF prefix before the VDJ, but the
# longest-ORF rule can begin at an upstream 5'-UTR ATG that reads in-frame into
# the VDJ — over-capturing the leader (observed: a 54-aa α leader on B1-2 where
# the real signal peptide is ~19 aa). The fix is to pick the START codon by
# Kozak context strength, not ORF length, and to bound the leader to a sane
# signal-peptide length.

# Length bands for a TCR signal peptide (aa). Real TCR leaders run ~15-25 (the
# curated DEFAULT_LEADERS span 18-21), but a few V genes have genuinely long
# leaders — e.g. IMGT TRBV13*01 is 29 aa (long polar n-region). So:
#   - len > SP_LENGTH_MAX (25): TRY to recover a better nested SP at the tail
#     (Kozak reselect / h-region trim). A real over-capture (TRAV1 → 17 aa)
#     gets trimmed regardless of how long it is.
#   - if NO better nested SP exists: 26-30 aa is ACCEPTED as a genuinely-long
#     leader (``long_leader``, not an error — e.g. TRBV13); only > SP_LENGTH_HARD
#     (30) is a ``too_long`` error (override-able via leader_fallback).
SP_LENGTH_MIN = 12
SP_LENGTH_MAX = 25
SP_LENGTH_HARD = 30

# Hydrophobic residues (Kyte-Doolittle positive side) for the signal-peptide
# h-region check — a signal peptide's defining feature is a hydrophobic core.
_HYDROPHOBIC_AA = frozenset("AILMFWVC")


def _has_h_region(leader: str, *, window: int = 6, min_hydrophobic: int = 4) -> bool:
    """True if the leader contains a hydrophobic h-region.

    Looks for any ``window``-residue stretch with at least ``min_hydrophobic``
    hydrophobic residues — the hydrophobic core every signal peptide has. The
    curated DEFAULT_LEADERS all pass comfortably.
    """
    if not leader or len(leader) < window:
        return False
    flags = [1 if aa in _HYDROPHOBIC_AA else 0 for aa in leader]
    return any(
        sum(flags[i : i + window]) >= min_hydrophobic
        for i in range(len(flags) - window + 1)
    )


def leader_qc(
    leader: str | None,
    *,
    lo: int = SP_LENGTH_MIN,
    hi: int = SP_LENGTH_MAX,
    require_h_region: bool = True,
) -> str:
    """Classify a signal peptide: ``ok`` / ``missing`` / ``no_met`` /
    ``internal_stop`` / ``too_short`` / ``too_long`` / ``no_h_region`` (#263).

    Order matters: a structurally broken leader (no Met, internal stop) is
    reported before a length/shape complaint.
    """
    if not leader or not isinstance(leader, str):
        return "missing"
    if not leader.startswith("M"):
        return "no_met"
    if "*" in leader:  # a leader should never carry a stop (see translate_dna)
        return "internal_stop"
    if len(leader) < lo:
        return "too_short"
    if len(leader) > hi:
        return "too_long"
    if require_h_region and not _has_h_region(leader):
        return "no_h_region"
    return "ok"


def _classify_leader(leader_aa: str, kozak_score: int, source: str) -> str:
    """Final QC verdict for a from_contig leader given its provenance (#263).

    Applies the length banding on top of :func:`leader_qc`:

    - A germline-anchored leader (``source == "contig_germline_anchored"``) is
      ``ok`` — the germline alignment validates its start AND length, so a
      genuinely long one (TRBV13's 29 aa) is NOT flagged ``too_long`` (#267).
    - An h-region structural trim → ``hregion_trimmed`` (usable for synthesis,
      resolved on the weaker h-region basis, kept visible for review).
    - A Kozak-reselected leader → judged purely on shape (in-range by
      construction).
    - A KEPT leader (``source == "contig"``): a 26-30 aa leader with no better
      nested SP is an accepted ``long_leader`` (a genuinely long native leader
      like TRBV13's 29 aa), NOT an error; only ``> SP_LENGTH_HARD`` is
      ``too_long`` (override-able). A well-shaped in-range leader whose start has
      weak Kozak context is ``weak_kozak_start`` (leaky-scan informs, doesn't act).
    """
    if source == "contig_germline_anchored":
        return "ok"
    if source.startswith("curated_fallback"):
        return "curated_fallback"  # substituted for an implausible from_contig SP
    shape = leader_qc(leader_aa)  # ok/too_short/too_long/no_met/no_h_region/…
    if source == "curated_default":
        # A hand-verified default leader (CD28/CD8A/…) or a user-supplied one —
        # still QC'd for shape so a broken custom leader is flagged, not trusted.
        return shape
    if source == "contig_hregion_trimmed":
        return "hregion_trimmed"
    if source == "contig_kozak_reselected":
        return shape
    # source == "contig": the ORF start was kept (no better nested SP found).
    if shape == "too_long":  # len > SP_LENGTH_MAX
        return "long_leader" if len(leader_aa) <= SP_LENGTH_HARD else "too_long"
    if shape == "ok" and kozak_score < KOZAK_ADEQUATE:
        return "weak_kozak_start"
    return shape


def _record_leader(
    result: dict,
    chain: str,
    leader_aa: str,
    leader_nt: str | None,
    source: str,
    *,
    v_gene: str | None = None,
    kozak_score: int | None = None,
    support: str | None = None,
    species: str = "human",
) -> None:
    """Central sink for EVERY signal-peptide derivation path (#263/#267)."""
    result[f"{chain}_leader_aa"] = leader_aa
    if leader_nt is not None:
        result[f"{chain}_leader_nt"] = leader_nt
    result[f"{chain}_leader_len"] = len(leader_aa) if isinstance(leader_aa, str) else None
    result[f"{chain}_leader_source"] = source
    result[f"{chain}_leader_qc"] = _classify_leader(leader_aa, kozak_score or 0, source)
    if kozak_score is not None:
        result[f"{chain}_leader_kozak_score"] = kozak_score
    if support is not None:
        result[f"{chain}_leader_support"] = support
    # Universal germline comparison (always present; NaN when gene absent).
    cmp = germline_compare_leader(leader_aa, v_gene, species=species)
    g_allele = g_aa = identity = diff = None
    if cmp is not None:
        allele, g_aa, identity, diff = cmp
        g_allele = f"{normalize_vgene(v_gene)}*{allele}"
    result[f"{chain}_germline_allele"] = g_allele
    result[f"{chain}_germline_leader_aa"] = g_aa
    result[f"{chain}_germline_identity"] = identity
    result[f"{chain}_leader_vs_germline"] = diff
    # Signal-peptide feature panel (von Heijne tripartite, #270): does the leader
    # have everything expected of an SP? Drives the keep-vs-switch policy.
    feats = signal_peptide_features(leader_aa)
    result[f"{chain}_sp_features_ok"] = feats["features_ok"]
    result[f"{chain}_sp_score"] = feats["score"]
    result[f"{chain}_sp_features"] = sp_features_summary(feats)


# Kozak context scoring (Kozak 1987/1991). The optimal context is gccRccATGG
# with the A of ATG counted as +1; of the ~10 positions only two carry most of
# the functional weight — a purine (A/G) at -3 (dominant) and a G at +4. We
# score just those: a -3 purine is "adequate", +4=G adds to it. The bar for
# initiation is ADEQUATE (-3 purine), NOT a textbook-consensus match — even a
# real start need not be optimal (the validated TRAV16 start M35 is -3=G/+4=A,
# adequate but not strong). Score: 0 weak · >=2 adequate · >=4 strong.
KOZAK_ADEQUATE = 2  # threshold: a -3 purine is present


def _kozak_score(contig_nt: str, atg_pos: int) -> int:
    """Kozak initiation-context score for the ATG whose A is at ``atg_pos``.

    ``-3`` purine (A/G) contributes +2 (the dominant position; +1 more if it's
    the preferred A), ``+4`` = G contributes +2. Returns ``-1`` for incomplete
    context (ATG hard against the contig 5' end) so such starts are never
    chosen over a real downstream one.
    """
    if atg_pos < 3 or atg_pos + 3 >= len(contig_nt):
        return -1
    seq = contig_nt.upper()
    m3, p4 = seq[atg_pos - 3], seq[atg_pos + 3]
    s = 0
    if m3 in "AG":
        s += 2  # -3 purine — the dominant position
    if m3 == "A":
        s += 1  # A preferred over G at -3
    if p4 == "G":
        s += 2  # +4 = G
    return s


# Standard constant region endings for QC
CONSTANT_REGION_ENDINGS = {
    "TRAC": "LLMTLRLWSS",
    "TRBC1": "VKRKDF",
    "TRBC2": "VKRKDSRG",
}

MOUSE_CONSTANT_REGION_ENDINGS = {
    "TRAC": "LLMTLRLWSS",
    "TRBC1": "VKRKNS",
    "TRBC2": "VKKKNS",
}

CONSTANT_REGION_ENDINGS_BY_SPECIES = {
    "human": CONSTANT_REGION_ENDINGS,
    "mouse": MOUSE_CONSTANT_REGION_ENDINGS,
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
_CANONICAL_CONSTANTS_MOUSE_FASTA = "canonical_constants_mouse.fasta"
_CANONICAL_ALLELES_MOUSE_FASTA = "canonical_constants_alleles_mouse.fasta"


def _load_canonical_constants_fasta(
    filename: str = _CANONICAL_CONSTANTS_FASTA,
) -> dict[str, str]:
    """Parse a packaged canonical-constants FASTA into a ``{gene: AA}`` dict."""
    fasta_text = (
        _pkg_files("tcrsift.refseqs").joinpath(filename).read_text()
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
            f"{filename} is missing canonical entries: {sorted(missing)}"
        )
    return out


_canonical = _load_canonical_constants_fasta(_CANONICAL_CONSTANTS_FASTA)
HUMAN_TRAC_AA: str = _canonical["TRAC"]
HUMAN_TRBC1_AA: str = _canonical["TRBC1"]
HUMAN_TRBC2_AA: str = _canonical["TRBC2"]
del _canonical

_mouse_canonical = _load_canonical_constants_fasta(_CANONICAL_CONSTANTS_MOUSE_FASTA)
MOUSE_TRAC_AA: str = _mouse_canonical["TRAC"]
MOUSE_TRBC1_AA: str = _mouse_canonical["TRBC1"]
MOUSE_TRBC2_AA: str = _mouse_canonical["TRBC2"]
del _mouse_canonical

# Fallback J→C junction residue used when no contig C-region NT is available.
_FALLBACK_JUNCTION_RESIDUE: dict[str, str] = {"alpha": "N", "beta": "E"}

# Human per-TRAJ α J→C junction residue (#242).
TRAJ_JUNCTION_RESIDUES: dict[str, str] = {
    "TRAJ3": "N", "TRAJ4": "Y", "TRAJ5": "N", "TRAJ6": "Y", "TRAJ7": "N", "TRAJ8": "N",
    "TRAJ9": "N", "TRAJ10": "N", "TRAJ11": "D", "TRAJ12": "D", "TRAJ13": "N", "TRAJ15": "N",
    "TRAJ16": "N", "TRAJ17": "N", "TRAJ18": "D", "TRAJ20": "N", "TRAJ21": "N", "TRAJ22": "D",
    "TRAJ23": "N", "TRAJ24": "D", "TRAJ25": "N", "TRAJ26": "Y", "TRAJ27": "N", "TRAJ28": "N",
    "TRAJ29": "N", "TRAJ30": "N", "TRAJ31": "N", "TRAJ32": "N", "TRAJ33": "D", "TRAJ34": "N",
    "TRAJ35": "H", "TRAJ36": "Y", "TRAJ37": "D", "TRAJ38": "N", "TRAJ39": "H", "TRAJ40": "N",
    "TRAJ41": "H", "TRAJ42": "N", "TRAJ43": "N", "TRAJ44": "D", "TRAJ45": "Y", "TRAJ46": "N",
    "TRAJ47": "Y", "TRAJ48": "N", "TRAJ49": "N", "TRAJ50": "N", "TRAJ52": "N", "TRAJ53": "N",
    "TRAJ54": "N", "TRAJ56": "D", "TRAJ57": "Y", "TRAJ58": "D",
}

# Mouse (GRCm39 Ensembl 112) per-TRAJ α J→C splice-junction residue: terminal
# J nucleotide + Trac exon 1 5' dinucleotide "AT" -> AAT/AAC (N), GAT/GAC (D),
# TAT/TAC (Y), or CAT/CAC (H).
MOUSE_TRAJ_JUNCTION_RESIDUES: dict[str, str] = {
    "TRAJ1": "D", "TRAJ2": "D", "TRAJ3": "H", "TRAJ4": "Y", "TRAJ5": "N",
    "TRAJ6": "Y", "TRAJ7": "N", "TRAJ8": "Y", "TRAJ9": "N", "TRAJ11": "D",
    "TRAJ12": "D", "TRAJ13": "N", "TRAJ14": "H", "TRAJ15": "N", "TRAJ16": "H",
    "TRAJ17": "D", "TRAJ18": "D", "TRAJ19": "N", "TRAJ20": "D", "TRAJ21": "N",
    "TRAJ22": "D", "TRAJ23": "N", "TRAJ24": "D", "TRAJ25": "N", "TRAJ26": "Y",
    "TRAJ27": "N", "TRAJ28": "N", "TRAJ29": "N", "TRAJ30": "N", "TRAJ31": "N",
    "TRAJ32": "N", "TRAJ33": "D", "TRAJ34": "N", "TRAJ35": "Y", "TRAJ36": "H",
    "TRAJ37": "D", "TRAJ38": "N", "TRAJ39": "D", "TRAJ40": "H", "TRAJ41": "D",
    "TRAJ42": "N", "TRAJ43": "N", "TRAJ44": "D", "TRAJ45": "Y", "TRAJ46": "N",
    "TRAJ47": "H", "TRAJ48": "N", "TRAJ49": "N", "TRAJ50": "N", "TRAJ51": "N",
    "TRAJ52": "N", "TRAJ53": "N", "TRAJ54": "N", "TRAJ55": "N", "TRAJ56": "D",
    "TRAJ57": "Y", "TRAJ58": "D", "TRAJ59": "Y", "TRAJ60": "Y", "TRAJ61": "N",
}


def _traj_junction_residue(j_gene: object, species: str = "human") -> str | None:
    """Look up the α J→C junction residue for a TRAJ gene (allele-insensitive),
    or None when the gene is unknown/unmapped."""
    if j_gene is None or (isinstance(j_gene, float) and pd.isna(j_gene)):
        return None
    key = re.sub(r"(?<=[A-Z/-])0+(\d)", r"\1", str(j_gene).strip().upper().split("*", 1)[0])
    table = MOUSE_TRAJ_JUNCTION_RESIDUES if species == "mouse" else TRAJ_JUNCTION_RESIDUES
    return table.get(key)


_VALID_JUNCTION_RESIDUES: dict[str, frozenset] = {
    "alpha": frozenset("NYDH"),
    "beta": frozenset("E"),
}


def _fallback_junction_residue(
    chain: str, j_gene: object = None, *, species: str = "human"
) -> tuple[str, str]:
    """J→C junction residue + provenance when no contig covers the seam (#242)."""
    if chain == "alpha":
        r = _traj_junction_residue(j_gene, species=species)
        if r:
            return r, "j_inferred"
        return "N", "canonical_fallback"
    return _FALLBACK_JUNCTION_RESIDUE.get(chain, ""), "canonical_fallback"


def _apply_cohort_alpha_junctions(
    df: pd.DataFrame, *, codon_species: str = "human"
) -> pd.DataFrame:
    """Replace the blanket-``N`` α junction fallback with a per-J consensus read
    from the donor's OWN contig-verified clones (#242)."""
    src_col, j_col = "alpha_junction_residue_source", "alpha_j_gene"
    if src_col not in df.columns or j_col not in df.columns or df.empty:
        return df
    by_j: dict = {}
    verified = df[df[src_col] == "contig"]
    for j_gene, grp in verified.groupby(j_col, dropna=True):
        res = grp["alpha_junction_residue"].dropna().astype(str)
        if len(res):
            by_j[str(j_gene)] = res.mode().iloc[0]
    if not by_j:
        return df
    for idx in df.index[df[src_col] == "canonical_fallback"]:
        j_gene = df.at[idx, j_col]
        r = by_j.get(str(j_gene))
        cur = df.at[idx, "alpha_junction_residue"]
        if not r or r == cur:
            continue
        const_aa = df.at[idx, "alpha_constant_aa"]
        if isinstance(const_aa, str) and const_aa:
            new_aa = r + const_aa[1:]
            df.at[idx, "alpha_constant_aa"] = new_aa
            full_aa = df.at[idx, "full_alpha_aa"]
            if isinstance(full_aa, str) and full_aa.endswith(const_aa):
                df.at[idx, "full_alpha_aa"] = full_aa[: -len(const_aa)] + new_aa
        new_codon = back_translate(r, species=codon_species)
        full_nt = df.at[idx, "full_alpha_nt"] if "full_alpha_nt" in df.columns else None
        old_const_nt = df.at[idx, "alpha_constant_nt"] if "alpha_constant_nt" in df.columns else None
        if isinstance(full_nt, str) and isinstance(old_const_nt, str) and full_nt.endswith(old_const_nt):
            df.at[idx, "full_alpha_nt"] = (
                full_nt[: -len(old_const_nt)] + new_codon + old_const_nt[3:]
            )
        full_opt = (
            df.at[idx, "full_alpha_nt_optimized"]
            if "full_alpha_nt_optimized" in df.columns
            else None
        )
        old_const_opt = (
            df.at[idx, "alpha_constant_nt_optimized"]
            if "alpha_constant_nt_optimized" in df.columns
            else None
        )
        if (
            isinstance(full_opt, str)
            and isinstance(old_const_opt, str)
            and full_opt.endswith(old_const_opt)
        ):
            df.at[idx, "full_alpha_nt_optimized"] = (
                full_opt[: -len(old_const_opt)] + new_codon + old_const_opt[3:]
            )
        for nt_col in ("alpha_constant_nt", "alpha_constant_nt_optimized"):
            nt = df.at[idx, nt_col] if nt_col in df.columns else None
            if isinstance(nt, str) and len(nt) >= 3:
                df.at[idx, nt_col] = new_codon + nt[3:]
        df.at[idx, "alpha_junction_residue"] = r
        df.at[idx, src_col] = "cohort_j_consensus"
    return df


def _load_canonical_alleles_fasta(
    filename: str = _CANONICAL_ALLELES_FASTA,
) -> dict[str, dict[str, str]]:
    """Parse a packaged multi-allele FASTA into a ``{gene: {allele: AA}}`` dict."""
    fasta_text = (
        _pkg_files("tcrsift.refseqs").joinpath(filename).read_text()
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
            f"{filename} is missing canonical gene entries: "
            f"{sorted(missing)}"
        )
    return out


HUMAN_CONSTANT_ALLELES: dict[str, dict[str, str]] = _load_canonical_alleles_fasta(
    _CANONICAL_ALLELES_FASTA
)
MOUSE_CONSTANT_ALLELES: dict[str, dict[str, str]] = _load_canonical_alleles_fasta(
    _CANONICAL_ALLELES_MOUSE_FASTA
)
CONSTANT_ALLELES_BY_SPECIES: dict[str, dict[str, dict[str, str]]] = {
    "human": HUMAN_CONSTANT_ALLELES,
    "mouse": MOUSE_CONSTANT_ALLELES,
}


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


def _count_divergence_positions(value) -> int:
    """Count AA divergences in a ``{chain}_allele_divergence_positions`` cell.

    The column is ``"3:N->K;7:V->I"`` (or None/empty when no contig was
    available or the donor matched the canonical allele). Returns the number of
    diverging positions (#187).
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return 0
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return 0
    return len([p for p in text.split(";") if p])


def _warn_constant_allele_divergence(df: pd.DataFrame, *, clone_col: str = "CDR3ab") -> None:
    """Loudly surface constant-allele divergences (#187).

    The donor's *gene* may be correct (e.g. the J-parity TRBC rule) while the
    *allele* diverges from canonical — a silent fidelity loss in the assembled
    constant. Emit a verbose-independent warning naming how many clones diverge
    and an example, so it's visible by default rather than buried in a column.
    No-op when no divergence columns are present or nothing diverges.
    """
    diverged: list[tuple[str, str, int, str]] = []  # (clone, chain, n, positions)
    for chain in ("alpha", "beta"):
        cnt_col = f"{chain}_constant_allele_divergence"
        pos_col = f"{chain}_allele_divergence_positions"
        if cnt_col not in df.columns:
            continue
        mask = df[cnt_col] > 0
        for idx in df.index[mask]:
            clone = str(df.at[idx, clone_col]) if clone_col in df.columns else str(idx)
            diverged.append(
                (clone, chain, int(df.at[idx, cnt_col]),
                 str(df.at[idx, pos_col]) if pos_col in df.columns else "")
            )
    if not diverged:
        return
    ex_clone, ex_chain, ex_n, ex_pos = diverged[0]
    logger.warning(
        "Constant-allele divergence: %d clone-chain(s) carry a sequenced "
        "constant that diverges from the assembled canonical allele "
        "(donor allele differs). e.g. %s %s: %d aa (%s). The gene call may be "
        "correct while the allele isn't — see {chain}_allele_divergence_positions / "
        "{chain}_constant_allele_divergence.",
        len(diverged), ex_clone, ex_chain, ex_n, ex_pos,
    )


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
    best_score = -1.0
    best_n_agree = -1
    for allele, allele_aa in pool.items():
        n_agree, n_compared = _score_allele_against_contig(contig_aa, allele_aa)
        if n_compared >= min_compared:
            score = n_agree / n_compared
            all_scores[allele] = round(score, 3)
            # Prefer higher fractional agreement (score = n_agree / n_compared),
            # not raw n_agree: alleles can differ in length, so a longer allele
            # compared over more positions could win on absolute count despite a
            # strictly worse fraction. Tie-break by more positions matched, then
            # FASTA order (dict iteration is insertion order in Python 3.7+).
            if score > best_score or (score == best_score and n_agree > best_n_agree):
                best_score = score
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
for _gene, _aa in (
    ("TRAC", MOUSE_TRAC_AA),
    ("TRBC1", MOUSE_TRBC1_AA),
    ("TRBC2", MOUSE_TRBC2_AA),
):
    _expected_tail = MOUSE_CONSTANT_REGION_ENDINGS[_gene]
    if not _aa.endswith(_expected_tail):
        raise RuntimeError(
            f"{_CANONICAL_CONSTANTS_MOUSE_FASTA}: {_gene} canonical AA tail {_aa[-15:]!r} "
            f"doesn't match MOUSE_CONSTANT_REGION_ENDINGS[{_gene!r}]={_expected_tail!r}"
        )
del _gene, _aa, _expected_tail

HUMAN_CONSTANT_REGIONS_AA: dict[str, str] = {
    "TRAC": HUMAN_TRAC_AA,
    "TRBC1": HUMAN_TRBC1_AA,
    "TRBC2": HUMAN_TRBC2_AA,
}

MOUSE_CONSTANT_REGIONS_AA: dict[str, str] = {
    "TRAC": MOUSE_TRAC_AA,
    "TRBC1": MOUSE_TRBC1_AA,
    "TRBC2": MOUSE_TRBC2_AA,
}

CONSTANT_REGIONS_AA_BY_SPECIES: dict[str, dict[str, str]] = {
    "human": HUMAN_CONSTANT_REGIONS_AA,
    "mouse": MOUSE_CONSTANT_REGIONS_AA,
}




def pick_canonical_constant(
    chain: str,
    c_gene: str | None = None,
    j_gene: str | None = None,
    *,
    species: str = "human",
) -> tuple[str, str]:
    """Choose the canonical constant-region AA for a chain and species."""
    regions = CONSTANT_REGIONS_AA_BY_SPECIES.get(species, HUMAN_CONSTANT_REGIONS_AA)
    if chain == "alpha":
        return "TRAC", regions["TRAC"]

    j_family = _beta_j_family(j_gene)
    expected_c = BETA_JC_PARITY.get(j_family) if j_family else None
    if expected_c == "TRBC1":
        return "TRBC1", regions["TRBC1"]
    if expected_c == "TRBC2":
        return "TRBC2", regions["TRBC2"]

    c_base = _beta_c_base(c_gene)
    if c_base == "TRBC2":
        return "TRBC2", regions["TRBC2"]
    if c_base == "TRBC1":
        return "TRBC1", regions["TRBC1"]
    return "TRBC1", regions["TRBC1"]


def _beta_j_family(j_gene: str | None) -> str:
    """Parse the J-family prefix (``"TRBJ1"`` / ``"TRBJ2"``) from a
    raw J-gene call. Tolerates lowercase, whitespace, zero-padded segment
    numbers (``"TRBJ02-01"``), and allele suffixes (``"TRBJ1-1*02"``). Returns
    an empty string when the input doesn't look like a TRBJ call."""
    if not isinstance(j_gene, str):
        return ""
    token = re.sub(r"(?<=[A-Z/-])0+(\d)", r"\1", j_gene.strip().upper().split("*", 1)[0])
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


def _ensure_qc_warnings_column(df: pd.DataFrame) -> pd.DataFrame:
    """Guarantee a ``qc_warnings`` column exists with a list per row (#278).

    The assembler stashes ``qc_warnings`` only on the clones that raise one, so
    the column is absent entirely on a clean cohort — leaving the assembled
    schema donor-dependent (144 vs 145 columns) and breaking naive
    ``df["qc_warnings"]`` access and cross-cohort ``concat``. This fills in an
    empty list for every row that has no warning (NaN where the column was added
    late, or simply absent), so downstream code can always read the column and
    every cell is a list.
    """
    if "qc_warnings" not in df.columns:
        df["qc_warnings"] = [[] for _ in range(len(df))]
        return df
    df["qc_warnings"] = [
        v if isinstance(v, list) else [] for v in df["qc_warnings"]
    ]
    return df


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

