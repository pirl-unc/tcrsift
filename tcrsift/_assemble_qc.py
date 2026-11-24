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

"""Post-assembly structural validation, synthesis-hazard QC, novel-allele audit, and FASTA export."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ._assemble_codons import (
    LINKERS,
    SYNTHESIS_RESTRICTION_SITES,
    translate_dna,
)
from ._assemble_constants import (
    _MAX_DIVERGENCE_SCAN_POSITIONS,
    _VALID_AA_CHARS,
    _VALID_JUNCTION_RESIDUES,
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_DIVERGENT_CONTIG,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_SPARSE_CONTIG,
    BETA_JC_PARITY,
    CONSTANT_AA_FLOOR,
    CONSTANT_ALLELES_BY_SPECIES,
    CONSTANT_REGION_ENDINGS,
    CONSTANT_REGION_ENDINGS_BY_SPECIES,
    CONSTANT_REGIONS_AA_BY_SPECIES,
    HUMAN_CONSTANT_ALLELES,
    HUMAN_CONSTANT_REGIONS_AA,
    MOUSE_CONSTANT_REGION_ENDINGS,
    MOUSE_CONSTANT_REGIONS_AA,
    _beta_c_base,
    _beta_j_family,
    _resolve_c_gene,
    _score_allele_against_contig,
    _traj_junction_residue,
    verify_canonical_constant_start,
)
from .validation import TCRsiftValidationError

logger = logging.getLogger(__name__)

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
    cand_positions = np.flatnonzero(expected_s.notna().to_numpy())
    if len(cand_positions) == 0:
        return messages

    canon_col_idx = df.columns.get_loc("beta_c_gene_canonical")
    for pos in cand_positions:
        idx = df.index[pos]
        row = df.iloc[pos]
        current = _resolve_c_gene(row, "beta")
        current_base = _beta_c_base(current)
        expected_c = expected_s.iloc[pos]
        if current_base and current_base != expected_c:
            df.iat[pos, canon_col_idx] = expected_c
            messages.append(
                ValidationMessage(
                    f"Clone {idx}: autocorrected β c_gene {current!r} → "
                    f"{expected_c} based on {row['beta_j_gene']} (J family rules)",
                    idx=idx,
                    severity="autocorrect",
                )
            )
    return messages


def _unverified_chains(row) -> list:
    """Chains present in ``row`` whose J→C junction residue was NOT read from
    the contig (or whose V-domain required IMGT germline rescue/reconstruction)."""
    out = []
    for chain in ("alpha", "beta"):
        if not isinstance(row.get(f"full_{chain}_aa"), str) or not row.get(f"full_{chain}_aa"):
            continue
        vdj_src = row.get(f"{chain}_vdj_source")
        if vdj_src in ("contig_5prime_germline_rescued", "germline_vj_reconstructed"):
            out.append(chain)
            continue
        src = row.get(f"{chain}_junction_residue_source")
        has_const = isinstance(row.get(f"{chain}_constant_aa"), str) and bool(row.get(f"{chain}_constant_aa"))
        if src is not None and src != "contig":
            out.append(chain)
        elif src is None and has_const:
            out.append(chain)
    return out


def assemble_fidelity_summary(df: pd.DataFrame) -> dict:
    """Summarize per-construct contig-verification across an assembled frame.

    Returns ``{n_total, n_verified, n_unverified, alpha_fallback, beta_fallback,
    unverified_idx}``. A construct is unverified when any present chain's J→C
    junction residue came from the canonical fallback rather than the donor's
    contig (#241/#243/#244) — meaning reference constants, a defaulted junction,
    and an unverified donor allele.
    """
    n_total = len(df)
    unverified_idx, a_fb, b_fb = [], 0, 0
    for idx, row in df.iterrows():
        chains = _unverified_chains(row)
        if chains:
            unverified_idx.append(idx)
        a_fb += "alpha" in chains
        b_fb += "beta" in chains
    return {
        "n_total": n_total,
        "n_verified": n_total - len(unverified_idx),
        "n_unverified": len(unverified_idx),
        "alpha_fallback": a_fb,
        "beta_fallback": b_fb,
        "unverified_idx": unverified_idx,
    }


def enforce_contig_fidelity(
    df: pd.DataFrame,
    *,
    allow_canonical_fallback: bool = False,
    context: str = "assembly",
) -> str:
    """Fail-closed gate for synthesis-ready constructs (#243/#244).

    Builds a cohort fidelity summary. When any construct is not contig-verified
    (canonical fallback: reference constants, defaulted junction residue,
    unverified donor allele) and ``allow_canonical_fallback`` is False, raises
    :class:`TCRsiftValidationError` — refusing to silently emit constructs whose
    junctions/alleles can't be verified. When the fallback IS allowed, returns a
    loud summary string for the caller to surface as a warning instead (#241).

    Returns the summary string (empty when everything is contig-verified).
    """
    s = assemble_fidelity_summary(df)
    if s["n_unverified"] == 0:
        return ""
    has_blanket_alpha = (
        "alpha_junction_residue_source" in df.columns
        and (df["alpha_junction_residue_source"] == "canonical_fallback").any()
    )
    junction_note = (
        "J→C junction residues defaulted (α blanket N is wrong for ~30% of J genes)"
        if has_blanket_alpha or "alpha_junction_residue_source" not in df.columns
        else "J→C junction residues germline-inferred (not contig-read)"
    )
    summary = (
        f"{s['n_unverified']}/{s['n_total']} constructs assembled with NO contig "
        f"verification (α fallback: {s['alpha_fallback']}, β fallback: "
        f"{s['beta_fallback']}) — reference constants, {junction_note}, donor alleles "
        f"unverified."
    )
    if not allow_canonical_fallback:
        raise TCRsiftValidationError(
            f"{context}: {summary}",
            hint="Pass --cellranger-dir / --contigs-dir so junctions and alleles "
            "are read from the donor's contigs, or pass --allow-canonical-fallback "
            "to consciously accept reference/approximate constants.",
        )
    return summary


def _infer_constant_species(
    const_aa: object, full_aa: object, c_base: str, default: str = "human"
) -> str:
    """Infer whether an assembled chain carries a human or mouse constant region."""
    mouse_bare = MOUSE_CONSTANT_REGIONS_AA.get(c_base)
    if mouse_bare:
        mb8 = mouse_bare[:8]
        if isinstance(const_aa, str) and (
            const_aa.startswith(mb8) or const_aa[1:].startswith(mb8)
        ):
            return "mouse"
    human_bare = HUMAN_CONSTANT_REGIONS_AA.get(c_base)
    if human_bare:
        hb8 = human_bare[:8]
        if isinstance(const_aa, str) and (
            const_aa.startswith(hb8) or const_aa[1:].startswith(hb8)
        ):
            return "human"
    if c_base in ("TRBC1", "TRBC2") and isinstance(full_aa, str) and full_aa:
        s = full_aa.rstrip("*")
        if s.endswith(MOUSE_CONSTANT_REGION_ENDINGS[c_base]):
            return "mouse"
        if s.endswith(CONSTANT_REGION_ENDINGS[c_base]):
            return "human"
    return default


def validate_sequences(
    df: pd.DataFrame,
    strict: bool = False,
    fix: bool = False,
    *,
    species: str | None = None,
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

    from .species import detect_vdj_species, normalize_species

    _sp = normalize_species(species, allow_auto=True) if species else None
    explicit_vdj_sp = None if _sp == "auto" else _sp
    detected_vdj_sp = explicit_vdj_sp or detect_vdj_species(df)

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
            elif row.get(f"{chain}_leader_source") == "missing_from_contig":
                _lb(idx,
                    f"full_{chain}_aa is missing its requested from_contig "
                    f"signal peptide (starts with {seq[:5]!r})")

            # Surface a from_contig signal-peptide QC flag (#263) as an
            # informational note — the extracted leader is kept (it's the best
            # available), but a weak-Kozak over-capture / out-of-range length /
            # missing h-region shouldn't pass silently. Never load-bearing here:
            # the assembler already substitutes a curated SP when leader_fallback
            # is set; this just keeps the flag visible in a validation pass.
            lqc = row.get(f"{chain}_leader_qc")
            if isinstance(lqc, str) and lqc not in ("ok", "curated_fallback"):
                support = row.get(f"{chain}_leader_support")
                _info(idx,
                    f"{chain} signal peptide QC={lqc} "
                    f"(len {row.get(f'{chain}_leader_len')}, "
                    f"kozak {row.get(f'{chain}_leader_kozak_score')}, "
                    f"support {support})")

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

            # J→C junction-seam guard (#235). The constant must begin with a
            # valid junction residue (β: E; α: N/Y/D/H) followed by the bare
            # canonical — NOT the bare canonical directly. A direct-canonical
            # start means the junction residue was dropped and the chain is 1 aa
            # short at the seam. This is the check the prior self-consistent
            # canonical-start comparison (observed vs the row's own constant_aa)
            # could not make, which is how #235 shipped silently.
            c_base = _resolve_c_gene(row, chain).strip().upper().split("*")[0]
            row_const_sp = _infer_constant_species(
                const, seq, c_base, default=explicit_vdj_sp or detected_vdj_sp or "human"
            )
            bare = CONSTANT_REGIONS_AA_BY_SPECIES[row_const_sp].get(c_base)
            if isinstance(const, str) and const and bare:
                bare8 = bare[:8]
                if const.startswith(bare8):
                    _lb(idx,
                        f"{chain}_constant_aa starts at the bare canonical "
                        f"{bare8!r} — J→C junction residue missing (chain 1 aa "
                        f"short at the seam)")
                elif not (
                    const[0] in _VALID_JUNCTION_RESIDUES.get(chain, frozenset())
                    and const[1:].startswith(bare8)
                ):
                    _lb(idx,
                        f"{chain}_constant_aa has an unexpected J→C seam: "
                        f"{const[:9]!r} (expected a valid junction residue "
                        f"{sorted(_VALID_JUNCTION_RESIDUES.get(chain, ()))} "
                        f"+ canonical {bare8!r})")
                elif chain == "alpha":
                    j_gene_val = row.get("alpha_j_gene")
                    if explicit_vdj_sp:
                        expected_set = {
                            _traj_junction_residue(j_gene_val, species=explicit_vdj_sp)
                        } - {None}
                    elif detected_vdj_sp and row_const_sp == detected_vdj_sp:
                        expected_set = {
                            _traj_junction_residue(j_gene_val, species=detected_vdj_sp)
                        } - {None}
                    else:
                        # Auto-detected species can be ambiguous on uppercased mouse
                        # TRAJ calls or chimeric constructs; accept either germline.
                        expected_set = {
                            _traj_junction_residue(j_gene_val, species="human"),
                            _traj_junction_residue(j_gene_val, species="mouse"),
                        } - {None}
                    if expected_set and const[0] not in expected_set:
                        expected_disp = sorted(expected_set)[0] if len(expected_set) == 1 else sorted(expected_set)
                        _lb(idx,
                            f"alpha J→C junction {const[0]!r} disagrees with the "
                            f"germline residue {expected_disp!r} for "
                            f"{str(j_gene_val)!r} — likely a "
                            f"defaulted/misread junction (#242)")

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

            c_gene = _resolve_c_gene(row, chain)
            c_gene_base = c_gene.strip().upper().split("*")[0]
            endings_map = CONSTANT_REGION_ENDINGS_BY_SPECIES[row_const_sp]
            if c_gene_base in endings_map:
                expected_end = endings_map[c_gene_base]
                seq_no_stop = seq.rstrip("*")
                if not seq_no_stop.endswith(expected_end):
                    _lb(idx,
                        f"{chain} doesn't end with canonical "
                        f"{c_gene_base} C-terminus ({expected_end!r}); "
                        f"got {seq_no_stop[-len(expected_end):]!r}")
                expected_const = row.get(f"{chain}_constant_aa")
                used_bare_fallback = False
                if not isinstance(expected_const, str) or not expected_const:
                    expected_const = CONSTANT_REGIONS_AA_BY_SPECIES[row_const_sp].get(c_gene_base)
                    used_bare_fallback = True
                if expected_const:
                    obs_start = _expected_constant_start_from_full(seq_no_stop, row, chain)
                    ok_start = verify_canonical_constant_start(
                        obs_start, expected_const, min_match=8
                    ) or (
                        used_bare_fallback
                        and len(obs_start) > 1
                        and verify_canonical_constant_start(
                            obs_start[1:], expected_const, min_match=8
                        )
                    )
                    if not ok_start:
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

    # 5. Surface any per-row qc_warnings the assembler stashed. These are
    # *informational, self-corrected* notes (fell back to canonical, allele
    # not called, contig diverged at the C-region boundary, …) — the
    # assembler already handled them and the assembled sequence is valid, so
    # they must NOT be load-bearing (that aborted `tcrsift run` on real
    # CellRanger data, #129). Genuine structural failures are caught by the
    # explicit checks above and remain load-bearing.
    if "qc_warnings" in df.columns:
        for idx, qcs in df["qc_warnings"].items():
            if isinstance(qcs, list):
                for msg in qcs:
                    _info(idx, str(msg))

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

        # Terminal residue (against canonical for the row's c_gene and constant species).
        def _ends_canonical(row, _chain=chain):
            s = row.get(f"full_{_chain}_aa", "")
            if not isinstance(s, str) or not s:
                return False
            c = _resolve_c_gene(row, _chain)
            base = c.strip().upper().split("*")[0]
            csp = _infer_constant_species(
                row.get(f"{_chain}_constant_aa"), s, base, default="human"
            )
            endings = CONSTANT_REGION_ENDINGS_BY_SPECIES[csp]
            s_clean = s.rstrip("*")
            if base in endings:
                return s_clean.endswith(endings[base])
            all_ends = set(CONSTANT_REGION_ENDINGS.values()) | set(
                MOUSE_CONSTANT_REGION_ENDINGS.values()
            )
            return any(s_clean.endswith(end) for end in all_ends)

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



def _gc_fraction(nt: str) -> float:
    s = str(nt).upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def _max_homopolymer(nt: str) -> int:
    s = str(nt).upper()
    best = cur = 0
    prev = ""
    for ch in s:
        cur = cur + 1 if ch == prev else 1
        prev = ch
        best = max(best, cur)
    return best


def _longest_repeat(nt: str, *, min_len: int = 8) -> int:
    """Longest substring (>= min_len) that occurs at least twice."""
    s = str(nt).upper()
    n = len(s)
    if n < min_len * 2:
        return 0

    def _has(length: int) -> bool:
        seen: set[str] = set()
        for i in range(n - length + 1):
            sub = s[i:i + length]
            if sub in seen:
                return True
            seen.add(sub)
        return False

    lo, hi, ans = min_len, n // 2, 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if _has(mid):
            ans = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return ans


def _restriction_hits(nt: str) -> str:
    # Scan BOTH strands. Non-palindromic Type IIS enzymes (BsaI GGTCTC,
    # BsmBI CGTCTC) cut from a site that may appear on the reverse strand —
    # i.e. as its reverse complement (GAGACC / GAGACG) in the given sequence.
    # A forward-only scan misses these; e.g. the P2A/F2A linker DNA carries a
    # reverse-strand BsmBI site that a Golden-Gate/Type IIS assembly would cut.
    s = str(nt).upper()
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return ";".join(
        nm for nm, site in SYNTHESIS_RESTRICTION_SITES.items()
        if site in s or site in rc
    )


def _terminal_stops(nt: str) -> list[str]:
    """Trailing stop codons of a CDS, outermost-last (e.g. ``["TAA", "TGA"]``)."""
    s = str(nt).upper()
    out: list[str] = []
    while len(s) >= 3 and s[-3:] in ("TAA", "TAG", "TGA"):
        out.insert(0, s[-3:])
        s = s[:-3]
    return out


def add_synthesis_qc(
    df: pd.DataFrame,
    *,
    nt_col: str | None = None,
    gc_range: tuple[float, float] = (0.35, 0.70),
) -> pd.DataFrame:
    """Flag synthesis hazards on the assembled construct NT (#206).

    Adds, per construct (scored on ``nt_col`` — default the codon-optimized
    single-chain CDS, else the unoptimized): ``synth_gc_fraction``,
    ``synth_max_homopolymer`` (longest single-base run), ``synth_max_repeat``
    (longest repeated substring), ``synth_restriction_sites`` (``;``-joined
    enzyme names hit), and ``synth_gc_ok`` (in ``gc_range``). Cross-construct:
    ``synth_duplicate_construct`` (identical ``single_chain_aa`` across picks)
    and ``synth_alpha_beta_swap`` (``CDR3_alpha == CDR3_beta``). No-op (returns
    a copy) when no construct-NT column is present. The codon optimizer already
    avoids these — this surfaces any that survive (esp. in the unoptimized CDS).
    """
    out = df.copy()
    col = nt_col or next(
        (c for c in ("single_chain_nt_optimized", "single_chain_nt") if c in out.columns),
        None,
    )
    if col is None:
        return out
    seqs = out[col].fillna("").astype(str)
    out["synth_gc_fraction"] = seqs.map(_gc_fraction).round(4)
    out["synth_max_homopolymer"] = seqs.map(_max_homopolymer)
    out["synth_max_repeat"] = seqs.map(_longest_repeat)
    out["synth_restriction_sites"] = seqs.map(_restriction_hits)
    # Terminal stops (#293 follow-up): make the construct's terminator visible
    # in the deliverable QC. The synthesis best-practice is two NON-REDUNDANT
    # stops (TAA;TGA — different release factors). `synth_terminal_stops` lists
    # them; `synth_dual_stop_ok` is True only with ≥2 distinct terminal stops.
    # (Internal stops are separately caught load-bearing by validate_sequences'
    # NT→AA roundtrip, so they don't need a column here.)
    out["synth_terminal_stops"] = seqs.map(lambda s: ";".join(_terminal_stops(s)))
    lo, hi = gc_range
    nonempty = seqs.str.len() > 0
    # Empty / missing constructs (e.g. an unpaired clone with no single_chain)
    # have nothing to terminate-check — treat them as OK, matching synth_gc_ok's
    # `| (~nonempty)` below so they don't inflate the "not ending in 2 stops" count.
    out["synth_dual_stop_ok"] = out["synth_terminal_stops"].map(
        lambda s: len(set(s.split(";"))) >= 2 if s else False
    ) | (~nonempty)
    out["synth_gc_ok"] = ((out["synth_gc_fraction"] >= lo) & (out["synth_gc_fraction"] <= hi)) | (~nonempty)
    if "single_chain_aa" in out.columns:
        sc = out["single_chain_aa"].fillna("").astype(str)
        out["synth_duplicate_construct"] = sc.duplicated(keep=False) & (sc != "")
    if {"CDR3_alpha", "CDR3_beta"} <= set(out.columns):
        a = out["CDR3_alpha"].astype(str)
        b = out["CDR3_beta"].astype(str)
        out["synth_alpha_beta_swap"] = (a == b) & (a != "") & (a.str.lower() != "nan")
    return out


def synthesis_qc_report(df: pd.DataFrame) -> str:
    """One-line-per-check synthesis-hazard tally (#206); '' if QC not run."""
    if "synth_gc_fraction" not in df.columns:
        return ""
    n = len(df)
    lines = ["[synthesis QC]"]
    n_gc = int((~df["synth_gc_ok"]).sum()) if "synth_gc_ok" in df.columns else 0
    lines.append(f"  GC outside window: {n_gc}/{n}")
    if "synth_max_homopolymer" in df.columns:
        lines.append(f"  homopolymer run >=9: {int((df['synth_max_homopolymer'] >= 9).sum())}/{n}")
    if "synth_max_repeat" in df.columns:
        lines.append(f"  repeat >=20 nt: {int((df['synth_max_repeat'] >= 20).sum())}/{n}")
    if "synth_restriction_sites" in df.columns:
        lines.append(f"  restriction site(s): {int((df['synth_restriction_sites'] != '').sum())}/{n}")
    if "synth_dual_stop_ok" in df.columns:
        lines.append(
            f"  NOT ending in 2 distinct stops: {int((~df['synth_dual_stop_ok']).sum())}/{n}"
        )
    if "synth_duplicate_construct" in df.columns:
        lines.append(f"  duplicate construct: {int(df['synth_duplicate_construct'].sum())}/{n}")
    if "synth_alpha_beta_swap" in df.columns:
        lines.append(f"  α/β CDR3 identical: {int(df['synth_alpha_beta_swap'].sum())}/{n}")
    return "\n".join(lines)


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


def _candidate_constant_alleles_for_chain(
    chain: str, *, species: str = "human"
) -> list[tuple[str, str, str]]:
    """Return ``(gene, allele_label, aa)`` references compatible with chain."""
    pools = CONSTANT_ALLELES_BY_SPECIES.get(species, HUMAN_CONSTANT_ALLELES)
    if chain == "alpha":
        genes = ("TRAC",)
    elif chain == "beta":
        genes = ("TRBC1", "TRBC2")
    else:
        genes = tuple(pools)

    refs: list[tuple[str, str, str]] = []
    for gene in genes:
        for allele, aa in pools.get(gene, {}).items():
            refs.append((gene, f"{gene}*{allele}", aa))
    return refs


def _default_reference_for_gene(
    canonical_name: object, *, species: str = "human"
) -> tuple[str | None, str | None]:
    """Return the default packaged allele label + AA for a known gene."""
    gene = str(canonical_name)
    pools = CONSTANT_ALLELES_BY_SPECIES.get(species, HUMAN_CONSTANT_ALLELES)
    regions = CONSTANT_REGIONS_AA_BY_SPECIES.get(species, HUMAN_CONSTANT_REGIONS_AA)
    pool = pools.get(gene, {})
    if pool:
        allele, aa = next(iter(pool.items()))
        return f"{gene}*{allele}", aa
    aa = regions.get(gene)
    if aa:
        return gene, aa
    return None, None


def _closest_reference_for_observations(
    chain: str,
    observations: pd.Series,
    *,
    species: str = "human",
) -> tuple[str | None, str | None]:
    """Pick the nearest packaged constant allele for a novel gene label."""
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
    for _, label, allele_aa in _candidate_constant_alleles_for_chain(
        chain, species=species
    ):
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
    *,
    species: str = "human",
) -> tuple[str | None, str | None]:
    """Resolve the reference allele used to describe divergences."""
    reference_label, reference_aa = _default_reference_for_gene(
        canonical_name, species=species
    )
    if reference_aa:
        return reference_label, reference_aa
    return _closest_reference_for_observations(
        chain, observations, species=species
    )


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
        keep_cols = [called_col, divs_col, gene_col]
        if v_col in df.columns:
            keep_cols.append(v_col)
        if "samples" in df.columns:
            keep_cols.append("samples")
        stored = df.loc[
            df[divs_col].notna() & df[called_col].notna(),
            keep_cols,
        ].copy()
        if "samples" not in stored.columns:
            stored["samples"] = ""
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
            const_col = f"{chain}_constant_aa"
            full_col = f"full_{chain}_aa"
            first_const = gene_slice[const_col].iloc[0] if const_col in gene_slice.columns else None
            first_full = gene_slice[full_col].iloc[0] if full_col in gene_slice.columns else None
            slice_sp = _infer_constant_species(
                first_const, first_full, str(canonical_name).upper(), default="human"
            )
            reference_label, reference_aa = _reference_for_gene_or_observations(
                chain, canonical_name, gene_slice[obs_col], species=slice_sp
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
    # NOTE: ``reference_allele`` is deliberately NOT a groupby key. The same
    # physical variant — identical (chain, gene, position, expected_aa,
    # observed_aa) — can arrive labelled against different reference alleles
    # (Source 1 tags it with the *called* allele, e.g. TRBC2*03; Source 2
    # recomputes against the gene *default*, TRBC2*01). ``expected_aa`` already
    # pins the reference residue at that position, so equal-expected rows are
    # the same variant and must aggregate together; keying on reference_allele
    # would split them and deflate n_clones/n_v_genes/n_samples below the
    # candidate thresholds. We pick a deterministic representative label below.
    grouped = long_df.groupby(
        [
            "chain", "gene", "position",
            "expected_aa", "observed_aa",
        ],
        dropna=False, sort=False, observed=True,
    )
    n_clones = grouped["clone_idx"].nunique().rename("n_clones")
    n_v_genes = grouped["v_gene"].nunique(dropna=True).rename("n_v_genes")
    n_samples = grouped["sample"].agg(_expand_samples_unique).rename("n_samples")
    reference_allele = grouped["reference_allele"].agg(
        lambda s: sorted(set(s.dropna().astype(str)))[0] if s.notna().any() else None
    ).rename("reference_allele")
    summary = pd.concat(
        [n_clones, n_v_genes, n_samples, reference_allele], axis=1
    ).reset_index()

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
    if ragged or len(translated) != n_codons or "*" in translated:
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
            if not isinstance(seq, str) or not seq.strip():
                continue

            # Build header
            cdr3ab = row.get("CDR3ab", idx)
            cdr3a = row.get("CDR3_alpha", "")
            cdr3b = row.get("CDR3_beta", "")

            header = f">{cdr3ab} CDR3a={cdr3a} CDR3b={cdr3b}"
            f.write(f"{header}\n{seq}\n")

    logger.info(f"Exported {len(df)} sequences to {output_path}")
