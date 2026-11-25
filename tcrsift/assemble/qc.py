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

from ..validation import TCRsiftValidationError
from .codons import (
    LINKERS,
    SYNTHESIS_RESTRICTION_SITES,
    translate_dna,
)
from .constants import (
    _VALID_AA_CHARS,
    _VALID_JUNCTION_RESIDUES,
    BETA_JC_PARITY,
    CONSTANT_AA_FLOOR,
    CONSTANT_REGION_ENDINGS,
    CONSTANT_REGION_ENDINGS_BY_SPECIES,
    CONSTANT_REGIONS_AA_BY_SPECIES,
    HUMAN_CONSTANT_REGIONS_AA,
    MOUSE_CONSTANT_REGION_ENDINGS,
    MOUSE_CONSTANT_REGIONS_AA,
    _beta_c_base,
    _beta_j_family,
    _resolve_c_gene,
    _traj_junction_residue,
    verify_canonical_constant_start,
)

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
    all_from_data = True
    for ch in ("alpha", "beta"):
        col = f"{ch}_junction_residue_source"
        if col in df.columns:
            non_contig = df[col].dropna()
            non_contig = non_contig[non_contig != "contig"]
            if len(non_contig) and not (non_contig == "from-data").all():
                all_from_data = False
        elif s[f"{ch}_fallback"] > 0:
            all_from_data = False
    if all_from_data:
        detail_note = "caller-supplied from-data constants (not verified against CellRanger contigs)"
    else:
        junction_note = (
            "J→C junction residues defaulted (α blanket N is wrong for ~30% of J genes)"
            if has_blanket_alpha or "alpha_junction_residue_source" not in df.columns
            else "J→C junction residues germline-inferred (not contig-read)"
        )
        detail_note = f"reference constants, {junction_note}, donor alleles unverified"
    summary = (
        f"{s['n_unverified']}/{s['n_total']} constructs assembled with NO contig "
        f"verification (α fallback: {s['alpha_fallback']}, β fallback: "
        f"{s['beta_fallback']}) — {detail_note}."
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

    from ..species import detect_vdj_species, normalize_species

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

            const_val = row.get(f"{chain}_constant_aa")
            has_const = isinstance(const_val, str) and bool(const_val)
            vdj_val = row.get(f"vdj_{chain}_aa") or row.get(f"VDJ_{chain}_aa")
            const_omitted = (
                (not has_const)
                and isinstance(vdj_val, str)
                and bool(vdj_val)
                and seq.rstrip("*").endswith(vdj_val.strip().rstrip("*"))
            )
            min_len = 80 if const_omitted else 200
            if len(seq) < min_len:
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

            if not const_omitted:
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

        # Detect VDJ-only assembly (include_constant=False): no constant column
        # populated and every full sequence ends with its VDJ_aa.
        const_col = f"{chain}_constant_aa"
        has_any_const = (
            const_col in df.columns
            and sub[const_col].apply(lambda s: isinstance(s, str) and bool(s)).any()
        )
        vdj_col_qc = f"vdj_{chain}_aa" if f"vdj_{chain}_aa" in sub.columns else f"VDJ_{chain}_aa"
        all_vdj_only = (
            not has_any_const
            and vdj_col_qc in sub.columns
            and sub.apply(
                lambda r: isinstance(r.get(vdj_col_qc), str)
                and bool(r.get(vdj_col_qc))
                and str(r.get(f"full_{chain}_aa", "")).rstrip("*").endswith(
                    str(r.get(vdj_col_qc, "")).strip().rstrip("*")
                ),
                axis=1,
            ).all()
        )

        # Terminal residue & 200-aa length floor apply when constants are assembled.
        if not all_vdj_only:
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
        else:
            lengths = sub_full.str.len()
            len_pass = lengths >= 80
            _add(
                f"{chain}_length_floor",
                f"{chain_label} VDJ length ≥ 80",
                chain,
                len_pass,
                total,
                median_value=float(lengths.median()),
                unit="aa",
            )

        # Constant length floor.
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
