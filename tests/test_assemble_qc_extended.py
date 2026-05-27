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

"""Tests for the extended assembly QC pass (#67).

Beyond the canonical-ending / canonical-start / qc_warnings checks
already covered by `test_assemble_validate_strict.py`, #67 asks for:

- separate constant length floors (alpha ≥ 100, beta ≥ 150),
- β-chain J→C in-cis parity (TRBJ1-* with TRBC1, TRBJ2-* with TRBC2),
- byte-for-byte `full == leader + vdj + constant`,
- no premature stop codon,
- methionine start when a leader is present,
- standard residue alphabet,
- single-chain construct integrity,
- an aggregate `[assemble QC] X / Y` summary string.

The fixtures here use the same `_ok_clone` shape as the existing
strict-mode tests; they extend it with leader + constant_aa columns
so the byte-for-byte and length-floor checks have data to exercise.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.assemble import (
    BETA_JC_PARITY,
    CONSTANT_AA_FLOOR,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    HUMAN_TRBC2_AA,
    ValidationMessage,
    assemble_qc_report,
    fix_jc_parity,
    validate_sequences,
)
from tcrsift.validation import TCRsiftValidationError


def _ok_clone(
    *,
    leader_alpha: str = "M" + "A" * 19,                   # 20 aa, starts with M
    leader_beta: str = "M" + "A" * 19,
    vdj_alpha: str = "CASS" + "A" * 60 + "VLPHA",         # ~69 aa, valid AAs only
    vdj_beta: str = "CASS" + "G" * 30 + "VETA",           # ~38 aa, valid AAs only
    c_alpha_canonical: str = "TRAC",
    c_beta_canonical: str = "TRBC1",
    beta_j_gene: str = "TRBJ1-2",
) -> dict:
    """Build a row that passes every new check.

    Critical: `full_*_aa` is built as `leader + vdj + constant` so the
    byte-for-byte check passes; `*_constant_aa` is the canonical
    sequence so length-floor and terminal-residue checks pass too.
    """
    const_alpha = HUMAN_TRAC_AA
    const_beta = HUMAN_TRBC1_AA if c_beta_canonical == "TRBC1" else HUMAN_TRBC2_AA
    return {
        "CDR3_alpha": vdj_alpha,
        "CDR3_beta": vdj_beta,
        "alpha_leader_aa": leader_alpha,
        "beta_leader_aa": leader_beta,
        "vdj_alpha_aa": vdj_alpha,
        "vdj_beta_aa": vdj_beta,
        "alpha_constant_aa": const_alpha,
        "beta_constant_aa": const_beta,
        "alpha_c_gene": c_alpha_canonical,
        "beta_c_gene": c_beta_canonical,
        "alpha_c_gene_canonical": c_alpha_canonical,
        "beta_c_gene_canonical": c_beta_canonical,
        "beta_j_gene": beta_j_gene,
        "full_alpha_aa": leader_alpha + vdj_alpha + const_alpha,
        "full_beta_aa": leader_beta + vdj_beta + const_beta,
    }


class TestConstantLengthFloor:
    """The pre-#67 length check only looked at the full chain. A short
    constant could hide behind a long leader + VDJ. Per-chain constant
    floor catches that."""

    def test_alpha_floor_violation_caught(self):
        row = _ok_clone()
        row["alpha_constant_aa"] = "M" * 50  # well under 100
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("alpha_constant_aa too short" in w for w in warnings)

    def test_beta_floor_violation_caught(self):
        row = _ok_clone()
        row["beta_constant_aa"] = "M" * 120  # under 150
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("beta_constant_aa too short" in w for w in warnings)

    def test_strict_raises_on_floor_violation(self):
        row = _ok_clone()
        row["alpha_constant_aa"] = "RT"  # the #66 truncation
        df = pd.DataFrame([row])
        with pytest.raises(TCRsiftValidationError):
            validate_sequences(df, strict=True)

    def test_floor_constants_match_expected_chain_lengths(self):
        """Defensive: alpha floor < TRAC length, beta floor <
        TRBC1/2 length. Otherwise the OK clone trivially fails."""
        assert CONSTANT_AA_FLOOR["alpha"] < len(HUMAN_TRAC_AA)
        assert CONSTANT_AA_FLOOR["beta"] < len(HUMAN_TRBC1_AA)
        assert CONSTANT_AA_FLOOR["beta"] < len(HUMAN_TRBC2_AA)


class TestBetaJCParity:
    """TRBJ1-* must pair with TRBC1; TRBJ2-* with TRBC2."""

    def test_trbj1_with_trbc2_flagged(self):
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-5")
        row["beta_c_gene_canonical"] = "TRBC2"
        row["beta_c_gene"] = "TRBC2"
        # full_beta_aa still uses TRBC1 constant from the helper -
        # rebuild so we isolate the J→C mismatch (otherwise the
        # canonical-end check also fires).
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "J→C parity mismatch" in w and "TRBJ1" in w for w in warnings
        )

    def test_trbj2_with_trbc1_flagged(self):
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ2-3")
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "J→C parity mismatch" in w and "TRBJ2" in w for w in warnings
        )

    def test_in_cis_pair_passes(self):
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-2")
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("J→C parity mismatch" in w for w in warnings)

    def test_allele_suffix_on_c_gene_handled(self):
        """`TRBC2*01` is the common allele form — strip before the
        parity comparison."""
        row = _ok_clone(c_beta_canonical="TRBC2", beta_j_gene="TRBJ2-1")
        row["beta_c_gene_canonical"] = "TRBC2*01"
        row["beta_c_gene"] = "TRBC2*01"
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        # TRBJ2 + TRBC2 is correct; no parity mismatch should fire.
        assert not any("J→C parity mismatch" in w for w in warnings)

    def test_parity_table_complete(self):
        """Defensive: every TRBJ family should map to one TRBC."""
        assert BETA_JC_PARITY == {"TRBJ1": "TRBC1", "TRBJ2": "TRBC2"}


class TestAssembleEndToEndJCOverride:
    """End-to-end: when CellRanger's β c_gene conflicts with the J
    family, assemble_full_sequences must propagate J-gene through the
    assembly path and produce a full_beta_aa with the J-correct
    canonical C-terminus. Without J-gene propagation in
    ``_assemble_clone``, ``pick_canonical_constant`` silently sees
    j_gene=None and the override is bypassed in the real pipeline (#90).
    """

    def test_aggregate_override_warning_fires_at_end_of_assembly(self, caplog):
        """Regression: the aggregate `Overrode CellRanger TRBC call on
        N / M β clones` warning must actually be logged when verbose=True
        and overrides occurred. Existing tests check side effects on
        the output frame; this checks the audit trail."""
        import logging

        from tcrsift.assemble import assemble_full_sequences

        vdj_a = "CASS" + "A" * 100 + "EQFF"
        vdj_b = "CASS" + "G" * 100 + "EYFF"
        df = pd.DataFrame([
            # Two clones with cross-parity, one OK.
            {
                "CDR3ab": "c1", "CDR3_alpha": vdj_a, "CDR3_beta": vdj_b,
                "VDJ_alpha_aa": vdj_a, "VDJ_beta_aa": vdj_b,
                "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
                "beta_j_gene": "TRBJ1-1", "samples": "S1",
            },
            {
                "CDR3ab": "c2", "CDR3_alpha": vdj_a + "X", "CDR3_beta": vdj_b + "X",
                "VDJ_alpha_aa": vdj_a + "X", "VDJ_beta_aa": vdj_b + "X",
                "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
                "beta_j_gene": "TRBJ2-5", "samples": "S1",  # also cross-parity
            },
            {
                "CDR3ab": "c3", "CDR3_alpha": vdj_a + "Y", "CDR3_beta": vdj_b + "Y",
                "VDJ_alpha_aa": vdj_a + "Y", "VDJ_beta_aa": vdj_b + "Y",
                "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC1",
                "beta_j_gene": "TRBJ1-1", "samples": "S1",  # consistent
            },
        ])
        with caplog.at_level(logging.WARNING, logger="tcrsift.assemble"):
            assemble_full_sequences(
                df, alpha_leader=None, beta_leader=None,
                verbose=True, show_progress=False,
            )
        msgs = [r.getMessage() for r in caplog.records]
        assert any(
            "Overrode CellRanger TRBC call on 2 / 3" in m for m in msgs
        ), f"expected aggregate override warning, got: {msgs}"

    def test_aggregate_override_warning_silenced_when_verbose_false(
        self, caplog
    ):
        """`verbose=False` should suppress the audit warning — useful
        for batch / notebook contexts. The override still happens; only
        the log line is gated."""
        import logging

        from tcrsift.assemble import assemble_full_sequences

        vdj_a = "CASS" + "A" * 100 + "EQFF"
        vdj_b = "CASS" + "G" * 100 + "EYFF"
        df = pd.DataFrame([{
            "CDR3ab": "c1", "CDR3_alpha": vdj_a, "CDR3_beta": vdj_b,
            "VDJ_alpha_aa": vdj_a, "VDJ_beta_aa": vdj_b,
            "alpha_c_gene": "TRAC", "beta_c_gene": "TRBC2",
            "beta_j_gene": "TRBJ1-1", "samples": "S1",
        }])
        with caplog.at_level(logging.WARNING, logger="tcrsift.assemble"):
            out = assemble_full_sequences(
                df, alpha_leader=None, beta_leader=None,
                verbose=False, show_progress=False,
            )
        # Override still applied.
        assert out["beta_c_gene_canonical"].iloc[0] == "TRBC1"
        # But the audit log line wasn't emitted.
        assert not any(
            "Overrode CellRanger TRBC call" in r.getMessage()
            for r in caplog.records
        )

    @pytest.mark.parametrize(
        "j_gene,cellranger_c,expected_c,expected_tail",
        [
            ("TRBJ1-1", "TRBC2", "TRBC1", "VKRKDF"),
            ("TRBJ1-4", "TRBC2", "TRBC1", "VKRKDF"),
            ("TRBJ2-5", "TRBC1", "TRBC2", "VKRKDSRG"),
        ],
    )
    def test_cross_parity_overridden_through_full_assembly(
        self, j_gene, cellranger_c, expected_c, expected_tail
    ):
        from tcrsift.assemble import assemble_full_sequences

        vdj_a = "M" + "A" * 19 + "CASS" + "V" * 100 + "EQFF"
        vdj_b = "M" + "A" * 19 + "CASS" + "G" * 100 + "EYFF"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": "CASS" + "V" * 100 + "EQFF",
            "CDR3_beta":  "CASS" + "G" * 100 + "EYFF",
            "VDJ_alpha_aa": vdj_a,
            "VDJ_beta_aa":  vdj_b,
            "alpha_c_gene": "TRAC",
            "beta_c_gene":  cellranger_c,
            "beta_j_gene":  j_gene,
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            verbose=False, show_progress=False,
        )
        assert out["beta_c_gene_canonical"].iloc[0] == expected_c
        assert out["full_beta_aa"].iloc[0].endswith(expected_tail)
        # Raw CellRanger call preserved as audit trail.
        assert out["beta_c_gene"].iloc[0] == cellranger_c
        # And validate_sequences should no longer fire parity mismatch.
        msgs = validate_sequences(out, strict=False)
        assert not any(
            m.severity == "load_bearing" and "parity mismatch" in m
            for m in msgs
        )


class TestFromDataConstantWritesCanonical:
    """Regression: ``constant_source='from-data'`` was leaving
    ``{chain}_c_gene_canonical`` unset, which silently degraded the
    autocorrect (#89) and aggregate-override warning (#90) on that
    branch. After the fix, the canonical name should be recorded with
    the J-family override applied."""

    def test_from_data_writes_canonical_with_jfamily_override(self):
        from tcrsift.assemble import (
            HUMAN_TRAC_AA,
            HUMAN_TRBC1_AA,
            assemble_full_sequences,
        )

        vdj_a = "CASS" + "A" * 100 + "EQFF"
        vdj_b = "CASS" + "G" * 100 + "EYFF"
        df = pd.DataFrame([{
            "CDR3ab": "c1",
            "CDR3_alpha": vdj_a,
            "CDR3_beta": vdj_b,
            "VDJ_alpha_aa": vdj_a,
            "VDJ_beta_aa": vdj_b,
            "alpha_c_gene": "TRAC",
            "beta_c_gene": "TRBC2",   # CellRanger
            "beta_j_gene": "TRBJ1-1", # J says TRBC1
            "alpha_constant_aa": HUMAN_TRAC_AA,
            "beta_constant_aa": HUMAN_TRBC1_AA,
            "samples": "S1",
        }])
        out = assemble_full_sequences(
            df, alpha_leader=None, beta_leader=None,
            constant_source="from-data",
            verbose=False, show_progress=False,
        )
        # The canonical column was written even though we took the
        # AA from the row, AND it reflects the J-family override.
        assert out["beta_c_gene_canonical"].iloc[0] == "TRBC1"
        assert out["alpha_c_gene_canonical"].iloc[0] == "TRAC"


class TestFixJCParity:
    """#89: autocorrect β J→C parity mismatches in-place via
    :func:`fix_jc_parity`, and via ``validate_sequences(fix=True)``."""

    def test_trbj1_with_trbc2_autocorrected(self):
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene_canonical"] = "TRBC2"
        row["beta_c_gene"] = "TRBC2"
        df = pd.DataFrame([row])
        messages = fix_jc_parity(df)
        assert df.loc[df.index[0], "beta_c_gene_canonical"] == "TRBC1"
        assert any("autocorrected" in m and "TRBC1" in m for m in messages)

    def test_trbj2_with_trbc1_autocorrected(self):
        row = _ok_clone(c_beta_canonical="TRBC2", beta_j_gene="TRBJ2-5")
        row["beta_c_gene_canonical"] = "TRBC1"
        row["beta_c_gene"] = "TRBC1"
        df = pd.DataFrame([row])
        messages = fix_jc_parity(df)
        assert df.loc[df.index[0], "beta_c_gene_canonical"] == "TRBC2"
        assert any("autocorrected" in m and "TRBC2" in m for m in messages)

    def test_fix_jc_parity_with_no_raw_beta_c_gene_column(self):
        """A frame produced by a path that wrote only
        `beta_c_gene_canonical` (no raw `beta_c_gene`) should still
        autocorrect when the canonical disagrees with J."""
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene_canonical"] = "TRBC2"
        # Explicitly NOT setting beta_c_gene.
        row.pop("beta_c_gene", None)
        df = pd.DataFrame([row])
        assert "beta_c_gene" not in df.columns
        messages = fix_jc_parity(df)
        assert df.loc[df.index[0], "beta_c_gene_canonical"] == "TRBC1"
        assert any("autocorrected" in m for m in messages)

    def test_no_correction_when_parity_already_consistent(self):
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-2")
        df = pd.DataFrame([row])
        messages = fix_jc_parity(df)
        assert messages == []
        assert df.loc[df.index[0], "beta_c_gene_canonical"] == "TRBC1"

    def test_validate_with_fix_clears_parity_mismatch(self):
        """``validate_sequences(fix=True)`` should run the autocorrect
        first so the J/C parity check no longer fires."""
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene_canonical"] = "TRBC2"
        row["beta_c_gene"] = "TRBC2"
        df = pd.DataFrame([row])
        warnings = validate_sequences(df, fix=True)
        assert not any("J→C parity mismatch" in w for w in warnings)
        # And the autocorrect note is reported.
        assert any("autocorrected" in w for w in warnings)

    def test_fix_jc_parity_recovers_when_canonical_is_nan(self):
        """Regression: ``float('nan')`` is truthy in Python, so the
        old ``row.get(canonical) or row.get(raw)`` short-circuited to
        NaN when canonical was NaN — fallback to raw never fired,
        autocorrect silently no-op'd."""
        import numpy as np

        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene"] = "TRBC2"
        row["beta_c_gene_canonical"] = np.nan  # the failure case
        df = pd.DataFrame([row])
        messages = fix_jc_parity(df)
        assert df.loc[df.index[0], "beta_c_gene_canonical"] == "TRBC1"
        assert any("autocorrected" in m for m in messages)

    def test_validate_jc_parity_fires_when_canonical_is_nan(self):
        """The parity check inside ``validate_sequences`` also has to
        fall back to raw ``beta_c_gene`` when canonical is NaN."""
        import numpy as np

        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene"] = "TRBC2"  # disagrees with TRBJ1
        row["beta_c_gene_canonical"] = np.nan
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("J→C parity mismatch" in w for w in warnings)

    def test_validate_with_fix_does_not_mask_other_failures(self):
        """Autocorrect only fixes J/C parity; other load-bearing
        failures (e.g. CDR3 not in full chain) must still fire."""
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene_canonical"] = "TRBC2"
        row["beta_c_gene"] = "TRBC2"
        row["CDR3_alpha"] = "NOT_IN_SEQUENCE_XYZ"
        df = pd.DataFrame([row])
        warnings = validate_sequences(df, fix=True)
        assert not any("J→C parity mismatch" in w for w in warnings)
        assert any("CDR3_alpha" in w and "not found" in w for w in warnings)


class TestValidationMessageStructure:
    """validate_sequences returns ValidationMessage(str) instances
    carrying .idx and .severity, so callers don't have to parse
    'Clone {idx}:' back out of the text."""

    def test_is_str_subclass_preserves_text(self):
        row = _ok_clone()
        df = pd.DataFrame([row])
        msgs = validate_sequences(df)
        for m in msgs:
            assert isinstance(m, ValidationMessage)
            assert isinstance(m, str)

    def test_idx_attribute_matches_row_index(self):
        """The idx attribute should be the DataFrame index label of
        the failing row, not a parsed-out string."""
        row = _ok_clone()
        row["CDR3_alpha"] = "NOT_IN_SEQ_XYZ"  # load-bearing failure
        df = pd.DataFrame([row], index=["custom-clone-id"])
        msgs = validate_sequences(df)
        load_bearing = [m for m in msgs if m.severity == "load_bearing"]
        assert load_bearing
        assert all(m.idx == "custom-clone-id" for m in load_bearing)

    def test_metadata_survives_pickle_and_copy(self):
        """`copy.copy` and `pickle.loads(pickle.dumps(...))` preserve
        idx/severity even though plain str operations don't.
        Documented behavior — locked in so a refactor to `__slots__`
        or `__reduce__` doesn't silently break it."""
        import copy
        import pickle

        m = ValidationMessage("Clone 7: chain too short", idx=7, severity="load_bearing")
        m2 = copy.copy(m)
        assert m2.idx == 7 and m2.severity == "load_bearing"
        m3 = pickle.loads(pickle.dumps(m))
        assert m3.idx == 7 and m3.severity == "load_bearing"

    def test_severity_classifies_messages(self):
        """Failures are load_bearing; 'unverifiable' notes are
        informational; autocorrect-from-fix is autocorrect."""
        row = _ok_clone(c_beta_canonical="TRBC1", beta_j_gene="TRBJ1-1")
        row["beta_c_gene_canonical"] = "TRBC2"
        row["beta_c_gene"] = "TRBC2"
        # Also drop both c_gene_canonical columns on alpha to trigger
        # the informational "unverifiable" note.
        row["alpha_c_gene"] = "MYSTERY"
        row["alpha_c_gene_canonical"] = "MYSTERY"
        df = pd.DataFrame([row])
        msgs = validate_sequences(df, fix=True)
        severities = {m.severity for m in msgs}
        assert "autocorrect" in severities
        assert "informational" in severities
        # The parity issue should be cleared by the autocorrect.
        assert not any(
            m.severity == "load_bearing" and "parity mismatch" in m
            for m in msgs
        )


class TestByteForByteEquality:
    """`full == leader + vdj + constant` exactly. The CDR3-substring
    check alone misses cases where the assembler drops/adds residues
    elsewhere."""

    def test_dropped_residue_caught(self):
        row = _ok_clone()
        # full_alpha_aa is one residue shorter than leader+vdj+const.
        full = row["alpha_leader_aa"] + row["vdj_alpha_aa"] + row["alpha_constant_aa"]
        row["full_alpha_aa"] = full[:-1]  # drop terminal S — also breaks C-term
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "full_alpha_aa != leader+vdj+constant" in w for w in warnings
        )

    def test_inserted_residue_caught(self):
        row = _ok_clone()
        full = row["full_alpha_aa"]
        # Insert an extra residue between leader and VDJ.
        leader_len = len(row["alpha_leader_aa"])
        row["full_alpha_aa"] = full[:leader_len] + "X" + full[leader_len:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "full_alpha_aa != leader+vdj+constant" in w for w in warnings
        )

    def test_exact_match_passes(self):
        row = _ok_clone()
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("!= leader+vdj+constant" in w for w in warnings)

    def test_skipped_when_parts_missing(self):
        """If any of leader/vdj/constant is missing, the byte-for-byte
        check can't run — must not falsely fire."""
        row = _ok_clone()
        del row["alpha_leader_aa"]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        # The byte-for-byte check on alpha is skipped; other checks
        # may still produce warnings but not this one.
        assert not any(
            "full_alpha_aa != leader+vdj+constant" in w for w in warnings
        )


class TestPrematureStop:
    """A `*` mid-chain is a clear frame error."""

    def test_premature_stop_caught(self):
        row = _ok_clone()
        # Insert a stop codon halfway through the alpha chain.
        full = row["full_alpha_aa"]
        midpoint = len(full) // 2
        row["full_alpha_aa"] = full[:midpoint] + "*" + full[midpoint:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("premature stop" in w for w in warnings)

    def test_trailing_stop_allowed(self):
        """A single trailing `*` from the stop codon is fine."""
        row = _ok_clone()
        row["full_alpha_aa"] = row["full_alpha_aa"] + "*"
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("premature stop" in w for w in warnings)


class TestMethionineStart:
    """A full chain with a leader should start with M (canonical
    translation start). Catches off-by-one frame errors in the leader."""

    def test_non_methionine_start_caught(self):
        row = _ok_clone()
        full = row["full_alpha_aa"]
        # Replace the leading M with R — clear assembly bug if a
        # leader was supplied.
        row["full_alpha_aa"] = "R" + full[1:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("doesn't start with M" in w for w in warnings)

    def test_no_leader_means_no_methionine_check(self):
        """If no leader column is present, M-start is not enforced —
        the assembly may have skipped the leader on purpose."""
        row = _ok_clone()
        del row["alpha_leader_aa"]
        full = row["full_alpha_aa"]
        row["full_alpha_aa"] = "R" + full[1:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("doesn't start with M" in w for w in warnings)


class TestStandardAlphabet:
    """Only ACDEFGHIKLMNPQRSTVWYX* are valid residues; anything else
    is corruption (e.g. lowercase letters, digits, punctuation)."""

    def test_lowercase_caught(self):
        row = _ok_clone()
        row["full_alpha_aa"] = row["full_alpha_aa"].replace("M", "m", 1)
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("invalid residues" in w for w in warnings)

    def test_digit_caught(self):
        row = _ok_clone()
        full = row["full_alpha_aa"]
        row["full_alpha_aa"] = full[:50] + "3" + full[51:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("invalid residues" in w for w in warnings)

    def test_X_is_allowed(self):
        """X (unknown residue) is valid output from translation when
        an NT triplet contains an N — should NOT trigger this check."""
        row = _ok_clone()
        full = row["full_alpha_aa"]
        row["full_alpha_aa"] = full[:50] + "X" + full[51:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("invalid residues" in w for w in warnings)


class TestSingleChainIntegrity:
    """`single_chain_aa` should contain the linker exactly once
    (β-linker-α construct order)."""

    # Canonical T2A AA (see LINKERS["T2A"]["aa"] in assemble.py).
    _CANONICAL_T2A = "EGRGSLLTCGDVEENPGP"

    def test_missing_linker_caught(self):
        row = _ok_clone()
        row["linker"] = self._CANONICAL_T2A
        row["single_chain_aa"] = row["full_beta_aa"][:-1] + row["full_alpha_aa"]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("missing linker" in w for w in warnings)

    def test_duplicated_linker_caught(self):
        row = _ok_clone()
        linker = self._CANONICAL_T2A
        row["linker"] = linker
        row["single_chain_aa"] = (
            row["full_beta_aa"][:-1] + linker + row["full_alpha_aa"] + linker
        )
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("linker" in w and "2 times" in w for w in warnings)

    def test_correct_construct_passes(self):
        row = _ok_clone()
        linker = self._CANONICAL_T2A
        row["linker"] = linker
        row["single_chain_aa"] = (
            row["full_beta_aa"].rstrip("*") + linker + row["full_alpha_aa"]
        )
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert not any("linker" in w for w in warnings)


class TestAssembleQCReport:
    """The aggregate `[assemble QC]` summary string that #67 specifies."""

    def test_clean_input_reports_pass(self):
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "[assemble QC] 1 / 1 rows checked" in report
        assert "PASS" in report
        assert "FAIL" not in report

    def test_truncated_constant_reports_fail(self):
        """The exact pre-#66 scenario — every row has an 11 aa alpha
        constant — should produce the FAIL banner the issue described."""
        row = _ok_clone()
        row["alpha_constant_aa"] = "RT"  # truncation
        row["full_alpha_aa"] = row["alpha_leader_aa"] + row["vdj_alpha_aa"] + "RT"
        df = pd.DataFrame([row])
        report = assemble_qc_report(df)
        assert "FAIL" in report
        assert "PASS" not in report
        # The constant-floor line should reflect 0/1 pass.
        assert "α constant" in report
        assert "0/1 pass" in report

    def test_terminal_residue_check_in_report(self):
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "terminal residue" in report

    def test_premature_stop_check_in_report(self):
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "no premature stop" in report

    def test_standard_alphabet_check_in_report(self):
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "standard alphabet" in report

    def test_beta_jc_parity_check_in_report(self):
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "J→C in-cis parity" in report

    def test_median_length_appears_in_summary(self):
        """The issue calls out median length in the failure example —
        verify it's surfaced."""
        df = pd.DataFrame([_ok_clone()])
        report = assemble_qc_report(df)
        assert "median" in report

    def test_skips_missing_columns(self):
        """An assembly without leaders shouldn't crash the report."""
        row = _ok_clone()
        for col in ("alpha_leader_aa", "beta_leader_aa"):
            del row[col]
        df = pd.DataFrame([row])
        report = assemble_qc_report(df)  # should not raise
        assert "[assemble QC]" in report

    def test_empty_frame_no_crash(self):
        """Calling on an empty frame is a no-op, not an exception."""
        df = pd.DataFrame({"full_alpha_aa": [], "full_beta_aa": []})
        report = assemble_qc_report(df)
        assert "[assemble QC] 0 / 0" in report


class TestAssemblyQCReportDataclass:
    """The structured :class:`AssemblyQCReport` (#67) replaces the
    string-only return for programmatic consumption."""

    def test_build_returns_structured_report(self):
        from tcrsift.assemble import (
            AssemblyQCCheck,
            AssemblyQCReport,
            build_assembly_qc_report,
        )

        df = pd.DataFrame([_ok_clone()])
        report = build_assembly_qc_report(df)
        assert isinstance(report, AssemblyQCReport)
        assert report.n_rows == 1
        assert report.passed is True
        assert all(isinstance(c, AssemblyQCCheck) for c in report.checks)
        assert len(report.checks) > 0

    def test_check_attributes_make_sense(self):
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([_ok_clone()])
        report = build_assembly_qc_report(df)
        for c in report.checks:
            assert 0 <= c.passed <= c.total
            assert c.failed == c.total - c.passed
            assert c.pass_rate == (1.0 if c.total == 0 else c.passed / c.total)
            assert c.is_passing == (c.failed == 0)

    def test_failing_check_recorded(self):
        from tcrsift.assemble import build_assembly_qc_report

        row = _ok_clone()
        row["alpha_constant_aa"] = "RT"
        row["full_alpha_aa"] = row["alpha_leader_aa"] + row["vdj_alpha_aa"] + "RT"
        df = pd.DataFrame([row])
        report = build_assembly_qc_report(df)
        assert report.passed is False
        bad = [c for c in report.checks if not c.is_passing]
        assert len(bad) >= 1
        # Constant floor check should be one of the failing ones.
        assert any(c.name == "alpha_constant_floor" for c in bad)

    def test_to_dataframe_has_expected_columns(self):
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([_ok_clone()])
        report = build_assembly_qc_report(df)
        rdf = report.to_dataframe()
        assert set(rdf.columns) == {
            "name", "label", "chain", "passed", "failed", "total",
            "pass_rate", "median", "unit",
        }
        assert len(rdf) == len(report.checks)

    def test_format_text_matches_string_helper(self):
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([_ok_clone()])
        assert build_assembly_qc_report(df).format_text() == assemble_qc_report(df)

    def test_str_calls_format_text(self):
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([_ok_clone()])
        report = build_assembly_qc_report(df)
        assert str(report) == report.format_text()


class TestT2ASingleChain:
    """Extra single-chain (β-linker-α) checks (#67):
    - byte-for-byte construct equality,
    - 2A peptide canonical AA when the linker matches a known name."""

    def _sc_clone(self, linker_aa: str = "EGRGSLLTCGDVEENPGP") -> dict:
        """β + linker + α, fully consistent."""
        row = _ok_clone()
        row["linker"] = linker_aa
        # The assembler strips trailing `*` from β before joining.
        row["single_chain_aa"] = (
            row["full_beta_aa"].rstrip("*") + linker_aa + row["full_alpha_aa"]
        )
        return row

    def test_byte_exact_construct_passes(self):
        df = pd.DataFrame([self._sc_clone()])
        warnings = validate_sequences(df)
        assert not any("single_chain_aa !=" in w for w in warnings)

    def test_byte_exact_construct_mismatch_caught(self):
        row = self._sc_clone()
        # Drop one residue between linker and α.
        sc = row["single_chain_aa"]
        linker_end = sc.index(row["linker"]) + len(row["linker"])
        row["single_chain_aa"] = sc[:linker_end] + sc[linker_end + 1:]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("single_chain_aa !=" in w for w in warnings)

    def test_t2a_canonical_passes(self):
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([self._sc_clone("EGRGSLLTCGDVEENPGP")])  # canonical T2A
        report = build_assembly_qc_report(df)
        canon = next(
            (c for c in report.checks if c.name == "single_chain_2a_canonical"),
            None,
        )
        assert canon is not None
        assert canon.is_passing

    def test_non_2a_custom_linker_accepted(self):
        """A custom non-2A linker (e.g. flexible GS) should NOT be
        flagged as non-canonical — the 2A canonical check only fires
        on linkers that look like 2A peptides (NPGP tail)."""
        from tcrsift.assemble import build_assembly_qc_report

        df = pd.DataFrame([self._sc_clone("GGGGSGGGGSGGGGS")])
        report = build_assembly_qc_report(df)
        canon = next(
            (c for c in report.checks if c.name == "single_chain_2a_canonical"),
            None,
        )
        assert canon is not None
        assert canon.is_passing


class TestPlotAssemblyQC:
    """The new :func:`plot_assembly_qc` visualization (#67)."""

    def test_returns_matplotlib_figure(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from tcrsift.assemble import build_assembly_qc_report
        from tcrsift.plots import plot_assembly_qc

        df = pd.DataFrame([_ok_clone()])
        report = build_assembly_qc_report(df)
        fig = plot_assembly_qc(report)
        assert fig is not None
        # Title should reflect PASS state.
        title = fig.axes[0].get_title()
        assert "PASS" in title and "1 clones" in title

    def test_fail_state_in_title(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from tcrsift.assemble import build_assembly_qc_report
        from tcrsift.plots import plot_assembly_qc

        row = _ok_clone()
        row["alpha_constant_aa"] = "RT"
        row["full_alpha_aa"] = row["alpha_leader_aa"] + row["vdj_alpha_aa"] + "RT"
        df = pd.DataFrame([row])
        fig = plot_assembly_qc(build_assembly_qc_report(df))
        title = fig.axes[0].get_title()
        assert "FAIL" in title

    def test_writes_file_when_path_given(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        from tcrsift.assemble import build_assembly_qc_report
        from tcrsift.plots import plot_assembly_qc

        df = pd.DataFrame([_ok_clone()])
        out = tmp_path / "qc.png"
        plot_assembly_qc(build_assembly_qc_report(df), out)
        assert out.exists() and out.stat().st_size > 0

    def test_empty_report_renders_placeholder(self):
        import matplotlib

        matplotlib.use("Agg")
        from tcrsift.assemble import build_assembly_qc_report
        from tcrsift.plots import plot_assembly_qc

        df = pd.DataFrame({"full_alpha_aa": [], "full_beta_aa": []})
        report = build_assembly_qc_report(df)
        fig = plot_assembly_qc(report)  # should not raise
        assert fig is not None
