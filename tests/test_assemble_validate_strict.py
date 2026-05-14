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

"""Tests for the strengthened ``validate_sequences`` (#68).

The pre-#68 validator silently passed when ``c_gene`` was missing —
that hole is what let the #66 truncated-constant bug ship undetected.
These tests lock the new behavior:

- Canonical-ending check fires per-gene with allele-suffix handling.
- Canonical-START check fires when the constant doesn't begin at the
  right canonical residues (the failure mode the old check couldn't
  see — the constant could end correctly but the splice could be
  totally wrong upstream).
- Missing ``c_gene`` produces a "didn't check" informational note —
  not a silent pass.
- Per-row ``qc_warnings`` populated by the assembler surface through.
- Strict mode raises on load-bearing failures but accepts
  informational notes.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.assemble import (
    CONSTANT_REGION_ENDINGS,
    HUMAN_CONSTANT_REGIONS_AA,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    validate_sequences,
)
from tcrsift.validation import TCRsiftValidationError


def _ok_clone(
    vdj_alpha: str = "CASS" + "A" * 60 + "VDJALPHA",     # ~72 aa
    vdj_beta: str = "CASS" + "B" * 30 + "VDJBETA",        # ~41 aa
    leader_aa: str = "MALELPLLLL" * 2,                    # 20 aa
    c_gene_alpha: str = "TRAC",
    c_gene_beta: str = "TRBC1",
) -> dict:
    """Build a single clonotype row whose full sequences pass every check.

    Uses the canonical TRAC / TRBC1 AA strings so the canonical-start
    and canonical-end checks both succeed.
    """
    return {
        "CDR3_alpha": vdj_alpha,
        "CDR3_beta": vdj_beta,
        "vdj_alpha_aa": vdj_alpha,
        "vdj_beta_aa": vdj_beta,
        "alpha_c_gene": c_gene_alpha,
        "beta_c_gene": c_gene_beta,
        "alpha_c_gene_canonical": c_gene_alpha,
        "beta_c_gene_canonical": c_gene_beta,
        "full_alpha_aa": leader_aa + vdj_alpha + HUMAN_TRAC_AA,
        "full_beta_aa": leader_aa + vdj_beta + HUMAN_TRBC1_AA,
    }


class TestValidateSequencesStrict:
    """Strict mode is the gate the PDF and the CLI ``assemble`` step
    need so structurally-broken TCRs don't reach downstream."""

    def test_clean_input_returns_empty_list(self):
        df = pd.DataFrame([_ok_clone()])
        warnings = validate_sequences(df, strict=False)
        assert warnings == []

    def test_strict_passes_on_clean_input(self):
        df = pd.DataFrame([_ok_clone()])
        assert validate_sequences(df, strict=True) == []

    def test_strict_raises_on_load_bearing_failure(self):
        """A row whose CDR3 isn't in the full sequence is load-bearing —
        cloning constructs would be wrong."""
        row = _ok_clone()
        row["full_alpha_aa"] = "GIBBERISHTHATDOESNOTCONTAINTHECDR3" + HUMAN_TRAC_AA
        df = pd.DataFrame([row])
        with pytest.raises(TCRsiftValidationError, match="validation failed"):
            validate_sequences(df, strict=True)


class TestCanonicalCTerminusCheck:
    """The pre-#68 check — the only existing test for the C-terminus —
    silently passed when c_gene was missing. Now we distinguish."""

    def test_truncated_constant_caught(self):
        """The #66 scenario: constant region truncated. The full
        sequence doesn't end with the canonical C-terminus."""
        row = _ok_clone()
        # Replace canonical TRAC tail with a stub (the #66 bug
        # produced 2-11 aa constants).
        row["full_alpha_aa"] = "M" * 20 + row["vdj_alpha_aa"] + "RT"
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "doesn't end with canonical TRAC C-terminus" in w for w in warnings
        )

    def test_allele_suffix_handled(self):
        """``TRBC2*01`` should still resolve to the TRBC2 canonical
        check (#68) — the old check did exact-string match against
        the dict, which missed the common allele-suffixed form."""
        row = _ok_clone(c_gene_beta="TRBC2*01")
        row["beta_c_gene_canonical"] = "TRBC2"  # set by _add_constant_regions
        df = pd.DataFrame([row])
        # The OK clone uses TRBC1 in the full sequence — with
        # c_gene_canonical=TRBC2 it now fails the C-terminus check.
        warnings = validate_sequences(df)
        assert any("TRBC2 C-terminus" in w for w in warnings)


class TestCanonicalStartCheck:
    """New in #68: the constant region must start at the right
    canonical residues, not just end with them. Catches splice errors
    that would have ended correctly by accident."""

    def test_wrong_start_caught(self):
        """A full sequence that ends with the canonical TRAC tail but
        starts with the wrong post-VDJ residues — what you'd get if
        the canonical splice picked the wrong gene or had an offset
        bug — is rejected."""
        row = _ok_clone()
        # End is canonical TRAC, but the chunk between VDJ and the
        # last 10 aa is replaced with garbage; the canonical-START
        # check should fire.
        row["full_alpha_aa"] = (
            "M" * 20
            + row["vdj_alpha_aa"]
            + "GARBAGESTARTXXX"
            + HUMAN_TRAC_AA[15:]  # tail is right, start is wrong
        )
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any(
            "doesn't match canonical TRAC" in w
            and "start" in w
            for w in warnings
        )


class TestCGeneUnverifiable:
    """The hole that hid #66: missing c_gene → silent pass. Fixed."""

    def test_missing_c_gene_emits_didnt_check_note(self):
        row = _ok_clone()
        row["alpha_c_gene"] = ""
        row["alpha_c_gene_canonical"] = ""
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        # The check should run on TRBC1 and emit a "didn't check"
        # note for alpha.
        assert any(
            "alpha c_gene" in w and "unverifiable" in w for w in warnings
        )

    def test_strict_does_not_raise_on_didnt_check_notes(self):
        """Informational notes are not load-bearing — strict mode
        should still pass when only the "didn't check" notes are
        present. The clone is fully valid except its c_gene call is
        absent, so we can verify TRBC1 (β) but not TRAC (α)."""
        row = _ok_clone()
        row["alpha_c_gene"] = ""
        row["alpha_c_gene_canonical"] = ""
        df = pd.DataFrame([row])
        # Should return the didn't-check note without raising.
        notes = validate_sequences(df, strict=True)
        assert any("unverifiable" in n for n in notes)


class TestQCWarningsSurfacing:
    """When ``_add_constant_regions`` detects a contig-vs-canonical
    start mismatch during assembly, it stashes a per-row
    ``qc_warnings`` list. ``validate_sequences`` must surface those
    so downstream gates see them (#68)."""

    def test_assembler_qc_warnings_become_validation_warnings(self):
        row = _ok_clone()
        row["qc_warnings"] = [
            "alpha constant start mismatch: observed 'XYZ' differs from canonical TRAC"
        ]
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("constant start mismatch" in w for w in warnings)

    def test_strict_raises_on_assembler_qc_warnings(self):
        row = _ok_clone()
        row["qc_warnings"] = ["alpha constant start mismatch: BAD"]
        df = pd.DataFrame([row])
        with pytest.raises(TCRsiftValidationError):
            validate_sequences(df, strict=True)


class TestLengthChecks:
    """The length window (200-450 aa) is now load-bearing."""

    def test_too_short_is_load_bearing(self):
        row = _ok_clone()
        row["full_alpha_aa"] = "M" * 50  # well under 200
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("too short" in w for w in warnings)
        with pytest.raises(TCRsiftValidationError):
            validate_sequences(df, strict=True)

    def test_too_long_is_load_bearing(self):
        row = _ok_clone()
        row["full_alpha_aa"] = "M" * 500
        df = pd.DataFrame([row])
        warnings = validate_sequences(df)
        assert any("too long" in w for w in warnings)
