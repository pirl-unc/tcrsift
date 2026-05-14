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

"""Tests for `_default_tcr_sequence_columns` (#65 — column order /
duplicates fix in the TCR sequence PDF auto-detect)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.plots import _default_tcr_sequence_columns  # noqa: E402


def test_columns_in_biological_assembly_order():
    """β leader → β VDJ → β constant → linker → α leader → α VDJ →
    α constant. Previously the α leader landed between the β leader
    and β VDJ, which is biologically meaningless (#65)."""
    df = pd.DataFrame(
        {
            "alpha_leader_aa": [""],
            "vdj_alpha_aa": [""],
            "alpha_constant_aa": [""],
            "linker": [""],
            "beta_leader_aa": [""],
            "vdj_beta_aa": [""],
            "beta_constant_aa": [""],
        }
    )
    cols = list(_default_tcr_sequence_columns(df))
    assert cols == [
        "beta_leader_aa",
        "vdj_beta_aa",
        "beta_constant_aa",
        "linker",
        "alpha_leader_aa",
        "vdj_alpha_aa",
        "alpha_constant_aa",
    ]


def test_full_aa_columns_excluded():
    """``full_alpha_aa`` and ``full_beta_aa`` are leader + VDJ + constant
    concatenations. Earlier versions rendered them alongside the parts,
    duplicating each chain. The auto-detect must skip them (#65)."""
    df = pd.DataFrame(
        {
            "beta_leader_aa": [""],
            "vdj_beta_aa": [""],
            "beta_constant_aa": [""],
            "full_beta_aa": [""],
            "alpha_leader_aa": [""],
            "vdj_alpha_aa": [""],
            "alpha_constant_aa": [""],
            "full_alpha_aa": [""],
        }
    )
    cols = list(_default_tcr_sequence_columns(df))
    assert "full_beta_aa" not in cols
    assert "full_alpha_aa" not in cols


def test_uppercase_vdj_alias_no_longer_matches():
    """Earlier versions carried both lowercase ``vdj_*_aa`` and
    uppercase ``VDJ_*_aa`` aliases — when both columns existed the
    chain rendered twice. Only the canonical lowercase form
    ``assemble.py`` writes to ``full_sequences.csv`` should match (#65).
    """
    # Frame with ONLY the uppercase alias present (pre-assemble shape).
    df = pd.DataFrame({"VDJ_beta_aa": [""], "VDJ_alpha_aa": [""]})
    cols = list(_default_tcr_sequence_columns(df))
    # Neither uppercase column should appear in the detected set.
    assert "VDJ_beta_aa" not in cols
    assert "VDJ_alpha_aa" not in cols
    assert cols == []  # nothing else to render


def test_both_cases_present_renders_chain_only_once():
    """When BOTH uppercase and lowercase columns happen to coexist on
    a frame, only the canonical lowercase one renders (no duplication)."""
    df = pd.DataFrame(
        {
            "vdj_beta_aa": [""],
            "VDJ_beta_aa": [""],  # alias — must not double-render
            "vdj_alpha_aa": [""],
            "VDJ_alpha_aa": [""],
        }
    )
    cols = list(_default_tcr_sequence_columns(df))
    assert cols.count("vdj_beta_aa") == 1
    assert cols.count("vdj_alpha_aa") == 1
    assert "VDJ_beta_aa" not in cols
    assert "VDJ_alpha_aa" not in cols


def test_partial_chain_subset():
    """A frame with only the β block columns yields just those, in
    order, with no placeholders or duplicates."""
    df = pd.DataFrame(
        {"beta_leader_aa": [""], "vdj_beta_aa": [""], "beta_constant_aa": [""]}
    )
    cols = list(_default_tcr_sequence_columns(df))
    assert cols == ["beta_leader_aa", "vdj_beta_aa", "beta_constant_aa"]


def test_single_chain_fallback():
    """When no canonical assembled columns are present, fall back to
    ``single_chain_aa`` if it exists."""
    df = pd.DataFrame({"single_chain_aa": ["MABCDE"]})
    cols = _default_tcr_sequence_columns(df)
    assert cols == {"single_chain_aa": "Single Chain"}


def test_empty_input_returns_empty_dict():
    """No matching columns + no single-chain fallback → empty result.
    Caller short-circuits and reports nothing to render."""
    df = pd.DataFrame({"unrelated_column": ["x"]})
    assert _default_tcr_sequence_columns(df) == {}


def test_display_labels_match_chains():
    """The display labels paired with each canonical column should
    name the chain ("Beta Leader", "Alpha VDJ", etc.) so the PDF
    legend isn't ambiguous."""
    df = pd.DataFrame(
        {
            "beta_leader_aa": [""],
            "vdj_beta_aa": [""],
            "beta_constant_aa": [""],
            "alpha_leader_aa": [""],
            "vdj_alpha_aa": [""],
            "alpha_constant_aa": [""],
        }
    )
    detected = _default_tcr_sequence_columns(df)
    assert detected["beta_leader_aa"] == "Beta Leader"
    assert detected["vdj_beta_aa"] == "Beta VDJ"
    assert detected["beta_constant_aa"] == "Beta Constant"
    assert detected["alpha_leader_aa"] == "Alpha Leader"
    assert detected["vdj_alpha_aa"] == "Alpha VDJ"
    assert detected["alpha_constant_aa"] == "Alpha Constant"


class TestCreateTcrSequencePdfStrict:
    """The new sanity gate on PDF rendering (#68). The pre-#68 PDF
    happily rendered nonsense — broken constants, missing CDR3,
    out-of-spec lengths — because nothing checked the input."""

    def test_strict_refuses_to_render_truncated_constant(self, tmp_path):
        """The #66 scenario: truncated constant region. PDF must
        refuse rather than emit a structurally-invalid figure."""
        from tcrsift.assemble import HUMAN_TRAC_AA
        from tcrsift.plots import create_tcr_sequence_pdf
        from tcrsift.validation import TCRsiftValidationError

        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        # The constant is truncated to "RT" — exactly what the #66
        # bug produced — but everything else is in order.
        df = pd.DataFrame(
            [
                {
                    "CDR3_alpha": vdj_alpha,
                    "CDR3_beta": "CASSGETAVDF",
                    "vdj_alpha_aa": vdj_alpha,
                    "vdj_beta_aa": "CASSGETAVDF",
                    "alpha_c_gene_canonical": "TRAC",
                    "beta_c_gene_canonical": "TRBC1",
                    "full_alpha_aa": "M" * 20 + vdj_alpha + "RT",  # truncated
                    "full_beta_aa": "M" * 20 + "CASSGETAVDF" + HUMAN_TRAC_AA,
                    "alpha_constant_aa": "RT",
                    "beta_constant_aa": HUMAN_TRAC_AA,
                }
            ]
        )
        out = tmp_path / "broken.pdf"
        with pytest.raises(TCRsiftValidationError, match="validation failures"):
            create_tcr_sequence_pdf(df, out, strict=True)
        # Nothing written.
        assert not out.exists()

    def test_strict_false_logs_warning_and_renders(self, tmp_path, caplog):
        """``strict=False`` is the escape hatch: render despite
        warnings, but tell the user."""
        import logging

        from tcrsift.assemble import HUMAN_TRAC_AA, HUMAN_TRBC1_AA
        from tcrsift.plots import create_tcr_sequence_pdf

        vdj_alpha = "CASS" + "A" * 60 + "VLPHA"
        df = pd.DataFrame(
            [
                {
                    "CDR3_alpha": vdj_alpha,
                    "CDR3_beta": "CASSGETAVDF",
                    "vdj_alpha_aa": vdj_alpha,
                    "vdj_beta_aa": "CASSGETAVDF",
                    "alpha_c_gene_canonical": "TRAC",
                    "beta_c_gene_canonical": "TRBC1",
                    "full_alpha_aa": "M" * 20 + vdj_alpha + "RT",  # truncated
                    "full_beta_aa": "M" * 20 + "CASSGETAVDF" + HUMAN_TRBC1_AA,
                    "alpha_constant_aa": "RT",
                    "beta_constant_aa": HUMAN_TRBC1_AA,
                }
            ]
        )
        out = tmp_path / "rendered_with_warning.pdf"
        with caplog.at_level(logging.WARNING, logger="tcrsift.plots"):
            try:
                create_tcr_sequence_pdf(df, out, strict=False)
            except ImportError:
                # reportlab not installed in CI — the warning still
                # got logged before the renderer tried to import.
                pass
        assert any(
            "validation failures" in r.message for r in caplog.records
        )
