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

"""Tests for ``tcrsift.format`` (#77) — pretty-printer for sample-sheet
enrichment-method strings."""

from __future__ import annotations

import pandas as pd

import tcrsift
from tcrsift.format import (
    pretty_method,
    pretty_methods,
    pretty_sample,
    pretty_samples,
)


class TestPrettyMethod:
    def test_pos_suffix(self):
        assert pretty_method("AIMpos") == "AIM⁺"
        assert pretty_method("tetpos") == "tet⁺"
        assert pretty_method("IFNpos") == "IFN⁺"

    def test_neg_suffix(self):
        assert pretty_method("CTYneg") == "CTY⁻"

    def test_compound_with_cty_already_second(self):
        # Already in display order — no reordering.
        assert pretty_method("AIMpos_CTYneg") == "AIM⁺CTY⁻"
        assert pretty_method("IFNpos_CTYneg") == "IFN⁺CTY⁻"

    def test_compound_with_cty_first_reorders(self):
        # CTY belongs after the positive marker.
        assert pretty_method("CTYneg_tetpos") == "tet⁺CTY⁻"

    def test_three_part_name_keeps_baseline_last(self):
        assert pretty_method("AIMpos_IFNpos_CTYneg") == "AIM⁺IFN⁺CTY⁻"

    def test_three_part_name_reorders_baseline_last(self):
        # Baseline (CTY) is pushed last regardless of input position — the
        # generalization of the old 2-part rule to N parts (#208).
        assert pretty_method("CTYneg_AIMpos_IFNpos") == "AIM⁺IFN⁺CTY⁻"
        assert pretty_method("IFNpos_CTYneg_AIMpos") == "AIM⁺IFN⁺CTY⁻"

    def test_baseline_override_per_call(self):
        # Make DMSO the baseline instead of CTY for this call.
        assert pretty_method("DMSOneg_AIMpos", last=("DMSO",)) == "AIM⁺DMSO⁻"

    def test_priority_pins_leading_order(self):
        # Without priority, alphabetical: AIM before tet. With priority, tet leads.
        assert pretty_method("AIMpos_tetpos") == "AIM⁺tet⁺"
        assert pretty_method("AIMpos_tetpos", priority=("tet",)) == "tet⁺AIM⁺"

    def test_unknown_suffix_strips_underscore_separator(self):
        # No ``pos``/``neg`` to translate, but the ``_`` separator
        # still drops out (the underscore is a name-encoding artifact,
        # not display content).
        assert pretty_method("custom_label") == "customlabel"

    def test_single_token_unknown_passes_through(self):
        assert pretty_method("custom") == "custom"

    def test_non_string_coerced(self):
        assert pretty_method(None) == "None"
        assert pretty_method(42) == "42"


class TestPrettyMethods:
    def test_list_mapping(self):
        out = pretty_methods(["AIMpos", "CTYneg", "AIMpos_CTYneg"])
        assert out == ["AIM⁺", "CTY⁻", "AIM⁺CTY⁻"]

    def test_pandas_series_mapping(self):
        s = pd.Series(["AIMpos", "tetpos"])
        assert pretty_methods(s) == ["AIM⁺", "tet⁺"]


class TestPrettySample:
    def test_strips_donor_suffix(self):
        # Sample = method + donor number suffix like "-2".
        assert pretty_sample("AIMpos_CTYneg-2") == "AIM⁺CTY⁻"
        assert pretty_sample("tetpos-1") == "tet⁺"

    def test_no_suffix_passes_through_method_rules(self):
        assert pretty_sample("AIMpos") == "AIM⁺"

    def test_empty_stem_falls_back_to_original(self):
        # "-2" alone strips to "" → return original.
        assert pretty_sample("-2") == "-2"

    def test_non_string_coerced(self):
        assert pretty_sample(None) == "None"


class TestPrettySamples:
    def test_list_mapping(self):
        out = pretty_samples(["AIMpos_CTYneg-2", "tetpos-1", "CTYneg-2"])
        assert out == ["AIM⁺CTY⁻", "tet⁺", "CTY⁻"]


class TestOrderConditions:
    def test_consistent_regardless_of_input_order(self):
        from tcrsift.format import order_conditions

        assert order_conditions(["CTYneg", "AIMpos"]) == ["AIMpos", "CTYneg"]
        assert order_conditions(["AIMpos", "CTYneg"]) == ["AIMpos", "CTYneg"]

    def test_baseline_always_last(self):
        from tcrsift.format import order_conditions

        out = order_conditions(["CTYneg", "tetpos", "AIMpos"])
        assert out[-1] == "CTYneg"

    def test_priority_then_baseline(self):
        from tcrsift.format import order_conditions

        out = order_conditions(["CTYneg", "AIMpos", "tetpos"], priority=("tet",))
        assert out == ["tetpos", "AIMpos", "CTYneg"]

    def test_override_last_marker(self):
        from tcrsift.format import order_conditions

        out = order_conditions(["DMSOneg", "AIMpos"], last=("DMSO",))
        assert out[-1] == "DMSOneg"


class TestSetBaselineMarkers:
    def test_global_override_and_restore(self):
        from tcrsift import format as fmt

        original = fmt.BASELINE_MARKERS
        try:
            fmt.set_baseline_markers("DMSO")
            assert fmt.pretty_method("DMSOneg_AIMpos") == "AIM⁺DMSO⁻"
            # CTY no longer special → alphabetical (AIM < CTY).
            assert fmt.pretty_method("CTYneg_AIMpos") == "AIM⁺CTY⁻"
        finally:
            fmt.set_baseline_markers(*original)
        assert fmt.BASELINE_MARKERS == original


class TestPdfSafe:
    def test_superscripts_to_ascii(self):
        from tcrsift.format import pdf_safe

        assert pdf_safe("AIM⁺CTY⁻") == "AIM+CTY-"

    def test_greek_chain_labels_spelled_out(self):
        from tcrsift.format import pdf_safe

        assert pdf_safe("β-T2A-α") == "beta-T2A-alpha"

    def test_set_symbols(self):
        from tcrsift.format import pdf_safe

        assert pdf_safe("AIM⁺ ∩ CTY⁻") == "AIM+ & CTY-"

    def test_winansi_glyphs_preserved(self):
        from tcrsift.format import pdf_safe

        # Em-dash + curly quote are in WinAnsi (cp1252) and must survive.
        assert pdf_safe("a — b") == "a — b"

    def test_winansi_multiplication_sign_preserved(self):
        from tcrsift.format import pdf_safe

        # × is in cp1252, so it must NOT be downgraded to 'x' (#202 review).
        assert pdf_safe("freq × PRISM") == "freq × PRISM"

    def test_unknown_glyph_replaced_not_box(self):
        from tcrsift.format import pdf_safe

        # An exotic glyph with no mapping becomes '?' (never a missing-glyph box).
        assert "★" not in pdf_safe("a★b")


class TestPrettyAntigen:
    """#317: antigen-agnostic condition labels — never a bare pool token when
    a label exists; the label is opaque free text (peptide / RNA / protein …)."""

    def test_explicit_label_wins(self):
        from tcrsift.format import pretty_antigen

        assert pretty_antigen("P2", label="P2 (KIF1C)") == "P2 (KIF1C)"

    def test_per_call_map(self):
        from tcrsift.format import pretty_antigen

        assert pretty_antigen("P2", labels={"P2": "P2 (KIF1C)"}) == "P2 (KIF1C)"

    def test_unmapped_falls_back_to_pretty_method(self):
        from tcrsift.format import pretty_antigen

        # No label for P9 → the raw token (pretty_method leaves it unchanged).
        assert pretty_antigen("P9", labels={"P2": "P2 (KIF1C)"}) == "P9"
        # A sort-style token still gets the pretty superscript.
        assert pretty_antigen("AIMpos") == "AIM⁺"

    def test_format_neutral_free_text(self):
        from tcrsift.format import pretty_antigen

        # Works for non-peptide antigens — the label is never parsed.
        for lbl in ("WT1 minigene pool B", "MART-1 26-35", "EBV lysate", "whole PRAME"):
            assert pretty_antigen("cond", label=lbl) == lbl

    def test_global_registry_set_and_clear(self):
        from tcrsift.format import pretty_antigen, set_antigen_labels

        set_antigen_labels({"WT1": "WT1 minigene pool B"})
        try:
            assert pretty_antigen("WT1") == "WT1 minigene pool B"
        finally:
            set_antigen_labels()  # clear so other tests aren't affected
        assert pretty_antigen("WT1") == "WT1"


class TestTopLevelExport:
    def test_importable_from_top_level(self):
        assert tcrsift.pretty_method("AIMpos") == "AIM⁺"
        assert tcrsift.pretty_sample("tetpos-3") == "tet⁺"
        assert tcrsift.order_conditions(["CTYneg", "AIMpos"]) == ["AIMpos", "CTYneg"]

    def test_pretty_antigen_importable_from_top_level(self):
        assert tcrsift.pretty_antigen("P2", label="P2 (KIF1C)") == "P2 (KIF1C)"
