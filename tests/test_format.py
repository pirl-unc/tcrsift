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

    def test_three_part_name_no_reorder(self):
        # The CTY-second rule is only applied for 2-part names.
        assert pretty_method("AIMpos_IFNpos_CTYneg") == "AIM⁺IFN⁺CTY⁻"

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


class TestTopLevelExport:
    def test_importable_from_top_level(self):
        assert tcrsift.pretty_method("AIMpos") == "AIM⁺"
        assert tcrsift.pretty_sample("tetpos-3") == "tet⁺"
