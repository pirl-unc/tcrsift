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

"""Tests for mnemonic TCR naming module."""

from tcrsift.mnemonic import (
    DEFAULT_PREFIXES,
    DEFAULT_SUFFIXES,
    _make_pronounceable,
    _split_into_words,
    tcr_name,
)


class TestMakePronounceable:
    """Tests for the _make_pronounceable helper function."""

    def test_amino_vowels_preserved(self):
        """Amino acid vowels (A, E, I, Y) should be preserved as-is."""
        # Input has amino acid vowels
        result = _make_pronounceable("laga")
        assert result == "laga"  # Original A preserved

    def test_inserted_vowels_use_o_u_or_diphthongs(self):
        """Inserted vowels should use O, U, or diphthongs - not A/E/I alone."""
        result = _make_pronounceable("lg")
        # Should have inserted vowel that's NOT a single a, e, or i
        # (since those are amino acid codes)
        inserted = result[1:-1]  # The inserted part between l and g
        # Either it's o/u, or a diphthong containing o/u
        assert any(c in "ou" for c in inserted) or inserted in {"ai", "ei", "oi", "au", "oo", "ee"}

    def test_distinguishes_lga_from_laga(self):
        """LGA (needs insertion) should differ from LAGA (has original A)."""
        lga = _make_pronounceable("lga")
        laga = _make_pronounceable("laga")
        assert lga != laga
        # LAGA should keep the original 'a'
        assert "a" in laga
        # LGA should have a diphthong or o/u inserted
        assert lga != "lga"  # Something was inserted

    def test_consonant_clusters_get_vowels(self):
        """Invalid consonant clusters should get vowels inserted."""
        result = _make_pronounceable("bcdf")
        # Should have vowels inserted between consonants
        assert len(result) > 4

    def test_valid_clusters_preserved(self):
        """Valid consonant clusters like 'st', 'tr' should be preserved."""
        result = _make_pronounceable("st")
        assert result == "st"

        result = _make_pronounceable("tr")
        assert result == "tr"

    def test_y_treated_as_vowel(self):
        """Y should be treated as a vowel."""
        result = _make_pronounceable("y")
        assert result == "y"

    def test_lowercase_output(self):
        """Output should be lowercase."""
        result = _make_pronounceable("ABCD")
        assert result == result.lower()


class TestSplitIntoWords:
    """Tests for the _split_into_words helper function."""

    def test_short_string_no_split(self):
        """Short strings should not be split."""
        result = _split_into_words("abc")
        assert len(result) == 1
        assert result[0] == "abc"

    def test_long_string_split(self):
        """Longer strings should be split at vowel boundaries."""
        # Use a string with vowels so split logic can find boundaries
        result = _split_into_words("abodefigohu", target_len=4)
        assert len(result) > 1

    def test_returns_list(self):
        """Should return a list of strings."""
        result = _split_into_words("abcdefgh")
        assert isinstance(result, list)
        for word in result:
            assert isinstance(word, str)


class TestTcrName:
    """Tests for tcr_name function."""

    def test_basic_name_generation(self):
        """Test basic TCR name generation."""
        result = tcr_name("CASSLGQAYEQYF")
        assert isinstance(result, str)
        assert len(result) > 0
        assert result[0].isupper()

    def test_empty_string_returns_anon(self):
        """Empty string should return 'Anon'."""
        result = tcr_name("")
        assert result == "Anon"

    def test_prefix_stripping_cass(self):
        """CASS prefix should be stripped."""
        # After stripping CASS and CAS, both should give same variable region
        result_with_cass = tcr_name("CASSLGQAY")
        result_with_cas = tcr_name("CASLGQAY")
        # Both strip to "LGQAY"
        assert result_with_cass == result_with_cas

    def test_prefix_stripping_cav(self):
        """CAV prefix should be stripped (alpha chains)."""
        result = tcr_name("CAVSDGGSF")
        # Should not start with "Cav"
        assert not result.lower().startswith("cav")

    def test_suffix_stripping_f(self):
        """F suffix should be stripped."""
        result_with_f = tcr_name("CASSLGQAYF")
        result_without_f = tcr_name("CASSLGQAY")
        assert result_with_f == result_without_f

    def test_semicolon_multiple_sequences(self):
        """Multiple sequences separated by semicolon should be joined with 'or'."""
        result = tcr_name("CASSLGQAYF;CASSIRASYF")
        assert " or " in result

    def test_deterministic_output(self):
        """Same input should always produce same output."""
        seq = "CASSLFGAGYEQYF"
        result1 = tcr_name(seq)
        result2 = tcr_name(seq)
        assert result1 == result2

    def test_different_sequences_different_names(self):
        """Different sequences should produce different names."""
        result1 = tcr_name("CASSLFGAG")
        result2 = tcr_name("CAVRSGYST")
        assert result1 != result2

    def test_similar_sequences_similar_names(self):
        """Similar sequences should produce similar names."""
        result1 = tcr_name("CASSIRASYEQYF")
        result2 = tcr_name("CASSIRANYEQYF")
        # Names should share common prefix (Iras/Iran)
        assert result1[:3].lower() == result2[:3].lower()

    def test_distinguishes_original_vs_inserted_vowels(self):
        """Should distinguish sequences with original vowels vs inserted ones."""
        # LAGA has original A, LGA needs insertion
        with_a = tcr_name("CASSLAGAF")  # Variable region: LAGA
        without_a = tcr_name("CASSLGAF")  # Variable region: LGA
        assert with_a != without_a

    def test_typical_cdr3_alpha(self):
        """Test with typical CDR3 alpha sequence."""
        result = tcr_name("CAVSDLEPNSSASKIIF")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_typical_cdr3_beta(self):
        """Test with typical CDR3 beta sequence."""
        result = tcr_name("CASSLAPGTGELFF")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_custom_prefixes(self):
        """Custom prefix list should be used."""
        result = tcr_name("XYZABCDE", prefixes=["XYZ"])
        assert not result.lower().startswith("x")

    def test_custom_suffixes(self):
        """Custom suffix list should be used."""
        result_with = tcr_name("ABCDEQQQ", suffixes=["QQQ"])
        result_without = tcr_name("ABCDE", suffixes=["QQQ"])
        assert result_with == result_without

    def test_no_prefix_stripping(self):
        """Empty prefix list should not strip prefixes."""
        result = tcr_name("CASSLGQAY", prefixes=[])
        assert result.lower().startswith("c")

    def test_no_suffix_stripping(self):
        """Empty suffix list should not strip suffixes."""
        result_with_f = tcr_name("CASSLGQAYF", suffixes=[])
        result_without_f = tcr_name("CASSLGQAY", suffixes=[])
        assert result_with_f != result_without_f

    def test_pronounceable_output(self):
        """Output should be pronounceable (no extremely long consonant runs)."""
        result = tcr_name("CASSRRRRRRRF")
        vowels = set("aeiouy")
        consonants_in_row = 0
        for c in result.lower():
            if c == " ":
                consonants_in_row = 0
            elif c in vowels:
                consonants_in_row = 0
            else:
                consonants_in_row += 1
            # Should never have more than 3 consonants in a row
            assert consonants_in_row <= 3

    def test_whitespace_handling(self):
        """Leading/trailing whitespace should be handled."""
        result = tcr_name("  CASSLGQAYF  ")
        assert result == tcr_name("CASSLGQAYF")

    def test_lowercase_input(self):
        """Lowercase input should work the same as uppercase."""
        result_upper = tcr_name("CASSLGQAYF")
        result_lower = tcr_name("casslgqayf")
        assert result_upper == result_lower

    def test_minimal_sequence_after_stripping(self):
        """Sequence that becomes empty after stripping should return Anon."""
        result = tcr_name("CASSF")
        assert result == "Anon"


class TestDefaultConstants:
    """Test that default constants are properly defined."""

    def test_default_prefixes_contains_cass(self):
        """Default prefixes should contain CASS."""
        assert "CASS" in DEFAULT_PREFIXES

    def test_default_prefixes_contains_cav(self):
        """Default prefixes should contain CAV."""
        assert "CAV" in DEFAULT_PREFIXES

    def test_default_suffixes_contains_f(self):
        """Default suffixes should contain F."""
        assert "F" in DEFAULT_SUFFIXES
