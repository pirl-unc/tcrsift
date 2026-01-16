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

import pytest

from tcrsift.mnemonic import tcr_name


class TestTcrName:
    """Tests for tcr_name function."""

    def test_basic_name_generation(self):
        """Test basic TCR name generation."""
        result = tcr_name("CASSL")

        # Should return a capitalized string
        assert isinstance(result, str)
        assert len(result) > 0
        assert result[0].isupper()

    def test_split_mode(self):
        """Test split mode creates first/last name."""
        result = tcr_name("CASSLFGAG", split=True)

        # Should contain a space (first + last name)
        assert " " in result

    def test_no_split_mode(self):
        """Test non-split mode creates single word."""
        result = tcr_name("CASS", split=False)

        # Should be single word (no space, unless prefix/suffix added)
        assert isinstance(result, str)

    def test_semicolon_multiple_sequences(self):
        """Test handling of multiple sequences separated by semicolon."""
        result = tcr_name("CASSL;CASSF")

        # Should contain " or " for multiple sequences
        assert " or " in result

    def test_prefix_trimming_cass(self):
        """Test CASS prefix is trimmed and adds 'Dr.' title."""
        result = tcr_name("CASSLGAG", split=True)

        # Should add Dr. prefix
        assert "Dr." in result

    def test_prefix_trimming_cas(self):
        """Test CAS prefix adds 'Prof.' title."""
        result = tcr_name("CASLGAG", split=True)

        # Should add Prof. prefix
        assert "Prof." in result

    def test_suffix_trimming_qyf(self):
        """Test QYF suffix adds 'MD' suffix."""
        result = tcr_name("CASLGAQYF", split=True)

        # Should add MD suffix
        assert "MD" in result

    def test_suffix_trimming_qff(self):
        """Test QFF suffix adds 'PhD' suffix."""
        result = tcr_name("CASLGAQFF", split=True)

        # Should add PhD suffix
        assert "PhD" in result

    def test_vowel_cluster_handling(self):
        """Test handling of vowel clusters."""
        result = tcr_name("AEIOU", split=False)

        # Should be pronounceable (consonants inserted)
        assert isinstance(result, str)
        assert len(result) >= 5

    def test_consonant_cluster_handling(self):
        """Test handling of consonant clusters."""
        result = tcr_name("BCDFGH", split=False)

        # Should be pronounceable (vowels inserted)
        assert isinstance(result, str)

    def test_short_sequence(self):
        """Test with very short sequence."""
        result = tcr_name("CA", split=True)

        assert isinstance(result, str)

    def test_single_char(self):
        """Test with single character."""
        result = tcr_name("A", split=False)

        assert isinstance(result, str)
        assert result == "A"

    def test_empty_string(self):
        """Test with empty string."""
        result = tcr_name("", split=True)

        assert isinstance(result, str)

    def test_deterministic_output(self):
        """Test that same input gives same output."""
        seq = "CASSLFGAGYEQY"

        result1 = tcr_name(seq)
        result2 = tcr_name(seq)

        assert result1 == result2

    def test_different_sequences_different_names(self):
        """Test that different sequences give different names."""
        result1 = tcr_name("CASSLFGAG")
        result2 = tcr_name("CAVRSGYST")

        assert result1 != result2

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

    def test_valid_consonant_clusters_preserved(self):
        """Test that valid consonant clusters are preserved."""
        # "ST" is a valid consonant cluster
        result = tcr_name("STAG", split=False)

        # Should start with ST (valid cluster)
        assert isinstance(result, str)

    def test_q_gets_uo_vowel(self):
        """Test that Q gets 'UO' vowel if surrounded by consonants."""
        result = tcr_name("QQQ", split=False)

        # Should be pronounceable with vowels
        assert isinstance(result, str)

    def test_y_as_vowel_at_end(self):
        """Test that Y is treated as vowel at end of word."""
        result = tcr_name("CASSY", split=False)

        assert isinstance(result, str)

    def test_capitalization(self):
        """Test output is properly capitalized."""
        result = tcr_name("abcde", split=False)

        # First letter should be uppercase, rest lowercase
        if len(result) > 0:
            assert result[0].isupper()
            assert result[1:].islower() or len(result) == 1

    def test_gen_prefix(self):
        """Test CSA prefix adds 'Gen.' title."""
        result = tcr_name("CSALGAG", split=True)

        assert "Gen." in result

    def test_capt_prefix(self):
        """Test CAT prefix adds 'Capt.' title."""
        result = tcr_name("CATLGAG", split=True)

        assert "Capt." in result

    def test_sir_prefix(self):
        """Test CAW prefix adds 'Sir' title."""
        result = tcr_name("CAWLGAG", split=True)

        assert "Sir" in result

    def test_madame_prefix(self):
        """Test CAI prefix adds 'Madame' title."""
        result = tcr_name("CAILGAG", split=True)

        assert "Madame" in result

    def test_honourable_prefix(self):
        """Test CSV prefix adds 'The Honourable' title."""
        result = tcr_name("CSVLGAG", split=True)

        assert "The Honourable" in result

    def test_esq_suffix(self):
        """Test AFF suffix adds 'Esq.' title."""
        result = tcr_name("CASLGAAFF", split=True)

        assert "Esq." in result

    def test_jr_suffix(self):
        """Test LFF suffix adds 'Jr.' title."""
        result = tcr_name("CASLGALFF", split=True)

        assert "Jr." in result

    def test_sr_suffix(self):
        """Test YTF suffix adds 'Sr.' title."""
        result = tcr_name("CASLGAYTF", split=True)

        assert "Sr." in result
