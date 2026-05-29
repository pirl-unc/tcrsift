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

"""Tests for #57: Levenshtein-1 fuzzy CDR3 matching in
``match_clonotypes``.

The fuzzy fallback fires only when:
1. ``match_mode="levenshtein"`` is requested explicitly
2. The clone doesn't hit any exact αβ or β-only match first
3. There exists a database β CDR3 within edit distance ≤ max_distance

αβ matching stays strict-exact in all modes — biologically fuzzy αβ
is too noisy (the paired-chain prior is the strong signal).
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.annotate import (
    _build_cdr3_neighbor_index,
    _canonical_variants,
    _find_lev_neighbors,
    _levenshtein_distance_at_most_1,
    match_clonotypes,
)


def _db(rows):
    df = pd.DataFrame(rows)
    if "database" not in df.columns:
        df["database"] = "TEST"
    if "is_viral" not in df.columns:
        df["is_viral"] = False
    if "species" not in df.columns:
        df["species"] = "Homo sapiens"
    if "antigen_gene" not in df.columns:
        df["antigen_gene"] = ""
    if "epitope" not in df.columns:
        df["epitope"] = ""
    return df


class TestLevenshteinDistanceAtMost1:
    """Direct unit tests for the inner distance check."""

    def test_identical_returns_0(self):
        assert _levenshtein_distance_at_most_1("CASS", "CASS") == 0

    def test_single_substitution_returns_1(self):
        assert _levenshtein_distance_at_most_1("CASS", "CXSS") == 1

    def test_substitution_at_end_returns_1(self):
        assert _levenshtein_distance_at_most_1("CASS", "CASR") == 1

    def test_single_deletion_returns_1(self):
        # Delete 'A' at index 1 → CSS
        assert _levenshtein_distance_at_most_1("CASS", "CSS") == 1

    def test_single_insertion_returns_1(self):
        # CSS → CASS (insert A at index 1)
        assert _levenshtein_distance_at_most_1("CSS", "CASS") == 1

    def test_two_substitutions_returns_none(self):
        assert _levenshtein_distance_at_most_1("CASS", "CXSR") is None

    def test_length_diff_more_than_one_returns_none(self):
        assert _levenshtein_distance_at_most_1("CASS", "CASSXY") is None

    def test_substitution_then_deletion_returns_none(self):
        # CASS → CXS (substitute A→X, delete final S) — Lev distance 2
        assert _levenshtein_distance_at_most_1("CASS", "CXS") is None

    def test_empty_strings(self):
        assert _levenshtein_distance_at_most_1("", "") == 0
        assert _levenshtein_distance_at_most_1("", "X") == 1


class TestCanonicalVariants:
    """The deletion-canonical variant generator is the core of the
    O(L) candidate-set lookup. Verify it generates exactly the
    expected variants."""

    def test_includes_self(self):
        assert "CASS" in set(_canonical_variants("CASS"))

    def test_includes_each_deletion(self):
        variants = set(_canonical_variants("CASS"))
        assert "ASS" in variants   # delete idx 0
        assert "CSS" in variants   # delete idx 1
        assert "CAS" in variants   # delete idx 2 OR idx 3 (dedup'd by set)

    def test_count_is_at_most_1_plus_len(self):
        # set dedup-collapses identical deletions ("CAS" twice for CASS).
        variants = set(_canonical_variants("CASS"))
        assert len(variants) <= 1 + len("CASS")

    def test_empty_string_yields_nothing(self):
        assert list(_canonical_variants("")) == []


class TestBuildCdr3NeighborIndex:
    def test_indexes_each_canonical_variant(self):
        series = pd.Series(["CASS", "CAST"])
        index = _build_cdr3_neighbor_index(series)
        # Both CDR3s share deletion-canonical "CAS" (delete final char).
        assert "CAS" in index
        assert index["CAS"] == {"CASS", "CAST"}

    def test_drops_nan_and_empty(self):
        series = pd.Series(["CASS", None, "", "CAST"])
        index = _build_cdr3_neighbor_index(series)
        # Only non-null, non-empty entries contribute.
        all_indexed = set().union(*index.values())
        assert all_indexed == {"CASS", "CAST"}


class TestFindLevNeighbors:
    """Integration of the index + distance check used by
    ``match_clonotypes``."""

    def test_finds_substitution_neighbor(self):
        index = _build_cdr3_neighbor_index(pd.Series(["CASSXEQYF", "CASSAEQYF"]))
        # Query has Y at idx 4 — neither DB entry matches exactly but
        # both are within Lev distance 1.
        results = _find_lev_neighbors("CASSYEQYF", index, max_distance=1)
        result_set = {cdr3 for cdr3, _ in results}
        assert result_set == {"CASSXEQYF", "CASSAEQYF"}
        # Each has distance 1.
        assert all(d == 1 for _, d in results)

    def test_finds_deletion_neighbor(self):
        # DB has CAASS; query CASS is one deletion away.
        index = _build_cdr3_neighbor_index(pd.Series(["CAASS"]))
        results = _find_lev_neighbors("CASS", index, max_distance=1)
        assert ("CAASS", 1) in results

    def test_excludes_exact_match(self):
        """``_find_lev_neighbors`` returns only *near* neighbors;
        the exact match is handled by the caller separately."""
        index = _build_cdr3_neighbor_index(pd.Series(["CASS"]))
        results = _find_lev_neighbors("CASS", index, max_distance=1)
        assert results == []

    def test_empty_query_returns_empty(self):
        index = _build_cdr3_neighbor_index(pd.Series(["CASS"]))
        assert _find_lev_neighbors("", index, max_distance=1) == []


class TestMatchClonotypesExactModeUnchanged:
    """``match_mode="exact"`` (default) preserves pre-2.1 behaviour."""

    def test_exact_mode_is_default(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAS"], "CDR3_beta": ["CASSYEQYF"],
        })
        # DB entry differs by one substitution from query.
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSAEQYF",
            "epitope": "EPI", "antigen_gene": "ANT",
        }])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        # Without fuzzy mode, no match.
        assert out["db_match"].iloc[0] is False or out["db_match"].iloc[0] == 0

    def test_exact_mode_explicit_no_fuzzy(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAS"], "CDR3_beta": ["CASSYEQYF"],
        })
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSAEQYF",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="exact",
            verbose=False, show_progress=False,
        )
        assert not out["db_match"].iloc[0]


class TestMatchClonotypesLevenshteinMode:
    """The new fuzzy fallback fires only when an exact match isn't found."""

    def test_exact_ab_still_wins_over_fuzzy(self):
        """When αβ exact matches, distance must be 0 and the strength
        is plain ``ab`` (no ``_near`` suffix)."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "CAS", "cdr3_beta": "CASS", "epitope": "E1"}])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein", max_distance=1,
            verbose=False, show_progress=False,
        )
        assert out["db_match"].iloc[0]
        assert out["db_match_strength"].iloc[0] == "ab"
        assert out["db_match_distance"].iloc[0] == 0

    def test_exact_beta_only_wins_over_fuzzy(self):
        """A β-only exact match still takes precedence over a fuzzy
        β-only neighbor at distance 1."""
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS", "epitope": "E1"},
            # Fuzzy neighbor — would be a Lev-1 hit but exact wins.
            {"cdr3_alpha": "OTHER", "cdr3_beta": "CAST", "epitope": "E2"},
        ])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        # Exact β match — distance 0, strength ``b_only``.
        assert out["db_match_strength"].iloc[0] == "b_only"
        assert out["db_match_distance"].iloc[0] == 0

    def test_fuzzy_match_fires_when_no_exact_hit(self):
        """The actual #57 scenario: no exact αβ or β-only, but a
        single-substitution β neighbor in the DB."""
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASSPTGLGQPQHF"],
        })
        # DB has a MART-1-like canonical β differing by 1 substitution
        # at position 8 (Q→N).
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGLGNPQHF",
            "epitope": "ELAGIGILTV", "antigen_gene": "MLANA",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein", max_distance=1,
            verbose=False, show_progress=False,
        )
        assert out["db_match"].iloc[0]
        assert out["db_match_strength"].iloc[0] == "b_only_near"
        assert out["db_match_distance"].iloc[0] == 1
        # And the antigen/epitope info from the near match is recorded.
        assert out["db_epitope"].iloc[0] == "ELAGIGILTV"
        assert out["db_protein"].iloc[0] == "MLANA"

    def test_indel_neighbor_also_matches(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASSPTGLGNPQHF"],
        })
        # DB has the same β with one extra residue (insertion) → Lev-1.
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGLGNXPQHF",
            "epitope": "EPI",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        assert out["db_match_strength"].iloc[0] == "b_only_near"
        assert out["db_match_distance"].iloc[0] == 1

    def test_no_fuzzy_match_when_distance_too_far(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASSPTGLGQPQHF"],
        })
        # DB β differs by 3 substitutions.
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGYGRPRHF",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        assert not out["db_match"].iloc[0]
        assert pd.isna(out["db_match_distance"].iloc[0]) or out["db_match_distance"].iloc[0] is None

    def test_fuzzy_records_cross_species_suffix(self):
        """When a fuzzy match's host is non-human, the strength gets
        the ``_cross`` suffix on top of ``_near``."""
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASSPTGLGQPQHF"],
        })
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGLGNPQHF",
            "host_species": "Mus musculus",
            "epitope": "EPI",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        assert out["db_match_strength"].iloc[0] == "b_only_near_cross"
        assert out["db_match_distance"].iloc[0] == 1


class TestMatchModeValidation:
    def test_invalid_mode_raises(self):
        from tcrsift.validation import TCRsiftValidationError

        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["Y"]})
        db = _db([{"cdr3_alpha": "X", "cdr3_beta": "Y"}])
        with pytest.raises(TCRsiftValidationError, match="match_mode"):
            match_clonotypes(
                clones, db,
                match_mode="tcrdist",  # not in MATCH_MODES
                verbose=False, show_progress=False,
            )

    def test_max_distance_above_1_warns_and_clamps(self, caplog):
        import logging

        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS"}])
        with caplog.at_level(logging.WARNING, logger="tcrsift.annotate"):
            match_clonotypes(
                clones, db,
                match_mode="levenshtein", max_distance=3,
                verbose=False, show_progress=False,
            )
        assert any(
            "clamping" in r.message and "1" in r.message
            for r in caplog.records
        )


class TestDbMatchDistanceColumn:
    def test_column_present_on_no_db_short_circuit(self):
        from tcrsift.annotate import annotate_clonotypes

        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["Y"]})
        out = annotate_clonotypes(clones)
        assert "db_match_distance" in out.columns
        assert out["db_match_distance"].iloc[0] is None

    def test_column_present_after_match(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS"}])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert "db_match_distance" in out.columns
        assert out["db_match_distance"].iloc[0] == 0
