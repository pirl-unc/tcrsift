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

"""Tests for #56: composite ``db_match_score`` with structured audit
columns.

New per-clone columns:
  - ``db_match_score`` ∈ [0, 1]
  - ``db_v_gene_agreement`` ∈ [0, 1] or None
  - ``db_n_studies`` int or None
  - ``db_category_votes`` "cat:frac;cat:frac;..." string or None
"""

from __future__ import annotations

import pandas as pd

from tcrsift.annotate import (
    _canonical_v_gene,
    _compute_match_score,
    annotate_clonotypes,
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


# ---------------------------------------------------------------------------
# Unit tests on the score helper


class TestCanonicalVGene:
    def test_strips_allele_suffix(self):
        assert _canonical_v_gene("TRBV6-1*01") == "TRBV6-1"

    def test_uppercases(self):
        assert _canonical_v_gene("trbv6-1") == "TRBV6-1"

    def test_strips_whitespace(self):
        assert _canonical_v_gene("  TRBV6-1  ") == "TRBV6-1"

    def test_non_string_returns_empty(self):
        assert _canonical_v_gene(None) == ""
        assert _canonical_v_gene(3) == ""
        assert _canonical_v_gene(float("nan")) == ""


class TestComputeMatchScore:
    """Unit tests for the formula. Documented in the helper's docstring;
    these tests pin the actual numbers across refactors."""

    def test_ab_exact_3_studies_v_agree(self):
        s = _compute_match_score("ab", distance=0, n_studies=3, v_gene_agreement=1.0)
        # base=1.0, v_factor=1.0, studies_factor=1-0.3/3=0.9 → 0.9
        assert s == 0.9

    def test_b_only_exact_1_study_v_agree(self):
        s = _compute_match_score(
            "b_only", distance=0, n_studies=1, v_gene_agreement=1.0,
        )
        # base=0.7, v_factor=1.0, studies_factor=0.7 → 0.49
        assert s == round(0.7 * 1.0 * 0.7, 3)

    def test_b_only_near_2_studies_v_disagree(self):
        s = _compute_match_score(
            "b_only_near", distance=1, n_studies=2, v_gene_agreement=0.0,
        )
        # base=0.7 * 0.7=0.49, v_factor=1-0.3=0.7, studies=1-0.15=0.85
        expected = round(0.49 * 0.7 * 0.85, 3)
        assert s == expected

    def test_ab_cross_penalty(self):
        s = _compute_match_score(
            "ab_cross", distance=0, n_studies=3, v_gene_agreement=1.0,
        )
        # base=1.0 * 0.5=0.5, v=1.0, studies=0.9 → 0.45
        assert s == round(0.5 * 1.0 * 0.9, 3)

    def test_b_only_near_cross_stacks_both_penalties(self):
        s = _compute_match_score(
            "b_only_near_cross", distance=1, n_studies=1, v_gene_agreement=1.0,
        )
        # base=0.7 * 0.7 * 0.5=0.245, v=1.0, studies=0.7 → 0.1715
        expected = round(0.7 * 0.7 * 0.5 * 1.0 * 0.7, 3)
        assert s == expected

    def test_unknown_v_gene_agreement_is_neutral(self):
        s_unknown = _compute_match_score(
            "ab", distance=0, n_studies=3, v_gene_agreement=None,
        )
        s_perfect = _compute_match_score(
            "ab", distance=0, n_studies=3, v_gene_agreement=1.0,
        )
        # Unknown v-gene is treated as 1.0 (neutral, not penalised).
        assert s_unknown == s_perfect

    def test_zero_or_missing_n_studies_falls_back_to_1(self):
        s_zero = _compute_match_score(
            "ab", distance=0, n_studies=0, v_gene_agreement=1.0,
        )
        s_one = _compute_match_score(
            "ab", distance=0, n_studies=1, v_gene_agreement=1.0,
        )
        s_none = _compute_match_score(
            "ab", distance=0, n_studies=None, v_gene_agreement=1.0,
        )
        assert s_zero == s_one == s_none

    def test_score_clamps_to_unit_interval(self):
        for strength in ("ab", "b_only", "ab_cross", "b_only_near_cross"):
            s = _compute_match_score(
                strength, distance=0, n_studies=10, v_gene_agreement=1.0,
            )
            assert 0.0 <= s <= 1.0


# ---------------------------------------------------------------------------
# Integration tests through match_clonotypes


class TestVGeneAgreementColumn:
    def test_agreement_when_clone_and_db_v_genes_match(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASS"],
            "beta_v_gene": ["TRBV6-1"],
        })
        db = _db([
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "v_beta": "TRBV6-1", "reference": "ref1"},
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "v_beta": "TRBV6-1*01", "reference": "ref2"},
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        # Both DB rows have V_beta that canonicalises to TRBV6-1 == clone's.
        assert out["db_v_gene_agreement"].iloc[0] == 1.0

    def test_partial_agreement(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASS"],
            "beta_v_gene": ["TRBV6-1"],
        })
        db = _db([
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "v_beta": "TRBV6-1", "reference": "r1"},
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "v_beta": "TRBV28", "reference": "r2"},
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert out["db_v_gene_agreement"].iloc[0] == 0.5

    def test_no_agreement(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASS"],
            "beta_v_gene": ["TRBV6-1"],
        })
        db = _db([
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "v_beta": "TRBV28", "reference": "r1"},
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert out["db_v_gene_agreement"].iloc[0] == 0.0

    def test_none_when_clone_has_no_v_gene(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASS"],
        })
        db = _db([{
            "cdr3_alpha": "Y", "cdr3_beta": "CASS",
            "v_beta": "TRBV6-1", "reference": "r1",
        }])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        # Column exists but value is None (couldn't compute).
        assert pd.isna(out["db_v_gene_agreement"].iloc[0])

    def test_none_when_db_has_no_v_gene_column(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASS"],
            "beta_v_gene": ["TRBV6-1"],
        })
        # DB has no v_beta column at all (e.g. IEDB-only annotation).
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": "r1"}])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert pd.isna(out["db_v_gene_agreement"].iloc[0])


class TestNStudiesColumn:
    def test_counts_unique_references(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": "ref1"},
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": "ref2"},
            # Duplicate of ref1 — same study, doesn't count twice.
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": "ref1"},
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert int(out["db_n_studies"].iloc[0]) == 2

    def test_zero_when_all_references_missing(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": None}])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert int(out["db_n_studies"].iloc[0]) == 0

    def test_none_when_db_has_no_reference_column(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS"}])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert pd.isna(out["db_n_studies"].iloc[0])


class TestCategoryVotesColumn:
    def test_per_category_fractions(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([
            # 2 viral + 1 self → 0.67 / 0.33 fraction.
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "species": "Epstein-Barr virus", "antigen_gene": "EBNA-3"},
            {"cdr3_alpha": "Y", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MART1"},
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        votes = out["db_category_votes"].iloc[0]
        # Expected ordering: viral first (most common).
        assert votes.startswith("viral:0.67")
        assert "tumor_self:0.33" in votes


class TestMatchScoreIntegration:
    def test_score_present_on_match(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{"cdr3_alpha": "Y", "cdr3_beta": "CASS", "reference": "ref1"}])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        score = out["db_match_score"].iloc[0]
        assert isinstance(score, float)
        assert 0.0 <= score <= 1.0

    def test_score_present_on_no_db_short_circuit(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["Y"]})
        out = annotate_clonotypes(clones)
        assert "db_match_score" in out.columns
        assert "db_v_gene_agreement" in out.columns
        assert "db_n_studies" in out.columns
        assert "db_category_votes" in out.columns
        assert out["db_match_score"].iloc[0] is None

    def test_ab_exact_high_quality_scores_high(self):
        """αβ exact + 3 references + V-gene agreement → highest score."""
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"],
            "beta_v_gene": ["TRBV6-1"],
        })
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "v_beta": "TRBV6-1", "reference": f"ref{i}"}
            for i in range(3)
        ])
        out = match_clonotypes(clones, db, verbose=False, show_progress=False)
        assert out["db_match_score"].iloc[0] >= 0.85

    def test_b_only_near_cross_low_quality_scores_low(self):
        """Fuzzy β-only + non-human host + single study + V disagree
        → score should be deep in 'low-confidence' territory."""
        clones = pd.DataFrame({
            "CDR3_alpha": ["X"], "CDR3_beta": ["CASSPTGLGQPQHF"],
            "beta_v_gene": ["TRBV6-1"],
        })
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGLGNPQHF",
            "v_beta": "TRBV28", "reference": "ref1",
            "host_species": "Mus musculus",
        }])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        assert out["db_match_strength"].iloc[0] == "b_only_near_cross"
        assert out["db_match_score"].iloc[0] < 0.3

    def test_score_high_passes_low_passes_filter_threshold(self):
        """The motivating use case from #56: ``score >= 0.6`` separates
        high-confidence αβ from low-confidence β-only-near."""
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAS", "X"],
            "CDR3_beta":  ["CASS", "CASSPTGLGQPQHF"],
            "beta_v_gene": ["TRBV6-1", "TRBV6-1"],
        })
        db = _db([
            # Clone 0: αβ exact, 3 studies, V agrees → high.
            *[
                {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
                 "v_beta": "TRBV6-1", "reference": f"r{i}"}
                for i in range(3)
            ],
            # Clone 1: β-only fuzzy, 1 study, V disagrees → low.
            {"cdr3_alpha": "OTHER", "cdr3_beta": "CASSPTGLGNPQHF",
             "v_beta": "TRBV28", "reference": "fuzzy_ref"},
        ])
        out = match_clonotypes(
            clones, db,
            match_mode="levenshtein",
            verbose=False, show_progress=False,
        )
        # High-quality αβ exact passes the score-≥-0.6 gate.
        assert out["db_match_score"].iloc[0] >= 0.6
        # Low-quality fuzzy + V-disagree does not.
        assert out["db_match_score"].iloc[1] < 0.6
