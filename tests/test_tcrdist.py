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

"""Unit tests for native zero-dependency TCRdist similarity and database matching."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift.annotate import match_clonotypes
from tcrsift.config import AnnotateConfig
from tcrsift.tcrdist import (
    DEFAULT_BETA_RADIUS,
    DEFAULT_PAIRED_RADIUS,
    TCRDistIndex,
    cdr3_distance,
    tcr_distance,
    trim_cdr3,
    v_gene_distance,
)
from tcrsift.validation import TCRsiftValidationError


class TestTCRDistFormula:
    """Exact Dash (2017) / tcrdist3 numerical parity tests."""

    def test_trim_cdr3_removes_anchors(self):
        assert trim_cdr3("CASSIRSSYEQYF") == "SIRSSYEQ"
        assert trim_cdr3("  cassirssyeqyf ") == "SIRSSYEQ"
        # Short strings (<= 5 aa) are preserved untouched
        assert trim_cdr3("CASSF") == "CASSF"
        assert trim_cdr3(None) == ""  # type: ignore[arg-type]

    def test_cdr3_distance_matches_tcrdist3(self):
        # Identical
        assert cdr3_distance("CASSIRSSYEQYF", "CASSIRSSYEQYF") == 0
        # Conservative I -> V (BLOSUM62=3 -> cost 1 -> *3 = 3)
        assert cdr3_distance("CASSIRSSYEQYF", "CASSVRSSYEQYF") == 3
        # Disruptive I -> D (BLOSUM62=-3 -> capped cost 4 -> *3 = 12)
        assert cdr3_distance("CASSIRSSYEQYF", "CASSDRSSYEQYF") == 12
        # Single-residue indel (gap_penalty=4 -> *3 = 12)
        assert cdr3_distance("CASSIRSSYEQYF", "CASSIRSYEQYF") == 12
        # Conserved 3 N-term / 2 C-term anchor differences are trimmed away
        assert cdr3_distance("CASSIRSSYEQYF", "CSASIRSSYEQFF") == 0

    def test_v_gene_distance_human_and_mouse(self):
        assert v_gene_distance("TRBV19*01", "TRBV19*02", "beta", species="human") == 0
        assert v_gene_distance("TRBV19*01", "TRBV28*01", "beta", species="human") == 31
        # Format / DV-slash canonicalization
        assert v_gene_distance("TRAV14/DV4*01", "TRAV14DV4", "alpha", species="human") == 0
        # Mouse V genes
        assert v_gene_distance("Trbv13-2*01", "TRBV13-2", "beta", species="mouse") == 0
        d_mouse = v_gene_distance("Trbv13-2", "Trbv19", "beta", species="mouse")
        assert isinstance(d_mouse, int) and d_mouse > 0
        # Unknown V gene returns None
        assert v_gene_distance("TRBV999", "TRBV19", "beta", species="human") is None

    def test_tcr_distance_joint_and_fallback(self):
        # Joint paired = d_alpha + d_beta
        d = tcr_distance(
            cdr3_a1="CAGPGGSQGNLIF",
            v_a1="TRAV27*01",
            cdr3_b1="CASSIRSSYEQYF",
            v_b1="TRBV19*01",
            cdr3_a2="CAGPGGSQGNLIF",
            v_a2="TRAV27*01",
            cdr3_b2="CASSVRSSYEQYF",
            v_b2="TRBV28*01",
        )
        assert d == 3 + 31  # conservative beta CDR3 swap (3) + TRBV19 vs TRBV28 (31)

        # Graceful CDR3-only when V genes omitted
        d_cdr3_only = tcr_distance(
            cdr3_b1="CASSIRSSYEQYF",
            cdr3_b2="CASSVRSSYEQYF",
        )
        assert d_cdr3_only == 3


class TestTCRDistIndexExactness:
    """Verify branch-and-bound TCRDistIndex matches brute-force tcr_distance exactly."""

    def test_paired_index_exact_vs_brute_force(self):
        ref = [
            ("CAGPGGSQGNLIF", "TRAV27", "CASSIRSSYEQYF", "TRBV19"),
            ("CAGPGGSQGNLIF", "TRAV27", "CASSVRSSYEQYF", "TRBV19"),
            ("CAGPGGSQGNLIF", "TRAV27", "CASSIRSYEQYF", "TRBV19"),
            ("CAVSDGGSQGNLIF", "TRAV12-2", "CASSIRSSYEQYF", "TRBV28"),
            ("CAASDTNTGKLTF", "TRAV14/DV4", "CASSLGQAYEQYF", "TRBV7-9"),
            # Duplicate entry to confirm all tied DB positions are returned
            ("CAGPGGSQGNLIF", "TRAV27*01", "CASSVRSSYEQYF", "TRBV19*01"),
        ]
        idx = TCRDistIndex(
            np.arange(len(ref), dtype=np.int64),
            [r[2] for r in ref],
            v_beta=[r[3] for r in ref],
            cdr3_alpha=[r[0] for r in ref],
            v_alpha=[r[1] for r in ref],
            species="human",
            paired=True,
        )

        queries = [
            ("CAGPGGSQGNLIF", "TRAV27", "CASSIRSSYEQYF", "TRBV19"),
            ("CAGPGGSQGNLIF", "TRAV27", "CASSVRSSYEQYF", "TRBV19"),
            ("CAGPGGSQGNLIF", "TRAV27", "CASSLRSSYEQYF", "TRBV19"),
            ("CAGPGGSQGNLIF", "TRAV27", "CASSIRSSYEQYF", "TRBV28"),
            ("CAASDTNTGKLTF", "TRAV14DV4", "CASSLGQVYEQYF", "TRBV7-9"),
        ]
        for radius in (0, 12, 24, 48):
            for qa, qva, qb, qvb in queries:
                dists = [
                    tcr_distance(
                        cdr3_a1=qa,
                        v_a1=qva,
                        cdr3_b1=qb,
                        v_b1=qvb,
                        cdr3_a2=ra,
                        v_a2=rva,
                        cdr3_b2=rb,
                        v_b2=rvb,
                    )
                    for ra, rva, rb, rvb in ref
                ]
                min_d = min(dists)
                expected = (
                    (min_d, np.array([i for i, d in enumerate(dists) if d == min_d], dtype=np.int64))
                    if min_d <= radius
                    else None
                )
                got = idx.query_nearest(
                    qb, v_beta=qvb, cdr3_alpha=qa, v_alpha=qva, max_distance=radius
                )
                if expected is None:
                    assert got is None
                else:
                    assert got is not None
                    assert got[0] == expected[0]
                    np.testing.assert_array_equal(got[1], expected[1])

    def test_beta_only_index_mouse(self):
        ref = [
            ("CASSDDRVNERLFF", "Trbv13-2"),
            ("CASSDERVNERLFF", "Trbv13-2"),
            ("CASSQDRGTEVFF", "Trbv19"),
        ]
        idx = TCRDistIndex(
            np.arange(len(ref), dtype=np.int64),
            [r[0] for r in ref],
            v_beta=[r[1] for r in ref],
            species="mouse",
            paired=False,
        )
        hit0 = idx.query_nearest("CASSDDRVNERLFF", v_beta="TRBV13-2*01", max_distance=0)
        assert hit0 is not None and hit0[0] == 0 and hit0[1].tolist() == [0]

        hit_near = idx.query_nearest("CASSDNRVNERLFF", v_beta="Trbv13-2", max_distance=12)
        assert hit_near is not None and hit_near[0] == 9 and hit_near[1].tolist() == [0]


def _make_db(rows: list[dict]) -> pd.DataFrame:
    defaults = {
        "cdr3_alpha": "",
        "cdr3_beta": "",
        "v_alpha": "",
        "v_beta": "",
        "epitope": "GILGFVFTL",
        "antigen_gene": "M1",
        "species": "InfluenzaA",
        "host_species": "Homo sapiens",
        "mhc_allele": "HLA-A*02:01",
        "database": "VDJdb",
        "reference": "PMID:12345",
        "is_viral": True,
    }
    return pd.DataFrame([{**defaults, **r} for r in rows])


class TestMatchClonotypesTCRDist:
    """End-to-end integration tests for match_clonotypes(..., match_mode='tcrdist')."""

    def test_paired_strict_ab_with_v_genes(self):
        db = _make_db(
            [
                {
                    "cdr3_alpha": "CAGPGGSQGNLIF",
                    "v_alpha": "TRAV27*01",
                    "cdr3_beta": "CASSIRSSYEQYF",
                    "v_beta": "TRBV19*01",
                    "epitope": "GILGFVFTL",
                },
                {
                    "cdr3_alpha": "CAASDTNTGKLTF",
                    "v_alpha": "TRAV14/DV4",
                    "cdr3_beta": "CASSLGQAYEQYF",
                    "v_beta": "TRBV7-9",
                    "epitope": "NLVPMVATV",
                    "antigen_gene": "pp65",
                    "species": "CMV",
                },
            ]
        )
        clones = pd.DataFrame(
            {
                "CDR3_alpha": ["CAGPGGSQGNLIF", "CAGPGGSQGNLIF", "CAGPGGSQGNLIF", "CAAXXXXXXKLTF"],
                "CDR3_beta": ["CASSIRSSYEQYF", "CASSVRSSYEQYF", "CASSIRSSYEQYF", "CASSYYYYYYEQYF"],
                "alpha_v_gene": ["TRAV27", "TRAV27", "TRAV27", "TRAV1-2"],
                "beta_v_gene": ["TRBV19", "TRBV19", "TRBV28", "TRBV20-1"],
            }
        )

        # Default paired radius = 48: row 0 exact (d=0), row 1 conservative swap (d=3),
        # row 2 V-gene shift TRBV19->TRBV28 (d=31 <= 48), row 3 unmatched
        out48 = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            match_mode="tcrdist",
            verbose=False,
            show_progress=False,
        )
        assert out48["db_match"].tolist() == [True, True, True, False]
        assert out48["db_match_strength"].tolist() == ["ab", "ab_near", "ab_near", None]
        assert out48["db_match_distance"].tolist() == [0, 3, 31, None]
        assert out48["db_epitope"].iloc[:3].tolist() == ["GILGFVFTL"] * 3

        # Tight radius = 24 excludes row 2 (d=31) while keeping rows 0 and 1
        out24 = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            match_mode="tcrdist",
            max_distance=24,
            verbose=False,
            show_progress=False,
        )
        assert out24["db_match"].tolist() == [True, True, False, False]

        # Zero radius keeps only exact V-gene + trimmed CDR3 match
        out0 = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            match_mode="tcrdist",
            max_distance=0,
            verbose=False,
            show_progress=False,
        )
        assert out0["db_match"].tolist() == [True, False, False, False]
        assert out0["db_match_distance"].iloc[0] == 0

    def test_ab_with_partial_fallback_and_b_only(self):
        db = _make_db(
            [
                {
                    "cdr3_alpha": "CAGPGGSQGNLIF",
                    "v_alpha": "TRAV27",
                    "cdr3_beta": "CASSIRSSYEQYF",
                    "v_beta": "TRBV19",
                    "epitope": "GILGFVFTL",
                }
            ]
        )
        # Clone has completely different alpha, but identical beta + V_beta -> falls back to b_only
        clones = pd.DataFrame(
            {
                "CDR3_alpha": ["CAVRRRRRRNKLIF"],
                "CDR3_beta": ["CASSVRSSYEQYF"],
                "alpha_v_gene": ["TRAV1-2"],
                "beta_v_gene": ["TRBV19"],
            }
        )
        out_partial = match_clonotypes(
            clones,
            db,
            match_strictness="ab_with_partial",
            match_mode="tcrdist",
            verbose=False,
            show_progress=False,
        )
        assert bool(out_partial["db_match"].iloc[0]) is True
        assert bool(out_partial["db_match_partial"].iloc[0]) is True
        assert out_partial["db_match_strength"].iloc[0] == "b_only_near"
        assert out_partial["db_match_distance"].iloc[0] == 3

        out_b = match_clonotypes(
            clones,
            db,
            match_strictness="b_only",
            match_mode="tcrdist",
            max_distance=12,
            verbose=False,
            show_progress=False,
        )
        assert bool(out_b["db_match"].iloc[0]) is True
        assert bool(out_b["db_match_partial"].iloc[0]) is False
        assert out_b["db_match_strength"].iloc[0] == "b_only_near"
        assert out_b["db_match_distance"].iloc[0] == 3

    def test_graceful_cdr3_only_fallback_when_v_genes_missing(self):
        db = _make_db(
            [
                {
                    "cdr3_alpha": "CAGPGGSQGNLIF",
                    "v_alpha": "TRAV27",
                    "cdr3_beta": "CASSIRSSYEQYF",
                    "v_beta": "TRBV19",
                }
            ]
        )
        # Query table lacks V-gene columns completely -> gracefully scores BLOSUM62 CDR3 only
        clones = pd.DataFrame(
            {
                "CDR3_alpha": ["CAGPGGSQGNLIF"],
                "CDR3_beta": ["CASSVRSSYEQYF"],
            }
        )
        out = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            match_mode="tcrdist",
            max_distance=12,
            verbose=False,
            show_progress=False,
        )
        assert bool(out["db_match"].iloc[0]) is True
        assert out["db_match_distance"].iloc[0] == 3

    def test_config_and_validation(self):
        cfg = AnnotateConfig(match_mode="tcrdist", max_distance=0)
        assert cfg.match_mode == "tcrdist"
        assert cfg.max_distance == 0
        cfg_default = AnnotateConfig(match_mode="tcrdist")
        assert cfg_default.max_distance is None
        with pytest.raises(ValueError, match="max_distance"):
            AnnotateConfig(match_mode="tcrdist", max_distance=-1)
        with pytest.raises(ValueError, match="max_distance"):
            AnnotateConfig(match_mode="levenshtein", max_distance=0)

        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["Y"]})
        db = _make_db([{"cdr3_alpha": "X", "cdr3_beta": "Y"}])
        with pytest.raises(TCRsiftValidationError, match="max_distance"):
            match_clonotypes(
                clones,
                db,
                match_mode="tcrdist",
                max_distance=-5,
                verbose=False,
                show_progress=False,
            )
        assert DEFAULT_PAIRED_RADIUS == 48
        assert DEFAULT_BETA_RADIUS == 12
