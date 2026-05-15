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

"""Tests for the #83 extensions to ``match_clonotypes``:

- ``db_host_species`` column populated from the DB's TCR donor species
- ``db_match_strength`` gains ``ab_cross`` / ``b_only_cross`` suffix
  when matched rows are non-human host
- ``db_category`` becomes ``"contradictory"`` when matches disagree on
  an informative category (rather than silently picking a mode)
"""

from __future__ import annotations

import pandas as pd

from tcrsift.annotate import (
    CATEGORY_CONTRADICTORY,
    CATEGORY_TUMOR_SELF,
    CATEGORY_VIRAL,
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
        df["species"] = "unknown"
    if "antigen_gene" not in df.columns:
        df["antigen_gene"] = ""
    if "epitope" not in df.columns:
        df["epitope"] = ""
    return df


class TestHostSpeciesColumn:
    def test_db_host_species_populated_when_host_present_in_db(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "host_species": "Homo sapiens", "species": "CMV",
        }])
        out = match_clonotypes(clones, db)
        assert out["db_host_species"].iloc[0] == "Homo sapiens"

    def test_db_host_species_is_none_when_absent_from_db(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS", "species": "CMV",
        }])
        out = match_clonotypes(clones, db)
        # Match succeeds, but host_species wasn't in the DB.
        assert out["db_match"].iloc[0]
        assert out["db_host_species"].iloc[0] is None

    def test_no_db_path_short_circuit_includes_host_column(self):
        """The no-DB short-circuit must still emit the new column so
        downstream code can rely on schema regardless of whether
        annotation actually ran."""
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["Y"]})
        out = annotate_clonotypes(clones)
        assert "db_host_species" in out.columns


class TestCrossSpeciesStrength:
    def test_mouse_host_yields_ab_cross_strength(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "host_species": "Mus musculus", "species": "OVA",
        }])
        out = match_clonotypes(clones, db)
        assert out["db_match_strength"].iloc[0] == "ab_cross"

    def test_human_host_yields_plain_ab_strength(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "host_species": "Homo sapiens", "species": "CMV",
        }])
        out = match_clonotypes(clones, db)
        assert out["db_match_strength"].iloc[0] == "ab"

    def test_mouse_host_with_b_only_fallback(self):
        clones = pd.DataFrame({"CDR3_alpha": ["X"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "OTHER", "cdr3_beta": "CASS",
            "host_species": "Mus musculus",
        }])
        out = match_clonotypes(clones, db, match_strictness="ab_with_partial")
        # β-only fallback + non-human host → "b_only_cross".
        assert out["db_match_strength"].iloc[0] == "b_only_cross"


class TestContradictoryCategory:
    def test_matches_with_disagreeing_categories_flagged(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        # Same clone matches both a viral and a tumor-self row.
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MART1"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_category"].iloc[0] == CATEGORY_CONTRADICTORY

    def test_agreeing_matches_keep_their_category(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        # Three matches, all viral with the same category.
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_category"].iloc[0] == CATEGORY_VIRAL

    def test_unknown_does_not_trigger_contradiction(self):
        """Mixing a real label with ``unknown`` is not a contradiction
        — pick the informative label."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MART1"},
            # A row with no species/antigen → ``unknown`` category.
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "", "antigen_gene": ""},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_category"].iloc[0] == CATEGORY_TUMOR_SELF
