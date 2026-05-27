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
  an informative category, with dominant/detail/count audit columns
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
        assert "db_category_dominant" in out.columns
        assert "db_category_detail" in out.columns
        assert "db_category_counts" in out.columns


class TestCrossSpeciesStrength:
    def test_vdjdb_homosapiens_host_yields_plain_ab_strength(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "host_species": "HomoSapiens", "species": "Human cytomegalovirus",
        }])
        out = match_clonotypes(clones, db, show_progress=False)
        assert out["db_host_species"].iloc[0] == "Homo sapiens"
        assert out["db_species"].iloc[0] == "CMV"
        assert out["db_match_strength"].iloc[0] == "ab"

    def test_mouse_host_yields_ab_cross_strength(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "host_species": "MusMusculus", "species": "OVA",
        }])
        out = match_clonotypes(clones, db, show_progress=False)
        assert out["db_host_species"].iloc[0] == "Mus musculus"
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
        assert out["db_category_dominant"].iloc[0] is None
        assert out["db_category_detail"].iloc[0] == CATEGORY_CONTRADICTORY
        assert out["db_category_counts"].iloc[0] == "tumor_self:1;viral:1"

    def test_disagreeing_categories_with_dominant_label_audited(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        # Viral is the dominant informative category, but there is also
        # a tumor-self match and an unknown row. Unknown rows are ignored
        # whenever informative category rows exist.
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Cytomegalovirus", "antigen_gene": "pp65"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Epstein-Barr virus", "antigen_gene": "EBNA-3"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MART1"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "", "antigen_gene": ""},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_category"].iloc[0] == CATEGORY_CONTRADICTORY
        assert out["db_category_dominant"].iloc[0] == CATEGORY_VIRAL
        assert out["db_category_detail"].iloc[0] == "contradictory_dominant_viral"
        assert out["db_category_counts"].iloc[0] == "viral:2;tumor_self:1"

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
        assert out["db_category_dominant"].iloc[0] == CATEGORY_VIRAL
        assert out["db_category_detail"].iloc[0] == CATEGORY_VIRAL
        assert out["db_category_counts"].iloc[0] == "viral:3"

    def test_unknown_does_not_trigger_contradiction(self):
        """Mixing a real label with ``unknown`` is not a contradiction
        — drop the unknown row and pick the informative label."""
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
        assert out["db_category_dominant"].iloc[0] == CATEGORY_TUMOR_SELF
        assert out["db_category_detail"].iloc[0] == CATEGORY_TUMOR_SELF
        assert out["db_category_counts"].iloc[0] == "tumor_self:1"

    def test_all_unknown_matches_stay_unknown(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "", "antigen_gene": ""},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": None, "antigen_gene": ""},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_category"].iloc[0] == "unknown"
        assert out["db_category_dominant"].iloc[0] is None
        assert out["db_category_detail"].iloc[0] == "unknown"
        assert out["db_category_counts"].iloc[0] == "unknown:2"


class TestAntigenIdentityCoherence:
    """#88: db_protein / db_protein_canonical / db_species / db_epitope
    must all describe the SAME antigen when matches span multiple DB
    rows. Previously each column took an independent mode, so the
    canonical label could be from one row, protein/species/epitope
    from another."""

    def test_mart1_majority_picks_consistent_mart1_metadata(self):
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        # Three MART1 matches, one SARS-CoV-2 ORF1ab match. Previously
        # mode-per-column could pick MART1 as the raw protein but
        # SARS-CoV-2 ORF1ab as the canonical if the latter had a
        # higher count there (or even on tiebreaker quirks). With the
        # group-by-canonical fix, all antigen-identity columns come
        # from the MART1 rows.
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV", "db_protein_canonical": "MART-1"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV", "db_protein_canonical": "MART-1"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV", "db_protein_canonical": "MART-1"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "ORF1ab",
             "epitope": "APKEIIFLEGETL",
             "db_protein_canonical": "SARS-CoV-2 ORF1ab"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_protein_canonical"].iloc[0] == "MART-1"
        assert out["db_protein"].iloc[0] == "MLANA"
        assert out["db_species"].iloc[0] == "Homo sapiens"
        assert out["db_epitope"].iloc[0] == "ALGIGILTV"

    def test_canonical_and_raw_protein_stay_coupled_on_tie(self):
        """Even when raw protein has more rows than the canonical, the
        canonical (user-facing label) is the grouping key — and once
        chosen, raw protein/species/epitope come from the rows tagged
        with that canonical."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            # 2 SARS-CoV-2 matches (canonical wins by majority)
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "ORF1ab",
             "epitope": "APKEIIFLEGETL",
             "db_protein_canonical": "SARS-CoV-2 ORF1ab"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "ORF1ab",
             "epitope": "APKEIIFLEGETL",
             "db_protein_canonical": "SARS-CoV-2 ORF1ab"},
            # 1 HIV match
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "HIV-1", "antigen_gene": "Gag",
             "epitope": "SLYNTVATL",
             "db_protein_canonical": "HIV-1 Gag"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_protein_canonical"].iloc[0] == "SARS-CoV-2 ORF1ab"
        assert out["db_protein"].iloc[0] == "ORF1ab"
        assert out["db_species"].iloc[0] == "SARS-CoV-2"
        assert out["db_epitope"].iloc[0] == "APKEIIFLEGETL"

    def test_single_row_match_passes_through_metadata_intact(self):
        """The smallest match-set: one DB row. The antigen-group
        machinery must reduce to a no-op restriction and pass all
        antigen-identity columns through verbatim."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([{
            "cdr3_alpha": "CAS", "cdr3_beta": "CASS",
            "species": "Homo sapiens", "antigen_gene": "MLANA",
            "epitope": "ALGIGILTV", "db_protein_canonical": "MART-1",
        }])
        out = match_clonotypes(clones, db)
        assert out["db_protein"].iloc[0] == "MLANA"
        assert out["db_protein_canonical"].iloc[0] == "MART-1"
        assert out["db_species"].iloc[0] == "Homo sapiens"
        assert out["db_epitope"].iloc[0] == "ALGIGILTV"

    def test_falls_back_to_antigen_gene_when_canonical_all_nan(self):
        """If `db_protein_canonical` column exists but is entirely
        empty/NaN for this match set, fall through to `antigen_gene`
        for grouping rather than silently un-grouping."""
        import numpy as np

        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV", "db_protein_canonical": np.nan},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV", "db_protein_canonical": np.nan},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "ORF1ab",
             "epitope": "APKEIIFLEGETL", "db_protein_canonical": np.nan},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_protein"].iloc[0] == "MLANA"
        assert out["db_species"].iloc[0] == "Homo sapiens"
        assert out["db_epitope"].iloc[0] == "ALGIGILTV"

    def test_host_species_and_cross_strength_come_from_antigen_matches(self):
        """Regression: previously host_species and the _cross strength
        suffix were derived from the full matches frame BEFORE the
        antigen group was picked. A mouse-model row in the match set
        could flip db_host_species and add _cross even when the
        winning canonical antigen came from a human-host row. Same
        archetype as #88 but on the cross-species columns."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            # Two human-host EBV rows (the winning antigen)
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Epstein-Barr virus", "antigen_gene": "EBNA-3",
             "host_species": "Homo sapiens",
             "db_protein_canonical": "EBV EBNA-3"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Epstein-Barr virus", "antigen_gene": "EBNA-3",
             "host_species": "Homo sapiens",
             "db_protein_canonical": "EBV EBNA-3"},
            # An unrelated mouse-host row tagged with a different antigen
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "OVA", "antigen_gene": "OVA",
             "host_species": "Mus musculus",
             "db_protein_canonical": "OVA"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_protein_canonical"].iloc[0] == "EBV EBNA-3"
        assert out["db_host_species"].iloc[0] == "Homo sapiens"
        assert out["db_match_strength"].iloc[0] == "ab"  # NOT ab_cross

    def test_canonical_drills_down_within_antigen_gene_group(self):
        """When grouped on antigen_gene (because db_protein_canonical
        was all-NaN or absent for the picked rows), the picked subset
        may still carry multiple canonical values. The fix drills down
        so the written canonical comes from the same rows as
        species/epitope."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            # Two rows share antigen_gene="Spike glycoprotein" but
            # differ on canonical (and on species/epitope).
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "Spike glycoprotein",
             "epitope": "VYAWNRKRI",
             "db_protein_canonical": "SARS-CoV-2 Spike"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "Spike glycoprotein",
             "epitope": "VYAWNRKRI",
             "db_protein_canonical": "SARS-CoV-2 Spike"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-1", "antigen_gene": "Spike glycoprotein",
             "epitope": "DIFFERENT_EPI",
             "db_protein_canonical": "SARS-CoV-1 Spike"},
        ])
        out = match_clonotypes(clones, db)
        # The dominant canonical (the SARS-CoV-2 pair) should drive
        # both protein_canonical AND species/epitope — they must
        # describe the same DB rows.
        assert out["db_protein_canonical"].iloc[0] == "SARS-CoV-2 Spike"
        assert out["db_species"].iloc[0] == "SARS-CoV-2"
        assert out["db_epitope"].iloc[0] == "VYAWNRKRI"

    def test_falls_back_to_antigen_gene_when_canonical_absent(self):
        """If no `db_protein_canonical` column exists, group on
        `antigen_gene` instead."""
        clones = pd.DataFrame({"CDR3_alpha": ["CAS"], "CDR3_beta": ["CASS"]})
        db = _db([
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "Homo sapiens", "antigen_gene": "MLANA",
             "epitope": "ALGIGILTV"},
            {"cdr3_alpha": "CAS", "cdr3_beta": "CASS",
             "species": "SARS-CoV-2", "antigen_gene": "ORF1ab",
             "epitope": "APKEIIFLEGETL"},
        ])
        out = match_clonotypes(clones, db)
        assert out["db_protein"].iloc[0] == "MLANA"
        # When grouped by antigen_gene, species/epitope stay coupled.
        assert out["db_species"].iloc[0] == "Homo sapiens"
        assert out["db_epitope"].iloc[0] == "ALGIGILTV"
