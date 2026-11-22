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

"""Tests for gene-vs-protein/proteoform separation, organism-aware antigen
canonicalization, host-species isolation, and IEDB v3 receptor parsing."""

from __future__ import annotations

import pandas as pd

from tcrsift.annotate import (
    CATEGORY_OTHER,
    CATEGORY_TUMOR_SELF,
    CATEGORY_VIRAL,
    canonicalize_antigen,
    classify_category,
    load_iedb,
    match_clonotypes,
    resolve_antigen_gene,
)


class TestAntigenCanonicalizationPrecision:
    """Word-boundary and organism-aware canonicalization without false collisions."""

    def test_host_snrpg_never_collides_with_flu_np(self):
        assert (
            canonicalize_antigen(
                "Small nuclear ribonucleoprotein G", species="Homo sapiens"
            )
            == "Small nuclear ribonucleoprotein G"
        )
        assert canonicalize_antigen("Nucleoprotein", species="Influenza A") == "Flu NP"
        assert canonicalize_antigen("Nucleoprotein", species="LCMV") == "LCMV NP"
        assert canonicalize_antigen("Nucleoprotein", species="Sendai virus") == "Sendai NP"
        assert canonicalize_antigen("Nucleoprotein", species="SARS-CoV-2") == "SARS-CoV-2 N"

    def test_shared_viral_protein_names_respect_organism(self):
        assert canonicalize_antigen("NS3", species="HCV") == "HCV NS3"
        assert canonicalize_antigen("NS3", species="Dengue virus") == "Dengue NS3"
        assert canonicalize_antigen("NS3", species="Zika virus") == "Zika NS3"
        assert canonicalize_antigen("Spike glycoprotein", species="SARS-CoV-1") == "SARS-CoV-1 Spike"
        assert canonicalize_antigen("Spike glycoprotein", species="SARS-CoV-2") == "SARS-CoV-2 Spike"
        assert canonicalize_antigen("IE1", species="Mouse CMV") == "MCMV IE1"
        assert canonicalize_antigen("M45", species="Mouse CMV") == "MCMV M45"
        assert canonicalize_antigen("GP33", species="LCMV") == "LCMV GP"
        assert canonicalize_antigen("p15E", species="MuLV") == "MuLV p15E"

    def test_paralog_and_substring_word_boundaries(self):
        # ERBB2 vs ERBB2IP
        assert canonicalize_antigen("ERBB2") == "HER2"
        assert canonicalize_antigen("ERBB2IP") == "ERBB2IP"

        # CEACAM5 (CEA) vs other CEACAM paralogs
        assert canonicalize_antigen("CEACAM5") == "CEA"
        assert canonicalize_antigen("Carcinoembryonic antigen-related cell adhesion molecule 5") == "CEA"
        assert canonicalize_antigen("CEACAM1") == "CEACAM1"
        assert canonicalize_antigen("CEACAM6") == "CEACAM6"
        assert (
            canonicalize_antigen("Carcinoembryonic antigen-related cell adhesion molecule 1")
            == "Carcinoembryonic antigen-related cell adhesion molecule 1"
        )

        # Numbered MAGE paralogs keep their exact numbers rather than collapsing into MAGE-A1
        assert canonicalize_antigen("MAGE-A1") == "MAGE-A1"
        assert canonicalize_antigen("MAGE-A3") == "MAGE-A3"
        assert canonicalize_antigen("MAGE-A6") == "MAGE-A6"
        assert canonicalize_antigen("MAGE-A10") == "MAGE-A10"
        assert canonicalize_antigen("MAGEA12") == "MAGE-A12"
        assert canonicalize_antigen("MAGEC2") == "MAGE-C2"

        # Tyrosinase-related protein 1 / 2 never collapse into Tyrosinase
        assert canonicalize_antigen("Tyrosinase-related protein 1") == "TRP-1"
        assert canonicalize_antigen("Tyrosinase-related protein 2") == "TRP-2"
        assert canonicalize_antigen("Tyrosinase") == "Tyrosinase"


class TestGeneVsProteinSeparation:
    """Explicit separation between db_protein_canonical and db_antigen_gene."""

    def test_ctag1a_ctag1b_ctag2_disambiguation(self):
        # CTAG1A and CTAG1B encode the identical NY-ESO-1 protein, while CTAG2 encodes LAGE-1
        assert canonicalize_antigen("CTAG1A") == "NY-ESO-1"
        assert canonicalize_antigen("CTAG1B") == "NY-ESO-1"
        assert canonicalize_antigen("CTAG2") == "LAGE-1"

        assert resolve_antigen_gene("CTAG1A", species="Homo sapiens") == "CTAG1A"
        assert resolve_antigen_gene("CTAG1B", species="Homo sapiens") == "CTAG1B"
        # When only the protein name NY-ESO-1 is given, paralog ambiguity is explicit
        assert resolve_antigen_gene("NY-ESO-1", species="Homo sapiens") == "CTAG1A;CTAG1B"
        assert resolve_antigen_gene("CTAG2", species="Homo sapiens") == "CTAG2"
        assert resolve_antigen_gene("LAGE-1", species="Homo sapiens") == "CTAG2"
        assert resolve_antigen_gene("LAGE-1", species="Mus musculus") == "Ctag2"
        # Human-specific CTAs/paralogs return None in mouse rather than fake symbols
        assert resolve_antigen_gene("NY-ESO-1", species="Mus musculus") is None
        assert resolve_antigen_gene("CEA", species="Mus musculus") is None
        assert resolve_antigen_gene("PRAME", species="Mus musculus") is None
        assert resolve_antigen_gene("PSA", species="Mus musculus") is None

    def test_endogenous_genes_vs_exogenous_foreign_antigens(self):
        assert resolve_antigen_gene("MART-1", species="Homo sapiens") == "MLANA"
        assert resolve_antigen_gene("MART-1", species="Mus musculus") == "Mlana"
        assert resolve_antigen_gene("Tyrosinase-related protein 1", species="Homo sapiens") == "TYRP1"
        assert resolve_antigen_gene("Tyrosinase-related protein 1", species="Mus musculus") == "Tyrp1"
        assert resolve_antigen_gene("TRP-2", species="Homo sapiens") == "DCT"
        assert resolve_antigen_gene("TRP-2", species="Mus musculus") == "Dct"
        assert resolve_antigen_gene("CEA", species="Homo sapiens") == "CEACAM5"
        assert resolve_antigen_gene("MAGE-A10", species="Homo sapiens") == "MAGEA10"

        # Foreign model antigens (OVA, HEL, KLH), hyphenated protein aliases, and viral proteins -> None
        assert resolve_antigen_gene("Ovalbumin", species="Gallus gallus") is None
        assert resolve_antigen_gene("OVA", species="Gallus gallus") is None
        assert resolve_antigen_gene("HEL", species="Mus musculus") is None
        assert resolve_antigen_gene("KLH", species="Homo sapiens") is None
        assert resolve_antigen_gene("IL-2", species="Homo sapiens") is None
        assert resolve_antigen_gene("pp65", species="CMV") is None

    def test_match_clonotypes_emits_db_antigen_gene_alongside_db_protein_canonical(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAGAGSQGNLIF", "CAVRDNYGQNFVF", "CAASDNYQLIW"],
            "CDR3_beta": ["CASSIRSSYEQYF", "CASSLGGTEAFF", "CASSRANYEQYF"],
        })
        db = pd.DataFrame([
            {
                "cdr3_alpha": "CAGAGSQGNLIF",
                "cdr3_beta": "CASSIRSSYEQYF",
                "epitope": "NLVPMVATV",
                "antigen_gene": "pp65",
                "species": "CMV",
                "host_species": "Homo sapiens",
                "is_viral": True,
                "database": "VDJdb",
            },
            {
                "cdr3_alpha": "CAVRDNYGQNFVF",
                "cdr3_beta": "CASSLGGTEAFF",
                "epitope": "SLLMWITQC",
                "antigen_gene": "NY-ESO-1",
                "species": "Homo sapiens",
                "host_species": "Homo sapiens",
                "is_viral": False,
                "database": "VDJdb",
            },
            {
                "cdr3_alpha": "CAASDNYQLIW",
                "cdr3_beta": "CASSRANYEQYF",
                "epitope": "SIINFEKL",
                "antigen_gene": "Ovalbumin",
                "species": "Gallus gallus",
                "host_species": "Mus musculus",
                "is_viral": False,
                "database": "IEDB",
            },
        ])
        out = match_clonotypes(clones, db, show_progress=False)

        # Row 0: Viral CMV pp65 -> db_protein_canonical="CMV pp65", db_antigen_gene=None
        assert out.loc[0, "db_protein_canonical"] == "CMV pp65"
        assert out.loc[0, "db_antigen_gene"] is None
        assert out.loc[0, "db_category"] == CATEGORY_VIRAL

        # Row 1: Human NY-ESO-1 -> db_protein_canonical="NY-ESO-1", db_antigen_gene="CTAG1A;CTAG1B"
        assert out.loc[1, "db_protein_canonical"] == "NY-ESO-1"
        assert out.loc[1, "db_antigen_gene"] == "CTAG1A;CTAG1B"
        assert out.loc[1, "db_category"] == CATEGORY_TUMOR_SELF

        # Row 2: Chicken OVA -> db_protein_canonical="OVA", db_antigen_gene=None, category="other"
        assert out.loc[2, "db_protein_canonical"] == "OVA"
        assert out.loc[2, "db_antigen_gene"] is None
        assert out.loc[2, "db_category"] == CATEGORY_OTHER


class TestHostSpeciesSymmetryAndIEDBV3:
    """Host-species cross-flagging, same-host filtering, and IEDB v3 anchor/host parsing."""

    def test_host_species_filter_same_vs_all_for_mouse_and_human(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAASDNYQLIW", "CAVDNNNDMRF"],
            "CDR3_beta": ["CASSRANYEQYF", "CASSLGQAYEQYF"],
        })
        db = pd.DataFrame([
            {
                "cdr3_alpha": "CAASDNYQLIW",
                "cdr3_beta": "CASSRANYEQYF",
                "epitope": "SVYDFFVWL",
                "antigen_gene": "Dct",
                "species": "Mus musculus",
                "host_species": "MusMusculus",
                "is_viral": False,
                "database": "VDJdb",
            },
            {
                "cdr3_alpha": "CAVDNNNDMRF",
                "cdr3_beta": "CASSLGQAYEQYF",
                "epitope": "EAAGIGILTV",
                "antigen_gene": "MLANA",
                "species": "Homo sapiens",
                "host_species": "HomoSapiens",
                "is_viral": False,
                "database": "VDJdb",
            },
        ])

        # Querying as mouse with host_species_filter="all": mouse row is plain "ab", human row is "ab_cross"
        mouse_all = match_clonotypes(
            clones, db, species="mouse", host_species_filter="all", show_progress=False
        )
        assert mouse_all.loc[0, "db_match_strength"] == "ab"
        assert mouse_all.loc[0, "db_category"] == CATEGORY_TUMOR_SELF
        assert mouse_all.loc[0, "db_antigen_gene"] == "Dct"
        assert mouse_all.loc[1, "db_match_strength"] == "ab_cross"

        # Querying as mouse with host_species_filter="same": human TCR row is excluded
        mouse_same = match_clonotypes(
            clones, db, species="mouse", host_species_filter="same", show_progress=False
        )
        assert bool(mouse_same.loc[0, "db_match"]) is True
        assert bool(mouse_same.loc[1, "db_match"]) is False

        # Chicken OVA is never misclassified as mouse self
        assert (
            classify_category(
                pd.Series(["Gallus gallus"]),
                pd.Series(["Ovalbumin"]),
                query_species="mouse",
            ).iloc[0]
            == CATEGORY_OTHER
        )

    def test_load_iedb_v3_upgrades_anchorless_cdr3_and_infers_mouse_host(self, tmp_path):
        iedb_csv = tmp_path / "tcr_full_v3.csv"
        content = (
            "Receptor,Receptor,Epitope,Epitope,Epitope,Assay,Host,Reference,"
            "Chain 1,Chain 1,Chain 1,Chain 1,Chain 1,"
            "Chain 2,Chain 2,Chain 2,Chain 2,Chain 2\n"
            "IEDB Receptor ID,Type,Name,Source Molecule,Source Organism,MHC Allele Names,IRI,IEDB IRI,"
            "Type,CDR3 Curated,CDR3 Calculated,Junction Calculated,Curated V Gene,"
            "Type,CDR3 Curated,CDR3 Calculated,Junction Calculated,Calculated V Gene\n"
            "1,alphabeta,SVYDFFVWL,Dopachrome tautomerase,Mus musculus,H2-Kb,"
            "http://ontology.iedb.org/ontology/ONTIE_0000462,http://www.iedb.org/reference/1029999,"
            "alpha,AASDNYQLI,AASDNYQLI,CAASDNYQLIW,TRAV14-1*01,"
            "beta,ASSRANYEQY,ASSRANYEQY,CASSRANYEQYF,TRBV12-1*01\n"
        )
        iedb_csv.write_text(content)
        loaded = load_iedb(iedb_csv)
        assert len(loaded) == 1
        row = loaded.iloc[0]
        assert row["cdr3_alpha"] == "CAASDNYQLIW"
        assert row["cdr3_beta"] == "CASSRANYEQYF"
        assert row["host_species"] == "Mus musculus"
        assert row["v_alpha"] == "TRAV14-1*01"
        assert row["v_beta"] == "TRBV12-1*01"
        assert row["reference"] == "http://www.iedb.org/reference/1029999"

    def test_iedb_junction_fgxg_guard_and_strain_host_and_nodb_schema(self, tmp_path):
        from tcrsift.annotate import (
            _upgrade_iedb_cdr3_with_junction,
            annotate_clonotypes,
            canonicalize_species,
        )

        # Strain-qualified mouse host labels canonicalize to Mus musculus
        assert canonicalize_species("Mus musculus C57BL/6") == "Mus musculus"
        assert canonicalize_species("BALB/c") == "Mus musculus"

        # Junction Calculated with trailing FGXG framework tail is NOT grafted
        upgraded = _upgrade_iedb_cdr3_with_junction(
            pd.Series(["ASSIRSSYEQY", "ASSIRSSYEQY"]),
            pd.Series(["CASSIRSSYEQYF", "CASSIRSSYEQYFGPGTRLTVT"]),
        )
        assert upgraded.iloc[0] == "CASSIRSSYEQYF"
        assert upgraded.iloc[1] == "ASSIRSSYEQY"

        # No-DB short-circuit produces identical annotation columns as DB matching
        clones = pd.DataFrame({"CDR3_alpha": ["CAASDNYQLIW"], "CDR3_beta": ["CASSRANYEQYF"]})
        no_db = annotate_clonotypes(clones)
        with_db = match_clonotypes(
            clones,
            pd.DataFrame([{"cdr3_alpha": "X", "cdr3_beta": "Y", "database": "VDJdb", "is_viral": False}]),
            show_progress=False,
        )
        assert set(no_db.columns) == set(with_db.columns)
