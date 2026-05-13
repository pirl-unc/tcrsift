"""
Tests for TCR annotation using public databases.
"""

import numpy as np
import pandas as pd
import pytest

from tcrsift.annotate import (
    _flag_viral,
    annotate_clonotypes,
    classify_category,
    get_annotation_summary,
    load_cedar,
    load_databases,
    load_iedb,
    load_vdjdb,
    match_clonotypes,
)


class TestFlagViral:
    """Tests for viral flagging."""

    def test_flag_cmv_viral(self):
        """CMV should be flagged as viral."""
        df = pd.DataFrame({"species": ["CMV", "Cytomegalovirus", "Human CMV"]})
        result = _flag_viral(df)
        assert all(result)

    def test_flag_ebv_viral(self):
        """EBV should be flagged as viral."""
        df = pd.DataFrame({"species": ["EBV", "Epstein-Barr virus"]})
        result = _flag_viral(df)
        assert all(result)

    def test_flag_hiv_viral(self):
        """HIV should be flagged as viral."""
        df = pd.DataFrame({"species": ["HIV", "Human Immunodeficiency Virus"]})
        result = _flag_viral(df)
        assert all(result)

    def test_non_viral_not_flagged(self):
        """Non-viral species should not be flagged."""
        df = pd.DataFrame({"species": ["Homo sapiens", "Mouse", "Self-antigen"]})
        result = _flag_viral(df)
        assert not any(result)

    def test_missing_species_not_flagged(self):
        """Missing species column should not be flagged."""
        df = pd.DataFrame({"epitope": ["PEPTIDE1", "PEPTIDE2"]})
        result = _flag_viral(df)
        assert not any(result)

    def test_null_species_not_flagged(self):
        """Null species should not be flagged."""
        df = pd.DataFrame({"species": [None, np.nan, ""]})
        result = _flag_viral(df)
        assert not any(result)


class TestClassifyCategory:
    """Tests for the species/antigen → category classifier."""

    def test_viral(self):
        cats = classify_category(
            pd.Series(["Human cytomegalovirus", "Epstein-Barr virus"]),
            pd.Series(["pp65", "EBNA-3"]),
        )
        assert (cats == "viral").all()

    def test_bacterial(self):
        cats = classify_category(
            pd.Series(["Mycobacterium tuberculosis", "Listeria monocytogenes"]),
            pd.Series(["ESAT-6", "LLO"]),
        )
        assert (cats == "bacterial").all()

    def test_self_homo_sapiens(self):
        cats = classify_category(
            pd.Series(["Homo sapiens (human)", "Homo sapiens"]),
            pd.Series(["beta-2-microglobulin", "insulin"]),
        )
        assert (cats == "self").all()

    def test_tumor_self_overrides_species(self):
        """MART-1 / NY-ESO-1 / MAGE peptides are Homo sapiens but should
        not bucket as plain self — they're tumor-associated self antigens."""
        cats = classify_category(
            pd.Series(["Homo sapiens", "Homo sapiens", "Homo sapiens"]),
            pd.Series(["MLANA (MART-1)", "NY-ESO-1", "MAGE-A3"]),
        )
        assert (cats == "tumor_self").all()

    def test_unknown_species(self):
        cats = classify_category(
            pd.Series(["", None]),
            pd.Series(["", None]),
        )
        assert (cats == "unknown").all()

    def test_other_known_but_unclassified(self):
        cats = classify_category(
            pd.Series(["Plasmodium falciparum"]),
            pd.Series(["CSP"]),
        )
        assert (cats == "other").all()


class TestLoadVdjdb:
    """Tests for loading VDJdb."""

    @pytest.fixture
    def mock_vdjdb_file(self, temp_dir):
        """Create a mock VDJdb file."""
        vdjdb_path = temp_dir / "vdjdb.tsv"
        df = pd.DataFrame(
            {
                "cdr3": ["CASSLGQAYEQYF", "CASSXYZAYEQYF"],
                "cdr3.alpha": ["CAVSDGGSQGNLIF", "CAVXYZQGNLIF"],
                "antigen.epitope": ["NLV", "GLC"],
                "antigen.species": ["CMV", "EBV"],
                "mhc.a": ["HLA-A*02:01", "HLA-B*08:01"],
            }
        )
        df.to_csv(vdjdb_path, sep="\t", index=False)
        return vdjdb_path

    @pytest.fixture
    def mock_vdjdb_dir(self, temp_dir, mock_vdjdb_file):
        """Create a mock VDJdb directory."""
        return temp_dir

    def test_load_from_file(self, mock_vdjdb_file):
        """Load VDJdb from file."""
        result = load_vdjdb(mock_vdjdb_file)

        assert len(result) == 2
        assert "cdr3_beta" in result.columns
        assert "cdr3_alpha" in result.columns
        assert "epitope" in result.columns
        assert "is_viral" in result.columns
        assert result["database"].iloc[0] == "VDJdb"

    def test_load_from_directory(self, mock_vdjdb_dir):
        """Load VDJdb from directory."""
        result = load_vdjdb(mock_vdjdb_dir)
        assert len(result) == 2

    def test_viral_flagging(self, mock_vdjdb_file):
        """VDJdb entries should have viral flag."""
        result = load_vdjdb(mock_vdjdb_file)
        # CMV and EBV are both viral
        assert all(result["is_viral"])


class TestLoadIedb:
    """Tests for loading IEDB."""

    @pytest.fixture
    def mock_iedb_file(self, temp_dir):
        """Create a mock IEDB file."""
        iedb_path = temp_dir / "iedb.tsv"
        df = pd.DataFrame(
            {
                "Chain 2 CDR3 Curated": ["CASSLGQAYEQYF"],
                "Chain 1 CDR3 Curated": ["CAVSDGGSQGNLIF"],
                "Epitope - Name": ["pp65"],
                "Epitope - Source Organism Name": ["CMV"],
            }
        )
        df.to_csv(iedb_path, sep="\t", index=False)
        return iedb_path

    def test_load_iedb(self, mock_iedb_file):
        """Load IEDB file."""
        result = load_iedb(mock_iedb_file)

        assert len(result) == 1
        assert "cdr3_beta" in result.columns
        assert result["database"].iloc[0] == "IEDB"

    def test_load_iedb_v3_hierarchical_header(self, temp_dir):
        """Parse the current v3 export shape: comma-separated, two-row header.

        The header starts ``Receptor,Receptor,Receptor,...`` and fields
        live one level down (``CDR3 Curated`` under ``Chain 1`` etc.).
        Synthesised here as a small CSV so the test doesn't depend on
        a 99MB cache file.
        """
        path = temp_dir / "tcr_full_v3.csv"
        content = (
            # Row 1 — section names (repeated per field)
            "Receptor,Receptor,Epitope,Epitope,Epitope,Assay,"
            "Chain 1,Chain 1,Chain 2,Chain 2\n"
            # Row 2 — field names
            "IEDB Receptor ID,Type,Name,Source Molecule,Source Organism,"
            "MHC Allele Names,Type,CDR3 Curated,Type,CDR3 Curated\n"
            # Data — one αβ row that should be picked up
            "1,alphabeta,NLVPMVATV,pp65,Human cytomegalovirus,HLA-A*02:01,"
            "alpha,CAVSDGGSQGNLIF,beta,CASSLGQAYEQYF\n"
            # A gammadelta row that should be dropped (not αβ-matchable)
            "2,gammadelta,SOMEPEP,foo,Homo sapiens,HLA-A*02:01,"
            "gamma,CAVGAMMA,delta,CASSDELTA\n"
        )
        path.write_text(content)

        result = load_iedb(path)

        assert len(result) == 1
        row = result.iloc[0]
        assert row["cdr3_alpha"] == "CAVSDGGSQGNLIF"
        assert row["cdr3_beta"] == "CASSLGQAYEQYF"
        assert row["epitope"] == "NLVPMVATV"
        assert row["antigen_gene"] == "pp65"
        assert row["species"] == "Human cytomegalovirus"
        assert row["mhc_allele"] == "HLA-A*02:01"
        assert bool(row["is_viral"]) is True
        assert row["database"] == "IEDB"

    def test_load_iedb_v3_falls_back_to_calculated_cdr3(self, temp_dir):
        """If a chain has no ``CDR3 Curated`` value but ``CDR3 Calculated``
        is present, the calculated value is used. Mirrors the priority
        VDJdb applies to its curated-over-inferred fields."""
        path = temp_dir / "tcr_full_v3.csv"
        content = (
            "Receptor,Receptor,Epitope,Chain 1,Chain 1,Chain 1,"
            "Chain 2,Chain 2,Chain 2\n"
            "IEDB Receptor ID,Type,Name,Type,CDR3 Curated,CDR3 Calculated,"
            "Type,CDR3 Curated,CDR3 Calculated\n"
            "1,alphabeta,EPI1,alpha,,CAVCALCALPHA,beta,CASSCURBETA,CASSCALCBETA\n"
        )
        path.write_text(content)

        result = load_iedb(path)

        assert len(result) == 1
        assert result.iloc[0]["cdr3_alpha"] == "CAVCALCALPHA"
        # Curated wins over Calculated when both are present.
        assert result.iloc[0]["cdr3_beta"] == "CASSCURBETA"


class TestLoadCedar:
    """Tests for loading CEDAR."""

    @pytest.fixture
    def mock_cedar_file(self, temp_dir):
        """Create a mock CEDAR file."""
        cedar_path = temp_dir / "cedar.tsv"
        df = pd.DataFrame(
            {
                "cdr3_b_aa": ["CASSLGQAYEQYF"],
                "cdr3_a_aa": ["CAVSDGGSQGNLIF"],
                "epitope_sequence": ["NLVPMVATV"],
                "organism": ["CMV"],
            }
        )
        df.to_csv(cedar_path, sep="\t", index=False)
        return cedar_path

    def test_load_cedar(self, mock_cedar_file):
        """Load CEDAR file."""
        result = load_cedar(mock_cedar_file)

        assert len(result) == 1
        assert "cdr3_beta" in result.columns
        assert result["database"].iloc[0] == "CEDAR"


class TestLoadDatabases:
    """Tests for combined database loading."""

    def test_load_combined(self, temp_dir):
        """Load multiple databases."""
        # Create mock files
        vdjdb_path = temp_dir / "vdjdb.tsv"
        iedb_path = temp_dir / "iedb.tsv"

        pd.DataFrame(
            {
                "cdr3": ["CASSTEST1"],
                "cdr3.alpha": ["CAVTEST1"],
                "antigen.species": ["CMV"],
            }
        ).to_csv(vdjdb_path, sep="\t", index=False)

        pd.DataFrame(
            {
                "Chain 2 CDR3 Curated": ["CASSTEST2"],
                "Chain 1 CDR3 Curated": ["CAVTEST2"],
                "Epitope - Source Organism Name": ["EBV"],
            }
        ).to_csv(iedb_path, sep="\t", index=False)

        result = load_databases(vdjdb_path=vdjdb_path, iedb_path=iedb_path)

        assert len(result) == 2
        assert set(result["database"].unique()) == {"VDJdb", "IEDB"}

    def test_load_raises_without_any_db(self):
        """Should raise if no database provided."""
        with pytest.raises(ValueError, match="At least one database"):
            load_databases()


class TestMatchClonotypes:
    """Tests for clonotype matching."""

    def test_match_by_cdr3ab(self, sample_clonotypes_df, sample_database_df):
        """Match by both alpha and beta chains."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3ab",
        )

        assert "db_match" in result.columns
        assert "db_epitope" in result.columns
        assert "is_viral" in result.columns

        # First clone should match (same CDR3ab)
        assert result["db_match"].iloc[0]

    def test_match_by_cdr3b_only(self, sample_clonotypes_df, sample_database_df):
        """Match by beta chain only."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3b_only",
        )

        # More matches expected with beta-only matching
        assert result["db_match"].sum() >= 1

    def test_no_match_annotation(self, sample_clonotypes_df):
        """Clonotypes without matches should have empty annotations."""
        # Database with non-matching sequences
        non_matching_db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSNOMATCH"],
                "cdr3_alpha": ["CAVNOMATCH"],
                "epitope": ["NOTFOUND"],
                "species": ["Unknown"],
                "database": ["TestDB"],
                "is_viral": [False],
            }
        )

        result = match_clonotypes(sample_clonotypes_df, non_matching_db, match_by="CDR3ab")

        assert result["db_match"].sum() == 0
        assert (~result["is_viral"]).all()

    def test_skips_with_warning_when_database_missing_cdr3_columns(
        self, sample_clonotypes_df, caplog
    ):
        """When a database loader silently produces a DataFrame without
        `cdr3_alpha`/`cdr3_beta` (e.g. format drift in IEDB's v3 CSV
        export — the column-mapping in load_iedb expects flat headers
        but v3 ships hierarchical ones), match_clonotypes should log a
        clear warning and return the clonotypes unannotated rather than
        crashing with KeyError (#46)."""
        import logging

        # Simulate a load_iedb output where the column-rename map missed
        # every header — the result is a DataFrame with raw IEDB headers
        # and no `cdr3_alpha`/`cdr3_beta`.
        malformed_db = pd.DataFrame(
            {
                "Chain 1": ["foo"],
                "Chain 2": ["bar"],
                "database": ["IEDB"],
                "is_viral": [False],
            }
        )

        with caplog.at_level(logging.WARNING, logger="tcrsift.annotate"):
            result = match_clonotypes(
                sample_clonotypes_df, malformed_db, match_by="CDR3ab"
            )

        # Annotation columns are present and empty (initialized then
        # returned without matching).
        assert "db_match" in result.columns
        assert not result["db_match"].any()
        # Warning names the missing columns and points at the cause.
        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("missing required columns" in m for m in warnings)
        assert any("cdr3_alpha" in m for m in warnings)

    def test_structured_columns_present(self, sample_clonotypes_df, sample_database_df):
        """match_clonotypes always emits the structured-annotation columns
        added in #48 (db_protein, db_mhc, db_category, db_match_strength)."""
        result = match_clonotypes(
            sample_clonotypes_df, sample_database_df, match_by="CDR3ab"
        )
        for col in ("db_protein", "db_mhc", "db_category", "db_match_strength"):
            assert col in result.columns

    def test_match_strictness_strict_ab_no_fallback(self, sample_clonotypes_df):
        """strict_ab should *not* fall back to β-only matching. A clone
        whose β matches but α doesn't should come back unmatched."""
        # DB has the β CDR3 but a different α CDR3 → only β matches.
        db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "cdr3_alpha": ["CAVDIFFERENT"],
                "epitope": ["NLV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )

        strict = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="strict_ab"
        )
        # Legacy ab_with_partial would mark this clone as matched (partial);
        # strict_ab must not.
        assert strict["db_match"].sum() == 0
        assert not strict["db_match_partial"].any()

        partial = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="ab_with_partial"
        )
        # The legacy mode picks up the β-only fallback.
        assert partial["db_match"].sum() >= 1
        assert partial["db_match_partial"].sum() >= 1
        # Partial matches are labeled b_only in db_match_strength.
        matched = partial[partial["db_match"]]
        assert (matched["db_match_strength"] == "b_only").all()

    def test_match_strictness_b_only_matches_strength(self, sample_clonotypes_df):
        """b_only strictness should mark all matches as db_match_strength=b_only,
        without setting db_match_partial (that flag is reserved for the αβ
        fallback path, not an explicit β-only run)."""
        db = pd.DataFrame(
            {
                "cdr3_beta": ["CASSLGQAYEQYF"],
                "cdr3_alpha": [None],
                "epitope": ["NLV"],
                "antigen_gene": ["pp65"],
                "species": ["Human cytomegalovirus"],
                "mhc_allele": ["HLA-A*02:01"],
                "database": ["VDJdb"],
                "is_viral": [True],
            }
        )

        result = match_clonotypes(
            sample_clonotypes_df, db, match_strictness="b_only"
        )
        matched = result[result["db_match"]]
        assert len(matched) >= 1
        assert (matched["db_match_strength"] == "b_only").all()
        # Not a fallback — explicit β-only mode.
        assert not matched["db_match_partial"].any()

    def test_match_strictness_invalid_raises(self, sample_clonotypes_df, sample_database_df):
        """An unknown strictness should produce a clear validation error."""
        from tcrsift.validation import TCRsiftValidationError

        with pytest.raises(TCRsiftValidationError, match="match_strictness"):
            match_clonotypes(
                sample_clonotypes_df,
                sample_database_df,
                match_strictness="bogus",
            )

    def test_db_category_populated_from_species_and_protein(self, sample_clonotypes_df):
        """db_category should derive from species + antigen, with tumor-self
        antigens (MART-1) overriding the species-derived bucket."""
        # Two DB entries: one viral (CMV/pp65), one tumor-self (Homo sapiens,
        # but MART-1 should override "self" → "tumor_self").
        sample_a = sample_clonotypes_df.iloc[0]
        sample_b = sample_clonotypes_df.iloc[1] if len(sample_clonotypes_df) > 1 else sample_a
        db = pd.DataFrame(
            {
                "cdr3_beta": [sample_a["CDR3_beta"], sample_b["CDR3_beta"]],
                "cdr3_alpha": [sample_a["CDR3_alpha"], sample_b["CDR3_alpha"]],
                "epitope": ["NLV", "ELAGIGILTV"],
                "antigen_gene": ["pp65", "MLANA"],
                "species": ["Human cytomegalovirus", "Homo sapiens"],
                "mhc_allele": ["HLA-A*02:01", "HLA-A*02:01"],
                "database": ["VDJdb", "IEDB"],
                "is_viral": [True, False],
            }
        )

        result = match_clonotypes(sample_clonotypes_df, db, match_by="CDR3ab")
        matched = result[result["db_match"]]
        assert len(matched) >= 1
        # Categories should be assigned per the curated tables.
        cats = set(matched["db_category"].dropna())
        assert cats.issubset({"viral", "tumor_self"})

    def test_viral_flag_propagation(self, sample_clonotypes_df, sample_database_df):
        """Viral flag should propagate from database."""
        result = match_clonotypes(
            sample_clonotypes_df,
            sample_database_df,
            match_by="CDR3ab",
        )

        # Matched clones should inherit viral status
        matched = result[result["db_match"]]
        if len(matched) > 0:
            # All matches in our test db are viral
            assert all(matched["is_viral"])


class TestAnnotateClonotypes:
    """Tests for main annotation function."""

    @pytest.fixture
    def mock_db_paths(self, temp_dir):
        """Create mock database files."""
        vdjdb_path = temp_dir / "vdjdb.tsv"
        pd.DataFrame(
            {
                "cdr3": ["CASSLGQAYEQYF"],
                "cdr3.alpha": ["CAVSDGGSQGNLIF"],
                "antigen.epitope": ["NLV"],
                "antigen.species": ["CMV"],
            }
        ).to_csv(vdjdb_path, sep="\t", index=False)

        return {"vdjdb_path": vdjdb_path}

    def test_annotate_clonotypes(self, sample_clonotypes_df, mock_db_paths):
        """Full annotation pipeline."""
        result = annotate_clonotypes(
            sample_clonotypes_df,
            vdjdb_path=mock_db_paths["vdjdb_path"],
        )

        assert "db_match" in result.columns
        assert len(result) == len(sample_clonotypes_df)

    def test_exclude_viral(self, sample_clonotypes_df, mock_db_paths):
        """Exclude viral clones option."""
        result = annotate_clonotypes(
            sample_clonotypes_df,
            vdjdb_path=mock_db_paths["vdjdb_path"],
            exclude_viral=True,
        )

        # Viral clones should be excluded
        assert not any(result.get("is_viral", pd.Series([False])))

    def test_no_database_paths_returns_default_annotation_columns(self, sample_clonotypes_df):
        """Annotation should be optional when no database paths are provided."""
        result = annotate_clonotypes(sample_clonotypes_df)

        assert len(result) == len(sample_clonotypes_df)
        assert "db_match" in result.columns
        assert "db_epitope" in result.columns
        assert "db_species" in result.columns
        assert "db_database" in result.columns
        assert "is_viral" in result.columns
        assert (~result["db_match"]).all()
        assert (~result["is_viral"]).all()


class TestGetAnnotationSummary:
    """Tests for annotation summary."""

    def test_summary_basic(self, sample_clonotypes_df, sample_database_df):
        """Summary should have basic stats."""
        annotated = match_clonotypes(sample_clonotypes_df, sample_database_df)
        summary = get_annotation_summary(annotated)

        assert "total" in summary
        assert "matched" in summary
        assert "viral" in summary

    def test_summary_database_breakdown(self, sample_clonotypes_df, sample_database_df):
        """Summary should have database breakdown."""
        annotated = match_clonotypes(sample_clonotypes_df, sample_database_df)
        summary = get_annotation_summary(annotated)

        if "database_breakdown" in summary:
            for db in ["VDJdb", "IEDB", "CEDAR"]:
                assert db in summary["database_breakdown"]
