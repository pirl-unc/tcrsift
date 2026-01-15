"""
Tests for TCR annotation using public databases.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from tcrsift.annotate import (
    load_vdjdb,
    load_iedb,
    load_cedar,
    load_databases,
    match_clonotypes,
    annotate_clonotypes,
    get_annotation_summary,
    VIRAL_SPECIES_PATTERNS,
    _flag_viral,
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


class TestLoadVdjdb:
    """Tests for loading VDJdb."""

    @pytest.fixture
    def mock_vdjdb_file(self, temp_dir):
        """Create a mock VDJdb file."""
        vdjdb_path = temp_dir / "vdjdb.tsv"
        df = pd.DataFrame({
            "cdr3": ["CASSLGQAYEQYF", "CASSXYZAYEQYF"],
            "cdr3.alpha": ["CAVSDGGSQGNLIF", "CAVXYZQGNLIF"],
            "antigen.epitope": ["NLV", "GLC"],
            "antigen.species": ["CMV", "EBV"],
            "mhc.a": ["HLA-A*02:01", "HLA-B*08:01"],
        })
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
        df = pd.DataFrame({
            "Chain 2 CDR3 Curated": ["CASSLGQAYEQYF"],
            "Chain 1 CDR3 Curated": ["CAVSDGGSQGNLIF"],
            "Epitope - Name": ["pp65"],
            "Epitope - Source Organism Name": ["CMV"],
        })
        df.to_csv(iedb_path, sep="\t", index=False)
        return iedb_path

    def test_load_iedb(self, mock_iedb_file):
        """Load IEDB file."""
        result = load_iedb(mock_iedb_file)

        assert len(result) == 1
        assert "cdr3_beta" in result.columns
        assert result["database"].iloc[0] == "IEDB"


class TestLoadCedar:
    """Tests for loading CEDAR."""

    @pytest.fixture
    def mock_cedar_file(self, temp_dir):
        """Create a mock CEDAR file."""
        cedar_path = temp_dir / "cedar.tsv"
        df = pd.DataFrame({
            "cdr3_b_aa": ["CASSLGQAYEQYF"],
            "cdr3_a_aa": ["CAVSDGGSQGNLIF"],
            "epitope_sequence": ["NLVPMVATV"],
            "organism": ["CMV"],
        })
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

        pd.DataFrame({
            "cdr3": ["CASSTEST1"],
            "cdr3.alpha": ["CAVTEST1"],
            "antigen.species": ["CMV"],
        }).to_csv(vdjdb_path, sep="\t", index=False)

        pd.DataFrame({
            "Chain 2 CDR3 Curated": ["CASSTEST2"],
            "Chain 1 CDR3 Curated": ["CAVTEST2"],
            "Epitope - Source Organism Name": ["EBV"],
        }).to_csv(iedb_path, sep="\t", index=False)

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
        assert result["db_match"].iloc[0] == True

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
        non_matching_db = pd.DataFrame({
            "cdr3_beta": ["CASSNOMATCH"],
            "cdr3_alpha": ["CAVNOMATCH"],
            "epitope": ["NOTFOUND"],
            "species": ["Unknown"],
            "database": ["TestDB"],
            "is_viral": [False],
        })

        result = match_clonotypes(sample_clonotypes_df, non_matching_db, match_by="CDR3ab")

        assert result["db_match"].sum() == 0
        assert all(result["is_viral"] == False)

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
        pd.DataFrame({
            "cdr3": ["CASSLGQAYEQYF"],
            "cdr3.alpha": ["CAVSDGGSQGNLIF"],
            "antigen.epitope": ["NLV"],
            "antigen.species": ["CMV"],
        }).to_csv(vdjdb_path, sep="\t", index=False)

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
