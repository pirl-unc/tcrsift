"""
Tests for SCT platform data loading and processing.
"""

import pandas as pd
import pytest

from tcrsift.sct import (
    _find_column,
    aggregate_sct,
    get_sct_specificities,
    load_sct,
)


class TestFindColumn:
    """Tests for _find_column helper function."""

    def test_find_first_matching_column(self):
        """Should return first matching column."""
        df = pd.DataFrame({"tra.CDR3": [1], "other": [2]})
        result = _find_column(df, ["tra.CDR3", "TRA_CDR3"])
        assert result == "tra.CDR3"

    def test_find_second_candidate(self):
        """Should return second candidate if first not found."""
        df = pd.DataFrame({"TRA_CDR3": [1], "other": [2]})
        result = _find_column(df, ["tra.CDR3", "TRA_CDR3"])
        assert result == "TRA_CDR3"

    def test_no_match_returns_none(self):
        """Should return None if no candidates found."""
        df = pd.DataFrame({"foo": [1], "bar": [2]})
        result = _find_column(df, ["tra.CDR3", "TRA_CDR3"])
        assert result is None


class TestLoadSct:
    """Tests for load_sct function."""

    @pytest.fixture
    def basic_sct_df(self):
        """Create a basic SCT-style DataFrame."""
        return pd.DataFrame({
            "tra.CDR3": ["CAVXXX", "CAVYYY", "CAVZZZ"],
            "trb.CDR3": ["CASSAAA", "CASSBBB", "CASSCCC"],
            "tra.readcount": [100, 50, 20],
            "trb.readcount": [150, 60, 30],
            "SNR": [5.0, 3.0, 1.5],
            "Top1.mutation.PE.eq.Top1.mutation.APC": ["Yes", "Yes", "No"],
        })

    @pytest.fixture
    def sct_excel_file(self, temp_dir, basic_sct_df):
        """Create a temporary Excel file with SCT data."""
        path = temp_dir / "sct_data.xlsx"
        basic_sct_df.to_excel(path, sheet_name="Cell", index=False)
        return path

    def test_load_sct_basic(self, sct_excel_file):
        """Test basic loading of SCT data."""
        df = load_sct(sct_excel_file, min_snr=0, min_reads_per_chain=0, require_mutation_match=False)
        assert len(df) == 3
        assert "CDR3_alpha" in df.columns
        assert "CDR3_beta" in df.columns
        assert "CDR3_pair" in df.columns

    def test_load_sct_standardized_columns(self, sct_excel_file):
        """Test that standardized columns are created."""
        df = load_sct(sct_excel_file, min_snr=0, min_reads_per_chain=0, require_mutation_match=False)
        assert df["CDR3_pair"].tolist() == ["CAVXXX/CASSAAA", "CAVYYY/CASSBBB", "CAVZZZ/CASSCCC"]

    def test_load_sct_snr_filter(self, sct_excel_file):
        """Test SNR filtering."""
        df = load_sct(sct_excel_file, min_snr=2.5, min_reads_per_chain=0, require_mutation_match=False)
        assert df["high_quality"].sum() == 2  # SNR 5.0 and 3.0 pass

    def test_load_sct_read_count_filter(self, sct_excel_file):
        """Test read count filtering."""
        df = load_sct(sct_excel_file, min_snr=0, min_reads_per_chain=40, require_mutation_match=False)
        assert df["high_quality"].sum() == 2  # rows with 100/150 and 50/60 pass

    def test_load_sct_mutation_match_filter(self, sct_excel_file):
        """Test mutation match filtering."""
        df = load_sct(sct_excel_file, min_snr=0, min_reads_per_chain=0, require_mutation_match=True)
        assert df["high_quality"].sum() == 2  # "Yes" rows pass

    def test_load_sct_chosen_criteria(self, sct_excel_file):
        """Test stricter 'chosen' criteria."""
        df = load_sct(sct_excel_file, min_snr=0, min_reads_per_chain=0, require_mutation_match=False)
        # chosen requires SNR >= 3.4, reads >= 50, etc.
        assert df["chosen"].sum() <= df["high_quality"].sum()

    def test_load_sct_nonexistent_file(self, temp_dir):
        """Test loading nonexistent file raises error."""
        from tcrsift.validation import TCRsiftValidationError
        with pytest.raises(TCRsiftValidationError):
            load_sct(temp_dir / "nonexistent.xlsx")

    def test_load_sct_missing_cdr3_columns(self, temp_dir):
        """Test loading file without CDR3 columns raises error."""
        from tcrsift.validation import TCRsiftValidationError
        df = pd.DataFrame({"foo": [1, 2], "bar": [3, 4]})
        path = temp_dir / "bad_data.xlsx"
        df.to_excel(path, sheet_name="Cell", index=False)
        with pytest.raises(TCRsiftValidationError, match="Could not find CDR3 columns"):
            load_sct(path)


class TestAggregateSct:
    """Tests for aggregate_sct function."""

    @pytest.fixture
    def sct_df_with_duplicates(self):
        """Create SCT DataFrame with duplicate CDR3 pairs."""
        return pd.DataFrame({
            "CDR3_pair": ["A/B", "A/B", "A/B", "C/D", "C/D"],
            "CDR3_alpha": ["A", "A", "A", "C", "C"],
            "CDR3_beta": ["B", "B", "B", "D", "D"],
            "SNR": [5.0, 4.0, 3.0, 2.0, 1.0],
            "high_quality": [True, True, False, True, False],
            "chosen": [True, False, False, False, False],
            "complete": [True, True, True, True, True],
        })

    def test_aggregate_sct_counts(self, sct_df_with_duplicates):
        """Test aggregation counts entries correctly."""
        result = aggregate_sct(sct_df_with_duplicates, verbose=False)
        assert len(result) == 2
        assert result[result["CDR3_pair"] == "A/B"]["num_original_entries"].values[0] == 3
        assert result[result["CDR3_pair"] == "C/D"]["num_original_entries"].values[0] == 2

    def test_aggregate_sct_boolean_any_all(self, sct_df_with_duplicates):
        """Test boolean aggregation with any/all."""
        result = aggregate_sct(sct_df_with_duplicates, verbose=False)
        ab_row = result[result["CDR3_pair"] == "A/B"].iloc[0]
        # high_quality: [True, True, False] -> any=True, all=False
        assert ab_row["high_quality.any"] == True
        # chosen: [True, False, False] -> any=True, all=False
        assert ab_row["chosen.any"] == True

    def test_aggregate_sct_preserves_key_columns(self, sct_df_with_duplicates):
        """Test that key columns are preserved."""
        result = aggregate_sct(sct_df_with_duplicates, verbose=False)
        assert "CDR3_pair" in result.columns
        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns


class TestGetSctSpecificities:
    """Tests for get_sct_specificities function."""

    def test_extract_specificities_from_mutation_column(self):
        """Test extracting specificities from mutation column."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D", "E/F"],
            "mutation": ["KRAS_G12D", "TP53_R175H", None],
        })
        result = get_sct_specificities(df)
        assert result == {"A/B": "KRAS_G12D", "C/D": "TP53_R175H"}

    def test_extract_specificities_alternate_column(self):
        """Test extracting from Top1.mutation.PE column."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
            "Top1.mutation.PE": ["BRAF_V600E", "NRAS_Q61R"],
        })
        result = get_sct_specificities(df)
        assert result == {"A/B": "BRAF_V600E", "C/D": "NRAS_Q61R"}

    def test_no_specificity_column_returns_empty(self):
        """Test that missing specificity column returns empty dict."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", "C/D"],
            "other_column": [1, 2],
        })
        result = get_sct_specificities(df)
        assert result == {}

    def test_handles_missing_cdr3_pair(self):
        """Test handling of missing CDR3_pair values."""
        df = pd.DataFrame({
            "CDR3_pair": ["A/B", None, "C/D"],
            "mutation": ["X", "Y", "Z"],
        })
        result = get_sct_specificities(df)
        assert "A/B" in result
        assert "C/D" in result
        assert None not in result
