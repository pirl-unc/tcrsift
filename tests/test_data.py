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

"""Tests for data annotation module."""

import numpy as np
import pandas as pd
import pytest

from tcrsift.data import annotate_combined_df


@pytest.fixture
def sample_combined_df():
    """Create sample combined gene expression + VDJ DataFrame."""
    return pd.DataFrame({
        "barcode": ["AAAA-1", "BBBB-1", "CCCC-1", "DDDD-1", "EEEE-1"],
        "sample": ["Sample_1", "Sample_1", "Sample_2", "Sample_2", "Sample_3"],
        "cdr3_aa1": ["CASSL", "CAVSD", None, "CATRD", "CAETF"],
        "cdr3_aa2": ["CASSF", "CASRG", "CASXX", None, "CASYY"],
        "ENSG00000167286.10": [10, 5, 8, 12, 0],  # CD3D
        "ENSG00000198851.10": [8, 4, 6, 10, 0],   # CD3E
        "ENSG00000160654.11": [5, 3, 4, 8, 0],    # CD3G
        "ENSG00000010610.10": [100, 5, 200, 10, 50],  # CD4
        "ENSG00000153563.16": [5, 100, 10, 200, 50],  # CD8A
        "ENSG00000172116.23": [3, 80, 8, 150, 40],    # CD8B
        "nCount_RNA": [5000, 6000, 7000, 8000, 9000],
        "nFeature_RNA": [1000, 1200, 1400, 1600, 1800],
        "percent.mt": [3.0, 4.0, 5.0, 6.0, 7.0],
    })


class TestAnnotateCombinedDf:
    """Tests for annotate_combined_df function."""

    def test_column_renaming(self, sample_combined_df):
        """Test that columns are properly renamed."""
        result = annotate_combined_df(sample_combined_df)

        assert "CD3D" in result.columns
        assert "CD3E" in result.columns
        assert "CD3G" in result.columns
        assert "CD4" in result.columns
        assert "CD8A" in result.columns
        assert "CD8B" in result.columns
        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns

    def test_combined_markers(self, sample_combined_df):
        """Test combined CD3 and CD8 markers."""
        result = annotate_combined_df(sample_combined_df)

        # CD3 = CD3D + CD3E + CD3G
        expected_cd3 = sample_combined_df["ENSG00000167286.10"] + \
                       sample_combined_df["ENSG00000198851.10"] + \
                       sample_combined_df["ENSG00000160654.11"]
        assert list(result["CD3"]) == list(expected_cd3)

        # CD8 = CD8A + CD8B
        expected_cd8 = sample_combined_df["ENSG00000153563.16"] + \
                       sample_combined_df["ENSG00000172116.23"]
        assert list(result["CD8"]) == list(expected_cd8)

    def test_tcr_completeness(self, sample_combined_df):
        """Test TCR completeness flags."""
        result = annotate_combined_df(sample_combined_df)

        # First row has both chains
        assert result.iloc[0]["TCR_complete"] == True
        assert result.iloc[0]["filter:has_alpha"] == True
        assert result.iloc[0]["filter:has_beta"] == True

        # Third row missing alpha
        assert result.iloc[2]["filter:has_alpha"] == False
        assert result.iloc[2]["TCR_complete"] == False

        # Fourth row missing beta
        assert result.iloc[3]["filter:has_beta"] == False
        assert result.iloc[3]["TCR_complete"] == False

    def test_ctaa_identifier(self, sample_combined_df):
        """Test CTaa clone identifier creation."""
        result = annotate_combined_df(sample_combined_df)

        assert result.iloc[0]["CTaa"] == "CASSL_CASSF"
        assert result.iloc[1]["CTaa"] == "CAVSD_CASRG"

    def test_tcell_type_annotation(self, sample_combined_df):
        """Test T cell type annotation."""
        result = annotate_combined_df(sample_combined_df, tcell_min_read_ratio_cd4_vs_cd8=3)

        # Row 0: CD4=100, CD8=8 -> Confident CD4+
        assert "CD4" in result.iloc[0]["Tcell_type"]

        # Row 1: CD4=5, CD8=180 -> Confident CD8+
        assert "CD8" in result.iloc[1]["Tcell_type"]

    def test_filter_columns(self, sample_combined_df):
        """Test filter columns are created."""
        result = annotate_combined_df(sample_combined_df)

        assert "filter:cd3_reads" in result.columns
        assert "filter:percent.mt" in result.columns
        assert "filter:num_reads" in result.columns
        assert "filter:num_genes" in result.columns
        assert "filter:all" in result.columns

    def test_filter_thresholds(self, sample_combined_df):
        """Test filter thresholds are applied correctly."""
        result = annotate_combined_df(
            sample_combined_df,
            min_percent_mt=2.0,
            max_percent_mt=6.0,
            min_genes=1000,
            max_genes=2000,
            min_reads=4000,
            max_reads=10000,
        )

        # Check percent.mt filter
        assert result.iloc[0]["filter:percent.mt"] == True  # 3.0 in range
        assert result.iloc[4]["filter:percent.mt"] == False  # 7.0 > 6.0

    def test_cell_id_creation(self, sample_combined_df):
        """Test Cell_ID is created correctly."""
        result = annotate_combined_df(sample_combined_df)

        assert "Cell_ID" in result.columns
        # Format: Peptide_Number-Barcode
        assert result.iloc[0]["Cell_ID"].endswith("AAAA-1")

    def test_peptide_number_extraction(self, sample_combined_df):
        """Test Peptide_Number is extracted from sample name."""
        result = annotate_combined_df(sample_combined_df)

        assert "Peptide_Number" in result.columns
        assert result.iloc[0]["Peptide_Number"] == "1"  # From Sample_1

    def test_seq_column(self, sample_combined_df):
        """Test Seq column is created."""
        result = annotate_combined_df(sample_combined_df)

        assert "Seq" in result.columns
        assert result.iloc[0]["Seq"] == "CASSL_CASSF"

    def test_cdr3ab_pairs(self, sample_combined_df):
        """Test CDR3a/b pairs column."""
        result = annotate_combined_df(sample_combined_df)

        assert "CDR3a/b" in result.columns
        assert result.iloc[0]["CDR3a/b"] == "CASSL/CASSF"

    def test_both_cd4_and_cd8(self, sample_combined_df):
        """Test Both_CD4_and_CD8 flag."""
        result = annotate_combined_df(sample_combined_df)

        assert "Both_CD4_and_CD8" in result.columns
        # Row 4: CD4=50, CD8=90 -> Both > 1
        assert result.iloc[4]["Both_CD4_and_CD8"] == True

    def test_confident_and_complete(self, sample_combined_df):
        """Test confident_and_complete flag."""
        result = annotate_combined_df(sample_combined_df)

        assert "confident_and_complete" in result.columns
        # Row 0: Confident CD4+ and complete TCR
        assert result.iloc[0]["confident_and_complete"] == True
        # Row 2: Missing alpha chain
        assert result.iloc[2]["confident_and_complete"] == False

    def test_filtered_confident_and_complete(self, sample_combined_df):
        """Test filtered_confident_and_complete flag."""
        result = annotate_combined_df(sample_combined_df)

        assert "filtered_confident_and_complete" in result.columns

    @pytest.mark.xfail(reason="annotate_combined_df doesn't handle empty DataFrames due to pandas type issues")
    def test_empty_dataframe(self):
        """Test with empty DataFrame."""
        # Use explicit dtypes to avoid pandas type conversion issues
        empty_df = pd.DataFrame({
            "barcode": pd.Series([], dtype=str),
            "sample": pd.Series([], dtype=str),
            "cdr3_aa1": pd.Series([], dtype=str),
            "cdr3_aa2": pd.Series([], dtype=str),
            "ENSG00000167286.10": pd.Series([], dtype=float),
            "ENSG00000198851.10": pd.Series([], dtype=float),
            "ENSG00000160654.11": pd.Series([], dtype=float),
            "ENSG00000010610.10": pd.Series([], dtype=float),
            "ENSG00000153563.16": pd.Series([], dtype=float),
            "ENSG00000172116.23": pd.Series([], dtype=float),
            "nCount_RNA": pd.Series([], dtype=float),
            "nFeature_RNA": pd.Series([], dtype=float),
            "percent.mt": pd.Series([], dtype=float),
        })

        result = annotate_combined_df(empty_df)

        assert len(result) == 0

    def test_semicolon_separated_cdr3(self):
        """Test handling of semicolon-separated CDR3 sequences (doublets)."""
        df = pd.DataFrame({
            "barcode": ["AAAA-1"],
            "sample": ["Sample_1"],
            "cdr3_aa1": ["CASSL;CAVSD"],  # Two alpha chains
            "cdr3_aa2": ["CASSF"],
            "ENSG00000167286.10": [10],
            "ENSG00000198851.10": [8],
            "ENSG00000160654.11": [5],
            "ENSG00000010610.10": [100],
            "ENSG00000153563.16": [5],
            "ENSG00000172116.23": [3],
            "nCount_RNA": [5000],
            "nFeature_RNA": [1000],
            "percent.mt": [3.0],
        })

        result = annotate_combined_df(df)

        # CTaa_pairs should contain all combinations
        assert "CTaa_pairs" in result.columns
        pairs = result.iloc[0]["CTaa_pairs"]
        assert "CASSL_CASSF" in pairs
        assert "CAVSD_CASSF" in pairs

    def test_custom_ratio_threshold(self, sample_combined_df):
        """Test custom T cell ratio threshold."""
        result_low = annotate_combined_df(sample_combined_df, tcell_min_read_ratio_cd4_vs_cd8=2)
        result_high = annotate_combined_df(sample_combined_df, tcell_min_read_ratio_cd4_vs_cd8=10)

        # Higher threshold should result in fewer "Confident" annotations
        confident_low = result_low["Tcell_type"].str.contains("Confident").sum()
        confident_high = result_high["Tcell_type"].str.contains("Confident").sum()

        assert confident_low >= confident_high

    def test_likely_tcell_type(self, sample_combined_df):
        """Test 'Likely' T cell type for borderline cases."""
        # Create a case with CD4 > 0 and CD8 = 0
        # For "Likely" instead of "Confident", CD4 must be > 0 and CD8 = 0
        # BUT ratio (1+CD4)/(1+CD8) must be <= threshold (default 3)
        # So we need CD4 = 1 or 2 to get ratio of 2/1 = 2 or 3/1 = 3
        df = pd.DataFrame({
            "barcode": ["AAAA-1"],
            "sample": ["Sample_1"],
            "cdr3_aa1": ["CASSL"],
            "cdr3_aa2": ["CASSF"],
            "ENSG00000167286.10": [10],
            "ENSG00000198851.10": [8],
            "ENSG00000160654.11": [5],
            "ENSG00000010610.10": [2],  # Low CD4: (1+2)/(1+0) = 3, not > 3
            "ENSG00000153563.16": [0],  # No CD8A
            "ENSG00000172116.23": [0],  # No CD8B
            "nCount_RNA": [5000],
            "nFeature_RNA": [1000],
            "percent.mt": [3.0],
        })

        result = annotate_combined_df(df)

        assert "Likely CD4+" in result.iloc[0]["Tcell_type"]
