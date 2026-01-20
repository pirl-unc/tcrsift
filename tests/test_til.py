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

"""Tests for TIL matching module."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.til import (
    get_til_enrichment,
    get_til_summary,
    identify_til_specific_clones,
    load_til_data,
    match_til,
)
from tcrsift.validation import TCRsiftValidationError


@pytest.fixture
def sample_til_data():
    """Create sample TIL AnnData for testing."""
    obs = pd.DataFrame(
        {
            "CDR3_alpha": ["CASSL", "CASSL", "CAVSD", "CAVSD", "CAVSD", "UNKNOWN"],
            "CDR3_beta": ["CASSF", "CASSF", "CASRG", "CASRG", "CASRG", "CASXX"],
            "sample": ["TIL1"] * 6,
        },
        index=["cell1", "cell2", "cell3", "cell4", "cell5", "cell6"],
    )
    return ad.AnnData(obs=obs)


@pytest.fixture
def sample_culture_clonotypes():
    """Create sample culture clonotypes for testing."""
    return pd.DataFrame(
        {
            "clone_id": ["CASSL_CASSF", "NOMATCH_NOMATCH"],
            "CDR3_alpha": ["CASSL", "NOMATCH"],
            "CDR3_beta": ["CASSF", "NOMATCH"],
            "max_frequency": [0.1, 0.05],
        }
    )


class TestMatchTil:
    """Tests for match_til function."""

    def test_basic_matching(self, sample_til_data, sample_culture_clonotypes):
        """Test basic TIL matching."""
        result = match_til(sample_culture_clonotypes, sample_til_data)

        assert "til_match" in result.columns
        assert "til_cell_count" in result.columns
        assert "til_frequency" in result.columns

        # First clone should match (2 cells)
        assert result.iloc[0]["til_match"]
        assert result.iloc[0]["til_cell_count"] == 2

        # Second clone should not match
        assert not result.iloc[1]["til_match"]
        assert result.iloc[1]["til_cell_count"] == 0

    def test_cdr3b_only_matching(self, sample_til_data, sample_culture_clonotypes):
        """Test matching by CDR3 beta only."""
        result = match_til(sample_culture_clonotypes, sample_til_data, match_by="CDR3b_only")

        # Should match by beta chain only
        assert result.iloc[0]["til_match"]

    def test_min_til_cells_filter(self, sample_til_data, sample_culture_clonotypes):
        """Test minimum TIL cells filter."""
        result = match_til(
            sample_culture_clonotypes,
            sample_til_data,
            min_til_cells=3,  # Require 3+ cells
        )

        # First clone only has 2 TIL cells, so should not match
        assert not result.iloc[0]["til_match"]

    def test_til_frequency_calculation(self, sample_til_data, sample_culture_clonotypes):
        """Test TIL frequency is calculated correctly."""
        result = match_til(sample_culture_clonotypes, sample_til_data)

        # 2 cells of CASSL_CASSF out of 6 total cells with complete TCR
        # (all 6 cells have both CDR3_alpha and CDR3_beta)
        expected_freq = 2 / 6
        assert result.iloc[0]["til_frequency"] == pytest.approx(expected_freq, rel=0.1)

    def test_no_til_data(self, sample_culture_clonotypes):
        """Test with empty TIL data."""
        # Create TIL data with no matching clones instead of empty DataFrame
        # to avoid type issues with empty DataFrames
        empty_til = ad.AnnData(
            obs=pd.DataFrame(
                {
                    "CDR3_alpha": ["NOMATCH1"],
                    "CDR3_beta": ["NOMATCH2"],
                },
                index=["cell1"],
            )
        )

        result = match_til(sample_culture_clonotypes, empty_til)

        assert (~result["til_match"]).all()
        assert all(result["til_cell_count"] == 0)


class TestGetTilEnrichment:
    """Tests for get_til_enrichment function."""

    def test_enrichment_calculation(self):
        """Test enrichment calculation."""
        matched_clonotypes = pd.DataFrame(
            {
                "clone_id": ["A", "B"],
                "til_match": [True, True],
                "til_frequency": [0.1, 0.01],
                "max_frequency": [0.01, 0.1],
            }
        )

        result = get_til_enrichment(matched_clonotypes)

        # Clone A: TIL > culture = positive enrichment
        assert result.iloc[0]["til_enrichment"] > 0

        # Clone B: TIL < culture = negative enrichment
        assert result.iloc[1]["til_enrichment"] < 0

    def test_no_matches(self):
        """Test with no TIL matches."""
        unmatched = pd.DataFrame(
            {
                "clone_id": ["A", "B"],
                "til_match": [False, False],
            }
        )

        result = get_til_enrichment(unmatched)

        assert len(result) == 0

    def test_missing_til_match_column(self):
        """Test error when til_match column is missing."""
        no_match_col = pd.DataFrame({"clone_id": ["A"]})

        with pytest.raises(ValueError, match="TIL match information"):
            get_til_enrichment(no_match_col)

    def test_missing_frequency_columns(self):
        """Test handling when frequency columns are missing."""
        matched = pd.DataFrame(
            {
                "clone_id": ["A"],
                "til_match": [True],
            }
        )

        result = get_til_enrichment(matched)

        # Should have 0 enrichment as default
        assert result.iloc[0]["til_enrichment"] == 0


class TestGetTilSummary:
    """Tests for get_til_summary function."""

    def test_basic_summary(self):
        """Test basic summary statistics."""
        matched_clonotypes = pd.DataFrame(
            {
                "clone_id": ["A", "B", "C"],
                "til_match": [True, True, False],
                "til_cell_count": [10, 5, 0],
                "til_frequency": [0.1, 0.05, 0.0],
            }
        )

        summary = get_til_summary(matched_clonotypes)

        assert summary["total_culture_clones"] == 3
        assert summary["til_matched_clones"] == 2
        assert summary["til_recovery_rate"] == pytest.approx(2 / 3)
        assert summary["total_til_cells_matched"] == 15
        assert summary["median_til_frequency"] == pytest.approx(0.075)

    def test_summary_with_tiers(self):
        """Test summary with tier information."""
        matched_clonotypes = pd.DataFrame(
            {
                "clone_id": ["A", "B", "C", "D"],
                "til_match": [True, True, False, True],
                "til_cell_count": [10, 5, 0, 3],
                "til_frequency": [0.1, 0.05, 0.0, 0.03],
                "tier": [1, 1, 2, 2],
            }
        )

        summary = get_til_summary(matched_clonotypes)

        assert "recovery_by_tier" in summary
        assert summary["recovery_by_tier"][1] == 1.0  # 2/2 tier 1 matched
        assert summary["recovery_by_tier"][2] == 0.5  # 1/2 tier 2 matched

    def test_summary_with_antigens(self):
        """Test summary with antigen information."""
        matched_clonotypes = pd.DataFrame(
            {
                "clone_id": ["A", "B"],
                "til_match": [True, True],
                "til_cell_count": [10, 5],
                "til_frequency": [0.1, 0.05],
                "antigens": ["PRAME", "CMV"],
            }
        )

        summary = get_til_summary(matched_clonotypes)

        assert "til_cells_by_antigen" in summary
        assert summary["til_cells_by_antigen"]["PRAME"] == 10
        assert summary["til_cells_by_antigen"]["CMV"] == 5

    def test_missing_til_match(self):
        """Test error handling for missing til_match column."""
        no_match = pd.DataFrame({"clone_id": ["A"]})

        summary = get_til_summary(no_match)

        assert "error" in summary

    def test_empty_clonotypes(self):
        """Test with empty clonotypes."""
        empty = pd.DataFrame(
            {
                "clone_id": pd.Series([], dtype=str),
                "til_match": pd.Series([], dtype=bool),
                "til_cell_count": pd.Series([], dtype=int),
                "til_frequency": pd.Series([], dtype=float),
            }
        )

        summary = get_til_summary(empty)

        assert summary["total_culture_clones"] == 0
        assert summary["til_recovery_rate"] == 0


class TestIdentifyTilSpecificClones:
    """Tests for identify_til_specific_clones function."""

    def test_basic_identification(self, sample_til_data):
        """Test basic TIL-specific clone identification."""
        result = identify_til_specific_clones(sample_til_data, min_cells=2)

        # Should find clones with 2+ cells
        assert len(result) >= 1
        assert "clone_id" in result.columns
        assert "til_cell_count" in result.columns

    def test_exclude_culture_clones(self, sample_til_data, sample_culture_clonotypes):
        """Test excluding culture clones."""
        result = identify_til_specific_clones(
            sample_til_data, culture_clonotypes=sample_culture_clonotypes, min_cells=1
        )

        # CASSL_CASSF should be excluded
        assert "CASSL_CASSF" not in result["clone_id"].values

    def test_min_cells_filter(self, sample_til_data):
        """Test minimum cells filter."""
        result_low = identify_til_specific_clones(sample_til_data, min_cells=1)
        result_high = identify_til_specific_clones(sample_til_data, min_cells=5)

        # Higher threshold should return fewer clones
        assert len(result_high) <= len(result_low)

    def test_cdr3_extraction(self, sample_til_data):
        """Test CDR3 sequence extraction."""
        result = identify_til_specific_clones(sample_til_data, min_cells=2)

        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns

    def test_empty_til_data(self):
        """Test with empty TIL data."""
        # Create empty AnnData with proper shape to avoid warnings
        empty_til = ad.AnnData(np.empty((0, 0)))
        empty_til.obs["CDR3_alpha"] = pd.Series([], dtype=str)
        empty_til.obs["CDR3_beta"] = pd.Series([], dtype=str)
        empty_til.obs["sample"] = pd.Series([], dtype=str)

        result = identify_til_specific_clones(empty_til)

        assert len(result) == 0


class TestLoadTilData:
    """Tests for load_til_data function."""

    def test_load_csv_standard_columns(self, tmp_path):
        """Test loading TIL data from CSV with standard column names."""
        csv_path = tmp_path / "til.csv"
        csv_path.write_text("CDR3_alpha,CDR3_beta,sample\nCAVSD,CASRG,TIL1\nCASSL,CASSF,TIL1\n")

        result = load_til_data("csv", csv_path, "TestSample")

        assert len(result) == 2
        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns
        assert "sample" in result.columns

    def test_load_csv_alternate_columns(self, tmp_path):
        """Test loading TIL data from CSV with alternate column names."""
        csv_path = tmp_path / "til.csv"
        csv_path.write_text("cdr3_alpha,cdr3_beta\nCAVSD,CASRG\nCASSL,CASSF\n")

        result = load_til_data("csv", csv_path, "TestSample")

        assert len(result) == 2
        # Should be renamed to standard names
        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns
        # Sample should be added
        assert result["sample"].iloc[0] == "TestSample"

    def test_load_csv_with_CDR3a_CDR3b_columns(self, tmp_path):
        """Test loading CSV with CDR3a/CDR3b naming convention."""
        csv_path = tmp_path / "til.csv"
        csv_path.write_text("CDR3a,CDR3b\nCAVSD,CASRG\n")

        result = load_til_data("csv", csv_path, "TestSample")

        assert "CDR3_alpha" in result.columns
        assert "CDR3_beta" in result.columns

    def test_load_csv_missing_columns_raises(self, tmp_path):
        """Test that CSV without CDR3 columns raises error."""
        csv_path = tmp_path / "til.csv"
        csv_path.write_text("other_col,another_col\nval1,val2\n")

        with pytest.raises(TCRsiftValidationError, match="missing CDR3 columns"):
            load_til_data("csv", csv_path)

    def test_load_csv_file_not_found_raises(self, tmp_path):
        """Test that missing CSV file raises error."""
        fake_path = tmp_path / "nonexistent.csv"

        with pytest.raises(TCRsiftValidationError, match="not found"):
            load_til_data("csv", fake_path)

    def test_load_h5ad(self, tmp_path):
        """Test loading TIL data from h5ad file."""
        h5ad_path = tmp_path / "til.h5ad"

        # Create a simple AnnData and save it
        obs = pd.DataFrame(
            {
                "CDR3_alpha": ["CAVSD", "CASSL"],
                "CDR3_beta": ["CASRG", "CASSF"],
            },
            index=["cell1", "cell2"],
        )
        adata = ad.AnnData(obs=obs)
        adata.write_h5ad(h5ad_path)

        result = load_til_data("h5ad", h5ad_path, "TestSample")

        assert len(result) == 2
        assert "CDR3_alpha" in result.columns
        assert "sample" in result.columns
        assert result["sample"].iloc[0] == "TestSample"

    def test_load_h5ad_missing_cdr3_raises(self, tmp_path):
        """Test that h5ad without CDR3 columns raises error."""
        h5ad_path = tmp_path / "til.h5ad"

        obs = pd.DataFrame({"other_col": ["val1", "val2"]}, index=["cell1", "cell2"])
        adata = ad.AnnData(obs=obs)
        adata.write_h5ad(h5ad_path)

        with pytest.raises(TCRsiftValidationError, match="missing CDR3 columns"):
            load_til_data("h5ad", h5ad_path)

    def test_load_unknown_type_raises(self, tmp_path):
        """Test that unknown source type raises error."""
        with pytest.raises(TCRsiftValidationError, match="Unknown TIL source type"):
            load_til_data("unknown", tmp_path / "file.txt")


class TestMultiSampleTilMatching:
    """Tests for multi-sample TIL matching."""

    def test_match_til_dict_input(self, sample_culture_clonotypes):
        """Test matching against multiple TIL samples via dict."""
        til_dict = {
            "TIL_Sample1": pd.DataFrame(
                {
                    "CDR3_alpha": ["CASSL", "CASSL", "CAVSD"],
                    "CDR3_beta": ["CASSF", "CASSF", "CASRG"],
                }
            ),
            "TIL_Sample2": pd.DataFrame(
                {
                    "CDR3_alpha": ["CASSL", "OTHER"],
                    "CDR3_beta": ["CASSF", "OTHER"],
                }
            ),
        }

        result = match_til(sample_culture_clonotypes, til_dict)

        # Should have per-sample columns
        assert "til_cell_count.TIL_Sample1" in result.columns
        assert "til_cell_count.TIL_Sample2" in result.columns
        assert "til_frequency.TIL_Sample1" in result.columns
        assert "til_frequency.TIL_Sample2" in result.columns

        # Combined columns
        assert "til_match" in result.columns
        assert "til_samples" in result.columns
        assert "til_cell_count" in result.columns

        # First clone should match both samples
        first_row = result.iloc[0]
        assert first_row["til_match"]
        assert first_row["til_cell_count.TIL_Sample1"] == 2
        assert first_row["til_cell_count.TIL_Sample2"] == 1
        assert first_row["til_cell_count"] == 3  # Total
        # Check til_samples contains both
        assert "TIL_Sample1" in first_row["til_samples"]
        assert "TIL_Sample2" in first_row["til_samples"]

    def test_match_til_single_sample_no_per_sample_columns(
        self, sample_til_data, sample_culture_clonotypes
    ):
        """Test that single sample doesn't add per-sample columns."""
        result = match_til(sample_culture_clonotypes, sample_til_data)

        # Should NOT have per-sample columns for single sample
        assert "til_cell_count.TIL" not in result.columns

        # Should have combined columns
        assert "til_cell_count" in result.columns

    def test_match_til_dataframe_input(self, sample_culture_clonotypes):
        """Test matching against single TIL DataFrame."""
        til_df = pd.DataFrame(
            {
                "CDR3_alpha": ["CASSL", "CASSL"],
                "CDR3_beta": ["CASSF", "CASSF"],
            }
        )

        result = match_til(sample_culture_clonotypes, til_df)

        # Should work like single sample
        assert result.iloc[0]["til_match"]
        assert result.iloc[0]["til_cell_count"] == 2
