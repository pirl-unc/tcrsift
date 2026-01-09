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

"""Tests for loader module."""

import pandas as pd
import pytest

from tcrsift.loader import (
    VDJ_SEGMENT_COLS,
    VDJ_SEGMENT_NT_COLS,
    _pivot_vdj_by_barcode,
)


class TestPivotVdjByBarcode:
    """Tests for _pivot_vdj_by_barcode function."""

    def test_basic_pivot(self):
        """Test basic pivoting of VDJ data."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA", "BBBB", "BBBB"],
            "chain": ["TRA", "TRB", "TRA", "TRB"],
            "cdr3": ["CASSL", "CASSF", "CAVSD", "CASRG"],
            "v_gene": ["TRAV1", "TRBV2", "TRAV3", "TRBV4"],
            "d_gene": [None, "TRBD1", None, "TRBD2"],
            "j_gene": ["TRAJ1", "TRBJ2", "TRAJ3", "TRBJ4"],
            "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC2"],
            "umis": [100, 200, 150, 250],
            "reads": [1000, 2000, 1500, 2500],
            "contig_id": ["c1", "c2", "c3", "c4"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        assert len(result) == 2
        assert "AAAA" in result.index
        assert "BBBB" in result.index
        assert result.loc["AAAA", "TRA_1_cdr3"] == "CASSL"
        assert result.loc["AAAA", "TRB_1_cdr3"] == "CASSF"
        assert result.loc["BBBB", "TRA_1_cdr3"] == "CAVSD"
        assert result.loc["BBBB", "TRB_1_cdr3"] == "CASRG"

    def test_doublet_handling(self):
        """Test handling of cells with multiple chains (doublets)."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA", "AAAA", "AAAA"],
            "chain": ["TRA", "TRA", "TRB", "TRB"],
            "cdr3": ["CASSL1", "CASSL2", "CASSF1", "CASSF2"],
            "v_gene": ["TRAV1", "TRAV2", "TRBV1", "TRBV2"],
            "d_gene": [None, None, "TRBD1", "TRBD2"],
            "j_gene": ["TRAJ1", "TRAJ2", "TRBJ1", "TRBJ2"],
            "c_gene": ["TRAC", "TRAC", "TRBC1", "TRBC2"],
            "umis": [100, 50, 200, 75],
            "reads": [1000, 500, 2000, 750],
            "contig_id": ["c1", "c2", "c3", "c4"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        assert len(result) == 1
        assert "AAAA" in result.index
        # Primary chain should be highest UMI
        assert result.loc["AAAA", "TRA_1_cdr3"] == "CASSL1"
        assert result.loc["AAAA", "TRA_2_cdr3"] == "CASSL2"
        assert result.loc["AAAA", "TRB_1_cdr3"] == "CASSF1"
        assert result.loc["AAAA", "TRB_2_cdr3"] == "CASSF2"
        # Doublet flags
        assert result.loc["AAAA", "multi_TRA"] == True
        assert result.loc["AAAA", "multi_TRB"] == True
        assert result.loc["AAAA", "multi_chain"] == True

    def test_umi_prioritization(self):
        """Test that chains are prioritized by UMI count."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA"],
            "chain": ["TRA", "TRA"],
            "cdr3": ["LOW_UMI", "HIGH_UMI"],
            "v_gene": ["TRAV1", "TRAV2"],
            "d_gene": [None, None],
            "j_gene": ["TRAJ1", "TRAJ2"],
            "c_gene": ["TRAC", "TRAC"],
            "umis": [10, 100],
            "reads": [100, 1000],
            "contig_id": ["c1", "c2"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        # Higher UMI should be primary chain
        assert result.loc["AAAA", "TRA_1_cdr3"] == "HIGH_UMI"
        assert result.loc["AAAA", "TRA_2_cdr3"] == "LOW_UMI"

    def test_chain_count_columns(self):
        """Test that chain count columns are properly created."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA", "BBBB"],
            "chain": ["TRA", "TRB", "TRB"],
            "cdr3": ["CASSL", "CASSF", "CASRG"],
            "v_gene": ["TRAV1", "TRBV2", "TRBV4"],
            "d_gene": [None, "TRBD1", "TRBD2"],
            "j_gene": ["TRAJ1", "TRBJ2", "TRBJ4"],
            "c_gene": ["TRAC", "TRBC1", "TRBC2"],
            "umis": [100, 200, 250],
            "reads": [1000, 2000, 2500],
            "contig_id": ["c1", "c2", "c3"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        # AAAA has both chains
        assert result.loc["AAAA", "TRA_count"] == 1
        assert result.loc["AAAA", "TRB_count"] == 1
        assert result.loc["AAAA", "has_TRA"] == True
        assert result.loc["AAAA", "has_TRB"] == True
        assert result.loc["AAAA", "has_both_chains"] == True

        # BBBB only has TRB
        assert result.loc["BBBB", "TRA_count"] == 0
        assert result.loc["BBBB", "TRB_count"] == 1
        assert result.loc["BBBB", "has_TRA"] == False
        assert result.loc["BBBB", "has_TRB"] == True
        assert result.loc["BBBB", "has_both_chains"] == False

    def test_cdr3ab_identifier(self):
        """Test creation of CDR3ab identifier."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA"],
            "chain": ["TRA", "TRB"],
            "cdr3": ["CASSL", "CASSF"],
            "v_gene": ["TRAV1", "TRBV2"],
            "d_gene": [None, "TRBD1"],
            "j_gene": ["TRAJ1", "TRBJ2"],
            "c_gene": ["TRAC", "TRBC1"],
            "umis": [100, 200],
            "reads": [1000, 2000],
            "contig_id": ["c1", "c2"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        assert result.loc["AAAA", "CDR3_alpha"] == "CASSL"
        assert result.loc["AAAA", "CDR3_beta"] == "CASSF"
        assert result.loc["AAAA", "CDR3ab"] == "CASSL_CASSF"

    def test_segment_preservation(self):
        """Test that VDJ segment columns are preserved in pivot."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA"],
            "chain": ["TRA", "TRB"],
            "cdr3": ["CASSL", "CASSF"],
            "v_gene": ["TRAV1", "TRBV2"],
            "d_gene": [None, "TRBD1"],
            "j_gene": ["TRAJ1", "TRBJ2"],
            "c_gene": ["TRAC", "TRBC1"],
            "umis": [100, 200],
            "reads": [1000, 2000],
            "contig_id": ["c1", "c2"],
            "fwr1": ["MRLV", "MGVT"],
            "cdr1": ["TSGF", "SGHD"],
            "fwr2": ["WYRQ", "WYQQ"],
            "cdr2": ["YSSG", "SNNE"],
            "fwr3": ["GKAP", "GKGP"],
            "fwr4": ["FGGG", "FGXG"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        # Check that segment columns are present
        assert "TRA_1_fwr1" in result.columns
        assert "TRA_1_cdr1" in result.columns
        assert "TRB_1_fwr1" in result.columns
        assert result.loc["AAAA", "TRA_1_fwr1"] == "MRLV"
        assert result.loc["AAAA", "TRB_1_cdr1"] == "SGHD"

    def test_vdj_sequence_preservation(self):
        """Test that combined VDJ sequences are preserved."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA", "AAAA"],
            "chain": ["TRA", "TRB"],
            "cdr3": ["CASSL", "CASSF"],
            "v_gene": ["TRAV1", "TRBV2"],
            "d_gene": [None, "TRBD1"],
            "j_gene": ["TRAJ1", "TRBJ2"],
            "c_gene": ["TRAC", "TRBC1"],
            "umis": [100, 200],
            "reads": [1000, 2000],
            "contig_id": ["c1", "c2"],
            "vdj_aa": ["MRLVTSGFWYRQYSSGCASSL", "MGVTSGHDFGGG"],
            "vdj_nt": ["ATGCGT...", "ATGGGT..."],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        assert "TRA_1_vdj_aa" in result.columns
        assert "TRB_1_vdj_nt" in result.columns
        assert result.loc["AAAA", "TRA_1_vdj_aa"] == "MRLVTSGFWYRQYSSGCASSL"

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        vdj_df = pd.DataFrame({
            "barcode": [],
            "chain": [],
            "cdr3": [],
            "v_gene": [],
            "d_gene": [],
            "j_gene": [],
            "c_gene": [],
            "umis": [],
            "reads": [],
            "contig_id": [],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        assert len(result) == 0

    def test_more_than_two_chains_filtered(self):
        """Test that more than 2 chains per type are filtered to top 2."""
        vdj_df = pd.DataFrame({
            "barcode": ["AAAA"] * 4,
            "chain": ["TRA", "TRA", "TRA", "TRA"],
            "cdr3": ["BEST", "SECOND", "THIRD", "FOURTH"],
            "v_gene": ["TRAV1", "TRAV2", "TRAV3", "TRAV4"],
            "d_gene": [None, None, None, None],
            "j_gene": ["TRAJ1", "TRAJ2", "TRAJ3", "TRAJ4"],
            "c_gene": ["TRAC", "TRAC", "TRAC", "TRAC"],
            "umis": [100, 80, 60, 40],
            "reads": [1000, 800, 600, 400],
            "contig_id": ["c1", "c2", "c3", "c4"],
        })

        result = _pivot_vdj_by_barcode(vdj_df)

        # Only top 2 should be kept
        assert result.loc["AAAA", "TRA_1_cdr3"] == "BEST"
        assert result.loc["AAAA", "TRA_2_cdr3"] == "SECOND"
        # No TRA_3 column
        assert "TRA_3_cdr3" not in result.columns


class TestVdjSegmentConstants:
    """Tests for VDJ segment column constants."""

    def test_segment_cols(self):
        """Test VDJ_SEGMENT_COLS list."""
        expected = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
        assert VDJ_SEGMENT_COLS == expected

    def test_segment_nt_cols(self):
        """Test VDJ_SEGMENT_NT_COLS list."""
        expected = ["fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt"]
        assert VDJ_SEGMENT_NT_COLS == expected
