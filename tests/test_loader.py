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

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from tcrsift.loader import (
    VDJ_SEGMENT_COLS,
    VDJ_SEGMENT_NT_COLS,
    _extract_tcell_markers,
    _pivot_vdj_by_barcode,
    combine_gex_and_vdj,
    load_cellranger_vdj,
    load_sample,
    load_samples,
)
from tcrsift.sample_sheet import Sample, SampleSheet
from tcrsift.validation import TCRsiftValidationError


class TestPivotVdjByBarcode:
    """Tests for _pivot_vdj_by_barcode function."""

    def test_basic_pivot(self):
        """Test basic pivoting of VDJ data."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

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
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        assert len(result) == 1
        assert "AAAA" in result.index
        # Primary chain should be highest UMI
        assert result.loc["AAAA", "TRA_1_cdr3"] == "CASSL1"
        assert result.loc["AAAA", "TRA_2_cdr3"] == "CASSL2"
        assert result.loc["AAAA", "TRB_1_cdr3"] == "CASSF1"
        assert result.loc["AAAA", "TRB_2_cdr3"] == "CASSF2"
        # Doublet flags
        assert result.loc["AAAA", "multi_TRA"]
        assert result.loc["AAAA", "multi_TRB"]
        assert result.loc["AAAA", "multi_chain"]

    def test_umi_prioritization(self):
        """Test that chains are prioritized by UMI count."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        # Higher UMI should be primary chain
        assert result.loc["AAAA", "TRA_1_cdr3"] == "HIGH_UMI"
        assert result.loc["AAAA", "TRA_2_cdr3"] == "LOW_UMI"

    def test_chain_count_columns(self):
        """Test that chain count columns are properly created."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        # AAAA has both chains
        assert result.loc["AAAA", "TRA_count"] == 1
        assert result.loc["AAAA", "TRB_count"] == 1
        assert result.loc["AAAA", "has_TRA"]
        assert result.loc["AAAA", "has_TRB"]
        assert result.loc["AAAA", "has_both_chains"]

        # BBBB only has TRB
        assert result.loc["BBBB", "TRA_count"] == 0
        assert result.loc["BBBB", "TRB_count"] == 1
        assert not result.loc["BBBB", "has_TRA"]
        assert result.loc["BBBB", "has_TRB"]
        assert not result.loc["BBBB", "has_both_chains"]

    def test_cdr3ab_identifier(self):
        """Test creation of CDR3ab identifier."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        assert result.loc["AAAA", "CDR3_alpha"] == "CASSL"
        assert result.loc["AAAA", "CDR3_beta"] == "CASSF"
        assert result.loc["AAAA", "CDR3ab"] == "CASSL_CASSF"

    def test_segment_preservation(self):
        """Test that VDJ segment columns are preserved in pivot."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        # Check that segment columns are present
        assert "TRA_1_fwr1" in result.columns
        assert "TRA_1_cdr1" in result.columns
        assert "TRB_1_fwr1" in result.columns
        assert result.loc["AAAA", "TRA_1_fwr1"] == "MRLV"
        assert result.loc["AAAA", "TRB_1_cdr1"] == "SGHD"

    def test_vdj_sequence_preservation(self):
        """Test that combined VDJ sequences are preserved."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        assert "TRA_1_vdj_aa" in result.columns
        assert "TRB_1_vdj_nt" in result.columns
        assert result.loc["AAAA", "TRA_1_vdj_aa"] == "MRLVTSGFWYRQYSSGCASSL"

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

        result = _pivot_vdj_by_barcode(vdj_df)

        assert len(result) == 0

    def test_more_than_two_chains_filtered(self):
        """Test that more than 2 chains per type are filtered to top 2."""
        vdj_df = pd.DataFrame(
            {
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
            }
        )

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


class TestRealisticVdjPivot:
    """Tests using realistic CellRanger-like VDJ data."""

    def test_pivot_with_full_segments(self, sample_vdj_df_with_segments):
        """Test pivoting VDJ data with full IMGT segment columns."""
        result = _pivot_vdj_by_barcode(sample_vdj_df_with_segments)

        assert len(result) == 1
        barcode = "AAACCTGAGAACTCGG-1"
        assert barcode in result.index

        # Check CDR3 sequences are preserved
        assert result.loc[barcode, "CDR3_alpha"] == "CAVSDGGSQGNLIF"
        assert result.loc[barcode, "CDR3_beta"] == "CASSLGQAYEQYF"

        # Check segment columns are pivoted
        assert "TRA_1_fwr1" in result.columns
        assert "TRB_1_fwr1" in result.columns
        assert result.loc[barcode, "TRA_1_fwr1"] == "MSLGLLCCVALSLLNAGTS"
        assert result.loc[barcode, "TRB_1_cdr1"] == "SGHATL"

    def test_realistic_umi_counts(self, sample_vdj_df):
        """Test that realistic UMI counts are handled correctly."""
        result = _pivot_vdj_by_barcode(sample_vdj_df)

        # Check UMI values are preserved in pivoted form
        assert "TRA_1_umis" in result.columns
        assert "TRB_1_umis" in result.columns

        # UMI values should be in realistic range (10-200 based on fixture)
        for col in result.columns:
            if col.endswith("_umis"):
                non_null_values = result[col].dropna()
                if len(non_null_values) > 0:
                    assert all(non_null_values >= 10)
                    assert all(non_null_values <= 200)

    def test_imgt_gene_naming(self, sample_vdj_df):
        """Test that IMGT gene names are preserved correctly."""
        result = _pivot_vdj_by_barcode(sample_vdj_df)

        # V genes should follow TRAV/TRBV naming
        assert "TRA_1_v_gene" in result.columns
        v_genes = result["TRA_1_v_gene"].dropna()
        assert all(vg.startswith("TRAV") for vg in v_genes)

        # J genes should follow TRAJ/TRBJ naming
        assert "TRB_1_j_gene" in result.columns
        j_genes = result["TRB_1_j_gene"].dropna()
        assert all(jg.startswith("TRBJ") for jg in j_genes)

    def test_cdr3_sequence_format(self, sample_vdj_df):
        """Test CDR3 sequences follow IMGT conventions."""
        result = _pivot_vdj_by_barcode(sample_vdj_df)

        # Alpha CDR3 typically starts with CAV or CAA
        alpha_cdr3 = result["CDR3_alpha"].dropna()
        assert all(cdr3.startswith("CA") for cdr3 in alpha_cdr3)

        # Beta CDR3 typically starts with CASS or CAS
        beta_cdr3 = result["CDR3_beta"].dropna()
        assert all(cdr3.startswith("CAS") for cdr3 in beta_cdr3)

    def test_productive_contigs_only(self, sample_vdj_df):
        """Test that all contigs in fixture are productive."""
        # Verify fixture has productive column
        assert "productive" in sample_vdj_df.columns
        assert all(sample_vdj_df["productive"])

        result = _pivot_vdj_by_barcode(sample_vdj_df)
        # All pivoted entries should come from productive contigs
        assert len(result) > 0

    def test_clonotype_grouping(self, sample_vdj_df):
        """Test that cells with same clonotype have same CDR3."""
        result = _pivot_vdj_by_barcode(sample_vdj_df)

        # AAAA and BBBB should have same clone (clonotype1)
        assert result.loc["AAAA", "CDR3ab"] == result.loc["BBBB", "CDR3ab"]

        # CCCC should have different clone
        assert result.loc["AAAA", "CDR3ab"] != result.loc["CCCC", "CDR3ab"]


class TestCombineGexAndVdjBarcodeMatching:
    """Tests for combine_gex_and_vdj barcode matching."""

    def test_vdj_suffix_mapping_with_gex_no_suffix(self):
        """VDJ barcodes with suffixes should map to GEX barcodes without suffix."""
        # GEX with no suffix
        adata = ad.AnnData(
            X=sp.csr_matrix((2, 1)),
            obs=pd.DataFrame(index=["AAAA", "BBBB"]),
            var=pd.DataFrame(index=["GENE1"]),
        )

        # VDJ with suffix
        vdj_df = pd.DataFrame(
            {
                "barcode": ["AAAA-1", "AAAA-1", "BBBB-1", "BBBB-1"],
                "chain": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3": ["CAVAAA", "CASSAAA", "CAVBBB", "CASSBBB"],
                "v_gene": ["TRAV1", "TRBV1", "TRAV2", "TRBV2"],
                "d_gene": [None, "TRBD1", None, "TRBD1"],
                "j_gene": ["TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2"],
                "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1"],
                "umis": [10, 12, 8, 9],
                "reads": [100, 120, 80, 90],
                "contig_id": ["c1", "c2", "c3", "c4"],
            }
        )

        combined = combine_gex_and_vdj(adata, vdj_df)

        assert combined.obs.loc["AAAA", "CDR3_alpha"] == "CAVAAA"
        assert combined.obs.loc["AAAA", "CDR3_beta"] == "CASSAAA"
        assert combined.obs.loc["BBBB", "CDR3_alpha"] == "CAVBBB"
        assert combined.obs.loc["BBBB", "CDR3_beta"] == "CASSBBB"


class TestLoadSamplesWithSampleSheet:
    """Tests for load_samples accepting SampleSheet objects."""

    def test_load_samples_with_samplesheet(self, mock_cellranger_vdj_dir):
        """load_samples should accept a SampleSheet instance."""
        sample = Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir))
        sheet = SampleSheet(samples=[sample])

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.n_obs > 0
        assert "sample" in adata.obs.columns
        assert adata.obs["sample"].iloc[0] == "S1"

    def test_load_samples_uns_sample_sheet_round_trips_h5ad(
        self, mock_cellranger_vdj_dir, tmp_path
    ):
        """uns['sample_sheet'] must survive write_h5ad round-trip even when
        optional fields are unset (regression for the int-keyed-dict bug)."""
        import json

        sample = Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir))
        sheet = SampleSheet(samples=[sample])

        adata = load_samples(sheet, show_progress=False, verbose=False)
        path = tmp_path / "loaded.h5ad"
        adata.write_h5ad(path)

        reloaded = ad.read_h5ad(path)
        assert "sample_sheet" in reloaded.uns
        records = json.loads(reloaded.uns["sample_sheet"])
        assert any(r.get("sample") == "S1" for r in records)

    def test_load_samples_concats_multiple_vdj_samples(self, mock_cellranger_vdj_dir):
        """Integration test: two VDJ-only samples loaded end-to-end via the
        real CellRanger loader merge into a single AnnData with the union of
        obs rows and the expected sample labels."""
        sheet = SampleSheet(
            samples=[
                Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir)),
                Sample(sample="S2", vdj_dir=str(mock_cellranger_vdj_dir)),
            ]
        )

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.n_obs > 0
        assert set(adata.obs["sample"]) == {"S1", "S2"}

    def test_load_samples_merges_gex_samples_via_concat_on_disk(
        self, mock_cellranger_vdj_dir, monkeypatch
    ):
        """Per-sample AnnData with a real GEX matrix should be spilled and
        merged via concat_on_disk — this is the load-bearing memory-bounded
        path (issue #30)."""
        from tcrsift import loader as loader_module

        n_cells, n_genes = 30, 50
        rng = np.random.default_rng(42)

        def fake_load_sample(sample, **_kwargs):
            X = sp.random(
                n_cells, n_genes, density=0.2, format="csr", random_state=rng
            )
            a = ad.AnnData(X=X)
            a.var_names = [f"G{i}" for i in range(n_genes)]
            a.obs_names = [f"{sample.sample}_cell{i}" for i in range(n_cells)]
            a.obs["sample"] = sample.sample
            return a

        monkeypatch.setattr(loader_module, "load_sample", fake_load_sample)

        sheet = SampleSheet(
            samples=[
                Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir)),
                Sample(sample="S2", vdj_dir=str(mock_cellranger_vdj_dir)),
            ]
        )

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.n_obs == 2 * n_cells
        assert adata.n_vars == n_genes
        assert adata.X is not None
        assert set(adata.obs["sample"]) == {"S1", "S2"}

    @pytest.mark.parametrize(
        "x_format",
        # Production data hits us in formats that synthetic CSR-only tests
        # never exercise: CellRanger writes CSC (#37), scipy 1.11+ ships a
        # separate `sparse_array` family that `isinstance(X, csc_matrix)`
        # would miss. concat_on_disk only accepts CSR or dense along obs;
        # _ensure_x normalizes everything else. CSC is the real production
        # shape; CSR is the synthetic happy path; dense is exercised when
        # dask is available (the dense concat code path uses dask).
        ["csr", "csc", "dense"],
    )
    def test_load_samples_handles_all_x_formats(
        self, mock_cellranger_vdj_dir, monkeypatch, x_format
    ):
        """Every input X format that anndata accepts should survive the
        spill → concat_on_disk → read-back round-trip."""
        if x_format == "dense":
            # anndata.experimental.concat_on_disk delegates dense concat to
            # dask, which is not a hard tcrsift dependency. CellRanger never
            # produces dense X so this is a synthetic-only path.
            pytest.importorskip("dask")
        from tcrsift import loader as loader_module

        n_cells, n_genes = 15, 12
        rng = np.random.default_rng(11)

        def fake_load_sample(sample, **_kwargs):
            if x_format == "dense":
                X = rng.poisson(2, (n_cells, n_genes)).astype(np.float32)
            else:
                X = sp.random(
                    n_cells, n_genes, density=0.25, format=x_format,
                    random_state=rng,
                )
            a = ad.AnnData(X=X)
            a.var_names = [f"G{i}" for i in range(n_genes)]
            a.obs_names = [f"{sample.sample}_c{i}" for i in range(n_cells)]
            a.obs["sample"] = sample.sample
            return a

        monkeypatch.setattr(loader_module, "load_sample", fake_load_sample)

        sheet = SampleSheet(
            samples=[
                Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir)),
                Sample(sample="S2", vdj_dir=str(mock_cellranger_vdj_dir)),
            ]
        )

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.n_obs == 2 * n_cells
        assert adata.n_vars == n_genes
        assert adata.X is not None

    def test_load_samples_mixed_modality_merges_uniformly(
        self, mock_cellranger_vdj_dir, monkeypatch
    ):
        """A sample sheet mixing one GEX-bearing and one VDJ-only sample should
        merge through the same concat_on_disk path. The combined AnnData carries
        the GEX vars; the VDJ-only sample's rows are present but their X entries
        are zero (outer-join fill)."""
        from tcrsift import loader as loader_module

        n_cells_per_sample = 10
        n_genes = 20

        def fake_load_sample(sample, **_kwargs):
            if sample.sample == "GEX":
                X = sp.random(
                    n_cells_per_sample, n_genes, density=0.3, format="csr",
                    random_state=0,
                )
                a = ad.AnnData(X=X)
                a.var_names = [f"G{i}" for i in range(n_genes)]
            else:
                a = ad.AnnData(shape=(n_cells_per_sample, 0))
            a.obs_names = [f"{sample.sample}_c{i}" for i in range(n_cells_per_sample)]
            a.obs["sample"] = sample.sample
            return a

        monkeypatch.setattr(loader_module, "load_sample", fake_load_sample)

        sheet = SampleSheet(
            samples=[
                Sample(sample="GEX", vdj_dir=str(mock_cellranger_vdj_dir)),
                Sample(sample="VDJ", vdj_dir=str(mock_cellranger_vdj_dir)),
            ]
        )

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.n_obs == 2 * n_cells_per_sample
        assert adata.n_vars == n_genes
        assert adata.X is not None
        assert set(adata.obs["sample"]) == {"GEX", "VDJ"}

        # Outer-join fill: the VDJ-only rows have zero X entries in the GEX
        # vars (no overlap means concat_on_disk fills with zeros, not NaN, for
        # the sparse case).
        vdj_rows = adata.obs["sample"] == "VDJ"
        vdj_block = adata[vdj_rows].X
        assert vdj_block.nnz == 0, "VDJ-only rows should have no expression entries"

    def test_load_samples_vdj_only_result_has_x_none(
        self, mock_cellranger_vdj_dir
    ):
        """When every sample is VDJ-only, combined.X should be None (the empty
        placeholders synthesized for concat_on_disk are stripped after merge)."""
        sheet = SampleSheet(
            samples=[
                Sample(sample="V1", vdj_dir=str(mock_cellranger_vdj_dir)),
                Sample(sample="V2", vdj_dir=str(mock_cellranger_vdj_dir)),
            ]
        )

        adata = load_samples(sheet, show_progress=False, verbose=False)

        assert adata.X is None
        assert adata.n_vars == 0

    def test_load_samples_respects_tmpdir_kwarg(
        self, mock_cellranger_vdj_dir, monkeypatch, tmp_path
    ):
        """The tmpdir kwarg should redirect the spill location off the system
        tempdir (relevant when $TMPDIR is on tmpfs)."""
        from tcrsift import loader as loader_module

        # Intercept concat_on_disk to capture the in_files paths it sees.
        seen_paths: list[Path] = []
        real_concat_on_disk = loader_module.concat_on_disk

        def tracking_concat_on_disk(in_files, out_file, **kwargs):
            seen_paths.extend(Path(p) for p in in_files)
            return real_concat_on_disk(in_files, out_file, **kwargs)

        monkeypatch.setattr(loader_module, "concat_on_disk", tracking_concat_on_disk)

        spill_root = tmp_path / "spill"
        spill_root.mkdir()
        sheet = SampleSheet(samples=[Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir))])

        adata = load_samples(
            sheet, show_progress=False, verbose=False, tmpdir=spill_root
        )

        assert adata.n_obs > 0
        assert len(seen_paths) == 1
        assert spill_root in seen_paths[0].parents


class TestLoadSampleMetadata:
    """Tests for sample metadata propagation."""

    def test_load_sample_adds_antigen_name(self, mock_cellranger_vdj_dir):
        """antigen_name should be propagated into adata.obs."""
        sample = Sample(
            sample="S1",
            vdj_dir=str(mock_cellranger_vdj_dir),
            antigen_name="CMV_pp65",
        )

        adata = load_sample(sample)

        assert "antigen_name" in adata.obs.columns
        assert adata.obs["antigen_name"].iloc[0] == "CMV_pp65"

    def test_load_sample_propagates_patient_id_and_enrichment_method(
        self, mock_cellranger_vdj_dir
    ):
        """patient_id and enrichment_method should land in adata.obs (issue #8)."""
        sample = Sample(
            sample="B1-2_AIMpos",
            vdj_dir=str(mock_cellranger_vdj_dir),
            patient_id="B1-2",
            enrichment_method="AIMpos",
        )

        adata = load_sample(sample)

        assert adata.obs["patient_id"].iloc[0] == "B1-2"
        assert adata.obs["enrichment_method"].iloc[0] == "AIMpos"

    def test_load_sample_omits_unset_donor_method(self, mock_cellranger_vdj_dir):
        """When donor/method are not set, columns should be absent (consistent
        with the #5 fix for other optional metadata)."""
        sample = Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir))

        adata = load_sample(sample)

        assert "patient_id" not in adata.obs.columns
        assert "enrichment_method" not in adata.obs.columns

    def test_load_sample_propagates_timepoint_apc_tissue(
        self, mock_cellranger_vdj_dir
    ):
        """timepoint, apc_type, and tissue should land in adata.obs (#9)."""
        sample = Sample(
            sample="D7_mDC",
            vdj_dir=str(mock_cellranger_vdj_dir),
            timepoint="D7",
            apc_type="mDC",
            tissue="blood",
        )

        adata = load_sample(sample)

        assert adata.obs["timepoint"].iloc[0] == "D7"
        assert adata.obs["apc_type"].iloc[0] == "mDC"
        assert adata.obs["tissue"].iloc[0] == "blood"

    def test_load_sample_omits_unset_axes(self, mock_cellranger_vdj_dir):
        """When axis fields aren't set, the columns should be absent."""
        sample = Sample(sample="S1", vdj_dir=str(mock_cellranger_vdj_dir))
        adata = load_sample(sample)
        assert "timepoint" not in adata.obs.columns
        assert "apc_type" not in adata.obs.columns
        assert "tissue" not in adata.obs.columns


class TestLoadCellrangerVdj:
    """Tests for load_cellranger_vdj function."""

    @pytest.fixture
    def mock_vdj_dir(self, tmp_path):
        """Create a mock CellRanger VDJ output directory."""
        vdj_dir = tmp_path / "vdj_outs"
        vdj_dir.mkdir()

        # Create filtered_contig_annotations.csv
        annotations = pd.DataFrame(
            {
                "barcode": [
                    "AAAA-1",
                    "AAAA-1",
                    "BBBB-1",
                    "BBBB-1",
                    "CCCC-1",
                    "CCCC-1",
                ],
                "is_cell": [True, True, True, True, True, True],
                "contig_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
                "high_confidence": [True, True, True, True, True, True],
                "length": [500, 600, 550, 620, 480, 590],
                "chain": ["TRA", "TRB", "TRA", "TRB", "TRA", "TRB"],
                "v_gene": [
                    "TRAV1-1",
                    "TRBV2-1",
                    "TRAV1-1",
                    "TRBV2-1",
                    "TRAV3-1",
                    "TRBV4-1",
                ],
                "d_gene": [None, "TRBD1", None, "TRBD1", None, "TRBD2"],
                "j_gene": [
                    "TRAJ10",
                    "TRBJ2-1",
                    "TRAJ10",
                    "TRBJ2-1",
                    "TRAJ20",
                    "TRBJ1-1",
                ],
                "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1", "TRAC", "TRBC2"],
                "full_length": [True, True, True, True, True, True],
                "productive": [True, True, True, True, True, True],
                "cdr3": [
                    "CAVSDLEPNSSASKIIF",
                    "CASSLGQAYEQYF",
                    "CAVSDLEPNSSASKIIF",
                    "CASSLGQAYEQYF",
                    "CAVRSGYSTLTF",
                    "CASSLAPGTGELFF",
                ],
                "cdr3_nt": ["TGTGCC", "TGTGCC", "TGTGCC", "TGTGCC", "TGTGCC", "TGTGCC"],
                "reads": [1000, 2000, 1500, 2500, 800, 1200],
                "umis": [100, 200, 150, 250, 80, 120],
                "raw_clonotype_id": [
                    "clonotype1",
                    "clonotype1",
                    "clonotype1",
                    "clonotype1",
                    "clonotype2",
                    "clonotype2",
                ],
            }
        )
        annotations.to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        return vdj_dir

    @pytest.fixture
    def mock_vdj_dir_with_clonotypes(self, mock_vdj_dir):
        """Create mock VDJ dir with clonotypes.csv containing MAIT/iNKT evidence."""
        clonotypes = pd.DataFrame(
            {
                "clonotype_id": ["clonotype1", "clonotype2"],
                "frequency": [2, 1],
                "proportion": [0.67, 0.33],
                "cdr3s_aa": [
                    "TRA:CAVSDLEPNSSASKIIF;TRB:CASSLGQAYEQYF",
                    "TRA:CAVRSGYSTLTF;TRB:CASSLAPGTGELFF",
                ],
                "mait_evidence": ["none", "none"],
                "inkt_evidence": ["none", "none"],
            }
        )
        clonotypes.to_csv(mock_vdj_dir / "clonotypes.csv", index=False)
        return mock_vdj_dir

    def test_load_basic(self, mock_vdj_dir):
        """Test basic VDJ loading."""
        df = load_cellranger_vdj(mock_vdj_dir, "test_sample", verbose=False)

        assert len(df) == 6
        assert "sample" in df.columns
        assert df["sample"].iloc[0] == "test_sample"
        assert "barcode" in df.columns
        assert "chain" in df.columns
        assert "cdr3" in df.columns

    def test_load_adds_sample_name(self, mock_vdj_dir):
        """Test that sample name is added to DataFrame."""
        df = load_cellranger_vdj(mock_vdj_dir, "my_sample", verbose=False)

        assert all(df["sample"] == "my_sample")
        assert all(df["vdj_dir"] == str(mock_vdj_dir))

    def test_load_with_clonotypes(self, mock_vdj_dir_with_clonotypes):
        """Test loading with clonotypes.csv containing MAIT/iNKT evidence."""
        df = load_cellranger_vdj(mock_vdj_dir_with_clonotypes, "test_sample", verbose=False)

        assert "mait_evidence" in df.columns
        assert "inkt_evidence" in df.columns

    def test_load_invalid_dir_raises(self, tmp_path):
        """Test that invalid directory raises error."""
        with pytest.raises(TCRsiftValidationError):
            load_cellranger_vdj(tmp_path / "nonexistent", "test", verbose=False)

    def test_load_empty_file_raises(self, tmp_path):
        """Test that empty annotations file raises error."""
        vdj_dir = tmp_path / "empty_vdj"
        vdj_dir.mkdir()

        # Create empty CSV with just headers
        empty_df = pd.DataFrame(columns=["barcode", "chain", "cdr3"])
        empty_df.to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        with pytest.raises(TCRsiftValidationError, match="empty"):
            load_cellranger_vdj(vdj_dir, "test", verbose=False)

    def test_load_missing_columns_raises(self, tmp_path):
        """Test that file with missing required columns raises error."""
        vdj_dir = tmp_path / "bad_vdj"
        vdj_dir.mkdir()

        # Create CSV missing required columns
        bad_df = pd.DataFrame({"some_col": [1, 2, 3], "other_col": ["a", "b", "c"]})
        bad_df.to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        with pytest.raises(TCRsiftValidationError, match="missing required columns"):
            load_cellranger_vdj(vdj_dir, "test", verbose=False)

    def test_load_uses_all_contig_fallback(self, tmp_path):
        """Test that loader falls back to all_contig_annotations.csv."""
        vdj_dir = tmp_path / "vdj_fallback"
        vdj_dir.mkdir()

        # Create only all_contig_annotations.csv (not filtered)
        annotations = pd.DataFrame(
            {
                "barcode": ["AAAA-1", "AAAA-1"],
                "chain": ["TRA", "TRB"],
                "cdr3": ["CAVSDL", "CASSL"],
                "v_gene": ["TRAV1", "TRBV2"],
                "j_gene": ["TRAJ1", "TRBJ2"],
                "umis": [100, 200],
                "reads": [1000, 2000],
            }
        )
        annotations.to_csv(vdj_dir / "all_contig_annotations.csv", index=False)

        df = load_cellranger_vdj(vdj_dir, "test_sample", verbose=False)
        assert len(df) == 2


class TestExtractTcellMarkers:
    """Tests for _extract_tcell_markers function."""

    @pytest.fixture
    def mock_adata_with_markers(self):
        """Create mock AnnData with T-cell marker genes."""
        n_cells = 50
        n_genes = 20

        # Create sparse expression matrix
        X = sp.random(n_cells, n_genes, density=0.3, format="csr")

        # Gene names including T-cell markers
        var_names = [f"Gene{i}" for i in range(n_genes - 6)]
        var_names.extend(["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"])

        adata = ad.AnnData(X)
        adata.var_names = var_names
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

        # Set CD8+ expression pattern for first half
        cd8a_idx = var_names.index("CD8A")
        cd8b_idx = var_names.index("CD8B")
        cd4_idx = var_names.index("CD4")

        # Modify expression: first 25 cells are CD8+, next 25 are CD4+
        X_dense = X.toarray()
        X_dense[:25, cd8a_idx] = np.random.poisson(20, 25)
        X_dense[:25, cd8b_idx] = np.random.poisson(15, 25)
        X_dense[:25, cd4_idx] = np.random.poisson(2, 25)
        X_dense[25:, cd4_idx] = np.random.poisson(20, 25)
        X_dense[25:, cd8a_idx] = np.random.poisson(2, 25)
        X_dense[25:, cd8b_idx] = np.random.poisson(2, 25)

        adata.X = sp.csr_matrix(X_dense)

        return adata

    def test_extract_markers_basic(self, mock_adata_with_markers):
        """Test basic T-cell marker extraction."""
        markers = _extract_tcell_markers(mock_adata_with_markers)

        assert isinstance(markers, pd.DataFrame)
        assert len(markers) == 50
        assert "CD3D" in markers.columns
        assert "CD8A" in markers.columns
        assert "CD4" in markers.columns

    def test_extract_markers_missing_genes(self):
        """Test extraction when some marker genes are missing."""
        # Create AnnData without CD8B
        n_cells = 10
        n_genes = 5
        X = sp.random(n_cells, n_genes, density=0.3, format="csr")

        adata = ad.AnnData(X)
        adata.var_names = ["Gene1", "Gene2", "CD3D", "CD4", "CD8A"]
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

        markers = _extract_tcell_markers(adata)

        assert "CD3D" in markers.columns
        assert "CD8A" in markers.columns
        # Missing genes are added with value 0
        assert "CD8B" in markers.columns
        assert (markers["CD8B"] == 0).all()

    def test_extract_markers_all_missing(self):
        """Exercise the no-markers-present branch (empty present_cols).

        Structurally distinct from the mixed case: skips the column slice +
        DataFrame construction from a 2D block entirely, falls through to
        the empty-DataFrame path with all columns added as zeros (#35)."""
        n_cells = 10
        X = sp.random(n_cells, 3, density=0.3, format="csr")
        adata = ad.AnnData(X)
        adata.var_names = ["GeneX", "GeneY", "GeneZ"]
        adata.obs_names = [f"c{i}" for i in range(n_cells)]

        markers = _extract_tcell_markers(adata)

        for sym in ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"]:
            assert sym in markers.columns
            assert (markers[sym] == 0).all()
        assert len(markers) == n_cells


class TestCombineGexAndVdj:
    """Tests for combine_gex_and_vdj function."""

    @pytest.fixture
    def mock_gex_adata(self):
        """Create mock GEX AnnData."""
        n_cells = 10
        n_genes = 15

        X = sp.random(n_cells, n_genes, density=0.3, format="csr")

        var_names = [f"Gene{i}" for i in range(n_genes - 4)]
        var_names.extend(["CD3D", "CD4", "CD8A", "CD8B"])

        adata = ad.AnnData(X)
        adata.var_names = var_names
        adata.obs_names = [f"CELL{i:04d}-1" for i in range(n_cells)]
        adata.obs["sample"] = "test_sample"

        return adata

    @pytest.fixture
    def mock_vdj_df(self):
        """Create mock VDJ DataFrame."""
        # Create VDJ data for 5 of the 10 cells (partial overlap)
        return pd.DataFrame(
            {
                "barcode": [
                    "CELL0000-1",
                    "CELL0000-1",
                    "CELL0001-1",
                    "CELL0001-1",
                    "CELL0002-1",
                    "CELL0002-1",
                    "CELL0003-1",
                    "CELL0003-1",
                    "CELL0004-1",
                    "CELL0004-1",
                ],
                "chain": [
                    "TRA",
                    "TRB",
                    "TRA",
                    "TRB",
                    "TRA",
                    "TRB",
                    "TRA",
                    "TRB",
                    "TRA",
                    "TRB",
                ],
                "cdr3": [
                    "CAV1",
                    "CASS1",
                    "CAV2",
                    "CASS2",
                    "CAV3",
                    "CASS3",
                    "CAV4",
                    "CASS4",
                    "CAV5",
                    "CASS5",
                ],
                "v_gene": ["TRAV1"] * 10,
                "j_gene": ["TRAJ1"] * 10,
                "d_gene": [None, "TRBD1"] * 5,
                "c_gene": ["TRAC", "TRBC1"] * 5,
                "umis": [100] * 10,
                "reads": [1000] * 10,
                "contig_id": [f"c{i}" for i in range(10)],
                "sample": ["test_sample"] * 10,
            }
        )

    def test_combine_basic(self, mock_gex_adata, mock_vdj_df):
        """Test basic combination of GEX and VDJ data."""
        result = combine_gex_and_vdj(mock_gex_adata, mock_vdj_df)

        assert isinstance(result, ad.AnnData)
        # Returns all GEX cells with VDJ columns merged (left join)
        assert len(result) == 10
        assert "CDR3_alpha" in result.obs.columns
        assert "CDR3_beta" in result.obs.columns
        # First 5 cells should have VDJ data, rest should be NaN
        assert result.obs["CDR3_alpha"].notna().sum() == 5

    def test_combine_preserves_gex(self, mock_gex_adata, mock_vdj_df):
        """Test that GEX expression data is preserved."""
        result = combine_gex_and_vdj(mock_gex_adata, mock_vdj_df)

        # Expression matrix should be preserved
        assert result.X is not None
        assert result.n_vars == mock_gex_adata.n_vars

    def test_combine_adds_tcell_markers(self, mock_gex_adata, mock_vdj_df):
        """Test that T-cell marker columns are added."""
        result = combine_gex_and_vdj(mock_gex_adata, mock_vdj_df)

        # Should have T-cell marker columns
        assert "CD3D" in result.obs.columns or "gex_CD3D" in result.obs.columns

    def test_combine_no_overlap_vdj_columns_nan(self, mock_gex_adata):
        """Test combining when there's no barcode overlap."""
        # VDJ data with completely different barcodes
        vdj_df = pd.DataFrame(
            {
                "barcode": ["XXXX-1", "XXXX-1"],
                "chain": ["TRA", "TRB"],
                "cdr3": ["CAV", "CASS"],
                "v_gene": ["TRAV1", "TRBV1"],
                "j_gene": ["TRAJ1", "TRBJ1"],
                "d_gene": [None, "TRBD1"],
                "c_gene": ["TRAC", "TRBC1"],
                "umis": [100, 100],
                "reads": [1000, 1000],
                "contig_id": ["c1", "c2"],
                "sample": ["test_sample", "test_sample"],
            }
        )

        result = combine_gex_and_vdj(mock_gex_adata, vdj_df)

        # Returns all GEX cells, VDJ columns are NaN for non-matching cells
        assert len(result) == 10
        assert "CDR3_alpha" in result.obs.columns
        # No cells should have VDJ data since barcodes don't overlap
        assert result.obs["CDR3_alpha"].isna().all()


class TestLoadCellrangerGex:
    """Tests for load_cellranger_gex function."""

    @pytest.fixture
    def mock_gex_dir(self, tmp_path):
        """Create a mock CellRanger GEX output directory with 10x matrix format."""
        import gzip

        from scipy.io import mmwrite

        gex_dir = tmp_path / "gex_outs"
        gex_dir.mkdir()
        matrix_dir = gex_dir / "filtered_feature_bc_matrix"
        matrix_dir.mkdir()

        n_cells = 20
        n_genes = 30

        # Create sparse expression matrix
        np.random.seed(42)
        X = sp.random(n_cells, n_genes, density=0.3, format="coo")
        X.data = np.random.poisson(10, len(X.data)).astype(float)

        # Write matrix.mtx then gzip it (scanpy expects .mtx.gz)
        mtx_path = matrix_dir / "matrix.mtx"
        mmwrite(mtx_path, X.T)  # 10x format is genes x cells
        with open(mtx_path, "rb") as f_in:
            with gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
                f_out.write(f_in.read())
        mtx_path.unlink()  # Remove uncompressed version

        # Write features.tsv.gz (gene_id, gene_name, feature_type)
        gene_names = [f"Gene{i}" for i in range(n_genes - 6)]
        gene_names.extend(["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"])
        # Add some MT genes for mito QC
        gene_names[0] = "MT-CO1"
        gene_names[1] = "MT-CO2"

        with gzip.open(matrix_dir / "features.tsv.gz", "wt") as f:
            for i, name in enumerate(gene_names):
                f.write(f"ENSG{i:011d}\t{name}\tGene Expression\n")

        # Write barcodes.tsv.gz
        with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt") as f:
            for i in range(n_cells):
                f.write(f"CELL{i:04d}-1\n")

        return gex_dir

    # Permissive QC bounds for the 20-cell × 30-gene mock fixture, which is
    # smaller than any real CellRanger sample and can't survive production
    # defaults (min_genes=250, min_counts=500). #39 made QC actually filter.
    _PERMISSIVE_QC = dict(
        min_genes=1,
        max_genes=100,
        min_counts=1,
        max_counts=10000,
        min_mito_pct=0,
        max_mito_pct=100,
    )

    def test_load_gex_basic(self, mock_gex_dir):
        """Test basic GEX loading."""
        from tcrsift.loader import load_cellranger_gex

        adata = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )

        assert isinstance(adata, ad.AnnData)
        assert adata.n_obs == 20
        assert adata.n_vars == 30
        assert "sample" in adata.obs.columns
        assert adata.obs["sample"].iloc[0] == "test_sample"

    def test_load_gex_adds_qc_metrics(self, mock_gex_dir):
        """Test that QC metrics are added."""
        from tcrsift.loader import load_cellranger_gex

        adata = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )

        # QC columns should be present
        assert "n_genes" in adata.obs.columns
        assert "n_counts" in adata.obs.columns
        assert "percent_mt" in adata.obs.columns

    def test_load_gex_adds_filter_flags(self, mock_gex_dir):
        """Test that filter flags are added."""
        from tcrsift.loader import load_cellranger_gex

        adata = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )

        # Filter flag columns should be present
        assert "filter:min_genes" in adata.obs.columns
        assert "filter:max_genes" in adata.obs.columns
        assert "filter:pass_qc" in adata.obs.columns

    def test_load_gex_applies_qc_filter(self, mock_gex_dir):
        """QC parameters should actually drop cells that fail the bounds,
        not just be advisory flags (#39 — regression cover)."""
        from tcrsift.loader import load_cellranger_gex

        # Read with permissive bounds to learn the unfiltered cell pool.
        unfiltered = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )
        assert unfiltered.n_obs == 20

        # Now tighten min_genes to a value that cuts most cells. The mock
        # fixture's cells have n_genes ≈ 9 (30 genes × 0.3 density), so
        # min_genes=12 will drop the lower half.
        filtered = load_cellranger_gex(
            mock_gex_dir,
            "test_sample",
            verbose=False,
            min_genes=12,
            max_genes=100,
            min_counts=1,
            max_counts=10000,
            min_mito_pct=0,
            max_mito_pct=100,
        )

        assert filtered.n_obs < unfiltered.n_obs, (
            "QC filter should drop cells with n_genes < 12; "
            f"got n_obs={filtered.n_obs} == unfiltered {unfiltered.n_obs}"
        )
        # Defense: empty result would let later assertions vacuously pass.
        assert filtered.n_obs > 0
        # Every surviving cell satisfies the bound.
        assert (filtered.obs["n_genes"] >= 12).all()
        # The pass_qc flag is all True on the survivors.
        assert filtered.obs["filter:pass_qc"].all()

    # Self-calibrating thresholds derived from the mock fixture's actual
    # per-cell QC distribution: each bound is set to that bound's median
    # over the unfiltered cells. This guarantees a non-trivial cut
    # (cells both above and below) without coupling the test to fixture
    # randomness — a fixture re-seeding still leaves a median to bisect.
    @pytest.mark.parametrize(
        "param,direction,obs_col",
        [
            ("min_genes", "drop_below", "n_genes"),
            ("max_genes", "drop_above", "n_genes"),
            ("min_counts", "drop_below", "n_counts"),
            ("max_counts", "drop_above", "n_counts"),
            ("min_mito_pct", "drop_below", "percent_mt"),
            ("max_mito_pct", "drop_above", "percent_mt"),
        ],
    )
    def test_each_qc_bound_actually_filters(
        self, mock_gex_dir, param, direction, obs_col
    ):
        """Contract test: every documented QC bound must produce a smaller
        cell pool when tightened past the fixture's median. Catches the
        #39 class of bug (parameter exists in signature, gets recorded as
        an obs flag, but never filters) for any future bound."""
        from tcrsift.loader import load_cellranger_gex

        unfiltered = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )
        median = float(unfiltered.obs[obs_col].median())

        # Set this one bound past the median; keep everything else
        # permissive. For "drop_below" bounds (min_*) we set the bound to
        # the median, which drops cells below it. For "drop_above" bounds
        # (max_*) we set the bound to the median, which drops cells above.
        kwargs = dict(self._PERMISSIVE_QC)
        kwargs[param] = median
        filtered = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **kwargs
        )

        assert filtered.n_obs < unfiltered.n_obs, (
            f"Bound {param}={median} should drop cells "
            f"({direction} median {obs_col}={median}); "
            f"got n_obs={filtered.n_obs} == unfiltered {unfiltered.n_obs}"
        )
        # Surviving cells satisfy the bound (inclusive on both sides per
        # load_cellranger_gex's `>=` / `<=` semantics).
        if direction == "drop_below":
            assert (filtered.obs[obs_col] >= median).all()
        else:
            assert (filtered.obs[obs_col] <= median).all()

    def test_load_gex_warns_on_low_pass_rate(self, mock_gex_dir, caplog):
        """When QC drops more than half of cells, emit a warning that names
        the relevant CLI flags (#41). The verbose-only info log isn't loud
        enough to surface a config typo like --min-mito 1.0 (meant 10)."""
        import logging

        from tcrsift.loader import load_cellranger_gex

        # min_genes well above the fixture's typical n_genes (~9) drops
        # nearly every cell.
        with caplog.at_level(logging.WARNING, logger="tcrsift.loader"):
            load_cellranger_gex(
                mock_gex_dir,
                "test_sample",
                min_genes=20,  # drops most of the 20-cell × 30-gene mock
                max_genes=100,
                min_counts=1,
                max_counts=10000,
                min_mito_pct=0,
                max_mito_pct=100,
                verbose=False,
            )

        warnings = [
            r.message for r in caplog.records
            if r.levelno == logging.WARNING and "QC dropped" in r.message
        ]
        assert warnings, "expected a QC-drop warning"
        # The warning should name the relevant flags so the user knows
        # what to inspect.
        assert any("--min-genes" in m for m in warnings)
        assert any("--min-mito" in m for m in warnings)

    def test_load_gex_no_warning_at_normal_pass_rate(self, mock_gex_dir, caplog):
        """Permissive bounds keep most cells — no warning should fire."""
        import logging

        from tcrsift.loader import load_cellranger_gex

        with caplog.at_level(logging.WARNING, logger="tcrsift.loader"):
            load_cellranger_gex(
                mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
            )

        qc_warnings = [
            r.message for r in caplog.records
            if r.levelno == logging.WARNING and "QC dropped" in r.message
        ]
        assert not qc_warnings, f"unexpected warning: {qc_warnings}"

    def test_min_mito_floor_warns_on_large_cull(self, mock_gex_dir, caplog):
        """A min_mito_pct FLOOR that culls an unusually large fraction (>10%)
        must warn loudly so it can't silently shrink the cohort (#168)."""
        import logging

        from tcrsift.loader import load_cellranger_gex

        # Set the floor at the 75th percentile so it drops ~75% of cells —
        # comfortably past the 10% warning threshold.
        unfiltered = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )
        high_floor = float(unfiltered.obs["percent_mt"].quantile(0.75))

        kwargs = dict(self._PERMISSIVE_QC)
        kwargs["min_mito_pct"] = high_floor if high_floor > 0 else 1.0
        with caplog.at_level(logging.WARNING, logger="tcrsift.loader"):
            load_cellranger_gex(
                mock_gex_dir, "test_sample", verbose=False, **kwargs
            )

        floor_warnings = [
            r.message for r in caplog.records
            if r.levelno == logging.WARNING and "FLOOR" in r.message
        ]
        assert floor_warnings, "expected a min_mito floor warning on a large cull"

    def test_min_mito_default_floor_no_warning(self, mock_gex_dir, caplog):
        """A small default floor trims a tail (info-level), not a warning — the
        warning is reserved for an unusually large (>10%) cull (#168)."""
        import logging

        from tcrsift.loader import load_cellranger_gex

        # Floor just above the minimum: drops at most the bottom sliver of cells,
        # which should stay under the 10% warn threshold.
        unfiltered = load_cellranger_gex(
            mock_gex_dir, "test_sample", verbose=False, **self._PERMISSIVE_QC
        )
        low_floor = float(unfiltered.obs["percent_mt"].quantile(0.05))

        kwargs = dict(self._PERMISSIVE_QC)
        kwargs["min_mito_pct"] = low_floor
        with caplog.at_level(logging.WARNING, logger="tcrsift.loader"):
            load_cellranger_gex(
                mock_gex_dir, "test_sample", verbose=False, **kwargs
            )

        floor_warnings = [
            r.message for r in caplog.records
            if r.levelno == logging.WARNING and "FLOOR" in r.message
        ]
        assert not floor_warnings, f"unexpected floor warning: {floor_warnings}"

    def test_min_mito_zero_no_floor_warning(self, mock_gex_dir, caplog):
        """Setting min_mito_pct=0 disables the floor, so no floor warning."""
        import logging

        from tcrsift.loader import load_cellranger_gex

        with caplog.at_level(logging.WARNING, logger="tcrsift.loader"):
            load_cellranger_gex(
                mock_gex_dir,
                "test_sample",
                verbose=False,
                min_genes=1,
                max_genes=100,
                min_counts=1,
                max_counts=10000,
                min_mito_pct=0,
                max_mito_pct=100,
            )

        floor_warnings = [
            r.message for r in caplog.records
            if r.levelno == logging.WARNING and "FLOOR" in r.message
        ]
        assert not floor_warnings, f"unexpected floor warning: {floor_warnings}"

    def test_load_gex_invalid_dir_raises(self, tmp_path):
        """Test that invalid directory raises error."""
        from tcrsift.loader import load_cellranger_gex

        with pytest.raises(TCRsiftValidationError):
            load_cellranger_gex(tmp_path / "nonexistent", "test", verbose=False)

    def test_load_gex_invalid_params_raises(self, mock_gex_dir):
        """Test that invalid QC parameters raise errors."""
        from tcrsift.loader import load_cellranger_gex

        # min_genes > max_genes
        with pytest.raises(TCRsiftValidationError, match="min_genes"):
            load_cellranger_gex(mock_gex_dir, "test", min_genes=1000, max_genes=100, verbose=False)

        # min_counts > max_counts
        with pytest.raises(TCRsiftValidationError, match="min_counts"):
            load_cellranger_gex(
                mock_gex_dir, "test", min_counts=1000, max_counts=100, verbose=False
            )

        # min_mito > max_mito
        with pytest.raises(TCRsiftValidationError, match="min_mito"):
            load_cellranger_gex(
                mock_gex_dir, "test", min_mito_pct=50, max_mito_pct=10, verbose=False
            )


class TestLoadSample:
    """Tests for load_sample function."""

    @pytest.fixture
    def mock_sample_dirs(self, tmp_path):
        """Create mock VDJ and GEX directories for a single sample."""
        import gzip

        from scipy.io import mmwrite

        # Create VDJ dir
        vdj_dir = tmp_path / "vdj_outs"
        vdj_dir.mkdir()

        annotations = pd.DataFrame(
            {
                "barcode": ["CELL0000-1", "CELL0000-1", "CELL0001-1", "CELL0001-1"],
                "chain": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3": ["CAVSDLEPNSSASKIIF", "CASSLGQAYEQYF", "CAVRSGYSTLTF", "CASSLAPGTGELFF"],
                "v_gene": ["TRAV1-1", "TRBV2-1", "TRAV3-1", "TRBV4-1"],
                "d_gene": [None, "TRBD1", None, "TRBD2"],
                "j_gene": ["TRAJ10", "TRBJ2-1", "TRAJ20", "TRBJ1-1"],
                "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC2"],
                "umis": [100, 200, 150, 250],
                "reads": [1000, 2000, 1500, 2500],
                "contig_id": ["c1", "c2", "c3", "c4"],
                "productive": [True, True, True, True],
            }
        )
        annotations.to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        # Create GEX dir
        gex_dir = tmp_path / "gex_outs"
        gex_dir.mkdir()
        matrix_dir = gex_dir / "filtered_feature_bc_matrix"
        matrix_dir.mkdir()

        n_cells = 5
        n_genes = 15

        np.random.seed(42)
        X = sp.random(n_cells, n_genes, density=0.4, format="coo")
        X.data = np.random.poisson(15, len(X.data)).astype(float)

        # Write and gzip matrix
        mtx_path = matrix_dir / "matrix.mtx"
        mmwrite(mtx_path, X.T)
        with open(mtx_path, "rb") as f_in:
            with gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
                f_out.write(f_in.read())
        mtx_path.unlink()

        gene_names = [f"Gene{i}" for i in range(n_genes - 6)]
        gene_names.extend(["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"])

        with gzip.open(matrix_dir / "features.tsv.gz", "wt") as f:
            for i, name in enumerate(gene_names):
                f.write(f"ENSG{i:011d}\t{name}\tGene Expression\n")

        # Barcodes match VDJ data
        with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt") as f:
            for i in range(n_cells):
                f.write(f"CELL{i:04d}-1\n")

        return {"vdj_dir": vdj_dir, "gex_dir": gex_dir}

    def test_load_sample_combines_gex_vdj(self, mock_sample_dirs):
        """Test loading a sample with both GEX and VDJ."""
        from tcrsift.loader import load_sample
        from tcrsift.sample_sheet import Sample

        sample = Sample(
            sample="test_sample",
            vdj_dir=str(mock_sample_dirs["vdj_dir"]),
            gex_dir=str(mock_sample_dirs["gex_dir"]),
        )

        adata = load_sample(sample, min_genes=1, min_counts=1, min_mito_pct=0)

        assert isinstance(adata, ad.AnnData)
        assert "CDR3_alpha" in adata.obs.columns
        assert "CDR3_beta" in adata.obs.columns
        assert "CD3D" in adata.obs.columns

    def test_load_sample_vdj_only(self, mock_sample_dirs):
        """Test loading a sample with only VDJ data."""
        from tcrsift.loader import load_sample
        from tcrsift.sample_sheet import Sample

        sample = Sample(
            sample="test_sample",
            vdj_dir=str(mock_sample_dirs["vdj_dir"]),
            gex_dir=None,
        )

        adata = load_sample(sample)

        assert isinstance(adata, ad.AnnData)
        assert "CDR3_alpha" in adata.obs.columns
        assert len(adata) == 2  # 2 unique barcodes in VDJ data

    def test_load_sample_adds_metadata(self, mock_sample_dirs):
        """Test that sample metadata is added."""
        from tcrsift.loader import load_sample
        from tcrsift.sample_sheet import Sample

        sample = Sample(
            sample="test_sample",
            vdj_dir=str(mock_sample_dirs["vdj_dir"]),
            gex_dir=str(mock_sample_dirs["gex_dir"]),
            antigen_type="short_peptide",
            antigen_description="Test antigen",
            source="culture",
        )

        adata = load_sample(sample, min_genes=1, min_counts=1, min_mito_pct=0)

        assert adata.obs["antigen_type"].iloc[0] == "short_peptide"
        assert adata.obs["antigen_description"].iloc[0] == "Test antigen"
        assert adata.obs["source"].iloc[0] == "culture"

    def test_load_sample_omits_unset_metadata_and_writes_h5ad(
        self, mock_sample_dirs, tmp_path
    ):
        """Unset metadata fields should be omitted entirely so the resulting
        AnnData can be written to h5ad (regression for issue #5)."""
        from tcrsift.loader import load_sample
        from tcrsift.sample_sheet import Sample

        sample = Sample(
            sample="test_sample",
            vdj_dir=str(mock_sample_dirs["vdj_dir"]),
            gex_dir=str(mock_sample_dirs["gex_dir"]),
            source="culture",
        )

        adata = load_sample(sample, min_genes=1, min_counts=1, min_mito_pct=0)

        for absent in (
            "antigen_type",
            "antigen_description",
            "antigen_name",
            "antigen_sequence",
            "epitope_sequence",
            "mhc_allele",
            "antigen_names",
            "antigen_sequences",
            "epitope_sequences",
        ):
            assert absent not in adata.obs.columns
        assert "source" in adata.obs.columns

        # The all-None columns from the old code path crashed h5py here.
        adata.write_h5ad(tmp_path / "out.h5ad")


class TestCombineGexAndVdjBugs:
    """Tests for bugs in combine_gex_and_vdj function."""

    @pytest.fixture
    def gex_adata_with_suffix(self):
        """Create GEX AnnData with -1 barcode suffix."""
        n_cells = 5
        n_genes = 10

        X = sp.random(n_cells, n_genes, density=0.3, format="csr")

        var_names = [f"Gene{i}" for i in range(n_genes - 4)]
        var_names.extend(["CD3D", "CD4", "CD8A", "CD8B"])

        adata = ad.AnnData(X)
        adata.var_names = var_names
        # Barcodes with -1 suffix (CellRanger format)
        adata.obs_names = [f"CELL{i:04d}-1" for i in range(n_cells)]
        adata.obs["sample"] = "test_sample"

        return adata

    @pytest.fixture
    def vdj_df_without_suffix(self):
        """Create VDJ DataFrame without barcode suffix."""
        return pd.DataFrame(
            {
                "barcode": ["CELL0000", "CELL0000", "CELL0001", "CELL0001"],
                "chain": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3": ["CAV1", "CASS1", "CAV2", "CASS2"],
                "v_gene": ["TRAV1", "TRBV1", "TRAV2", "TRBV2"],
                "j_gene": ["TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2"],
                "d_gene": [None, "TRBD1", None, "TRBD1"],
                "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1"],
                "umis": [100, 200, 150, 250],
                "reads": [1000, 2000, 1500, 2500],
                "contig_id": ["c1", "c2", "c3", "c4"],
                "sample": ["test_sample"] * 4,
            }
        )

    def test_combine_does_not_modify_original_barcodes(
        self, gex_adata_with_suffix, vdj_df_without_suffix
    ):
        """Test that combine_gex_and_vdj does not destructively modify barcodes.

        BUG: Currently the function modifies adata.obs_names in place when
        barcodes don't match, which is a side effect that callers don't expect.
        """
        original_barcodes = list(gex_adata_with_suffix.obs_names)

        combine_gex_and_vdj(gex_adata_with_suffix, vdj_df_without_suffix)

        # After combining, original AnnData's barcodes should be unchanged
        # (or we should work on a copy)
        assert list(gex_adata_with_suffix.obs_names) == original_barcodes

    def test_combine_handles_different_gem_groups(self):
        """Test that combine handles -2, -3 gem group suffixes, not just -1.

        BUG: Currently only handles -1 suffix, but CellRanger can produce
        -2, -3, etc. for multiplexed samples.
        """
        n_cells = 3
        n_genes = 8

        X = sp.random(n_cells, n_genes, density=0.3, format="csr")
        var_names = ["Gene1", "Gene2", "Gene3", "Gene4", "CD3D", "CD4", "CD8A", "CD8B"]

        adata = ad.AnnData(X)
        adata.var_names = var_names
        # Mix of gem groups
        adata.obs_names = ["CELL0000-1", "CELL0001-2", "CELL0002-3"]
        adata.obs["sample"] = "test_sample"

        vdj_df = pd.DataFrame(
            {
                "barcode": ["CELL0000", "CELL0000", "CELL0001", "CELL0001"],
                "chain": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3": ["CAV1", "CASS1", "CAV2", "CASS2"],
                "v_gene": ["TRAV1", "TRBV1", "TRAV2", "TRBV2"],
                "j_gene": ["TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2"],
                "d_gene": [None, "TRBD1", None, "TRBD1"],
                "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1"],
                "umis": [100, 200, 150, 250],
                "reads": [1000, 2000, 1500, 2500],
                "contig_id": ["c1", "c2", "c3", "c4"],
                "sample": ["test_sample"] * 4,
            }
        )

        result = combine_gex_and_vdj(adata, vdj_df)

        # Should match CELL0000-1 to CELL0000 and CELL0001-2 to CELL0001
        assert result.obs.loc["CELL0000-1", "CDR3_alpha"] == "CAV1"
        assert result.obs.loc["CELL0001-2", "CDR3_alpha"] == "CAV2"


class TestMitochondrialGeneDetection:
    """Test mitochondrial gene detection in load_cellranger_gex."""

    @pytest.fixture
    def mock_gex_dir_with_mt_genes(self, tmp_path):
        """Create mock GEX directory with MT genes using ENSEMBL IDs as var_names."""
        import gzip

        from scipy.io import mmwrite

        gex_dir = tmp_path / "gex_outs"
        gex_dir.mkdir()
        matrix_dir = gex_dir / "filtered_feature_bc_matrix"
        matrix_dir.mkdir()

        n_cells = 10
        n_genes = 6

        # Create expression matrix with higher MT expression
        np.random.seed(42)
        X = sp.random(n_cells, n_genes, density=0.5, format="coo")
        X = X.toarray()
        # Make MT genes (first 2) have higher expression
        X[:, 0] = np.random.poisson(50, n_cells)  # MT-ND1
        X[:, 1] = np.random.poisson(40, n_cells)  # MT-CO1
        X = sp.coo_matrix(X)

        # Write matrix.mtx.gz
        mtx_path = matrix_dir / "matrix.mtx"
        mmwrite(mtx_path, X.T)
        with open(mtx_path, "rb") as f_in:
            with gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
                f_out.write(f_in.read())
        mtx_path.unlink()

        # Write features.tsv.gz with ENSEMBL IDs and gene symbols
        # Include MT genes to test detection
        with gzip.open(matrix_dir / "features.tsv.gz", "wt") as f:
            f.write("ENSG00000198888\tMT-ND1\tGene Expression\n")
            f.write("ENSG00000198804\tMT-CO1\tGene Expression\n")
            f.write("ENSG00000167286\tCD3D\tGene Expression\n")
            f.write("ENSG00000010610\tCD4\tGene Expression\n")
            f.write("ENSG00000153563\tCD8A\tGene Expression\n")
            f.write("ENSG00000172116\tCD8B\tGene Expression\n")

        # Write barcodes.tsv.gz
        with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt") as f:
            for i in range(n_cells):
                f.write(f"CELL{i:04d}-1\n")

        return gex_dir

    def test_mt_detection_with_ensembl_ids(self, mock_gex_dir_with_mt_genes):
        """Test that MT genes are detected when var_names are ENSEMBL IDs."""
        from tcrsift.loader import load_cellranger_gex

        adata = load_cellranger_gex(
            mock_gex_dir_with_mt_genes,
            "test_sample",
            min_genes=1,
            min_counts=1,
            min_mito_pct=0,
            max_mito_pct=100,
            verbose=False,
        )

        # Check that MT genes were detected
        assert "mt" in adata.var.columns
        mt_genes = adata.var[adata.var["mt"]].index.tolist()

        # Should detect 2 MT genes (MT-ND1 and MT-CO1 via gene_symbols column)
        assert len(mt_genes) == 2

        # percent_mt should be calculated and non-zero
        assert "percent_mt" in adata.obs.columns
        assert adata.obs["percent_mt"].mean() > 0


class TestExtractTcellMarkersWithVersionedIds:
    """Test T-cell marker extraction with versioned ENSEMBL IDs."""

    def test_extract_markers_versioned_ensembl(self):
        """Test marker extraction when genes use versioned ENSEMBL IDs."""
        n_cells = 10
        n_genes = 8

        X = sp.random(n_cells, n_genes, density=0.3, format="csr").toarray()

        # Use versioned ENSEMBL IDs (like real CellRanger output)
        var_names = [
            "ENSG00000167286.10",  # CD3D
            "ENSG00000198851.5",  # CD3E
            "ENSG00000160654.11",  # CD3G
            "ENSG00000010610.8",  # CD4
            "ENSG00000153563.16",  # CD8A
            "ENSG00000172116.23",  # CD8B
            "OtherGene1",
            "OtherGene2",
        ]

        adata = ad.AnnData(sp.csr_matrix(X))
        adata.var_names = var_names
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

        markers = _extract_tcell_markers(adata)

        # Should find all markers even with versioned IDs
        assert "CD3D" in markers.columns
        assert "CD3E" in markers.columns
        assert "CD3G" in markers.columns
        assert "CD4" in markers.columns
        assert "CD8A" in markers.columns
        assert "CD8B" in markers.columns

        # Values should be extracted (not all zeros)
        # Since we used random sparse matrix, check structure is correct
        assert len(markers) == n_cells


class TestLoaderBugs:
    """Tests for specific bugs in loader.py."""

    def test_combine_gex_and_vdj_no_sample_name_parameter(self):
        """Test that sample_name was removed from combine_gex_and_vdj.

        The sample_name parameter was unused, so it was removed.
        """
        import inspect

        from tcrsift.loader import combine_gex_and_vdj

        # Get the function signature
        sig = inspect.signature(combine_gex_and_vdj)
        param_names = list(sig.parameters.keys())

        # sample_name should NOT be a parameter anymore
        assert "sample_name" not in param_names
        # Only adata and vdj_df should be parameters
        assert param_names == ["adata", "vdj_df"]
