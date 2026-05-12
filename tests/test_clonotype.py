"""
Tests for clonotype aggregation.
"""

import pandas as pd
import pytest

from tcrsift.clonotype import (
    aggregate_clonotypes,
    build_clone_sample_long,
    calculate_clone_frequencies,
    export_clonotypes_airr,
    get_clonotype_summary,
)


class TestAggregateClonotypes:
    """Tests for aggregate_clonotypes function."""

    @pytest.fixture
    def adata_with_tcr(self, sample_adata):
        """Create AnnData with TCR info for clonotype testing."""
        adata = sample_adata.copy()
        # Add required columns for clonotyping
        adata.obs["sample"] = ["S1"] * 50 + ["S2"] * 50

        # Create clone patterns:
        # Clone A: cells 0-19 (20 cells, same CDR3ab)
        # Clone B: cells 20-34 (15 cells, different CDR3ab)
        # Clone C: cells 35-44 (10 cells, different CDR3ab)
        # Rest: unique CDR3s

        cdr3_alpha = []
        cdr3_beta = []

        for i in range(100):
            if i < 20:
                cdr3_alpha.append("CAVSDGGSQGNLIF")
                cdr3_beta.append("CASSLGQAYEQYF")
            elif i < 35:
                cdr3_alpha.append("CAVSAGGSQGNLIF")
                cdr3_beta.append("CASSLAGAYEQYF")
            elif i < 45:
                cdr3_alpha.append("CAVNAGGSQGNLIF")
                cdr3_beta.append("CASSNAGAYEQYF")
            else:
                cdr3_alpha.append(f"CAVUNIQUE{i}QGNLIF")
                cdr3_beta.append(f"CASSUNIQUE{i}YEQYF")

        adata.obs["CDR3_alpha"] = cdr3_alpha
        adata.obs["CDR3_beta"] = cdr3_beta

        # Add Tcell_type for testing consensus
        adata.obs["Tcell_type"] = ["Confident CD8+"] * 50 + ["Confident CD4+"] * 50
        adata.obs["is_CD8"] = [True] * 50 + [False] * 50
        adata.obs["is_CD4"] = [False] * 50 + [True] * 50

        return adata

    def test_aggregate_by_cdr3ab(self, adata_with_tcr):
        """Test aggregation by both alpha and beta chains."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert len(clonotypes) > 0
        assert "CDR3ab" in clonotypes.columns
        assert "CDR3_alpha" in clonotypes.columns
        assert "CDR3_beta" in clonotypes.columns
        assert "cell_count" in clonotypes.columns

    def test_aggregate_by_cdr3b_only(self, adata_with_tcr):
        """Test aggregation by beta chain only."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3b_only")

        assert len(clonotypes) > 0
        assert "CDR3ab" in clonotypes.columns
        assert "CDR3_beta" in clonotypes.columns

    def test_aggregate_invalid_group_by(self, adata_with_tcr):
        """Invalid group_by should raise."""
        with pytest.raises(ValueError, match="Invalid group_by"):
            aggregate_clonotypes(adata_with_tcr, group_by="invalid")

    def test_aggregate_cell_counts(self, adata_with_tcr):
        """Cell counts should be correct."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        # Find Clone A (the largest)
        clone_a = clonotypes[clonotypes["CDR3_alpha"] == "CAVSDGGSQGNLIF"]
        assert len(clone_a) == 1
        assert clone_a["cell_count"].iloc[0] == 20

    def test_aggregate_sample_info(self, adata_with_tcr):
        """Sample information should be aggregated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "samples" in clonotypes.columns
        assert "n_samples" in clonotypes.columns

    def test_aggregate_tcell_consensus(self, adata_with_tcr):
        """T cell type consensus should be calculated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "Tcell_type_consensus" in clonotypes.columns

    def test_aggregate_frequency(self, adata_with_tcr):
        """Frequency should be calculated."""
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")

        assert "max_frequency" in clonotypes.columns
        # All frequencies should be between 0 and 1
        assert all(clonotypes["max_frequency"] >= 0)
        assert all(clonotypes["max_frequency"] <= 1)

    def test_handle_doublets_flag(self, adata_with_tcr):
        """Doublet handling with flag mode."""
        adata_with_tcr.obs["multi_chain"] = [i % 10 == 0 for i in range(100)]

        clonotypes = aggregate_clonotypes(adata_with_tcr, handle_doublets="flag")

        assert "n_doublet_cells" in clonotypes.columns

    def test_handle_doublets_remove(self, adata_with_tcr):
        """Doublet handling with remove mode."""
        adata_with_tcr.obs["multi_chain"] = [i % 10 == 0 for i in range(100)]

        clonotypes = aggregate_clonotypes(adata_with_tcr, handle_doublets="remove")

        # Total cells should be less due to doublet removal
        total_cells = clonotypes["cell_count"].sum()
        assert total_cells < 100

    def test_empty_result_no_complete_clones(self, sample_adata):
        """Raises error when no complete clones."""
        from tcrsift.validation import TCRsiftValidationError

        sample_adata.obs["sample"] = "S1"
        sample_adata.obs["CDR3_alpha"] = None
        sample_adata.obs["CDR3_beta"] = None

        with pytest.raises(TCRsiftValidationError, match="No complete clones found"):
            aggregate_clonotypes(sample_adata, group_by="CDR3ab")

    def test_aggregate_handles_categorical_cdr3_columns(self, adata_with_tcr):
        """After an h5ad round-trip CDR3 columns return as Categorical;
        aggregate_clonotypes must still work (regression for issue #11)."""
        # Simulate the post-roundtrip state by converting to Categorical
        # (no "" in categories — that's exactly the failing condition).
        adata_with_tcr.obs["CDR3_alpha"] = pd.Categorical(
            adata_with_tcr.obs["CDR3_alpha"].astype(object)
        )
        adata_with_tcr.obs["CDR3_beta"] = pd.Categorical(
            adata_with_tcr.obs["CDR3_beta"].astype(object)
        )

        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")
        assert len(clonotypes) > 0
        assert "CDR3ab" in clonotypes.columns

    def test_aggregate_handles_categorical_cdr3_with_nan(self, adata_with_tcr):
        """Same as above but with some NaN entries — `.fillna('')` must not
        raise on the resulting Categorical (issue #11)."""
        # Drop alpha for the last 10 cells, then categoricalize.
        alpha = adata_with_tcr.obs["CDR3_alpha"].astype(object).copy()
        alpha.iloc[-10:] = None
        adata_with_tcr.obs["CDR3_alpha"] = pd.Categorical(alpha)
        adata_with_tcr.obs["CDR3_beta"] = pd.Categorical(
            adata_with_tcr.obs["CDR3_beta"].astype(object)
        )

        # Without the fix this raised:
        #   TypeError: Cannot setitem on a Categorical with a new category ()
        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3ab")
        assert len(clonotypes) > 0

    def test_aggregate_b_only_handles_categorical_cdr3(self, adata_with_tcr):
        """CDR3b_only mode — same Categorical-fillna concern at clonotype.py:148."""
        adata_with_tcr.obs["CDR3_beta"] = pd.Categorical(
            adata_with_tcr.obs["CDR3_beta"].astype(object)
        )

        clonotypes = aggregate_clonotypes(adata_with_tcr, group_by="CDR3b_only")
        assert len(clonotypes) > 0


class TestAggregateDonorMethodColumns:
    """#8 chunk 2 — per-clone donor/method aggregations."""

    def _make_donor_method_adata(self):
        """A clone shared across two donors and several methods within each:
        - Clone A: B1-2 (AIMpos, IFNpos, tetpos), B1-3 (AIMpos)
        - Clone B: B1-2 (tetpos)
        """
        import anndata as ad
        import numpy as np

        rows = (
            [("B1-2", "AIMpos", "CAVA", "CASS_A")] * 3
            + [("B1-2", "IFNpos", "CAVA", "CASS_A")] * 2
            + [("B1-2", "tetpos", "CAVA", "CASS_A")] * 2
            + [("B1-3", "AIMpos", "CAVA", "CASS_A")] * 2
            + [("B1-2", "tetpos", "CAVB", "CASS_B")] * 1
        )
        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [f"{r[0]}_{r[1]}" for r in rows],
                    "patient_id": [r[0] for r in rows],
                    "enrichment_method": [r[1] for r in rows],
                    "CDR3_alpha": [r[2] for r in rows],
                    "CDR3_beta": [r[3] for r in rows],
                }
            ),
        )
        return adata

    def test_donor_columns_present_when_patient_id_set(self):
        adata = self._make_donor_method_adata()
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")

        assert "n_donors" in clonotypes.columns
        assert "donors" in clonotypes.columns

        clone_a = clonotypes[clonotypes["CDR3ab"] == "CAVA_CASS_A"].iloc[0]
        assert clone_a["n_donors"] == 2
        assert set(clone_a["donors"].split(";")) == {"B1-2", "B1-3"}

        clone_b = clonotypes[clonotypes["CDR3ab"] == "CAVB_CASS_B"].iloc[0]
        assert clone_b["n_donors"] == 1

    def test_method_columns_present_when_enrichment_method_set(self):
        adata = self._make_donor_method_adata()
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")

        assert "n_methods" in clonotypes.columns
        clone_a = clonotypes[clonotypes["CDR3ab"] == "CAVA_CASS_A"].iloc[0]
        # AIMpos seen in both donors, IFNpos+tetpos in B1-2 only -> 3 distinct
        assert clone_a["n_methods"] == 3
        assert set(clone_a["methods"].split(";")) == {"AIMpos", "IFNpos", "tetpos"}

    def test_methods_per_donor_dict_and_max(self):
        import json

        adata = self._make_donor_method_adata()
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")

        assert "methods_per_donor" in clonotypes.columns
        assert "max_methods_per_donor" in clonotypes.columns

        clone_a = clonotypes[clonotypes["CDR3ab"] == "CAVA_CASS_A"].iloc[0]
        d = json.loads(clone_a["methods_per_donor"])
        assert sorted(d["B1-2"]) == ["AIMpos", "IFNpos", "tetpos"]
        assert d["B1-3"] == ["AIMpos"]
        # Filter knob (planned --min-methods-per-donor) reads this:
        assert clone_a["max_methods_per_donor"] == 3

        clone_b = clonotypes[clonotypes["CDR3ab"] == "CAVB_CASS_B"].iloc[0]
        d_b = json.loads(clone_b["methods_per_donor"])
        assert d_b == {"B1-2": ["tetpos"]}
        assert clone_b["max_methods_per_donor"] == 1

    def _make_minimal_adata_no_donor_method(self):
        import anndata as ad
        import numpy as np

        n = 4
        return ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": ["A", "A", "B", "B"],
                    "CDR3_alpha": ["CAVA"] * n,
                    "CDR3_beta": ["CASS_A"] * n,
                }
            ),
        )

    def test_donor_method_columns_absent_when_neither_set(self):
        adata = self._make_minimal_adata_no_donor_method()
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")
        for col in (
            "donors",
            "n_donors",
            "methods",
            "n_methods",
            "methods_per_donor",
            "max_methods_per_donor",
        ):
            assert col not in clonotypes.columns

    def test_timepoint_apc_tissue_aggregations(self):
        """#9 chunk 2 — per-axis timepoint/apc/tissue aggregations."""
        import anndata as ad
        import numpy as np

        # Clone A: B1-2 across (D7,mDC), (D14,mDC), (D14,B-LCL); B1-3 at (D7,mDC).
        # Clone B: B1-2 only at (D7,mDC).
        rows = (
            [("B1-2", "D7", "mDC", "blood", "CAVA", "CASS_A")] * 3
            + [("B1-2", "D14", "mDC", "blood", "CAVA", "CASS_A")] * 2
            + [("B1-2", "D14", "B-LCL", "blood", "CAVA", "CASS_A")] * 2
            + [("B1-3", "D7", "mDC", "blood", "CAVA", "CASS_A")] * 2
            + [("B1-2", "D7", "mDC", "blood", "CAVB", "CASS_B")] * 1
        )
        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [f"{r[0]}_{r[1]}_{r[2]}" for r in rows],
                    "patient_id": [r[0] for r in rows],
                    "timepoint": [r[1] for r in rows],
                    "apc_type": [r[2] for r in rows],
                    "tissue": [r[3] for r in rows],
                    "CDR3_alpha": [r[4] for r in rows],
                    "CDR3_beta": [r[5] for r in rows],
                }
            ),
        )

        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")

        # Single-axis lists / counts
        for col in (
            "timepoints", "n_timepoints",
            "apcs", "n_apcs",
            "tissues", "n_tissues",
        ):
            assert col in clonotypes.columns, col

        clone_a = clonotypes[clonotypes["CDR3ab"] == "CAVA_CASS_A"].iloc[0]
        assert set(clone_a["timepoints"].split(";")) == {"D7", "D14"}
        assert clone_a["n_timepoints"] == 2
        assert set(clone_a["apcs"].split(";")) == {"mDC", "B-LCL"}
        assert clone_a["n_apcs"] == 2
        assert clone_a["n_tissues"] == 1

        # Nested per-donor
        import json
        for col in ("timepoints_per_donor", "max_timepoints_per_donor",
                    "apcs_per_donor", "max_apcs_per_donor"):
            assert col in clonotypes.columns, col

        tp_pd = json.loads(clone_a["timepoints_per_donor"])
        assert sorted(tp_pd["B1-2"]) == ["D14", "D7"]
        assert tp_pd["B1-3"] == ["D7"]
        assert clone_a["max_timepoints_per_donor"] == 2

        apc_pd = json.loads(clone_a["apcs_per_donor"])
        assert sorted(apc_pd["B1-2"]) == ["B-LCL", "mDC"]
        assert apc_pd["B1-3"] == ["mDC"]
        assert clone_a["max_apcs_per_donor"] == 2

        clone_b = clonotypes[clonotypes["CDR3ab"] == "CAVB_CASS_B"].iloc[0]
        assert clone_b["max_timepoints_per_donor"] == 1
        assert clone_b["max_apcs_per_donor"] == 1

    def test_per_donor_axis_columns_only_when_both_axes_set(self):
        """timepoints_per_donor only appears if both patient_id and timepoint
        are populated; same for apcs_per_donor."""
        import anndata as ad
        import numpy as np

        # patient_id + apc_type but no timepoint
        n = 4
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": ["S1", "S1", "S2", "S2"],
                    "patient_id": ["B1-2", "B1-2", "B1-3", "B1-3"],
                    "apc_type": ["mDC", "mDC", "B-LCL", "B-LCL"],
                    "CDR3_alpha": ["CAVA"] * n,
                    "CDR3_beta": ["CASS_A"] * n,
                }
            ),
        )
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")
        assert "apcs_per_donor" in clonotypes.columns
        assert "timepoints_per_donor" not in clonotypes.columns
        assert "max_timepoints_per_donor" not in clonotypes.columns

    def test_axis_columns_absent_when_unset(self):
        """No timepoint/apc/tissue columns when neither sample-sheet field set."""
        import anndata as ad
        import numpy as np

        n = 4
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": ["S1"] * n,
                    "CDR3_alpha": ["CAVA"] * n,
                    "CDR3_beta": ["CASS_A"] * n,
                }
            ),
        )
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")
        for col in (
            "timepoints", "n_timepoints",
            "apcs", "n_apcs",
            "tissues", "n_tissues",
            "timepoints_per_donor", "max_timepoints_per_donor",
            "apcs_per_donor", "max_apcs_per_donor",
        ):
            assert col not in clonotypes.columns, col

    def test_methods_per_donor_dropped_when_only_one_axis_set(self):
        """methods_per_donor / max_methods_per_donor only appear when BOTH
        patient_id and enrichment_method are present."""
        import anndata as ad
        import numpy as np

        # patient_id present, enrichment_method missing
        n = 4
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": ["A", "A", "B", "B"],
                    "patient_id": ["B1-2", "B1-2", "B1-3", "B1-3"],
                    "CDR3_alpha": ["CAVA"] * n,
                    "CDR3_beta": ["CASS_A"] * n,
                }
            ),
        )
        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab")
        assert "n_donors" in clonotypes.columns
        assert "methods_per_donor" not in clonotypes.columns
        assert "max_methods_per_donor" not in clonotypes.columns


class TestGetClonotypeSummary:
    """Tests for get_clonotype_summary function."""

    def test_summary_statistics(self, sample_clonotypes_df):
        """Summary should have correct statistics."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        assert "n_clonotypes" in summary
        assert summary["n_clonotypes"] == 5

        assert "n_cells" in summary
        assert summary["n_cells"] == 29  # 15 + 8 + 3 + 2 + 1

        assert "n_singletons" in summary
        assert summary["n_singletons"] == 1  # only last one with count=1

        assert "n_expanded" in summary
        assert summary["n_expanded"] == 4  # 4 clones with count > 1

    def test_multi_sample_clones(self, sample_clonotypes_df):
        """Should count multi-sample clones."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        assert "n_multi_sample" in summary
        # First clone has n_samples=2, others have n_samples=1
        assert summary["n_multi_sample"] == 1


class TestExportClonotypesAirr:
    """Tests for AIRR format export."""

    def test_export_airr(self, sample_clonotypes_df, temp_dir):
        """Export to AIRR format should work."""
        output_path = temp_dir / "clonotypes.tsv"

        export_clonotypes_airr(sample_clonotypes_df, str(output_path))

        assert output_path.exists()

        # Load and check format
        airr_df = pd.read_csv(output_path, sep="\t")
        assert "junction_aa_tra" in airr_df.columns or len(airr_df) > 0
        assert "junction_aa_trb" in airr_df.columns or len(airr_df) > 0
        assert "clone_count" in airr_df.columns
        assert "productive" in airr_df.columns

    def test_airr_column_mapping(self, sample_clonotypes_df, temp_dir):
        """Check AIRR column mapping."""
        output_path = temp_dir / "clonotypes.tsv"

        export_clonotypes_airr(sample_clonotypes_df, str(output_path))

        airr_df = pd.read_csv(output_path, sep="\t")

        # Check that CDR3 sequences are mapped correctly
        if "junction_aa_trb" in airr_df.columns:
            assert airr_df["junction_aa_trb"].iloc[0] == "CASSLGQAYEQYF"


class TestCalculateCloneFrequencies:
    """Tests for calculate_clone_frequencies function."""

    @pytest.fixture
    def adata_with_clone_info(self, sample_adata):
        """Create AnnData with clone information."""
        adata = sample_adata.copy()
        adata.obs["sample"] = ["S1"] * 50 + ["S2"] * 50
        adata.obs["CDR3_alpha"] = ["CASSL"] * 30 + ["CAVSD"] * 70
        adata.obs["CDR3_beta"] = ["CASSF"] * 30 + ["CASRG"] * 70
        adata.obs["CDR3ab"] = adata.obs["CDR3_alpha"] + "_" + adata.obs["CDR3_beta"]
        adata.obs["is_complete_clone"] = True
        return adata

    def test_calculate_frequencies_basic(self, adata_with_clone_info, sample_clonotypes_df):
        """Test basic frequency calculation."""
        result = calculate_clone_frequencies(sample_clonotypes_df, adata_with_clone_info)

        assert "sample_frequencies" in result.columns
        assert "n_conditions_present" in result.columns

    def test_frequency_values(self, adata_with_clone_info):
        """Test frequency values are calculated correctly."""
        # Create clonotypes from this adata
        clonotypes = aggregate_clonotypes(adata_with_clone_info, group_by="CDR3ab")

        result = calculate_clone_frequencies(clonotypes, adata_with_clone_info)

        # All frequencies should be between 0 and 1
        if "max_frequency" in result.columns:
            assert all(result["max_frequency"] >= 0)
            assert all(result["max_frequency"] <= 1)


class TestAggregateClonotypesExtended:
    """Extended tests for aggregate_clonotypes edge cases."""

    @pytest.fixture
    def adata_with_vdj_genes(self, sample_adata):
        """Create AnnData with VDJ gene columns."""
        adata = sample_adata.copy()
        adata.obs["sample"] = "S1"
        adata.obs["CDR3_alpha"] = "CASSL"
        adata.obs["CDR3_beta"] = "CASSF"
        adata.obs["TRA_1_v_gene"] = "TRAV1"
        adata.obs["TRA_1_j_gene"] = "TRAJ1"
        adata.obs["TRA_1_c_gene"] = "TRAC"
        adata.obs["TRB_1_v_gene"] = "TRBV2"
        adata.obs["TRB_1_j_gene"] = "TRBJ2"
        adata.obs["TRB_1_c_gene"] = "TRBC1"
        adata.obs["TRA_1_vdj_aa"] = "MRLVTSGF"
        adata.obs["TRA_1_vdj_nt"] = "ATGCGT"
        adata.obs["TRB_1_vdj_aa"] = "MGVTSGHD"
        adata.obs["TRB_1_vdj_nt"] = "ATGGGT"
        adata.obs["TRA_1_umis"] = 100
        adata.obs["TRB_1_umis"] = 200
        adata.obs["TRA_1_reads"] = 1000
        adata.obs["TRB_1_reads"] = 2000
        adata.obs["TRA_1_contig_id"] = "contig_1"
        adata.obs["TRB_1_contig_id"] = "contig_2"
        adata.obs["antigen_description"] = "TestAntigen"
        adata.obs["source"] = "culture"
        return adata

    def test_aggregate_with_vdj_genes(self, adata_with_vdj_genes):
        """Test aggregation includes VDJ gene information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert len(clonotypes) > 0
        assert "alpha_v_gene" in clonotypes.columns
        assert "beta_v_gene" in clonotypes.columns
        assert clonotypes["alpha_v_gene"].iloc[0] == "TRAV1"

    def test_aggregate_with_vdj_sequences(self, adata_with_vdj_genes):
        """Test aggregation includes VDJ sequences."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "VDJ_alpha_aa" in clonotypes.columns
        assert "VDJ_beta_aa" in clonotypes.columns
        assert clonotypes["VDJ_alpha_aa"].iloc[0] == "MRLVTSGF"

    def test_aggregate_with_umi_metrics(self, adata_with_vdj_genes):
        """Test aggregation includes UMI metrics."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "alpha_umis_mean" in clonotypes.columns
        assert "beta_reads_sum" in clonotypes.columns

    def test_aggregate_with_contig_ids(self, adata_with_vdj_genes):
        """Test aggregation includes contig IDs."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "alpha_contig_ids" in clonotypes.columns
        assert "beta_contig_ids" in clonotypes.columns

    def test_aggregate_with_antigen_info(self, adata_with_vdj_genes):
        """Test aggregation includes antigen information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "antigens" in clonotypes.columns
        assert "n_antigens" in clonotypes.columns
        assert "n_conditions" in clonotypes.columns
        assert "TestAntigen" in clonotypes["antigens"].iloc[0]
        assert clonotypes["n_conditions"].iloc[0] == clonotypes["n_antigens"].iloc[0]

    def test_aggregate_with_source_info(self, adata_with_vdj_genes):
        """Test aggregation includes source information."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3ab")

        assert "sources" in clonotypes.columns
        assert clonotypes["sources"].iloc[0] == "culture"

    def test_aggregate_cdr3b_only_with_alpha(self, adata_with_vdj_genes):
        """Test CDR3b_only mode still captures alpha when available."""
        clonotypes = aggregate_clonotypes(adata_with_vdj_genes, group_by="CDR3b_only")

        assert len(clonotypes) > 0
        assert "CDR3_beta" in clonotypes.columns
        # CDR3_alpha should also be captured when available
        assert "CDR3_alpha" in clonotypes.columns


class TestGetClonotypeSummaryExtended:
    """Extended tests for get_clonotype_summary edge cases."""

    def test_summary_without_n_samples(self):
        """Test summary when n_samples column is missing."""
        df = pd.DataFrame(
            {
                "CDR3ab": ["A", "B"],
                "cell_count": [5, 3],
                # No n_samples column
            }
        )

        summary = get_clonotype_summary(df)

        assert "n_multi_sample" in summary
        assert summary["n_multi_sample"] == 0  # Default when column missing


class TestRealisticClonotyping:
    """Tests using realistic clonotype data matching real CellRanger outputs."""

    def test_clonotype_cdr3_format(self, sample_clonotypes_df):
        """Test that clonotype CDR3 sequences follow IMGT conventions."""
        df = sample_clonotypes_df

        # All alpha CDR3s should start with CA (conserved cysteine + alanine)
        for cdr3 in df["CDR3_alpha"]:
            assert cdr3.startswith("CA"), f"Alpha CDR3 {cdr3} doesn't start with CA"

        # All beta CDR3s should start with CAS (conserved cysteine + alanine + serine)
        for cdr3 in df["CDR3_beta"]:
            assert cdr3.startswith("CAS"), f"Beta CDR3 {cdr3} doesn't start with CAS"

    def test_clonotype_gene_naming(self, sample_clonotypes_df):
        """Test that V/J/C gene names follow IMGT nomenclature."""
        df = sample_clonotypes_df

        # V genes should be TRAV/TRBV with numeric suffix
        for v_gene in df["alpha_v_gene"]:
            assert v_gene.startswith("TRAV"), f"Alpha V gene {v_gene} doesn't start with TRAV"
        for v_gene in df["beta_v_gene"]:
            assert v_gene.startswith("TRBV"), f"Beta V gene {v_gene} doesn't start with TRBV"

        # J genes should be TRAJ/TRBJ
        for j_gene in df["alpha_j_gene"]:
            assert j_gene.startswith("TRAJ"), f"Alpha J gene {j_gene} doesn't start with TRAJ"

        # C genes should be TRAC/TRBC
        for c_gene in df["alpha_c_gene"]:
            assert c_gene == "TRAC", f"Alpha C gene {c_gene} should be TRAC"
        for c_gene in df["beta_c_gene"]:
            assert c_gene.startswith("TRBC"), f"Beta C gene {c_gene} doesn't start with TRBC"

    def test_realistic_cell_counts(self, sample_clonotypes_df):
        """Test that cell counts represent typical clonal expansion patterns."""
        df = sample_clonotypes_df

        # Most clones should be small (1-10 cells)
        small_clones = df[df["cell_count"] <= 10]
        assert len(small_clones) > 0, "Should have small clones"

        # Some clones can be expanded
        expanded = df[df["cell_count"] > 1]
        assert len(expanded) > 0, "Should have expanded clones"

        # Cell counts should be positive integers
        assert all(df["cell_count"] > 0)
        assert all(df["cell_count"] == df["cell_count"].astype(int))

    def test_frequency_range(self, sample_clonotypes_df):
        """Test that clone frequencies are in valid range."""
        df = sample_clonotypes_df

        # All frequencies should be between 0 and 1
        assert all(df["max_frequency"] >= 0)
        assert all(df["max_frequency"] <= 1)

        # Larger clones should have higher frequencies
        sorted_by_count = df.sort_values("cell_count", ascending=False)
        sorted_by_freq = df.sort_values("max_frequency", ascending=False)

        # Top clone by count should also be among top by frequency
        top_by_count = sorted_by_count.iloc[0]["CDR3ab"]
        top_3_by_freq = sorted_by_freq.head(3)["CDR3ab"].tolist()
        assert top_by_count in top_3_by_freq

    def test_tcell_type_annotations(self, sample_clonotypes_df):
        """Test T cell type consensus annotations."""
        df = sample_clonotypes_df

        # Valid T cell types
        valid_types = [
            "Confident CD4+",
            "Confident CD8+",
            "Likely CD4+",
            "Likely CD8+",
            "Unknown",
            "Mixed",
        ]

        for tcell_type in df["Tcell_type_consensus"]:
            assert tcell_type in valid_types, f"Unknown T cell type: {tcell_type}"

    def test_CDR3ab_format(self, sample_clonotypes_df):
        """Test that clone IDs are formatted as CDR3_alpha_CDR3_beta."""
        df = sample_clonotypes_df

        for idx, row in df.iterrows():
            expected_id = f"{row['CDR3_alpha']}_{row['CDR3_beta']}"
            assert row["CDR3ab"] == expected_id, (
                f"Clone ID {row['CDR3ab']} doesn't match expected {expected_id}"
            )

    def test_antigen_name_fallback(self):
        """antigen_name should be used when antigen_description is missing."""
        import anndata as ad

        obs = pd.DataFrame(
            {
                "CDR3_alpha": ["CAVAAA", "CAVAAA"],
                "CDR3_beta": ["CASSBBB", "CASSBBB"],
                "sample": ["S1", "S1"],
                "antigen_name": ["CMV", "CMV"],
            },
            index=["cell1", "cell2"],
        )
        adata = ad.AnnData(obs=obs)

        clonotypes = aggregate_clonotypes(adata, group_by="CDR3ab", show_progress=False)

        assert "antigens" in clonotypes.columns
        assert clonotypes["antigens"].iloc[0] == "CMV"


class TestClonotypeSummaryRealistic:
    """Tests for clonotype summary with realistic data."""

    def test_summary_statistics_realistic(self, sample_clonotypes_df):
        """Test summary statistics with realistic clonotype data."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        # Number of clonotypes
        assert summary["n_clonotypes"] == len(sample_clonotypes_df)

        # Total cells should be sum of cell_count
        assert summary["n_cells"] == sample_clonotypes_df["cell_count"].sum()

        # Singletons are clones with count == 1
        expected_singletons = (sample_clonotypes_df["cell_count"] == 1).sum()
        assert summary["n_singletons"] == expected_singletons

        # Expanded are clones with count > 1
        expected_expanded = (sample_clonotypes_df["cell_count"] > 1).sum()
        assert summary["n_expanded"] == expected_expanded

    def test_multi_sample_detection(self, sample_clonotypes_df):
        """Test detection of clones present in multiple samples."""
        summary = get_clonotype_summary(sample_clonotypes_df)

        # Count clones with n_samples > 1
        expected_multi = (sample_clonotypes_df["n_samples"] > 1).sum()
        assert summary["n_multi_sample"] == expected_multi


class TestBuildCloneSampleLong:
    """#20 chunk 1 — long-format (clone, sample) companion table."""

    def _make_adata(self, with_donor=True, with_method=True, with_umis=True):
        """Two clones across two samples (different donors + methods)."""
        import anndata as ad
        import numpy as np

        rows = (
            # (sample, donor, method, alpha, beta, tra_umis, trb_umis)
            [("B1-2_AIM", "B1-2", "AIMpos", "CAVA", "CASS_A", 5, 7)] * 4
            + [("B1-2_tet", "B1-2", "tetpos", "CAVA", "CASS_A", 6, 8)] * 2
            + [("B1-3_AIM", "B1-3", "AIMpos", "CAVB", "CASS_B", 4, 6)] * 3
            + [("B1-2_AIM", "B1-2", "AIMpos", "CAVB", "CASS_B", 5, 7)] * 1
        )
        n = len(rows)
        obs_data = {
            "sample": [r[0] for r in rows],
            "CDR3_alpha": [r[3] for r in rows],
            "CDR3_beta": [r[4] for r in rows],
        }
        if with_donor:
            obs_data["patient_id"] = [r[1] for r in rows]
        if with_method:
            obs_data["enrichment_method"] = [r[2] for r in rows]
        if with_umis:
            obs_data["TRA_1_umis"] = [r[5] for r in rows]
            obs_data["TRB_1_umis"] = [r[6] for r in rows]
        return ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(obs_data),
        )

    def test_basic_long_table_shape(self):
        adata = self._make_adata()
        long = build_clone_sample_long(adata)
        # 4 (clone, sample) pairs:
        # CAVA_CASS_A in B1-2_AIM, CAVA_CASS_A in B1-2_tet,
        # CAVB_CASS_B in B1-3_AIM, CAVB_CASS_B in B1-2_AIM
        assert len(long) == 4
        assert set(long.columns) >= {
            "CDR3ab", "sample", "donor", "method",
            "cells", "frequency", "n_alpha_umis", "n_beta_umis",
        }

    def test_cell_counts_match_input(self):
        adata = self._make_adata()
        long = build_clone_sample_long(adata)
        row_a_aim = long[
            (long["CDR3ab"] == "CAVA_CASS_A") & (long["sample"] == "B1-2_AIM")
        ].iloc[0]
        assert row_a_aim["cells"] == 4
        row_a_tet = long[
            (long["CDR3ab"] == "CAVA_CASS_A") & (long["sample"] == "B1-2_tet")
        ].iloc[0]
        assert row_a_tet["cells"] == 2

    def test_donor_method_propagation(self):
        adata = self._make_adata()
        long = build_clone_sample_long(adata)
        for _, row in long.iterrows():
            sample = row["sample"]
            if sample.startswith("B1-2"):
                assert row["donor"] == "B1-2"
            else:
                assert row["donor"] == "B1-3"
            if "AIM" in sample:
                assert row["method"] == "AIMpos"
            elif "tet" in sample:
                assert row["method"] == "tetpos"

    def test_umi_sums(self):
        adata = self._make_adata()
        long = build_clone_sample_long(adata)
        row = long[
            (long["CDR3ab"] == "CAVA_CASS_A") & (long["sample"] == "B1-2_AIM")
        ].iloc[0]
        # 4 cells x 5 alpha umis each = 20
        assert row["n_alpha_umis"] == 20
        assert row["n_beta_umis"] == 28

    def test_frequency_uses_complete_clones_when_available(self):
        adata = self._make_adata()
        # Mark some cells incomplete; frequency denominator should drop them
        adata.obs["is_complete_clone"] = True
        long = build_clone_sample_long(adata)
        # B1-2_AIM has 5 complete cells (4 of clone A + 1 of clone B)
        # Clone A in B1-2_AIM: 4 / 5 = 0.8
        row = long[
            (long["CDR3ab"] == "CAVA_CASS_A") & (long["sample"] == "B1-2_AIM")
        ].iloc[0]
        assert abs(row["frequency"] - 4 / 5) < 1e-9

    def test_axis_columns_optional(self):
        """donor/method should not appear when axis fields aren't on obs."""
        adata = self._make_adata(with_donor=False, with_method=False)
        long = build_clone_sample_long(adata)
        assert "donor" not in long.columns
        assert "method" not in long.columns

    def test_skips_empty_cdr3(self):
        """Cells with empty CDR3ab are excluded from the long table."""
        import anndata as ad
        import numpy as np

        adata = ad.AnnData(
            X=np.zeros((4, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": ["S1", "S1", "S2", "S2"],
                    "CDR3_alpha": ["CAVA", "", None, "CAVB"],
                    "CDR3_beta": ["CASS_A", "", None, "CASS_B"],
                }
            ),
        )
        long = build_clone_sample_long(adata)
        # Only the two non-empty cells survive.
        assert len(long) == 2
        assert set(long["CDR3ab"]) == {"CAVA_CASS_A", "CAVB_CASS_B"}

    def test_raises_when_no_cdr3_columns(self):
        from tcrsift.validation import TCRsiftValidationError
        import anndata as ad
        import numpy as np

        adata = ad.AnnData(
            X=np.zeros((2, 1), dtype=np.float32),
            obs=pd.DataFrame({"sample": ["S1", "S2"]}),
        )
        with pytest.raises(TCRsiftValidationError, match="CDR3"):
            build_clone_sample_long(adata)


class TestBuildPerMethodRankings:
    """#20 chunk 2 — per-(donor, method) ranked CSVs."""

    def _setup(self):
        from tcrsift.clonotype import build_clone_sample_long
        import anndata as ad
        import numpy as np

        # Two donors, three methods. Clone A is shared across both donors;
        # clone B is private to B1-2; clone C is private to B1-3.
        rows = (
            # (sample, donor, method, alpha, beta)
            [("B1-2_AIM", "B1-2", "AIMpos", "CAVA", "CASS_A")] * 10
            + [("B1-2_tet", "B1-2", "tetpos", "CAVA", "CASS_A")] * 4
            + [("B1-2_IFN", "B1-2", "IFNpos", "CAVA", "CASS_A")] * 2
            + [("B1-3_AIM", "B1-3", "AIMpos", "CAVA", "CASS_A")] * 6
            + [("B1-2_AIM", "B1-2", "AIMpos", "CAVB", "CASS_B")] * 8
            + [("B1-3_AIM", "B1-3", "AIMpos", "CAVC", "CASS_C")] * 5
        )
        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [r[0] for r in rows],
                    "patient_id": [r[1] for r in rows],
                    "enrichment_method": [r[2] for r in rows],
                    "CDR3_alpha": [r[3] for r in rows],
                    "CDR3_beta": [r[4] for r in rows],
                }
            ),
        )
        long_df = build_clone_sample_long(adata)
        # Synthetic filtered table — A and B pass; C does not.
        filtered = pd.DataFrame(
            {
                "CDR3ab": ["CAVA_CASS_A", "CAVB_CASS_B"],
                "tier": ["tier1", "tier2"],
                "n_donors": [2, 1],
                "cell_count": [22, 8],
                "max_frequency": [0.5, 0.4],
            }
        )
        return filtered, long_df

    def test_emits_one_table_per_donor_method(self):
        from tcrsift.clonotype import build_per_method_rankings

        filtered, long_df = self._setup()
        rankings = build_per_method_rankings(filtered, long_df, top_n=100)

        # Filtered has clones A and B. A appears in (B1-2,AIM/tet/IFN) and
        # (B1-3, AIM); B appears in (B1-2, AIM). So expected keys:
        expected_keys = {
            ("B1-2", "AIMpos"),
            ("B1-2", "tetpos"),
            ("B1-2", "IFNpos"),
            ("B1-3", "AIMpos"),
        }
        assert set(rankings.keys()) == expected_keys

    def test_per_table_content(self):
        from tcrsift.clonotype import build_per_method_rankings

        filtered, long_df = self._setup()
        rankings = build_per_method_rankings(filtered, long_df, top_n=100)

        b12_aim = rankings[("B1-2", "AIMpos")]
        # Clones A and B both appear in (B1-2, AIMpos)
        assert set(b12_aim["CDR3ab"]) == {"CAVA_CASS_A", "CAVB_CASS_B"}
        # Tier and sharing should be carried through
        assert "tier" in b12_aim.columns
        assert "sharing" in b12_aim.columns
        a_row = b12_aim[b12_aim["CDR3ab"] == "CAVA_CASS_A"].iloc[0]
        assert a_row["tier"] == "tier1"
        assert a_row["sharing"] == "public"  # n_donors == 2
        b_row = b12_aim[b12_aim["CDR3ab"] == "CAVB_CASS_B"].iloc[0]
        assert b_row["sharing"] == "private"  # n_donors == 1
        # Donor / method shouldn't appear in row data — encoded in dict key
        assert "donor" not in b12_aim.columns
        assert "method" not in b12_aim.columns

    def test_only_filtered_clones_appear(self):
        from tcrsift.clonotype import build_per_method_rankings

        filtered, long_df = self._setup()
        rankings = build_per_method_rankings(filtered, long_df, top_n=100)

        # Clone C wasn't in `filtered`, so should not appear in any output
        for key, df in rankings.items():
            assert "CAVC_CASS_C" not in set(df["CDR3ab"]), key

    def test_top_n_truncation(self):
        from tcrsift.clonotype import build_per_method_rankings

        filtered, long_df = self._setup()
        rankings = build_per_method_rankings(filtered, long_df, top_n=1)
        for key, df in rankings.items():
            assert len(df) <= 1

    def test_returns_empty_when_method_axis_absent(self):
        """No method column on long_df -> empty dict."""
        from tcrsift.clonotype import build_per_method_rankings

        long_df = pd.DataFrame(
            {"CDR3ab": ["CAVA_CASS_A"], "sample": ["S1"], "cells": [10], "frequency": [1.0]}
        )
        filtered = pd.DataFrame({"CDR3ab": ["CAVA_CASS_A"]})
        rankings = build_per_method_rankings(filtered, long_df)
        assert rankings == {}

    def test_returns_empty_when_filtered_is_empty(self):
        from tcrsift.clonotype import build_per_method_rankings

        _, long_df = self._setup()
        rankings = build_per_method_rankings(pd.DataFrame(), long_df)
        assert rankings == {}

    def test_no_donor_axis_uses_synthetic_all(self):
        """When method is on long_df but donor isn't, fold under
        donor='all' so file naming stays uniform."""
        from tcrsift.clonotype import build_per_method_rankings

        long_df = pd.DataFrame(
            {
                "CDR3ab": ["CAVA_CASS_A", "CAVB_CASS_B"],
                "sample": ["S1", "S2"],
                "method": ["AIMpos", "tetpos"],
                "cells": [10, 5],
                "frequency": [0.1, 0.05],
            }
        )
        filtered = pd.DataFrame({"CDR3ab": ["CAVA_CASS_A", "CAVB_CASS_B"]})
        rankings = build_per_method_rankings(filtered, long_df)
        assert ("all", "AIMpos") in rankings
        assert ("all", "tetpos") in rankings


class TestBuildMethodOverlapMatrices:
    """#27 chunk 3 — method × method overlap matrices per donor."""

    def _setup(self):
        from tcrsift.clonotype import build_clone_sample_long
        import anndata as ad
        import numpy as np

        # B1-2: clone A in (AIM, tet, IFN); clone B in (AIM, tet);
        # clone C in (IFN). B1-3: clone A in (AIM); clone D in (AIM, tet).
        rows = (
            [("B1-2_AIM", "B1-2", "AIMpos", "CAVA", "CASS_A")] * 4
            + [("B1-2_tet", "B1-2", "tetpos", "CAVA", "CASS_A")] * 2
            + [("B1-2_IFN", "B1-2", "IFNpos", "CAVA", "CASS_A")] * 2
            + [("B1-2_AIM", "B1-2", "AIMpos", "CAVB", "CASS_B")] * 3
            + [("B1-2_tet", "B1-2", "tetpos", "CAVB", "CASS_B")] * 2
            + [("B1-2_IFN", "B1-2", "IFNpos", "CAVC", "CASS_C")] * 2
            + [("B1-3_AIM", "B1-3", "AIMpos", "CAVA", "CASS_A")] * 4
            + [("B1-3_AIM", "B1-3", "AIMpos", "CAVD", "CASS_D")] * 3
            + [("B1-3_tet", "B1-3", "tetpos", "CAVD", "CASS_D")] * 2
        )
        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [r[0] for r in rows],
                    "patient_id": [r[1] for r in rows],
                    "enrichment_method": [r[2] for r in rows],
                    "CDR3_alpha": [r[3] for r in rows],
                    "CDR3_beta": [r[4] for r in rows],
                }
            ),
        )
        long_df = build_clone_sample_long(adata)
        # All four clones pass.
        filtered = pd.DataFrame(
            {
                "CDR3ab": [
                    "CAVA_CASS_A",
                    "CAVB_CASS_B",
                    "CAVC_CASS_C",
                    "CAVD_CASS_D",
                ]
            }
        )
        return filtered, long_df

    def test_jaccard_default(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        filtered, long_df = self._setup()
        out = build_method_overlap_matrices(filtered, long_df, similarity="jaccard")
        assert set(out.keys()) == {"B1-2", "B1-3"}

        b12 = out["B1-2"]
        # B1-2 methods: AIMpos, IFNpos, tetpos -> sorted lex order
        assert list(b12.index) == ["AIMpos", "IFNpos", "tetpos"]
        # AIMpos has {A, B}; tetpos has {A, B}; IFNpos has {A, C}.
        # AIM ∩ tet = {A, B} = 2 ; AIM ∪ tet = {A, B} = 2 -> 1.0
        assert abs(b12.loc["AIMpos", "tetpos"] - 1.0) < 1e-9
        # AIM ∩ IFN = {A} = 1 ; AIM ∪ IFN = {A, B, C} = 3 -> 1/3
        assert abs(b12.loc["AIMpos", "IFNpos"] - 1 / 3) < 1e-9
        # tet ∩ IFN = {A} = 1 ; tet ∪ IFN = {A, B, C} = 3 -> 1/3
        assert abs(b12.loc["tetpos", "IFNpos"] - 1 / 3) < 1e-9
        # diagonals are 1.0
        for m in b12.index:
            assert b12.loc[m, m] == 1.0
        # symmetric
        for i in b12.index:
            for j in b12.columns:
                assert b12.loc[i, j] == b12.loc[j, i]

    def test_count_metric(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        filtered, long_df = self._setup()
        out = build_method_overlap_matrices(filtered, long_df, similarity="count")
        b12 = out["B1-2"]
        # AIM has {A, B} -> diagonal = 2
        assert int(b12.loc["AIMpos", "AIMpos"]) == 2
        # AIM ∩ IFN = 1
        assert int(b12.loc["AIMpos", "IFNpos"]) == 1

    def test_dice_metric(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        filtered, long_df = self._setup()
        out = build_method_overlap_matrices(filtered, long_df, similarity="dice")
        b12 = out["B1-2"]
        # AIM = {A, B}, IFN = {A, C}; dice = 2*1 / (2+2) = 0.5
        assert abs(b12.loc["AIMpos", "IFNpos"] - 0.5) < 1e-9

    def test_invalid_similarity_raises(self):
        from tcrsift.clonotype import build_method_overlap_matrices
        from tcrsift.validation import TCRsiftValidationError

        with pytest.raises(TCRsiftValidationError, match="similarity"):
            build_method_overlap_matrices(
                pd.DataFrame({"CDR3ab": ["A"]}),
                pd.DataFrame({"CDR3ab": ["A"], "method": ["X"], "donor": ["D"]}),
                similarity="cosine",
            )

    def test_returns_empty_when_method_axis_absent(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        long_df = pd.DataFrame({"CDR3ab": ["A"], "sample": ["S1"]})
        filtered = pd.DataFrame({"CDR3ab": ["A"]})
        assert build_method_overlap_matrices(filtered, long_df) == {}

    def test_returns_empty_when_filtered_empty(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        _, long_df = self._setup()
        assert build_method_overlap_matrices(pd.DataFrame(), long_df) == {}

    def test_no_donor_axis_uses_synthetic_all(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        long_df = pd.DataFrame(
            {
                "CDR3ab": ["A", "A", "B"],
                "sample": ["S1", "S2", "S1"],
                "method": ["AIM", "tet", "AIM"],
                "cells": [1, 1, 1],
                "frequency": [0.1, 0.1, 0.1],
            }
        )
        filtered = pd.DataFrame({"CDR3ab": ["A", "B"]})
        out = build_method_overlap_matrices(filtered, long_df)
        assert "all" in out


class TestBuildMethodRecoveryTable:
    """#27 chunk 4 — per-(donor, method) recovery of a target tier."""

    def _setup(self):
        from tcrsift.clonotype import build_clone_sample_long
        import anndata as ad
        import numpy as np

        # 5 tier-1 target clones; B1-2's AIMpos catches all 5; B1-2's tetpos
        # catches only 2; B1-3's AIMpos catches 4; B1-3's tetpos catches 0.
        rows = []
        for cdr3 in ["A", "B", "C", "D", "E"]:
            rows.append(("B1-2_AIM", "B1-2", "AIMpos", f"CAV{cdr3}", f"CASS_{cdr3}"))
        for cdr3 in ["A", "B"]:
            rows.append(("B1-2_tet", "B1-2", "tetpos", f"CAV{cdr3}", f"CASS_{cdr3}"))
        for cdr3 in ["A", "B", "C", "D"]:
            rows.append(("B1-3_AIM", "B1-3", "AIMpos", f"CAV{cdr3}", f"CASS_{cdr3}"))
        # B1-3 tetpos sees only a non-target clone (Z); so recovery == 0.
        rows.append(("B1-3_tet", "B1-3", "tetpos", "CAVZ", "CASS_Z"))

        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [r[0] for r in rows],
                    "patient_id": [r[1] for r in rows],
                    "enrichment_method": [r[2] for r in rows],
                    "CDR3_alpha": [r[3] for r in rows],
                    "CDR3_beta": [r[4] for r in rows],
                }
            ),
        )
        long_df = build_clone_sample_long(adata)
        filtered = pd.DataFrame(
            {
                "CDR3ab": [f"CAV{c}_CASS_{c}" for c in "ABCDE"],
                "tier": ["tier1"] * 5,
            }
        )
        return filtered, long_df

    def test_recovery_counts(self):
        from tcrsift.clonotype import build_method_recovery_table

        filtered, long_df = self._setup()
        out = build_method_recovery_table(filtered, long_df, tier="tier1")
        assert set(out["donor"].unique()) == {"B1-2", "B1-3"}
        # B1-2 / AIMpos: all 5 of B1-2's target clones recovered.
        b12_aim = out[(out["donor"] == "B1-2") & (out["method"] == "AIMpos")].iloc[0]
        assert b12_aim["recovered"] == 5
        assert b12_aim["total"] == 5
        # B1-2 / tetpos: 2 recovered out of 5.
        b12_tet = out[(out["donor"] == "B1-2") & (out["method"] == "tetpos")].iloc[0]
        assert b12_tet["recovered"] == 2
        assert abs(b12_tet["fraction"] - 2 / 5) < 1e-9
        # B1-3 / AIMpos: 4 recovered.
        b13_aim = out[(out["donor"] == "B1-3") & (out["method"] == "AIMpos")].iloc[0]
        assert b13_aim["recovered"] == 4
        # B1-3 / tetpos: 0 recovered.
        b13_tet = out[(out["donor"] == "B1-3") & (out["method"] == "tetpos")].iloc[0]
        assert b13_tet["recovered"] == 0
        assert b13_tet["fraction"] == 0.0

    def test_returns_empty_when_method_absent(self):
        from tcrsift.clonotype import build_method_recovery_table

        long_df = pd.DataFrame({"CDR3ab": ["A"], "sample": ["S1"]})
        filtered = pd.DataFrame({"CDR3ab": ["A"], "tier": ["tier1"]})
        out = build_method_recovery_table(filtered, long_df)
        assert out.empty

    def test_falls_back_to_all_filtered_when_no_tier_column(self):
        from tcrsift.clonotype import build_method_recovery_table

        filtered, long_df = self._setup()
        # Strip tier; recovery should still target the full filtered set.
        filtered_no_tier = filtered.drop(columns=["tier"])
        out = build_method_recovery_table(
            filtered_no_tier, long_df, tier="tier1"
        )
        # B1-2 / AIMpos still recovers all 5
        b12_aim = out[(out["donor"] == "B1-2") & (out["method"] == "AIMpos")].iloc[0]
        assert b12_aim["recovered"] == 5

    def test_donor_total_does_not_fall_back_to_cohort_target(self):
        """Review fix: when a donor has zero target clones in its long_df,
        `total` is 0 (not the cohort-wide target count). Prevents
        misleading "0/N" denominators that imply the donor was supposed
        to recover targets it never had access to."""
        from tcrsift.clonotype import build_method_recovery_table

        # Two donors. B1-2 has all 5 target clones; B1-3 has none of them
        # (only a non-target clone Z).
        long_df = pd.DataFrame(
            {
                "CDR3ab": [
                    "CAVA_CASS_A", "CAVB_CASS_B", "CAVC_CASS_C",
                    "CAVD_CASS_D", "CAVE_CASS_E",
                    "CAVZ_CASS_Z",
                ],
                "sample": ["S1"] * 5 + ["S2"],
                "donor": ["B1-2"] * 5 + ["B1-3"],
                "method": ["AIMpos"] * 5 + ["AIMpos"],
                "cells": [1] * 6,
                "frequency": [0.1] * 6,
            }
        )
        filtered = pd.DataFrame(
            {
                "CDR3ab": [f"CAV{c}_CASS_{c}" for c in "ABCDE"],
                "tier": ["tier1"] * 5,
            }
        )
        out = build_method_recovery_table(filtered, long_df, tier="tier1")
        b13 = out[out["donor"] == "B1-3"]
        # B1-3 has zero target clones in its long table -> total should be 0,
        # not 5.
        assert all(b13["total"] == 0)
        assert all(b13["recovered"] == 0)
        assert all(b13["fraction"] == 0.0)


class TestSingleMethodMethodOverlap:
    """Review polish: cover n_methods == 1 case for build_method_overlap_matrices."""

    def test_single_method_returns_one_by_one_diagonal_one(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        long_df = pd.DataFrame(
            {
                "CDR3ab": ["A", "B"],
                "sample": ["S1", "S1"],
                "donor": ["B1-2", "B1-2"],
                "method": ["AIMpos", "AIMpos"],
                "cells": [1, 1],
                "frequency": [0.1, 0.1],
            }
        )
        filtered = pd.DataFrame({"CDR3ab": ["A", "B"]})
        out = build_method_overlap_matrices(filtered, long_df, similarity="jaccard")
        assert "B1-2" in out
        m = out["B1-2"]
        assert m.shape == (1, 1)
        assert m.iat[0, 0] == 1.0

    def test_single_method_count_metric_returns_clone_count(self):
        from tcrsift.clonotype import build_method_overlap_matrices

        long_df = pd.DataFrame(
            {
                "CDR3ab": ["A", "B", "C"],
                "sample": ["S1", "S1", "S1"],
                "donor": ["B1-2", "B1-2", "B1-2"],
                "method": ["AIMpos", "AIMpos", "AIMpos"],
                "cells": [1, 1, 1],
                "frequency": [0.1, 0.1, 0.1],
            }
        )
        filtered = pd.DataFrame({"CDR3ab": ["A", "B", "C"]})
        out = build_method_overlap_matrices(filtered, long_df, similarity="count")
        assert int(out["B1-2"].iat[0, 0]) == 3
