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

"""End-to-end tests for TCRsift pipeline.

These tests exercise the complete workflow from raw VDJ data through
filtering and annotation, using realistic synthetic data.
"""

import anndata as ad
import numpy as np
import pandas as pd
import pytest


class TestEndToEndPipeline:
    """End-to-end pipeline tests with realistic synthetic data."""

    @pytest.fixture
    def synthetic_cellranger_data(self):
        """Generate realistic CellRanger VDJ output for 500 cells.

        Creates:
        - filtered_contig_annotations.csv with paired alpha/beta chains
        - Realistic distribution of clone sizes (few expanded, many rare)
        - Some cells with only single chain
        """
        np.random.seed(42)
        n_cells = 500

        # Generate clone distribution: few expanded, many rare
        # ~10 expanded clones (10-50 cells each)
        # ~50 medium clones (3-10 cells each)
        # ~remaining are singletons/doublets
        clone_sizes = []
        clone_cdr3s = []

        # Expanded clones
        for i in range(10):
            size = np.random.randint(10, 50)
            clone_sizes.extend([i] * size)
            clone_cdr3s.append((f"CAV{i:03d}QGNLIF", f"CASS{i:03d}YEQYF"))

        # Medium clones
        for i in range(10, 60):
            size = np.random.randint(3, 10)
            clone_sizes.extend([i] * size)
            clone_cdr3s.append((f"CAV{i:03d}SGTYKYIF", f"CASS{i:03d}DTQYF"))

        # Fill remaining with singletons/doublets
        next_clone = 60
        while len(clone_sizes) < n_cells:
            size = np.random.choice([1, 2], p=[0.7, 0.3])
            clone_sizes.extend([next_clone] * min(size, n_cells - len(clone_sizes)))
            clone_cdr3s.append((f"CAV{next_clone:03d}NTGNQFYF", f"CASS{next_clone:03d}NEQFF"))
            next_clone += 1

        clone_sizes = clone_sizes[:n_cells]
        np.random.shuffle(clone_sizes)

        # Build contigs (2 per cell - alpha and beta)
        rows = []
        for cell_idx in range(n_cells):
            clone_id = clone_sizes[cell_idx]
            cdr3_alpha, cdr3_beta = clone_cdr3s[clone_id]
            barcode = f"CELL{cell_idx:04d}-1"

            # Alpha chain
            rows.append(
                {
                    "barcode": barcode,
                    "contig_id": f"{barcode}_contig_1",
                    "chain": "TRA",
                    "v_gene": f"TRAV{np.random.randint(1, 30)}-1",
                    "j_gene": f"TRAJ{np.random.randint(1, 60)}",
                    "c_gene": "TRAC",
                    "d_gene": None,
                    "cdr3": cdr3_alpha,
                    "cdr3_nt": "N" * (len(cdr3_alpha) * 3),
                    "umis": np.random.randint(5, 200),
                    "reads": np.random.randint(50, 2000),
                    "productive": True,
                    "full_length": True,
                    "is_cell": True,
                    "high_confidence": True,
                    "raw_clonotype_id": f"clonotype{clone_id}",
                }
            )

            # Beta chain (95% have paired beta)
            if np.random.random() < 0.95:
                rows.append(
                    {
                        "barcode": barcode,
                        "contig_id": f"{barcode}_contig_2",
                        "chain": "TRB",
                        "v_gene": f"TRBV{np.random.randint(1, 30)}-1",
                        "j_gene": f"TRBJ{np.random.choice(['1-1', '1-2', '2-1', '2-7'])}",
                        "c_gene": np.random.choice(["TRBC1", "TRBC2"]),
                        "d_gene": np.random.choice(["TRBD1", "TRBD2"]),
                        "cdr3": cdr3_beta,
                        "cdr3_nt": "N" * (len(cdr3_beta) * 3),
                        "umis": np.random.randint(5, 200),
                        "reads": np.random.randint(50, 2000),
                        "productive": True,
                        "full_length": True,
                        "is_cell": True,
                        "high_confidence": True,
                        "raw_clonotype_id": f"clonotype{clone_id}",
                    }
                )

        return pd.DataFrame(rows)

    @pytest.fixture
    def synthetic_gex_data(self, synthetic_cellranger_data):
        """Generate matching GEX data with CD4/CD8 expression."""
        np.random.seed(42)
        barcodes = synthetic_cellranger_data["barcode"].unique()
        n_cells = len(barcodes)
        n_genes = 100

        # Random expression
        X = np.random.poisson(5, (n_cells, n_genes)).astype(float)

        # Gene names
        var_names = [f"Gene_{i}" for i in range(n_genes - 6)]
        var_names.extend(["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"])

        # 60% CD8+, 30% CD4+, 10% ambiguous
        cd8_cells = int(n_cells * 0.6)
        cd4_cells = int(n_cells * 0.3)

        # CD8+ cells
        X[:cd8_cells, -2] = np.random.poisson(25, cd8_cells)  # CD8A
        X[:cd8_cells, -1] = np.random.poisson(20, cd8_cells)  # CD8B
        X[:cd8_cells, -3] = np.random.poisson(2, cd8_cells)  # CD4

        # CD4+ cells
        X[cd8_cells : cd8_cells + cd4_cells, -3] = np.random.poisson(25, cd4_cells)  # CD4
        X[cd8_cells : cd8_cells + cd4_cells, -2] = np.random.poisson(2, cd4_cells)  # CD8A
        X[cd8_cells : cd8_cells + cd4_cells, -1] = np.random.poisson(2, cd4_cells)  # CD8B

        # CD3 for all
        X[:, -6:-3] = np.random.poisson(15, (n_cells, 3))

        adata = ad.AnnData(X)
        adata.var_names = var_names
        adata.obs_names = list(barcodes)
        adata.obs["sample"] = "test_sample"

        return adata

    @pytest.fixture
    def synthetic_til_data(self):
        """Generate TIL data with some overlapping clones."""
        np.random.seed(123)
        n_til_cells = 1000

        # Some clones overlap with culture, some are TIL-specific
        rows = []
        for i in range(n_til_cells):
            # 20% overlap with culture clones (using same CDR3 pattern)
            if np.random.random() < 0.2:
                clone_idx = np.random.randint(0, 20)  # Overlap with expanded culture clones
                cdr3_alpha = f"CAV{clone_idx:03d}QGNLIF"
                cdr3_beta = f"CASS{clone_idx:03d}YEQYF"
            else:
                # TIL-specific clones
                clone_idx = np.random.randint(500, 1000)
                cdr3_alpha = f"CAV{clone_idx:03d}TILSPEC"
                cdr3_beta = f"CASS{clone_idx:03d}TILSEQ"

            rows.append(
                {
                    "CDR3_alpha": cdr3_alpha,
                    "CDR3_beta": cdr3_beta,
                    "sample": "TIL_sample",
                }
            )

        df = pd.DataFrame(rows)
        df.index = [f"TIL_cell_{i}" for i in range(len(df))]
        return ad.AnnData(obs=df)

    @pytest.fixture
    def synthetic_database(self):
        """Generate a mock annotation database."""
        return pd.DataFrame(
            {
                "cdr3_beta": [
                    "CASS000YEQYF",  # Match first clone
                    "CASS001YEQYF",  # Match second clone
                    "CASSVIRALSEQ",  # Viral sequence
                ],
                "cdr3_alpha": [
                    "CAV000QGNLIF",
                    "CAV001QGNLIF",
                    "",
                ],
                "epitope": ["pp65", "MART-1", "EBV_BMLF1"],
                "species": ["CMV", "self", "EBV"],
                "database": ["VDJdb", "IEDB", "VDJdb"],
                "is_viral": [True, False, True],
            }
        )

    def test_load_aggregate_filter_pipeline(self, synthetic_cellranger_data, tmp_path):
        """Test loading VDJ data, aggregating to clonotypes, and filtering."""
        from tcrsift.clonotype import aggregate_clonotypes
        from tcrsift.filter import filter_clonotypes
        from tcrsift.loader import _pivot_vdj_by_barcode

        # Step 1: Pivot VDJ to cell-level
        vdj_df = synthetic_cellranger_data.copy()
        cells_df = _pivot_vdj_by_barcode(vdj_df)
        cells_df["sample"] = "test_sample"

        assert len(cells_df) > 0
        assert "CDR3_alpha" in cells_df.columns
        assert "CDR3_beta" in cells_df.columns

        # Step 2: Convert to AnnData and aggregate to clonotypes
        adata = ad.AnnData(obs=cells_df)
        clonotypes = aggregate_clonotypes(adata)

        assert len(clonotypes) > 0
        assert "CDR3ab" in clonotypes.columns
        assert "cell_count" in clonotypes.columns

        # Should have mix of clone sizes
        assert clonotypes["cell_count"].max() >= 10  # Some expanded
        assert clonotypes["cell_count"].min() <= 2  # Some rare

        # Step 3: Filter clonotypes (threshold method)
        filtered = filter_clonotypes(
            clonotypes,
            method="threshold",
            min_cells=2,
            tcell_type="both",
        )

        assert len(filtered) > 0
        assert "tier" in filtered.columns
        assert filtered["cell_count"].min() >= 2

    def test_logistic_filtering_pipeline(self, clonotypes_for_logistic):
        """Test logistic regression filtering with proper synthetic data."""
        from tcrsift.filter import filter_clonotypes, get_filter_summary

        # Filter with logistic method
        result = filter_clonotypes(
            clonotypes_for_logistic,
            method="logistic",
            min_cells=1,
            tcell_type="both",
        )

        assert "tier" in result.columns

        # Should have multiple tiers assigned
        tier_counts = result["tier"].value_counts()
        assert len(tier_counts) >= 2

        # Get summary
        summary = get_filter_summary(result)
        assert "total_clonotypes" in summary
        assert "tier_counts" in summary

    def test_til_matching_pipeline(self, synthetic_cellranger_data, synthetic_til_data, tmp_path):
        """Test TIL matching with culture clonotypes."""
        from tcrsift.clonotype import aggregate_clonotypes
        from tcrsift.loader import _pivot_vdj_by_barcode
        from tcrsift.til import get_til_summary, match_til

        # Load and aggregate culture data
        vdj_df = synthetic_cellranger_data.copy()
        cells_df = _pivot_vdj_by_barcode(vdj_df)
        cells_df["sample"] = "culture"
        adata = ad.AnnData(obs=cells_df)
        clonotypes = aggregate_clonotypes(adata)

        # Match against TIL
        matched = match_til(
            clonotypes,
            synthetic_til_data,
            match_by="CDR3ab",
            min_til_cells=1,
        )

        assert "til_match" in matched.columns
        assert "til_cell_count" in matched.columns

        # Should have some matches (we designed overlap)
        assert matched["til_match"].sum() > 0

        # Get summary
        summary = get_til_summary(matched)
        assert summary["til_matched_clones"] > 0
        assert summary["til_recovery_rate"] > 0

    def test_annotation_pipeline(self, synthetic_cellranger_data, synthetic_database, tmp_path):
        """Test database annotation workflow."""
        from tcrsift.annotate import match_clonotypes
        from tcrsift.clonotype import aggregate_clonotypes
        from tcrsift.loader import _pivot_vdj_by_barcode

        # Load and aggregate
        vdj_df = synthetic_cellranger_data.copy()
        cells_df = _pivot_vdj_by_barcode(vdj_df)
        cells_df["sample"] = "test"
        adata = ad.AnnData(obs=cells_df)
        clonotypes = aggregate_clonotypes(adata)

        # Annotate with synthetic database
        annotated = match_clonotypes(
            clonotypes,
            synthetic_database,
            match_by="CDR3ab",
        )

        assert "is_viral" in annotated.columns
        assert "db_match" in annotated.columns

    def test_full_pipeline_with_filtering_and_til(
        self, synthetic_cellranger_data, synthetic_til_data, clonotypes_for_logistic
    ):
        """Test complete pipeline: load -> aggregate -> filter -> TIL match."""
        from tcrsift.clonotype import aggregate_clonotypes
        from tcrsift.filter import filter_clonotypes, split_by_tier
        from tcrsift.loader import _pivot_vdj_by_barcode
        from tcrsift.til import get_til_summary, match_til

        # Load and aggregate
        vdj_df = synthetic_cellranger_data.copy()
        cells_df = _pivot_vdj_by_barcode(vdj_df)
        cells_df["sample"] = "culture"
        adata = ad.AnnData(obs=cells_df)
        clonotypes = aggregate_clonotypes(adata)

        # Filter with threshold method
        filtered = filter_clonotypes(
            clonotypes,
            method="threshold",
            min_cells=2,
            tcell_type="both",
        )

        # Split by tier
        tier_dfs = split_by_tier(filtered)
        assert len(tier_dfs) > 0

        # Match TIL for each tier
        for tier_name, tier_df in tier_dfs.items():
            matched = match_til(tier_df, synthetic_til_data, min_til_cells=1)
            summary = get_til_summary(matched)
            assert "til_matched_clones" in summary


class TestCLIEndToEnd:
    """Test CLI commands end-to-end."""

    def test_cli_help(self):
        """Test that CLI help works."""
        import subprocess

        result = subprocess.run(
            ["python", "-m", "tcrsift.cli", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "tcrsift" in result.stdout.lower() or "usage" in result.stdout.lower()

    def test_cli_version(self):
        """Test that CLI version works."""
        import subprocess

        result = subprocess.run(
            ["python", "-m", "tcrsift.cli", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        result = subprocess.run(
            ["python", "-m", "tcrsift.cli", "-v"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
