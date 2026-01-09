"""
Pytest fixtures for TCRsift tests.
"""

import pytest
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import tempfile
import os


@pytest.fixture
def sample_vdj_df():
    """Create a sample VDJ DataFrame for testing."""
    return pd.DataFrame({
        "barcode": ["AAAA", "AAAA", "BBBB", "BBBB", "CCCC", "CCCC", "DDDD", "DDDD"],
        "contig_id": ["AAAA_1", "AAAA_2", "BBBB_1", "BBBB_2", "CCCC_1", "CCCC_2", "DDDD_1", "DDDD_2"],
        "chain": ["TRA", "TRB", "TRA", "TRB", "TRA", "TRB", "TRA", "TRB"],
        "cdr3": [
            "CAVSDGGSQGNLIF", "CASSLGQAYEQYF",
            "CAVSDGGSQGNLIF", "CASSLGQAYEQYF",  # Same clone as AAAA
            "CAVSAGGSQGNLIF", "CASSLGQAYEQYF",  # Different alpha
            "CAVNAGGSQGNLIF", "CASSLAGAYEQYF",  # Different clone
        ],
        "v_gene": ["TRAV1", "TRBV1", "TRAV1", "TRBV1", "TRAV2", "TRBV1", "TRAV3", "TRBV2"],
        "j_gene": ["TRAJ1", "TRBJ1", "TRAJ1", "TRBJ1", "TRAJ1", "TRBJ1", "TRAJ2", "TRBJ2"],
        "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1", "TRAC", "TRBC2", "TRAC", "TRBC1"],
        "umis": [10, 15, 8, 12, 5, 8, 20, 25],
        "reads": [100, 150, 80, 120, 50, 80, 200, 250],
        "sample": ["S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2"],
    })


@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for testing."""
    n_cells = 100
    n_genes = 50

    # Random expression matrix
    X = np.random.poisson(5, (n_cells, n_genes)).astype(float)

    # Gene names including T cell markers
    var_names = [f"Gene_{i}" for i in range(n_genes - 6)]
    var_names.extend(["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"])

    # Set some cells to be CD8+ and some CD4+
    # CD8+ cells: high CD8A/CD8B, low CD4
    X[:30, -2] = np.random.poisson(20, 30)  # CD8A
    X[:30, -1] = np.random.poisson(15, 30)  # CD8B
    X[:30, -3] = np.random.poisson(1, 30)   # CD4

    # CD4+ cells: high CD4, low CD8
    X[30:60, -3] = np.random.poisson(20, 30)  # CD4
    X[30:60, -2] = np.random.poisson(1, 30)   # CD8A
    X[30:60, -1] = np.random.poisson(1, 30)   # CD8B

    # Unknown cells: similar levels
    X[60:, -3] = np.random.poisson(5, 40)   # CD4
    X[60:, -2] = np.random.poisson(5, 40)   # CD8A
    X[60:, -1] = np.random.poisson(5, 40)   # CD8B

    # CD3 for all T cells
    X[:, -6] = np.random.poisson(10, n_cells)  # CD3D
    X[:, -5] = np.random.poisson(10, n_cells)  # CD3E
    X[:, -4] = np.random.poisson(10, n_cells)  # CD3G

    adata = ad.AnnData(X)
    adata.var_names = var_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
    adata.obs["sample"] = ["S1"] * 50 + ["S2"] * 50

    # Add some TCR info
    adata.obs["CDR3_alpha"] = ["CAVSDGGSQGNLIF"] * 20 + ["CAVSAGGSQGNLIF"] * 30 + [None] * 50
    adata.obs["CDR3_beta"] = ["CASSLGQAYEQYF"] * 20 + ["CASSLAGAYEQYF"] * 30 + [None] * 50
    adata.obs["has_TRA"] = [True] * 50 + [False] * 50
    adata.obs["has_TRB"] = [True] * 50 + [False] * 50
    adata.obs["has_both_chains"] = [True] * 50 + [False] * 50

    return adata


@pytest.fixture
def sample_clonotypes_df():
    """Create a sample clonotypes DataFrame for testing."""
    return pd.DataFrame({
        "clone_id": [
            "CAVSDGGSQGNLIF_CASSLGQAYEQYF",
            "CAVSAGGSQGNLIF_CASSLGQAYEQYF",
            "CAVNAGGSQGNLIF_CASSLAGAYEQYF",
            "CAVXAGGSQGNLIF_CASSXAGAYEQYF",
            "CAVYAGGSQGNLIF_CASSYAGAYEQYF",
        ],
        "CDR3_alpha": [
            "CAVSDGGSQGNLIF", "CAVSAGGSQGNLIF", "CAVNAGGSQGNLIF",
            "CAVXAGGSQGNLIF", "CAVYAGGSQGNLIF",
        ],
        "CDR3_beta": [
            "CASSLGQAYEQYF", "CASSLGQAYEQYF", "CASSLAGAYEQYF",
            "CASSXAGAYEQYF", "CASSYAGAYEQYF",
        ],
        "cell_count": [15, 8, 3, 2, 1],
        "samples": ["S1;S2", "S1", "S2", "S1", "S2"],
        "n_samples": [2, 1, 1, 1, 1],
        "max_frequency": [0.15, 0.08, 0.03, 0.02, 0.01],
        "Tcell_type_consensus": [
            "Confident CD8+", "Confident CD8+", "Confident CD4+",
            "Likely CD8+", "Unknown",
        ],
        "alpha_v_gene": ["TRAV1", "TRAV2", "TRAV3", "TRAV4", "TRAV5"],
        "beta_v_gene": ["TRBV1", "TRBV1", "TRBV2", "TRBV3", "TRBV4"],
    })


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_csv_sample_sheet(temp_dir):
    """Create a sample CSV sample sheet."""
    csv_path = temp_dir / "samples.csv"
    df = pd.DataFrame({
        "sample": ["S1", "S2"],
        "vdj_dir": ["/fake/path/S1/vdj", "/fake/path/S2/vdj"],
        "gex_dir": ["/fake/path/S1/gex", "/fake/path/S2/gex"],
        "antigen_type": ["short_peptide", "long_peptide"],
        "antigen_description": ["CMV pp65", "KRAS G12D"],
        "source": ["culture", "culture"],
    })
    df.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def sample_yaml_sample_sheet(temp_dir):
    """Create a sample YAML sample sheet."""
    yaml_path = temp_dir / "samples.yaml"
    yaml_content = """
samples:
  - sample: "S1"
    vdj_dir: "/fake/path/S1/vdj"
    gex_dir: "/fake/path/S1/gex"
    antigen_type: "short_peptide"
    antigen_description: "CMV pp65"
    source: "culture"
  - sample: "S2"
    vdj_dir: "/fake/path/S2/vdj"
    gex_dir: "/fake/path/S2/gex"
    antigen_type: "long_peptide"
    antigen_description: "KRAS G12D"
    source: "culture"
"""
    yaml_path.write_text(yaml_content)
    return yaml_path


@pytest.fixture
def sample_database_df():
    """Create a sample database DataFrame for testing annotation."""
    return pd.DataFrame({
        "cdr3_beta": [
            "CASSLGQAYEQYF",  # Matches our test clone
            "CASSXYZAYEQYF",  # Different
            "CASSLGQAYEQYF",  # Duplicate with different epitope
        ],
        "cdr3_alpha": [
            "CAVSDGGSQGNLIF",
            "CAVXYZQGNLIF",
            "CAVSDGGSQGNLIF",
        ],
        "epitope": ["NLV", "GLC", "pp65"],
        "species": ["CMV", "EBV", "CMV"],
        "database": ["VDJdb", "VDJdb", "IEDB"],
        "is_viral": [True, True, True],
    })
