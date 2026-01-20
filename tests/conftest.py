"""
Pytest fixtures for TCRsift tests.

Fixtures use realistic data patterns based on CellRanger VDJ output:
- CDR3 sequences follow IMGT conventions (alpha starts "CA", beta starts "CASS")
- VDJ segment columns match CellRanger filtered_contig_annotations.csv format
- UMI/read counts reflect typical single-cell sequencing ranges
"""

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

# Realistic VDJ segment sequences from IMGT reference
# These are truncated examples based on typical TCR sequences
REALISTIC_TRA_FWR1 = "MSLGLLCCVALSLLNAGTS"
REALISTIC_TRA_CDR1 = "TSINN"
REALISTIC_TRA_FWR2 = "WYRQQKPGQPPE"
REALISTIC_TRA_CDR2 = "LLSN"
REALISTIC_TRA_FWR3 = "NKGQRQPQRLNASLDKSGV"
REALISTIC_TRA_CDR3 = "CAVSDGGSQGNLIF"  # Complete CDR3 alpha
REALISTIC_TRA_FWR4 = "FGKGTHLIIQP"

REALISTIC_TRB_FWR1 = "MGSRLLCWVLLCLLGAGPVKA"
REALISTIC_TRB_CDR1 = "SGHATL"
REALISTIC_TRB_FWR2 = "WYRQQKAGGLQP"
REALISTIC_TRB_CDR2 = "FNNFY"
REALISTIC_TRB_FWR3 = "ISNKLSLNVSDGSQLFLNIK"
REALISTIC_TRB_CDR3 = "CASSLGQAYEQYF"  # Complete CDR3 beta
REALISTIC_TRB_FWR4 = "FGPGTRLLVL"

# Nucleotide sequences (simplified - actual sequences are longer)
REALISTIC_TRA_CDR3_NT = "TGTGCTGTGTCAGATGGAGGAAGCCAGGGAAATCTCATCTTT"
REALISTIC_TRB_CDR3_NT = "TGTGCCAGCAGTTTGGGACAGGCTTACGAGCAGTACTTC"


@pytest.fixture
def sample_vdj_df():
    """Create a sample VDJ DataFrame for testing.

    Matches CellRanger filtered_contig_annotations.csv format with:
    - Realistic CDR3 sequences following IMGT convention
    - Realistic UMI/read counts (35-185 UMIs, 350-1850 reads)
    - CellRanger-like contig_id format

    Note: Does not include VDJ segment columns (fwr1, cdr1, etc.) to avoid
    conflicts with the loader's pivot logic. Use sample_vdj_df_with_segments
    for tests requiring full segment data.
    """
    return pd.DataFrame(
        {
            "barcode": ["AAAA", "AAAA", "BBBB", "BBBB", "CCCC", "CCCC", "DDDD", "DDDD"],
            "contig_id": [
                "AAAA_contig_1",
                "AAAA_contig_2",
                "BBBB_contig_1",
                "BBBB_contig_2",
                "CCCC_contig_1",
                "CCCC_contig_2",
                "DDDD_contig_1",
                "DDDD_contig_2",
            ],
            "chain": ["TRA", "TRB", "TRA", "TRB", "TRA", "TRB", "TRA", "TRB"],
            "cdr3": [
                "CAVSDGGSQGNLIF",
                "CASSLGQAYEQYF",  # Clone A
                "CAVSDGGSQGNLIF",
                "CASSLGQAYEQYF",  # Clone A (same)
                "CAVSAGGSQGNLIF",
                "CASSLGQAYEQYF",  # Clone B (different alpha)
                "CAVNAGGSQGNLIF",
                "CASSLAGAYEQYF",  # Clone C (different clone)
            ],
            "v_gene": [
                "TRAV12-2",
                "TRBV28",
                "TRAV12-2",
                "TRBV28",
                "TRAV12-1",
                "TRBV28",
                "TRAV12-3",
                "TRBV7-9",
            ],
            "d_gene": [None, "TRBD1", None, "TRBD1", None, "TRBD2", None, "TRBD1"],
            "j_gene": [
                "TRAJ33",
                "TRBJ2-7",
                "TRAJ33",
                "TRBJ2-7",
                "TRAJ33",
                "TRBJ2-7",
                "TRAJ42",
                "TRBJ2-1",
            ],
            "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1", "TRAC", "TRBC2", "TRAC", "TRBC1"],
            "umis": [85, 120, 42, 78, 35, 52, 145, 185],
            "reads": [850, 1200, 420, 780, 350, 520, 1450, 1850],
            "sample": ["S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2"],
            "productive": [True, True, True, True, True, True, True, True],
            "full_length": [True, True, True, True, True, True, True, True],
            "raw_clonotype_id": [
                "clonotype1",
                "clonotype1",
                "clonotype1",
                "clonotype1",
                "clonotype2",
                "clonotype2",
                "clonotype3",
                "clonotype3",
            ],
        }
    )


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
    X[:30, -3] = np.random.poisson(1, 30)  # CD4

    # CD4+ cells: high CD4, low CD8
    X[30:60, -3] = np.random.poisson(20, 30)  # CD4
    X[30:60, -2] = np.random.poisson(1, 30)  # CD8A
    X[30:60, -1] = np.random.poisson(1, 30)  # CD8B

    # Unknown cells: similar levels
    X[60:, -3] = np.random.poisson(5, 40)  # CD4
    X[60:, -2] = np.random.poisson(5, 40)  # CD8A
    X[60:, -1] = np.random.poisson(5, 40)  # CD8B

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
    """Create a sample clonotypes DataFrame for testing.

    Represents aggregated clonotype data as output from aggregate_clonotypes().
    Includes realistic:
    - CDR3 sequences following IMGT conventions
    - V/J/C gene annotations matching IMGT nomenclature
    - Cell counts and sample distributions
    - T cell type consensus annotations
    """
    return pd.DataFrame(
        {
            "clone_id": [
                "CAVSDGGSQGNLIF_CASSLGQAYEQYF",
                "CAVSAGGSQGNLIF_CASSLGQAYEQYF",
                "CAVNAGGSQGNLIF_CASSLAGAYEQYF",
                "CAVTGGSNTGKLIF_CASSFSGANVLTF",
                "CAVRSGGSYLTF_CASSLGGETQYF",
            ],
            "CDR3_alpha": [
                "CAVSDGGSQGNLIF",
                "CAVSAGGSQGNLIF",
                "CAVNAGGSQGNLIF",
                "CAVTGGSNTGKLIF",
                "CAVRSGGSYLTF",
            ],
            "CDR3_beta": [
                "CASSLGQAYEQYF",
                "CASSLGQAYEQYF",
                "CASSLAGAYEQYF",
                "CASSFSGANVLTF",
                "CASSLGGETQYF",
            ],
            "cell_count": [15, 8, 3, 2, 1],
            "samples": ["S1;S2", "S1", "S2", "S1", "S2"],
            "n_samples": [2, 1, 1, 1, 1],
            "max_frequency": [0.15, 0.08, 0.03, 0.02, 0.01],
            "Tcell_type_consensus": [
                "Confident CD8+",
                "Confident CD8+",
                "Confident CD4+",
                "Likely CD8+",
                "Unknown",
            ],
            "alpha_v_gene": ["TRAV12-2", "TRAV12-1", "TRAV12-3", "TRAV21", "TRAV8-6"],
            "alpha_j_gene": ["TRAJ33", "TRAJ33", "TRAJ42", "TRAJ37", "TRAJ52"],
            "alpha_c_gene": ["TRAC", "TRAC", "TRAC", "TRAC", "TRAC"],
            "beta_v_gene": ["TRBV28", "TRBV28", "TRBV7-9", "TRBV5-1", "TRBV28"],
            "beta_j_gene": ["TRBJ2-7", "TRBJ2-7", "TRBJ2-1", "TRBJ2-6", "TRBJ2-5"],
            "beta_c_gene": ["TRBC1", "TRBC1", "TRBC1", "TRBC2", "TRBC1"],
        }
    )


@pytest.fixture
def sample_vdj_df_with_segments():
    """Create VDJ DataFrame with full VDJ segment columns.

    Mimics CellRanger output with all IMGT segments:
    - fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4 (amino acid)
    - fwr1_nt, cdr1_nt, etc. (nucleotide)
    """
    return pd.DataFrame(
        {
            "barcode": ["AAACCTGAGAACTCGG-1", "AAACCTGAGAACTCGG-1"],
            "contig_id": ["AAACCTGAGAACTCGG-1_contig_1", "AAACCTGAGAACTCGG-1_contig_2"],
            "chain": ["TRA", "TRB"],
            "cdr3": [REALISTIC_TRA_CDR3, REALISTIC_TRB_CDR3],
            "cdr3_nt": [REALISTIC_TRA_CDR3_NT, REALISTIC_TRB_CDR3_NT],
            "v_gene": ["TRAV12-2", "TRBV28"],
            "d_gene": [None, "TRBD1"],
            "j_gene": ["TRAJ33", "TRBJ2-7"],
            "c_gene": ["TRAC", "TRBC1"],
            "fwr1": [REALISTIC_TRA_FWR1, REALISTIC_TRB_FWR1],
            "cdr1": [REALISTIC_TRA_CDR1, REALISTIC_TRB_CDR1],
            "fwr2": [REALISTIC_TRA_FWR2, REALISTIC_TRB_FWR2],
            "cdr2": [REALISTIC_TRA_CDR2, REALISTIC_TRB_CDR2],
            "fwr3": [REALISTIC_TRA_FWR3, REALISTIC_TRB_FWR3],
            "fwr4": [REALISTIC_TRA_FWR4, REALISTIC_TRB_FWR4],
            "umis": [95, 142],
            "reads": [950, 1420],
            "sample": ["S1", "S1"],
            "productive": [True, True],
            "full_length": [True, True],
            "raw_clonotype_id": ["clonotype1", "clonotype1"],
        }
    )


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_csv_sample_sheet(temp_dir):
    """Create a sample CSV sample sheet."""
    csv_path = temp_dir / "samples.csv"
    df = pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "vdj_dir": ["/fake/path/S1/vdj", "/fake/path/S2/vdj"],
            "gex_dir": ["/fake/path/S1/gex", "/fake/path/S2/gex"],
            "antigen_type": ["short_peptide", "long_peptide"],
            "antigen_description": ["CMV pp65", "KRAS G12D"],
            "source": ["culture", "culture"],
        }
    )
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
    return pd.DataFrame(
        {
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
        }
    )


@pytest.fixture
def sample_full_length_clonotypes():
    """Create clonotypes DataFrame with full-length TCR sequences.

    Based on realistic TCR sequences as they would appear after
    full-length assembly. Includes:
    - VDJ sequences (without leader/constant)
    - Full sequences (leader + VDJ + constant)
    - Single-chain constructs (beta-T2A-alpha)

    Lengths are realistic:
    - VDJ region: ~100-120 AA
    - Full chain: ~250-280 AA
    """
    # Realistic VDJ sequence (~110 AA each)
    vdj_alpha = (
        REALISTIC_TRA_FWR1
        + REALISTIC_TRA_CDR1
        + REALISTIC_TRA_FWR2
        + REALISTIC_TRA_CDR2
        + REALISTIC_TRA_FWR3
        + REALISTIC_TRA_CDR3
        + REALISTIC_TRA_FWR4
    )
    vdj_beta = (
        REALISTIC_TRB_FWR1
        + REALISTIC_TRB_CDR1
        + REALISTIC_TRB_FWR2
        + REALISTIC_TRB_CDR2
        + REALISTIC_TRB_FWR3
        + REALISTIC_TRB_CDR3
        + REALISTIC_TRB_FWR4
    )

    # Leader sequences (typically 20-25 AA)
    leader_alpha = "MKSLRVLLVILWLQLDS"
    leader_beta = "MGSRLLCCVVLCLLGAGPVKA"

    # Constant region endings (truncated for testing)
    constant_alpha = "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS"
    constant_beta = "DLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"

    return pd.DataFrame(
        {
            "clone_id": ["clone1", "clone2"],
            "CDR3_alpha": [REALISTIC_TRA_CDR3, "CAVSAGGSQGNLIF"],
            "CDR3_beta": [REALISTIC_TRB_CDR3, "CASSLAGAYEQYF"],
            "VDJ_alpha_aa": [vdj_alpha, vdj_alpha.replace(REALISTIC_TRA_CDR3, "CAVSAGGSQGNLIF")],
            "VDJ_beta_aa": [vdj_beta, vdj_beta.replace(REALISTIC_TRB_CDR3, "CASSLAGAYEQYF")],
            "alpha_leader_aa": [leader_alpha, leader_alpha],
            "beta_leader_aa": [leader_beta, leader_beta],
            "alpha_constant_aa": [constant_alpha, constant_alpha],
            "beta_constant_aa": [constant_beta, constant_beta],
            "full_alpha_aa": [
                leader_alpha + vdj_alpha + constant_alpha,
                leader_alpha
                + vdj_alpha.replace(REALISTIC_TRA_CDR3, "CAVSAGGSQGNLIF")
                + constant_alpha,
            ],
            "full_beta_aa": [
                leader_beta + vdj_beta + constant_beta,
                leader_beta + vdj_beta.replace(REALISTIC_TRB_CDR3, "CASSLAGAYEQYF") + constant_beta,
            ],
            "alpha_c_gene": ["TRAC", "TRAC"],
            "beta_c_gene": ["TRBC1", "TRBC1"],
            "alpha_v_gene": ["TRAV12-2", "TRAV12-1"],
            "beta_v_gene": ["TRBV28", "TRBV7-9"],
            "samples": ["S1", "S2"],
            "cell_count": [10, 5],
        }
    )


@pytest.fixture
def clonotypes_for_logistic():
    """Create clonotypes suitable for logistic regression model fitting.

    Generates 200 clones with realistic variation to avoid perfect separation:
    - Varied max_frequency spanning 0.001 to 0.5
    - Mix of viral (20%) and non-viral (80%) clones
    - Varied n_samples (1-5)
    - Realistic cell counts

    The data is designed so the logistic model can fit without warnings:
    - Some high-frequency viral clones (model sees high freq != always good)
    - Some low-frequency non-viral clones (noise in target)
    """
    np.random.seed(42)
    n = 200

    # Generate frequencies with realistic distribution (log-normal-ish)
    # Mix of rare and expanded clones
    frequencies = np.concatenate(
        [
            np.random.beta(1, 20, n // 2),  # Many rare clones (0-0.1)
            np.random.beta(2, 5, n // 4),  # Some expanded (0.1-0.3)
            np.random.beta(5, 3, n // 4),  # Some highly expanded (0.3-0.6)
        ]
    )
    np.random.shuffle(frequencies)
    frequencies = frequencies[:n]

    # Viral status: 20% viral, but NOT perfectly correlated with frequency
    # Some high-freq clones are viral, some low-freq are non-viral
    is_viral = np.random.random(n) < 0.20

    # Cell counts correlated with frequency but with noise
    cell_counts = (frequencies * 500 + np.random.poisson(3, n)).astype(int)
    cell_counts = np.maximum(cell_counts, 1)

    # n_samples: higher frequency clones tend to appear in more samples
    n_samples = np.ones(n, dtype=int)
    n_samples[frequencies > 0.1] = np.random.choice([1, 2], size=(frequencies > 0.1).sum())
    n_samples[frequencies > 0.3] = np.random.choice([2, 3, 4], size=(frequencies > 0.3).sum())

    return pd.DataFrame(
        {
            "clone_id": [f"clone_{i}" for i in range(n)],
            "CDR3_alpha": [f"CAV{i:03d}QGNLIF" for i in range(n)],
            "CDR3_beta": [f"CASS{i:03d}YEQYF" for i in range(n)],
            "cell_count": cell_counts,
            "max_frequency": frequencies,
            "n_samples": n_samples,
            "is_viral": is_viral,
            "Tcell_type_consensus": ["Confident CD8+"] * (n // 2) + ["Confident CD4+"] * (n // 2),
        }
    )


@pytest.fixture
def mock_cellranger_vdj_dir(temp_dir):
    """Create a mock CellRanger VDJ output directory.

    Structure matches CellRanger vdj output:
    - filtered_contig_annotations.csv
    - clonotypes.csv
    - all_contig.fasta
    """
    vdj_dir = temp_dir / "vdj_output"
    vdj_dir.mkdir()

    # Create filtered_contig_annotations.csv
    annotations_df = pd.DataFrame(
        {
            "barcode": [
                "AAACCTGAGAACTCGG-1",
                "AAACCTGAGAACTCGG-1",
                "AAACCTGCATCGGGTC-1",
                "AAACCTGCATCGGGTC-1",
            ],
            "is_cell": [True, True, True, True],
            "contig_id": [
                "AAACCTGAGAACTCGG-1_contig_1",
                "AAACCTGAGAACTCGG-1_contig_2",
                "AAACCTGCATCGGGTC-1_contig_1",
                "AAACCTGCATCGGGTC-1_contig_2",
            ],
            "high_confidence": [True, True, True, True],
            "length": [556, 714, 548, 702],
            "chain": ["TRA", "TRB", "TRA", "TRB"],
            "v_gene": ["TRAV12-2", "TRBV28", "TRAV12-1", "TRBV7-9"],
            "d_gene": [None, "TRBD1", None, "TRBD2"],
            "j_gene": ["TRAJ33", "TRBJ2-7", "TRAJ42", "TRBJ2-1"],
            "c_gene": ["TRAC", "TRBC1", "TRAC", "TRBC1"],
            "full_length": [True, True, True, True],
            "productive": [True, True, True, True],
            "cdr3": [REALISTIC_TRA_CDR3, REALISTIC_TRB_CDR3, "CAVSAGGSQGNLIF", "CASSLAGAYEQYF"],
            "cdr3_nt": [
                REALISTIC_TRA_CDR3_NT,
                REALISTIC_TRB_CDR3_NT,
                "TGTGCTGTGTCAGCTGGAGGAAGCCAGGGAAATCTCATCTTT",
                "TGTGCCAGCAGTTTGGCTGGAGCTTACGAGCAGTACTTC",
            ],
            "reads": [850, 1200, 420, 780],
            "umis": [85, 120, 42, 78],
            "raw_clonotype_id": ["clonotype1", "clonotype1", "clonotype2", "clonotype2"],
        }
    )
    annotations_df.to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

    # Create clonotypes.csv
    clonotypes_df = pd.DataFrame(
        {
            "clonotype_id": ["clonotype1", "clonotype2"],
            "frequency": [2, 2],
            "proportion": [0.5, 0.5],
            "cdr3s_aa": [
                f"TRA:{REALISTIC_TRA_CDR3};TRB:{REALISTIC_TRB_CDR3}",
                "TRA:CAVSAGGSQGNLIF;TRB:CASSLAGAYEQYF",
            ],
            "cdr3s_nt": [
                f"TRA:{REALISTIC_TRA_CDR3_NT};TRB:{REALISTIC_TRB_CDR3_NT}",
                "TRA:TGTGCTGTGTCAGCTGGAGGAAGCCAGGGAAATCTCATCTTT;TRB:TGTGCCAGCAGTTTGGCTGGAGCTTACGAGCAGTACTTC",
            ],
        }
    )
    clonotypes_df.to_csv(vdj_dir / "clonotypes.csv", index=False)

    # Create all_contig.fasta
    fasta_content = f""">AAACCTGAGAACTCGG-1_contig_1
ATGGCTAAGCTGAGGATGCTGAGCCTGCTGCTGTGTCTGCTGGCACTGCTGACCAGCGCAGGCTCC{REALISTIC_TRA_CDR3_NT}TTCGGCAAGGGA
>AAACCTGAGAACTCGG-1_contig_2
ATGGGCTCCAGGCTCCTCCTGTGCGTGCTGCTGCTCTGCCTGGGGGCAGGCTCC{REALISTIC_TRB_CDR3_NT}TTCGGACCGGGCACC
>AAACCTGCATCGGGTC-1_contig_1
ATGGCTAAGCTGAGGATGCTGAGCCTGCTGCTGTGTCTGCTGGCACTGCTGACCAGCGCAGGCTCCTGTGCTGTGTCAGCTGGAGGAAGCCAGGGAAATCTCATCTTTTTCGGCAAGGGA
>AAACCTGCATCGGGTC-1_contig_2
ATGGGCTCCAGGCTCCTCCTGTGCGTGCTGCTGCTCTGCCTGGGGGCAGGCTCCTGTGCCAGCAGTTTGGCTGGAGCTTACGAGCAGTACTTCTTCGGACCGGGCACC
"""
    (vdj_dir / "all_contig.fasta").write_text(fasta_content)

    return vdj_dir
