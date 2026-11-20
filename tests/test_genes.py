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

"""Tests for genes module."""

import numpy as np
import pytest

from tcrsift.genes import (
    CD3D,
    TCELL_MARKERS,
    Gene,
    build_versioned_rename_map,
    find_column_for_gene,
    find_gene,
    find_gene_by_ensembl,
    find_gene_by_symbol,
    strip_ensembl_version,
)


class TestGene:
    """Tests for Gene dataclass."""

    def test_gene_creation(self):
        """Test creating a Gene."""
        gene = Gene("TEST", "ENSG00000000001", "Test gene")
        assert gene.symbol == "TEST"
        assert gene.ensembl_id == "ENSG00000000001"
        assert gene.description == "Test gene"

    def test_matches_ensembl_exact(self):
        """Test exact ENSEMBL ID matching."""
        assert CD3D.matches_ensembl("ENSG00000167286")

    def test_matches_ensembl_with_version(self):
        """Test matching versioned ENSEMBL ID."""
        assert CD3D.matches_ensembl("ENSG00000167286.10")
        assert CD3D.matches_ensembl("ENSG00000167286.15")

    def test_matches_ensembl_wrong_id(self):
        """Test non-matching ENSEMBL ID."""
        assert not CD3D.matches_ensembl("ENSG00000010610")  # This is CD4
        assert not CD3D.matches_ensembl("ENSG00000010610.5")

    def test_gene_is_frozen(self):
        """Test that Gene is immutable."""
        with pytest.raises(AttributeError):
            CD3D.symbol = "CHANGED"


class TestStripEnsemblVersion:
    """Tests for strip_ensembl_version function."""

    def test_strips_version(self):
        """Test stripping version suffix."""
        assert strip_ensembl_version("ENSG00000167286.10") == "ENSG00000167286"

    def test_no_version(self):
        """Test ID without version."""
        assert strip_ensembl_version("ENSG00000167286") == "ENSG00000167286"

    def test_multiple_dots(self):
        """Test handling multiple dots (should only split on first)."""
        # This shouldn't happen in practice, but test the behavior
        assert strip_ensembl_version("ENSG00000167286.10.extra") == "ENSG00000167286"


class TestFindGeneBySymbol:
    """Tests for find_gene_by_symbol function."""

    def test_find_cd3d(self):
        """Test finding CD3D by symbol."""
        gene = find_gene_by_symbol("CD3D")
        assert gene is not None
        assert gene.symbol == "CD3D"

    def test_find_case_insensitive(self):
        """Test case-insensitive symbol matching."""
        gene = find_gene_by_symbol("cd3d")
        assert gene is not None
        assert gene.symbol == "CD3D"

    def test_not_found(self):
        """Test gene not found."""
        gene = find_gene_by_symbol("NOTREAL")
        assert gene is None


class TestFindGeneByEnsembl:
    """Tests for find_gene_by_ensembl function."""

    def test_find_by_base_id(self):
        """Test finding by base ENSEMBL ID."""
        gene = find_gene_by_ensembl("ENSG00000167286")
        assert gene is not None
        assert gene.symbol == "CD3D"

    def test_find_by_versioned_id(self):
        """Test finding by versioned ENSEMBL ID."""
        gene = find_gene_by_ensembl("ENSG00000167286.10")
        assert gene is not None
        assert gene.symbol == "CD3D"

    def test_find_different_version(self):
        """Test finding with different version number."""
        gene = find_gene_by_ensembl("ENSG00000167286.99")
        assert gene is not None
        assert gene.symbol == "CD3D"

    def test_not_found(self):
        """Test ENSEMBL ID not found."""
        gene = find_gene_by_ensembl("ENSG99999999999")
        assert gene is None


class TestFindGene:
    """Tests for find_gene function (combined lookup)."""

    def test_find_by_symbol(self):
        """Test finding by symbol."""
        gene = find_gene("CD4")
        assert gene is not None
        assert gene.symbol == "CD4"

    def test_find_by_ensembl(self):
        """Test finding by ENSEMBL ID."""
        gene = find_gene("ENSG00000010610")
        assert gene is not None
        assert gene.symbol == "CD4"

    def test_find_by_versioned_ensembl(self):
        """Test finding by versioned ENSEMBL ID."""
        gene = find_gene("ENSG00000010610.5")
        assert gene is not None
        assert gene.symbol == "CD4"


class TestFindColumnForGene:
    """Tests for find_column_for_gene function."""

    def test_find_by_symbol(self):
        """Test finding column by gene symbol."""
        columns = ["other", "CD3D", "more"]
        col = find_column_for_gene(CD3D, columns)
        assert col == "CD3D"

    def test_find_by_symbol_case_insensitive(self):
        """Test finding column by symbol (case-insensitive)."""
        columns = ["other", "cd3d", "more"]
        col = find_column_for_gene(CD3D, columns)
        assert col == "cd3d"

    def test_find_by_ensembl_base(self):
        """Test finding column by base ENSEMBL ID."""
        columns = ["other", "ENSG00000167286", "more"]
        col = find_column_for_gene(CD3D, columns)
        assert col == "ENSG00000167286"

    def test_find_by_ensembl_versioned(self):
        """Test finding column by versioned ENSEMBL ID."""
        columns = ["other", "ENSG00000167286.10", "more"]
        col = find_column_for_gene(CD3D, columns)
        assert col == "ENSG00000167286.10"

    def test_prefers_symbol_over_ensembl(self):
        """Test that symbol match is preferred over ENSEMBL."""
        columns = ["ENSG00000167286.10", "CD3D", "other"]
        col = find_column_for_gene(CD3D, columns)
        assert col == "CD3D"

    def test_not_found(self):
        """Test column not found."""
        columns = ["other", "stuff", "more"]
        col = find_column_for_gene(CD3D, columns)
        assert col is None


class TestBuildVersionedRenameMap:
    """Tests for build_versioned_rename_map function."""

    def test_renames_versioned_ids(self):
        """Test building rename map for versioned ENSEMBL IDs."""
        columns = ["ENSG00000167286.10", "ENSG00000010610.5", "other_col"]
        rename_map = build_versioned_rename_map(columns)

        assert rename_map == {
            "ENSG00000167286.10": "CD3D",
            "ENSG00000010610.5": "CD4",
        }

    def test_skips_non_ensembl_columns(self):
        """Test that non-ENSEMBL columns are skipped."""
        columns = ["CD3D", "gene_name", "other"]
        rename_map = build_versioned_rename_map(columns)
        assert rename_map == {}

    def test_skips_unknown_ensembl_ids(self):
        """Test that unknown ENSEMBL IDs are skipped."""
        columns = ["ENSG99999999999.1", "ENSG00000167286.10"]
        rename_map = build_versioned_rename_map(columns)
        assert rename_map == {"ENSG00000167286.10": "CD3D"}

    def test_empty_columns(self):
        """Test with empty column list."""
        rename_map = build_versioned_rename_map([])
        assert rename_map == {}


class TestTcellMarkers:
    """Tests for TCELL_MARKERS constant."""

    def test_contains_expected_markers(self):
        """Test that TCELL_MARKERS contains all expected genes."""
        symbols = {g.symbol for g in TCELL_MARKERS}
        assert symbols == {"CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"}

    def test_all_have_ensembl_ids(self):
        """Test that all markers have ENSEMBL IDs."""
        for gene in TCELL_MARKERS:
            assert gene.ensembl_id.startswith("ENSG")


class TestCanonicalizeGene:
    """One shared IMGT canonicalizer for all V/J matching (the 'real gap')."""

    def test_allele_stripped_by_default(self):
        from tcrsift.genes import canonicalize_gene
        assert canonicalize_gene("TRBV20-1*01") == "TRBV20-1"

    def test_keep_allele_normalizes_to_two_digits(self):
        from tcrsift.genes import canonicalize_gene
        assert canonicalize_gene("TCRBV20-1*1", keep_allele=True) == "TRBV20-1*01"

    def test_adaptive_prefix_and_leading_zeros(self):
        from tcrsift.genes import canonicalize_gene
        assert canonicalize_gene("TCRBV20-01*01") == "TRBV20-1"
        assert canonicalize_gene("TRBV09") == "TRBV9"
        assert canonicalize_gene("TRAJ08") == "TRAJ8"

    def test_dv_slash_collapses(self):
        from tcrsift.genes import canonicalize_gene
        # the TRAV14/DV4 ↔ TRAV14DV4 mismatch the spec flagged
        assert canonicalize_gene("TRAV14/DV4") == canonicalize_gene("TRAV14DV4")
        assert canonicalize_gene("TRAV38-2/DV8") == "TRAV38-2DV8"

    def test_case_and_whitespace(self):
        from tcrsift.genes import canonicalize_gene
        assert canonicalize_gene(" trbv20-1 ") == "TRBV20-1"

    def test_empty_and_none(self):
        from tcrsift.genes import canonicalize_gene
        assert canonicalize_gene("") is None
        assert canonicalize_gene(None) is None


class TestHarmonizeGeneSymbols:
    """Signature genes are matched across reference-annotation renames."""

    def test_newer_symbol_is_matched_to_published_name(self):
        from tcrsift.genes import harmonize_gene_symbols

        out = harmonize_gene_symbols(np.array(["CD8A", "PLAAT3", "CARS1"]), ["PLA2G16", "CARS", "CD8A"])
        assert list(out) == ["CD8A", "PLA2G16", "CARS"]

    def test_older_symbol_is_matched_to_newer_name(self):
        from tcrsift.genes import harmonize_gene_symbols

        assert list(harmonize_gene_symbols(np.array(["FAM129A"]), ["NIBAN1"])) == ["NIBAN1"]

    def test_gene_without_symbol_is_matched_by_ensembl_id(self):
        from tcrsift.genes import harmonize_gene_symbols

        out = harmonize_gene_symbols(np.array(["ENSG00000267364"]), ["AC022706.1"])
        assert list(out) == ["AC022706.1"]

    def test_present_and_unknown_genes_are_untouched(self):
        from tcrsift.genes import harmonize_gene_symbols

        symbols = np.array(["PLA2G16", "PLAAT3", "GZMB"])
        out = harmonize_gene_symbols(symbols, ["PLA2G16", "NOT_A_GENE"])
        assert list(out) == ["PLA2G16", "PLAAT3", "GZMB"]
        assert out is not symbols

    def test_every_synonym_names_one_gene(self):
        from tcrsift.genes import GENE_SYNONYMS

        names = [name for group in GENE_SYNONYMS for name in group]
        assert len(names) == len(set(names))
        assert all(group[0].startswith("ENSG") and len(group) >= 2 for group in GENE_SYNONYMS)


def test_neotcr8_scores_on_reference_with_renamed_genes():
    """NeoTCR8 lists PLA2G16/CARS; a 2020-A or later matrix has PLAAT3/CARS1."""
    import anndata as ad
    import scipy.sparse as sp

    from tcrsift.signature_methods import NEOTCR8_GENES_HGNC, score_by_name

    renamed = {"PLA2G16": "PLAAT3", "CARS": "CARS1", "AC022706.1": "ENSG00000267364"}
    genes = [renamed.get(g, g) for g in NEOTCR8_GENES_HGNC] + [f"FILL{i}" for i in range(400)]
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=sp.csr_matrix(np.log1p(rng.poisson(rng.integers(1, 30, len(genes)), (60, len(genes))))))
    adata.var_names = genes

    scores = score_by_name(adata, "NeoTCR8", log1p=False, on_missing="error")

    assert len(scores) == 60 and np.isfinite(scores).all()
