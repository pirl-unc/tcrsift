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

"""Tests for unified datacache coverage, mouse germline specs, offline pgen_models,
and offline pyensembl static Ensembl-ID fallback."""

from __future__ import annotations

import io
from unittest.mock import MagicMock

import pandas as pd
import pytest

from tcrsift import leaders, vregion
from tcrsift.datacache import (
    DATABASES,
    GROUPS,
    ensure_cached,
    resolve_db_selection,
)
from tcrsift.genes import ensembl_to_gene_name
from tcrsift.pgen_models import (
    cached_model_file,
    ensure_model,
    fetch_repertoire,
    train_model,
)
from tcrsift.til_select import load_iedb, load_vdjdb


class TestMouseGermlineCacheAndIsolation:
    """Mouse IMGT V-REGION and leader database specs + species isolation."""

    def test_germline_groups_preserve_human_and_expose_mouse(self):
        assert GROUPS["germline"] == ["imgt_trav_vregion", "imgt_trbv_vregion"]
        assert GROUPS["germline_mouse"] == [
            "imgt_trav_vregion_mouse",
            "imgt_trbv_vregion_mouse",
            "imgt_trav_leader_mouse",
            "imgt_trbv_leader_mouse",
        ]
        assert set(resolve_db_selection(None, ["germline_mouse"])) == set(
            GROUPS["germline_mouse"]
        )
        for name in GROUPS["germline_mouse"]:
            assert DATABASES[name].download_url is not None

    def test_mouse_leaders_and_vregions_do_not_leak_human_when_uncached(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("TCRSIFT_DATA_DIR", str(tmp_path))
        leaders.clear_reference_cache()
        vregion.clear_reference_cache()

        # Mouse returns empty/None when uncached rather than human TRAV12-2
        assert leaders.germline_vgene_leaders("TRAV12-2", species="mouse") == []
        assert leaders.germline_anchor_leader("METLLGLLILWLQLQWVSS", "TRAV12-2", species="mouse") is None
        assert vregion.reference_is_available(species="mouse") is False

        # Human bundled leaders remain available
        assert len(leaders.germline_vgene_leaders("TRAV12-2", species="human")) >= 1

    def test_mouse_leaders_and_vregions_load_from_cache(self, tmp_path, monkeypatch):
        monkeypatch.setenv("TCRSIFT_DATA_DIR", str(tmp_path))
        leaders.clear_reference_cache()
        vregion.clear_reference_cache()

        # Query BEFORE populating cache to verify empty {} is not permanently memoized
        assert leaders.germline_vgene_leaders("TRAV1", species="mouse") == []
        assert vregion.reference_is_available(species="mouse") is False

        # Populate mock IMGT HTML/FASTA files in cache using canonical spec filenames
        ldr_dir = tmp_path / "imgt_trav_leader_mouse"
        ldr_dir.mkdir(parents=True)
        (ldr_dir / DATABASES["imgt_trav_leader_mouse"].filename).write_text(
            "<html><body><pre>\n"
            ">J00660|TRAV1*01|Mus musculus|F|L-PART1+L-PART2|1..63|63 nt|1| | | | |63+0=63| | |\n"
            "atgaagaggctgctgtgttctctgctggggcttctgtgcacccaggtttgctgggtgaaagga\n"
            "</pre></body></html>\n"
        )

        vreg_dir = tmp_path / "imgt_trav_vregion_mouse"
        vreg_dir.mkdir(parents=True)
        (vreg_dir / DATABASES["imgt_trav_vregion_mouse"].filename).write_text(
            "<html><body><pre>\n"
            ">J00660|TRAV1*01|Mus musculus|F|V-REGION|1..282|282 nt|1| | | |94 AA|94+0=94| | |\n"
            "cagcagcaggtgagacagagcccacagagcctgaccgtgtgggagggcgagacagccatc"
            "ctgaactgcagctacgaggacagcaccttcgactacttcccctggtaccagcagttccca"
            "ggcgagggcccagccctgctgatcagcatcagaagcgtgagcgacaagaaggaggacggc"
            "agattcaccatcttcttcaacaagagagagaagaagctgtccctgcacatcgccgacagc"
            "cagccaggcgacagcgccacctacttctgcgccgccagc\n"
            "</pre></body></html>\n"
        )

        mouse_ldrs = leaders.germline_vgene_leaders("TRAV1", species="mouse")
        assert len(mouse_ldrs) == 1
        assert mouse_ldrs[0][2].startswith("M")
        assert vregion.reference_is_available(species="mouse") is True
        leaders.clear_reference_cache()
        vregion.clear_reference_cache()


class TestEnsureCachedAndTilSelectFormats:
    """`ensure_cached` auto-download behavior and `til_select` v3 format sniffing."""

    def test_ensure_cached_respects_auto_download_flag(self, tmp_path):
        assert ensure_cached("imgt_trav_leader_mouse", tmp_path, auto_download=False) is None

        fasta_payload = (
            b"<html><body><pre>\n"
            b">J00660|TRAV1*01|Mus musculus|F|L-PART1+L-PART2|1..63|63 nt|1|\n"
            b"atgaagaggctgctgtgttctctgctggggcttctgtgcacccaggtttgctgggtgaaagga\n"
            b"</pre></body></html>\n"
        )

        def fake_opener(req, timeout=None):
            resp = MagicMock()
            resp.read = io.BytesIO(fasta_payload).read
            resp.__enter__ = lambda s: s
            resp.__exit__ = lambda s, *exc: False
            return resp

        resolved = ensure_cached(
            "imgt_trav_leader_mouse",
            tmp_path,
            auto_download=True,
            opener=fake_opener,
        )
        assert resolved is not None and resolved.exists()

    def test_til_select_loads_vdjdb_full_and_iedb_csv_v3(self, tmp_path):
        vdjdb_file = tmp_path / "vdjdb_full.txt"
        pd.DataFrame([
            {
                "cdr3.alpha": "CAGAGSQGNLIF",
                "cdr3.beta": "CASSIRSSYEQYF",
                "species": "HomoSapiens",
                "antigen.epitope": "NLVPMVATV",
                "antigen.gene": "pp65",
                "antigen.species": "CMV",
                "mhc.a": "HLA-A*02:01",
                "vdjdb.score": 3,
            }
        ]).to_csv(vdjdb_file, sep="\t", index=False)
        vdf = load_vdjdb(vdjdb_file)
        assert len(vdf) == 1
        assert vdf.iloc[0]["cdr3_beta"] == "CASSIRSSYEQYF"

        # Comma-separated IEDB v3 export (stored by datacache as tcr_full_v3.csv)
        iedb_csv = tmp_path / "tcr_full_v3.csv"
        header_row_1 = [
            "Receptor", "Receptor", "Epitope", "Epitope", "Epitope", "Assay",
            "Chain 1", "Chain 1", "Chain 1", "Chain 2", "Chain 2", "Chain 2",
        ]
        header_row_2 = [
            "Receptor ID", "Type", "Name", "Source Molecule", "Source Organism", "MHC Allele Names",
            "Type", "CDR3 Curated", "CDR3 Calculated", "Type", "CDR3 Curated", "CDR3 Calculated",
        ]
        data_row = [
            "101", "alphabeta", "GILGFVFTL", "Matrix protein 1", "Influenza A virus", "HLA-A*02:01",
            "alpha", "CAGAGSQGNLIF", "", "beta", "CASSIRSSYEQYF", "",
        ]
        iedb_csv.write_text(
            ",".join(header_row_1) + "\n"
            + ",".join(header_row_2) + "\n"
            + ",".join(data_row) + "\n"
        )
        idf = load_iedb(iedb_csv)
        assert len(idf) == 1
        assert idf.iloc[0]["cdr3_beta"] == "CASSIRSSYEQYF"
        assert bool(idf.iloc[0]["is_viral"]) is True


class TestOfflinePgenModelsAndPyensembl:
    """Offline α/β repertoire training in pgen_models and offline pyensembl fallback."""

    def test_pgen_models_offline_alpha_and_kmer_gene_and_species_namespace(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("TCRSIFT_DATA_DIR", str(tmp_path))
        # Species namespacing in model cache paths
        human_path = cached_model_file("kmer", "ppost", "beta", species="human")
        mouse_path = cached_model_file("kmer", "ppost", "beta", species="mouse")
        assert human_path.name == "kmer_ppost_beta.npz"
        assert mouse_path.name == "kmer_ppost_mouse_beta.npz"

        # Alpha repertoire resolves offline from bundled observed_pbmc_10x.csv.gz
        rep_path = fetch_repertoire("alpha", species="human")
        rep = pd.read_csv(rep_path)
        assert len(rep) > 1000
        assert {"seq", "v", "j"} <= set(rep.columns)

        # kmer_gene training works offline from the extracted repertoire
        kg_path = train_model(
            "alpha",
            role="ppost",
            backend="kmer_gene",
            url=str(rep_path),
            species="human",
        )
        assert kg_path.exists()

        # ensure_model("alpha", role="ppost", backend="kmer") succeeds offline
        ensured = ensure_model(
            "alpha",
            role="ppost",
            backend="kmer",
            species="human",
        )
        assert ensured is not None

        # Mouse fetch_repertoire without local file raises clear offline guidance
        with pytest.raises(FileNotFoundError, match="mouse"):
            fetch_repertoire("beta", species="mouse")

    def test_offline_ensembl_to_gene_name_static_fallback(self):
        assert ensembl_to_gene_name("ENSG00000172116.8") == "CD8B"
        assert ensembl_to_gene_name("ENSG00000198888") == "MT-ND1"
        assert ensembl_to_gene_name("ENSMUSG00000053044") == "Cd8b1"
        assert ensembl_to_gene_name("ENSMUSG00000064351") == "mt-Co1"

    def test_load_10x_h5_marker_counts_resolves_mouse_cd8b1_and_ensembl(self, tmp_path):
        import h5py
        import numpy as np
        import scipy.sparse as sp

        from tcrsift.til_select import _load_10x_h5_marker_counts

        h5_path = tmp_path / "filtered_feature_bc_matrix.h5"
        # 2 cells x 3 features (Cd8b1, ENSMUSG00000015437=Gzmb, Cd4) -> CSC shape (3, 2)
        csc = sp.csc_matrix(np.array([[15, 0], [25, 5], [0, 20]], dtype=np.int32))
        with h5py.File(h5_path, "w") as f:
            grp = f.create_group("matrix")
            grp.create_dataset("barcodes", data=np.array([b"BC1-1", b"BC2-1"]))
            grp.create_dataset("data", data=csc.data)
            grp.create_dataset("indices", data=csc.indices)
            grp.create_dataset("indptr", data=csc.indptr)
            grp.create_dataset("shape", data=np.array(csc.shape, dtype=np.int64))
            feat = grp.create_group("features")
            feat.create_dataset(
                "name", data=np.array([b"Cd8b1", b"ENSMUSG00000015437", b"Cd4"])
            )
            feat.create_dataset(
                "id",
                data=np.array(
                    [b"ENSMUSG00000053044", b"ENSMUSG00000015437", b"ENSMUSG00000023274"]
                ),
            )
            feat.create_dataset(
                "feature_type",
                data=np.array([b"Gene Expression", b"Gene Expression", b"Gene Expression"]),
            )

        out = _load_10x_h5_marker_counts(h5_path, ["CD8B", "GZMB", "CD4"])
        assert list(out["count_CD8B"]) == [15.0, 0.0]
        assert list(out["count_GZMB"]) == [25.0, 5.0]
        assert list(out["count_CD4"]) == [0.0, 20.0]
