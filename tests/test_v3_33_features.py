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

"""Unit tests for v3.33.0 features:
- Static IMGT V-QUEST gapped FASTA parsing + datacache IPv4 fallback
- GENCODE v19 / Cell Ranger hg19 -> GRCh38 NeoTCR4/NeoTCR8/NeoTCR_PBL gene synonyms
- SampleSheet pre_sorted propagation + coreceptor-dropout rescue in phenotype_cells
- FOXP3+ Treg audit (`foxp3_pos_fraction`, `is_treg_like`, `treg_FOXP3`) & `--exclude-treg`
- Fast host_species_filter='same' pre-filtering in match_clonotypes
"""

from __future__ import annotations

import io
import socket
import urllib.error
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from tcrsift import datacache, vregion
from tcrsift.annotate import match_clonotypes
from tcrsift.cli import main as cli_main
from tcrsift.genes import GENE_SYNONYMS, harmonize_gene_symbols
from tcrsift.phenotype import phenotype_cells
from tcrsift.sample_sheet import Sample
from tcrsift.signature_methods import NEOTCR4_GENES_HGNC, score_by_name


class TestDatacacheIpv4FallbackAndGappedFasta:
    """Verify IPv4 fallback in _default_urlopen and gapped IMGT V-QUEST FASTA parsing."""

    def test_default_urlopen_falls_back_to_ipv4_on_403_or_timeout(self, monkeypatch):
        calls = []

        def fake_urlopen(req, timeout=20.0):
            families = {r[0] for r in socket.getaddrinfo("example.org", 443)}
            calls.append(families)
            if socket.AF_INET6 in families:
                raise urllib.error.HTTPError(
                    req.full_url, 403, "Forbidden", hdrs=None, fp=None
                )
            return io.BytesIO(b">ok\nACGT\n")

        def fake_getaddrinfo(host, port, *args, **kwargs):
            return [
                (socket.AF_INET6, socket.SOCK_STREAM, 6, "", ("::1", 443, 0, 0)),
                (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("127.0.0.1", 443)),
            ]

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        monkeypatch.setattr(socket, "getaddrinfo", fake_getaddrinfo)

        resp = datacache._default_urlopen("https://example.org/test.fasta", timeout=5.0)
        assert resp.read() == b">ok\nACGT\n"
        assert len(calls) == 2
        assert socket.AF_INET6 in calls[0]
        assert calls[1] == {socket.AF_INET}

    def test_default_urlopen_reraises_404_without_retry(self, monkeypatch):
        n_calls = 0

        def fake_urlopen(req, timeout=20.0):
            nonlocal n_calls
            n_calls += 1
            raise urllib.error.HTTPError(
                req.full_url, 404, "Not Found", hdrs=None, fp=None
            )

        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        with pytest.raises(urllib.error.HTTPError) as exc:
            datacache._default_urlopen("https://example.org/missing.fasta", timeout=5.0)
        assert exc.value.code == 404
        assert n_calls == 1

    def test_vregion_parses_gapped_imgt_vquest_fasta(self, tmp_path):
        vregion.clear_reference_cache()
        vdir = tmp_path / "imgt_trbv_vregion_mouse"
        vdir.mkdir(parents=True)
        # Include IMGT '.' alignment gaps inside V-REGION
        (vdir / datacache.DATABASES["imgt_trbv_vregion_mouse"].filename).write_text(
            ">X01643|TRBV13-2*01|Mus musculus|F|V-REGION|1..285|285 nt|1| | | | |285+27=312| | |\n"
            "gaggctgcagtcacccaaagcccaagaaacaaggtggcagtaacaggaggaaaggtgaca"
            "ttgagctgt...aatcagact...aataaccacaacaacatgtactggtatcggcaggac"
            "acggggcatggg...ctgaggctgatccattactcatatggtgctggcagcactgagaaa"
            "ggagatatccctgatggatacaaggcctccagaccaagccaagagaacttctccctcatt"
            "ctggagttggctaccccctctcagacatcagtgtacttctgtgccagcggtgat\n"
        )
        seqs = vregion.germline_vregion_sequences(
            "Trbv13-2*01", species="mouse", data_dir=tmp_path
        )
        assert len(seqs) == 1
        allele, func, aa, nt = seqs[0]
        assert allele == "01" and func == "F"
        assert "." not in nt and "*" not in aa and aa.endswith("CASGD")
        vregion.clear_reference_cache()

    def test_cellranger_mouse_hyphen_dv_genes_normalize_to_imgt_slash(self):
        from tcrsift.genes import canonicalize_gene
        from tcrsift.leaders import normalize_vgene

        assert normalize_vgene("TRAV13-4-dv7") == "TRAV13-4/DV7"
        assert normalize_vgene("TRAV6-7-dv9") == "TRAV6-7/DV9"
        assert normalize_vgene("TRAV16d-dv11") == "TRAV16D/DV11"
        assert normalize_vgene("TRAV15d-1-dv6d-1") == "TRAV15D-1/DV6D-1"
        assert normalize_vgene("TRAV15-1-dv6-1") == "TRAV15-1/DV6-1"
        assert canonicalize_gene("TRAV13-4-dv7") == canonicalize_gene("TRAV13-4/DV7")
        assert canonicalize_gene("TRAV15d-1-dv6d-1") == canonicalize_gene("TRAV15D-1/DV6D-1")


class TestHg19GeneSynonymsAndPreSortedPhenotype:
    """Verify hg19 Cell Ranger gene synonyms and FACS pre_sorted coreceptor-dropout rescue."""

    def test_neotcr4_scores_on_hg19_gpr56_and_rp11_symbols(self):
        renamed = {"ADGRG1": "GPR56", "LINC01943": "RP11-291B21.2"}
        genes = [renamed.get(g, g) for g in NEOTCR4_GENES_HGNC] + [
            f"BG_{i}" for i in range(300)
        ]
        rng = np.random.default_rng(7)
        X = sp.csr_matrix(np.log1p(rng.poisson(3, size=(40, len(genes)))))
        adata = ad.AnnData(X=X, var=pd.DataFrame(index=genes))
        scores = score_by_name(adata, "NeoTCR4", log1p=False, on_missing="error")
        assert len(scores) == 40 and np.isfinite(scores).all()

    def test_gene_synonyms_Remain_unique_and_ensembl_keyed(self):
        tokens = [tok for grp in GENE_SYNONYMS for tok in grp]
        assert len(tokens) == len(set(tokens))
        assert harmonize_gene_symbols(["GPR56", "PPAP2A"], ["ADGRG1", "PLPP1"]).tolist() == [
            "ADGRG1",
            "PLPP1",
        ]

    def test_sample_pre_sorted_normalizes_case_and_rescues_cd4_dropout(self):
        s = Sample(sample="S_cd4", vdj_dir="/tmp/v", gex_dir="/tmp/g", pre_sorted="cd4")
        assert s.pre_sorted == "CD4"

        obs = pd.DataFrame(
            {
                "CD3": [10, 10, 10, 1, 10],
                "CD4": [5, 0, 0, 0, 0],
                "CD8": [0, 0, 6, 0, 0],
                "pre_sorted": ["CD4", "CD4", "CD4", "CD4", pd.NA],
            },
            index=[f"c{i}" for i in range(5)],
        )
        adata = ad.AnnData(X=np.zeros((5, 1)), obs=obs)
        out = phenotype_cells(adata, min_cd3_reads=3, verbose=False)
        # c0: CD4=5 -> Confident CD4+
        # c1: CD4=0, CD8=0, CD3=10, pre_sorted=CD4 -> rescued Likely CD4+
        # c2: CD4=0, CD8=6, pre_sorted=CD4 -> stays Confident CD8+ (sorting impurity rejected)
        # c3: CD4=0, CD8=0, CD3=1 (<3), pre_sorted=CD4 -> stays Unknown
        # c4: CD4=0, CD8=0, CD3=10, pre_sorted=NA -> stays Unknown
        assert out.obs["Tcell_type"].tolist() == [
            "Confident CD4+",
            "Likely CD4+",
            "Confident CD8+",
            "Unknown",
            "Unknown",
        ]
        assert out.obs["is_CD4"].tolist() == [True, True, False, False, False]

    def test_rehydrate_obs_preserves_mixed_pre_sorted_categorical(self):
        from tcrsift._dtypes import rehydrate_obs

        obs = pd.DataFrame(
            {"pre_sorted": pd.Categorical(["CD4", None]), "species": pd.Categorical(["mouse", "mouse"])},
            index=["c0", "c1"],
        )
        adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
        rehydrate_obs(adata)
        assert adata.obs["pre_sorted"].tolist() == ["CD4", None]
        assert adata.obs["species"].tolist() == ["mouse", "mouse"]


class TestFoxp3TregAuditAndExclusion:
    """Verify FOXP3/Foxp3 Treg fraction, risk_flags, and --exclude-treg in prioritize."""

    def test_treg_audit_and_opt_in_exclusion_human_and_mouse(self, tmp_path, monkeypatch):
        from tcrsift import loader
        from tcrsift.signature_methods import SIGNATURES

        genes = sorted({g for sig in SIGNATURES.values() for g in sig.all_genes}) + [
            "FOXP3",
            "FOXP3",  # duplicate var_name regression check
            *[f"BG_{i}" for i in range(200)],
        ]
        rng = np.random.default_rng(11)
        X = rng.poisson(2, size=(14, len(genes))).astype(float)
        # Clone 1 (4 cells): effector + FOXP3-high Treg-like (3/4 = 0.75 FOXP3+)
        # Clone 2 (4 cells): effector + FOXP3-negative (0/4 = 0.0 FOXP3+)
        # Clone 3 (6 cells): quiescent background so both Clone 1 and Clone 2 have z > 0
        for g in ("PRF1", "GZMB", "IFNG", "TNFRSF9"):
            X[:8, genes.index(g)] = 60.0
        X[:, genes.index("FOXP3")] = [8, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        X[:, genes.index("FOXP3") + 1] = 0.0

        obs = pd.DataFrame(
            {
                "sample": ["s1"] * 14,
                "patient_id": ["p1"] * 14,
                "CDR3_alpha": ["CAVTREGQGNLIF"] * 4 + ["CAVEFFGQGNLIF"] * 4 + ["CAVAGQGNLIF"] * 6,
                "CDR3_beta": ["CASSTREGYEQYF"] * 4 + ["CASSEFFGYEQYF"] * 4 + ["CASSAGYEQYF"] * 6,
                "TRA_1_umis": 5,
                "TRB_1_umis": 5,
                "TRA_1_reads": 50,
                "TRB_1_reads": 50,
                "CD3": 20,
                "CD4": 15,
                "CD8": 0,
            },
            index=[f"cell_{i}" for i in range(14)],
        )
        cells = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
        sheet = tmp_path / "samples.yaml"
        sheet.write_text("samples:\n  - sample: s1\n    vdj_dir: v\n    gex_dir: g\n")
        monkeypatch.setattr(loader, "load_samples", lambda *a, **kw: cells.copy())

        # Default (--no-exclude-treg): Treg clone is audited in risk_flags but retained
        out_keep = tmp_path / "keep_treg"
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_keep),
                "--tcell-type",
                "cd4",
                "--signatures",
                "Cytolytic",
                "--signature-cutoff",
                "legacy",
                "--signature-quantile",
                "0",
            ]
        )
        audit_keep = pd.read_csv(out_keep / "all_scored_clones.csv").set_index("CDR3_alpha")
        assert audit_keep.loc["CAVTREGQGNLIF", "foxp3_pos_fraction"] == pytest.approx(0.75)
        assert bool(audit_keep.loc["CAVTREGQGNLIF", "is_treg_like"]) is True
        assert "treg_FOXP3" in str(audit_keep.loc["CAVTREGQGNLIF", "risk_flags"])
        assert bool(audit_keep.loc["CAVTREGQGNLIF", "selected_for_review"]) is True
        assert audit_keep.loc["CAVEFFGQGNLIF", "foxp3_pos_fraction"] == pytest.approx(0.0)
        assert bool(audit_keep.loc["CAVEFFGQGNLIF", "is_treg_like"]) is False

        # Opt-in --exclude-treg: Treg clone is excluded with reason treg_FOXP3, effector survives
        out_drop = tmp_path / "drop_treg"
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_drop),
                "--tcell-type",
                "cd4",
                "--signatures",
                "Cytolytic",
                "--signature-cutoff",
                "legacy",
                "--signature-quantile",
                "0",
                "--exclude-treg",
            ]
        )
        audit_drop = pd.read_csv(out_drop / "all_scored_clones.csv").set_index("CDR3_alpha")
        cands_drop = pd.read_csv(out_drop / "candidate_clones.csv")
        assert audit_drop.loc["CAVTREGQGNLIF", "excluded_reason"] == "treg_FOXP3"
        assert cands_drop["CDR3_alpha"].tolist() == ["CAVEFFGQGNLIF"]

    @pytest.mark.parametrize("bad_val", ["0", "-0.2", "1.1", "nan"])
    def test_invalid_treg_foxp3_fraction_fails_before_loading(self, tmp_path, bad_val):
        with pytest.raises(SystemExit) as exc:
            cli_main(
                [
                    "prioritize",
                    "missing.yaml",
                    "-o",
                    str(tmp_path),
                    "--treg-foxp3-fraction",
                    bad_val,
                ]
            )
        assert exc.value.code == 1


class TestHostSpeciesPreFilter:
    """Verify host_species_filter='same' isolates mouse vs human database records."""

    def test_mouse_query_matches_only_mus_musculus_records(self):
        clones = pd.DataFrame(
            {
                "CDR3_alpha": ["CAASDTNTGKLTF", "CAVSDGGSQGNLIF"],
                "CDR3_beta": ["CASSDDRVNERLFF", "CASSIRSSYEQYF"],
                "CDR3ab": [
                    "CAASDTNTGKLTF_CASSDDRVNERLFF",
                    "CAVSDGGSQGNLIF_CASSIRSSYEQYF",
                ],
            }
        )
        db = pd.DataFrame(
            {
                "cdr3_alpha": ["CAASDTNTGKLTF", "CAVSDGGSQGNLIF"],
                "cdr3_beta": ["CASSDDRVNERLFF", "CASSIRSSYEQYF"],
                "epitope": ["KYNKANVFL", "GILGFVFTL"],
                "antigen_gene": ["IGRP", "M1"],
                "species": ["Mus musculus", "Influenza A"],
                "host_species": ["MusMusculus", "HomoSapiens"],
                "is_viral": [False, True],
                "database": ["VDJdb", "VDJdb"],
            }
        )
        out = match_clonotypes(
            clones,
            db,
            match_strictness="strict_ab",
            species="mouse",
            host_species_filter="same",
            verbose=False,
            show_progress=False,
        )
        assert out["db_match"].tolist() == [True, False]
        assert out.loc[0, "db_epitope"] == "KYNKANVFL"
        assert out.loc[0, "db_host_species"] == "Mus musculus"
