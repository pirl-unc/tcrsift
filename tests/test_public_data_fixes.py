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

"""Regression and feature tests for public-dataset compatibility fixes (v3.30.0)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.io import mmwrite

from tcrsift.annotate import annotate_clonotypes, load_vdjdb
from tcrsift.assemble import assemble_full_sequences, load_contigs, parse_fasta
from tcrsift.cli import main as tcrsift_main
from tcrsift.clonotype import aggregate_clonotypes
from tcrsift.datacache import cached_path
from tcrsift.filter import score_stimulation_response
from tcrsift.loader import (
    combine_gex_and_vdj,
    load_cellranger_gex,
    load_cellranger_vdj,
    load_sample,
)
from tcrsift.multimer import (
    call_multimer_specificity,
    classify_feature_columns,
    extract_feature_barcodes,
    summarize_clone_multimer_specificity,
)
from tcrsift.sample_sheet import Sample, SampleSheet
from tcrsift.validation import (
    TCRsiftValidationError,
    resolve_cellranger_gex_source,
    resolve_cellranger_vdj_file,
)


def _write_mtx_dir(
    matrix_dir: Path,
    X_cells_by_feats: np.ndarray,
    barcodes: list[str],
    features: list[tuple[str, str, str]],
) -> None:
    matrix_dir.mkdir(parents=True, exist_ok=True)
    mtx_tmp = matrix_dir / "matrix.mtx"
    mmwrite(mtx_tmp, sp.coo_matrix(X_cells_by_feats.T))
    with open(mtx_tmp, "rb") as fin, gzip.open(matrix_dir / "matrix.mtx.gz", "wb") as fout:
        fout.write(fin.read())
    mtx_tmp.unlink()
    with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt", encoding="utf-8") as fb:
        for bc in barcodes:
            fb.write(f"{bc}\n")
    with gzip.open(matrix_dir / "features.tsv.gz", "wt", encoding="utf-8") as ff:
        for fid, fname, ftype in features:
            ff.write(f"{fid}\t{fname}\t{ftype}\n")


class TestVdjdbCacheAndPairedResolution:
    """Ensure cached VDJdb resolves vdjdb_full.txt and long-format recovers complex.id > 0 pairs."""

    def test_cached_path_prefers_directory_when_vdjdb_txt_and_full_coexist(self, tmp_path):
        vdir = tmp_path / "vdjdb"
        vdir.mkdir()
        # Real VDJdb release zips contain BOTH vdjdb.txt and vdjdb_full.txt.
        pd.DataFrame(
            {
                "complex.id": [0],
                "gene": ["TRB"],
                "cdr3": ["CASSIRSSYEQYF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
            }
        ).to_csv(vdir / "vdjdb.txt", sep="\t", index=False)
        pd.DataFrame(
            {
                "cdr3.alpha": ["CAAGAGSQGNLIF"],
                "cdr3.beta": ["CASSIRSSYEQYF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
                "species": ["HomoSapiens"],
            }
        ).to_csv(vdir / "vdjdb_full.txt", sep="\t", index=False)

        resolved = cached_path("vdjdb", tmp_path)
        assert resolved == vdir

        loaded = load_vdjdb(resolved, verbose=False)
        assert len(loaded) == 1
        assert loaded.iloc[0]["cdr3_alpha"] == "CAAGAGSQGNLIF"
        assert loaded.iloc[0]["cdr3_beta"] == "CASSIRSSYEQYF"

        clones = pd.DataFrame(
            {"CDR3_alpha": ["CAAGAGSQGNLIF"], "CDR3_beta": ["CASSIRSSYEQYF"]}
        )
        matched = annotate_clonotypes(
            clones,
            vdjdb_path=resolved,
            match_strictness="strict_ab",
            verbose=False,
            show_progress=False,
        )
        assert bool(matched.iloc[0]["db_match"]) is True
        assert matched.iloc[0]["db_match_strength"] == "ab"
        assert matched.iloc[0]["db_epitope"] == "GILGFVFTL"

    def test_long_format_recovers_paired_alpha_via_positive_complex_id(self, tmp_path):
        vfile = tmp_path / "vdjdb.txt"
        pd.DataFrame(
            {
                "complex.id": [101, 101, 0, 0],
                "gene": ["TRA", "TRB", "TRA", "TRB"],
                "cdr3": [
                    "CAAGAGSQGNLIF",
                    "CASSIRSSYEQYF",
                    "CAVUNPAIREDALPHA",
                    "CASSSINGLEBETAF",
                ],
                "v.segm": ["TRAV27*01", "TRBV19*01", "TRAV12-2*01", "TRBV7-9*01"],
                "antigen.epitope": ["GILGFVFTL", "GILGFVFTL", "ELAGIGILTV", "NLVPMVATV"],
                "antigen.species": ["InfluenzaA", "InfluenzaA", "HomoSapiens", "CMV"],
                "species": ["HomoSapiens"] * 4,
            }
        ).to_csv(vfile, sep="\t", index=False)

        loaded = load_vdjdb(vfile, verbose=False)
        assert len(loaded) == 2
        paired = loaded[loaded["cdr3_beta"] == "CASSIRSSYEQYF"].iloc[0]
        assert paired["cdr3_alpha"] == "CAAGAGSQGNLIF"
        assert paired["v_alpha"] == "TRAV27*01"
        unpaired = loaded[loaded["cdr3_beta"] == "CASSSINGLEBETAF"].iloc[0]
        assert pd.isna(unpaired["cdr3_alpha"])
        assert "CAVUNPAIREDALPHA" not in set(loaded["cdr3_beta"].dropna())

    def test_pick_vdjdb_dated_filename_fallback(self, tmp_path):
        vdir = tmp_path / "vdjdb"
        vdir.mkdir()
        pd.DataFrame(
            {
                "complex.id": [0],
                "gene": ["TRB"],
                "cdr3": ["CASSIRSSYEQYF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
            }
        ).to_csv(vdir / "vdjdb-2024-08-15.txt", sep="\t", index=False)
        loaded = load_vdjdb(vdir, verbose=False)
        assert len(loaded) == 1
        assert loaded.iloc[0]["cdr3_beta"] == "CASSIRSSYEQYF"

    def test_cli_annotate_database_match_and_levenshtein(self, tmp_path):
        vfile = tmp_path / "vdjdb_full.txt"
        pd.DataFrame(
            {
                "cdr3.alpha": ["CAAGAGSQGNLIF"],
                "cdr3.beta": ["CASSIRSSYEQYF"],
                "antigen.epitope": ["GILGFVFTL"],
                "antigen.species": ["InfluenzaA"],
                "species": ["HomoSapiens"],
            }
        ).to_csv(vfile, sep="\t", index=False)
        clones_csv = tmp_path / "clones.csv"
        pd.DataFrame(
            {
                "CDR3ab": ["CAVWRONG_CASSIRSSYEQYF", "CAVOTHER_CASSIRSSYEQFF"],
                "CDR3_alpha": ["CAVWRONG", "CAVOTHER"],
                "CDR3_beta": ["CASSIRSSYEQYF", "CASSIRSSYEQFF"],
            }
        ).to_csv(clones_csv, index=False)

        out_strict = tmp_path / "strict.csv"
        tcrsift_main(
            [
                "annotate",
                "-i",
                str(clones_csv),
                "-o",
                str(out_strict),
                "--vdjdb",
                str(vfile),
                "--database-match",
                "strict_ab",
                "--flag-only",
            ]
        )
        df_strict = pd.read_csv(out_strict)
        assert df_strict["db_match"].sum() == 0

        out_lev = tmp_path / "lev.csv"
        tcrsift_main(
            [
                "annotate",
                "-i",
                str(clones_csv),
                "-o",
                str(out_lev),
                "--vdjdb",
                str(vfile),
                "--database-match",
                "b_only",
                "--match-mode",
                "levenshtein",
                "--max-distance",
                "1",
                "--flag-only",
            ]
        )
        df_lev = pd.read_csv(out_lev)
        assert df_lev["db_match"].tolist() == [True, True]
        assert df_lev["db_match_distance"].tolist() == [0, 1]


class TestFlexibleCellRangerLoadingAndJsonEnrichment:
    """Test gzipped/prefixed GEO & 10x inputs, gex_path, raw matrix, and JSON contig assembly."""

    def test_prefixed_gzipped_vdj_and_partial_5p_rescue(self, tmp_path):
        sdir = tmp_path / "GSM5216976_sample"
        sdir.mkdir()
        df = pd.DataFrame(
            {
                "barcode": ["BC1-1", "BC1-1", "BC2-1", "BC3-1"],
                "is_cell": ["True", "True", "True", "True"],
                "contig_id": ["BC1-1_c1", "BC1-1_c2", "BC2-1_c1", "BC3-1_c1"],
                "high_confidence": ["True", "True", "True", "True"],
                "chain": ["TRA", "TRB", "TRB", "TRA"],
                "v_gene": ["TRAV27", "TRBV19", "TRBV6-6", "TRAV12-2"],
                "d_gene": ["None", "TRBD1", "None", "None"],
                "j_gene": ["TRAJ42", "TRBJ2-7", "None", "TRAJ20"],
                "c_gene": ["TRAC", "TRBC2", "TRBC1", "TRAC"],
                "productive": ["None", "True", "None", "False"],
                "cdr3": [
                    "CAAGAGSQGNLIF",       # 5'-partial (productive=None) with complete V+J+CDR3 -> rescue
                    "CASSIRSSYEQYF",       # productive=True -> keep
                    "CVPVDCLAAPSPTPRFLW",  # j_gene=None -> must drop
                    "CAVR*STOP*KLVF",      # productive=False -> must drop
                ],
                "umis": [12, 18, 5, 9],
                "reads": [120, 180, 50, 90],
            }
        )
        with gzip.open(
            sdir / "GSM5216976_Sample_filtered_contig_annotations.csv.gz", "wt", encoding="utf-8"
        ) as fh:
            df.to_csv(fh, index=False)

        resolved, unfiltered = resolve_cellranger_vdj_file(sdir)
        assert not unfiltered
        assert resolved.name.endswith(".csv.gz")

        loaded = load_cellranger_vdj(sdir, "GT6_stim", verbose=False)
        assert sorted(loaded["contig_id"].tolist()) == ["BC1-1_c1", "BC1-1_c2"]

    def test_ambiguous_prefixed_vdj_files_raise(self, tmp_path):
        (tmp_path / "A_filtered_contig_annotations.csv").write_text("barcode,chain\n")
        (tmp_path / "B_filtered_contig_annotations.csv").write_text("barcode,chain\n")
        with pytest.raises(TCRsiftValidationError, match="Multiple VDJ"):
            resolve_cellranger_vdj_file(tmp_path)

    def test_raw_feature_bc_matrix_and_gex_path_loading(self, tmp_path, monkeypatch):
        import scanpy as sc

        # 1. raw_feature_bc_matrix fallback
        gdir = tmp_path / "raw_only_outs"
        raw_mtx = gdir / "raw_feature_bc_matrix"
        X = np.array([[10, 5, 4, 0, 8, 6], [8, 4, 3, 7, 0, 0]], dtype=float)
        feats = [
            ("ENSG00000167286", "CD3D", "Gene Expression"),
            ("ENSG00000198851", "CD3E", "Gene Expression"),
            ("ENSG00000160654", "CD3G", "Gene Expression"),
            ("ENSG00000010610", "CD4", "Gene Expression"),
            ("ENSG00000153563", "CD8A", "Gene Expression"),
            ("ENSG00000172116", "CD8B", "Gene Expression"),
        ]
        _write_mtx_dir(raw_mtx, X, ["BC1-1", "BC2-1"], feats)
        kind, src = resolve_cellranger_gex_source(gdir)
        assert kind == "mtx" and src == raw_mtx
        adata_raw = load_cellranger_gex(
            gdir, "S_raw", min_genes=1, min_counts=1, min_mito_pct=0.0, verbose=False
        )
        assert adata_raw.n_obs == 2

        # 2. Sample.gex_path with sample-prefixed .h5 file
        h5_file = tmp_path / "donor1_filtered_feature_bc_matrix.h5"
        h5_file.write_bytes(b"placeholder")
        mock_h5_adata = ad.AnnData(
            X=sp.csr_matrix(X),
            obs=pd.DataFrame(index=["BC1-1", "BC2-1"]),
            var=pd.DataFrame(
                {"gene_ids": [f[0] for f in feats]},
                index=[f[1] for f in feats],
            ),
        )
        monkeypatch.setattr(sc, "read_10x_h5", lambda *_a, **_kw: mock_h5_adata.copy())

        vdj_dir = tmp_path / "vdj_s1"
        vdj_dir.mkdir()
        pd.DataFrame(
            {
                "barcode": ["BC1-1", "BC1-1"],
                "chain": ["TRA", "TRB"],
                "v_gene": ["TRAV27", "TRBV19"],
                "d_gene": [None, "TRBD1"],
                "j_gene": ["TRAJ42", "TRBJ2-7"],
                "c_gene": ["TRAC", "TRBC2"],
                "productive": [True, True],
                "cdr3": ["CAAGAGSQGNLIF", "CASSIRSSYEQYF"],
                "umis": [10, 15],
                "reads": [100, 150],
                "contig_id": ["c1", "c2"],
            }
        ).to_csv(vdj_dir / "filtered_contig_annotations.csv", index=False)

        sample = Sample(
            sample="S_h5",
            gex_path=str(h5_file),
            vdj_dir=str(vdj_dir),
            source="control",
        )
        assert sample.is_control()
        sheet = SampleSheet(samples=[sample])
        assert len(sheet.get_control_samples()) == 1
        loaded = load_sample(sample, min_genes=1, min_counts=1, min_mito_pct=0.0)
        assert loaded is not None and loaded.n_obs == 2
        assert loaded.obs.loc["BC1-1", "CDR3ab"] == "CAAGAGSQGNLIF_CASSIRSSYEQYF"
        assert loaded.obs["source"].iloc[0] == "control"

    def test_json_vdj_enrichment_and_contig_loading_for_assembly(self, tmp_path):
        """Pre-CR6 CSV without fwr1..fwr4 + all_contig_annotations.json.gz enables full assembly."""
        from tcrsift.assemble import HUMAN_CONSTANT_REGIONS_AA, back_translate
        from tcrsift.leaders import germline_vgene_leaders

        sdir = tmp_path / "GT6_stim"
        sdir.mkdir()

        # Construct realistic full-length alpha (TRAV27 + TRAJ42 + TRAC) and beta (TRBV5-1 + TRBJ2-1 + TRBC2) contigs
        alpha_leader_aa = germline_vgene_leaders("TRAV27", species="human")[0][2]
        beta_leader_aa = germline_vgene_leaders("TRBV5-1", species="human")[0][2]
        alpha_vdj_aa = (
            "TQLLEQSPQFLSIQEGENLTVYCNSSSVFSSLQWYRQEPGEGPVLLVTVVTGGEVKKLKRLTFQFGDARKDSSLHITAAQPGDTGLYLC"
            "CAAGAGSQGNLIFGKGTKLSVKP"
        )
        beta_vdj_aa = (
            "KAGVTQTPRYLIKTRGQQVTLSCSPISGHRSVSWYQQTPGQGLQFLFEYFSETQRNKGNFPGRFSGRQFSNSRSEMNVSTLELGDSALYL"
            "CASSPFGGGINEQFFGPGTRLTVL"
        )
        utr5 = "TTCCCAGGAGGACCAAGCCC"
        alpha_cds_nt = (
            back_translate(alpha_leader_aa + alpha_vdj_aa, species="human")
            + back_translate("N" + HUMAN_CONSTANT_REGIONS_AA["TRAC"][:18], species="human")
        )
        beta_cds_nt = (
            back_translate(beta_leader_aa + beta_vdj_aa, species="human")
            + back_translate("E" + HUMAN_CONSTANT_REGIONS_AA["TRBC2"][:18], species="human")
        )
        alpha_seq = utr5 + alpha_cds_nt
        beta_seq = utr5 + beta_cds_nt
        a_j_end = len(utr5) + 3 * len(alpha_leader_aa + alpha_vdj_aa) + 1  # +1 overhang like real CellRanger
        b_j_end = len(utr5) + 3 * len(beta_leader_aa + beta_vdj_aa) + 1

        json_records = [
            {
                "contig_name": "BC1-1_contig_1",
                "barcode": "BC1-1",
                "is_cell": True,
                "high_confidence": True,
                "productive": True,
                "cdr3": "CAAGAGSQGNLIF",
                "start_codon_pos": len(utr5),
                "sequence": alpha_seq,
                "annotations": [
                    {"feature": {"region_type": "L-REGION+V-REGION", "gene_name": "TRAV27"}, "contig_match_start": len(utr5), "contig_match_end": a_j_end - 45},
                    {"feature": {"region_type": "J-REGION", "gene_name": "TRAJ42"}, "contig_match_start": a_j_end - 45, "contig_match_end": a_j_end},
                    {"feature": {"region_type": "C-REGION", "gene_name": "TRAC"}, "contig_match_start": a_j_end, "contig_match_end": len(alpha_seq)},
                ],
            },
            {
                "contig_name": "BC1-1_contig_2",
                "barcode": "BC1-1",
                "is_cell": True,
                "high_confidence": True,
                "productive": True,
                "cdr3": "CASSPFGGGINEQFF",
                "start_codon_pos": len(utr5),
                "sequence": beta_seq,
                "annotations": [
                    {"feature": {"region_type": "L-REGION+V-REGION", "gene_name": "TRBV5-1"}, "contig_match_start": len(utr5), "contig_match_end": b_j_end - 45},
                    {"feature": {"region_type": "J-REGION", "gene_name": "TRBJ2-1"}, "contig_match_start": b_j_end - 45, "contig_match_end": b_j_end},
                    {"feature": {"region_type": "C-REGION", "gene_name": "TRBC2"}, "contig_match_start": b_j_end, "contig_match_end": len(beta_seq)},
                ],
            },
        ]
        with gzip.open(sdir / "GSM1_all_contig_annotations.json.gz", "wt", encoding="utf-8") as fj:
            json.dump(json_records, fj)

        pd.DataFrame(
            {
                "barcode": ["BC1-1", "BC1-1"],
                "is_cell": [True, True],
                "contig_id": ["BC1-1_contig_1", "BC1-1_contig_2"],
                "high_confidence": [True, True],
                "chain": ["TRA", "TRB"],
                "v_gene": ["TRAV27", "TRBV5-1"],
                "d_gene": [None, "TRBD1"],
                "j_gene": ["TRAJ42", "TRBJ2-1"],
                "c_gene": ["TRAC", "TRBC2"],
                "productive": [True, True],
                "cdr3": ["CAAGAGSQGNLIF", "CASSPFGGGINEQFF"],
                "umis": [20, 30],
                "reads": [200, 300],
            }
        ).to_csv(sdir / "GSM1_filtered_contig_annotations.csv.gz", index=False, compression="gzip")

        vdj = load_cellranger_vdj(sdir, "GT6_stim", verbose=False)
        assert "vdj_aa" in vdj.columns and vdj["vdj_aa"].notna().all()
        assert vdj.loc[vdj["chain"] == "TRA", "vdj_aa"].iloc[0] == alpha_vdj_aa
        assert vdj.loc[vdj["chain"] == "TRB", "vdj_aa"].iloc[0] == beta_vdj_aa

        sample = Sample(sample="GT6_stim", vdj_dir=str(sdir))
        adata = load_sample(sample)
        clones = aggregate_clonotypes(adata, show_progress=False, verbose=False)
        assembled = assemble_full_sequences(
            clones,
            contigs_dir=tmp_path,
            alpha_leader="from_contig",
            beta_leader="from_contig",
            include_constant=True,
            linker="T2A",
            show_progress=False,
            verbose=False,
        )
        assert bool(assembled.iloc[0]["construct_contig_verified"]) is True
        assert assembled.iloc[0]["alpha_leader_aa"] == alpha_leader_aa
        assert assembled.iloc[0]["beta_leader_aa"] == beta_leader_aa

    def test_parse_fasta_gzipped(self, tmp_path):
        gz_fasta = tmp_path / "filtered_contig.fasta.gz"
        with gzip.open(gz_fasta, "wt", encoding="utf-8") as fh:
            fh.write(">c1 extra\nATGCGT\nTTAA\n>c2\nGGCC\n")
        parsed = parse_fasta(gz_fasta)
        assert parsed == {"c1": "ATGCGTTTAA", "c2": "GGCC"}
        contigs = load_contigs(tmp_path.parent, sample_name_from="parent")
        assert contigs[tmp_path.name]["c1"] == "ATGCGTTTAA"


class TestClonotypeSamplelessAndMultimerAndStimulation:
    """Test sampleless aggregate_clonotypes, multimer specificity calling, and paired stimulation delta."""

    def test_aggregate_clonotypes_without_sample_column(self):
        adata = ad.AnnData(
            X=sp.csr_matrix(np.ones((2, 6))),
            obs=pd.DataFrame(
                {
                    "CDR3_alpha": ["CAAGAGSQGNLIF", "CAAGAGSQGNLIF"],
                    "CDR3_beta": ["CASSIRSSYEQYF", "CASSIRSSYEQYF"],
                    "TRA_1_umis": [10, 12],
                    "TRB_1_umis": [15, 18],
                },
                index=["BC1-1", "BC2-1"],
            ),
            var=pd.DataFrame(index=["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"]),
        )
        clones = aggregate_clonotypes(adata, show_progress=False, verbose=False)
        assert len(clones) == 1
        assert clones.iloc[0]["cell_count"] == 2
        assert clones.iloc[0]["n_samples"] == 1

        # Also verify combine_gex_and_vdj propagates single sample from vdj_df
        adata2 = ad.AnnData(
            X=sp.csr_matrix(np.ones((1, 6))),
            obs=pd.DataFrame(index=["BC1-1"]),
            var=pd.DataFrame(index=["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"]),
        )
        vdj = pd.DataFrame(
            {
                "barcode": ["BC1-1", "BC1-1"],
                "chain": ["TRA", "TRB"],
                "cdr3": ["CAAGAGSQGNLIF", "CASSIRSSYEQYF"],
                "v_gene": ["TRAV27", "TRBV19"],
                "d_gene": [None, "TRBD1"],
                "j_gene": ["TRAJ42", "TRBJ2-7"],
                "c_gene": ["TRAC", "TRBC2"],
                "umis": [10, 15],
                "reads": [100, 150],
                "contig_id": ["c1", "c2"],
                "sample": ["donor1", "donor1"],
            }
        )
        comb = combine_gex_and_vdj(adata2, vdj)
        assert comb.obs.loc["BC1-1", "sample"] == "donor1"

    def test_multimer_extraction_classification_and_calling(self, tmp_path):
        mtx_dir = tmp_path / "outs" / "filtered_feature_bc_matrix"
        X = np.array(
            [
                [5, 120, 2, 1, 80],   # BC1: Flu specific
                [6, 95, 1, 2, 75],    # BC2: Flu specific (same clone)
                [4, 60, 55, 1, 50],   # BC3: poly-reactive (Flu ~ EBV)
                [7, 25, 2, 12, 60],   # BC4: high NC background (25 / 12 < 5)
                [5, 2, 1, 0, 90],     # BC5: negative
            ],
            dtype=float,
        )
        feats = [
            ("ENSG00000167286", "CD3D", "Gene Expression"),
            ("DEX_FLU", "A0201_GILGFVFTL_Flu", "Antibody Capture"),
            ("DEX_EBV", "A0201_GLCTLVAML_EBV", "Antibody Capture"),
            ("DEX_NC1", "A0201_ALIAPVHAV_NC", "Antibody Capture"),
            ("AB_CD8", "CD8a_TotalSeqC", "Antibody Capture"),
        ]
        barcodes = [f"BC{i}-1" for i in range(1, 6)]
        _write_mtx_dir(mtx_dir, X, barcodes, feats)

        fb = extract_feature_barcodes(tmp_path / "outs")
        assert list(fb.columns) == [
            "A0201_GILGFVFTL_Flu",
            "A0201_GLCTLVAML_EBV",
            "A0201_ALIAPVHAV_NC",
            "CD8a_TotalSeqC",
        ]
        classes = classify_feature_columns(fb.columns)
        assert classes["multimers"] == ["A0201_GILGFVFTL_Flu", "A0201_GLCTLVAML_EBV"]
        assert classes["negative_controls"] == ["A0201_ALIAPVHAV_NC"]
        assert classes["surface_proteins"] == ["CD8a_TotalSeqC"]

        calls = call_multimer_specificity(fb, min_umi=10.0, min_nc_ratio=5.0, min_dominance_ratio=2.0)
        assert calls.loc["BC1-1", "multimer_call"] == "A0201_GILGFVFTL_Flu"
        assert calls.loc["BC2-1", "multimer_call"] == "A0201_GILGFVFTL_Flu"
        assert bool(calls.loc["BC3-1", "is_multi_reactive"]) is True
        assert calls.loc["BC3-1", "multimer_call"] is None
        assert calls.loc["BC4-1", "multimer_call"] is None
        assert calls.loc["BC5-1", "multimer_call"] is None

        obs = pd.DataFrame(
            {
                "CDR3ab": [
                    "CAAGAGSQGNLIF_CASSIRSSYEQYF",
                    "CAAGAGSQGNLIF_CASSIRSSYEQYF",
                    "CAVPOLY_CASSPOLY",
                    "CAVSTICKY_CASSSTICKY",
                    "CAVSTICKY_CASSSTICKY",
                ]
            },
            index=barcodes,
        )
        clone_summary = summarize_clone_multimer_specificity(obs, calls)
        flu_row = clone_summary[clone_summary["CDR3ab"] == "CAAGAGSQGNLIF_CASSIRSSYEQYF"].iloc[0]
        assert bool(flu_row["multimer_confident"]) is True
        assert flu_row["multimer_concordance"] == pytest.approx(1.0)
        sticky_row = clone_summary[clone_summary["CDR3ab"] == "CAVSTICKY_CASSSTICKY"].iloc[0]
        assert bool(sticky_row["multimer_confident"]) is False

    def test_score_stimulation_response_separates_activated_from_bystander(self):
        # Simulate overnight peptide stimulation:
        # - CLONE_REACTIVE has 4 cells in unstim (AIM=-0.2) and 5 cells in stim (AIM=+1.8) -> strong delta_AIM
        # - CLONE_EXPANDED has 0 cells in unstim and 12 cells in stim (AIM=+1.5) -> strong Fisher enrichment
        # - CLONE_BYSTANDER has 20 cells in unstim (AIM=0.0) and 20 cells in stim (AIM=0.05) -> unresponsive bystander
        rows = []
        for i in range(4):
            rows.append(("unstim", "control", "CLONE_REACTIVE", -0.2))
        for i in range(5):
            rows.append(("pep_stim", "culture", "CLONE_REACTIVE", 1.8))
        for i in range(12):
            rows.append(("pep_stim", "culture", "CLONE_EXPANDED", 1.5))
        for i in range(20):
            rows.append(("unstim", "control", "CLONE_BYSTANDER", 0.0))
        for i in range(20):
            rows.append(("pep_stim", "culture", "CLONE_BYSTANDER", 0.05))
        for i in range(60):
            rows.append(("unstim", "control", f"SINGLE_U_{i}", -0.1))
            rows.append(("pep_stim", "culture", f"SINGLE_S_{i}", -0.1))

        obs = pd.DataFrame(rows, columns=["sample", "source", "CDR3ab", "signature_AIM"])
        resp = score_stimulation_response(obs, score_cols=["signature_AIM"]).set_index("CDR3ab")

        assert resp.loc["CLONE_EXPANDED", "log2_fold_change"] > 3.5
        assert resp.loc["CLONE_EXPANDED", "fisher_fdr"] < 0.05
        assert resp.loc["CLONE_REACTIVE", "delta_signature_AIM"] == pytest.approx(2.0)
        assert abs(resp.loc["CLONE_BYSTANDER", "delta_signature_AIM"]) < 0.1
        assert resp.loc["CLONE_BYSTANDER", "fisher_pvalue"] > 0.4

    def test_audit_regressions_partial_segments_multimer_and_control_exclusion(self, tmp_path):
        # 1. CR6+ 5'-partial contig with missing fwr1 must NOT emit truncated vdj_aa,
        #    and full-length sibling cell / contig must take priority in clonotyping.
        vdj_csv = tmp_path / "filtered_contig_annotations.csv"
        pd.DataFrame(
            [
                # BC1 has higher UMI on a 5'-partial contig (productive=None, fwr1=None)
                {
                    "barcode": "BC1-1",
                    "is_cell": True,
                    "contig_id": "BC1_c1",
                    "high_confidence": True,
                    "chain": "TRA",
                    "v_gene": "TRAV27",
                    "d_gene": None,
                    "j_gene": "TRAJ42",
                    "c_gene": "TRAC",
                    "full_length": False,
                    "productive": "None",
                    "fwr1": None,
                    "cdr1": "SVFSS",
                    "fwr2": "WYRQ",
                    "cdr2": "VVTG",
                    "fwr3": "RFTI",
                    "cdr3": "CAAGAGSQGNLIF",
                    "fwr4": "FGKG",
                    "umis": 30,
                    "reads": 300,
                },
                {
                    "barcode": "BC1-1",
                    "is_cell": True,
                    "contig_id": "BC1_c2",
                    "high_confidence": True,
                    "chain": "TRB",
                    "v_gene": "TRBV19",
                    "d_gene": "TRBD1",
                    "j_gene": "TRBJ2-7",
                    "c_gene": "TRBC2",
                    "full_length": True,
                    "productive": True,
                    "fwr1": "NAGVT",
                    "cdr1": "LNHDA",
                    "fwr2": "WYQQ",
                    "cdr2": "SQIV",
                    "fwr3": "GDSV",
                    "cdr3": "CASSIRSSYEQYF",
                    "fwr4": "FGPG",
                    "umis": 25,
                    "reads": 250,
                },
                # BC2 has lower UMI on TRA, but complete fwr1..fwr4!
                {
                    "barcode": "BC2-1",
                    "is_cell": True,
                    "contig_id": "BC2_c1",
                    "high_confidence": True,
                    "chain": "TRA",
                    "v_gene": "TRAV27",
                    "d_gene": None,
                    "j_gene": "TRAJ42",
                    "c_gene": "TRAC",
                    "full_length": True,
                    "productive": True,
                    "fwr1": "TQLLE",
                    "cdr1": "SVFSS",
                    "fwr2": "WYRQ",
                    "cdr2": "VVTG",
                    "fwr3": "RFTI",
                    "cdr3": "CAAGAGSQGNLIF",
                    "fwr4": "FGKG",
                    "umis": 10,
                    "reads": 100,
                },
                {
                    "barcode": "BC2-1",
                    "is_cell": True,
                    "contig_id": "BC2_c2",
                    "high_confidence": True,
                    "chain": "TRB",
                    "v_gene": "TRBV19",
                    "d_gene": "TRBD1",
                    "j_gene": "TRBJ2-7",
                    "c_gene": "TRBC2",
                    "full_length": True,
                    "productive": True,
                    "fwr1": "NAGVT",
                    "cdr1": "LNHDA",
                    "fwr2": "WYQQ",
                    "cdr2": "SQIV",
                    "fwr3": "GDSV",
                    "cdr3": "CASSIRSSYEQYF",
                    "fwr4": "FGPG",
                    "umis": 12,
                    "reads": 120,
                },
            ]
        ).to_csv(vdj_csv, index=False)

        loaded_vdj = load_cellranger_vdj(vdj_csv, "s1", verbose=False).set_index("contig_id")
        assert pd.isna(loaded_vdj.loc["BC1_c1", "vdj_aa"])
        assert loaded_vdj.loc["BC2_c1", "vdj_aa"] == "TQLLESVFSSWYRQVVTGRFTICAAGAGSQGNLIFFGKG"

        adata_mix = ad.AnnData(
            X=sp.csr_matrix(np.ones((2, 3))),
            obs=pd.DataFrame({"sample": ["s1", "s1"]}, index=["BC1-1", "BC2-1"]),
            var=pd.DataFrame(index=["CD3D", "CD8A", "CD8B"]),
        )
        comb_mix = combine_gex_and_vdj(adata_mix, loaded_vdj.reset_index())
        agg_mix = aggregate_clonotypes(comb_mix, show_progress=False, verbose=False)
        assert agg_mix.iloc[0]["VDJ_alpha_aa"] == "TQLLESVFSSWYRQVVTGRFTICAAGAGSQGNLIFFGKG"

        # 2. classify_feature_columns handles Hashtag1, CD45/HLA-DR, BEAM-T Antigen Capture, & RefSeq NC_
        fb_test = pd.DataFrame(
            [[0.0, 20.0, 500.0, 30.0, 15.0, 1.0]],
            index=["BC1-1"],
            columns=[
                "SARS2_NC_045512_Spike",
                "BEAM_Pep1_TotalSeqC",
                "Hashtag1",
                "CD45",
                "HLA-DR",
                "Blank_HLA_A02",
            ],
        )
        fb_test.attrs["feature_types"] = {
            "SARS2_NC_045512_Spike": "Antigen Capture",
            "BEAM_Pep1_TotalSeqC": "Antigen Capture",
            "Hashtag1": "Antibody Capture",
            "CD45": "Antibody Capture",
            "HLA-DR": "Antibody Capture",
            "Blank_HLA_A02": "Antigen Capture",
        }
        fc = classify_feature_columns(fb_test)
        assert fc["multimers"] == ["SARS2_NC_045512_Spike", "BEAM_Pep1_TotalSeqC"]
        assert fc["negative_controls"] == ["Blank_HLA_A02"]
        assert fc["hashtags"] == ["Hashtag1"]
        assert set(fc["surface_proteins"]) == {"CD45", "HLA-DR"}

        # Duplicate barcodes across unaligned obs_df and calls_df must raise TCRsiftValidationError
        dup_obs = pd.DataFrame({"CDR3ab": ["CA_CB", "CA2_CB2"]}, index=["BC1-1", "BC1-1"])
        dup_calls = pd.DataFrame({"multimer_call": ["M1"]}, index=["BC1-1"])
        with pytest.raises(TCRsiftValidationError, match="Duplicate cell barcodes"):
            summarize_clone_multimer_specificity(dup_obs, dup_calls)

        # Partial single-chain "_CASS..." entries must be excluded when clone_col == "CDR3ab"
        part_obs = pd.DataFrame(
            {
                "sample": ["ctrl", "stim", "ctrl", "stim", "stim"],
                "source": ["control", "culture", "control", "culture", "culture"],
                "CDR3ab": ["_CASSONLY", "_CASSONLY", "CAV_CASS", "CAV_CASS", "CAV_CASS"],
            }
        )
        part_resp = score_stimulation_response(part_obs)
        assert list(part_resp["CDR3ab"]) == ["CAV_CASS"]
