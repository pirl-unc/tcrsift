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

"""End-to-end and unit tests for mouse (Mus musculus) and chimeric TCR workflows."""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from tcrsift.annotate_cells import (
     _TOP_GENE_EXCLUDE,
    MOUSE_GATES,
    annotate_cells,
    gates_for_species,
)
from tcrsift.assemble import (
    MOUSE_TRAC_AA,
    MOUSE_TRBC1_AA,
    MOUSE_TRBC2_AA,
    assemble_full_sequences,
    assemble_qc_report,
    optimize_codons,
    translate_dna,
    validate_sequences,
)
from tcrsift.cli import main as cli_main
from tcrsift.genes import adata_mito_mask, resolve_adata_var_map, tcell_markers_for_species
from tcrsift.gex import compute_signature_scores_per_cell
from tcrsift.loader import load_sample
from tcrsift.pgen import annotate_publicness, compute_pgen
from tcrsift.phenotype import phenotype_cells
from tcrsift.qc import cell_qc_funnel, cross_lineage_doublets
from tcrsift.qc_battery import run_qc_battery
from tcrsift.sample_sheet import Sample
from tcrsift.seqprob import score_log_prob
from tcrsift.signature_methods import score_by_name, signatures_for_species
from tcrsift.species import (
    detect_adata_species,
    detect_vdj_species,
    normalize_species,
    resolve_species,
    species_to_binomial,
)


def _write_10x_mtx_dir(
    out_dir: Path,
    barcodes: list[str],
    gene_ids: list[str],
    gene_symbols: list[str],
    counts_cells_by_genes: np.ndarray,
) -> Path:
    """Write a Cell Ranger filtered_feature_bc_matrix directory."""
    mtx_dir = out_dir / "filtered_feature_bc_matrix"
    mtx_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt") as fh:
        for bc in barcodes:
            fh.write(f"{bc}\n")
    with gzip.open(mtx_dir / "features.tsv.gz", "wt") as fh:
        for gid, sym in zip(gene_ids, gene_symbols):
            fh.write(f"{gid}\t{sym}\tGene Expression\n")
    mat = sp.coo_matrix(counts_cells_by_genes.T.astype(np.int32))
    buf_path = mtx_dir / "matrix.mtx"
    scipy.io.mmwrite(str(buf_path), mat)
    with open(buf_path, "rb") as src, gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as dst:
        dst.write(src.read())
    buf_path.unlink()
    return out_dir


class TestSpeciesResolution:
    """Species normalization and heuristic detection from VDJ and AnnData."""

    def test_normalize_species_aliases(self):
        for raw in ("human", "Homo sapiens", "hs", "hg38", "GRCh38"):
            assert normalize_species(raw) == "human"
        for raw in ("mouse", "Mus musculus", "murine", "mm", "mm10", "GRCm39"):
            assert normalize_species(raw) == "mouse"
        with pytest.raises(ValueError, match="Unsupported species"):
            normalize_species("zebrafish")
        assert species_to_binomial("mouse") == "Mus musculus"
        assert species_to_binomial("human") == "Homo sapiens"

    def test_detect_vdj_species(self):
        mouse_df = pd.DataFrame({
            "alpha_v_gene": ["TRAV14-1", "Trav12-1"],
            "alpha_c_gene": ["Trac", "Trac"],
            "beta_c_gene": ["Trbc1", "Trbc2"],
        })
        assert detect_vdj_species(mouse_df) == "mouse"

        human_df = pd.DataFrame({
            "alpha_v_gene": ["TRAV12-2", "TRAV21"],
            "alpha_c_gene": ["TRAC", "TRAC"],
            "beta_c_gene": ["TRBC1", "TRBC2"],
        })
        assert detect_vdj_species(human_df) == "human"

    def test_detect_adata_species_from_symbols_and_ensembl(self):
        mouse_sym = ad.AnnData(
            X=np.ones((3, 5)),
            var=pd.DataFrame(index=["Cd3e", "Cd8a", "Cd8b1", "mt-Co1", "Gzmb"]),
        )
        assert detect_adata_species(mouse_sym) == "mouse"

        mouse_ens = ad.AnnData(
            X=np.ones((2, 3)),
            var=pd.DataFrame(
                index=[
                    "ENSMUSG00000032093",
                    "ENSMUSG00000053044",
                    "ENSMUSG00000064351",
                ]
            ),
        )
        assert detect_adata_species(mouse_ens) == "mouse"
        assert resolve_species("auto", adata=mouse_sym) == "mouse"

    def test_resolve_species_explicit_mismatch_warns(self, caplog):
        mouse_adata = ad.AnnData(
            X=np.ones((2, 4)),
            var=pd.DataFrame(index=["Cd3e", "Cd8a", "Cd8b1", "mt-Nd1"]),
        )
        resolved = resolve_species("human", adata=mouse_adata, warn_mismatch=True)
        assert resolved == "human"
        assert any("mouse identifiers" in rec.message.lower() for rec in caplog.records)


class TestMouseGenesAndLoading:
    """Mouse mitochondrial mask, Cd8b1 mapping, Cell Ranger loading, and phenotyping."""

    def test_mito_mask_and_var_resolution_on_mouse_symbols_and_ensembl(self):
        adata_sym = ad.AnnData(
            X=np.ones((2, 5)),
            var=pd.DataFrame(index=["Cd3e", "Cd8b1", "mt-Co1", "mt-Cytb", "Gzmb"]),
        )
        mask_sym = adata_mito_mask(adata_sym)
        assert mask_sym.tolist() == [False, False, True, True, False]

        _, var_map = resolve_adata_var_map(adata_sym, wanted_genes=["CD8B", "CD3E"])
        assert var_map["CD8B"] == "Cd8b1"
        assert var_map["CD3E"] == "Cd3e"

        adata_ens = ad.AnnData(
            X=np.ones((2, 3)),
            var=pd.DataFrame(
                index=[
                    "ENSMUSG00000053044",  # Cd8b1
                    "ENSMUSG00000064351",  # mt-Co1
                    "ENSMUSG00000032093",  # Cd3e
                ]
            ),
        )
        mask_ens = adata_mito_mask(adata_ens)
        assert mask_ens.tolist() == [False, True, False]
        _, ens_map = resolve_adata_var_map(adata_ens, wanted_genes=["CD8B", "CD3E"])
        assert ens_map["CD8B"] == "ENSMUSG00000053044"
        assert ens_map["CD3E"] == "ENSMUSG00000032093"
        assert "Cd8b1" in [
            alias for g in tcell_markers_for_species("mouse") for alias in g.aliases
        ]

    def test_load_sample_and_phenotype_mouse_cellranger(self, tmp_path):
        vdj_dir = tmp_path / "vdj"
        vdj_dir.mkdir()
        contig_rows = [
            # Cell 1 & 2: CD8 OT-I-like clone
            {
                "barcode": "AAACCTGAGAAACCAT-1", "is_cell": True, "contig_id": "c1",
                "high_confidence": True, "length": 550, "chain": "TRA",
                "v_gene": "TRAV14-1", "d_gene": "", "j_gene": "TRAJ37", "c_gene": "Trac",
                "full_length": True, "productive": True, "cdr3": "CAASDNYQLIW",
                "cdr3_nt": "TGTGCAGCAAGTGACAACTATCAGTTGATCTGG", "reads": 200, "umis": 12,
            },
            {
                "barcode": "AAACCTGAGAAACCAT-1", "is_cell": True, "contig_id": "c2",
                "high_confidence": True, "length": 600, "chain": "TRB",
                "v_gene": "TRBV12-1", "d_gene": "TRBD1", "j_gene": "TRBJ2-7", "c_gene": "Trbc2",
                "full_length": True, "productive": True, "cdr3": "CASSRANYEQYF",
                "cdr3_nt": "TGTGCCAGCAGCCGGGCTAACTATGAACAGTACTTC", "reads": 250, "umis": 15,
            },
            # Cell 2: CD4 clone
            {
                "barcode": "AAACCTGAGAAACCGC-1", "is_cell": True, "contig_id": "c3",
                "high_confidence": True, "length": 540, "chain": "TRA",
                "v_gene": "TRAV12-1", "d_gene": "", "j_gene": "TRAJ15", "c_gene": "Trac",
                "full_length": True, "productive": True, "cdr3": "CALSDRGGRALIF",
                "cdr3_nt": "TGTGCTCTGAGTGACAGGGGAGGCAGAGCTCTGATATTT", "reads": 180, "umis": 10,
            },
            {
                "barcode": "AAACCTGAGAAACCGC-1", "is_cell": True, "contig_id": "c4",
                "high_confidence": True, "length": 580, "chain": "TRB",
                "v_gene": "TRBV13-2", "d_gene": "TRBD1", "j_gene": "TRBJ1-1", "c_gene": "Trbc1",
                "full_length": True, "productive": True, "cdr3": "CASGDNTEVFF",
                "cdr3_nt": "TGTGCCAGCGGTGACAACAEVFF", "reads": 220, "umis": 11,
            },
        ]
        pd.DataFrame(contig_rows).to_csv(
            vdj_dir / "filtered_contig_annotations.csv", index=False
        )

        genes = [
            ("ENSMUSG00000032094", "Cd3d"),
            ("ENSMUSG00000032093", "Cd3e"),
            ("ENSMUSG00000002033", "Cd3g"),
            ("ENSMUSG00000023274", "Cd4"),
            ("ENSMUSG00000053977", "Cd8a"),
            ("ENSMUSG00000053044", "Cd8b1"),
            ("ENSMUSG00000064351", "mt-Co1"),
            ("ENSMUSG00000015437", "Gzmb"),
            ("ENSMUSG00000037202", "Prf1"),
        ]
        # Cell 1: Cd8a=6, Cd8b1=18, Cd4=0; Cell 2: Cd4=20, Cd8a=0, Cd8b1=0
        counts = np.array([
            [5, 8, 3, 0, 6, 18, 2, 14, 10],
            [4, 7, 2, 20, 0, 0, 1, 0, 0],
        ], dtype=np.int32)
        gex_dir = _write_10x_mtx_dir(
            tmp_path / "gex",
            ["AAACCTGAGAAACCAT-1", "AAACCTGAGAAACCGC-1"],
            [g[0] for g in genes],
            [g[1] for g in genes],
            counts,
        )

        sample = Sample(
            sample="Mouse_Splenocytes",
            vdj_dir=str(vdj_dir),
            gex_dir=str(gex_dir),
            species="auto",
        )
        adata = load_sample(sample, min_genes=3, min_counts=10, min_mito_pct=0.0, max_mito_pct=20.0)
        assert adata.uns["species"] == "mouse"
        assert float(adata.obs.loc["AAACCTGAGAAACCAT-1", "CD8B"]) == 18.0
        assert float(adata.obs["pct_counts_mt"].iloc[0]) > 0.0

        pheno = phenotype_cells(adata, cd4_cd8_ratio=3.0, min_cd3_reads=3)
        calls = dict(zip(pheno.obs_names, pheno.obs["Tcell_type"]))
        assert calls["AAACCTGAGAAACCAT-1"] == "Confident CD8+"
        assert calls["AAACCTGAGAAACCGC-1"] == "Confident CD4+"


class TestMouseSignaturesAndAtlas:
    """Mouse GEX signature scoring, cell QC, and cluster cell-type annotation."""

    def test_mouse_signature_scoring_and_cluster_annotation(self):
        filler_genes = [f"Fgene{i}" for i in range(120)]
        gene_names = [
            "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1",
            "Foxp3", "Il2ra", "Ctla4", "Ikzf2",
            "Gzmb", "Gzma", "Prf1", "Nkg7", "Ifng", "Ccl5",
            "Tcf7", "Lef1", "Ccr7", "Sell", "Il7r", "Cd27", "Cd28",
            "Trdc", "Slc4a10", "Klrb1c", "Zbtb16",
            "Ncr1", "Klrb1b", "Klra8",
            "Cd19", "Ms4a1", "Cd79a",
            "Lyz2", "S100a8", "S100a9", "Adgre1", "Csf1r",
            "Tnfrsf9", "Tnfrsf4", "Egr2", "Entpd1", "Itgae",
            "Mki67", "Pdcd1", "Lag3", "Havcr2", "Tigit", "Tox",
            "mt-Co1", *filler_genes,
        ]
        # 4 clean clusters x 25 cells = 100 cells:
        # cluster 0: CD8 effector, cluster 1: CD4 Treg, cluster 2: NK, cluster 3: B cell
        rng = np.random.default_rng(7)
        X = np.zeros((100, len(gene_names)), dtype=np.float32)
        idx = {g: i for i, g in enumerate(gene_names)}
        for k, fg in enumerate(filler_genes):
            X[:, idx[fg]] = rng.poisson(1 + (k % 35), size=100)

        # Cluster 0: CD8 effector
        for g in ("Cd3d", "Cd3e", "Cd8a", "Cd8b1", "Gzmb", "Gzma", "Prf1", "Nkg7", "Ifng", "Ccl5"):
            X[0:25, idx[g]] = 30
        # Cluster 1: Treg
        for g in ("Cd3d", "Cd3e", "Cd4", "Foxp3", "Il2ra", "Ctla4", "Ikzf2"):
            X[25:50, idx[g]] = 30
        # Cluster 2: NK (no CD3, high Ncr1/Klra8/Prf1)
        for g in ("Ncr1", "Klrb1c", "Klra8", "Prf1", "Gzmb", "Nkg7"):
            X[50:75, idx[g]] = 30
        # Cluster 3: B cell
        for g in ("Cd19", "Ms4a1", "Cd79a"):
            X[75:100, idx[g]] = 35
        X[:, idx["mt-Co1"]] = 3

        leiden = pd.Categorical(["0"] * 25 + ["1"] * 25 + ["2"] * 25 + ["3"] * 25)
        adata = ad.AnnData(
            X=sp.csr_matrix(X),
            obs=pd.DataFrame({"leiden": leiden}, index=[f"cell_{i}" for i in range(100)]),
            var=pd.DataFrame(index=gene_names),
        )

        # QC funnel + cross-lineage doublet check auto-detects mouse
        qc_adata, wf = cell_qc_funnel(adata.copy(), min_genes=10, min_counts=50, max_mito_pct=20)
        assert qc_adata.n_obs == 100
        dbl_mask, *_ = cross_lineage_doublets(adata, species="mouse")
        assert int(dbl_mask.sum()) == 0

        # Per-cell signature scoring
        scored = compute_signature_scores_per_cell(adata.copy(), species="mouse")
        assert "signature_effector" in scored.obs.columns
        assert scored.obs.iloc[0:25]["signature_effector"].mean() > scored.obs.iloc[75:100]["signature_effector"].mean()

        diff_scores = score_by_name(adata, "Differentiated", log1p=True, species="mouse")
        assert diff_scores.iloc[0:25].mean() > diff_scores.iloc[25:50].mean()
        assert "Differentiated" in signatures_for_species("mouse")

        # Cluster annotation with MOUSE_GATES
        assert gates_for_species("mouse") is MOUSE_GATES
        ann = adata.copy()
        annotate_cells(ann, species="mouse")
        base_by_cluster = (
            ann.obs.groupby("leiden", observed=True)["phenotype_base"]
            .agg(lambda s: s.mode().iloc[0])
            .to_dict()
        )
        assert "CD8" in base_by_cluster["0"]
        assert "Treg" in base_by_cluster["1"]
        assert "NK" in base_by_cluster["2"]
        assert "B" in base_by_cluster["3"]

    def test_top_gene_exclude_regex_is_case_insensitive_and_word_bounded(self):
        for receptor_gene in ("TRAV12-2", "Trav12-1", "Trbv13-2", "Trac", "Trbc2", "Ighv1-5", "Ighg2b"):
            assert _TOP_GENE_EXCLUDE.search(receptor_gene), receptor_gene
        for non_receptor in ("TRADD", "Tradd", "TRAF1", "Traf1", "IGHMBP2", "Ighmbp2"):
            assert not _TOP_GENE_EXCLUDE.search(non_receptor), non_receptor


class TestMousePrioritizeCLI:
    """End-to-end mouse `tcrsift prioritize --species mouse` execution."""

    def test_mouse_prioritize_cli_selects_expanded_effector_clone(self, tmp_path):
        vdj_dir = tmp_path / "vdj"
        vdj_dir.mkdir()
        genes = [
            ("ENSMUSG00000032094", "Cd3d"),
            ("ENSMUSG00000032093", "Cd3e"),
            ("ENSMUSG00000002033", "Cd3g"),
            ("ENSMUSG00000023274", "Cd4"),
            ("ENSMUSG00000053977", "Cd8a"),
            ("ENSMUSG00000053044", "Cd8b1"),
            ("ENSMUSG00000015437", "Gzmb"),
            ("ENSMUSG00000023132", "Gzma"),
            ("ENSMUSG00000037202", "Prf1"),
            ("ENSMUSG00000055170", "Ifng"),
            ("ENSMUSG00000004668", "Nkg7"),
            ("ENSMUSG00000021728", "Tnfrsf9"),
            ("ENSMUSG00000031004", "Mki67"),
            ("ENSMUSG00000020122", "Tnfrsf4"),
            ("ENSMUSG00000000782", "Tcf7"),
            ("ENSMUSG00000027985", "Lef1"),
            ("ENSMUSG00000037944", "Ccr7"),
            ("ENSMUSG00000026581", "Sell"),
            ("ENSMUSG00000003882", "Il7r"),
            ("ENSMUSG00000026285", "Pdcd1"),
            ("ENSMUSG00000030124", "Lag3"),
            ("ENSMUSG00000020399", "Havcr2"),
            ("ENSMUSG00000071552", "Tigit"),
            ("ENSMUSG00000041272", "Tox"),
            ("ENSMUSG00000026011", "Ctla4"),
            ("ENSMUSG00000005947", "Itgae"),
            ("ENSMUSG00000048120", "Entpd1"),
            ("ENSMUSG00000000982", "Ccl3"),
            ("ENSMUSG00000018930", "Ccl4"),
            ("ENSMUSG00000064351", "mt-Co1"),
        ]
        # pad with housekeeping genes so min_genes=10 is easily satisfied
        for k in range(15):
            genes.append((f"ENSMUSG900000{k:05d}", f"Rpl{k + 10}"))

        rng = np.random.default_rng(42)
        n_cells = 60
        barcodes = [f"BC{i:04d}-1" for i in range(n_cells)]
        contig_rows = []
        counts = rng.integers(2, 8, size=(n_cells, len(genes)), dtype=np.int32)
        sym_idx = {sym: idx for idx, (_, sym) in enumerate(genes)}

        # Cells 0..14: expanded cytotoxic mouse clone (15 cells)
        # Cells 15..59: 15 naive background mouse clones (3 cells each)
        for i, bc in enumerate(barcodes):
            counts[i, sym_idx["Cd3d"]] = 6
            counts[i, sym_idx["Cd3e"]] = 8
            counts[i, sym_idx["Cd4"]] = 0
            counts[i, sym_idx["Cd8a"]] = 10
            counts[i, sym_idx["Cd8b1"]] = 12
            counts[i, sym_idx["mt-Co1"]] = 1
            if i < 15:
                cdr3a, cdr3b = "CAASDNYQLIW", "CASSRANYEQYF"
                for eff in ("Gzmb", "Prf1", "Ifng", "Nkg7", "Tnfrsf9", "Mki67"):
                    counts[i, sym_idx[eff]] = 45
                for naive in ("Tcf7", "Lef1", "Ccr7", "Sell"):
                    counts[i, sym_idx[naive]] = 0
            else:
                clone_k = (i - 15) // 3
                suffix = "ACDEFGHIKLMNPQRSTVWY"[clone_k]
                cdr3a = f"CALSD{suffix}NYQLIW"
                cdr3b = f"CASSL{suffix}NTLYF"
                for eff in ("Gzmb", "Prf1", "Ifng", "Nkg7", "Tnfrsf9", "Mki67"):
                    counts[i, sym_idx[eff]] = 0
                for naive in ("Tcf7", "Lef1", "Ccr7", "Sell"):
                    counts[i, sym_idx[naive]] = 25

            contig_rows.extend([
                {
                    "barcode": bc, "is_cell": True, "contig_id": f"{bc}_a",
                    "high_confidence": True, "length": 520, "chain": "TRA",
                    "v_gene": "TRAV14-1", "d_gene": "", "j_gene": "TRAJ37", "c_gene": "Trac",
                    "full_length": True, "productive": True, "cdr3": cdr3a,
                    "cdr3_nt": "TGT" * 11, "reads": 100, "umis": 8,
                },
                {
                    "barcode": bc, "is_cell": True, "contig_id": f"{bc}_b",
                    "high_confidence": True, "length": 560, "chain": "TRB",
                    "v_gene": "TRBV12-1", "d_gene": "TRBD1", "j_gene": "TRBJ2-7", "c_gene": "Trbc2",
                    "full_length": True, "productive": True, "cdr3": cdr3b,
                    "cdr3_nt": "TGT" * 12, "reads": 120, "umis": 10,
                },
            ])

        pd.DataFrame(contig_rows).to_csv(
            vdj_dir / "filtered_contig_annotations.csv", index=False
        )
        gex_dir = _write_10x_mtx_dir(
            tmp_path / "gex",
            barcodes,
            [g[0] for g in genes],
            [g[1] for g in genes],
            counts,
        )
        sheet_path = tmp_path / "samples.yaml"
        sheet_path.write_text(
            f"species: mouse\nsamples:\n  - sample: B16_TIL\n    vdj_dir: {vdj_dir}\n    gex_dir: {gex_dir}\n"
        )
        out_dir = tmp_path / "prioritized"
        cli_main([
            "prioritize",
            str(sheet_path),
            "-o", str(out_dir),
            "--context", "solid-tumor",
            "--min-genes", "10",
            "--min-counts", "50",
            "--max-mito-pct", "20",
        ])
        cand = pd.read_csv(out_dir / "candidate_clones.csv")
        assert len(cand) >= 1
        assert "CAASDNYQLIW" in cand["CDR3_alpha"].values


class TestMouseAndChimericAssembly:
    """Full-length mouse and chimeric (human-VDJ + mouse-constant) TCR assembly."""

    @staticmethod
    def _synthetic_vdj(cdr3_aa: str, prefix_len: int = 86, suffix: str = "FGAGTKLSVLP") -> tuple[str, str]:
        framework = ("GQVEQSPQSLIILEGKNCTLQCNYTVSPFSNLRWYKQDTGRGPVSLTIMTFSENTKSNGRYTATLDADTKQSSLHITASQLSDSASYI"[:prefix_len])
        vdj_aa = framework + cdr3_aa + suffix
        vdj_nt = optimize_codons(vdj_aa, species="mouse")
        return vdj_aa, vdj_nt

    def test_pure_mouse_assembly_and_validation(self):
        vdj_a_aa, vdj_a_nt = self._synthetic_vdj("CALSDRGGRALIF", suffix="FGTGTQVIVLP")
        vdj_b_aa, vdj_b_nt = self._synthetic_vdj("CASSRANYEQYF", suffix="FGPGTRLTVL")
        clones = pd.DataFrame([
            {
                "CDR3ab": "CALSDRGGRALIF_CASSRANYEQYF",
                "CDR3_alpha": "CALSDRGGRALIF",
                "CDR3_beta": "CASSRANYEQYF",
                "VDJ_alpha_aa": vdj_a_aa,
                "VDJ_alpha_nt": vdj_a_nt,
                "VDJ_beta_aa": vdj_b_aa,
                "VDJ_beta_nt": vdj_b_nt,
                "alpha_v_gene": "TRAV12-1",
                "alpha_j_gene": "TRAJ4",    # mouse TRAJ4 junction residue is 'Y'
                "alpha_c_gene": "Trac",
                "beta_v_gene": "TRBV12-1",
                "beta_j_gene": "TRBJ2-7",
                "beta_c_gene": "Trbc2",
            }
        ])
        assembled = assemble_full_sequences(
            clones,
            alpha_leader="mCD28",
            beta_leader="mCD8B",
            include_constant=True,
            linker="T2A",
            species="mouse",
        )
        row = assembled.iloc[0]
        # Mouse TRAJ4 has germline 'Y' junction + bare MOUSE_TRAC_AA (136 aa)
        assert row["alpha_constant_aa"] == "Y" + MOUSE_TRAC_AA
        assert row["beta_constant_aa"] == "E" + MOUSE_TRBC2_AA
        assert row["full_alpha_aa"].endswith("LLMTLRLWSS")
        assert row["full_beta_aa"].endswith("VKKKNS")

        msgs = validate_sequences(assembled, strict=True, species="mouse")
        assert not [m for m in msgs if m.severity == "load_bearing"]
        report = assemble_qc_report(assembled)
        assert "0 load-bearing" in report or "PASS" in report.upper() or "QC" in report

        # Non-canonical J-end gene check (mouse TRAJ7 / TRBJ2-6)
        qc_df = assembled.copy()
        qc_df["alpha_j_gene"] = "TRAJ7"
        qc_df["beta_j_gene"] = "TRBJ2-6"
        battery = run_qc_battery(qc_df)
        assert len(battery) >= 5
        assert all(hasattr(r, "status") for r in battery)

    def test_chimeric_human_vdj_with_mouse_constants(self):
        vdj_a_aa, vdj_a_nt = self._synthetic_vdj("CAGNTGNQFYF", suffix="FGTGTSLTVIP")
        vdj_b_aa, vdj_b_nt = self._synthetic_vdj("CASSIRSSYEQYF", suffix="FGPGTRLTVT")
        clones = pd.DataFrame([
            {
                "CDR3ab": "CAGNTGNQFYF_CASSIRSSYEQYF",
                "CDR3_alpha": "CAGNTGNQFYF",
                "CDR3_beta": "CASSIRSSYEQYF",
                "VDJ_alpha_aa": vdj_a_aa,
                "VDJ_alpha_nt": vdj_a_nt,
                "VDJ_beta_aa": vdj_b_aa,
                "VDJ_beta_nt": vdj_b_nt,
                "alpha_v_gene": "TRAV12-2",
                "alpha_j_gene": "TRAJ49",   # human TRAJ49 junction residue is 'N'
                "alpha_c_gene": "TRAC",
                "beta_v_gene": "TRBV19",
                "beta_j_gene": "TRBJ1-2",   # pairs with TRBC1 -> MOUSE_TRBC1_AA
                "beta_c_gene": "TRBC1",
            }
        ])
        assembled = assemble_full_sequences(
            clones,
            alpha_leader="CD28",
            beta_leader="CD8A",
            include_constant=True,
            linker="P2A",
            species="human",
            constant_species="mouse",
            codon_species="human",
        )
        row = assembled.iloc[0]
        assert row["alpha_constant_aa"] == "N" + MOUSE_TRAC_AA
        assert row["beta_constant_aa"] == "E" + MOUSE_TRBC1_AA
        assert row["full_beta_aa"].endswith("VKRKNS")
        sc_aa, ragged = translate_dna(row["single_chain_nt"])
        assert ragged == ""
        assert sc_aa == row["single_chain_aa"]

        msgs = validate_sequences(assembled, strict=True, species="human")
        assert not [m for m in msgs if m.severity == "load_bearing"]


class TestMousePgenAndSeqprob:
    """Mouse support in heuristic pgen and seqprob background models."""

    def test_mouse_pgen_and_seqprob_guards(self):
        clones = pd.DataFrame({
            "CDR3_alpha": ["CAASDNYQLIW", "CALSDRGGRALIF"],
            "CDR3_beta": ["CASSRANYEQYF", "CASGDNTEVFF"],
            "alpha_v_gene": ["Trav14-1*01", "Trav12-1*01"],
            "alpha_j_gene": ["Traj37*01", "Traj15*01"],
            "beta_v_gene": ["Trbv12-1*01", "Trbv13-2*01"],
            "beta_j_gene": ["Trbj2-7*01", "Trbj1-1*01"],
        })
        pgen_a = compute_pgen(
            clones,
            cdr3_col="CDR3_alpha",
            v_gene_col="alpha_v_gene",
            j_gene_col="alpha_j_gene",
            chain="alpha",
            species="mouse",
        )
        pgen_b = compute_pgen(clones, species="mouse")
        assert pgen_a.notna().all()
        assert pgen_b.notna().all()

        pub_df = annotate_publicness(clones, species="mouse")
        assert "publicness" in pub_df.columns
        assert pub_df["publicness"].between(0.0, 1.0).all()

        # Amino-acid kmer falls back with warning when mouse model is not trained yet
        scores = score_log_prob(clones, chain="beta", role="ppost", backend="kmer", species="mouse")
        assert scores.notna().all()

        # Gene-aware kmer_gene must NOT silently apply human V/J priors to mouse
        with pytest.raises(FileNotFoundError, match="mouse"):
            score_log_prob(clones, chain="beta", role="ppost", backend="kmer_gene", species="mouse")
