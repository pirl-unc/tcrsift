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

"""Unit tests for v3.34.0 workflow streamlining:
- Tiered `tcrsift --help` formatting (Primary workflows / Exploration & annotation / Modular steps)
- End-to-end `tcrsift prioritize --assemble` (with and without `--germline-vdj-fallback` and empty shortlists)
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift import datacache, loader, vregion
from tcrsift.cli import create_parser
from tcrsift.cli import main as cli_main
from tcrsift.prioritize import PrioritizeConfig
from tcrsift.signature_methods import SIGNATURES
from tcrsift.til_prioritize import run_til_prioritize


class TestTieredCliHelp:
    """Verify root CLI help groups workflows cleanly without breaking subcommands."""

    def test_root_help_shows_three_workflow_tiers_and_clean_usage(self):
        parser = create_parser()
        help_text = parser.format_help()
        assert "usage: tcrsift [-h] [-v] <command> ..." in help_text
        assert "{load,cells," not in help_text
        assert "Primary workflows:" in help_text
        assert "Exploration & annotation:" in help_text
        assert "Modular pipeline steps & legacy commands:" in help_text
        # Primary commands appear before legacy aliases
        assert help_text.index("prioritize") < help_text.index("til-prioritize")
        assert help_text.index("run") < help_text.index("til-select")

    def test_nested_subcommand_help_remains_standard(self):
        parser = create_parser()
        # Every subcommand still parses cleanly
        args = parser.parse_args(
            ["prioritize", "samples.yaml", "-o", "out", "--context", "solid-tumor", "--assemble"]
        )
        assert args.command == "prioritize"
        assert args.assemble is True
        assert args.single_chain_linker == "T2A"


class TestPrioritizeAssembleIntegration:
    """Verify `prioritize --assemble` produces validated full-length constructs in one pass."""

    def _mock_prioritize_adata(self) -> ad.AnnData:
        genes = sorted({g for sig in SIGNATURES.values() for g in sig.all_genes}) + [
            f"BG_{i}" for i in range(200)
        ]
        rng = np.random.default_rng(19)
        X = rng.poisson(2, size=(10, len(genes))).astype(float)
        for g in ("PRF1", "GZMB", "IFNG", "TNFRSF9"):
            X[:4, genes.index(g)] = 65.0

        obs = pd.DataFrame(
            {
                "sample": ["s1"] * 10,
                "patient_id": ["p1"] * 10,
                "CDR3_alpha": ["CAVSDGGSQGNLIF"] * 4 + ["CAASDTNTGKLTF"] * 6,
                "CDR3_beta": ["CASSIRSSYEQYF"] * 4 + ["CASSDDRVNERLFF"] * 6,
                "TRA_1_v_gene": ["TRAV12-2"] * 4 + ["TRAV1-2"] * 6,
                "TRA_1_j_gene": ["TRAJ42"] * 4 + ["TRAJ33"] * 6,
                "TRA_1_c_gene": ["TRAC"] * 10,
                "TRB_1_v_gene": ["TRBV19"] * 4 + ["TRBV6-5"] * 6,
                "TRB_1_j_gene": ["TRBJ2-7"] * 4 + ["TRBJ1-1"] * 6,
                "TRB_1_c_gene": ["TRBC2"] * 4 + ["TRBC1"] * 6,
                "TRA_1_umis": 5,
                "TRB_1_umis": 5,
                "TRA_1_reads": 50,
                "TRB_1_reads": 50,
                "CD3": 20,
                "CD4": 0,
                "CD8": 15,
            },
            index=[f"cell_{i}" for i in range(10)],
        )
        return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))

    def _seed_human_germline_cache(self, cache_dir):
        vregion.clear_reference_cache()
        specs = {
            "imgt_trav_vregion": (
                ">X01|TRAV12-2*01|Homo sapiens|F|V-REGION|1..282|282 nt|1| | | | |282+0=282| | |\n"
                "caaaaggaagaagazcagcaagagagactcagggccagtcagtgacgcagcccgatgct"  # replaced below
            ),
        }
        del specs
        # Minimal valid IMGT gapped/ungapped FASTAs for TRAV12-2, TRBV19, TRAJ42, TRBJ2-7
        entries = {
            "imgt_trav_vregion": (
                ">M1|TRAV12-2*01|Homo sapiens|F|V-REGION|1..315|315 nt|1| | | | |315+0=315| | |\n"
                "aaaggaagaagagcagcaagagagactcagggccagtcagtgacgcagcccgatgctc"
                "gcgtcactgtctctgaaggagcctctctgcagctgagatgcaagtattcctactctggga"
                "caccttatcagttatggtacgtccagtacccgcggcaggggctgcagctgctcctcaagt"
                "actattcaggagacccagtggttcaaggagtgaatggcttcgaggctgagttcagcaaga"
                "gcaactcttccttccacctgcggaaagcctctgtgcactggagcgactctgctgtgtact"
                "tctgtgctgtgagtgac\n"
            ),
            "imgt_trbv_vregion": (
                ">M2|TRBV19*01|Homo sapiens|F|V-REGION|1..285|285 nt|1| | | | |285+0=285| | |\n"
                "gatggtggaattactcaatccccaaagtacctgttcagaaaggaaggacagaatgtgacc"
                "ctgagttgtgaacagaatttgaaccacgatgccatgtactggtaccgacaggacccaggg"
                "caagggctgagattgatctactactcacagatagtaaatgactttcagaaaggagatatag"
                "ctgaagggtacagcgtctctcgggagaagaaggaatcctttcctctcactgtgacatcgg"
                "cccaaaagaacccgacagctttctatctctgtgccagtagtata\n"
            ),
            "imgt_traj_jregion": (
                ">M3|TRAJ42*01|Homo sapiens|F|J-REGION|1..54|54 nt|1| | | | |54+0=54| | |\n"
                "ggaagccaaggaaatctcatctttggaaaaggcactaaactctctgttaaacca\n"
            ),
            "imgt_trbj_jregion": (
                ">M4|TRBJ2-7*01|Homo sapiens|F|J-REGION|1..45|45 nt|1| | | | |45+0=45| | |\n"
                "tcctacgagcagtacttcgggccgggcaccaggctcacggtcaca\n"
            ),
        }
        for db_name, content in entries.items():
            d = cache_dir / db_name
            d.mkdir(parents=True, exist_ok=True)
            (d / datacache.DATABASES[db_name].filename).write_text(content)

    def test_prioritize_with_assemble_writes_constructs_and_fasta(self, tmp_path, monkeypatch):
        cache_dir = tmp_path / "cache"
        self._seed_human_germline_cache(cache_dir)
        cells = self._mock_prioritize_adata()

        sheet = tmp_path / "samples.yaml"
        sheet.write_text("samples:\n  - sample: s1\n    vdj_dir: v\n    gex_dir: g\n")
        monkeypatch.setattr(loader, "load_samples", lambda *a, **kw: cells.copy())

        out_dir = tmp_path / "out_asm"
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_dir),
                "--cache-dir",
                str(cache_dir),
                "--signatures",
                "Cytolytic",
                "--signature-cutoff",
                "legacy",
                "--signature-quantile",
                "0",
                "--no-exclude-known-mart1",
                "--assemble",
                "--germline-vdj-fallback",
            ]
        )

        constructs_csv = out_dir / "candidate_constructs.csv"
        constructs_fasta = out_dir / "candidate_constructs.fasta"
        meta = json.loads((out_dir / "prioritization.json").read_text())

        assert constructs_csv.exists()
        assert constructs_fasta.exists()
        df_c = pd.read_csv(constructs_csv)
        assert len(df_c) == 1
        assert df_c.loc[0, "CDR3_alpha"] == "CAVSDGGSQGNLIF"
        assert isinstance(df_c.loc[0, "single_chain_aa"], str)
        assert len(df_c.loc[0, "single_chain_aa"]) > 500
        assert ">" in constructs_fasta.read_text()
        assert meta["assembly"]["enabled"] is True
        assert meta["assembly"]["complete_constructs"] == 1
        assert meta["assembly"]["load_bearing_validation_errors"] == 0
        vregion.clear_reference_cache()

    def test_prioritize_invalid_single_chain_linker_fails_fast(self, tmp_path):
        cfg = PrioritizeConfig(
            sample_sheet=tmp_path / "missing.yaml",
            output_dir=tmp_path / "out",
            single_chain_linker="BAD_LINKER",
        )
        with pytest.raises(ValueError, match="--single-chain-linker"):
            run_til_prioritize(cfg)

    def test_prioritize_assemble_contig_auto_detect_and_truncated_leader_fallback(
        self, tmp_path, monkeypatch
    ):
        """Verify 5'-truncated from_contig clones receive germline M... leaders and 0 frame errors."""
        cache_dir = tmp_path / "cache"
        self._seed_human_germline_cache(cache_dir)
        cells = self._mock_prioritize_adata()
        # Provide 5'-truncated partial VDJ (plus +1 trailing nt at J->C boundary) + contig JSON
        vdj_dir = tmp_path / "s1_vdj"
        vdj_dir.mkdir(parents=True, exist_ok=True)
        partial_a_aa = "GSQGNLIFGKGTKLSVKP"  # suffix of CAVSDGGSQGNLIF + TRAJ42 FR4
        # Wait: include full CDR3_alpha inside partial_VDJ_alpha_aa as Cell Ranger does:
        partial_a_aa = "DSAVYFCAVSDGGSQGNLIFGKGTKLSVKP"
        partial_a_nt = (
            "a"  # +1 5' truncation offset
            + "gactctgctgtgtacttctgtgctgtgagtgacggaagccaaggaaatctcatctttggaaaaggcactaaactctctgttaaacca"
            + "a"  # +1 3' J->C splice nucleotide
        )
        # Contig continues with TRAC junction codon 'aac' (N) + IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS
        trac_prefix_nt = "aacatccagaacccagacccagccgtgtaccagctgagagactctaaatccagtgacaagtctgtctgcctattcaccgattttgattctcaaacaaatgtgtcacaaagtaaggattctgatgtgtatatcacagacaaaactgtgctagacatgaggtctatggacttcaagagcaacagtgctgtggcctggagcaacaaatctgactttgcatgtgcaaacgccttcaacaacagcattattccAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTAA"
        contig_seq = (partial_a_nt[:-1] + trac_prefix_nt).upper()
        (vdj_dir / "filtered_contig.fasta").write_text(f">c_alpha_1\n{contig_seq}\n")

        cells.obs["TRA_1_partial_vdj_aa"] = [partial_a_aa] * 4 + [None] * 6
        cells.obs["TRA_1_partial_vdj_nt"] = [partial_a_nt] * 4 + [None] * 6
        cells.obs["TRA_1_contig_id"] = ["c_alpha_1"] * 4 + ["c_none"] * 6

        sheet = tmp_path / "samples.yaml"
        sheet.write_text(f"samples:\n  - sample: s1\n    vdj_dir: {vdj_dir}\n    gex_dir: g\n")
        monkeypatch.setattr(loader, "load_samples", lambda *a, **kw: cells.copy())

        out_dir = tmp_path / "out_trunc"
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_dir),
                "--cache-dir",
                str(cache_dir),
                "--signatures",
                "Cytolytic",
                "--signature-cutoff",
                "legacy",
                "--signature-quantile",
                "0",
                "--no-exclude-known-mart1",
                "--assemble",
                "--germline-vdj-fallback",
                "--single-chain-linker",
                "NONE",
            ]
        )
        df_c = pd.read_csv(out_dir / "candidate_constructs.csv")
        fasta_txt = (out_dir / "candidate_constructs.fasta").read_text()
        meta = json.loads((out_dir / "prioritization.json").read_text())

        assert len(df_c) == 1
        assert str(df_c.loc[0, "full_alpha_aa"]).startswith("M")
        assert str(df_c.loc[0, "full_beta_aa"]).startswith("M")
        assert "|beta " in fasta_txt and "|alpha " in fasta_txt
        assert meta["assembly"]["complete_constructs"] == 1
        assert meta["assembly"]["load_bearing_validation_errors"] == 0

        # Rerun WITHOUT --assemble cleans up stale construct files
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_dir),
                "--signatures",
                "Cytolytic",
                "--signature-cutoff",
                "legacy",
                "--signature-quantile",
                "0",
                "--no-exclude-known-mart1",
            ]
        )
        assert not (out_dir / "candidate_constructs.csv").exists()
        assert not (out_dir / "candidate_constructs.fasta").exists()
        vregion.clear_reference_cache()

    def test_prioritize_assemble_empty_shortlist_writes_schema_headers(
        self, tmp_path, monkeypatch
    ):
        cells = self._mock_prioritize_adata()
        sheet = tmp_path / "samples.yaml"
        sheet.write_text("samples:\n  - sample: s1\n    vdj_dir: v\n    gex_dir: g\n")
        monkeypatch.setattr(loader, "load_samples", lambda *a, **kw: cells.copy())

        out_dir = tmp_path / "out_empty"
        cli_main(
            [
                "prioritize",
                str(sheet),
                "-o",
                str(out_dir),
                "--signatures",
                "Cytolytic",
                "--min-cells",
                "99",
                "--assemble",
            ]
        )
        df_c = pd.read_csv(out_dir / "candidate_constructs.csv")
        meta = json.loads((out_dir / "prioritization.json").read_text())
        assert len(df_c) == 0
        assert {"full_alpha_aa", "full_beta_aa", "single_chain_aa", "qc_warnings"} <= set(df_c.columns)
        assert (out_dir / "candidate_constructs.fasta").read_text() == ""
        assert meta["assembly"]["complete_constructs"] == 0
