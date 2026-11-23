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

"""Comprehensive tests for v3.32.0 architecture, IMGT germline rescue, and prioritization upgrades."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tcrsift
import tcrsift.vregion as vr
from tcrsift.assemble import (
    HUMAN_CONSTANT_REGIONS_AA,
    VDJ_SOURCE_5PRIME_RESCUED,
    VDJ_SOURCE_CONTIG_FULL,
    VDJ_SOURCE_VJ_RECONSTRUCTED,
    assemble_full_sequences,
    back_translate,
    enforce_contig_fidelity,
    reconstruct_germline_vdj,
    translate_dna,
)
from tcrsift.cli import create_parser, main
from tcrsift.clonotype import aggregate_clonotypes
from tcrsift.filter import filter_clonotypes, filter_clonotypes_stimulation
from tcrsift.leaders import germline_vgene_leaders
from tcrsift.loader import load_cellranger_vdj, load_sample
from tcrsift.multimer import call_multimers_for_adata
from tcrsift.prioritize import (
    PRIORITY_CONTEXT_SIGNATURES,
    PrioritizeConfig,
    select_round_robin,
)
from tcrsift.sample_sheet import Sample, SampleSheet
from tcrsift.validation import TCRsiftValidationError


def _write_synthetic_germline_cache(cache_dir: Path) -> tuple[str, str, str, str]:
    """Seed synthetic human TRAV27, TRBV5-1, TRAJ42, and TRBJ2-1 IMGT GENElect files."""
    trav_fr = (
        "TQLLEQSPQFLSIQEGENLTVYCNSSSVFSSLQWYRQEPGEGPVLLVTVVTGGEVKKLKRLTFQFGDARKDSSLHITAAQPGDTGLYL"
    )
    trav_full_aa = trav_fr + "CAG"
    trbv_fr = (
        "KAGVTQTPRYLIKTRGQQVTLSCSPISGHRSVSWYQQTPGQGLQFLFEYFSETQRNKGNFPGRFSGRQFSNSRSEMNVSTLELGDSALYL"
    )
    trbv_full_aa = trbv_fr + "CASS"

    traj_aa = "NYGGSQGNLIFGKGTKLSVKP"
    trbj_aa = "SYNEQFFGPGTRLTVL"

    trav_nt = back_translate(trav_full_aa, species="human")
    trbv_nt = back_translate(trbv_full_aa, species="human")
    # Prepend 2-nt heptamer offset and append 1-nt 3' splice overhang like real IMGT J-REGIONs!
    traj_nt = "GA" + back_translate(traj_aa, species="human") + "A"
    trbj_nt = "CT" + back_translate(trbj_aa, species="human") + "G"

    files = {
        ("imgt_trav_vregion", "TRAV_V-REGION.html"): (
            f"<pre>\n>X1|TRAV27*01|Homo sapiens|F|V-REGION|1..{len(trav_nt)}|\n{trav_nt}\n</pre>\n"
        ),
        ("imgt_trbv_vregion", "TRBV_V-REGION.html"): (
            f"<pre>\n>X2|TRBV5-1*01|Homo sapiens|F|V-REGION|1..{len(trbv_nt)}|\n{trbv_nt}\n</pre>\n"
        ),
        ("imgt_traj_jregion", "TRAJ_J-REGION.html"): (
            f"<pre>\n>X3|TRAJ42*01|Homo sapiens|F|J-REGION|1..{len(traj_nt)}|\n{traj_nt}\n</pre>\n"
        ),
        ("imgt_trbj_jregion", "TRBJ_J-REGION.html"): (
            f"<pre>\n>X4|TRBJ2-1*01|Homo sapiens|F|J-REGION|1..{len(trbj_nt)}|\n{trbj_nt}\n</pre>\n"
        ),
    }
    for (db, fname), payload in files.items():
        d = cache_dir / db
        d.mkdir(parents=True, exist_ok=True)
        (d / fname).write_text(payload)
    vr.clear_reference_cache()
    return trav_fr, trbv_fr, traj_aa, trbj_aa


class TestGermlineVDJReconstruction:
    """IMGT V+J germline parsing, 5'-partial contig rescue, and bare V+CDR3+J assembly."""

    def test_jregion_parser_handles_frame_offset_and_splice_overhang(self, tmp_path):
        _trav_fr, _trbv_fr, traj_aa, trbj_aa = _write_synthetic_germline_cache(tmp_path)

        j_alpha = vr.germline_jregion_sequences("TRAJ42", species="human", data_dir=tmp_path)
        assert len(j_alpha) == 1
        allele, func, aa, nt = j_alpha[0]
        assert (allele, func, aa) == ("01", "F", traj_aa)
        assert translate_dna(nt)[0] == traj_aa

        # Zero-padded ImmunoSEQ/VDJdb synonym lookup works transparently.
        j_beta = vr.germline_jregion_sequences("TRBJ02-01*01", species="human", data_dir=tmp_path)
        assert len(j_beta) == 1 and j_beta[0][2] == trbj_aa
        v_beta = vr.germline_vregion_sequences("TRBV05-01", species="human", data_dir=tmp_path)
        assert len(v_beta) == 1

    def test_bare_vj_reconstruction_and_provenance(self, tmp_path):
        trav_fr, trbv_fr, _traj_aa, _trbj_aa = _write_synthetic_germline_cache(tmp_path)
        clones = pd.DataFrame(
            [
                {
                    "CDR3ab": "CAAGAGSQGNLIF_CASSPFGGGINEQFF",
                    "CDR3_alpha": "CAAGAGSQGNLIF",
                    "CDR3_beta": "CASSPFGGGINEQFF",
                    "alpha_v_gene": "TRAV27",
                    "alpha_j_gene": "TRAJ42",
                    "beta_v_gene": "TRBV05-01",
                    "beta_j_gene": "TRBJ02-01",
                }
            ]
        )
        rec = reconstruct_germline_vdj(clones, species="human", data_dir=tmp_path)
        assert rec["alpha_vdj_source"].iloc[0] == VDJ_SOURCE_VJ_RECONSTRUCTED
        assert rec["beta_vdj_source"].iloc[0] == VDJ_SOURCE_VJ_RECONSTRUCTED
        assert rec["VDJ_alpha_aa"].iloc[0] == trav_fr + "CAAGAGSQGNLIF" + "GKGTKLSVKP"
        assert rec["VDJ_beta_aa"].iloc[0] == trbv_fr + "CASSPFGGGINEQFF" + "GPGTRLTVL"
        assert translate_dna(rec["VDJ_alpha_nt"].iloc[0].upper())[0] == rec["VDJ_alpha_aa"].iloc[0]
        assert translate_dna(rec["VDJ_beta_nt"].iloc[0].upper())[0] == rec["VDJ_beta_aa"].iloc[0]
        assert "alpha_vdj_germline_reconstructed" in rec["qc_warnings"].iloc[0]
        assert "beta_vdj_germline_reconstructed" in rec["qc_warnings"].iloc[0]

        assembled = assemble_full_sequences(
            clones,
            germline_vdj_fallback=True,
            data_dir=tmp_path,
            verbose=False,
            show_progress=False,
        )
        assert assembled["single_chain_aa"].notna().all()
        assert assembled["beta_c_gene_canonical"].iloc[0] == "TRBC2"
        assert not bool(assembled["construct_contig_verified"].iloc[0])
        assert "alpha_vdj_germline_reconstructed" in assembled["qc_warnings"].iloc[0]
        with pytest.raises(TCRsiftValidationError):
            enforce_contig_fidelity(assembled, allow_canonical_fallback=False)
        msg = enforce_contig_fidelity(assembled, allow_canonical_fallback=True)
        assert "1/1 constructs assembled with NO contig verification" in msg

    def test_noncanonical_j_motifs_and_truncated_bracketed_v_alleles(self, tmp_path):
        """Real IMGT exceptions: TRAJ35*01 (CGSG), TRAJ16*01 (FARG), and 5'-partial [F]/(F) V-alleles."""
        trav_fr, trbv_fr, _traj_aa, _trbj_aa = _write_synthetic_germline_cache(tmp_path)
        # Add TRAJ35*01 (CGSG) and TRAJ16*01 (FARG) to TRAJ cache, plus a 5'-truncated [F] allele on TRAV27!
        traj35_aa = "IGFGNVLHCGSGTQVIVLP"
        traj16_aa = "FSDGQKLLFARGTMLKVDL"
        traj35_nt = "A" + back_translate(traj35_aa, species="human") + "C"
        traj16_nt = "GA" + back_translate(traj16_aa, species="human") + "A"
        (tmp_path / "imgt_traj_jregion" / "TRAJ_J-REGION.html").write_text(
            f"<pre>\n"
            f">J35|TRAJ35*01|Homo sapiens|F|J-REGION|1..{len(traj35_nt)}|\n{traj35_nt}\n"
            f">J16|TRAJ16*01|Homo sapiens|F|J-REGION|1..{len(traj16_nt)}|\n{traj16_nt}\n"
            f"</pre>\n"
        )
        # TRAV27*02 is a 58-aa 5'-truncated [F] allele with a longer CDR3 germline prefix ("CAGAG")!
        trunc_trav_aa = trav_fr[35:] + "CAGAG"
        trunc_trav_nt = back_translate(trunc_trav_aa, species="human")
        full_trav_nt = back_translate(trav_fr + "CAG", species="human")
        (tmp_path / "imgt_trav_vregion" / "TRAV_V-REGION.html").write_text(
            f"<pre>\n"
            f">V1|TRAV27*01|Homo sapiens|F|V-REGION|1..{len(full_trav_nt)}|\n{full_trav_nt}\n"
            f">V2|TRAV27*02|Homo sapiens|[F]|V-REGION|1..{len(trunc_trav_nt)}|\n{trunc_trav_nt}\n"
            f"</pre>\n"
        )
        vr.clear_reference_cache()

        clones = pd.DataFrame(
            [
                {
                    "CDR3ab": "CAGIGFGNVLHC_CASSPFGGGINEQFF",
                    "CDR3_alpha": "CAGIGFGNVLHC",
                    "CDR3_beta": "CASSPFGGGINEQFF",
                    "alpha_v_gene": "TRAV27",
                    "alpha_j_gene": "TRAJ35",
                    "beta_v_gene": "TRBV5-1",
                    "beta_j_gene": "TRBJ02-01",
                },
                {
                    "CDR3ab": "CAGARFSDGQKLLF_CASSPFGGGINEQFF",
                    "CDR3_alpha": "CAGARFSDGQKLLF",
                    "CDR3_beta": "CASSPFGGGINEQFF",
                    "alpha_v_gene": "TRAV27",
                    "alpha_j_gene": "TRAJ016",
                    "beta_v_gene": "TRBV5-1",
                    "beta_j_gene": "TRBJ2-1",
                },
            ]
        )
        rec = reconstruct_germline_vdj(clones, species="human", data_dir=tmp_path)
        # Must pick full-length TRAV27*01 (len 89 FR), NOT 5'-truncated TRAV27*02 [F]!
        assert rec["VDJ_alpha_aa"].iloc[0] == trav_fr + "CAGIGFGNVLHC" + "GSGTQVIVLP"
        assert rec["VDJ_alpha_aa"].iloc[1] == trav_fr + "CAGARFSDGQKLLF" + "ARGTMLKVDL"
        assert rec["alpha_vdj_germline_allele"].iloc[0] == "TRAV27*01;TRAJ35*01"
        assert rec["alpha_vdj_germline_allele"].iloc[1] == "TRAV27*01;TRAJ16*01"

        # Also test ragged 5' nucleotide preservation on partial_VDJ_beta_nt!
        partial_beta_aa = trbv_fr[25:] + "CASSPFGGGINEQFFGPGTRLTVL"
        synonymous_donor_nt = "TTT" + back_translate(partial_beta_aa[1:], species="human")
        if partial_beta_aa[0] == "F":
            donor_codons = synonymous_donor_nt
        else:
            # Flip last Leucine ("L") codon to synonymous "TTA" so it differs from canonical back_translate ("CTG")!
            assert partial_beta_aa.endswith("L")
            donor_codons = back_translate(partial_beta_aa[:-1], species="human") + "TTA"
        ragged_partial_nt = "AC" + donor_codons + "G"  # 2 leading ragged bases + 1 trailing overhang
        partial_df = pd.DataFrame(
            [
                {
                    "CDR3ab": "CAAGAGSQGNLIF_CASSPFGGGINEQFF",
                    "CDR3_alpha": "CAAGAGSQGNLIF",
                    "CDR3_beta": "CASSPFGGGINEQFF",
                    "alpha_v_gene": "TRAV27",
                    "alpha_j_gene": "TRAJ35",
                    "beta_v_gene": "TRBV5-1",
                    "beta_j_gene": "TRBJ2-1",
                    "partial_VDJ_beta_aa": partial_beta_aa,
                    "partial_VDJ_beta_nt": ragged_partial_nt,
                }
            ]
        )
        rec_p = reconstruct_germline_vdj(partial_df, species="human", data_dir=tmp_path)
        assert rec_p["beta_vdj_source"].iloc[0] == VDJ_SOURCE_5PRIME_RESCUED
        assert rec_p["VDJ_beta_nt"].iloc[0].upper().endswith("TTA")

    def test_five_prime_truncated_json_contigs_rescue_missing_fr1_and_preserve_donor_snps(
        self, tmp_path
    ):
        """Simulate Fischer 2021 GSE171037 5'-truncated contigs (start_codon_pos=None)."""
        trav_fr, trbv_fr, _traj_aa, _trbj_aa = _write_synthetic_germline_cache(tmp_path)

        # Introduce a real donor framework SNP in FR3 of alpha (pos -10: A->V) that MUST survive rescue!
        donor_trav_fr = trav_fr[:-10] + "V" + trav_fr[-9:]
        alpha_vdj_aa = donor_trav_fr + "CAAGAGSQGNLIFGKGTKLSVKP"
        beta_leader_aa = germline_vgene_leaders("TRBV5-1", species="human")[0][2]
        beta_vdj_aa = trbv_fr + "CASSPFGGGINEQFFGPGTRLTVL"

        # Alpha contig is 5'-truncated 22 codons into V-REGION ( plus 2 ragged 5' bases )!
        trunc_codons = 22
        alpha_partial_aa = alpha_vdj_aa[trunc_codons:]
        alpha_c_tail_nt = back_translate(
            "N" + HUMAN_CONSTANT_REGIONS_AA["TRAC"][:18], species="human"
        )
        alpha_seq = (
            "GT"  # 2 ragged non-frame bases at 5' read start
            + back_translate(alpha_partial_aa, species="human")
            + alpha_c_tail_nt
        )
        a_j_end = 2 + 3 * len(alpha_partial_aa) + 1  # +1 J->C overhang

        # Beta contig is full-length with leader start codon.
        utr5 = "TTCCCAGGAGGACCAAGCCC"
        beta_seq = (
            utr5
            + back_translate(beta_leader_aa + beta_vdj_aa, species="human")
            + back_translate("E" + HUMAN_CONSTANT_REGIONS_AA["TRBC2"][:18], species="human")
        )
        b_j_end = len(utr5) + 3 * len(beta_leader_aa + beta_vdj_aa) + 1

        sdir = tmp_path / "fischer_sample"
        sdir.mkdir()
        json_records = [
            {
                "contig_name": "BC1-1_contig_1",
                "barcode": "BC1-1",
                "is_cell": True,
                "high_confidence": True,
                "full_length": False,
                "productive": None,
                "cdr3": "CAAGAGSQGNLIF",
                "start_codon_pos": None,
                "sequence": alpha_seq,
                "annotations": [
                    {
                        "feature": {"region_type": "V-REGION", "gene_name": "TRAV27"},
                        "contig_match_start": 0,
                        "contig_match_end": a_j_end - 42,
                        "annotation_match_start": 3 * trunc_codons - 2,
                    },
                    {
                        "feature": {"region_type": "J-REGION", "gene_name": "TRAJ42"},
                        "contig_match_start": a_j_end - 42,
                        "contig_match_end": a_j_end,
                    },
                    {
                        "feature": {"region_type": "C-REGION", "gene_name": "TRAC"},
                        "contig_match_start": a_j_end,
                        "contig_match_end": len(alpha_seq),
                    },
                ],
            },
            {
                "contig_name": "BC1-1_contig_2",
                "barcode": "BC1-1",
                "is_cell": True,
                "high_confidence": True,
                "full_length": True,
                "productive": True,
                "cdr3": "CASSPFGGGINEQFF",
                "start_codon_pos": len(utr5),
                "sequence": beta_seq,
                "annotations": [
                    {
                        "feature": {"region_type": "L-REGION+V-REGION", "gene_name": "TRBV5-1"},
                        "contig_match_start": len(utr5),
                        "contig_match_end": b_j_end - 42,
                    },
                    {
                        "feature": {"region_type": "J-REGION", "gene_name": "TRBJ2-1"},
                        "contig_match_start": b_j_end - 42,
                        "contig_match_end": b_j_end,
                    },
                    {
                        "feature": {"region_type": "C-REGION", "gene_name": "TRBC2"},
                        "contig_match_start": b_j_end,
                        "contig_match_end": len(beta_seq),
                    },
                ],
            },
        ]
        with gzip.open(sdir / "all_contig_annotations.json.gz", "wt", encoding="utf-8") as fh:
            json.dump(json_records, fh)

        pd.DataFrame(
            {
                "barcode": ["BC1-1", "BC1-1"],
                "is_cell": [True, True],
                "contig_id": ["BC1-1_contig_1", "BC1-1_contig_2"],
                "high_confidence": [True, True],
                "chain": ["TRA", "TRB"],
                "v_gene": ["TRAV27", "TRBV5-1"],
                "d_gene": ["None", "TRBD1"],
                "j_gene": ["TRAJ42", "TRBJ2-1"],
                "c_gene": ["TRAC", "TRBC2"],
                "full_length": [False, True],
                "productive": ["None", "True"],
                "cdr3": ["CAAGAGSQGNLIF", "CASSPFGGGINEQFF"],
                "cdr3_nt": [
                    back_translate("CAAGAGSQGNLIF", species="human"),
                    back_translate("CASSPFGGGINEQFF", species="human"),
                ],
                "reads": [400, 500],
                "umis": [12, 15],
            }
        ).to_csv(sdir / "filtered_contig_annotations.csv.gz", index=False)

        vdj_df = load_cellranger_vdj(sdir, "fischer_sample", verbose=False)
        tra_row = vdj_df[vdj_df["chain"] == "TRA"].iloc[0]
        assert pd.isna(tra_row["vdj_aa"])
        assert tra_row["partial_vdj_aa"] == alpha_partial_aa

        adata = load_sample(Sample(sample="fischer_sample", vdj_dir=str(sdir)))
        clones = aggregate_clonotypes(adata, verbose=False, show_progress=False)
        assert "VDJ_alpha_aa" not in clones.columns or pd.isna(clones["VDJ_alpha_aa"].iloc[0])
        assert clones["partial_VDJ_alpha_aa"].iloc[0] == alpha_partial_aa

        # Default assembly (germline_vdj_fallback=False) leaves alpha unassembled (never a VDJ-less leader+C stub).
        asm_default = assemble_full_sequences(
            clones, contigs_dir=tmp_path, verbose=False, show_progress=False
        )
        assert (
            "full_alpha_aa" not in asm_default.columns
            or pd.isna(asm_default["full_alpha_aa"].iloc[0])
        )

        # With germline_vdj_fallback=True, missing 5' FR1 is rescued AND donor FR3 SNP + donor J->C junction survive!
        asm_rescued = assemble_full_sequences(
            clones,
            contigs_dir=tmp_path,
            germline_vdj_fallback=True,
            data_dir=tmp_path,
            verbose=False,
            show_progress=False,
        )
        assert asm_rescued["alpha_vdj_source"].iloc[0] == VDJ_SOURCE_5PRIME_RESCUED
        assert asm_rescued["beta_vdj_source"].iloc[0] == VDJ_SOURCE_CONTIG_FULL
        assert asm_rescued["vdj_alpha_aa"].iloc[0] == alpha_vdj_aa  # exact including donor V SNP!
        assert asm_rescued["alpha_junction_residue_source"].iloc[0] == "contig"
        assert asm_rescued["beta_junction_residue_source"].iloc[0] == "contig"
        assert asm_rescued["single_chain_aa"].notna().all()
        assert not bool(asm_rescued["construct_contig_verified"].iloc[0])


class TestPrioritizeMethodologicalUpgrades:
    """Priority-tiered round-robin, safe strict_ab exclusions, and --call-multimers."""

    def test_priority_tiered_round_robin_ranks_tumor_signatures_before_backfill(self):
        clones = pd.DataFrame(
            {
                "donor": ["D1", "D1", "D1"],
                "CDR3ab": ["CA1_CB1", "CA2_CB2", "CA3_CB3"],
                "eligible_for_review": [True, True, True],
            }
        )
        scores = pd.DataFrame(
            {
                "donor": ["D1", "D1", "D1"],
                "sample": ["S1", "S1", "S1"],
                "lineage": ["cd8", "cd8", "cd8"],
                "CDR3ab": ["CA1_CB1", "CA2_CB2", "CA3_CB3"],
                "cells": [10, 10, 10],
                "frequency": [0.05, 0.04, 0.03],
                "signature_Cytolytic": [3.0, 0.1, 0.2],
                "signature_Cytolytic_percentile": [1.0, 0.3, 0.6],
                "signature_TumorReactive": [0.1, 2.5, 0.2],
                "signature_TumorReactive_percentile": [0.3, 1.0, 0.6],
            }
        )
        # Flat mode gives Cytolytic turn 1 (listed first).
        flat = select_round_robin(
            clones,
            scores,
            ["Cytolytic", "TumorReactive"],
            max_clones=1,
            quantile=0.8,
            cutoff="legacy",
        )
        assert flat.loc[flat["selection_rank"] == 1, "CDR3ab"].iloc[0] == "CA1_CB1"

        # Tiered mode gives TumorReactive turn 1 and backfills Cytolytic into rank 2 when budget=2!
        tiered_1 = select_round_robin(
            clones,
            scores,
            ["Cytolytic", "TumorReactive"],
            max_clones=1,
            quantile=0.8,
            cutoff="legacy",
            priority_signatures=PRIORITY_CONTEXT_SIGNATURES["solid-tumor"],
        )
        assert tiered_1.loc[tiered_1["selection_rank"] == 1, "CDR3ab"].iloc[0] == "CA2_CB2"

        tiered_2 = select_round_robin(
            clones,
            scores,
            ["Cytolytic", "TumorReactive"],
            max_clones=2,
            quantile=0.8,
            cutoff="legacy",
            priority_signatures=PRIORITY_CONTEXT_SIGNATURES["solid-tumor"],
        )
        ranked = tiered_2[tiered_2["selected_for_review"]]["CDR3ab"].tolist()
        assert ranked == ["CA2_CB2", "CA1_CB1"]

    def test_database_match_defaults_to_strict_ab_preserving_beta_partial_tumor_clones(
        self, tmp_path, monkeypatch
    ):
        # Default CLI parse: database_match == strict_ab for paired α/β prioritize runs.
        args_default = create_parser().parse_args(
            ["til-prioritize", "samples.yaml", "-o", str(tmp_path), "--exclude-known-viral"]
        )
        cfg_default = PrioritizeConfig.from_namespace(args_default)
        assert cfg_default.database_match == "strict_ab"

        # Explicit --database-match b_only overrides cleanly.
        args_explicit = create_parser().parse_args(
            [
                "til-prioritize",
                "samples.yaml",
                "-o",
                str(tmp_path),
                "--database-match",
                "b_only",
            ]
        )
        cfg_explicit = PrioritizeConfig.from_namespace(args_explicit)
        assert cfg_explicit.database_match == "b_only"

        # End-to-end check: under default strict_ab, a tumor clone sharing only CDR3b
        # with a paired viral DB entry is neither mis-flagged viral nor discarded.
        vdjdb = tmp_path / "vdjdb.tsv"
        pd.DataFrame(
            [
                {
                    "gene": "TRA",
                    "cdr3": "CAVDIFFERENTF",
                    "species": "HomoSapiens",
                    "antigen.epitope": "GILGFVFTL",
                    "antigen.gene": "M1",
                    "antigen.species": "InfluenzaA",
                    "complex.id": 1,
                },
                {
                    "gene": "TRB",
                    "cdr3": "CASSLGQAYEQYF",  # shared public beta with our tumor clone!
                    "species": "HomoSapiens",
                    "antigen.epitope": "GILGFVFTL",
                    "antigen.gene": "M1",
                    "antigen.species": "InfluenzaA",
                    "complex.id": 1,
                },
            ]
        ).to_csv(vdjdb, sep="\t", index=False)

        from tests.test_til_prioritize import _mock_samples, _til_cells

        cells = _til_cells()[:6].copy()
        sheet = _mock_samples(tmp_path, monkeypatch, cells)
        main(
            [
                "til-prioritize",
                sheet,
                "-o",
                str(tmp_path / "out"),
                "--vdjdb",
                str(vdjdb),
                "--exclude-known-viral",
                "--signature-cutoff",
                "legacy",
                "--signatures",
                "Cytolytic",
                "--signature-quantile",
                "0",
                "--min-signature-score",
                "-100",
            ]
        )
        audit = pd.read_csv(tmp_path / "out" / "all_scored_clones.csv")
        cands = pd.read_csv(tmp_path / "out" / "candidate_clones.csv")
        target = audit[audit["CDR3_beta"] == "CASSLGQAYEQYF"].iloc[0]
        assert bool(target["is_viral"]) is False
        assert bool(target["known_viral_match"]) is False
        assert len(cands) == 2

    def test_call_multimers_for_adata_handles_shared_barcodes_across_samples(
        self, tmp_path, monkeypatch
    ):
        from tcrsift import multimer as mm_mod

        # Two samples reuse identical 10x barcodes ("BC1", "BC2", "BC3"), binding DIFFERENT multimers!
        obs = pd.DataFrame(
            {
                "sample": ["S1", "S1", "S1", "S2", "S2", "S2"],
                "donor": ["D1", "D1", "D1", "D2", "D2", "D2"],
                "CDR3ab": ["CA1_CB1", "CA1_CB1", "CA2_CB2", "CA1_CB1", "CA1_CB1", "CA3_CB3"],
            },
            index=["BC1", "BC2", "BC3", "BC1", "BC2", "BC3"],
        )
        import warnings

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Observation names are not unique", category=UserWarning)
            adata = ad.AnnData(X=np.zeros((6, 1)), obs=obs)
        sheet = SampleSheet(
            samples=[
                Sample(sample="S1", gex_dir=str(tmp_path / "gex1")),
                Sample(sample="S2", gex_dir=str(tmp_path / "gex2")),
            ]
        )

        def fake_extract(path):
            if str(path).endswith("gex1"):
                df = pd.DataFrame(
                    {
                        "A0201_GILGFVFTL": [80.0, 95.0, 1.0],
                        "A0201_NLVPMVATV": [1.0, 2.0, 0.0],
                        "NEGATIVE_CONTROL": [1.0, 1.0, 1.0],
                    },
                    index=["BC1", "BC2", "BC3"],
                )
            else:
                df = pd.DataFrame(
                    {
                        "A0201_GILGFVFTL": [0.0, 1.0, 1.0],
                        "A0201_NLVPMVATV": [120.0, 110.0, 2.0],
                        "NEGATIVE_CONTROL": [1.0, 1.0, 1.0],
                    },
                    index=["BC1", "BC2", "BC3"],
                )
            df.attrs["feature_types"] = {c: "Antigen Capture" for c in df.columns}
            return df

        monkeypatch.setattr(mm_mod, "extract_feature_barcodes", fake_extract)
        res = call_multimers_for_adata(
            adata, sheet, clone_col="CDR3ab", group_col="donor", min_umi=10.0, min_nc_ratio=5.0
        )
        d1 = res[(res["donor"] == "D1") & (res["CDR3ab"] == "CA1_CB1")].iloc[0]
        d2 = res[(res["donor"] == "D2") & (res["CDR3ab"] == "CA1_CB1")].iloc[0]
        assert d1["multimer_call"] == "A0201_GILGFVFTL" and bool(d1["multimer_confident"])
        assert d2["multimer_call"] == "A0201_NLVPMVATV" and bool(d2["multimer_confident"])


class TestStimulationFilteringMethod:
    """First-class method='stimulation' for both expansion and 16 h AIM upregulation."""

    def test_stimulation_filter_keeps_expanded_and_16h_aim_upregulated_clones(self):
        clones = pd.DataFrame(
            {
                "CDR3ab": ["CA_EXP_CB", "CA_AIM16H_CB", "CA_BYSTANDER_CB"],
                "CDR3_alpha": ["CAEXP", "CAAIM", "CABYS"],
                "CDR3_beta": ["CBEXP", "CBAIM", "CBBYS"],
                "cell_count": [25, 4, 30],
                "max_frequency": [0.05, 0.008, 0.06],
                "Tcell_type_consensus": ["CD8+ T cell"] * 3,
                # Clone 1 expanded in 10d culture; Clone 2 unexpanded at 16h (FDR=1.0) but strong delta_AIM;
                # Clone 3 abundant viral bystander present equally in control with zero AIM shift!
                "stim_fisher_fdr": [0.0005, 1.0, 0.95],
                "stim_log2_fc": [3.2, 0.1, 0.05],
                "delta_signature_AIM": [0.1, 1.2, -0.05],
            }
        )
        out = filter_clonotypes(
            clones,
            method="stimulation",
            stim_fdr_threshold=0.05,
            stim_log2_fc_threshold=1.0,
            stim_delta_threshold=0.5,
            verbose=False,
        )
        assert set(out["CDR3ab"]) == {"CA_EXP_CB", "CA_AIM16H_CB"}
        tiers = dict(zip(out["CDR3ab"], out["tier"]))
        assert tiers["CA_EXP_CB"] == "tier1"
        assert tiers["CA_AIM16H_CB"] == "tier1"

    def test_stimulation_filter_errors_clearly_without_control_columns(self):
        clones = pd.DataFrame(
            {
                "CDR3ab": ["CA1_CB1"],
                "CDR3_alpha": ["CA1"],
                "CDR3_beta": ["CB1"],
                "cell_count": [5],
                "max_frequency": [0.02],
            }
        )
        with pytest.raises(TCRsiftValidationError, match="stimulation"):
            filter_clonotypes_stimulation(clones)
        with pytest.raises(TCRsiftValidationError, match="stim_fdr_threshold"):
            filter_clonotypes_stimulation(clones, stim_fdr_threshold=0.0)
        with pytest.raises(TCRsiftValidationError, match="stim_delta_threshold"):
            filter_clonotypes_stimulation(clones, stim_delta_threshold=-0.1)


def test_core_import_and_kmer_pgen_without_torch_or_tcrpeg():
    """Verify tcrsift imports and scores default k-mer Pgen/Ppost when torch & tcrpeg are absent."""
    import subprocess
    import sys

    code = (
        "import importlib.abc, sys\n"
        "class _Block(importlib.abc.MetaPathFinder):\n"
        "    def find_spec(self, fullname, path, target=None):\n"
        "        if fullname.split('.')[0] in {'torch', 'tcrpeg'}:\n"
        "            return None\n"
        "        return None\n"
        "sys.modules['torch'] = None\n"
        "sys.modules['tcrpeg'] = None\n"
        "import pandas as pd, tcrsift\n"
        "from tcrsift.annotate_tcrs import add_pgen_ppost\n"
        "from tcrsift.seqprob import TCRpegProbabilityModel\n"
        "assert not TCRpegProbabilityModel.available()\n"
        "df = pd.DataFrame({'CDR3_alpha': ['CAVSDGGSQGNLIF'], 'CDR3_beta': ['CASSLGQAYEQYF']})\n"
        "out = add_pgen_ppost(df, backend='kmer', species='human')\n"
        "assert out['pgen_alpha'].notna().all() and out['ppost_beta'].notna().all()\n"
    )
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert res.returncode == 0, res.stderr


def test_public_reexports_intact():
    assert callable(tcrsift.reconstruct_germline_vdj)
    assert callable(tcrsift.filter_clonotypes_stimulation)
    assert callable(tcrsift.call_multimers_for_adata)
    assert "solid-tumor" in tcrsift.PRIORITY_CONTEXT_SIGNATURES
