"""Tests for v3.40.0 default germline VDJ fallback with automatic truncation detection."""

from __future__ import annotations

import logging

import pandas as pd

from tcrsift import cli
from tcrsift.assemble import (
    VDJ_SOURCE_5PRIME_RESCUED,
    VDJ_SOURCE_CONTIG_FULL,
    VDJ_SOURCE_VJ_RECONSTRUCTED,
    assemble_full_sequences,
    reconstruct_germline_vdj,
)
from tcrsift.assemble import vregion as vr
from tcrsift.assemble.codons import back_translate
from tcrsift.config import AssembleConfig, TCRsiftConfig
from tcrsift.prioritize import PrioritizeConfig
from tests.test_v3_32_features import _write_synthetic_germline_cache


def test_defaults_enable_germline_vdj_fallback(tmp_path):
    """Config dataclasses and CLI parsers enable germline_vdj_fallback by default."""
    assert AssembleConfig().germline_vdj_fallback is True
    assert (
        PrioritizeConfig(
            sample_sheet=tmp_path / "samples.yaml",
            output_dir=tmp_path / "out",
        ).germline_vdj_fallback
        is True
    )

    parser = cli.create_parser()

    # tcrsift assemble: default True, supports --no-germline-vdj-fallback & --germline-vdj-fallback
    a_def = parser.parse_args(["assemble", "-i", "in.csv", "-o", "out.csv"])
    assert a_def.germline_vdj_fallback is True
    a_no = parser.parse_args(
        ["assemble", "-i", "in.csv", "-o", "out.csv", "--no-germline-vdj-fallback"]
    )
    assert a_no.germline_vdj_fallback is False
    a_yes = parser.parse_args(
        ["assemble", "-i", "in.csv", "-o", "out.csv", "--germline-vdj-fallback"]
    )
    assert a_yes.germline_vdj_fallback is True

    # tcrsift prioritize / til-prioritize: default True, supports --no-germline-vdj-fallback
    p_def = parser.parse_args(["prioritize", "s.yaml", "-o", "out"])
    assert p_def.germline_vdj_fallback is True
    p_no = parser.parse_args(
        ["prioritize", "s.yaml", "-o", "out", "--no-germline-vdj-fallback"]
    )
    assert p_no.germline_vdj_fallback is False

    # tcrsift run: default None so YAML AssembleConfig default (True) is preserved,
    # and --no-germline-vdj-fallback overrides YAML to False.
    r_def = parser.parse_args(["run", "-s", "s.yaml", "-o", "out"])
    assert r_def.germline_vdj_fallback is None
    assert TCRsiftConfig().merge_with_args(r_def).assemble.germline_vdj_fallback is True

    r_no = parser.parse_args(
        ["run", "-s", "s.yaml", "-o", "out", "--no-germline-vdj-fallback"]
    )
    assert r_no.germline_vdj_fallback is False
    assert TCRsiftConfig().merge_with_args(r_no).assemble.germline_vdj_fallback is False

    # tcrsift report selected: default None, accepts both positive and negative flags
    rs_def = parser.parse_args(
        ["report", "selected", "--selection", "s.csv", "--clonotypes", "c.csv", "-o", "out"]
    )
    assert rs_def.germline_vdj_fallback is None
    rs_no = parser.parse_args(
        [
            "report",
            "selected",
            "--selection",
            "s.csv",
            "--clonotypes",
            "c.csv",
            "-o",
            "out",
            "--no-germline-vdj-fallback",
        ]
    )
    assert rs_no.germline_vdj_fallback is False


def test_default_assembly_preserves_intact_donor_vdj_and_does_not_fabricate_nt(tmp_path):
    """Intact donor VDJ sequences are untouched by default fallback (even AA-only inputs)."""
    trav_fr, trbv_fr, _traj_aa, _trbj_aa = _write_synthetic_germline_cache(tmp_path)
    # Introduce a donor polymorphism so we verify germline *01 never overwrites intact VDJ
    donor_alpha = trav_fr[:-8] + "W" + trav_fr[-7:] + "CAAGAGSQGNLIFGKGTKLSVKP"
    donor_beta = trbv_fr[:-6] + "Y" + trbv_fr[-5:] + "CASSPFGGGINEQFFGPGTRLTVL"

    # AA-only intact input (no VDJ_*_nt column supplied)
    clones = pd.DataFrame(
        [
            {
                "CDR3ab": "CAAGAGSQGNLIF_CASSPFGGGINEQFF",
                "CDR3_alpha": "CAAGAGSQGNLIF",
                "CDR3_beta": "CASSPFGGGINEQFF",
                "alpha_v_gene": "TRAV27",
                "alpha_j_gene": "TRAJ42",
                "beta_v_gene": "TRBV5-1",
                "beta_j_gene": "TRBJ2-1",
                "VDJ_alpha_aa": donor_alpha,
                "VDJ_beta_aa": donor_beta,
            }
        ]
    )

    rec = reconstruct_germline_vdj(clones, species="human", data_dir=tmp_path)
    assert rec["alpha_vdj_source"].iloc[0] == VDJ_SOURCE_CONTIG_FULL
    assert rec["beta_vdj_source"].iloc[0] == VDJ_SOURCE_CONTIG_FULL
    assert "VDJ_alpha_nt" not in rec.columns
    assert "VDJ_beta_nt" not in rec.columns
    assert rec["VDJ_alpha_aa"].iloc[0] == donor_alpha
    assert rec["VDJ_beta_aa"].iloc[0] == donor_beta
    assert rec["qc_warnings"].iloc[0] == []

    asm = assemble_full_sequences(
        clones, data_dir=tmp_path, verbose=False, show_progress=False
    )
    assert asm["vdj_alpha_aa"].iloc[0] == donor_alpha
    assert asm["vdj_beta_aa"].iloc[0] == donor_beta
    assert asm["alpha_vdj_source"].iloc[0] == VDJ_SOURCE_CONTIG_FULL
    assert asm["beta_vdj_source"].iloc[0] == VDJ_SOURCE_CONTIG_FULL
    assert pd.isna(asm["full_alpha_nt_contig"].iloc[0])
    assert pd.isna(asm["full_beta_nt_contig"].iloc[0])
    assert asm["single_chain_aa"].notna().all()


def test_default_assembly_auto_rescues_truncated_and_bare_vj_and_warns_if_uncached(
    tmp_path, caplog
):
    """Default assembly rescues truncated/bare chains when cached, and warns cleanly when uncached."""
    cache = tmp_path / "cache"
    trav_fr, trbv_fr, _traj_aa, _trbj_aa = _write_synthetic_germline_cache(cache)

    partial_alpha = trav_fr[20:] + "CAAGAGSQGNLIFGKGTKLSVKP"
    clones = pd.DataFrame(
        [
            {
                "CDR3ab": "CAAGAGSQGNLIF_CASSPFGGGINEQFF",
                "CDR3_alpha": "CAAGAGSQGNLIF",
                "CDR3_beta": "CASSPFGGGINEQFF",
                "alpha_v_gene": "TRAV27",
                "alpha_j_gene": "TRAJ42",
                "beta_v_gene": "TRBV5-1",
                "beta_j_gene": "TRBJ2-1",
                "partial_VDJ_alpha_aa": partial_alpha,
                "partial_VDJ_alpha_nt": back_translate(partial_alpha, species="human").lower(),
            }
        ]
    )

    # 1. With cache present: default rescues alpha (5'-partial) and beta (bare V+J) automatically
    asm = assemble_full_sequences(
        clones, data_dir=cache, verbose=False, show_progress=False
    )
    assert asm["alpha_vdj_source"].iloc[0] == VDJ_SOURCE_5PRIME_RESCUED
    assert asm["beta_vdj_source"].iloc[0] == VDJ_SOURCE_VJ_RECONSTRUCTED
    assert asm["single_chain_aa"].notna().all()
    assert not bool(asm["construct_contig_verified"].iloc[0])
    assert "alpha_vdj_germline_reconstructed" in asm["qc_warnings"].iloc[0]
    assert "beta_vdj_germline_reconstructed" in asm["qc_warnings"].iloc[0]

    # 2. Strict opt-out leaves both chains unassembled
    asm_strict = assemble_full_sequences(
        clones,
        germline_vdj_fallback=False,
        data_dir=cache,
        verbose=False,
        show_progress=False,
    )
    assert (
        "single_chain_aa" not in asm_strict.columns
        or pd.isna(asm_strict["single_chain_aa"].iloc[0])
    )

    # 3. With empty cache: does not crash, logs actionable warning
    empty_cache = tmp_path / "empty_cache"
    empty_cache.mkdir()
    vr.clear_reference_cache()
    with caplog.at_level(logging.WARNING, logger="tcrsift.assemble.construct"):
        asm_uncached = assemble_full_sequences(
            clones, data_dir=empty_cache, verbose=False, show_progress=False
        )
    assert (
        "single_chain_aa" not in asm_uncached.columns
        or pd.isna(asm_uncached["single_chain_aa"].iloc[0])
    )
    assert any("tcrsift data download --group germline" in rec.message for rec in caplog.records)
    vr.clear_reference_cache()
