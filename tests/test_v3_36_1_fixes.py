"""Regression tests for v3.36.1 holistic audit quick fixes."""

from __future__ import annotations

import argparse
import inspect
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tcrsift
from tcrsift import annotate, cli, mnemonic
from tcrsift.assemble import codons, germline
from tcrsift.config import AssembleConfig, TCRsiftConfig
from tcrsift.provenance import collect_provenance
from tcrsift.validation import TCRsiftValidationError, validate_assemble_args, validate_run_args


def test_top_level_match_clonotypes_is_annotate_not_til_select():
    """tcrsift.match_clonotypes resolves to tcrsift.annotate.match_clonotypes."""
    assert tcrsift.match_clonotypes is annotate.match_clonotypes
    sig = inspect.signature(tcrsift.match_clonotypes)
    assert "match_strictness" in sig.parameters
    assert "host_species_filter" in sig.parameters


def test_mnemonic_module_has_docstring():
    """tcrsift.mnemonic.__doc__ is populated (docstring precedes __future__ import)."""
    assert isinstance(mnemonic.__doc__, str)
    assert "mnemonic" in mnemonic.__doc__.lower()


def test_python_m_tcrsift_version():
    """python -m tcrsift --version works cleanly."""
    proc = subprocess.run(
        [sys.executable, "-m", "tcrsift", "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert f"tcrsift {tcrsift.__version__}" in proc.stdout


def test_assemble_config_normalizes_string_none_leaders():
    """YAML assemble.alpha_leader: none normalizes to None."""
    cfg = TCRsiftConfig._from_dict({"assemble": {"alpha_leader": "none", "beta_leader": "NONE"}})
    assert cfg.assemble.alpha_leader is None
    assert cfg.assemble.beta_leader is None
    direct = AssembleConfig(alpha_leader="none", beta_leader="CD8A")
    assert direct.alpha_leader is None
    assert direct.beta_leader == "CD8A"


def test_cmd_phenotype_applies_cd3_gate(tmp_path):
    """Standalone `tcrsift phenotype` filters cells failing --min-cd3-reads."""
    X = np.array(
        [
            [0.0, 10.0, 0.0, 5.0, 0.0, 0.0],  # CD3=5 (pass), CD8
            [0.0, 10.0, 0.0, 1.0, 0.0, 0.0],  # CD3=1 (fail min_cd3_reads=3)
            [10.0, 0.0, 0.0, 4.0, 2.0, 0.0],  # CD3=6 (pass), CD4
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {
                "sample": ["S1", "S1", "S1"],
                "CD4": [0.0, 0.0, 10.0],
                "CD8": [10.0, 10.0, 0.0],
                "CD3": [5.0, 1.0, 6.0],
            },
            index=["bc1", "bc2", "bc3"],
        ),
        var=pd.DataFrame(index=["CD4", "CD8A", "CD8B", "CD3D", "CD3E", "CD3G"]),
    )
    in_path = tmp_path / "in.h5ad"
    out_path = tmp_path / "out.h5ad"
    adata.write_h5ad(in_path)

    cli.main(["phenotype", "-i", str(in_path), "-o", str(out_path), "--min-cd3-reads", "3"])
    out = ad.read_h5ad(out_path)
    assert list(out.obs_names) == ["bc1", "bc3"]


def test_report_selected_registers_germline_fallback_flags():
    """`tcrsift report selected` accepts --germline-vdj-fallback, --auto-download-germline, --cache-dir."""
    parser = cli.create_parser()
    args = parser.parse_args(
        [
            "report",
            "selected",
            "--selection",
            "sel.csv",
            "--clonotypes",
            "cl.csv",
            "-o",
            "out",
            "--germline-vdj-fallback",
            "--auto-download-germline",
            "--cache-dir",
            "/tmp/cache",
        ]
    )
    assert args.germline_vdj_fallback is True
    assert args.auto_download_germline is True
    assert args.cache_dir == "/tmp/cache"


def test_validate_args_respects_no_assemble_no_leaders_and_yaml_from_contig(tmp_path):
    """validate_assemble_args and validate_run_args handle disabled leaders/assembly and YAML from_contig."""
    ss = tmp_path / "samples.yaml"
    ss.write_text("samples:\n  - sample: S1\n")

    # --no-leaders disables from_contig requirement in assemble
    validate_assemble_args(
        argparse.Namespace(
            alpha_leader="from_contig",
            beta_leader="CD8A",
            leaders_from_contigs=False,
            no_leaders=True,
            contigs_dir=None,
            cellranger_dir=None,
        )
    )

    # --no-assemble disables from_contig requirement in run
    validate_run_args(
        argparse.Namespace(
            sample_sheet=str(ss),
            alpha_leader=None,
            beta_leader=None,
            leaders_from_contigs=True,
            no_assemble=True,
            contigs_dir=None,
            cellranger_dir=None,
        )
    )

    # YAML alpha_leader=from_contig without contigs_dir fails validation
    cfg = TCRsiftConfig._from_dict({"assemble": {"alpha_leader": "from_contig"}})
    with pytest.raises(TCRsiftValidationError, match="--contigs-dir is required"):
        validate_run_args(
            argparse.Namespace(
                sample_sheet=str(ss),
                alpha_leader=None,
                beta_leader=None,
                leaders_from_contigs=False,
                no_assemble=False,
                contigs_dir=None,
                cellranger_dir=None,
            ),
            config=cfg,
        )


def test_reconstruct_germline_vdj_chain_inherits_mouse_codon_species():
    """reconstruct_germline_vdj_chain(species='mouse') uses mouse codon table by default."""
    # Amino acid sequence whose preferred human vs mouse codon differs (e.g. G/S/T/P/R/I)
    row = {"VDJ_beta_aa": "CASSQDRTGELFF", "VDJ_beta_nt": None}
    aa, nt_mouse, src, _ = germline.reconstruct_germline_vdj_chain(row, "beta", species="mouse")
    _, nt_human, _, _ = germline.reconstruct_germline_vdj_chain(row, "beta", species="human")
    assert aa == "CASSQDRTGELFF"
    assert src == germline.VDJ_SOURCE_CONTIG_FULL
    assert nt_mouse == codons.back_translate("CASSQDRTGELFF", species="mouse").lower()
    assert nt_human == codons.back_translate("CASSQDRTGELFF", species="human").lower()


def test_provenance_records_oncoref():
    """collect_provenance includes oncoref in dependency versions."""
    prov = collect_provenance()
    assert "oncoref" in prov["dependencies"]
