"""Regression tests for v3.37.0 dead-code removal and structural cleanup."""

from __future__ import annotations

import importlib
import subprocess
import sys

import pytest

from tcrsift import assemble, cli, qc, qc_battery
from tcrsift.config import TCRsiftConfig
from tcrsift.validation import TCRsiftValidationError


def test_dead_modules_removed_and_assemble_helpers_moved():
    """Dead tcrsift.data is gone; leaders/sigpep/vregion live in tcrsift.assemble."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("tcrsift.data")
    assert importlib.import_module("tcrsift.assemble.leaders") is not None
    assert importlib.import_module("tcrsift.assemble.sigpep") is not None
    assert importlib.import_module("tcrsift.assemble.vregion") is not None
    assert hasattr(assemble, "germline_vgene_leaders")
    assert hasattr(assemble, "signal_peptide_features")


def test_model_module_has_no_eager_matplotlib_import():
    """Importing tcrsift.model in a fresh interpreter does not import matplotlib."""
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, tcrsift.model; assert 'matplotlib' not in sys.modules",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert proc.returncode == 0


def test_qc_battery_result_no_longer_collides_with_qc_result():
    """qc_battery uses BatteryCheckResult instead of shadowing tcrsift.qc.QCResult."""
    assert hasattr(qc_battery, "BatteryCheckResult")
    assert not hasattr(qc_battery, "QCResult")
    assert hasattr(qc, "QCResult")


def test_config_handles_empty_yaml_section_and_invalid_nested_key(caplog):
    """Empty YAML section does not warn; typo in nested section raises TCRsiftValidationError."""
    with caplog.at_level("WARNING", logger="tcrsift.config"):
        cfg = TCRsiftConfig._from_dict({"load": None, "assemble": {}})
    assert cfg.load.min_genes == 250
    assert "Unknown configuration key" not in caplog.text

    with pytest.raises(TCRsiftValidationError, match="Unknown configuration key.*min_gene"):
        TCRsiftConfig._from_dict({"load": {"min_gene": 200}})


def test_report_and_pgen_without_subcommand_exit_nonzero():
    """Invoking `tcrsift report` or `tcrsift pgen` with no subcommand exits non-zero."""
    for cmd in (["report"], ["pgen"]):
        with pytest.raises(SystemExit) as exc_info:
            cli.main(cmd)
        assert exc_info.value.code == 1


def test_prioritize_exposes_allow_canonical_fallback_toggle():
    """`tcrsift prioritize` supports --allow-canonical-fallback / --no-allow-canonical-fallback."""
    parser = cli.create_parser()
    default_args = parser.parse_args(["prioritize", "samples.yaml", "-o", "out"])
    strict_args = parser.parse_args(
        ["prioritize", "samples.yaml", "-o", "out", "--no-allow-canonical-fallback"]
    )
    assert default_args.allow_canonical_fallback is True
    assert strict_args.allow_canonical_fallback is False
