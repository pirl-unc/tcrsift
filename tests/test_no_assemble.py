# Licensed under the Apache License, Version 2.0 (the "License").
"""#246: `run --no-assemble` / assemble.enabled=false — select-then-assemble."""
from __future__ import annotations

from tcrsift.cli import create_parser
from tcrsift.config import AssembleConfig, TCRsiftConfig


class TestNoAssembleFlag:
    def test_flag_parses(self):
        a = create_parser().parse_args(["run", "-s", "s.csv", "-o", "out", "--no-assemble"])
        assert a.no_assemble is True

    def test_default_is_off(self):
        a = create_parser().parse_args(["run", "-s", "s.csv", "-o", "out"])
        assert getattr(a, "no_assemble", False) is False


class TestAssembleEnabledConfig:
    def test_default_enabled(self):
        assert AssembleConfig().enabled is True

    def test_from_yaml_can_disable(self, tmp_path):
        p = tmp_path / "c.yaml"
        p.write_text("assemble:\n  enabled: false\n")
        cfg = TCRsiftConfig.from_yaml(str(p))
        assert cfg.assemble.enabled is False

    # (removed test_gate_logic: it asserted on a closure re-defining the
    # `enabled and not no_assemble` gate *inside the test*, so it exercised no
    # product code. The real gate lives in cli.cmd_run; its two inputs are
    # already covered by test_flag_parses / test_default_enabled / the yaml test.)
