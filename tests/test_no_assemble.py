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

    def test_gate_logic(self):
        # The run gate: assemble only when enabled AND not --no-assemble.
        def gate(enabled, no_assemble):
            return enabled and not no_assemble
        assert gate(True, False) is True       # default: assemble
        assert gate(True, True) is False        # --no-assemble overrides
        assert gate(False, False) is False      # config disabled
        assert gate(False, True) is False
