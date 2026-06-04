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

"""Tests for the standalone `tcrsift select` command (#122)."""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.cli import create_parser
from tcrsift.validation import TCRsiftValidationError


def _long_csv(path):
    # Two clones; A1_B1 is a strong tier1 in AIMpos, A2_B2 is weak.
    pd.DataFrame({
        "CDR3ab": ["A1_B1", "A1_B1", "A2_B2"],
        "sample": ["S1", "S2", "S1"],
        "method": ["AIMpos", "AIMpos", "CTYneg"],
        "cells": [15, 12, 1],
        "frequency": [0.15, 0.13, 0.001],
    }).to_csv(path, index=False)


def _config_yaml(path):
    path.write_text(
        "selection:\n"
        "  rules:\n"
        "    shared:\n"
        "      include_tiers: [tier1, tier2]\n"
        "  global_rank:\n"
        "    rule_order: [shared, method_pair, private]\n"
    )


def _run(argv):
    args = create_parser().parse_args(argv)
    args.func(args)


def test_select_emits_selected_clones(tmp_path):
    long_csv = tmp_path / "clone_sample_long.csv"
    cfg = tmp_path / "config.yaml"
    out = tmp_path / "selected.csv"
    _long_csv(long_csv)
    _config_yaml(cfg)

    _run(["select", "-i", str(long_csv), "-c", str(cfg), "-o", str(out)])

    assert out.exists()
    sel = pd.read_csv(out)
    assert "selection_rule" in sel.columns and "global_rank" in sel.columns
    # The strong tier1 clone is selected by the 'shared' rule.
    assert "A1_B1" in set(sel["CDR3ab"])
    assert (sel["selection_rule"] == "shared").any()


def test_select_joins_clonotype_metadata(tmp_path):
    long_csv = tmp_path / "clone_sample_long.csv"
    cfg = tmp_path / "config.yaml"
    clo = tmp_path / "clonotypes.csv"
    out = tmp_path / "selected.csv"
    _long_csv(long_csv)
    _config_yaml(cfg)
    pd.DataFrame({"CDR3ab": ["A1_B1"], "alpha_v_gene": ["TRAV12-2"]}).to_csv(clo, index=False)

    _run(["select", "-i", str(long_csv), "-c", str(cfg), "-o", str(out),
          "--clonotypes", str(clo)])

    sel = pd.read_csv(out)
    row = sel[sel["CDR3ab"] == "A1_B1"].iloc[0]
    assert row["alpha_v_gene"] == "TRAV12-2"


def test_select_emits_per_sample_tiers(tmp_path):
    long_csv = tmp_path / "clone_sample_long.csv"
    cfg = tmp_path / "config.yaml"
    out = tmp_path / "selected.csv"
    _long_csv(long_csv)
    _config_yaml(cfg)

    _run(["select", "-i", str(long_csv), "-c", str(cfg), "-o", str(out),
          "--emit-per-sample-tiers"])

    tier_file = tmp_path / "clone_sample_per_sample_tiers.csv"
    assert tier_file.exists()
    tiers = pd.read_csv(tier_file)
    assert "per_sample_tier" in tiers.columns


def test_select_without_rules_errors(tmp_path):
    long_csv = tmp_path / "clone_sample_long.csv"
    cfg = tmp_path / "empty.yaml"
    out = tmp_path / "selected.csv"
    _long_csv(long_csv)
    cfg.write_text("selection: {}\n")

    with pytest.raises(TCRsiftValidationError):
        _run(["select", "-i", str(long_csv), "-c", str(cfg), "-o", str(out)])
