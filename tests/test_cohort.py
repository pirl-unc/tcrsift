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

"""Tests for multi-donor cohort overlap analysis (#125 increment 7)."""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift import cohort


def _write_donor(root, name, clones, cells=None, *, selected=False, beta=None):
    d = root / name / "data"
    d.mkdir(parents=True, exist_ok=True)
    cols = {"CDR3ab": clones}
    if cells is not None:
        cols["cell_count"] = cells
    if beta is not None:
        cols["CDR3_beta"] = beta
    fname = "selected_clones.csv" if selected else "clonotypes.csv"
    pd.DataFrame(cols).to_csv(d / fname, index=False)
    return root / name


class TestReadDonorClones:
    def test_reads_clonotypes_with_cells(self, tmp_path):
        ddir = _write_donor(tmp_path, "B1-2", ["A_X", "B_Y"], [10, 5])
        out = cohort.read_donor_clones(ddir)
        assert list(out["CDR3ab"]) == ["A_X", "B_Y"]
        assert list(out["cells"]) == [10.0, 5.0]
        assert list(out["CDR3_beta"]) == ["X", "Y"]  # parsed from CDR3ab

    def test_cell_count_defaults_to_one(self, tmp_path):
        ddir = _write_donor(tmp_path, "D", ["A_X"])  # no cell column
        assert list(cohort.read_donor_clones(ddir)["cells"]) == [1.0]

    def test_selected_only_reads_selected_file(self, tmp_path):
        _write_donor(tmp_path, "D", ["A_X", "B_Y"], [10, 5])  # clonotypes
        _write_donor(tmp_path, "D", ["A_X"], selected=True)   # selected subset
        out = cohort.read_donor_clones(tmp_path / "D", selected_only=True)
        assert list(out["CDR3ab"]) == ["A_X"]

    def test_explicit_beta_column_preferred(self, tmp_path):
        ddir = _write_donor(tmp_path, "D", ["A_X"], [3], beta=["BETA1"])
        assert list(cohort.read_donor_clones(ddir)["CDR3_beta"]) == ["BETA1"]

    def test_missing_file_raises(self, tmp_path):
        (tmp_path / "empty").mkdir()
        with pytest.raises(FileNotFoundError):
            cohort.read_donor_clones(tmp_path / "empty")

    def test_missing_cdr3ab_raises(self, tmp_path):
        d = tmp_path / "bad" / "data"
        d.mkdir(parents=True)
        pd.DataFrame({"foo": [1]}).to_csv(d / "clonotypes.csv", index=False)
        with pytest.raises(ValueError, match="CDR3ab"):
            cohort.read_donor_clones(tmp_path / "bad")


class TestCohortOverlap:
    def _long(self):
        return pd.DataFrame({
            "donor": ["B1-2", "B1-2", "B1-2", "B1-3", "B1-3"],
            "CDR3ab": ["A_X", "B_Y", "C_Z", "A_X", "D_W"],
            "cells": [10, 5, 2, 8, 3],
            "CDR3_beta": ["X", "Y", "Z", "X", "W"],
        })

    def test_paired_jaccard(self):
        m = cohort.cohort_overlap(self._long())
        # union {A_X,B_Y,C_Z,D_W}=4, intersection {A_X}=1 -> 0.25
        assert m["jaccard"].loc["B1-2", "B1-3"] == 0.25
        assert m["jaccard"].loc["B1-2", "B1-2"] == 1.0

    def test_beta_only_matrices_present_by_default(self):
        m = cohort.cohort_overlap(self._long())
        assert "jaccard_beta_only" in m
        # beta sets: {X,Y,Z} vs {X,W} -> inter {X}=1, union 4 -> 0.25
        assert m["jaccard_beta_only"].loc["B1-2", "B1-3"] == 0.25

    def test_beta_only_can_be_disabled(self):
        m = cohort.cohort_overlap(self._long(), include_beta_only=False)
        assert "jaccard_beta_only" not in m


class TestRunCohortAnalysis:
    def test_writes_matrices_csv(self, tmp_path):
        _write_donor(tmp_path, "B1-2", ["A_X", "B_Y", "C_Z"], [10, 5, 2])
        _write_donor(tmp_path, "B1-3", ["A_X", "D_W"], [8, 3])
        out = tmp_path / "cohort"
        m = cohort.run_cohort_analysis(
            {"B1-2": tmp_path / "B1-2", "B1-3": tmp_path / "B1-3"}, out,
        )
        assert (out / "cohort_jaccard.csv").exists()
        assert (out / "cohort_cell_weighted_jaccard.csv").exists()
        assert m["jaccard"].loc["B1-2", "B1-3"] == 0.25

    def test_no_tables_skips_csv(self, tmp_path):
        _write_donor(tmp_path, "B1-2", ["A_X"], [10])
        _write_donor(tmp_path, "B1-3", ["A_X"], [8])
        out = tmp_path / "cohort"
        cohort.run_cohort_analysis(
            {"B1-2": tmp_path / "B1-2", "B1-3": tmp_path / "B1-3"}, out,
            emit_tables=False,
        )
        assert not (out / "cohort_jaccard.csv").exists()


class TestCohortCli:
    def test_donor_specs_parse(self):
        from tcrsift.cli import create_parser

        args = create_parser().parse_args([
            "cohort", "--donor", "B1-2=out/B1-2", "--donor", "B1-3=out/B1-3",
            "-o", "cohort_out", "--selected-only",
        ])
        assert args.donor == ["B1-2=out/B1-2", "B1-3=out/B1-3"]
        assert args.selected_only is True

    def test_fewer_than_two_donors_exits(self, tmp_path):
        import argparse

        from tcrsift.cli import cmd_cohort

        args = argparse.Namespace(
            donor=["only=dir"], output_dir=str(tmp_path), selected_only=False,
            no_beta_only=False, no_tables=False, emit_plots=False,
        )
        with pytest.raises(SystemExit):
            cmd_cohort(args)

    def test_cmd_cohort_end_to_end(self, tmp_path, capsys):
        import argparse

        from tcrsift.cli import cmd_cohort

        _write_donor(tmp_path, "B1-2", ["A_X", "B_Y"], [10, 5])
        _write_donor(tmp_path, "B1-3", ["A_X"], [8])
        out = tmp_path / "cohort"
        args = argparse.Namespace(
            donor=[f"B1-2={tmp_path / 'B1-2'}", f"B1-3={tmp_path / 'B1-3'}"],
            output_dir=str(out), selected_only=False, no_beta_only=False,
            no_tables=False, emit_plots=False,
        )
        cmd_cohort(args)
        assert (out / "cohort_jaccard.csv").exists()
        assert "Cohort overlap across 2 donors" in capsys.readouterr().out
