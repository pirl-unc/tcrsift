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

"""#195 generate-config overwrite guard + #196 report-selected contig overrides."""

from __future__ import annotations

import pytest

from tcrsift.cli import create_parser


def _run(argv):
    args = create_parser().parse_args(argv)
    args.func(args)


# --- #195: generate-config no silent overwrite ------------------------------
class TestGenerateConfigOverwrite:
    def test_stdout_by_default(self, capsys, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        _run(["generate-config"])
        out = capsys.readouterr().out
        assert "load:" in out and "attribution:" in out
        # Nothing written to disk.
        assert not (tmp_path / "tcrsift_config.yaml").exists()

    def test_refuses_to_overwrite_existing(self, tmp_path):
        target = tmp_path / "cfg.yaml"
        target.write_text("min_genes: 999  # my customized config\n")
        with pytest.raises(SystemExit):
            _run(["generate-config", "-o", str(target)])
        # The original is untouched.
        assert "my customized config" in target.read_text()

    def test_force_overwrites(self, tmp_path):
        target = tmp_path / "cfg.yaml"
        target.write_text("min_genes: 999\n")
        _run(["generate-config", "-o", str(target), "--force"])
        assert "load:" in target.read_text()
        assert "min_genes: 999\n" != target.read_text()

    def test_writes_when_absent(self, tmp_path):
        target = tmp_path / "new.yaml"
        _run(["generate-config", "-o", str(target)])
        assert target.exists() and "load:" in target.read_text()


# --- #196: report selected contig overrides ---------------------------------
class TestReportSelectedContigArgs:
    def test_contig_flags_parse_and_reach_args(self):
        a = create_parser().parse_args([
            "report", "selected", "--selection", "s.csv", "--clonotypes", "c.csv",
            "-o", "out", "--cellranger-dir", "B1-2/per_sample_outs",
            "--contigs-dir", "ctg", "--sample-name-from", "grandparent",
        ])
        assert a.cellranger_dir == "B1-2/per_sample_outs"
        assert a.contigs_dir == "ctg"
        assert a.sample_name_from == "grandparent"

    def test_overrides_applied(self, tmp_path, monkeypatch):
        # The CLI override should land in assemble_kwargs even without a config.
        import pandas as pd

        import tcrsift.cli as cli

        captured = {}

        def fake_build(selected, clonotypes, output_dir, *, obs=None,
                       assemble_kwargs=None, provenance_cols=None, **_kw):
            captured.update(assemble_kwargs or {})
            return pd.DataFrame()

        monkeypatch.setattr("tcrsift.report.build_selected_report", fake_build)
        sel = tmp_path / "s.csv"
        pd.DataFrame({"CDR3ab": ["A_B"]}).to_csv(sel, index=False)
        clo = tmp_path / "c.csv"
        pd.DataFrame({"CDR3ab": ["A_B"]}).to_csv(clo, index=False)
        args = create_parser().parse_args([
            "report", "selected", "--selection", str(sel), "--clonotypes", str(clo),
            "-o", str(tmp_path / "out"), "--cellranger-dir", "B1-2/pso",
            "--sample-name-from", "grandparent",
        ])
        cli.cmd_report_selected(args)
        assert captured["cellranger_dir"] == "B1-2/pso"
        assert captured["sample_name_from"] == "grandparent"
