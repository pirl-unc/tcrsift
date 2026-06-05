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

"""Finishing #122 (filter per-sample tier) + #123 (report sequence/cohort/polished)."""

from __future__ import annotations

import matplotlib
import pandas as pd

matplotlib.use("Agg")

from tcrsift.cli import create_parser  # noqa: E402
from tcrsift.selection import pivot_per_sample_tiers  # noqa: E402


def _run(argv):
    args = create_parser().parse_args(argv)
    args.func(args)


# --- #122: filter --emit-per-sample-tier ------------------------------------
class TestPerSampleTier:
    def test_pivot_per_sample_tiers(self):
        long = pd.DataFrame({
            "CDR3ab": ["A", "A", "B"],
            "sample": ["S1", "S2", "S1"],
            "cells": [15, 2, 1],
            "frequency": [0.15, 0.02, 0.01],
        })
        wide = pivot_per_sample_tiers(long)
        assert set(wide.columns) >= {"CDR3ab", "tier_in_S1", "tier_in_S2"}
        a = wide[wide["CDR3ab"] == "A"].iloc[0]
        assert a["tier_in_S1"] == "tier1"  # 15 cells / 15% -> tier1

    def test_filter_emit_per_sample_tier_cli(self, tmp_path):
        clones = pd.DataFrame({
            "CDR3ab": ["A", "B"], "cell_count": [15, 8], "max_frequency": [0.15, 0.05],
            "n_samples": [2, 1], "CDR3_alpha": ["x", "y"], "CDR3_beta": ["p", "q"],
        })
        ci = tmp_path / "clonotypes.csv"
        clones.to_csv(ci, index=False)
        long = pd.DataFrame({
            "CDR3ab": ["A", "A", "B"], "sample": ["S1", "S2", "S1"],
            "cells": [15, 12, 8], "frequency": [0.15, 0.13, 0.08],
        })
        li = tmp_path / "long.csv"
        long.to_csv(li, index=False)
        out = tmp_path / "filt"
        _run(["filter", "-i", str(ci), "-o", str(out), "--min-cells", "1",
              "--tcell-type", "both", "--emit-per-sample-tier", "--clone-sample-long", str(li)])
        written = pd.read_csv(out / "all_filtered.csv")
        assert any(c.startswith("tier_in_") for c in written.columns)


# --- #123: report sequence --------------------------------------------------
class TestReportSequence:
    def test_emits_pdf_with_overlay(self, tmp_path):
        sel = pd.DataFrame([{
            "CDR3ab": "A_B", "selection_rule": "shared", "global_rank": 1,
            "full_alpha_aa": "M" + "A" * 30, "full_beta_aa": "M" + "G" * 30,
            "single_chain_aa": "M" + "A" * 60,
        }])
        si = tmp_path / "sel.csv"
        sel.to_csv(si, index=False)
        long = pd.DataFrame({
            "CDR3ab": ["A_B", "A_B"], "sample": ["S1", "S2"], "method": ["AIMpos", "CTYneg"],
            "cells": [12, 3], "frequency": [0.12, 0.03],
        })
        li = tmp_path / "long.csv"
        long.to_csv(li, index=False)
        out = tmp_path / "seq.pdf"
        _run(["report", "sequence", "-i", str(si), "-o", str(out),
              "--per-method-source", str(li), "--sort-by", "global_rank"])
        assert out.exists() and out.stat().st_size > 0


# --- #123: report cohort ----------------------------------------------------
class TestReportCohort:
    def _donor(self, base, name, clones):
        d = base / name
        (d / "data").mkdir(parents=True)
        pd.DataFrame(clones).to_csv(d / "data" / "clonotypes.csv", index=False)
        return d

    def test_emits_overlap_tables(self, tmp_path):
        d1 = self._donor(tmp_path, "B1-2", {
            "CDR3ab": ["A_B", "C_D"], "cell_count": [10, 5],
            "CDR3_alpha": ["A", "C"], "CDR3_beta": ["B", "D"],
        })
        d2 = self._donor(tmp_path, "B1-3", {
            "CDR3ab": ["A_B", "E_F"], "cell_count": [8, 4],
            "CDR3_alpha": ["A", "E"], "CDR3_beta": ["B", "F"],
        })
        out = tmp_path / "cohort"
        _run(["report", "cohort", "--donor", f"B1-2={d1}", "--donor", f"B1-3={d2}",
              "-o", str(out), "--no-plots"])
        csvs = list(out.glob("*.csv"))
        assert csvs, "expected cohort overlap CSV exports"


# --- #123: report polished --------------------------------------------------
class TestReportPolished:
    def test_rerenders_with_style(self, tmp_path):
        run_dir = tmp_path / "run"
        (run_dir / "data").mkdir(parents=True)
        pd.DataFrame({
            "CDR3ab": ["A_B", "C_D", "E_F"], "cell_count": [10, 5, 1],
            "alpha_v_gene": ["TRAV1", "TRAV2", "TRAV1"],
            "beta_v_gene": ["TRBV1", "TRBV2", "TRBV1"],
        }).to_csv(run_dir / "data" / "clonotypes.csv", index=False)
        out = tmp_path / "polished"
        _run(["report", "polished", "-i", str(run_dir), "-o", str(out), "--style", "paper"])
        # Polished render uses plot_format=pdf -> at least one figure emitted.
        assert out.exists()
        assert any(out.iterdir())

    def test_set_polished_style_applies_font(self):
        import matplotlib as mpl

        from tcrsift.plots import set_polished_style

        set_polished_style("clean-white")
        assert "DejaVu Sans" in mpl.rcParams["font.family"]
