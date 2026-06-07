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

"""N-way set overlap (UpSet) over selected-clone sets (#208)."""

from __future__ import annotations

import pandas as pd

from tcrsift.clonotype import build_selection_sets, set_overlap_table


class TestBuildSelectionSets:
    def test_groups_clones_by_set_col(self):
        long_df = pd.DataFrame({
            "CDR3ab": ["A", "B", "A", "C"],
            "method": ["m1", "m1", "m2", "m2"],
        })
        sets = build_selection_sets(long_df, set_col="method")
        assert sets == {"m1": {"A", "B"}, "m2": {"A", "C"}}

    def test_clone_universe_restriction(self):
        long_df = pd.DataFrame({
            "CDR3ab": ["A", "B", "C"],
            "method": ["m1", "m1", "m2"],
        })
        sets = build_selection_sets(long_df, set_col="method", clones={"A", "C"})
        assert sets == {"m1": {"A"}, "m2": {"C"}}

    def test_missing_col_returns_empty(self):
        long_df = pd.DataFrame({"CDR3ab": ["A"]})
        assert build_selection_sets(long_df, set_col="method") == {}

    def test_empty_sets_dropped(self):
        long_df = pd.DataFrame({"CDR3ab": ["A"], "method": ["m1"]})
        sets = build_selection_sets(long_df, set_col="method", clones={"Z"})
        assert sets == {}


class TestSetOverlapTable:
    def test_intersection_patterns(self):
        sets = {"m1": {"A", "B"}, "m2": {"A", "C"}}
        table = set_overlap_table(sets)
        assert list(table.columns) == ["sets", "degree", "n_clones"]
        rows = {r["sets"]: r["n_clones"] for _, r in table.iterrows()}
        assert rows["m1+m2"] == 1   # A
        assert rows["m1"] == 1      # B
        assert rows["m2"] == 1      # C
        # Every clone counted exactly once.
        assert table["n_clones"].sum() == 3

    def test_degree_column(self):
        sets = {"a": {"X"}, "b": {"X"}, "c": {"X"}}
        table = set_overlap_table(sets)
        assert len(table) == 1
        assert table.loc[0, "degree"] == 3
        assert table.loc[0, "sets"] == "a+b+c"
        assert table.loc[0, "n_clones"] == 1

    def test_sorted_by_size(self):
        sets = {"m1": {"A", "B", "C"}, "m2": {"A"}}
        table = set_overlap_table(sets)
        # m1-only (B, C => 2) should rank above the shared (A => 1).
        assert table.iloc[0]["n_clones"] >= table.iloc[-1]["n_clones"]

    def test_empty(self):
        table = set_overlap_table({})
        assert table.empty
        assert list(table.columns) == ["sets", "degree", "n_clones"]


class TestPlotSetOverlap:
    def test_writes_png(self, tmp_path):
        from tcrsift.plots import plot_set_overlap

        sets = {"m1": {"A", "B"}, "m2": {"A", "C"}, "m3": {"B", "C", "D"}}
        out = tmp_path / "upset.png"
        plot_set_overlap(sets, out)
        assert out.exists() and out.stat().st_size > 0

    def test_noop_single_set(self, tmp_path):
        from tcrsift.plots import plot_set_overlap

        out = tmp_path / "upset.png"
        plot_set_overlap({"m1": {"A", "B"}}, out)
        assert not out.exists()

    def test_fallback_bar_without_upsetplot(self, tmp_path, monkeypatch):
        import builtins

        from tcrsift.plots import plot_set_overlap

        real_import = builtins.__import__

        def no_upsetplot(name, *a, **k):
            if name == "upsetplot":
                raise ImportError("forced")
            return real_import(name, *a, **k)

        monkeypatch.setattr(builtins, "__import__", no_upsetplot)
        out = tmp_path / "bar.png"
        plot_set_overlap({"m1": {"A", "B"}, "m2": {"A", "C"}}, out)
        assert out.exists() and out.stat().st_size > 0
