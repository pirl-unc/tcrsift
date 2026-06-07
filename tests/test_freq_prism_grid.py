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

"""freq × PRISM trade-off grid helper + plot (#207)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from tcrsift.selection import freq_prism_grid, select_freq_prism_per_condition


def _feat(n_per_cond=12, n_cond=2, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for ci in range(n_cond):
        for i in range(n_per_cond):
            rows.append({
                "CDR3ab": f"clone_{ci}_{i}",
                "method": f"cond{ci}",
                "frequency": float(rng.uniform(0.01, 0.5)),
                "ppost_alpha": float(rng.uniform(0, 1)),
                "ppost_beta": float(rng.uniform(0, 1)),
                "antigen_response_score": float(rng.uniform(0, 1)),
                "naive_score": float(rng.uniform(0, 1)),
            })
    return pd.DataFrame(rows)


class TestFreqPrismGrid:
    def test_shape_and_columns(self):
        feat = _feat()
        grid = freq_prism_grid(
            feat, condition_col="method",
            freq_values=(0, 5, 10), prism_values=(0, 5, 10),
        )
        assert list(grid.columns) == ["top_freq", "top_prism", "n_clones"]
        assert len(grid) == 9  # 3 × 3
        assert (grid["n_clones"] >= 0).all()

    def test_monotone_nondecreasing_in_each_axis(self):
        feat = _feat(seed=3)
        grid = freq_prism_grid(
            feat, condition_col="method",
            freq_values=(0, 5, 10, 15), prism_values=(0, 5, 10, 15),
        )
        piv = grid.pivot(index="top_prism", columns="top_freq", values="n_clones")
        # More picks (either axis) can only add clones to the union.
        for _, row in piv.iterrows():
            assert list(row.values) == sorted(row.values)
        for _, col in piv.items():
            assert list(col.values) == sorted(col.values)

    def test_origin_is_zero(self):
        feat = _feat()
        grid = freq_prism_grid(
            feat, condition_col="method", freq_values=(0,), prism_values=(0,),
        )
        assert grid.loc[0, "n_clones"] == 0

    def test_cell_matches_direct_selection(self):
        feat = _feat(seed=7)
        grid = freq_prism_grid(
            feat, condition_col="method", freq_values=(10,), prism_values=(5,),
        )
        direct = select_freq_prism_per_condition(
            feat, condition_col="method", top_freq=10, top_prism=5,
        )
        assert grid.loc[0, "n_clones"] == direct["CDR3ab"].nunique()

    def test_empty_input(self):
        grid = freq_prism_grid(
            pd.DataFrame(columns=["CDR3ab", "method", "frequency"]),
            condition_col="method", freq_values=(0, 5), prism_values=(0, 5),
        )
        assert (grid["n_clones"] == 0).all()


class TestPlotFreqPrismGrid:
    def test_writes_png(self, tmp_path):
        from tcrsift.plots import plot_freq_prism_grid

        feat = _feat()
        grid = freq_prism_grid(feat, condition_col="method")
        out = tmp_path / "grid.png"
        plot_freq_prism_grid(grid, out, chosen=(10, 5))
        assert out.exists() and out.stat().st_size > 0

    def test_empty_grid_noop(self, tmp_path):
        from tcrsift.plots import plot_freq_prism_grid

        out = tmp_path / "grid.png"
        plot_freq_prism_grid(
            pd.DataFrame(columns=["top_freq", "top_prism", "n_clones"]), out
        )
        assert not out.exists()
