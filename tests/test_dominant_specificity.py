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

"""Dominant-condition-fraction antigen-specific selection (#316).

A clone is antigen-specific when it is essentially confined to ONE condition
(peptide pool / antigen), clears an absolute per-condition cell floor, and is
not a viral/public bystander.
"""

from __future__ import annotations

import pandas as pd

from tcrsift.selection import select_by_dominant_specificity


def _long():
    # B1: 19/20 in P2 (spec .95, dominant 19) — specific
    # B2: 5/5/5 across P1/P2/P3 (spec .33) — bystander
    # B3: 3/3 in P3, CD8 (culture-rare)
    rows = [
        {"CDR3_beta": "B1", "condition": "P2", "Num_Combined_TCRs": 19, "lineage": "CD4+"},
        {"CDR3_beta": "B1", "condition": "P1", "Num_Combined_TCRs": 1, "lineage": "CD4+"},
        {"CDR3_beta": "B2", "condition": "P1", "Num_Combined_TCRs": 5, "lineage": "CD4+"},
        {"CDR3_beta": "B2", "condition": "P2", "Num_Combined_TCRs": 5, "lineage": "CD4+"},
        {"CDR3_beta": "B2", "condition": "P3", "Num_Combined_TCRs": 5, "lineage": "CD4+"},
        {"CDR3_beta": "B3", "condition": "P3", "Num_Combined_TCRs": 3, "lineage": "CD8+"},
    ]
    return pd.DataFrame(rows)


class TestDominantSpecificity:
    def test_confinement_and_absolute_floor(self):
        res = select_by_dominant_specificity(
            _long(), condition_col="condition", min_specificity=0.95,
            min_dominant_cells=5,
        ).set_index("CDR3_beta")
        # B1 confined to P2 with 19 cells → specific.
        assert bool(res.loc["B1", "is_dominant_specific"]) is True
        assert res.loc["B1", "dominant_condition"] == "P2"
        assert abs(res.loc["B1", "specificity"] - 0.95) < 1e-9
        # B2 spread across three pools → bystander.
        assert bool(res.loc["B2", "is_dominant_specific"]) is False
        # B3 confined but only 3 cells < default floor 5 → not specific.
        assert bool(res.loc["B3", "is_dominant_specific"]) is False

    def test_per_lineage_min(self):
        # CD8 floor 2 admits B3 (3 cells); CD4 floor 5 still excludes it.
        res = select_by_dominant_specificity(
            _long(), condition_col="condition", lineage_col="lineage",
            per_lineage_min={"CD8": 2, "CD4": 5},
        ).set_index("CDR3_beta")
        assert res.loc["B3", "min_dominant_cells"] == 2
        assert bool(res.loc["B3", "is_dominant_specific"]) is True
        assert res.loc["B1", "min_dominant_cells"] == 5
        assert bool(res.loc["B1", "is_dominant_specific"]) is True

    def test_exclude_clones(self):
        res = select_by_dominant_specificity(
            _long(), condition_col="condition", min_dominant_cells=1,
            exclude_clones={"B1"},
        ).set_index("CDR3_beta")
        # B1 would qualify but is excluded (viral/public bystander).
        assert bool(res.loc["B1", "is_dominant_specific"]) is False

    def test_antigen_label_column(self):
        res = select_by_dominant_specificity(
            _long(), condition_col="condition",
            label_map={"P2": "P2 (KIF1C)", "P3": "P3 (EXOC4)"},
        ).set_index("CDR3_beta")
        assert res.loc["B1", "dominant_antigen"] == "P2 (KIF1C)"
        # No dominant_antigen column when no label_map is given.
        res2 = select_by_dominant_specificity(_long(), condition_col="condition")
        assert "dominant_antigen" not in res2.columns

    def test_count_column_choosable(self):
        df = _long().rename(columns={"Num_Combined_TCRs": "Num_Complete_TCRs"})
        res = select_by_dominant_specificity(
            df, condition_col="condition", count_column="Num_Complete_TCRs",
            min_dominant_cells=5,
        ).set_index("CDR3_beta")
        assert bool(res.loc["B1", "is_dominant_specific"]) is True

    def test_missing_count_column_raises(self):
        df = _long().drop(columns=["Num_Combined_TCRs"])
        try:
            select_by_dominant_specificity(df, condition_col="condition")
        except ValueError as e:
            assert "count column" in str(e).lower()
        else:
            raise AssertionError("expected ValueError for missing count column")

    def test_beta_keyable_via_clone_col(self):
        df = _long().rename(columns={"CDR3_beta": "beta"})
        res = select_by_dominant_specificity(
            df, condition_col="condition", clone_col="beta", min_dominant_cells=5,
        )
        assert "beta" in res.columns
        assert set(res["beta"]) == {"B1", "B2", "B3"}
