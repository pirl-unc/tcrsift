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

"""Tests for the externally-chosen clone helpers (#302).

- ``tcrsift.gex.annotate_chosen`` — αβ-pair / CDR3β key matching of a
  scored clone table against an externally-selected clone list.
- ``tcrsift.plots.plot_signature_map`` — two-score signature scatter with
  a highlighted subset and optional categorical coloring.

The plot tests assert axis/file emission and layer behavior rather than
visual content.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from tcrsift.gex import annotate_chosen  # noqa: E402
from tcrsift.plots import plot_signature_map  # noqa: E402
from tcrsift.validation import TCRsiftValidationError  # noqa: E402


def _table(index=None):
    return pd.DataFrame(
        {
            "CDR3_alpha": ["CAA", "CAB", "cac", "CAD ", None, "CAA"],
            "CDR3_beta": ["CBA", "CBB", "CBC", "CBD", "CBE", "CBZ"],
        },
        index=index,
    )


class TestAnnotateChosenPair:
    def test_exact_pair_with_normalization(self):
        tbl = _table()
        # Mixed case / trailing space on both the table and the chosen list.
        chosen = [("caa", "cba"), ("CAC", "cbc"), ("cad ", "CBD")]
        mask = annotate_chosen(tbl, chosen, match="pair")
        assert mask.tolist() == [True, False, True, True, False, False]
        assert mask.name == "chosen"

    def test_pair_requires_both_chains(self):
        tbl = _table()
        # Right beta, wrong alpha -> no pair match (would match under beta).
        mask = annotate_chosen(tbl, [("WRONG", "CBA")], match="pair")
        assert not mask.any()

    def test_dataframe_input(self):
        tbl = _table()
        chosen = pd.DataFrame(
            {"CDR3_alpha": ["CAA", "CAB"], "CDR3_beta": ["CBA", "CBB"]}
        )
        mask = annotate_chosen(tbl, chosen, match="pair")
        assert mask.tolist() == [True, True, False, False, False, False]

    def test_empty_chains_never_match(self):
        tbl = pd.DataFrame({"CDR3_alpha": ["", "CAA"], "CDR3_beta": ["", "CBA"]})
        # An empty (α, β) in the chosen list must not match the empty row.
        mask = annotate_chosen(tbl, [("", ""), ("CAA", "CBA")], match="pair")
        assert mask.tolist() == [False, True]

    def test_index_alignment_preserved(self):
        idx = [100, 200, 300, 400, 500, 600]
        tbl = _table(index=idx)
        mask = annotate_chosen(tbl, [("CAA", "CBA")], match="pair")
        assert list(mask.index) == idx
        assert mask.loc[100]
        assert not mask.loc[600]

    def test_custom_column_names(self):
        tbl = pd.DataFrame({"a": ["CAA"], "b": ["CBA"]})
        chosen = pd.DataFrame({"a": ["CAA"], "b": ["CBA"]})
        mask = annotate_chosen(tbl, chosen, alpha_col="a", beta_col="b", match="pair")
        assert mask.tolist() == [True]

    def test_custom_name(self):
        mask = annotate_chosen(_table(), [("CAA", "CBA")], match="pair", name="is_pick")
        assert mask.name == "is_pick"


class TestAnnotateChosenBeta:
    def test_bare_beta_strings(self):
        tbl = _table()
        mask = annotate_chosen(tbl, ["cbb", "CBE"], match="beta")
        assert mask.tolist() == [False, True, False, False, True, False]

    def test_beta_ignores_alpha(self):
        tbl = _table()
        # CBA appears with alpha CAA (row 0); beta-only matches regardless.
        mask = annotate_chosen(tbl, ["CBA"], match="beta")
        assert mask.tolist() == [True, False, False, False, False, False]

    def test_beta_from_pair_tuples(self):
        tbl = _table()
        # Beta match also accepts (α, β) tuples, using only β.
        mask = annotate_chosen(tbl, [("ignored", "CBB")], match="beta")
        assert mask.tolist() == [False, True, False, False, False, False]

    def test_beta_dataframe_does_not_need_alpha(self):
        tbl = _table()
        chosen = pd.DataFrame({"CDR3_beta": ["CBC"]})  # no alpha column
        mask = annotate_chosen(tbl, chosen, match="beta")
        assert mask.tolist() == [False, False, True, False, False, False]

    def test_empty_beta_never_matches(self):
        tbl = pd.DataFrame({"CDR3_beta": ["", "CBX"]})
        mask = annotate_chosen(tbl, ["", "CBX"], match="beta")
        assert mask.tolist() == [False, True]


class TestAnnotateChosenErrors:
    def test_invalid_match(self):
        with pytest.raises(TCRsiftValidationError):
            annotate_chosen(_table(), [("CAA", "CBA")], match="fuzzy")

    def test_missing_beta_column(self):
        tbl = pd.DataFrame({"CDR3_alpha": ["CAA"]})
        with pytest.raises(TCRsiftValidationError):
            annotate_chosen(tbl, ["CBA"], match="beta")

    def test_pair_missing_alpha_column(self):
        tbl = pd.DataFrame({"CDR3_beta": ["CBA"]})
        with pytest.raises(TCRsiftValidationError):
            annotate_chosen(tbl, [("CAA", "CBA")], match="pair")

    def test_pair_input_not_tuples(self):
        with pytest.raises(TCRsiftValidationError):
            annotate_chosen(_table(), ["CBA"], match="pair")

    def test_dataframe_missing_alpha_for_pair(self):
        tbl = _table()
        chosen = pd.DataFrame({"CDR3_beta": ["CBA"]})
        with pytest.raises(TCRsiftValidationError):
            annotate_chosen(tbl, chosen, match="pair")


def _score_table(n=40, seed=0):
    rng = np.random.RandomState(seed)
    return pd.DataFrame(
        {
            "CDR3_alpha": [f"CAA{i}" for i in range(n)],
            "CDR3_beta": [f"CBB{i}" for i in range(n)],
            "antigen_response_score": rng.rand(n),
            "cytolytic_score": rng.rand(n),
            "specificity": rng.choice(["MART1", "viral", "NeoAg"], size=n),
        }
    )


class TestPlotSignatureMap:
    def teardown_method(self):
        plt.close("all")

    def test_returns_axes_chosen_vs_all(self):
        tbl = _score_table()
        hl = annotate_chosen(
            tbl, [(f"CAA{i}", f"CBB{i}") for i in (1, 4, 9)], match="pair"
        )
        ax = plot_signature_map(
            tbl, "antigen_response_score", "cytolytic_score", highlight_mask=hl
        )
        assert isinstance(ax, plt.Axes)
        # x/y labels default to the column name with underscores → spaces.
        assert ax.get_xlabel() == "antigen response score"
        assert ax.get_ylabel() == "cytolytic score"

    def test_background_mask_layer(self):
        tbl = _score_table()
        hl = annotate_chosen(tbl, [("CAA1", "CBB1")], match="pair")
        ax = plot_signature_map(
            tbl,
            "antigen_response_score",
            "cytolytic_score",
            highlight_mask=hl,
            background_mask=~hl,
        )
        assert isinstance(ax, plt.Axes)

    def test_label_by_renders_legend(self):
        tbl = _score_table()
        hl = annotate_chosen(tbl, [("CAA2", "CBB2")], match="pair")
        ax = plot_signature_map(
            tbl,
            "antigen_response_score",
            "cytolytic_score",
            highlight_mask=hl,
            label_by="specificity",
        )
        legend = ax.get_legend()
        assert legend is not None
        labels = {t.get_text() for t in legend.get_texts()}
        # Category entries plus the emphasis note.
        assert {"MART1", "viral", "NeoAg"}.issubset(labels)
        assert "highlighted" in labels

    def test_compose_into_provided_axis_does_not_save(self, tmp_path):
        tbl = _score_table()
        fig, ax0 = plt.subplots()
        out = plot_signature_map(
            tbl,
            "antigen_response_score",
            "cytolytic_score",
            label_by="specificity",
            ax=ax0,
            output_path=tmp_path / "should_not_exist.png",
        )
        assert out is ax0
        assert not (tmp_path / "should_not_exist.png").exists()

    def test_output_path_saves_and_returns_path(self, tmp_path):
        tbl = _score_table()
        hl = annotate_chosen(tbl, [("CAA1", "CBB1")], match="pair")
        out = plot_signature_map(
            tbl,
            "antigen_response_score",
            "cytolytic_score",
            highlight_mask=hl,
            output_path=tmp_path / "sigmap.png",
        )
        from pathlib import Path

        assert Path(out).exists()

    def test_highlight_none_emphasizes_all(self):
        tbl = _score_table()
        ax = plot_signature_map(tbl, "antigen_response_score", "cytolytic_score")
        assert isinstance(ax, plt.Axes)

    def test_array_mask_positional(self):
        tbl = _score_table(n=5)
        mask = np.array([True, False, False, True, False])
        ax = plot_signature_map(
            tbl, "antigen_response_score", "cytolytic_score", highlight_mask=mask
        )
        assert isinstance(ax, plt.Axes)

    def test_missing_score_column_raises(self):
        with pytest.raises(ValueError):
            plot_signature_map(_score_table(), "nope", "cytolytic_score")

    def test_missing_label_by_raises(self):
        with pytest.raises(ValueError):
            plot_signature_map(
                _score_table(),
                "antigen_response_score",
                "cytolytic_score",
                label_by="absent",
            )

    def test_mask_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            plot_signature_map(
                _score_table(n=5),
                "antigen_response_score",
                "cytolytic_score",
                highlight_mask=np.array([True, False]),
            )
