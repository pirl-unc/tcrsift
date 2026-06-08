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

"""Tests for clonal-expansion QC metrics + unexpanded warning (#161)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift.qc import cdr3_anchor_integrity, clonal_expansion_metrics


def _long(rows):
    return pd.DataFrame(rows, columns=["sample", "CDR3ab", "cells"])


class TestClonalExpansionMetrics:
    def test_expanded_sample_not_flagged(self):
        # One dominant clone (B1-2-like): top1 ~24%, low effective fraction.
        rows = [("S", f"c{i}", 1543 if i == 0 else 3) for i in range(2000)]
        m = clonal_expansion_metrics(_long(rows))
        r = m.iloc[0]
        assert r["top1_clone_fraction"] > 0.2
        assert not r["unexpanded"]
        assert r["warning"] == ""

    def test_unexpanded_sample_flagged(self):
        # Near-uniform (B1-4-like): top clone 7 cells, rest singletons.
        rows = [("S", f"c{i}", 7 if i < 3 else 1) for i in range(7000)]
        m = clonal_expansion_metrics(_long(rows))
        r = m.iloc[0]
        assert r["top1_clone_fraction"] < 0.01
        assert r["effective_clones"] / r["observed_clones"] > 0.8
        assert r["unexpanded"]
        assert "near-polyclonal" in r["warning"]

    def test_metrics_values(self):
        # Two clones, 75/25 split → known Shannon.
        m = clonal_expansion_metrics(_long([("S", "a", 75), ("S", "b", 25)]))
        r = m.iloc[0]
        assert r["n_cells"] == 100
        assert r["observed_clones"] == 2
        assert r["top1_clone_fraction"] == 0.75
        assert r["gini_simpson"] == pytest.approx(1 - (0.75**2 + 0.25**2))
        h = -(0.75 * np.log(0.75) + 0.25 * np.log(0.25))
        assert r["effective_clones"] == pytest.approx(np.exp(h))

    def test_fraction_cells_in_clones_ge_n(self):
        rows = [("S", "big", 50), ("S", "mid", 10)] + [
            ("S", f"s{i}", 1) for i in range(40)
        ]
        m = clonal_expansion_metrics(_long(rows), ge_n=10)
        r = m.iloc[0]
        # 50 + 10 of 100 cells are in clones with >= 10 cells
        assert r["fraction_cells_in_clones_ge_10"] == pytest.approx(0.6)

    def test_per_sample_rows(self):
        rows = [("A", "x", 100), ("B", "y", 50), ("B", "z", 50)]
        m = clonal_expansion_metrics(_long(rows))
        assert set(m["sample"]) == {"A", "B"}
        # A is monoclonal → clonality 1
        assert m.set_index("sample").loc["A", "clonality"] == pytest.approx(1.0)

    def test_aggregates_clone_across_rows(self):
        # Same clone appearing twice in a sample sums its cells.
        m = clonal_expansion_metrics(_long([("S", "a", 30), ("S", "a", 70)]))
        assert m.iloc[0]["observed_clones"] == 1
        assert m.iloc[0]["n_cells"] == 100

    def test_missing_group_col_raises(self):
        with pytest.raises(ValueError, match="missing"):
            clonal_expansion_metrics(_long([("S", "a", 1)]), group_col="donor")

    def test_donor_grouping(self):
        df = pd.DataFrame({
            "donor": ["d1", "d1", "d2"],
            "CDR3ab": ["a", "b", "c"],
            "cells": [10, 5, 20],
        })
        m = clonal_expansion_metrics(df, group_col="donor")
        assert set(m["donor"]) == {"d1", "d2"}


class TestCdr3AnchorIntegrity:
    def test_fraction_ending_fw(self):
        df = pd.DataFrame({
            "CDR3_beta": ["CASSLGTGELFF", "CASSPGQGAYEQYF", "CASSENDSINC"],
        })
        # 2 of 3 end in F/W
        out = cdr3_anchor_integrity(df, cdr3_cols=("CDR3_beta",))
        assert out["CDR3_beta"] == pytest.approx(2 / 3)

    def test_missing_column_skipped(self):
        out = cdr3_anchor_integrity(pd.DataFrame({"x": [1]}),
                                    cdr3_cols=("CDR3_beta",))
        assert out == {}

    def test_empty_cdr3s_nan(self):
        out = cdr3_anchor_integrity(pd.DataFrame({"CDR3_beta": ["", None]}),
                                    cdr3_cols=("CDR3_beta",))
        assert np.isnan(out["CDR3_beta"])


def _adata(obs, n_var=3):
    import anndata as ad
    return ad.AnnData(np.zeros((len(obs), n_var)), obs=obs)


def _barcodes(n, suffix, seed):
    rng = np.random.default_rng(seed)
    return [f"{''.join(rng.choice(list('ACGT'), 16))}-{suffix}" for _ in range(n)]


class TestSampleIntegrityQC:
    def _data(self):
        bc_a = _barcodes(100, 1, 0)
        bc_b = _barcodes(100, 2, 1)
        bc_c = [b.rsplit("-", 1)[0] + "-3" for b in bc_a]  # C = copy of A
        obs = pd.DataFrame({
            "sample": ["A"] * 100 + ["B"] * 100 + ["C"] * 100,
            "CDR3_beta": (["CASSLF"] * 80 + [None] * 20)
            + (["CASSLF"] * 95 + [None] * 5)
            + (["CASSLF"] * 10 + [None] * 90),
        }, index=bc_a + bc_b + bc_c)
        return _adata(obs)

    def test_gex_vdj_overlap_fraction(self):
        from tcrsift.qc import gex_vdj_overlap
        out = gex_vdj_overlap(self._data()).set_index("sample")
        assert out.loc["A", "gex_vdj_overlap_fraction"] == pytest.approx(0.80)
        assert out.loc["C", "n_with_vdj"] == 10

    def test_gex_vdj_overlap_duplicate_index_across_samples(self):
        # 10x reuses one barcode whitelist and load_samples does NOT uniquify
        # obs_names, so the same label appears in multiple samples. Per-sample
        # counts must not bleed across samples and the fraction must stay <= 1.0
        # (the old label-based .loc pulled rows from every sample, F7).
        from tcrsift.qc import gex_vdj_overlap

        obs = pd.DataFrame(
            {
                "sample": ["A", "A", "B", "B"],
                "CDR3_beta": ["CASSLF", None, "CASSLF", "CASSLF"],
            },
            index=["BC1-1", "BC2-1", "BC1-1", "BC3-1"],  # BC1-1 in both A and B
        )
        out = gex_vdj_overlap(_adata(obs)).set_index("sample")
        assert out.loc["A", "n_cells"] == 2
        assert out.loc["A", "n_with_vdj"] == 1  # only A's own row, not B's BC1-1
        assert out.loc["A", "gex_vdj_overlap_fraction"] == pytest.approx(0.5)
        assert out.loc["B", "n_with_vdj"] == 2
        assert (out["gex_vdj_overlap_fraction"] <= 1.0).all()

    def test_duplicate_input_high_jaccard(self):
        from tcrsift.qc import inter_sample_barcode_jaccard
        jac = inter_sample_barcode_jaccard(self._data())
        # A and C share all base barcodes; A and B are independent.
        assert jac.loc["A", "C"] == pytest.approx(1.0)
        assert jac.loc["A", "B"] == pytest.approx(0.0)
        assert (np.diag(jac.values) == 1.0).all()

    def test_integrity_warnings(self):
        from tcrsift.qc import sample_integrity_qc
        q = sample_integrity_qc(self._data()).set_index("sample")
        # C: low GEX↔VDJ overlap AND duplicate-of-A
        assert "low GEX↔VDJ overlap" in q.loc["C", "warning"]
        assert "duplicate" in q.loc["A", "warning"]  # A flagged as overlapping C
        assert q.loc["B", "warning"] == ""           # B is clean

    def test_no_sample_col_jaccard_empty(self):
        from tcrsift.qc import inter_sample_barcode_jaccard
        obs = pd.DataFrame({"x": [1, 2]}, index=["AA-1", "BB-1"])
        assert inter_sample_barcode_jaccard(_adata(obs)).empty

    def test_missing_vdj_cols_zero_overlap(self):
        from tcrsift.qc import gex_vdj_overlap
        obs = pd.DataFrame({"sample": ["A", "A"]}, index=["AA-1", "BB-1"])
        out = gex_vdj_overlap(_adata(obs))
        assert out.iloc[0]["gex_vdj_overlap_fraction"] == 0.0

    def test_read_cellranger_metrics_missing_returns_empty(self):
        from tcrsift.qc import read_cellranger_metrics
        assert read_cellranger_metrics("/nonexistent/path") == {}
        assert read_cellranger_metrics(None) == {}
