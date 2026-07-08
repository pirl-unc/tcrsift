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

"""γδ / ambient-αβ CD8 relabel_fn logic (#351).

Pure per-sub-cluster decision — no re-embedding — so it runs in a base install
without the `atlas` extra. The h37 rule: call γδ iff mean TRDC >= 0.9 (log1p
CP10K) AND mean αβ-capture < 0.25.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from tcrsift.embed import gd_cd8_relabel


def _subcluster(trdc_level, ab_fraction, n=40, *, sparse_x=False, extra_gene=True):
    """A sub-cluster with a uniform TRDC level and a given αβ-capture fraction."""
    genes = ["TRDC"] + (["CD8A"] if extra_gene else [])
    X = np.zeros((n, len(genes)), dtype=float)
    X[:, 0] = trdc_level
    if sparse_x:
        X = sparse.csr_matrix(X)
    n_true = int(round(ab_fraction * n))
    ab = np.array([True] * n_true + [False] * (n - n_true))
    obs = pd.DataFrame(
        {"has_ab_contig": ab}, index=[f"c{i}" for i in range(n)]
    )
    return ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=obs)


class TestGdCd8Relabel:
    def test_gd_when_high_trdc_and_low_capture(self):
        fn = gd_cd8_relabel()
        assert fn(_subcluster(1.0, 0.10)) == ("gd T", "gdt")

    def test_cd8_when_capture_high(self):
        # TRDC clears the bar but αβ-capture doesn't → not γδ.
        fn = gd_cd8_relabel()
        assert fn(_subcluster(1.0, 0.80)) == ("CD8 effector/cytotoxic", "gdt_cd8")

    def test_cd8_when_trdc_low(self):
        # αβ-capture is low but TRDC is below 0.9 → not γδ (BOTH required).
        fn = gd_cd8_relabel()
        assert fn(_subcluster(0.5, 0.05)) == ("CD8 effector/cytotoxic", "gdt_cd8")

    def test_both_conditions_are_required(self):
        fn = gd_cd8_relabel()
        assert fn(_subcluster(0.5, 0.80))[0] == "CD8 effector/cytotoxic"

    def test_boundary_is_inclusive_on_trdc_exclusive_on_capture(self):
        fn = gd_cd8_relabel()
        # TRDC exactly 0.9 (>=) with capture just under 0.25 → γδ.
        assert fn(_subcluster(0.9, 0.20))[0] == "gd T"
        # capture exactly 0.25 (not < 0.25) → not γδ.
        assert fn(_subcluster(1.0, 0.25))[0] == "CD8 effector/cytotoxic"

    def test_sparse_matrix_supported(self):
        fn = gd_cd8_relabel()
        assert fn(_subcluster(1.0, 0.10, sparse_x=True)) == ("gd T", "gdt")

    def test_ambiguous_capture_warns_but_still_calls(self, caplog):
        fn = gd_cd8_relabel()
        with caplog.at_level(logging.WARNING, logger="tcrsift.embed"):
            label = fn(_subcluster(1.0, 0.40))  # 0.25 < 0.40 < 0.55
        # Ambiguous band → call is made (capture >= 0.25 so αβ-CD8), warning logged.
        assert label[0] == "CD8 effector/cytotoxic"
        assert "between the" in caplog.text

    def test_no_warning_outside_ambiguous_band(self, caplog):
        fn = gd_cd8_relabel()
        with caplog.at_level(logging.WARNING, logger="tcrsift.embed"):
            fn(_subcluster(1.0, 0.05))
        assert caplog.text == ""

    def test_missing_ab_col_raises(self):
        fn = gd_cd8_relabel()
        sub = _subcluster(1.0, 0.10)
        del sub.obs["has_ab_contig"]
        with pytest.raises(KeyError, match="has_ab_contig"):
            fn(sub)

    def test_custom_thresholds_and_labels(self):
        fn = gd_cd8_relabel(
            gene="TRGC1",
            ab_col="paired",
            trdc_min=0.5,
            ab_max=0.5,
            gd_label="γδ",
            cd8_label="CD8",
            gd_cluster="g",
            cd8_cluster="c",
        )
        X = np.full((10, 1), 0.6)
        obs = pd.DataFrame({"paired": [True] * 2 + [False] * 8},  # 0.2 < 0.5
                           index=[f"c{i}" for i in range(10)])
        sub = ad.AnnData(X=X, var=pd.DataFrame(index=["TRGC1"]), obs=obs)
        assert fn(sub) == ("γδ", "g")
