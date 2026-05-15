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

"""Tests for the Pgen / publicness visualizations (#58)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from tcrsift.pgen import annotate_publicness  # noqa: E402
from tcrsift.plots import (  # noqa: E402
    plot_pgen_distribution,
    plot_publicness_vs_match_score,
)


def _sample_clonotypes(n=20, seed=0):
    rng = np.random.default_rng(seed)
    cdr3_lens = rng.integers(11, 18, size=n)
    cdr3s = [
        "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=L)) for L in cdr3_lens
    ]
    df = pd.DataFrame({
        "CDR3_beta": cdr3s,
        "beta_v_gene": rng.choice(
            ["TRBV20-1", "TRBV9", "TRBV28", "TRBV6-5", "TRBV15"],
            size=n,
        ),
        "beta_j_gene": rng.choice(
            ["TRBJ2-1", "TRBJ1-1", "TRBJ2-7"], size=n
        ),
        "tier": rng.choice(["tier1", "tier2", "tier3"], size=n),
        "n_db_matches": rng.integers(0, 50, size=n),
    })
    return annotate_publicness(df)


class TestPlotPgenDistribution:
    def test_emits_file(self, tmp_path):
        df = _sample_clonotypes()
        out = plot_pgen_distribution(df, tmp_path / "pgen.png")
        assert out is not None and out.exists() and out.stat().st_size > 0

    def test_skips_when_pgen_column_missing(self, tmp_path):
        df = pd.DataFrame({"CDR3_beta": ["CASS"]})
        out = plot_pgen_distribution(df, tmp_path / "pgen.png")
        assert out is None

    def test_skips_on_all_nan(self, tmp_path):
        df = pd.DataFrame({"log10_pgen": [np.nan, np.nan]})
        out = plot_pgen_distribution(df, tmp_path / "pgen.png")
        assert out is None

    def test_group_facet_works(self, tmp_path):
        df = _sample_clonotypes()
        out = plot_pgen_distribution(
            df, tmp_path / "pgen_by_tier.png", group_col="tier"
        )
        assert out.exists()


class TestPlotPublicnessVsMatchScore:
    def test_emits_file(self, tmp_path):
        df = _sample_clonotypes()
        out = plot_publicness_vs_match_score(
            df, tmp_path / "publicness_scatter.png"
        )
        assert out is not None and out.exists()

    def test_skips_when_columns_missing(self, tmp_path):
        df = pd.DataFrame({"CDR3_beta": ["CASS"]})
        out = plot_publicness_vs_match_score(
            df, tmp_path / "scatter.png"
        )
        assert out is None

    def test_works_without_tier_column(self, tmp_path):
        df = _sample_clonotypes().drop(columns=["tier"])
        out = plot_publicness_vs_match_score(
            df, tmp_path / "no_tier.png"
        )
        assert out.exists()
