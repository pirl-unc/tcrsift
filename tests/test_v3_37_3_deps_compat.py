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

"""Regression tests for v3.37.3 dependency, Pandas 3, and harmonypy 2.0 fixes."""

from __future__ import annotations

import builtins
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.annotate import _match_viral_species_lower
from tcrsift.config import EmbedConfig
from tcrsift.embed import embed_cells, gd_cd8_relabel
from tcrsift.pgen_models import _bundled_observed_beta, train_model


def test_viral_species_matcher_works_on_pyarrow_string_dtype():
    pytest.importorskip("pyarrow")
    series = pd.Series(
        [
            "influenza a virus",
            "flu a",
            "haemophilus influenzae",
            "parainfluenza",
            "homo sapiens",
            "cytomegalovirus",
        ],
        dtype="string[pyarrow]",
    )
    matched = _match_viral_species_lower(series)
    assert list(matched) == [True, True, False, False, False, True]


def test_embed_cells_harmony_small_n_obs_under_60():
    pytest.importorskip("harmonypy")
    pytest.importorskip("igraph")
    rng = np.random.default_rng(7)
    n_obs, n_vars = 40, 25
    counts = rng.poisson(3.0, size=(n_obs, n_vars)).astype(np.float32)
    obs = pd.DataFrame(
        {"sample": ["S1"] * 20 + ["S2"] * 20},
        index=[f"c{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_vars)])
    adata = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    out = embed_cells(
        adata,
        EmbedConfig(enabled=True, n_pcs=5, n_neighbors=5, batch_key="sample", seed=0),
    )
    assert "X_pca_harmony" in out.obsm
    assert out.obsm["X_pca_harmony"].shape[0] == n_obs
    assert "leiden" in out.obs.columns


def test_gd_cd8_relabel_tolerates_na_in_ab_col():
    obs = pd.DataFrame(
        {"has_ab": pd.Series([False, False, pd.NA], dtype="boolean")},
        index=["c0", "c1", "c2"],
    )
    var = pd.DataFrame(index=["TRDC", "CD8A"])
    sub = ad.AnnData(X=np.array([[2.0, 0.1], [2.5, 0.0], [1.8, 0.2]], dtype=np.float32), obs=obs, var=var)
    relabel = gd_cd8_relabel(ab_col="has_ab")
    label, cluster = relabel(sub)
    assert "γδ" in label or "gd" in cluster.lower()


def test_bundled_observed_beta_returns_none_when_tcrpeg_missing(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "tcrpeg":
            raise ImportError("No module named 'tcrpeg'")
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "tcrpeg", raising=False)
    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert _bundled_observed_beta() is None


def test_train_model_kmer_gene_drops_nan_seq_rows_before_reading_vj(tmp_path, monkeypatch):
    monkeypatch.setenv("TCRSIFT_DATA_DIR", str(tmp_path / "cache"))
    rep_csv = tmp_path / "rep.csv"
    pd.DataFrame(
        [
            {"seq": "CASSLGQAYEQYF", "v": "TRBV7-2", "j": "TRBJ2-7"},
            {"seq": "CASSIRSSYEQYF", "v": "TRBV19", "j": "TRBJ2-7"},
            {"seq": np.nan, "v": "TRBV99-FAKE", "j": "TRBJ99-FAKE"},
        ]
    ).to_csv(rep_csv, index=False)

    out_path = train_model("beta", backend="kmer_gene", role="ppost", url=str(rep_csv))
    assert out_path.exists()
    data = np.load(out_path, allow_pickle=True)
    v_keys = set(str(x) for x in data["v_keys"]) if "v_keys" in data else set()
    assert "TRBV99-FAKE" not in v_keys
