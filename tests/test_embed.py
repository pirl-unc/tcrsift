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

"""Single-cell embedding (#311).

Pearson residuals → PCA → optional Harmony → UMAP → Leiden, producing the
atlas the annotator (#312) and atlas plots (#315) consume.
"""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.config import EmbedConfig, TCRsiftConfig
from tcrsift.embed import embed_cells

# The embedding needs the `atlas` extra (Leiden clustering). Skip the whole
# module cleanly when it isn't installed so a base `pip install -e .[dev]` still
# runs; CI installs `.[dev,atlas]` so these actually execute.
pytest.importorskip("igraph")
pytest.importorskip("leidenalg")

# harmonypy logs verbosely at INFO; keep the test output clean.
logging.getLogger("harmonypy").setLevel(logging.WARNING)


@pytest.fixture
def two_pop_adata():
    """120 cells across two transcriptional populations and two batches, with a
    receptor gene (TRBV20-1) and a mito gene (MT-CO1) that must be excluded from
    the clustering panel."""
    rng = np.random.default_rng(0)
    n, g = 120, 200
    X = rng.poisson(0.5, size=(n, g)).astype(float)
    X[:60, 0:20] += rng.poisson(8, size=(60, 20))    # population A
    X[60:, 20:40] += rng.poisson(8, size=(60, 20))   # population B
    genes = [f"GENE{i}" for i in range(g)]
    genes[0] = "TRBV20-1"
    genes[1] = "MT-CO1"
    obs = pd.DataFrame(
        {"sample": (["s1", "s2"] * (n // 2))[:n]},
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=obs)


class TestEmbedConfig:
    def test_composed_on_root_config(self):
        c = TCRsiftConfig()
        assert hasattr(c, "embed")
        assert c.embed.enabled is False  # opt-in

    def test_round_trip_through_dict(self):
        c = TCRsiftConfig._from_dict(
            {"embed": {"enabled": True, "n_pcs": 20, "batch_key": "sample"}}
        )
        assert c.embed.enabled and c.embed.n_pcs == 20
        assert c.embed.batch_key == "sample"
        assert "embed" in c.to_dict()

    def test_flat_keys_route_to_embed(self):
        c = TCRsiftConfig._from_dict({"embed": True, "leiden_resolution": 0.5})
        assert c.embed.enabled is True
        assert c.embed.leiden_resolution == 0.5


class TestEmbedCells:
    def test_writes_embedding_and_clusters(self, two_pop_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = embed_cells(two_pop_adata, EmbedConfig(n_top_genes=80, n_pcs=15))
        assert out.obsm["X_umap"].shape == (two_pop_adata.n_obs, 2)
        assert "X_pca" in out.obsm
        assert "leiden" in out.obs
        # Two clean populations → at least two clusters.
        assert out.obs["leiden"].nunique() >= 2

    def test_input_not_mutated(self, two_pop_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            embed_cells(two_pop_adata, EmbedConfig(n_top_genes=80, n_pcs=15))
        assert "X_umap" not in two_pop_adata.obsm
        assert "leiden" not in two_pop_adata.obs

    def test_receptor_and_mito_excluded_from_panel(self, two_pop_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = embed_cells(two_pop_adata, EmbedConfig(n_top_genes=80, n_pcs=15))
        panel = out.uns["embed_panel"]
        assert "TRBV20-1" not in panel
        assert "MT-CO1" not in panel

    def test_harmony_orientation_correct(self, two_pop_adata):
        # The orientation guard must yield (cells × features), not the
        # transposed (features × cells) the broken scanpy wrapper produces.
        pytest.importorskip("harmonypy")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = embed_cells(
                two_pop_adata, EmbedConfig(n_top_genes=80, n_pcs=15, batch_key="sample")
            )
        assert "X_pca_harmony" in out.obsm
        assert out.obsm["X_pca_harmony"].shape[0] == two_pop_adata.n_obs

    def test_informative_genes_panel_override(self, two_pop_adata):
        panel = [f"GENE{i}" for i in range(2, 40)]  # skip receptor/MT at 0,1
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = embed_cells(
                two_pop_adata, EmbedConfig(informative_genes=panel, n_pcs=15)
            )
        assert set(out.uns["embed_panel"]).issubset(set(panel))
        assert out.obs["leiden"].nunique() >= 2

    def test_missing_batch_key_raises(self, two_pop_adata):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with pytest.raises(ValueError, match="batch_key"):
                embed_cells(
                    two_pop_adata, EmbedConfig(n_top_genes=80, batch_key="not_a_col")
                )
