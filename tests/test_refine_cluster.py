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

"""General recovery-pass primitive ``refine_cluster`` (#352).

Re-embed one compartment on its own Pearson residuals, sub-cluster it, and hand
each sub-cluster to a caller callback for relabeling — the shared engine behind
hand-rolled γδ/αβ-CD8, DC-recovery and tumor sub-state passes.
"""

from __future__ import annotations

import logging
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from tcrsift.embed import refine_cluster

# Re-embedding needs the `atlas` extra (Leiden clustering).
pytest.importorskip("igraph")
pytest.importorskip("leidenalg")

logging.getLogger("harmonypy").setLevel(logging.WARNING)


@pytest.fixture
def mixed_atlas():
    """200 cells, 200 genes. Cells 0-99 are one first-pass compartment
    ("CD8 T") that actually holds two clean transcriptional sub-populations —
    A (marker GENE0..19) and B (marker GENE20..39) — the way a cytotoxic
    cluster co-embeds genuine γδ with ambient-αβ CD8. Cells 100-199 are a
    separate "Myeloid" compartment that must be left untouched."""
    rng = np.random.default_rng(0)
    n, g = 200, 200
    X = rng.poisson(0.5, size=(n, g)).astype(float)
    X[0:50, 0:20] += rng.poisson(8, size=(50, 20))      # CD8 sub-pop A
    X[50:100, 20:40] += rng.poisson(8, size=(50, 20))   # CD8 sub-pop B
    X[100:200, 40:60] += rng.poisson(8, size=(100, 20))  # Myeloid
    genes = [f"GENE{i}" for i in range(g)]
    obs = pd.DataFrame(
        {"phenotype": ["CD8 T"] * 100 + ["Myeloid"] * 100},
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=obs)


def _dominant_marker_label(sub):
    """Toy relabel_fn: name a sub-cluster by whichever marker set dominates."""
    a = float(np.asarray(sub[:, "GENE0"].X).mean())
    b = float(np.asarray(sub[:, "GENE20"].X).mean())
    return "CD8 alpha-beta" if a >= b else "gamma-delta T"


def _refine(adata, mask, relabel_fn, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return refine_cluster(adata, mask, relabel_fn=relabel_fn, **kw)


class TestRefineCluster:
    def test_splits_compartment_into_subpopulations(self, mixed_atlas):
        _refine(mixed_atlas, "CD8 T", _dominant_marker_label)
        pheno = mixed_atlas.obs["phenotype"]
        # Sub-pop A → one label, sub-pop B → the other, cleanly split at 50.
        assert set(pheno.iloc[0:50]) == {"CD8 alpha-beta"}
        assert set(pheno.iloc[50:100]) == {"gamma-delta T"}

    def test_leaves_other_compartments_untouched(self, mixed_atlas):
        _refine(mixed_atlas, "CD8 T", _dominant_marker_label)
        assert set(mixed_atlas.obs["phenotype"].iloc[100:200]) == {"Myeloid"}

    def test_mutates_in_place_and_returns_same_object(self, mixed_atlas):
        out = _refine(mixed_atlas, "CD8 T", _dominant_marker_label)
        assert out is mixed_atlas

    def test_boolean_mask(self, mixed_atlas):
        mask = np.zeros(mixed_atlas.n_obs, dtype=bool)
        mask[0:100] = True
        _refine(mixed_atlas, mask, _dominant_marker_label)
        assert set(mixed_atlas.obs["phenotype"].iloc[100:200]) == {"Myeloid"}
        assert "gamma-delta T" in set(mixed_atlas.obs["phenotype"])

    def test_callable_mask(self, mixed_atlas):
        _refine(
            mixed_atlas,
            lambda a: (a.obs["phenotype"] == "CD8 T").to_numpy(),
            _dominant_marker_label,
        )
        assert set(mixed_atlas.obs["phenotype"].iloc[50:100]) == {"gamma-delta T"}

    def test_none_return_leaves_labels_untouched(self, mixed_atlas):
        _refine(mixed_atlas, "CD8 T", lambda sub: None)
        assert set(mixed_atlas.obs["phenotype"].iloc[0:100]) == {"CD8 T"}

    def test_final_cluster_tag_stamped(self, mixed_atlas):
        def relabel(sub):
            label = _dominant_marker_label(sub)
            return label, f"cd8::{label}"

        _refine(
            mixed_atlas, "CD8 T", relabel, final_cluster_col="final_cluster"
        )
        finals = mixed_atlas.obs["final_cluster"]
        assert set(finals.iloc[0:50]) == {"cd8::CD8 alpha-beta"}
        assert set(finals.iloc[50:100]) == {"cd8::gamma-delta T"}
        # Untouched compartment keeps no tag.
        assert finals.iloc[100:200].isna().all()

    def test_label_only_pair_still_stamps_final_cluster(self, mixed_atlas):
        # (None, tag) → leave the label, still record the final-cluster tag.
        _refine(
            mixed_atlas,
            "CD8 T",
            lambda sub: (None, "cd8-refined"),
            final_cluster_col="final_cluster",
        )
        assert set(mixed_atlas.obs["phenotype"].iloc[0:100]) == {"CD8 T"}
        assert set(mixed_atlas.obs["final_cluster"].iloc[0:100]) == {"cd8-refined"}

    def test_too_few_cells_skips_and_warns(self, mixed_atlas, caplog):
        before = mixed_atlas.obs["phenotype"].copy()
        with caplog.at_level(logging.WARNING, logger="tcrsift.embed"):
            _refine(mixed_atlas, "CD8 T", _dominant_marker_label, min_cells=500)
        assert (mixed_atlas.obs["phenotype"] == before).all()
        assert "skipping re-embedding" in caplog.text

    def test_string_mask_requires_label_col(self, mixed_atlas):
        with pytest.raises(ValueError, match="not in adata.obs"):
            _refine(mixed_atlas, "CD8 T", _dominant_marker_label, label_col="missing")

    def test_mask_length_mismatch_raises(self, mixed_atlas):
        with pytest.raises(ValueError, match="mask length"):
            _refine(mixed_atlas, np.ones(5, dtype=bool), _dominant_marker_label)

    def test_creates_label_col_when_absent(self, mixed_atlas):
        del mixed_atlas.obs["phenotype"]
        mask = np.zeros(mixed_atlas.n_obs, dtype=bool)
        mask[0:100] = True
        _refine(mixed_atlas, mask, _dominant_marker_label, label_col="phenotype")
        assert set(mixed_atlas.obs["phenotype"].iloc[0:50]) == {"CD8 alpha-beta"}
        # Cells outside the mask were created but never relabeled.
        assert mixed_atlas.obs["phenotype"].iloc[100:200].isna().all()
