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
    genes[150] = "MT-CO1"       # a mito + a receptor gene the verbatim panel
    genes[151] = "TRBV20-1"     # must keep (embed_cells would have dropped both)
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


class TestFaithfulReembed:
    """The re-embed must run the caller's recipe on the panel VERBATIM — the
    byte-for-byte property a hand-rolled pass needs to be deletable."""

    def test_informative_genes_used_verbatim_not_curated(self, mixed_atlas):
        # A panel of ONLY a mito + a receptor gene: embed_cells drops both and
        # would raise "none present"; the faithful re-embed keeps them, so the
        # call succeeds and the compartment is relabeled.
        _refine(
            mixed_atlas,
            "CD8 T",
            lambda sub: "refined",
            informative_genes=["MT-CO1", "TRBV20-1"],
        )
        assert set(mixed_atlas.obs["phenotype"].iloc[0:100]) == {"refined"}

    def test_deterministic_partition(self, mixed_atlas):
        a = mixed_atlas.copy()
        b = mixed_atlas.copy()
        _refine(a, "CD8 T", _dominant_marker_label, random_state=0)
        _refine(b, "CD8 T", _dominant_marker_label, random_state=0)
        assert list(a.obs["phenotype"]) == list(b.obs["phenotype"])


@pytest.fixture
def gd_compartment():
    """120-cell "CD8 T" compartment holding a genuine γδ sub-pop (high TRDC,
    low αβ-capture, marker set A) and an ambient-αβ CD8 sub-pop (no TRDC, high
    αβ-capture, marker set B); plus a 80-cell untouched myeloid block."""
    rng = np.random.default_rng(0)
    n, g = 200, 200
    X = rng.poisson(0.5, size=(n, g)).astype(float)
    X[0:60, 0:20] += rng.poisson(8, size=(60, 20))       # γδ marker set
    X[60:120, 20:40] += rng.poisson(8, size=(60, 20))    # αβ-CD8 marker set
    X[120:200, 40:60] += rng.poisson(8, size=(80, 20))   # myeloid (excluded)
    genes = [f"GENE{i}" for i in range(g)]
    genes[100] = "TRDC"
    X[:, 100] = 0.0
    X[0:60, 100] = 3.0                                   # TRDC only in γδ cells
    ab = np.zeros(n, dtype=bool)
    ab[0:60] = rng.random(60) < 0.10                     # γδ: sparse αβ-capture
    ab[60:120] = rng.random(60) < 0.85                   # αβ-CD8: high capture
    obs = pd.DataFrame(
        {
            "phenotype": ["CD8 T"] * 120 + ["Myeloid"] * 80,
            "has_ab_contig": ab,
        },
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, var=pd.DataFrame(index=genes), obs=obs)


class TestGdCd8ThroughRefineCluster:
    def test_disaggregates_gd_from_ab_cd8(self, gd_compartment):
        from tcrsift.embed import gd_cd8_relabel

        _refine(
            gd_compartment,
            "CD8 T",
            gd_cd8_relabel(),
            resolution=0.5,
            final_cluster_col="final_cluster",
        )
        pheno = gd_compartment.obs["phenotype"]
        fc = gd_compartment.obs["final_cluster"]
        assert set(pheno.iloc[0:60]) == {"gd T"}
        assert set(fc.iloc[0:60]) == {"gdt"}
        assert set(pheno.iloc[60:120]) == {"CD8 effector/cytotoxic"}
        assert set(fc.iloc[60:120]) == {"gdt_cd8"}
        # Myeloid untouched, no final-cluster tag.
        assert set(pheno.iloc[120:200]) == {"Myeloid"}
        assert fc.iloc[120:200].isna().all()
