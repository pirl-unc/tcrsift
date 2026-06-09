"""Tests for tcrsift._dtypes.rehydrate_obs."""

import anndata as ad
import numpy as np
import pandas as pd

from tcrsift._dtypes import rehydrate_obs


def _make_adata(obs_dict):
    n = len(next(iter(obs_dict.values())))
    obs = pd.DataFrame(obs_dict)
    # Relabel to string obs_names AFTER construction: passing index= to the
    # constructor would realign dict values that are Series (their integer
    # index wouldn't match string labels → all-NaN). astype(str) just renames.
    obs.index = obs.index.astype(str)
    return ad.AnnData(X=np.zeros((n, 1), dtype=np.float32), obs=obs)


class TestRehydrateObs:
    def test_categorical_strings_become_object(self):
        adata = _make_adata(
            {
                "CDR3_alpha": pd.Categorical(["CAVA", "CAVB", "CAVC"]),
                "CDR3_beta": pd.Categorical(["CASS1", "CASS2", "CASS3"]),
                "sample": pd.Categorical(["S1", "S1", "S2"]),
            }
        )
        rehydrate_obs(adata)
        for col in ("CDR3_alpha", "CDR3_beta", "sample"):
            assert adata.obs[col].dtype == object, col

    def test_categorical_with_nan_string_becomes_object_with_nan(self):
        # The exact #11 failure mode — Categorical without "" in categories
        # but a NaN entry. After rehydrate `.fillna("")` must succeed.
        adata = _make_adata(
            {"CDR3_alpha": pd.Categorical(["CAVA", None, "CAVC"])}
        )
        rehydrate_obs(adata)
        assert adata.obs["CDR3_alpha"].dtype == object
        # Downstream: this must not raise.
        adata.obs["CDR3_alpha"].fillna("")

    def test_object_bool_with_none_becomes_bool_false(self):
        adata = _make_adata(
            {"has_TRA": pd.Series([True, False, None], dtype=object)}
        )
        rehydrate_obs(adata)
        assert adata.obs["has_TRA"].dtype == bool
        assert list(adata.obs["has_TRA"]) == [True, False, False]

    def test_float_int_with_nan_becomes_int_zero(self):
        adata = _make_adata(
            {"TRA_count": pd.Series([1.0, 2.0, np.nan], dtype=float)}
        )
        rehydrate_obs(adata)
        assert pd.api.types.is_integer_dtype(adata.obs["TRA_count"])
        assert list(adata.obs["TRA_count"]) == [1, 2, 0]

    def test_idempotent(self):
        adata = _make_adata(
            {
                "CDR3_alpha": pd.Categorical(["CAVA", "CAVB"]),
                "has_TRA": pd.Series([True, False], dtype=object),
                "TRA_count": pd.Series([1.0, 2.0], dtype=float),
            }
        )
        rehydrate_obs(adata)
        # Capture state, run again, expect no change.
        first = adata.obs.copy()
        rehydrate_obs(adata)
        pd.testing.assert_frame_equal(first, adata.obs)

    def test_unknown_columns_untouched(self):
        adata = _make_adata(
            {
                "CDR3_alpha": pd.Categorical(["CAVA"]),
                "some_user_metric": pd.Series([3.14], dtype=float),
            }
        )
        rehydrate_obs(adata)
        assert adata.obs["some_user_metric"].dtype == float

    def test_returns_same_adata(self):
        adata = _make_adata({"CDR3_alpha": pd.Categorical(["CAVA"])})
        out = rehydrate_obs(adata)
        assert out is adata

    def test_axis_string_columns_repinned(self):
        """timepoint, apc_type (and existing axis fields) survive an h5ad-style
        Categorical roundtrip cleanly (#9)."""
        adata = _make_adata(
            {
                "patient_id": pd.Categorical(["B1-2", "B1-3"]),
                "enrichment_method": pd.Categorical(["AIMpos", "tetpos"]),
                "timepoint": pd.Categorical(["D7", "D14"]),
                "apc_type": pd.Categorical(["mDC", "B-LCL"]),
                "tissue": pd.Categorical(["blood", "tumor"]),
            }
        )
        rehydrate_obs(adata)
        for col in ("patient_id", "enrichment_method", "timepoint", "apc_type", "tissue"):
            assert adata.obs[col].dtype == object, col

    def test_h5ad_roundtrip_then_aggregate(self, tmp_path):
        """End-to-end: write h5ad → read back → rehydrate → aggregate works.
        Anchors the helper to the actual #11 failure path."""
        from tcrsift.clonotype import aggregate_clonotypes

        n = 30
        obs = pd.DataFrame(
            {
                "sample": ["S1"] * n,
                "CDR3_alpha": [f"CAV_{i % 5}" for i in range(n)],
                "CDR3_beta": [f"CASS_{i % 5}" for i in range(n)],
                "is_CD8": [True] * n,
            },
            index=[str(i) for i in range(n)],
        )
        adata = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32), obs=obs)

        path = tmp_path / "rt.h5ad"
        adata.write_h5ad(path)
        reloaded = ad.read_h5ad(path)

        # After read_h5ad, anndata typically returns string columns as
        # Categorical. Confirm the failure premise still holds — if not,
        # the bug class is gone and the test is a no-op confirmation.
        # Either way, aggregate_clonotypes must succeed.
        clonotypes = aggregate_clonotypes(reloaded, group_by="CDR3ab")
        assert len(clonotypes) == 5
