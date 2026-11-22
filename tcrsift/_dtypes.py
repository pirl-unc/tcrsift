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

"""
Internal helpers for normalizing AnnData obs dtypes after h5ad round-trips.

anndata serialization is not dtype-preserving for object columns: string
columns are returned as ``Categorical``, bool columns as object/Categorical,
and integer columns as ``float64`` whenever NaN values were introduced. This
breaks downstream code that assumes the in-memory original dtype:

  - Categorical CDR3 columns reject ``.fillna("")`` and string concatenation
    (issue #11)
  - object/None bool columns can't be serialized back to h5ad (issue #3)
  - object columns of all-None values trip h5py during write (issue #5)

``rehydrate_obs(adata)`` re-pins the dtypes for the columns tcrsift owns,
to be called at pipeline-stage entry. Idempotent and cheap.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

# Columns tcrsift writes as plain Python strings. After an h5ad round-trip
# anndata returns them as pandas.Categorical, which rejects setitem with
# new categories.
_STRING_COLS: frozenset[str] = frozenset(
    {
        "CDR3_alpha",
        "CDR3_beta",
        "CDR3ab",
        "sample",
        "patient_id",
        "enrichment_method",
        "timepoint",
        "apc_type",
        "antigen_type",
        "antigen_description",
        "antigen_name",
        "antigen_sequence",
        "epitope_sequence",
        "mhc_allele",
        "antigen_names",
        "antigen_sequences",
        "epitope_sequences",
        "source",
        "expected_tcell_type",
        "Tcell_type",
        "tissue",
        "TRA_1_v_gene",
        "TRA_1_d_gene",
        "TRA_1_j_gene",
        "TRA_1_c_gene",
        "TRB_1_v_gene",
        "TRB_1_d_gene",
        "TRB_1_j_gene",
        "TRB_1_c_gene",
        "TRA_1_contig_id",
        "TRB_1_contig_id",
    }
)

# Columns tcrsift writes as bool.
_BOOL_COLS: frozenset[str] = frozenset(
    {
        "has_TRA",
        "has_TRB",
        "has_both_chains",
        "multi_TRA",
        "multi_TRB",
        "multi_chain",
        "is_complete_clone",
        "is_doublet",
        "is_CD4",
        "is_CD8",
        "TRA_pass_umi",
        "TRB_pass_umi",
    }
)

# Columns tcrsift writes as int.
_INT_COLS: frozenset[str] = frozenset(
    {
        "TRA_count",
        "TRB_count",
    }
)


def _to_clean_object_series(s: pd.Series) -> pd.Series:
    """Convert Categorical / Arrow StringDtype / nullable object columns to plain
    ``object`` with ``pd.NA`` / ``np.nan`` normalized to ``None`` for h5py."""
    arr = s.astype(object)
    return arr.where(pd.notna(arr), None)


def rehydrate_obs(adata: ad.AnnData) -> ad.AnnData:
    """Normalize known obs columns (and any StringDtype columns) to the dtypes
    tcrsift and h5py expect.

    Mutates ``adata.obs`` in place; returns ``adata`` for chaining. Apply at
    every pipeline-stage entry and prior to ``write_h5ad`` to immunize against
    h5ad-roundtrip and Pandas 2.x/3.x ``StringDtype`` / ``pd.NA`` drift.

    Idempotent: re-running on already-normalized data is a no-op.
    """
    obs = adata.obs
    columns = set(obs.columns)

    for col in list(obs.columns):
        s = obs[col]
        if col in _STRING_COLS or isinstance(s.dtype, pd.StringDtype):
            if not pd.api.types.is_object_dtype(s) or s.isna().any():
                obs[col] = _to_clean_object_series(s)
        elif isinstance(s.dtype, pd.CategoricalDtype) and isinstance(
            s.cat.categories.dtype, pd.StringDtype
        ):
            obs[col] = s.cat.rename_categories(s.cat.categories.astype(object))
        elif pd.api.types.is_object_dtype(s) and s.isna().any():
            # Ensure pd.NA inside arbitrary object columns becomes None/NaN
            obs[col] = s.where(pd.notna(s), None)

    if isinstance(obs.index.dtype, pd.StringDtype):
        obs.index = obs.index.astype(object)

    var = getattr(adata, "var", None)
    if var is not None:
        for col in list(var.columns):
            if isinstance(var[col].dtype, pd.StringDtype):
                var[col] = _to_clean_object_series(var[col])
        if isinstance(var.index.dtype, pd.StringDtype):
            var.index = var.index.astype(object)

    for col in _BOOL_COLS & columns:
        if obs[col].dtype != np.dtype(bool):
            values = obs[col].astype(object).to_numpy()
            mask = pd.isna(values)
            obs[col] = np.where(mask, False, values).astype(bool)

    for col in _INT_COLS & columns:
        if obs[col].dtype != np.dtype(np.int64):
            values = pd.to_numeric(obs[col], errors="coerce").to_numpy()
            mask = pd.isna(values)
            obs[col] = np.where(mask, 0, values).astype(np.int64)

    return adata
