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

"""Centralized DataFrame/Series string, CDR3, clone-key, and boolean normalization.

Provides dtype-safe helpers across Pandas 2 ``object``, ``Categorical``, and
Pandas 3 Arrow/nullable ``StringDtype`` columns so pipeline stages do not
maintain diverging missing-token sets or partial single-chain ``CDR3ab`` guards.
"""

from __future__ import annotations

import pandas as pd

MISSING_SCALAR_TOKENS: frozenset[str] = frozenset({"", "nan", "none", "<na>", "null"})
MISSING_CLONE_TOKENS: frozenset[str] = frozenset({"", "_", "nan", "none", "<na>", "null"})
TRUTHY_STR_TOKENS: frozenset[str] = frozenset({"true", "t", "1", "yes", "y"})


def clean_str_series(
    series: pd.Series,
    *,
    missing_tokens: frozenset[str] = MISSING_SCALAR_TOKENS,
) -> pd.Series:
    """Return a clean ``object``-dtype string Series with missing tokens normalized to ``""``."""
    s = series.astype(object) if isinstance(series.dtype, pd.CategoricalDtype) else series
    out = s.astype("string").fillna("").str.strip()
    if missing_tokens:
        out = out.mask(out.str.lower().isin(missing_tokens), "")
    return out.astype(object)


def clean_cdr3_series(
    df_or_series: pd.DataFrame | pd.Series,
    col: str | None = None,
) -> pd.Series:
    """Extract and normalize a single-chain CDR3 Series to ``object`` strings (``""`` if missing)."""
    if isinstance(df_or_series, pd.DataFrame):
        if col is None or col not in df_or_series.columns:
            return pd.Series("", index=df_or_series.index, dtype=object)
        series = df_or_series[col]
    else:
        series = df_or_series
    return clean_str_series(series, missing_tokens=MISSING_CLONE_TOKENS)


def paired_cdr3ab_key(
    df: pd.DataFrame,
    alpha_col: str = "CDR3_alpha",
    beta_col: str = "CDR3_beta",
) -> tuple[pd.Series, pd.Series]:
    """Build ``(cdr3ab_key, valid_mask)`` requiring both α and β CDR3s to be non-empty."""
    alpha = clean_cdr3_series(df, alpha_col)
    beta = clean_cdr3_series(df, beta_col)
    valid = (alpha != "") & (beta != "")
    key = (alpha + "_" + beta).where(valid, "")
    return key.astype(object), valid


def valid_clone_mask(
    series: pd.Series,
    *,
    clone_col: str = "CDR3ab",
) -> pd.Series:
    """Boolean mask of non-empty clone identifiers, rejecting single-chain ``_beta`` / ``alpha_`` keys."""
    cleaned = clean_str_series(series, missing_tokens=MISSING_CLONE_TOKENS)
    valid = series.notna() & (cleaned != "")
    if clone_col == "CDR3ab":
        s_str = cleaned.astype("string")
        valid &= ~s_str.str.startswith("_", na=True) & ~s_str.str.endswith("_", na=True)
    return valid.astype(bool)


def truthy_flag_series(series: pd.Series) -> pd.Series:
    """Coerce bool, nullable boolean, numeric (1/0/1.0/0.0), or stringified flag Series to ``bool``."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    s = series.astype(object) if isinstance(series.dtype, pd.CategoricalDtype) else series
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce").fillna(0).ne(0).astype(bool)
    lowered = s.astype("string").fillna("").str.strip().str.lower()
    num_truthy = pd.to_numeric(lowered, errors="coerce").fillna(0).ne(0)
    return pd.Series(
        (lowered.isin(TRUTHY_STR_TOKENS) | num_truthy).to_numpy(dtype=bool),
        index=series.index,
        dtype=bool,
    )
