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

"""Signature-based candidate shortlisting (#75).

After the tier cascade has labeled clones at varying FDR thresholds,
the typical follow-up is a shortlist:

  Selected = tier1 ∪ tier2 ∪ top-N-per-signature(tier3 ∪ tier4 ∪ tier5)

For each named signature (antigen-response, cytolytic, tumor-reactive
by default), take the top-N clones from the tier3+ pool by descending
signature score. This pulls in clones that didn't make the strict-FDR
cohort cut but have phenotype evidence we don't want to discard.

The function here is pure: it augments the input frame with boolean
tracking columns describing how each clone made the shortlist, so a
reviewer can audit *why* a clone is in the final list. Callers
filter/export downstream.
"""

from __future__ import annotations

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def select_candidates(
    clonotypes: pd.DataFrame,
    *,
    signatures: tuple[str, ...] = (
        "antigen_response", "cytolytic", "tumor_reactive",
    ),
    top_n: int = 3,
    tier_col: str = "tier",
    score_col_template: str = "signature_{name}",
    tier1_label: str = "tier1",
    tier2_label: str = "tier2",
    tier3plus_labels: tuple[str, ...] = ("tier3", "tier4", "tier5"),
) -> pd.DataFrame:
    """Augment ``clonotypes`` with shortlist tracking columns.

    For each clone, computes:

    - ``tier1_or_tier2`` — bool. True if ``tier`` ∈ {tier1, tier2}.
    - ``tier3plus_top_by_{name}_signature`` — bool, one per signature.
      True if the clone is in the top ``top_n`` clones of the tier3+
      pool when ranked by ``signature_{name}`` descending.
    - ``is_selected`` — bool. The OR of all the above.

    Does **not** filter the frame. The caller can export the
    ``is_selected == True`` subset to ``candidate_shortlist.csv``.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotypes frame from the standard tcrsift pipeline. Must
        contain ``tier_col`` plus the per-signature score columns
        referenced by ``score_col_template``.
    signatures : tuple[str, ...]
        Names of signatures to rank against. Each must have a column
        named ``score_col_template.format(name=name)``.
    top_n : int
        Picks per signature from the tier3+ pool.
    tier_col : str
        Column holding the per-clone tier label.
    score_col_template : str
        ``str.format``-compatible template that maps a signature name
        to its score column.
    tier1_label, tier2_label, tier3plus_labels : str
        Tier labels in the frame. Defaults match
        :data:`tcrsift.filter.DEFAULT_TIER_DEFINITIONS`.

    Returns
    -------
    pd.DataFrame
        Copy of ``clonotypes`` with the new boolean columns appended.
        Skipped signatures (missing score column) produce ``False`` in
        their column and emit a warning.
    """
    df = clonotypes.copy()

    if tier_col not in df.columns:
        raise ValueError(
            f"select_candidates: missing tier column {tier_col!r}; "
            "run the filter step first"
        )

    tier_str = df[tier_col].astype(str)
    df["tier1_or_tier2"] = tier_str.isin({tier1_label, tier2_label})
    pool_mask = tier_str.isin(set(tier3plus_labels))

    sig_columns: list[str] = []
    for sig in signatures:
        col = score_col_template.format(name=sig)
        track_col = f"tier3plus_top_by_{sig}_signature"
        sig_columns.append(track_col)
        if col not in df.columns:
            logger.warning(
                f"select_candidates: missing score column {col!r} for "
                f"signature {sig!r}; setting {track_col!r} to False"
            )
            df[track_col] = False
            continue
        # Rank the tier3+ subset by descending score; ties broken by
        # original order. NaN scores sort to the bottom.
        pool = df.loc[pool_mask, col]
        if pool.empty or pool.dropna().empty:
            df[track_col] = False
            continue
        top_indices = (
            pool.rank(ascending=False, method="first", na_option="bottom")
            .nsmallest(top_n)
            .index
        )
        flags = pd.Series(False, index=df.index)
        flags.loc[top_indices] = True
        df[track_col] = flags

    df["is_selected"] = df["tier1_or_tier2"] | df[sig_columns].any(axis=1)

    n_selected = int(df["is_selected"].sum())
    logger.info(
        f"select_candidates: selected {n_selected:,} / {len(df):,} clones "
        f"(top_n={top_n} per signature; signatures={list(signatures)})"
    )
    return df
