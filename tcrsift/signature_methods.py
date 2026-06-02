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

"""Gene-expression signatures as synthetic selection methods.

The selection language (:mod:`tcrsift.selection`) is method-agnostic — a
"method" is just a value in the per-``(clone, sample)`` table carrying a
tier + frequency. Physical sorts (``AIMpos``, ``tetpos``…) come from
samples, but a method need not be a physical sort. This module scores
cells on a gene set, calls signature-positive cells, and emits them as a
**synthetic method** so the rules can use a signature like a sort —
include it in ``shared``, spawn ``private_<signature>``, or pair it with a
sort via ``method_pair`` ("tetpos AND signature-positive").

Two selection signatures targeting different axes, composed from the
canonical gene sets in :mod:`tcrsift.signatures` so they can be tried
head-to-head:

- ``TumorReactive`` — chronic-exposure / tumor-reactive / exhaustion axis
  (CXCL13, CD39, CD103, plus the exhaustion panel). Most *complementary*
  to an activation sort like AIMpos.
- ``AntigenExperienced`` — antigen-experienced / activated-in-culture
  (recent-antigen 4-1BB+Ki-67 plus the effector-activation panel). Shares
  more with activation sorts; offered so the two can be compared.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .signatures import (
    ACTIVATION_GENES_HGNC,
    ANTIGEN_RESPONSE_GENES_HGNC,
    EXHAUSTION_GENES_HGNC,
    TUMOR_REACTIVE_GENES_HGNC,
)

logger = logging.getLogger(__name__)


def _dedup(*groups) -> tuple[str, ...]:
    seen: dict[str, None] = {}
    for group in groups:
        for g in group:
            seen.setdefault(g, None)
    return tuple(seen)


# Chronic-exposure / tumor-reactive axis: CXCL13 + CD39 (tumor_reactive),
# the exhaustion panel, and CD103 (ITGAE, tissue-resident).
TUMOR_REACTIVE_SIGNATURE: tuple[str, ...] = _dedup(
    TUMOR_REACTIVE_GENES_HGNC, EXHAUSTION_GENES_HGNC, ("ITGAE",),
)

# Antigen-experienced / activated-in-culture axis: recent-antigen
# (4-1BB + Ki-67) plus the broad effector-activation panel.
ANTIGEN_EXPERIENCED_SIGNATURE: tuple[str, ...] = _dedup(
    ANTIGEN_RESPONSE_GENES_HGNC, ACTIVATION_GENES_HGNC,
)

SELECTION_SIGNATURES: dict[str, tuple[str, ...]] = {
    "TumorReactive": TUMOR_REACTIVE_SIGNATURE,
    "AntigenExperienced": ANTIGEN_EXPERIENCED_SIGNATURE,
}


def score_signature(
    expr: pd.DataFrame,
    genes,
    *,
    z_score: bool = True,
    min_genes_present: int = 2,
) -> pd.Series:
    """Per-cell signature score from a cells × genes expression frame.

    Mean expression across the signature genes present in
    ``expr.columns``; with ``z_score`` (default) each gene is standardized
    across cells first so highly-expressed genes don't dominate. Missing
    genes are dropped (with a warning); fewer than ``min_genes_present``
    present → all-zero score (can't score this panel).
    """
    present = [g for g in genes if g in expr.columns]
    missing = [g for g in genes if g not in expr.columns]
    if missing:
        logger.warning(
            "score_signature: %d/%d genes absent from panel: %s",
            len(missing), len(list(genes)), missing,
        )
    if len(present) < min_genes_present:
        logger.warning(
            "score_signature: only %d gene(s) present (<%d); returning zeros.",
            len(present), min_genes_present,
        )
        return pd.Series(0.0, index=expr.index)
    sub = expr[present].astype(float)
    if z_score:
        mu = sub.mean(axis=0)
        sd = sub.std(axis=0, ddof=0).replace(0.0, 1.0)
        sub = (sub - mu) / sd
    return sub.mean(axis=1)


def _largest_gap_cutoff(values: np.ndarray, *, search_top: float = 0.5) -> float:
    """Cutoff at the biggest drop in the top of the sorted distribution.

    Sorts descending and, within the top ``search_top`` fraction, finds
    the largest consecutive drop — the visual gap separating a high
    cluster from the bulk — and returns the midpoint. The count above it
    is whatever clusters there: small when a cluster is clearly
    separated, near-empty when the scores are smooth.
    """
    v = np.sort(values)[::-1]  # descending
    n = len(v)
    if n < 2:
        return float(v[0]) if n else float("inf")
    # Only look for the gap among the top candidates, so we isolate a
    # small high cluster rather than splitting the whole distribution.
    k = max(2, int(np.ceil(n * search_top)))
    top = v[:k]
    diffs = top[:-1] - top[1:]
    i = int(np.argmax(diffs))  # gap between top[i] and top[i+1]
    return float((top[i] + top[i + 1]) / 2.0)


def _otsu_cutoff(values: np.ndarray) -> float:
    """1-D Otsu threshold: the split maximizing between-group variance."""
    v = np.sort(values)
    n = len(v)
    if n < 2:
        return float(v[0]) if n else float("inf")
    best_t, best_var = (v[0] + v[-1]) / 2.0, -1.0
    csum = np.cumsum(v)
    total = csum[-1]
    for i in range(1, n):
        w0 = i / n
        w1 = 1.0 - w0
        m0 = csum[i - 1] / i
        m1 = (total - csum[i - 1]) / (n - i)
        var = w0 * w1 * (m0 - m1) ** 2
        if var > best_var:
            best_var = var
            best_t = float((v[i - 1] + v[i]) / 2.0)
    return best_t


def call_positive(
    scores: pd.Series,
    *,
    method: str = "quantile",
    quantile: float | None = 0.75,
    threshold: float | None = None,
    search_top: float = 0.5,
) -> pd.Series:
    """Boolean positive call over a score Series (cells *or* clones).

    Works on any per-item score, so the same call selects positive cells
    or — given a per-clone signature score — the small set of clones that
    cluster above the rest.

    Methods:

    - ``"quantile"`` (default) — fixed top fraction (``quantile=0.75`` →
      top quartile). Predictable count.
    - ``"gap"`` — **adaptive**: cut at the largest gap in the top
      ``search_top`` of the sorted scores, taking only the separated high
      cluster. Variable, usually-small count; the data picks it.
    - ``"otsu"`` — adaptive bimodal split (between-group variance).

    An explicit ``threshold`` overrides ``method``.
    """
    arr = scores.to_numpy(dtype=float)
    if threshold is not None:
        cut = float(threshold)
    elif method == "quantile":
        if quantile is None:
            raise ValueError("call_positive(method='quantile') needs quantile")
        cut = float(np.quantile(arr, quantile))
    elif method == "gap":
        cut = _largest_gap_cutoff(arr, search_top=search_top)
    elif method == "otsu":
        cut = _otsu_cutoff(arr)
    else:
        raise ValueError(
            f"call_positive: unknown method {method!r} "
            "(expected 'quantile', 'gap', or 'otsu')"
        )
    return scores > cut


def signature_methods_long(
    obs: pd.DataFrame,
    positive_by_signature: dict[str, pd.Series],
    *,
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Emit signature-positive cells as synthetic clone-sample-long rows.

    Per signature, positive cells are aggregated per ``(clone, sample)``
    into a row shaped like :func:`tcrsift.clonotype.build_clone_sample_long`
    output — ``(CDR3ab, sample, method, cells, frequency)`` — where
    ``method`` is the signature name and ``frequency`` is the clone's
    within-sample share of that signature's positive cells. Concatenate
    with the real clone-sample-long table and the rules treat the
    signature as just another method.
    """
    cols = [clone_col, sample_col, "method", "cells", "frequency"]
    frames: list[pd.DataFrame] = []
    for sig, positive in positive_by_signature.items():
        mask = positive.reindex(obs.index, fill_value=False).to_numpy()
        pos = obs.loc[mask, [clone_col, sample_col]].dropna()
        if pos.empty:
            continue
        counts = (
            pos.groupby([clone_col, sample_col], observed=True)
            .size().rename("cells").reset_index()
        )
        totals = counts.groupby(sample_col, observed=True)["cells"].transform("sum")
        counts["frequency"] = counts["cells"] / totals
        counts["method"] = sig
        frames.append(counts[cols])
    if not frames:
        return pd.DataFrame(columns=cols)
    return pd.concat(frames, ignore_index=True)


def build_signature_methods(
    expr: pd.DataFrame,
    obs: pd.DataFrame,
    *,
    signatures: dict[str, tuple[str, ...]] | None = None,
    quantile: float = 0.75,
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Score → call-positive → emit synthetic-method rows for each signature.

    Convenience wrapper: ``expr`` is a cells × genes frame aligned to
    ``obs`` (which carries ``clone_col`` / ``sample_col``). Returns a
    clone-sample-long frame ready to concat with the real one before
    :func:`tcrsift.selection.select_from_clone_sample_long`.
    """
    sigs = signatures if signatures is not None else SELECTION_SIGNATURES
    positive_by_signature = {
        name: call_positive(score_signature(expr, genes), quantile=quantile)
        for name, genes in sigs.items()
    }
    return signature_methods_long(
        obs, positive_by_signature, clone_col=clone_col, sample_col=sample_col,
    )
