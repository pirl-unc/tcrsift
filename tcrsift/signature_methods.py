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
tier + frequency. A method need not be a physical sort: this module scores
cells on a gene set, calls signature-positive cells, and emits them as a
**synthetic method** so the rules can use a signature like a sort.

Signatures are **signed** (``genes_up`` − ``genes_down``, so loss-of-naive
axes work). Canonical scoring (defaults): per gene, the z-score across
cells of ``log(1 + TPM)``, then the mean of those z-scores across the
signature (``combine="mean"`` instead averages the values directly).
Small "focal" and larger "broad" panels are offered per axis.

NOTE: a signature's *meaning* still depends on the sample source — e.g.
``TumorReactive`` (CD39/CD103/TOX) only reads "chronic in-situ antigen
exposure" in fresh tissue TILs, not in PBMC culture. That caveat lives in
each signature's ``description``; it is documentation, not an enforced
guard, so you can still score any signature on any data (and should, when
exploring). Gene memberships are literature-backed defaults for tuning.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Signature:
    """A named, optionally-signed gene signature.

    ``genes_down`` subtracts (loss-of-naive etc.). ``panel`` is a
    "focal"/"broad" label for the small-vs-large variants. ``description``
    carries the biological caveat (e.g. which sample source it's meaningful
    in) as prose.
    """

    name: str
    genes_up: tuple[str, ...]
    genes_down: tuple[str, ...] = ()
    panel: str = "focal"
    description: str = ""

    @property
    def all_genes(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(self.genes_up + self.genes_down))


_INVARIANT_SIGNATURES = [
    Signature(
        "Proliferation", ("MKI67",), panel="focal",
        description="Recently divided (Ki-67).",
    ),
    Signature(
        "ProliferationBroad",
        ("MKI67", "TOP2A", "PCNA", "TYMS", "CCNB1", "CDK1", "BIRC5"),
        panel="broad", description="Cell-cycle / proliferation program.",
    ),
    Signature(
        "Cytolytic", ("PRF1", "GZMB"), panel="focal",
        description="Canonical cytotoxic effector pair (Caushi/Krishna/Hanada).",
    ),
    Signature(
        "CytolyticBroad",
        ("PRF1", "GZMB", "GNLY", "NKG7", "GZMH", "GZMA", "KLRD1"),
        panel="broad", description="Broad cytotoxic-effector panel.",
    ),
    Signature(
        "Differentiated", (), ("CCR7", "TCF7"), panel="focal",
        description="Loss of naive (down CCR7/TCF7) — antigen-experienced.",
    ),
    Signature(
        "DifferentiatedBroad", (), ("CCR7", "SELL", "TCF7", "LEF1"),
        panel="broad", description="Loss of the full naive program.",
    ),
]

# The two headline composites the pilot uses. Their descriptions note the
# sample source they're meaningful in — read with that in mind.

_COMPOSITE_SIGNATURES = [
    Signature(
        "TumorReactive",
        ("CXCL13", "ENTPD1", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX",
         "CTLA4", "ITGAE"),
        panel="broad",
        description="Chronic in-situ antigen engagement / tumor-reactive "
        "exhaustion (CD39/CXCL13 + exhaustion panel + CD103; Simoni/Duhen "
        "2018). Meaningful in fresh tumor TILs; reads as noise in PBMC culture.",
    ),
    Signature(
        "AntigenExperienced",
        ("TNFRSF9", "MKI67", "IFNG", "GZMB", "PRF1", "GNLY", "NKG7"),
        panel="broad",
        description="Recent cognate engagement + effector program "
        "(4-1BB/Ki-67 + activation). In-culture activation axis; correlated "
        "with an AIM sort by construction.",
    ),
    Signature(
        "AIM", ("TNFRSF9", "TNFRSF4"), panel="focal",
        description="Transcriptional read-out of the AIM (activation-induced "
        "marker) assay: 4-1BB (CD137) + OX40 (CD134). The GEX analogue of an "
        "AIMpos sort — use it to check sort/transcriptome concordance.",
    ),
    Signature(
        "AIMBroad",
        ("TNFRSF9", "TNFRSF4", "CD69", "IL2RA", "CD40LG"),
        panel="broad",
        description="Broad AIM panel: 4-1BB, OX40, CD69, CD25 (IL2RA), "
        "CD40L. Surface activation markers used across AIM-assay variants.",
    ),
    Signature(
        "CirculatingMemory", ("GZMK", "IL7R"), ("CCR7",), panel="focal",
        description="Resting memory differentiation in blood "
        "(effector-memory up, naive down).",
    ),
]

SIGNATURES: dict[str, Signature] = {
    s.name: s for s in (_INVARIANT_SIGNATURES + _COMPOSITE_SIGNATURES)
}

# Convenience grouping of the two headline selection composites (the names
# the pilot uses).
SELECTION_SIGNATURES: dict[str, Signature] = {
    "TumorReactive": SIGNATURES["TumorReactive"],
    "AntigenExperienced": SIGNATURES["AntigenExperienced"],
}


def expression_frame_from_adata(
    adata, genes, *, layer: str | None = None, on_missing: str = "error",
) -> pd.DataFrame:
    """Extract a cells × genes (symbol-keyed) expression frame from AnnData.

    Resolves HGNC symbols against ``adata.var_names`` or, when those are
    Ensembl IDs (the ``load_samples`` default), a gene-symbol column in
    ``adata.var`` (``gene_symbols`` / ``feature_name`` / …). Densifies the
    requested columns only; duplicate symbols collapse to the highest-mean
    copy.

    ``on_missing`` governs requested genes absent from the matrix:
    ``"error"`` (default) raises — scoring a signature on a panel that
    can't measure all its genes silently changes the signature's meaning,
    so we fail loud; ``"warn"`` logs and omits; ``"ignore"`` omits silently.
    """
    from scipy.sparse import issparse

    var = adata.var
    symbol_series = None
    for col in ("gene_symbols", "feature_name", "symbol", "gene_symbol", "gene_name"):
        if col in getattr(var, "columns", []):
            symbol_series = var[col].astype(str)
            break
    if symbol_series is None:
        symbol_series = pd.Series(adata.var_names.astype(str))
    sym_upper = symbol_series.str.upper().to_numpy()

    X = adata.layers[layer] if layer is not None else adata.X
    cols: dict[str, np.ndarray] = {}
    missing: list[str] = []
    for g in genes:
        idx = np.where(sym_upper == str(g).upper())[0]
        if idx.size == 0:
            missing.append(g)
            continue
        if idx.size > 1:
            means = np.asarray(X[:, idx].mean(axis=0)).ravel()
            idx = idx[[int(np.argmax(means))]]
        col = X[:, idx[0]]
        cols[g] = (
            np.asarray(col.todense()).ravel() if issparse(X)
            else np.asarray(col).ravel()
        )
    if missing:
        msg = (
            f"expression_frame_from_adata: {len(missing)} gene(s) not found "
            f"in the expression matrix: {missing}"
        )
        if on_missing == "error":
            raise KeyError(msg)
        if on_missing == "warn":
            logger.warning("%s", msg)
    return pd.DataFrame(cols, index=adata.obs_names.astype(str))


def _combine_block(
    expr: pd.DataFrame,
    genes: list[str],
    combine: str,
    background: pd.DataFrame | None = None,
    log1p: bool = True,
) -> pd.Series:
    if not genes:
        return pd.Series(0.0, index=expr.index)
    sub = expr[genes].astype(float)
    if log1p:
        sub = np.log1p(sub)  # log(1 + TPM)
    if combine == "zscore":
        # Standardize each gene across cells against the BACKGROUND
        # population's mean/sd (default: the input cells themselves), then
        # average the per-gene z-scores — the "mean of z-scores". Pass a
        # background frame to z-score against a defined reference (full
        # dataset, a control population, or per-sample).
        ref = (background if background is not None else expr)[genes].astype(float)
        if log1p:
            ref = np.log1p(ref)
        mu = ref.mean(axis=0)
        sd = ref.std(axis=0, ddof=0).replace(0.0, 1.0)
        sub = (sub - mu) / sd
    elif combine == "mean":
        pass  # mean of log1p(TPM), or raw TPM when log1p=False
    else:
        raise ValueError(
            f"score_signature: unknown combine={combine!r} (expected 'zscore' or 'mean')"
        )
    return sub.mean(axis=1)


def score_signature(
    expr: pd.DataFrame,
    signature: Signature,
    *,
    combine: str = "zscore",
    log1p: bool = True,
    background: pd.DataFrame | None = None,
    on_missing: str = "error",
) -> pd.Series:
    """Per-cell signed signature score from a cells × genes (TPM) frame.

    Canonical scoring (defaults): for each gene, take the **z-score across
    cells of ``log(1 + TPM)``**, then average those z-scores across the
    signature's genes — ``score = mean_z(genes_up) − mean_z(genes_down)``.

    ``combine``:
    - ``"zscore"`` (default) — mean of per-gene z-scores; genes contribute
      equally regardless of baseline.
    - ``"mean"`` — mean of the per-gene values (``log1p`` still applies
      unless disabled); absolute scale, high-expression genes dominate.

    ``log1p`` (default True) applies ``log(1+x)`` before scoring — pass
    TPM/normalized counts. Set False if the input is already log-space.
    ``background`` is the population whose per-gene mean/sd define the
    z-score (``"zscore"`` only); defaults to the input cells — pass the
    full dataset / a control population / a per-sample slice to set the
    reference explicitly.

    ``on_missing`` governs signature genes absent from ``expr`` (and
    ``background``): ``"error"`` (default) **raises** — a partial panel is a
    different signature, so we never silently degrade it; ``"warn"`` scores
    on the present genes with a warning; ``"ignore"`` scores silently.
    """
    cols = set(expr.columns)
    if background is not None:
        cols &= set(background.columns)
    up = [g for g in signature.genes_up if g in cols]
    down = [g for g in signature.genes_down if g in cols]
    missing = [g for g in signature.all_genes if g not in cols]
    if missing:
        msg = (
            f"score_signature[{signature.name}]: required gene(s) have no "
            f"expression in the input: {missing}"
        )
        if on_missing == "error":
            raise KeyError(msg)
        if on_missing == "warn":
            logger.warning("%s", msg)
        # "ignore": score on what's present
    if not up and not down:
        return pd.Series(0.0, index=expr.index)
    return (
        _combine_block(expr, up, combine, background, log1p)
        - _combine_block(expr, down, combine, background, log1p)
    )


def _largest_gap_cutoff(
    values: np.ndarray, *, search_top: float = 0.5, min_gap_ratio: float = 3.0,
) -> float:
    # The gap is sought only within the top ``search_top`` fraction so the
    # cut lands high (a minority "above the rest"), not at a low outlier.
    # Consequence: a cluster larger than ``search_top`` of cells has its
    # separating gap below the window and won't be found — raise
    # ``search_top`` for a signature expected to be broadly positive.
    v = np.sort(values)[::-1]  # descending
    n = len(v)
    if n < 2:
        return float(v[0]) if n else float("inf")
    k = max(2, int(np.ceil(n * search_top)))
    top = v[:k]
    gaps = top[:-1] - top[1:]
    i = int(np.argmax(gaps))
    largest = float(gaps[i])
    # "Clusters above the others" only if the top gap stands out vs the
    # typical spacing; on a smooth distribution there is no distinct high
    # cluster, so return a cutoff above the max → empty positive set
    # (rather than slicing off an arbitrary top few).
    # Typical spacing = median of ALL consecutive drops (the bulk gap). Using
    # only positive drops would be skewed by ties (e.g. many zero-expression
    # cells) into a too-strict floor; the bulk median is ~0 then, so a real
    # cluster above the zeros is still picked.
    drops = -np.diff(v)  # consecutive drops over the full sorted array (≥0)
    typical = float(np.median(drops)) if drops.size else 0.0
    if typical > 0 and largest < min_gap_ratio * typical:
        return float("inf")
    return float((top[i] + top[i + 1]) / 2.0)


def _otsu_cutoff(values: np.ndarray) -> float:
    v = np.sort(values)
    n = len(v)
    if n < 2:
        return float(v[0]) if n else float("inf")
    csum = np.cumsum(v)
    total = csum[-1]
    best_t, best_var = float((v[0] + v[-1]) / 2.0), -1.0
    for i in range(1, n):
        w0 = i / n
        m0 = csum[i - 1] / i
        m1 = (total - csum[i - 1]) / (n - i)
        var = w0 * (1.0 - w0) * (m0 - m1) ** 2
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
    min_gap_ratio: float = 3.0,
) -> pd.Series:
    """Boolean positive call over any score Series (cells or clones).

    ``method``: ``"quantile"`` (fixed top fraction), ``"gap"`` (adaptive —
    largest gap in the top of the sorted scores → the separated high
    cluster; returns **none** when no gap stands out by ``min_gap_ratio``×
    the typical spacing, i.e. no distinct cluster), or ``"otsu"`` (adaptive
    bimodal split). Explicit ``threshold`` overrides ``method``.

    ``search_top`` (``"gap"`` only) bounds where the gap is sought to the
    top fraction of cells, so the positive set is a minority "above the
    rest". A positive cluster larger than ``search_top`` won't be found —
    raise it (e.g. 0.8) for a signature expected to be broadly positive.
    """
    arr = scores.to_numpy(dtype=float)
    if threshold is not None:
        cut = float(threshold)
    elif method == "quantile":
        if quantile is None:
            raise ValueError("call_positive(method='quantile') needs quantile")
        cut = float(np.quantile(arr, quantile))
    elif method == "gap":
        cut = _largest_gap_cutoff(arr, search_top=search_top, min_gap_ratio=min_gap_ratio)
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

    Per signature, positive cells aggregate per ``(clone, sample)`` into a
    row shaped like :func:`tcrsift.clonotype.build_clone_sample_long`
    output — ``(CDR3ab, sample, method=<signature>, cells, frequency)`` —
    ready to concat with the real table.
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
    signatures: Iterable[Signature | str] | None = None,
    combine: str = "zscore",
    log1p: bool = True,
    background: pd.DataFrame | None = None,
    positive_method: str = "gap",
    quantile: float = 0.75,
    on_missing: str = "error",
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Score → call-positive → emit synthetic-method rows for each signature.

    ``signatures`` is an iterable of :class:`Signature` (or names in
    :data:`SIGNATURES`), defaulting to :data:`SELECTION_SIGNATURES`. Each is
    scored (z-score of log1p TPM, mean across genes), thresholded
    (``positive_method``: ``gap``/``otsu``/``quantile``), and emitted as a
    synthetic method ready to concat with the real clone-sample-long table.
    ``on_missing`` (default ``"error"``) forbids scoring a signature whose
    genes aren't all in ``expr`` — see :func:`score_signature`.
    """
    if signatures is None:
        signatures = SELECTION_SIGNATURES.values()
    resolved = [
        SIGNATURES[s] if isinstance(s, str) else s for s in signatures
    ]
    positive_by_signature: dict[str, pd.Series] = {}
    for sig in resolved:
        scores = score_signature(
            expr, sig, combine=combine, log1p=log1p, background=background,
            on_missing=on_missing,
        )
        positive_by_signature[sig.name] = call_positive(
            scores, method=positive_method, quantile=quantile,
        )
    return signature_methods_long(
        obs, positive_by_signature, clone_col=clone_col, sample_col=sample_col,
    )
