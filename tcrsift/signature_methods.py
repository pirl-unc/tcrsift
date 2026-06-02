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

"""Context-aware gene-expression signatures as synthetic selection methods.

The selection language (:mod:`tcrsift.selection`) is method-agnostic — a
"method" is just a value in the per-``(clone, sample)`` table carrying a
tier + frequency. A method need not be a physical sort: this module scores
cells on a gene set, calls signature-positive cells, and emits them as a
**synthetic method** so the rules can use a signature like a sort.

**A signature's meaning is context-dependent.** CD39/CD103/TOX reads
"chronically engaged antigen in situ" *only* in fresh tissue TILs — in a
PBMC peptide-stim culture, CD103 needs TGF-β/residency it never sees and
TOX/PD-1 needs chronic stim a short culture doesn't provide. So each
:class:`Signature` declares the ``contexts`` it is interpretable in, and
:func:`build_signature_methods` refuses to score it against mismatched
data (e.g. a ``tissue``-only signature on ``source=culture`` samples).
Signatures with no declared context are *invariant* (valid anywhere).

Signatures are **signed** (``genes_up`` − ``genes_down``) so loss-of-naive
axes work, and scored by either z-scored mean (default; weights genes
equally) or raw-TPM mean (``combine="mean"``; absolute expression, so
high-TPM genes dominate). Small "focal" and larger "broad" panels are
offered per axis so they can be compared.

NOTE: gene memberships and context assignments below are sensible
literature-backed defaults intended for review/tuning, not gospel.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Sample-source contexts a signature can be interpreted in.
CULTURE = "culture"        # in-vitro peptide/antigen stim of PBMC
CIRCULATING = "circulating"  # resting blood / tetramer-sorted PBMC
TISSUE = "tissue"          # fresh tumor / tissue TILs


@dataclass(frozen=True)
class Signature:
    """A named, context-tagged, optionally-signed gene signature.

    ``contexts`` empty = invariant (interpretable in any sample source).
    ``genes_down`` subtracts (loss-of-naive etc.). ``panel`` is just a
    "focal"/"broad" label for the small-vs-large variants.
    """

    name: str
    genes_up: tuple[str, ...]
    genes_down: tuple[str, ...] = ()
    contexts: frozenset[str] = field(default_factory=frozenset)
    panel: str = "focal"
    description: str = ""

    @property
    def all_genes(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(self.genes_up + self.genes_down))

    def valid_in(self, context: str | None) -> bool:
        """True if interpretable in ``context`` (invariant ⇒ always)."""
        return not self.contexts or context is None or context in self.contexts


# --- Context-invariant axes (same meaning in blood, culture, tissue) ------

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

# --- Context-specific axes (only interpretable in the right source) -------

_CONTEXTUAL_SIGNATURES = [
    Signature(
        "ActivatedInCulture", ("NR4A1", "TNFRSF9"),
        contexts=frozenset({CULTURE}), panel="focal",
        description="Recent cognate TCR engagement in vitro (Nur77, 4-1BB). "
        "Correlated with an AIM sort by construction.",
    ),
    Signature(
        "ActivatedInCultureBroad",
        ("NR4A1", "NR4A2", "NR4A3", "EGR2", "EGR3", "FOS", "JUN",
         "CD69", "TNFRSF9", "TNFRSF4", "IL2RA", "CD40LG", "IFNG", "TNF"),
        contexts=frozenset({CULTURE}), panel="broad",
        description="Broad acute-activation program (in-vitro stim).",
    ),
    Signature(
        "NeoantigenExperiencedTIL", ("ENTPD1", "CXCL13"),
        contexts=frozenset({TISSUE}), panel="focal",
        description="Chronic in-situ antigen engagement (CD39/CXCL13; "
        "Simoni/Duhen 2018). FRESH TUMOR TILs ONLY.",
    ),
    Signature(
        "NeoantigenExperiencedTILBroad",
        ("ENTPD1", "ITGAE", "TOX", "PDCD1", "LAG3", "HAVCR2", "CXCL13",
         "TIGIT", "CTLA4"),
        contexts=frozenset({TISSUE}), panel="broad",
        description="Tumor-reactive/exhausted TIL program. FRESH TILs ONLY.",
    ),
    Signature(
        "CirculatingMemory", ("GZMK", "IL7R"), ("CCR7",),
        contexts=frozenset({CIRCULATING}), panel="focal",
        description="Resting memory differentiation in blood "
        "(effector-memory up, naive down).",
    ),
]

SIGNATURES: dict[str, Signature] = {
    s.name: s for s in (_INVARIANT_SIGNATURES + _CONTEXTUAL_SIGNATURES)
}

# Map a sample-sheet `source` value to a signature context.
_SOURCE_TO_CONTEXT: dict[str, str] = {
    "culture": CULTURE,
    "tetramer": CIRCULATING,
    "sct": CIRCULATING,
    "pbmc": CIRCULATING,
    "til": TISSUE,
    "tumor": TISSUE,
    "tissue": TISSUE,
}


def infer_context(source: str | None) -> str | None:
    """Map a sample-sheet ``source`` to a signature context (or None)."""
    if source is None:
        return None
    return _SOURCE_TO_CONTEXT.get(str(source).strip().lower())


def _combine_block(expr: pd.DataFrame, genes: list[str], combine: str) -> pd.Series:
    if not genes:
        return pd.Series(0.0, index=expr.index)
    sub = expr[genes].astype(float)
    if combine == "zscore":
        mu = sub.mean(axis=0)
        sd = sub.std(axis=0, ddof=0).replace(0.0, 1.0)
        sub = (sub - mu) / sd
    elif combine == "mean":
        pass  # raw TPM/normalized mean — absolute expression
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
    min_genes_present: int = 2,
) -> pd.Series:
    """Per-cell signed signature score from a cells × genes frame.

    ``score = combine(genes_up) − combine(genes_down)``. ``combine`` is
    ``"zscore"`` (per-gene standardize, equal weight) or ``"mean"`` (raw
    TPM mean, so high-expression genes dominate). Missing genes are
    dropped with a warning; fewer than ``min_genes_present`` present →
    all-zero score.
    """
    up = [g for g in signature.genes_up if g in expr.columns]
    down = [g for g in signature.genes_down if g in expr.columns]
    missing = [g for g in signature.all_genes if g not in expr.columns]
    if missing:
        logger.warning(
            "score_signature[%s]: %d gene(s) absent from panel: %s",
            signature.name, len(missing), missing,
        )
    if len(up) + len(down) < min_genes_present:
        logger.warning(
            "score_signature[%s]: <%d genes present; returning zeros.",
            signature.name, min_genes_present,
        )
        return pd.Series(0.0, index=expr.index)
    return _combine_block(expr, up, combine) - _combine_block(expr, down, combine)


def _largest_gap_cutoff(values: np.ndarray, *, search_top: float = 0.5) -> float:
    v = np.sort(values)[::-1]
    n = len(v)
    if n < 2:
        return float(v[0]) if n else float("inf")
    k = max(2, int(np.ceil(n * search_top)))
    top = v[:k]
    i = int(np.argmax(top[:-1] - top[1:]))
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
) -> pd.Series:
    """Boolean positive call over any score Series (cells or clones).

    ``method``: ``"quantile"`` (fixed top fraction), ``"gap"`` (adaptive —
    largest gap in the top of the sorted scores → the separated high
    cluster), or ``"otsu"`` (adaptive bimodal split). Explicit
    ``threshold`` overrides ``method``.
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


class SignatureContextError(ValueError):
    """Raised when a signature is applied to a context it can't be read in."""


def build_signature_methods(
    expr: pd.DataFrame,
    obs: pd.DataFrame,
    *,
    signatures,
    context: str | None = None,
    combine: str = "zscore",
    positive_method: str = "gap",
    quantile: float = 0.75,
    on_context_mismatch: str = "raise",
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Score → call-positive → emit synthetic-method rows for each signature.

    ``signatures`` is an iterable of :class:`Signature` (or names in
    :data:`SIGNATURES`). ``context`` is the sample source's context (see
    :func:`infer_context`); a signature not valid in ``context`` is, per
    ``on_context_mismatch``, ``"raise"`` (default), ``"warn"`` (skip with a
    warning), or ``"ignore"`` (score anyway). This is the guard that stops
    a ``tissue``-only signature being scored on ``culture`` data.
    """
    resolved = [
        SIGNATURES[s] if isinstance(s, str) else s for s in signatures
    ]
    positive_by_signature: dict[str, pd.Series] = {}
    for sig in resolved:
        if not sig.valid_in(context):
            msg = (
                f"signature {sig.name!r} is only interpretable in "
                f"{sorted(sig.contexts)}, not context={context!r}"
            )
            if on_context_mismatch == "raise":
                raise SignatureContextError(msg)
            if on_context_mismatch == "warn":
                logger.warning("Skipping %s", msg)
                continue
            # "ignore": fall through and score anyway
        scores = score_signature(expr, sig, combine=combine)
        positive_by_signature[sig.name] = call_positive(
            scores, method=positive_method, quantile=quantile,
        )
    return signature_methods_long(
        obs, positive_by_signature, clone_col=clone_col, sample_col=sample_col,
    )
