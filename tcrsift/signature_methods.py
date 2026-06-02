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
    min_genes_present: int = 2,
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
    reference explicitly. Missing genes are dropped; fewer than
    ``min_genes_present`` present → all-zero score.
    """
    cols = set(expr.columns)
    if background is not None:
        cols &= set(background.columns)
    up = [g for g in signature.genes_up if g in cols]
    down = [g for g in signature.genes_down if g in cols]
    missing = [g for g in signature.all_genes if g not in cols]
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
    return (
        _combine_block(expr, up, combine, background, log1p)
        - _combine_block(expr, down, combine, background, log1p)
    )


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
    log1p: bool = True,
    background: pd.DataFrame | None = None,
    positive_method: str = "gap",
    quantile: float = 0.75,
    on_context_mismatch: str = "warn",
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Score → call-positive → emit synthetic-method rows for each signature.

    ``signatures`` is an iterable of :class:`Signature` (or names in
    :data:`SIGNATURES`). ``context`` is the sample source's context (see
    :func:`infer_context`).

    A signature is always *computable* — the score is just an off-context
    program you may genuinely want to inspect (e.g. "do any culture cells
    carry an exhaustion-like program?"). So a context mismatch does not
    block by default; it only flags interpretation. ``on_context_mismatch``:

    - ``"warn"`` (default) — log a warning and score anyway.
    - ``"ignore"`` — score silently.
    - ``"raise"`` — block (opt-in strict, e.g. for a locked pipeline).
    """
    resolved = [
        SIGNATURES[s] if isinstance(s, str) else s for s in signatures
    ]
    positive_by_signature: dict[str, pd.Series] = {}
    for sig in resolved:
        if not sig.valid_in(context):
            msg = (
                f"signature {sig.name!r} is only interpretable in "
                f"{sorted(sig.contexts)}, not context={context!r} — "
                "scoring anyway; read the result with that in mind"
            )
            if on_context_mismatch == "raise":
                raise SignatureContextError(msg)
            if on_context_mismatch == "warn":
                logger.warning("%s", msg)
            # "warn" and "ignore" both fall through and score
        scores = score_signature(
            expr, sig, combine=combine, log1p=log1p, background=background,
        )
        positive_by_signature[sig.name] = call_positive(
            scores, method=positive_method, quantile=quantile,
        )
    return signature_methods_long(
        obs, positive_by_signature, clone_col=clone_col, sample_col=sample_col,
    )
