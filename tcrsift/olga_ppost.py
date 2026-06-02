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

"""OLGA/SONIA-backed Pgen + Ppost per-clonotype selection features (#143).

The reproducible **precursor-frequency / publicness** axis. To select
high-specificity rare-precursor clones and de-prioritize high-precursor,
broadly cross-reactive ones, the discriminating quantity is set by V(D)J
recombination + thymic selection — the TCR *sequence*, not the
transcriptome. On the B1-2/B1-3 MART-1 pilot, ``log10 Ppost`` reproducibly
separates the public (TRAV12-2) repertoire from private clones across both
donors (AUROC 0.72 / 0.72), while RNA signatures flip direction.

- **Pgen** (OLGA, Sethna 2019): generation probability of the CDR3 under
  the recombination model.
- **Ppost** (SONIA, Sethna 2020): post-selection probability ``Pgen · Q``,
  where ``Q`` is the learned selection factor. ``Ppost > Pgen`` for public
  clones — selection amplifies their publicness.

License boundary
----------------
OLGA and SONIA are **GPL-3.0**. To keep tcrsift Apache-2.0 they are a
**runtime-optional dependency, not vendored**: this module imports them
lazily and only when a Pgen/Ppost computation is actually requested.
Install with ``pip install tcrsift[pgen]``.

Caveats (also surfaced in docs)
-------------------------------
- Pgen/Ppost estimate **naive precursor frequency**, a *proxy* for
  cross-reactivity risk — not a direct measure.
- They say nothing about **functional avidity** (tetramer MFI / koff) —
  that needs a wet-lab assay and is out of scope.

The lightweight numpy-only proxy in :mod:`tcrsift.pgen` remains the default
for ranking when the optional models aren't installed; this module is the
calibrated gold-standard path.
"""

from __future__ import annotations

import contextlib
import io
import logging
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# chain -> (SONIA chain_type, default-model directory name)
_CHAIN_DEFS: dict[str, tuple[str, str]] = {
    "alpha": ("humanTRA", "human_T_alpha"),
    "beta": ("humanTRB", "human_T_beta"),
}

_INSTALL_HINT = (
    "OLGA + SONIA are required for gold-standard Pgen/Ppost but are not "
    "installed. They are a GPL-3.0 optional extra (tcrsift stays Apache-2.0 "
    "by not vendoring them). Install with:\n\n    pip install tcrsift[pgen]\n\n"
    "Or use the numpy-only proxy in tcrsift.pgen.compute_pgen for ranking."
)


def olga_sonia_available() -> bool:
    """True if both ``olga`` and ``sonia`` can be imported."""
    import importlib.util

    return (
        importlib.util.find_spec("olga") is not None
        and importlib.util.find_spec("sonia") is not None
    )


def _require_olga_sonia() -> None:
    if not olga_sonia_available():
        raise ImportError(_INSTALL_HINT)


@dataclass
class ChainModel:
    """Loaded OLGA pgen model + SONIA selection model for one chain.

    Holds everything needed to compute Pgen (OLGA) and Q/Ppost (SONIA)
    without spawning multiprocessing pools — ``Q`` is evaluated through
    SONIA's serial ``evaluate_selection_factors`` and ``Pgen`` through
    OLGA's per-sequence ``compute_aa_CDR3_pgen``.
    """

    chain: str
    evaluator: object  # sonia.evaluate_model.EvaluateModel
    norm_productive: float
    v_alleles: frozenset[str]
    j_alleles: frozenset[str]

    @property
    def pgen_model(self):
        return self.evaluator.pgen_model


# Loaded models are heavy (TensorFlow graph + IGoR tables); cache per chain.
_MODEL_CACHE: dict[str, ChainModel] = {}


def load_chain_model(chain: str) -> ChainModel:
    """Load (and cache) the OLGA+SONIA default human-T model for a chain.

    ``chain`` is ``"alpha"`` or ``"beta"``. Raises :class:`ImportError`
    with an install hint when the optional extra is missing.
    """
    chain = chain.lower()
    if chain not in _CHAIN_DEFS:
        raise ValueError(f"chain must be 'alpha' or 'beta', got {chain!r}")
    if chain in _MODEL_CACHE:
        return _MODEL_CACHE[chain]

    _require_olga_sonia()
    import sonia
    from sonia.evaluate_model import EvaluateModel
    from sonia.sonia_leftpos_rightpos import SoniaLeftposRightpos

    chain_type, model_dir_name = _CHAIN_DEFS[chain]
    model_dir = os.path.join(
        os.path.dirname(sonia.__file__), "default_models", model_dir_name
    )
    # chain_type MUST be passed explicitly: load_dir only loads the SONIA
    # weights/features, not the recombination topology — without it an alpha
    # (VJ) model silently loads as the default beta (VDJ) one.
    # processes=1 keeps EvaluateModel from constructing a multiprocessing
    # Pool; we never call the parallel compute_all_pgens path anyway. Both
    # constructors print load chatter ("Cannot find data_seqs.tsv", etc.) —
    # mute it.
    with _suppress_stdout():
        qm = SoniaLeftposRightpos(chain_type=chain_type, load_dir=model_dir)
        ev = EvaluateModel(sonia_model=qm, processes=1)

    v_alleles = frozenset(g[0] for g in ev.genomic_data.genV)
    j_alleles = frozenset(g[0] for g in ev.genomic_data.genJ)
    model = ChainModel(
        chain=chain,
        evaluator=ev,
        norm_productive=float(qm.norm_productive),
        v_alleles=v_alleles,
        j_alleles=j_alleles,
    )
    _MODEL_CACHE[chain] = model
    return model


def supported_alleles(chain: str, segment: str) -> frozenset[str]:
    """Return the V or J alleles the loaded model recognizes (#150 hook).

    ``segment`` is ``"v"`` or ``"j"``. This is the explicit "supported
    alleles" set that #150 maps unsupported/novel alleles onto before
    Pgen/Ppost, instead of letting OLGA silently fall back to a default
    usage mask.
    """
    model = load_chain_model(chain)
    seg = segment.lower()
    if seg == "v":
        return model.v_alleles
    if seg == "j":
        return model.j_alleles
    raise ValueError(f"segment must be 'v' or 'j', got {segment!r}")


@contextlib.contextmanager
def _suppress_stdout():
    """Silence OLGA's per-sequence 'Unfamiliar gene/allele' chatter on stdout.

    OLGA prints directly to stdout when it doesn't recognize a gene; that's
    the silent-default-mask path #150 handles explicitly. We mute the raw
    chatter here so library use stays quiet — substitution logging is done
    by the caller against :func:`supported_alleles`, not by parsing stdout.
    """
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        yield


def normalize_gene_name(gene: str | None) -> str | None:
    """Map a CellRanger gene name to OLGA/SONIA's allele-suffixed form.

    OLGA/SONIA anchors are allele-suffixed (``TRAV12-2*01``); CellRanger
    gene calls often lack the allele. Append ``*01`` when missing (#143,
    validated). Returns ``None`` for empty / non-string input.
    """
    if not isinstance(gene, str):
        return None
    g = gene.strip()
    if not g:
        return None
    if "*" not in g:
        g = f"{g}*01"
    return g


def _log10_or_nan(p: float) -> float:
    """log10(p) with 0 / negative / non-finite mapped to NaN."""
    if p is None or not np.isfinite(p) or p <= 0.0:
        return float("nan")
    return float(np.log10(p))


def compute_pgen_ppost(
    df: pd.DataFrame,
    *,
    chain: str = "beta",
    cdr3_col: str = "CDR3_beta",
    v_gene_col: str = "beta_v_gene",
    j_gene_col: str = "beta_j_gene",
    with_ppost: bool = True,
    pgen_col: str = "log10_pgen_olga",
    ppost_col: str = "log10_ppost",
    q_col: str = "sonia_q",
) -> pd.DataFrame:
    """Add OLGA ``log10 Pgen`` and SONIA ``log10 Ppost`` columns to ``df``.

    For each row the CDR3 amino-acid sequence + V/J gene calls are scored
    against the default human-T model for ``chain``. V/J names are mapped to
    OLGA's allele-suffixed form via :func:`normalize_gene_name` (append
    ``*01`` when missing). Pgen is normalized by SONIA's productive
    normalization so it matches SONIA's Ppost = Pgen·Q convention.

    Rows whose CDR3 is empty/non-string, or whose Pgen the model evaluates
    to 0 (out-of-frame, or a gene not recognized even after ``*01``), get
    NaN — never a silently-wrong default. The count of such rows is logged.

    Returns a copy of ``df`` with ``pgen_col`` and (when ``with_ppost``)
    ``q_col`` + ``ppost_col`` added. Does not mutate the input.
    """
    chain = chain.lower()
    if chain not in _CHAIN_DEFS:
        raise ValueError(f"chain must be 'alpha' or 'beta', got {chain!r}")
    if cdr3_col not in df.columns:
        raise ValueError(f"compute_pgen_ppost: missing {cdr3_col!r} column")

    model = load_chain_model(chain)
    pgen_model = model.pgen_model

    out = df.copy()
    n = len(out)
    cdr3 = out[cdr3_col]
    v_series = out[v_gene_col] if v_gene_col in out.columns else pd.Series([None] * n, index=out.index)
    j_series = out[j_gene_col] if j_gene_col in out.columns else pd.Series([None] * n, index=out.index)

    log10_pgen = np.full(n, np.nan)
    pgen_normed = np.full(n, np.nan)
    # SONIA seqs for the Q pass: [cdr3, V, J]; only valid rows are evaluated.
    valid_idx: list[int] = []
    sonia_seqs: list[list[str]] = []

    with _suppress_stdout():
        for i, (seq, v, j) in enumerate(zip(cdr3, v_series, j_series)):
            if not isinstance(seq, str) or not seq:
                continue
            v_norm = normalize_gene_name(v)
            j_norm = normalize_gene_name(j)
            if v_norm is None or j_norm is None:
                continue
            try:
                p = pgen_model.compute_aa_CDR3_pgen(seq, v_norm, j_norm)
            except Exception as exc:  # OLGA raises on malformed CDR3 chars
                logger.debug("Pgen failed for %r (%s/%s): %s", seq, v_norm, j_norm, exc)
                continue
            p_norm = p / model.norm_productive if model.norm_productive else p
            lp = _log10_or_nan(p_norm)
            log10_pgen[i] = lp
            if np.isfinite(lp):
                pgen_normed[i] = p_norm
                valid_idx.append(i)
                sonia_seqs.append([seq, v_norm, j_norm])

    out[pgen_col] = log10_pgen
    n_nan = int(np.isnan(log10_pgen).sum())
    if n_nan:
        logger.info(
            "compute_pgen_ppost(%s): %d/%d rows had no computable Pgen "
            "(empty CDR3, missing V/J, or Pgen=0) → NaN", chain, n_nan, n,
        )

    if with_ppost:
        q_arr = np.full(n, np.nan)
        ppost = np.full(n, np.nan)
        if sonia_seqs:
            with _suppress_stdout():
                q_vals = np.asarray(
                    model.evaluator.evaluate_selection_factors(sonia_seqs)
                ).ravel()
            for k, i in enumerate(valid_idx):
                q_arr[i] = q_vals[k]
                ppost[i] = _log10_or_nan(pgen_normed[i] * q_vals[k])
        out[q_col] = q_arr
        out[ppost_col] = ppost

    return out


def flag_private_candidates(
    df: pd.DataFrame,
    *,
    score_col: str = "log10_ppost",
    freq_col: str | None = None,
    score_quantile: float = 0.25,
    min_freq: float = 0.01,
    out_col: str = "private_candidate",
) -> pd.Series:
    """Flag low-Pgen/Ppost **and** expanded clones — the private candidates.

    A private candidate is rare-precursor (``score_col`` in the bottom
    ``score_quantile`` of the non-NaN distribution) and yet expanded
    (``freq_col`` ≥ ``min_freq``) — exactly the high-specificity clones the
    publicness axis is meant to surface. When ``freq_col`` is None the
    expansion gate is skipped (flag = rare-precursor alone).

    Returns a boolean Series aligned to ``df`` (NaN scores → False).
    """
    if score_col not in df.columns:
        raise ValueError(f"flag_private_candidates: missing {score_col!r}")
    scores = df[score_col].astype(float)
    non_nan = scores.dropna()
    if non_nan.empty:
        return pd.Series(False, index=df.index, name=out_col)
    cutoff = float(non_nan.quantile(score_quantile))
    rare = scores <= cutoff
    if freq_col is not None and freq_col in df.columns:
        expanded = pd.to_numeric(df[freq_col], errors="coerce").fillna(0.0) >= min_freq
        flag = rare & expanded
    else:
        flag = rare
    flag = flag.fillna(False)
    flag.name = out_col
    return flag
