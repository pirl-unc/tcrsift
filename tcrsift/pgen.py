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

"""Lightweight Pgen estimator for TCR β / α CDR3 sequences (#58).

The gold-standard generation-probability tool is **OLGA**
(Sethna et al. 2019), which fits V/D/J/insertion models and computes
Pgen by dynamic programming. OLGA gives the right answer but ships
~30 MB of pre-trained model files plus a runtime dependency chain.

Most tcrsift use cases (publicness discounting in DB annotation,
ranking convergent vs. private sequences) just need a coarse Pgen
*proxy* — values that rank sequences from "highly generatable" to
"clearly private", not necessarily numerically calibrated. This
module provides that proxy in <200 lines, with no runtime dependency
beyond numpy.

Decomposition:

    log Pgen(CDR3, V, J, chain)
      ≈ log P(V)             # gene-usage marginal
      + log P(J)             # gene-usage marginal
      + log P(length)        # CDR3 AA length, per chain
      + log P(n_inserted)    # estimated total non-templated nt
      + log P(CDR3 | length) # length-normalized AA-composition term

For an observed sequence, `n_inserted` is estimated from the CDR3
nucleotide length minus a typical templated boundary contribution
(see :data:`_TEMPLATED_NT_TYPICAL`). When the caller doesn't have
nucleotide-level data, we approximate from AA length.

Marginal parameters bundled inline (not as external data files) —
they're small (<3 KB total), and inlining keeps the module
self-contained and easy to audit. Sources are documented per-table.

For real Pgens use OLGA; for ranking/discounting :func:`compute_pgen`
is sufficient.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Iterable

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# Bundled marginal parameters
# =============================================================================
#
# Numbers chosen from order-of-magnitude estimates in the published
# OLGA pre-trained marginals (Sethna 2019), Britanova 2014, Pogorelyy
# 2018, and Robins 2009. They're approximate by design — tcrsift's
# publicness discount only needs the relative ordering to be sensible,
# not numerically calibrated against OLGA.


# Per-junction N-insertion geometric parameter (mean = (1 - p) / p).
# β has two junctions (VβDβ + DβJβ); α has one (VαJα).
_N_INSERT_GEOM_P = {
    "beta": 0.18,   # mean ~4.5 nt per junction → ~9 nt total across two junctions
    "alpha": 0.22,  # mean ~3.5 nt for the single junction
}

# Typical templated nucleotide contribution to CDR3 from V + J + D
# combined. Used to back out an estimate of total inserted nt from
# observed CDR3 length. Values from IMGT mean V/J templated lengths
# after the conserved C / before the conserved F.
_TEMPLATED_NT_TYPICAL = {
    "beta":  21,   # ~5 aa V + ~2 aa D (mean of 5 D-trims) + ~2 aa J = 9 aa = 27 nt; minus some exo trimming → ~21
    "alpha": 24,   # ~6 aa V + ~2 aa J = 8 aa = 24 nt
}

# CDR3 AA length distribution (mean, std) per chain. From Robins 2009,
# Madi 2014. Gaussian approximation — actual distribution is right-skewed.
_CDR3_LENGTH_PARAMS = {
    "beta":  (14.0, 2.5),
    "alpha": (13.5, 2.5),
}

# Marginal V/J usage frequencies for the most-common V/J genes (sums
# to ~0.6 for V, ~1.0 for J). Anything unlisted falls back to the
# uniform-mass complement spread across the remaining ~30 V or ~5 J.
# Values are approximate, derived from Britanova 2014 / Pogorelyy 2018
# orderings — the goal is rank-correct, not numerically exact.
_TRBV_USAGE: dict[str, float] = {
    "TRBV20-1": 0.054, "TRBV28": 0.046, "TRBV5-1": 0.040,
    "TRBV29-1": 0.040, "TRBV9": 0.038, "TRBV6-5": 0.035,
    "TRBV19": 0.032, "TRBV6-1": 0.029, "TRBV4-1": 0.028,
    "TRBV2": 0.027, "TRBV30": 0.026, "TRBV7-9": 0.025,
    "TRBV4-3": 0.024, "TRBV12-3": 0.022, "TRBV6-6": 0.020,
    "TRBV5-4": 0.018, "TRBV3-1": 0.016, "TRBV27": 0.015,
    "TRBV18": 0.014, "TRBV10-3": 0.013,
}

_TRBJ_USAGE: dict[str, float] = {
    "TRBJ2-1": 0.165, "TRBJ2-7": 0.143, "TRBJ2-3": 0.115,
    "TRBJ1-1": 0.103, "TRBJ2-5": 0.099, "TRBJ1-2": 0.083,
    "TRBJ1-5": 0.073, "TRBJ2-2": 0.059, "TRBJ1-6": 0.051,
    "TRBJ1-4": 0.048, "TRBJ1-3": 0.034, "TRBJ2-4": 0.018,
    "TRBJ2-6": 0.009,
}

_TRAV_USAGE: dict[str, float] = {
    # Top ~15 TRAV by published usage. Less peaked than TRBV.
    "TRAV1-2": 0.085, "TRAV12-1": 0.060, "TRAV21": 0.055,
    "TRAV13-1": 0.050, "TRAV8-1": 0.045, "TRAV17": 0.040,
    "TRAV41": 0.040, "TRAV12-2": 0.038, "TRAV26-1": 0.035,
    "TRAV14": 0.032, "TRAV24": 0.030, "TRAV38-2": 0.028,
    "TRAV29": 0.025, "TRAV35": 0.023, "TRAV5": 0.020,
}

_TRAJ_USAGE: dict[str, float] = {
    "TRAJ49": 0.040, "TRAJ48": 0.035, "TRAJ33": 0.033,
    "TRAJ31": 0.030, "TRAJ29": 0.028, "TRAJ12": 0.026,
    "TRAJ27": 0.024, "TRAJ20": 0.022, "TRAJ43": 0.020,
    "TRAJ34": 0.020,
    # ~50 more TRAJ genes; the listed ones cover ~30% of usage.
}


def _strip_allele(gene: str | None) -> str:
    """Drop the IMGT allele suffix (``TRBV20-1*01`` → ``TRBV20-1``)."""
    if not isinstance(gene, str):
        return ""
    return gene.split("*", 1)[0].strip()


def _v_usage_log(gene: str, chain: str) -> float:
    """log P(V) for the given gene; falls back to a uniform tail mass."""
    table = _TRBV_USAGE if chain == "beta" else _TRAV_USAGE
    base = _strip_allele(gene)
    if base in table:
        p = table[base]
    else:
        # Total uncovered tail mass spread over an estimated 40 genes.
        listed_mass = sum(table.values())
        p = max(0.0005, (1.0 - listed_mass) / 40)
    return math.log(p)


def _j_usage_log(gene: str, chain: str) -> float:
    """log P(J) for the given gene; falls back to a uniform tail mass."""
    table = _TRBJ_USAGE if chain == "beta" else _TRAJ_USAGE
    base = _strip_allele(gene)
    if base in table:
        p = table[base]
    else:
        listed_mass = sum(table.values())
        # β has ~13 J; α has ~50. Tail spread accordingly.
        n_tail = 5 if chain == "beta" else 50
        p = max(0.0005, (1.0 - listed_mass) / n_tail)
    return math.log(p)


def _length_log(cdr3_aa_len: int, chain: str) -> float:
    """log P(CDR3 AA length) — Gaussian approximation per chain."""
    mu, sigma = _CDR3_LENGTH_PARAMS[chain]
    z = (cdr3_aa_len - mu) / sigma
    # Density of Normal(mu, sigma) at this length, normalized to ~AA
    # length scale rather than continuous probability — that's fine
    # for ranking. Use unit-width bin: p ≈ pdf(z) / sigma.
    pdf = math.exp(-0.5 * z * z) / (sigma * math.sqrt(2 * math.pi))
    return math.log(max(pdf, 1e-30))


def _n_inserted_log(n_inserted: int, chain: str) -> float:
    """log P(n_inserted total nucleotides). β: sum of two geometrics;
    α: single geometric."""
    n_inserted = max(0, int(n_inserted))
    p = _N_INSERT_GEOM_P[chain]
    if chain == "beta":
        # Sum of two iid geometric(p) → negative binomial(r=2, p).
        # P(N=n) = C(n+1, 1) * (1-p)^n * p^2
        return math.log((n_inserted + 1) * ((1 - p) ** n_inserted) * (p ** 2))
    else:
        # Geometric(p): P(N=n) = (1-p)^n * p
        return math.log(((1 - p) ** n_inserted) * p)


def _estimate_n_inserted(cdr3_aa: str, chain: str) -> int:
    """Estimate total non-templated nt from observed CDR3 AA length."""
    if not isinstance(cdr3_aa, str) or not cdr3_aa:
        return 0
    return max(0, 3 * len(cdr3_aa) - _TEMPLATED_NT_TYPICAL[chain])


# AA composition baseline — heuristic uniform-over-20-AA log prob per
# position. Better than nothing; OLGA uses a position-specific matrix.
_LOG_AA_BASELINE = math.log(1.0 / 20.0)


def _composition_log(cdr3_aa: str) -> float:
    """log P(CDR3 AA | length) — coarse uniform-AA baseline. Equivalent
    to ``length × log(1/20)``, which makes longer sequences less
    likely a priori (matches the OLGA-style decomposition)."""
    if not isinstance(cdr3_aa, str):
        return 0.0
    return len(cdr3_aa) * _LOG_AA_BASELINE


def pgen_components(
    cdr3_aa: str,
    v_gene: str | None,
    j_gene: str | None,
    *,
    chain: str = "beta",
    n_inserted: int | None = None,
) -> dict[str, float]:
    """Decompose log10 Pgen into its five terms.

    Useful for debugging: shows which term is driving a sequence
    toward "public" or "private".
    """
    chain = chain.lower()
    if chain not in ("alpha", "beta"):
        raise ValueError(f"chain must be 'alpha' or 'beta', got {chain!r}")
    if n_inserted is None:
        n_inserted = _estimate_n_inserted(cdr3_aa, chain)
    L = len(cdr3_aa) if isinstance(cdr3_aa, str) else 0
    LN10 = math.log(10.0)
    return {
        "log10_p_v":          _v_usage_log(v_gene or "", chain) / LN10,
        "log10_p_j":          _j_usage_log(j_gene or "", chain) / LN10,
        "log10_p_length":     _length_log(L, chain) / LN10,
        "log10_p_n_inserted": _n_inserted_log(n_inserted, chain) / LN10,
        "log10_p_composition": _composition_log(cdr3_aa) / LN10,
    }


def pgen_single(
    cdr3_aa: str,
    v_gene: str | None,
    j_gene: str | None,
    *,
    chain: str = "beta",
    n_inserted: int | None = None,
) -> float:
    """Return log10 Pgen estimate for one sequence."""
    parts = pgen_components(
        cdr3_aa, v_gene, j_gene, chain=chain, n_inserted=n_inserted
    )
    return float(sum(parts.values()))


def compute_pgen(
    df: pd.DataFrame,
    *,
    cdr3_col: str = "CDR3_beta",
    v_gene_col: str = "beta_v_gene",
    j_gene_col: str = "beta_j_gene",
    chain: str = "beta",
    n_inserted_col: str | None = None,
    output_col: str = "log10_pgen",
) -> pd.Series:
    """Compute log10 Pgen estimate for each row in ``df``.

    Parameters
    ----------
    df : pd.DataFrame
        Per-clone frame with at least the ``cdr3_col``.
    cdr3_col : str
        Column with CDR3 amino-acid sequence.
    v_gene_col, j_gene_col : str
        Columns with V / J gene calls (allele suffix tolerated).
    chain : str
        ``"alpha"`` or ``"beta"``.
    n_inserted_col : str | None
        Optional column with pre-computed total N-insertion length
        (nt). When None, estimated from CDR3 AA length and a typical
        templated contribution from :data:`_TEMPLATED_NT_TYPICAL`.
    output_col : str
        Not currently used (the returned Series is unnamed).

    Returns
    -------
    pd.Series
        log10 Pgen estimate per row, indexed like ``df``.
    """
    chain = chain.lower()
    if chain not in ("alpha", "beta"):
        raise ValueError(f"chain must be 'alpha' or 'beta', got {chain!r}")
    if cdr3_col not in df.columns:
        raise ValueError(f"compute_pgen: missing {cdr3_col!r} column")

    cdr3 = df[cdr3_col]
    v = df[v_gene_col] if v_gene_col in df.columns else pd.Series([None] * len(df), index=df.index)
    j = df[j_gene_col] if j_gene_col in df.columns else pd.Series([None] * len(df), index=df.index)
    if n_inserted_col and n_inserted_col in df.columns:
        n_ins = df[n_inserted_col]
    else:
        n_ins = pd.Series([None] * len(df), index=df.index)

    out = []
    for c, vg, jg, ni in zip(cdr3, v, j, n_ins):
        if not isinstance(c, str) or not c:
            out.append(np.nan)
            continue
        ni_val = int(ni) if ni is not None and not pd.isna(ni) else None
        out.append(pgen_single(c, vg, jg, chain=chain, n_inserted=ni_val))
    return pd.Series(out, index=df.index, name=output_col)


def publicness_score(
    log10_pgen: pd.Series | Iterable[float],
    *,
    low_pgen_cutoff: float = -30.0,
    high_pgen_cutoff: float = -18.0,
) -> pd.Series:
    """Map log10 Pgen → [0, 1] publicness in a monotone-decreasing way.

    - ``log10 Pgen ≥ high_pgen_cutoff`` (very generatable) → 1.0
    - ``log10 Pgen ≤ low_pgen_cutoff`` (very private) → 0.0
    - linear interpolation in between.

    A high publicness score means the sequence is likely to be public
    and DB matches against it should be discounted.

    Defaults are calibrated to *this module's* coarse Pgen estimator,
    not OLGA. Typical β CDR3s land around log10 ≈ −25; very public
    short sequences with canonical V/J around log10 ≈ −18. The
    estimator is uniformly more negative than OLGA because its AA
    composition term is uniform-over-20 rather than position-specific.
    Override the cutoffs if you're feeding in real OLGA values
    (use ~−20 / −8 for OLGA).
    """
    arr = pd.Series(np.asarray(list(log10_pgen), dtype=float))
    width = high_pgen_cutoff - low_pgen_cutoff
    if width <= 0:
        raise ValueError("high_pgen_cutoff must be > low_pgen_cutoff")
    clipped = arr.clip(lower=low_pgen_cutoff, upper=high_pgen_cutoff)
    score = (clipped - low_pgen_cutoff) / width
    score = score.where(arr.notna(), other=np.nan)
    return score


def annotate_publicness(
    df: pd.DataFrame,
    *,
    cdr3_col: str = "CDR3_beta",
    v_gene_col: str = "beta_v_gene",
    j_gene_col: str = "beta_j_gene",
    chain: str = "beta",
    pgen_col: str = "log10_pgen",
    publicness_col: str = "publicness",
) -> pd.DataFrame:
    """Augment ``df`` with ``log10_pgen`` + ``publicness`` columns.

    Convenience wrapper combining :func:`compute_pgen` and
    :func:`publicness_score`. Does not modify the input frame.
    """
    out = df.copy()
    out[pgen_col] = compute_pgen(
        out,
        cdr3_col=cdr3_col,
        v_gene_col=v_gene_col,
        j_gene_col=j_gene_col,
        chain=chain,
    )
    out[publicness_col] = publicness_score(out[pgen_col])
    return out
