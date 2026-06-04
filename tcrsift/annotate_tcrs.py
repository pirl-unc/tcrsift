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

"""Uniform per-clone TCR annotation API + PRISM selection (#158).

Adds optional, **always-defined** columns to any table of TCRs — generated
sequences, ``clonotypes.csv``, filtered/selected clones — and a named
selection method, **PRISM**, built on them. Every annotation is
independently toggleable; the sequence-probability columns never return
0/NaN for a valid junction (#157), unlike the OLGA/SONIA runtime path.

Annotations
-----------
- **Pgen / Ppost** per chain (``pgen_alpha``/``pgen_beta``,
  ``ppost_alpha``/``ppost_beta``) via :mod:`tcrsift.seqprob`. **Ppost
  (observed-repertoire background) is the default publicness measure**;
  ``log_q_<chain> = log Ppost − log Pgen`` is the data-driven selection
  factor.
- **GEX signature scores** (``antigen_response_score``, ``naive_score``)
  from per-cell expression, **z-scored within each sample or donor group**
  (#144/#145), then averaged per clone.

PRISM (Percentile-Rank In-Silico Multi-criterion)
--------------------------------------------------
Rank every clone by the **mean of percentile ranks** of:
low ``ppost_alpha``, low ``ppost_beta``, high ``antigen_response_score``,
low ``naive_score`` — then take the top-K. On the B1-2/B1-3 pilot PRISM
selected ~33% / 45% condition-private clones vs 6% / 17% for frequency.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import seqprob
from .insilico_filter import FilterPredicate, percentile_rank

logger = logging.getLogger(__name__)

# PRISM's default criteria (#158): low α/β Ppost, high antigen-response,
# low naive. Each is a (column, direction) lower-is-better percentile dim.
PRISM_DEFAULT_PREDICATES = [
    FilterPredicate("ppost_alpha", "low", 1.0),
    FilterPredicate("ppost_beta", "low", 1.0),
    FilterPredicate("antigen_response_score", "high", 1.0),
    FilterPredicate("naive_score", "low", 1.0),
]


def naive_signature():
    """The naïve/stem signature used for ``naive_score`` (#141/#145).

    The down-pole of the Differentiated contrast: TCF7, LEF1, CCR7, SELL,
    IL7R, CD27, CD28.
    """
    from .signature_methods import Signature

    return Signature(
        "Naive",
        ("TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28"),
        panel="broad",
        description="Naïve/stem-memory program (#141/#145) for naive_score.",
    )


def add_pgen_ppost(
    df: pd.DataFrame,
    *,
    chains: tuple[str, ...] = ("alpha", "beta"),
    backend: str = "kmer",
    cdr3_cols: dict[str, str] | None = None,
    with_pgen: bool = True,
    with_ppost: bool = True,
    with_q: bool = True,
    auto_train: bool = True,
) -> pd.DataFrame:
    """Add ``pgen_<chain>`` / ``ppost_<chain>`` / ``log_q_<chain>`` columns.

    Ppost (observed-repertoire background) is the publicness default; Pgen is
    the OLGA-generated background; ``log_q`` is their difference (the
    selection factor). ``cdr3_cols`` maps chain → CDR3 column (default
    ``CDR3_<chain>``). Chains whose CDR3 column is absent are skipped.
    Returns a copy.
    """
    import numpy as np

    from .pgen_models import ensure_model

    out = df.copy()
    cdr3_cols = cdr3_cols or {}

    def _resolve(chain: str, role: str):
        """(model, estimator) honestly, or (None, None) if unavailable.

        Pgen may fall back across *estimators* (tcrpeg→k-mer) — still genuine
        Pgen. Ppost is role-pure: never substituted by Pgen, so an absent
        Ppost stays unavailable (→ NaN) rather than masquerading as Pgen.
        """
        try:
            return ensure_model(chain, backend=backend, role=role,
                                auto_train=auto_train), backend
        except ImportError as exc:  # pgen training needs OLGA
            if role == "pgen":
                logger.warning("%s; using shipped k-mer Pgen for %s", exc, chain)
                from .seqprob import load_background_model
                return load_background_model(chain, "kmer", "pgen"), "kmer"
            return None, None
        except FileNotFoundError as exc:
            if role == "pgen" and backend != "kmer":
                from .seqprob import load_background_model
                try:
                    return load_background_model(chain, "kmer", "pgen"), "kmer"
                except FileNotFoundError:
                    pass
            logger.warning("add_pgen_ppost: %s %s unavailable for %s (%s) → NaN",
                           backend, role, chain, exc)
            return None, None

    for chain in chains:
        cdr3_col = cdr3_cols.get(chain, f"CDR3_{chain}")
        if cdr3_col not in out.columns:
            logger.info("add_pgen_ppost: no %r column; skipping %s",
                        cdr3_col, chain)
            continue
        pgen_est = ppost_est = None
        if with_pgen:
            model, pgen_est = _resolve(chain, "pgen")
            out[f"pgen_{chain}"] = (
                seqprob.score_log_prob(out, chain=chain, cdr3_col=cdr3_col,
                                       model=model, out_col=f"pgen_{chain}")
                if model is not None else np.nan
            )
        if with_ppost:
            model, ppost_est = _resolve(chain, "ppost")
            out[f"ppost_{chain}"] = (
                seqprob.score_log_prob(out, chain=chain, cdr3_col=cdr3_col,
                                       model=model, out_col=f"ppost_{chain}")
                if model is not None else np.nan
            )
        # Q = log Ppost − log Pgen is only meaningful within one estimator;
        # don't emit a cross-estimator ratio. (Rank on Ppost, not Q — Q alone
        # is the weaker signal; publicness lives in the observed frequency.)
        if with_q and with_pgen and with_ppost:
            if pgen_est is not None and pgen_est == ppost_est:
                out[f"log_q_{chain}"] = out[f"ppost_{chain}"] - out[f"pgen_{chain}"]
            else:
                out[f"log_q_{chain}"] = np.nan
    return out


def score_gex_signature_per_clone(
    per_cell: pd.DataFrame,
    signature,
    *,
    clone_col: str = "CDR3ab",
    group_col: str | None = None,
    gex_prefix: str = "gex",
    log1p: bool = True,
    out_col: str | None = None,
) -> pd.Series:
    """Per-clone GEX signature score, z-scored within sample/donor groups.

    ``per_cell`` is a per-cell frame with ``{gex_prefix}.<SYMBOL>`` columns,
    a ``clone_col`` and (optionally) a ``group_col`` (sample or donor). Each
    signature gene is z-scored across cells **within each group** (so a
    high-baseline sample doesn't dominate, #144/#145), averaged across the
    signature's genes per cell, then averaged across each clone's cells.

    ``signature`` is a :class:`tcrsift.signature_methods.Signature` or a name
    in its registry. Returns a per-clone Series indexed by clone id.
    """
    from .signature_methods import SIGNATURES, score_signature

    sig = SIGNATURES[signature] if isinstance(signature, str) else signature
    cols = {g: f"{gex_prefix}.{g}" for g in sig.all_genes}
    present = {g: c for g, c in cols.items() if c in per_cell.columns}
    if not present:
        raise ValueError(
            f"score_gex_signature_per_clone: none of {sig.name}'s genes "
            f"({list(sig.all_genes)}) found as {gex_prefix}.<gene> columns"
        )
    expr = per_cell[list(present.values())].copy()
    expr.columns = list(present.keys())  # rename to bare symbols
    groups = per_cell[group_col] if group_col and group_col in per_cell.columns else None
    per_cell_score = score_signature(
        expr, sig, combine="zscore", log1p=log1p, groups=groups,
        on_missing="ignore",
    )
    tmp = pd.DataFrame({"clone": per_cell[clone_col].values,
                        "score": per_cell_score.values})
    per_clone = tmp.groupby("clone", observed=True)["score"].mean()
    per_clone.name = out_col or f"{sig.name}_score"
    return per_clone


def add_gex_signature_scores(
    df: pd.DataFrame,
    per_cell: pd.DataFrame,
    *,
    signatures: dict[str, object] | None = None,
    clone_col: str = "CDR3ab",
    group_col: str | None = None,
    gex_prefix: str = "gex",
) -> pd.DataFrame:
    """Join per-clone, group-z-scored signature scores onto ``df``.

    ``signatures`` maps output column → signature (name or object); defaults
    to ``{"antigen_response_score": "AcuteActivation", "naive_score":
    "DifferentiatedNaive"}`` analogues. ``group_col`` is the sample/donor
    column in ``per_cell`` to z-score within. Returns a copy of ``df`` with
    the score columns joined on ``clone_col``.
    """
    if signatures is None:
        signatures = {
            "antigen_response_score": "AcuteActivation",
            "naive_score": naive_signature(),
        }
    out = df.copy()
    for out_col, sig in signatures.items():
        per_clone = score_gex_signature_per_clone(
            per_cell, sig, clone_col=clone_col, group_col=group_col,
            gex_prefix=gex_prefix, out_col=out_col,
        )
        out[out_col] = out[clone_col].map(per_clone)
    return out


def add_paired_ppost(
    df: pd.DataFrame,
    *,
    alpha_col: str = "ppost_alpha",
    beta_col: str = "ppost_beta",
    paired_col: str = "ppost_paired",
    either_col: str = "ppost_either",
    pmi: callable | None = None,
    pmi_weight: float = 0.0,
) -> pd.DataFrame:
    """Add integrated αβ publicness columns — two distinct, both exposed.

    - ``paired_col`` = ``log Pα + log Pβ`` (+ optional ``pmi_weight·PMI(Vα,Vβ)``):
      the **joint** generation log-prob, since α/β recombine independently
      (sum-of-logs). NaN if either chain is missing.
    - ``either_col`` = ``max(log Pα, log Pβ)``: the value to drive a
      **depletion decision** — a clone is public if *either* chain is common,
      so don't hide a public α behind a private β. Uses whichever chain is
      present when one is missing.

    ``pmi`` is an optional ``(v_alpha, v_beta) → float`` empirical pairing
    term (log PMI); with ``pmi_weight=0`` (default) pairing is ignored and
    rare/unseen pairs contribute nothing. Returns a copy.
    """
    out = df.copy()
    a = out[alpha_col] if alpha_col in out.columns else np.nan
    b = out[beta_col] if beta_col in out.columns else np.nan
    out[paired_col] = a + b
    if pmi is not None and pmi_weight:
        va = out.get("alpha_v_gene")
        vb = out.get("beta_v_gene")
        if va is not None and vb is not None:
            out[paired_col] = out[paired_col] + pmi_weight * pd.Series(
                [pmi(x, y) for x, y in zip(va, vb)], index=out.index)
    out[either_col] = pd.concat([a, b], axis=1).max(axis=1) \
        if (alpha_col in out.columns and beta_col in out.columns) \
        else (a if beta_col not in out.columns else b)
    return out


def add_pairing_promiscuity(
    df: pd.DataFrame,
    *,
    alpha_col: str = "CDR3_alpha",
    beta_col: str = "CDR3_beta",
    promiscuous_min_partners: int = 3,
) -> pd.DataFrame:
    """Add α-CDR3 β-pairing promiscuity as a publicness / cross-reactivity flag (#148).

    For each clonotype, ``alpha_beta_promiscuity`` is the number of *distinct*
    CDR3β chains its CDR3α pairs with across the table — a directly
    interpretable publicness signal (germline-like, degenerate public αs like
    Melan-A/TRAV12-2 pair with many βs) to **down-weight** in selection.
    ``alpha_promiscuous`` flags clones at/above ``promiscuous_min_partners``
    (default 3). ``alpha_cdr3_length`` is exposed too — α (not β) CDR3 length is
    the informative sequence feature here.

    Promiscuity is a tail phenomenon (median 1 partner): it flags the public
    degenerate set, not all clones. It is a publicness / cross-reactivity
    *proxy*, not a specificity or avidity measurement (cf. #143/#146). Returns
    a copy; a no-op (no columns added) when the CDR3 columns are absent.
    """
    out = df.copy()
    if alpha_col not in out.columns or beta_col not in out.columns:
        return out

    a = out[alpha_col].astype("string")
    b = out[beta_col].astype("string")

    def _valid(s):
        return s.notna() & (s.str.len() > 0) & (s.str.lower() != "nan")

    valid_a = _valid(a)
    pair_ok = valid_a & _valid(b)
    # Distinct β partners per α, counted over confidently paired clones only
    # (multi-β ambiguity would otherwise inflate counts, #148).
    partners = b[pair_ok].groupby(a[pair_ok], observed=True).nunique()
    out["alpha_beta_promiscuity"] = (
        a.map(partners).fillna(0).astype(int).where(valid_a, 0)
    )
    out["alpha_promiscuous"] = out["alpha_beta_promiscuity"] >= promiscuous_min_partners
    out["alpha_cdr3_length"] = a.str.len().where(valid_a).astype("Int64")
    return out


def prism_score(
    df: pd.DataFrame,
    *,
    predicates: list[FilterPredicate] | None = None,
    weights: list[float] | None = None,
    group_col: str | None = None,
    score_col: str = "prism_score",
    rank_col: str = "prism_rank",
) -> pd.DataFrame:
    """Compute the PRISM score: (weighted) mean of per-dimension percentile ranks.

    Each predicate contributes a lower-is-better percentile rank in [0, 1]
    (0 = best); PRISM averages them (optionally weighted) so a low score =
    a strong multi-criterion candidate. ``rank_col`` is the 1-based ordering
    (1 = best). When ``group_col`` is given, percentile ranks are computed
    within each group (e.g. per assay condition). Missing any dimension →
    NaN PRISM score (ranked last). Returns a copy with the two columns added.
    """
    predicates = predicates or PRISM_DEFAULT_PREDICATES
    missing = [p.score for p in predicates if p.score not in df.columns]
    if missing:
        raise ValueError(f"prism_score: missing predicate columns {missing}")
    if weights is None:
        weights = [1.0] * len(predicates)
    if len(weights) != len(predicates):
        raise ValueError("prism_score: weights length != predicates length")

    out = df.copy()
    grouped = group_col is not None and group_col in out.columns

    parts = []
    for pred in predicates:
        lower, method = pred.lower_is_better, pred.rank_method
        if grouped:
            pr = out.groupby(group_col, observed=True)[pred.score].transform(
                lambda s, _l=lower, _m=method: percentile_rank(
                    s, lower_is_better=_l, method=_m,
                )
            )
        else:
            pr = percentile_rank(out[pred.score], lower_is_better=lower,
                                 method=method)
        parts.append(pr)
    mat = pd.concat(parts, axis=1)
    w = np.asarray(weights, dtype=float)
    scores = (mat * w).sum(axis=1, skipna=False) / w.sum()

    out[score_col] = scores
    out[rank_col] = scores.rank(method="min", ascending=True, na_option="bottom")
    return out


def select_prism(
    df: pd.DataFrame,
    *,
    k: int,
    group_col: str | None = None,
    score_col: str = "prism_score",
    selected_col: str = "prism_selected",
    **prism_kwargs,
) -> pd.DataFrame:
    """Compute PRISM and flag the top-``k`` clones (per group if given).

    Returns a copy with the PRISM columns plus a boolean ``selected_col``.
    """
    out = prism_score(df, group_col=group_col, score_col=score_col,
                      **prism_kwargs)
    if group_col is not None and group_col in out.columns:
        rank = out.groupby(group_col, observed=True)[score_col].rank(
            method="first", ascending=True, na_option="bottom",
        )
    else:
        rank = out[score_col].rank(method="first", ascending=True,
                                   na_option="bottom")
    out[selected_col] = out[score_col].notna() & (rank <= k)
    return out


def annotate_tcrs(
    df: pd.DataFrame,
    *,
    methods: list[str] | None = None,
    per_cell: pd.DataFrame | None = None,
    backend: str = "kmer",
    chains: tuple[str, ...] = ("alpha", "beta"),
    clone_col: str = "CDR3ab",
    group_col: str | None = None,
    gex_prefix: str = "gex",
) -> pd.DataFrame:
    """Add the requested annotation columns to a TCR table (#158).

    ``methods`` subset of ``{"pgen", "ppost", "gex_signatures", "prism"}``
    (default all that are computable). ``per_cell`` (with ``{gex_prefix}.``
    columns + ``clone_col`` + ``group_col``) is required for
    ``gex_signatures``/``prism``. Returns a copy with the columns added.
    """
    methods = methods or ["pgen", "ppost", "gex_signatures", "prism"]
    out = df.copy()
    if "pgen" in methods or "ppost" in methods:
        out = add_pgen_ppost(
            out, chains=chains, backend=backend,
            with_pgen=("pgen" in methods or "prism" in methods),
            with_ppost=("ppost" in methods or "prism" in methods),
        )
    if ("gex_signatures" in methods or "prism" in methods) and per_cell is not None:
        out = add_gex_signature_scores(
            out, per_cell, clone_col=clone_col, group_col=group_col,
            gex_prefix=gex_prefix,
        )
    if "prism" in methods:
        have = [p.score for p in PRISM_DEFAULT_PREDICATES if p.score in out.columns]
        if len(have) == len(PRISM_DEFAULT_PREDICATES):
            out = prism_score(out, group_col=group_col)
        else:
            logger.warning(
                "annotate_tcrs: PRISM needs %s; have only %s — skipping PRISM",
                [p.score for p in PRISM_DEFAULT_PREDICATES], have,
            )
    return out


