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

"""Data-driven CDR3 sequence-probability models — the publicness axis.

Replaces the brittle OLGA/SONIA runtime path (:mod:`tcrsift.olga_ppost`)
for the **precursor-frequency / publicness** axis. Instead of a fixed,
allele-masked GPL prior, a background generation/occurrence model is **fit
once on an external reference repertoire** and reused; ``log_pgen(seq)`` is
then a fast, dependency-light, calibrated score for "how generatable /
common is this CDR3" — lower = more private / rarer precursor.

Two interchangeable backends behind one :class:`SequenceProbabilityModel`
interface:

- :class:`KmerProbabilityModel` — an order-``k`` Markov model over CDR3
  amino acids (numpy-only, no GPL, the default). The shipped default models
  in :mod:`tcrsift.refseqs` are fit offline on OLGA-generated synthetic
  repertoires (OLGA used *once at build time* to produce training
  sequences — never at runtime, so tcrsift stays Apache-2.0).
- :class:`TCRpegProbabilityModel` — wraps **TCRpeg** (Jiang & Li 2023), an
  autoregressive deep model. Optional extra: ``pip install tcrsift[tcrpeg]``.

Both are trained on an external reference (not the experiment's own clones)
so the probability is a genuine background, not circular with the selection
target.
"""

from __future__ import annotations

import abc
import logging
from collections.abc import Iterable

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# The 20 canonical amino acids, plus begin-of-sequence / end-of-sequence
# sentinels. CDR3 length is modelled implicitly via the EOS emission, so no
# separate length term is needed.
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_AA_INDEX = {aa: i for i, aa in enumerate(AA_ALPHABET)}
BOS = len(AA_ALPHABET)        # 20
EOS = len(AA_ALPHABET) + 1    # 21
N_SYM = len(AA_ALPHABET) + 2  # 22


class SequenceProbabilityModel(abc.ABC):
    """A fittable per-sequence log-probability model over CDR3 strings."""

    @abc.abstractmethod
    def fit(self, sequences: Iterable[str]) -> SequenceProbabilityModel:
        """Train on an iterable of CDR3 amino-acid strings. Returns self."""

    @abc.abstractmethod
    def log_prob(self, sequences: Iterable[str]) -> np.ndarray:
        """Natural-log probability per sequence (NaN for unscorable input)."""

    @abc.abstractmethod
    def save(self, path) -> None:
        """Persist the fitted model to ``path``."""

    @classmethod
    @abc.abstractmethod
    def load(cls, path) -> SequenceProbabilityModel:
        """Load a model previously written by :meth:`save`."""


def _encode(seq: str) -> list[int] | None:
    """CDR3 string → list of AA indices, or None if it has a non-AA char."""
    if not isinstance(seq, str) or not seq:
        return None
    out = []
    for ch in seq:
        idx = _AA_INDEX.get(ch)
        if idx is None:
            return None
        out.append(idx)
    return out


class KmerProbabilityModel(SequenceProbabilityModel):
    """Order-``k`` Markov model over CDR3 amino acids (numpy-only).

    ``log P(CDR3) = Σ_i log P(a_i | a_{i-k} … a_{i-1})`` with the sequence
    padded by ``order`` BOS sentinels and terminated by EOS, so both the
    composition *and* the length are captured. Add-``alpha`` (Laplace)
    smoothing keeps unseen contexts from giving ``-inf``.

    Parameters are a dense ``(N_SYM**order, N_SYM)`` log-probability table,
    compact enough to ship: the shipped defaults are order 2 (~20-32 KB
    float32 per chain); order 3 would be ~1 MB.
    """

    def __init__(self, *, order: int = 2, alpha: float = 1.0, chain: str = ""):
        if order < 1:
            raise ValueError(f"order must be >= 1, got {order}")
        self.order = int(order)
        self.alpha = float(alpha)
        self.chain = chain
        self.n_train = 0
        self.n_contexts = N_SYM**self.order
        self._logp: np.ndarray | None = None  # (n_contexts, N_SYM)

    # -- context id helpers ------------------------------------------------
    def _context_id(self, ctx: list[int]) -> int:
        cid = 0
        for s in ctx:
            cid = cid * N_SYM + s
        return cid

    def fit(self, sequences: Iterable[str]) -> KmerProbabilityModel:
        counts = np.zeros((self.n_contexts, N_SYM), dtype=np.float64)
        n = 0
        skipped = 0
        for seq in sequences:
            ids = _encode(seq)
            if ids is None:
                skipped += 1
                continue
            padded = [BOS] * self.order + ids + [EOS]
            for i in range(self.order, len(padded)):
                cid = self._context_id(padded[i - self.order:i])
                counts[cid, padded[i]] += 1.0
            n += 1
        if n == 0:
            raise ValueError("KmerProbabilityModel.fit: no scorable sequences")
        counts += self.alpha
        totals = counts.sum(axis=1, keepdims=True)
        self._logp = np.log(counts / totals).astype(np.float32)
        self.n_train = n
        if skipped:
            logger.info(
                "KmerProbabilityModel.fit: skipped %d non-AA sequences", skipped
            )
        return self

    def _log_prob_one(self, seq: str) -> float:
        ids = _encode(seq)
        if ids is None:
            return float("nan")
        padded = [BOS] * self.order + ids + [EOS]
        lp = 0.0
        for i in range(self.order, len(padded)):
            cid = self._context_id(padded[i - self.order:i])
            lp += float(self._logp[cid, padded[i]])
        return lp

    def log_prob(self, sequences: Iterable[str]) -> np.ndarray:
        if self._logp is None:
            raise RuntimeError("KmerProbabilityModel is not fitted")
        return np.array([self._log_prob_one(s) for s in sequences], dtype=float)

    def save(self, path) -> None:
        if self._logp is None:
            raise RuntimeError("KmerProbabilityModel is not fitted")
        np.savez_compressed(
            path,
            logp=self._logp,
            order=np.int64(self.order),
            alpha=np.float64(self.alpha),
            n_train=np.int64(self.n_train),
            chain=np.array(self.chain),
        )

    @classmethod
    def load(cls, path) -> KmerProbabilityModel:
        with np.load(path, allow_pickle=False) as data:
            model = cls(
                order=int(data["order"]),
                alpha=float(data["alpha"]),
                chain=str(data["chain"]),
            )
            model._logp = data["logp"].astype(np.float32)
            model.n_train = int(data["n_train"])
        if model._logp.shape != (model.n_contexts, N_SYM):
            raise ValueError(
                f"loaded k-mer table shape {model._logp.shape} != expected "
                f"{(model.n_contexts, N_SYM)} for order {model.order}"
            )
        return model


def _fit_gene_marginal(genes, alpha: float = 1.0):
    """Laplace-smoothed log gene marginal → (dict, unseen-tail logprob).

    Genes are canonicalized first (shared :func:`tcrsift.genes.canonicalize_gene`),
    so format variants collapse. The tail is a single reserved back-off slot used
    as the prior for ANY unseen (or missing) gene at scoring time — so
    ``P(gene) > 0`` always, never a "gene not found" zero. Note it is a shared
    back-off prior, not a per-gene mass: the seen genes plus one tail slot sum to
    1, but assigning ``tail`` to several distinct unseen genes does not (this term
    is an additive ranking offset, not a normalized distribution over novel genes).
    """
    from collections import Counter

    from .genes import canonicalize_gene

    counts: Counter = Counter()
    for g in genes:
        cg = canonicalize_gene(g)
        if cg:
            counts[cg] += 1
    total = sum(counts.values())
    # +1 distinct slot for the unseen-gene tail.
    denom = total + alpha * (len(counts) + 1)
    if denom == 0:
        return {}, 0.0
    logp = {g: float(np.log((c + alpha) / denom)) for g, c in counts.items()}
    tail = float(np.log(alpha / denom))
    return logp, tail


class GeneAwareKmerModel(SequenceProbabilityModel):
    """Gene-aware Ppost: ``logP(V) + logP(J) + logP_kmer(CDR3)``.

    All three terms come from one reference repertoire: V and J from
    Laplace-smoothed gene marginals (so an unseen gene gets a finite tail
    probability, never zero), and the CDR3 from an order-``k`` Markov model
    (default **order 2** — most data-efficient below ~10⁵ reference seqs).
    Gene names are canonicalized via :func:`tcrsift.genes.canonicalize_gene`
    so format variants (alleles, Adaptive, ``TRAV14/DV4`` vs ``TRAV14DV4``)
    all resolve.

    ``log_prob`` scores gene-aware when V/J are supplied and degrades to the
    CDR3-only k-mer score when they aren't. NB (the caveat we agreed): don't
    benchmark a gene-aware score with TRAV12-2 AUROC — the label *is* the
    V gene; validate against a V-gene-independent publicness label.
    """

    def __init__(self, *, order: int = 2, alpha: float = 1.0, chain: str = ""):
        self.order = int(order)
        self.alpha = float(alpha)
        self.chain = chain
        self.n_train = 0
        self._cdr3 = KmerProbabilityModel(order=order, alpha=alpha, chain=chain)
        self._v_logp: dict[str, float] = {}
        self._j_logp: dict[str, float] = {}
        self._v_tail = 0.0
        self._j_tail = 0.0

    def fit(self, sequences, v_genes=None, j_genes=None) -> GeneAwareKmerModel:
        seqs = list(sequences)
        self._cdr3.fit(seqs)
        if v_genes is not None:
            self._v_logp, self._v_tail = _fit_gene_marginal(v_genes, self.alpha)
        if j_genes is not None:
            self._j_logp, self._j_tail = _fit_gene_marginal(j_genes, self.alpha)
        self.n_train = self._cdr3.n_train
        return self

    def _gene_term(self, genes, logp: dict[str, float], tail: float) -> np.ndarray:
        from .genes import canonicalize_gene

        # A missing/blank gene backs off to the unseen-gene tail prior, exactly
        # like an unrecognized gene — returning NaN here would poison the whole
        # gene-aware score (CDR3 + other gene included) for that clone.
        return np.array(
            [logp.get(canonicalize_gene(g), tail) for g in genes],
            dtype=float,
        )

    def log_prob(self, sequences, v_genes=None, j_genes=None) -> np.ndarray:
        lp = self._cdr3.log_prob(sequences)
        if v_genes is not None and self._v_logp:
            lp = lp + self._gene_term(v_genes, self._v_logp, self._v_tail)
        if j_genes is not None and self._j_logp:
            lp = lp + self._gene_term(j_genes, self._j_logp, self._j_tail)
        return lp

    @property
    def gene_aware(self) -> bool:
        return bool(self._v_logp or self._j_logp)

    def save(self, path) -> None:
        if self._cdr3._logp is None:
            raise RuntimeError("GeneAwareKmerModel is not fitted")
        np.savez_compressed(
            path,
            kind=np.array("gene_aware_kmer"),
            logp=self._cdr3._logp,
            order=np.int64(self.order),
            alpha=np.float64(self.alpha),
            n_train=np.int64(self.n_train),
            chain=np.array(self.chain),
            v_names=np.array(list(self._v_logp.keys())),
            v_logp=np.array(list(self._v_logp.values()), dtype=float),
            v_tail=np.float64(self._v_tail),
            j_names=np.array(list(self._j_logp.keys())),
            j_logp=np.array(list(self._j_logp.values()), dtype=float),
            j_tail=np.float64(self._j_tail),
        )

    @classmethod
    def load(cls, path) -> GeneAwareKmerModel:
        with np.load(path, allow_pickle=False) as data:
            model = cls(order=int(data["order"]), alpha=float(data["alpha"]),
                        chain=str(data["chain"]))
            model._cdr3._logp = data["logp"].astype(np.float32)
            model._cdr3.n_train = int(data["n_train"])
            model.n_train = int(data["n_train"])
            # Tolerate older CDR3-only k-mer files (no gene marginals): they
            # load and score CDR3-only until refit with V/J.
            if "v_names" in data:
                model._v_logp = dict(zip(data["v_names"].tolist(),
                                         data["v_logp"].tolist()))
                model._v_tail = float(data["v_tail"])
            if "j_names" in data:
                model._j_logp = dict(zip(data["j_names"].tolist(),
                                         data["j_logp"].tolist()))
                model._j_tail = float(data["j_tail"])
        return model


class TCRpegProbabilityModel(SequenceProbabilityModel):
    """TCRpeg-backed CDR3 probability (optional ``[tcrpeg]`` extra).

    Wraps the autoregressive TCRpeg model (Jiang & Li 2023). Heavier
    (PyTorch) but better-calibrated than the k-mer Markov model. Trained on
    the same external reference. Lazy import; raises :class:`ImportError`
    with an install hint when the extra is missing.
    """

    _INSTALL_HINT = (
        "TCRpeg (+ torch) is required for the TCRpeg sequence-probability "
        "backend but is not installed. Install with:\n\n"
        "    pip install tcrsift[tcrpeg]\n\n"
        "Or use the numpy-only KmerProbabilityModel (the default backend)."
    )

    def __init__(
        self,
        *,
        max_length: int = 30,
        embedding_size: int = 32,
        hidden_size: int = 64,
        num_layers: int = 1,
        device: str = "cpu",
        epochs: int = 20,
        batch_size: int = 1000,
        lr: float = 1e-3,
        chain: str = "",
        embedding_path: str | None = None,
    ):
        self.max_length = max_length
        self.embedding_size = embedding_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.device = device
        self.epochs = epochs
        self.batch_size = batch_size
        self.lr = lr
        self.chain = chain
        self.embedding_path = embedding_path
        self.n_train = 0
        self._model = None

    @staticmethod
    def available() -> bool:
        import importlib.util

        return (
            importlib.util.find_spec("tcrpeg") is not None
            and importlib.util.find_spec("torch") is not None
        )

    def _require(self) -> None:
        if not self.available():
            raise ImportError(self._INSTALL_HINT)

    def _resolve_embedding_path(self) -> str:
        """Absolute path to the AA embedding TCRpeg needs.

        TCRpeg's default ``embedding_path`` is relative to the CWD; the file
        actually ships inside the installed ``tcrpeg/data/`` dir, so resolve
        to that bundled copy unless the caller gave an explicit path.
        """
        if self.embedding_path is not None:
            return self.embedding_path
        import os as _os

        import tcrpeg  # pylint: disable=import-error

        return _os.path.join(
            _os.path.dirname(tcrpeg.__file__),
            "data", f"embedding_{self.embedding_size}.txt",
        )

    def _new_model(self, sequences: list[str] | None = None):
        from tcrpeg.TCRpeg import TCRpeg  # pylint: disable=import-error

        model = TCRpeg(
            max_length=self.max_length,
            embedding_size=self.embedding_size,
            hidden_size=self.hidden_size,
            num_layers=self.num_layers,
            device=self.device,
            load_data=sequences is not None,
            path_train=sequences,
            embedding_path=self._resolve_embedding_path(),
        )
        model.create_model()
        return model

    def fit(self, sequences: Iterable[str]) -> TCRpegProbabilityModel:
        self._require()
        seqs = [s for s in sequences if isinstance(s, str) and s]
        if not seqs:
            raise ValueError("TCRpegProbabilityModel.fit: no scorable sequences")
        self._model = self._new_model(seqs)
        self._model.train_tcrpeg(
            epochs=self.epochs, batch_size=self.batch_size, lr=self.lr,
        )
        self.n_train = len(seqs)
        return self

    def log_prob(self, sequences: Iterable[str]) -> np.ndarray:
        self._require()
        if self._model is None:
            raise RuntimeError("TCRpegProbabilityModel is not fitted")
        seqs = list(sequences)
        scorable = [isinstance(s, str) and bool(s) and _encode(s) is not None
                    for s in seqs]
        out = np.full(len(seqs), np.nan)
        good = [s for s, ok in zip(seqs, scorable) if ok]
        if good:
            # sampling_tcrpeg already returns natural-log probabilities.
            logs = np.asarray(self._model.sampling_tcrpeg(good), dtype=float)
            j = 0
            for i, ok in enumerate(scorable):
                if ok:
                    out[i] = logs[j]
                    j += 1
        return out

    def save(self, path) -> None:
        self._require()
        if self._model is None:
            raise RuntimeError("TCRpegProbabilityModel is not fitted")
        self._model.save(str(path))

    @classmethod
    def load(cls, path, **kwargs) -> TCRpegProbabilityModel:
        obj = cls(**kwargs)
        obj._require()
        from tcrpeg.TCRpeg import TCRpeg  # pylint: disable=import-error

        model = TCRpeg(
            max_length=obj.max_length,
            embedding_size=obj.embedding_size,
            hidden_size=obj.hidden_size,
            num_layers=obj.num_layers,
            device=obj.device,
            embedding_path=obj._resolve_embedding_path(),
        )
        model.create_model(load=True, path=str(path))
        obj._model = model
        return obj


# Default is the **gene-agnostic** order-2 CDR3 k-mer ("kmer"). Leave-one-
# donor-out CV on a V-gene-independent label (cross-donor recurrence) shows
# adding V/J marginals doesn't help (α 0.795 vs 0.798, β 0.828 vs 0.838) and
# V-only is ≤ chance — so gene-aware ("kmer_gene") is available but OFF by
# default, and the V gene is better used as a separate specificity feature.
# Both views read the same shipped npz (it carries CDR3 table + V/J marginals).
BACKENDS: dict[str, type[SequenceProbabilityModel]] = {
    "kmer": KmerProbabilityModel,
    "kmer_gene": GeneAwareKmerModel,
    "tcrpeg": TCRpegProbabilityModel,
}

# Below this many unique reference CDR3s the order-2 k-mer is more data-efficient
# than the TCRpeg GRU (which undertrains under ~10^5 seqs); only switch to
# tcrpeg once an OTS-scale reference is cached (#160).
KMER_REFERENCE_MAX = 30_000


def select_backend_for_reference_size(
    n_unique_cdr3: int, *, kmer_max: int = KMER_REFERENCE_MAX
) -> str:
    """Pick the Pgen/Ppost backend by reference size (#160).

    Returns ``"kmer"`` for small references (``< kmer_max`` unique CDR3s) — the
    data-efficient, always-positive default — and ``"tcrpeg"`` only once the
    reference is OTS-scale, where the GRU starts to win. Don't default
    ``pgen train`` to TCRpeg unconditionally.
    """
    return "kmer" if int(n_unique_cdr3) < int(kmer_max) else "tcrpeg"


def _default_model_path(chain: str, backend: str, role: str = "ppost"):
    """Path to a shipped default model under :mod:`tcrsift.refseqs`."""
    from importlib.resources import files

    chain = chain.lower()
    if chain not in ("alpha", "beta"):
        raise ValueError(f"chain must be 'alpha' or 'beta', got {chain!r}")
    if role not in ("pgen", "ppost"):
        raise ValueError(f"role must be 'pgen' or 'ppost', got {role!r}")
    if backend != "kmer":
        raise ValueError(f"no shipped default model for backend {backend!r}")
    return files("tcrsift.refseqs").joinpath(f"kmer_{role}_{chain}.npz")


_MODEL_CACHE: dict[tuple[str, str, str], SequenceProbabilityModel] = {}


def load_background_model(
    chain: str = "beta", backend: str = "kmer", role: str = "ppost",
) -> SequenceProbabilityModel:
    """Load (and cache) a shipped default background model.

    ``role`` is ``"ppost"`` (default — fit on an *observed* repertoire, the
    post-selection publicness measure) or ``"pgen"`` (fit on an
    OLGA-generated reference, pre-selection generation probability). Only the
    ``"kmer"`` backend ships defaults. **Role-pure**: raises
    :class:`FileNotFoundError` when the requested role isn't shipped for that
    chain (e.g. no observed-α ppost) — it never silently returns the Pgen
    model in place of Ppost. Callers decide how to degrade.
    """
    chain = chain.lower()
    key = (chain, backend, role)
    if key in _MODEL_CACHE:
        return _MODEL_CACHE[key]
    path = _default_model_path(chain, backend, role)
    if not path.is_file():
        raise FileNotFoundError(
            f"no shipped {backend} {role} model for chain {chain!r} at {path}"
        )
    model = BACKENDS[backend].load(str(path))
    _MODEL_CACHE[key] = model
    return model


def score_log_prob(
    df: pd.DataFrame,
    *,
    chain: str = "beta",
    cdr3_col: str | None = None,
    v_gene_col: str | None = None,
    j_gene_col: str | None = None,
    backend: str = "kmer",
    role: str = "ppost",
    model: SequenceProbabilityModel | None = None,
    out_col: str | None = None,
) -> pd.Series:
    """Per-clone natural-log probability under a background model.

    ``role="ppost"`` (default) scores against the observed-repertoire model
    (post-selection publicness); ``role="pgen"`` against the generated model.
    Uses ``model`` if given, else the shipped default for ``(chain, backend,
    role)``. ``cdr3_col`` defaults to ``CDR3_<chain>``.

    When the model is **gene-aware** (a :class:`GeneAwareKmerModel` with V/J
    marginals) and the V/J columns are present (default ``<chain>_v_gene`` /
    ``<chain>_j_gene``), the score includes ``logP(V) + logP(J)`` — gene names
    are canonicalized inside the model. Returns a Series aligned to ``df``.
    """
    cdr3_col = cdr3_col or f"CDR3_{chain}"
    if cdr3_col not in df.columns:
        raise ValueError(f"score_log_prob: missing {cdr3_col!r} column")
    if model is None:
        model = load_background_model(chain, backend, role)
    cdr3 = df[cdr3_col].tolist()
    if getattr(model, "gene_aware", False):
        v_col = v_gene_col or f"{chain}_v_gene"
        j_col = j_gene_col or f"{chain}_j_gene"
        v = df[v_col].tolist() if v_col in df.columns else None
        j = df[j_col].tolist() if j_col in df.columns else None
        scores = model.log_prob(cdr3, v_genes=v, j_genes=j)
    else:
        scores = model.log_prob(cdr3)
    return pd.Series(scores, index=df.index, name=out_col or f"log_{role}")


# Back-compat alias: log_pgen scoring is score_log_prob(role="pgen").
def score_log_pgen(
    df: pd.DataFrame,
    *,
    chain: str = "beta",
    cdr3_col: str | None = None,
    backend: str = "kmer",
    model: SequenceProbabilityModel | None = None,
    out_col: str = "log_pgen",
) -> pd.Series:
    """Per-clone log Pgen (generated-repertoire background). See
    :func:`score_log_prob`."""
    return score_log_prob(
        df, chain=chain, cdr3_col=cdr3_col, backend=backend, role="pgen",
        model=model, out_col=out_col,
    )


def score_log_q(
    df: pd.DataFrame,
    *,
    chain: str = "beta",
    cdr3_col: str | None = None,
    backend: str = "kmer",
    out_col: str | None = None,
) -> pd.Series:
    """Per-clone log selection factor ``log Q = log Ppost − log Pgen``.

    The data-driven selection factor (route 1): the log-ratio of the
    observed-repertoire model to the generated model. Positive ⇒ enriched by
    selection relative to pure generation. Unnormalized (a per-clone
    constant offset), which is irrelevant for ranking.
    """
    lp_post = score_log_prob(df, chain=chain, cdr3_col=cdr3_col,
                             backend=backend, role="ppost")
    lp_gen = score_log_prob(df, chain=chain, cdr3_col=cdr3_col,
                            backend=backend, role="pgen")
    q = lp_post - lp_gen
    q.name = out_col or f"log_q_{chain}"
    return q
