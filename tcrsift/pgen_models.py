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

"""Manage Pgen/Ppost background models: fetch data, train, cache, auto-load.

The lifecycle layer over :mod:`tcrsift.seqprob`. Trained models and their
training repertoires are cached under ``~/.cache/tcrsift/seqprob`` (the
shared :func:`tcrsift.datacache.resolve_cache_dir`). The default backend is
the **gene-agnostic order-2 k-mer** (cross-donor CV shows it matches
gene-aware / TCRpeg, and it's correct at every reference size we can get);
TCRpeg is opt-in (worth it only at OTS scale). :func:`ensure_model`
auto-trains and caches on first use when a model isn't shipped/cached.

Training data, by role (kept strictly non-circular — never the experiment's
own clones):

- **pgen** — pre-selection generation probability. The shipped numpy k-mer
  Pgen model is used; it was generated offline (OLGA, no longer a dependency)
  and cannot be retrained at runtime — see refseqs/PROVENANCE.md.
- **ppost** — an **external observed neutral** repertoire (post-selection
  publicness). The shipped α/β defaults are fit on the **5 pooled public
  10x healthy-PBMC VDJ-T sets** (:func:`fetch_healthy_pbmc`, ≈16k α / 18k β
  unique CDR3). Scale up with **OTS** (Observed TCR Space, OPIG; paired-αβ,
  CC-BY-4.0) via ``fetch_repertoire(url=...)`` — that's what unlocks TCRpeg.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd

from . import seqprob

logger = logging.getLogger(__name__)

# Default is the gene-agnostic order-2 k-mer — correct at every reference
# size we can get, and cross-donor CV shows it matches gene-aware/TCRpeg.
# TCRpeg is opt-in (auto-select to it only makes sense at OTS scale).
DEFAULT_BACKEND = "kmer"

# Recommended neutral (observed, non-antigen-sorted) repertoire sources for
# the Ppost background. OTS is paired-αβ; download a per-study/filtered CSV
# from the search page and pass its path/URL to `fetch_repertoire`.
REPERTOIRE_SOURCES = {
    "healthy-pbmc": "5 pooled public 10x healthy-PBMC VDJ-T sets (auto-download)",
    "ots": "https://opig.stats.ox.ac.uk/webapps/ots",          # paired αβ, CC-BY-4.0
    "ireceptor": "https://gateway.ireceptor.org",              # healthy PBMC AIRR-seq
}

# The 5 public 10x healthy-PBMC paired VDJ-T sets (CC-licensed, direct
# download) pooled into the neutral α+β Ppost reference — single-cell paired
# VDJ, so both chains come from one non-antigen-enriched source (#160 §2).
HEALTHY_PBMC_10X = [
    ("sc5p_v2_hs_PBMC_10k",
     "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/"
     "sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/"
     "sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_filtered_contig_annotations.csv"),
    ("sc5p_v2_hs_PBMC_1k",
     "https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/"
     "sc5p_v2_hs_PBMC_1k_multi_5gex_5fb_b_t/"
     "sc5p_v2_hs_PBMC_1k_multi_5gex_5fb_b_t_vdj_t_filtered_contig_annotations.csv"),
    ("vdj_v1_hs_pbmc3",
     "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/"
     "vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv"),
    ("vdj_nextgem_hs_pbmc3",
     "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/"
     "vdj_nextgem_hs_pbmc3/vdj_nextgem_hs_pbmc3_t_filtered_contig_annotations.csv"),
    ("vdj_v1_hs_pbmc_t",
     "https://cf.10xgenomics.com/samples/cell-vdj/2.2.0/"
     "vdj_v1_hs_pbmc_t/vdj_v1_hs_pbmc_t_filtered_contig_annotations.csv"),
]

_STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# 10x's CDN 403s urllib's default User-Agent; send a browser-like one.
_HTTP_HEADERS = {"User-Agent": "Mozilla/5.0 (tcrsift)"}


def _download(url: str, dest: Path) -> None:
    """Download ``url`` → ``dest`` with a browser User-Agent (CDN-friendly)."""
    import shutil
    import urllib.request

    req = urllib.request.Request(url, headers=_HTTP_HEADERS)  # noqa: S310
    with urllib.request.urlopen(req) as resp, open(dest, "wb") as fh:  # noqa: S310
        shutil.copyfileobj(resp, fh)


def fetch_healthy_pbmc(*, force: bool = False) -> dict[str, Path]:
    """Download + filter + pool the 5 public 10x healthy-PBMC VDJ-T sets.

    The reproducible neutral reference for Ppost (#160 §2): per dataset, keep
    ``productive & high_confidence``, split by ``chain`` (TRA→α, TRB→β), take
    the ``cdr3`` amino-acid sequence (length 5–30, standard-20 AA only),
    capture ``v_gene``/``j_gene``, then pool across the 5 and dedup by CDR3.
    Caches per-chain ``seq,v,j`` CSVs (≈16k α / 18k β unique CDR3). Returns
    ``{chain: path}``. Never the experiment's own clones — healthy donors only.
    """
    self_dir = _data_dir()
    self_dir.mkdir(parents=True, exist_ok=True)
    buckets = {"alpha": [], "beta": []}
    for name, url in HEALTHY_PBMC_10X:
        raw = self_dir / f"_10x_{name}.csv"
        if force or not raw.is_file():
            logger.info("fetch_healthy_pbmc: downloading %s", name)
            _download(url, raw)
        frame = pd.read_csv(raw)
        prod = frame["productive"].astype(str).str.lower() == "true"
        hi = frame["high_confidence"].astype(str).str.lower() == "true"
        frame = frame[prod & hi].dropna(subset=["cdr3"])
        frame = frame[frame["cdr3"].str.len().between(5, 30)]
        frame = frame[frame["cdr3"].apply(lambda s: set(s) <= _STD_AA)]
        for chain, locus in (("alpha", "TRA"), ("beta", "TRB")):
            sub = frame[frame["chain"] == locus][["cdr3", "v_gene", "j_gene"]]
            buckets[chain].append(
                sub.rename(columns={"cdr3": "seq", "v_gene": "v", "j_gene": "j"}))
    out: dict[str, Path] = {}
    for chain, parts in buckets.items():
        pooled = pd.concat(parts).drop_duplicates(subset=["seq"])
        dest = cached_repertoire_file(chain)
        pooled.to_csv(dest, index=False)
        out[chain] = dest
        logger.info("fetch_healthy_pbmc: %d unique %s CDR3 → %s",
                    len(pooled), chain, dest)
    return out


def cache_dir() -> Path:
    from .datacache import resolve_cache_dir

    return resolve_cache_dir() / "seqprob"


def _models_dir() -> Path:
    return cache_dir() / "models"


def _data_dir() -> Path:
    return cache_dir() / "data"


def cached_model_file(backend: str, role: str, chain: str) -> Path:
    ext = "npz" if backend == "kmer" else "pt"
    return _models_dir() / f"{backend}_{role}_{chain}.{ext}"


def cached_repertoire_file(chain: str) -> Path:
    return _data_dir() / f"observed_{chain}.csv"


# ---------------------------------------------------------------------------
# Training data
# ---------------------------------------------------------------------------

def _bundled_observed_beta() -> list[str]:
    """Healthy β CDR3s bundled with TCRpeg (external, non-circular)."""
    import tcrpeg  # pylint: disable=import-error

    path = Path(tcrpeg.__file__).parent / "data" / "TCRs_train.csv"
    try:
        return pd.read_csv(path, compression="gzip")["seq"].dropna().tolist()
    except Exception:
        return pd.read_csv(path)["seq"].dropna().tolist()


def fetch_repertoire(
    chain: str,
    *,
    url: str | None = None,
    source: str | None = None,
    seq_col: str = "seq",
    v_col: str | None = None,
    chain_prefix: str | None = None,
) -> Path:
    """Fetch + cache an observed repertoire (CDR3s) for Ppost training.

    With ``url`` (a CSV/TSV, optionally gzipped, with a CDR3 column), the
    file is downloaded, the ``seq_col`` extracted (optionally filtered to the
    chain via a ``v_col`` whose values start with ``chain_prefix``, e.g.
    ``"TRA"``/``"TRB"``), and written to the per-chain cache. With no ``url``
    and ``chain == "beta"``, the bundled healthy β repertoire is cached as-is.
    Returns the cached CSV path. ``source`` is an informational preset key
    (see :data:`REPERTOIRE_SOURCES`).
    """
    dest = cached_repertoire_file(chain)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if url is None:
        if chain == "beta":
            seqs = _bundled_observed_beta()
            logger.info("fetch_repertoire: using bundled healthy β (%d seqs)",
                        len(seqs))
            pd.DataFrame({"seq": seqs}).to_csv(dest, index=False)
            return dest
        raise ValueError(
            f"fetch_repertoire: no bundled observed {chain} repertoire — pass "
            f"url= a neutral repertoire (e.g. OTS: {REPERTOIRE_SOURCES['ots']})"
        )

    logger.info("fetch_repertoire: downloading %s (%s)", url, source or "url")
    raw = _data_dir() / f"_download_{chain}{_suffix(url)}"
    raw.parent.mkdir(parents=True, exist_ok=True)
    _download(url, raw)
    opener = gzip.open if str(raw).endswith(".gz") else open
    # Key the separator off the resolved file extension, not a substring of the
    # whole URL (a query string or path segment containing ".tsv" on a real CSV
    # would otherwise pick the wrong separator).
    sep = "\t" if _suffix(url) in (".tsv", ".tsv.gz") else ","
    with opener(raw, "rt") as fh:
        frame = pd.read_csv(fh, sep=sep)
    if seq_col not in frame.columns:
        raise ValueError(
            f"fetch_repertoire: {seq_col!r} not in downloaded columns "
            f"{list(frame.columns)}"
        )
    if v_col and chain_prefix and v_col in frame.columns:
        frame = frame[frame[v_col].astype(str).str.startswith(chain_prefix)]
    seqs = frame[seq_col].dropna().astype(str)
    seqs = seqs[seqs.str.len() > 0]
    pd.DataFrame({"seq": seqs}).to_csv(dest, index=False)
    logger.info("fetch_repertoire: cached %d %s seqs → %s", len(seqs), chain, dest)
    return dest


def _suffix(url: str) -> str:
    for s in (".csv.gz", ".tsv.gz", ".csv", ".tsv", ".gz"):
        if url.endswith(s):
            return s
    return ".csv"


def _observed_seqs(chain: str, *, url: str | None = None) -> list[str]:
    """CDR3s for Ppost training: cached file → fetch (url/bundled)."""
    cached = cached_repertoire_file(chain)
    if cached.is_file() and url is None:
        return pd.read_csv(cached)["seq"].dropna().astype(str).tolist()
    path = fetch_repertoire(chain, url=url)
    return pd.read_csv(path)["seq"].dropna().astype(str).tolist()


# ---------------------------------------------------------------------------
# Train + cache
# ---------------------------------------------------------------------------

def train_model(
    chain: str,
    *,
    backend: str = DEFAULT_BACKEND,
    role: str = "ppost",
    url: str | None = None,
    order: int = 3,
    epochs: int = 30,
    batch_size: int = 2000,
    device: str | None = None,
) -> Path:
    """Train a (backend, role, chain) model and cache it; return its path.

    Only Ppost (observed-repertoire) models are trainable here. Pgen models are
    shipped pre-built (``refseqs/kmer_pgen_*.npz``); they were generated offline
    with OLGA, which is no longer a tcrsift dependency, so Pgen can't be
    (re)trained at runtime — use the shipped model.
    """
    if role == "ppost":
        seqs = _observed_seqs(chain, url=url)
    elif role == "pgen":
        raise ValueError(
            "Pgen models are shipped pre-built and cannot be trained at "
            "runtime (OLGA, the synthetic-sequence generator, is no longer a "
            "dependency). Use the shipped k-mer Pgen model, or regenerate "
            "offline — see refseqs/PROVENANCE.md."
        )
    else:
        raise ValueError(f"role must be 'pgen'/'ppost', got {role!r}")

    out = cached_model_file(backend, role, chain)
    out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("train_model: fitting %s %s %s on %d seqs → %s",
                backend, role, chain, len(seqs), out)
    if backend == "kmer":
        model = seqprob.KmerProbabilityModel(order=order, chain=chain).fit(seqs)
    elif backend == "tcrpeg":
        dev = device or _default_device()
        model = seqprob.TCRpegProbabilityModel(
            epochs=epochs, batch_size=batch_size, device=dev, chain=chain,
        ).fit(seqs)
    else:
        raise ValueError(f"unknown backend {backend!r}")
    model.save(str(out))
    return out


def _default_device() -> str:
    # CUDA when present, else CPU. NOT Apple MPS: TCRpeg's embedding is float64
    # and MPS rejects float64 (and MPS-device checkpoints don't reload), so
    # MPS is unsupported here — CPU is the portable, save/load-clean default.
    try:
        import torch  # pylint: disable=import-error

        if torch.cuda.is_available():
            return "cuda:0"
    except Exception:
        pass
    return "cpu"


# ---------------------------------------------------------------------------
# Resolve a usable model: cache → shipped → auto-train (→ fallback)
# ---------------------------------------------------------------------------

_CACHE: dict[tuple, seqprob.SequenceProbabilityModel] = {}


def ensure_model(
    chain: str,
    *,
    backend: str = DEFAULT_BACKEND,
    role: str = "ppost",
    auto_train: bool = True,
    url: str | None = None,
) -> seqprob.SequenceProbabilityModel:
    """Return the (backend, role, chain) model, training it on first use.

    Strictly role-pure: precedence is in-memory cache → on-disk cached model
    → shipped default (k-mer only) → auto-train (when ``auto_train``). It
    **never substitutes one role for another** — a Ppost request never falls
    back to a Pgen model (that would privilege Pgen and fake selection
    correction). Raises :class:`FileNotFoundError` (no model / no training
    data) or :class:`ValueError` (Pgen can't be trained at runtime — use the
    shipped k-mer model) so the caller can decide how to degrade honestly.
    Cross-*estimator* fallback within a role (e.g. TCRpeg-Pgen → k-mer-Pgen) is
    the caller's call, not this one's.
    """
    key = (chain.lower(), backend, role)
    if key in _CACHE:
        return _CACHE[key]

    path = cached_model_file(backend, role, chain)
    if path.is_file():
        model = seqprob.BACKENDS[backend].load(str(path))
        _CACHE[key] = model
        return model

    # Shipped k-mer defaults (instant, no training). Raises FileNotFoundError
    # for a role/chain with no shipped model (e.g. k-mer ppost alpha).
    if backend == "kmer":
        model = seqprob.load_background_model(chain, "kmer", role)
        _CACHE[key] = model
        return model

    if not auto_train:
        raise FileNotFoundError(
            f"no cached {backend} {role} model for {chain}; run "
            f"`tcrsift pgen train` or set auto_train=True"
        )

    # Auto-train (tcrpeg, ppost only). Propagates FileNotFoundError (no observed
    # data) and ValueError (Pgen isn't trainable at runtime) to the caller.
    if role == "ppost" and chain == "alpha" and url is None \
            and not cached_repertoire_file("alpha").is_file():
        raise FileNotFoundError(
            "no observed alpha repertoire for ppost; fetch a neutral paired "
            "repertoire (e.g. OTS) via `tcrsift pgen fetch --url ...`"
        )
    train_model(chain, backend=backend, role=role, url=url)
    model = seqprob.BACKENDS[backend].load(
        str(cached_model_file(backend, role, chain)))
    _CACHE[key] = model
    return model
