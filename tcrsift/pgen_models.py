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


def _download(url: str, dest: Path, *, timeout: float = 30.0) -> None:
    """Download ``url`` → ``dest`` atomically with scheme validation & timeout."""
    import os
    import shutil
    import tempfile
    import urllib.parse
    import urllib.request

    parsed = urllib.parse.urlparse(url)
    if parsed.scheme not in ("http", "https"):
        raise ValueError(
            f"Unsupported URL scheme {parsed.scheme!r} for {url!r}; only http/https allowed"
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers=_HTTP_HEADERS)  # noqa: S310
    with tempfile.NamedTemporaryFile(
        prefix=f".{dest.name}.", suffix=".tmp", delete=False, dir=dest.parent
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp_path, "wb") as fh:  # noqa: S310
            shutil.copyfileobj(resp, fh)
        os.replace(tmp_path, dest)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


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


def cached_model_file(
    backend: str, role: str, chain: str, *, species: str = "human"
) -> Path:
    from .species import normalize_species

    sp = normalize_species(species)
    ext = "npz" if backend in ("kmer", "kmer_gene") else "pt"
    tag = f"_{sp}" if sp != "human" else ""
    return _models_dir() / f"{backend}_{role}{tag}_{chain}.{ext}"


def cached_repertoire_file(chain: str, *, species: str = "human") -> Path:
    from .species import normalize_species

    sp = normalize_species(species)
    tag = f"_{sp}" if sp != "human" else ""
    return _data_dir() / f"observed{tag}_{chain}.csv"


# ---------------------------------------------------------------------------
# Training data
# ---------------------------------------------------------------------------

def _bundled_observed_pbmc(chain: str) -> pd.DataFrame | None:
    """Read bundled healthy-PBMC ``observed_pbmc_10x.csv.gz`` (TRA + TRB with V/J)."""
    from importlib.resources import files

    locus = "TRA" if chain.lower() == "alpha" else "TRB"
    try:
        ref = files("tcrsift.refseqs").joinpath("observed_pbmc_10x.csv.gz")
        if not ref.is_file():
            return None
        with gzip.open(str(ref), "rt") as fh:
            df = pd.read_csv(fh)
    except Exception:
        return None
    if not {"chain", "seq"} <= set(df.columns):
        return None
    sub = df[df["chain"] == locus].dropna(subset=["seq"]).copy()
    cols = [c for c in ("seq", "v", "j") if c in sub.columns]
    return sub[cols].reset_index(drop=True) if not sub.empty else None


def _bundled_observed_beta() -> list[str] | None:
    """Healthy β CDR3s bundled with TCRpeg (external, non-circular)."""
    try:
        import tcrpeg  # pylint: disable=import-error
    except ImportError:
        return None

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
    j_col: str | None = None,
    chain_prefix: str | None = None,
    species: str = "human",
) -> Path:
    """Fetch + cache an observed repertoire (CDR3s) for Ppost training.

    With ``url`` (a CSV/TSV or local file path, optionally gzipped, with a CDR3
    column), the file is loaded, ``seq_col`` (plus optional ``v_col``/``j_col``)
    extracted, and written to the per-chain/per-species cache. With no ``url``
    and ``species == "human"``, the bundled healthy 10x PBMC α/β repertoire
    (``tcrsift/refseqs/observed_pbmc_10x.csv.gz``) is used completely offline
    (falling back to TCRpeg's bundled β CDR3s if absent).
    """
    from .species import normalize_species

    sp = normalize_species(species)
    chain = chain.lower()
    dest = cached_repertoire_file(chain, species=sp)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if url is None:
        if sp == "human":
            bundled = _bundled_observed_pbmc(chain)
            if bundled is not None:
                logger.info(
                    "fetch_repertoire: using bundled healthy 10x PBMC %s (%d seqs)",
                    chain, len(bundled),
                )
                bundled.to_csv(dest, index=False)
                return dest
            if chain == "beta":
                seqs = _bundled_observed_beta()
                if seqs:
                    logger.info(
                        "fetch_repertoire: using bundled healthy β (%d seqs)",
                        len(seqs),
                    )
                    pd.DataFrame({"seq": seqs}).to_csv(dest, index=False)
                    return dest
        raise FileNotFoundError(
            f"fetch_repertoire: no bundled observed {sp} {chain} repertoire — pass "
            f"url= a neutral {sp} repertoire CSV/TSV (e.g. OTS: {REPERTOIRE_SOURCES['ots']})"
        )

    url_path = Path(url)
    if url_path.is_file():
        raw = url_path
    else:
        logger.info("fetch_repertoire: downloading %s (%s)", url, source or "url")
        raw = _data_dir() / f"_download_{sp}_{chain}{_suffix(url)}"
        raw.parent.mkdir(parents=True, exist_ok=True)
        _download(url, raw)
    opener = gzip.open if str(raw).endswith(".gz") else open
    sep = "\t" if _suffix(str(raw)) in (".tsv", ".tsv.gz") else ","
    with opener(raw, "rt") as fh:
        frame = pd.read_csv(fh, sep=sep)
    if seq_col not in frame.columns:
        raise ValueError(
            f"fetch_repertoire: {seq_col!r} not in downloaded columns "
            f"{list(frame.columns)}"
        )
    if v_col and chain_prefix and v_col in frame.columns:
        frame = frame[frame[v_col].astype(str).str.upper().str.startswith(chain_prefix.upper())]
    frame = frame[frame[seq_col].notna()].copy()
    frame["seq"] = frame[seq_col].astype(str)
    frame = frame[frame["seq"].str.len() > 0]
    out_cols = {"seq": frame["seq"]}
    resolved_v = v_col if (v_col and v_col in frame.columns) else ("v" if "v" in frame.columns else ("v_gene" if "v_gene" in frame.columns else None))
    resolved_j = j_col if (j_col and j_col in frame.columns) else ("j" if "j" in frame.columns else ("j_gene" if "j_gene" in frame.columns else None))
    if resolved_v:
        out_cols["v"] = frame[resolved_v]
    if resolved_j:
        out_cols["j"] = frame[resolved_j]
    pd.DataFrame(out_cols).to_csv(dest, index=False)
    logger.info("fetch_repertoire: cached %d %s (%s) seqs → %s", len(frame), chain, sp, dest)
    return dest


def _suffix(url: str) -> str:
    import urllib.parse

    path = urllib.parse.urlparse(url).path or url
    for s in (".csv.gz", ".tsv.gz", ".csv", ".tsv", ".gz"):
        if path.endswith(s):
            return s
    return ".csv"


def _observed_repertoire_df(
    chain: str, *, url: str | None = None, species: str = "human"
) -> pd.DataFrame:
    """Observed repertoire frame (``seq`` and optional ``v``, ``j``) for Ppost training."""
    cached = cached_repertoire_file(chain, species=species)
    if cached.is_file() and url is None:
        return pd.read_csv(cached)
    path = fetch_repertoire(chain, url=url, species=species)
    return pd.read_csv(path)


def _observed_seqs(
    chain: str, *, url: str | None = None, species: str = "human"
) -> list[str]:
    """CDR3s for Ppost training: cached file → fetch (url/bundled)."""
    df = _observed_repertoire_df(chain, url=url, species=species)
    return df["seq"].dropna().astype(str).tolist()


# ---------------------------------------------------------------------------
# Train + cache
# ---------------------------------------------------------------------------

def train_model(
    chain: str,
    *,
    backend: str = DEFAULT_BACKEND,
    role: str = "ppost",
    url: str | None = None,
    order: int = 2,
    epochs: int = 30,
    batch_size: int = 2000,
    device: str | None = None,
    species: str = "human",
) -> Path:
    """Train a (backend, role, chain) model and cache it; return its path."""
    from .species import normalize_species

    sp = normalize_species(species)
    if role == "ppost":
        rep_df = _observed_repertoire_df(chain, url=url, species=sp)
        rep_df = rep_df.dropna(subset=["seq"]).reset_index(drop=True)
        seqs = rep_df["seq"].astype(str).tolist()
    elif role == "pgen":
        raise ValueError(
            "Pgen models are shipped pre-built and cannot be trained at "
            "runtime (OLGA, the synthetic-sequence generator, is no longer a "
            "dependency). Use the shipped k-mer Pgen model, or regenerate "
            "offline — see refseqs/PROVENANCE.md."
        )
    else:
        raise ValueError(f"role must be 'pgen'/'ppost', got {role!r}")

    out = cached_model_file(backend, role, chain, species=sp)
    out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("train_model: fitting %s %s %s (%s) on %d seqs → %s",
                backend, role, chain, sp, len(seqs), out)
    if backend == "kmer":
        model = seqprob.KmerProbabilityModel(order=order, chain=chain).fit(seqs)
    elif backend == "kmer_gene":
        kmer_order = min(order, 2)
        v_genes = rep_df["v"].tolist() if "v" in rep_df.columns else None
        j_genes = rep_df["j"].tolist() if "j" in rep_df.columns else None
        model = seqprob.GeneAwareKmerModel(order=kmer_order, chain=chain).fit(
            seqs, v_genes=v_genes, j_genes=j_genes
        )
    elif backend == "tcrpeg":
        dev = device or _default_device()
        model = seqprob.TCRpegProbabilityModel(
            epochs=epochs, batch_size=batch_size, device=dev, chain=chain,
        ).fit(seqs)
    else:
        raise ValueError(f"unknown backend {backend!r}")
    model.save(str(out))
    # Evict any cached model for this (chain, backend, role, species, cache_dir)
    cdir = str(cache_dir())
    for k in [k for k in _CACHE if k[:4] == (chain.lower(), backend, role, sp) and (len(k) < 5 or k[4] == cdir)]:
        _CACHE.pop(k, None)
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
    species: str = "human",
) -> seqprob.SequenceProbabilityModel:
    """Return the (backend, role, chain, species) model, training it on first use."""
    from .species import normalize_species

    sp = normalize_species(species)
    key = (chain.lower(), backend, role, sp, str(cache_dir()))
    if url is None:
        if key in _CACHE:
            return _CACHE[key]

        path = cached_model_file(backend, role, chain, species=sp)
        if path.is_file():
            model = seqprob.BACKENDS[backend].load(str(path))
            _CACHE[key] = model
            return model

    # Shipped k-mer / kmer_gene defaults (instant, no training when present).
    # For mouse with auto_train=True and an available mouse repertoire (or url),
    # or whenever the caller explicitly provides a training `url`, train instead
    # of returning a shipped default.
    can_train_mouse_ppost = (
        sp == "mouse"
        and role == "ppost"
        and auto_train
        and (url is not None or cached_repertoire_file(chain, species=sp).is_file())
    )
    if backend in ("kmer", "kmer_gene") and url is None and not can_train_mouse_ppost:
        try:
            model = seqprob.load_background_model(chain, backend, role, species=sp)
            _CACHE[key] = model
            return model
        except FileNotFoundError:
            if not auto_train:
                raise

    if not auto_train:
        raise FileNotFoundError(
            f"no cached {backend} {role} model for {sp} {chain}; run "
            f"`tcrsift pgen train` or set auto_train=True"
        )

    train_model(chain, backend=backend, role=role, url=url, species=sp)
    model = seqprob.BACKENDS[backend].load(
        str(cached_model_file(backend, role, chain, species=sp))
    )
    _CACHE[key] = model
    return model
