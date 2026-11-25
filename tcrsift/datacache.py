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

"""
Cache management for IEDB / VDJdb / CEDAR public databases.

The cache holds copies of TCR-annotation reference databases so users don't
have to pass `--vdjdb`/`--iedb`/`--cedar` paths on every invocation, and so
fresh installs have a clear discoverable location to drop the files.

Resolution order for the cache directory (first that's set wins):
  1. explicit `data_dir=` passed to the API
  2. `--cache-dir` CLI flag
  3. `TCRSIFT_DATA_DIR` environment variable
  4. `$XDG_CACHE_HOME/tcrsift`
  5. `~/.cache/tcrsift`

Within the cache, each database lives at `<cache_dir>/<db>/` with the
canonical filename declared in `DATABASES`, plus a `_meta.json` sidecar
recording the source URL, download timestamp, sha256, and file size for
every automated download. The annotate path consults the cache when no
explicit `--vdjdb` / `--iedb` / `--cedar` is provided.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import logging
import os
import shutil
import tempfile
import threading
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)

_META_FILENAME = "_meta.json"
_IPV4_FALLBACK_LOCK = threading.Lock()


@dataclass(frozen=True)
class DatabaseSpec:
    """Descriptor for a managed reference database.

    Attributes
    ----------
    name : str
        Cache key; also the per-db subdirectory name under the cache root.
    filename : str
        Canonical filename the annotate loader expects. For VDJdb the loader
        globs `vdjdb*.{txt,tsv}` against the parent dir, so this is just a
        sentinel — `cached_path` returns the directory when versioned files
        are present.
    source : str
        Human-readable URL surfaced in `list` output and warnings.
    download_url : str | None
        Stable download URL, or None when no automated download is
        supported. VDJdb's URL is resolved dynamically via GitHub's
        latest-release API (see `_resolve_vdjdb_url`); pass the sentinel
        ``"<vdjdb-latest>"`` to trigger that branch.
    archive_format : str
        ``"raw"`` (no extraction) or ``"zip"``.
    archive_member : str | None
        When `archive_format == "zip"`, the path inside the archive that
        should be extracted as the canonical filename. None means extract
        the whole archive in place.
    """

    name: str
    filename: str
    source: str
    download_url: str | None = None
    archive_format: str = "raw"
    archive_member: str | None = None


# Sentinel value: resolve via GitHub releases API at download time.
_VDJDB_LATEST_SENTINEL = "<vdjdb-latest>"


# Order matters only for stable `list` output.
DATABASES: dict[str, DatabaseSpec] = {
    "vdjdb": DatabaseSpec(
        name="vdjdb",
        filename="vdjdb.txt",
        source="https://vdjdb.cdr3.net/",
        download_url=_VDJDB_LATEST_SENTINEL,
        archive_format="zip",
        # The release zip contains multiple files; the annotate loader
        # accepts a directory and globs vdjdb*.{txt,tsv}, so we let
        # extraction drop every file into the per-db dir.
        archive_member=None,
    ),
    "iedb": DatabaseSpec(
        name="iedb",
        filename="tcr_full_v3.csv",
        source="https://www.iedb.org/database_export_v3.php",
        download_url=(
            "https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
        ),
        archive_format="zip",
        archive_member="tcr_full_v3.csv",
    ),
    # Companion to ``iedb``: the per-epitope authoritative table.
    # ``load_iedb`` uses it to override Source Molecule / Source Organism
    # names with the shorter, more canonical strings the epitope table
    # ships (e.g. ``Protein Tax-1`` instead of the receptor file's
    # ``transcriptional activator Tax``). Optional — when absent,
    # ``load_iedb`` just returns the receptor values unchanged.
    "iedb_epitope": DatabaseSpec(
        name="iedb_epitope",
        filename="epitope_full_v3.csv",
        source="https://www.iedb.org/database_export_v3.php",
        download_url=(
            "https://www.iedb.org/downloader.php?file_name=doc/epitope_full_v3.zip"
        ),
        archive_format="zip",
        archive_member="epitope_full_v3.csv",
    ),
    # CEDAR (cancer-scoped IEDB receptors). Its bulk export mirrors IEDB's:
    # ``doc/receptor_full_v3.zip`` contains ``tcr_full_v3.csv`` in the same v3
    # two-row-header schema, so it auto-downloads and ``load_cedar`` reuses the
    # IEDB v3 parser. (A legacy manual web-export ``cedar.tsv`` is still honored
    # — see ``cached_path`` and ``load_cedar``'s format sniff.)
    "cedar": DatabaseSpec(
        name="cedar",
        filename="tcr_full_v3.csv",
        source="https://cedar.iedb.org/",
        download_url=(
            "https://cedar.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip"
        ),
        archive_format="zip",
        archive_member="tcr_full_v3.csv",
    ),
    # IMGT V-QUEST reference directory FASTAs (one per locus and species). Backs
    # framework (FR1–FR3) and J-region (FR4) germline comparison and 5'-partial
    # contig rescue in tcrsift.assemble.vregion and tcrsift.assemble.germline. Cached
    # rather than vendored; tcrsift.assemble.vregion parses both raw gapped FASTA and
    # legacy GENElect <pre> HTML.
    "imgt_trav_vregion": DatabaseSpec(
        name="imgt_trav_vregion",
        filename="TRAV_V-REGION.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAV.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Homo_sapiens/TR/TRAV.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_trbv_vregion": DatabaseSpec(
        name="imgt_trbv_vregion",
        filename="TRBV_V-REGION.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBV.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Homo_sapiens/TR/TRBV.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_traj_jregion": DatabaseSpec(
        name="imgt_traj_jregion",
        filename="TRAJ_J-REGION.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAJ.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Homo_sapiens/TR/TRAJ.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_trbj_jregion": DatabaseSpec(
        name="imgt_trbj_jregion",
        filename="TRBJ_J-REGION.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBJ.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Homo_sapiens/TR/TRBJ.fasta"
        ),
        archive_format="raw",
    ),
    # Mouse (Mus musculus) IMGT V-REGION, J-REGION, and L-PART1+L-PART2
    # references for murine framework/J germline comparison and 5'-partial rescue.
    "imgt_trav_vregion_mouse": DatabaseSpec(
        name="imgt_trav_vregion_mouse",
        filename="TRAV_V-REGION_mouse.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRAV.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Mus_musculus/TR/TRAV.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_trbv_vregion_mouse": DatabaseSpec(
        name="imgt_trbv_vregion_mouse",
        filename="TRBV_V-REGION_mouse.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRBV.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Mus_musculus/TR/TRBV.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_traj_jregion_mouse": DatabaseSpec(
        name="imgt_traj_jregion_mouse",
        filename="TRAJ_J-REGION_mouse.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRAJ.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Mus_musculus/TR/TRAJ.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_trbj_jregion_mouse": DatabaseSpec(
        name="imgt_trbj_jregion_mouse",
        filename="TRBJ_J-REGION_mouse.html",
        source="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRBJ.fasta",
        download_url=(
            "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/"
            "Mus_musculus/TR/TRBJ.fasta"
        ),
        archive_format="raw",
    ),
    "imgt_trav_leader_mouse": DatabaseSpec(
        name="imgt_trav_leader_mouse",
        filename="TRAV_LEADER_mouse.html",
        source="https://www.imgt.org/genedb/ (GENElect 8.1 TRAV, species=Mus musculus, IMGTlabel=L-PART1+L-PART2)",
        download_url=(
            "https://www.imgt.org/genedb/GENElect?query=8.1+TRAV"
            "&species=Mus+musculus&IMGTlabel=L-PART1%2BL-PART2"
        ),
        archive_format="raw",
    ),
    "imgt_trbv_leader_mouse": DatabaseSpec(
        name="imgt_trbv_leader_mouse",
        filename="TRBV_LEADER_mouse.html",
        source="https://www.imgt.org/genedb/ (GENElect 8.1 TRBV, species=Mus musculus, IMGTlabel=L-PART1+L-PART2)",
        download_url=(
            "https://www.imgt.org/genedb/GENElect?query=8.1+TRBV"
            "&species=Mus+musculus&IMGTlabel=L-PART1%2BL-PART2"
        ),
        archive_format="raw",
    ),
}


# Logical groupings for bulk `data download` / `data clear` by workflow, so a
# user can grab everything one pipeline needs without naming each DB.
#   germline       → human IMGT V-REGION + J-REGION references (assembly germline QC & rescue)
#   germline_mouse → mouse IMGT V-REGION + J-REGION + leader references
#   annotation     → the epitope-annotation databases (the `annotate` command)
#   all            → every registered database
GROUPS: dict[str, list[str]] = {
    "germline": [
        "imgt_trav_vregion",
        "imgt_trbv_vregion",
        "imgt_traj_jregion",
        "imgt_trbj_jregion",
    ],
    "germline_mouse": [
        "imgt_trav_vregion_mouse",
        "imgt_trbv_vregion_mouse",
        "imgt_traj_jregion_mouse",
        "imgt_trbj_jregion_mouse",
        "imgt_trav_leader_mouse",
        "imgt_trbv_leader_mouse",
    ],
    "annotation": ["vdjdb", "iedb", "iedb_epitope", "cedar"],
    "all": list(DATABASES),
}


def resolve_db_selection(
    dbs: list[str] | None,
    groups: list[str] | None,
    *,
    downloadable_only: bool = False,
) -> list[str]:
    """Expand ``--db`` names + ``--group`` aliases into a de-duplicated db list.

    When both are empty the selection defaults to every database (the historical
    no-argument behaviour). With ``downloadable_only=True`` (the download path),
    any database without an automated endpoint is dropped from a *group/default*
    expansion — a user asking for a whole category shouldn't fail on a manual-only
    member — while an explicit ``--db <name>`` is left in so the command still
    reports that DB's manual instruction. (Every DB is auto-downloadable today,
    so this only matters if a manual-only DB is added later.) Raises
    ``ValueError`` on an unknown name or group.
    """
    selected: list[str] = []
    explicit: set[str] = set()
    if dbs:
        for name in dbs:
            if name not in DATABASES:
                raise ValueError(f"Unknown database {name!r}; known: {sorted(DATABASES)}")
            selected.append(name)
            explicit.add(name)
    if groups:
        for g in groups:
            if g not in GROUPS:
                raise ValueError(f"Unknown group {g!r}; known: {sorted(GROUPS)}")
            selected.extend(GROUPS[g])
    if not dbs and not groups:
        selected = list(DATABASES)
    # de-dupe, preserve order
    seen: set[str] = set()
    ordered = [n for n in selected if not (n in seen or seen.add(n))]
    if downloadable_only:
        ordered = [
            n for n in ordered
            if DATABASES[n].download_url is not None or n in explicit
        ]
    return ordered


def resolve_cache_dir(data_dir: str | Path | None = None) -> Path:
    """Resolve the cache directory using the documented precedence chain.

    Does not create the directory; callers that write should call
    `path.mkdir(parents=True, exist_ok=True)` themselves so read-only
    operations (`list`) don't side-effect the filesystem.
    """
    if data_dir is not None:
        return Path(data_dir).expanduser()
    env_dir = os.environ.get("TCRSIFT_DATA_DIR")
    if env_dir:
        return Path(env_dir).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg).expanduser() / "tcrsift"
    return Path.home() / ".cache" / "tcrsift"


def _dir_has_vdjdb_files(directory: Path) -> bool:
    """True when ``directory`` directly contains legacy TSV or 2026+ relational VDJdb files."""
    if not directory.is_dir():
        return False
    sidecars = (".meta.txt", ".scored.txt", "_broken.txt", ".meta.tsv", ".scored.tsv")
    canonical_names = ("vdjdb_full.txt", "vdjdb.txt", "vdjdb.slim.txt", "vdjdb.parquet")
    if any((directory / n).is_file() for n in canonical_names):
        return True
    if (directory / "chains.tsv.gz").is_file() and (directory / "records.tsv.gz").is_file():
        return True
    if (directory / "chains.tsv").is_file() and (directory / "records.tsv").is_file():
        return True
    return any(
        p.is_file() and not p.name.lower().endswith(sidecars)
        for p in (*directory.glob("vdjdb-*.txt"), *directory.glob("vdjdb-*.tsv"))
    )


def cached_path(db: str, data_dir: str | Path | None = None) -> Path | None:
    """Return the cached path for `db` if present, else None.

    For VDJdb the annotate loader accepts either a directory containing
    legacy `vdjdb*.txt` files, a 2026+ relational bundle (`chains.tsv.gz` +
    `records.tsv.gz` / `vdjdb.parquet`), or a versioned `vdjdb-YYYY-MM-DD/`
    subdirectory extracted from the GitHub release zip.
    """
    spec = DATABASES.get(db)
    if spec is None:
        return None
    base = resolve_cache_dir(data_dir) / db
    if not base.exists():
        return None
    if db == "vdjdb" and base.is_dir():
        if _dir_has_vdjdb_files(base) or any(
            _dir_has_vdjdb_files(sub) for sub in base.glob("vdjdb-*") if sub.is_dir()
        ):
            return base
        return None
    canonical = base / spec.filename
    if canonical.exists():
        return canonical
    # CEDAR back-compat: the canonical is now the bulk `tcr_full_v3.csv`, but a
    # user who hand-exported the older `cedar.tsv` (or `.csv`) from the web UI
    # should still be found. `load_cedar` sniffs the format either way.
    if db == "cedar":
        for legacy in ("cedar.tsv", "cedar.csv"):
            if (base / legacy).exists():
                return base / legacy
    return None


def ensure_cached(
    db: str,
    data_dir: str | Path | None = None,
    *,
    auto_download: bool = False,
    opener: Callable[[str], object] | None = None,
) -> Path | None:
    """Return the cached path for ``db``, optionally downloading if absent.

    When ``auto_download=False`` (default), equivalent to :func:`cached_path`.
    When ``auto_download=True`` and ``db`` is missing from the cache, fetches
    it via :func:`download_database` and returns the populated path.
    """
    existing = cached_path(db, data_dir)
    if existing is not None:
        return existing
    if not auto_download:
        return None
    return download_database(db, data_dir=data_dir, opener=opener)


@dataclass(frozen=True)
class CacheEntry:
    """State of a single managed database in the cache."""

    name: str
    present: bool
    path: Path  # canonical path (may not exist when present=False)
    size_bytes: int  # 0 when absent or directory entry
    source: str


def inspect_cache(data_dir: str | Path | None = None) -> list[CacheEntry]:
    """Inspect the cache and return one CacheEntry per known database.

    Always returns every entry in `DATABASES` so `list` output is uniform
    regardless of which DBs the user has populated.
    """
    base = resolve_cache_dir(data_dir)
    entries: list[CacheEntry] = []
    for db, spec in DATABASES.items():
        resolved = cached_path(db, data_dir)
        canonical = base / db / spec.filename
        if resolved is None:
            entries.append(
                CacheEntry(
                    name=db, present=False, path=canonical,
                    size_bytes=0, source=spec.source,
                )
            )
            continue
        size = _path_size(resolved)
        entries.append(
            CacheEntry(
                name=db, present=True, path=resolved,
                size_bytes=size, source=spec.source,
            )
        )
    return entries


def _path_size(path: Path) -> int:
    """Total bytes at `path` — single file or sum-of-files when directory."""
    if path.is_file():
        return path.stat().st_size
    if path.is_dir():
        return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
    return 0


def _invalidate_germline_memo_caches() -> None:
    """Reset in-memory leader/V-REGION caches after germline cache updates."""
    try:
        from .assemble import leaders, vregion

        leaders.clear_reference_cache()
        vregion.clear_reference_cache()
    except Exception:
        pass


def clear_cache(
    db: str | None = None, data_dir: str | Path | None = None,
) -> list[Path]:
    """Remove cached files for `db` (or all known DBs when None).

    Returns the list of paths that were removed. Missing DBs and missing
    cache directories are silently treated as no-ops so the command is
    idempotent.
    """
    base = resolve_cache_dir(data_dir)
    targets = [db] if db else list(DATABASES)
    removed: list[Path] = []
    for name in targets:
        if name not in DATABASES:
            raise ValueError(
                f"Unknown database {name!r}; known: {sorted(DATABASES)}"
            )
        db_dir = base / name
        if db_dir.exists():
            shutil.rmtree(db_dir)
            removed.append(db_dir)
    if removed:
        _invalidate_germline_memo_caches()
    return removed


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


class DownloadError(RuntimeError):
    """Raised when an automated database download fails."""


def _default_urlopen(url: str, *, timeout: float = 20.0):
    """Open ``url`` with a User-Agent, cookie jar, and IPv4 fallback on hang/403.

    Several upstream scientific mirrors (including ``www.imgt.org`` and
    ``ftp.ncbi.nlm.nih.gov``) advertise IPv6 ``AAAA`` records that time out or
    return HTTP 403 on dual-stack hosts while succeeding immediately over IPv4,
    and IMGT ``GENElect`` leader queries redirect via a ``JSESSIONID`` cookie.
    """
    import socket
    import urllib.error
    import urllib.request

    from .version import __version__

    req = urllib.request.Request(url, headers={"User-Agent": f"tcrsift/{__version__}"})

    def _open_once():
        # Respect monkeypatched urllib.request.urlopen in tests; otherwise attach
        # an ephemeral cookie jar so Struts/Tomcat 302 redirects succeed.
        if getattr(urllib.request.urlopen, "__module__", "") == "urllib.request":
            opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor())
            return opener.open(req, timeout=timeout)
        return urllib.request.urlopen(req, timeout=timeout)

    try:
        return _open_once()
    except (urllib.error.URLError, TimeoutError, OSError) as first_err:
        if isinstance(first_err, urllib.error.HTTPError) and first_err.code not in (403, 502, 503, 504):
            raise
        with _IPV4_FALLBACK_LOCK:
            orig_getaddrinfo = socket.getaddrinfo

            def _ipv4_only(*args, **kwargs):
                res = orig_getaddrinfo(*args, **kwargs)
                v4 = [r for r in res if r[0] == socket.AF_INET]
                return v4 or res

            socket.getaddrinfo = _ipv4_only
            try:
                return _open_once()
            finally:
                socket.getaddrinfo = orig_getaddrinfo


def _resolve_vdjdb_url(opener: Callable[[str], object]) -> str:
    """Resolve the latest VDJdb release asset URL via the GitHub releases API.

    `opener` is injected for testability — defaults to urllib.request.urlopen
    in `download_database`.
    """
    api = "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest"
    try:
        with opener(api) as resp:  # type: ignore[arg-type]
            payload = json.load(resp)
    except Exception as e:
        raise DownloadError(
            f"Could not query GitHub for latest VDJdb release: {e}"
        ) from e
    for asset in payload.get("assets", []):
        name = asset.get("name", "")
        if name.startswith("vdjdb-") and name.endswith(".zip"):
            return asset["browser_download_url"]
    raise DownloadError(
        "No vdjdb-*.zip asset found in latest GitHub release "
        f"(tag={payload.get('tag_name')!r}). Either the release format "
        "changed or the API response was unexpected."
    )


def _sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_metadata(db_dir: Path, url: str, downloaded_at: str, size: int, sha256: str) -> None:
    meta = {
        "source_url": url,
        "downloaded_at": downloaded_at,
        "size_bytes": size,
        "sha256": sha256,
    }
    (db_dir / _META_FILENAME).write_text(json.dumps(meta, indent=2) + "\n")


def read_metadata(db: str, data_dir: str | Path | None = None) -> dict | None:
    """Return the metadata sidecar for `db` if present, else None."""
    base = resolve_cache_dir(data_dir) / db / _META_FILENAME
    if not base.exists():
        return None
    return json.loads(base.read_text())


def download_database(
    db: str,
    data_dir: str | Path | None = None,
    force: bool = False,
    opener: Callable[[str], object] | None = None,
) -> Path:
    """Download `db` into the cache and return its canonical path.

    Skips the download when the canonical path already exists unless
    `force=True`. Raises `DownloadError` for missing/unsupported DBs and
    network/extraction failures.

    `opener` is injected for testability — defaults to
    urllib.request.urlopen.
    """
    spec = DATABASES.get(db)
    if spec is None:
        raise ValueError(f"Unknown database {db!r}; known: {sorted(DATABASES)}")
    if spec.download_url is None:
        raise DownloadError(
            f"{db!r} has no automated download endpoint. "
            f"Download manually from {spec.source} and place files under the cache."
        )

    if opener is None:
        opener = _default_urlopen

    existing = cached_path(db, data_dir)
    if existing is not None and not force:
        logger.info("%s already cached at %s (use force=True to refresh)", db, existing)
        return existing

    base = resolve_cache_dir(data_dir)
    base.mkdir(parents=True, exist_ok=True)
    db_dir = base / db

    url = spec.download_url
    if url == _VDJDB_LATEST_SENTINEL:
        url = _resolve_vdjdb_url(opener)

    # Stage download and extraction into a sibling temp dir so a failed
    # `force=True` refresh never destroys an existing working cache.
    staging_dir = Path(tempfile.mkdtemp(prefix=f".{db}_staging_", dir=base))
    tmp_path: Path | None = None
    try:
        logger.info("Downloading %s from %s", db, url)
        with tempfile.NamedTemporaryFile(
            prefix=f"{db}_", suffix=".tmp", delete=False, dir=base,
        ) as tmp:
            tmp_path = Path(tmp.name)
            try:
                with opener(url) as resp:  # type: ignore[arg-type]
                    shutil.copyfileobj(resp, tmp)
            except Exception as e:
                raise DownloadError(f"Failed to download {db} from {url}: {e}") from e

        if spec.archive_format == "zip":
            _extract_zip(tmp_path, staging_dir, spec.archive_member)
        elif spec.archive_format == "raw":
            shutil.move(str(tmp_path), str(staging_dir / spec.filename))
            tmp_path = None
        else:
            raise DownloadError(
                f"Unknown archive_format {spec.archive_format!r} for {db}"
            )

        if db == "vdjdb":
            has_vdjdb = _dir_has_vdjdb_files(staging_dir) or any(
                _dir_has_vdjdb_files(sub)
                for sub in staging_dir.glob("vdjdb-*")
                if sub.is_dir()
            )
            if not has_vdjdb:
                raise DownloadError(
                    "Downloaded vdjdb but no VDJdb files found after extraction."
                )
            staged_target = staging_dir
        else:
            staged_target = staging_dir / spec.filename
            if not staged_target.exists():
                raise DownloadError(
                    f"Downloaded {db} but canonical path is still missing — "
                    f"archive_member={spec.archive_member!r} likely incorrect."
                )

        size = _path_size(staged_target)
        sha256 = _sha256_of(staged_target) if staged_target.is_file() else ""
        _write_metadata(
            staging_dir,
            url=url,
            downloaded_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            size=size,
            sha256=sha256,
        )

        # Atomically swap staging_dir into db_dir (keeping a backup until move succeeds).
        backup_dir: Path | None = None
        if db_dir.exists():
            backup_dir = Path(tempfile.mkdtemp(prefix=f".{db}_old_", dir=base))
            backup_dir.rmdir()
            shutil.move(str(db_dir), str(backup_dir))
        try:
            shutil.move(str(staging_dir), str(db_dir))
        except Exception:
            if backup_dir is not None and backup_dir.exists() and not db_dir.exists():
                shutil.move(str(backup_dir), str(db_dir))
            raise
        if backup_dir is not None and backup_dir.exists():
            shutil.rmtree(backup_dir, ignore_errors=True)
    except Exception:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink(missing_ok=True)

    resolved = cached_path(db, data_dir)
    if resolved is None:
        raise DownloadError(
            f"Downloaded {db} but canonical path is still missing — "
            f"archive_member={spec.archive_member!r} likely incorrect."
        )
    _invalidate_germline_memo_caches()
    return resolved


def _extract_zip(archive: Path, dest: Path, member: str | None) -> None:
    """Extract `member` (or the whole archive when member is None) into `dest`.

    Validates every archive entry path against ``dest.resolve()`` to block
    Zip Slip path-traversal entries.
    """
    dest_resolved = dest.resolve()
    try:
        zf_ctx = zipfile.ZipFile(archive)
    except zipfile.BadZipFile as exc:
        raise DownloadError(
            f"Downloaded archive {archive} is not a valid zip file: {exc}"
        ) from exc
    with zf_ctx as zf:
        if member is None:
            for info in zf.infolist():
                target = (dest_resolved / info.filename).resolve()
                try:
                    target.relative_to(dest_resolved)
                except ValueError as exc:
                    raise DownloadError(
                        f"Refusing to extract unsafe archive member {info.filename!r} "
                        f"outside destination {dest_resolved}"
                    ) from exc
            zf.extractall(dest_resolved)
            return
        # Locate the member case-insensitively against the archive's basenames
        # so a release that wraps files in a top-level dir doesn't trip us.
        member_lower = member.lower()
        candidates = [
            n for n in zf.namelist() if n.rsplit("/", 1)[-1].lower() == member_lower
        ]
        if not candidates:
            raise DownloadError(
                f"Archive member {member!r} not found in {archive.name}. "
                f"Contents: {zf.namelist()[:10]}{'...' if len(zf.namelist()) > 10 else ''}"
            )
        target = (dest_resolved / Path(member).name).resolve()
        try:
            target.relative_to(dest_resolved)
        except ValueError as exc:
            raise DownloadError(f"Unsafe target member name {member!r}") from exc
        with zf.open(candidates[0]) as src, target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
