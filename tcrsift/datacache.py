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
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)

_META_FILENAME = "_meta.json"


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
    "cedar": DatabaseSpec(
        name="cedar",
        filename="cedar.tsv",
        source="https://cedar.iedb.org/",
        # CEDAR has no stable single-file download endpoint; users must
        # populate this DB manually.
        download_url=None,
    ),
    # IMGT/GENE-DB V-REGION germline reference (one HTML per locus). Backs the
    # framework (FR1–FR3) germline comparison in tcrsift.vregion — large and
    # rarely-changing, so cached rather than vendored. The raw GENElect HTML is
    # stored as-is; tcrsift.vregion parses the <pre> FASTA and translates it.
    "imgt_trav_vregion": DatabaseSpec(
        name="imgt_trav_vregion",
        filename="TRAV_V-REGION.html",
        source="https://www.imgt.org/genedb/ (GENElect 8.1 TRAV, IMGTlabel=V-REGION)",
        download_url=(
            "https://www.imgt.org/genedb/GENElect?query=8.1+TRAV"
            "&species=Homo+sapiens&IMGTlabel=V-REGION"
        ),
        archive_format="raw",
    ),
    "imgt_trbv_vregion": DatabaseSpec(
        name="imgt_trbv_vregion",
        filename="TRBV_V-REGION.html",
        source="https://www.imgt.org/genedb/ (GENElect 8.1 TRBV, IMGTlabel=V-REGION)",
        download_url=(
            "https://www.imgt.org/genedb/GENElect?query=8.1+TRBV"
            "&species=Homo+sapiens&IMGTlabel=V-REGION"
        ),
        archive_format="raw",
    ),
}


# Logical groupings for bulk `data download` / `data clear` by workflow, so a
# user can grab everything one pipeline needs without naming each DB.
#   germline    → the IMGT V-REGION framework reference (assembly germline QC)
#   annotation  → the epitope-annotation databases (the `annotate` command)
#   all         → every registered database
GROUPS: dict[str, list[str]] = {
    "germline": ["imgt_trav_vregion", "imgt_trbv_vregion"],
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
    databases without an automated endpoint (e.g. CEDAR) are dropped from a
    *group/default* expansion — a user asking for a whole category shouldn't fail
    on its manual-only member — while an explicit ``--db cedar`` is left in so the
    command still reports the manual instruction. Raises ``ValueError`` on an
    unknown name or group.
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


def cached_path(db: str, data_dir: str | Path | None = None) -> Path | None:
    """Return the cached path for `db` if present, else None.

    For VDJdb the annotate loader accepts either a directory containing
    `vdjdb*.txt` files or a direct file path. We return the directory if it
    holds at least one matching file, otherwise the canonical filename if it
    exists, otherwise None.
    """
    spec = DATABASES.get(db)
    if spec is None:
        return None
    base = resolve_cache_dir(data_dir) / db
    if not base.exists():
        return None
    canonical = base / spec.filename
    if canonical.exists():
        return canonical
    # VDJdb releases ship with versioned filenames; the annotate loader
    # globs `vdjdb*.{txt,tsv}` against the parent dir, so pass the
    # directory when at least one match exists.
    if db == "vdjdb":
        matches = list(base.glob("vdjdb*.txt")) + list(base.glob("vdjdb*.tsv"))
        if matches:
            return base
    return None


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
    return removed


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


class DownloadError(RuntimeError):
    """Raised when an automated database download fails."""


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
        import urllib.request
        opener = urllib.request.urlopen

    existing = cached_path(db, data_dir)
    if existing is not None and not force:
        logger.info("%s already cached at %s (use force=True to refresh)", db, existing)
        return existing

    base = resolve_cache_dir(data_dir)
    db_dir = base / db
    if db_dir.exists():
        # `force=True` path: clear before re-downloading so partial state
        # from a failed previous run can't shadow the new download.
        shutil.rmtree(db_dir)
    db_dir.mkdir(parents=True)

    url = spec.download_url
    if url == _VDJDB_LATEST_SENTINEL:
        url = _resolve_vdjdb_url(opener)

    logger.info("Downloading %s from %s", db, url)
    with tempfile.NamedTemporaryFile(
        prefix=f"{db}_", suffix=".tmp", delete=False, dir=base,
    ) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with opener(url) as resp:  # type: ignore[arg-type]
                shutil.copyfileobj(resp, tmp)
        except Exception as e:
            tmp_path.unlink(missing_ok=True)
            shutil.rmtree(db_dir, ignore_errors=True)
            raise DownloadError(f"Failed to download {db} from {url}: {e}") from e

    try:
        if spec.archive_format == "zip":
            _extract_zip(tmp_path, db_dir, spec.archive_member)
        elif spec.archive_format == "raw":
            shutil.move(str(tmp_path), str(db_dir / spec.filename))
            tmp_path = db_dir / spec.filename
        else:
            raise DownloadError(
                f"Unknown archive_format {spec.archive_format!r} for {db}"
            )
    except Exception:
        shutil.rmtree(db_dir, ignore_errors=True)
        raise
    finally:
        if tmp_path.exists() and spec.archive_format == "zip":
            tmp_path.unlink()

    resolved = cached_path(db, data_dir)
    if resolved is None:
        raise DownloadError(
            f"Downloaded {db} but canonical path is still missing — "
            f"archive_member={spec.archive_member!r} likely incorrect."
        )

    size = _path_size(resolved)
    sha256 = _sha256_of(resolved) if resolved.is_file() else ""
    _write_metadata(
        db_dir,
        url=url,
        downloaded_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        size=size,
        sha256=sha256,
    )
    return resolved


def _extract_zip(archive: Path, dest: Path, member: str | None) -> None:
    """Extract `member` (or the whole archive when member is None) into `dest`."""
    with zipfile.ZipFile(archive) as zf:
        if member is None:
            zf.extractall(dest)
            return
        # Locate the member case-insensitively against the archive's basenames
        # so a release that wraps files in a top-level dir doesn't trip us.
        candidates = [n for n in zf.namelist() if n.rsplit("/", 1)[-1] == member]
        if not candidates:
            raise DownloadError(
                f"Archive member {member!r} not found in {archive.name}. "
                f"Contents: {zf.namelist()[:10]}{'...' if len(zf.namelist()) > 10 else ''}"
            )
        with zf.open(candidates[0]) as src, (dest / member).open("wb") as dst:
            shutil.copyfileobj(src, dst)
