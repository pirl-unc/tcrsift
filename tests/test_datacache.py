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

"""Tests for tcrsift.datacache — reference-database cache management (#32)."""

import io
import json
import zipfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tcrsift.datacache import (
    DATABASES,
    GROUPS,
    DownloadError,
    cached_path,
    clear_cache,
    download_database,
    inspect_cache,
    read_metadata,
    resolve_cache_dir,
    resolve_db_selection,
)


class TestResolveDbSelection:
    """`resolve_db_selection` expands --db names + --group aliases."""

    def test_default_is_all(self):
        assert resolve_db_selection(None, None) == list(DATABASES)

    def test_germline_group(self):
        assert resolve_db_selection(None, ["germline"]) == [
            "imgt_trav_vregion", "imgt_trbv_vregion"
        ]

    def test_annotation_group_includes_cedar_now_downloadable(self):
        # CEDAR auto-downloads now, so the annotation group keeps it.
        sel = resolve_db_selection(None, ["annotation"], downloadable_only=True)
        assert "cedar" in sel and "vdjdb" in sel

    def test_downloadable_only_drops_manual_in_group(self, monkeypatch):
        # Mechanism check: a manual-only DB (no download_url) is dropped from a
        # group/default expansion under downloadable_only...
        import dataclasses

        from tcrsift import datacache as dc
        manual = dataclasses.replace(DATABASES["cedar"], download_url=None)
        monkeypatch.setitem(dc.DATABASES, "cedar", manual)
        sel = resolve_db_selection(None, ["annotation"], downloadable_only=True)
        assert "cedar" not in sel and "vdjdb" in sel
        # ...but an explicit --db <name> is kept so its manual note still surfaces.
        assert resolve_db_selection(["cedar"], None, downloadable_only=True) == ["cedar"]

    def test_db_and_group_union_dedup(self):
        sel = resolve_db_selection(["vdjdb"], ["germline"])
        assert sel == ["vdjdb", "imgt_trav_vregion", "imgt_trbv_vregion"]

    def test_unknown_group_raises(self):
        with pytest.raises(ValueError):
            resolve_db_selection(None, ["nope"])

    def test_unknown_db_raises(self):
        with pytest.raises(ValueError):
            resolve_db_selection(["nope"], None)

    def test_groups_reference_real_dbs(self):
        for members in GROUPS.values():
            assert all(m in DATABASES for m in members)


class TestResolveCacheDir:
    """Resolution precedence: explicit > env > XDG > ~/.cache/tcrsift."""

    def test_explicit_wins(self, tmp_path, monkeypatch):
        monkeypatch.setenv("TCRSIFT_DATA_DIR", "/should/not/be/used")
        monkeypatch.setenv("XDG_CACHE_HOME", "/also/ignored")
        assert resolve_cache_dir(tmp_path) == tmp_path

    def test_env_var_used_when_no_explicit(self, tmp_path, monkeypatch):
        monkeypatch.setenv("TCRSIFT_DATA_DIR", str(tmp_path))
        monkeypatch.setenv("XDG_CACHE_HOME", "/ignored")
        assert resolve_cache_dir() == tmp_path

    def test_xdg_fallback(self, tmp_path, monkeypatch):
        monkeypatch.delenv("TCRSIFT_DATA_DIR", raising=False)
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
        assert resolve_cache_dir() == tmp_path / "tcrsift"

    def test_home_cache_fallback(self, monkeypatch):
        monkeypatch.delenv("TCRSIFT_DATA_DIR", raising=False)
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        assert resolve_cache_dir() == Path.home() / ".cache" / "tcrsift"

    def test_does_not_create_directory(self, tmp_path):
        target = tmp_path / "doesnotexist"
        resolved = resolve_cache_dir(target)
        assert resolved == target
        assert not target.exists()


class TestCachedPath:
    """`cached_path` returns the resolved path or None per db."""

    def test_returns_none_when_cache_missing(self, tmp_path):
        assert cached_path("vdjdb", tmp_path) is None
        assert cached_path("iedb", tmp_path) is None

    def test_returns_canonical_file_when_present(self, tmp_path):
        (tmp_path / "iedb").mkdir()
        canonical = tmp_path / "iedb" / DATABASES["iedb"].filename
        canonical.write_text("dummy")
        assert cached_path("iedb", tmp_path) == canonical

    def test_vdjdb_returns_directory_for_versioned_files(self, tmp_path):
        """VDJdb releases ship versioned filenames; cached_path should return
        the parent dir so the annotate loader's glob picks them up."""
        vdjdb_dir = tmp_path / "vdjdb"
        vdjdb_dir.mkdir()
        (vdjdb_dir / "vdjdb-2024-08-15.txt").write_text("dummy")
        assert cached_path("vdjdb", tmp_path) == vdjdb_dir

    def test_returns_none_for_unknown_db(self, tmp_path):
        assert cached_path("unknown_db", tmp_path) is None


class TestInspectCache:
    """`inspect_cache` returns one entry per known DB, present or not."""

    def test_all_missing_when_empty(self, tmp_path):
        entries = inspect_cache(tmp_path)
        names = {e.name for e in entries}
        assert names == set(DATABASES)
        assert all(not e.present for e in entries)
        assert all(e.size_bytes == 0 for e in entries)

    def test_present_entry_carries_size(self, tmp_path):
        (tmp_path / "iedb").mkdir()
        f = tmp_path / "iedb" / DATABASES["iedb"].filename
        f.write_bytes(b"x" * 1024)
        entries = inspect_cache(tmp_path)
        iedb = next(e for e in entries if e.name == "iedb")
        assert iedb.present
        assert iedb.size_bytes == 1024

    def test_mixed_state(self, tmp_path):
        (tmp_path / "vdjdb").mkdir()
        (tmp_path / "vdjdb" / "vdjdb-2024-01-01.txt").write_text("x")
        entries = {e.name: e for e in inspect_cache(tmp_path)}
        assert entries["vdjdb"].present
        assert not entries["iedb"].present
        assert not entries["cedar"].present


class TestClearCache:
    """`clear_cache` is idempotent and reports what it removed."""

    def test_clear_specific_db(self, tmp_path):
        (tmp_path / "iedb").mkdir()
        (tmp_path / "iedb" / "x.tsv").write_text("dummy")
        removed = clear_cache("iedb", data_dir=tmp_path)
        assert removed == [tmp_path / "iedb"]
        assert not (tmp_path / "iedb").exists()

    def test_clear_all_when_db_none(self, tmp_path):
        for db in DATABASES:
            (tmp_path / db).mkdir()
            (tmp_path / db / "dummy.txt").write_text("x")
        removed = clear_cache(data_dir=tmp_path)
        assert {p.name for p in removed} == set(DATABASES)
        for db in DATABASES:
            assert not (tmp_path / db).exists()

    def test_clear_is_idempotent(self, tmp_path):
        # No directory exists — should not raise.
        assert clear_cache(data_dir=tmp_path) == []
        assert clear_cache("iedb", data_dir=tmp_path) == []

    def test_unknown_db_raises(self, tmp_path):
        with pytest.raises(ValueError, match="Unknown database"):
            clear_cache("not_a_db", data_dir=tmp_path)


# ---------------------------------------------------------------------------
# Downloads — all HTTP I/O is injected via `opener=` for hermetic tests.
# ---------------------------------------------------------------------------


def _mock_opener(url_to_bytes: dict[str, bytes]):
    """Build a fake urlopen that returns the bytes for each known URL.

    Returned object supports the context-manager protocol used by callers
    (`with opener(url) as resp: ...`).
    """
    def opener(url):
        if url not in url_to_bytes:
            raise AssertionError(f"Unexpected URL in test: {url}")
        return _FakeResponse(url_to_bytes[url])
    return opener


class _FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _make_iedb_zip(payload: bytes = b"cdr3,epitope\nCASS,FLY\n") -> bytes:
    """Build an in-memory zip mirroring IEDB's release shape."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(DATABASES["iedb"].archive_member, payload)
    return buf.getvalue()


def _make_vdjdb_zip(filename: str = "vdjdb-2024-01-01.txt") -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr(filename, b"vdjdb data\n")
    return buf.getvalue()


class TestDownloadIEDB:
    """IEDB uses a stable download URL; the simpler test case."""

    def test_downloads_and_extracts(self, tmp_path):
        zip_bytes = _make_iedb_zip()
        opener = _mock_opener({DATABASES["iedb"].download_url: zip_bytes})

        result = download_database("iedb", data_dir=tmp_path, opener=opener)

        canonical = tmp_path / "iedb" / DATABASES["iedb"].filename
        assert result == canonical
        assert canonical.exists()
        assert canonical.read_bytes().startswith(b"cdr3,epitope")

    def test_writes_metadata_sidecar(self, tmp_path):
        zip_bytes = _make_iedb_zip()
        opener = _mock_opener({DATABASES["iedb"].download_url: zip_bytes})

        download_database("iedb", data_dir=tmp_path, opener=opener)

        meta = read_metadata("iedb", data_dir=tmp_path)
        assert meta is not None
        assert meta["source_url"] == DATABASES["iedb"].download_url
        assert meta["size_bytes"] > 0
        assert len(meta["sha256"]) == 64  # hex sha256
        assert "downloaded_at" in meta

    def test_skips_when_already_cached(self, tmp_path):
        opener = MagicMock(side_effect=AssertionError("should not be called"))
        # Seed the cache.
        (tmp_path / "iedb").mkdir()
        (tmp_path / "iedb" / DATABASES["iedb"].filename).write_text("seed")

        result = download_database("iedb", data_dir=tmp_path, opener=opener)

        assert result == tmp_path / "iedb" / DATABASES["iedb"].filename
        assert result.read_text() == "seed"
        opener.assert_not_called()

    def test_force_redownloads(self, tmp_path):
        # Pre-seed with stale content, then refresh with new content.
        (tmp_path / "iedb").mkdir()
        canonical = tmp_path / "iedb" / DATABASES["iedb"].filename
        canonical.write_text("stale")

        fresh = _make_iedb_zip(payload=b"fresh content\n")
        opener = _mock_opener({DATABASES["iedb"].download_url: fresh})

        result = download_database(
            "iedb", data_dir=tmp_path, force=True, opener=opener
        )

        assert result.read_bytes() == b"fresh content\n"

    def test_missing_archive_member_raises(self, tmp_path):
        bogus = io.BytesIO()
        with zipfile.ZipFile(bogus, "w") as zf:
            zf.writestr("wrong_name.csv", b"x")
        opener = _mock_opener({DATABASES["iedb"].download_url: bogus.getvalue()})

        with pytest.raises(DownloadError, match="not found"):
            download_database("iedb", data_dir=tmp_path, opener=opener)
        # Failed download shouldn't leave a partial cache dir.
        assert not (tmp_path / "iedb").exists()


class TestDownloadVDJdb:
    """VDJdb resolves its URL dynamically via the GitHub releases API."""

    def test_resolves_via_github_api(self, tmp_path):
        api_payload = json.dumps({
            "tag_name": "2024-08-15",
            "assets": [
                {"name": "vdjdb-2024-08-15.zip",
                 "browser_download_url": "https://example.com/vdjdb-2024-08-15.zip"},
                {"name": "Source code (zip)",
                 "browser_download_url": "https://example.com/source.zip"},
            ],
        }).encode()
        zip_bytes = _make_vdjdb_zip("vdjdb-2024-08-15.txt")
        opener = _mock_opener({
            "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest": api_payload,
            "https://example.com/vdjdb-2024-08-15.zip": zip_bytes,
        })

        result = download_database("vdjdb", data_dir=tmp_path, opener=opener)

        # cached_path returns the directory for VDJdb when versioned files
        # are present (the annotate loader globs).
        assert result == tmp_path / "vdjdb"
        assert (tmp_path / "vdjdb" / "vdjdb-2024-08-15.txt").exists()

    def test_missing_release_asset_raises(self, tmp_path):
        api_payload = json.dumps({
            "tag_name": "2024-08-15",
            "assets": [{"name": "irrelevant.txt", "browser_download_url": "x"}],
        }).encode()
        opener = _mock_opener({
            "https://api.github.com/repos/antigenomics/vdjdb-db/releases/latest": api_payload,
        })

        with pytest.raises(DownloadError, match="No vdjdb-.*zip asset"):
            download_database("vdjdb", data_dir=tmp_path, opener=opener)


class TestDownloadErrors:
    def test_manual_only_db_has_no_download(self, tmp_path, monkeypatch):
        # A DB without a download_url raises with the manual instruction. (CEDAR
        # is now auto-downloadable, so synthesize a manual spec to exercise it.)
        import dataclasses

        from tcrsift import datacache as dc
        manual = dataclasses.replace(DATABASES["cedar"], download_url=None)
        monkeypatch.setitem(dc.DATABASES, "cedar", manual)
        with pytest.raises(DownloadError, match="no automated download"):
            download_database("cedar", data_dir=tmp_path)

    def test_unknown_db(self, tmp_path):
        with pytest.raises(ValueError, match="Unknown database"):
            download_database("nope", data_dir=tmp_path)
