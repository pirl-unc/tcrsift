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

"""Reference-database cache CLI subcommands: ``data list|clear|download``."""

from __future__ import annotations

from .common import _format_size, setup_logging


def cmd_data_list(args):
    """List managed reference databases and their cache state."""
    from ..datacache import DATABASES, inspect_cache, resolve_cache_dir

    cache_dir = resolve_cache_dir(args.cache_dir)
    entries = inspect_cache(args.cache_dir)

    print(f"Cache directory: {cache_dir}")
    print()
    print(f"{'database':<10} {'status':<8} {'size':>10}  path")
    print("-" * 80)
    total = 0
    for e in entries:
        status = "present" if e.present else "missing"
        size = _format_size(e.size_bytes)
        print(f"{e.name:<10} {status:<8} {size:>10}  {e.path}")
        total += e.size_bytes
    print("-" * 80)
    print(f"{'total':<10} {'':<8} {_format_size(total):>10}")

    missing = [e for e in entries if not e.present]
    auto = [e for e in missing if DATABASES[e.name].download_url is not None]
    manual = [e for e in missing if DATABASES[e.name].download_url is None]
    if auto:
        print()
        names = " ".join(e.name for e in auto)
        print(f"Missing — fetch with: tcrsift data download --db {names}")
    if manual:
        print()
        print("Missing (manual only) — download and place under the cache:")
        for e in manual:
            print(f"  {e.name}: {e.source} → {e.path}")


def cmd_data_clear(args):
    """Remove cached files for one or more managed databases / groups."""
    from ..datacache import clear_cache, resolve_db_selection

    # Expand --db names + --group aliases (default = all). Order preserved so
    # the user's first mention determines report order.
    targets = resolve_db_selection(args.db, getattr(args, "group", None))
    removed_total: list = []
    for name in targets:
        removed_total.extend(clear_cache(db=name, data_dir=args.cache_dir))

    if not removed_total:
        print("Nothing to remove.")
        return
    for path in removed_total:
        print(f"Removed {path}")


def cmd_data_download(args):
    """Download one or more supported reference databases / groups into the cache."""
    from ..datacache import DownloadError, download_database, resolve_db_selection

    setup_logging(bool(getattr(args, "verbose", False)))
    # Expand --db names + --group aliases. Default (neither) = all auto-
    # downloadable DBs; group/default expansions skip manual-only DBs (CEDAR),
    # while an explicit `--db cedar` is kept so its manual instruction surfaces.
    targets = resolve_db_selection(
        args.db, getattr(args, "group", None), downloadable_only=True
    )
    failures: list[tuple[str, str]] = []
    for db in targets:
        try:
            path = download_database(db, data_dir=args.cache_dir, force=args.force)
            print(f"  {db}: cached at {path}")
        except DownloadError as e:
            failures.append((db, str(e)))
            print(f"  {db}: FAILED — {e}")
    if failures:
        raise SystemExit(1)


def register(subparsers) -> None:
    """Register the ``data`` command group (``list``, ``clear``, ``download``)."""
    from ..datacache import DATABASES, GROUPS

    p_data = subparsers.add_parser(
        "data",
        help="Manage cached reference databases (VDJdb, IEDB, CEDAR, IMGT V/J-REGION)",
    )
    p_data_sub = p_data.add_subparsers(dest="subcommand", required=True)

    p_data_list = p_data_sub.add_parser(
        "list", help="Show cached database status and disk usage"
    )
    p_data_list.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_list.set_defaults(func=cmd_data_list)

    p_data_clear = p_data_sub.add_parser(
        "clear", help="Remove cached files for one or all databases"
    )
    p_data_clear.add_argument(
        "--db",
        choices=sorted(DATABASES),
        action="extend",
        nargs="+",
        help=(
            "Database(s) to clear (default: clear all). "
            "Repeatable and/or space-separated: "
            "`--db vdjdb iedb` and `--db vdjdb --db iedb` both work."
        ),
    )
    p_data_clear.add_argument(
        "--group",
        choices=sorted(GROUPS),
        action="extend",
        nargs="+",
        help="Clear a whole workflow group (germline / germline_mouse / annotation / all).",
    )
    p_data_clear.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_clear.set_defaults(func=cmd_data_clear)

    p_data_download = p_data_sub.add_parser(
        "download",
        help="Download reference databases (VDJdb, IEDB, CEDAR, IMGT human/mouse germline).",
    )
    p_data_download.add_argument(
        "--db",
        choices=sorted(DATABASES),
        action="extend",
        nargs="+",
        help=(
            "Database(s) to download (default: all registered DBs). "
            "Repeatable and/or space-separated: "
            "`--db vdjdb iedb` and `--db vdjdb --db iedb` both work."
        ),
    )
    p_data_download.add_argument(
        "--group",
        choices=sorted(GROUPS),
        action="extend",
        nargs="+",
        help=(
            "Download a whole workflow group: `germline` (human IMGT V-REGION + J-REGION), "
            "`germline_mouse` (mouse IMGT V-REGION + J-REGION + leaders), "
            "`annotation` (VDJdb/IEDB/CEDAR), or `all`. Combines with --db."
        ),
    )
    p_data_download.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the database is already cached",
    )
    p_data_download.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_download.add_argument("--verbose", action="store_true", help="Verbose logging")
    p_data_download.set_defaults(func=cmd_data_download)
