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

"""Tests for shared datacache.cli_help integration, argcomplete, and lazy CLI imports."""

from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path

import pytest
from datacache.cli_help import (
    ANSI_BOLD,
    ANSI_BOLD_CYAN,
    ANSI_CYAN,
    ANSI_RESET,
    HelpArgumentParser,
    HelpFormatter,
)

import tcrsift.cli as cli_pkg
from tcrsift.cli import create_parser, main
from tcrsift.cli.til import _PRIORITIZE_CONTEXTS, _PRIORITIZE_QC_DEFAULTS


class _Stream(io.StringIO):
    def __init__(self, tty: bool):
        super().__init__()
        self._tty = tty

    def isatty(self) -> bool:
        return self._tty


def test_create_parser_uses_shared_help_argument_parser():
    parser = create_parser()
    assert isinstance(parser, HelpArgumentParser)
    assert issubclass(parser.formatter_class, HelpFormatter)

    subparsers_action = next(a for a in parser._actions if hasattr(a, "choices") and isinstance(a.choices, dict))
    for sub in subparsers_action.choices.values():
        assert isinstance(sub, HelpArgumentParser)
        assert issubclass(sub.formatter_class, HelpFormatter)


def test_help_plain_vs_tty_and_no_color(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.delenv("TERM", raising=False)

    parser = create_parser()
    plain = parser.format_help()
    assert "\033[" not in plain
    assert "Primary workflows:" in plain

    tty_stream = _Stream(True)
    parser.print_help(file=tty_stream)
    styled = tty_stream.getvalue()
    assert ANSI_BOLD + "usage:" + ANSI_RESET in styled
    assert ANSI_CYAN + "--version" + ANSI_RESET in styled
    assert ANSI_BOLD_CYAN + "Primary workflows:" + ANSI_RESET in styled

    non_tty_stream = _Stream(False)
    parser.print_help(file=non_tty_stream)
    assert non_tty_stream.getvalue() == plain

    monkeypatch.setenv("NO_COLOR", "1")
    no_color_stream = _Stream(True)
    parser.print_help(file=no_color_stream)
    assert no_color_stream.getvalue() == plain


def test_argcomplete_marker_and_registration(monkeypatch):
    cli_dir = Path(cli_pkg.__file__).parent
    for fname in ("__init__.py", "__main__.py"):
        head = (cli_dir / fname).read_text(encoding="utf-8")[:1024]
        assert "PYTHON_ARGCOMPLETE_OK" in head, f"Missing PYTHON_ARGCOMPLETE_OK in {fname}"

    called: list[object] = []
    monkeypatch.setattr(cli_pkg, "enable_argcomplete", lambda p, **kw: called.append(p) or True)
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert len(called) == 1
    assert isinstance(called[0], HelpArgumentParser)


def test_cli_help_stays_import_light():
    code = (
        "import sys; "
        "import tcrsift.cli; "
        "p = tcrsift.cli.create_parser(); "
        "p.format_help(); "
        "heavy = [m for m in ('numpy', 'pandas', 'scipy', 'anndata', 'scanpy', "
        "'matplotlib', 'seaborn', 'torch', 'h5py', 'yaml') if m in sys.modules]; "
        "assert not heavy, f'Unexpected heavy imports on CLI help: {heavy}'"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr


def test_cli_prioritize_constants_match_prioritize_module():
    from tcrsift.prioritize import CONTEXT_SIGNATURES, QC_DEFAULTS

    assert _PRIORITIZE_CONTEXTS == tuple(CONTEXT_SIGNATURES.keys())
    assert _PRIORITIZE_QC_DEFAULTS == QC_DEFAULTS


def test_cli_lazy_pipeline_exports_resolve():
    from tcrsift.pipeline import _common_contig_root, run_pipeline

    assert cli_pkg.run_pipeline is run_pipeline
    assert cli_pkg._common_contig_root is _common_contig_root
