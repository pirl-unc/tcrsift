"""Tests for package import behavior."""

from __future__ import annotations

import subprocess
import sys


def test_root_lazy_import_does_not_pull_scanpy():
    """Importing a light-weight symbol from the package root should stay light-weight."""
    script = """
import builtins
import importlib

real_import = builtins.__import__

def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
    if name.startswith("scanpy"):
        raise AssertionError("scanpy should not be imported for load_sample_sheet")
    return real_import(name, globals, locals, fromlist, level)

builtins.__import__ = guarded_import
package = importlib.import_module("tcrsift")
load_sample_sheet = package.load_sample_sheet
assert callable(load_sample_sheet)
"""

    result = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
