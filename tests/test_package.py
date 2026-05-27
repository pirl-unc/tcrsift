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


def test_root_export_table_drives_public_api():
    """The lazy export table should stay aligned with __all__ and dir()."""
    import tcrsift

    expected_exports = list(tcrsift._EXPORT_TO_MODULE)

    assert tcrsift.__all__ == expected_exports
    assert set(expected_exports).issubset(set(dir(tcrsift)))


def test_most_common_replaces_safe_mode_in_1_0():
    """1.0 breaking rename: ``safe_mode`` is gone; the same function is
    now exposed as ``most_common``. Lock both halves so a future
    accidental restoration of the old name doesn't sneak back in."""
    import pandas as pd

    import tcrsift

    assert hasattr(tcrsift, "most_common")
    assert tcrsift.most_common(pd.Series(["A", "A", "B"])) == "A"
    assert tcrsift.most_common(pd.Series([], dtype=object), default="x") == "x"

    # The old name must not be importable from the package root.
    assert not hasattr(tcrsift, "safe_mode")

    # And not from the submodule either — the rename should be in the
    # function definition itself, not just an export shim.
    from tcrsift import validation

    assert hasattr(validation, "most_common")
    assert not hasattr(validation, "safe_mode")
