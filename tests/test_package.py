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


def test_validation_module_imports_on_python_3_9():
    """Regression for 1.0.0: ``tcrsift.validation`` was importing-broken
    on Python 3.9 because ``pick_representative_cell`` used PEP 604
    ``pd.Series | None`` in its annotation without
    ``from __future__ import annotations``. CI (3.9) failed at module
    collection. The fix is the future-import; this test just exercises
    the import path so the regression can't sneak back in via the
    full-suite-passes-on-3.12-only blind spot."""
    # Importing the module is the test — if a 3.9-incompatible
    # annotation slips in, this will raise TypeError at collection.
    from tcrsift import validation

    assert hasattr(validation, "pick_representative_cell")
    assert hasattr(validation, "most_common")


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
