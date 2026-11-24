"""Structural and contract tests for tcrsift.cli, tcrsift.assemble, and tcrsift.plots subpackages."""

from __future__ import annotations

import ast
import importlib
import subprocess
import sys
from pathlib import Path

import tcrsift
from tcrsift import assemble, cli, plots


def test_subpackage_all_exports_resolve():
    """Every symbol listed in __all__ on cli/assemble/plots resolves to a real attribute."""
    for pkg in (cli, assemble, plots):
        assert hasattr(pkg, "__all__")
        assert len(pkg.__all__) > 0
        assert len(pkg.__all__) == len(set(pkg.__all__)), f"Duplicate entries in {pkg.__name__}.__all__"
        for name in pkg.__all__:
            assert not name.startswith("_"), f"Private symbol {name!r} in {pkg.__name__}.__all__"
            assert hasattr(pkg, name), f"{pkg.__name__}.__all__ missing attribute {name!r}"


def test_cli_leader_choices_match_assemble_default_leaders():
    """CLI leader choices stay synchronized with assemble.DEFAULT_LEADERS."""
    assert cli.LEADER_CHOICES == tuple(assemble.DEFAULT_LEADERS.keys())


def test_assemble_intra_package_imports_form_dag():
    """tcrsift.assemble submodules have no intra-package import cycles (even lazy ones)."""
    pkg_dir = Path(assemble.__file__).parent
    graph: dict[str, set[str]] = {}
    for py_file in sorted(pkg_dir.glob("*.py")):
        if py_file.name == "__init__.py":
            continue
        mod = py_file.stem
        deps: set[str] = set()
        tree = ast.parse(py_file.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                deps.add(node.module.split(".")[0])
        graph[mod] = deps

    # Standard DFS cycle check
    visiting: set[str] = set()
    visited: set[str] = set()

    def _visit(node: str) -> None:
        assert node not in visiting, f"Circular intra-package import involving tcrsift.assemble.{node}"
        if node in visited:
            return
        visiting.add(node)
        for dep in graph.get(node, ()):
            _visit(dep)
        visiting.remove(node)
        visited.add(node)

    for mod in graph:
        _visit(mod)


def test_python_m_tcrsift_cli_version():
    """python -m tcrsift.cli --version prints the package version cleanly."""
    proc = subprocess.run(
        [sys.executable, "-m", "tcrsift.cli", "--version"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert f"tcrsift {tcrsift.__version__}" in proc.stdout


def test_top_level_lazy_exports_all_resolve():
    """Every symbol in tcrsift.__all__ resolves cleanly across the refactored subpackages."""
    for name in tcrsift.__all__:
        assert getattr(tcrsift, name) is not None
    # Verify submodules are directly importable
    for mod_name in (
        "tcrsift.cli.pipeline",
        "tcrsift.cli.til",
        "tcrsift.assemble.construct",
        "tcrsift.assemble.contigs",
        "tcrsift.assemble.leader_policy",
        "tcrsift.plots.annotations",
        "tcrsift.plots.repertoire",
    ):
        assert importlib.import_module(mod_name) is not None
