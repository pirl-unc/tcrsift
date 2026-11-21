# Contributing

Bug fixes, documentation improvements, benchmarks, and new features are welcome
via GitHub pull requests.

## Development setup

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
pip install -e ".[dev,docs]"
```

TCRsift targets Python 3.9–3.12 on Linux and macOS.

## Running tests and linters

```bash
# Full test suite (add -n auto for parallel execution via pytest-xdist)
pytest tests/

# With coverage report
pytest tests/ --cov=tcrsift --cov-report=term-missing

# Ruff and Pylint checks
./lint.sh
```

## Building documentation locally

```bash
mkdocs serve          # Live preview at http://127.0.0.1:8000
mkdocs build --strict # CI check: fails on broken internal links or warnings
```

## Pull request checklist

- Document public functions and classes with
  [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
- Include regression tests alongside bug fixes and behavioral tests for new
  features.
- Ensure `./lint.sh`, `pytest`, and `mkdocs build --strict` pass cleanly.
- Follow the [Code of Conduct](https://github.com/pirl-unc/tcrsift/blob/main/code-of-conduct.md);
  contributions are licensed under Apache 2.0.
