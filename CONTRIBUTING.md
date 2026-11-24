# Contributing to TCRsift

Thank you for your interest in contributing to TCRsift.

## Development setup

TCRsift supports Python 3.9–3.12. Install the package in editable mode with
development and documentation dependencies:

```bash
git clone https://github.com/pirl-unc/tcrsift.git
cd tcrsift
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,docs]"
```

## Checks before opening a pull request

```bash
./lint.sh
pytest tests/ -n auto
mkdocs build --strict
```

Set `TCRSIFT_LINT_INCLUDE_PYLINT=1 ./lint.sh` to also run the `pylint --errors-only`
check enforced by CI.

See [docs/contributing.md](docs/contributing.md) for full details and
[code-of-conduct.md](code-of-conduct.md) for community standards.
