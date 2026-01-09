# Contributing

We welcome contributions to TCRsift!

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/pirl-unc/tcrsift.git
   cd tcrsift
   ```

2. Install in development mode:
   ```bash
   pip install -e ".[dev]"
   ```

3. Install pre-commit hooks:
   ```bash
   pre-commit install
   ```

## Running Tests

```bash
pytest tests/ -v
```

With coverage:

```bash
pytest tests/ --cov=tcrsift --cov-report=term-missing
```

## Code Style

We use:

- **ruff** for linting
- **pylint** for additional checks

Run linting:

```bash
./lint.sh
```

## Documentation

Build documentation locally:

```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs serve
```

Then open http://localhost:8000.

## Submitting Changes

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests and linting
5. Submit a pull request

See [CONTRIBUTING.md](https://github.com/pirl-unc/tcrsift/blob/main/CONTRIBUTING.md) for detailed guidelines.
