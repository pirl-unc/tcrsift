#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: ./deploy.sh [--dry-run] [version]" >&2
}

DRY_RUN=0
VERSION_INPUT=""
UV_AVAILABLE=0

if command -v uv &> /dev/null; then
  UV_AVAILABLE=1
  export UV_CACHE_DIR="${UV_CACHE_DIR:-/tmp/uv_cache}"
  mkdir -p "${UV_CACHE_DIR}"
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "${VERSION_INPUT}" ]]; then
        VERSION_INPUT="$1"
        shift
      else
        echo "Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

CURRENT_VERSION="$(
  python3 - <<'PY'
import re
from pathlib import Path

text = Path("tcrsift/version.py").read_text()
match = re.search(r'__version__ = "([^"]+)"', text)
if not match:
    raise SystemExit("Failed to read tcrsift/version.py; unexpected format.")
print(match.group(1))
PY
)"

VERSION="${CURRENT_VERSION}"
if [[ -n "${VERSION_INPUT}" ]]; then
  VERSION="${VERSION_INPUT#v}"
fi

# Ensure a clean tree before deploy.
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree not clean; commit or stash changes before deploy." >&2
  exit 1
fi

./lint.sh
./test.sh

# Bump version, commit, and push if version was provided.
if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
  VERSION="${VERSION}" python3 - <<'PY'
import os
import re
from pathlib import Path

version = os.environ["VERSION"]
path = Path("tcrsift/version.py")
text = path.read_text()
new_text, n = re.subn(
    r'__version__ = "([^"]+)"',
    f'__version__ = "{version}"',
    text,
    count=1,
)
if n != 1:
    raise SystemExit("Failed to update tcrsift/version.py; unexpected format.")
path.write_text(new_text)
PY
  git add tcrsift/version.py
  git commit -m "Bump version to ${VERSION}"
  git push
fi

if [[ "${UV_AVAILABLE}" -eq 1 ]]; then
  uv pip install --upgrade build twine
else
  python3 -m pip install --upgrade build
  python3 -m pip install --upgrade twine
fi
rm -rf dist
python3 -m build
git --version

if [[ "${DRY_RUN}" -eq 1 ]]; then
  if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
    echo "Dry run: would bump version to ${VERSION}"
  fi
  echo "Dry run: ready to release v${VERSION}"
  exit 0
fi

python3 -m twine upload dist/*
git tag "v${VERSION}"
git push --tags
