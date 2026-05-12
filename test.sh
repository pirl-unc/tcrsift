#!/bin/bash
set -e

# Run lint first
bash lint.sh

export NUMBA_DISABLE_JIT="${NUMBA_DISABLE_JIT:-1}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp/numba_cache}"

pytest tests/ -v --cov=tcrsift --cov-report=term-missing "$@"
