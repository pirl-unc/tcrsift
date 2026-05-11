#!/bin/bash
set -e

# ruff is the default linter — it's fast and has stable Python 3.12 support,
# so it's safe to gate `test.sh` on it. pylint is opt-in via
# TCRSIFT_LINT_INCLUDE_PYLINT=1 because its astroid dep has lagged Python
# language features in the past (PEP 695 type-alias parsing broke in
# older releases), making it a fragile blocker on shared dev envs.

echo "Running ruff..."
ruff check tcrsift/
echo "Ruff check passed!"

if [[ "${TCRSIFT_LINT_INCLUDE_PYLINT:-0}" == "1" ]]; then
  echo "Running pylint..."
  pylint tcrsift/ \
    --errors-only \
    --disable=unsubscriptable-object,not-an-iterable,no-member,invalid-unary-operand-type
  echo "Pylint check passed!"
else
  echo "Skipping pylint (set TCRSIFT_LINT_INCLUDE_PYLINT=1 to include)."
fi
