#!/bin/bash
set -e

echo "Running ruff..."
ruff check tcrsift/
echo "Ruff check passed!"

echo "Running pylint..."
pylint tcrsift/ \
  --errors-only \
  --disable=unsubscriptable-object,not-an-iterable,no-member,invalid-unary-operand-type
echo "Pylint check passed!"
