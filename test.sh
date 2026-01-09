#!/bin/bash
set -e

pytest tests/ -v --cov=tcrsift --cov-report=term-missing "$@"
