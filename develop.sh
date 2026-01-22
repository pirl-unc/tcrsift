#!/bin/bash
# Development setup script - source this file: `source develop.sh` or `. develop.sh`
#
# This script:
# 1. Activates the virtual environment if not already active
# 2. Installs the package in editable mode using uv

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if being sourced or executed
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${YELLOW}Warning: This script should be sourced to activate the venv in your current shell.${NC}"
    echo -e "${YELLOW}Usage: source develop.sh${NC}"
    echo ""
    echo "Continuing anyway (venv will only be active within this script)..."
fi

# Function to find and activate venv
activate_venv() {
    # Already in a venv?
    if [[ -n "$VIRTUAL_ENV" ]]; then
        echo -e "${GREEN}Virtual environment already active: $VIRTUAL_ENV${NC}"
        return 0
    fi

    # Check for local .venv
    if [[ -d ".venv" ]]; then
        echo -e "${GREEN}Activating local .venv...${NC}"
        source .venv/bin/activate
        return 0
    fi

    # Check for local venv
    if [[ -d "venv" ]]; then
        echo -e "${GREEN}Activating local venv...${NC}"
        source venv/bin/activate
        return 0
    fi

    # No venv found - create one with uv
    if command -v uv &> /dev/null; then
        echo -e "${YELLOW}No virtual environment found. Creating .venv with uv...${NC}"
        uv venv .venv
        source .venv/bin/activate
        return 0
    fi

    # Fallback to python -m venv
    echo -e "${YELLOW}No virtual environment found. Creating .venv with python...${NC}"
    python3 -m venv .venv
    source .venv/bin/activate
    return 0
}

# Activate venv
activate_venv

# Install in editable mode
echo -e "${GREEN}Installing package in editable mode...${NC}"
if command -v uv &> /dev/null; then
    uv pip install -e .
else
    pip install -e .
fi

echo -e "${GREEN}Done! Development environment ready.${NC}"
