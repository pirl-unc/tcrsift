#!/usr/bin/env python3
"""Compatibility wrapper for ``tcrsift til-prioritize``.

Prefer the installed command:
    tcrsift til-prioritize samples.yaml -o candidates/
"""

import sys

from tcrsift.cli import main

if __name__ == "__main__":
    main(["til-prioritize", *sys.argv[1:]])
