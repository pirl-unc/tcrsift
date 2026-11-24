# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Shared helpers and constants for TCRsift CLI subcommands."""

from __future__ import annotations

import logging

LEADER_CHOICES: tuple[str, ...] = (
    "CD8A",
    "CD28",
    "mCD8B",
    "mCD28",
    "H2Kb",
    "IgK",
    "TRAC",
    "TRBC",
)


def setup_logging(verbose: bool = False) -> None:
    """Configure root logging for CLI invocations."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _format_size(n: int) -> str:
    """Human-friendly byte count, e.g. '12 MB', '4.2 GB'."""
    if n == 0:
        return "-"
    units = ("B", "KB", "MB", "GB", "TB")
    size = float(n)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.0f} {unit}" if size >= 10 or unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{n} B"
