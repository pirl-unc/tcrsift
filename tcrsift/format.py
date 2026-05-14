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

"""Display formatters for sample-sheet method names (#77).

The raw enrichment-method strings in tcrsift sample sheets are
pipeline-friendly but ugly in figures: ``AIMpos``, ``CTYneg``,
``AIMpos_CTYneg``, ``tetpos``, ``IFNpos``. Built-in plots route every
method/sample axis label through :func:`pretty_method` so figures
show the readable ``AIM⁺``, ``CTY⁻``, ``AIM⁺CTY⁻`` forms.

Rules:

- Suffix ``pos`` → Unicode superscript ``⁺``.
- Suffix ``neg`` → Unicode superscript ``⁻``.
- For two-part names joined by ``_``, reorder so a ``CTY*`` marker
  appears second (``CTYneg_tetpos`` → ``tet⁺CTY⁻``). ``CTY`` is a
  baseline/exclusion gate in our cohort, so it reads more naturally
  after the positive marker.

Pass-through behavior: any name that doesn't match the rules is
returned unchanged. Non-string input is coerced via ``str()``.
"""

from __future__ import annotations

from collections.abc import Iterable


def _format_part(part: str) -> str:
    """One ``pos``/``neg`` translation for a single marker token."""
    if part.endswith("pos"):
        return part[:-3] + "⁺"
    if part.endswith("neg"):
        return part[:-3] + "⁻"
    return part


def pretty_method(name: str) -> str:
    """Format a tcrsift enrichment-method string for display.

    Examples:
        >>> pretty_method("AIMpos")
        'AIM⁺'
        >>> pretty_method("CTYneg")
        'CTY⁻'
        >>> pretty_method("AIMpos_CTYneg")
        'AIM⁺CTY⁻'
        >>> pretty_method("CTYneg_tetpos")
        'tet⁺CTY⁻'
    """
    if not isinstance(name, str):
        return str(name)
    parts = name.split("_")
    formatted = [_format_part(p) for p in parts]
    # CTY is the baseline gate — read it after the positive marker.
    if (
        len(formatted) == 2
        and "CTY" in formatted[0]
        and "CTY" not in formatted[1]
    ):
        formatted = [formatted[1], formatted[0]]
    return "".join(formatted)


def pretty_methods(names: Iterable[str]) -> list[str]:
    """Map :func:`pretty_method` over a list/Series of method names."""
    return [pretty_method(str(n)) for n in names]


def pretty_sample(name: str) -> str:
    """Format a sample name like ``AIMpos_CTYneg-2`` as the prettified
    method (``AIM⁺CTY⁻``), dropping the trailing ``-N`` donor suffix.

    Use in plots that label by sample-of-origin in a single-donor
    context (the donor is named elsewhere in the figure title).
    """
    if not isinstance(name, str):
        return str(name)
    stem = name.rsplit("-", 1)[0] if "-" in name else name
    if not stem:
        return name
    return pretty_method(stem)


def pretty_samples(names: Iterable[str]) -> list[str]:
    """Map :func:`pretty_sample` over a list/Series of sample names."""
    return [pretty_sample(str(n)) for n in names]
