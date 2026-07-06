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
- A name combining several markers (joined by ``_``, or an N-way set
  combination joined by ``+``) is reordered so any *baseline* marker
  (default ``CTY``, an exclusion gate in our cohort) always reads last
  (``CTYneg_tetpos`` → ``tet⁺CTY⁻``). This gives a consistent label no
  matter which order the markers happened to appear in the data.

Baseline markers default to :data:`BASELINE_MARKERS` (``("CTY",)``) but
can be overridden globally via :func:`set_baseline_markers` or per call
via the ``last=`` argument. An explicit ``priority=`` list pins the
leading order of non-baseline markers.

Pass-through behavior: any name that doesn't match the rules is
returned unchanged. Non-string input is coerced via ``str()``.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

# Markers that should always read *last* in a multi-marker label — baseline /
# exclusion gates that read more naturally after the positive markers. Override
# globally with set_baseline_markers(), or per-call with the ``last=`` argument.
BASELINE_MARKERS: tuple[str, ...] = ("CTY",)


def set_baseline_markers(*markers: str) -> None:
    """Override the global baseline markers pushed last in combined labels.

    e.g. ``set_baseline_markers("CTY", "DMSO")`` makes both gates read last in
    any combined method/condition label. Call with no arguments to clear.
    """
    global BASELINE_MARKERS
    BASELINE_MARKERS = tuple(markers)


def _is_baseline(token: str, last: Sequence[str]) -> bool:
    return any(m and m in token for m in last)


# Antigen-agnostic display labels (#317). A condition / stimulation "pool" is a
# grouping token in the data (a sample, a peptide-pool number, an antigen
# code); its human-readable ANTIGEN label is opaque free text supplied by the
# user (via the sample sheet's ``antigen_name``, or a config map) — tcrsift
# never parses it, so it works for peptides, RNA/DNA-encoded antigens, whole
# protein, lysate, tetramer, etc. The rule everywhere: never render a bare
# condition token when an antigen label exists.
_ANTIGEN_LABELS: dict[str, str] = {}


def set_antigen_labels(mapping: dict | None = None, /, **kwargs: str) -> None:
    """Register a global condition-token → antigen-label map (#317).

    e.g. ``set_antigen_labels({"P2": "P2 (KIF1C)", "sampleA": "MART-1 26-35"})``.
    Values are opaque display strings (peptide pool, minigene, protein, …).
    Call with no arguments to clear. See :func:`pretty_antigen`.
    """
    global _ANTIGEN_LABELS
    merged = dict(mapping or {})
    merged.update(kwargs)
    _ANTIGEN_LABELS = {str(k): str(v) for k, v in merged.items()}


def condition_sort_key(
    name: str,
    *,
    priority: Sequence[str] | None = None,
    last: Sequence[str] | None = None,
) -> tuple:
    """Stable ordering key for a single marker/condition token (#208).

    Sorts so that: explicit ``priority`` markers lead (in the given order),
    then everything else, then ``last`` (baseline) markers — alphabetical
    within each band for determinism. Used to give a consistent order to the
    markers inside a combined label and to N-way overlap set combinations.
    """
    last = list(last) if last is not None else list(BASELINE_MARKERS)
    priority = list(priority or ())
    base = 1 if _is_baseline(name, last) else 0
    pidx = next((i for i, p in enumerate(priority) if p and p in name), len(priority))
    return (base, pidx, name)


def order_conditions(
    names: Iterable[str],
    *,
    priority: Sequence[str] | None = None,
    last: Sequence[str] | None = None,
) -> list[str]:
    """Order condition/marker names consistently (baseline markers last).

    Deterministic (alphabetical within each band) so the same set of conditions
    always renders in the same order regardless of input order — e.g.
    ``["CTYneg", "AIMpos"]`` and ``["AIMpos", "CTYneg"]`` both order as
    ``["AIMpos", "CTYneg"]``.
    """
    return sorted(
        (str(n) for n in names),
        key=lambda n: condition_sort_key(n, priority=priority, last=last),
    )


# Glyphs we emit for screen/matplotlib (TrueType fonts) that the reportlab
# base-14 PDF fonts (Helvetica/Courier, WinAnsi) can't render — they show as
# missing-glyph boxes. Mapped to ASCII before any PDF drawString (#202).
_PDF_TEXT_SUBS = {
    "⁺": "+",
    "⁻": "-",
    "β": "beta",
    "α": "alpha",
    "∩": "&",
    "∪": "|",
    "≥": ">=",
    "≤": "<=",
    "→": "->",
}


def pdf_safe(text: str) -> str:
    """Make text safe for the reportlab base-14 PDF fonts (#202).

    The figures use Unicode (``AIM⁺``, ``β-T2A-α``, ``∩``) that the standard
    PDF fonts lack, so they render as black boxes. This maps the known glyphs to
    ASCII (``AIM+``, ``beta-T2A-alpha``, ``&``) and drops anything else outside
    WinAnsi, so PDF text never shows a missing-glyph box. (Matplotlib output is
    unaffected — it keeps the nicer Unicode.)
    """
    s = str(text)
    for k, v in _PDF_TEXT_SUBS.items():
        s = s.replace(k, v)
    # Backstop: WinAnsi (cp1252) is reportlab's default base-14 encoding; it
    # keeps em-dash / curly quotes but replaces any remaining exotic glyph.
    return s.encode("cp1252", "replace").decode("cp1252")


def _format_part(part: str) -> str:
    """One ``pos``/``neg`` translation for a single marker token."""
    if part.endswith("pos"):
        return part[:-3] + "⁺"
    if part.endswith("neg"):
        return part[:-3] + "⁻"
    return part


def pretty_method(
    name: str,
    *,
    priority: Sequence[str] | None = None,
    last: Sequence[str] | None = None,
) -> str:
    """Format a tcrsift enrichment-method string for display.

    Markers joined by ``_`` are reordered so baseline markers (default
    ``CTY``) read last, for any number of parts (not just two). Override the
    baseline set with ``last=`` (or globally via :func:`set_baseline_markers`)
    and pin leading order with ``priority=``.

    Examples:
        >>> pretty_method("AIMpos")
        'AIM⁺'
        >>> pretty_method("CTYneg")
        'CTY⁻'
        >>> pretty_method("AIMpos_CTYneg")
        'AIM⁺CTY⁻'
        >>> pretty_method("CTYneg_tetpos")
        'tet⁺CTY⁻'
        >>> pretty_method("CTYneg_AIMpos_IFNpos")
        'AIM⁺IFN⁺CTY⁻'
    """
    if not isinstance(name, str):
        return str(name)
    parts = name.split("_")
    # Reorder on the RAW tokens (baseline detection is independent of the
    # superscript translation), then format each surviving part.
    ordered = order_conditions(parts, priority=priority, last=last)
    return "".join(_format_part(p) for p in ordered)


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


def pretty_antigen(
    condition: object, *, label: str | None = None, labels: dict | None = None
) -> str:
    """Human-readable antigen label for a condition / stimulation token (#317).

    Antigen-agnostic: the label is opaque free text (a peptide pool like
    ``"P2 (KIF1C)"``, an epitope, an RNA/DNA minigene, whole protein, a
    lysate…) — tcrsift never parses it. Resolution order (never a bare
    condition token when a label exists):

    1. an explicit ``label=`` argument,
    2. the per-call ``labels`` map (``condition → antigen``),
    3. the global map set via :func:`set_antigen_labels`,
    4. fallback to :func:`pretty_method` on the raw token.

    Route reactive-pool legends, per-culture axes, and candidate tables through
    this so a pool is never shown by number alone.
    """
    if label is not None and str(label).strip():
        return str(label)
    registry = labels if labels is not None else _ANTIGEN_LABELS
    key = str(condition)
    mapped = registry.get(key)
    if mapped is not None and str(mapped).strip():
        return str(mapped)
    return pretty_method(key)


def pretty_antigens(
    conditions: Iterable, *, labels: dict | None = None
) -> list[str]:
    """Map :func:`pretty_antigen` over a list/Series of condition tokens."""
    return [pretty_antigen(c, labels=labels) for c in conditions]
