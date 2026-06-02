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

"""Canonical T-cell gene-set signatures.

Single source of truth for the gene sets used across TIL selection
(``til_select.py``) and the per-sample signature scatter
(``plots.py``). Pulled out of ``til_select.py`` so they can be reused
for selecting other kinds of T cells — antigen-response screens,
exhaustion-state phenotyping, healthy-donor panels — without inheriting
the TIL-specific scoring code.

All gene symbols are HGNC (human). When working with Ensembl IDs,
translate via the 10x ``features.tsv.gz`` mapping.

Signatures grouped by intent:

- :data:`EFFECTOR_GENES_HGNC` — cytotoxic-effector panel
  (IFNG, GZMB, PRF1, GNLY, NKG7). Used in cytotoxicity readouts.
  ``ACTIVATION_GENES_HGNC`` is a deprecated alias (#142): this set is
  effector differentiation, *not* immediate-early activation (CD69,
  NR4A1-3, EGR1/2, FOS, JUN) — the old name invited that conflation.
- :data:`NAIVE_STEM_GENES_HGNC` — naïve / stem-memory program
  (TCF7, LEF1, CCR7, SELL, IL7R, CD27, CD28). The "down" pole of the
  effector−naïve differentiation contrast (#141), the best axis for
  separating antigen-expanded clones from naïve-like bystanders.
- :data:`ANTIGEN_RESPONSE_GENES_HGNC` — focal 2-gene recent-antigen
  marker: TNFRSF9 (4-1BB / CD137, the AIM-assay marker per Wölfl
  2007, Frentsch 2005, Bacher 2013) + MKI67 (Ki-67, proliferation).
  Note (#142): this acute/proliferation axis runs *inverse* to
  in-vitro clonal expansion at snapshot — don't use it to score
  "antigen-expanded".
- :data:`CYTOLYTIC_GENES_HGNC` — minimal canonical cytotoxic
  effector pair (PRF1 perforin, GZMB granzyme B). Used as the
  effector readout in Caushi 2021, Krishna 2021, Hanada 2022.
- :data:`EXHAUSTION_GENES_HGNC` — canonical exhausted-T-cell
  surface markers (PDCD1, LAG3, HAVCR2, TIGIT, TOX, CTLA4).
- :data:`TUMOR_REACTIVE_GENES_HGNC` — TIL-resident tumor-reactive
  phenotype: CXCL13 (Workel 2019, Cohen 2022, Veatch 2022) + ENTPD1
  / CD39 (Duhen 2018, Simoni 2018, Thommen 2018).

:data:`T_CELL_SIGNATURES` is a snake-case name → tuple dict for
convenient iteration.
"""

from __future__ import annotations

EFFECTOR_GENES_HGNC: tuple[str, ...] = ("IFNG", "GZMB", "PRF1", "GNLY", "NKG7")
# Deprecated alias (#142): mis-named — this is cytotoxic effector, not
# immediate-early activation. Kept for back-compat; prefer EFFECTOR_GENES_HGNC.
ACTIVATION_GENES_HGNC: tuple[str, ...] = EFFECTOR_GENES_HGNC
NAIVE_STEM_GENES_HGNC: tuple[str, ...] = (
    "TCF7", "LEF1", "CCR7", "SELL", "IL7R", "CD27", "CD28",
)
ANTIGEN_RESPONSE_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "MKI67")
CYTOLYTIC_GENES_HGNC: tuple[str, ...] = ("PRF1", "GZMB")
EXHAUSTION_GENES_HGNC: tuple[str, ...] = (
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4",
)
TUMOR_REACTIVE_GENES_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")

T_CELL_SIGNATURES: dict[str, tuple[str, ...]] = {
    "effector":          EFFECTOR_GENES_HGNC,
    "activation":        ACTIVATION_GENES_HGNC,  # deprecated alias of effector
    "naive_stem":        NAIVE_STEM_GENES_HGNC,
    "antigen_response":  ANTIGEN_RESPONSE_GENES_HGNC,
    "cytolytic":         CYTOLYTIC_GENES_HGNC,
    "exhaustion":        EXHAUSTION_GENES_HGNC,
    "tumor_reactive":    TUMOR_REACTIVE_GENES_HGNC,
}
