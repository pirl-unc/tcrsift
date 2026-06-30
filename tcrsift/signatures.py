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
  Kept as the back-compat subset of :data:`AIM_GENES_HGNC`. Per the
  expansion-correlation analysis (#303), both members positively track
  *culture* clonal expansion (TNFRSF9 ρ≈0.28, MKI67 ρ≈0.46) despite
  TNFRSF9's sparsity — sparse-but-specific, not a reason to demote it.
- :data:`AIM_GENES_HGNC` — activation-induced-marker panel for the
  in-vitro restimulation / antigen-response axis: TNFRSF9 (4-1BB),
  TNFRSF4 (OX40), IL2RA (CD25) + MKI67 (proliferation). These track
  *culture* expansion but sit near zero in TIL (#303), so this is the
  culture-compartment activation readout, not an in-vivo state.
- :data:`CYTOLYTIC_GENES_HGNC` — minimal canonical cytotoxic
  effector pair (PRF1 perforin, GZMB granzyme B). Used as the
  effector readout in Caushi 2021, Krishna 2021, Hanada 2022. Tracks
  *TIL* expansion (ρ≈0.29–0.34) but is flat/negative in culture (#303).
- :data:`EXHAUSTION_GENES_HGNC` — canonical co-inhibitory /
  chronic-antigen exhausted-T-cell set (PDCD1, TOX, LAG3, HAVCR2,
  TIGIT). CTLA4 was dropped (#303): it is an activation/Treg marker,
  not specific to the exhausted state. Every member positively tracks
  *TIL* expansion in the H37 tumor — the expanded intratumoral clones
  are the chronically-stimulated ones. Residency/tumor-reactive
  extension: ENTPD1 (CD39), CXCL13, ITGAE (CD103).
- :data:`TUMOR_REACTIVE_GENES_HGNC` — TIL-resident tumor-reactive
  phenotype: CXCL13 (Workel 2019, Cohen 2022, Veatch 2022) + ENTPD1
  / CD39 (Duhen 2018, Simoni 2018, Thommen 2018).
- :data:`EXPANSION_CORE_GENES_HGNC` — cross-compartment expansion
  core (MKI67, TNFRSF9, EGR2, IFNG, CXCL13, HAVCR2): the markers whose
  clone-level expression tracks expansion in *both* the peptide-culture
  and TIL compartments (#303), so the one panel usable to score
  "expanding & antigen-engaged" regardless of compartment.
- :data:`MARKER_PANEL_HGNC` — the default per-clone GEX scoring panel.
  Widened (#303) to the union of the functional panels above (plus
  CD4/CD8 lineage) so every panel is computable from a single
  extraction pass. Not a single-intent signature (so it is
  intentionally absent from :data:`T_CELL_SIGNATURES`); it is the union
  display panel shown per clone. ``til_select`` re-exports it as
  ``MARKER_GENES_DEFAULT`` for back-compat.

Compartment note (#303): the panels mean different things by
compartment. AIM markers track *in-vitro culture* expansion; cytolytic
and exhaustion track *in-vivo TIL* state/expansion; the expansion core
tracks both. Detectability (% of cells expressing) is *not* the right
yardstick — CD69 is the most detectable activation marker yet
anti-correlates with expansion, while sparse TNFRSF9 works.

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
# Activation-induced markers (#303): the culture-compartment antigen-response
# axis. ANTIGEN_RESPONSE_GENES_HGNC above is its back-compat 2-gene subset.
AIM_GENES_HGNC: tuple[str, ...] = ("TNFRSF9", "TNFRSF4", "IL2RA", "MKI67")
CYTOLYTIC_GENES_HGNC: tuple[str, ...] = ("PRF1", "GZMB")
# Co-inhibitory / chronic-antigen exhaustion set (#303). CTLA4 dropped:
# it is an activation/Treg marker, not specific to the exhausted state.
EXHAUSTION_GENES_HGNC: tuple[str, ...] = (
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT",
)
TUMOR_REACTIVE_GENES_HGNC: tuple[str, ...] = ("CXCL13", "ENTPD1")
# Cross-compartment expansion core (#303): tracks clonal expansion in BOTH
# the peptide-culture and TIL compartments.
EXPANSION_CORE_GENES_HGNC: tuple[str, ...] = (
    "MKI67", "TNFRSF9", "EGR2", "IFNG", "CXCL13", "HAVCR2",
)

# Default per-clone GEX scoring panel. Widened (#303) to the union of the
# functional panels (plus CD4/CD8 lineage and the immediate-early NR4A/GITR
# context markers) so every panel is computable from one extraction pass. A
# union display panel, not a single intent signature — kept out of
# T_CELL_SIGNATURES below on purpose.
MARKER_PANEL_HGNC: tuple[str, ...] = (
    # lineage
    "CD4", "CD8A", "CD8B",
    # activation / AIM (+ immediate-early context: GITR, NR4A, EGR2)
    "MKI67", "TNFRSF9", "TNFRSF4", "IL2RA", "TNFRSF18", "NR4A1", "NR4A2",
    "EGR2", "CTLA4",
    # cytolytic / effector
    "PRF1", "GZMB", "GZMK", "NKG7", "GNLY", "IFNG",
    # exhaustion / chronic antigen (+ residency / tumor-reactive)
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "CXCL13", "ITGAE",
)

T_CELL_SIGNATURES: dict[str, tuple[str, ...]] = {
    "effector":          EFFECTOR_GENES_HGNC,
    "activation":        ACTIVATION_GENES_HGNC,  # deprecated alias of effector
    "naive_stem":        NAIVE_STEM_GENES_HGNC,
    "antigen_response":  ANTIGEN_RESPONSE_GENES_HGNC,
    "aim":               AIM_GENES_HGNC,
    "cytolytic":         CYTOLYTIC_GENES_HGNC,
    "exhaustion":        EXHAUSTION_GENES_HGNC,
    "tumor_reactive":    TUMOR_REACTIVE_GENES_HGNC,
    "expansion_core":    EXPANSION_CORE_GENES_HGNC,
}
