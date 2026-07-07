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

"""Tests for the canonical T-cell signatures module.

Previously the focal signature gene sets lived in ``til_select.py`` and
``plots.py``. They've moved to ``tcrsift.signatures`` so they can drive
non-TIL selections too (antigen-response screens, exhaustion panels,
etc.) — these tests lock the shape of the new module and verify the
back-compat aliases on both call sites still resolve.
"""

from __future__ import annotations

import tcrsift
from tcrsift import signatures


class TestModuleShape:
    """The five canonical signatures must be exposed as tuples plus a
    convenience dict."""

    def test_all_five_signatures_exposed(self):
        for name in (
            "ACTIVATION_GENES_HGNC",
            "ANTIGEN_RESPONSE_GENES_HGNC",
            "CYTOLYTIC_GENES_HGNC",
            "EXHAUSTION_GENES_HGNC",
            "TUMOR_REACTIVE_GENES_HGNC",
        ):
            assert hasattr(signatures, name), f"missing {name}"
            assert isinstance(getattr(signatures, name), tuple)

    def test_t_cell_signatures_dict_complete(self):
        assert set(signatures.T_CELL_SIGNATURES) == {
            "effector",
            "activation",  # deprecated alias of effector (#142)
            "naive_stem",
            "antigen_response",
            "aim",  # #303 activation-induced markers
            "cytolytic",
            "exhaustion",
            "tumor_reactive",
            "expansion_core",  # #303 cross-compartment core
        }
        # Each value is a non-empty tuple of strings.
        for name, genes in signatures.T_CELL_SIGNATURES.items():
            assert isinstance(genes, tuple) and genes, name
            assert all(isinstance(g, str) and g.isupper() for g in genes), name

    def test_focal_signatures_have_canonical_gene_sets(self):
        """Lock the specific gene-set membership so downstream
        consumers can rely on it without re-checking."""
        assert signatures.ANTIGEN_RESPONSE_GENES_HGNC == ("TNFRSF9", "MKI67")
        assert signatures.CYTOLYTIC_GENES_HGNC == ("PRF1", "GZMB")
        assert signatures.TUMOR_REACTIVE_GENES_HGNC == ("CXCL13", "ENTPD1")


class TestExpansionCorrelationPanels:
    """#303: panels re-grounded on clone-level expansion correlation."""

    def test_aim_panel(self):
        # Activation-induced markers (culture compartment).
        assert signatures.AIM_GENES_HGNC == ("TNFRSF9", "TNFRSF4", "IL2RA", "MKI67")

    def test_antigen_response_is_aim_subset(self):
        # Kept as the 2-gene back-compat subset of AIM.
        assert set(signatures.ANTIGEN_RESPONSE_GENES_HGNC).issubset(
            set(signatures.AIM_GENES_HGNC)
        )

    def test_exhaustion_dropped_ctla4(self):
        assert signatures.EXHAUSTION_GENES_HGNC == (
            "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT",
        )
        assert "CTLA4" not in signatures.EXHAUSTION_GENES_HGNC

    def test_expansion_core_panel(self):
        assert signatures.EXPANSION_CORE_GENES_HGNC == (
            "MKI67", "TNFRSF9", "EGR2", "IFNG", "CXCL13", "HAVCR2",
        )

    def test_activation_alias_untouched(self):
        # The new AIM panel must NOT collide with the deprecated
        # ACTIVATION alias, which still means effector (#142).
        assert signatures.ACTIVATION_GENES_HGNC == signatures.EFFECTOR_GENES_HGNC
        assert signatures.AIM_GENES_HGNC != signatures.ACTIVATION_GENES_HGNC

    def test_new_panels_top_level_importable(self):
        assert tcrsift.AIM_GENES_HGNC == signatures.AIM_GENES_HGNC
        assert tcrsift.EXPANSION_CORE_GENES_HGNC == signatures.EXPANSION_CORE_GENES_HGNC
        assert tcrsift.MARKER_PANEL_HGNC == signatures.MARKER_PANEL_HGNC


class TestTopLevelExports:
    """The signatures should be importable from the top-level
    ``tcrsift`` package via the lazy __getattr__."""

    def test_top_level_import_works(self):
        # These should all resolve through ``tcrsift/__init__.py``.
        assert tcrsift.ANTIGEN_RESPONSE_GENES_HGNC == ("TNFRSF9", "MKI67")
        assert tcrsift.CYTOLYTIC_GENES_HGNC == ("PRF1", "GZMB")
        assert tcrsift.TUMOR_REACTIVE_GENES_HGNC == ("CXCL13", "ENTPD1")
        assert tcrsift.ACTIVATION_GENES_HGNC[:2] == ("IFNG", "GZMB")
        assert tcrsift.EXHAUSTION_GENES_HGNC[0] == "PDCD1"

    def test_t_cell_signatures_dict_top_level(self):
        assert "antigen_response" in tcrsift.T_CELL_SIGNATURES


class TestSingleSourceOfTruth:
    """``plots`` and ``til_select`` both re-export from ``signatures``;
    the constants they expose must be the *same object* (not a
    duplicate tuple), so any change to the canonical definition can't
    silently fail to propagate."""

    def test_plots_re_exports_same_object(self):
        from tcrsift.plots import ANTIGEN_RESPONSE_GENES_HGNC

        assert ANTIGEN_RESPONSE_GENES_HGNC is signatures.ANTIGEN_RESPONSE_GENES_HGNC

    def test_til_select_alias_same_object(self):
        from tcrsift.til_select import ANTIGEN_RESPONSE_GENES_DEFAULT

        assert ANTIGEN_RESPONSE_GENES_DEFAULT is signatures.ANTIGEN_RESPONSE_GENES_HGNC

    def test_til_select_marker_panel_alias_same_object(self):
        from tcrsift.til_select import MARKER_GENES_DEFAULT

        assert MARKER_GENES_DEFAULT is signatures.MARKER_PANEL_HGNC


class TestMarkerPanel:
    """The default per-clone GEX scoring panel lives in ``signatures``
    and is re-exported from ``til_select`` as ``MARKER_GENES_DEFAULT``."""

    def test_panel_has_lineage_and_no_dupes(self):
        panel = signatures.MARKER_PANEL_HGNC
        assert {"CD4", "CD8A", "CD8B"}.issubset(panel)  # lineage markers
        assert len(panel) == len(set(panel))  # no duplicates
        assert all(g.isupper() for g in panel)

    def test_panel_is_union_of_functional_panels(self):
        # #303: the panel was widened so every functional panel is
        # computable from one extraction pass — i.e. it must be a
        # superset of each. This guards against a panel gene being added
        # without also widening the extraction.
        panel = set(signatures.MARKER_PANEL_HGNC)
        for genes in (
            signatures.AIM_GENES_HGNC,
            signatures.ANTIGEN_RESPONSE_GENES_HGNC,
            signatures.CYTOLYTIC_GENES_HGNC,
            signatures.EXHAUSTION_GENES_HGNC,
            signatures.TUMOR_REACTIVE_GENES_HGNC,
            signatures.EXPANSION_CORE_GENES_HGNC,
        ):
            assert set(genes).issubset(panel), set(genes) - panel

    def test_panel_is_not_an_intent_signature(self):
        # The union display panel is intentionally excluded from the
        # snake-case intent-signature dict.
        assert signatures.MARKER_PANEL_HGNC not in signatures.T_CELL_SIGNATURES.values()

    def test_cli_default_matches_panel(self):
        # The ``til-select --marker-genes`` default is derived from the
        # canonical panel, so the two can't drift apart.
        from tcrsift.cli import create_parser

        parser = create_parser()
        ns = parser.parse_args(["til-select"])
        assert ns.marker_genes == ",".join(signatures.MARKER_PANEL_HGNC)


class TestSolidTumorCellTypePreset:
    """Blood/PBMC default flagged + a solid-tumor cell-type preset parallel to
    the SOLID_TUMOR_LINEAGE_PROGRAMS QC preset (#340)."""

    def test_pbmc_alias_is_the_default(self):
        assert signatures.PBMC_CELL_TYPE_SIGNATURES is signatures.CELL_TYPE_SIGNATURES

    def test_solid_tumor_preset_extends_default(self):
        base = set(signatures.CELL_TYPE_SIGNATURES)
        solid = set(signatures.SOLID_TUMOR_CELL_TYPE_SIGNATURES)
        assert base < solid  # strict superset
        added = solid - base
        assert added == {
            "Mesothelial", "Osteoclast", "Skeletal muscle", "Adipocyte",
            "Schwann/nerve",
        }

    def test_solid_tumor_preset_preserves_shared_gene_lists(self):
        for t, genes in signatures.CELL_TYPE_SIGNATURES.items():
            assert signatures.SOLID_TUMOR_CELL_TYPE_SIGNATURES[t] == genes

    def test_tumor_not_a_signature(self):
        # Malignant cells stay MarkerCountOverride territory, not a gene set.
        for t in signatures.SOLID_TUMOR_CELL_TYPE_SIGNATURES:
            assert "tumor" not in t.lower()

    def test_all_gene_lists_nonempty(self):
        for genes in signatures.SOLID_TUMOR_CELL_TYPE_SIGNATURES.values():
            assert genes and all(isinstance(g, str) for g in genes)

    def test_exposed_at_package_top_level(self):
        assert tcrsift.SOLID_TUMOR_CELL_TYPE_SIGNATURES is (
            signatures.SOLID_TUMOR_CELL_TYPE_SIGNATURES
        )
        assert tcrsift.PBMC_CELL_TYPE_SIGNATURES is signatures.CELL_TYPE_SIGNATURES
