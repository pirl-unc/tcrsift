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
            "activation",
            "antigen_response",
            "cytolytic",
            "exhaustion",
            "tumor_reactive",
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
