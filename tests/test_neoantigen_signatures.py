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

"""Neoantigen-reactivity signature registry + dispatcher (#309).

The published signatures are NOT interchangeable weighted gene sums, so the
registry records each one's true structure (genes / sign / units / method /
citation) and the dispatcher routes each to the matching scorer. These tests
lock the registry shape and the MANAscore ``weighted_z`` formula.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from tcrsift import signature_methods as sm
from tcrsift import signatures


class TestNeoantigenConstants:
    def test_manascore_directions(self):
        # +CXCL13, +ENTPD1, -IL7R (log-normalized) — the reproducible part.
        assert signatures.MANASCORE_UP_HGNC == ("CXCL13", "ENTPD1")
        assert signatures.MANASCORE_DOWN_HGNC == ("IL7R",)
        assert signatures.MANASCORE_WEIGHTS_HGNC == {
            "CXCL13": 1.0, "ENTPD1": 1.0, "IL7R": -1.0,
        }

    def test_neotcr_set_sizes(self):
        # 243-gene CD8 / 40-gene CD4 sets (Lowery Table S10), deduped.
        assert len(signatures.NEOTCR8_GENES_HGNC) == 243
        assert len(signatures.NEOTCR4_GENES_HGNC) == 40
        assert len(set(signatures.NEOTCR8_GENES_HGNC)) == 243
        # A hyphenated HLA symbol survived embedding intact.
        assert "HLA-DPA1" in signatures.NEOTCR8_GENES_HGNC
        assert "CXCL13" in signatures.NEOTCR4_GENES_HGNC

    def test_neotcr_pbl_set(self):
        # 151-gene circulating (PBL) neoantigen-reactive CD8 set (Yossef Cancer
        # Cell 2023, Table S2D, cluster C9 avg_log2FC>=0.5), deduped and verbatim.
        pbl = signatures.NEOTCRPBL_GENES_HGNC
        assert len(pbl) == 151
        assert len(set(pbl)) == 151  # no duplicates
        assert pbl[0] == "PASK"  # top gene by log2FC
        # Hallmark blood-memory / inhibitory members present verbatim.
        for g in ("SELL", "LEF1", "KLF2", "CTLA4", "TIGIT", "ITGAE", "HLA-DRB1"):
            assert g in pbl
        # Distinct from the TIL-derived NeoTCR8 (different signature).
        assert set(pbl) != set(signatures.NEOTCR8_GENES_HGNC)


class TestRegistryStructure:
    def test_registry_records_method_units_citation(self):
        reg = sm.NEOANTIGEN_SIGNATURES
        assert set(reg) == {"MANAscore", "NeoTCR8", "NeoTCR4", "NeoTCR_PBL"}
        assert reg["MANAscore"].method == "weighted_z"
        assert reg["MANAscore"].units == "log1p"
        assert "39900903" in reg["MANAscore"].citation  # PMID
        for name in ("NeoTCR8", "NeoTCR4"):
            assert reg[name].method == "geneset_enrichment"
            assert reg[name].units == "ranks"
            assert "35113651" in reg[name].citation
        # NeoTCR_PBL is a geneset too, but from the Yossef Cancer Cell 2023 paper.
        assert reg["NeoTCR_PBL"].method == "geneset_enrichment"
        assert reg["NeoTCR_PBL"].units == "ranks"
        assert "38039963" in reg["NeoTCR_PBL"].citation

    def test_signature_defaults_backward_compatible(self):
        # Pre-existing signatures stay plain zscore/log1p (no behaviour change).
        s = sm.SIGNATURES["Cytolytic"]
        assert s.method == "zscore" and s.units == "log1p" and s.citation == ""


class TestScoreByName:
    def _mana_frame(self):
        return pd.DataFrame({
            "CXCL13": [0.0, 10.0, 0.0, 5.0],
            "ENTPD1": [0.0, 8.0, 0.0, 4.0],
            "IL7R":   [10.0, 0.0, 9.0, 0.0],
        })

    def test_manascore_matches_reference_formula(self):
        expr = self._mana_frame()
        got = sm.score_by_name(expr, "manascore")  # case-insensitive lookup

        def z(col):
            v = np.log1p(expr[col].to_numpy())
            return (v - v.mean()) / v.std(ddof=0)

        expected = (z("CXCL13") + z("ENTPD1") - z("IL7R")) / np.sqrt(3)
        assert np.allclose(got.to_numpy(), expected)

    def test_manascore_high_low_ordering(self):
        s = sm.score_by_name(self._mana_frame(), "MANAscore")
        # CXCL13/ENTPD1-high IL7R-low cell (1) scores above the inverse (0).
        assert s.iloc[1] > 0 > s.iloc[0]

    def test_unknown_signature_raises(self):
        try:
            sm.score_by_name(self._mana_frame(), "not_a_signature")
        except KeyError as e:
            assert "not_a_signature" in str(e)
        else:
            raise AssertionError("expected KeyError")

    def test_geneset_enrichment_frame_fallback(self, caplog):
        # On a bare frame (no full gene universe) geneset_enrichment falls back
        # to a mean-z proxy and warns — it must not raise.
        expr = pd.DataFrame({"CXCL13": [0.0, 5.0], "GZMB": [1.0, 9.0]})
        with caplog.at_level(logging.WARNING):
            got = sm.score_by_name(expr, "NeoTCR8", on_missing="ignore")
        assert len(got) == 2
        assert any("mean-z proxy" in r.message for r in caplog.records)

    def test_weighted_z_ignores_missing_down_gene(self):
        # IL7R absent → score from up genes only, still finite.
        expr = pd.DataFrame({"CXCL13": [0.0, 4.0], "ENTPD1": [0.0, 3.0]})
        got = sm.score_by_name(expr, "manascore", on_missing="ignore")
        assert np.isfinite(got.to_numpy()).all()
        assert got.iloc[1] > got.iloc[0]
