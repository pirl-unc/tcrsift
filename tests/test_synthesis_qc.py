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

"""Synthesis-hazard QC for assembled constructs (#206).

Covers the primitives (GC, homopolymer, repeat, restriction sites), the
`add_synthesis_qc` column pass (including cross-construct duplicate / α-β swap
flags and nt_col fallback), and the `synthesis_qc_report` tally text.
"""

from __future__ import annotations

import pandas as pd

from tcrsift.assemble import (
    SYNTHESIS_RESTRICTION_SITES,
    _gc_fraction,
    _longest_repeat,
    _max_homopolymer,
    _restriction_hits,
    add_synthesis_qc,
    synthesis_qc_report,
)


class TestPrimitives:
    def test_gc_fraction(self):
        assert _gc_fraction("GCGC") == 1.0
        assert _gc_fraction("ATAT") == 0.0
        assert _gc_fraction("ATGC") == 0.5
        assert _gc_fraction("") == 0.0

    def test_max_homopolymer(self):
        assert _max_homopolymer("AAAA") == 4
        assert _max_homopolymer("ATGCAAATG") == 3
        assert _max_homopolymer("ATGC") == 1
        assert _max_homopolymer("") == 0

    def test_longest_repeat(self):
        # "ACGTACGT" -> "ACGT" repeats but len 4 < min_len 8 -> 0
        assert _longest_repeat("ACGTACGT") == 0
        block = "ACGTACGTAC"  # 10 nt
        seq = block + "TTT" + block  # the 10-mer recurs
        assert _longest_repeat(seq) >= 10
        assert _longest_repeat("") == 0

    def test_restriction_hits(self):
        assert _restriction_hits("AAAGAATTCAAA") == "EcoRI"  # GAATTC
        assert _restriction_hits("ACGTACGT") == ""
        # Multiple sites are ;-joined.
        hit = _restriction_hits("GAATTCGGATCC")
        assert "EcoRI" in hit and "BamHI" in hit

    def test_restriction_hits_reverse_strand(self):
        # Non-palindromic Type IIS sites must be caught on the reverse strand:
        # BsmBI CGTCTC -> reverse complement GAGACG; BsaI GGTCTC -> GAGACC (F1).
        assert "BsmBI" in _restriction_hits("AAA" + "GAGACG" + "AAA")
        assert "BsaI" in _restriction_hits("TTT" + "GAGACC" + "TTT")
        # Palindromic sites already caught either way (EcoRI is its own RC).
        assert _restriction_hits("GAGACG").count(";") == 0  # only BsmBI


class TestAddSynthesisQc:
    def test_columns_added(self):
        df = pd.DataFrame(
            {"single_chain_nt_optimized": ["ATGC" * 20, "GAATTC" + "AAAAAAAAAA"]}
        )
        out = add_synthesis_qc(df)
        for col in (
            "synth_gc_fraction",
            "synth_max_homopolymer",
            "synth_max_repeat",
            "synth_restriction_sites",
            "synth_gc_ok",
        ):
            assert col in out.columns
        assert out.loc[1, "synth_restriction_sites"] == "EcoRI"
        assert out.loc[1, "synth_max_homopolymer"] >= 10

    def test_gc_window(self):
        df = pd.DataFrame(
            {"single_chain_nt": ["ATATATATAT", "GCGCGCGCGC", "ATGCATGCAT"]}
        )
        out = add_synthesis_qc(df, gc_range=(0.35, 0.70))
        assert not out.loc[0, "synth_gc_ok"]  # all-AT, GC=0
        assert not out.loc[1, "synth_gc_ok"]  # all-GC, GC=1
        assert out.loc[2, "synth_gc_ok"]      # GC=0.5

    def test_nt_col_fallback_prefers_optimized(self):
        df = pd.DataFrame(
            {
                "single_chain_nt": ["AAAAAAAA"],
                "single_chain_nt_optimized": ["ATGCATGC"],
            }
        )
        out = add_synthesis_qc(df)
        # Scored the optimized column (homopolymer 1, not 8).
        assert out.loc[0, "synth_max_homopolymer"] == 1

    def test_explicit_nt_col(self):
        df = pd.DataFrame({"my_nt": ["GAATTC"]})
        out = add_synthesis_qc(df, nt_col="my_nt")
        assert out.loc[0, "synth_restriction_sites"] == "EcoRI"

    def test_noop_without_nt_column(self):
        df = pd.DataFrame({"CDR3ab": ["CASSX_CAVX"]})
        out = add_synthesis_qc(df)
        assert "synth_gc_fraction" not in out.columns
        assert list(out.columns) == ["CDR3ab"]

    def test_empty_sequence_gc_ok(self):
        df = pd.DataFrame({"single_chain_nt": ["", "ATGCATGCAT"]})
        out = add_synthesis_qc(df)
        # Empty string isn't flagged as a GC failure.
        assert bool(out.loc[0, "synth_gc_ok"])

    def test_duplicate_construct(self):
        df = pd.DataFrame(
            {
                "single_chain_nt": ["ATGC", "ATGC", "GGCC"],
                "single_chain_aa": ["MX", "MX", "GA"],
            }
        )
        out = add_synthesis_qc(df)
        assert bool(out.loc[0, "synth_duplicate_construct"])
        assert bool(out.loc[1, "synth_duplicate_construct"])
        assert not bool(out.loc[2, "synth_duplicate_construct"])

    def test_alpha_beta_swap(self):
        df = pd.DataFrame(
            {
                "single_chain_nt": ["ATGC", "ATGC"],
                "CDR3_alpha": ["CAVRDX", "CAVRDX"],
                "CDR3_beta": ["CAVRDX", "CASSPX"],
            }
        )
        out = add_synthesis_qc(df)
        assert bool(out.loc[0, "synth_alpha_beta_swap"])
        assert not bool(out.loc[1, "synth_alpha_beta_swap"])


class TestSynthesisQcReport:
    def test_empty_when_not_run(self):
        assert synthesis_qc_report(pd.DataFrame({"x": [1]})) == ""

    def test_tally_text(self):
        df = pd.DataFrame(
            {"single_chain_nt": ["ATATATATAT", "GAATTC" + "AAAAAAAAA"]}
        )
        out = add_synthesis_qc(df)
        text = synthesis_qc_report(out)
        assert "[synthesis QC]" in text
        assert "GC outside window" in text
        assert "restriction site" in text


class TestRestrictionSiteCatalog:
    def test_known_enzymes_present(self):
        assert SYNTHESIS_RESTRICTION_SITES["EcoRI"] == "GAATTC"
        assert "BsaI" in SYNTHESIS_RESTRICTION_SITES
