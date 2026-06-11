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

"""Tests for the germline V-REGION framework comparison (#286 follow-up).

All tests inject a small synthetic reference (no network / no data cache), so
the comparison logic is fully covered without the bundled IMGT reference.
"""

from __future__ import annotations

import pandas as pd
import pytest

import tcrsift.vregion as vr
from tcrsift.assemble import annotate_vregion_divergence

# A synthetic germline framework for one gene, two alleles.
_GERM = "GAVVSQHPSWVICKSGTSVKIECRSLDFQATTMFWYRQFP"  # 40 aa "framework"


@pytest.fixture
def synthetic_reference(monkeypatch):
    """Inject a one-gene reference and restore the memo afterwards."""
    ref = {"TRBV20-1": [("01", "F", _GERM), ("02", "F", _GERM[:-1] + "L")]}
    monkeypatch.setattr(vr, "_REFERENCE", ref)
    yield ref


def test_parse_vregion_fasta_translates_and_skips_stops():
    # frame-0 clean record kept; a record with an internal stop dropped.
    fasta = (
        "<pre>\n"
        ">X|TRBV9*01|Homo sapiens|F|V-REGION|1..6|\n"
        "GCTGGT\n"  # -> AG
        ">X|TRBVbad*01|Homo sapiens|F|V-REGION|1..6|\n"
        "TGAGGT\n"  # -> *G (internal stop) -> skipped
        "</pre>"
    )
    ref = vr._parse_vregion_fasta(fasta)
    assert "TRBV9" in ref and ref["TRBV9"][0] == ("01", "F", "AG")
    assert "TRBVbad" not in ref


def test_germline_vregion_lookup(synthetic_reference):
    alleles = [a for a, _f, _aa in vr.germline_vregion("TRBV20-1")]
    assert alleles == ["01", "02"]
    assert vr.germline_vregion("TRBV99") == []  # absent gene


def test_compare_identical_framework(synthetic_reference):
    vdj = _GERM + "CASSXYZF"  # framework + CDR3
    res = vr.germline_compare_vregion(vdj, "CASSXYZF", "TRBV20-1")
    allele, g_aa, identity, diff = res
    assert allele == "01" and identity == 1.0 and diff == "identical"


def test_compare_single_substitution(synthetic_reference):
    # one framework SNP at position 3 (V->N)
    fr = _GERM[:2] + "N" + _GERM[3:]
    res = vr.germline_compare_vregion(fr + "CASSF", "CASSF", "TRBV20-1")
    allele, g_aa, identity, diff = res
    assert diff == f"{_GERM[2]}3N"  # germline residue + pos + donor residue
    assert identity < 1.0


def test_compare_returns_none_without_cdr3_locatable(synthetic_reference):
    # CDR3 not a substring of VDJ → framework can't be isolated → None
    assert vr.germline_compare_vregion(_GERM, "ZZZZ", "TRBV20-1") is None


def test_compare_returns_none_for_absent_gene(synthetic_reference):
    assert vr.germline_compare_vregion(_GERM + "CF", "CF", "TRBV99") is None


def test_compare_picks_closest_allele(synthetic_reference):
    # donor matches allele *02 (last residue L) exactly → identity 1.0 vs *02
    fr = _GERM[:-1] + "L"
    res = vr.germline_compare_vregion(fr + "CF", "CF", "TRBV20-1")
    allele, _g, identity, diff = res
    assert allele == "02" and identity == 1.0 and diff == "identical"


def test_annotate_noop_without_reference(monkeypatch):
    # reference empty → no columns added, assembly stays inert.
    monkeypatch.setattr(vr, "_REFERENCE", {})
    df = pd.DataFrame([{"VDJ_beta_aa": _GERM + "CF", "CDR3_beta": "CF",
                        "beta_v_gene": "TRBV20-1"}])
    out = annotate_vregion_divergence(df)
    assert "beta_vregion_vs_germline" not in out.columns


def test_annotate_adds_columns_with_reference(synthetic_reference):
    fr = _GERM[:2] + "N" + _GERM[3:]
    df = pd.DataFrame([
        {"VDJ_beta_aa": _GERM + "CASSF", "CDR3_beta": "CASSF", "beta_v_gene": "TRBV20-1"},
        {"VDJ_beta_aa": fr + "CASSF", "CDR3_beta": "CASSF", "beta_v_gene": "TRBV20-1"},
    ])
    out = annotate_vregion_divergence(df)
    assert list(out["beta_vregion_vs_germline"]) == ["identical", f"{_GERM[2]}3N"]
    assert out["beta_vregion_allele"].iloc[0] == "TRBV20-1*01"
    assert out["beta_vregion_donor_aa"].iloc[1] == fr


def test_annotate_then_collect_framework_rows(synthetic_reference):
    from tcrsift.assemble import collect_germline_variants
    fr = _GERM[:2] + "N" + _GERM[3:]
    df = pd.DataFrame([
        {"VDJ_beta_aa": fr + "CASSF", "CDR3_beta": "CASSF", "beta_v_gene": "TRBV20-1"},
    ])
    out = collect_germline_variants(annotate_vregion_divergence(df))
    fw = out[out["region"] == "framework"]
    assert len(fw) == 1
    assert fw.iloc[0]["variant"] == f"{_GERM[2]}3N"
    assert fw.iloc[0]["gene"] == "TRBV20-1"
