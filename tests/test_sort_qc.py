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

"""Tests for per-condition signature-consistency QC (#161)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tcrsift.sort_qc import (
    DEFAULT_SORT_SIGNATURE_MAP,
    parse_sort_label,
    sort_signature_consistency,
)


class TestParseSortLabel:
    def test_strips_donor_suffix(self):
        assert parse_sort_label("AIMpos-2") == ["AIMpos"]

    def test_compound_split(self):
        assert parse_sort_label("AIMpos_CTYneg-2") == ["AIMpos", "CTYneg"]

    def test_empty(self):
        assert parse_sort_label("") == []
        assert parse_sort_label(None) == []


class TestDefaultMap:
    def test_ctyneg_is_proliferation_not_cytolytic(self):
        # CTY = Cell Trace Yellow dilution = divided/expanded, NOT low cytolytic.
        assert DEFAULT_SORT_SIGNATURE_MAP["CTYneg"] == [("Proliferation", "high")]
        assert DEFAULT_SORT_SIGNATURE_MAP["AIMpos"] == [("AcuteActivation", "high")]


def _cohort(seed=0):
    rng = np.random.default_rng(seed)
    rows = []

    def cells(donor, method, sigs, n=200):
        for _ in range(n):
            rows.append({
                "donor": donor, "enrichment_method": method,
                "signature_AcuteActivation": rng.normal(sigs.get("AcuteActivation", 0), 1),
                "signature_Proliferation": rng.normal(sigs.get("Proliferation", 0), 1),
                "signature_AntigenExperienced": rng.normal(sigs.get("AntigenExperienced", 0), 1),
            })
    for d in ("B1-2", "B1-3"):
        s = d[-1]
        cells(d, f"AIMpos-{s}", {"AcuteActivation": 1.5})
        cells(d, f"CTYneg-{s}", {"Proliferation": 1.5})
        cells(d, f"IFNpos-{s}", {"AntigenExperienced": 1.5})
    # B1-4: nothing enriched (failed/weird donor)
    for m in ("AIMpos-4", "CTYneg-4", "IFNpos-4"):
        cells("B1-4", m, {})
    return pd.DataFrame(rows)


class TestSortSignatureConsistency:
    def test_consistent_donors_pass(self):
        out = sort_signature_consistency(_cohort()).set_index(
            ["donor", "enrichment_method"])
        assert out.loc[("B1-2", "CTYneg-2"), "within_donor_consistent"]
        assert out.loc[("B1-3", "AIMpos-3"), "within_donor_consistent"]
        assert out.loc[("B1-2", "AIMpos-2"), "warning"] == ""

    def test_unenriched_sort_flagged(self):
        out = sort_signature_consistency(_cohort()).set_index(
            ["donor", "enrichment_method"])
        r = out.loc[("B1-4", "CTYneg-4")]
        assert not r["within_donor_consistent"]
        assert "not enriched as labelled" in r["warning"]
        assert "Proliferation" in r["warning"]

    def test_cross_donor_outlier_flagged(self):
        out = sort_signature_consistency(_cohort())
        by_donor = out.groupby("donor")["donor_outlier"].any()
        assert not by_donor["B1-2"]
        assert not by_donor["B1-3"]
        assert by_donor["B1-4"]            # the weird donor
        corr = out.groupby("donor")["donor_pattern_corr"].first()
        assert corr["B1-2"] > 0.5
        assert corr["B1-4"] < 0.5

    def test_under_three_donors_no_cross_donor(self):
        df = _cohort()
        df = df[df["donor"].isin(["B1-2", "B1-3"])]
        out = sort_signature_consistency(df)
        assert out["donor_pattern_corr"].isna().all()
        assert not out["donor_outlier"].any()

    def test_unmapped_sort_not_flagged(self):
        df = _cohort()
        # tetpos has no signature in the map → reported but not consistency-flagged
        df = pd.concat([df, pd.DataFrame([{
            "donor": "B1-2", "enrichment_method": "tetpos-2",
            "signature_AcuteActivation": 0.0, "signature_Proliferation": 0.0,
            "signature_AntigenExperienced": 0.0,
        }])], ignore_index=True)
        out = sort_signature_consistency(df).set_index(
            ["donor", "enrichment_method"])
        r = out.loc[("B1-2", "tetpos-2")]
        assert r["within_donor_consistent"]   # no expectation → not flagged
        assert "not enriched" not in r["warning"]

    def test_missing_method_col_raises(self):
        with pytest.raises(ValueError, match="enrichment_method"):
            sort_signature_consistency(pd.DataFrame({"x": [1]}))

    def test_no_signature_cols_raises(self):
        with pytest.raises(ValueError, match="no signature score columns"):
            sort_signature_consistency(pd.DataFrame({
                "enrichment_method": ["AIMpos-2"], "donor": ["B1-2"]}))

    def test_custom_map(self):
        df = _cohort()
        custom = {"AIMpos": [("Proliferation", "high")]}  # wrong-on-purpose
        out = sort_signature_consistency(df, sort_signature_map=custom).set_index(
            ["donor", "enrichment_method"])
        # AIMpos was enriched for AcuteActivation, not Proliferation → flagged
        assert not out.loc[("B1-2", "AIMpos-2"), "within_donor_consistent"]
