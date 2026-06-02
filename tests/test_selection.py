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

"""Tests for the tcrsift.selection consolidation module (#125 step 1)."""

from __future__ import annotations

import pandas as pd

from tcrsift import selection


class TestPublicSurface:
    def test_reexports_are_importable(self):
        # The whole point of the module is one import surface for the
        # shared primitives that previously lived across filter/clonotype.
        assert callable(selection.build_clone_sample_long)
        assert callable(selection.build_clone_method_long)
        assert callable(selection.per_sample_tier)
        assert callable(selection.attach_per_sample_tiers)
        assert "tier1" in selection.DEFAULT_THRESHOLD_TIERS

    def test_reexported_per_sample_tier_is_the_filter_impl(self):
        from tcrsift.filter import per_sample_tier as filter_impl

        assert selection.per_sample_tier is filter_impl


class TestAttachPerSampleTiers:
    def _long(self):
        # One row per (clone, sample) with within-sample cells + frequency,
        # shaped like build_clone_sample_long output.
        return pd.DataFrame(
            {
                "CDR3ab": ["A", "B", "C", "D"],
                "sample": ["S1", "S1", "S2", "S2"],
                "cells": [12, 3, 2, 1],
                "frequency": [0.02, 0.002, 0.0006, 0.5],
            }
        )

    def test_assigns_expected_abundance_tiers(self):
        out = selection.attach_per_sample_tiers(self._long())
        tiers = dict(zip(out["CDR3ab"], out["per_sample_tier"]))
        assert tiers["A"] == "tier1"   # cells>=10, freq>=0.01
        assert tiers["B"] == "tier3"   # cells>=3, freq>=0.001
        assert tiers["C"] == "tier4"   # cells>=2, freq>=0.0005
        assert tiers["D"] is None      # cells<2 clears no tier

    def test_does_not_mutate_input(self):
        long_df = self._long()
        selection.attach_per_sample_tiers(long_df)
        assert "per_sample_tier" not in long_df.columns

    def test_custom_out_column_name(self):
        out = selection.attach_per_sample_tiers(self._long(), out_col="ps_tier")
        assert "ps_tier" in out.columns
        assert "per_sample_tier" not in out.columns

    def test_empty_frame_yields_empty_tier_column(self):
        empty = pd.DataFrame(
            {
                "CDR3ab": pd.Series([], dtype=str),
                "sample": pd.Series([], dtype=str),
                "cells": pd.Series([], dtype=int),
                "frequency": pd.Series([], dtype=float),
            }
        )
        out = selection.attach_per_sample_tiers(empty)
        assert "per_sample_tier" in out.columns
        assert len(out) == 0

    def test_missing_required_column_raises_clear_error(self):
        bad = pd.DataFrame({"CDR3ab": ["A"], "sample": ["S1"], "cells": [5]})
        try:
            selection.attach_per_sample_tiers(bad)
        except KeyError as e:
            assert "frequency" in str(e)
        else:
            raise AssertionError("expected KeyError for missing 'frequency'")

    def test_round_trips_with_build_clone_sample_long(self):
        import anndata as ad
        import numpy as np

        rows = (
            [("S1", "CAVA", "CASS_A")] * 12
            + [("S1", "CAVB", "CASS_B")] * 1
        )
        n = len(rows)
        adata = ad.AnnData(
            X=np.zeros((n, 1), dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "sample": [r[0] for r in rows],
                    "CDR3_alpha": [r[1] for r in rows],
                    "CDR3_beta": [r[2] for r in rows],
                }
            ),
        )
        long_df = selection.build_clone_sample_long(adata)
        out = selection.attach_per_sample_tiers(long_df)
        assert "per_sample_tier" in out.columns
        assert len(out) == len(long_df)


def _cml(rows):
    """Build a clone-method-long frame: rows of (clone, method, tier, freq)."""
    return pd.DataFrame(
        rows, columns=["CDR3ab", "method", "tier", "max_freq_in_method"]
    )


class TestAttachMethodTiers:
    def test_attaches_per_method_tier(self):
        cml = pd.DataFrame(
            {
                "CDR3ab": ["A", "B"],
                "method": ["M1", "M1"],
                "cells_in_method": [12, 1],
                "max_freq_in_method": [0.02, 0.5],
            }
        )
        out = selection.attach_method_tiers(cml)
        tiers = dict(zip(out["CDR3ab"], out["tier"]))
        assert tiers["A"] == "tier1"   # 12 cells, 2% freq
        assert tiers["B"] is None      # 1 cell clears nothing

    def test_empty_frame(self):
        out = selection.attach_method_tiers(
            pd.DataFrame(columns=["CDR3ab", "method", "cells_in_method", "max_freq_in_method"])
        )
        assert "tier" in out.columns and len(out) == 0


class TestBuildSelectionRoutesShared:
    def test_shared_includes_tier1_2_ranked_by_frequency(self):
        cml = _cml([
            ("A", "M1", "tier1", 0.2),
            ("B", "M1", "tier2", 0.05),
            ("C", "M1", "tier3", 0.005),
        ])
        config = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"],
                                        "rank_by": "max_frequency"}}}
        out = selection.build_selection_rules(cml, config)
        assert list(out["CDR3ab"]) == ["A", "B"]          # tier3 C excluded
        assert list(out["selection_rule"]) == ["shared", "shared"]
        assert list(out["rank_within_rule"]) == [1, 2]   # 0.2 before 0.05
        assert list(out["global_rank"]) == [1, 2]


class TestBuildSelectionRoutesPrivate:
    def test_private_excludes_tier3plus_in_other_methods(self):
        cml = _cml([
            ("X", "AIMpos", "tier3", 0.01),   # tier3 in AIMpos only
            ("Y", "AIMpos", "tier3", 0.005),  # tier3 in BOTH -> disqualified
            ("Y", "IFNpos", "tier3", 0.005),
        ])
        config = {"rules": {"private": {
            "include_tier": "tier3",
            "exclude_tier3plus_in_other_methods": True,
            "top_n": 3, "rank_by": "frequency_in_method",
            "apply_to_methods": ["AIMpos", "IFNpos"],
        }}}
        out = selection.build_selection_rules(cml, config)
        private_aim = out[out["selection_rule"] == "private_AIMpos"]
        assert list(private_aim["CDR3ab"]) == ["X"]
        # Y is tier3 in both methods, so it is private to neither.
        assert "Y" not in set(out["CDR3ab"])

    def test_private_respects_top_n(self):
        cml = _cml([
            ("X", "AIMpos", "tier3", 0.03),
            ("Y", "AIMpos", "tier3", 0.02),
            ("Z", "AIMpos", "tier3", 0.01),
        ])
        config = {"rules": {"private": {
            "include_tier": "tier3", "top_n": 2,
            "apply_to_methods": ["AIMpos"],
        }}}
        out = selection.build_selection_rules(cml, config)
        assert list(out["CDR3ab"]) == ["X", "Y"]   # top 2 by frequency


class TestBuildSelectionRoutesPair:
    def test_pair_requires_tier_in_all_members(self):
        cml = _cml([
            ("P", "AIMpos", "tier3", 0.02),
            ("P", "AIMpos_CTYneg", "tier3", 0.02),
            ("Q", "AIMpos", "tier3", 0.01),  # missing from AIMpos_CTYneg
        ])
        config = {"rules": {"method_pair": {
            "pairs": {"AIM": ["AIMpos", "AIMpos_CTYneg"]},
            "require_tier_in_all_members": "tier3",
            "exclude_tier3plus_outside_pair": True,
            "top_n": 3, "rank_by": "mean_freq_across_pair",
        }}}
        out = selection.build_selection_rules(cml, config)
        cty = out[out["selection_rule"] == "method_pair_AIM"]
        assert list(cty["CDR3ab"]) == ["P"]
        assert "Q" not in set(out["CDR3ab"])


class TestBuildSelectionRoutesGlobalRank:
    def test_rule_order_interleave_and_dedup(self):
        cml = _cml([
            ("S1", "M_shared", "tier1", 0.3),
            ("S2", "M_shared", "tier2", 0.2),
            ("P", "AIMpos", "tier3", 0.05),
            ("P", "AIMpos_CTYneg", "tier3", 0.05),
            ("R", "CTYneg", "tier3", 0.01),
            # D is tier1 in CTYneg: shared must claim it, private must NOT.
            ("D", "CTYneg", "tier1", 0.4),
        ])
        config = {
            "rules": {
                "shared": {"include_tiers": ["tier1", "tier2"],
                           "rank_by": "max_frequency"},
                "method_pair": {"pairs": {"AIM": ["AIMpos", "AIMpos_CTYneg"]},
                             "require_tier_in_all_members": "tier3",
                             "exclude_tier3plus_outside_pair": True, "top_n": 3},
                "private": {"include_tier": "tier3", "top_n": 3,
                            "exclude_tier3plus_in_other_methods": True,
                            "apply_to_methods": ["CTYneg"]},
            },
            "global_rank": {"rule_order": ["shared", "method_pair", "private"]},
        }
        out = selection.build_selection_rules(cml, config)
        # D (0.4) and S1 (0.3) and S2 (0.2) shared; P pair; R private.
        assert list(out["selection_rule"]) == [
            "shared", "shared", "shared", "method_pair_AIM", "private_CTYneg",
        ]
        assert list(out["CDR3ab"]) == ["D", "S1", "S2", "P", "R"]
        assert list(out["global_rank"]) == [1, 2, 3, 4, 5]
        # D was tier1 in CTYneg but shared claimed it first — no dupes.
        assert out["CDR3ab"].is_unique

    def test_empty_input_yields_empty_frame_with_columns(self):
        out = selection.build_selection_rules(
            _cml([]), {"rules": {"shared": {"include_tiers": ["tier1"]}}}
        )
        assert len(out) == 0
        for col in ("selection_rule", "rank_within_rule", "global_rank"):
            assert col in out.columns


# Pilot-shaped config mirroring the B1-2 enrichment design (#122): 7
# methods, three A/B method pairs, shared/private/pair routes.
_PILOT_CONFIG = {
    "rules": {
        "shared": {"include_tiers": ["tier1", "tier2"], "rank_by": "max_frequency"},
        "method_pair": {
            "pairs": {
                "AIM": ["AIMpos", "AIMpos_CTYneg"],
                "IFN": ["IFNpos", "IFNpos_CTYneg"],
                "tet": ["tetpos", "CTYneg_tetpos"],
            },
            "require_tier_in_all_members": "tier3",
            "exclude_tier3plus_outside_pair": True, "top_n": 3,
        },
        "private": {
            "include_tier": "tier3", "exclude_tier3plus_in_other_methods": True,
            "top_n": 3,
            "apply_to_methods": ["AIMpos", "AIMpos_CTYneg", "CTYneg",
                                 "CTYneg_tetpos", "IFNpos", "IFNpos_CTYneg", "tetpos"],
        },
    },
    "global_rank": {"rule_order": ["shared", "method_pair", "private"]},
}


class TestSelectionRoutesPilotShape:
    """Schema checkpoint on B1-2-shaped data before selected_clones.csv
    bakes the route definitions in (#122/#125 increment 5)."""

    def _pilot_cml(self):
        rows = [
            ("MART1_a", "AIMpos", "tier1", 0.20), ("MART1_a", "tetpos", "tier1", 0.18),
            ("MART1_a", "IFNpos", "tier2", 0.05),
            ("MART1_b", "AIMpos", "tier2", 0.04), ("MART1_b", "CTYneg", "tier2", 0.03),
            ("pair_AIM", "AIMpos", "tier3", 0.012),
            ("pair_AIM", "AIMpos_CTYneg", "tier3", 0.011),
            ("pair_tet", "tetpos", "tier3", 0.010),
            ("pair_tet", "CTYneg_tetpos", "tier3", 0.009),
            ("priv_CTY", "CTYneg", "tier3", 0.008),
            ("priv_IFN", "IFNpos", "tier3", 0.007),
            ("noise", "AIMpos", "tier4", 0.0006),
        ]
        return _cml(rows)

    def test_routes_match_pilot_expectations(self):
        out = selection.build_selection_rules(self._pilot_cml(), _PILOT_CONFIG)
        route = dict(zip(out["CDR3ab"], out["selection_rule"]))
        assert route["MART1_a"] == "shared"
        assert route["MART1_b"] == "shared"
        assert route["pair_AIM"] == "method_pair_AIM"
        assert route["pair_tet"] == "method_pair_tet"
        assert route["priv_CTY"] == "private_CTYneg"
        assert route["priv_IFN"] == "private_IFNpos"
        # tier4-only noise selected by nothing
        assert "noise" not in route

    def test_global_rank_is_contiguous_rule_order(self):
        out = selection.build_selection_rules(self._pilot_cml(), _PILOT_CONFIG)
        assert list(out["global_rank"]) == list(range(1, len(out) + 1))
        # shared block precedes pair block precedes private block
        routes = list(out["selection_rule"])
        assert routes[0] == "shared" and routes[1] == "shared"
        assert routes[2].startswith("method_pair_") and routes[3].startswith("method_pair_")
        assert routes[4].startswith("private_") and routes[5].startswith("private_")

    def test_pair_clone_not_double_claimed_by_private(self):
        # pair_AIM is tier3 in AIMpos; private_AIMpos must NOT re-select it.
        out = selection.build_selection_rules(self._pilot_cml(), _PILOT_CONFIG)
        assert out["CDR3ab"].is_unique
        assert (out["selection_rule"] == "private_AIMpos").sum() == 0


class TestSelectFromCloneSampleLong:
    def test_end_to_end_from_long_table(self):
        # clone_sample_long shape: per (clone, sample) cells + frequency + method.
        csl = pd.DataFrame({
            "CDR3ab": ["A", "A", "B"],
            "sample": ["S1", "S2", "S1"],
            "method": ["AIMpos", "tetpos", "CTYneg"],
            "cells": [12, 11, 1],
            "frequency": [0.2, 0.18, 0.0006],
        })
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"],
                                     "rank_by": "max_frequency"}}}
        out = selection.select_from_clone_sample_long(csl, cfg)
        assert list(out["CDR3ab"]) == ["A"]   # B is tier-less noise
        assert out["selection_rule"].iloc[0] == "shared"

    def test_no_method_axis_yields_empty(self):
        csl = pd.DataFrame({
            "CDR3ab": ["A"], "sample": ["S1"], "cells": [12], "frequency": [0.2],
        })
        out = selection.select_from_clone_sample_long(
            csl, {"rules": {"shared": {"include_tiers": ["tier1"]}}}
        )
        assert out.empty

    def test_no_routes_yields_empty(self):
        csl = pd.DataFrame({
            "CDR3ab": ["A"], "sample": ["S1"], "method": ["AIMpos"],
            "cells": [12], "frequency": [0.2],
        })
        assert selection.select_from_clone_sample_long(csl, {}).empty


class TestPerMethodEvidence:
    def _cml(self):
        return pd.DataFrame({
            "CDR3ab": ["A", "A", "B"],
            "method": ["AIMpos", "tetpos", "CTYneg"],
            "cells_in_method": [12, 8, 3],
            "max_freq_in_method": [0.2, 0.1, 0.005],
            "tier": ["tier1", "tier2", "tier3"],
        })

    def test_extract_per_method_evidence_for_clone(self):
        ev = selection.extract_per_method_evidence(self._cml(), "A")
        assert list(ev["method"]) == ["AIMpos", "tetpos"]
        assert list(ev["tier"]) == ["tier1", "tier2"]

    def test_extract_restricts_and_orders_by_methods(self):
        ev = selection.extract_per_method_evidence(
            self._cml(), "A", methods=["tetpos", "AIMpos"]
        )
        assert list(ev["method"]) == ["tetpos", "AIMpos"]  # follows methods order

    def test_build_pdf_annotations_includes_route_header_and_evidence(self):
        sel = pd.DataFrame({
            "CDR3ab": ["A"], "selection_rule": ["shared"], "global_rank": [1],
        })
        ann = selection.build_pdf_annotations(self._cml(), selection_df=sel)
        assert any("selection: shared" in line for line in ann["A"])
        assert any("AIMpos" in line and "tier1" in line for line in ann["A"])
        # B has no route header (not in selection_df) but still has evidence.
        assert all("selection:" not in line for line in ann["B"])
        assert any("CTYneg" in line for line in ann["B"])

    def test_build_pdf_annotations_without_selection(self):
        ann = selection.build_pdf_annotations(self._cml())
        assert set(ann) == {"A", "B"}
        assert all("selection:" not in line for line in ann["A"])


class TestSharedBoundsAndExclusion:
    def _cml(self):
        return _cml([
            ("A", "M1", "tier1", 0.40),
            ("B", "M1", "tier1", 0.30),
            ("C", "M1", "tier2", 0.20),
            ("D", "M1", "tier2", 0.05),
        ])

    def test_shared_uncapped_by_default(self):
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"]}}}
        out = selection.build_selection_rules(self._cml(), cfg)
        assert list(out["CDR3ab"]) == ["A", "B", "C", "D"]

    def test_shared_top_n_caps_block(self):
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"], "top_n": 2}}}
        out = selection.build_selection_rules(self._cml(), cfg)
        assert list(out["CDR3ab"]) == ["A", "B"]  # top 2 by frequency

    def test_shared_min_frequency_floor(self):
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"],
                                     "min_frequency": 0.1}}}
        out = selection.build_selection_rules(self._cml(), cfg)
        assert list(out["CDR3ab"]) == ["A", "B", "C"]  # D (0.05) dropped

    def test_exclude_clones_removed_from_all_routes(self):
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"]}}}
        out = selection.build_selection_rules(
            self._cml(), cfg, exclude_clones={"A", "B"}
        )
        assert list(out["CDR3ab"]) == ["C", "D"]
        # excluded clones never appear, even though they were tier1
        assert "A" not in set(out["CDR3ab"])

    def test_select_from_long_passes_exclude_clones(self):
        csl = pd.DataFrame({
            "CDR3ab": ["A", "B"], "sample": ["S1", "S1"], "method": ["M1", "M1"],
            "cells": [12, 11], "frequency": [0.2, 0.18],
        })
        cfg = {"rules": {"shared": {"include_tiers": ["tier1", "tier2"]}}}
        out = selection.select_from_clone_sample_long(
            csl, cfg, exclude_clones={"A"}
        )
        assert list(out["CDR3ab"]) == ["B"]


class TestSelectionWithoutMethodPairs:
    """Not every study has paired methods — a config with only shared /
    private routes must work with no method_pair block (#122)."""

    def test_shared_and_private_only_no_pair_route(self):
        cml = _cml([
            ("S", "AIMpos", "tier1", 0.30),   # shared
            ("P", "CTYneg", "tier3", 0.01),   # private to CTYneg
            ("N", "AIMpos", "tier4", 0.0006),  # noise
        ])
        config = {
            "rules": {
                "shared": {"include_tiers": ["tier1", "tier2"]},
                "private": {"include_tier": "tier3", "top_n": 3,
                            "apply_to_methods": ["AIMpos", "CTYneg"]},
            },
            "global_rank": {"rule_order": ["shared", "private"]},
        }
        out = selection.build_selection_rules(cml, config)
        routes = dict(zip(out["CDR3ab"], out["selection_rule"]))
        assert routes["S"] == "shared"
        assert routes["P"] == "private_CTYneg"
        assert "N" not in routes
        # No method_pair route exists anywhere.
        assert not any(r.startswith("method_pair") for r in out["selection_rule"])

    def test_default_rule_order_skips_absent_method_pair(self):
        # Default rule_order includes method_pair, but with none configured
        # it's simply skipped — no error, no rows.
        cml = _cml([("S", "AIMpos", "tier1", 0.3)])
        out = selection.build_selection_rules(
            cml, {"rules": {"shared": {"include_tiers": ["tier1"]}}}
        )
        assert list(out["selection_rule"]) == ["shared"]
