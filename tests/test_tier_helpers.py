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

"""Tests for the #85 tier-eligibility helpers."""

from __future__ import annotations

from tcrsift.filter import (
    DEFAULT_THRESHOLD_TIERS,
    clone_clears_tier,
    strictest_tier_met,
)


class TestCloneClearsTier:
    def test_clears_when_cells_and_freq_above(self):
        tier = {"min_cells": 5, "min_frequency": 0.005}
        assert clone_clears_tier(10, 0.01, tier)

    def test_fails_when_cells_below(self):
        tier = {"min_cells": 5, "min_frequency": 0.005}
        assert not clone_clears_tier(4, 0.01, tier)

    def test_fails_when_freq_below(self):
        tier = {"min_cells": 5, "min_frequency": 0.005}
        assert not clone_clears_tier(10, 0.001, tier)

    def test_exact_boundary_clears(self):
        tier = {"min_cells": 5, "min_frequency": 0.005}
        assert clone_clears_tier(5, 0.005, tier)

    def test_n_conditions_check_optional(self):
        tier = {"min_cells": 5, "min_frequency": 0.005, "max_conditions": 2}
        assert clone_clears_tier(10, 0.01, tier, n_conditions=2)
        assert not clone_clears_tier(10, 0.01, tier, n_conditions=3)
        # When n_conditions is None, the max_conditions key is ignored.
        assert clone_clears_tier(10, 0.01, tier, n_conditions=None)

    def test_non_numeric_input_fails_safely(self):
        tier = {"min_cells": 5, "min_frequency": 0.005}
        assert not clone_clears_tier(None, 0.01, tier)
        assert not clone_clears_tier("not-a-number", 0.01, tier)


class TestStrictestTierMet:
    def test_returns_tier1_for_high_clone(self):
        # tier1 default: min_cells=10, min_frequency=0.01.
        assert strictest_tier_met(20, 0.02) == "tier1"

    def test_returns_lower_tier_when_strict_fails(self):
        # Clears tier2 (5 cells / 0.5%) but not tier1.
        assert strictest_tier_met(6, 0.006) == "tier2"

    def test_returns_none_when_all_fail(self):
        assert strictest_tier_met(0, 0.0) is None

    def test_respects_custom_tier_defs(self):
        my_tiers = {
            "elite": {"min_cells": 100, "min_frequency": 0.1},
            "ok": {"min_cells": 1, "min_frequency": 0.0},
        }
        assert strictest_tier_met(
            5, 0.005,
            tier_defs=my_tiers,
            tier_order=("elite", "ok"),
        ) == "ok"

    def test_default_tier_defs_is_public(self):
        """Per #85, ``DEFAULT_THRESHOLD_TIERS`` is exposed as public
        API for downstream callers."""
        assert "tier1" in DEFAULT_THRESHOLD_TIERS
        assert "min_cells" in DEFAULT_THRESHOLD_TIERS["tier1"]
