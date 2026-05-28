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
    per_sample_tier,
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


class TestAbundanceOnlyDefaultTiers:
    """#99: the bundled defaults dropped ``max_conditions`` to make the
    tier label an abundance-only signal. Specificity moved to explicit
    downstream filtering via ``n_conditions``."""

    def test_defaults_carry_no_max_conditions(self):
        for name, tier_def in DEFAULT_THRESHOLD_TIERS.items():
            assert "max_conditions" not in tier_def, (
                f"{name} still carries max_conditions — that defeats "
                "the abundance-only tier semantics from #99"
            )
            # And the abundance keys are present.
            assert "min_cells" in tier_def
            assert "min_frequency" in tier_def

    def test_n_conditions_ignored_on_default_tiers(self):
        """Passing ``n_conditions`` to ``strictest_tier_met`` against
        the bundled defaults should be a no-op (no max_conditions key
        means the specificity check is skipped). Same clone should
        get the same tier regardless of how many conditions it spans."""
        without_cond = strictest_tier_met(20, 0.02)
        with_cond_1 = strictest_tier_met(20, 0.02, n_conditions=1)
        with_cond_999 = strictest_tier_met(20, 0.02, n_conditions=999)
        assert without_cond == with_cond_1 == with_cond_999 == "tier1"

    def test_user_supplied_max_conditions_still_honoured(self):
        """The clone_clears_tier predicate still honours
        ``max_conditions`` when present in a *user-supplied* tier dict
        — only the bundled defaults dropped it."""
        custom = {
            "tier1": {"min_cells": 10, "min_frequency": 0.01, "max_conditions": 2},
        }
        # Without n_conditions the cap is ignored.
        assert strictest_tier_met(20, 0.02, tier_defs=custom,
                                  tier_order=("tier1",)) == "tier1"
        # With n_conditions over the cap, the clone is demoted.
        assert strictest_tier_met(20, 0.02, n_conditions=3, tier_defs=custom,
                                  tier_order=("tier1",)) is None


class TestPerSampleTier:
    """#99: per-sample tier helper for downstream report code so PDF /
    heatmap consumers stop re-implementing the threshold table."""

    def test_tier1_for_high_per_sample_signal(self):
        # tier1 default: 10 cells, 0.01 freq.
        assert per_sample_tier(20, 0.02) == "tier1"

    def test_tier3_for_borderline(self):
        # tier3 default: 3 cells, 0.001 freq.
        assert per_sample_tier(4, 0.002) == "tier3"

    def test_none_when_below_all_tiers(self):
        # min_cells everywhere is >= 2; this clears nothing.
        assert per_sample_tier(0, 0.0) is None

    def test_specificity_not_consulted_per_sample(self):
        """Per-sample tier is by definition single-sample; there's no
        ``n_conditions`` concept to apply. The helper exposes no
        ``n_conditions`` parameter, so callers can't accidentally pass
        one and produce a confusing mixed-semantics result."""
        import inspect

        sig = inspect.signature(per_sample_tier)
        assert "n_conditions" not in sig.parameters

    def test_respects_custom_tier_defs(self):
        my_tiers = {
            "elite": {"min_cells": 100, "min_frequency": 0.1},
            "ok": {"min_cells": 1, "min_frequency": 0.0},
        }
        assert per_sample_tier(5, 0.005,
                               tier_defs=my_tiers,
                               tier_order=("elite", "ok")) == "ok"
        assert per_sample_tier(150, 0.2,
                               tier_defs=my_tiers,
                               tier_order=("elite", "ok")) == "elite"
