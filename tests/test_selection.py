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
