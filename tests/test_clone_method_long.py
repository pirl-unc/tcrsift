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

"""Tests for #81 (per-clone, per-method long table) and #82
(sample × sample overlap matrices)."""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.clonotype import (
    build_clone_method_long,
    compute_sample_overlap_matrices,
)


class TestBuildCloneMethodLong:
    def test_aggregates_cells_sum_across_samples_within_method(self):
        long_df = pd.DataFrame([
            {"CDR3ab": "A", "sample": "S1", "method": "AIMpos", "cells": 5, "frequency": 0.1},
            {"CDR3ab": "A", "sample": "S2", "method": "AIMpos", "cells": 3, "frequency": 0.06},
            {"CDR3ab": "A", "sample": "S3", "method": "CTYneg", "cells": 2, "frequency": 0.04},
        ])
        out = build_clone_method_long(long_df).set_index(["CDR3ab", "method"])
        # cells.sum across the two AIMpos samples = 8.
        assert out.loc[("A", "AIMpos"), "cells_in_method"] == 8
        # freq.max across the two AIMpos samples = 0.1.
        assert abs(out.loc[("A", "AIMpos"), "max_freq_in_method"] - 0.1) < 1e-9
        # n_samples = 2 for AIMpos, 1 for CTYneg.
        assert out.loc[("A", "AIMpos"), "n_samples_in_method"] == 2
        assert out.loc[("A", "CTYneg"), "n_samples_in_method"] == 1

    def test_empty_input_returns_empty_frame_with_expected_columns(self):
        out = build_clone_method_long(pd.DataFrame(columns=["CDR3ab", "method"]))
        assert list(out.columns) == [
            "CDR3ab", "method", "cells_in_method",
            "max_freq_in_method", "n_samples_in_method",
        ]
        assert len(out) == 0

    def test_missing_method_column_raises(self):
        long_df = pd.DataFrame({"CDR3ab": ["A"], "sample": ["S1"], "cells": [1], "frequency": [0.1]})
        with pytest.raises(ValueError, match="method"):
            build_clone_method_long(long_df)


class TestComputeSampleOverlapMatrices:
    def _build(self, rows):
        return pd.DataFrame(rows)

    def test_jaccard_perfect_overlap(self):
        long_df = self._build([
            {"CDR3ab": "A", "sample": "S1", "cells": 5},
            {"CDR3ab": "B", "sample": "S1", "cells": 3},
            {"CDR3ab": "A", "sample": "S2", "cells": 4},
            {"CDR3ab": "B", "sample": "S2", "cells": 2},
        ])
        m = compute_sample_overlap_matrices(long_df)
        # Same clone set → Jaccard = 1.0 off-diagonal too.
        assert abs(m["jaccard"].loc["S1", "S2"] - 1.0) < 1e-9

    def test_jaccard_no_overlap(self):
        long_df = self._build([
            {"CDR3ab": "A", "sample": "S1", "cells": 5},
            {"CDR3ab": "B", "sample": "S2", "cells": 4},
        ])
        m = compute_sample_overlap_matrices(long_df)
        assert m["jaccard"].loc["S1", "S2"] == 0.0
        assert m["jaccard"].loc["S1", "S1"] == 1.0

    def test_cell_weighted_jaccard_abundance_aware(self):
        """Same clone set but different abundances → Jaccard=1.0 but
        cell-weighted Jaccard < 1.0."""
        long_df = self._build([
            {"CDR3ab": "A", "sample": "S1", "cells": 10},
            {"CDR3ab": "A", "sample": "S2", "cells": 1},
        ])
        m = compute_sample_overlap_matrices(long_df)
        assert m["jaccard"].loc["S1", "S2"] == 1.0
        # min/max = 1/10 = 0.1
        assert abs(m["cell_weighted_jaccard"].loc["S1", "S2"] - 0.1) < 1e-9

    def test_restrict_to_clones(self):
        long_df = self._build([
            {"CDR3ab": "A", "sample": "S1", "cells": 5},
            {"CDR3ab": "B", "sample": "S1", "cells": 5},
            {"CDR3ab": "A", "sample": "S2", "cells": 5},
        ])
        # Restrict to {A} → both samples share clone set {A} → Jaccard = 1.0.
        m = compute_sample_overlap_matrices(long_df, restrict_clones={"A"})
        assert m["jaccard"].loc["S1", "S2"] == 1.0
        # Without restriction, S1 has {A, B} and S2 has {A} → 1/2.
        m2 = compute_sample_overlap_matrices(long_df)
        assert m2["jaccard"].loc["S1", "S2"] == 0.5

    def test_diagonal_is_one(self):
        long_df = self._build([
            {"CDR3ab": "A", "sample": "S1", "cells": 5},
            {"CDR3ab": "A", "sample": "S2", "cells": 3},
        ])
        m = compute_sample_overlap_matrices(long_df)
        for s in ("S1", "S2"):
            assert m["jaccard"].loc[s, s] == 1.0
            assert m["cell_weighted_jaccard"].loc[s, s] == 1.0

    def test_empty_long_df_returns_empty_matrix(self):
        m = compute_sample_overlap_matrices(
            pd.DataFrame(columns=["CDR3ab", "sample", "cells"])
        )
        assert m["jaccard"].shape == (0, 0)
