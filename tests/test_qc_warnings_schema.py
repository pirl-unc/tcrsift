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

"""Tests for #278 — `qc_warnings` is always emitted (empty list per row when
no warning fired) so the assembled schema is stable across cohorts."""

from __future__ import annotations

import pandas as pd

from tcrsift.assemble import _ensure_qc_warnings_column


def test_adds_column_when_absent():
    df = pd.DataFrame({"clone_id": [1, 2, 3]})
    out = _ensure_qc_warnings_column(df)
    assert "qc_warnings" in out.columns
    assert all(v == [] for v in out["qc_warnings"])


def test_fills_missing_rows_with_empty_list():
    # one clone raised a warning, the others didn't -> NaN where absent
    df = pd.DataFrame({
        "clone_id": [1, 2, 3],
        "qc_warnings": [["fell back to canonical"], None, float("nan")],
    })
    out = _ensure_qc_warnings_column(df)
    assert out["qc_warnings"].iloc[0] == ["fell back to canonical"]
    assert out["qc_warnings"].iloc[1] == []
    assert out["qc_warnings"].iloc[2] == []
    # every cell is a list -> `isinstance(x, list)` consumers stay correct
    assert all(isinstance(v, list) for v in out["qc_warnings"])


def test_empty_frame_gets_column():
    df = pd.DataFrame({"clone_id": []})
    out = _ensure_qc_warnings_column(df)
    assert "qc_warnings" in out.columns
    assert len(out) == 0
