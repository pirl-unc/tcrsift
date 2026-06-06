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

"""#201: `run` defaults assembly to the sample sheet's vdj_dir contigs."""

from __future__ import annotations

import os

from tcrsift.cli import _common_contig_root, create_parser


class TestCommonContigRoot:
    def test_common_ancestor_of_per_sample_dirs(self, tmp_path):
        base = tmp_path / "B1-2" / "per_sample_outs"
        dirs = [str(base / s / "vdj_t") for s in ("AIMpos", "CTYneg", "tetpos")]
        root = _common_contig_root(dirs)
        assert root == os.path.abspath(str(base))

    def test_single_dir_returns_itself(self, tmp_path):
        d = str(tmp_path / "x" / "vdj_t")
        assert _common_contig_root([d]) == os.path.abspath(d)

    def test_empty_returns_none(self):
        assert _common_contig_root([]) is None
        assert _common_contig_root([None, ""]) is None

    def test_root_level_common_path_refused(self):
        # Two dirs whose only common ancestor is "/" must be refused (rglob /
        # the whole filesystem is unsafe).
        assert _common_contig_root(["/aaa/vdj", "/bbb/vdj"]) is None


class TestNoSheetContigsFlag:
    def test_flag_parses(self):
        a = create_parser().parse_args(
            ["run", "-s", "sheet.csv", "-o", "out", "--no-sheet-contigs"]
        )
        assert a.no_sheet_contigs is True

    def test_default_is_off(self):
        a = create_parser().parse_args(["run", "-s", "sheet.csv", "-o", "out"])
        assert getattr(a, "no_sheet_contigs", False) is False
