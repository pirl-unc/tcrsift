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

"""Tests for load_contigs sample-name policies (#124).

Covers the CellRanger ``per_sample_outs/{sample}/vdj_t/`` layout
(grandparent), the legacy symlinked layout (parent), sample-sheet-driven
naming, the ``--cellranger-dir`` shorthand, and the error paths.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tcrsift.assemble import load_contigs


def _write(path, text=">c\nACGT\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _build_cellranger_layout(root):
    """per_sample_outs/{sample}/vdj_t/filtered_contig.fasta."""
    for sample in ("AIMpos-2", "CTYneg-2"):
        _write(
            root / sample / "vdj_t" / "filtered_contig.fasta",
            f">c_{sample}\nACGT\n",
        )


class TestSampleNameFrom:
    def test_grandparent_layout_matches_per_sample_contigs(self, tmp_path):
        _build_cellranger_layout(tmp_path)
        contigs = load_contigs(tmp_path, sample_name_from="grandparent")
        assert set(contigs) == {"AIMpos-2", "CTYneg-2"}
        assert "c_AIMpos-2" in contigs["AIMpos-2"]
        assert "c_CTYneg-2" in contigs["CTYneg-2"]

    def test_parent_layout_backward_compat(self, tmp_path):
        for sample in ("AIMpos-2", "CTYneg-2"):
            _write(tmp_path / sample / "filtered_contig.fasta", f">c_{sample}\nACGT\n")
        contigs = load_contigs(tmp_path)  # default parent
        assert set(contigs) == {"AIMpos-2", "CTYneg-2"}

    def test_sample_name_from_sheet(self, tmp_path):
        sheet = pd.DataFrame([
            {"sample": "AIMpos-2", "vdj_dir": str(tmp_path / "weird-name-X")},
            {"sample": "CTYneg-2", "vdj_dir": str(tmp_path / "weird-name-Y")},
        ])
        _write(tmp_path / "weird-name-X" / "filtered_contig.fasta", ">cX\nACGT\n")
        _write(tmp_path / "weird-name-Y" / "filtered_contig.fasta", ">cY\nACGT\n")
        contigs = load_contigs(tmp_path, sample_name_from="sheet", sample_sheet=sheet)
        assert "cX" in contigs["AIMpos-2"]
        assert "cY" in contigs["CTYneg-2"]

    def test_cellranger_dir_is_shorthand_for_grandparent(self, tmp_path):
        _build_cellranger_layout(tmp_path)
        c1 = load_contigs(cellranger_dir=tmp_path)
        c2 = load_contigs(tmp_path, sample_name_from="grandparent")
        assert c1 == c2

    def test_grandparent_does_not_collapse_to_vdj_t(self, tmp_path):
        # The bug #124 fixes: parent naming collapses every sample to "vdj_t".
        _build_cellranger_layout(tmp_path)
        parent = load_contigs(tmp_path, sample_name_from="parent")
        assert set(parent) == {"vdj_t"}  # the broken behavior, documented
        grand = load_contigs(tmp_path, sample_name_from="grandparent")
        assert "vdj_t" not in grand


class TestSheetEdgeCases:
    def test_unmatched_file_is_skipped(self, tmp_path):
        sheet = pd.DataFrame([
            {"sample": "S1", "vdj_dir": str(tmp_path / "known")},
        ])
        _write(tmp_path / "known" / "filtered_contig.fasta", ">c1\nACGT\n")
        _write(tmp_path / "orphan" / "filtered_contig.fasta", ">c2\nACGT\n")
        contigs = load_contigs(tmp_path, sample_name_from="sheet", sample_sheet=sheet)
        assert set(contigs) == {"S1"}
        assert "c1" in contigs["S1"]

    def test_sheet_missing_required_columns_raises(self, tmp_path):
        bad = pd.DataFrame([{"sample": "S1"}])  # no vdj_dir
        _write(tmp_path / "d" / "filtered_contig.fasta")
        with pytest.raises(ValueError, match="vdj_dir"):
            load_contigs(tmp_path, sample_name_from="sheet", sample_sheet=bad)


class TestErrorPaths:
    def test_both_contig_dir_and_cellranger_dir_raises(self, tmp_path):
        with pytest.raises(ValueError, match="exactly one"):
            load_contigs(tmp_path, cellranger_dir=tmp_path)

    def test_neither_raises(self):
        with pytest.raises(ValueError, match="either contig_dir or cellranger_dir"):
            load_contigs()

    def test_unknown_sample_name_from_raises(self, tmp_path):
        with pytest.raises(ValueError, match="sample_name_from"):
            load_contigs(tmp_path, sample_name_from="bogus")

    def test_sheet_mode_requires_sample_sheet(self, tmp_path):
        with pytest.raises(ValueError, match="requires sample_sheet"):
            load_contigs(tmp_path, sample_name_from="sheet")


class TestAssembleFullSequencesThreading:
    """assemble_full_sequences exposes the same #124 controls and forwards
    them to load_contigs."""

    def test_cellranger_and_contigs_dir_mutually_exclusive(self, tmp_path):
        from tcrsift.assemble import assemble_full_sequences

        _build_cellranger_layout(tmp_path)
        with pytest.raises(ValueError, match="exactly one"):
            assemble_full_sequences(
                pd.DataFrame({"CDR3ab": ["A_B"], "CDR3_alpha": ["A"], "CDR3_beta": ["B"]}),
                contigs_dir=str(tmp_path),
                cellranger_dir=str(tmp_path),
                verbose=False,
                show_progress=False,
            )
