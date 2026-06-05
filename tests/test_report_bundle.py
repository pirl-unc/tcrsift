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

"""Tests for the `tcrsift report bundle` cohort figure PDF (#123)."""

from __future__ import annotations

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pypdf import PdfReader  # noqa: E402

from tcrsift.cli import create_parser  # noqa: E402
from tcrsift.report import CATEGORY_ORDER, _categorize, bundle_figure_pdf  # noqa: E402


def _fig(path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [1, 0])
    fig.savefig(path)
    plt.close(fig)


def _run_dir(base, name, *, fmt="pdf", stems=("qc_a", "clonotype_b", "method_c")):
    d = base / name
    plots = d / "plots"
    plots.mkdir(parents=True)
    for s in stems:
        _fig(plots / f"{s}.{fmt}")
    # An excluded non-figure PDF that must NOT be concatenated.
    _fig(plots / "tcr_sequences.pdf")
    return d


def test_categorize_orders_by_step():
    names = ["method_x.pdf", "qc_y.pdf", "clonotype_z.pdf", "weird_w.pdf"]
    cats = _categorize(names)
    order = [n for _c, n in cats]
    assert order[0] == "qc_y.pdf"  # QC first
    assert order.index("clonotype_z.pdf") < order.index("method_x.pdf")
    assert cats[-1] == ("Other", "weird_w.pdf")  # unmatched last


def test_bundle_vector_pdfs(tmp_path):
    d1 = _run_dir(tmp_path, "B1-2")
    d2 = _run_dir(tmp_path, "B1-3")
    out = tmp_path / "bundle.pdf"
    bundle_figure_pdf([d1, d2], out, title="Cohort")
    assert out.exists()
    # cover(1) + per-run [title(1) + 3 figures] x2 = 1 + 8 = 9 pages.
    assert len(PdfReader(str(out)).pages) == 9


def test_excludes_non_figure_pdfs(tmp_path):
    d1 = _run_dir(tmp_path, "B1-2")
    out = tmp_path / "bundle.pdf"
    bundle_figure_pdf([d1], out)
    # cover + title + 3 figures = 5 (tcr_sequences.pdf excluded).
    assert len(PdfReader(str(out)).pages) == 5


def test_png_fallback(tmp_path):
    d1 = _run_dir(tmp_path, "B1-4", fmt="png")
    out = tmp_path / "bundle.pdf"
    bundle_figure_pdf([d1], out)
    assert len(PdfReader(str(out)).pages) == 5  # raster pages embedded


def test_prefers_pdf_over_png_same_stem(tmp_path):
    d = tmp_path / "B1-2"
    plots = d / "plots"
    plots.mkdir(parents=True)
    _fig(plots / "qc_a.pdf")
    _fig(plots / "qc_a.png")  # same stem -> not double-counted
    out = tmp_path / "bundle.pdf"
    bundle_figure_pdf([d], out)
    assert len(PdfReader(str(out)).pages) == 3  # cover + title + 1 figure


def test_empty_raises(tmp_path):
    empty = tmp_path / "empty"
    (empty / "plots").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        bundle_figure_pdf([empty], tmp_path / "b.pdf")


def test_cli_report_bundle(tmp_path):
    d1 = _run_dir(tmp_path, "B1-2")
    d2 = _run_dir(tmp_path, "B1-3")
    out = tmp_path / "figures.pdf"
    args = create_parser().parse_args([
        "report", "bundle", str(d1), str(d2), "-o", str(out), "--labels", "Donor2,Donor3",
    ])
    args.func(args)
    assert out.exists() and len(PdfReader(str(out)).pages) == 9


def test_category_order_is_nonempty():
    assert CATEGORY_ORDER and all(isinstance(c[0], str) for c in CATEGORY_ORDER)
