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

"""Unified reporting: assemble a cohort figure bundle (#123).

Folds the hand-written ``make_figure_set.py`` from downstream analyses into the
tool: concatenate each run's per-figure plots into one categorized,
multi-sample PDF with section/title pages. Vector figures (``plot_format: pdf``)
are concatenated as true vector pages so text stays selectable; PNG-only runs
fall back to embedded raster pages.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Category label -> filename-stem prefixes, in report order. Mirrors the
# pipeline step order (QC -> phenotype -> clonotypes -> filtering -> annotation
# -> assembly -> method recovery); anything unmatched lands under "Other".
CATEGORY_ORDER: list[tuple[str, tuple[str, ...]]] = [
    ("QC", ("qc_",)),
    ("Phenotype", ("phenotype_",)),
    ("Clonotypes", ("clonotype_",)),
    ("Filtering & selection", ("filter_", "funnel_")),
    ("Annotation (public DBs)", ("annotation_",)),
    ("Assembly", ("assembly_",)),
    ("Method recovery / overlap", ("method_",)),
]

# Per-run PDFs that aren't figures (the report itself, the sequence sheet).
_EXCLUDE_STEMS = {"tcrsift_report", "tcr_sequences", "selected_clones_sequences"}

_LETTER = (612.0, 792.0)  # reportlab letter, in points


def _categorize(names: list[str]) -> list[tuple[str, str]]:
    """Order figure filenames by category; unmatched go last under 'Other'."""
    used: set[str] = set()
    ordered: list[tuple[str, str]] = []
    for cat, prefixes in CATEGORY_ORDER:
        for n in names:
            if n in used:
                continue
            if any(Path(n).name.startswith(p) for p in prefixes):
                used.add(n)
                ordered.append((cat, n))
    for n in names:
        if n not in used:
            ordered.append(("Other", n))
    return ordered


def _pretty(stem: str) -> str:
    return Path(stem).stem.replace("_", " ")


def _title_reader(title: str, subtitle: str = ""):
    from pypdf import PdfReader
    from reportlab.pdfgen import canvas

    w, h = _LETTER
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=_LETTER)
    c.setFont("Helvetica-Bold", 30)
    c.drawCentredString(w / 2, h / 2 + 10, title)
    if subtitle:
        c.setFont("Helvetica", 14)
        c.setFillGray(0.35)
        c.drawCentredString(w / 2, h / 2 - 25, subtitle)
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf)


def _image_reader(img_path: Path, header: str):
    from pypdf import PdfReader
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfgen import canvas

    w, h = _LETTER
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=_LETTER)
    c.setFont("Helvetica-Bold", 12)
    c.drawString(40, h - 40, header)
    im = ImageReader(str(img_path))
    iw, ih = im.getSize()
    scale = min((w - 80) / iw, (h - 110) / ih)
    c.drawImage(im, (w - iw * scale) / 2, (h - 90 - ih * scale) / 2, iw * scale, ih * scale)
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf)


def _find_plots_dir(run_dir: Path) -> Path:
    """A run's figures live in ``<run_dir>/plots`` (the run layout), but accept
    a plots dir passed directly too."""
    plots = run_dir / "plots"
    return plots if plots.is_dir() else run_dir


def bundle_figure_pdf(
    run_dirs: list[str | Path],
    out_path: str | Path,
    *,
    labels: list[str] | None = None,
    title: str = "TCRsift figure set",
    subtitle: str = "",
    extra_images: list[tuple[str, str | Path]] | None = None,
) -> Path:
    """Concatenate per-run figures into one categorized cohort PDF (#123).

    For each run directory, figures from ``<run_dir>/plots`` are ordered by
    :data:`CATEGORY_ORDER` behind a per-run title page. Vector ``.pdf`` figures
    are preferred (concatenated as vector pages); a run with only ``.png``
    figures falls back to embedded raster pages. ``extra_images`` (``(header,
    path)`` pairs) are appended as raster pages — for cross-run montages that
    have no vector source.

    Returns the output path. Raises FileNotFoundError if no figures are found.
    """
    from pypdf import PdfWriter

    run_dirs = [Path(d) for d in run_dirs]
    labels = labels or [d.name for d in run_dirs]
    if len(labels) != len(run_dirs):
        raise ValueError("labels must match run_dirs length")

    writer = PdfWriter()

    def _add(reader):
        for page in reader.pages:
            writer.add_page(page)

    _add(_title_reader(title, subtitle))

    total_figures = 0
    for run_dir, label in zip(run_dirs, labels):
        plots = _find_plots_dir(run_dir)
        pdfs = sorted(
            p.name for p in plots.glob("*.pdf") if p.stem not in _EXCLUDE_STEMS
        )
        pngs = sorted(
            p.name for p in plots.glob("*.png") if p.stem not in _EXCLUDE_STEMS
        )
        # Prefer vector; only use a PNG when there's no PDF of the same stem.
        pdf_stems = {Path(n).stem for n in pdfs}
        png_only = [n for n in pngs if Path(n).stem not in pdf_stems]
        names = pdfs + png_only
        if not names:
            logger.warning("  report bundle: no figures in %s — skipping", plots)
            continue
        _add(_title_reader(label, f"{len(names)} figures"))
        for _cat, name in _categorize(names):
            fpath = plots / name
            if fpath.suffix == ".pdf":
                from pypdf import PdfReader

                _add(PdfReader(str(fpath)))
            else:
                _add(_image_reader(fpath, f"{label} · {_pretty(name)}  [raster]"))
            total_figures += 1

    for header, img in extra_images or []:
        img = Path(img)
        if img.is_file():
            _add(_image_reader(img, header))

    if total_figures == 0 and not (extra_images or []):
        raise FileNotFoundError(
            f"No figures found under {[str(d) for d in run_dirs]}. Run the "
            "pipeline with generate_plots first (plot_format: pdf for vector)."
        )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as fh:
        writer.write(fh)
    logger.info("Wrote figure bundle %s (%d pages)", out_path, len(writer.pages))
    return out_path
