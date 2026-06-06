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

import pandas as pd

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


# Clonotype column -> cell per-chain column suffix, for rebuilding an alpha row
# from a specific cell+slot (#188 dual-alpha variants / #184 repair).
_ALPHA_COL_TO_CELL = {
    "CDR3_alpha": "cdr3", "alpha_v_gene": "v_gene", "alpha_j_gene": "j_gene",
    "alpha_c_gene": "c_gene", "VDJ_alpha_aa": "vdj_aa", "VDJ_alpha_nt": "vdj_nt",
}


def _rep_cell_both_alphas(obs, beta: str, a_lo: str, a_hi: str):
    """Return ``(alpha_cdr3 -> slot_prefix, cell_row)`` for the highest-UMI cell
    that genuinely co-expresses both alphas with this beta, else None."""
    s = obs.astype({c: str for c in ("TRA_1_cdr3", "TRA_2_cdr3", "TRB_1_cdr3") if c in obs.columns})
    if not {"TRA_1_cdr3", "TRA_2_cdr3", "TRB_1_cdr3"} <= set(obs.columns):
        return None
    m = (s["TRB_1_cdr3"] == beta) & (
        ((s["TRA_1_cdr3"] == a_lo) & (s["TRA_2_cdr3"] == a_hi))
        | ((s["TRA_1_cdr3"] == a_hi) & (s["TRA_2_cdr3"] == a_lo))
    )
    cells = obs[m.values]
    if not len(cells):
        return None
    tot = sum(
        cells[c].fillna(0) for c in ("TRA_1_umis", "TRA_2_umis", "TRB_1_umis")
        if c in cells.columns
    )
    cell = cells.loc[tot.idxmax()]
    return {str(cell["TRA_1_cdr3"]): "TRA_1", str(cell["TRA_2_cdr3"]): "TRA_2"}, cell


def expand_dual_alpha_variants(
    selected, obs, *, partners_col: str = "merged_alpha_partners"
):
    """Emit BOTH alpha-variants of each merged dual-alpha clone (#188).

    A merged dual-alpha clone carries two alphas (``merged_alpha_partners``) but
    a single-chain construct holds one. For each such clone, find a cell that
    co-expresses both alphas + the beta and emit two rows — one per alpha, with
    that alpha's VDJ/genes pulled from the cell's matching slot so each
    construct's CDR3α matches its own VDJ. Non-merged clones pass through.

    Returns ``(expanded_df, variant_of)`` where ``variant_of`` maps each variant
    clone id to its canonical selected clone. A no-op (returns a copy) when
    ``obs`` is None or there's no ``merged_alpha_partners`` column.
    """
    out = selected.copy()
    if partners_col not in out.columns or obs is None:
        out["dual_alpha_variant"] = None
        return out.reset_index(drop=True), {}

    rows = []
    variant_of: dict[str, str] = {}
    for _, r in out.iterrows():
        partners = str(r.get(partners_col) or "")
        alphas = [a for a in partners.split(";") if a]
        beta = str(r.get("CDR3_beta", "") or "")
        rc = _rep_cell_both_alphas(obs, beta, *sorted(alphas)) if len(alphas) == 2 else None
        if rc:
            amap, cell = rc
            for a in alphas:
                prefix = amap.get(a)
                if prefix is None:
                    continue
                vr = r.copy()
                for col, suf in _ALPHA_COL_TO_CELL.items():
                    cellcol = f"{prefix}_{suf}"
                    if cellcol in cell.index:
                        vr[col] = cell[cellcol]
                vr["CDR3ab"] = f"{a}_{beta}"
                vr["dual_alpha_variant"] = a
                vr["selected_clone"] = r["CDR3ab"]
                variant_of[vr["CDR3ab"]] = r["CDR3ab"]
                rows.append(vr)
        else:
            rr = r.copy()
            rr["dual_alpha_variant"] = None
            rr["selected_clone"] = r["CDR3ab"]
            rows.append(rr)
    return pd.DataFrame(rows).reset_index(drop=True), variant_of


_CONSTRUCT_LEGEND = [
    "Each page is one synthesis-ready single-chain construct, in N->C order:",
    "",
    "   [β leader] — VDJ-β — [β constant] — [T2A] — [α leader] — VDJ-α — [α constant]",
    "",
    "• Architecture: beta-T2A-alpha. The T2A is a 2A 'ribosomal skip' peptide —",
    "  the ribosome releases the beta chain and re-initiates the alpha, so one",
    "  transcript yields two separate chains in ~1:1 ratio.",
    "• Leaders (signal peptides): default CD8A (beta) / CD28 (alpha) — configurable.",
    "• Constants: J-family parity sets the TRBC allele (TRBJ1->TRBC1, TRBJ2->TRBC2);",
    "  alpha uses TRAC. With contigs, the donor's observed allele is verified and any",
    "  divergence from canonical is flagged (constant_allele_divergence).",
    "• Two CDS forms per construct in selected_clones.csv: full_*_nt (unoptimized,",
    "  the observed/native codons) and *_optimized (codon-optimized for expression).",
    "• Dual-alpha (allelic-inclusion) clones are emitted as TWO labeled construct",
    "  variants of the same clone (same beta, each alpha) — pick one to synthesize.",
    "",
    "Per-clone header shows selection provenance (route, rank, frequency); the",
    "color legend on each page maps the V / CDR3 / J / constant / leader / linker segments.",
]


def _legend_reader(title: str, lines: list[str]):
    """A reportlab text page (PdfReader) — used for the construct cover (#202)."""
    from pypdf import PdfReader
    from reportlab.pdfgen import canvas

    w, h = _LETTER
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=_LETTER)
    c.setFont("Helvetica-Bold", 20)
    c.drawString(54, h - 70, title)
    c.setFont("Helvetica", 10.5)
    y = h - 110
    for line in lines:
        c.drawString(54, y, line)
        y -= 16
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf)


def combine_selected_pdfs(
    pdfs: list[str | Path],
    out_path: str | Path,
    *,
    labels: list[str] | None = None,
    title: str = "Selected clones — cohort",
    include_legend: bool = True,
) -> Path:
    """Combine per-donor ``selected_clones_sequences.pdf`` into one cohort PDF (#202).

    Cover (optional construct legend) + per-donor section title + that donor's
    pages, concatenated as vector. Mirrors ``report bundle`` but for the
    selected-clone reports. Returns the output path.
    """
    from pypdf import PdfReader, PdfWriter

    pdfs = [Path(p) for p in pdfs]
    labels = labels or [p.parent.name or p.stem for p in pdfs]
    if len(labels) != len(pdfs):
        raise ValueError("labels must match pdfs length")

    writer = PdfWriter()

    def _add(reader):
        for page in reader.pages:
            writer.add_page(page)

    _add(_title_reader(title, f"{len(pdfs)} donors"))
    if include_legend:
        _add(_legend_reader("Construct legend (β-T2A-α)", _CONSTRUCT_LEGEND))
    for pdf, label in zip(pdfs, labels):
        if not pdf.is_file():
            logger.warning("  selected-combine: missing %s — skipping", pdf)
            continue
        _add(_title_reader(label))
        _add(PdfReader(str(pdf)))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as fh:
        writer.write(fh)
    logger.info("Wrote combined selected report %s (%d pages)", out_path, len(writer.pages))
    return out_path


def build_selected_report(
    selected,
    clonotypes,
    output_dir,
    *,
    obs=None,
    assemble_kwargs: dict | None = None,
    provenance_cols: list[str] | None = None,
    title: str = "Selected clones",
    cover: bool = True,
):
    """Assemble a selected clone set into a synthesis-ready deliverable (#188).

    Joins the selection (clone ids + provenance) to ``clonotypes`` (the full
    VDJ/genes, now #184-correct), emits both alpha-variants for merged dual-alpha
    clones (when ``obs`` is given), assembles β-T2A-α constructs, runs the
    assembly QC, and writes ``selected_clones.csv`` (sequences + provenance),
    ``selected_clones_sequences.pdf``, and ``selected_clones_qc.txt`` under
    ``output_dir``. Returns the assembled DataFrame.
    """
    from .assemble import assemble_full_sequences, assemble_qc_report, validate_sequences
    from .plots import create_tcr_sequence_pdf

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prov_cols = [c for c in (provenance_cols or []) if c in selected.columns]
    keep = ["CDR3ab", *prov_cols]
    sel = selected[[c for c in keep if c in selected.columns]].drop_duplicates("CDR3ab")
    merged = sel.merge(clonotypes, on="CDR3ab", how="left", suffixes=("", "_clo"))

    expanded, variant_of = expand_dual_alpha_variants(merged, obs)
    n_variants = len(variant_of)

    assembled = assemble_full_sequences(
        expanded, verbose=False, show_progress=False, **(assemble_kwargs or {})
    )
    validate_sequences(assembled, strict=False)
    qc_text = assemble_qc_report(assembled)
    (out_dir / "selected_clones_qc.txt").write_text(qc_text + "\n")

    # Per-clone provenance lines for the sequence PDF.
    annotations = None
    if prov_cols:
        annotations = {}
        for _, r in expanded.iterrows():
            key = r["CDR3ab"]
            lines = [f"{c}: {r[c]}" for c in prov_cols if pd.notna(r.get(c))]
            if r.get("dual_alpha_variant"):
                lines.append(f"dual-alpha variant of {r.get('selected_clone')}")
            if lines:
                annotations[key] = lines

    seq_pdf = out_dir / "selected_clones_sequences.pdf"
    create_tcr_sequence_pdf(
        assembled, seq_pdf, strict=False, annotations=annotations,
    )
    # Prepend a self-documenting construct cover/legend page (#202).
    if cover:
        from pypdf import PdfReader, PdfWriter

        writer = PdfWriter()
        for page in _legend_reader(f"{title} — β-T2A-α constructs", _CONSTRUCT_LEGEND).pages:
            writer.add_page(page)
        for page in PdfReader(str(seq_pdf)).pages:
            writer.add_page(page)
        with seq_pdf.open("wb") as fh:
            writer.write(fh)

    assembled.to_csv(out_dir / "selected_clones.csv", index=False)
    logger.info(
        "Selected-clones report: %d constructs (%d dual-alpha variants) → %s",
        len(assembled), n_variants, out_dir,
    )
    return assembled
