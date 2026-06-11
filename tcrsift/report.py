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
# Combined per-run bundles + sequence sheets: never embed these into the cohort
# figure bundle (they'd double-include every figure). "tcrsift_report" kept for
# back-compat with pre-2.78 runs; "all-figures" is the current per-donor bundle.
_EXCLUDE_STEMS = {
    "all-figures", "tcrsift_report", "tcr_sequences", "selected_clones_sequences",
}

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


def _fit_font_size(
    text: str, font: str, start: float, floor: float, max_width: float
) -> float:
    """Largest size in ``[floor, start]`` whose ``text`` fits ``max_width`` (#combined-pdf).

    ``drawCentredString`` doesn't wrap, so a long title at a fixed size runs off
    both margins and clips. Shrink to fit instead, bottoming out at ``floor``.
    """
    from reportlab.pdfbase.pdfmetrics import stringWidth

    size = start
    while size > floor and stringWidth(text, font, size) > max_width:
        size -= 1
    return size


def _title_reader(title: str, subtitle: str = ""):
    from pypdf import PdfReader
    from reportlab.pdfgen import canvas

    from .format import pdf_safe

    w, h = _LETTER
    max_w = w - 2 * 54  # 54pt margins
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=_LETTER)

    title = pdf_safe(title)
    c.setFont("Helvetica-Bold", _fit_font_size(title, "Helvetica-Bold", 30, 10, max_w))
    c.drawCentredString(w / 2, h / 2 + 10, title)
    if subtitle:
        subtitle = pdf_safe(subtitle)
        c.setFont("Helvetica", _fit_font_size(subtitle, "Helvetica", 14, 8, max_w))
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


# Output dir names that don't identify a donor — skip up to the grandparent.
_GENERIC_DIR_NAMES = frozenset({"data", "plots", "output", "out", "results", "."})
_DONOR_COLUMNS = ("donor", "patient_id", "patient", "subject")


def resolve_report_name(
    output_dir: str | Path,
    *,
    clones_df=None,
    cli_name: str | None = None,
    default: str | None = None,
) -> str:
    """Resolve a donor/run label for a report cover + filename (#262 follow-up).

    Cascade: explicit ``cli_name`` → a unanimous donor/patient field in
    ``clones_df`` → the output dir name (parent, or grandparent when the parent
    is a generic name like ``data``/``plots``) → ``default`` → the literal
    ``"run"``. So per-donor runs get their donor name automatically while an
    explicit ``--report-name`` always wins.
    """
    if cli_name and str(cli_name).strip():
        return str(cli_name).strip()
    if clones_df is not None and len(getattr(clones_df, "columns", [])):
        for col in _DONOR_COLUMNS:
            if col in clones_df.columns:
                vals = clones_df[col].dropna().astype(str).unique()
                if len(vals) == 1 and vals[0].strip():
                    return vals[0].strip()
    p = Path(output_dir).resolve()
    for cand in (p, p.parent):
        if cand.name and cand.name not in _GENERIC_DIR_NAMES:
            return cand.name
    return (default or "run").strip() or "run"


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
# ``alpha_contig_ids`` is NOT refreshed from the rep cell here — it's set
# per-variant from the co-expressing cells' contigs (see expand_dual_alpha_
# variants). The merged clone's list only carries the DOMINANT α's contigs, so
# without a per-variant override the non-dominant α variant searches the wrong
# contigs in _extract_c_region_nt_from_contig → no_contig → blanket-N fallback,
# even when that α is fully sequenced through the C region (#235 dual-α footgun).


def _rep_cell_both_alphas(obs, beta: str, a_lo: str, a_hi: str):
    """Return ``(alpha_cdr3 -> slot_prefix, cell_row, alpha_contigs)`` for the
    highest-UMI cell that genuinely co-expresses both alphas with this beta,
    else None.

    ``alpha_contigs`` maps each α CDR3 to a ``;``-joined list of EVERY contig id
    for that α across ALL co-expressing cells (not just the rep cell) — so the
    per-α C-region extraction keeps its multi-contig consensus and maximal
    coverage, mirroring how the merged clone's contig list is built. The rep
    cell still supplies the VDJ/gene bytes (one coherent source)."""
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
    # Collect each α's contig ids across all co-expressing cells (the α may sit
    # in TRA_1 in one cell, TRA_2 in another), de-duplicated, order-preserving.
    s_cells = s.loc[cells.index]
    alpha_contigs: dict[str, str] = {}
    for a in (a_lo, a_hi):
        ids: list[str] = []
        for slot in ("TRA_1", "TRA_2"):
            cid_col = f"{slot}_contig_id"
            if cid_col not in cells.columns:
                continue
            hit = cells[cid_col][(s_cells[f"{slot}_cdr3"] == a).values]
            for cid in hit.dropna().astype(str):
                if cid and cid.lower() != "nan" and cid not in ids:
                    ids.append(cid)
        alpha_contigs[a] = ";".join(ids)
    return {str(cell["TRA_1_cdr3"]): "TRA_1", str(cell["TRA_2_cdr3"]): "TRA_2"}, cell, alpha_contigs


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
            amap, cell, alpha_contigs = rc
            for a in alphas:
                prefix = amap.get(a)
                if prefix is None:
                    continue
                vr = r.copy()
                for col, suf in _ALPHA_COL_TO_CELL.items():
                    cellcol = f"{prefix}_{suf}"
                    if cellcol in cell.index:
                        vr[col] = cell[cellcol]
                # Point this variant at ITS OWN α's contigs (across the
                # co-expressing cells), never the merged clone's dominant-α
                # list. None when no contig id is available, so extraction
                # falls back cleanly rather than searching the wrong α (#235).
                vr["alpha_contig_ids"] = alpha_contigs.get(a) or None
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


# Construct legend / cover page (#202). Entries are (kind, text):
#   "h1"   section heading (bold)
#   "p"    body line
#   "mono" fixed-width line (the architecture diagram / column names)
#   ""     blank spacer
# Rendered by _legend_reader; kept structured so the layout stays tidy and the
# text is easy to edit. ASCII-only (beta/alpha spelled out, plain +/-) so it is
# safe in the base-14 PDF fonts — see format.pdf_safe (#202).
_CONSTRUCT_LEGEND = [
    ("p", "Each page is one synthesis-ready single-chain construct."),
    ("p", "Segments run N -> C in this order:"),
    ("", ""),
    ("mono", "  [beta leader]-VDJ_beta-[beta C]-[T2A]-[alpha leader]-VDJ_alpha-[alpha C]"),
    ("", ""),
    ("h1", "Architecture: beta-T2A-alpha"),
    ("p", "T2A is a self-cleaving 2A 'ribosomal-skip' peptide: the ribosome"),
    ("p", "releases the beta chain and re-initiates the alpha, so one transcript"),
    ("p", "yields two separate chains in ~1:1 ratio."),
    ("", ""),
    ("h1", "Leaders (signal peptides)"),
    ("p", "Default CD8A (beta) / CD28 (alpha). Configurable."),
    ("", ""),
    ("h1", "Constants"),
    ("p", "J-family parity sets the TRBC allele (TRBJ1 -> TRBC1, TRBJ2 -> TRBC2);"),
    ("p", "alpha uses TRAC. With contigs the donor's observed allele is verified,"),
    ("p", "and any divergence from canonical is flagged (constant_allele_divergence)."),
    ("", ""),
    ("h1", "Sequence forms (per construct, in selected_clones.csv)"),
    ("mono", "  full_*_nt        unoptimized - observed / native codons"),
    ("mono", "  *_nt_optimized   codon-optimized for expression"),
    ("", ""),
    ("h1", "Dual-alpha (allelic inclusion)"),
    ("p", "Emitted as TWO labeled variants of one clone (shared beta, one per"),
    ("p", "alpha). Pick one to synthesize."),
    ("", ""),
    ("h1", "On each page"),
    ("p", "The header shows selection provenance (route, rank, frequency); the"),
    ("p", "color key maps the V / CDR3 / J / constant / leader / linker segments."),
]


def _legend_reader(title: str, entries: list):
    """A reportlab text page (PdfReader) — the construct cover (#202).

    ``entries`` is a list of ``(kind, text)`` tuples (see _CONSTRUCT_LEGEND) or
    plain strings (rendered as body lines). All text is routed through
    ``pdf_safe`` so non-WinAnsi glyphs don't render as boxes (#202).
    """
    from pypdf import PdfReader
    from reportlab.pdfgen import canvas

    from .format import pdf_safe

    w, h = _LETTER
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=_LETTER)
    c.setFont("Helvetica-Bold", 20)
    c.drawString(54, h - 70, pdf_safe(title))
    y = h - 108
    for entry in entries:
        kind, text = entry if isinstance(entry, tuple) else ("p", entry)
        if not text:
            y -= 9
            continue
        if kind == "h1":
            y -= 4
            c.setFont("Helvetica-Bold", 12)
            c.drawString(54, y, pdf_safe(text))
        elif kind == "mono":
            c.setFont("Courier", 9.5)
            c.drawString(60, y, pdf_safe(text))
        else:
            c.setFont("Helvetica", 10.5)
            c.drawString(60, y, pdf_safe(text))
        y -= 15
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf)


def _is_legend_page(page) -> bool:
    """True when a PDF page is a construct-legend cover (#combined-pdf).

    Detected by prose unique to the legend body (``ribosomal``-skip etc.) that
    never appears on a construct sequence page — so a cover=False donor PDF or a
    real construct page is never mistaken for the legend and stripped.
    """
    try:
        text = page.extract_text() or ""
    except Exception:
        return False
    return "ribosomal" in text and "T2A" in text


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
        reader = PdfReader(str(pdf))
        pages = list(reader.pages)
        # Each donor's selected_clones_sequences.pdf already begins with its own
        # construct-legend cover. When the cohort legend is shown once up front,
        # that per-donor cover is redundant — drop it so it doesn't repeat before
        # every section (#combined-pdf).
        if include_legend and pages and _is_legend_page(pages[0]):
            pages = pages[1:]
        for page in pages:
            writer.add_page(page)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as fh:
        writer.write(fh)
    logger.info("Wrote combined selected report %s (%d pages)", out_path, len(writer.pages))
    return out_path


# Per-condition selection breakdown columns (#selection-cols) — emitted into
# selected_clones.csv but kept OUT of the sequence PDF, where they'd just
# duplicate the formatted ``selection_detail`` block and clutter the page.
_SELECTION_BREAKDOWN_COLS = frozenset({
    "selection_conditions",
    "frequency_per_condition",
    "prism_per_condition",
    "freq_rank_per_condition",
    "prism_rank_per_condition",
})


def _format_selection_item(item: str) -> str:
    """Format one compact selection token for the sequence-PDF footer.

    ``"AIM+=freq#6(0.90%)"`` -> ``"AIM+ freq#6 (0.90%)"`` — a space after the
    condition (no ``=``) and a space before the parenthesised frequency. Tokens
    without an ``=`` are returned unchanged.
    """
    item = item.strip()
    if not item:
        return ""
    cond, sep, rest = item.partition("=")
    if not sep:
        return item
    rest = rest.replace("(", " (", 1)
    return f"{cond} {rest}"


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
    allow_canonical_fallback: bool = False,
    annotations=None,
):
    """Assemble a selected clone set into a synthesis-ready deliverable (#188).

    Joins the selection (clone ids + provenance) to ``clonotypes`` (the full
    VDJ/genes, now #184-correct), emits both alpha-variants for merged dual-alpha
    clones (when ``obs`` is given), assembles β-T2A-α constructs, runs the
    assembly QC, and writes ``selected_clones.csv`` (sequences + provenance),
    ``selected_clones_sequences.pdf``, and ``selected_clones_qc.txt`` under
    ``output_dir``. Returns the assembled DataFrame.
    """
    from .assemble import (
        assemble_full_sequences,
        assemble_qc_report,
        synthesis_qc_report,
        validate_sequences,
    )
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
    # Fail closed on canonical-fallback constructs unless explicitly allowed
    # (#241/#243/#244): this is the command that shipped the all-canonical
    # deliverable when --cellranger-dir was omitted. The QC text is written
    # below regardless; this gate refuses to call the output synthesis-ready
    # when junctions/alleles weren't contig-verified.
    from .assemble import enforce_contig_fidelity

    fidelity = enforce_contig_fidelity(
        assembled,
        allow_canonical_fallback=allow_canonical_fallback,
        context="report selected",
    )
    if fidelity:
        logger.warning(fidelity)
    validate_sequences(assembled, strict=False)
    qc_text = assemble_qc_report(assembled)
    synth_text = synthesis_qc_report(assembled)
    if synth_text:
        qc_text = f"{qc_text}\n\n{synth_text}"
    # Dual-α rollup (#237): make the clone→construct expansion explicit so it's
    # impossible to overlook that some clones carry two α and become 2 constructs.
    n_dual = len(set(variant_of.values())) if variant_of else 0
    if n_dual:
        dual_line = (
            f"Dual-α clones: {n_dual} (each carries two α → 2 single-chain "
            f"constructs; {n_variants} dual-α construct rows total)."
        )
        logger.info(dual_line)
        qc_text = f"{dual_line}\n\n{qc_text}"
    (out_dir / "selected_clones_qc.txt").write_text(qc_text + "\n")

    # Leader QC outputs (#275): emit from the select-then-assemble path too, not
    # only the monolithic `run` assembly — the assembled shortlist carries the
    # leader/germline columns these need.
    from .assemble import collect_germline_variants

    _gv = collect_germline_variants(assembled)
    if not _gv.empty:
        _gv.to_csv(out_dir / "germline_variants.csv", index=False)
        logger.info("  Germline variants: %d distinct (by region) → germline_variants.csv", len(_gv))
    try:
        from .plots import plot_leader_summary

        plot_leader_summary(assembled, out_dir / "leader_summary.png")
    except Exception as e:  # plotting is best-effort; don't fail the report
        logger.warning("leader_summary plot failed: %s", e, exc_info=True)

    # Per-clone provenance lines for the sequence PDF. NB: keep this distinct
    # from the ``annotations`` parameter (the public-DB frame) — they collided in
    # 2.92.0 and silently emptied every db_* column (#296).
    pdf_annotations = None
    if prov_cols:
        # Stable variant numbering: within each parent clone, number its
        # dual-alpha variants #1, #2, … by sorted variant id, so the PDF label
        # is reproducible (#188 emits them in partner-list order, not sorted).
        variant_rank: dict[str, int] = {}
        if "dual_alpha_variant" in expanded.columns and "selected_clone" in expanded.columns:
            dav = expanded[expanded["dual_alpha_variant"].notna()]
            for _parent, grp in dav.groupby("selected_clone", sort=True):
                for i, cdr in enumerate(sorted(grp["CDR3ab"].astype(str)), 1):
                    variant_rank[cdr] = i
        pdf_annotations = {}
        for _, r in expanded.iterrows():
            key = r["CDR3ab"]
            lines: list[str] = []
            for c in prov_cols:
                if not pd.notna(r.get(c)):
                    continue
                cl = c.lower()
                if cl in ("selection", "selection_detail"):
                    # One condition per indented line, reformatted; no sub-label
                    # (the bold "Selection:" header already titles the block).
                    # ``selection`` is the pre-2.90 name; ``selection_detail`` the
                    # current one — accept both.
                    for item in str(r[c]).split(";"):
                        formatted = _format_selection_item(item)
                        if formatted:
                            lines.append(f"    {formatted}")
                elif cl in _SELECTION_BREAKDOWN_COLS:
                    # Per-condition freq/PRISM/rank breakdowns are CSV-only — they
                    # duplicate selection_detail and would clutter the page.
                    continue
                else:
                    lines.append(f"{c}: {r[c]}")
            if r.get("dual_alpha_variant"):
                n = variant_rank.get(str(key))
                parent = r.get("selected_clone")
                lines.append(
                    f"dual-alpha variant #{n} of {parent}" if n
                    else f"dual-alpha variant of {parent}"
                )
            if lines:
                pdf_annotations[key] = lines

    seq_pdf = out_dir / "selected_clones_sequences.pdf"
    create_tcr_sequence_pdf(
        assembled, seq_pdf, strict=False, annotations=pdf_annotations,
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

    # Public-DB (IEDB/CEDAR/VDJdb) match annotation onto the selected clones
    # (#selected-anno). report selected joins to the un-annotated clonotypes, so
    # without this the deliverable can't show which TCR-T candidates match a
    # viral / public epitope. Curated columns only; dual-α variants inherit their
    # parent clone's annotation. No-op (columns still emitted, empty) when no
    # annotations frame is supplied.
    assembled = _attach_public_db_annotation(assembled, annotations, variant_of)

    # Per-donor selected-clones × conditions frequency table + heatmap
    # (#selected-freq): rows = selected CDR3ab, columns = conditions, cells =
    # within-condition frequency. Parsed from the frequency_per_condition
    # provenance; skipped (with a note) when that column is absent.
    try:
        _emit_frequency_by_condition(assembled, out_dir, name=title)
    except Exception as e:  # best-effort; never fail the build on the side table
        logger.warning("selected frequency table/heatmap failed: %s", e, exc_info=True)

    assembled.to_csv(out_dir / "selected_clones.csv", index=False)

    # Independent QC battery (#279). Re-derives every construct fact from first
    # principles and cross-checks the assembled frame against the source
    # clonotypes and the raw cellranger contigs — NOT trusting tcrsift's own QC
    # columns. Writes a PASS/FAIL/SKIP `qc-summary.md` so every deliverable
    # self-validates. Checks whose inputs are absent (no contig dir) report SKIP.
    try:
        _emit_qc_summary(
            assembled,
            out_dir,
            clonotypes=clonotypes,
            assemble_kwargs=assemble_kwargs,
            cohort_name=title,
        )
    except Exception as e:  # the battery is a gate, not a hard failure of the build
        logger.warning("qc-summary battery failed: %s", e, exc_info=True)

    logger.info(
        "Selected-clones report: %d constructs (%d dual-alpha variants) → %s",
        len(assembled), n_variants, out_dir,
    )
    return assembled


# Curated public-DB (IEDB/CEDAR/VDJdb) annotation columns surfaced on the
# selected clones so the deliverable shows which TCR-T candidates match a viral /
# public epitope (#selected-anno).
_PUBLIC_DB_COLS = [
    "is_viral", "db_match", "db_match_strength", "db_category",
    "db_epitope", "db_protein_canonical", "db_species", "db_database",
]


def _attach_public_db_annotation(assembled, annotations, variant_of):
    """Merge curated public-DB match columns onto the assembled selected clones.

    Keyed by ``CDR3ab``; a dual-α variant (whose CDR3ab differs from its parent)
    inherits the parent clone's annotation via ``variant_of``. Always leaves the
    curated columns present (NaN when no annotations / no match), so the schema is
    stable whether or not an annotations frame was supplied.
    """
    import pandas as pd

    have = (
        annotations is not None
        and not getattr(annotations, "empty", True)
        and "CDR3ab" in getattr(annotations, "columns", [])
    )
    if have:
        present = [c for c in _PUBLIC_DB_COLS if c in annotations.columns]
        src = annotations[["CDR3ab", *present]].copy()
        # Deterministic dedup: real annotation frames are one row per CDR3ab, but
        # if a multi-match file is passed, keep the STRONGEST match per clone
        # (a real db_match / viral flag over an empty row) reproducibly, rather
        # than an arbitrary first row.
        sort_keys = [c for c in ("db_match", "is_viral") if c in present]
        if sort_keys:
            src = src.sort_values(
                ["CDR3ab", *sort_keys],
                ascending=[True, *([False] * len(sort_keys))],
                na_position="last", kind="stable",
            )
        anno = src.drop_duplicates("CDR3ab", keep="first").set_index("CDR3ab")

        def _key(cdr):
            if cdr in anno.index:
                return cdr
            parent = variant_of.get(cdr)
            return parent if parent in anno.index else None

        keys = [_key(str(c)) for c in assembled["CDR3ab"].astype(str)]
        for c in present:
            col_map = anno[c].to_dict()
            assembled[c] = [col_map.get(k, pd.NA) for k in keys]
    for c in _PUBLIC_DB_COLS:
        if c not in assembled.columns:
            assembled[c] = pd.NA
    return assembled


def _parse_freq_per_condition(s: str) -> dict:
    """``"AIM⁺=0.90%;CTY⁻=0.34%"`` → ``{"AIM⁺": 0.90, "CTY⁻": 0.34}`` (percent)."""
    out: dict = {}
    for item in str(s).split(";"):
        item = item.strip()
        if "=" not in item:
            continue
        cond, _, val = item.rpartition("=")
        val = val.strip().rstrip("%")
        try:
            out[cond.strip()] = float(val)
        except ValueError:
            pass
    return out


def _emit_frequency_by_condition(assembled, out_dir, *, name: str = "selected"):
    """Write the per-donor selected-clones × conditions frequency table + heatmap.

    Rows = selected CDR3ab, columns = conditions, cells = within-condition
    frequency (%), parsed from the ``frequency_per_condition`` provenance. Skips
    (returns None) when that column is absent or carries no parseable values.
    """
    import pandas as pd

    if "frequency_per_condition" not in assembled.columns:
        return None
    rows: dict = {}
    for _, r in assembled.iterrows():
        cdr = r.get("CDR3ab")
        fpc = r.get("frequency_per_condition")
        if isinstance(fpc, str) and fpc:
            parsed = _parse_freq_per_condition(fpc)
            if parsed:
                rows[cdr] = parsed
    if not rows:
        return None
    pivot = pd.DataFrame.from_dict(rows, orient="index").sort_index(axis=1)
    pivot.index.name = "CDR3ab"
    pivot.to_csv(out_dir / "selected_frequency_by_condition.csv")
    logger.info(
        "  Selected freq × condition: %d clones × %d conditions → "
        "selected_frequency_by_condition.csv",
        pivot.shape[0], pivot.shape[1],
    )
    # Public-DB category per clone for the heatmap's annotation strip (#selected-anno):
    # so viral/public TCR-T candidates are readable off the frequency heatmap.
    anno = None
    if "db_category" in assembled.columns and "CDR3ab" in assembled.columns:
        anno = (
            assembled[["CDR3ab", "db_category"]]
            .drop_duplicates("CDR3ab")
            .set_index("CDR3ab")["db_category"]
        )
    try:
        from .plots import plot_selected_frequency_heatmap

        plot_selected_frequency_heatmap(
            pivot, out_dir / "selected_frequency_heatmap.png", title=name,
            annotations=anno,
        )
    except Exception as e:  # plot is best-effort
        logger.warning("selected frequency heatmap failed: %s", e, exc_info=True)
    return pivot


def _emit_qc_summary(
    assembled,
    out_dir: Path,
    *,
    clonotypes=None,
    assemble_kwargs: dict | None = None,
    cohort_name: str = "deliverable",
) -> bool:
    """Run the #279 QC battery on an assembled frame and write ``qc-summary.md``.

    Returns ``True`` when every check passed (SKIP does not fail). Logs a loud
    warning when any check FAILed so a real integrity regression isn't silent.
    """
    from .qc_battery import contig_cdr3_map, qc_summary_markdown, run_qc_battery

    ak = assemble_kwargs or {}
    cmap = contig_cdr3_map(
        contig_dir=ak.get("contigs_dir"),
        cellranger_dir=ak.get("cellranger_dir"),
    )
    # Prefer the resolved donor/dir name as the cohort header (consistent with
    # the rest of the reporting); fall back to the passed title.
    try:
        cohort_name = resolve_report_name(
            out_dir, clones_df=assembled, default=cohort_name
        )
    except Exception:
        pass
    results = run_qc_battery(assembled, clonotypes=clonotypes, contig_map=cmap)
    markdown, all_ok = qc_summary_markdown(
        {cohort_name: results}, deliverable=str(out_dir)
    )
    (out_dir / "qc-summary.md").write_text(markdown)
    failed = [r.name for r in results if r.status == "FAIL"]
    if failed:
        logger.warning(
            "  QC battery FAILED %d/%d checks (%s) → see qc-summary.md",
            len(failed), len(results), ", ".join(failed),
        )
    else:
        logger.info(
            "  QC battery: %d/%d checks pass → qc-summary.md",
            sum(1 for r in results if r.status == "PASS"), len(results),
        )
    return all_ok
