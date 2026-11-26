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

"""Reporting CLI subcommands: ``report bundle|cohort-figures|selected|...``."""

from __future__ import annotations

import argparse
from pathlib import Path

from .common import setup_logging


def cmd_report_bundle(args):
    """Assemble a categorized, multi-run figure bundle PDF (#123)."""
    from ..report import bundle_figure_pdf

    setup_logging(args.verbose)
    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    out = bundle_figure_pdf(
        args.run_dirs,
        args.output,
        labels=labels,
        title=args.title,
        subtitle=args.subtitle or "",
    )
    print(f"Wrote figure bundle: {out}")


def cmd_report_cohort_figures(args):
    """Cross-donor all-figures bundle with one section per donor (#262 follow-up)."""
    from ..report import bundle_figure_pdf, resolve_report_name

    setup_logging(args.verbose)
    labels = [resolve_report_name(d) for d in args.run_dirs]
    out = bundle_figure_pdf(
        args.run_dirs,
        args.output,
        labels=labels,
        title="TCRsift — cross-donor figures",
        subtitle=f"{len(args.run_dirs)} donors",
    )
    print(f"Wrote cross-donor figure bundle: {out}")


def cmd_report_selected(args):
    """Assemble a selected clone set into a synthesis-ready PDF + CSV (#188)."""
    import anndata as ad
    import pandas as pd

    from ..config import TCRsiftConfig
    from ..report import build_selected_report

    setup_logging(args.verbose)
    selected = pd.read_csv(args.selection)
    clonotypes = pd.read_csv(args.clonotypes)
    obs = ad.read_h5ad(args.h5ad).obs if args.h5ad else None
    print(f"Selecting {selected['CDR3ab'].nunique()} clones from {args.clonotypes}")

    annotations = None
    if getattr(args, "annotations", None):
        annotations = pd.read_csv(args.annotations)
    elif "is_viral" in clonotypes.columns or "db_database" in clonotypes.columns:
        annotations = clonotypes

    assemble_kwargs = {}
    if args.config:
        cfg = TCRsiftConfig.from_yaml(args.config)
        asm = cfg.assemble
        assemble_kwargs = dict(
            alpha_leader=asm.alpha_leader,
            beta_leader=asm.beta_leader,
            include_constant=asm.include_constant,
            constant_source=asm.constant_source,
            linker=asm.linker if asm.single_chain else None,
            contigs_dir=asm.contigs_dir,
            sample_name_from=asm.sample_name_from,
            cellranger_dir=asm.cellranger_dir,
            trac_allele=getattr(asm, "trac_allele", "auto"),
            trbc1_allele=getattr(asm, "trbc1_allele", "auto"),
            trbc2_allele=getattr(asm, "trbc2_allele", "auto"),
            stop_codons=tuple(getattr(asm, "stop_codons", ("TAA", "TGA"))),
            leader_fallback=getattr(asm, "leader_fallback", None) or "germline",
            force_alpha_leader=getattr(asm, "force_alpha_leader", None),
            force_beta_leader=getattr(asm, "force_beta_leader", None),
            secondary_alpha_leader=getattr(asm, "secondary_alpha_leader", None),
            secondary_beta_leader=getattr(asm, "secondary_beta_leader", None),
            codon_optimization=getattr(asm, "codon_optimization", "focal"),
            max_codon_repeats=getattr(asm, "max_codon_repeats", 1),
            avoid_enzymes=tuple(getattr(asm, "avoid_enzymes", []) or []),
            species=getattr(cfg, "species", "auto") or "auto",
            constant_species=getattr(asm, "constant_species", None),
            codon_species=getattr(asm, "codon_species", None),
            germline_vdj_fallback=bool(getattr(asm, "germline_vdj_fallback", True)),
            auto_download_germline=bool(
                getattr(asm, "auto_download_germline", False)
                or getattr(cfg.annotate, "auto_download", False)
            ),
        )
    if args.cellranger_dir:
        assemble_kwargs["cellranger_dir"] = args.cellranger_dir
        if not args.contigs_dir:
            assemble_kwargs["contigs_dir"] = None
    if args.contigs_dir:
        assemble_kwargs["contigs_dir"] = args.contigs_dir
        if not args.cellranger_dir:
            assemble_kwargs["cellranger_dir"] = None
    if args.sample_name_from:
        assemble_kwargs["sample_name_from"] = args.sample_name_from
    if getattr(args, "sample_sheet", None):
        from ..validation import validate_file_exists

        validate_file_exists(args.sample_sheet, "sample sheet (--sample-sheet)")
        assemble_kwargs["sample_sheet"] = args.sample_sheet
    if (
        assemble_kwargs.get("sample_name_from") == "sheet"
        and not assemble_kwargs.get("sample_sheet")
    ):
        from ..validation import TCRsiftValidationError

        raise TCRsiftValidationError(
            "--sample-sheet is required when --sample-name-from sheet is used",
            hint="Provide --sample-sheet <path> or use --sample-name-from parent/grandparent.",
        )
    if getattr(args, "species", None) is not None:
        assemble_kwargs["species"] = args.species
    assemble_kwargs.setdefault("species", "auto")
    if getattr(args, "constant_species", None) is not None:
        assemble_kwargs["constant_species"] = args.constant_species
    if getattr(args, "codon_species", None) is not None:
        assemble_kwargs["codon_species"] = args.codon_species
    if getattr(args, "codon_optimization", None) is not None:
        assemble_kwargs["codon_optimization"] = args.codon_optimization
    if getattr(args, "max_codon_repeats", None) is not None:
        assemble_kwargs["max_codon_repeats"] = args.max_codon_repeats
    if getattr(args, "avoid_enzymes", None) is not None:
        assemble_kwargs["avoid_enzymes"] = tuple(args.avoid_enzymes)
    if getattr(args, "germline_vdj_fallback", None) is not None:
        assemble_kwargs["germline_vdj_fallback"] = bool(args.germline_vdj_fallback)
    if getattr(args, "auto_download_germline", False):
        assemble_kwargs["auto_download_germline"] = True
    if getattr(args, "cache_dir", None):
        assemble_kwargs["data_dir"] = args.cache_dir

    prov_cols = [c for c in selected.columns if c != "CDR3ab"]
    build_selected_report(
        selected,
        clonotypes,
        args.output_dir,
        obs=obs,
        assemble_kwargs=assemble_kwargs,
        provenance_cols=prov_cols,
        cover=not getattr(args, "no_cover", False),
        allow_canonical_fallback=getattr(args, "allow_canonical_fallback", False),
        annotations=annotations,
    )
    print(f"Wrote selected-clones report to {args.output_dir}")


def cmd_report_selected_combine(args):
    """Combine per-donor selected-clone PDFs into one cohort PDF (#202)."""
    from ..report import combine_selected_pdfs

    setup_logging(args.verbose)
    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    out = combine_selected_pdfs(
        args.pdfs,
        args.output,
        labels=labels,
        title=args.title,
        include_legend=not args.no_legend,
    )
    print(f"Wrote combined selected report: {out}")


def cmd_report_sequence(args):
    """Render a sequence PDF with per-method evidence + selection overlay (#123)."""
    import pandas as pd

    from ..plots import create_tcr_sequence_pdf

    setup_logging(args.verbose)
    df = pd.read_csv(args.input)
    if args.sort_by and args.sort_by in df.columns:
        df = df.sort_values(args.sort_by).reset_index(drop=True)

    annotations = None
    if args.per_method_source:
        from ..clonotype import build_clone_method_long
        from ..selection import attach_method_tiers, build_pdf_annotations

        cml = attach_method_tiers(build_clone_method_long(pd.read_csv(args.per_method_source)))
        sel = df if "selection_rule" in df.columns else None
        annotations = build_pdf_annotations(cml, selection_df=sel)

    create_tcr_sequence_pdf(df, args.output, strict=False, annotations=annotations)
    print(f"Wrote sequence PDF ({len(df)} clones): {args.output}")


def cmd_report_cohort(args):
    """Cross-donor overlap report (Jaccard / cell-weighted / beta-only) (#123)."""
    from ..cohort import run_cohort_analysis
    from ..validation import TCRsiftValidationError

    setup_logging(args.verbose)
    donor_dirs: dict[str, str] = {}
    for spec in args.donor or []:
        if "=" not in spec:
            raise TCRsiftValidationError(f"--donor must be NAME=DIR, got {spec!r}")
        name, ddir = spec.split("=", 1)
        donor_dirs[name.strip()] = ddir.strip()
    if len(donor_dirs) < 2:
        raise TCRsiftValidationError(
            "report cohort needs >=2 donors (--donor NAME=DIR, repeatable)"
        )
    matrices = run_cohort_analysis(
        donor_dirs,
        args.output_dir,
        selected_only=args.selected_only,
        include_beta_only=not args.no_beta_only,
        emit_tables=not args.no_tables,
        emit_plots=not args.no_plots,
    )
    print(f"Cohort overlap across {len(donor_dirs)} donors → {args.output_dir}")
    for name, mat in matrices.items():
        print(f"  cohort_{name}: {mat.shape[0]}x{mat.shape[1]}")


def cmd_report_polished(args):
    """Re-render a run's figures with a publication style profile (#123)."""
    import pandas as pd

    from ..plots import (
        get_plot_format,
        plot_clonotypes,
        plot_filter,
        set_plot_format,
        set_polished_style,
    )

    setup_logging(args.verbose)
    prev_fmt = get_plot_format()
    set_polished_style(args.style)
    set_plot_format("pdf")
    try:
        run_dir = Path(args.input)
        data = run_dir / "data" if (run_dir / "data").is_dir() else run_dir
        out = Path(args.output_dir)
        out.mkdir(parents=True, exist_ok=True)

        n = 0
        clo = data / "clonotypes.csv"
        if clo.exists():
            plot_clonotypes(pd.read_csv(clo), out)
            n += 1
        tier_csvs = sorted(data.glob("filtered_tier*.csv"))
        if tier_csvs:
            filt = pd.concat([pd.read_csv(t) for t in tier_csvs], ignore_index=True)
            plot_filter(filt, out)
            n += 1
        if n == 0:
            print(f"  No re-renderable data (clonotypes.csv / filtered_tier*.csv) under {data}")
        print(f"Wrote polished figures ({args.style}) → {out}")
    finally:
        set_plot_format(prev_fmt)


def register(subparsers) -> None:
    """Register the ``report`` command group."""
    p_report = subparsers.add_parser("report", help="Unified reporting (figure bundles, selected clones, cohorts)")
    p_report.set_defaults(func=lambda _a: (p_report.print_help(), 1)[1])
    p_report_sub = p_report.add_subparsers(dest="report_command")

    pb = p_report_sub.add_parser(
        "bundle",
        help="Concatenate per-run figures into one categorized cohort PDF",
    )
    pb.add_argument(
        "run_dirs",
        nargs="+",
        help="Run output dirs (each with a plots/ subdir, or a plots dir directly)",
    )
    pb.add_argument("--output", "-o", required=True, help="Output bundle PDF")
    pb.add_argument(
        "--labels",
        help="Comma-separated section labels per dir (default: directory names)",
    )
    pb.add_argument("--title", default="TCRsift figure set", help="Cover-page title")
    pb.add_argument("--subtitle", default="", help="Cover-page subtitle")
    pb.add_argument("--verbose", action="store_true", help="Verbose output")
    pb.set_defaults(func=cmd_report_bundle)

    pcf = p_report_sub.add_parser(
        "cohort-figures",
        help="Cross-donor figure bundle → all-figures-cross-donor.pdf (one section "
        "per donor, donor name resolved from each run dir)",
    )
    pcf.add_argument(
        "run_dirs",
        nargs="+",
        help="Per-donor run output dirs (each with a plots/ subdir)",
    )
    pcf.add_argument(
        "--output",
        "-o",
        default="all-figures-cross-donor.pdf",
        help="Output PDF (default: all-figures-cross-donor.pdf)",
    )
    pcf.add_argument("--verbose", action="store_true", help="Verbose output")
    pcf.set_defaults(func=cmd_report_cohort_figures)

    ps = p_report_sub.add_parser(
        "selected",
        help="Assemble a selected clone set into a synthesis-ready PDF + CSV",
    )
    ps.add_argument(
        "--selection",
        required=True,
        help="Selected-clones CSV (CDR3ab + provenance columns, e.g. from `tcrsift select`)",
    )
    ps.add_argument(
        "--clonotypes",
        required=True,
        help="clonotypes.csv (all clones, for VDJ/genes — incl. rescued/low-tier picks)",
    )
    ps.add_argument(
        "--annotations",
        help="annotated.csv (public-DB IEDB/CEDAR/VDJdb matches) to surface "
        "is_viral / db_category / db_epitope / db_species on selected_clones.csv. "
        "Optional — auto-used if --clonotypes already carries the db_* columns.",
    )
    ps.add_argument(
        "--h5ad",
        help="phenotyped.h5ad — needed to emit both alpha-variants of merged dual-alpha clones",
    )
    ps.add_argument("--config", "-c", help="YAML config (assemble options)")
    ps.add_argument(
        "--cellranger-dir",
        help="Per-donor CellRanger per_sample_outs dir for allele-faithful "
        "constants (overrides config; mirrors `run`, #196)",
    )
    ps.add_argument(
        "--contigs-dir",
        help="Per-donor contigs dir (overrides config, #196)",
    )
    ps.add_argument(
        "--sample-name-from",
        choices=["parent", "grandparent", "sheet"],
        help="Contig sample-name policy (overrides config, #196); 'sheet' requires --sample-sheet",
    )
    ps.add_argument(
        "--sample-sheet",
        metavar="FILE",
        default=None,
        help="Sample sheet CSV/YAML mapping vdj_dir paths to sample names "
        "(required when --sample-name-from sheet)",
    )
    ps.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species override ('auto', 'human', or 'mouse'; overrides config)",
    )
    ps.add_argument(
        "--constant-species",
        choices=["human", "mouse"],
        default=None,
        help="Constant-region species override ('human' or 'mouse'; overrides config)",
    )
    ps.add_argument(
        "--codon-species",
        choices=["human", "mouse"],
        default=None,
        help="Codon-optimization table override ('human' or 'mouse'; overrides config)",
    )
    ps.add_argument(
        "--codon-optimization",
        choices=["focal", "full"],
        default=None,
        help="How full_*_nt_optimized is generated ('focal' or 'full'; overrides config)",
    )
    ps.add_argument(
        "--max-codon-repeats",
        type=int,
        default=None,
        help="Longest allowed run of identical consecutive codons in focal/full recoding",
    )
    ps.add_argument(
        "--avoid-enzymes",
        action="extend",
        nargs="+",
        default=None,
        metavar="ENZYME",
        help="Restriction enzyme name(s) or ACGT site(s) to recode out of constructs",
    )
    ps.add_argument(
        "--germline-vdj-fallback",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Auto-detect and rescue missing or 5'-truncated VDJ chains from cached IMGT V+J germlines "
        "(default: enabled; pass --no-germline-vdj-fallback for strict contig-only assembly)",
    )
    ps.add_argument(
        "--auto-download-germline",
        action="store_true",
        help="Automatically fetch missing IMGT germline V/J references when germline VDJ fallback is enabled",
    )
    ps.add_argument(
        "--cache-dir",
        default=None,
        help="Override local reference-data cache directory (default: ~/.cache/tcrsift)",
    )
    ps.add_argument("--output-dir", "-o", required=True, help="Output directory")
    ps.add_argument(
        "--no-cover",
        action="store_true",
        help="Skip the construct cover/legend page on the sequence PDF (#202)",
    )
    ps.add_argument("--verbose", action="store_true", help="Verbose output")
    ps.add_argument(
        "--allow-canonical-fallback",
        action="store_true",
        help="Permit constructs assembled from canonical/reference constants when "
        "no contig verifies the donor's junction/allele. Without it, report "
        "selected fails closed if any construct is unverified (#241/#243).",
    )
    ps.set_defaults(func=cmd_report_selected)

    pscomb = p_report_sub.add_parser(
        "selected-combine",
        help="Combine per-donor selected_clones_sequences.pdf into one cohort PDF",
    )
    pscomb.add_argument(
        "pdfs",
        nargs="+",
        help="Per-donor selected_clones_sequences.pdf files",
    )
    pscomb.add_argument("--output", "-o", required=True, help="Output combined PDF")
    pscomb.add_argument("--labels", help="Comma-separated donor labels (default: dir names)")
    pscomb.add_argument(
        "--title",
        default="Selected clones — cohort",
        help="Cover-page title",
    )
    pscomb.add_argument(
        "--no-legend",
        action="store_true",
        help="Skip the construct legend page",
    )
    pscomb.add_argument("--verbose", action="store_true", help="Verbose output")
    pscomb.set_defaults(func=cmd_report_selected_combine)

    pseq = p_report_sub.add_parser(
        "sequence",
        help="Render a sequence PDF with per-method evidence overlay",
    )
    pseq.add_argument("--input", "-i", required=True, help="Clones CSV (assembled full_* columns)")
    pseq.add_argument("--output", "-o", required=True, help="Output PDF")
    pseq.add_argument(
        "--per-method-source",
        help="clone_sample_long.csv to overlay per-method evidence (cells/freq/tier)",
    )
    pseq.add_argument("--sort-by", help="Column to sort pages by (e.g. global_rank)")
    pseq.add_argument("--verbose", action="store_true", help="Verbose output")
    pseq.set_defaults(func=cmd_report_sequence)

    pcoh = p_report_sub.add_parser(
        "cohort",
        help="Cross-donor overlap report (Jaccard / cell-weighted / beta-only)",
    )
    pcoh.add_argument(
        "--donor",
        action="append",
        metavar="NAME=DIR",
        help="Donor run dir (repeatable; >=2 required)",
    )
    pcoh.add_argument("--output-dir", "-o", required=True, help="Output directory")
    pcoh.add_argument("--selected-only", action="store_true", help="Restrict to selected clones")
    pcoh.add_argument("--no-beta-only", action="store_true", help="Skip the beta-only overlap")
    pcoh.add_argument("--no-tables", action="store_true", help="Skip CSV matrix exports")
    pcoh.add_argument("--no-plots", action="store_true", help="Skip heatmap PDFs")
    pcoh.add_argument("--verbose", action="store_true", help="Verbose output")
    pcoh.set_defaults(func=cmd_report_cohort)

    ppol = p_report_sub.add_parser(
        "polished",
        help="Re-render a run's figures with a publication style profile",
    )
    ppol.add_argument("--input", "-i", required=True, help="Run output dir (with data/ + plots)")
    ppol.add_argument("--output-dir", "-o", required=True, help="Output directory")
    ppol.add_argument(
        "--style",
        choices=["clean-white", "paper"],
        default="clean-white",
        help="Style profile (default: clean-white)",
    )
    ppol.add_argument("--verbose", action="store_true", help="Verbose output")
    ppol.set_defaults(func=cmd_report_polished)
