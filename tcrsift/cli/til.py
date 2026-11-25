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

"""TIL matching, clonotyping, and ex vivo prioritization CLI subcommands.

Commands: ``match-til``, ``til-clonotype``, ``prioritize``, ``til-prioritize``,
``til-select``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..validation import TCRsiftValidationError
from .common import setup_logging


def _parse_per_lineage_min(spec: str | None) -> dict | None:
    """Parse a ``CD8=2,CD4=5`` string into ``{"CD8": 2, "CD4": 5}`` (#328)."""
    if not spec:
        return None
    out: dict = {}
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise argparse.ArgumentTypeError(
                f"--per-lineage-min entry {part!r} must be LINEAGE=INT (e.g. CD8=2)"
            )
        key, val = part.split("=", 1)
        try:
            out[key.strip()] = int(val)
        except ValueError as e:
            raise argparse.ArgumentTypeError(
                f"--per-lineage-min value in {part!r} must be an integer"
            ) from e
    return out or None


def _load_til_source(args):
    """Load TIL data from whichever CLI source flag was provided."""
    from ..sample_sheet import load_sample_sheet
    from ..til import load_til_data, load_til_samples, load_til_specs

    if args.til_sample:
        til_data = load_til_specs(args.til_sample)
        total_cells = sum(len(df) for df in til_data.values())
        print(
            f"Loaded {total_cells} TIL cells from {len(til_data)} sample(s): {list(til_data.keys())}"
        )
        return til_data

    if args.sample_sheet:
        sample_sheet = load_sample_sheet(args.sample_sheet)
        til_samples = sample_sheet.get_til_samples()
        if not til_samples:
            raise TCRsiftValidationError(
                "No TIL samples found in sample sheet",
                hint="Make sure your sample sheet has samples with source='til'",
            )
        til_data = load_til_samples(til_samples)
        total_cells = sum(len(df) for df in til_data.values())
        print(
            f"Loaded {total_cells} TIL cells from {len(til_data)} sample(s): {list(til_data.keys())}"
        )
        return til_data

    if args.til_h5ad:
        til_data = load_til_data("h5ad", args.til_h5ad)
        print(f"Loaded {len(til_data)} TIL cells from h5ad: {args.til_h5ad}")
        return til_data

    if args.til_csv:
        til_data = load_til_data("csv", args.til_csv)
        print(f"Loaded {len(til_data)} TIL cells from CSV: {args.til_csv}")
        return til_data

    if args.til_vdj_dir:
        til_data = load_til_data("vdj_dir", args.til_vdj_dir)
        print(f"Loaded {len(til_data)} TIL cells from VDJ directory: {args.til_vdj_dir}")
        return til_data

    raise TCRsiftValidationError(
        "No TIL data source specified",
        hint="Provide one of: --sample-sheet, --til-h5ad, --til-csv, --til-vdj-dir, or --til-sample",
    )


def cmd_match_til(args):
    """Match culture clonotypes against TIL data."""
    import pandas as pd

    from ..plots import plot_til
    from ..til import get_til_summary, match_til
    from ..validation import validate_match_til_args

    setup_logging(args.verbose)
    validate_match_til_args(args)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} culture clonotypes from {args.input}")

    til_data = _load_til_source(args)

    matched = match_til(
        clonotypes,
        til_data,
        match_by=args.match_by,
        min_til_cells=args.min_til_cells,
    )

    summary = get_til_summary(matched)
    print("\nTIL Matching Summary:")
    for key, value in summary.items():
        if isinstance(value, dict):
            print(f"  {key}:")
            for k, v in value.items():
                print(f"    {k}: {v}")
        else:
            print(f"  {key}: {value}")

    matched.to_csv(args.output, index=False)
    print(f"\nSaved TIL-matched clonotypes to {args.output}")

    if args.plot_til:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_til(matched, output_dir)
        print(f"TIL plots saved to {output_dir}")


def cmd_til_clonotype(args):
    """Aggregate TIL-only data into clonotype counts across samples."""
    from ..til import summarize_til_clonotypes
    from ..validation import validate_match_til_args

    setup_logging(args.verbose)
    validate_match_til_args(args)

    til_data = _load_til_source(args)
    clonotypes = summarize_til_clonotypes(
        til_data,
        match_by=args.match_by,
        min_cells=args.min_cells,
    )

    clonotypes.to_csv(args.output, index=False)
    print(f"Saved {len(clonotypes)} TIL clonotypes to {args.output}")


def cmd_til_prioritize(args):
    """Prioritize clones from a standard VDJ+GEX sample sheet."""
    from ..til_prioritize import run_til_prioritize

    setup_logging(args.verbose)
    _, candidates = run_til_prioritize(args)
    print(f"Wrote {len(candidates)} candidates and the full audit table to {args.output_dir}")


def cmd_til_select(args):
    """Prioritize promising TIL clonotypes from multi-timepoint VDJ+GEX data."""
    from ..til_select import run_til_select

    setup_logging(args.verbose)
    master_df = run_til_select(args)
    n_final = int(master_df["is_candidate_tumor_reactive"].sum())
    print(f"Saved harmonized table: {args.out_table}")
    print(f"Final candidate clonotypes: {n_final}")
    if getattr(args, "dominant_specificity", False) and "is_dominant_specific" in master_df.columns:
        n_dom = int(master_df["is_dominant_specific"].sum())
        print(f"Dominant-specificity (condition-concentrated) clones: {n_dom}")
    print(f"Figures/subsets directory: {args.fig_dir}")


def register(subparsers) -> None:
    """Register ``match-til``, ``til-clonotype``, ``prioritize``, ``til-prioritize``, and ``til-select``."""
    from ..signatures import AIM_GENES_HGNC, EXHAUSTION_GENES_HGNC, MARKER_PANEL_HGNC
    from ..til_prioritize import add_cli_args as add_til_prioritize_args

    # match-til
    p_til = subparsers.add_parser(
        "match-til",
        help="Match against TIL data",
        description="""
Match culture-expanded clonotypes against TIL (tumor-infiltrating lymphocyte) data.

REQUIRED INPUTS:
  --input/-i     Clonotypes CSV from culture expansion
  --output/-o    Output CSV path

TIL DATA SOURCE (provide ONE of the following):
  --sample-sheet   YAML/CSV sample sheet with TIL samples (source='til')
                   Supports multiple TIL samples with different input types
  --til-h5ad       Single TIL AnnData h5ad file
  --til-csv        Single TIL CSV file (must have CDR3_alpha/CDR3_beta columns)
  --til-vdj-dir    Single CellRanger VDJ output directory
  --til-sample     Repeatable direct TIL sample spec:
                   NAME=TYPE:PATH (TYPE is csv|h5ad|vdj)
                   Example: --til-sample T1=csv:/path/t1.csv
                            --til-sample T2=vdj:/path/vdj_outs

SAMPLE SHEET FORMAT (for multiple TIL samples):
  samples:
    - sample: "TIL_Sample1"
      source: til
      vdj_dir: "/path/to/vdj"         # CellRanger VDJ directory
    - sample: "TIL_Sample2"
      source: til
      til_csv: "/path/to/til.csv"     # CSV file
    - sample: "TIL_Sample3"
      source: til
      til_h5ad: "/path/to/til.h5ad"   # Pre-processed h5ad

OUTPUT:
  For multiple TIL samples, per-sample columns are added:
    til_cell_count.{sample}   - Cell count in specific TIL sample
    til_frequency.{sample}    - Frequency in specific TIL sample
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    til_required = p_til.add_argument_group("required arguments")
    til_required.add_argument("--input", "-i", required=True, help="Input culture clonotypes CSV")
    til_required.add_argument(
        "--output", "-o", required=True, help="Output CSV with TIL match annotations"
    )

    til_source = p_til.add_argument_group(
        "TIL data source (provide one)",
        "Multiple input formats supported. For multi-sample TIL matching, use --sample-sheet or repeat --til-sample.",
    )
    til_source.add_argument(
        "--sample-sheet", "-s", metavar="PATH", help="Sample sheet with TIL samples (YAML or CSV)"
    )
    til_source.add_argument(
        "--til-h5ad",
        metavar="PATH",
        help="Single TIL h5ad file (AnnData with CDR3 columns in .obs)",
    )
    til_source.add_argument(
        "--til-csv", metavar="PATH", help="Single TIL CSV file (must have CDR3_alpha/CDR3_beta)"
    )
    til_source.add_argument(
        "--til-vdj-dir", metavar="PATH", help="Single CellRanger VDJ output directory"
    )
    til_source.add_argument(
        "--til-sample",
        action="append",
        metavar="NAME=TYPE:PATH",
        help="Repeatable direct TIL sample spec. TYPE is csv, h5ad, or vdj (alias for vdj_dir).",
    )

    til_opts = p_til.add_argument_group("matching options")
    til_opts.add_argument(
        "--match-by",
        choices=["CDR3ab", "CDR3b_only"],
        default="CDR3ab",
        help="Match by both chains or beta only (default: CDR3ab)",
    )
    til_opts.add_argument(
        "--min-til-cells", type=int, default=1, help="Min TIL cells to count as match (default: 1)"
    )

    til_out = p_til.add_argument_group("output options")
    til_out.add_argument("--plot-til", action="store_true", help="Generate TIL matching plots")
    til_out.add_argument(
        "--output-dir", metavar="DIR", help="Output directory for plots (default: next to output)"
    )
    til_out.add_argument("--verbose", action="store_true", help="Verbose output")
    p_til.set_defaults(func=cmd_match_til)

    # til-clonotype
    p_til_clono = subparsers.add_parser(
        "til-clonotype",
        help="Aggregate TIL-only data into clonotype counts",
        description="""
Aggregate TIL-only data into clonotype counts/frequencies across one or more TIL samples.

TIL DATA SOURCE (provide ONE of the following):
  --sample-sheet   YAML/CSV sample sheet with TIL samples (source='til')
  --til-h5ad       Single TIL AnnData h5ad file
  --til-csv        Single TIL CSV file
  --til-vdj-dir    Single CellRanger VDJ output directory
  --til-sample     Repeatable direct TIL sample spec:
                   NAME=TYPE:PATH (TYPE is csv|h5ad|vdj)
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_til_clono.add_argument("--output", "-o", required=True, help="Output CSV path")
    p_til_clono.add_argument("--sample-sheet", "-s", metavar="PATH", help="TIL sample sheet")
    p_til_clono.add_argument("--til-h5ad", metavar="PATH", help="Single TIL h5ad file")
    p_til_clono.add_argument("--til-csv", metavar="PATH", help="Single TIL CSV file")
    p_til_clono.add_argument(
        "--til-vdj-dir", metavar="PATH", help="Single CellRanger VDJ output directory"
    )
    p_til_clono.add_argument(
        "--til-sample",
        action="append",
        metavar="NAME=TYPE:PATH",
        help="Repeatable direct TIL sample spec. TYPE is csv, h5ad, or vdj (alias for vdj_dir).",
    )
    p_til_clono.add_argument(
        "--match-by",
        choices=["CDR3ab", "CDR3b_only"],
        default="CDR3ab",
        help="Aggregate by both chains or beta-only (default: CDR3ab)",
    )
    p_til_clono.add_argument(
        "--min-cells",
        type=int,
        default=1,
        help="Minimum total TIL cells per clonotype (default: 1)",
    )
    p_til_clono.add_argument("--verbose", action="store_true", help="Verbose output")
    p_til_clono.set_defaults(func=cmd_til_clonotype)

    # prioritize / til-prioritize
    for command, context in (("prioritize", "generic"), ("til-prioritize", "solid-tumor")):
        cmd_help = (
            "Prioritize ex vivo tissue or blood clonotypes across expression signatures"
            if command == "prioritize"
            else "Alias for 'prioritize --context solid-tumor'"
        )
        p_prioritize = subparsers.add_parser(
            command,
            help=cmd_help,
            description=(
                "Prioritize paired CellRanger VDJ + GEX clones. Score within sample and "
                "CD4/CD8 lineage, then select round-robin from each signature/sample list "
                "with one total clone budget and per-list score cutoffs. Writes candidate_clones.csv, "
                "all_scored_clones.csv, clone_sample_scores.csv, prioritization.json, and "
                "(with --assemble) candidate_constructs.csv + candidate_constructs.fasta."
            ),
        )
        add_til_prioritize_args(p_prioritize, context=context)
        p_prioritize.set_defaults(func=cmd_til_prioritize)

    # til-select
    p_til_select = subparsers.add_parser(
        "til-select",
        help="Select promising TIL clonotypes from 2+ ordered VDJ+GEX samples",
        description="""
Prioritize promising TIL clonotypes from two or more ordered tumor samples or
timepoints. Independent samples can be ordered for harmonization, but frequency
increase flags are meaningful only for a genuine longitudinal series.

Input model is compatible with legacy v2 scripts:
  - consensus_annotations.<TP>.csv
  - clonotypes.<TP>.csv
  - filtered_contig_annotations.<TP>.csv
  - sample_filtered_feature_bc_matrix.<TP>.h5

Timepoint mappings can be provided explicitly with:
  --samples T1=consensus_annotations.T1.csv,clonotypes.T1.csv ...
or through:
  --config config.yaml
or auto-discovered from --data-dir.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_til_select.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional YAML mapping timepoints to consensus/clonotypes paths.",
    )
    p_til_select.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Directory containing per-timepoint input files (default: ./data).",
    )
    p_til_select.add_argument(
        "--samples",
        "--inputs",
        dest="samples",
        nargs="+",
        default=None,
        help=(
            "Explicit inputs as LABEL=CONSENSUS_PATH,CLONOTYPES_PATH. "
            "Example: T1=consensus_annotations.T1.csv,clonotypes.T1.csv"
        ),
    )
    p_til_select.add_argument(
        "--count-column",
        type=str,
        default=None,
        help="Optional clonotype count column; auto-detected when omitted.",
    )
    p_til_select.add_argument("--verbose", action="store_true", help="Verbose logging")
    p_til_select.add_argument("--top-k", type=int, default=20, help="Top-k clones to output/plot")
    p_til_select.add_argument(
        "--min-cells-per-clone",
        type=int,
        default=2,
        help="Minimum total cells across timepoints for base filtering (default: 2).",
    )
    p_til_select.add_argument(
        "--min-cd8-cp10k",
        type=float,
        default=0.0,
        help="Minimum CD8 CP10K threshold for base filtering (strict >, default: 0).",
    )
    p_til_select.add_argument(
        "--max-cd4-to-cd8-ratio",
        type=float,
        default=1.0,
        help="Maximum CD4/CD8 ratio threshold for base filtering (strict <, default: 1.0).",
    )
    p_til_select.add_argument(
        "--increase-ratio-nonzero-min",
        type=float,
        default=1.5,
        help="Minimum last-vs-prior nonzero frequency ratio for increasing branch (default: 1.5).",
    )
    p_til_select.add_argument(
        "--increase-ratio-all-timepoints-min",
        type=float,
        default=1.5,
        help="Minimum last-vs-all-timepoints ratio for all-positive branch (default: 1.5).",
    )
    p_til_select.add_argument(
        "--immunogenic-percentile",
        type=float,
        default=0.90,
        help="Percentile threshold for per-gene/panel immunogenic selection (default: 0.90).",
    )
    p_til_select.add_argument(
        "--immunogenic-percentile-slack-frac",
        type=float,
        default=0.01,
        help="Relative slack below percentile threshold for immunogenic selection (default: 0.01).",
    )
    p_til_select.add_argument(
        "--immunogenic-min-cp10k",
        type=float,
        default=0.0,
        help="Minimum CP10K required for immunogenic eligibility (default: 0).",
    )
    p_til_select.add_argument(
        "--immunogenic-require-above-median",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Require per-gene/panel score above base-selected median (default: false).",
    )
    p_til_select.add_argument(
        "--cytotoxic-last-min-z",
        type=float,
        default=0.25,
        help="Minimum terminal cytotoxic z-score (default: 0.25).",
    )
    p_til_select.add_argument(
        "--cytotoxic-last-min-cp10k",
        type=float,
        default=0.05,
        help="Minimum terminal cytotoxic CP10K score (default: 0.05).",
    )
    p_til_select.add_argument(
        "--became-cytotoxic-min-delta-z",
        type=float,
        default=0.5,
        help="Minimum cytotoxic z-score increase (last-first) for became-cytotoxic flag.",
    )
    p_til_select.add_argument(
        "--trend-increase-ratio-min",
        type=float,
        default=1.5,
        help=argparse.SUPPRESS,
    )
    p_til_select.add_argument(
        "--trend-decrease-ratio-max",
        type=float,
        default=0.5,
        help=argparse.SUPPRESS,
    )
    p_til_select.add_argument(
        "--marker-genes",
        type=str,
        default=",".join(MARKER_PANEL_HGNC),
        help="Comma-separated marker genes for per-clone GEX scoring.",
    )
    p_til_select.add_argument(
        "--immunogenic-genes",
        type=str,
        default="GZMB,PRF1,IFNG,MKI67,TNFRSF9",
        help="Comma-separated genes for immunogenic branch ranking.",
    )
    p_til_select.add_argument(
        "--cytotoxic-genes",
        type=str,
        default="GZMB,PRF1,IFNG,MKI67,TNFRSF9",
        help="Comma-separated genes for cytotoxic aggregate scores.",
    )
    p_til_select.add_argument(
        "--cytolytic-genes",
        type=str,
        default="PRF1,GZMB",
        help="Comma-separated genes for cytolytic score branch.",
    )
    p_til_select.add_argument(
        "--antigen-response-genes",
        type=str,
        default="TNFRSF9,MKI67",
        help="Comma-separated genes for antigen-response score branch.",
    )
    p_til_select.add_argument(
        "--activation-genes",
        type=str,
        default=",".join(AIM_GENES_HGNC),
        help="Comma-separated genes for the activation (AIM) score branch (#303).",
    )
    p_til_select.add_argument(
        "--exhaustion-genes",
        type=str,
        default=None,
        help="Comma-separated genes for the exhaustion score branch (#303; "
        f"default: {','.join(EXHAUSTION_GENES_HGNC)}).",
    )
    p_til_select.add_argument(
        "--enrichment-genes",
        type=str,
        default="",
        help="DEPRECATED (#303): alias for --exhaustion-genes; seeds the "
        "exhaustion branch when --exhaustion-genes is not given.",
    )
    p_til_select.add_argument(
        "--pyensembl-release",
        type=int,
        default=110,
        help=argparse.SUPPRESS,
    )
    p_til_select.add_argument(
        "--rank-by",
        type=str,
        default="mean_frequency",
        choices=[
            "mean_frequency",
            "max_frequency",
            "total_cells",
            "marker_score_cp10k_mean",
            "marker_score_z_mean",
        ],
        help="Ranking metric for top-k outputs.",
    )
    p_til_select.add_argument(
        "--fig-dir",
        type=Path,
        default=Path("figures"),
        help="Output directory for figures/subsets (default: ./figures).",
    )
    p_til_select.add_argument(
        "--out-table",
        type=Path,
        default=Path("abTCR_harmonized.csv"),
        help="Output CSV for harmonized table (default: ./abTCR_harmonized.csv).",
    )
    p_til_select.add_argument(
        "--out-heatmap",
        type=Path,
        default=None,
        help="Output PNG for top-k heatmap (default: FIG_DIR/abTCR_topk_heatmap.png).",
    )
    p_til_select.add_argument(
        "--out-top",
        type=Path,
        default=None,
        help="Output CSV for top-k table (default: FIG_DIR/abTCR_topk.csv).",
    )
    p_til_select.add_argument(
        "--vdjdb", type=Path, default=None, help="Path to VDJdb TSV (overrides cache)"
    )
    p_til_select.add_argument(
        "--iedb", type=Path, default=None, help="Path to IEDB CSV/TSV (overrides cache)"
    )
    p_til_select.add_argument(
        "--cedar", type=Path, default=None, help="Path to CEDAR CSV/TSV (overrides cache)"
    )
    p_til_select.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Reference-database cache directory (default: ~/.cache/tcrsift).",
    )
    p_til_select.add_argument(
        "--auto-download",
        action="store_true",
        help="Automatically download missing public TCR databases into the cache.",
    )
    p_til_select.add_argument(
        "--match-by",
        type=str,
        default="CDR3b_only",
        choices=["CDR3ab", "CDR3b_only"],
        help="Matching strategy for public DB annotation (default: CDR3b_only).",
    )
    p_til_select.add_argument(
        "--out-annotated",
        type=Path,
        default=None,
        help="Output CSV for annotated table (default: FIG_DIR/abTCR_annotated.csv).",
    )
    p_til_select.add_argument(
        "--out-annotated-heatmap",
        type=Path,
        default=None,
        help="Compatibility output path for annotated heatmap (reserved).",
    )
    p_til_select.add_argument(
        "--out-selected-report",
        type=Path,
        default=None,
        help="Output PDF for selected clone report (default: FIG_DIR/selected_clones_report.pdf).",
    )

    dom_group = p_til_select.add_argument_group("dominant-specificity selection (#328)")
    dom_group.add_argument(
        "--dominant-specificity",
        action="store_true",
        help="Flag clones confined to one sample/antigen pool (dominant-fraction "
        "rule): adds is_dominant_specific + specificity + dominant_antigen to the "
        "harmonized table and writes abTCR_dominant_specificity.csv. The condition "
        "axis is the per-sample cell counts, so name samples by antigen pool.",
    )
    dom_group.add_argument(
        "--min-specificity",
        type=float,
        default=0.95,
        help="Min dominant-condition fraction (dominant/total) to call a clone "
        "specific (default: 0.95).",
    )
    dom_group.add_argument(
        "--dominant-min-cells",
        type=int,
        default=5,
        help="Absolute dominant-condition cell floor when no --per-lineage-min "
        "entry applies (default: 5).",
    )
    dom_group.add_argument(
        "--per-lineage-min",
        type=_parse_per_lineage_min,
        default=None,
        metavar="CD8=2,CD4=5",
        help="Per-lineage dominant-cell floors, e.g. CD8=2,CD4=5. Lineage comes "
        "from --dominant-lineage-col or is auto-derived from cd4_to_cd8_ratio.",
    )
    dom_group.add_argument(
        "--dominant-lineage-col",
        type=str,
        default=None,
        help="Harmonized-table column holding the per-clone lineage for "
        "--per-lineage-min (default: auto-derive from cd4_to_cd8_ratio).",
    )
    dom_group.add_argument(
        "--antigen-labels",
        type=Path,
        default=None,
        help="CSV mapping condition -> antigen label (e.g. P15 -> 'P15 (EXOC4)'), "
        "columns condition,label. Prettifies dominant_antigen without a Python call.",
    )
    p_til_select.set_defaults(func=cmd_til_select)
