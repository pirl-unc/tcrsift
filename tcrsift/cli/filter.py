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

"""Filtering, selection, cross-experiment unification, and cohort CLI subcommands.

Commands: ``filter``, ``select``, ``unify``, ``cohort``.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from ..validation import TCRsiftValidationError
from .common import setup_logging

logger = logging.getLogger(__name__)


def cmd_filter(args):
    """Filter clonotypes with tiered confidence levels."""
    import pandas as pd

    from ..filter import (
        DEFAULT_FDR_TIERS,
        bucket_by_donor_sharing,
        filter_clonotypes,
        get_filter_summary,
        resolve_fdr_tiers_for_method,
        resolve_filter_mode_kwargs,
        split_by_tier,
    )
    from ..plots import plot_filter

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    fdr_tiers = None
    if args.fdr_tiers:
        fdr_tiers = [float(x) for x in args.fdr_tiers.split(",")]
    fdr_tiers, warn_fdr_inert = resolve_fdr_tiers_for_method(
        args.method, fdr_tiers, DEFAULT_FDR_TIERS
    )
    if warn_fdr_inert:
        logger.warning(
            "fdr_tiers customized but method=%r; fdr_tiers only applies to "
            "method='logistic' — using fixed abundance tiers, list ignored.",
            args.method,
        )

    mode = getattr(args, "filter_mode", "fdr")
    user_kwargs = {
        "min_donors": args.min_donors,
        "min_methods_per_donor": args.min_methods_per_donor,
        "min_cells_per_method": args.min_cells_per_method,
        "min_frequency_per_method": args.min_frequency_per_method,
        "min_timepoints": args.min_timepoints,
        "min_timepoints_per_donor": args.min_timepoints_per_donor,
        "min_apcs": args.min_apcs,
        "min_apcs_per_donor": args.min_apcs_per_donor,
        "min_til_cells_per_donor": args.min_til_cells_per_donor,
    }
    mode_kwargs = resolve_filter_mode_kwargs(mode, user_kwargs)

    if mode == "cross-donor-public" and "n_donors" in clonotypes.columns:
        if (clonotypes["n_donors"] >= 2).sum() == 0:
            logger.warning(
                "filter_mode='cross-donor-public' requested but no clone has "
                "n_donors >= 2. Result will be empty. Did you mean to set "
                "donors_share_antigen / supply patient_id across donors?"
            )

    filtered = filter_clonotypes(
        clonotypes,
        method=args.method,
        tcell_type=args.tcell_type,
        min_cells=args.min_cells,
        min_frequency=args.min_frequency,
        require_complete=args.require_complete,
        exclude_viral=args.exclude_viral,
        fdr_tiers=fdr_tiers,
        stim_fdr_threshold=getattr(args, "stim_fdr_threshold", 0.1),
        stim_log2_fc_threshold=getattr(args, "stim_log2_fc_threshold", 1.0),
        stim_delta_threshold=getattr(args, "stim_delta_threshold", 0.5),
        **mode_kwargs,
    )

    summary = get_filter_summary(filtered)
    print("\nFilter Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    tier_dfs = split_by_tier(filtered)
    for tier, tier_df in tier_dfs.items():
        tier_path = output_dir / f"{tier}.csv"
        tier_df.to_csv(tier_path, index=False)
        print(f"Saved {len(tier_df)} {tier} clonotypes to {tier_path}")

    if getattr(args, "emit_per_sample_tier", False):
        if not getattr(args, "clone_sample_long", None):
            print(
                "  --emit-per-sample-tier needs --clone-sample-long PATH "
                "(the clone_sample_long.csv from `tcrsift run`); skipping."
            )
        else:
            from ..selection import pivot_per_sample_tiers

            long_df = pd.read_csv(args.clone_sample_long)
            wide = pivot_per_sample_tiers(long_df)
            filtered = filtered.merge(wide, on="CDR3ab", how="left")
            n_tier_cols = sum(c.startswith("tier_in_") for c in filtered.columns)
            print(f"  Added {n_tier_cols} per-sample tier column(s)")

    filtered.to_csv(output_dir / "all_filtered.csv", index=False)

    if mode == "shared-high-freq":
        filtered.to_csv(output_dir / "all_shared_high_freq.csv", index=False)
    elif mode == "cross-donor-public":
        filtered.to_csv(output_dir / "all_cross_donor_public.csv", index=False)

    sharing_buckets = bucket_by_donor_sharing(filtered)
    for bucket, bucket_df in sharing_buckets.items():
        bpath = output_dir / f"{bucket}.csv"
        bucket_df.to_csv(bpath, index=False)
        print(f"Saved {len(bucket_df)} {bucket} clonotypes to {bpath}")

    if args.plot_filter:
        plot_dir = Path(args.output_dir) if args.output_dir else output_dir / "plots"
        plot_filter(filtered, plot_dir)
        print(f"Filter plots saved to {plot_dir}")


def cmd_select(args):
    """Apply the config-driven multi-route selection language to a long table.

    Standalone front-end (#122) to the selection engine that ``run`` uses
    internally: reads a clone-sample-long CSV (the ``clone_sample_long.csv``
    that ``run`` emits), runs the configured ``selection.rules``, and writes
    the selected clones (selection columns + ``global_rank``). Optionally joins
    per-clone metadata from a clonotypes CSV and emits per-sample tiers.
    """
    import pandas as pd

    from ..config import TCRsiftConfig
    from ..selection import attach_per_sample_tiers, select_from_clone_sample_long

    setup_logging(args.verbose)

    long_df = pd.read_csv(args.input)
    print(f"Loaded {len(long_df)} (clone, sample) rows from {args.input}")

    clone_scores = None
    scores_path = getattr(args, "clone_scores", None) or getattr(args, "clonotypes", None)
    if scores_path:
        clone_scores = pd.read_csv(scores_path)

    if getattr(args, "prism", False):
        from ..selection import select_freq_prism_per_condition

        if clone_scores is None:
            raise TCRsiftValidationError(
                "--prism needs --clone-scores (per-clone ppost/GEX scores).",
                hint="Run `tcrsift annotate` / add_pgen_ppost first to produce "
                "ppost_alpha/ppost_beta/antigen_response_score/naive_score.",
            )
        if args.condition_col:
            if args.condition_col not in long_df.columns:
                raise TCRsiftValidationError(
                    f"--condition-col {args.condition_col!r} is not a column in "
                    "the long table.",
                    hint=f"Available columns: {list(long_df.columns)}",
                )
            cond_col = args.condition_col
        else:
            cond_col = next(
                (c for c in ("method", "sample") if c in long_df.columns),
                None,
            )
            if cond_col is None:
                raise TCRsiftValidationError(
                    "No condition column in the long table.",
                    hint="Pass --condition-col (e.g. method or sample).",
                )
        score_cols = [c for c in clone_scores.columns if c != "CDR3ab"]
        feat = long_df.merge(clone_scores[["CDR3ab", *score_cols]], on="CDR3ab", how="left")
        tie_break = None
        if args.config:
            tie_break = (TCRsiftConfig.from_yaml(args.config).selection or {}).get("tie_break")
        if getattr(args, "tie_break", None) is not None:
            tie_break = [k.strip() for k in args.tie_break.split(",") if k.strip()]
        rescue_target = getattr(args, "rescue_target", 0) or 0
        rescue_rank_col = getattr(args, "rescue_rank_col", None)
        if rescue_target > 0 and not rescue_rank_col:
            raise TCRsiftValidationError(
                "--rescue-target > 0 requires --rescue-rank-col.",
                hint="Pass a coverage column (from --clone-scores) to rank "
                "sub-threshold clones for rescue.",
            )
        if rescue_rank_col and rescue_rank_col not in feat.columns:
            raise TCRsiftValidationError(
                f"--rescue-rank-col {rescue_rank_col!r} is not a column in the "
                "feature table.",
                hint="It must be present in --clone-scores or the long table.",
            )
        sel = select_freq_prism_per_condition(
            feat,
            condition_col=cond_col,
            freq_col="frequency",
            gate=args.gate,
            top_freq=args.top_freq,
            top_prism=args.top_prism,
            tie_break=tie_break,
            rescue_target=rescue_target,
            rescue_rank_col=rescue_rank_col,
        )
        if rescue_target > 0 and "rescued" in sel.columns:
            n_rescued = int(sel["rescued"].sum())
            if n_rescued:
                print(
                    f"  Coverage rescue: {n_rescued} sub-threshold clone×condition "
                    f"pick(s) added (rescue_target={rescue_target}, "
                    f"by {rescue_rank_col!r}) — low-confidence, flagged rescued=True."
                )
        sel.to_csv(args.output, index=False)
        print(
            f"Wrote {sel['CDR3ab'].nunique()} distinct clones "
            f"({len(sel)} clone×condition rows) by per-condition freq∪PRISM → {args.output}"
        )
        try:
            from ..plots import plot_freq_prism_scatter, plot_prism_components
            from ..selection import prism_candidates

            cand = prism_candidates(
                feat,
                condition_col=cond_col,
                freq_col="frequency",
                gate=args.gate,
                top_freq=args.top_freq,
                top_prism=args.top_prism,
                tie_break=tie_break,
            )
            scatter_png = Path(args.output).with_name("freq_prism_scatter.png")
            if plot_freq_prism_scatter(
                cand, scatter_png, condition_col=cond_col, gate=args.gate
            ):
                print(f"Wrote freq×PRISM scatter → {scatter_png}")
            selected_ids = set(
                sel.loc[sel["selection_route"].isin(["prism", "both"]), "CDR3ab"]
            )
            comp_png = Path(args.output).with_name("prism_components.png")
            if plot_prism_components(clone_scores, selected_ids, comp_png):
                print(f"Wrote PRISM-component distributions → {comp_png}")
        except Exception as exc:  # pragma: no cover - plotting is best-effort
            print(f"  (skipped selection context plots: {exc})")
        if getattr(args, "prism_grid", False):
            from ..selection import freq_prism_grid

            grid = freq_prism_grid(
                feat, condition_col=cond_col, freq_col="frequency", gate=args.gate
            )
            grid_csv = Path(args.output).with_name("freq_prism_grid.csv")
            grid.to_csv(grid_csv, index=False)
            print(f"Wrote freq×PRISM trade-off grid → {grid_csv}")
            try:
                from ..plots import plot_freq_prism_grid

                grid_png = Path(args.output).with_name("freq_prism_grid.png")
                plot_freq_prism_grid(
                    grid, grid_png, chosen=(args.top_freq, args.top_prism)
                )
                print(f"Wrote freq×PRISM heatmap → {grid_png}")
            except Exception as exc:  # pragma: no cover - plotting is best-effort
                print(f"  (skipped heatmap: {exc})")
        return

    selection = {}
    if args.config:
        selection = TCRsiftConfig.from_yaml(args.config).selection or {}
    if not selection.get("rules"):
        raise TCRsiftValidationError(
            "No selection rules configured.",
            hint="Pass --config with a 'selection:' block ('rules'), or use "
            "--prism for the per-condition freq∪PRISM recipe.",
        )

    rules = select_from_clone_sample_long(long_df, selection, clone_scores=clone_scores)
    if rules.empty:
        print("  No clones selected (no method axis, no rules matched, or empty input).")
    elif args.clonotypes:
        clo = pd.read_csv(args.clonotypes)
        if "CDR3ab" in clo.columns:
            rules = (
                rules.merge(clo, on="CDR3ab", how="left")
                .sort_values("global_rank")
                .reset_index(drop=True)
            )

    rules.to_csv(args.output, index=False)
    n_rules = rules["selection_rule"].nunique() if "selection_rule" in rules.columns else 0
    print(f"Wrote {len(rules)} selected clones across {n_rules} rules → {args.output}")

    if getattr(args, "emit_per_sample_tiers", False):
        tiered = attach_per_sample_tiers(long_df)
        tier_path = Path(args.output).with_name("clone_sample_per_sample_tiers.csv")
        tiered.to_csv(tier_path, index=False)
        print(f"Wrote per-sample tiers → {tier_path}")


def cmd_unify(args):
    """Unify clonotype data from multiple experiments."""
    import pandas as pd

    from ..unify import (
        add_phenotype_confidence,
        find_top_condition,
        get_unify_summary,
        merge_experiments,
    )

    setup_logging(args.verbose)

    experiments = []
    for path in args.inputs:
        df = pd.read_csv(path)
        name = Path(path).stem
        experiments.append((df, name))
        print(f"Loaded {len(df):,} clonotypes from {path}")

    merged = merge_experiments(
        experiments,
        add_occurrence_flags=args.add_occurrence_flags,
        add_combined_stats=args.add_combined_stats,
        verbose=True,
    )

    if args.add_phenotype_confidence:
        merged = add_phenotype_confidence(
            merged,
            ratio_threshold=args.phenotype_ratio_threshold,
            verbose=True,
        )

    if args.conditions:
        conditions = [c.strip() for c in args.conditions.split(",")]
        merged = find_top_condition(merged, conditions, verbose=True)

    summary = get_unify_summary(merged)
    print("\nUnify Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value:,}" if isinstance(value, (int, float)) else f"  {key}: {value}")

    merged.to_csv(args.output, index=False)
    print(f"\nSaved {len(merged):,} unified clonotypes to {args.output}")


def cmd_cohort(args):
    """Run multi-donor clone-overlap analysis (#125)."""
    from ..cohort import run_cohort_analysis

    donor_dirs: dict[str, str] = {}
    for spec in args.donor or []:
        if "=" not in spec:
            print(f"ERROR: --donor must be NAME=DIR, got {spec!r}")
            sys.exit(1)
        name, ddir = spec.split("=", 1)
        donor_dirs[name.strip()] = ddir.strip()
    if len(donor_dirs) < 2:
        print("ERROR: cohort needs >=2 donors (--donor NAME=DIR, repeatable)")
        sys.exit(1)

    matrices = run_cohort_analysis(
        donor_dirs,
        args.output_dir,
        selected_only=args.selected_only,
        include_beta_only=not args.no_beta_only,
        emit_tables=not args.no_tables,
        emit_plots=args.emit_plots,
    )
    print(
        f"Cohort overlap across {len(donor_dirs)} donors "
        f"({', '.join(donor_dirs)}) -> {args.output_dir}"
    )
    for name, mat in matrices.items():
        print(f"  cohort_{name}: {mat.shape[0]}x{mat.shape[1]}")


def register_filter(subparsers) -> None:
    """Register the ``filter`` subcommand."""
    p_filter = subparsers.add_parser("filter", help="Filter clonotypes with tiered confidence")
    p_filter.add_argument("--input", "-i", required=True, help="Input clonotypes CSV")
    p_filter.add_argument("--output", "-o", required=True, help="Output directory for tier CSVs")
    p_filter.add_argument(
        "--tcell-type", choices=["cd8", "cd4", "both"], default="cd8", help="T cell type filter"
    )
    p_filter.add_argument(
        "--method",
        choices=["threshold", "logistic", "stimulation"],
        default="threshold",
        help="Filtering method ('threshold', 'logistic', or 'stimulation')",
    )
    p_filter.add_argument(
        "--stim-fdr-threshold",
        type=float,
        default=0.1,
        help="Max Fisher exact FDR for method=stimulation expansion branch (default: 0.1)",
    )
    p_filter.add_argument(
        "--stim-log2-fc-threshold",
        type=float,
        default=1.0,
        help="Min log2 fold-change over control for method=stimulation expansion branch (default: 1.0)",
    )
    p_filter.add_argument(
        "--stim-delta-threshold",
        type=float,
        default=0.5,
        help="Min within-clone delta signature score (e.g. delta_signature_AIM) for 16 h AIM stimulation branch (default: 0.5)",
    )
    p_filter.add_argument("--min-cells", type=int, default=2, help="Min cells per clone")
    p_filter.add_argument("--min-frequency", type=float, default=0.0, help="Min frequency")
    p_filter.add_argument(
        "--require-complete", action="store_true", default=True, help="Require complete TCR"
    )
    p_filter.add_argument("--no-require-complete", dest="require_complete", action="store_false")
    p_filter.add_argument(
        "--fdr-tiers", default="0.15,0.1,0.01,0.001,0.0001", help="FDR tiers (comma-separated)"
    )
    p_filter.add_argument("--exclude-viral", action="store_true", help="Exclude viral clones")
    p_filter.add_argument(
        "--emit-per-sample-tier",
        action="store_true",
        help="Add tier_in_{sample} columns to all_filtered.csv (needs "
        "--clone-sample-long); feeds the `tcrsift select` rule language (#122)",
    )
    p_filter.add_argument(
        "--clone-sample-long",
        help="clone_sample_long.csv to source per-sample tiers from (#122)",
    )
    p_filter.add_argument(
        "--filter-mode",
        choices=["fdr", "shared-high-freq", "cross-donor-public"],
        default="fdr",
        help="Named filter preset (default: fdr). 'shared-high-freq' applies "
        "min-methods-per-donor=2 and min-frequency-per-method=0.01. "
        "'cross-donor-public' adds min-donors=2.",
    )
    p_filter.add_argument(
        "--min-donors",
        type=int,
        default=0,
        help="Min distinct donors clone must appear in (#8/#15)",
    )
    p_filter.add_argument(
        "--min-methods-per-donor",
        type=int,
        default=0,
        help="Min distinct enrichment methods within at least one donor (#8/#15)",
    )
    p_filter.add_argument(
        "--min-cells-per-method",
        type=int,
        default=0,
        help="Min cells in at least one enrichment method (#15)",
    )
    p_filter.add_argument(
        "--min-frequency-per-method",
        type=float,
        default=0.0,
        help="Min frequency within a single enrichment method (#15)",
    )
    p_filter.add_argument(
        "--min-timepoints",
        type=int,
        default=0,
        help="Min distinct timepoints clone must appear in (#9)",
    )
    p_filter.add_argument(
        "--min-timepoints-per-donor",
        type=int,
        default=0,
        help="Min distinct timepoints within at least one donor (#9)",
    )
    p_filter.add_argument(
        "--min-apcs",
        type=int,
        default=0,
        help="Min distinct APCs clone must appear across (#9)",
    )
    p_filter.add_argument(
        "--min-apcs-per-donor",
        type=int,
        default=0,
        help="Min distinct APCs within at least one donor (#9)",
    )
    p_filter.add_argument(
        "--min-til-cells-per-donor",
        type=int,
        default=0,
        help="Min TIL cells of the same donor in which the clone appears (#9)",
    )
    p_filter.add_argument("--plot-filter", action="store_true", help="Generate filter plots")
    p_filter.add_argument("--output-dir", help="Output directory for plots")
    p_filter.add_argument("--verbose", action="store_true", help="Verbose output")
    p_filter.set_defaults(func=cmd_filter)


def register_select(subparsers) -> None:
    """Register the ``select`` subcommand."""
    p_select = subparsers.add_parser(
        "select",
        help="Apply config-driven selection rules to a clone-sample-long table",
    )
    p_select.add_argument(
        "--input",
        "-i",
        required=True,
        help="clone_sample_long.csv (one row per clone×sample with method/freq)",
    )
    p_select.add_argument("--output", "-o", required=True, help="Output selected-clones CSV")
    p_select.add_argument(
        "--config",
        "-c",
        help="YAML config with a 'selection:' block (rules + global_rank). "
        "Not needed with --prism.",
    )
    p_select.add_argument(
        "--prism",
        action="store_true",
        help="Per-condition top-K freq ∪ top-K PRISM recipe (needs --clone-scores), "
        "matching the validated selection_analysis (#193)",
    )
    p_select.add_argument(
        "--condition-col",
        help="Condition column for --prism (default: method, else sample)",
    )
    p_select.add_argument(
        "--gate",
        type=float,
        default=0.001,
        help="--prism frequency gate (default: 0.001)",
    )
    p_select.add_argument(
        "--top-freq",
        type=int,
        default=10,
        help="--prism top-K by frequency (default: 10)",
    )
    p_select.add_argument(
        "--top-prism",
        type=int,
        default=5,
        help="--prism top-K by PRISM (default: 5)",
    )
    p_select.add_argument(
        "--tie-break",
        help="--prism: comma-separated tie-break keys for boundary ties "
        "(default: single_alpha,umi,cd8_purity,cdr3ab; use 'cdr3ab' for "
        "CDR3ab-only). Precedence: this flag > config selection.tie_break > "
        "default (#229).",
    )
    p_select.add_argument(
        "--rescue-target",
        type=int,
        default=0,
        help="--prism: when a condition has fewer than N gated clones, rescue "
        "the best sub-threshold clones by --rescue-rank-col so thin conditions "
        "still get a meaningful top-K (default: 0 = off) (#228).",
    )
    p_select.add_argument(
        "--rescue-rank-col",
        help="--prism: column (from --clone-scores, e.g. a coverage/UMI metric) "
        "ranking sub-threshold clones for rescue (descending). Required when "
        "--rescue-target > 0.",
    )
    p_select.add_argument(
        "--prism-grid",
        action="store_true",
        help="With --prism: also sweep a (top_freq × top_prism) grid and write "
        "the distinct-clone trade-off CSV + heatmap PNG (#207)",
    )
    p_select.add_argument(
        "--clonotypes",
        help="Optional clonotypes CSV to join per-clone metadata onto the output "
        "(also used as per-clone scores for a percentile_rank route if it has them)",
    )
    p_select.add_argument(
        "--clone-scores",
        help="Per-clone scores CSV (ppost_alpha/beta, GEX gene-set scores) for a "
        "route's rank_by: {percentile_rank: [...]} — reproduces PRISM (#193)",
    )
    p_select.add_argument(
        "--emit-per-sample-tiers",
        action="store_true",
        help="Also write per-sample abundance tiers for the long table",
    )
    p_select.add_argument("--verbose", action="store_true", help="Verbose output")
    p_select.set_defaults(func=cmd_select)


def register_unify(subparsers) -> None:
    """Register the ``unify`` subcommand."""
    p_unify = subparsers.add_parser("unify", help="Unify clonotypes from multiple experiments")
    p_unify.add_argument(
        "--inputs", "-i", nargs="+", required=True, help="Input CSV files to merge"
    )
    p_unify.add_argument("--output", "-o", required=True, help="Output unified CSV")
    p_unify.add_argument(
        "--add-occurrence-flags",
        action="store_true",
        default=True,
        help="Add 'occurs_in_*' columns",
    )
    p_unify.add_argument("--no-occurrence-flags", dest="add_occurrence_flags", action="store_false")
    p_unify.add_argument(
        "--add-combined-stats", action="store_true", default=True, help="Add combined statistics"
    )
    p_unify.add_argument("--no-combined-stats", dest="add_combined_stats", action="store_false")
    p_unify.add_argument(
        "--add-phenotype-confidence",
        action="store_true",
        default=True,
        help="Add phenotype confidence columns",
    )
    p_unify.add_argument(
        "--no-phenotype-confidence", dest="add_phenotype_confidence", action="store_false"
    )
    p_unify.add_argument(
        "--phenotype-ratio-threshold",
        type=float,
        default=10.0,
        help="CD4/CD8 ratio for confident classification (default: 10.0)",
    )
    p_unify.add_argument(
        "--conditions", help="Comma-separated condition names for top-condition analysis"
    )
    p_unify.add_argument("--verbose", action="store_true", help="Verbose output")
    p_unify.set_defaults(func=cmd_unify)


def register_cohort(subparsers) -> None:
    """Register the ``cohort`` subcommand."""
    p_cohort = subparsers.add_parser(
        "cohort", help="Multi-donor clone-overlap analysis across run outputs"
    )
    p_cohort.add_argument(
        "--donor",
        action="append",
        metavar="NAME=DIR",
        help="Donor as NAME=run_output_dir (repeatable; need >=2).",
    )
    p_cohort.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Directory for cohort overlap tables / plots.",
    )
    p_cohort.add_argument(
        "--selected-only",
        action="store_true",
        help="Overlap selected_clones.csv instead of clonotypes.csv.",
    )
    p_cohort.add_argument(
        "--no-beta-only",
        action="store_true",
        help="Skip the beta-only overlap matrices.",
    )
    p_cohort.add_argument(
        "--no-tables",
        action="store_true",
        help="Don't write CSV matrices.",
    )
    p_cohort.add_argument(
        "--emit-plots",
        action="store_true",
        help="Render heatmap PDFs.",
    )
    p_cohort.set_defaults(func=cmd_cohort)
