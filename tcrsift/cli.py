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

"""
Command-line interface for TCRsift.

TCRsift: TCR selection from antigen-specific culture and scRNA/VDJ sequencing data.
"""

import argparse
import logging
import sys
from pathlib import Path

from .config import add_config_args, load_config_with_args
from .validation import (
    TCRsiftValidationError,
    validate_annotate_gex_args,
    validate_assemble_args,
    validate_file_exists,
    validate_run_args,
)
from .version import __version__


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# =============================================================================
# Load Command
# =============================================================================


def cmd_load(args):
    """Load CellRanger outputs from sample sheet."""
    from .loader import load_samples
    from .plots import plot_qc

    setup_logging(args.verbose)

    adata = load_samples(
        args.sample_sheet,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        min_counts=args.min_counts,
        max_counts=args.max_counts,
        max_mito_pct=args.max_mito,
        min_mito_pct=args.min_mito,
    )

    # Save
    adata.write_h5ad(args.output)
    print(f"Saved {len(adata)} cells to {args.output}")

    # Generate plots if requested
    if args.plot_qc:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_qc(adata, output_dir)
        print(f"QC plots saved to {output_dir}")


# =============================================================================
# Phenotype Command
# =============================================================================


def cmd_phenotype(args):
    """Classify T cells as CD4/CD8."""
    import anndata as ad

    from .phenotype import get_phenotype_summary, phenotype_cells, validate_phenotype_vs_expected
    from .plots import plot_phenotype

    setup_logging(args.verbose)

    adata = ad.read_h5ad(args.input)
    print(f"Loaded {len(adata)} cells from {args.input}")

    adata = phenotype_cells(
        adata,
        cd4_cd8_ratio=args.cd4_cd8_ratio,
        min_cd3_reads=args.min_cd3_reads,
    )

    # Validate against expected
    warnings = validate_phenotype_vs_expected(adata)
    for w in warnings:
        print(f"WARNING: {w}")

    # Print summary
    summary = get_phenotype_summary(adata)
    print("\nPhenotype Summary:")
    print(summary.to_string())

    # Save
    adata.write_h5ad(args.output)
    print(f"\nSaved phenotyped data to {args.output}")

    # Generate plots if requested
    if args.plot_phenotype:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_phenotype(adata, output_dir)
        print(f"Phenotype plots saved to {output_dir}")


# =============================================================================
# Clonotype Command
# =============================================================================


def cmd_clonotype(args):
    """Aggregate cells into clonotypes."""
    import anndata as ad

    from .clonotype import aggregate_clonotypes, export_clonotypes_airr, get_clonotype_summary
    from .plots import plot_clonotypes

    setup_logging(args.verbose)

    adata = ad.read_h5ad(args.input)
    print(f"Loaded {len(adata)} cells from {args.input}")

    clonotypes = aggregate_clonotypes(
        adata,
        group_by=args.group_by,
        min_umi=args.min_umi,
        handle_doublets=args.handle_doublets,
    )

    # Print summary
    summary = get_clonotype_summary(clonotypes)
    print("\nClonotype Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    # Save
    clonotypes.to_csv(args.output, index=False)
    print(f"\nSaved {len(clonotypes)} clonotypes to {args.output}")

    # AIRR format if requested
    if args.airr:
        export_clonotypes_airr(clonotypes, args.airr)
        print(f"Saved AIRR format to {args.airr}")

    # Generate plots if requested
    if args.plot_clonotypes:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_clonotypes(clonotypes, output_dir)
        print(f"Clonotype plots saved to {output_dir}")


# =============================================================================
# Filter Command
# =============================================================================


def cmd_filter(args):
    """Filter clonotypes with tiered confidence levels."""
    import logging

    import pandas as pd

    from .filter import (
        filter_clonotypes,
        get_filter_summary,
        resolve_filter_mode_kwargs,
        split_by_tier,
    )
    from .plots import plot_filter

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    # Parse FDR tiers
    fdr_tiers = None
    if args.fdr_tiers:
        fdr_tiers = [float(x) for x in args.fdr_tiers.split(",")]

    # Resolve named filter mode preset, then layer user-supplied knobs.
    mode = getattr(args, "filter_mode", "fdr")
    user_kwargs = {
        "min_donors": args.min_donors,
        "min_methods_per_donor": args.min_methods_per_donor,
        "min_cells_per_method": args.min_cells_per_method,
        "min_frequency_per_method": args.min_frequency_per_method,
    }
    mode_kwargs = resolve_filter_mode_kwargs(mode, user_kwargs)

    # `cross-donor-public` only makes sense for multi-donor cohorts. Per
    # #15: warn rather than error when n_donors == 1 across all clones.
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
        **mode_kwargs,
    )

    # Print summary
    summary = get_filter_summary(filtered)
    print("\nFilter Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    # Save by tier
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    tier_dfs = split_by_tier(filtered)
    for tier, tier_df in tier_dfs.items():
        tier_path = output_dir / f"{tier}.csv"
        tier_df.to_csv(tier_path, index=False)
        print(f"Saved {len(tier_df)} {tier} clonotypes to {tier_path}")

    # Also save combined
    filtered.to_csv(output_dir / "all_filtered.csv", index=False)

    # Mode-named output (#15 chunk 4)
    if mode == "shared-high-freq":
        filtered.to_csv(output_dir / "all_shared_high_freq.csv", index=False)
    elif mode == "cross-donor-public":
        filtered.to_csv(output_dir / "all_cross_donor_public.csv", index=False)

    # Public/private bucketing (#15 chunk 4)
    from .filter import bucket_by_donor_sharing

    sharing_buckets = bucket_by_donor_sharing(filtered)
    for bucket, bucket_df in sharing_buckets.items():
        bpath = output_dir / f"{bucket}.csv"
        bucket_df.to_csv(bpath, index=False)
        print(f"Saved {len(bucket_df)} {bucket} clonotypes to {bpath}")

    # Generate plots if requested
    if args.plot_filter:
        plot_dir = Path(args.output_dir) if args.output_dir else output_dir / "plots"
        plot_filter(filtered, plot_dir)
        print(f"Filter plots saved to {plot_dir}")


# =============================================================================
# Annotate Command
# =============================================================================


def cmd_annotate(args):
    """Annotate clonotypes with public database matches."""
    import pandas as pd

    from .annotate import annotate_clonotypes, get_annotation_summary
    from .plots import plot_annotations

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    annotated = annotate_clonotypes(
        clonotypes,
        vdjdb_path=args.vdjdb,
        iedb_path=args.iedb,
        cedar_path=args.cedar,
        match_by=args.match_by if hasattr(args, "match_by") else "CDR3ab",
        exclude_viral=args.exclude_viral,
        flag_only=args.flag_only,
    )

    # Print summary
    summary = get_annotation_summary(annotated)
    print("\nAnnotation Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    # Save
    annotated.to_csv(args.output, index=False)
    print(f"\nSaved annotated clonotypes to {args.output}")

    # Generate plots if requested
    if args.plot_annotations:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_annotations(annotated, output_dir)
        print(f"Annotation plots saved to {output_dir}")


# =============================================================================
# Match-TIL Command
# =============================================================================


def cmd_match_til(args):
    """Match culture clonotypes against TIL data."""
    import pandas as pd

    from .plots import plot_til
    from .sample_sheet import load_sample_sheet
    from .til import get_til_summary, load_til_data, load_til_samples, load_til_specs, match_til
    from .validation import TCRsiftValidationError, validate_match_til_args

    setup_logging(args.verbose)
    validate_match_til_args(args)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} culture clonotypes from {args.input}")

    # Determine TIL data source - multiple options supported
    til_data = None

    if args.til_sample:
        til_data = load_til_specs(args.til_sample)
        total_cells = sum(len(df) for df in til_data.values())
        print(
            f"Loaded {total_cells} TIL cells from {len(til_data)} sample(s): {list(til_data.keys())}"
        )

    elif args.sample_sheet:
        # Load TIL samples from sample sheet
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

    elif args.til_h5ad:
        # Single h5ad file
        til_data = load_til_data("h5ad", args.til_h5ad)
        print(f"Loaded {len(til_data)} TIL cells from h5ad: {args.til_h5ad}")

    elif args.til_csv:
        # Single CSV file
        til_data = load_til_data("csv", args.til_csv)
        print(f"Loaded {len(til_data)} TIL cells from CSV: {args.til_csv}")

    elif args.til_vdj_dir:
        # Single VDJ directory
        til_data = load_til_data("vdj_dir", args.til_vdj_dir)
        print(f"Loaded {len(til_data)} TIL cells from VDJ directory: {args.til_vdj_dir}")

    else:
        raise TCRsiftValidationError(
            "No TIL data source specified",
            hint="Provide one of: --sample-sheet, --til-h5ad, --til-csv, --til-vdj-dir, or --til-sample",
        )

    matched = match_til(
        clonotypes,
        til_data,
        match_by=args.match_by,
        min_til_cells=args.min_til_cells,
    )

    # Print summary
    summary = get_til_summary(matched)
    print("\nTIL Matching Summary:")
    for key, value in summary.items():
        if isinstance(value, dict):
            print(f"  {key}:")
            for k, v in value.items():
                print(f"    {k}: {v}")
        else:
            print(f"  {key}: {value}")

    # Save
    matched.to_csv(args.output, index=False)
    print(f"\nSaved TIL-matched clonotypes to {args.output}")

    # Generate plots if requested
    if args.plot_til:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_til(matched, output_dir)
        print(f"TIL plots saved to {output_dir}")


# =============================================================================
# TIL-Clonotype Command
# =============================================================================


def cmd_til_clonotype(args):
    """Aggregate TIL-only data into clonotype counts across samples."""
    from .sample_sheet import load_sample_sheet
    from .til import load_til_data, load_til_samples, load_til_specs, summarize_til_clonotypes
    from .validation import TCRsiftValidationError, validate_match_til_args

    setup_logging(args.verbose)
    validate_match_til_args(args)

    til_data = None

    if args.til_sample:
        til_data = load_til_specs(args.til_sample)
        total_cells = sum(len(df) for df in til_data.values())
        print(
            f"Loaded {total_cells} TIL cells from {len(til_data)} sample(s): {list(til_data.keys())}"
        )
    elif args.sample_sheet:
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
    elif args.til_h5ad:
        til_data = load_til_data("h5ad", args.til_h5ad)
        print(f"Loaded {len(til_data)} TIL cells from h5ad: {args.til_h5ad}")
    elif args.til_csv:
        til_data = load_til_data("csv", args.til_csv)
        print(f"Loaded {len(til_data)} TIL cells from CSV: {args.til_csv}")
    elif args.til_vdj_dir:
        til_data = load_til_data("vdj_dir", args.til_vdj_dir)
        print(f"Loaded {len(til_data)} TIL cells from VDJ directory: {args.til_vdj_dir}")
    else:
        raise TCRsiftValidationError(
            "No TIL data source specified",
            hint="Provide one of: --sample-sheet, --til-h5ad, --til-csv, --til-vdj-dir, or --til-sample",
        )

    clonotypes = summarize_til_clonotypes(
        til_data,
        match_by=args.match_by,
        min_cells=args.min_cells,
    )

    clonotypes.to_csv(args.output, index=False)
    print(f"Saved {len(clonotypes)} TIL clonotypes to {args.output}")


# =============================================================================
# TIL-Select Command
# =============================================================================


def cmd_til_select(args):
    """Prioritize promising TIL clonotypes from multi-timepoint VDJ+GEX data."""
    from .til_select import run_til_select

    setup_logging(args.verbose)
    master_df = run_til_select(args)
    n_final = int(master_df["is_candidate_tumor_reactive"].sum())
    print(f"Saved harmonized table: {args.out_table}")
    print(f"Final candidate clonotypes: {n_final}")
    print(f"Figures/subsets directory: {args.fig_dir}")


# =============================================================================
# Assemble Command
# =============================================================================


def cmd_assemble(args):
    """Assemble full-length TCR sequences."""
    import pandas as pd

    from .assemble import assemble_full_sequences, export_fasta, validate_sequences
    from .plots import plot_assembly

    setup_logging(args.verbose)

    # Early validation of conditional requirements
    validate_assemble_args(args)
    validate_file_exists(args.input, "input clonotypes CSV")

    # Handle leader shortcuts
    alpha_leader = args.alpha_leader
    beta_leader = args.beta_leader

    if getattr(args, "no_leaders", False):
        alpha_leader = None
        beta_leader = None
    elif getattr(args, "leaders_from_contigs", False):
        alpha_leader = "from_contig"
        beta_leader = "from_contig"
    else:
        # Convert "none" string to None
        if alpha_leader == "none":
            alpha_leader = None
        if beta_leader == "none":
            beta_leader = None

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    assembled = assemble_full_sequences(
        clonotypes,
        contigs_dir=args.contigs_dir,
        alpha_leader=alpha_leader,
        beta_leader=beta_leader,
        include_constant=args.include_constant,
        constant_source=args.constant_source,
        linker=args.linker if args.single_chain else None,
    )

    # Validate
    warnings = validate_sequences(assembled)
    for w in warnings:
        print(f"WARNING: {w}")

    # Save
    assembled.to_csv(args.output, index=False)
    print(f"\nSaved assembled sequences to {args.output}")

    # AIRR format if requested
    if args.airr:
        from .clonotype import export_clonotypes_airr

        export_clonotypes_airr(assembled, args.airr)
        print(f"Saved AIRR format to {args.airr}")

    # FASTA if requested
    if args.fasta:
        seq_col = "single_chain_aa" if args.single_chain else "full_beta_aa"
        export_fasta(assembled, args.fasta, sequence_col=seq_col)
        print(f"Saved FASTA to {args.fasta}")

    # Generate plots if requested
    if args.plot_assembly:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_assembly(assembled, output_dir)
        print(f"Assembly plots saved to {output_dir}")


# =============================================================================
# Run Command (Unified Pipeline)
# =============================================================================


def cmd_run(args):
    """Run the complete TCRsift pipeline."""
    from pathlib import Path

    from .annotate import annotate_clonotypes
    from .assemble import assemble_full_sequences, export_fasta
    from .clonotype import aggregate_clonotypes, export_clonotypes_airr
    from .filter import filter_clonotypes, split_by_tier
    from .loader import load_samples
    from .phenotype import filter_by_tcell_type, phenotype_cells
    from .plots import (
        create_pipeline_funnel,
        create_tcr_sequence_pdf,
        generate_report,
        plot_annotations,
        plot_assembly,
        plot_clonotypes,
        plot_filter,
        plot_phenotype,
        plot_qc,
        plot_til,
    )
    from .sample_sheet import load_sample_sheet
    from .til import load_til_specs, match_til

    # Early validation of arguments
    validate_run_args(args)

    # Load config with CLI overrides
    config = load_config_with_args(args)

    # Load sample sheet once (used for auto-detect + TIL loading)
    sample_sheet = load_sample_sheet(args.sample_sheet)

    # Auto-detect TIL samples from sample sheet if not explicitly specified
    if not config.til.til_samples:
        til_samples = sample_sheet.get_til_samples()
        if til_samples:
            til_sample_names = [s.sample for s in til_samples]
            config.til.til_samples = til_sample_names
            print(f"Auto-detected TIL samples from sample sheet: {til_sample_names}")
    setup_logging(config.verbose)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Save config for reproducibility
    config.to_yaml(output_dir / "config.yaml")

    print("=" * 60)
    print("TCRsift Pipeline")
    print("=" * 60)

    # Track counts for funnel plot
    funnel_counts = {}

    # Step 1: Load
    print("\n[1/7] Loading samples...")
    adata = load_samples(
        args.sample_sheet,
        min_genes=config.load.min_genes,
        max_genes=config.load.max_genes,
        min_counts=config.load.min_counts,
        max_counts=config.load.max_counts,
        min_mito_pct=config.load.min_mito_pct,
        max_mito_pct=config.load.max_mito_pct,
    )
    adata.write_h5ad(data_dir / "loaded.h5ad")
    funnel_counts["Raw Cells"] = len(adata)
    print(f"  Loaded {len(adata)} cells")

    # Exclude TIL samples from the main culture pipeline
    if "source" in adata.obs.columns:
        non_til_mask = adata.obs["source"] != "til"
        n_til = (~non_til_mask).sum()
        if n_til > 0:
            if non_til_mask.sum() == 0:
                raise TCRsiftValidationError(
                    "No non-TIL samples found for culture pipeline",
                    hint="The run command expects at least one non-TIL sample. "
                    "Use match-til or load/clonotype commands for TIL-only analysis.",
                )
            adata = adata[non_til_mask].copy()
            funnel_counts["Raw Cells"] = len(adata)
            print(f"  Excluding {n_til} TIL cells from culture pipeline")

    # Count cells with VDJ
    if "CDR3_beta" in adata.obs.columns:
        n_with_vdj = adata.obs["CDR3_beta"].notna().sum()
        funnel_counts["With VDJ"] = n_with_vdj
        print(f"  With VDJ data: {n_with_vdj} cells")

    if config.output.generate_plots:
        plot_qc(adata, plots_dir)

    # Step 2: Phenotype
    print("\n[2/7] Phenotyping T cells...")
    adata = phenotype_cells(
        adata,
        cd4_cd8_ratio=config.phenotype.cd4_cd8_ratio,
        min_cd3_reads=config.phenotype.min_cd3_reads,
    )
    adata.write_h5ad(data_dir / "phenotyped.h5ad")
    if config.output.generate_plots:
        plot_phenotype(adata, plots_dir)

    # Filter by T cell type if specified
    tcell_type = config.filter.tcell_type
    if tcell_type != "both":
        adata = filter_by_tcell_type(adata, tcell_type)
        print(f"  Filtered to {tcell_type.upper()}+: {len(adata)} cells")

    funnel_counts["Phenotyped"] = len(adata)

    # Step 3: Clonotype
    print("\n[3/7] Aggregating clonotypes...")
    clonotypes = aggregate_clonotypes(
        adata,
        group_by=config.clonotype.group_by,
        min_umi=config.clonotype.min_umi,
        handle_doublets=config.clonotype.handle_doublets,
    )
    clonotypes.to_csv(data_dir / "clonotypes.csv", index=False)
    funnel_counts["Clonotypes"] = len(clonotypes)
    print(f"  Found {len(clonotypes)} clonotypes")
    if config.output.generate_plots:
        plot_clonotypes(clonotypes, plots_dir)

    # Step 4: Filter
    print("\n[4/7] Filtering clonotypes...")
    # Resolve named filter mode preset, then layer user-supplied knobs.
    from .filter import resolve_filter_mode_kwargs

    user_filter_kwargs = {
        "min_donors": config.filter.min_donors,
        "min_methods_per_donor": config.filter.min_methods_per_donor,
        "min_cells_per_method": config.filter.min_cells_per_method,
        "min_frequency_per_method": config.filter.min_frequency_per_method,
    }
    mode_kwargs = resolve_filter_mode_kwargs(
        config.filter.filter_mode, user_filter_kwargs
    )

    if config.filter.filter_mode == "cross-donor-public" and "n_donors" in clonotypes.columns:
        if (clonotypes["n_donors"] >= 2).sum() == 0:
            print(
                "  WARNING: filter_mode='cross-donor-public' requested but no clone "
                "has n_donors >= 2. Result will likely be empty."
            )

    filtered = filter_clonotypes(
        clonotypes,
        method=config.filter.method,
        tcell_type=tcell_type,
        min_cells=config.filter.min_cells,
        min_frequency=config.filter.min_frequency,
        require_complete=config.filter.require_complete,
        fdr_tiers=config.filter.fdr_tiers,
        **mode_kwargs,
    )

    # Save by tier
    tier_dfs = split_by_tier(filtered)
    tier_counts = {}
    for tier, tier_df in tier_dfs.items():
        tier_df.to_csv(data_dir / f"filtered_{tier}.csv", index=False)
        tier_counts[tier] = len(tier_df)
        print(f"  {tier}: {len(tier_df)} clonotypes")

    # Mode-named output for non-FDR modes (#15 chunk 4). Lets users find the
    # "shared-high-freq passes" set without having to read tier CSVs that
    # don't really mean tiers under that mode.
    if config.filter.filter_mode == "shared-high-freq":
        path = data_dir / "filtered_shared_high_freq.csv"
        filtered.to_csv(path, index=False)
        print(f"  shared-high-freq: {len(filtered)} clonotypes -> {path.name}")
    elif config.filter.filter_mode == "cross-donor-public":
        path = data_dir / "filtered_cross_donor_public.csv"
        filtered.to_csv(path, index=False)
        print(f"  cross-donor-public: {len(filtered)} clonotypes -> {path.name}")

    # Public/private bucketing (#15 chunk 4). Written when n_donors info is
    # present. Empty buckets are skipped so single-donor cohorts don't get
    # a 0-row public_across_donors.csv.
    from .filter import bucket_by_donor_sharing

    sharing_buckets = bucket_by_donor_sharing(filtered)
    for bucket, bucket_df in sharing_buckets.items():
        path = data_dir / f"filtered_{bucket}.csv"
        bucket_df.to_csv(path, index=False)
        print(f"  {bucket}: {len(bucket_df)} clonotypes -> {path.name}")

    funnel_counts["Filtered"] = len(filtered)

    if config.output.generate_plots:
        plot_filter(filtered, plots_dir)

    # Step 5: Annotate (if databases provided)
    annotated = filtered
    has_annotation = (
        config.annotate.vdjdb_path or config.annotate.iedb_path or config.annotate.cedar_path
    )
    if has_annotation:
        print("\n[5/7] Annotating with public databases...")
        annotated = annotate_clonotypes(
            filtered,
            vdjdb_path=config.annotate.vdjdb_path,
            iedb_path=config.annotate.iedb_path,
            cedar_path=config.annotate.cedar_path,
            match_by=config.annotate.match_by,
            exclude_viral=config.annotate.exclude_viral,
            flag_only=config.annotate.flag_only,
        )
        annotated.to_csv(data_dir / "annotated.csv", index=False)
        n_viral = annotated["is_viral"].sum() if "is_viral" in annotated.columns else 0
        print(f"  Flagged {n_viral} viral clonotypes")
        if config.output.generate_plots:
            plot_annotations(annotated, plots_dir)
    else:
        print("\n[5/7] Skipping annotation (no databases provided)")

    # Step 6: TIL matching (if TIL samples specified)
    til_matched = annotated
    if getattr(args, "til_sample", None):
        print("\n[6/7] Matching against TIL samples from --til-sample...")
        til_data = load_til_specs(args.til_sample)
        total_cells = sum(len(df) for df in til_data.values())
        print(f"  Loaded {total_cells} TIL cells from {len(til_data)} sample(s)")

        if total_cells > 0:
            til_matched = match_til(
                annotated,
                til_data,
                match_by=config.til.match_by,
                min_til_cells=config.til.min_til_cells,
            )
            til_matched.to_csv(data_dir / "til_matched.csv", index=False)
            n_til = til_matched["til_match"].sum() if "til_match" in til_matched.columns else 0
            print(f"  Found {n_til} clonotypes in TILs")
            if config.output.generate_plots:
                plot_til(til_matched, plots_dir)
        else:
            print("  No TIL cells loaded")
    elif config.til.til_samples:
        print("\n[6/7] Matching against TIL samples...")
        from .til import load_til_samples as _load_til_samples

        # Get TIL samples from sample sheet
        til_sample_objs = sample_sheet.get_til_samples()
        til_sample_names = [s.sample for s in til_sample_objs]

        # Filter to only requested samples if config specifies specific names
        if config.til.til_samples != til_sample_names:
            til_sample_objs = [s for s in til_sample_objs if s.sample in config.til.til_samples]

        if til_sample_objs:
            # Use new flexible TIL loading (supports h5ad, CSV, VDJ directory)
            til_data = _load_til_samples(til_sample_objs)
            total_cells = sum(len(df) for df in til_data.values())
            print(f"  Loaded {total_cells} TIL cells from {len(til_data)} sample(s)")

            if total_cells > 0:
                til_matched = match_til(
                    annotated,
                    til_data,
                    match_by=config.til.match_by,
                    min_til_cells=config.til.min_til_cells,
                )
                til_matched.to_csv(data_dir / "til_matched.csv", index=False)
                n_til = til_matched["til_match"].sum() if "til_match" in til_matched.columns else 0
                print(f"  Found {n_til} clonotypes in TILs")
                if config.output.generate_plots:
                    plot_til(til_matched, plots_dir)
            else:
                print("  No TIL cells loaded")
        else:
            print("  No TIL samples found matching configuration")
    else:
        print("\n[6/7] Skipping TIL matching (no TIL samples specified)")

    # Step 7: Assemble (if requested)
    assembled = None
    has_leaders = (
        config.assemble.alpha_leader is not None or config.assemble.beta_leader is not None
    )
    if config.assemble.single_chain or has_leaders or config.assemble.include_constant:
        print("\n[7/7] Assembling full-length sequences...")
        assembled = assemble_full_sequences(
            til_matched,
            contigs_dir=config.assemble.contigs_dir,
            alpha_leader=config.assemble.alpha_leader,
            beta_leader=config.assemble.beta_leader,
            include_constant=config.assemble.include_constant,
            constant_source=config.assemble.constant_source,
            linker=config.assemble.linker if config.assemble.single_chain else None,
        )
        assembled.to_csv(data_dir / "full_sequences.csv", index=False)
        print(f"  Assembled {len(assembled)} sequences")
        if config.output.generate_plots:
            plot_assembly(assembled, plots_dir)

        # Generate sequence PDF
        if config.output.generate_report:
            create_tcr_sequence_pdf(assembled, output_dir / "tcr_sequences.pdf")

        if config.output.output_airr:
            airr_path = data_dir / "full_sequences.airr.tsv"
            export_clonotypes_airr(assembled, str(airr_path))
            print(f"  Exported AIRR table: {airr_path}")

        if config.output.output_fasta:
            fasta_path = data_dir / "full_sequences.fasta"
            seq_col = "single_chain_aa" if config.assemble.single_chain else "full_beta_aa"
            export_fasta(assembled, fasta_path, sequence_col=seq_col)
            print(f"  Exported FASTA: {fasta_path}")
    else:
        print("\n[7/7] Skipping assembly")
        if config.output.output_airr or config.output.output_fasta:
            print("  Skipping AIRR/FASTA export because assembly was skipped")

    # Generate funnel plot
    if config.output.generate_plots:
        create_pipeline_funnel(
            raw_cells=funnel_counts.get("Raw Cells", 0),
            with_vdj=funnel_counts.get("With VDJ", funnel_counts.get("Raw Cells", 0)),
            phenotyped=funnel_counts.get("Phenotyped", 0),
            clonotypes=funnel_counts.get("Clonotypes", 0),
            filtered=funnel_counts.get("Filtered", 0),
            tier_counts=tier_counts,
            output_dir=plots_dir,
        )

    # Generate report
    if config.output.generate_report:
        print("\nGenerating report...")
        report_format = str(config.output.report_format).lower()
        if report_format == "both":
            generate_report(plots_dir, format="html")
            generate_report(plots_dir, format="pdf")
        elif report_format in {"html", "pdf"}:
            generate_report(plots_dir, format=report_format)
        else:
            raise TCRsiftValidationError(
                f"Invalid report format: '{config.output.report_format}'",
                hint="Use report_format: pdf, html, or both in config.",
            )

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)


# =============================================================================
# Load-SCT Command
# =============================================================================


def cmd_load_sct(args):
    """Load TCR data from SCT platform."""
    from .sct import aggregate_sct, load_sct

    setup_logging(args.verbose)

    df = load_sct(
        args.input,
        sheet_name=args.sheet_name,
        min_snr=args.min_snr,
        min_reads_per_chain=args.min_reads,
        require_mutation_match=args.require_mutation_match,
        require_compact_match=args.require_compact_match,
        verbose=True,
    )

    # Aggregate if requested
    if args.aggregate:
        df = aggregate_sct(df, verbose=True)

    # Save
    df.to_csv(args.output, index=False)
    print(f"\nSaved {len(df)} {'clonotypes' if args.aggregate else 'cells'} to {args.output}")


# =============================================================================
# Annotate-GEX Command
# =============================================================================


def cmd_annotate_gex(args):
    """Annotate TCR data with gene expression from a 10x HDF5 file."""
    import pandas as pd

    from .gex import (
        DEFAULT_GENE_LIST,
        aggregate_gex_by_clonotype,
        augment_with_gex,
        compute_cd4_cd8_counts,
    )

    setup_logging(args.verbose)

    # Early validation
    validate_annotate_gex_args(args)
    validate_file_exists(args.input, "input CSV")

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df):,} rows from {args.input}")

    # Parse custom gene list if provided
    gene_list = None
    if args.genes:
        gene_list = [g.strip() for g in args.genes.split(",")]
        print(f"Using custom gene list: {gene_list}")
    else:
        print(f"Using default T cell gene list ({len(DEFAULT_GENE_LIST)} genes)")

    # Step 1: Augment with GEX if barcode column exists
    augmented_df = None
    if args.barcode_col in df.columns:
        print(f"\nAugmenting with gene expression from {args.gex_file}...")
        df = augment_with_gex(
            df,
            args.gex_file,
            barcode_col=args.barcode_col,
            gene_list=gene_list,
            col_prefix=args.prefix,
            include_qc=not args.no_qc,
            verbose=args.verbose,
        )
        # Save augmented per-cell data for CD4/CD8 counts
        augmented_df = df.copy()
    else:
        print(
            f"Warning: Barcode column '{args.barcode_col}' not found - skipping per-cell augmentation"
        )

    # Step 2: Aggregate by clonotype if requested
    if args.aggregate:
        print(f"\nAggregating by {args.group_col}...")
        df = aggregate_gex_by_clonotype(
            df,
            group_col=args.group_col,
            gex_prefix=args.prefix,
            operations=["sum", "mean"],
            verbose=args.verbose,
        )

    # Step 3: Compute CD4/CD8 counts if requested
    if args.cd4_cd8_counts:
        if augmented_df is None:
            print("Warning: Cannot compute CD4/CD8 counts without GEX augmentation - skipping")
        else:
            print("\nComputing CD4/CD8 cell counts...")
            cd4_cd8_df = compute_cd4_cd8_counts(
                augmented_df,  # Use augmented per-cell data (has gex.CD4, gex.CD8)
                group_col=args.group_col,
                gex_prefix=args.prefix,
                verbose=args.verbose,
            )
            # Merge CD4/CD8 counts into result
            if args.aggregate:
                df = df.merge(
                    cd4_cd8_df[[args.group_col, "CD4_only.count", "CD8_only.count"]],
                    on=args.group_col,
                    how="left",
                )
            else:
                # For non-aggregated output, add per-clonotype counts as additional columns
                df = df.merge(
                    cd4_cd8_df[[args.group_col, "CD4_only.count", "CD8_only.count"]],
                    on=args.group_col,
                    how="left",
                )

    # Save
    df.to_csv(args.output, index=False)
    print(f"\nSaved {len(df):,} rows to {args.output}")


# =============================================================================
# Unify Command
# =============================================================================


def cmd_unify(args):
    """Unify clonotype data from multiple experiments."""
    import pandas as pd

    from .unify import (
        add_phenotype_confidence,
        find_top_condition,
        get_unify_summary,
        merge_experiments,
    )

    setup_logging(args.verbose)

    # Load input files
    experiments = []
    for path in args.inputs:
        df = pd.read_csv(path)
        # Extract name from filename if not specified
        name = Path(path).stem
        experiments.append((df, name))
        print(f"Loaded {len(df):,} clonotypes from {path}")

    # Merge experiments
    merged = merge_experiments(
        experiments,
        add_occurrence_flags=args.add_occurrence_flags,
        add_combined_stats=args.add_combined_stats,
        verbose=True,
    )

    # Add phenotype confidence
    if args.add_phenotype_confidence:
        merged = add_phenotype_confidence(
            merged,
            ratio_threshold=args.phenotype_ratio_threshold,
            verbose=True,
        )

    # Find top condition if specified
    if args.conditions:
        conditions = [c.strip() for c in args.conditions.split(",")]
        merged = find_top_condition(merged, conditions, verbose=True)

    # Print summary
    summary = get_unify_summary(merged)
    print("\nUnify Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value:,}" if isinstance(value, (int, float)) else f"  {key}: {value}")

    # Save
    merged.to_csv(args.output, index=False)
    print(f"\nSaved {len(merged):,} unified clonotypes to {args.output}")


# =============================================================================
# Mnemonic Command
# =============================================================================


def cmd_mnemonic(args):
    """Generate memorable names for TCR sequences."""
    import pandas as pd

    from .mnemonic import tcr_name

    setup_logging(args.verbose)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} entries from {args.input}")

    # Auto-detect CDR3 column
    cdr3_col = args.cdr3_col
    if not cdr3_col:
        for candidate in ["CDR3_beta", "CDR3_alpha", "cdr3", "cdr3_beta", "cdr3_alpha"]:
            if candidate in df.columns:
                cdr3_col = candidate
                break

    if not cdr3_col or cdr3_col not in df.columns:
        print(f"ERROR: Could not find CDR3 column. Available columns: {list(df.columns)}")
        return

    # Generate names
    df[args.name_col] = df[cdr3_col].apply(lambda x: tcr_name(x) if pd.notna(x) else None)

    # Save
    df.to_csv(args.output, index=False)
    print(f"Saved with mnemonic names to {args.output}")


# =============================================================================
# Main Parser
# =============================================================================


def create_parser():
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        prog="tcrsift",
        description="TCRsift: TCR selection from antigen-specific culture and scRNA/VDJ sequencing data",
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # -------------------------------------------------------------------------
    # Load command
    # -------------------------------------------------------------------------
    p_load = subparsers.add_parser("load", help="Load CellRanger outputs from sample sheet")
    p_load.add_argument("--sample-sheet", "-s", required=True, help="Sample sheet (CSV or YAML)")
    p_load.add_argument("--output", "-o", required=True, help="Output h5ad file")
    p_load.add_argument(
        "--min-genes", type=int, default=250, help="Min genes per cell (default: 250)"
    )
    p_load.add_argument(
        "--max-genes", type=int, default=15000, help="Max genes per cell (default: 15000)"
    )
    p_load.add_argument("--min-counts", type=int, default=500, help="Min UMI counts (default: 500)")
    p_load.add_argument(
        "--max-counts", type=int, default=100000, help="Max UMI counts (default: 100000)"
    )
    p_load.add_argument("--min-mito", type=float, default=2.0, help="Min mito %% (default: 2)")
    p_load.add_argument("--max-mito", type=float, default=8.0, help="Max mito %% (default: 8)")
    p_load.add_argument("--plot-qc", action="store_true", help="Generate QC plots")
    p_load.add_argument("--output-dir", help="Output directory for plots")
    p_load.add_argument("--verbose", action="store_true", help="Verbose output")
    p_load.set_defaults(func=cmd_load)

    # -------------------------------------------------------------------------
    # Phenotype command
    # -------------------------------------------------------------------------
    p_pheno = subparsers.add_parser("phenotype", help="Classify T cells as CD4/CD8")
    p_pheno.add_argument("--input", "-i", required=True, help="Input h5ad from load")
    p_pheno.add_argument("--output", "-o", required=True, help="Output h5ad")
    p_pheno.add_argument(
        "--cd4-cd8-ratio", type=float, default=3.0, help="Ratio for confident calls (default: 3.0)"
    )
    p_pheno.add_argument(
        "--min-cd3-reads", type=int, default=10, help="Min CD3 reads (default: 10)"
    )
    p_pheno.add_argument("--plot-phenotype", action="store_true", help="Generate phenotype plots")
    p_pheno.add_argument("--output-dir", help="Output directory for plots")
    p_pheno.add_argument("--verbose", action="store_true", help="Verbose output")
    p_pheno.set_defaults(func=cmd_phenotype)

    # -------------------------------------------------------------------------
    # Clonotype command
    # -------------------------------------------------------------------------
    p_clone = subparsers.add_parser("clonotype", help="Aggregate cells into clonotypes")
    p_clone.add_argument("--input", "-i", required=True, help="Input h5ad with phenotypes")
    p_clone.add_argument("--output", "-o", required=True, help="Output CSV")
    p_clone.add_argument(
        "--group-by", choices=["CDR3ab", "CDR3b_only"], default="CDR3ab", help="Grouping strategy"
    )
    p_clone.add_argument(
        "--handle-doublets", choices=["flag", "remove", "keep-primary"], default="flag"
    )
    p_clone.add_argument("--min-umi", type=int, default=2, help="Min UMIs per chain (default: 2)")
    p_clone.add_argument("--airr", help="Also output AIRR format to this path")
    p_clone.add_argument("--plot-clonotypes", action="store_true", help="Generate clonotype plots")
    p_clone.add_argument("--output-dir", help="Output directory for plots")
    p_clone.add_argument("--verbose", action="store_true", help="Verbose output")
    p_clone.set_defaults(func=cmd_clonotype)

    # -------------------------------------------------------------------------
    # Filter command
    # -------------------------------------------------------------------------
    p_filter = subparsers.add_parser("filter", help="Filter clonotypes with tiered confidence")
    p_filter.add_argument("--input", "-i", required=True, help="Input clonotypes CSV")
    p_filter.add_argument("--output", "-o", required=True, help="Output directory for tier CSVs")
    p_filter.add_argument(
        "--tcell-type", choices=["cd8", "cd4", "both"], default="cd8", help="T cell type filter"
    )
    p_filter.add_argument(
        "--method", choices=["threshold", "logistic"], default="threshold", help="Filtering method"
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
        "--filter-mode",
        choices=["fdr", "shared-high-freq", "cross-donor-public"],
        default="fdr",
        help="Named filter preset (default: fdr). 'shared-high-freq' applies "
        "min-methods-per-donor=2 and min-frequency-per-method=0.01. "
        "'cross-donor-public' adds min-donors=2.",
    )
    p_filter.add_argument(
        "--min-donors", type=int, default=0,
        help="Min distinct donors clone must appear in (#8/#15)",
    )
    p_filter.add_argument(
        "--min-methods-per-donor", type=int, default=0,
        help="Min distinct enrichment methods within at least one donor (#8/#15)",
    )
    p_filter.add_argument(
        "--min-cells-per-method", type=int, default=0,
        help="Min cells in at least one enrichment method (#15)",
    )
    p_filter.add_argument(
        "--min-frequency-per-method", type=float, default=0.0,
        help="Min frequency within a single enrichment method (#15)",
    )
    p_filter.add_argument("--plot-filter", action="store_true", help="Generate filter plots")
    p_filter.add_argument("--output-dir", help="Output directory for plots")
    p_filter.add_argument("--verbose", action="store_true", help="Verbose output")
    p_filter.set_defaults(func=cmd_filter)

    # -------------------------------------------------------------------------
    # Annotate command
    # -------------------------------------------------------------------------
    p_annot = subparsers.add_parser("annotate", help="Annotate with public databases")
    p_annot.add_argument("--input", "-i", required=True, help="Input filtered CSV")
    p_annot.add_argument("--output", "-o", required=True, help="Output annotated CSV")
    p_annot.add_argument("--vdjdb", help="Path to VDJdb")
    p_annot.add_argument("--iedb", help="Path to IEDB")
    p_annot.add_argument("--cedar", help="Path to CEDAR")
    p_annot.add_argument(
        "--match-by",
        choices=["CDR3ab", "CDR3b_only"],
        default="CDR3ab",
        help="Matching strategy (default: CDR3ab)",
    )
    p_annot.add_argument("--exclude-viral", action="store_true", help="Remove viral clones")
    p_annot.add_argument("--flag-only", action="store_true", help="Just flag, don't remove")
    p_annot.add_argument(
        "--plot-annotations", action="store_true", help="Generate annotation plots"
    )
    p_annot.add_argument("--output-dir", help="Output directory for plots")
    p_annot.add_argument("--verbose", action="store_true", help="Verbose output")
    p_annot.set_defaults(func=cmd_annotate)

    # -------------------------------------------------------------------------
    # Match-TIL command
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # TIL-Clonotype command
    # -------------------------------------------------------------------------
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
        "--min-cells", type=int, default=1, help="Minimum total TIL cells per clonotype (default: 1)"
    )
    p_til_clono.add_argument("--verbose", action="store_true", help="Verbose output")
    p_til_clono.set_defaults(func=cmd_til_clonotype)

    # -------------------------------------------------------------------------
    # TIL-Select command
    # -------------------------------------------------------------------------
    p_til_select = subparsers.add_parser(
        "til-select",
        help="Select promising TIL clonotypes from multi-timepoint VDJ+GEX inputs",
        description="""
Prioritize promising TIL clonotypes from one or more tumor timepoints.

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
        help="Compatibility option for trend reporting (default: 1.5).",
    )
    p_til_select.add_argument(
        "--trend-decrease-ratio-max",
        type=float,
        default=0.5,
        help="Compatibility option for trend reporting (default: 0.5).",
    )
    p_til_select.add_argument(
        "--marker-genes",
        type=str,
        default="CD4,CD8A,CD8B,GZMB,PRF1,IFNG,MKI67,TNFRSF9,CXCL13,ENTPD1",
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
        "--pyensembl-release",
        type=int,
        default=110,
        help="Compatibility option used in selected-clone PDF report metadata.",
    )
    p_til_select.add_argument(
        "--rank-by",
        type=str,
        default="mean_frequency",
        choices=["mean_frequency", "max_frequency", "total_cells", "marker_score_cp10k_mean", "marker_score_z_mean"],
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
    p_til_select.add_argument("--vdjdb", type=Path, default=None, help="Path to VDJdb TSV")
    p_til_select.add_argument("--iedb", type=Path, default=None, help="Path to IEDB TSV")
    p_til_select.add_argument("--cedar", type=Path, default=None, help="Path to CEDAR TSV")
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
    p_til_select.set_defaults(func=cmd_til_select)

    # -------------------------------------------------------------------------
    # Assemble command
    # -------------------------------------------------------------------------
    p_asm = subparsers.add_parser(
        "assemble",
        help="Assemble full-length sequences",
        description="""
Assemble full-length TCR sequences from clonotypes.

REQUIRED INPUTS:
  --input/-i     Clonotypes CSV (must have CDR3_alpha and CDR3_beta columns)
  --output/-o    Output CSV path

CONDITIONALLY REQUIRED:
  --contigs-dir  Required ONLY when using --leaders-from-contigs or
                 --alpha-leader=from_contig or --beta-leader=from_contig
                 (for extracting native leader peptides from CellRanger FASTAs)
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required inputs group
    asm_required = p_asm.add_argument_group("required arguments")
    asm_required.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input clonotypes CSV (must have CDR3_alpha and CDR3_beta)",
    )
    asm_required.add_argument(
        "--output", "-o", required=True, help="Output CSV with assembled sequences"
    )

    # Leader options group
    asm_leaders = p_asm.add_argument_group(
        "leader peptide options",
        "Leader peptides are added to the N-terminus. Use default leaders (CD8A/CD28) "
        "or extract native leaders from CellRanger contigs.",
    )
    asm_leaders.add_argument(
        "--alpha-leader",
        choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
        default="CD28",
        help="Alpha chain leader (default: CD28)",
    )
    asm_leaders.add_argument(
        "--beta-leader",
        choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
        default="CD8A",
        help="Beta chain leader (default: CD8A)",
    )
    asm_leaders.add_argument(
        "--no-leaders", action="store_true", help="Disable leaders on both chains"
    )
    asm_leaders.add_argument(
        "--leaders-from-contigs",
        action="store_true",
        help="Extract native leaders from CellRanger FASTAs (REQUIRES --contigs-dir)",
    )
    asm_leaders.add_argument(
        "--contigs-dir",
        metavar="DIR",
        help="Directory with CellRanger contig FASTAs. "
        "REQUIRED when using --leaders-from-contigs or "
        "--alpha/beta-leader=from_contig",
    )

    # Sequence options group
    asm_seq = p_asm.add_argument_group("sequence options")
    asm_seq.add_argument(
        "--include-constant",
        action="store_true",
        help="Include constant region in assembled sequences",
    )
    asm_seq.add_argument(
        "--constant-source",
        choices=["ensembl", "from-data"],
        default="ensembl",
        help="Source for constant regions (default: ensembl)",
    )
    asm_seq.add_argument(
        "--single-chain",
        action="store_true",
        help="Generate single-chain constructs (beta-linker-alpha)",
    )
    asm_seq.add_argument(
        "--linker", default="T2A", help="Linker peptide for single-chain (default: T2A)"
    )

    # Output options group
    asm_out = p_asm.add_argument_group("additional output options")
    asm_out.add_argument("--airr", metavar="PATH", help="Also output in AIRR format to this path")
    asm_out.add_argument("--fasta", metavar="PATH", help="Also output in FASTA format to this path")
    asm_out.add_argument("--plot-assembly", action="store_true", help="Generate assembly QC plots")
    asm_out.add_argument(
        "--output-dir", metavar="DIR", help="Output directory for plots (default: next to output)"
    )
    asm_out.add_argument("--verbose", action="store_true", help="Verbose output")
    p_asm.set_defaults(func=cmd_assemble)

    # -------------------------------------------------------------------------
    # Run command (unified pipeline)
    # -------------------------------------------------------------------------
    p_run = subparsers.add_parser("run", help="Run complete pipeline")
    p_run.add_argument("--sample-sheet", "-s", required=True, help="Sample sheet (CSV or YAML)")
    p_run.add_argument("--output-dir", "-o", required=True, help="Output directory")
    add_config_args(p_run)  # Add --config option

    # Load step parameters
    load_group = p_run.add_argument_group("Load options")
    load_group.add_argument("--min-genes", type=int, help="Min genes per cell (default: 250)")
    load_group.add_argument("--max-genes", type=int, help="Max genes per cell (default: 15000)")
    load_group.add_argument("--min-counts", type=int, help="Min UMI counts (default: 500)")
    load_group.add_argument("--max-counts", type=int, help="Max UMI counts (default: 100000)")
    load_group.add_argument(
        "--min-mito", type=float, dest="min_mito_pct", help="Min mito %% (default: 2)"
    )
    load_group.add_argument(
        "--max-mito", type=float, dest="max_mito_pct", help="Max mito %% (default: 8)"
    )

    # Phenotype step parameters
    pheno_group = p_run.add_argument_group("Phenotype options")
    pheno_group.add_argument(
        "--cd4-cd8-ratio", type=float, help="Ratio for confident calls (default: 3.0)"
    )
    pheno_group.add_argument("--min-cd3-reads", type=int, help="Min CD3 reads (default: 10)")

    # Clonotype step parameters
    clone_group = p_run.add_argument_group("Clonotype options")
    clone_group.add_argument(
        "--group-by", choices=["CDR3ab", "CDR3b_only"], help="Grouping strategy (default: CDR3ab)"
    )
    clone_group.add_argument(
        "--handle-doublets",
        choices=["flag", "remove", "keep-primary"],
        help="Doublet handling (default: flag)",
    )
    clone_group.add_argument("--min-umi", type=int, help="Min UMIs per chain (default: 2)")

    # Filter step parameters
    filter_group = p_run.add_argument_group("Filter options")
    filter_group.add_argument(
        "--tcell-type", choices=["cd8", "cd4", "both"], help="T cell type filter (default: cd8)"
    )
    filter_group.add_argument(
        "--method", choices=["threshold", "logistic"], help="Filtering method (default: threshold)"
    )
    filter_group.add_argument("--min-cells", type=int, help="Min cells per clone (default: 2)")
    filter_group.add_argument("--min-frequency", type=float, help="Min frequency (default: 0.0)")
    filter_group.add_argument(
        "--require-complete", action="store_true", default=None, help="Require complete TCR"
    )
    filter_group.add_argument(
        "--no-require-complete", dest="require_complete", action="store_false"
    )
    filter_group.add_argument(
        "--fdr-tiers", help="FDR tiers comma-separated (default: 0.15,0.1,0.01,0.001,0.0001)"
    )
    filter_group.add_argument(
        "--filter-mode",
        choices=["fdr", "shared-high-freq", "cross-donor-public"],
        help="Named filter preset (default: fdr). 'shared-high-freq' applies "
        "min-methods-per-donor=2 and min-frequency-per-method=0.01.",
    )
    filter_group.add_argument(
        "--min-donors", type=int, help="Min distinct donors per clone",
    )
    filter_group.add_argument(
        "--min-methods-per-donor", type=int,
        help="Min distinct enrichment methods within at least one donor",
    )
    filter_group.add_argument(
        "--min-cells-per-method", type=int,
        help="Min cells in at least one enrichment method",
    )
    filter_group.add_argument(
        "--min-frequency-per-method", type=float,
        help="Min frequency within a single enrichment method",
    )

    # Annotate step parameters
    annot_group = p_run.add_argument_group("Annotation options")
    annot_group.add_argument("--vdjdb", dest="vdjdb_path", help="Path to VDJdb")
    annot_group.add_argument("--iedb", dest="iedb_path", help="Path to IEDB")
    annot_group.add_argument("--cedar", dest="cedar_path", help="Path to CEDAR")
    annot_group.add_argument(
        "--match-by", choices=["CDR3ab", "CDR3b_only"], help="Matching strategy (default: CDR3ab)"
    )
    annot_group.add_argument(
        "--exclude-viral", action="store_true", default=None, help="Remove viral clones"
    )
    annot_group.add_argument(
        "--flag-only", action="store_true", default=None, help="Flag but don't remove viral"
    )

    # TIL step parameters
    til_group = p_run.add_argument_group("TIL matching options")
    til_group.add_argument("--til-samples", help="Comma-separated TIL sample names")
    til_group.add_argument(
        "--til-sample",
        action="append",
        metavar="NAME=TYPE:PATH",
        help="Repeatable direct TIL sample spec. TYPE is csv, h5ad, or vdj (alias for vdj_dir).",
    )
    til_group.add_argument(
        "--til-match-by", choices=["CDR3ab", "CDR3b_only"], help="TIL matching strategy"
    )
    til_group.add_argument("--min-til-cells", type=int, help="Min TIL cells to count (default: 1)")

    # Assemble step parameters
    asm_group = p_run.add_argument_group(
        "Assembly options",
        "Options for assembling full-length TCR sequences. "
        "NOTE: --contigs-dir is REQUIRED when using --leaders-from-contigs or leader=from_contig.",
    )
    asm_group.add_argument(
        "--alpha-leader",
        choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
        help="Alpha chain leader (default: CD28). Use 'none' for no leader.",
    )
    asm_group.add_argument(
        "--beta-leader",
        choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
        help="Beta chain leader (default: CD8A). Use 'none' for no leader.",
    )
    asm_group.add_argument(
        "--no-leaders", action="store_true", help="Disable leaders on both chains"
    )
    asm_group.add_argument(
        "--leaders-from-contigs",
        action="store_true",
        help="Extract native leaders from CellRanger FASTAs (REQUIRES --contigs-dir)",
    )
    asm_group.add_argument(
        "--include-constant",
        action="store_true",
        default=None,
        help="Include constant region (default: True)",
    )
    asm_group.add_argument("--no-include-constant", dest="include_constant", action="store_false")
    asm_group.add_argument(
        "--constant-source",
        choices=["ensembl", "from-data"],
        help="Constant region source (default: ensembl)",
    )
    asm_group.add_argument(
        "--linker", choices=["T2A", "P2A", "E2A", "F2A"], help="Linker peptide (default: T2A)"
    )
    asm_group.add_argument(
        "--contigs-dir",
        metavar="DIR",
        help="Directory with CellRanger contig FASTAs. "
        "REQUIRED when using --leaders-from-contigs or leader=from_contig",
    )
    asm_group.add_argument(
        "--single-chain",
        action="store_true",
        default=None,
        help="Generate single-chain constructs (default: True)",
    )
    asm_group.add_argument("--no-single-chain", dest="single_chain", action="store_false")

    # Output options
    out_group = p_run.add_argument_group("Output options")
    out_group.add_argument(
        "--skip-plots",
        dest="generate_plots",
        action="store_false",
        default=None,
        help="Skip plot generation",
    )
    out_group.add_argument(
        "--no-report",
        dest="generate_report",
        action="store_false",
        default=None,
        help="Skip report generation",
    )
    out_group.add_argument("--verbose", action="store_true", help="Verbose output")

    p_run.set_defaults(func=cmd_run)

    # -------------------------------------------------------------------------
    # Load-SCT command
    # -------------------------------------------------------------------------
    p_sct = subparsers.add_parser("load-sct", help="Load SCT platform data")
    p_sct.add_argument("--input", "-i", required=True, help="Input SCT Excel file")
    p_sct.add_argument("--output", "-o", required=True, help="Output CSV")
    p_sct.add_argument("--sheet-name", default="Cell", help="Excel sheet name (default: Cell)")
    p_sct.add_argument(
        "--min-snr", type=float, default=2.0, help="Min signal-to-noise ratio (default: 2.0)"
    )
    p_sct.add_argument(
        "--min-reads", type=int, default=10, help="Min reads per chain (default: 10)"
    )
    p_sct.add_argument(
        "--require-mutation-match",
        action="store_true",
        default=True,
        help="Require PE/APC mutation match",
    )
    p_sct.add_argument(
        "--no-require-mutation-match", dest="require_mutation_match", action="store_false"
    )
    p_sct.add_argument(
        "--require-compact-match", action="store_true", help="Require comPACT ID match"
    )
    p_sct.add_argument("--aggregate", action="store_true", help="Aggregate to unique clonotypes")
    p_sct.add_argument("--verbose", action="store_true", help="Verbose output")
    p_sct.set_defaults(func=cmd_load_sct)

    # -------------------------------------------------------------------------
    # Annotate-GEX command
    # -------------------------------------------------------------------------
    p_gex = subparsers.add_parser(
        "annotate-gex",
        help="Annotate TCR data with gene expression",
        description="""
Annotate TCR data with gene expression values from 10x Genomics data.

REQUIRED INPUTS:
  --input/-i     CSV file with cell barcodes (or clonotypes with barcodes)
  --output/-o    Output CSV path
  --gex-file     10x filtered_feature_bc_matrix.h5 file (always required)
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required inputs group
    gex_required = p_gex.add_argument_group("required arguments")
    gex_required.add_argument(
        "--input", "-i", required=True, help="Input CSV (cells or clonotypes with barcode column)"
    )
    gex_required.add_argument(
        "--output", "-o", required=True, help="Output CSV with added GEX columns"
    )
    gex_required.add_argument(
        "--gex-file", required=True, help="10x filtered_feature_bc_matrix.h5 file (REQUIRED)"
    )

    # GEX options group
    gex_opts = p_gex.add_argument_group("gene expression options")
    gex_opts.add_argument(
        "--barcode-col",
        default="barcode",
        help="Column containing cell barcodes (default: barcode)",
    )
    gex_opts.add_argument(
        "--genes",
        metavar="GENE1,GENE2,...",
        help="Comma-separated genes to extract (default: T cell markers)",
    )
    gex_opts.add_argument("--prefix", default="gex", help="Prefix for GEX columns (default: gex)")
    gex_opts.add_argument(
        "--no-qc", action="store_true", help="Skip QC metrics (n_reads, n_genes, pct_mito)"
    )

    # Aggregation options group
    gex_agg = p_gex.add_argument_group("aggregation options")
    gex_agg.add_argument(
        "--aggregate", action="store_true", help="Aggregate expression by clonotype (sum, mean)"
    )
    gex_agg.add_argument(
        "--group-col",
        default="CDR3_pair",
        help="Column to group by when aggregating (default: CDR3_pair)",
    )
    gex_agg.add_argument(
        "--cd4-cd8-counts",
        action="store_true",
        help="Compute CD4-only and CD8-only cell counts per clonotype",
    )
    gex_agg.add_argument("--verbose", action="store_true", help="Verbose output")
    p_gex.set_defaults(func=cmd_annotate_gex)

    # -------------------------------------------------------------------------
    # Unify command
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Mnemonic command
    # -------------------------------------------------------------------------
    p_mnem = subparsers.add_parser("mnemonic", help="Generate memorable TCR names")
    p_mnem.add_argument("--input", "-i", required=True, help="Input CSV")
    p_mnem.add_argument("--output", "-o", required=True, help="Output CSV")
    p_mnem.add_argument(
        "--cdr3-col", help="Column with CDR3 sequences (auto-detected if not specified)"
    )
    p_mnem.add_argument("--name-col", default="tcr_name", help="Output column name")
    p_mnem.add_argument("--verbose", action="store_true")
    p_mnem.set_defaults(func=cmd_mnemonic)

    # -------------------------------------------------------------------------
    # Generate config command
    # -------------------------------------------------------------------------
    p_config = subparsers.add_parser("generate-config", help="Generate example config file")
    p_config.add_argument("--output", "-o", default="tcrsift_config.yaml", help="Output YAML file")
    p_config.set_defaults(func=cmd_generate_config)

    return parser


def cmd_generate_config(args):
    """Generate an example configuration file."""
    from .config import generate_example_config

    generate_example_config(args.output)
    print(f"Generated example config: {args.output}")
    print("\nYou can customize this file and use it with:")
    print(f"  tcrsift run --config {args.output} --sample-sheet samples.csv -o output/")


def main(args=None):
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args(args)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        args.func(args)
    except Exception as e:
        logging.error(f"Error: {e}")
        if hasattr(args, "verbose") and args.verbose:
            raise
        sys.exit(1)


if __name__ == "__main__":
    main()
