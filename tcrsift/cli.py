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
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
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
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
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
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        plot_clonotypes(clonotypes, output_dir)
        print(f"Clonotype plots saved to {output_dir}")


# =============================================================================
# Filter Command
# =============================================================================

def cmd_filter(args):
    """Filter clonotypes with tiered confidence levels."""
    import pandas as pd

    from .filter import filter_clonotypes, get_filter_summary, split_by_tier
    from .plots import plot_filter

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    # Parse FDR tiers
    fdr_tiers = None
    if args.fdr_tiers:
        fdr_tiers = [float(x) for x in args.fdr_tiers.split(",")]

    filtered = filter_clonotypes(
        clonotypes,
        method=args.method,
        tcell_type=args.tcell_type,
        min_cells=args.min_cells,
        min_frequency=args.min_frequency,
        require_complete=args.require_complete,
        exclude_viral=args.exclude_viral,
        fdr_tiers=fdr_tiers,
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
        match_by=args.match_by if hasattr(args, 'match_by') else "CDR3ab",
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
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        plot_annotations(annotated, output_dir)
        print(f"Annotation plots saved to {output_dir}")


# =============================================================================
# Match-TIL Command
# =============================================================================

def cmd_match_til(args):
    """Match culture clonotypes against TIL data."""
    import anndata as ad
    import pandas as pd

    from .plots import plot_til
    from .til import get_til_summary, match_til

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} culture clonotypes from {args.input}")

    til_data = ad.read_h5ad(args.til_data)
    print(f"Loaded {len(til_data)} TIL cells from {args.til_data}")

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
        print(f"  {key}: {value}")

    # Save
    matched.to_csv(args.output, index=False)
    print(f"\nSaved TIL-matched clonotypes to {args.output}")

    # Generate plots if requested
    if args.plot_til:
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        plot_til(matched, output_dir)
        print(f"TIL plots saved to {output_dir}")


# =============================================================================
# Assemble Command
# =============================================================================

def cmd_assemble(args):
    """Assemble full-length TCR sequences."""
    import pandas as pd

    from .assemble import assemble_full_sequences, export_fasta, validate_sequences
    from .plots import plot_assembly

    setup_logging(args.verbose)

    # Handle leader shortcuts
    alpha_leader = args.alpha_leader
    beta_leader = args.beta_leader

    if getattr(args, 'no_leaders', False):
        alpha_leader = None
        beta_leader = None
    elif getattr(args, 'leaders_from_contigs', False):
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
        output_dir = Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        plot_assembly(assembled, output_dir)
        print(f"Assembly plots saved to {output_dir}")


# =============================================================================
# Run Command (Unified Pipeline)
# =============================================================================

def cmd_run(args):
    """Run the complete TCRsift pipeline."""
    from pathlib import Path

    import anndata as ad

    from .annotate import annotate_clonotypes
    from .assemble import assemble_full_sequences
    from .clonotype import aggregate_clonotypes
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
    from .til import match_til

    # Load config with CLI overrides
    config = load_config_with_args(args)

    # Auto-detect TIL samples from sample sheet if not explicitly specified
    if not config.til.til_samples:
        sample_sheet = load_sample_sheet(args.sample_sheet)
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
    filtered = filter_clonotypes(
        clonotypes,
        method=config.filter.method,
        tcell_type=tcell_type,
        min_cells=config.filter.min_cells,
        min_frequency=config.filter.min_frequency,
        require_complete=config.filter.require_complete,
        fdr_tiers=config.filter.fdr_tiers,
    )

    # Save by tier
    tier_dfs = split_by_tier(filtered)
    tier_counts = {}
    for tier, tier_df in tier_dfs.items():
        tier_df.to_csv(data_dir / f"filtered_{tier}.csv", index=False)
        tier_counts[tier] = len(tier_df)
        print(f"  {tier}: {len(tier_df)} clonotypes")

    funnel_counts["Filtered"] = len(filtered)

    if config.output.generate_plots:
        plot_filter(filtered, plots_dir)

    # Step 5: Annotate (if databases provided)
    annotated = filtered
    has_annotation = config.annotate.vdjdb_path or config.annotate.iedb_path or config.annotate.cedar_path
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
    if config.til.til_samples:
        print("\n[6/7] Matching against TIL samples...")
        # Load TIL data from the loaded h5ad
        til_adata = ad.read_h5ad(data_dir / "phenotyped.h5ad")
        til_samples = config.til.til_samples
        til_adata = til_adata[til_adata.obs["sample"].isin(til_samples)]

        if len(til_adata) > 0:
            til_matched = match_til(
                annotated,
                til_adata,
                match_by=config.til.match_by,
                min_til_cells=config.til.min_til_cells,
            )
            til_matched.to_csv(data_dir / "til_matched.csv", index=False)
            n_til = til_matched["til_match"].sum() if "til_match" in til_matched.columns else 0
            print(f"  Found {n_til} clonotypes in TILs")
            if config.output.generate_plots:
                plot_til(til_matched, plots_dir)
        else:
            print("  No TIL samples found")
    else:
        print("\n[6/7] Skipping TIL matching (no TIL samples specified)")

    # Step 7: Assemble (if requested)
    assembled = None
    has_leaders = config.assemble.alpha_leader is not None or config.assemble.beta_leader is not None
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
    else:
        print("\n[7/7] Skipping assembly")

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
        generate_report(plots_dir, format="html")
        generate_report(plots_dir, format="pdf")

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
        print(f"Warning: Barcode column '{args.barcode_col}' not found - skipping per-cell augmentation")

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
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # -------------------------------------------------------------------------
    # Load command
    # -------------------------------------------------------------------------
    p_load = subparsers.add_parser("load", help="Load CellRanger outputs from sample sheet")
    p_load.add_argument("--sample-sheet", "-s", required=True, help="Sample sheet (CSV or YAML)")
    p_load.add_argument("--output", "-o", required=True, help="Output h5ad file")
    p_load.add_argument("--min-genes", type=int, default=250, help="Min genes per cell (default: 250)")
    p_load.add_argument("--max-genes", type=int, default=15000, help="Max genes per cell (default: 15000)")
    p_load.add_argument("--min-counts", type=int, default=500, help="Min UMI counts (default: 500)")
    p_load.add_argument("--max-counts", type=int, default=100000, help="Max UMI counts (default: 100000)")
    p_load.add_argument("--min-mito", type=float, default=2.0, help="Min mito %% (default: 2)")
    p_load.add_argument("--max-mito", type=float, default=8.0, help="Max mito %% (default: 8)")
    p_load.add_argument("--plot-qc", action="store_true", help="Generate QC plots")
    p_load.add_argument("--output-dir", help="Output directory for plots")
    p_load.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_load.set_defaults(func=cmd_load)

    # -------------------------------------------------------------------------
    # Phenotype command
    # -------------------------------------------------------------------------
    p_pheno = subparsers.add_parser("phenotype", help="Classify T cells as CD4/CD8")
    p_pheno.add_argument("--input", "-i", required=True, help="Input h5ad from load")
    p_pheno.add_argument("--output", "-o", required=True, help="Output h5ad")
    p_pheno.add_argument("--cd4-cd8-ratio", type=float, default=3.0, help="Ratio for confident calls (default: 3.0)")
    p_pheno.add_argument("--min-cd3-reads", type=int, default=10, help="Min CD3 reads (default: 10)")
    p_pheno.add_argument("--plot-phenotype", action="store_true", help="Generate phenotype plots")
    p_pheno.add_argument("--output-dir", help="Output directory for plots")
    p_pheno.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_pheno.set_defaults(func=cmd_phenotype)

    # -------------------------------------------------------------------------
    # Clonotype command
    # -------------------------------------------------------------------------
    p_clone = subparsers.add_parser("clonotype", help="Aggregate cells into clonotypes")
    p_clone.add_argument("--input", "-i", required=True, help="Input h5ad with phenotypes")
    p_clone.add_argument("--output", "-o", required=True, help="Output CSV")
    p_clone.add_argument("--group-by", choices=["CDR3ab", "CDR3b_only"], default="CDR3ab", help="Grouping strategy")
    p_clone.add_argument("--handle-doublets", choices=["flag", "remove", "keep-primary"], default="flag")
    p_clone.add_argument("--min-umi", type=int, default=2, help="Min UMIs per chain (default: 2)")
    p_clone.add_argument("--airr", help="Also output AIRR format to this path")
    p_clone.add_argument("--plot-clonotypes", action="store_true", help="Generate clonotype plots")
    p_clone.add_argument("--output-dir", help="Output directory for plots")
    p_clone.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_clone.set_defaults(func=cmd_clonotype)

    # -------------------------------------------------------------------------
    # Filter command
    # -------------------------------------------------------------------------
    p_filter = subparsers.add_parser("filter", help="Filter clonotypes with tiered confidence")
    p_filter.add_argument("--input", "-i", required=True, help="Input clonotypes CSV")
    p_filter.add_argument("--output", "-o", required=True, help="Output directory for tier CSVs")
    p_filter.add_argument("--tcell-type", choices=["cd8", "cd4", "both"], default="cd8", help="T cell type filter")
    p_filter.add_argument("--method", choices=["threshold", "logistic"], default="threshold", help="Filtering method")
    p_filter.add_argument("--min-cells", type=int, default=2, help="Min cells per clone")
    p_filter.add_argument("--min-frequency", type=float, default=0.0, help="Min frequency")
    p_filter.add_argument("--require-complete", action="store_true", default=True, help="Require complete TCR")
    p_filter.add_argument("--no-require-complete", dest="require_complete", action="store_false")
    p_filter.add_argument("--fdr-tiers", default="0.15,0.1,0.01,0.001,0.0001", help="FDR tiers (comma-separated)")
    p_filter.add_argument("--exclude-viral", action="store_true", help="Exclude viral clones")
    p_filter.add_argument("--plot-filter", action="store_true", help="Generate filter plots")
    p_filter.add_argument("--output-dir", help="Output directory for plots")
    p_filter.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
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
    p_annot.add_argument("--exclude-viral", action="store_true", help="Remove viral clones")
    p_annot.add_argument("--flag-only", action="store_true", help="Just flag, don't remove")
    p_annot.add_argument("--plot-annotations", action="store_true", help="Generate annotation plots")
    p_annot.add_argument("--output-dir", help="Output directory for plots")
    p_annot.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_annot.set_defaults(func=cmd_annotate)

    # -------------------------------------------------------------------------
    # Match-TIL command
    # -------------------------------------------------------------------------
    p_til = subparsers.add_parser("match-til", help="Match against TIL data")
    p_til.add_argument("--input", "-i", required=True, help="Input culture clonotypes CSV")
    p_til.add_argument("--output", "-o", required=True, help="Output CSV with TIL matches")
    p_til.add_argument("--til-data", required=True, help="TIL data h5ad")
    p_til.add_argument("--match-by", choices=["CDR3ab", "CDR3b_only"], default="CDR3ab")
    p_til.add_argument("--min-til-cells", type=int, default=1, help="Min TIL cells to count")
    p_til.add_argument("--plot-til", action="store_true", help="Generate TIL plots")
    p_til.add_argument("--output-dir", help="Output directory for plots")
    p_til.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_til.set_defaults(func=cmd_match_til)

    # -------------------------------------------------------------------------
    # Assemble command
    # -------------------------------------------------------------------------
    p_asm = subparsers.add_parser("assemble", help="Assemble full-length sequences")
    p_asm.add_argument("--input", "-i", required=True, help="Input clonotypes CSV")
    p_asm.add_argument("--output", "-o", required=True, help="Output CSV with sequences")
    p_asm.add_argument("--contigs-dir", help="Directory with CellRanger contig FASTAs")
    p_asm.add_argument("--alpha-leader", choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
                      default="CD28", help="Alpha chain leader (default: CD28)")
    p_asm.add_argument("--beta-leader", choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
                      default="CD8A", help="Beta chain leader (default: CD8A)")
    p_asm.add_argument("--no-leaders", action="store_true", help="Disable leaders on both chains")
    p_asm.add_argument("--leaders-from-contigs", action="store_true",
                      help="Extract native leaders from contig FASTAs (requires --contigs-dir)")
    p_asm.add_argument("--include-constant", action="store_true", help="Include constant region")
    p_asm.add_argument("--constant-source", choices=["ensembl", "from-data"], default="ensembl")
    p_asm.add_argument("--single-chain", action="store_true", help="Generate single-chain constructs")
    p_asm.add_argument("--linker", default="T2A", help="Linker for single-chain (default: T2A)")
    p_asm.add_argument("--airr", help="Also output AIRR format")
    p_asm.add_argument("--fasta", help="Also output FASTA format")
    p_asm.add_argument("--plot-assembly", action="store_true", help="Generate assembly plots")
    p_asm.add_argument("--output-dir", help="Output directory for plots")
    p_asm.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
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
    load_group.add_argument("--min-mito", type=float, dest="min_mito_pct", help="Min mito %% (default: 2)")
    load_group.add_argument("--max-mito", type=float, dest="max_mito_pct", help="Max mito %% (default: 8)")

    # Phenotype step parameters
    pheno_group = p_run.add_argument_group("Phenotype options")
    pheno_group.add_argument("--cd4-cd8-ratio", type=float, help="Ratio for confident calls (default: 3.0)")
    pheno_group.add_argument("--min-cd3-reads", type=int, help="Min CD3 reads (default: 10)")

    # Clonotype step parameters
    clone_group = p_run.add_argument_group("Clonotype options")
    clone_group.add_argument("--group-by", choices=["CDR3ab", "CDR3b_only"], help="Grouping strategy (default: CDR3ab)")
    clone_group.add_argument("--handle-doublets", choices=["flag", "remove", "keep-primary"], help="Doublet handling (default: flag)")
    clone_group.add_argument("--min-umi", type=int, help="Min UMIs per chain (default: 2)")

    # Filter step parameters
    filter_group = p_run.add_argument_group("Filter options")
    filter_group.add_argument("--tcell-type", choices=["cd8", "cd4", "both"], help="T cell type filter (default: cd8)")
    filter_group.add_argument("--method", choices=["threshold", "logistic"], help="Filtering method (default: threshold)")
    filter_group.add_argument("--min-cells", type=int, help="Min cells per clone (default: 2)")
    filter_group.add_argument("--min-frequency", type=float, help="Min frequency (default: 0.0)")
    filter_group.add_argument("--require-complete", action="store_true", default=None, help="Require complete TCR")
    filter_group.add_argument("--no-require-complete", dest="require_complete", action="store_false")
    filter_group.add_argument("--fdr-tiers", help="FDR tiers comma-separated (default: 0.15,0.1,0.01,0.001,0.0001)")

    # Annotate step parameters
    annot_group = p_run.add_argument_group("Annotation options")
    annot_group.add_argument("--vdjdb", dest="vdjdb_path", help="Path to VDJdb")
    annot_group.add_argument("--iedb", dest="iedb_path", help="Path to IEDB")
    annot_group.add_argument("--cedar", dest="cedar_path", help="Path to CEDAR")
    annot_group.add_argument("--match-by", choices=["CDR3ab", "CDR3b_only"], help="Matching strategy (default: CDR3ab)")
    annot_group.add_argument("--exclude-viral", action="store_true", default=None, help="Remove viral clones")
    annot_group.add_argument("--flag-only", action="store_true", default=None, help="Flag but don't remove viral")

    # TIL step parameters
    til_group = p_run.add_argument_group("TIL matching options")
    til_group.add_argument("--til-samples", help="Comma-separated TIL sample names")
    til_group.add_argument("--til-match-by", choices=["CDR3ab", "CDR3b_only"], help="TIL matching strategy")
    til_group.add_argument("--min-til-cells", type=int, help="Min TIL cells to count (default: 1)")

    # Assemble step parameters
    asm_group = p_run.add_argument_group("Assembly options")
    asm_group.add_argument("--alpha-leader", choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
                          help="Alpha chain leader (default: CD28). Use 'none' for no leader.")
    asm_group.add_argument("--beta-leader", choices=["CD8A", "CD28", "IgK", "TRAC", "TRBC", "from_contig", "none"],
                          help="Beta chain leader (default: CD8A). Use 'none' for no leader.")
    asm_group.add_argument("--no-leaders", action="store_true", help="Disable leaders on both chains")
    asm_group.add_argument("--leaders-from-contigs", action="store_true",
                          help="Extract native leaders from contig FASTAs (requires --contigs-dir)")
    asm_group.add_argument("--include-constant", action="store_true", default=None, help="Include constant region (default: True)")
    asm_group.add_argument("--no-include-constant", dest="include_constant", action="store_false")
    asm_group.add_argument("--constant-source", choices=["ensembl", "from-data"], help="Constant region source (default: ensembl)")
    asm_group.add_argument("--linker", choices=["T2A", "P2A", "E2A", "F2A"], help="Linker peptide (default: T2A)")
    asm_group.add_argument("--contigs-dir", help="Directory with CellRanger contig FASTAs (for native leaders)")
    asm_group.add_argument("--single-chain", action="store_true", default=None, help="Generate single-chain constructs (default: True)")
    asm_group.add_argument("--no-single-chain", dest="single_chain", action="store_false")

    # Output options
    out_group = p_run.add_argument_group("Output options")
    out_group.add_argument("--skip-plots", dest="generate_plots", action="store_false", default=None, help="Skip plot generation")
    out_group.add_argument("--no-report", dest="generate_report", action="store_false", default=None, help="Skip report generation")
    out_group.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    p_run.set_defaults(func=cmd_run)

    # -------------------------------------------------------------------------
    # Load-SCT command
    # -------------------------------------------------------------------------
    p_sct = subparsers.add_parser("load-sct", help="Load SCT platform data")
    p_sct.add_argument("--input", "-i", required=True, help="Input SCT Excel file")
    p_sct.add_argument("--output", "-o", required=True, help="Output CSV")
    p_sct.add_argument("--sheet-name", default="Cell", help="Excel sheet name (default: Cell)")
    p_sct.add_argument("--min-snr", type=float, default=2.0, help="Min signal-to-noise ratio (default: 2.0)")
    p_sct.add_argument("--min-reads", type=int, default=10, help="Min reads per chain (default: 10)")
    p_sct.add_argument("--require-mutation-match", action="store_true", default=True,
                          help="Require PE/APC mutation match")
    p_sct.add_argument("--no-require-mutation-match", dest="require_mutation_match", action="store_false")
    p_sct.add_argument("--require-compact-match", action="store_true", help="Require comPACT ID match")
    p_sct.add_argument("--aggregate", action="store_true", help="Aggregate to unique clonotypes")
    p_sct.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_sct.set_defaults(func=cmd_load_sct)

    # -------------------------------------------------------------------------
    # Annotate-GEX command
    # -------------------------------------------------------------------------
    p_gex = subparsers.add_parser("annotate-gex", help="Annotate TCR data with gene expression")
    p_gex.add_argument("--input", "-i", required=True, help="Input CSV (cells or clonotypes)")
    p_gex.add_argument("--output", "-o", required=True, help="Output CSV with GEX columns")
    p_gex.add_argument("--gex-file", required=True,
                       help="Path to 10x filtered_feature_bc_matrix.h5 file")
    p_gex.add_argument("--barcode-col", default="barcode",
                       help="Column containing cell barcodes (default: barcode)")
    p_gex.add_argument("--genes", help="Comma-separated list of genes to extract (default: T cell markers)")
    p_gex.add_argument("--prefix", default="gex",
                       help="Prefix for GEX columns (default: gex)")
    p_gex.add_argument("--no-qc", action="store_true",
                       help="Skip QC metrics (n_reads, n_genes, pct_mito)")
    p_gex.add_argument("--aggregate", action="store_true",
                       help="Aggregate expression by clonotype (sum, mean)")
    p_gex.add_argument("--group-col", default="CDR3_pair",
                       help="Column to group by when aggregating (default: CDR3_pair)")
    p_gex.add_argument("--cd4-cd8-counts", action="store_true",
                       help="Compute CD4-only and CD8-only cell counts per clonotype")
    p_gex.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_gex.set_defaults(func=cmd_annotate_gex)

    # -------------------------------------------------------------------------
    # Unify command
    # -------------------------------------------------------------------------
    p_unify = subparsers.add_parser("unify", help="Unify clonotypes from multiple experiments")
    p_unify.add_argument("--inputs", "-i", nargs="+", required=True, help="Input CSV files to merge")
    p_unify.add_argument("--output", "-o", required=True, help="Output unified CSV")
    p_unify.add_argument("--add-occurrence-flags", action="store_true", default=True,
                        help="Add 'occurs_in_*' columns")
    p_unify.add_argument("--no-occurrence-flags", dest="add_occurrence_flags", action="store_false")
    p_unify.add_argument("--add-combined-stats", action="store_true", default=True,
                        help="Add combined statistics")
    p_unify.add_argument("--no-combined-stats", dest="add_combined_stats", action="store_false")
    p_unify.add_argument("--add-phenotype-confidence", action="store_true", default=True,
                        help="Add phenotype confidence columns")
    p_unify.add_argument("--no-phenotype-confidence", dest="add_phenotype_confidence", action="store_false")
    p_unify.add_argument("--phenotype-ratio-threshold", type=float, default=10.0,
                        help="CD4/CD8 ratio for confident classification (default: 10.0)")
    p_unify.add_argument("--conditions", help="Comma-separated condition names for top-condition analysis")
    p_unify.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    p_unify.set_defaults(func=cmd_unify)

    # -------------------------------------------------------------------------
    # Mnemonic command
    # -------------------------------------------------------------------------
    p_mnem = subparsers.add_parser("mnemonic", help="Generate memorable TCR names")
    p_mnem.add_argument("--input", "-i", required=True, help="Input CSV")
    p_mnem.add_argument("--output", "-o", required=True, help="Output CSV")
    p_mnem.add_argument("--cdr3-col", help="Column with CDR3 sequences (auto-detected if not specified)")
    p_mnem.add_argument("--name-col", default="tcr_name", help="Output column name")
    p_mnem.add_argument("--verbose", "-v", action="store_true")
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
        if hasattr(args, 'verbose') and args.verbose:
            raise
        sys.exit(1)


if __name__ == "__main__":
    main()
