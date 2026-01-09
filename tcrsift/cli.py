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

import sys
import logging
from pathlib import Path
from typing import Optional

import argparse


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
    from .phenotype import phenotype_cells, validate_phenotype_vs_expected, get_phenotype_summary
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
    from .clonotype import aggregate_clonotypes, get_clonotype_summary, export_clonotypes_airr
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
    from .filter import filter_clonotypes, split_by_tier, get_filter_summary
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
    import pandas as pd
    import anndata as ad
    from .til import match_til, get_til_summary
    from .plots import plot_til

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
    from .assemble import assemble_full_sequences, validate_sequences, export_fasta
    from .plots import plot_assembly

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    assembled = assemble_full_sequences(
        clonotypes,
        contigs_dir=args.contigs_dir,
        include_leader=args.include_leader,
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
    import anndata as ad
    import pandas as pd
    from pathlib import Path

    from .loader import load_samples
    from .phenotype import phenotype_cells, filter_by_tcell_type
    from .clonotype import aggregate_clonotypes
    from .filter import filter_clonotypes, split_by_tier
    from .annotate import annotate_clonotypes
    from .til import match_til
    from .assemble import assemble_full_sequences
    from .plots import (
        plot_qc, plot_phenotype, plot_clonotypes,
        plot_filter, plot_annotations, plot_til, plot_assembly,
        generate_report,
    )

    setup_logging(args.verbose)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    print("=" * 60)
    print("TCRsift Pipeline")
    print("=" * 60)

    # Step 1: Load
    print("\n[1/7] Loading samples...")
    adata = load_samples(args.sample_sheet)
    adata.write_h5ad(data_dir / "loaded.h5ad")
    print(f"  Loaded {len(adata)} cells")
    if not args.skip_plots:
        plot_qc(adata, plots_dir)

    # Step 2: Phenotype
    print("\n[2/7] Phenotyping T cells...")
    adata = phenotype_cells(adata)
    adata.write_h5ad(data_dir / "phenotyped.h5ad")
    if not args.skip_plots:
        plot_phenotype(adata, plots_dir)

    # Filter by T cell type if specified
    if args.tcell_type != "both":
        adata = filter_by_tcell_type(adata, args.tcell_type)
        print(f"  Filtered to {args.tcell_type.upper()}+: {len(adata)} cells")

    # Step 3: Clonotype
    print("\n[3/7] Aggregating clonotypes...")
    clonotypes = aggregate_clonotypes(adata)
    clonotypes.to_csv(data_dir / "clonotypes.csv", index=False)
    print(f"  Found {len(clonotypes)} clonotypes")
    if not args.skip_plots:
        plot_clonotypes(clonotypes, plots_dir)

    # Step 4: Filter
    print("\n[4/7] Filtering clonotypes...")
    fdr_tiers = [float(x) for x in args.fdr_tiers.split(",")] if args.fdr_tiers else None
    filtered = filter_clonotypes(
        clonotypes,
        method=args.method,
        min_cells=args.min_cells,
        fdr_tiers=fdr_tiers,
    )

    # Save by tier
    tier_dfs = split_by_tier(filtered)
    for tier, tier_df in tier_dfs.items():
        tier_df.to_csv(data_dir / f"filtered_{tier}.csv", index=False)
        print(f"  {tier}: {len(tier_df)} clonotypes")

    if not args.skip_plots:
        plot_filter(filtered, plots_dir)

    # Step 5: Annotate (if databases provided)
    annotated = filtered
    if args.vdjdb or args.iedb or args.cedar:
        print("\n[5/7] Annotating with public databases...")
        annotated = annotate_clonotypes(
            filtered,
            vdjdb_path=args.vdjdb,
            iedb_path=args.iedb,
            cedar_path=args.cedar,
            exclude_viral=args.exclude_viral,
        )
        annotated.to_csv(data_dir / "annotated.csv", index=False)
        n_viral = annotated["is_viral"].sum() if "is_viral" in annotated.columns else 0
        print(f"  Flagged {n_viral} viral clonotypes")
        if not args.skip_plots:
            plot_annotations(annotated, plots_dir)
    else:
        print("\n[5/7] Skipping annotation (no databases provided)")

    # Step 6: TIL matching (if TIL samples specified)
    til_matched = annotated
    if args.til_samples:
        print("\n[6/7] Matching against TIL samples...")
        # Load TIL data from the loaded h5ad
        til_adata = ad.read_h5ad(data_dir / "phenotyped.h5ad")
        til_samples = args.til_samples.split(",")
        til_adata = til_adata[til_adata.obs["sample"].isin(til_samples)]

        if len(til_adata) > 0:
            til_matched = match_til(annotated, til_adata)
            til_matched.to_csv(data_dir / "til_matched.csv", index=False)
            n_til = til_matched["til_match"].sum() if "til_match" in til_matched.columns else 0
            print(f"  Found {n_til} clonotypes in TILs")
            if not args.skip_plots:
                plot_til(til_matched, plots_dir)
        else:
            print("  No TIL samples found")
    else:
        print("\n[6/7] Skipping TIL matching (no TIL samples specified)")

    # Step 7: Assemble (if requested)
    if args.assemble:
        print("\n[7/7] Assembling full-length sequences...")
        assembled = assemble_full_sequences(
            til_matched,
            include_leader=args.include_leader,
            include_constant=args.include_constant,
        )
        assembled.to_csv(data_dir / "full_sequences.csv", index=False)
        print(f"  Assembled {len(assembled)} sequences")
        if not args.skip_plots:
            plot_assembly(assembled, plots_dir)
    else:
        print("\n[7/7] Skipping assembly")

    # Generate report
    if args.report:
        print("\nGenerating report...")
        generate_report(plots_dir, format="html")
        generate_report(plots_dir, format="pdf")

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)


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
    p_asm.add_argument("--include-leader", action="store_true", help="Include leader peptide")
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
    p_run.add_argument("--tcell-type", choices=["cd8", "cd4", "both"], default="cd8")
    p_run.add_argument("--method", choices=["threshold", "logistic"], default="threshold")
    p_run.add_argument("--min-cells", type=int, default=2)
    p_run.add_argument("--fdr-tiers", default="0.15,0.1,0.01,0.001,0.0001")
    p_run.add_argument("--vdjdb", help="Path to VDJdb")
    p_run.add_argument("--iedb", help="Path to IEDB")
    p_run.add_argument("--cedar", help="Path to CEDAR")
    p_run.add_argument("--exclude-viral", action="store_true")
    p_run.add_argument("--til-samples", help="Comma-separated TIL sample names")
    p_run.add_argument("--assemble", action="store_true", help="Run assembly step")
    p_run.add_argument("--include-leader", action="store_true")
    p_run.add_argument("--include-constant", action="store_true")
    p_run.add_argument("--report", action="store_true", help="Generate combined report")
    p_run.add_argument("--skip-plots", action="store_true", help="Skip plot generation")
    p_run.add_argument("--verbose", "-v", action="store_true")
    p_run.set_defaults(func=cmd_run)

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

    return parser


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
