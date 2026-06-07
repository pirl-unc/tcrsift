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
    from .config import AttributionConfig
    from .plots import plot_clonotypes

    setup_logging(args.verbose)

    adata = ad.read_h5ad(args.input)
    print(f"Loaded {len(adata)} cells from {args.input}")

    attribution = AttributionConfig(
        enabled=getattr(args, "attribution", False),
        dual_alpha_min_cells=getattr(args, "dual_alpha_min_cells", 2),
        share_weighting=getattr(args, "share_weighting", "proportional"),
        **{
            k: v
            for k, v in (
                ("dual_alpha_split", getattr(args, "dual_alpha_split", None)),
                ("dual_beta_split", getattr(args, "dual_beta_split", None)),
            )
            if v is not None
        },
    )
    clonotypes = aggregate_clonotypes(
        adata,
        group_by=args.group_by,
        min_umi=args.min_umi,
        handle_doublets=args.handle_doublets,
        attribution=attribution,
        doublet_warn_rate=getattr(args, "doublet_warn_rate", 0.1),
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
        "min_timepoints": args.min_timepoints,
        "min_timepoints_per_donor": args.min_timepoints_per_donor,
        "min_apcs": args.min_apcs,
        "min_apcs_per_donor": args.min_apcs_per_donor,
        "min_til_cells_per_donor": args.min_til_cells_per_donor,
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
    # Per-sample tier exposure (#122): merge `tier_in_{sample}` columns onto the
    # filtered table from a clone_sample_long, so the selection language can
    # evaluate per-method include_tier / exclude_tier3plus without recomputing.
    if getattr(args, "emit_per_sample_tier", False):
        if not getattr(args, "clone_sample_long", None):
            print(
                "  --emit-per-sample-tier needs --clone-sample-long PATH "
                "(the clone_sample_long.csv from `tcrsift run`); skipping."
            )
        else:
            from .selection import pivot_per_sample_tiers

            long_df = pd.read_csv(args.clone_sample_long)
            wide = pivot_per_sample_tiers(long_df)
            filtered = filtered.merge(wide, on="CDR3ab", how="left")
            n_tier_cols = sum(c.startswith("tier_in_") for c in filtered.columns)
            print(f"  Added {n_tier_cols} per-sample tier column(s)")

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
    from .datacache import cached_path
    from .plots import plot_annotations

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    cache_dir = getattr(args, "cache_dir", None)
    annotated = annotate_clonotypes(
        clonotypes,
        vdjdb_path=args.vdjdb or cached_path("vdjdb", cache_dir),
        iedb_path=args.iedb or cached_path("iedb", cache_dir),
        cedar_path=args.cedar or cached_path("cedar", cache_dir),
        iedb_epitope_path=(
            getattr(args, "iedb_epitope", None)
            or cached_path("iedb_epitope", cache_dir)
        ),
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

    from .assemble import (
        assemble_full_sequences,
        assemble_qc_report,
        export_fasta,
        synthesis_qc_report,
        validate_sequences,
    )
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

    # Parse --stop-codons (#116). Empty string → no stops.
    stop_codons = tuple(
        c.strip().upper() for c in args.stop_codons.split(",") if c.strip()
    )

    assembled = assemble_full_sequences(
        clonotypes,
        contigs_dir=args.contigs_dir,
        cellranger_dir=getattr(args, "cellranger_dir", None),
        sample_name_from=getattr(args, "sample_name_from", None) or "parent",
        alpha_leader=alpha_leader,
        beta_leader=beta_leader,
        include_constant=args.include_constant,
        constant_source=args.constant_source,
        linker=args.linker if args.single_chain else None,
        trac_allele=args.trac_allele,
        trbc1_allele=args.trbc1_allele,
        trbc2_allele=args.trbc2_allele,
        stop_codons=stop_codons,
    )

    # Validate (per-row warnings + aggregate QC report).
    warnings = validate_sequences(assembled)
    for w in warnings:
        print(f"WARNING: {w}")
    print(assemble_qc_report(assembled))
    synth_text = synthesis_qc_report(assembled)
    if synth_text:
        print(synth_text)

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
# Data Command (Reference-database cache)
# =============================================================================


def _format_size(n: int) -> str:
    """Human-friendly byte count, e.g. '12 MB', '4.2 GB'."""
    if n == 0:
        return "-"
    units = ("B", "KB", "MB", "GB", "TB")
    size = float(n)
    for unit in units:
        if size < 1024 or unit == units[-1]:
            return f"{size:.0f} {unit}" if size >= 10 or unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{n} B"


def cmd_audit_alleles(args):
    """Run cohort-level allele audit + novel-allele detection (#119, #120).

    Reads a ``full_sequences.csv`` produced by ``tcrsift assemble`` or
    ``tcrsift run`` and prints the audit report. When ``--output`` is
    given, also writes the long-form novel-allele candidates CSV.
    """
    import pandas as pd

    from .assemble import allele_audit_report, detect_novel_alleles

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")
    print()

    report = allele_audit_report(
        df,
        min_pct=args.min_pct,
        min_v_spread=args.min_v_spread,
        min_samples=args.min_samples,
    )
    print(report)

    if args.output:
        candidates = detect_novel_alleles(
            df,
            min_pct=args.min_pct,
            min_v_spread=args.min_v_spread,
            min_samples=args.min_samples,
        )
        candidates.to_csv(args.output, index=False)
        print()
        print(f"Wrote novel-allele candidate table to {args.output}")


def cmd_ppost(args):
    """Compute OLGA Pgen + SONIA Ppost per clonotype (#143).

    Reads a clonotype CSV with CDR3 + V/J columns, adds ``log10_pgen_olga``,
    ``sonia_q`` and ``log10_ppost`` for the chosen chain(s), and writes the
    augmented table. OLGA/SONIA are a GPL-3.0 optional extra — install with
    ``pip install tcrsift[pgen]``.
    """
    import pandas as pd

    from .olga_ppost import (
        compute_pgen_ppost,
        flag_private_candidates,
        olga_sonia_available,
    )

    if not olga_sonia_available():
        from .olga_ppost import _INSTALL_HINT
        print(_INSTALL_HINT)
        return 1

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")

    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    for chain in chains:
        cdr3_col = args.cdr3_col or f"CDR3_{chain}"
        v_col = args.v_gene_col or f"{chain}_v_gene"
        j_col = args.j_gene_col or f"{chain}_j_gene"
        if cdr3_col not in df.columns:
            print(f"  [{chain}] skipped: no {cdr3_col!r} column")
            continue
        suffix = f"_{chain}"
        df = compute_pgen_ppost(
            df, chain=chain, cdr3_col=cdr3_col, v_gene_col=v_col,
            j_gene_col=j_col,
            pgen_col=f"log10_pgen_olga{suffix}",
            ppost_col=f"log10_ppost{suffix}",
            q_col=f"sonia_q{suffix}",
        )
        df[f"private_candidate{suffix}"] = flag_private_candidates(
            df, score_col=f"log10_ppost{suffix}",
            freq_col=args.freq_col if args.freq_col in df.columns else None,
        )
        n_pc = int(df[f"private_candidate{suffix}"].sum())
        print(f"  [{chain}] computed Pgen/Ppost; {n_pc} private candidates")

    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")


def cmd_log_pgen(args):
    """Score per-clone log Pgen with a data-driven background (seqprob).

    Reads a clonotype CSV and adds ``log_pgen_<chain>`` using the shipped
    OLGA-generated k-mer background (default) or a TCRpeg model
    (``--backend tcrpeg``, needs ``pip install tcrsift[tcrpeg]``). This is
    the data-driven replacement for the OLGA/SONIA runtime path: no GPL deps
    for the default k-mer backend.
    """
    import pandas as pd

    from .seqprob import score_log_pgen

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")

    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    for chain in chains:
        cdr3_col = args.cdr3_col or f"CDR3_{chain}"
        if cdr3_col not in df.columns:
            print(f"  [{chain}] skipped: no {cdr3_col!r} column")
            continue
        try:
            df[f"log_pgen_{chain}"] = score_log_pgen(
                df, chain=chain, cdr3_col=cdr3_col, backend=args.backend,
            )
        except (FileNotFoundError, ImportError) as exc:
            print(f"  [{chain}] {exc}")
            return 1
        n_scored = int(df[f"log_pgen_{chain}"].notna().sum())
        print(f"  [{chain}] scored {n_scored} clones (backend={args.backend})")

    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output}")


def cmd_pgen_fetch(args):
    """Fetch + cache a neutral observed repertoire for Ppost training."""
    from .pgen_models import (
        REPERTOIRE_SOURCES,
        fetch_healthy_pbmc,
        fetch_repertoire,
    )

    # The pooled 5×10x healthy-PBMC reference downloads both chains at once.
    if args.source == "healthy-pbmc" and not args.url:
        print("Fetching the 5 pooled 10x healthy-PBMC VDJ-T sets...")
        for chain, path in fetch_healthy_pbmc().items():
            print(f"  [{chain}] pooled healthy-PBMC reference → {path}")
        return

    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    if args.source and args.source in REPERTOIRE_SOURCES and not args.url:
        print(f"Source {args.source!r}: {REPERTOIRE_SOURCES[args.source]}")
        print("Download a healthy/neutral (not antigen-sorted) repertoire CSV "
              "from there and pass it via --url, or use --source healthy-pbmc.")
    for chain in chains:
        try:
            path = fetch_repertoire(
                chain, url=args.url, source=args.source,
                seq_col=args.seq_col, v_col=args.v_col,
                chain_prefix=("TRA" if chain == "alpha" else "TRB")
                if args.v_col else None,
            )
            print(f"  [{chain}] cached → {path}")
        except ValueError as exc:
            print(f"  [{chain}] {exc}")


def cmd_pgen_train(args):
    """Train + cache Pgen/Ppost background models (kmer or tcrpeg)."""
    from .pgen_models import train_model

    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    roles = ["pgen", "ppost"] if args.role == "both" else [args.role]
    for chain in chains:
        for role in roles:
            try:
                path = train_model(
                    chain, backend=args.backend, role=role, url=args.url,
                    n_olga=args.n_olga, epochs=args.epochs, device=args.device,
                )
                print(f"  trained {args.backend} {role} {chain} → {path}")
            except (ImportError, ValueError) as exc:
                print(f"  [{role} {chain}] skipped: {exc}")


def cmd_pgen_annotate(args):
    """Annotate a TCR CSV with pgen/ppost (α/β); auto-trains on first use."""
    import pandas as pd

    from .annotate_tcrs import add_pgen_ppost

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")
    chains = ("alpha", "beta") if args.chain == "both" else (args.chain,)
    df = add_pgen_ppost(df, chains=chains, backend=args.backend,
                        auto_train=not args.no_auto_train)
    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output} (backend={args.backend})")


def cmd_pgen_status(args):
    """Show cached Pgen/Ppost training data + trained models."""
    from .pgen_models import cache_dir, cached_model_file, cached_repertoire_file

    print(f"seqprob cache: {cache_dir()}")
    print("Training repertoires:")
    for chain in ("alpha", "beta"):
        p = cached_repertoire_file(chain)
        print(f"  {chain}: {'present' if p.is_file() else 'absent'}  ({p})")
    print("Models:")
    for backend in ("kmer", "tcrpeg"):
        for role in ("pgen", "ppost"):
            for chain in ("alpha", "beta"):
                p = cached_model_file(backend, role, chain)
                if p.is_file():
                    print(f"  {backend} {role} {chain}: {p}")


def cmd_data_list(args):
    """List managed reference databases and their cache state."""
    from .datacache import inspect_cache, resolve_cache_dir

    cache_dir = resolve_cache_dir(args.cache_dir)
    entries = inspect_cache(args.cache_dir)

    print(f"Cache directory: {cache_dir}")
    print()
    print(f"{'database':<10} {'status':<8} {'size':>10}  path")
    print("-" * 80)
    total = 0
    for e in entries:
        status = "present" if e.present else "missing"
        size = _format_size(e.size_bytes)
        print(f"{e.name:<10} {status:<8} {size:>10}  {e.path}")
        total += e.size_bytes
    print("-" * 80)
    print(f"{'total':<10} {'':<8} {_format_size(total):>10}")

    missing = [e for e in entries if not e.present]
    if missing:
        print()
        print("Missing databases — download manually and place under the cache:")
        for e in missing:
            print(f"  {e.name}: {e.source} → {e.path}")


def cmd_data_clear(args):
    """Remove cached files for one or more managed databases."""
    from .datacache import clear_cache

    # `args.db` is None (clear all) or a list[str] of names. Deduplicate
    # while preserving order so the user's first mention determines
    # report order.
    targets: list[str] | None = None
    if args.db:
        seen: set[str] = set()
        targets = [n for n in args.db if not (n in seen or seen.add(n))]

    removed_total: list = []
    if targets is None:
        removed_total = clear_cache(db=None, data_dir=args.cache_dir)
    else:
        for name in targets:
            removed_total.extend(clear_cache(db=name, data_dir=args.cache_dir))

    if not removed_total:
        print("Nothing to remove.")
        return
    for path in removed_total:
        print(f"Removed {path}")


def cmd_data_download(args):
    """Download one or more supported reference databases into the cache."""
    from .datacache import DATABASES, DownloadError, download_database

    setup_logging(verbose=True)
    # `args.db` is None (download all) or a list[str]. Deduplicate while
    # preserving order — argparse with action=extend doesn't dedupe, and
    # repeating a name shouldn't trigger a second download.
    if args.db:
        seen: set[str] = set()
        targets = [n for n in args.db if not (n in seen or seen.add(n))]
    else:
        targets = [
            name for name, spec in DATABASES.items() if spec.download_url is not None
        ]
    failures: list[tuple[str, str]] = []
    for db in targets:
        try:
            path = download_database(db, data_dir=args.cache_dir, force=args.force)
            print(f"  {db}: cached at {path}")
        except DownloadError as e:
            failures.append((db, str(e)))
            print(f"  {db}: FAILED — {e}")
    if failures:
        raise SystemExit(1)


# =============================================================================
# Run Command (Unified Pipeline)
# =============================================================================


def _common_contig_root(vdj_dirs):
    """Common-ancestor directory of the sample sheet's vdj_dirs, to scan for
    contigs in 'sheet' mode (#201). Returns ``None`` when there's no shared
    ancestor or it resolves to the filesystem root (too broad to rglob safely)."""
    import os

    paths = [os.path.abspath(str(d)) for d in vdj_dirs if d]
    if not paths:
        return None
    try:
        root = os.path.commonpath(paths)
    except ValueError:
        return None  # mixed drives / no common path
    p = Path(root)
    if p == Path(p.anchor):  # filesystem root — refuse
        return None
    return root


def cmd_run(args):
    """Run the complete TCRsift pipeline."""
    from pathlib import Path

    import pandas as pd

    from .annotate import annotate_clonotypes
    from .assemble import assemble_full_sequences, export_fasta
    from .clonotype import aggregate_clonotypes, export_clonotypes_airr
    from .filter import (
        filter_clonotypes,
        pick_strictest_populated_tier,
        resolve_fdr_tiers_for_method,
        split_by_tier,
    )
    from .loader import load_samples
    from .phenotype import apply_cd3_gate, filter_by_tcell_type, phenotype_cells
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
        set_plot_format,
    )
    from .sample_sheet import load_sample_sheet
    from .til import load_til_specs, match_til

    # Early validation of arguments
    validate_run_args(args)

    # Load config with CLI overrides
    config = load_config_with_args(args)

    # Figure output format (#169): png (default) or a vector format (pdf/svg)
    # emitted alongside a PNG. Applies to every save_figure in this run.
    set_plot_format(config.output.plot_format)

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

    # CD3 read gate on the cells entering clonotype aggregation. Previously
    # min_cd3_reads only set a column (filter:min_cd3) and never gated selection
    # (#172). Default 3 (a gentle contamination floor); 0 disables. A raw CD3 UMI
    # count tracks sequencing depth, so this is a low floor, not a large cull.
    min_cd3 = config.phenotype.min_cd3_reads
    n_before_cd3 = len(adata)
    adata = apply_cd3_gate(adata, min_cd3)
    if len(adata) != n_before_cd3:
        print(f"  CD3 gate (>={min_cd3} reads): {len(adata)}/{n_before_cd3} cells pass")

    funnel_counts["Phenotyped"] = len(adata)

    # Step 3: Clonotype
    print("\n[3/7] Aggregating clonotypes...")
    clonotypes = aggregate_clonotypes(
        adata,
        group_by=config.clonotype.group_by,
        min_umi=config.clonotype.min_umi,
        handle_doublets=config.clonotype.handle_doublets,
        attribution=config.attribution,
        doublet_warn_rate=config.clonotype.doublet_warn_rate,
    )
    if config.attribution.enabled:
        print("  Attribution ON: weighted fractional clone counts (#176)")
    clonotypes.to_csv(data_dir / "clonotypes.csv", index=False)
    funnel_counts["Clonotypes"] = len(clonotypes)
    print(f"  Found {len(clonotypes)} clonotypes")
    if config.output.generate_plots:
        plot_clonotypes(clonotypes, plots_dir)

    # Long-format per-(clone, sample) companion table (#20 chunk 1). Built
    # whenever it'll be consumed downstream — the chunk-1 CSV emission OR
    # the chunk-2 per-method rankings (which need long_df even when the
    # CSV is suppressed). Disk emission is gated by --emit-clone-sample-long.
    emit_long = config.output.emit_clone_sample_long
    n_samples_in_run = adata.obs["sample"].nunique() if "sample" in adata.obs.columns else 0
    should_emit_long = (
        emit_long == "always"
        or (emit_long == "auto" and n_samples_in_run >= 2)
    )
    have_method_axis = "enrichment_method" in adata.obs.columns

    long_df = None
    if should_emit_long or have_method_axis:
        from .clonotype import build_clone_sample_long
        long_df = build_clone_sample_long(adata, attribution=config.attribution)

    # Per-sample clonal-expansion QC (#161): a one-line expansion fingerprint
    # that flags an unexpanded/near-polyclonal sample (failed expansion, swap,
    # or baseline sort) at a glance instead of via a forensic dig.
    if long_df is not None and "sample" in long_df.columns:
        from .qc import clonal_expansion_metrics
        exp_qc = clonal_expansion_metrics(long_df)
        if not exp_qc.empty:
            exp_qc.to_csv(data_dir / "expansion_qc.csv", index=False)
            print(f"  Wrote expansion_qc.csv: {len(exp_qc)} samples "
                  f"(top1% / effective-clones / clonality per sample)")
            for _, r in exp_qc[exp_qc["unexpanded"]].iterrows():
                print(f"  ⚠️  sample {r['sample']!r}: {r['warning']}")

    # Per-sample integrity QC (#161): GEX↔VDJ overlap, inter-sample barcode
    # Jaccard (duplicate-input/swap), CDR3 anchor integrity, cellranger
    # Sample ID — surfaces a mislabel/swap/duplicate without a forensic dig.
    if "sample" in adata.obs.columns:
        from .qc import inter_sample_barcode_jaccard, sample_integrity_qc
        sample_qc = sample_integrity_qc(adata)
        if not sample_qc.empty:
            sample_qc.to_csv(data_dir / "sample_qc.csv", index=False)
            print(f"  Wrote sample_qc.csv: {len(sample_qc)} samples "
                  "(GEX↔VDJ overlap / barcode Jaccard / anchor integrity)")
            jac = inter_sample_barcode_jaccard(adata)
            if not jac.empty:
                jac.to_csv(data_dir / "sample_barcode_jaccard.csv")
            for _, r in sample_qc[sample_qc["warning"] != ""].iterrows():
                print(f"  ⚠️  sample {r['sample']!r}: {r['warning']}")

    # Per-condition signature-consistency QC (#161): does each sort behave
    # like its label (its GEX signature enriched), consistently across donors?
    if "enrichment_method" in adata.obs.columns:
        try:
            from .sort_qc import sort_signature_consistency_from_adata
            sort_qc = sort_signature_consistency_from_adata(adata)
        except Exception as exc:  # GEX / signature genes may be absent
            logging.debug("sort-consistency QC skipped: %s", exc)
            sort_qc = pd.DataFrame()
        if not sort_qc.empty:
            sort_qc.to_csv(data_dir / "sort_consistency_qc.csv", index=False)
            print(f"  Wrote sort_consistency_qc.csv: {len(sort_qc)} "
                  "(donor, sort) rows (signature enrichment vs label)")
            for _, r in sort_qc[sort_qc["warning"] != ""].iterrows():
                print(f"  ⚠️  {r['donor']}/{r['enrichment_method']}: {r['warning']}")

    if should_emit_long and long_df is not None:
        # Enrich the emitted table with a per-sample abundance tier (the
        # per-(clone, sample) analogue of the global clonotype tier) via
        # the shared selection helper, so downstream report/selection
        # consumers don't have to re-derive it. Backward-compatible: adds
        # one column to the long table.
        from .selection import attach_per_sample_tiers
        emitted_long = attach_per_sample_tiers(long_df)
        emitted_long.to_csv(data_dir / "clone_sample_long.csv", index=False)
        print(
            f"  Wrote clone_sample_long.csv: {len(emitted_long)} rows "
            f"({emitted_long['CDR3ab'].nunique()} clones x {n_samples_in_run} samples)"
        )

    # Step 4: Filter
    print("\n[4/7] Filtering clonotypes...")
    # Resolve named filter mode preset, then layer user-supplied knobs.
    from .filter import resolve_filter_mode_kwargs

    user_filter_kwargs = {
        "min_donors": config.filter.min_donors,
        "min_methods_per_donor": config.filter.min_methods_per_donor,
        "min_cells_per_method": config.filter.min_cells_per_method,
        "min_frequency_per_method": config.filter.min_frequency_per_method,
        "min_timepoints": config.filter.min_timepoints,
        "min_timepoints_per_donor": config.filter.min_timepoints_per_donor,
        "min_apcs": config.filter.min_apcs,
        "min_apcs_per_donor": config.filter.min_apcs_per_donor,
        "min_til_cells_per_donor": config.filter.min_til_cells_per_donor,
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

    # FDR scope resolution (#26). 'auto' resolves to 'per-donor' for
    # multi-donor cohorts that don't share antigen, else 'global'.
    from .filter import (
        resolve_fdr_scope,
        split_clonotypes_by_donor,
        split_clonotypes_by_sample,
    )

    n_donors_in_run = (
        int(clonotypes["n_donors"].max())
        if "n_donors" in clonotypes.columns and len(clonotypes)
        else 1
    )
    sheet_shares_antigen = getattr(sample_sheet, "donors_share_antigen", False)
    donors_share_antigen = (
        config.filter.donors_share_antigen or sheet_shares_antigen
    )
    resolved_scope = resolve_fdr_scope(
        config.filter.fdr_scope,
        n_donors=n_donors_in_run,
        donors_share_antigen=donors_share_antigen,
    )
    print(f"  FDR scope: {resolved_scope}")

    # fdr_tiers only drives logistic (FDR) tiering. Under method='threshold' it
    # is inert, so a config that sets fdr_tiers (or filter_mode='fdr') alongside
    # method='threshold' would silently get fixed abundance tiers (#170). Pass it
    # through only when it's actually used (logistic); when a *customized* list
    # can't take effect, warn loudly once (here — passing None below keeps
    # filter_clonotypes from emitting a duplicate warning).
    from .config import FilterConfig

    fdr_tiers_arg, fdr_inert = resolve_fdr_tiers_for_method(
        config.filter.method, config.filter.fdr_tiers, FilterConfig().fdr_tiers
    )
    if fdr_inert:
        print(
            f"  WARNING: fdr_tiers is set but filter.method={config.filter.method!r}; "
            "fdr_tiers only applies to method='logistic'. Using fixed abundance "
            "tiers — the fdr_tiers list is ignored. Set method: logistic for "
            "FDR-based tiering."
        )

    def _filter_one(df_in):
        return filter_clonotypes(
            df_in,
            method=config.filter.method,
            tcell_type=tcell_type,
            min_cells=config.filter.min_cells,
            min_frequency=config.filter.min_frequency,
            require_complete=config.filter.require_complete,
            fdr_tiers=fdr_tiers_arg,
            **mode_kwargs,
        )

    if resolved_scope == "global":
        filtered = _filter_one(clonotypes)
        # Save by tier (existing behavior)
        tier_dfs = split_by_tier(filtered)
        tier_counts = {}
        for tier, tier_df in tier_dfs.items():
            tier_df.to_csv(data_dir / f"filtered_{tier}.csv", index=False)
            tier_counts[tier] = len(tier_df)
            print(f"  {tier}: {len(tier_df)} clonotypes")
    else:
        # Per-donor or per-sample scope. Build long_df if not already.
        if long_df is None:
            from .clonotype import build_clone_sample_long
            long_df = build_clone_sample_long(adata, attribution=config.attribution)
        if resolved_scope == "per-donor":
            subsets = split_clonotypes_by_donor(clonotypes, long_df)
            label = "donor"
        else:
            subsets = split_clonotypes_by_sample(clonotypes, long_df)
            label = "sample"

        per_subset_filtered = []
        tier_counts = {}
        for key, sub_clonotypes in subsets.items():
            if len(sub_clonotypes) == 0:
                continue
            try:
                sub_filtered = _filter_one(sub_clonotypes)
            except TCRsiftValidationError as exc:
                # Narrow exception: filter_clonotypes_threshold raises this
                # when the subset is empty after filtering. Other exceptions
                # (real bugs, OOM, dtype mismatches) propagate so they're not
                # masked under a benign "no clones" message.
                print(
                    f"  WARNING: no clones survived filtering for "
                    f"{label}={key}: {exc}"
                )
                continue
            sub_filtered = sub_filtered.copy()
            sub_filtered[f"scope_{label}"] = key
            per_subset_filtered.append(sub_filtered)

            # Per-subset tier files alongside the existing global naming
            # convention (e.g. filtered_tier1_B1-2.csv for per-donor scope).
            sub_tier_dfs = split_by_tier(sub_filtered)
            safe_key = str(key).replace("/", "_").replace(" ", "_")
            for tier, tier_df in sub_tier_dfs.items():
                tier_df.to_csv(
                    data_dir / f"filtered_{tier}_{safe_key}.csv",
                    index=False,
                )
                tier_counts[f"{tier}_{safe_key}"] = len(tier_df)
                print(f"  {label}={key} {tier}: {len(tier_df)} clonotypes")

        filtered = (
            pd.concat(per_subset_filtered, ignore_index=True)
            if per_subset_filtered
            else pd.DataFrame()
        )

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

    # Per-method ranked CSVs (#20 chunk 2). One file per (donor, method)
    # under data/filtered_by_method/, top-N clones ranked by within-bucket
    # frequency. Skipped when enrichment_method axis isn't populated.
    if long_df is not None:
        from .clonotype import (
            build_method_overlap_matrices,
            build_method_recovery_table,
            build_per_method_rankings,
        )

        rankings = build_per_method_rankings(
            filtered, long_df, top_n=config.output.per_method_top_n
        )
        if rankings:
            method_dir = data_dir / "filtered_by_method"
            method_dir.mkdir(parents=True, exist_ok=True)
            for (donor, method), ranked in rankings.items():
                # Sanitize the filename — encoding/method labels may have
                # path-unfriendly characters (slashes etc.); replace with
                # underscores defensively.
                safe_donor = donor.replace("/", "_").replace(" ", "_")
                safe_method = method.replace("/", "_").replace(" ", "_")
                ranked.to_csv(
                    method_dir / f"{safe_donor}__{safe_method}.csv",
                    index=False,
                )
            print(
                f"  Wrote {len(rankings)} per-method ranked CSVs to "
                f"{method_dir.name}/"
            )

        # Method × method overlap matrix per donor (#27 chunk 3). CSV +
        # heatmap PNG. Skipped when method axis isn't populated.
        similarity = config.output.method_overlap_similarity
        overlap_matrices = build_method_overlap_matrices(
            filtered, long_df, similarity=similarity
        )
        if overlap_matrices:
            for donor, mat in overlap_matrices.items():
                safe_donor = donor.replace("/", "_").replace(" ", "_")
                mat.to_csv(data_dir / f"method_overlap_{safe_donor}.csv")
                if config.output.generate_plots:
                    from .plots import plot_method_overlap

                    plot_method_overlap(
                        mat,
                        plots_dir / f"method_overlap_{safe_donor}.png",
                        similarity=similarity,
                        donor=donor,
                    )
            print(
                f"  Method overlap ({similarity}): "
                f"{len(overlap_matrices)} donor matrices"
            )

        # Method-recovery table + bar plot (#27 chunk 4). Targets the strictest
        # *populated* tier on `filtered` (not a hardcoded tier1, which empties
        # the table for unexpanded cohorts and silently drops the panel, #167);
        # falls back to '*' (all filtered clones) when no tier is populated.
        target_tier = pick_strictest_populated_tier(filtered)
        recovery = build_method_recovery_table(
            filtered, long_df, tier=target_tier
        )
        if not recovery.empty:
            recovery.to_csv(data_dir / "method_recovery.csv", index=False)
            if config.output.generate_plots:
                from .plots import plot_method_recovery

                plot_method_recovery(
                    recovery,
                    plots_dir / "method_recovery.png",
                    tier_label=target_tier,
                )
            print(
                f"  Method recovery ({target_tier}): "
                f"{recovery['method'].nunique()} methods × "
                f"{recovery['donor'].nunique()} donors"
            )

    funnel_counts["Filtered"] = len(filtered)

    if config.output.generate_plots:
        plot_filter(filtered, plots_dir)

    # Step 5: Annotate (if databases provided or cached).
    # Explicit config paths win; otherwise fall back to the data cache
    # (`tcrsift data list` shows what's available).
    from .datacache import cached_path

    cache_dir = getattr(args, "cache_dir", None)
    vdjdb_path = config.annotate.vdjdb_path or cached_path("vdjdb", cache_dir)
    iedb_path = config.annotate.iedb_path or cached_path("iedb", cache_dir)
    cedar_path = config.annotate.cedar_path or cached_path("cedar", cache_dir)
    # No config entry yet — fall back to the cache. Users who downloaded
    # via ``tcrsift data download --db iedb_epitope`` get the override
    # automatically.
    iedb_epitope_path = cached_path("iedb_epitope", cache_dir)
    annotated = filtered
    has_annotation = vdjdb_path or iedb_path or cedar_path
    if has_annotation:
        print("\n[5/7] Annotating with public databases...")
        annotated = annotate_clonotypes(
            filtered,
            vdjdb_path=vdjdb_path,
            iedb_path=iedb_path,
            cedar_path=cedar_path,
            iedb_epitope_path=iedb_epitope_path,
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
    # Selection-shortlist funnel counts (populated when a selection config
    # runs below; consumed by the funnel plot).
    n_selected_clones: int | None = None
    n_non_viral: int | None = None
    has_leaders = (
        config.assemble.alpha_leader is not None or config.assemble.beta_leader is not None
    )
    if config.assemble.single_chain or has_leaders or config.assemble.include_constant:
        print("\n[7/7] Assembling full-length sequences...")
        # Contig sample-name policy (#124). For "sheet" mode, hand the
        # loaded sample sheet to the assembler as a (sample, vdj_dir) frame.
        _sample_name_from = getattr(config.assemble, "sample_name_from", "parent")
        _contigs_dir = config.assemble.contigs_dir

        # Allele fidelity on by default (#201): when nothing explicit was
        # configured but the sample sheet has vdj_dirs, assemble against those
        # contigs (the run already loaded them) so allele-faithful constants +
        # the #187 divergence check are active — instead of silently assembling
        # canonical constants. --contigs-dir/--cellranger-dir remain overrides.
        _explicit_contigs = config.assemble.contigs_dir or getattr(
            config.assemble, "cellranger_dir", None
        )
        _sheet_vdj_dirs = [s.vdj_dir for s in sample_sheet.samples if s.vdj_dir]
        if (
            _sample_name_from == "parent"
            and not _explicit_contigs
            and _sheet_vdj_dirs
            and not getattr(args, "no_sheet_contigs", False)
        ):
            _root = _common_contig_root(_sheet_vdj_dirs)
            if _root:
                _contigs_dir = _root
                _sample_name_from = "sheet"
                print(
                    f"  Allele fidelity: assembling against sample-sheet contigs "
                    f"under {_root} (#201; pass --no-sheet-contigs / --contigs-dir to change)"
                )

        _contig_sheet = (
            pd.DataFrame(
                [
                    {"sample": s.sample, "vdj_dir": s.vdj_dir}
                    for s in sample_sheet.samples
                    if s.vdj_dir
                ]
            )
            if _sample_name_from == "sheet"
            else None
        )
        assembled = assemble_full_sequences(
            til_matched,
            contigs_dir=_contigs_dir,
            cellranger_dir=getattr(config.assemble, "cellranger_dir", None),
            sample_name_from=_sample_name_from,
            sample_sheet=_contig_sheet,
            alpha_leader=config.assemble.alpha_leader,
            beta_leader=config.assemble.beta_leader,
            include_constant=config.assemble.include_constant,
            constant_source=config.assemble.constant_source,
            linker=config.assemble.linker if config.assemble.single_chain else None,
            # #113 allele controls: ``tcrsift run`` uses defaults
            # (``auto``); fine-tune via ``tcrsift assemble`` if needed.
            trac_allele=getattr(config.assemble, "trac_allele", "auto"),
            trbc1_allele=getattr(config.assemble, "trbc1_allele", "auto"),
            trbc2_allele=getattr(config.assemble, "trbc2_allele", "auto"),
            # #116 stop codons: default to dual-stop. Config-driven
            # callers can pass a list via config.assemble.stop_codons.
            stop_codons=tuple(
                getattr(config.assemble, "stop_codons", ("TAA", "TGA"))
            ),
        )
        assembled.to_csv(data_dir / "full_sequences.csv", index=False)
        print(f"  Assembled {len(assembled)} sequences")

        # Config-driven clone selection (#122/#125). When a `selection`
        # block is present, match/rank the assembled clones against the
        # selection rules and emit selected_clones.csv (selection columns
        # first, then full seqs), ordered by global_rank.
        if config.selection.get("rules") and long_df is not None:
            from .selection import select_from_clone_sample_long
            # Gene-expression signatures as synthetic methods (#sig). When
            # selection.signatures is set, score those signatures on the GEX
            # matrix and inject them as extra (clone, sample, method) rows so
            # the rules can use a signature like a sort. Augments a *copy* of
            # the long table — the emitted clone_sample_long.csv is untouched.
            selection_long = long_df
            sig_names = config.selection.get("signatures")
            if sig_names:
                from .signature_methods import (
                    SIGNATURES,
                    build_signature_methods,
                    expression_frame_from_adata,
                )
                genes = sorted({
                    g for n in sig_names for g in SIGNATURES[n].all_genes
                })
                # Extract leniently so build_signature_methods can raise a
                # per-signature error naming the missing genes — a signature
                # is never scored on a partial panel (it'd be a different
                # signature). A missing gene is a hard error, not a skip.
                expr = expression_frame_from_adata(adata, genes, on_missing="ignore")
                obs_sig = adata.obs
                if "CDR3ab" not in obs_sig.columns and {
                    "CDR3_alpha", "CDR3_beta"
                } <= set(obs_sig.columns):
                    obs_sig = obs_sig.assign(
                        CDR3ab=obs_sig["CDR3_alpha"].fillna("")
                        + "_" + obs_sig["CDR3_beta"].fillna("")
                    )
                if not ({"CDR3ab", "sample"} <= set(obs_sig.columns)):
                    raise ValueError(
                        f"selection.signatures={sig_names} requested but "
                        "adata.obs lacks CDR3ab/sample — cannot build "
                        "signature methods."
                    )
                # Raises (on_missing='error') if any signature gene has no
                # expression — never silently degrades the panel.
                sig_long = build_signature_methods(
                    expr, obs_sig, signatures=sig_names,
                    positive_method=config.selection.get(
                        "signature_positive_method", "gap"
                    ),
                )
                # If the cohort has no enrichment-method axis, long_df has no
                # 'method' column; the concat introduces it (real rows get
                # NaN method, dropped downstream) so the signatures become the
                # selection methods. With a method axis, signatures are added
                # alongside the real sorts.
                selection_long = pd.concat(
                    [long_df, sig_long], ignore_index=True
                )
                print(
                    f"  Added {len(sig_long)} signature-method rows "
                    f"({', '.join(sig_names)}) to the selection input"
                )
            # exclude_viral: drop public-DB viral bystanders (is_viral from
            # the annotate step) so they don't enrich into the selection.
            exclude_clones = None
            # Defaults ON: viral bystanders are dropped from selection
            # unless the config explicitly sets exclude_viral: false.
            if config.selection.get("exclude_viral", True) and "is_viral" in til_matched.columns:
                viral_mask = til_matched["is_viral"].fillna(False).astype(bool)
                exclude_clones = set(
                    til_matched.loc[viral_mask, "CDR3ab"].astype(str)
                )
                if exclude_clones:
                    # Surface what was dropped, for audit — don't silently
                    # remove them. Keep the DB-annotation columns when present.
                    audit_cols = [
                        c for c in (
                            "CDR3ab", "is_viral", "antigen", "antigen_species",
                            "epitope", "database", "annotation_category",
                        )
                        if c in til_matched.columns
                    ]
                    til_matched.loc[viral_mask, audit_cols].to_csv(
                        data_dir / "excluded_viral_clones.csv", index=False
                    )
                    n_non_viral = int((~viral_mask).sum())
                    print(
                        f"  Selection excluding {len(exclude_clones)} viral clones "
                        f"-> data/excluded_viral_clones.csv"
                    )
            rules = select_from_clone_sample_long(
                selection_long, config.selection, exclude_clones=exclude_clones,
            )
            if not rules.empty and "CDR3ab" in assembled.columns:
                selected = (
                    rules.merge(assembled, on="CDR3ab", how="left")
                    .sort_values("global_rank")
                    .reset_index(drop=True)
                )
                selected.to_csv(data_dir / "selected_clones.csv", index=False)
                n_selected_clones = len(selected)
                print(
                    f"  Wrote selected_clones.csv: {len(selected)} clones "
                    f"across {rules['selection_rule'].nunique()} rules"
                )

        if config.output.generate_plots:
            plot_assembly(assembled, plots_dir)

        # Generate sequence PDF. When per-method data is available, overlay
        # each clone's per-method evidence + selection rule (#123/#125).
        if config.output.generate_report:
            pdf_annotations = None
            if long_df is not None and "method" in long_df.columns:
                from .selection import (
                    attach_method_tiers,
                    build_clone_method_long,
                    build_pdf_annotations,
                    select_from_clone_sample_long,
                )
                cml = attach_method_tiers(build_clone_method_long(long_df))
                sel_df = (
                    select_from_clone_sample_long(long_df, config.selection)
                    if config.selection.get("rules") else None
                )
                pdf_annotations = build_pdf_annotations(cml, selection_df=sel_df)
            # strict="skip": render the valid clones and drop+warn any that
            # genuinely fail validation, rather than aborting the whole run
            # at this final (cosmetic) step — full_sequences.csv is already
            # written above (#129).
            create_tcr_sequence_pdf(
                assembled,
                output_dir / "tcr_sequences.pdf",
                annotations=pdf_annotations,
                strict="skip",
            )

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
            # Selection shortlist as a sibling funnel; viral exclusion (if
            # any) shows as a "Non-viral" stage before "Selected".
            selected_count=n_selected_clones,
            emit_selected_variant=n_selected_clones is not None,
            non_viral=n_non_viral,
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
        compute_signature_scores_per_clonotype,
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

    # Step 2b: Per-clonotype signature scores (#74). Always runs when
    # the per-cell augmented frame is available — scores are NaN for
    # signatures whose genes weren't in the augmented set.
    if augmented_df is not None:
        print("\nComputing per-clonotype signature scores...")
        sig_df = compute_signature_scores_per_clonotype(
            augmented_df,
            group_col=args.group_col,
            gex_prefix=args.prefix,
            cd8_only=True,
            verbose=args.verbose,
        )
        sig_cols = [c for c in sig_df.columns if c.startswith("signature_")]
        if sig_cols:
            df = df.merge(
                sig_df[[args.group_col, *sig_cols]],
                on=args.group_col,
                how="left",
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


def cmd_select(args):
    """Apply the config-driven multi-route selection language to a long table.

    Standalone front-end (#122) to the selection engine that ``run`` uses
    internally: reads a clone-sample-long CSV (the ``clone_sample_long.csv``
    that ``run`` emits), runs the configured ``selection.rules``, and writes
    the selected clones (selection columns + ``global_rank``). Optionally joins
    per-clone metadata from a clonotypes CSV and emits per-sample tiers.
    """
    import pandas as pd

    from .config import TCRsiftConfig
    from .selection import attach_per_sample_tiers, select_from_clone_sample_long

    setup_logging(args.verbose)

    long_df = pd.read_csv(args.input)
    print(f"Loaded {len(long_df)} (clone, sample) rows from {args.input}")

    # Per-clone scores (ppost_alpha/beta, GEX gene-set scores, …).
    clone_scores = None
    scores_path = getattr(args, "clone_scores", None) or getattr(args, "clonotypes", None)
    if scores_path:
        clone_scores = pd.read_csv(scores_path)

    # PRISM mode (#193): the validated per-condition `top-K freq ∪ top-K PRISM`
    # recipe, matching the downstream selection_analysis exactly. Condition =
    # the per-(clone, sample) method/sample column.
    if getattr(args, "prism", False):
        from .selection import select_freq_prism_per_condition

        if clone_scores is None:
            raise TCRsiftValidationError(
                "--prism needs --clone-scores (per-clone ppost/GEX scores).",
                hint="Run `tcrsift annotate` / add_pgen_ppost first to produce "
                "ppost_alpha/ppost_beta/antigen_response_score/naive_score.",
            )
        cond_col = next(
            (c for c in (args.condition_col, "method", "sample") if c and c in long_df.columns),
            None,
        )
        if cond_col is None:
            raise TCRsiftValidationError(
                "No condition column in the long table.",
                hint="Pass --condition-col (e.g. method or sample).",
            )
        score_cols = [c for c in clone_scores.columns if c != "CDR3ab"]
        feat = long_df.merge(clone_scores[["CDR3ab", *score_cols]], on="CDR3ab", how="left")
        sel = select_freq_prism_per_condition(
            feat, condition_col=cond_col, freq_col="frequency",
            gate=args.gate, top_freq=args.top_freq, top_prism=args.top_prism,
        )
        sel.to_csv(args.output, index=False)
        print(
            f"Wrote {sel['CDR3ab'].nunique()} distinct clones "
            f"({len(sel)} clone×condition rows) by per-condition freq∪PRISM → {args.output}"
        )
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

    # Per-sample tier exposure (#122): label each (clone, sample) row with its
    # within-sample abundance tier — the per-sample analogue of the global tier.
    if getattr(args, "emit_per_sample_tiers", False):
        tiered = attach_per_sample_tiers(long_df)
        tier_path = Path(args.output).with_name("clone_sample_per_sample_tiers.csv")
        tiered.to_csv(tier_path, index=False)
        print(f"Wrote per-sample tiers → {tier_path}")


def cmd_report_bundle(args):
    """Assemble a categorized, multi-run figure bundle PDF (#123)."""
    from .report import bundle_figure_pdf

    setup_logging(args.verbose)
    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    out = bundle_figure_pdf(
        args.run_dirs, args.output, labels=labels,
        title=args.title, subtitle=args.subtitle or "",
    )
    print(f"Wrote figure bundle: {out}")


def cmd_report_selected(args):
    """Assemble a selected clone set into a synthesis-ready PDF + CSV (#188)."""
    import anndata as ad
    import pandas as pd

    from .config import TCRsiftConfig
    from .report import build_selected_report

    setup_logging(args.verbose)
    selected = pd.read_csv(args.selection)
    clonotypes = pd.read_csv(args.clonotypes)
    obs = ad.read_h5ad(args.h5ad).obs if args.h5ad else None
    print(f"Selecting {selected['CDR3ab'].nunique()} clones from {args.clonotypes}")

    assemble_kwargs = {}
    if args.config:
        asm = TCRsiftConfig.from_yaml(args.config).assemble
        assemble_kwargs = dict(
            alpha_leader=asm.alpha_leader, beta_leader=asm.beta_leader,
            include_constant=asm.include_constant, constant_source=asm.constant_source,
            linker=asm.linker, contigs_dir=asm.contigs_dir,
            sample_name_from=asm.sample_name_from, cellranger_dir=asm.cellranger_dir,
        )
    # Per-donor contig overrides (#196). The config's contigs path is per-donor
    # but the config is typically shared, so without these the selected report
    # re-assembles canonical constants with NO contigs — losing donor-allele
    # fidelity + the #187 divergence flags that the contig-wired `run` produced.
    if args.cellranger_dir:
        assemble_kwargs["cellranger_dir"] = args.cellranger_dir
    if args.contigs_dir:
        assemble_kwargs["contigs_dir"] = args.contigs_dir
    if args.sample_name_from:
        assemble_kwargs["sample_name_from"] = args.sample_name_from

    prov_cols = [c for c in selected.columns if c != "CDR3ab"]
    build_selected_report(
        selected, clonotypes, args.output_dir, obs=obs,
        assemble_kwargs=assemble_kwargs, provenance_cols=prov_cols,
        cover=not getattr(args, "no_cover", False),
    )
    print(f"Wrote selected-clones report to {args.output_dir}")


def cmd_report_selected_combine(args):
    """Combine per-donor selected-clone PDFs into one cohort PDF (#202)."""
    from .report import combine_selected_pdfs

    setup_logging(args.verbose)
    labels = [s.strip() for s in args.labels.split(",")] if args.labels else None
    out = combine_selected_pdfs(
        args.pdfs, args.output, labels=labels, title=args.title,
        include_legend=not args.no_legend,
    )
    print(f"Wrote combined selected report: {out}")


def cmd_report_sequence(args):
    """Render a sequence PDF with per-method evidence + selection overlay (#123)."""
    import pandas as pd

    from .plots import create_tcr_sequence_pdf

    setup_logging(args.verbose)
    df = pd.read_csv(args.input)
    if args.sort_by and args.sort_by in df.columns:
        df = df.sort_values(args.sort_by).reset_index(drop=True)

    annotations = None
    if args.per_method_source:
        from .clonotype import build_clone_method_long
        from .selection import attach_method_tiers, build_pdf_annotations

        cml = attach_method_tiers(build_clone_method_long(pd.read_csv(args.per_method_source)))
        sel = df if "selection_rule" in df.columns else None
        annotations = build_pdf_annotations(cml, selection_df=sel)

    create_tcr_sequence_pdf(df, args.output, strict=False, annotations=annotations)
    print(f"Wrote sequence PDF ({len(df)} clones): {args.output}")


def cmd_report_cohort(args):
    """Cross-donor overlap report (Jaccard / cell-weighted / beta-only) (#123)."""
    from .cohort import run_cohort_analysis

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
        donor_dirs, args.output_dir,
        selected_only=args.selected_only,
        include_beta_only=not args.no_beta_only,
        emit_tables=not args.no_tables, emit_plots=not args.no_plots,
    )
    print(f"Cohort overlap across {len(donor_dirs)} donors → {args.output_dir}")
    for name, mat in matrices.items():
        print(f"  cohort_{name}: {mat.shape[0]}x{mat.shape[1]}")


def cmd_report_polished(args):
    """Re-render a run's figures with a publication style profile (#123)."""
    import pandas as pd

    from .plots import plot_clonotypes, plot_filter, set_plot_format, set_polished_style

    setup_logging(args.verbose)
    set_polished_style(args.style)
    set_plot_format("pdf")

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
    p_load.add_argument(
        "--min-mito",
        type=float,
        default=2.0,
        help="Min mito %% FLOOR — drops cells below it (default: 2; set 0 to disable, #168)",
    )
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
        "--min-cd3-reads", type=int, default=3, help="Min CD3 reads (default: 3; set 0 to disable, #172)"
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
    p_clone.add_argument(
        "--doublet-warn-rate",
        type=float,
        default=0.1,
        help="Warn when the multi-chain rate meets/exceeds this fraction "
        "(default: 0.1; 0 disables, #165)",
    )
    p_clone.add_argument(
        "--attribution",
        action="store_true",
        help="Partial-information attribution: weighted fractional clone counts "
        "(alpha-dropout sharing, doublet split, dual-alpha merge). #176",
    )
    p_clone.add_argument(
        "--dual-alpha-min-cells",
        type=int,
        default=2,
        help="Cells sharing a dual-alpha triple to merge as one clone (default: 2)",
    )
    p_clone.add_argument(
        "--share-weighting",
        choices=["proportional", "uniform"],
        default="proportional",
        help="How shared/split weight is apportioned (default: proportional)",
    )
    p_clone.add_argument(
        "--dual-alpha-split",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Split non-merged dual-alpha cells as doublets (default: inherit "
        "doublet_split; --no-dual-alpha-split collapses to primary, #181)",
    )
    p_clone.add_argument(
        "--dual-beta-split",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Split dual-beta cells as doublets (default: inherit doublet_split, #181)",
    )
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
    p_filter.add_argument(
        "--min-timepoints", type=int, default=0,
        help="Min distinct timepoints clone must appear in (#9)",
    )
    p_filter.add_argument(
        "--min-timepoints-per-donor", type=int, default=0,
        help="Min distinct timepoints within at least one donor (#9)",
    )
    p_filter.add_argument(
        "--min-apcs", type=int, default=0,
        help="Min distinct APCs clone must appear across (#9)",
    )
    p_filter.add_argument(
        "--min-apcs-per-donor", type=int, default=0,
        help="Min distinct APCs within at least one donor (#9)",
    )
    p_filter.add_argument(
        "--min-til-cells-per-donor", type=int, default=0,
        help="Min TIL cells of the same donor in which the clone appears (#9)",
    )
    p_filter.add_argument("--plot-filter", action="store_true", help="Generate filter plots")
    p_filter.add_argument("--output-dir", help="Output directory for plots")
    p_filter.add_argument("--verbose", action="store_true", help="Verbose output")
    p_filter.set_defaults(func=cmd_filter)

    # -------------------------------------------------------------------------
    # Select command (#122): config-driven multi-route selection language
    # -------------------------------------------------------------------------
    p_select = subparsers.add_parser(
        "select",
        help="Apply config-driven selection rules to a clone-sample-long table",
    )
    p_select.add_argument(
        "--input", "-i", required=True,
        help="clone_sample_long.csv (one row per clone×sample with method/freq)",
    )
    p_select.add_argument("--output", "-o", required=True, help="Output selected-clones CSV")
    p_select.add_argument(
        "--config", "-c",
        help="YAML config with a 'selection:' block (rules + global_rank). "
        "Not needed with --prism.",
    )
    p_select.add_argument(
        "--prism", action="store_true",
        help="Per-condition top-K freq ∪ top-K PRISM recipe (needs --clone-scores), "
        "matching the validated selection_analysis (#193)",
    )
    p_select.add_argument(
        "--condition-col",
        help="Condition column for --prism (default: method, else sample)",
    )
    p_select.add_argument(
        "--gate", type=float, default=0.001, help="--prism frequency gate (default: 0.001)",
    )
    p_select.add_argument(
        "--top-freq", type=int, default=10, help="--prism top-K by frequency (default: 10)",
    )
    p_select.add_argument(
        "--top-prism", type=int, default=5, help="--prism top-K by PRISM (default: 5)",
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

    # -------------------------------------------------------------------------
    # Report command (#123): unified reporting (cohort figure bundle)
    # -------------------------------------------------------------------------
    p_report = subparsers.add_parser("report", help="Unified reporting (figure bundles)")
    p_report_sub = p_report.add_subparsers(dest="report_command")
    pb = p_report_sub.add_parser(
        "bundle",
        help="Concatenate per-run figures into one categorized cohort PDF",
    )
    pb.add_argument(
        "run_dirs", nargs="+",
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

    ps = p_report_sub.add_parser(
        "selected",
        help="Assemble a selected clone set into a synthesis-ready PDF + CSV",
    )
    ps.add_argument(
        "--selection", required=True,
        help="Selected-clones CSV (CDR3ab + provenance columns, e.g. from `tcrsift select`)",
    )
    ps.add_argument(
        "--clonotypes", required=True,
        help="clonotypes.csv (all clones, for VDJ/genes — incl. rescued/low-tier picks)",
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
        "--contigs-dir", help="Per-donor contigs dir (overrides config, #196)",
    )
    ps.add_argument(
        "--sample-name-from", choices=["parent", "grandparent", "sheet"],
        help="Contig sample-name policy (overrides config, #196)",
    )
    ps.add_argument("--output-dir", "-o", required=True, help="Output directory")
    ps.add_argument(
        "--no-cover", action="store_true",
        help="Skip the construct cover/legend page on the sequence PDF (#202)",
    )
    ps.add_argument("--verbose", action="store_true", help="Verbose output")
    ps.set_defaults(func=cmd_report_selected)

    pscomb = p_report_sub.add_parser(
        "selected-combine",
        help="Combine per-donor selected_clones_sequences.pdf into one cohort PDF",
    )
    pscomb.add_argument(
        "pdfs", nargs="+", help="Per-donor selected_clones_sequences.pdf files",
    )
    pscomb.add_argument("--output", "-o", required=True, help="Output combined PDF")
    pscomb.add_argument("--labels", help="Comma-separated donor labels (default: dir names)")
    pscomb.add_argument(
        "--title", default="Selected clones — cohort", help="Cover-page title",
    )
    pscomb.add_argument(
        "--no-legend", action="store_true", help="Skip the construct legend page",
    )
    pscomb.add_argument("--verbose", action="store_true", help="Verbose output")
    pscomb.set_defaults(func=cmd_report_selected_combine)

    pseq = p_report_sub.add_parser(
        "sequence", help="Render a sequence PDF with per-method evidence overlay",
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
        "cohort", help="Cross-donor overlap report (Jaccard / cell-weighted / beta-only)",
    )
    pcoh.add_argument(
        "--donor", action="append", metavar="NAME=DIR",
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
        "polished", help="Re-render a run's figures with a publication style profile",
    )
    ppol.add_argument("--input", "-i", required=True, help="Run output dir (with data/ + plots)")
    ppol.add_argument("--output-dir", "-o", required=True, help="Output directory")
    ppol.add_argument(
        "--style", choices=["clean-white", "paper"], default="clean-white",
        help="Style profile (default: clean-white)",
    )
    ppol.add_argument("--verbose", action="store_true", help="Verbose output")
    ppol.set_defaults(func=cmd_report_polished)

    # -------------------------------------------------------------------------
    # Annotate command
    # -------------------------------------------------------------------------
    p_annot = subparsers.add_parser("annotate", help="Annotate with public databases")
    p_annot.add_argument("--input", "-i", required=True, help="Input filtered CSV")
    p_annot.add_argument("--output", "-o", required=True, help="Output annotated CSV")
    p_annot.add_argument("--vdjdb", help="Path to VDJdb (overrides cache)")
    p_annot.add_argument("--iedb", help="Path to IEDB receptor table (overrides cache)")
    p_annot.add_argument(
        "--iedb-epitope",
        dest="iedb_epitope",
        help="Path to IEDB epitope_full_v3.csv (overrides cache). "
        "When provided, IEDB receptor antigen/species fields are "
        "replaced by the epitope-table's shorter canonical names.",
    )
    p_annot.add_argument("--cedar", help="Path to CEDAR (overrides cache)")
    p_annot.add_argument(
        "--cache-dir",
        help="Reference-database cache directory (default: $TCRSIFT_DATA_DIR or "
        "$XDG_CACHE_HOME/tcrsift or ~/.cache/tcrsift). See `tcrsift data list`.",
    )
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
    # Data command (reference-database cache management)
    # -------------------------------------------------------------------------
    from .datacache import DATABASES

    p_data = subparsers.add_parser(
        "data", help="Manage cached reference databases (VDJdb, IEDB, CEDAR)"
    )
    p_data_sub = p_data.add_subparsers(dest="subcommand", required=True)

    p_data_list = p_data_sub.add_parser(
        "list", help="Show cached database status and disk usage"
    )
    p_data_list.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_list.set_defaults(func=cmd_data_list)

    p_data_clear = p_data_sub.add_parser(
        "clear", help="Remove cached files for one or all databases"
    )
    p_data_clear.add_argument(
        "--db",
        choices=sorted(DATABASES),
        action="extend",
        nargs="+",
        help=(
            "Database(s) to clear (default: clear all). "
            "Repeatable and/or space-separated: "
            "`--db vdjdb iedb` and `--db vdjdb --db iedb` both work."
        ),
    )
    p_data_clear.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_clear.set_defaults(func=cmd_data_clear)

    p_data_download = p_data_sub.add_parser(
        "download",
        help="Download supported databases (VDJdb, IEDB). CEDAR is manual-only.",
    )
    p_data_download.add_argument(
        "--db",
        choices=sorted(DATABASES),
        action="extend",
        nargs="+",
        help=(
            "Database(s) to download (default: all DBs with a download URL). "
            "Repeatable and/or space-separated: "
            "`--db vdjdb iedb` and `--db vdjdb --db iedb` both work."
        ),
    )
    p_data_download.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the database is already cached",
    )
    p_data_download.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_download.set_defaults(func=cmd_data_download)

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
        "--enrichment-genes",
        type=str,
        default="CXCL13,ENTPD1",
        help="Comma-separated genes for exhaustion/tumor-reactivity enrichment branch.",
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
    asm_leaders.add_argument(
        "--cellranger-dir",
        metavar="DIR",
        default=None,
        help="Raw CellRanger per_sample_outs/ dir; shorthand for "
        "--contigs-dir DIR --sample-name-from grandparent (#124).",
    )
    asm_leaders.add_argument(
        "--sample-name-from",
        choices=["parent", "grandparent", "sheet"],
        default=None,
        help="How to derive each contig's sample name from its FASTA path "
        "(default: parent). 'grandparent' for per_sample_outs/{sample}/vdj_t/.",
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
        "--trac-allele",
        default="auto",
        help="TRAC allele override (#113). 'auto' (default) scores "
        "packaged alleles against contig translation; pass an allele "
        "label like '01' to force a specific canonical.",
    )
    asm_seq.add_argument(
        "--trbc1-allele",
        default="auto",
        help="TRBC1 allele override (#113). 'auto' (default) or an "
        "allele label like '01'.",
    )
    asm_seq.add_argument(
        "--trbc2-allele",
        default="auto",
        help="TRBC2 allele override (#113). 'auto' (default) picks "
        "between *01-protein (E at mature pos 9 — major-allele) and "
        "*03-protein (K at pos 9) based on the donor's contig. Pass "
        "'01' or '03' to force.",
    )
    asm_seq.add_argument(
        "--stop-codons",
        default="TAA,TGA",
        help="Comma-separated stop codons to append to codon-"
        "optimized constants (#116). Default 'TAA,TGA' = two non-"
        "redundant stops (different release factors, reduces "
        "read-through). Pass 'TAA' for single-stop (pre-2.4 "
        "behavior) or '' to omit. Each entry must be one of "
        "TAA/TAG/TGA.",
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
    # Audit alleles subcommand (#119, #120)
    # -------------------------------------------------------------------------
    p_audit = subparsers.add_parser(
        "audit-alleles",
        help="Cohort-level allele audit + novel-allele detection",
    )
    p_audit.add_argument(
        "input",
        help="Path to full_sequences.csv from tcrsift assemble / run",
    )
    p_audit.add_argument(
        "-o", "--output",
        metavar="PATH",
        help="Write novel-allele candidates CSV to this path "
             "(in addition to the text report on stdout)",
    )
    p_audit.add_argument(
        "--min-pct", type=float, default=0.05,
        help="Minimum fraction of in-cohort clones for a "
             "novel-allele candidate (default: 0.05 = 5%%)",
    )
    p_audit.add_argument(
        "--min-v-spread", type=int, default=3,
        help="Minimum distinct V genes the variant must appear "
             "in (default: 3)",
    )
    p_audit.add_argument(
        "--min-samples", type=int, default=2,
        help="Minimum distinct samples the variant must "
             "appear in (default: 2)",
    )
    p_audit.set_defaults(func=cmd_audit_alleles)

    # -------------------------------------------------------------------------
    # Pgen/Ppost command (OLGA + SONIA, optional extra)
    # -------------------------------------------------------------------------
    p_ppost = subparsers.add_parser(
        "ppost",
        help="Compute OLGA Pgen + SONIA Ppost per clonotype "
             "(needs: pip install tcrsift[pgen])",
    )
    p_ppost.add_argument("input", help="Path to a clonotype CSV")
    p_ppost.add_argument(
        "-o", "--output", required=True, metavar="PATH",
        help="Write the Pgen/Ppost-annotated CSV here",
    )
    p_ppost.add_argument(
        "--chain", choices=["alpha", "beta", "both"], default="beta",
        help="Which chain(s) to score (default: beta)",
    )
    p_ppost.add_argument(
        "--cdr3-col", default=None,
        help="CDR3 AA column (default: CDR3_<chain>)",
    )
    p_ppost.add_argument(
        "--v-gene-col", default=None,
        help="V-gene column (default: <chain>_v_gene)",
    )
    p_ppost.add_argument(
        "--j-gene-col", default=None,
        help="J-gene column (default: <chain>_j_gene)",
    )
    p_ppost.add_argument(
        "--freq-col", default="frequency",
        help="Clone-frequency column for the private-candidate flag "
             "(default: frequency; skipped if absent)",
    )
    p_ppost.set_defaults(func=cmd_ppost)

    # -------------------------------------------------------------------------
    # log-pgen command (data-driven k-mer / TCRpeg background)
    # -------------------------------------------------------------------------
    p_lpg = subparsers.add_parser(
        "log-pgen",
        help="Score per-clone log Pgen with a data-driven background "
             "(k-mer default, no extra deps; or --backend tcrpeg)",
    )
    p_lpg.add_argument("input", help="Path to a clonotype CSV")
    p_lpg.add_argument(
        "-o", "--output", required=True, metavar="PATH",
        help="Write the log-Pgen-annotated CSV here",
    )
    p_lpg.add_argument(
        "--chain", choices=["alpha", "beta", "both"], default="beta",
        help="Which chain(s) to score (default: beta)",
    )
    p_lpg.add_argument(
        "--backend", choices=["kmer", "tcrpeg"], default="kmer",
        help="Probability backend (default: kmer, numpy-only)",
    )
    p_lpg.add_argument(
        "--cdr3-col", default=None,
        help="CDR3 AA column (default: CDR3_<chain>)",
    )
    p_lpg.set_defaults(func=cmd_log_pgen)

    # -------------------------------------------------------------------------
    # pgen group: manage Pgen/Ppost background models (fetch / train / annotate)
    # -------------------------------------------------------------------------
    p_pgen = subparsers.add_parser(
        "pgen",
        help="Manage Pgen/Ppost background models (fetch data, train, "
             "annotate); TCRpeg default, auto-trains on first use",
    )
    p_pgen_sub = p_pgen.add_subparsers(dest="pgen_command")

    pf = p_pgen_sub.add_parser("fetch", help="Fetch a neutral observed "
                               "repertoire for Ppost training")
    pf.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pf.add_argument("--url", default=None,
                    help="CSV/TSV (optionally .gz) with a CDR3 column")
    pf.add_argument("--source", default=None,
                    choices=["healthy-pbmc", "ots", "ireceptor"],
                    help="'healthy-pbmc' auto-downloads the 5 pooled 10x sets; "
                         "'ots'/'ireceptor' are informational pointers")
    pf.add_argument("--seq-col", default="seq", help="CDR3 column (default seq)")
    pf.add_argument("--v-col", default=None,
                    help="V-gene column to filter chain (TRA*/TRB*)")
    pf.set_defaults(func=cmd_pgen_fetch)

    pt = p_pgen_sub.add_parser("train", help="Train + cache background models")
    pt.add_argument("--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer")
    pt.add_argument("--role", choices=["pgen", "ppost", "both"], default="both")
    pt.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pt.add_argument("--url", default=None, help="observed repertoire for ppost")
    pt.add_argument("--n-olga", type=int, default=200_000,
                    help="synthetic seqs for pgen training (default 200k)")
    pt.add_argument("--epochs", type=int, default=30, help="TCRpeg epochs")
    pt.add_argument("--device", default=None, help="torch device (auto: mps/cuda/cpu)")
    pt.set_defaults(func=cmd_pgen_train)

    pa = p_pgen_sub.add_parser("annotate", help="Annotate a TCR CSV with "
                               "pgen/ppost (α/β); auto-trains on first use")
    pa.add_argument("input", help="Path to a clonotype CSV")
    pa.add_argument("-o", "--output", required=True, metavar="PATH")
    pa.add_argument("--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer")
    pa.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pa.add_argument("--no-auto-train", action="store_true",
                    help="error instead of training a missing TCRpeg model")
    pa.set_defaults(func=cmd_pgen_annotate)

    ps = p_pgen_sub.add_parser("status", help="Show cached training data + models")
    ps.set_defaults(func=cmd_pgen_status)

    p_pgen.set_defaults(func=lambda _a: p_pgen.print_help())

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
        "--min-mito",
        type=float,
        dest="min_mito_pct",
        help="Min mito %% FLOOR — drops cells below it (default: 2; set 0 to disable, #168)",
    )
    load_group.add_argument(
        "--max-mito", type=float, dest="max_mito_pct", help="Max mito %% (default: 8)"
    )

    # Phenotype step parameters
    pheno_group = p_run.add_argument_group("Phenotype options")
    pheno_group.add_argument(
        "--cd4-cd8-ratio", type=float, help="Ratio for confident calls (default: 3.0)"
    )
    pheno_group.add_argument(
        "--min-cd3-reads", type=int, help="Min CD3 reads (default: 3; set 0 to disable, #172)"
    )

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
    clone_group.add_argument(
        "--doublet-warn-rate",
        type=float,
        help="Warn when the multi-chain rate meets/exceeds this fraction "
        "(default: 0.1; 0 disables, #165)",
    )

    # Attribution options (#176). Opt-in; flags default to None so they don't
    # clobber config values via merge_with_args when left unset.
    attr_group = p_run.add_argument_group("Attribution options")
    attr_group.add_argument(
        "--attribution",
        action="store_true",
        default=None,
        help="Partial-information attribution: weighted fractional clone counts "
        "(off by default, #176)",
    )
    attr_group.add_argument(
        "--dual-alpha-min-cells",
        type=int,
        help="Cells sharing a dual-alpha triple to merge as one clone (default: 2)",
    )
    attr_group.add_argument(
        "--share-weighting",
        choices=["proportional", "uniform"],
        help="How shared/split weight is apportioned (default: proportional)",
    )
    attr_group.add_argument(
        "--dual-alpha-split",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Split non-merged dual-alpha cells as doublets (default: inherit "
        "doublet_split). --no-dual-alpha-split collapses them to the primary pair (#181)",
    )
    attr_group.add_argument(
        "--dual-beta-split",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Split dual-beta cells as doublets (default: inherit doublet_split, #181)",
    )

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
    filter_group.add_argument(
        "--min-timepoints", type=int,
        help="Min distinct timepoints clone must appear in",
    )
    filter_group.add_argument(
        "--min-timepoints-per-donor", type=int,
        help="Min distinct timepoints within at least one donor",
    )
    filter_group.add_argument(
        "--min-apcs", type=int,
        help="Min distinct APCs clone must appear across",
    )
    filter_group.add_argument(
        "--min-apcs-per-donor", type=int,
        help="Min distinct APCs within at least one donor",
    )
    filter_group.add_argument(
        "--min-til-cells-per-donor", type=int,
        help="Min TIL cells of the same donor in which the clone appears",
    )
    filter_group.add_argument(
        "--fdr-scope",
        choices=["auto", "global", "per-donor", "per-sample"],
        help="FDR null scope (default: auto). 'auto' = per-donor when "
        "n_donors>1 and the sheet doesn't set donors_share_antigen, "
        "else global. (#26)",
    )
    filter_group.add_argument(
        "--donors-share-antigen",
        action="store_true",
        default=None,
        help="Mark the cohort as having donors that share antigen + MHC + "
        "cohort. Locks fdr_scope='auto' resolution to 'global'.",
    )

    # Annotate step parameters
    annot_group = p_run.add_argument_group("Annotation options")
    annot_group.add_argument("--vdjdb", dest="vdjdb_path", help="Path to VDJdb (overrides cache)")
    annot_group.add_argument("--iedb", dest="iedb_path", help="Path to IEDB (overrides cache)")
    annot_group.add_argument("--cedar", dest="cedar_path", help="Path to CEDAR (overrides cache)")
    annot_group.add_argument(
        "--cache-dir",
        help="Reference-database cache dir (default: $TCRSIFT_DATA_DIR or "
        "$XDG_CACHE_HOME/tcrsift or ~/.cache/tcrsift). See `tcrsift data list`.",
    )
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
        "--cellranger-dir",
        metavar="DIR",
        default=None,
        help="Raw CellRanger per_sample_outs/ dir; shorthand for "
        "--contigs-dir DIR --sample-name-from grandparent (#124).",
    )
    asm_group.add_argument(
        "--sample-name-from",
        choices=["parent", "grandparent", "sheet"],
        default=None,
        help="How to derive each contig's sample name from its FASTA path "
        "(default: parent). 'grandparent' for per_sample_outs/{sample}/vdj_t/.",
    )
    asm_group.add_argument(
        "--no-sheet-contigs",
        action="store_true",
        help="Don't auto-assemble against the sample sheet's vdj_dir contigs; "
        "use canonical constants instead (allele fidelity off, #201).",
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
    out_group.add_argument(
        "--plot-format",
        choices=["png", "pdf", "svg"],
        help="Figure format (default: png); pdf/svg emit a vector file "
        "alongside a PNG for publication figures (#169)",
    )
    out_group.add_argument(
        "--emit-clone-sample-long",
        choices=["auto", "always", "never"],
        help="Emit data/clone_sample_long.csv: 'auto' (default) writes it "
        "when the sheet has >=2 samples; 'always' / 'never' force.",
    )
    out_group.add_argument(
        "--per-method-top-n",
        type=int,
        help="Per-(donor, method) ranked CSVs: top-N clones each. "
        "Default 100. Skipped when enrichment_method axis isn't populated.",
    )
    out_group.add_argument(
        "--method-overlap-similarity",
        choices=["jaccard", "dice", "count"],
        help="Similarity metric for the method × method overlap matrix "
        "(default: jaccard).",
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
    # Cohort command (#125) — multi-donor clone-overlap analysis
    # -------------------------------------------------------------------------
    p_cohort = subparsers.add_parser(
        "cohort", help="Multi-donor clone-overlap analysis across run outputs"
    )
    p_cohort.add_argument(
        "--donor", action="append", metavar="NAME=DIR",
        help="Donor as NAME=run_output_dir (repeatable; need >=2).",
    )
    p_cohort.add_argument(
        "--output-dir", "-o", required=True,
        help="Directory for cohort overlap tables / plots.",
    )
    p_cohort.add_argument(
        "--selected-only", action="store_true",
        help="Overlap selected_clones.csv instead of clonotypes.csv.",
    )
    p_cohort.add_argument(
        "--no-beta-only", action="store_true",
        help="Skip the beta-only overlap matrices.",
    )
    p_cohort.add_argument(
        "--no-tables", action="store_true", help="Don't write CSV matrices.",
    )
    p_cohort.add_argument(
        "--emit-plots", action="store_true", help="Render heatmap PDFs.",
    )
    p_cohort.set_defaults(func=cmd_cohort)

    # Generate config command
    # -------------------------------------------------------------------------
    p_config = subparsers.add_parser("generate-config", help="Generate example config file")
    p_config.add_argument(
        "--output", "-o", default=None,
        help="Output YAML file (default: print to stdout; refuses to overwrite "
        "an existing file without --force, #195)",
    )
    p_config.add_argument(
        "--force", action="store_true", help="Overwrite an existing output file",
    )
    p_config.set_defaults(func=cmd_generate_config)

    return parser


def cmd_cohort(args):
    """Run multi-donor clone-overlap analysis (#125)."""
    from .cohort import run_cohort_analysis

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


def cmd_generate_config(args):
    """Generate an example configuration file (#195: no silent overwrite)."""
    import tempfile

    from .config import generate_example_config

    # Default to stdout so `tcrsift generate-config | less` is safe and never
    # clobbers a working config (#195).
    if args.output is None:
        with tempfile.NamedTemporaryFile("w+", suffix=".yaml", delete=False) as tf:
            tmp = Path(tf.name)
        generate_example_config(tmp)
        print(tmp.read_text(), end="")
        tmp.unlink()
        return

    out_path = Path(args.output)
    if out_path.exists() and not args.force:
        print(
            f"ERROR: {out_path} already exists — refusing to overwrite. "
            "Use --force to overwrite, -o PATH for a different file, or omit -o "
            "to print the example to stdout."
        )
        sys.exit(1)
    generate_example_config(out_path)
    print(f"Generated example config: {out_path}")
    print("\nYou can customize this file and use it with:")
    print(f"  tcrsift run --config {out_path} --sample-sheet samples.csv -o output/")


def main(args=None):
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args(args)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        rc = args.func(args)
    except Exception as e:
        logging.error(f"Error: {e}")
        if hasattr(args, "verbose") and args.verbose:
            raise
        sys.exit(1)
    if rc:
        sys.exit(rc)


if __name__ == "__main__":
    main()
