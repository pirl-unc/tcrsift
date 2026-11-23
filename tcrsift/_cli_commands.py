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

"""Subcommand handler implementations for the TCRsift CLI."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .config import load_config_with_args
from .pipeline import run_pipeline
from .validation import (
    TCRsiftValidationError,
    validate_annotate_gex_args,
    validate_assemble_args,
    validate_file_exists,
    validate_run_args,
)


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


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
        species=getattr(args, "species", "auto"),
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
    import pandas as pd

    from .filter import (
        DEFAULT_FDR_TIERS,
        filter_clonotypes,
        get_filter_summary,
        resolve_fdr_tiers_for_method,
        resolve_filter_mode_kwargs,
        split_by_tier,
    )
    from .plots import plot_filter

    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    # Parse FDR tiers. fdr_tiers only drives logistic tiering; under
    # method='threshold' it is inert (#170). Resolve so the default tiers left
    # untouched don't trigger a spurious "ignored" warning on every run — only
    # a genuine customization warns. Mirrors cmd_run.
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
        stim_fdr_threshold=getattr(args, "stim_fdr_threshold", 0.1),
        stim_log2_fc_threshold=getattr(args, "stim_log2_fc_threshold", 1.0),
        stim_delta_threshold=getattr(args, "stim_delta_threshold", 0.5),
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
    from .datacache import ensure_cached
    from .plots import plot_annotations

    setup_logging(args.verbose)

    clonotypes = pd.read_csv(args.input)
    print(f"Loaded {len(clonotypes)} clonotypes from {args.input}")

    cache_dir = getattr(args, "cache_dir", None)
    auto_dl = bool(getattr(args, "auto_download", False))
    annotated = annotate_clonotypes(
        clonotypes,
        vdjdb_path=args.vdjdb or ensure_cached("vdjdb", cache_dir, auto_download=auto_dl),
        iedb_path=args.iedb or ensure_cached("iedb", cache_dir, auto_download=auto_dl),
        cedar_path=args.cedar or ensure_cached("cedar", cache_dir, auto_download=auto_dl),
        iedb_epitope_path=(
            getattr(args, "iedb_epitope", None)
            or ensure_cached("iedb_epitope", cache_dir, auto_download=auto_dl)
        ),
        match_by=args.match_by if hasattr(args, "match_by") else "CDR3ab",
        match_strictness=getattr(args, "database_match", None),
        match_mode=getattr(args, "match_mode", None) or "exact",
        max_distance=(
            args.max_distance if getattr(args, "max_distance", None) is not None else 1
        ),
        exclude_viral=args.exclude_viral,
        flag_only=args.flag_only,
        species=getattr(args, "species", "auto") or "auto",
        host_species_filter=getattr(args, "host_species_filter", "all") or "all",
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
    from .validation import validate_match_til_args  # TCRsiftValidationError is module-level

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
    from .validation import validate_match_til_args  # TCRsiftValidationError is module-level

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
# TIL Prioritization Commands
# =============================================================================


def cmd_til_prioritize(args):
    """Prioritize clones from a standard VDJ+GEX sample sheet."""
    from .til_prioritize import run_til_prioritize

    setup_logging(args.verbose)
    _, candidates = run_til_prioritize(args)
    print(f"Wrote {len(candidates)} candidates and the full audit table to {args.output_dir}")


def cmd_til_select(args):
    """Prioritize promising TIL clonotypes from multi-timepoint VDJ+GEX data."""
    from .til_select import run_til_select

    setup_logging(args.verbose)
    master_df = run_til_select(args)
    n_final = int(master_df["is_candidate_tumor_reactive"].sum())
    print(f"Saved harmonized table: {args.out_table}")
    print(f"Final candidate clonotypes: {n_final}")
    if getattr(args, "dominant_specificity", False) and "is_dominant_specific" in master_df.columns:
        n_dom = int(master_df["is_dominant_specific"].sum())
        print(f"Dominant-specificity (condition-concentrated) clones: {n_dom}")
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
        leader_fallback=getattr(args, "leader_fallback", None) or "germline",
        force_alpha_leader=getattr(args, "force_alpha_leader", None)
        or getattr(args, "force_leader", None),
        force_beta_leader=getattr(args, "force_beta_leader", None)
        or getattr(args, "force_leader", None),
        secondary_alpha_leader=getattr(args, "secondary_alpha_leader", None)
        or getattr(args, "leader_secondary", None),
        secondary_beta_leader=getattr(args, "secondary_beta_leader", None)
        or getattr(args, "leader_secondary", None),
        codon_optimization=getattr(args, "codon_optimization", None) or "focal",
        max_codon_repeats=(
            1 if getattr(args, "max_codon_repeats", None) is None
            else args.max_codon_repeats
        ),
        avoid_enzymes=tuple(getattr(args, "avoid_enzymes", None) or ()),
        species=getattr(args, "species", "auto") or "auto",
        constant_species=getattr(args, "constant_species", None),
        codon_species=getattr(args, "codon_species", None),
        germline_vdj_fallback=bool(getattr(args, "germline_vdj_fallback", False)),
        data_dir=getattr(args, "cache_dir", None),
        auto_download_germline=bool(getattr(args, "auto_download_germline", False)),
    )

    # Fail closed on unverified (canonical-fallback) constructs unless the
    # user explicitly opted in (#243/#244). Raises TCRsiftValidationError when
    # any construct's junction/allele couldn't be contig-verified.
    from .assemble import enforce_contig_fidelity

    fidelity = enforce_contig_fidelity(
        assembled,
        allow_canonical_fallback=getattr(args, "allow_canonical_fallback", False),
        context="assemble",
    )
    if fidelity:
        print(f"WARNING: {fidelity}")

    # Validate (per-row warnings + aggregate QC report).
    warnings = validate_sequences(
        assembled, species=getattr(args, "species", "auto") or "auto"
    )
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


def cmd_log_pgen(args):
    """Score per-clone log Pgen with a data-driven background (seqprob).

    Reads a clonotype CSV and adds ``log_pgen_<chain>`` using the shipped
    OLGA-generated k-mer background (default) or a TCRpeg model
    (``--backend tcrpeg``). This is the data-driven replacement for the
    OLGA/SONIA runtime path: no GPL deps for the default k-mer backend.
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
                df,
                chain=chain,
                cdr3_col=cdr3_col,
                backend=args.backend,
                species=getattr(args, "species", "human") or "human",
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

    sp = getattr(args, "species", "human") or "human"
    # The pooled 5×10x healthy-PBMC reference downloads both chains at once (human).
    if args.source == "healthy-pbmc" and not args.url and sp == "human":
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
                species=sp,
            )
            print(f"  [{chain}] cached → {path}")
        except (FileNotFoundError, ValueError) as exc:
            print(f"  [{chain}] {exc}")


def cmd_pgen_train(args):
    """Train + cache a custom Ppost model; packaged Pgen is not trainable."""
    from .pgen_models import train_model

    sp = getattr(args, "species", "human") or "human"
    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    roles = ["pgen", "ppost"] if args.role == "both" else [args.role]
    for chain in chains:
        for role in roles:
            try:
                path = train_model(
                    chain, backend=args.backend, role=role, url=args.url,
                    epochs=args.epochs, device=args.device, species=sp,
                )
                print(f"  trained {args.backend} {role} {chain} ({sp}) → {path}")
            except (FileNotFoundError, ImportError, ValueError) as exc:
                print(f"  [{role} {chain}] skipped: {exc}")


def cmd_pgen_annotate(args):
    """Annotate a TCR CSV with pgen/ppost (α/β); auto-trains on first use."""
    import pandas as pd

    from .annotate_tcrs import add_pgen_ppost

    sp = getattr(args, "species", "auto") or "auto"
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")
    chains = ("alpha", "beta") if args.chain == "both" else (args.chain,)
    df = add_pgen_ppost(df, chains=chains, backend=args.backend,
                        auto_train=not args.no_auto_train, species=sp)
    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output} (backend={args.backend})")


def cmd_pgen_status(args):
    """Show cached Pgen/Ppost training data + trained models."""
    from .pgen_models import cache_dir, cached_model_file, cached_repertoire_file

    species_list = (
        [args.species]
        if getattr(args, "species", None) in ("human", "mouse")
        else ["human", "mouse"]
    )
    print(f"seqprob cache: {cache_dir()}")
    print("Training repertoires:")
    for sp in species_list:
        for chain in ("alpha", "beta"):
            p = cached_repertoire_file(chain, species=sp)
            if sp == "human" or p.is_file():
                print(f"  {sp} {chain}: {'present' if p.is_file() else 'absent'}  ({p})")
    print("Models:")
    for sp in species_list:
        for backend in ("kmer", "kmer_gene", "tcrpeg"):
            for role in ("pgen", "ppost"):
                for chain in ("alpha", "beta"):
                    p = cached_model_file(backend, role, chain, species=sp)
                    if p.is_file():
                        print(f"  {sp} {backend} {role} {chain}: {p}")


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

    from .datacache import DATABASES

    missing = [e for e in entries if not e.present]
    auto = [e for e in missing if DATABASES[e.name].download_url is not None]
    manual = [e for e in missing if DATABASES[e.name].download_url is None]
    if auto:
        print()
        names = " ".join(e.name for e in auto)
        print(f"Missing — fetch with: tcrsift data download --db {names}")
    if manual:
        print()
        print("Missing (manual only) — download and place under the cache:")
        for e in manual:
            print(f"  {e.name}: {e.source} → {e.path}")


def cmd_data_clear(args):
    """Remove cached files for one or more managed databases / groups."""
    from .datacache import clear_cache, resolve_db_selection

    # Expand --db names + --group aliases (default = all). Order preserved so
    # the user's first mention determines report order.
    targets = resolve_db_selection(args.db, getattr(args, "group", None))
    removed_total: list = []
    for name in targets:
        removed_total.extend(clear_cache(db=name, data_dir=args.cache_dir))

    if not removed_total:
        print("Nothing to remove.")
        return
    for path in removed_total:
        print(f"Removed {path}")


def cmd_data_download(args):
    """Download one or more supported reference databases / groups into the cache."""
    from .datacache import DownloadError, download_database, resolve_db_selection

    setup_logging(verbose=True)
    # Expand --db names + --group aliases. Default (neither) = all auto-
    # downloadable DBs; group/default expansions skip manual-only DBs (CEDAR),
    # while an explicit `--db cedar` is kept so its manual instruction surfaces.
    targets = resolve_db_selection(
        args.db, getattr(args, "group", None), downloadable_only=True
    )
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




def cmd_run(args):
    """Run the complete TCRsift pipeline."""
    validate_run_args(args)
    config = load_config_with_args(args)
    if getattr(args, "no_pgen_ppost", False):
        config.clonotype.add_pgen_ppost = False
    setup_logging(config.verbose)

    run_pipeline(
        args.sample_sheet,
        args.output_dir,
        config=config,
        til_sample_specs=getattr(args, "til_sample", None),
        cache_dir=getattr(args, "cache_dir", None),
        no_assemble=bool(getattr(args, "no_assemble", False)),
        no_sheet_contigs=bool(getattr(args, "no_sheet_contigs", False)),
        allow_canonical_fallback=bool(getattr(args, "allow_canonical_fallback", False)),
        report_name=getattr(args, "report_name", None),
    )


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
    gex_species = getattr(args, "species", "auto") or "auto"
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

    # Step 2b: Per-clonotype signature scores (#74).
    if augmented_df is not None:
        from .species import resolve_species

        resolved_gex_sp = gex_species
        if resolved_gex_sp == "auto":
            try:
                import scanpy as sc

                from .species import detect_adata_species

                resolved_gex_sp = (
                    detect_adata_species(sc.read_10x_h5(str(args.gex_file)))
                    or resolve_species("auto", df=augmented_df, warn_mismatch=False)
                )
            except Exception:
                resolved_gex_sp = resolve_species("auto", df=augmented_df, warn_mismatch=False)
        else:
            resolved_gex_sp = resolve_species(gex_species, df=augmented_df, warn_mismatch=False)
        print("\nComputing per-clonotype signature scores...")
        sig_df = compute_signature_scores_per_clonotype(
            augmented_df,
            group_col=args.group_col,
            gex_prefix=args.prefix,
            cd8_only=True,
            species=resolved_gex_sp,
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
        return 1

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
        # If --condition-col was given explicitly, require it — don't silently
        # fall through to method/sample, which would produce a plausible but
        # wrong per-axis selection on a typo or wrong file.
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
                (c for c in ("method", "sample") if c in long_df.columns), None,
            )
            if cond_col is None:
                raise TCRsiftValidationError(
                    "No condition column in the long table.",
                    hint="Pass --condition-col (e.g. method or sample).",
                )
        score_cols = [c for c in clone_scores.columns if c != "CDR3ab"]
        feat = long_df.merge(clone_scores[["CDR3ab", *score_cols]], on="CDR3ab", how="left")
        # Non-PRISM tie-break for boundary ties (#221/#223). Precedence:
        # --tie-break flag > config selection.tie_break > DEFAULT_TIE_BREAK.
        # Its columns (merged_alpha_partners/TRA_1_umis/TRB_1_umis/
        # Tcell_type_purity) come in via --clonotypes; absent ones are skipped.
        tie_break = None
        if args.config:
            tie_break = (TCRsiftConfig.from_yaml(args.config).selection or {}).get("tie_break")
        if getattr(args, "tie_break", None) is not None:
            tie_break = [k.strip() for k in args.tie_break.split(",") if k.strip()]
        # Coverage rescue for thin conditions (#228).
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
            feat, condition_col=cond_col, freq_col="frequency",
            gate=args.gate, top_freq=args.top_freq, top_prism=args.top_prism,
            tie_break=tie_break, rescue_target=rescue_target,
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
        # Context plots next to the selection output (#248/#249): the 2-D
        # freq×PRISM selection space and the per-component distributions —
        # both need the PRISM score for ALL gated candidates, which only the
        # select step has.
        try:
            from .plots import plot_freq_prism_scatter, plot_prism_components
            from .selection import prism_candidates

            cand = prism_candidates(
                feat, condition_col=cond_col, freq_col="frequency",
                gate=args.gate, top_freq=args.top_freq, top_prism=args.top_prism,
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
            from .selection import freq_prism_grid

            grid = freq_prism_grid(
                feat, condition_col=cond_col, freq_col="frequency", gate=args.gate,
            )
            grid_csv = Path(args.output).with_name("freq_prism_grid.csv")
            grid.to_csv(grid_csv, index=False)
            print(f"Wrote freq×PRISM trade-off grid → {grid_csv}")
            try:
                from .plots import plot_freq_prism_grid

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


def cmd_report_cohort_figures(args):
    """Cross-donor all-figures bundle with one section per donor (#262 follow-up)."""
    from .report import bundle_figure_pdf, resolve_report_name

    setup_logging(args.verbose)
    labels = [resolve_report_name(d) for d in args.run_dirs]
    out = bundle_figure_pdf(
        args.run_dirs, args.output, labels=labels,
        title="TCRsift — cross-donor figures",
        subtitle=f"{len(args.run_dirs)} donors",
    )
    print(f"Wrote cross-donor figure bundle: {out}")


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

    # Public-DB (IEDB/CEDAR/VDJdb) match annotation to surface on the deliverable
    # (#selected-anno): an explicit --annotations file, else the clonotypes frame
    # itself when it already carries the db_* columns (e.g. annotated.csv passed
    # as --clonotypes). None → the columns are emitted empty.
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
        )
    # Per-donor contig + species overrides (#196).
    if args.cellranger_dir:
        assemble_kwargs["cellranger_dir"] = args.cellranger_dir
    if args.contigs_dir:
        assemble_kwargs["contigs_dir"] = args.contigs_dir
    if args.sample_name_from:
        assemble_kwargs["sample_name_from"] = args.sample_name_from
    if getattr(args, "species", None) is not None:
        assemble_kwargs["species"] = args.species
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

    prov_cols = [c for c in selected.columns if c != "CDR3ab"]
    build_selected_report(
        selected, clonotypes, args.output_dir, obs=obs,
        assemble_kwargs=assemble_kwargs, provenance_cols=prov_cols,
        cover=not getattr(args, "no_cover", False),
        allow_canonical_fallback=getattr(args, "allow_canonical_fallback", False),
        annotations=annotations,
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

    from .plots import (
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


# =============================================================================
# Cells command group (qc / embed / annotate) — the single-cell GEX path (#342)
# =============================================================================

# Sentinel distinguishing "threshold not passed on the CLI" (let cell_qc_funnel
# use its own default) from an explicit value — including ``none`` (disable the
# bound). argparse never applies ``type=`` to a default, so this survives untouched.
_UNSET_CLI = object()

_CELLS_QC_THRESHOLDS = (
    ("min_genes", "min_genes", "Min genes per cell"),
    ("max_genes", "max_genes", "Max genes per cell"),
    ("min_counts", "min_counts", "Min UMI counts"),
    ("max_counts", "max_counts", "Max UMI counts"),
    ("min_mito", "min_mito_pct", "Min mito %% FLOOR (drops low-mito ambient)"),
    ("max_mito", "max_mito_pct", "Max mito %% (drops dying cells)"),
)


def _cells_threshold(value):
    """Parse a qc-threshold CLI value: ``none``/``off``/``disable`` → ``None``
    (disable the bound, per cell_qc_funnel's #333 semantics); else a number
    (int when integer-valued, so waterfall reasons read cleanly)."""
    s = str(value).strip().lower()
    if s in ("none", "off", "disable", "disabled"):
        return None
    try:
        f = float(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"{value!r} must be a number or 'none' (to disable the bound)"
        ) from e
    return int(f) if f.is_integer() else f


def cmd_cells_qc(args):
    """Per-cell QC + doublet funnel on an AnnData (cell_qc_funnel, #314/#342)."""
    import anndata as ad

    from .qc import cell_qc_funnel, lineage_programs_for_species
    from .species import resolve_species

    setup_logging(args.verbose)
    adata = ad.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells from {args.input}")

    explicit_sp = getattr(args, "species", None)
    resolved_sp = resolve_species(
        explicit_sp or getattr(adata, "uns", {}).get("species", "auto"),
        adata=adata,
        warn_mismatch=bool(explicit_sp not in (None, "auto")),
        context="cells qc",
    )
    kwargs = {
        param: getattr(args, cli_name)
        for cli_name, param, _ in _CELLS_QC_THRESHOLDS
        if getattr(args, cli_name) is not _UNSET_CLI
    }
    kwargs["species"] = resolved_sp
    config = None
    if args.pbmc_preset:
        from .config import LoadConfig

        config = LoadConfig()
    if args.solid_tumor:
        kwargs["lineage_sets"] = lineage_programs_for_species(
            resolved_sp, solid_tumor=True
        )

    out, waterfall = cell_qc_funnel(adata, config=config, **kwargs)
    out.uns["species"] = resolved_sp
    print("\nQC waterfall:")
    print(waterfall.to_string(index=False))
    out.write_h5ad(args.output)
    print(f"\nSaved {out.n_obs} cells to {args.output}")
    if args.waterfall_csv:
        waterfall.to_csv(args.waterfall_csv, index=False)
        print(f"Waterfall written to {args.waterfall_csv}")


def cmd_cells_embed(args):
    """Variance-stabilized embedding → Leiden clusters (embed_cells, #311/#342)."""
    import anndata as ad

    from .config import EmbedConfig
    from .embed import embed_cells
    from .species import resolve_species

    setup_logging(args.verbose)
    adata = ad.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} cells from {args.input}")
    explicit_sp = getattr(args, "species", None)
    resolved_sp = resolve_species(
        explicit_sp or getattr(adata, "uns", {}).get("species", "auto"),
        adata=adata,
        warn_mismatch=bool(explicit_sp not in (None, "auto")),
        context="cells embed",
    )
    adata.uns["species"] = resolved_sp

    config = EmbedConfig(
        enabled=True,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        batch_key=args.batch_key,
        leiden_resolution=args.resolution,
        n_neighbors=args.n_neighbors,
        seed=args.seed,
    )
    atlas = embed_cells(adata, config)
    atlas.uns["species"] = resolved_sp
    atlas.write_h5ad(args.output)
    print(
        f"Embedded {atlas.n_obs} cells: "
        f"{atlas.obs['leiden'].nunique()} Leiden clusters → {args.output}"
    )


def cmd_cells_annotate(args):
    """Per-cluster cell-type + T/B-state annotation (annotate_cells, #312/#342)."""
    import anndata as ad

    from .annotate_cells import annotate_cells
    from .signatures import PBMC_CULTURE_TYPES, cell_type_signatures_for_species
    from .species import resolve_species

    setup_logging(args.verbose)
    atlas = ad.read_h5ad(args.input)
    print(f"Loaded {atlas.n_obs} cells from {args.input}")

    explicit_sp = getattr(args, "species", None)
    resolved_sp = resolve_species(
        explicit_sp or getattr(atlas, "uns", {}).get("species", "auto"),
        adata=atlas,
        warn_mismatch=bool(explicit_sp not in (None, "auto")),
        context="cells annotate",
    )
    reference = (
        cell_type_signatures_for_species(resolved_sp, solid_tumor=True)
        if args.solid_tumor
        else None
    )
    allowed_types = PBMC_CULTURE_TYPES if args.pbmc_culture else None
    annotate_cells(
        atlas,
        leiden_col=args.leiden_col,
        reference=reference,
        allowed_types=allowed_types,
        min_subtype_frac=args.min_subtype_frac,
        add_suffix=not args.no_suffix,
        species=resolved_sp,
    )
    atlas.uns["species"] = resolved_sp
    atlas.write_h5ad(args.output)
    print("\nPhenotype counts:")
    print(atlas.obs["phenotype"].value_counts().to_string())
    print(f"\nSaved annotated atlas to {args.output}")





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


