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

"""Data ingestion, phenotyping, and clonotype aggregation CLI subcommands.

Commands: ``load``, ``load-sct``, ``phenotype``, ``clonotype``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .common import setup_logging


def cmd_load(args):
    """Load CellRanger outputs from sample sheet."""
    from ..loader import load_samples
    from ..plots import plot_qc

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

    adata.write_h5ad(args.output)
    print(f"Saved {len(adata)} cells to {args.output}")

    if args.plot_qc:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_qc(adata, output_dir)
        print(f"QC plots saved to {output_dir}")


def cmd_load_sct(args):
    """Load TCR data from SCT platform."""
    from ..sct import aggregate_sct, load_sct

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

    if args.aggregate:
        df = aggregate_sct(df, verbose=True)

    df.to_csv(args.output, index=False)
    print(f"\nSaved {len(df)} {'clonotypes' if args.aggregate else 'cells'} to {args.output}")


def cmd_phenotype(args):
    """Classify T cells as CD4/CD8."""
    import anndata as ad

    from ..phenotype import (
        apply_cd3_gate,
        get_phenotype_summary,
        phenotype_cells,
        validate_phenotype_vs_expected,
    )
    from ..plots import plot_phenotype

    setup_logging(args.verbose)

    adata = ad.read_h5ad(args.input)
    print(f"Loaded {len(adata)} cells from {args.input}")

    adata = phenotype_cells(
        adata,
        cd4_cd8_ratio=args.cd4_cd8_ratio,
        min_cd3_reads=args.min_cd3_reads,
    )

    n_before_cd3 = len(adata)
    adata = apply_cd3_gate(adata, args.min_cd3_reads)
    if len(adata) != n_before_cd3:
        print(f"CD3 gate (>={args.min_cd3_reads} reads): {len(adata)}/{n_before_cd3} cells pass")

    warnings = validate_phenotype_vs_expected(adata)
    for w in warnings:
        print(f"WARNING: {w}")

    summary = get_phenotype_summary(adata)
    print("\nPhenotype Summary:")
    print(summary.to_string())

    adata.write_h5ad(args.output)
    print(f"\nSaved phenotyped data to {args.output}")

    if args.plot_phenotype:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_phenotype(adata, output_dir)
        print(f"Phenotype plots saved to {output_dir}")


def cmd_clonotype(args):
    """Aggregate cells into clonotypes."""
    import anndata as ad

    from ..clonotype import aggregate_clonotypes, export_clonotypes_airr, get_clonotype_summary
    from ..config import AttributionConfig
    from ..plots import plot_clonotypes

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

    summary = get_clonotype_summary(clonotypes)
    print("\nClonotype Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    clonotypes.to_csv(args.output, index=False)
    print(f"\nSaved {len(clonotypes)} clonotypes to {args.output}")

    if args.airr:
        export_clonotypes_airr(clonotypes, args.airr)
        print(f"Saved AIRR format to {args.airr}")

    if args.plot_clonotypes:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_clonotypes(clonotypes, output_dir)
        print(f"Clonotype plots saved to {output_dir}")


def register_load(subparsers) -> None:
    """Register the ``load`` subcommand."""
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
    p_load.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default="auto",
        help="Organism species: 'human', 'mouse', or 'auto' (default: auto-detect)",
    )
    p_load.add_argument("--plot-qc", action="store_true", help="Generate QC plots")
    p_load.add_argument("--output-dir", help="Output directory for plots")
    p_load.add_argument("--verbose", action="store_true", help="Verbose output")
    p_load.set_defaults(func=cmd_load)


def register_load_sct(subparsers) -> None:
    """Register the ``load-sct`` subcommand."""
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


def register_phenotype(subparsers) -> None:
    """Register the ``phenotype`` subcommand."""
    p_pheno = subparsers.add_parser("phenotype", help="Classify T cells as CD4/CD8")
    p_pheno.add_argument("--input", "-i", required=True, help="Input h5ad from load")
    p_pheno.add_argument("--output", "-o", required=True, help="Output h5ad")
    p_pheno.add_argument(
        "--cd4-cd8-ratio", type=float, default=3.0, help="Ratio for confident calls (default: 3.0)"
    )
    p_pheno.add_argument(
        "--min-cd3-reads",
        type=int,
        default=3,
        help="Min CD3 reads (default: 3; set 0 to disable, #172)",
    )
    p_pheno.add_argument("--plot-phenotype", action="store_true", help="Generate phenotype plots")
    p_pheno.add_argument("--output-dir", help="Output directory for plots")
    p_pheno.add_argument("--verbose", action="store_true", help="Verbose output")
    p_pheno.set_defaults(func=cmd_phenotype)


def register_clonotype(subparsers) -> None:
    """Register the ``clonotype`` subcommand."""
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
