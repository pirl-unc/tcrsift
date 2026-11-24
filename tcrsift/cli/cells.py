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

"""Single-cell GEX atlas CLI subcommands: ``cells qc|embed|annotate`` (#342)."""

from __future__ import annotations

import argparse

from .common import setup_logging

# Sentinel distinguishing "threshold not passed on the CLI" (let cell_qc_funnel
# use its own default) from an explicit value — including ``none`` (disable the
# bound). argparse never applies ``type=`` to a default, so this survives untouched.
_UNSET_CLI = object()

_CELLS_QC_THRESHOLDS = (
    ("min_genes", "min_genes", "Min genes per cell (default: 250; pass 'none' to disable)"),
    ("max_genes", "max_genes", "Max genes per cell (default: disabled; pass 'none' to disable)"),
    ("min_counts", "min_counts", "Min UMI counts (default: 500; pass 'none' to disable)"),
    ("max_counts", "max_counts", "Max UMI counts (default: disabled; pass 'none' to disable)"),
    ("min_mito", "min_mito_pct", "Min mito %% floor (default: disabled; pass 'none' to disable)"),
    ("max_mito", "max_mito_pct", "Max mito %% (default: 8.0; pass 'none' to disable)"),
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

    from ..qc import cell_qc_funnel, lineage_programs_for_species
    from ..species import resolve_species

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
        from ..config import LoadConfig

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

    from ..config import EmbedConfig
    from ..embed import embed_cells
    from ..species import resolve_species

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

    from ..annotate_cells import annotate_cells
    from ..signatures import PBMC_CULTURE_TYPES, cell_type_signatures_for_species
    from ..species import resolve_species

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


def register(subparsers) -> None:
    """Register the ``cells`` command group (``qc``, ``embed``, ``annotate``)."""
    p_cells = subparsers.add_parser(
        "cells", help="Single-cell GEX path: QC → embed → annotate (.h5ad in/out)"
    )
    p_cells_sub = p_cells.add_subparsers(dest="cells_command", required=True)

    # cells qc
    p_cells_qc = p_cells_sub.add_parser(
        "qc", help="Per-cell QC + doublet funnel (cell_qc_funnel)"
    )
    p_cells_qc.add_argument("--input", "-i", required=True, help="Input h5ad (raw counts)")
    p_cells_qc.add_argument("--output", "-o", required=True, help="Output h5ad (filtered)")
    for cli_name, _param, helptext in _CELLS_QC_THRESHOLDS:
        p_cells_qc.add_argument(
            f"--{cli_name.replace('_', '-')}",
            type=_cells_threshold,
            default=_UNSET_CLI,
            metavar="N|none",
            help=helptext,
        )
    p_cells_qc.add_argument(
        "--pbmc-preset",
        action="store_true",
        help="Use the opinionated PBMC LoadConfig() thresholds as the baseline",
    )
    p_cells_qc.add_argument(
        "--solid-tumor",
        action="store_true",
        help="Use SOLID_TUMOR_LINEAGE_PROGRAMS for the cross-lineage doublet gate",
    )
    p_cells_qc.add_argument("--waterfall-csv", help="Also write the QC waterfall to this CSV")
    p_cells_qc.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species: 'human', 'mouse', or 'auto' (default: auto-detect)",
    )
    p_cells_qc.add_argument("--verbose", action="store_true", help="Verbose output")
    p_cells_qc.set_defaults(func=cmd_cells_qc)

    # cells embed
    p_cells_embed = p_cells_sub.add_parser(
        "embed", help="Pearson-residual PCA → Harmony → UMAP → Leiden (embed_cells)"
    )
    p_cells_embed.add_argument("--input", "-i", required=True, help="Input h5ad (raw counts)")
    p_cells_embed.add_argument("--output", "-o", required=True, help="Output h5ad (embedded)")
    p_cells_embed.add_argument(
        "--batch-key",
        default=None,
        help="obs column to integrate over with Harmony (default: no integration)",
    )
    p_cells_embed.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Leiden resolution (default: 1.0)",
    )
    p_cells_embed.add_argument("--n-pcs", type=int, default=40, help="PCA components (default: 40)")
    p_cells_embed.add_argument(
        "--n-neighbors", type=int, default=15, help="kNN neighbors (default: 15)"
    )
    p_cells_embed.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Highly-variable genes for the clustering panel (default: 2000)",
    )
    p_cells_embed.add_argument("--seed", type=int, default=0, help="Random seed (default: 0)")
    p_cells_embed.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species: 'human', 'mouse', or 'auto' (default: auto-detect)",
    )
    p_cells_embed.add_argument("--verbose", action="store_true", help="Verbose output")
    p_cells_embed.set_defaults(func=cmd_cells_embed)

    # cells annotate
    p_cells_annotate = p_cells_sub.add_parser(
        "annotate", help="Per-cluster cell-type + T/B-state annotation (annotate_cells)"
    )
    p_cells_annotate.add_argument(
        "--input", "-i", required=True, help="Input h5ad (embedded, has obs['leiden'])"
    )
    p_cells_annotate.add_argument("--output", "-o", required=True, help="Output h5ad (annotated)")
    p_cells_annotate.add_argument(
        "--leiden-col", default="leiden", help="Cluster column in obs (default: leiden)"
    )
    p_cells_annotate.add_argument(
        "--solid-tumor",
        action="store_true",
        help="Use solid-tumor cell-type signatures instead of the PBMC default",
    )
    p_cells_annotate.add_argument(
        "--pbmc-culture",
        action="store_true",
        help="Restrict the argmax to PBMC_CULTURE_TYPES (Ficoll/culture)",
    )
    p_cells_annotate.add_argument(
        "--min-subtype-frac",
        type=float,
        default=0.01,
        help="Sub-clusters below this fraction of their type lose the gene suffix (default: 0.01)",
    )
    p_cells_annotate.add_argument(
        "--no-suffix",
        action="store_true",
        help="Skip the distinguishing-gene suffix (phenotype == phenotype_base)",
    )
    p_cells_annotate.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species: 'human', 'mouse', or 'auto' (default: auto-detect)",
    )
    p_cells_annotate.add_argument("--verbose", action="store_true", help="Verbose output")
    p_cells_annotate.set_defaults(func=cmd_cells_annotate)
