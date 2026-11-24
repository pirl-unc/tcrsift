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

TCRsift prioritizes TCR clonotypes from antigen-enrichment cultures and
single-cell VDJ and gene-expression data.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from ._cli_commands import (  # noqa: F401
    _CELLS_QC_THRESHOLDS,
    _UNSET_CLI,
    _cells_threshold,
    _format_size,
    _parse_per_lineage_min,
    cmd_annotate,
    cmd_annotate_gex,
    cmd_assemble,
    cmd_audit_alleles,
    cmd_cells_annotate,
    cmd_cells_embed,
    cmd_cells_qc,
    cmd_clonotype,
    cmd_cohort,
    cmd_data_clear,
    cmd_data_download,
    cmd_data_list,
    cmd_filter,
    cmd_generate_config,
    cmd_load,
    cmd_load_sct,
    cmd_log_pgen,
    cmd_match_til,
    cmd_mnemonic,
    cmd_pgen_annotate,
    cmd_pgen_fetch,
    cmd_pgen_status,
    cmd_pgen_train,
    cmd_phenotype,
    cmd_report_bundle,
    cmd_report_cohort,
    cmd_report_cohort_figures,
    cmd_report_polished,
    cmd_report_selected,
    cmd_report_selected_combine,
    cmd_report_sequence,
    cmd_run,
    cmd_select,
    cmd_til_clonotype,
    cmd_til_prioritize,
    cmd_til_select,
    cmd_unify,
    setup_logging,
)
from .config import add_config_args
from .pipeline import _common_contig_root, run_pipeline  # noqa: F401
from .version import __version__

# =============================================================================
# Main Parser
# =============================================================================

_PRIMARY_WORKFLOW_COMMANDS = ("prioritize", "run", "assemble", "data")
_EXPLORE_ANNOTATE_COMMANDS = ("cells", "annotate", "report")


class _TieredHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Group top-level ``tcrsift --help`` subcommands by workflow role."""

    def _format_action(self, action):
        if not (isinstance(action, argparse._SubParsersAction) and action.dest == "command"):
            return super()._format_action(action)

        help_by_cmd = {
            sub.dest.split(" (", 1)[0]: (sub.help or "")
            for sub in action._get_subactions()
        }
        shown: set[str] = set()
        sections = [
            ("Primary workflows", [c for c in _PRIMARY_WORKFLOW_COMMANDS if c in help_by_cmd]),
            ("Exploration & annotation", [c for c in _EXPLORE_ANNOTATE_COMMANDS if c in help_by_cmd]),
        ]
        shown.update(c for _, cmds in sections for c in cmds)
        plumbing = [c for c in help_by_cmd if c not in shown]
        if plumbing:
            sections.append(("Modular pipeline steps & legacy commands", plumbing))

        lines: list[str] = []
        for idx, (heading, cmds) in enumerate(sections):
            if idx > 0:
                lines.append("")
            lines.append(f"{heading}:")
            for cmd in cmds:
                desc = " ".join(help_by_cmd[cmd].split())
                lines.append(f"  {cmd:<18} {desc}")
        lines.append("")
        return "\n".join(lines)


def create_parser():
    """Create the argument parser."""
    from .signatures import AIM_GENES_HGNC, EXHAUSTION_GENES_HGNC, MARKER_PANEL_HGNC
    from .til_prioritize import add_cli_args as add_til_prioritize_args

    parser = argparse.ArgumentParser(
        prog="tcrsift",
        formatter_class=_TieredHelpFormatter,
        description=(
            "TCRsift: prioritize TCR clonotypes from fresh tissue/blood or "
            "antigen-enrichment cultures and assemble synthesis-ready constructs."
        ),
        epilog=(
            "Run 'tcrsift <command> --help' for command-specific flags.\n"
            "Quick start:\n"
            "  Ex vivo TIL / blood:  tcrsift prioritize samples.yaml -o out/ --context solid-tumor --assemble\n"
            "  Peptide / culture:    tcrsift run --sample-sheet samples.yaml -o out/\n"
            "  Reference cache:      tcrsift data download --db vdjdb --group germline"
        ),
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser._positionals.title = None

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

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

    # -------------------------------------------------------------------------
    # Cells command group: single-cell GEX path (qc → embed → annotate) (#342)
    # -------------------------------------------------------------------------
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
            help=f"{helptext} (default: cell_qc_funnel's; pass 'none' to disable)",
        )
    p_cells_qc.add_argument(
        "--pbmc-preset", action="store_true",
        help="Use the opinionated PBMC LoadConfig() thresholds as the baseline",
    )
    p_cells_qc.add_argument(
        "--solid-tumor", action="store_true",
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
        "--batch-key", default=None,
        help="obs column to integrate over with Harmony (default: no integration)",
    )
    p_cells_embed.add_argument(
        "--resolution", type=float, default=1.0,
        help="Leiden resolution (default: 1.0)",
    )
    p_cells_embed.add_argument("--n-pcs", type=int, default=40, help="PCA components (default: 40)")
    p_cells_embed.add_argument(
        "--n-neighbors", type=int, default=15, help="kNN neighbors (default: 15)"
    )
    p_cells_embed.add_argument(
        "--n-top-genes", type=int, default=2000,
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
        "--solid-tumor", action="store_true",
        help="Use solid-tumor cell-type signatures instead of the PBMC default",
    )
    p_cells_annotate.add_argument(
        "--pbmc-culture", action="store_true",
        help="Restrict the argmax to PBMC_CULTURE_TYPES (Ficoll/culture)",
    )
    p_cells_annotate.add_argument(
        "--min-subtype-frac", type=float, default=0.01,
        help="Sub-clusters below this fraction of their type lose the gene suffix (default: 0.01)",
    )
    p_cells_annotate.add_argument(
        "--no-suffix", action="store_true",
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
        "--tie-break",
        help="--prism: comma-separated tie-break keys for boundary ties "
        "(default: single_alpha,umi,cd8_purity,cdr3ab; use 'cdr3ab' for "
        "CDR3ab-only). Precedence: this flag > config selection.tie_break > "
        "default (#229).",
    )
    p_select.add_argument(
        "--rescue-target", type=int, default=0,
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
        "--prism-grid", action="store_true",
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

    # -------------------------------------------------------------------------
    # Report command (#123): unified reporting (cohort figure bundle)
    # -------------------------------------------------------------------------
    p_report = subparsers.add_parser("report", help="Unified reporting (figure bundles)")
    p_report.set_defaults(func=lambda _a: p_report.print_help())
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

    pcf = p_report_sub.add_parser(
        "cohort-figures",
        help="Cross-donor figure bundle → all-figures-cross-donor.pdf (one section "
        "per donor, donor name resolved from each run dir)",
    )
    pcf.add_argument(
        "run_dirs", nargs="+",
        help="Per-donor run output dirs (each with a plots/ subdir)",
    )
    pcf.add_argument(
        "--output", "-o", default="all-figures-cross-donor.pdf",
        help="Output PDF (default: all-figures-cross-donor.pdf)",
    )
    pcf.add_argument("--verbose", action="store_true", help="Verbose output")
    pcf.set_defaults(func=cmd_report_cohort_figures)

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
        "--contigs-dir", help="Per-donor contigs dir (overrides config, #196)",
    )
    ps.add_argument(
        "--sample-name-from", choices=["parent", "grandparent", "sheet"],
        help="Contig sample-name policy (overrides config, #196)",
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
    ps.add_argument("--output-dir", "-o", required=True, help="Output directory")
    ps.add_argument(
        "--no-cover", action="store_true",
        help="Skip the construct cover/legend page on the sequence PDF (#202)",
    )
    ps.add_argument("--verbose", action="store_true", help="Verbose output")
    ps.add_argument(
        "--allow-canonical-fallback", action="store_true",
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
    p_annot.add_argument(
        "--database-match",
        dest="database_match",
        choices=["strict_ab", "ab_with_partial", "b_only"],
        default=None,
        help="Explicit chain-matching strictness (overrides --match-by when set): "
        "'strict_ab' (paired α+β exact only), 'ab_with_partial' (paired α+β with β-only fallback), "
        "or 'b_only' (CDR3β only).",
    )
    p_annot.add_argument(
        "--match-mode",
        choices=["exact", "levenshtein"],
        default="exact",
        help="CDR3β distance mode: 'exact' (default) or 'levenshtein' (edit distance <= --max-distance).",
    )
    p_annot.add_argument(
        "--max-distance",
        type=int,
        default=1,
        help="Maximum Levenshtein edit distance for --match-mode levenshtein (default: 1).",
    )
    p_annot.add_argument("--exclude-viral", action="store_true", help="Remove viral clones")
    p_annot.add_argument("--flag-only", action="store_true", help="Just flag, don't remove")
    p_annot.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default="auto",
        help="Query TCR organism ('human', 'mouse', or 'auto' default)",
    )
    p_annot.add_argument(
        "--host-species-filter",
        choices=["all", "same"],
        default="all",
        help="Filter DB entries by TCR donor host species: 'all' (default, tags "
        "heterologous hosts with *_cross) or 'same' (restrict to same-species hosts).",
    )
    p_annot.add_argument(
        "--auto-download",
        action="store_true",
        help="Automatically download missing public TCR databases into the cache",
    )
    p_annot.add_argument(
        "--plot-annotations", action="store_true", help="Generate annotation plots"
    )
    p_annot.add_argument("--output-dir", help="Output directory for plots")
    p_annot.add_argument("--verbose", action="store_true", help="Verbose output")
    p_annot.set_defaults(func=cmd_annotate)

    # -------------------------------------------------------------------------
    # Data command (reference-database cache management)
    # -------------------------------------------------------------------------
    from .datacache import DATABASES, GROUPS

    p_data = subparsers.add_parser(
        "data",
        help="Manage cached reference databases (VDJdb, IEDB, CEDAR, IMGT V/J-REGION)",
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
        "--group",
        choices=sorted(GROUPS),
        action="extend",
        nargs="+",
        help="Clear a whole workflow group (germline / germline_mouse / annotation / all).",
    )
    p_data_clear.add_argument(
        "--cache-dir", help="Override the default cache directory"
    )
    p_data_clear.set_defaults(func=cmd_data_clear)

    p_data_download = p_data_sub.add_parser(
        "download",
        help="Download reference databases (VDJdb, IEDB, CEDAR, IMGT human/mouse germline).",
    )
    p_data_download.add_argument(
        "--db",
        choices=sorted(DATABASES),
        action="extend",
        nargs="+",
        help=(
            "Database(s) to download (default: all registered DBs). "
            "Repeatable and/or space-separated: "
            "`--db vdjdb iedb` and `--db vdjdb --db iedb` both work."
        ),
    )
    p_data_download.add_argument(
        "--group",
        choices=sorted(GROUPS),
        action="extend",
        nargs="+",
        help=(
            "Download a whole workflow group: `germline` (human IMGT V-REGION + J-REGION), "
            "`germline_mouse` (mouse IMGT V-REGION + J-REGION + leaders), "
            "`annotation` (VDJdb/IEDB/CEDAR), or `all`. Combines with --db."
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
    # Sample-sheet clone prioritization commands share one implementation.
    # -------------------------------------------------------------------------
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

    # Legacy ordered/timepoint TIL selection command
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
        default=",".join(EXHAUSTION_GENES_HGNC),
        help="Comma-separated genes for the exhaustion score branch (#303).",
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
    p_til_select.add_argument("--vdjdb", type=Path, default=None, help="Path to VDJdb TSV (overrides cache)")
    p_til_select.add_argument("--iedb", type=Path, default=None, help="Path to IEDB CSV/TSV (overrides cache)")
    p_til_select.add_argument("--cedar", type=Path, default=None, help="Path to CEDAR CSV/TSV (overrides cache)")
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
    # Dominant-fraction antigen-specificity selection (#328 / #316).
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
    _LEADER_CHOICES = ["CD8A", "CD28", "mCD8B", "mCD28", "H2Kb", "IgK", "TRAC", "TRBC"]
    asm_leaders.add_argument(
        "--alpha-leader",
        choices=[*_LEADER_CHOICES, "from_contig", "none"],
        default="CD28",
        help="Alpha chain leader (default: CD28)",
    )
    asm_leaders.add_argument(
        "--beta-leader",
        choices=[*_LEADER_CHOICES, "from_contig", "none"],
        default="CD8A",
        help="Beta chain leader (default: CD8A)",
    )
    asm_leaders.add_argument(
        "--no-leaders", action="store_true", help="Disable leaders on both chains"
    )
    asm_leaders.add_argument(
        "--leader-fallback",
        choices=["germline", "keep", *_LEADER_CHOICES],
        default=None,
        metavar="TARGET",
        help="With --*-leader=from_contig: what a REJECTED leader (bad SP, or "
        "diverges from germline AND inconsistent across the gene's clones) "
        "switches to — 'germline' (default: the gene's germline V-leader), a "
        "curated key (CD8A/mCD8B/H2Kb/…), or 'keep' to never switch (just flag) (#263/#270). "
        "An SP-sound, consistent contig leader is always kept.",
    )
    asm_leaders.add_argument(
        "--force-leader", metavar="TARGET",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Force this leader (germline/curated) on BOTH chains, ignoring the "
        "contig. Per-chain: --force-alpha-leader / --force-beta-leader (#270).",
    )
    asm_leaders.add_argument("--force-alpha-leader", metavar="TARGET",
                             choices=["germline", *_LEADER_CHOICES])
    asm_leaders.add_argument("--force-beta-leader", metavar="TARGET",
                             choices=["germline", *_LEADER_CHOICES])
    asm_leaders.add_argument(
        "--leader-secondary", metavar="TARGET",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Also emit a SECONDARY construct (twin row, leader_variant="
        "secondary:<leader>) using this leader on BOTH chains, even for good "
        "calls. Per-chain: --secondary-alpha-leader / --secondary-beta-leader (#270).",
    )
    asm_leaders.add_argument("--secondary-alpha-leader", metavar="TARGET",
                             choices=["germline", *_LEADER_CHOICES])
    asm_leaders.add_argument("--secondary-beta-leader", metavar="TARGET",
                             choices=["germline", *_LEADER_CHOICES])
    asm_leaders.add_argument(
        "--codon-optimization", choices=["focal", "full"], default=None, metavar="MODE",
        help="How full_*_nt_optimized is generated: 'focal' (default) keeps the "
        "donor's native leader+VDJ codons and makes only the minimal synonymous "
        "swaps to clear constraints; 'full' recodes the whole ORF (#293).",
    )
    asm_leaders.add_argument(
        "--max-codon-repeats", type=int, default=None, metavar="N",
        help="Focal recode: longest allowed run of identical consecutive codons "
        "(default 1 = no immediate repeat — breaks native poly-codon runs like a "
        "triple-Gly poly-G). 0 disables.",
    )
    asm_leaders.add_argument(
        "--avoid-enzymes", action="extend", nargs="+", default=None, metavar="ENZYME",
        help="Focal recode: restriction enzyme name(s) (e.g. BsaI EcoRI) or raw "
        "ACGT site(s) to synonymously recode out of the construct (default none).",
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
        choices=["canonical", "ensembl", "from-data"],
        default="ensembl",
        help="Constant source: packaged canonical (canonical or legacy alias "
             "ensembl, the default) or explicit from-data columns",
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
    asm_seq.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default="auto",
        help="VDJ organism species ('human', 'mouse', or 'auto' default).",
    )
    asm_seq.add_argument(
        "--constant-species",
        choices=["human", "mouse"],
        default=None,
        help="Constant-region species override (e.g. '--constant-species mouse' "
        "on human VDJ for murinized constant domains avoiding endogenous mispairing).",
    )
    asm_seq.add_argument(
        "--codon-species",
        choices=["human", "mouse"],
        default=None,
        help="Host codon-optimization table override ('human' or 'mouse'; "
        "defaults to constant-species).",
    )
    asm_seq.add_argument(
        "--germline-vdj-fallback",
        action="store_true",
        help="Rescue 5'-truncated contigs and bare (v_gene, CDR3, j_gene) rows "
        "from cached IMGT V-REGION + J-REGION germline references (`tcrsift data "
        "download --group germline` / `--group germline_mouse`). Stamps "
        "{chain}_vdj_source and qc_warnings, and sets construct_contig_verified=False.",
    )
    asm_seq.add_argument(
        "--auto-download-germline",
        action="store_true",
        help="With --germline-vdj-fallback, automatically download missing IMGT "
        "V-REGION and J-REGION germline references into the cache.",
    )
    asm_seq.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Reference-database cache directory (default: ~/.cache/tcrsift).",
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
    p_asm.add_argument(
        "--allow-canonical-fallback", action="store_true",
        help="Permit constructs assembled from canonical/reference constants when "
        "no contig verifies the donor's junction/allele (fail-closed otherwise) "
        "(#243).",
    )
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
        "--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer",
        help="Probability backend (default: kmer, numpy-only)",
    )
    p_lpg.add_argument(
        "--cdr3-col", default=None,
        help="CDR3 AA column (default: CDR3_<chain>)",
    )
    p_lpg.add_argument(
        "--species", choices=["human", "mouse"], default="human",
        help="Organism species (default: human)",
    )
    p_lpg.set_defaults(func=cmd_log_pgen)

    # -------------------------------------------------------------------------
    # pgen group: manage Pgen/Ppost background models (fetch / train / annotate)
    # -------------------------------------------------------------------------
    p_pgen = subparsers.add_parser(
        "pgen",
        help="Annotate Pgen/Ppost or manage custom Ppost backgrounds "
             "(packaged k-mer models are the default)",
    )
    p_pgen_sub = p_pgen.add_subparsers(dest="pgen_command")

    pf = p_pgen_sub.add_parser("fetch", help="Fetch a neutral observed "
                               "repertoire for Ppost training")
    pf.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pf.add_argument("--url", "--repertoire", dest="url", default=None,
                    help="CSV/TSV (optionally .gz) with a CDR3 column")
    pf.add_argument("--source", default=None,
                    choices=["healthy-pbmc", "ots", "ireceptor"],
                    help="'healthy-pbmc' auto-downloads the 5 pooled 10x sets; "
                         "'ots'/'ireceptor' are informational pointers")
    pf.add_argument("--seq-col", default="seq", help="CDR3 column (default seq)")
    pf.add_argument("--v-col", default=None,
                    help="V-gene column to filter chain (TRA*/TRB*)")
    pf.add_argument("--species", choices=["human", "mouse"], default="human",
                    help="Organism species (default: human)")
    pf.set_defaults(func=cmd_pgen_fetch)

    pt = p_pgen_sub.add_parser(
        "train",
        help="Train + cache a custom Ppost model (Pgen models are packaged)",
    )
    pt.add_argument("--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer")
    pt.add_argument(
        "--role",
        choices=["pgen", "ppost", "both"],
        default="ppost",
        help="Role to request (default: ppost); packaged Pgen cannot be retrained",
    )
    pt.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pt.add_argument("--url", "--repertoire", dest="url", default=None,
                    help="observed repertoire for ppost")
    pt.add_argument("--epochs", type=int, default=30, help="TCRpeg epochs")
    pt.add_argument("--device", default=None, help="torch device (auto: cuda/cpu)")
    pt.add_argument("--species", choices=["human", "mouse"], default="human",
                    help="Organism species (default: human)")
    pt.set_defaults(func=cmd_pgen_train)

    pa = p_pgen_sub.add_parser("annotate", help="Annotate a TCR CSV with "
                               "pgen/ppost (α/β); auto-trains on first use")
    pa.add_argument("input", help="Path to a clonotype CSV")
    pa.add_argument("-o", "--output", required=True, metavar="PATH")
    pa.add_argument("--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer")
    pa.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pa.add_argument("--no-auto-train", action="store_true",
                    help="error instead of training a missing TCRpeg model")
    pa.add_argument("--species", choices=["auto", "human", "mouse"], default="auto",
                    help="Organism species (default: auto)")
    pa.set_defaults(func=cmd_pgen_annotate)

    ps = p_pgen_sub.add_parser("status", help="Show cached training data + models")
    ps.add_argument("--species", choices=["human", "mouse"], default=None,
                    help="Filter status to a single species")
    ps.set_defaults(func=cmd_pgen_status)

    p_pgen.set_defaults(func=lambda _a: p_pgen.print_help())

    # -------------------------------------------------------------------------
    # Run command (unified pipeline)
    # -------------------------------------------------------------------------
    p_run = subparsers.add_parser("run", help="Run complete pipeline")
    p_run.add_argument("--sample-sheet", "-s", required=True, help="Sample sheet (CSV or YAML)")
    p_run.add_argument("--output-dir", "-o", required=True, help="Output directory")
    p_run.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species ('human', 'mouse', or 'auto' default)",
    )
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
    clone_group.add_argument(
        "--no-pgen-ppost",
        action="store_true",
        help="Skip the default k-mer Pgen/Ppost annotation of clonotypes "
        "(otherwise ppost_alpha/ppost_beta are emitted for PRISM)",
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
        "--method",
        choices=["threshold", "logistic", "stimulation"],
        help="Filtering method (default: threshold; 'stimulation' uses paired control expansion or AIM delta)",
    )
    filter_group.add_argument(
        "--stim-fdr-threshold",
        type=float,
        default=None,
        help="Max Fisher exact FDR for --method stimulation (default: 0.1)",
    )
    filter_group.add_argument(
        "--stim-log2-fc-threshold",
        type=float,
        default=None,
        help="Min log2 fold-change over control for --method stimulation (default: 1.0)",
    )
    filter_group.add_argument(
        "--stim-delta-threshold",
        type=float,
        default=None,
        help="Min within-clone AIM/activation delta over control for --method stimulation (default: 0.5)",
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
    annot_group.add_argument(
        "--iedb-epitope",
        dest="iedb_epitope_path",
        default=None,
        help="Path to IEDB epitope_full_v3.csv (overrides cache)",
    )
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
        "--database-match",
        choices=["strict_ab", "ab_with_partial", "b_only"],
        default=None,
        help="Chain-matching policy: 'ab_with_partial' (default), 'strict_ab', or 'b_only'.",
    )
    annot_group.add_argument(
        "--match-mode",
        choices=["exact", "levenshtein"],
        default=None,
        help="CDR3 comparison mode: 'exact' (default) or 'levenshtein'.",
    )
    annot_group.add_argument(
        "--max-distance",
        type=int,
        default=None,
        help="Maximum Levenshtein edit distance when --match-mode=levenshtein (default: 1).",
    )
    annot_group.add_argument(
        "--exclude-viral", action="store_true", default=None, help="Remove viral clones"
    )
    annot_group.add_argument(
        "--flag-only", action="store_true", default=None, help="Flag but don't remove viral"
    )
    annot_group.add_argument(
        "--host-species-filter",
        choices=["all", "same"],
        default=None,
        help="Filter DB entries by TCR donor host species: 'all' (default) or 'same'.",
    )
    annot_group.add_argument(
        "--auto-download",
        action="store_true",
        default=None,
        help="Automatically download missing public TCR databases into the cache.",
    )
    annot_group.add_argument(
        "--call-multimers",
        action="store_true",
        default=None,
        help="Extract 10x Feature Barcode pMHC multimers and summarize per-clone specificity.",
    )
    annot_group.add_argument(
        "--multimer-min-umi",
        type=float,
        default=None,
        help="Min top-multimer UMI count for --call-multimers (default: 10.0).",
    )
    annot_group.add_argument(
        "--multimer-min-nc-ratio",
        type=float,
        default=None,
        help="Min fold-change over negative-control multimers (default: 5.0).",
    )
    annot_group.add_argument(
        "--multimer-min-concordance",
        type=float,
        default=None,
        help="Min clonal concordance for multimer_confident (default: 0.5).",
    )
    annot_group.add_argument(
        "--multimer-min-specific-cells",
        type=int,
        default=None,
        help="Min cells calling modal multimer for multimer_confident (default: 2).",
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
        choices=[*_LEADER_CHOICES, "from_contig", "none"],
        help="Alpha chain leader (default: CD28). Use 'none' for no leader.",
    )
    asm_group.add_argument(
        "--beta-leader",
        choices=[*_LEADER_CHOICES, "from_contig", "none"],
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
        choices=["canonical", "ensembl", "from-data"],
        help="Constant source: packaged canonical (canonical or legacy alias "
             "ensembl, the default) or explicit from-data columns",
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
    asm_group.add_argument(
        "--constant-species",
        choices=["human", "mouse"],
        default=None,
        help="Constant-region species override ('human' or 'mouse', e.g. murinized constants).",
    )
    asm_group.add_argument(
        "--codon-species",
        choices=["human", "mouse"],
        default=None,
        help="Codon-optimization table override ('human' or 'mouse').",
    )
    asm_group.add_argument(
        "--leader-fallback",
        choices=["germline", "keep", *_LEADER_CHOICES],
        default=None,
        help="Fallback leader when a native from_contig leader is rejected (default: germline).",
    )
    asm_group.add_argument(
        "--force-alpha-leader",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Force this leader on alpha chains.",
    )
    asm_group.add_argument(
        "--force-beta-leader",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Force this leader on beta chains.",
    )
    asm_group.add_argument(
        "--secondary-alpha-leader",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Emit secondary construct twin using this leader on alpha.",
    )
    asm_group.add_argument(
        "--secondary-beta-leader",
        choices=["germline", *_LEADER_CHOICES],
        default=None,
        help="Emit secondary construct twin using this leader on beta.",
    )
    asm_group.add_argument(
        "--codon-optimization",
        choices=["focal", "full"],
        default=None,
        help="Codon optimization mode: 'focal' (default) or 'full'.",
    )
    asm_group.add_argument(
        "--max-codon-repeats",
        type=int,
        default=None,
        help="Longest allowed run of identical consecutive codons (default: 1; 0 disables).",
    )
    asm_group.add_argument(
        "--avoid-enzymes",
        action="extend",
        nargs="+",
        default=None,
        metavar="ENZYME",
        help="Restriction enzyme name(s) or ACGT site(s) to recode out of constructs.",
    )
    asm_group.add_argument(
        "--germline-vdj-fallback",
        action="store_true",
        default=None,
        help="Rescue 5'-truncated contigs and bare (v_gene, CDR3, j_gene) rows "
        "from cached IMGT V-REGION + J-REGION germline references.",
    )
    asm_group.add_argument(
        "--auto-download-germline",
        action="store_true",
        default=None,
        help="With --germline-vdj-fallback, automatically download missing IMGT "
        "V-REGION and J-REGION germline references into the cache.",
    )

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
        choices=["png", "pdf", "svg", "both"],
        help="Figure format (default: png); pdf/svg emit that vector file "
        "alongside a PNG, 'both' emits PDF + SVG + PNG, for publication "
        "figures (#169, #258)",
    )
    out_group.add_argument(
        "--emit-clone-sample-long",
        choices=["auto", "always", "never"],
        help="Emit data/clone_sample_long.csv: 'auto' (default) writes it "
        "when the sheet has >=2 samples; 'always' / 'never' force.",
    )
    out_group.add_argument(
        "--report-name", metavar="NAME",
        help="Donor/run label for the report cover (all-figures.pdf). Overrides "
        "the auto-detected name (sample-sheet donor field, else output dir name).",
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
    out_group.add_argument(
        "--allow-canonical-fallback", action="store_true",
        help="Permit constructs assembled from canonical/reference constants when "
        "no contig verifies the donor's junction/allele (fail-closed otherwise) "
        "(#243).",
    )
    out_group.add_argument(
        "--no-assemble", action="store_true",
        help="Skip whole-repertoire assembly (select-then-assemble): write the "
        "analysis outputs (clonotypes, clone_sample_long, phenotyped h5ad) but "
        "build synthesis sequences only for the selected basket via `select` + "
        "`report selected`. ~15–40× less assembly; the fail-closed gate (#243) "
        "then applies to the basket, not non-selected clones (#246).",
    )

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
    gex_opts.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default="auto",
        help="Organism species for signature scoring ('human', 'mouse', or 'auto')",
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
    if isinstance(rc, int) and not isinstance(rc, bool) and rc != 0:
        sys.exit(rc)


if __name__ == "__main__":
    main()
