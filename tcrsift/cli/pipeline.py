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

"""Unified end-to-end pipeline and configuration CLI subcommands.

Commands: ``run``, ``generate-config``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..config import TCRsiftConfig, add_config_args, generate_example_config, load_config_with_args
from ..pipeline import run_pipeline
from ..validation import validate_run_args
from .common import LEADER_CHOICES, setup_logging


def cmd_run(args):
    """Run the complete TCRsift pipeline."""
    setup_logging(bool(getattr(args, "verbose", False)))
    config = load_config_with_args(args)
    validate_run_args(args, config=config)
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


def cmd_generate_config(args):
    """Generate an example configuration file (prints to stdout when -o is omitted)."""
    if args.output is None:
        import yaml

        print(
            yaml.dump(TCRsiftConfig().to_dict(), default_flow_style=False, sort_keys=False),
            end="",
        )
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


def register_run(subparsers) -> None:
    """Register the ``run`` subcommand."""
    p_run = subparsers.add_parser("run", help="Run complete pipeline")
    p_run.add_argument("--sample-sheet", "-s", required=True, help="Sample sheet (CSV or YAML)")
    p_run.add_argument("--output-dir", "-o", required=True, help="Output directory")
    p_run.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default=None,
        help="Organism species ('human', 'mouse', or 'auto' default)",
    )
    add_config_args(p_run)

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

    pheno_group = p_run.add_argument_group("Phenotype options")
    pheno_group.add_argument(
        "--cd4-cd8-ratio", type=float, help="Ratio for confident calls (default: 3.0)"
    )
    pheno_group.add_argument(
        "--min-cd3-reads", type=int, help="Min CD3 reads (default: 3; set 0 to disable, #172)"
    )

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
    filter_group.add_argument("--min-donors", type=int, help="Min distinct donors per clone")
    filter_group.add_argument(
        "--min-methods-per-donor",
        type=int,
        help="Min distinct enrichment methods within at least one donor",
    )
    filter_group.add_argument(
        "--min-cells-per-method",
        type=int,
        help="Min cells in at least one enrichment method",
    )
    filter_group.add_argument(
        "--min-frequency-per-method",
        type=float,
        help="Min frequency within a single enrichment method",
    )
    filter_group.add_argument(
        "--min-timepoints",
        type=int,
        help="Min distinct timepoints clone must appear in",
    )
    filter_group.add_argument(
        "--min-timepoints-per-donor",
        type=int,
        help="Min distinct timepoints within at least one donor",
    )
    filter_group.add_argument(
        "--min-apcs",
        type=int,
        help="Min distinct APCs clone must appear across",
    )
    filter_group.add_argument(
        "--min-apcs-per-donor",
        type=int,
        help="Min distinct APCs within at least one donor",
    )
    filter_group.add_argument(
        "--min-til-cells-per-donor",
        type=int,
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

    asm_group = p_run.add_argument_group(
        "Assembly options",
        "Options for assembling full-length TCR sequences. "
        "NOTE: --contigs-dir is REQUIRED when using --leaders-from-contigs or leader=from_contig.",
    )
    asm_group.add_argument(
        "--alpha-leader",
        choices=[*LEADER_CHOICES, "from_contig", "none"],
        help="Alpha chain leader (default: CD28). Use 'none' for no leader.",
    )
    asm_group.add_argument(
        "--beta-leader",
        choices=[*LEADER_CHOICES, "from_contig", "none"],
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
        choices=["germline", "keep", *LEADER_CHOICES],
        default=None,
        help="Fallback leader when a native from_contig leader is rejected (default: germline).",
    )
    asm_group.add_argument(
        "--force-alpha-leader",
        choices=["germline", *LEADER_CHOICES],
        default=None,
        help="Force this leader on alpha chains.",
    )
    asm_group.add_argument(
        "--force-beta-leader",
        choices=["germline", *LEADER_CHOICES],
        default=None,
        help="Force this leader on beta chains.",
    )
    asm_group.add_argument(
        "--secondary-alpha-leader",
        choices=["germline", *LEADER_CHOICES],
        default=None,
        help="Emit secondary construct twin using this leader on alpha.",
    )
    asm_group.add_argument(
        "--secondary-beta-leader",
        choices=["germline", *LEADER_CHOICES],
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
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Auto-detect and rescue 5'-truncated contigs and bare (v_gene, CDR3, j_gene) "
        "rows from cached IMGT V-REGION + J-REGION germline references "
        "(default: enabled; pass --no-germline-vdj-fallback for strict contig-only assembly).",
    )
    asm_group.add_argument(
        "--auto-download-germline",
        action="store_true",
        default=None,
        help="Automatically download missing IMGT V-REGION and J-REGION germline "
        "references into the cache when germline VDJ fallback is enabled.",
    )

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
        "--report-name",
        metavar="NAME",
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
        help="Similarity metric for the method × method overlap matrix (default: jaccard).",
    )
    out_group.add_argument("--verbose", action="store_true", help="Verbose output")
    out_group.add_argument(
        "--allow-canonical-fallback",
        action="store_true",
        help="Permit constructs assembled from canonical/reference constants when "
        "no contig verifies the donor's junction/allele (fail-closed otherwise) "
        "(#243).",
    )
    out_group.add_argument(
        "--no-assemble",
        action="store_true",
        help="Skip whole-repertoire assembly (select-then-assemble): write the "
        "analysis outputs (clonotypes, clone_sample_long, phenotyped h5ad) but "
        "build synthesis sequences only for the selected basket via `select` + "
        "`report selected`. ~15–40× less assembly; the fail-closed gate (#243) "
        "then applies to the basket, not non-selected clones (#246).",
    )
    p_run.set_defaults(func=cmd_run)


def register_generate_config(subparsers) -> None:
    """Register the ``generate-config`` subcommand."""
    p_config = subparsers.add_parser("generate-config", help="Generate example config file")
    p_config.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output YAML file (default: print to stdout; refuses to overwrite "
        "an existing file without --force, #195)",
    )
    p_config.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing output file",
    )
    p_config.set_defaults(func=cmd_generate_config)
