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

"""Full-length TCR construct assembly and allele-audit CLI subcommands.

Commands: ``assemble``, ``audit-alleles``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..validation import validate_assemble_args, validate_file_exists
from .common import LEADER_CHOICES, setup_logging


def cmd_assemble(args):
    """Assemble full-length TCR sequences."""
    import pandas as pd

    from ..assemble import (
        assemble_full_sequences,
        assemble_qc_report,
        enforce_contig_fidelity,
        export_fasta,
        synthesis_qc_report,
        validate_sequences,
    )
    from ..plots import plot_assembly

    setup_logging(args.verbose)

    validate_assemble_args(args)
    validate_file_exists(args.input, "input clonotypes CSV")

    alpha_leader = args.alpha_leader
    beta_leader = args.beta_leader

    if getattr(args, "no_leaders", False):
        alpha_leader = None
        beta_leader = None
    elif getattr(args, "leaders_from_contigs", False):
        alpha_leader = "from_contig"
        beta_leader = "from_contig"
    else:
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
        sample_sheet=getattr(args, "sample_sheet", None),
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
            1
            if getattr(args, "max_codon_repeats", None) is None
            else args.max_codon_repeats
        ),
        avoid_enzymes=tuple(getattr(args, "avoid_enzymes", None) or ()),
        species=getattr(args, "species", "auto") or "auto",
        constant_species=getattr(args, "constant_species", None),
        codon_species=getattr(args, "codon_species", None),
        germline_vdj_fallback=bool(getattr(args, "germline_vdj_fallback", True)),
        data_dir=getattr(args, "cache_dir", None),
        auto_download_germline=bool(getattr(args, "auto_download_germline", False)),
    )

    # Fail closed on unverified (canonical-fallback) constructs unless the
    # user explicitly opted in (#243/#244). Raises TCRsiftValidationError when
    # any construct's junction/allele couldn't be contig-verified.
    fidelity = enforce_contig_fidelity(
        assembled,
        allow_canonical_fallback=getattr(args, "allow_canonical_fallback", False),
        context="assemble",
    )
    if fidelity:
        print(f"WARNING: {fidelity}")

    warnings = validate_sequences(
        assembled, species=getattr(args, "species", "auto") or "auto"
    )
    for w in warnings:
        print(f"WARNING: {w}")
    print(assemble_qc_report(assembled))
    synth_text = synthesis_qc_report(assembled)
    if synth_text:
        print(synth_text)

    assembled.to_csv(args.output, index=False)
    print(f"\nSaved assembled sequences to {args.output}")

    if args.airr:
        from ..clonotype import export_clonotypes_airr

        export_clonotypes_airr(assembled, args.airr)
        print(f"Saved AIRR format to {args.airr}")

    if args.fasta:
        seq_col = "single_chain_aa" if args.single_chain else "full_beta_aa"
        export_fasta(assembled, args.fasta, sequence_col=seq_col)
        print(f"Saved FASTA to {args.fasta}")

    if args.plot_assembly:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_assembly(assembled, output_dir)
        print(f"Assembly plots saved to {output_dir}")


def cmd_audit_alleles(args):
    """Run cohort-level allele audit + novel-allele detection (#119, #120).

    Reads a ``full_sequences.csv`` produced by ``tcrsift assemble`` or
    ``tcrsift run`` and prints the audit report. When ``--output`` is
    given, also writes the long-form novel-allele candidates CSV.
    """
    import pandas as pd

    from ..assemble import allele_audit_report, detect_novel_alleles

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


def register_assemble(subparsers) -> None:
    """Register the ``assemble`` subcommand."""
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

    asm_leaders = p_asm.add_argument_group(
        "leader peptide options",
        "Leader peptides are added to the N-terminus. Use default leaders (CD8A/CD28) "
        "or extract native leaders from CellRanger contigs.",
    )
    asm_leaders.add_argument(
        "--alpha-leader",
        choices=[*LEADER_CHOICES, "from_contig", "none"],
        default="CD28",
        help="Alpha chain leader (default: CD28)",
    )
    asm_leaders.add_argument(
        "--beta-leader",
        choices=[*LEADER_CHOICES, "from_contig", "none"],
        default="CD8A",
        help="Beta chain leader (default: CD8A)",
    )
    asm_leaders.add_argument(
        "--no-leaders", action="store_true", help="Disable leaders on both chains"
    )
    asm_leaders.add_argument(
        "--leader-fallback",
        choices=["germline", "keep", *LEADER_CHOICES],
        default=None,
        metavar="TARGET",
        help="With --*-leader=from_contig: what a REJECTED leader (bad SP, or "
        "diverges from germline AND inconsistent across the gene's clones) "
        "switches to — 'germline' (default: the gene's germline V-leader), a "
        "curated key (CD8A/mCD8B/H2Kb/…), or 'keep' to never switch (just flag) (#263/#270). "
        "An SP-sound, consistent contig leader is always kept.",
    )
    asm_leaders.add_argument(
        "--force-leader",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
        default=None,
        help="Force this leader (germline/curated) on BOTH chains, ignoring the "
        "contig. Per-chain: --force-alpha-leader / --force-beta-leader (#270).",
    )
    asm_leaders.add_argument(
        "--force-alpha-leader",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
    )
    asm_leaders.add_argument(
        "--force-beta-leader",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
    )
    asm_leaders.add_argument(
        "--leader-secondary",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
        default=None,
        help="Also emit a SECONDARY construct (twin row, leader_variant="
        "secondary:<leader>) using this leader on BOTH chains, even for good "
        "calls. Per-chain: --secondary-alpha-leader / --secondary-beta-leader (#270).",
    )
    asm_leaders.add_argument(
        "--secondary-alpha-leader",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
    )
    asm_leaders.add_argument(
        "--secondary-beta-leader",
        metavar="TARGET",
        choices=["germline", *LEADER_CHOICES],
    )
    asm_leaders.add_argument(
        "--codon-optimization",
        choices=["focal", "full"],
        default=None,
        metavar="MODE",
        help="How full_*_nt_optimized is generated: 'focal' (default) keeps the "
        "donor's native leader+VDJ codons and makes only the minimal synonymous "
        "swaps to clear constraints; 'full' recodes the whole ORF (#293).",
    )
    asm_leaders.add_argument(
        "--max-codon-repeats",
        type=int,
        default=None,
        metavar="N",
        help="Focal recode: longest allowed run of identical consecutive codons "
        "(default 1 = no immediate repeat — breaks native poly-codon runs like a "
        "triple-Gly poly-G). 0 disables.",
    )
    asm_leaders.add_argument(
        "--avoid-enzymes",
        action="extend",
        nargs="+",
        default=None,
        metavar="ENZYME",
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
        "(default: parent). 'grandparent' for per_sample_outs/{sample}/vdj_t/; "
        "'sheet' requires --sample-sheet.",
    )
    asm_leaders.add_argument(
        "--sample-sheet",
        metavar="FILE",
        default=None,
        help="Sample sheet CSV/YAML mapping vdj_dir paths to sample names "
        "(required when --sample-name-from sheet).",
    )

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
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Auto-detect and rescue 5'-truncated contigs and bare (v_gene, CDR3, j_gene) "
        "rows from cached IMGT V-REGION + J-REGION germline references (`tcrsift data "
        "download --group germline` / `--group germline_mouse`). Intact donor contigs are "
        "preserved verbatim; imputed chains stamp {chain}_vdj_source and qc_warnings and "
        "set construct_contig_verified=False. Pass --no-germline-vdj-fallback for strict "
        "contig-only assembly (default: enabled).",
    )
    asm_seq.add_argument(
        "--auto-download-germline",
        action="store_true",
        help="Automatically download missing IMGT V-REGION and J-REGION germline "
        "references into the cache when germline VDJ fallback is enabled.",
    )
    asm_seq.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Reference-database cache directory (default: ~/.cache/tcrsift).",
    )

    asm_out = p_asm.add_argument_group("additional output options")
    asm_out.add_argument("--airr", metavar="PATH", help="Also output in AIRR format to this path")
    asm_out.add_argument("--fasta", metavar="PATH", help="Also output in FASTA format to this path")
    asm_out.add_argument("--plot-assembly", action="store_true", help="Generate assembly QC plots")
    asm_out.add_argument(
        "--output-dir", metavar="DIR", help="Output directory for plots (default: next to output)"
    )
    asm_out.add_argument(
        "--allow-canonical-fallback",
        action="store_true",
        help="Permit constructs assembled from canonical/reference constants when "
        "no contig verifies the donor's junction/allele (fail-closed otherwise).",
    )
    asm_out.add_argument("--verbose", action="store_true", help="Verbose output")
    p_asm.set_defaults(func=cmd_assemble)


def register_audit_alleles(subparsers) -> None:
    """Register the ``audit-alleles`` subcommand."""
    p_audit = subparsers.add_parser(
        "audit-alleles",
        help="Cohort-level allele audit + novel-allele detection",
    )
    p_audit.add_argument(
        "input",
        help="Path to full_sequences.csv from tcrsift assemble / run",
    )
    p_audit.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="Write novel-allele candidates CSV to this path "
        "(in addition to the text report on stdout)",
    )
    p_audit.add_argument(
        "--min-pct",
        type=float,
        default=0.05,
        help="Minimum fraction of in-cohort clones for a "
        "novel-allele candidate (default: 0.05 = 5%%)",
    )
    p_audit.add_argument(
        "--min-v-spread",
        type=int,
        default=3,
        help="Minimum distinct V genes the variant must appear in (default: 3)",
    )
    p_audit.add_argument(
        "--min-samples",
        type=int,
        default=2,
        help="Minimum distinct samples the variant must appear in (default: 2)",
    )
    p_audit.set_defaults(func=cmd_audit_alleles)
