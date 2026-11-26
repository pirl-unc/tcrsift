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

"""Database annotation, GEX augmentation, naming, and Pgen/Ppost CLI subcommands.

Commands: ``annotate``, ``annotate-gex``, ``mnemonic``, ``log-pgen``,
``pgen fetch|train|annotate|status``.
"""

from __future__ import annotations

from pathlib import Path

from .common import setup_logging


def cmd_annotate(args):
    """Annotate clonotypes with public database matches."""
    import pandas as pd

    from ..annotate import annotate_clonotypes, get_annotation_summary
    from ..datacache import ensure_cached
    from ..plots import plot_annotations

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
        max_distance=getattr(args, "max_distance", None),
        exclude_viral=args.exclude_viral,
        flag_only=args.flag_only,
        species=getattr(args, "species", "auto") or "auto",
        host_species_filter=getattr(args, "host_species_filter", "all") or "all",
    )

    summary = get_annotation_summary(annotated)
    print("\nAnnotation Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    annotated.to_csv(args.output, index=False)
    print(f"\nSaved annotated clonotypes to {args.output}")

    if args.plot_annotations:
        output_dir = (
            Path(args.output_dir) if args.output_dir else Path(args.output).parent / "plots"
        )
        plot_annotations(annotated, output_dir)
        print(f"Annotation plots saved to {output_dir}")


def cmd_annotate_gex(args):
    """Annotate TCR data with gene expression from a 10x HDF5 file."""
    import pandas as pd

    from ..gex import (
        DEFAULT_GENE_LIST,
        aggregate_gex_by_clonotype,
        augment_with_gex,
        compute_cd4_cd8_counts,
        compute_signature_scores_per_clonotype,
    )
    from ..validation import validate_annotate_gex_args, validate_file_exists

    setup_logging(args.verbose)

    validate_annotate_gex_args(args)
    validate_file_exists(args.input, "input CSV")

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df):,} rows from {args.input}")

    gene_list = None
    if args.genes:
        gene_list = [g.strip() for g in args.genes.split(",")]
        print(f"Using custom gene list: {gene_list}")
    else:
        print(f"Using default T cell gene list ({len(DEFAULT_GENE_LIST)} genes)")

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
        augmented_df = df.copy()
    else:
        print(
            f"Warning: Barcode column '{args.barcode_col}' not found - skipping per-cell augmentation"
        )

    if args.aggregate:
        print(f"\nAggregating by {args.group_col}...")
        df = aggregate_gex_by_clonotype(
            df,
            group_col=args.group_col,
            gex_prefix=args.prefix,
            operations=["sum", "mean"],
            verbose=args.verbose,
        )

    if augmented_df is not None:
        sig_group_col = (
            args.group_col
            if args.group_col in augmented_df.columns
            else ("CDR3ab" if "CDR3ab" in augmented_df.columns else None)
        )
        if sig_group_col is None:
            print(
                f"Warning: '{args.group_col}' not in input table — "
                "skipping per-clonotype signature aggregation (pass --group-col to specify)."
            )
        else:
            from ..species import resolve_species

            resolved_gex_sp = gex_species
            if resolved_gex_sp == "auto":
                try:
                    import scanpy as sc

                    from ..species import detect_adata_species

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
                group_col=sig_group_col,
                gex_prefix=args.prefix,
                cd8_only=True,
                species=resolved_gex_sp,
                verbose=args.verbose,
            )
            sig_cols = [c for c in sig_df.columns if c.startswith("signature_")]
            if sig_cols and sig_group_col in df.columns:
                df = df.merge(
                    sig_df[[sig_group_col, *sig_cols]],
                    on=sig_group_col,
                    how="left",
                )

    if args.cd4_cd8_counts:
        if augmented_df is None:
            print("Warning: Cannot compute CD4/CD8 counts without GEX augmentation - skipping")
        else:
            print("\nComputing CD4/CD8 cell counts...")
            cd4_cd8_df = compute_cd4_cd8_counts(
                augmented_df,
                group_col=args.group_col,
                gex_prefix=args.prefix,
                verbose=args.verbose,
            )
            df = df.merge(
                cd4_cd8_df[[args.group_col, "CD4_only.count", "CD8_only.count"]],
                on=args.group_col,
                how="left",
            )

    df.to_csv(args.output, index=False)
    print(f"\nSaved {len(df):,} rows to {args.output}")


def cmd_mnemonic(args):
    """Generate memorable names for TCR sequences."""
    import pandas as pd

    from ..mnemonic import tcr_name

    setup_logging(args.verbose)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} entries from {args.input}")

    cdr3_col = args.cdr3_col
    if not cdr3_col:
        for candidate in ["CDR3_beta", "CDR3_alpha", "cdr3", "cdr3_beta", "cdr3_alpha"]:
            if candidate in df.columns:
                cdr3_col = candidate
                break

    if not cdr3_col or cdr3_col not in df.columns:
        print(f"ERROR: Could not find CDR3 column. Available columns: {list(df.columns)}")
        return 1

    df[args.name_col] = df[cdr3_col].apply(lambda x: tcr_name(x) if pd.notna(x) else None)

    df.to_csv(args.output, index=False)
    print(f"Saved with mnemonic names to {args.output}")


def cmd_log_pgen(args):
    """Score per-clone log Pgen with a data-driven background (seqprob)."""
    import pandas as pd

    from ..seqprob import score_log_pgen

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
    from ..pgen_models import (
        REPERTOIRE_SOURCES,
        fetch_healthy_pbmc,
        fetch_repertoire,
    )

    sp = getattr(args, "species", "human") or "human"
    if args.source == "healthy-pbmc" and not args.url and sp == "human":
        print("Fetching the 5 pooled 10x healthy-PBMC VDJ-T sets...")
        for chain, path in fetch_healthy_pbmc().items():
            print(f"  [{chain}] pooled healthy-PBMC reference → {path}")
        return

    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    if args.source and args.source in REPERTOIRE_SOURCES and not args.url:
        print(f"Source {args.source!r}: {REPERTOIRE_SOURCES[args.source]}")
        print(
            "Download a healthy/neutral (not antigen-sorted) repertoire CSV "
            "from there and pass it via --url, or use --source healthy-pbmc."
        )
    failures: list[str] = []
    for chain in chains:
        try:
            path = fetch_repertoire(
                chain,
                url=args.url,
                source=args.source,
                seq_col=args.seq_col,
                v_col=args.v_col,
                chain_prefix=("TRA" if chain == "alpha" else "TRB") if args.v_col else None,
                species=sp,
            )
            print(f"  [{chain}] cached → {path}")
        except (FileNotFoundError, ValueError, OSError, RuntimeError) as exc:
            failures.append(chain)
            print(f"  [{chain}] {exc}")
    if failures:
        raise SystemExit(1)


def cmd_pgen_train(args):
    """Train + cache a custom Ppost model; packaged Pgen is not trainable."""
    from ..pgen_models import train_model

    sp = getattr(args, "species", "human") or "human"
    chains = ["alpha", "beta"] if args.chain == "both" else [args.chain]
    roles = ["pgen", "ppost"] if args.role == "both" else [args.role]
    failures: list[str] = []
    for chain in chains:
        for role in roles:
            try:
                path = train_model(
                    chain,
                    backend=args.backend,
                    role=role,
                    url=args.url,
                    epochs=args.epochs,
                    device=args.device,
                    species=sp,
                )
                print(f"  trained {args.backend} {role} {chain} ({sp}) → {path}")
            except (FileNotFoundError, ImportError, ValueError, OSError, RuntimeError) as exc:
                failures.append(f"{role}/{chain}")
                print(f"  [{role} {chain}] skipped: {exc}")
    if failures:
        raise SystemExit(1)


def cmd_pgen_annotate(args):
    """Annotate a TCR CSV with pgen/ppost (α/β); auto-trains on first use."""
    import pandas as pd

    from ..annotate_tcrs import add_pgen_ppost

    sp = getattr(args, "species", "auto") or "auto"
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} clones from {args.input}")
    chains = ("alpha", "beta") if args.chain == "both" else (args.chain,)
    df = add_pgen_ppost(
        df,
        chains=chains,
        backend=args.backend,
        auto_train=not args.no_auto_train,
        species=sp,
    )
    df.to_csv(args.output, index=False)
    print(f"Wrote {args.output} (backend={args.backend})")


def cmd_pgen_status(args):
    """Show cached Pgen/Ppost training data + trained models."""
    from ..pgen_models import cache_dir, cached_model_file, cached_repertoire_file

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


def register_annotate(subparsers) -> None:
    """Register the ``annotate`` subcommand."""
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
        choices=["exact", "levenshtein", "tcrdist"],
        default="exact",
        help="CDR3/TCR distance mode: 'exact' (default), 'levenshtein' (edit distance <= 1), "
        "or 'tcrdist' (V-gene CDR1/2/2.5 + BLOSUM62 trimmed-CDR3 distance).",
    )
    p_annot.add_argument(
        "--max-distance",
        type=int,
        default=None,
        help="Maximum match distance (defaults: 1 for levenshtein; 48 paired / 12 b_only for tcrdist).",
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


def register_annotate_gex(subparsers) -> None:
    """Register the ``annotate-gex`` subcommand."""
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
    )

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


def register_mnemonic(subparsers) -> None:
    """Register the ``mnemonic`` subcommand."""
    p_mnem = subparsers.add_parser("mnemonic", help="Generate memorable TCR names")
    p_mnem.add_argument("--input", "-i", required=True, help="Input CSV")
    p_mnem.add_argument("--output", "-o", required=True, help="Output CSV")
    p_mnem.add_argument(
        "--cdr3-col", help="Column with CDR3 sequences (auto-detected if not specified)"
    )
    p_mnem.add_argument("--name-col", default="tcr_name", help="Output column name")
    p_mnem.add_argument("--verbose", action="store_true")
    p_mnem.set_defaults(func=cmd_mnemonic)


def register_pgen(subparsers) -> None:
    """Register ``log-pgen`` and the ``pgen`` command group."""
    p_lpg = subparsers.add_parser(
        "log-pgen",
        help="Score per-clone log Pgen with a data-driven background "
        "(k-mer default, no extra deps; or --backend tcrpeg)",
    )
    p_lpg.add_argument("input", help="Path to a clonotype CSV")
    p_lpg.add_argument(
        "-o",
        "--output",
        required=True,
        metavar="PATH",
        help="Write the log-Pgen-annotated CSV here",
    )
    p_lpg.add_argument(
        "--chain",
        choices=["alpha", "beta", "both"],
        default="beta",
        help="Which chain(s) to score (default: beta)",
    )
    p_lpg.add_argument(
        "--backend",
        choices=["kmer", "kmer_gene", "tcrpeg"],
        default="kmer",
        help="Probability backend (default: kmer, numpy-only)",
    )
    p_lpg.add_argument(
        "--cdr3-col",
        default=None,
        help="CDR3 AA column (default: CDR3_<chain>)",
    )
    p_lpg.add_argument(
        "--species",
        choices=["human", "mouse"],
        default="human",
        help="Organism species (default: human)",
    )
    p_lpg.set_defaults(func=cmd_log_pgen)

    p_pgen = subparsers.add_parser(
        "pgen",
        help="Annotate Pgen/Ppost or manage custom Ppost backgrounds "
        "(packaged k-mer models are the default)",
    )
    p_pgen_sub = p_pgen.add_subparsers(dest="pgen_command")

    pf = p_pgen_sub.add_parser(
        "fetch", help="Fetch a neutral observed repertoire for Ppost training"
    )
    pf.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pf.add_argument(
        "--url",
        "--repertoire",
        dest="url",
        default=None,
        help="CSV/TSV (optionally .gz) with a CDR3 column",
    )
    pf.add_argument(
        "--source",
        default=None,
        choices=["healthy-pbmc", "ots", "ireceptor"],
        help="'healthy-pbmc' auto-downloads the 5 pooled 10x sets; "
        "'ots'/'ireceptor' are informational pointers",
    )
    pf.add_argument("--seq-col", default="seq", help="CDR3 column (default seq)")
    pf.add_argument(
        "--v-col",
        default=None,
        help="V-gene column to filter chain (TRA*/TRB*)",
    )
    pf.add_argument(
        "--species",
        choices=["human", "mouse"],
        default="human",
        help="Organism species (default: human)",
    )
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
    pt.add_argument(
        "--url",
        "--repertoire",
        dest="url",
        default=None,
        help="observed repertoire for ppost",
    )
    pt.add_argument("--epochs", type=int, default=30, help="TCRpeg epochs")
    pt.add_argument("--device", default=None, help="torch device (auto: cuda/cpu)")
    pt.add_argument(
        "--species",
        choices=["human", "mouse"],
        default="human",
        help="Organism species (default: human)",
    )
    pt.set_defaults(func=cmd_pgen_train)

    pa = p_pgen_sub.add_parser(
        "annotate",
        help="Annotate a TCR CSV with pgen/ppost (α/β); auto-trains on first use",
    )
    pa.add_argument("input", help="Path to a clonotype CSV")
    pa.add_argument("-o", "--output", required=True, metavar="PATH")
    pa.add_argument("--backend", choices=["kmer", "kmer_gene", "tcrpeg"], default="kmer")
    pa.add_argument("--chain", choices=["alpha", "beta", "both"], default="both")
    pa.add_argument(
        "--no-auto-train",
        action="store_true",
        help="error instead of training a missing TCRpeg model",
    )
    pa.add_argument(
        "--species",
        choices=["auto", "human", "mouse"],
        default="auto",
        help="Organism species (default: auto)",
    )
    pa.set_defaults(func=cmd_pgen_annotate)

    ps = p_pgen_sub.add_parser("status", help="Show cached training data + models")
    ps.add_argument(
        "--species",
        choices=["human", "mouse"],
        default=None,
        help="Filter status to a single species",
    )
    ps.set_defaults(func=cmd_pgen_status)

    p_pgen.set_defaults(func=lambda _a: (p_pgen.print_help(), 1)[1])
