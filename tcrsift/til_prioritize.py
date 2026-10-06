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

"""Prioritize clonotypes from paired 10x VDJ + GEX samples.

This workflow produces an auditable shortlist for experimental validation. It combines:

- observed expansion (cells and within-sample frequency),
- context-appropriate expression signatures with round-robin selection,
- exact public-database annotation,
- sequence-background publicness and alpha-chain promiscuity flags.

It writes every scored clone before writing the filtered shortlist, so no
exclusion is hidden.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .prioritize import (
    CONTEXT_SIGNATURES,
    QC_DEFAULTS,
    SIGNATURE_LINEAGE,
    filter_cells,
    gene_filter_mask,
    parse_expression_limits,
    resolve_signatures,
    select_round_robin,
)

logger = logging.getLogger(__name__)

MART1_PATTERN = re.compile(
    r"MART[- ]?1|MELAN[- ]?A|MLANA|E(?:AA|LA)GIGILTV|AAGIGILTV",
    flags=re.IGNORECASE,
)


def _score_within_samples(adata, name: str) -> np.ndarray:
    """Score one signature separately within each sample and compatible lineage."""
    from .signature_methods import score_by_name

    scores = np.full(adata.n_obs, np.nan, dtype=float)
    for (_, lineage), positions in adata.obs.groupby(["sample", "lineage"], observed=True).indices.items():
        if SIGNATURE_LINEAGE.get(name, lineage) != lineage:
            continue
        sample_adata = adata[positions].copy()
        values = score_by_name(
            sample_adata,
            name,
            log1p=False,  # the caller has already made log1p(CP10K)
            on_missing="error",
        )
        scores[positions] = values.to_numpy(dtype=float)
    return scores


def _analysis_units(adata) -> pd.Series:
    """Return a complete per-cell patient key, or one cohort-wide unit."""
    if "patient_id" not in adata.obs:
        return pd.Series("cohort", index=adata.obs_names, dtype="string")
    units = adata.obs["patient_id"].astype("string").str.strip()
    missing = units.isna() | units.eq("")
    if missing.any():
        missing_samples = sorted(adata.obs.loc[missing, "sample"].dropna().astype(str).unique())
        raise ValueError(
            "patient_id must be populated for every sample when the column is "
            f"used; missing for {missing_samples}"
        )
    return units


def _aggregate_within_patients(
    adata,
    aggregate_clonotypes,
    build_clone_sample_long,
    min_umi=2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate shared samples together without merging clones across patients."""
    units = _analysis_units(adata)
    clonotype_tables = []
    clone_sample_tables = []
    for patient in pd.unique(units):
        patient_adata = adata[units.eq(patient).to_numpy()].copy()
        patient_clones = aggregate_clonotypes(patient_adata, min_umi=min_umi)
        patient_clones.insert(0, "donor", patient)
        clonotype_tables.append(patient_clones)

        lineages = patient_adata.obs.get("lineage", pd.Series("both", index=patient_adata.obs_names))
        for lineage in pd.unique(lineages):
            subset = patient_adata[lineages.eq(lineage).to_numpy()].copy()
            patient_long = build_clone_sample_long(subset, min_umi=min_umi)
            patient_long["donor"] = patient
            patient_long["lineage"] = lineage
            clone_sample_tables.append(patient_long)

    return (
        pd.concat(clonotype_tables, ignore_index=True),
        pd.concat(clone_sample_tables, ignore_index=True),
    )


def _ensure_clone_key(adata) -> None:
    if "CDR3ab" in adata.obs:
        return
    required = {"CDR3_alpha", "CDR3_beta"}
    if not required.issubset(adata.obs):
        raise ValueError("Loaded cells have no CDR3ab or paired CDR3_alpha/CDR3_beta columns")
    alpha = adata.obs["CDR3_alpha"].fillna("").astype(str).str.strip()
    beta = adata.obs["CDR3_beta"].fillna("").astype(str).str.strip()
    adata.obs["CDR3ab"] = (alpha + "_" + beta).where(
        alpha.ne("") & beta.ne(""),
        other=pd.NA,
    )


def _known_mart1_mask(df: pd.DataFrame) -> pd.Series:
    annotation_cols = [
        col
        for col in (
            "db_epitope",
            "db_protein",
            "db_protein_canonical",
            "db_species",
            "db_category_detail",
        )
        if col in df
    ]
    if not annotation_cols:
        return pd.Series(False, index=df.index)
    text = df[annotation_cols].fillna("").astype(str).agg(" ".join, axis=1)
    return text.str.contains(MART1_PATTERN, na=False)


def _trav12_2_mask(df: pd.DataFrame) -> pd.Series:
    if "alpha_v_gene" not in df:
        return pd.Series(False, index=df.index)
    gene = df["alpha_v_gene"].fillna("").astype(str).str.upper().str.split("*", n=1).str[0]
    return gene.eq("TRAV12-2")


def _join_flags(df: pd.DataFrame, columns: list[tuple[str, str]]) -> pd.Series:
    values = []
    for _, row in df.iterrows():
        values.append(";".join(label for col, label in columns if bool(row[col])))
    return pd.Series(values, index=df.index, dtype="string")


def run_til_prioritize(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score and prioritize clones from a standard CellRanger sample sheet."""
    names = resolve_signatures(args.context, args.tcell_type, args.signatures, args.exclude_signatures)
    if not 0 <= args.signature_quantile <= 1:
        raise ValueError("--signature-quantile must be in [0, 1]")
    if args.exclude_public_quantile is not None and not (0 < args.exclude_public_quantile <= 1):
        raise ValueError("--exclude-public-quantile must be in (0, 1]")
    if not 0 <= args.min_frequency <= 1:
        raise ValueError("--min-frequency must be in [0, 1]")
    if args.min_cells < 1:
        raise ValueError("--min-cells must be >= 1")
    if not 1 <= args.min_signature_support <= len(names):
        raise ValueError(f"--min-signature-support must be between 1 and {len(names)}")
    for key in ("max_clones", "min_vdj_umis", "min_vdj_reads", "min_cd3",
                "min_alpha_cdr3_length", "min_beta_cdr3_length", *QC_DEFAULTS):
        if not np.isfinite(getattr(args, key)) or getattr(args, key) < 0:
            raise ValueError(f"--{key.replace('_', '-')} must be finite and >= 0")
    for low, high in (("min_genes", "max_genes"), ("min_counts", "max_counts"),
                      ("min_mito_pct", "max_mito_pct")):
        if getattr(args, low) > getattr(args, high):
            raise ValueError(f"--{low.replace('_', '-')} must not exceed --{high.replace('_', '-')}")
    if args.max_mito_pct > 100:
        raise ValueError("--max-mito-pct must not exceed 100")
    for option in ("min_expression", "max_expression"):
        try:
            parse_expression_limits(getattr(args, option))
        except ValueError as exc:
            raise ValueError(f"--{option.replace('_', '-')}: {exc}") from exc

    import scanpy as sc

    from .annotate import annotate_clonotypes
    from .annotate_tcrs import add_paired_ppost, add_pairing_promiscuity, add_pgen_ppost
    from .clonotype import aggregate_clonotypes, build_clone_sample_long
    from .loader import load_samples
    from .phenotype import phenotype_cells
    from .sample_sheet import load_sample_sheet

    args.output_dir.mkdir(parents=True, exist_ok=True)

    sheet = load_sample_sheet(args.sample_sheet)
    sample_names = [s.sample for s in sheet]
    if len(set(sample_names)) != len(sample_names):
        raise ValueError("Each sample must have a unique name")
    if any(not s.gex_dir or not s.vdj_dir for s in sheet):
        raise ValueError("Every sample requires paired gex_dir and vdj_dir inputs")
    if any(s.patient_id is not None for s in sheet) and any(
        s.patient_id is None or not str(s.patient_id).strip() for s in sheet
    ):
        raise ValueError("patient_id must be populated for every sample when the column is used")
    logger.info("Context: %s; T cells: %s; signatures: %s", args.context, args.tcell_type, ", ".join(names))
    adata = load_samples(sheet, **{key: getattr(args, key) for key in QC_DEFAULTS})
    if "sample" not in adata.obs or adata.obs["sample"].isna().any() or adata.obs["sample"].nunique() < 1:
        raise ValueError("prioritize requires at least one named sample")
    _analysis_units(adata)  # Validate before a cell filter can hide missing donors.
    _ensure_clone_key(adata)
    adata = phenotype_cells(adata, min_cd3_reads=args.min_cd3)
    adata, cell_filters, sample_counts = filter_cells(adata, args)

    clonotypes, clone_sample = _aggregate_within_patients(
        adata,
        aggregate_clonotypes,
        build_clone_sample_long,
        min_umi=args.min_vdj_umis,
    )

    # score_genes expects log-normalized expression. Using a copy preserves the
    # raw-count matrix used by the rest of the workflow.
    scored_adata = adata.copy()
    sc.pp.normalize_total(scored_adata, target_sum=10_000)
    sc.pp.log1p(scored_adata)

    signature_cols = []
    for name in names:
        col = f"signature_{name}"
        scored_adata.obs[col] = _score_within_samples(scored_adata, name)
        signature_cols.append(col)

    per_cell_scores = scored_adata.obs[["CDR3ab", "sample", "lineage", *signature_cols]].dropna(
        subset=["CDR3ab"]
    )
    per_sample_scores = (
        per_cell_scores.groupby(["CDR3ab", "sample", "lineage"], observed=True)[signature_cols]
        .mean()
        .reset_index()
    )
    clone_sample = clone_sample.merge(
        per_sample_scores,
        on=["CDR3ab", "sample", "lineage"],
        how="left",
        validate="one_to_one",
    )
    signature_percentile_cols = []
    for col in signature_cols:
        percentile_col = f"{col}_percentile"
        clone_sample[percentile_col] = clone_sample.groupby(
            ["sample", "lineage"],
            observed=True,
        )[col].rank(method="average", pct=True)
        signature_percentile_cols.append(percentile_col)
    clone_sample.to_csv(args.output_dir / "clone_sample_scores.csv", index=False)

    # A strong signal in any sample is retained in the clone-level review table;
    # the long table above shows which sample supplied it.
    signature_max = (
        clone_sample.groupby(["donor", "CDR3ab"], observed=True)[
            [*signature_cols, *signature_percentile_cols]
        ]
        .max()
        .add_suffix("_max")
        .reset_index()
    )
    clonotypes = clonotypes.merge(
        signature_max,
        on=["donor", "CDR3ab"],
        how="left",
        validate="one_to_one",
    )

    clonotypes = add_pgen_ppost(clonotypes, backend="kmer")
    clonotypes = add_paired_ppost(clonotypes)
    clonotypes = add_pairing_promiscuity(clonotypes)
    clonotypes = annotate_clonotypes(
        clonotypes,
        vdjdb_path=args.vdjdb,
        iedb_path=args.iedb,
        cedar_path=args.cedar,
        match_strictness=args.database_match,
        flag_only=True,
    )

    clonotypes["known_viral_match"] = clonotypes["is_viral"].fillna(False).astype(bool)
    clonotypes["known_mart1_match"] = _known_mart1_mask(clonotypes)
    clonotypes["uses_trav12_2"] = _trav12_2_mask(clonotypes)

    # ppost_either is ln(P) for the more common of the two chains. Higher
    # (less negative) values are more public; the percentile is cohort-relative.
    clonotypes["publicness_percentile"] = clonotypes.groupby("donor", observed=True)[
        "ppost_either"
    ].rank(method="average", pct=True)
    clonotypes["high_publicness"] = False
    if args.exclude_public_quantile is not None:
        clonotypes["high_publicness"] = (
            clonotypes["publicness_percentile"] >= args.exclude_public_quantile
        ).fillna(False)

    percentile_cols = [f"{col}_max" for col in signature_percentile_cols]
    clonotypes["signature_support_count"] = (
        clonotypes[percentile_cols].ge(args.signature_quantile).sum(axis=1)
    )
    clonotypes["meets_abundance"] = clonotypes["cell_count"].ge(args.min_cells) & clonotypes[
        "max_frequency"
    ].ge(args.min_frequency)
    for chain in ("alpha", "beta"):
        clonotypes[f"{chain}_cdr3_length"] = clonotypes[f"CDR3_{chain}"].str.len()
        clonotypes[f"short_{chain}_cdr3"] = clonotypes[f"{chain}_cdr3_length"].lt(
            getattr(args, f"min_{chain}_cdr3_length")
        )
    clonotypes["gene_filter_failed"] = ~(
        gene_filter_mask(clonotypes, args.include_v_gene, args.exclude_v_gene, "v")
        & gene_filter_mask(clonotypes, args.include_j_gene, args.exclude_j_gene, "j")
    )

    risk_flags = [
        ("known_viral_match", "known_viral"),
        ("known_mart1_match", "known_MART1"),
        ("uses_trav12_2", "TRAV12-2"),
        ("high_publicness", "high_publicness"),
        ("alpha_promiscuous", "promiscuous_alpha"),
    ]
    clonotypes["risk_flags"] = _join_flags(clonotypes, risk_flags)

    excluded = pd.Series(False, index=clonotypes.index)
    exclusion_flags = [("short_alpha_cdr3", "short_alpha_cdr3"),
                       ("short_beta_cdr3", "short_beta_cdr3"),
                       ("gene_filter_failed", "gene_filter")]
    for col, _ in exclusion_flags:
        excluded |= clonotypes[col]
    if args.exclude_known_viral:
        excluded |= clonotypes["known_viral_match"]
        exclusion_flags.append(("known_viral_match", "known_viral"))
    if args.exclude_known_mart1:
        excluded |= clonotypes["known_mart1_match"]
        exclusion_flags.append(("known_mart1_match", "known_MART1"))
    if args.exclude_trav12_2:
        excluded |= clonotypes["uses_trav12_2"]
        exclusion_flags.append(("uses_trav12_2", "TRAV12-2"))
    if args.exclude_public_quantile is not None:
        excluded |= clonotypes["high_publicness"]
        exclusion_flags.append(("high_publicness", "high_publicness"))
    clonotypes["excluded_reason"] = _join_flags(clonotypes, exclusion_flags)

    clonotypes["eligible_for_review"] = (
        clonotypes["meets_abundance"]
        & clonotypes["signature_support_count"].ge(args.min_signature_support)
        & ~excluded
    )
    clonotypes = select_round_robin(
        clonotypes, clone_sample, names, args.max_clones, args.signature_quantile,
    )
    candidates = clonotypes[clonotypes["selected_for_review"]].copy()

    clonotypes.to_csv(args.output_dir / "all_scored_clones.csv", index=False)
    candidates.to_csv(args.output_dir / "candidate_clones.csv", index=False)
    from .version import __version__

    config = {
        "version": __version__, "arguments": vars(args) | {"func": args.func.__name__},
        "resolved_signatures": names, "signature_lineages": SIGNATURE_LINEAGE,
        "cell_filters": cell_filters, "sample_cell_counts": sample_counts,
        "selection": "round-robin per patient over signature/sample/lineage lists",
        "database_annotation_available": bool(args.vdjdb or args.iedb or args.cedar),
        "candidate_count": len(candidates),
    }
    (args.output_dir / "prioritization.json").write_text(json.dumps(config, indent=2, default=str) + "\n")
    if args.exclude_known_mart1 and not config["database_annotation_available"]:
        logger.warning("MART-1 exclusion is enabled, but no reference database was supplied; "
                       "unannotated clones cannot be identified as MART-1 matches")
    return clonotypes, candidates


def add_cli_args(parser: argparse.ArgumentParser, *, context="generic") -> None:
    """Register the sample-sheet workflow options on a CLI parser."""
    parser.add_argument("sample_sheet", type=Path, help="Paired CellRanger VDJ + GEX YAML/CSV (one or more samples)")
    parser.add_argument("-o", "--output-dir", type=Path, required=True,
                        help="Directory for the candidate, audit, and per-sample score CSVs")
    parser.add_argument("--vdjdb", type=Path, help="Optional VDJdb annotation file")
    parser.add_argument("--iedb", type=Path, help="Optional IEDB annotation file")
    parser.add_argument("--cedar", type=Path, help="Optional CEDAR annotation file")
    parser.add_argument(
        "--database-match",
        choices=("strict_ab", "ab_with_partial", "b_only"),
        default="ab_with_partial",
        help="Public-database match strictness (default: ab_with_partial)",
    )
    parser.add_argument("--min-cells", type=int, default=2,
                        help="Minimum clone cell count across a patient's samples (default: 2)")
    parser.add_argument("--min-frequency", type=float, default=0.001,
                        help="Minimum clone frequency in at least one sample (default: 0.001)")
    parser.add_argument("--signature-quantile", type=float, default=0.0,
                        help="Optional within-sample/lineage percentile floor in [0, 1] (default: 0, disabled)")
    parser.add_argument("--min-signature-support", type=int, default=1,
                        help="Minimum signatures meeting the percentile floor (default: 1)")
    parser.add_argument(
        "--exclude-known-viral",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Exclude known viral database matches (default: disabled)",
    )
    parser.add_argument(
        "--exclude-known-mart1",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude known MART-1/Melan-A database matches (default: enabled)",
    )
    parser.add_argument(
        "--exclude-trav12-2",
        action="store_true",
        help="Aggressive heuristic; TRAV12-2 alone does not establish MART-1 specificity",
    )
    parser.add_argument(
        "--exclude-public-quantile",
        type=float,
        default=None,
        metavar="Q",
        help="Optionally remove the most public cohort fraction, e.g. 0.90",
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose logging and error tracebacks")
    parser.add_argument("--context", choices=CONTEXT_SIGNATURES, default=context,
                        help=f"Signature preset (default: {context}); blood includes vaccine studies, heme includes AML")
    parser.add_argument("--tcell-type", choices=("cd8", "cd4", "both"), default="cd8",
                        help="Retain classified CD8, CD4, or both (default: cd8); unknown cells are excluded")
    parser.add_argument("--signatures", nargs="+", help="Replace the preset with named signatures in this order")
    parser.add_argument("--exclude-signatures", nargs="+", help="Remove named signatures from the selected set")
    parser.add_argument("--max-clones", type=int, default=100,
                        help="Maximum unique clones per patient, selected round-robin (default: 100; 0 = unlimited)")
    for key, default in QC_DEFAULTS.items():
        parser.add_argument(f"--{key.replace('_', '-')}", type=type(default), default=default,
                            help=f"Per-cell GEX QC (default: {default})")
    for option, default, help_text in (
        ("min-vdj-umis", 2, "Minimum UMIs in EACH primary VDJ chain per cell"),
        ("min-vdj-reads", 0, "Minimum reads in EACH primary VDJ chain per cell"),
        ("min-cd3", 10, "Minimum summed CD3D/E/G raw GEX counts per cell"),
        ("min-alpha-cdr3-length", 0, "Minimum alpha CDR3 amino-acid length, including anchors"),
        ("min-beta-cdr3-length", 0, "Minimum beta CDR3 amino-acid length, including anchors"),
    ):
        parser.add_argument(f"--{option}", type=int, default=default,
                            help=f"{help_text} (default: {default}; 0 = disabled)")
    for direction in ("min", "max"):
        parser.add_argument(f"--{direction}-expression", action="append", default=[], metavar="GENE=VALUE",
                            help="Per-cell log1p(CP10K) marker threshold; repeat for AND conditions")
    for action in ("include", "exclude"):
        for segment in ("v", "j"):
            parser.add_argument(f"--{action}-{segment}-gene", action="append", default=[], metavar="GENE",
                                help=f"{action.capitalize()} exact {segment.upper()} genes in either chain; repeatable, allele-insensitive")
