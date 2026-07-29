#!/usr/bin/env python3
"""Prioritize clonotypes from multiple paired 10x VDJ + GEX TIL samples.

This is an auditable example, not a specificity predictor. It combines:

- observed expansion (cells and within-sample frequency),
- several TIL/neoantigen-reactivity expression signatures,
- exact public-database annotation,
- sequence-background publicness and alpha-chain promiscuity flags.

It writes every scored clone before writing the filtered shortlist, so no
exclusion is hidden.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Tumor-local signatures are useful here. NeoTCR_PBL is intentionally omitted:
# it was derived for circulating cells, not TILs.
SIGNATURE_NAMES = (
    "TumorReactive",
    "Cytolytic",
    "Differentiated",
    "MANAscore",
    "NeoTCR8",
    "NeoTCR4",
)
MART1_PATTERN = re.compile(
    r"MART[- ]?1|MELAN[- ]?A|MLANA|E(?:AA|LA)GIGILTV|AAGIGILTV",
    flags=re.IGNORECASE,
)


def _score_within_samples(adata, name: str) -> np.ndarray:
    """Score one signature separately within each sample."""
    from tcrsift import score_by_name

    scores = np.full(adata.n_obs, np.nan, dtype=float)
    samples = adata.obs["sample"].astype(str).to_numpy()
    for sample in pd.unique(samples):
        positions = np.flatnonzero(samples == sample)
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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate shared samples together without merging clones across patients."""
    units = _analysis_units(adata)
    clonotype_tables = []
    clone_sample_tables = []
    for patient in pd.unique(units):
        patient_adata = adata[units.eq(patient).to_numpy()].copy()
        patient_clones = aggregate_clonotypes(patient_adata)
        patient_clones.insert(0, "donor", patient)
        clonotype_tables.append(patient_clones)

        patient_long = build_clone_sample_long(patient_adata)
        patient_long["donor"] = patient
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


def run(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not 0 < args.signature_quantile <= 1:
        raise ValueError("--signature-quantile must be in (0, 1]")
    if args.exclude_public_quantile is not None and not (0 < args.exclude_public_quantile <= 1):
        raise ValueError("--exclude-public-quantile must be in (0, 1]")
    if not 0 <= args.min_frequency <= 1:
        raise ValueError("--min-frequency must be in [0, 1]")
    if args.min_cells < 1:
        raise ValueError("--min-cells must be >= 1")
    if not 1 <= args.min_signature_support <= len(SIGNATURE_NAMES):
        raise ValueError(f"--min-signature-support must be between 1 and {len(SIGNATURE_NAMES)}")

    import scanpy as sc

    from tcrsift import (
        add_paired_ppost,
        add_pgen_ppost,
        aggregate_clonotypes,
        annotate_clonotypes,
        build_clone_sample_long,
        load_samples,
        phenotype_cells,
    )
    from tcrsift.annotate_tcrs import add_pairing_promiscuity

    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_samples(args.sample_sheet)
    if "sample" not in adata.obs or adata.obs["sample"].nunique() < 2:
        raise ValueError("This example requires at least two named samples")
    _ensure_clone_key(adata)
    adata = phenotype_cells(adata)

    clonotypes, clone_sample = _aggregate_within_patients(
        adata,
        aggregate_clonotypes,
        build_clone_sample_long,
    )

    # score_genes expects log-normalized expression. Using a copy preserves the
    # raw-count matrix used by the rest of the workflow.
    scored_adata = adata.copy()
    sc.pp.normalize_total(scored_adata, target_sum=10_000)
    sc.pp.log1p(scored_adata)

    signature_cols = []
    for name in SIGNATURE_NAMES:
        col = f"signature_{name}"
        scored_adata.obs[col] = _score_within_samples(scored_adata, name)
        signature_cols.append(col)

    per_cell_scores = scored_adata.obs[["CDR3ab", "sample", *signature_cols]].dropna(
        subset=["CDR3ab"]
    )
    per_sample_scores = (
        per_cell_scores.groupby(["CDR3ab", "sample"], observed=True)[signature_cols]
        .mean()
        .reset_index()
    )
    clone_sample = clone_sample.merge(
        per_sample_scores,
        on=["CDR3ab", "sample"],
        how="left",
        validate="one_to_one",
    )
    signature_percentile_cols = []
    for col in signature_cols:
        percentile_col = f"{col}_percentile"
        clone_sample[percentile_col] = clone_sample.groupby(
            "sample",
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
    clonotypes["publicness_percentile"] = pd.to_numeric(
        clonotypes.get("ppost_either"), errors="coerce"
    ).rank(method="average", pct=True)
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

    risk_flags = [
        ("known_viral_match", "known_viral"),
        ("known_mart1_match", "known_MART1"),
        ("uses_trav12_2", "TRAV12-2"),
        ("high_publicness", "high_publicness"),
        ("alpha_promiscuous", "promiscuous_alpha"),
    ]
    clonotypes["risk_flags"] = _join_flags(clonotypes, risk_flags)

    excluded = pd.Series(False, index=clonotypes.index)
    exclusion_flags = []
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

    clonotypes["selected_for_review"] = (
        clonotypes["meets_abundance"]
        & clonotypes["signature_support_count"].ge(args.min_signature_support)
        & ~excluded
    )
    clonotypes = clonotypes.sort_values(
        ["selected_for_review", "signature_support_count", "max_frequency", "cell_count"],
        ascending=[False, False, False, False],
    )
    candidates = clonotypes[clonotypes["selected_for_review"]].copy()

    clonotypes.to_csv(args.output_dir / "all_scored_clones.csv", index=False)
    candidates.to_csv(args.output_dir / "candidate_clones.csv", index=False)
    return clonotypes, candidates


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path, help="Multi-sample TIL YAML/CSV")
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument("--vdjdb", type=Path)
    parser.add_argument("--iedb", type=Path)
    parser.add_argument("--cedar", type=Path)
    parser.add_argument(
        "--database-match",
        choices=("strict_ab", "ab_with_partial", "b_only"),
        default="ab_with_partial",
        help="Public-database match strictness (default: ab_with_partial)",
    )
    parser.add_argument("--min-cells", type=int, default=2)
    parser.add_argument("--min-frequency", type=float, default=0.001)
    parser.add_argument("--signature-quantile", type=float, default=0.90)
    parser.add_argument("--min-signature-support", type=int, default=1)
    parser.add_argument(
        "--exclude-known-viral",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--exclude-known-mart1",
        action=argparse.BooleanOptionalAction,
        default=True,
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
    parser.add_argument("--verbose", action="store_true")
    return parser


def main() -> None:
    args = create_parser().parse_args()
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)
    _, candidates = run(args)
    print(f"Wrote {len(candidates)} candidates and the full audit table to {args.output_dir}")


if __name__ == "__main__":
    main()
