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

"""End-to-end TCRsift antigen-enrichment pipeline (:func:`run_pipeline`).

Provides the Python API backing ``tcrsift run``, orchestrating sample loading,
T-cell phenotyping, optional single-cell atlas embedding, paired stimulation
scoring, clonotype aggregation, statistical/threshold filtering, public TCR
database annotation, TIL matching, full-length construct assembly, and QC/PDF
reporting.
"""

from __future__ import annotations

import copy
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._schema import paired_cdr3ab_key
from .config import FilterConfig, TCRsiftConfig
from .sample_sheet import SampleSheet, load_sample_sheet
from .validation import TCRsiftValidationError

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Structured output returned by :func:`run_pipeline`."""

    output_dir: Path
    data_dir: Path
    plots_dir: Path
    species: str
    clonotypes: pd.DataFrame
    filtered: pd.DataFrame
    annotated: pd.DataFrame
    til_matched: pd.DataFrame
    assembled: pd.DataFrame | None
    stimulation_response: pd.DataFrame | None
    funnel_counts: dict[str, int]

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)


def _resolve_sheet_path(path_val: str | Path, base_dir: Path | None) -> Path:
    """Resolve a sample-sheet path relative to ``base_dir`` when not absolute."""
    p = Path(str(path_val))
    if not p.is_absolute() and base_dir is not None:
        return (base_dir / p).resolve()
    return p.resolve()


def _common_contig_root(vdj_dirs, *, base_dir: Path | None = None) -> str | None:
    """Common-ancestor directory of the sample sheet's ``vdj_dirs``, to scan for
    contigs in ``'sheet'`` mode (#201). Returns ``None`` when there is no shared
    ancestor or it resolves to the filesystem root (too broad to rglob safely)."""
    dir_paths = []
    for d in vdj_dirs:
        if not d:
            continue
        p_raw = _resolve_sheet_path(d, base_dir)
        dir_paths.append(str(p_raw.parent if p_raw.is_file() else p_raw))
    if not dir_paths:
        return None
    try:
        root = os.path.commonpath(dir_paths)
    except ValueError:
        return None  # mixed drives / no common path
    p = Path(root)
    if p == Path(p.anchor):  # filesystem root — refuse
        return None
    return root


def _prepare_signature_expression(adata, genes: list[str]) -> tuple[pd.DataFrame, bool]:
    """Return ``(expr_df, apply_log1p)`` normalized to CP10K when ``adata.X`` holds raw counts."""
    from .signature_methods import expression_frame_from_adata

    expr = expression_frame_from_adata(adata, genes, on_missing="ignore")
    already_log1p = "log1p" in getattr(adata, "uns", {})
    apply_log1p = not already_log1p
    if (
        not already_log1p
        and getattr(adata, "X", None) is not None
        and adata.n_obs > 0
        and adata.n_vars > 0
    ):
        x_mat = adata.X
        x_min, x_max = float(x_mat.min()), float(x_mat.max())
        x_vals = x_mat.data if hasattr(x_mat, "data") else np.asarray(x_mat).ravel()
        sample_vals = x_vals[: min(len(x_vals), 10000)]
        is_int_like = len(sample_vals) == 0 or bool(
            np.allclose(sample_vals, np.round(sample_vals), atol=1e-4)
        )
        if x_min >= 0 and (x_max > 25 or is_int_like):
            raw_sums = np.asarray(x_mat.sum(axis=1)).ravel().astype(float)
            if np.nanmax(raw_sums) > 50:
                lib_safe = np.where(raw_sums > 0, raw_sums, np.nan)
                expr = expr.div(lib_safe, axis=0).fillna(0.0) * 1e4
        elif x_min >= 0 and x_max <= 25 and not is_int_like:
            apply_log1p = False
    return expr, apply_log1p


def run_pipeline(
    sample_sheet: str | Path | SampleSheet,
    output_dir: str | Path,
    config: TCRsiftConfig | None = None,
    *,
    til_sample_specs: list[str] | None = None,
    cache_dir: str | Path | None = None,
    no_assemble: bool = False,
    no_sheet_contigs: bool = False,
    allow_canonical_fallback: bool = False,
    report_name: str | None = None,
) -> PipelineResult:
    """Run the end-to-end TCRsift culture/stimulation pipeline.

    Parameters
    ----------
    sample_sheet : str, Path, or SampleSheet
        Path to a YAML/CSV sample sheet or a pre-built :class:`SampleSheet`.
    output_dir : str or Path
        Output directory for ``data/``, ``plots/``, ``config.yaml``, and reports.
    config : TCRsiftConfig, optional
        Pipeline configuration (defaults to ``TCRsiftConfig()``).
    til_sample_specs : list of str, optional
        Repeatable ``NAME=TYPE:PATH`` TIL sample specifications (overrides sheet TILs).
    cache_dir : str or Path, optional
        Reference-database cache directory override.
    no_assemble : bool
        Skip Step 7 whole-repertoire construct assembly when ``True``.
    no_sheet_contigs : bool
        Disable automatic sample-sheet contig assembly root detection (#201).
    allow_canonical_fallback : bool
        Permit canonical constant/leader fallback when donor contigs cannot verify fidelity.
    report_name : str, optional
        Explicit report filename stem override.

    Returns
    -------
    PipelineResult
        Structured summary of output paths, intermediate/final DataFrames, and funnel counts.
    """
    from ._schema import truthy_flag_series
    from .annotate import annotate_clonotypes
    from .assemble import (
        assemble_full_sequences,
        collect_germline_variants,
        enforce_contig_fidelity,
        export_fasta,
    )
    from .clonotype import aggregate_clonotypes, export_clonotypes_airr
    from .datacache import ensure_cached
    from .filter import (
        bucket_by_donor_sharing,
        filter_clonotypes,
        pick_strictest_populated_tier,
        resolve_fdr_scope,
        resolve_fdr_tiers_for_method,
        resolve_filter_mode_kwargs,
        split_by_tier,
        split_clonotypes_by_donor,
        split_clonotypes_by_sample,
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
    from .provenance import write_provenance
    from .species import resolve_species
    from .til import load_til_specs, match_til

    if config is None:
        config = TCRsiftConfig()
    else:
        config = copy.deepcopy(config)

    set_plot_format(config.output.plot_format)

    if config.output.baseline_markers is not None:
        from .format import set_baseline_markers

        set_baseline_markers(*config.output.baseline_markers)

    sheet_obj = (
        sample_sheet
        if isinstance(sample_sheet, SampleSheet)
        else load_sample_sheet(sample_sheet)
    )

    if not config.til.til_samples:
        til_samples = sheet_obj.get_til_samples()
        if til_samples:
            til_sample_names = [s.sample for s in til_samples]
            config.til.til_samples = til_sample_names
            print(f"Auto-detected TIL samples from sample sheet: {til_sample_names}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    config.to_yaml(output_dir / "config.yaml")
    write_provenance(output_dir)

    print("=" * 60)
    print("TCRsift Pipeline")
    print("=" * 60)

    funnel_counts: dict[str, int] = {}

    # Step 1: Load
    print("\n[1/7] Loading samples...")
    requested_species = (
        config.species
        if config.species != "auto"
        else (getattr(sheet_obj, "species", None) or "auto")
    )
    adata = load_samples(
        sheet_obj if isinstance(sample_sheet, SampleSheet) else sample_sheet,
        min_genes=config.load.min_genes,
        max_genes=config.load.max_genes,
        min_counts=config.load.min_counts,
        max_counts=config.load.max_counts,
        min_mito_pct=config.load.min_mito_pct,
        max_mito_pct=config.load.max_mito_pct,
        species=requested_species,
    )
    run_species = resolve_species(requested_species, adata=adata, warn_mismatch=False)
    adata.write_h5ad(data_dir / "loaded.h5ad")
    funnel_counts["Raw Cells"] = len(adata)
    print(f"  Loaded {len(adata)} cells")

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

    if "CDR3_beta" in adata.obs.columns:
        n_with_vdj = int(adata.obs["CDR3_beta"].notna().sum())
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

    tcell_type = config.filter.tcell_type
    if tcell_type != "both":
        adata = filter_by_tcell_type(adata, tcell_type)
        print(f"  Filtered to {tcell_type.upper()}+: {len(adata)} cells")

    min_cd3 = config.phenotype.min_cd3_reads
    n_before_cd3 = len(adata)
    adata = apply_cd3_gate(adata, min_cd3)
    if len(adata) != n_before_cd3:
        print(f"  CD3 gate (>={min_cd3} reads): {len(adata)}/{n_before_cd3} cells pass")

    funnel_counts["Phenotyped"] = len(adata)

    if config.embed.enabled:
        print("\n  Embedding cells (Pearson residuals → PCA → Harmony → UMAP → Leiden)...")
        try:
            from .annotate_cells import annotate_cells, gates_from_phenotype_config
            from .embed import embed_cells
            from .plots import plot_umap

            atlas = embed_cells(adata, config.embed)
            annotate_cells(
                atlas,
                allowed_types=config.phenotype.allowed_types,
                gates=gates_from_phenotype_config(config.phenotype, species=run_species),
                min_subtype_frac=config.phenotype.min_subtype_frac,
                species=run_species,
            )
            atlas.write_h5ad(data_dir / "atlas.h5ad")
            print(
                f"  Atlas: {atlas.obs['leiden'].nunique()} Leiden clusters, "
                f"{atlas.obs['phenotype_base'].nunique()} cell types "
                f"→ {data_dir / 'atlas.h5ad'}"
            )
            if config.output.generate_plots:
                plot_umap(atlas, "phenotype", plots_dir / "atlas_umap_phenotype.png")
                plot_umap(atlas, "leiden", plots_dir / "atlas_umap_leiden.png")
                print(f"  Atlas UMAPs (phenotype, leiden) → {plots_dir}")
        except ImportError as e:
            print(f"  [embed] skipped — install tcrsift[atlas] to enable ({e})")

    # Step 3: Clonotype (with paired control comparison when source == 'control' is present)
    stim_resp_df = None
    stim_join_col = "CDR3ab" if config.clonotype.group_by == "CDR3ab" else "CDR3_beta"
    if "source" in adata.obs.columns and (adata.obs["source"].astype(str) == "control").any():
        non_ctrl_mask = adata.obs["source"].astype(str) != "control"
        n_ctrl_cells = int((~non_ctrl_mask).sum())
        if int(non_ctrl_mask.sum()) == 0:
            raise TCRsiftValidationError(
                "No non-control stimulated/culture cells remain for clonotype aggregation",
                hint="Provide at least one positive sample (e.g. source='culture') alongside source='control'.",
            )
        try:
            from .filter import score_stimulation_response

            stim_score_cols: list[str] = []
            if getattr(adata, "X", None) is not None and adata.n_vars > 0:
                import scanpy as sc

                from .signature_methods import score_by_name

                stim_adata = adata.copy()
                already_log1p = "log1p" in getattr(stim_adata, "uns", {})
                if not already_log1p and stim_adata.n_obs and stim_adata.n_vars:
                    x_min = float(stim_adata.X.min())
                    x_max = float(stim_adata.X.max())
                    if x_min >= 0 and x_max > 25:
                        sc.pp.normalize_total(stim_adata, target_sum=10_000)
                        sc.pp.log1p(stim_adata)
                for sig_name in ("AIM", "AcuteActivation"):
                    try:
                        vals = score_by_name(
                            stim_adata,
                            sig_name,
                            log1p=False,
                            on_missing="warn",
                            species=run_species,
                        )
                        col_name = f"signature_{sig_name}"
                        adata.obs[col_name] = vals.to_numpy(dtype=float)
                        stim_score_cols.append(col_name)
                    except Exception as sig_exc:
                        logger.debug("Skipping stimulation %s score: %s", sig_name, sig_exc)
                del stim_adata

            if (
                "patient_id" in adata.obs.columns
                and adata.obs["patient_id"].dropna().nunique() > 1
            ):
                per_donor_resp = []
                for _, donor_idx in adata.obs.groupby("patient_id", observed=True).indices.items():
                    donor_adata = adata[donor_idx]
                    donor_src = donor_adata.obs["source"].astype(str)
                    if (donor_src == "control").any() and (donor_src != "control").any():
                        per_donor_resp.append(
                            score_stimulation_response(
                                donor_adata,
                                clone_col=stim_join_col,
                                score_cols=stim_score_cols or None,
                            )
                        )
                stim_resp_df = (
                    pd.concat(per_donor_resp, ignore_index=True)
                    .sort_values(["fisher_pvalue", "stim_cells"], ascending=[True, False], kind="mergesort")
                    .drop_duplicates(stim_join_col)
                    if per_donor_resp
                    else None
                )
            else:
                stim_resp_df = score_stimulation_response(
                    adata,
                    clone_col=stim_join_col,
                    score_cols=stim_score_cols or None,
                )
            if stim_resp_df is not None:
                stim_resp_df.to_csv(data_dir / "stimulation_response.csv", index=False)
                print(
                    f"  Scored stimulation vs control ({n_ctrl_cells} control cells) -> stimulation_response.csv"
                )
        except TCRsiftValidationError as exc:
            print(f"  [stimulation] skipped control comparison: {exc}")
        adata = adata[non_ctrl_mask].copy()
        print(f"  Excluded {n_ctrl_cells} control cells before positive clonotype aggregation")

    print("\n[3/7] Aggregating clonotypes...")
    clonotypes = aggregate_clonotypes(
        adata,
        group_by=config.clonotype.group_by,
        min_umi=config.clonotype.min_umi,
        handle_doublets=config.clonotype.handle_doublets,
        attribution=config.attribution,
        doublet_warn_rate=config.clonotype.doublet_warn_rate,
    )
    if stim_resp_df is not None and stim_join_col in clonotypes.columns:
        extra_stim_cols = [
            c
            for c in stim_resp_df.columns
            if c.startswith(("stim_mean_", "control_mean_", "delta_"))
        ]
        stim_cols = (
            stim_resp_df[
                [
                    stim_join_col,
                    "control_cells",
                    "control_frequency",
                    "log2_fold_change",
                    "fisher_pvalue",
                    "fisher_fdr",
                    *extra_stim_cols,
                ]
            ]
            .rename(
                columns={
                    "log2_fold_change": "stim_log2_fc",
                    "fisher_pvalue": "stim_fisher_pvalue",
                    "fisher_fdr": "stim_fisher_fdr",
                }
            )
            .drop_duplicates(stim_join_col)
        )
        clonotypes = clonotypes.merge(stim_cols, on=stim_join_col, how="left")
        clonotypes["control_cells"] = clonotypes["control_cells"].fillna(0).astype(int)
        clonotypes["control_frequency"] = clonotypes["control_frequency"].fillna(0.0)

    if getattr(config.annotate, "call_multimers", False):
        from .multimer import call_multimers_for_adata

        mm_df = call_multimers_for_adata(
            adata,
            sheet_obj,
            clone_col=stim_join_col,
            min_umi=config.annotate.multimer_min_umi,
            min_nc_ratio=config.annotate.multimer_min_nc_ratio,
            min_specific_cells=config.annotate.multimer_min_specific_cells,
            min_concordance=config.annotate.multimer_min_concordance,
        )
        if not mm_df.empty and stim_join_col in clonotypes.columns:
            mm_df.to_csv(data_dir / "multimer_specificity.csv", index=False)
            clonotypes = clonotypes.merge(
                mm_df.drop(columns=["clone_cells"], errors="ignore"),
                on=stim_join_col,
                how="left",
            )
            clonotypes["multimer_confident"] = truthy_flag_series(
                clonotypes.get("multimer_confident", pd.Series(False, index=clonotypes.index))
            )
            print(
                f"  Called pMHC multimers across {int(clonotypes['multimer_confident'].sum())} confident clones -> multimer_specificity.csv"
            )
    if config.attribution.enabled:
        print("  Attribution ON: weighted fractional clone counts (#176)")

    if config.clonotype.add_pgen_ppost:
        from .annotate_tcrs import add_pgen_ppost

        clonotypes = add_pgen_ppost(clonotypes, backend="kmer", species=run_species)
        scored = [c for c in ("ppost_alpha", "ppost_beta") if c in clonotypes.columns]
        if scored:
            print(f"  Added k-mer Pgen/Ppost → {', '.join(scored)} (for PRISM)")

    clonotypes.to_csv(data_dir / "clonotypes.csv", index=False)
    funnel_counts["Clonotypes"] = len(clonotypes)
    print(f"  Found {len(clonotypes)} clonotypes")
    if config.output.generate_plots:
        plot_clonotypes(clonotypes, plots_dir)

    emit_long = config.output.emit_clone_sample_long
    n_samples_in_run = adata.obs["sample"].nunique() if "sample" in adata.obs.columns else 0
    should_emit_long = (
        emit_long == "always"
        or (emit_long == "auto" and n_samples_in_run >= 2)
    )
    have_method_axis = "enrichment_method" in adata.obs.columns
    has_selection_rules = bool(config.selection.get("rules"))

    long_df = None
    if should_emit_long or have_method_axis or has_selection_rules:
        from .clonotype import build_clone_sample_long

        long_df = build_clone_sample_long(adata, attribution=config.attribution)

    if long_df is not None and "sample" in long_df.columns:
        from .qc import clonal_expansion_metrics

        exp_qc = clonal_expansion_metrics(long_df)
        if not exp_qc.empty:
            exp_qc.to_csv(data_dir / "expansion_qc.csv", index=False)
            print(
                f"  Wrote expansion_qc.csv: {len(exp_qc)} samples "
                f"(top1% / effective-clones / clonality per sample)"
            )
            for _, r in exp_qc[exp_qc["unexpanded"]].iterrows():
                print(f"  ⚠️  sample {r['sample']!r}: {r['warning']}")

    if "sample" in adata.obs.columns:
        from .qc import inter_sample_barcode_jaccard, sample_integrity_qc

        sample_qc = sample_integrity_qc(adata)
        if not sample_qc.empty:
            sample_qc.to_csv(data_dir / "sample_qc.csv", index=False)
            print(
                f"  Wrote sample_qc.csv: {len(sample_qc)} samples "
                "(GEX↔VDJ overlap / barcode Jaccard / anchor integrity)"
            )
            jac = inter_sample_barcode_jaccard(adata)
            if not jac.empty:
                jac.to_csv(data_dir / "sample_barcode_jaccard.csv")
            for _, r in sample_qc[sample_qc["warning"] != ""].iterrows():
                print(f"  ⚠️  sample {r['sample']!r}: {r['warning']}")

    if "enrichment_method" in adata.obs.columns:
        try:
            from .sort_qc import sort_signature_consistency_from_adata

            sort_qc = sort_signature_consistency_from_adata(adata, species=run_species)
        except Exception as exc:
            logger.debug("sort-consistency QC skipped: %s", exc)
            sort_qc = pd.DataFrame()
        if not sort_qc.empty:
            sort_qc.to_csv(data_dir / "sort_consistency_qc.csv", index=False)
            print(
                f"  Wrote sort_consistency_qc.csv: {len(sort_qc)} "
                "(donor, sort) rows (signature enrichment vs label)"
            )
            for _, r in sort_qc[sort_qc["warning"] != ""].iterrows():
                print(f"  ⚠️  {r['donor']}/{r['enrichment_method']}: {r['warning']}")

    if should_emit_long and long_df is not None:
        from .selection import attach_per_sample_tiers

        emitted_long = attach_per_sample_tiers(long_df)
        emitted_long.to_csv(data_dir / "clone_sample_long.csv", index=False)
        print(
            f"  Wrote clone_sample_long.csv: {len(emitted_long)} rows "
            f"({emitted_long['CDR3ab'].nunique()} clones x {n_samples_in_run} samples)"
        )

    # Step 4: Filter
    print("\n[4/7] Filtering clonotypes...")
    user_filter_kwargs = {
        "min_donors": config.filter.min_donors,
        "min_methods_per_donor": config.filter.min_methods_per_donor,
        "min_cells_per_method": config.filter.min_cells_per_method,
        "min_frequency_per_method": config.filter.min_frequency_per_method,
        "min_timepoints": config.filter.min_timepoints,
        "min_timepoints_per_donor": config.filter.min_timepoints_per_donor,
        "min_apcs": config.filter.min_apcs,
        "min_apcs_per_donor": config.filter.min_apcs_per_donor,
        # min_til_cells_per_donor is applied in Step 6 after match_til populates
        # max_til_cells_per_donor (applying it here would silently no-op).
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

    _cohort_donors: set[str] = set()
    if "patient_id" in adata.obs.columns:
        _cohort_donors = {
            str(v).strip()
            for v in adata.obs["patient_id"].dropna().unique()
            if str(v).strip()
        }
    if not _cohort_donors and "donors" in clonotypes.columns and len(clonotypes):
        for _dval in clonotypes["donors"].dropna().astype(str):
            _cohort_donors.update(t.strip() for t in _dval.split(";") if t.strip())
    n_donors_in_run = (
        len(_cohort_donors)
        if _cohort_donors
        else (
            int(clonotypes["n_donors"].max())
            if "n_donors" in clonotypes.columns and len(clonotypes)
            else 1
        )
    )
    sheet_shares_antigen = getattr(sheet_obj, "donors_share_antigen", False)
    donors_share_antigen = config.filter.donors_share_antigen or sheet_shares_antigen
    resolved_scope = resolve_fdr_scope(
        config.filter.fdr_scope,
        n_donors=n_donors_in_run,
        donors_share_antigen=donors_share_antigen,
    )
    print(f"  FDR scope: {resolved_scope}")

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

    def _filter_one(df_in: pd.DataFrame) -> pd.DataFrame:
        return filter_clonotypes(
            df_in,
            method=config.filter.method,
            tcell_type=tcell_type,
            min_cells=config.filter.min_cells,
            min_frequency=config.filter.min_frequency,
            require_complete=config.filter.require_complete,
            fdr_tiers=fdr_tiers_arg,
            stim_fdr_threshold=config.filter.stim_fdr_threshold,
            stim_log2_fc_threshold=config.filter.stim_log2_fc_threshold,
            stim_delta_threshold=config.filter.stim_delta_threshold,
            **mode_kwargs,
        )

    if resolved_scope == "global":
        filtered = _filter_one(clonotypes)
        tier_dfs = split_by_tier(filtered)
        tier_counts = {}
        for tier, tier_df in tier_dfs.items():
            tier_df.to_csv(data_dir / f"filtered_{tier}.csv", index=False)
            tier_counts[tier] = len(tier_df)
            print(f"  {tier}: {len(tier_df)} clonotypes")
    else:
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
                print(
                    f"  WARNING: no clones survived filtering for "
                    f"{label}={key}: {exc}"
                )
                continue
            sub_filtered = sub_filtered.copy()
            sub_filtered[f"scope_{label}"] = key
            per_subset_filtered.append(sub_filtered)

            sub_tier_dfs = split_by_tier(sub_filtered)
            safe_key = str(key).replace("/", "_").replace(" ", "_")
            for tier, tier_df in sub_tier_dfs.items():
                tier_df.to_csv(
                    data_dir / f"filtered_{tier}_{safe_key}.csv",
                    index=False,
                )
                tier_counts[f"{tier}_{safe_key}"] = len(tier_df)
                print(f"  {label}={key} {tier}: {len(tier_df)} clonotypes")

        if per_subset_filtered:
            combined_sub = pd.concat(per_subset_filtered, ignore_index=True)
            scope_col = f"scope_{label}"
            tier_rank_map = {"tier1": 1, "tier2": 2, "tier3": 3, "tier4": 4, "tier5": 5}
            if "CDR3ab" in combined_sub.columns and combined_sub["CDR3ab"].duplicated().any():
                combined_sub = combined_sub.assign(
                    _trank=combined_sub.get("tier", pd.Series(dtype=object)).map(
                        lambda t: tier_rank_map.get(t, 99)
                    )
                ).sort_values(["_trank", "max_frequency"], ascending=[True, False], kind="mergesort")
                scopes_by_clone = (
                    combined_sub.groupby("CDR3ab", observed=True)[scope_col]
                    .apply(lambda s: ";".join(dict.fromkeys(str(x) for x in s.dropna())))
                    .to_dict()
                )
                deduped = combined_sub.drop_duplicates("CDR3ab", keep="first").drop(columns=["_trank"]).copy()
                deduped[scope_col] = deduped["CDR3ab"].map(scopes_by_clone)
            else:
                deduped = combined_sub.copy()
            if "CDR3ab" in deduped.columns and "CDR3ab" in clonotypes.columns:
                cohort_cols = [
                    c for c in (
                        "cell_count", "max_frequency", "mean_frequency", "n_samples",
                        "n_conditions", "n_methods", "max_methods_per_donor",
                        "max_cells_per_method", "max_frequency_per_method",
                        "n_timepoints", "max_timepoints_per_donor",
                        "n_apcs", "max_apcs_per_donor",
                    )
                    if c in clonotypes.columns and c in deduped.columns
                ]
                if cohort_cols:
                    base_lookup = clonotypes.drop_duplicates("CDR3ab").set_index("CDR3ab")[cohort_cols]
                    for c in cohort_cols:
                        deduped[c] = deduped["CDR3ab"].map(base_lookup[c]).values
            filtered = deduped.reset_index(drop=True)
        else:
            filtered = pd.DataFrame()

    if config.filter.filter_mode == "shared-high-freq":
        path = data_dir / "filtered_shared_high_freq.csv"
        filtered.to_csv(path, index=False)
        print(f"  shared-high-freq: {len(filtered)} clonotypes -> {path.name}")
    elif config.filter.filter_mode == "cross-donor-public":
        path = data_dir / "filtered_cross_donor_public.csv"
        filtered.to_csv(path, index=False)
        print(f"  cross-donor-public: {len(filtered)} clonotypes -> {path.name}")

    sharing_buckets = bucket_by_donor_sharing(filtered)
    for bucket, bucket_df in sharing_buckets.items():
        path = data_dir / f"filtered_{bucket}.csv"
        bucket_df.to_csv(path, index=False)
        print(f"  {bucket}: {len(bucket_df)} clonotypes -> {path.name}")

    if long_df is not None:
        from .clonotype import (
            build_method_overlap_matrices,
            build_method_recovery_table,
            build_per_method_rankings,
            build_selection_sets,
            set_overlap_table,
        )

        rankings = build_per_method_rankings(
            filtered, long_df, top_n=config.output.per_method_top_n
        )
        if rankings:
            method_dir = data_dir / "filtered_by_method"
            method_dir.mkdir(parents=True, exist_ok=True)
            for (donor, method), ranked in rankings.items():
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

        method_sets = build_selection_sets(
            long_df, set_col="method", clones=set(filtered["CDR3ab"]) if "CDR3ab" in filtered.columns else set()
        )
        if len(method_sets) >= 2:
            set_overlap_table(method_sets).to_csv(
                data_dir / "method_overlap_sets.csv", index=False
            )
            if config.output.generate_plots:
                from .plots import plot_set_overlap

                plot_set_overlap(
                    method_sets,
                    plots_dir / "method_overlap_upset.png",
                    title="Selected-clone overlap across methods",
                )
            print(f"  Method set overlap (UpSet): {len(method_sets)} methods")

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

    # Step 5: Annotate
    auto_dl = bool(getattr(config.annotate, "auto_download", False))
    vdjdb_path = config.annotate.vdjdb_path or ensure_cached("vdjdb", cache_dir, auto_download=auto_dl)
    iedb_path = config.annotate.iedb_path or ensure_cached("iedb", cache_dir, auto_download=auto_dl)
    cedar_path = config.annotate.cedar_path or ensure_cached("cedar", cache_dir, auto_download=auto_dl)
    iedb_epitope_path = getattr(config.annotate, "iedb_epitope_path", None) or ensure_cached(
        "iedb_epitope", cache_dir, auto_download=auto_dl
    )
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
            match_strictness=getattr(config.annotate, "match_strictness", None),
            match_mode=getattr(config.annotate, "match_mode", "exact"),
            max_distance=getattr(config.annotate, "max_distance", 1),
            exclude_viral=config.annotate.exclude_viral,
            flag_only=config.annotate.flag_only,
            species=run_species,
            host_species_filter=getattr(config.annotate, "host_species_filter", "all"),
        )
        annotated.to_csv(data_dir / "annotated.csv", index=False)
        n_viral = annotated["is_viral"].sum() if "is_viral" in annotated.columns else 0
        print(f"  Flagged {n_viral} viral clonotypes")
        if config.output.generate_plots:
            plot_annotations(annotated, plots_dir)
    else:
        print("\n[5/7] Skipping annotation (no databases provided)")

    # Step 6: TIL matching
    til_matched = annotated
    if til_sample_specs:
        print("\n[6/7] Matching against TIL samples from --til-sample...")
        til_data = load_til_specs(til_sample_specs)
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

        til_sample_objs = sheet_obj.get_til_samples()
        til_sample_names = [s.sample for s in til_sample_objs]
        if config.til.til_samples != til_sample_names:
            til_sample_objs = [s for s in til_sample_objs if s.sample in config.til.til_samples]

        if til_sample_objs:
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

    # Apply min_til_cells_per_donor now that match_til has populated max_til_cells_per_donor.
    min_til_per_donor = config.filter.min_til_cells_per_donor
    if min_til_per_donor is not None and min_til_per_donor > 0:
        if "max_til_cells_per_donor" in til_matched.columns:
            before_til = len(til_matched)
            til_matched = til_matched[
                pd.to_numeric(til_matched["max_til_cells_per_donor"], errors="coerce").fillna(0)
                >= min_til_per_donor
            ].copy()
            til_matched.to_csv(data_dir / "til_matched.csv", index=False)
            print(
                f"  Applied min_til_cells_per_donor >= {min_til_per_donor}: "
                f"{before_til:,} -> {len(til_matched):,} clonotypes"
            )
        else:
            logger.warning(
                "min_til_cells_per_donor=%s requested but no TIL matching ran "
                "(column 'max_til_cells_per_donor' absent); skipping.",
                min_til_per_donor,
            )

    # Step 7: Assemble
    assembled = None
    n_selected_clones: int | None = None
    n_non_viral: int | None = None
    rules = None
    if config.selection.get("rules") and long_df is not None:
        from .selection import select_from_clone_sample_long

        selection_long = long_df
        sig_names = config.selection.get("signatures")
        if sig_names:
            from .signature_methods import (
                build_signature_methods,
                signatures_for_species,
            )

            sig_registry = signatures_for_species(run_species)
            genes = sorted({
                g for n in sig_names for g in sig_registry[n].all_genes
            })
            expr, apply_sig_log1p = _prepare_signature_expression(adata, genes)
            obs_sig = adata.obs
            if "CDR3ab" not in obs_sig.columns and {
                "CDR3_alpha", "CDR3_beta"
            } <= set(obs_sig.columns):
                obs_sig = obs_sig.assign(CDR3ab=paired_cdr3ab_key(obs_sig)[0])
            if not ({"CDR3ab", "sample"} <= set(obs_sig.columns)):
                raise ValueError(
                    f"selection.signatures={sig_names} requested but "
                    "adata.obs lacks CDR3ab/sample — cannot build "
                    "signature methods."
                )
            sig_long = build_signature_methods(
                expr, obs_sig, signatures=sig_names,
                log1p=apply_sig_log1p,
                positive_method=config.selection.get(
                    "signature_positive_method", "gap"
                ),
                species=run_species,
            )
            selection_long = pd.concat(
                [long_df, sig_long], ignore_index=True
            )
            print(
                f"  Added {len(sig_long)} signature-method rows "
                f"({', '.join(sig_names)}) to the selection input"
            )
        exclude_clones = None
        if config.selection.get("exclude_viral", True) and "is_viral" in til_matched.columns:
            viral_mask = truthy_flag_series(til_matched["is_viral"])
            exclude_clones = set(
                til_matched.loc[viral_mask, "CDR3ab"].astype(str)
            )
            if exclude_clones:
                audit_cols = [
                    c for c in (
                        "CDR3ab", "CDR3_alpha", "CDR3_beta", "is_viral",
                        "db_protein_canonical", "db_antigen_gene", "db_protein",
                        "db_species", "db_epitope", "db_database", "db_category",
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
            selection_long,
            config.selection,
            exclude_clones=exclude_clones,
            clone_scores=til_matched,
        )

    has_leaders = (
        config.assemble.alpha_leader is not None or config.assemble.beta_leader is not None
    )
    assemble_enabled = getattr(config.assemble, "enabled", True) and not no_assemble
    if assemble_enabled and (
        config.assemble.single_chain or has_leaders or config.assemble.include_constant
    ):
        print("\n[7/7] Assembling full-length sequences...")
        _sample_name_from = getattr(config.assemble, "sample_name_from", "parent")
        _contigs_dir = config.assemble.contigs_dir

        _explicit_contigs = config.assemble.contigs_dir or getattr(
            config.assemble, "cellranger_dir", None
        )
        _sheet_base_dir = Path(sample_sheet).resolve().parent if sample_sheet else None
        _sheet_vdj_dirs = [
            str(_resolve_sheet_path(s.vdj_dir, _sheet_base_dir))
            for s in sheet_obj.samples
            if s.vdj_dir
        ]
        if (
            _sample_name_from in ("parent", "sheet")
            and not _explicit_contigs
            and _sheet_vdj_dirs
            and not no_sheet_contigs
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
                    {
                        "sample": s.sample,
                        "vdj_dir": str(_resolve_sheet_path(s.vdj_dir, _sheet_base_dir)),
                    }
                    for s in sheet_obj.samples
                    if s.vdj_dir
                ]
            )
            if _sample_name_from == "sheet"
            else None
        )
        _needs_leader_contigs = (
            config.assemble.alpha_leader == "from_contig"
            or config.assemble.beta_leader == "from_contig"
        )
        if _needs_leader_contigs and not (
            _contigs_dir or getattr(config.assemble, "cellranger_dir", None)
        ):
            raise TCRsiftValidationError(
                "Contigs directory required when extracting native leaders from contigs",
                hint="Provide --contigs-dir, --cellranger-dir, set assemble.contigs_dir "
                "in your config, or populate vdj_dir in your sample sheet.",
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
            trac_allele=getattr(config.assemble, "trac_allele", "auto"),
            trbc1_allele=getattr(config.assemble, "trbc1_allele", "auto"),
            trbc2_allele=getattr(config.assemble, "trbc2_allele", "auto"),
            stop_codons=tuple(
                getattr(config.assemble, "stop_codons", ("TAA", "TGA"))
            ),
            leader_fallback=getattr(config.assemble, "leader_fallback", None) or "germline",
            force_alpha_leader=getattr(config.assemble, "force_alpha_leader", None),
            force_beta_leader=getattr(config.assemble, "force_beta_leader", None),
            secondary_alpha_leader=getattr(config.assemble, "secondary_alpha_leader", None),
            secondary_beta_leader=getattr(config.assemble, "secondary_beta_leader", None),
            codon_optimization=getattr(config.assemble, "codon_optimization", "focal"),
            max_codon_repeats=getattr(config.assemble, "max_codon_repeats", 1),
            avoid_enzymes=tuple(getattr(config.assemble, "avoid_enzymes", []) or []),
            species=run_species,
            constant_species=getattr(config.assemble, "constant_species", None),
            codon_species=getattr(config.assemble, "codon_species", None),
            germline_vdj_fallback=bool(
                getattr(config.assemble, "germline_vdj_fallback", True)
            ),
            data_dir=cache_dir,
            auto_download_germline=bool(
                getattr(config.assemble, "auto_download_germline", False)
                or getattr(config.annotate, "auto_download", False)
            ),
        )
        assembled.to_csv(data_dir / "full_sequences.csv", index=False)
        print(f"  Assembled {len(assembled)} sequences")

        _gv = collect_germline_variants(assembled)
        if not _gv.empty:
            _gv.to_csv(data_dir / "germline_variants.csv", index=False)
            print(f"  Germline variants: {len(_gv)} distinct (by region) → germline_variants.csv")

        _fidelity = enforce_contig_fidelity(
            assembled,
            allow_canonical_fallback=allow_canonical_fallback,
            context="run",
        )
        if _fidelity:
            print(f"  WARNING: {_fidelity}")

        if config.output.generate_plots:
            plot_assembly(assembled, plots_dir)

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
                    rules
                    if rules is not None
                    else (
                        select_from_clone_sample_long(long_df, config.selection)
                        if config.selection.get("rules") else None
                    )
                )
                pdf_annotations = build_pdf_annotations(cml, selection_df=sel_df)
            create_tcr_sequence_pdf(
                assembled,
                output_dir / "tcr_sequences.pdf",
                annotations=pdf_annotations,
                strict="skip",
                species=run_species,
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
        if not assemble_enabled:
            print(
                "\n[7/7] Skipping whole-repertoire assembly (--no-assemble). "
                "Analysis outputs written; assemble the selected basket with "
                "`tcrsift select` + `tcrsift report selected`."
            )
        else:
            print("\n[7/7] Skipping assembly")
        if config.output.output_airr or config.output.output_fasta:
            print("  Skipping AIRR/FASTA export because assembly was skipped")

    if rules is not None and not rules.empty:
        _selection_base = assembled if assembled is not None else til_matched
        if "CDR3ab" in _selection_base.columns:
            selected = (
                rules.merge(_selection_base, on="CDR3ab", how="left")
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
        create_pipeline_funnel(
            raw_cells=funnel_counts.get("Raw Cells", 0),
            with_vdj=funnel_counts.get("With VDJ", funnel_counts.get("Raw Cells", 0)),
            phenotyped=funnel_counts.get("Phenotyped", 0),
            clonotypes=funnel_counts.get("Clonotypes", 0),
            filtered=funnel_counts.get("Filtered", 0),
            tier_counts=tier_counts,
            output_dir=plots_dir,
            selected_count=n_selected_clones,
            emit_selected_variant=n_selected_clones is not None,
            non_viral=n_non_viral,
        )

        if long_df is not None and "method" in getattr(long_df, "columns", []):
            from .plots import plot_vgene_usage_by_method

            plot_vgene_usage_by_method(
                long_df, clonotypes, plots_dir / "vgene_usage_by_method.png",
            )

        if assembled is not None and any(
            f"{c}_leader_aa" in getattr(assembled, "columns", []) for c in ("alpha", "beta")
        ):
            from .plots import plot_leader_summary

            plot_leader_summary(assembled, plots_dir / "leader_summary.png")

    if config.output.generate_report:
        if not list(plots_dir.glob("*.png")):
            print("\nSkipping figure report generation (no plot PNGs produced)")
        else:
            print("\nGenerating report...")
            from .report import resolve_report_name

            resolved_report_name = resolve_report_name(
                output_dir,
                clones_df=clonotypes,
                cli_name=report_name,
            )
            report_format = str(config.output.report_format).lower()
            if report_format == "both":
                generate_report(plots_dir, format="html", report_name=resolved_report_name)
                generate_report(plots_dir, format="pdf", report_name=resolved_report_name)
            elif report_format in {"html", "pdf"}:
                generate_report(plots_dir, format=report_format, report_name=resolved_report_name)
            else:
                raise TCRsiftValidationError(
                    f"Invalid report format: '{config.output.report_format}'",
                    hint="Use report_format: pdf, html, or both in config.",
                )

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"Results saved to: {output_dir}")
    print("=" * 60)

    return PipelineResult(
        output_dir=output_dir,
        data_dir=data_dir,
        plots_dir=plots_dir,
        species=run_species,
        clonotypes=clonotypes,
        filtered=filtered,
        annotated=annotated,
        til_matched=til_matched,
        assembled=assembled,
        stimulation_response=stim_resp_df,
        funnel_counts=funnel_counts,
    )
