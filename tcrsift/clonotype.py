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
Clonotype aggregation for TCRsift.

Groups cells by TCR CDR3 sequences to identify clonal populations.
"""

from __future__ import annotations

import json
import logging

import anndata as ad
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from ._dtypes import rehydrate_obs
from .validation import (
    TCRsiftValidationError,
    pick_representative_cell,
    safe_divide,
    validate_anndata,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


def aggregate_clonotypes(
    adata: ad.AnnData,
    group_by: str = "CDR3ab",
    min_umi: int = 2,
    handle_doublets: str = "flag",
    verbose: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Aggregate cells into clonotypes based on CDR3 sequences.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with VDJ and phenotype information
    group_by : str
        How to group clones: "CDR3ab" (alpha+beta) or "CDR3b_only" (beta only)
    min_umi : int
        Minimum UMI count for a chain to be considered
    handle_doublets : str
        How to handle cells with multiple chains: "flag", "remove", "keep-primary"
    verbose : bool
        Print detailed progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        DataFrame with one row per unique clonotype
    """
    # Validate inputs
    adata = validate_anndata(adata, "input AnnData", min_cells=1)
    validate_numeric_param(min_umi, "min_umi", min_value=0)

    # Re-pin obs dtypes — h5ad round-trips can return string columns as
    # Categorical, bool as object, etc., which breaks downstream `.fillna`
    # and string concat. Idempotent.
    rehydrate_obs(adata)

    valid_group_by = ["CDR3ab", "CDR3b_only"]
    if group_by not in valid_group_by:
        raise TCRsiftValidationError(
            f"Invalid group_by: '{group_by}'",
            hint=f"Valid options are: {valid_group_by}",
        )

    valid_doublet_handling = ["flag", "remove", "keep-primary"]
    if handle_doublets not in valid_doublet_handling:
        raise TCRsiftValidationError(
            f"Invalid handle_doublets: '{handle_doublets}'",
            hint=f"Valid options are: {valid_doublet_handling}",
        )

    # Check for required columns
    if group_by == "CDR3ab":
        required = ["CDR3_alpha", "CDR3_beta"]
    else:
        required = ["CDR3_beta"]

    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        available = [c for c in adata.obs.columns if "CDR3" in c or "cdr3" in c.lower()]
        raise TCRsiftValidationError(
            f"Missing required CDR3 columns for {group_by} grouping: {missing}",
            hint=f"Available CDR3-related columns: {available}. "
            "Make sure VDJ data was loaded correctly.",
        )

    logger.info(f"Aggregating clonotypes by {group_by} from {len(adata):,} cells")

    df = adata.obs.copy()

    # Handle doublets
    if "multi_chain" in df.columns:
        n_doublets = df["multi_chain"].sum()
        if n_doublets > 0:
            if verbose:
                logger.info(
                    f"  Found {n_doublets:,} cells with multiple chains ({n_doublets / len(df) * 100:.1f}%)"
                )

            if handle_doublets == "remove":
                df = df[~df["multi_chain"]]
                if verbose:
                    logger.info(f"  Removed doublets, {len(df):,} cells remaining")
            elif handle_doublets == "flag":
                # Keep but flag
                df["is_doublet"] = df["multi_chain"]
                if verbose:
                    logger.info("  Flagging doublets (keeping all cells)")
            # keep-primary: use primary chain (already sorted by UMI)
            elif verbose:
                logger.info("  Using primary chain for doublets")

    # Apply UMI filter
    if "TRA_1_umis" in df.columns and min_umi > 0:
        df["TRA_pass_umi"] = df["TRA_1_umis"].fillna(0) >= min_umi
        df["TRB_pass_umi"] = df["TRB_1_umis"].fillna(0) >= min_umi
    else:
        df["TRA_pass_umi"] = True
        df["TRB_pass_umi"] = True

    # Build clone identifier. CDR3_alpha/CDR3_beta have been re-pinned to
    # object by rehydrate_obs above, so .fillna("") is safe.
    if group_by == "CDR3ab":
        df["CDR3ab"] = df["CDR3_alpha"].fillna("") + "_" + df["CDR3_beta"].fillna("")
        df["is_complete_clone"] = (
            df["CDR3_alpha"].notna()
            & (df["CDR3_alpha"] != "")
            & df["CDR3_beta"].notna()
            & (df["CDR3_beta"] != "")
            & df["TRA_pass_umi"]
            & df["TRB_pass_umi"]
        )
    elif group_by == "CDR3b_only":
        df["CDR3ab"] = df["CDR3_beta"].fillna("")
        df["is_complete_clone"] = (
            df["CDR3_beta"].notna() & (df["CDR3_beta"] != "") & df["TRB_pass_umi"]
        )
    else:
        raise ValueError(f"Invalid group_by: {group_by}. Use 'CDR3ab' or 'CDR3b_only'")

    # Filter to complete clones
    df_complete = df[df["is_complete_clone"]].copy()
    if verbose:
        logger.info(
            f"  Found {len(df_complete):,} cells with complete TCR ({len(df_complete) / len(df) * 100:.1f}%)"
        )

    if len(df_complete) == 0:
        raise TCRsiftValidationError(
            "No complete clones found after filtering",
            hint=f"Check that cells have valid CDR3 sequences. "
            f"Grouping by: {group_by}, min_umi: {min_umi}. "
            f"Total cells: {len(df)}, cells with CDR3_beta: {df['CDR3_beta'].notna().sum() if 'CDR3_beta' in df.columns else 'N/A'}",
        )

    # Count unique clones for progress bar
    n_unique_clones = df_complete["CDR3ab"].nunique()
    if verbose:
        logger.info(f"  Aggregating {n_unique_clones:,} unique clone IDs...")

    # Aggregate by clone
    clonotypes = _aggregate_clone_data(df_complete, group_by, show_progress=show_progress)

    if verbose:
        logger.info(f"  Found {len(clonotypes):,} unique clonotypes")
        # Log clone size distribution
        if len(clonotypes) > 0:
            n_singletons = (clonotypes["cell_count"] == 1).sum()
            n_expanded = (clonotypes["cell_count"] > 1).sum()
            max_size = clonotypes["cell_count"].max()
            logger.info(
                f"    Singletons: {n_singletons:,}, Expanded clones: {n_expanded:,}, Max clone size: {max_size:,}"
            )

    return clonotypes


def _aggregate_clone_data(
    df: pd.DataFrame, group_by: str, show_progress: bool = True
) -> pd.DataFrame:
    """Aggregate cell-level data to clone-level."""

    clone_data = []

    # Per-method total cell counts (denominator for max_frequency_per_method).
    # Computed once outside the per-clone loop. Only meaningful when
    # enrichment_method was populated in the sample sheet (#8).
    method_cell_totals: pd.Series | None = None
    if "enrichment_method" in df.columns:
        method_cell_totals = df["enrichment_method"].value_counts()

    # Create iterator with optional progress bar
    grouped = df.groupby("CDR3ab", observed=True)
    if show_progress:
        grouped = tqdm(
            grouped,
            desc="Aggregating clones",
            unit="clone",
            total=df["CDR3ab"].nunique(),
        )

    for cdr3ab, clone_df in grouped:
        if cdr3ab == "_" or cdr3ab == "":
            continue

        record = {
            "CDR3ab": cdr3ab,
            "cell_count": len(clone_df),
            "cell_barcodes": ";".join(clone_df.index.tolist()),
        }

        # CDR3 sequences
        if group_by == "CDR3ab":
            # Split on first underscore only to handle edge cases
            parts = cdr3ab.split("_", 1)
            record["CDR3_alpha"] = parts[0] if len(parts) > 0 else ""
            record["CDR3_beta"] = parts[1] if len(parts) > 1 else ""
        else:
            # For CDR3b_only grouping, pull the α-chain identity from
            # the cell with the strongest α-UMI evidence (#94). Using
            # per-column ``most_common`` here would risk decoupling
            # CDR3_alpha from V/J/AA/NT picked elsewhere if a future
            # change ever reads more α columns at this site.
            record["CDR3_beta"] = cdr3ab
            if "CDR3_alpha" in clone_df.columns:
                rep_alpha = pick_representative_cell(
                    clone_df, umi_cols=("TRA_1_umis",)
                )
                record["CDR3_alpha"] = (
                    rep_alpha.get("CDR3_alpha", "")
                    if rep_alpha is not None else ""
                )
            else:
                record["CDR3_alpha"] = ""

        # Sample and condition information
        record["samples"] = ";".join(clone_df["sample"].unique())
        record["n_samples"] = clone_df["sample"].nunique()
        record["n_conditions"] = record["n_samples"]

        # Donor / enrichment-method aggregations (#8). Only populated when
        # the corresponding sample-sheet fields were supplied.
        has_donor = "patient_id" in clone_df.columns
        has_method = "enrichment_method" in clone_df.columns

        if has_donor:
            donors = sorted(clone_df["patient_id"].dropna().astype(str).unique())
            record["donors"] = ";".join(donors)
            record["n_donors"] = len(donors)

        if has_method:
            methods = sorted(
                clone_df["enrichment_method"].dropna().astype(str).unique()
            )
            record["methods"] = ";".join(methods)
            record["n_methods"] = len(methods)

        if has_donor and has_method:
            # methods_per_donor is the dict-shape view from #8/#9:
            # {donor: sorted list of distinct methods seen for this clone in
            # that donor}. Stored as a JSON string for clean h5ad/CSV
            # serialization. max_methods_per_donor is the derived scalar
            # used by the planned --min-methods-per-donor filter.
            methods_per_donor: dict[str, list[str]] = {}
            grouped_dm = clone_df.dropna(subset=["patient_id", "enrichment_method"])
            for donor, group in grouped_dm.groupby("patient_id", observed=True):
                donor_methods = sorted(
                    group["enrichment_method"].astype(str).unique()
                )
                if donor_methods:
                    methods_per_donor[str(donor)] = donor_methods
            record["methods_per_donor"] = json.dumps(
                methods_per_donor, sort_keys=True
            )
            record["max_methods_per_donor"] = max(
                (len(m) for m in methods_per_donor.values()), default=0
            )

        # Per-method cell counts and frequencies (#15 chunk 3).
        # max_cells_per_method backs --min-cells-per-method;
        # max_frequency_per_method backs --min-frequency-per-method.
        if has_method:
            method_cell_counts = clone_df["enrichment_method"].value_counts()
            record["max_cells_per_method"] = (
                int(method_cell_counts.max()) if len(method_cell_counts) else 0
            )
            if method_cell_totals is not None and len(method_cell_counts):
                # Align by method, divide. Methods present in this clone are
                # guaranteed to exist in method_cell_totals (same source df).
                freqs = method_cell_counts / method_cell_totals.reindex(
                    method_cell_counts.index
                )
                record["max_frequency_per_method"] = float(freqs.max())
            else:
                record["max_frequency_per_method"] = 0.0

        # Timepoint / APC / tissue axis aggregations (#9 chunk 2). Each block
        # only fires when the corresponding sample-sheet field was supplied
        # — same conditional-emission convention as the donor/method block
        # above so designs that don't use these axes don't get NaN columns.
        has_timepoint = "timepoint" in clone_df.columns
        has_apc = "apc_type" in clone_df.columns
        has_tissue = "tissue" in clone_df.columns

        if has_timepoint:
            timepoints = sorted(
                clone_df["timepoint"].dropna().astype(str).unique()
            )
            record["timepoints"] = ";".join(timepoints)
            record["n_timepoints"] = len(timepoints)

        if has_apc:
            apcs = sorted(clone_df["apc_type"].dropna().astype(str).unique())
            record["apcs"] = ";".join(apcs)
            record["n_apcs"] = len(apcs)

        if has_tissue:
            tissues = sorted(clone_df["tissue"].dropna().astype(str).unique())
            record["tissues"] = ";".join(tissues)
            record["n_tissues"] = len(tissues)

        # Nested timepoint-per-donor aggregation (mirrors methods_per_donor).
        # Backs the planned --min-timepoints-per-donor filter knob.
        if has_donor and has_timepoint:
            timepoints_per_donor: dict[str, list[str]] = {}
            grouped_td = clone_df.dropna(subset=["patient_id", "timepoint"])
            for donor, group in grouped_td.groupby("patient_id", observed=True):
                donor_tps = sorted(group["timepoint"].astype(str).unique())
                if donor_tps:
                    timepoints_per_donor[str(donor)] = donor_tps
            record["timepoints_per_donor"] = json.dumps(
                timepoints_per_donor, sort_keys=True
            )
            record["max_timepoints_per_donor"] = max(
                (len(t) for t in timepoints_per_donor.values()), default=0
            )

        # Nested apc-per-donor aggregation. Backs --min-apcs-per-donor.
        if has_donor and has_apc:
            apcs_per_donor: dict[str, list[str]] = {}
            grouped_ad = clone_df.dropna(subset=["patient_id", "apc_type"])
            for donor, group in grouped_ad.groupby("patient_id", observed=True):
                donor_apcs = sorted(group["apc_type"].astype(str).unique())
                if donor_apcs:
                    apcs_per_donor[str(donor)] = donor_apcs
            record["apcs_per_donor"] = json.dumps(
                apcs_per_donor, sort_keys=True
            )
            record["max_apcs_per_donor"] = max(
                (len(a) for a in apcs_per_donor.values()), default=0
            )

        # Antigen information if available
        if "antigen_description" in clone_df.columns or "antigen_name" in clone_df.columns:
            antigen_desc = (
                clone_df["antigen_description"]
                if "antigen_description" in clone_df.columns
                else pd.Series(index=clone_df.index, dtype=object)
            )
            antigen_name = (
                clone_df["antigen_name"]
                if "antigen_name" in clone_df.columns
                else pd.Series(index=clone_df.index, dtype=object)
            )
            antigens_series = antigen_desc.where(
                antigen_desc.notna() & (antigen_desc != ""), antigen_name
            )
            antigens = antigens_series.dropna().unique()
            record["antigens"] = ";".join(str(a) for a in antigens)
            record["n_antigens"] = len(antigens)
            record["n_conditions"] = record["n_antigens"]

        # Source information
        if "source" in clone_df.columns:
            record["sources"] = ";".join(clone_df["source"].unique())

        # T cell type consensus
        if "Tcell_type" in clone_df.columns:
            type_counts = clone_df["Tcell_type"].value_counts()
            if len(type_counts) > 0:
                record["Tcell_type_consensus"] = type_counts.index[0]
                record["Tcell_type_purity"] = safe_divide(
                    type_counts.iloc[0], len(clone_df), default=0.0
                )
            else:
                record["Tcell_type_consensus"] = "Unknown"
                record["Tcell_type_purity"] = 0.0

            # CD4/CD8 counts
            record["n_CD8"] = clone_df["is_CD8"].sum() if "is_CD8" in clone_df.columns else 0
            record["n_CD4"] = clone_df["is_CD4"].sum() if "is_CD4" in clone_df.columns else 0

        # Gene calls + VDJ AA/NT (#94). All coupled per-chain columns
        # must come from a SINGLE representative cell, not from
        # independent per-column modes. Independent modes break ties
        # alphabetically per column, so AA and NT can end up sourced
        # from different cells and produce internally inconsistent
        # (AA, NT) pairs (translate(NT) != AA). Per-chain
        # representatives are picked by per-chain UMI rank so each
        # chain's data comes from the cell with the most read support
        # for THAT chain — even when α and β were captured in
        # different cells.
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            rep = pick_representative_cell(
                clone_df, umi_cols=(f"{prefix}_umis",)
            )
            for gene in ["v_gene", "j_gene", "c_gene"]:
                col = f"{prefix}_{gene}"
                if col in clone_df.columns:
                    record[f"{chain}_{gene}"] = (
                        rep.get(col) if rep is not None else None
                    )
            vdj_col = f"{prefix}_vdj_aa"
            if vdj_col in clone_df.columns:
                record[f"VDJ_{chain}_aa"] = (
                    rep.get(vdj_col) if rep is not None else None
                )
            vdj_nt_col = f"{prefix}_vdj_nt"
            if vdj_nt_col in clone_df.columns:
                record[f"VDJ_{chain}_nt"] = (
                    rep.get(vdj_nt_col) if rep is not None else None
                )

        # Quality metrics
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            for metric in ["umis", "reads"]:
                col = f"{prefix}_{metric}"
                if col in clone_df.columns:
                    record[f"{chain}_{metric}_mean"] = clone_df[col].mean()
                    record[f"{chain}_{metric}_sum"] = clone_df[col].sum()

        # Contig IDs
        for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
            contig_col = f"{prefix}_contig_id"
            if contig_col in clone_df.columns:
                record[f"{chain}_contig_ids"] = ";".join(
                    clone_df[contig_col].dropna().astype(str).tolist()
                )

        # Doublet information
        if "is_doublet" in clone_df.columns:
            record["n_doublet_cells"] = clone_df["is_doublet"].sum()

        # Calculate clone frequency within each sample
        # Note: Uses outer 'df' (not df_complete) to get total count of complete
        # clones per sample, which is the correct denominator for frequency
        sample_freqs = []
        for sample in clone_df["sample"].unique():
            sample_total = df[df["sample"] == sample]["is_complete_clone"].sum()
            sample_count = (clone_df["sample"] == sample).sum()
            freq = safe_divide(sample_count, sample_total, default=0.0)
            sample_freqs.append(freq)
        record["max_frequency"] = max(sample_freqs) if sample_freqs else 0
        record["mean_frequency"] = np.mean(sample_freqs) if sample_freqs else 0

        clone_data.append(record)

    return pd.DataFrame(clone_data)


def calculate_clone_frequencies(
    clonotypes: pd.DataFrame,
    adata: ad.AnnData,
) -> pd.DataFrame:
    """
    Calculate detailed frequency information for each clone.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame from aggregate_clonotypes
    adata : ad.AnnData
        Original AnnData with cell-level data

    Returns
    -------
    pd.DataFrame
        Clonotypes with additional frequency columns
    """
    df = adata.obs.copy()

    # Calculate total complete TCRs per sample
    sample_totals = df.groupby(
        "sample", observed=True,
    )["is_complete_clone"].sum().to_dict()

    freq_data = []
    for _, clone_row in clonotypes.iterrows():
        cdr3ab = clone_row["CDR3ab"]
        clone_cells = df[df["CDR3ab"] == cdr3ab]

        sample_freqs = {}
        for sample in clone_cells["sample"].unique():
            sample_count = (clone_cells["sample"] == sample).sum()
            sample_total = sample_totals.get(sample, 0)
            sample_freqs[sample] = safe_divide(sample_count, sample_total, default=0.0)

        freq_data.append(
            {
                "CDR3ab": cdr3ab,
                "sample_frequencies": sample_freqs,
                "max_frequency": max(sample_freqs.values()) if sample_freqs else 0,
                "n_conditions_present": len(sample_freqs),
            }
        )

    freq_df = pd.DataFrame(freq_data)

    # Merge back
    clonotypes = clonotypes.merge(
        freq_df[["CDR3ab", "sample_frequencies", "n_conditions_present"]],
        on="CDR3ab",
        how="left",
        suffixes=("", "_new"),
    )

    return clonotypes


def build_clone_sample_long(adata: ad.AnnData) -> pd.DataFrame:
    """Build a long-format (clone, sample) DataFrame from adata.obs.

    One row per (CDR3ab, sample) pair where the clone has at least one
    cell in that sample. Includes ``donor`` and ``method`` columns when
    ``patient_id`` and ``enrichment_method`` are populated; same for
    ``timepoint`` and ``apc_type``. UMI sums (``n_alpha_umis`` /
    ``n_beta_umis``) are emitted when the per-chain UMI columns are
    present.

    Both the per-clone cell count (numerator) and the frequency
    denominator are restricted to complete-clone cells (both CDR3s
    present and both chains >=2 UMI), mirroring the convention used by
    ``max_frequency`` in ``aggregate_clonotypes`` so the two are directly
    comparable. Consequently per-sample frequencies sum to 1.0 and
    single-chain clones do not appear in the table (#175).

    Implements #20 chunk 1 — surfaces the per-(clone, sample) view that
    users were previously reconstructing by parsing the semicolon-
    delimited ``samples`` string on ``clonotypes.csv``.
    """
    rehydrate_obs(adata)
    obs = adata.obs.copy()

    if "sample" not in obs.columns:
        raise TCRsiftValidationError(
            "adata.obs missing 'sample' column",
            hint="Run load_samples first.",
        )

    # Build CDR3ab if not already present (aggregate_clonotypes works on a
    # copy and doesn't write CDR3ab back to adata.obs).
    if "CDR3ab" not in obs.columns:
        if "CDR3_alpha" in obs.columns and "CDR3_beta" in obs.columns:
            obs["CDR3ab"] = (
                obs["CDR3_alpha"].fillna("") + "_" + obs["CDR3_beta"].fillna("")
            )
        elif "CDR3_beta" in obs.columns:
            obs["CDR3ab"] = obs["CDR3_beta"].fillna("")
        else:
            raise TCRsiftValidationError(
                "adata.obs missing CDR3 columns; cannot build long table",
                hint="Make sure VDJ data was loaded.",
            )

    # Filter to rows with a real CDR3ab. NaN != "" evaluates to True
    # in pandas, so a plain ``!= ""`` check would let chainless cells
    # through and silently inflate the denominator (#80).
    cdr3ab_str = obs["CDR3ab"].astype(str)
    valid = obs[
        obs["CDR3ab"].notna()
        & (cdr3ab_str != "")
        & (cdr3ab_str != "_")
        & (cdr3ab_str != "nan")
    ].copy()

    # Frequency denominator: re-derive ``is_complete_clone`` with the
    # same UMI-gated rule ``aggregate_clonotypes`` uses (#80). The
    # column already present on ``adata.obs`` from an earlier pipeline
    # step often lacks the UMI gate, which made the denominator
    # diverge from ``clonotypes.csv:max_frequency``. Doing it locally
    # here keeps the two paths in lockstep.
    def _has(col):
        return col in valid.columns

    tra_pass_umi = (
        valid["TRA_1_umis"].fillna(0).astype(float) >= 2
        if _has("TRA_1_umis") else pd.Series(True, index=valid.index)
    )
    trb_pass_umi = (
        valid["TRB_1_umis"].fillna(0).astype(float) >= 2
        if _has("TRB_1_umis") else pd.Series(True, index=valid.index)
    )
    if _has("CDR3_alpha") and _has("CDR3_beta"):
        is_complete = (
            valid["CDR3_alpha"].notna()
            & (valid["CDR3_alpha"].astype(str) != "")
            & valid["CDR3_beta"].notna()
            & (valid["CDR3_beta"].astype(str) != "")
            & tra_pass_umi
            & trb_pass_umi
        )
    elif _has("CDR3_beta"):
        is_complete = (
            valid["CDR3_beta"].notna()
            & (valid["CDR3_beta"].astype(str) != "")
            & trb_pass_umi
        )
    else:
        # No chain columns at all — fall back to the raw row count.
        is_complete = pd.Series(True, index=valid.index)
    valid["is_complete_clone"] = is_complete
    complete = valid[is_complete]
    sample_totals = complete.groupby("sample", observed=True).size()

    # Cell counts and UMI sums per (clone, sample). The numerator must use the
    # SAME complete-cell rule as ``sample_totals`` (the denominator). Grouping
    # over ``valid`` instead admitted single-chain clones and both-chain cells
    # failing the UMI>=2 gate — cells the denominator drops — so per-sample
    # frequencies summed to >1 and un-synthesizable single-chain "clones" leaked
    # into the table (#175). Restricting to ``complete`` restores the invariant
    # (per-sample freqs sum to 1.0) and drops the single-chain rows.
    #
    # Keep the (CDR3ab, sample) MultiIndex while assigning the UMI sums so pandas
    # aligns them by key rather than by position, then flatten to columns.
    grouped = complete.groupby(["CDR3ab", "sample"], observed=True)
    out = grouped.size().rename("cells").to_frame()

    if "TRA_1_umis" in complete.columns:
        out["n_alpha_umis"] = grouped["TRA_1_umis"].sum().fillna(0).astype(int)
    if "TRB_1_umis" in complete.columns:
        out["n_beta_umis"] = grouped["TRB_1_umis"].sum().fillna(0).astype(int)

    out = out.reset_index()

    out["frequency"] = out.apply(
        lambda r: safe_divide(r["cells"], sample_totals.get(r["sample"], 0), default=0.0),
        axis=1,
    )

    # Propagate per-sample axis metadata (one value per sample).
    axis_to_short = {
        "patient_id": "donor",
        "enrichment_method": "method",
        "timepoint": "timepoint",
        "apc_type": "apc",
    }
    sample_meta_cols = [c for c in axis_to_short if c in valid.columns]
    if sample_meta_cols:
        sample_meta = valid.groupby("sample", observed=True)[sample_meta_cols].first()
        for col in sample_meta_cols:
            out[axis_to_short[col]] = out["sample"].map(sample_meta[col])

    # Column order matches the issue spec, with optional axes inserted.
    leading = ["CDR3ab", "sample"]
    for col in ("donor", "method", "timepoint", "apc"):
        if col in out.columns:
            leading.append(col)
    leading.extend(["cells", "frequency"])
    for col in ("n_alpha_umis", "n_beta_umis"):
        if col in out.columns:
            leading.append(col)
    extra = [c for c in out.columns if c not in leading]
    out = out[leading + extra]

    return out.sort_values(["CDR3ab", "sample"]).reset_index(drop=True)


def build_clone_method_long(
    long_df: pd.DataFrame,
    *,
    method_col: str = "method",
    clone_col: str = "CDR3ab",
    cells_col: str = "cells",
    freq_col: str = "frequency",
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Per-(clone, method) cell-count + frequency table (#81).

    Aggregates ``long_df`` (one row per (clone, sample)) into one row
    per (clone, method) using the same denominator the clone_sample
    table uses. This is the authoritative input for per-method picks,
    selection-route heatmaps, per-method tier-eligibility queries.

    Columns emitted:

    - ``CDR3ab``, ``method``
    - ``cells_in_method`` — sum of cells across samples within the method
    - ``max_freq_in_method`` — max per-sample frequency observed in the method
    - ``n_samples_in_method`` — count of samples in which the clone appears

    Downstream code was rolling this aggregation by hand in multiple
    places, drifting on ``cells=max`` vs ``cells=sum``. This function
    is the single source of truth (see #81).
    """
    if long_df.empty:
        return pd.DataFrame(
            columns=[
                clone_col, method_col, "cells_in_method",
                "max_freq_in_method", "n_samples_in_method",
            ]
        )
    if method_col not in long_df.columns:
        raise ValueError(
            f"build_clone_method_long: missing {method_col!r} column "
            "(populate enrichment_method on the sample sheet)"
        )

    agg = (
        long_df.groupby([clone_col, method_col], observed=True)
        .agg(
            cells_in_method=(cells_col, "sum"),
            max_freq_in_method=(freq_col, "max"),
            n_samples_in_method=(sample_col, "nunique"),
        )
        .reset_index()
    )
    return agg.sort_values([clone_col, method_col]).reset_index(drop=True)


def compute_sample_overlap_matrices(
    long_df: pd.DataFrame,
    *,
    clone_col: str = "CDR3ab",
    sample_col: str = "sample",
    cells_col: str = "cells",
    restrict_clones: set | None = None,
) -> dict[str, pd.DataFrame]:
    """Sample × sample overlap matrices (#82).

    Returns a dict with two square DataFrames indexed by sample:

    - ``"jaccard"`` — clone-set Jaccard. |A ∩ B| / |A ∪ B| over the
      sets of clones observed in each sample.
    - ``"cell_weighted_jaccard"`` — cell-weighted Jaccard.
      Σ min(cells_A[c], cells_B[c]) / Σ max(cells_A[c], cells_B[c])
      over clones in A ∪ B. Captures abundance similarity, not just
      presence/absence.

    Parameters
    ----------
    long_df : pd.DataFrame
        Output of :func:`build_clone_sample_long`.
    clone_col, sample_col, cells_col : str
        Column names.
    restrict_clones : set | None
        When set, restrict the matrices to clones in this set. Useful
        for "selected clones only" / "tier1 only" variants.

    Returns
    -------
    dict[str, pd.DataFrame]
        ``{"jaccard": …, "cell_weighted_jaccard": …}``. Diagonal is
        always 1.0.
    """
    sub = long_df
    if restrict_clones is not None:
        sub = sub[sub[clone_col].isin(restrict_clones)]

    samples = sorted(long_df[sample_col].astype(str).unique())
    if not samples or sub.empty:
        empty = pd.DataFrame(0.0, index=samples, columns=samples)
        for s in samples:
            empty.loc[s, s] = 1.0
        return {"jaccard": empty, "cell_weighted_jaccard": empty}

    # Pivot to (sample × clone) matrix of cell counts.
    pivot = (
        sub.pivot_table(
            index=sample_col, columns=clone_col,
            values=cells_col, aggfunc="sum", fill_value=0,
        )
        .reindex(samples, fill_value=0)
    )
    sets = {s: set(pivot.columns[pivot.loc[s] > 0]) for s in samples}

    jaccard = pd.DataFrame(0.0, index=samples, columns=samples)
    cw_jaccard = pd.DataFrame(0.0, index=samples, columns=samples)
    for a in samples:
        a_vec = pivot.loc[a].astype(float).to_numpy()
        for b in samples:
            if a == b:
                jaccard.loc[a, b] = 1.0
                cw_jaccard.loc[a, b] = 1.0
                continue
            inter = sets[a] & sets[b]
            union = sets[a] | sets[b]
            jaccard.loc[a, b] = (
                len(inter) / len(union) if union else 0.0
            )
            b_vec = pivot.loc[b].astype(float).to_numpy()
            min_sum = float(np.minimum(a_vec, b_vec).sum())
            max_sum = float(np.maximum(a_vec, b_vec).sum())
            cw_jaccard.loc[a, b] = min_sum / max_sum if max_sum > 0 else 0.0

    return {"jaccard": jaccard, "cell_weighted_jaccard": cw_jaccard}


def build_per_method_rankings(
    filtered: pd.DataFrame,
    long_df: pd.DataFrame,
    top_n: int = 100,
) -> dict[tuple[str, str], pd.DataFrame]:
    """Build per-(donor, method) ranked clone tables.

    For each populated ``(donor, method)`` pair in ``long_df``, returns up
    to ``top_n`` clones (those that survived filtering) ranked by their
    max within-(donor, method) frequency descending. Each row is annotated
    with the clone's tier (when ``filter_mode='fdr'``) and a derived
    ``sharing`` label ("private" / "public" based on ``n_donors``) so
    output CSVs are self-describing regardless of which filter mode
    produced ``filtered``.

    Returns ``{}`` when the long table doesn't carry a ``method`` axis —
    designs that don't supply ``enrichment_method`` get no per-method
    output. Implements #20 chunk 2.
    """
    # The long table uses the short-name 'donor' / 'method' columns from
    # build_clone_sample_long, not the upstream patient_id /
    # enrichment_method.
    if "method" not in long_df.columns:
        return {}
    if "donor" not in long_df.columns:
        # Method axis present but no donor info — fold a synthetic "all"
        # donor key so file naming stays uniform and downstream consumers
        # don't need a special case.
        long_df = long_df.copy()
        long_df["donor"] = "all"

    if filtered is None or len(filtered) == 0:
        return {}

    # Aggregate to one row per (donor, method, CDR3ab): max freq across
    # any sample replicates within that (donor, method) bucket; sum cells.
    grouped = (
        long_df.groupby(["donor", "method", "CDR3ab"], as_index=False, observed=True)
        .agg(cells=("cells", "sum"), frequency=("frequency", "max"))
    )

    # Pull useful annotation columns from `filtered`.
    annot_cols = ["CDR3ab"]
    for col in ("tier", "max_frequency", "cell_count", "n_donors", "n_methods"):
        if col in filtered.columns:
            annot_cols.append(col)
    annot = filtered[annot_cols].drop_duplicates("CDR3ab").copy()

    # Derive a private/public sharing label when n_donors is on the
    # filtered table; mode-agnostic.
    if "n_donors" in annot.columns:
        annot["sharing"] = (
            annot["n_donors"].fillna(1).astype(int).apply(
                lambda n: "public" if n >= 2 else "private"
            )
        )

    grouped = grouped.merge(annot, on="CDR3ab", how="inner")

    rankings: dict[tuple[str, str], pd.DataFrame] = {}
    for (donor, method), group in grouped.groupby(["donor", "method"], observed=True):
        ranked = (
            group.sort_values("frequency", ascending=False)
            .head(top_n)
            .reset_index(drop=True)
        )
        # Drop donor/method from the row schema since they're encoded in
        # the file name.
        ranked = ranked.drop(columns=["donor", "method"])
        rankings[(str(donor), str(method))] = ranked

    return rankings


def build_method_overlap_matrices(
    filtered: pd.DataFrame,
    long_df: pd.DataFrame,
    similarity: str = "jaccard",
) -> dict[str, pd.DataFrame]:
    """Per-donor method × method overlap matrices over filter-passing clones.

    For each donor that has at least two distinct methods in ``long_df``,
    builds a symmetric (n_methods × n_methods) DataFrame whose off-diagonal
    cells are pairwise overlap of selected-clone sets and whose diagonal is
    the per-method count.

    ``similarity`` ∈ {``"jaccard"``, ``"dice"``, ``"count"``}:
      - ``jaccard``: ``|A ∩ B| / |A ∪ B|`` — defaults to 0 when both empty.
        Diagonal = 1.0.
      - ``dice``:    ``2|A ∩ B| / (|A| + |B|)`` — defaults to 0 when both
        empty. Diagonal = 1.0.
      - ``count``:   raw intersection counts. Diagonal = per-method clone
        count for clarity.

    Returns ``{}`` when ``long_df`` lacks the method axis or ``filtered``
    is empty. Implements #27 chunk 3.
    """
    if similarity not in ("jaccard", "dice", "count"):
        raise TCRsiftValidationError(
            f"Invalid similarity: {similarity!r}",
            hint="Valid options: 'jaccard', 'dice', 'count'.",
        )

    if "method" not in long_df.columns:
        return {}
    if filtered is None or len(filtered) == 0:
        return {}

    # Restrict long_df to filter-passing clones.
    selected = long_df[long_df["CDR3ab"].isin(set(filtered["CDR3ab"]))]
    if "donor" not in selected.columns:
        # Fold under a synthetic 'all' donor so output naming stays uniform.
        selected = selected.assign(donor="all")

    matrices: dict[str, pd.DataFrame] = {}
    for donor, donor_df in selected.groupby("donor", observed=True):
        method_to_clones: dict[str, set[str]] = {
            str(m): set(g["CDR3ab"].unique())
            for m, g in donor_df.groupby("method", observed=True)
        }
        methods = sorted(method_to_clones)
        n = len(methods)
        if n == 0:
            continue

        dtype = int if similarity == "count" else float
        mat = pd.DataFrame(
            np.zeros((n, n), dtype=dtype), index=methods, columns=methods
        )
        for i, mi in enumerate(methods):
            ci = method_to_clones[mi]
            for j, mj in enumerate(methods):
                cj = method_to_clones[mj]
                if i == j:
                    mat.iat[i, j] = len(ci) if similarity == "count" else 1.0
                    continue
                inter = len(ci & cj)
                if similarity == "count":
                    mat.iat[i, j] = inter
                elif similarity == "jaccard":
                    union = len(ci | cj)
                    mat.iat[i, j] = inter / union if union else 0.0
                else:  # dice
                    denom = len(ci) + len(cj)
                    mat.iat[i, j] = (2 * inter) / denom if denom else 0.0
        matrices[str(donor)] = mat
    return matrices


def build_method_recovery_table(
    filtered: pd.DataFrame,
    long_df: pd.DataFrame,
    tier: str = "tier1",
) -> pd.DataFrame:
    """Per-(donor, method) recovery of ``tier``-level filtered clones.

    Returns a long DataFrame ``[donor, method, recovered, total, fraction]``
    where ``total`` is the number of clones in ``filtered`` carrying
    ``tier`` for that donor (or all of ``filtered`` when ``tier=='*'``),
    and ``recovered`` is how many of those clones appear in ``long_df`` for
    that ``(donor, method)`` bucket.

    Implements #27 chunk 4 — backs the method-recovery report panel.
    """
    if "method" not in long_df.columns or filtered is None or len(filtered) == 0:
        return pd.DataFrame(
            columns=["donor", "method", "recovered", "total", "fraction"]
        )
    if "donor" not in long_df.columns:
        long_df = long_df.assign(donor="all")

    if tier == "*":
        target = filtered
    elif "tier" in filtered.columns:
        target = filtered[filtered["tier"] == tier]
    else:
        # No tier annotation; fall back to the full filtered set (so the
        # panel still works under non-FDR filter modes).
        target = filtered

    if len(target) == 0:
        return pd.DataFrame(
            columns=["donor", "method", "recovered", "total", "fraction"]
        )

    target_clones = set(target["CDR3ab"])

    # Restrict long_df to those clones; group by (donor, method).
    selected = long_df[long_df["CDR3ab"].isin(target_clones)]
    rows = []
    for donor in sorted(long_df["donor"].dropna().astype(str).unique()):
        # Total per donor: target clones that appear at all in this donor.
        # If a donor has zero target clones, total is 0 and fraction is 0 —
        # don't fall back to the cohort-wide target count, which would
        # produce misleading "0/N" denominators that imply the donor was
        # supposed to recover a target it never had access to.
        donor_clones_total = set(
            long_df[long_df["donor"].astype(str) == donor]["CDR3ab"]
        )
        donor_target = target_clones & donor_clones_total
        donor_total = len(donor_target)
        for method in sorted(
            long_df["method"].dropna().astype(str).unique()
        ):
            sub = selected[
                (selected["donor"].astype(str) == donor)
                & (selected["method"].astype(str) == method)
            ]
            recovered = sub["CDR3ab"].nunique()
            rows.append(
                {
                    "donor": donor,
                    "method": method,
                    "recovered": int(recovered),
                    "total": int(donor_total),
                    "fraction": (
                        recovered / donor_total if donor_total else 0.0
                    ),
                }
            )
    return pd.DataFrame(rows)


def get_clonotype_summary(clonotypes: pd.DataFrame) -> dict:
    """
    Get summary statistics for clonotypes.

    Returns
    -------
    dict
        Summary statistics
    """
    return {
        "n_clonotypes": len(clonotypes),
        "n_cells": clonotypes["cell_count"].sum(),
        "median_clone_size": clonotypes["cell_count"].median(),
        "max_clone_size": clonotypes["cell_count"].max(),
        "n_singletons": (clonotypes["cell_count"] == 1).sum(),
        "n_expanded": (clonotypes["cell_count"] > 1).sum(),
        "n_multi_sample": (clonotypes["n_samples"] > 1).sum()
        if "n_samples" in clonotypes.columns
        else 0,
    }


def export_clonotypes_airr(clonotypes: pd.DataFrame, output_path: str):
    """
    Export clonotypes in AIRR format.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame
    output_path : str
        Output file path (.tsv)
    """
    # Map to AIRR standard columns
    airr_mapping = {
        "CDR3ab": "clone_id",
        "CDR3_alpha": "junction_aa_tra",
        "CDR3_beta": "junction_aa_trb",
        "alpha_v_gene": "v_call_tra",
        "alpha_j_gene": "j_call_tra",
        "beta_v_gene": "v_call_trb",
        "beta_d_gene": "d_call_trb",
        "beta_j_gene": "j_call_trb",
        "cell_count": "clone_count",
    }

    airr_df = pd.DataFrame()
    for src_col, dst_col in airr_mapping.items():
        if src_col in clonotypes.columns:
            airr_df[dst_col] = clonotypes[src_col]

    # Add required AIRR fields with defaults
    if "sequence_id" not in airr_df.columns:
        airr_df["sequence_id"] = [f"clone_{i}" for i in range(len(airr_df))]
    if "productive" not in airr_df.columns:
        airr_df["productive"] = "T"

    airr_df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Exported {len(airr_df)} clonotypes to AIRR format: {output_path}")
