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
Data loading functions for TCRsift.

Handles loading CellRanger VDJ and GEX outputs into unified data structures.
"""

from __future__ import annotations

import logging
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from anndata.experimental import concat_on_disk
from scipy.sparse import csr_matrix, issparse
from tqdm.auto import tqdm

from ._dtypes import rehydrate_obs
from .genes import TCELL_MARKERS, adata_mito_mask, find_column_for_gene, resolve_adata_var_map
from .sample_sheet import Sample, SampleSheet, load_sample_sheet
from .species import resolve_species
from .validation import (
    TCRsiftValidationError,
    resolve_cellranger_gex_source,
    resolve_cellranger_vdj_file,
    validate_file_exists,
    validate_numeric_param,
)

logger = logging.getLogger(__name__)


_LIST_METADATA_COLS = frozenset(
    {"antigen_names", "antigen_sequences", "epitope_sequences"}
)


def _add_obs_columns(adata: ad.AnnData, new_cols: dict) -> None:
    """Add/overwrite multiple ``adata.obs`` columns in a single concat.

    Assigning columns one at a time in a loop (``adata.obs[c] = …``) triggers
    pandas' "DataFrame is highly fragmented" PerformanceWarning once obs grows
    past ~100 columns. Building one frame and concatenating it in one shot
    avoids the repeated inserts (and the warnings). Scalar values broadcast to
    ``n_obs``; sample-level list metadata fields are ``;``-joined so they
    broadcast cleanly across cells and serialize to h5ad.
    """
    if not new_cols:
        return
    built: dict = {}
    for key, val in new_cols.items():
        if key in _LIST_METADATA_COLS and isinstance(val, (list, tuple)):
            joined = ";".join(str(x) for x in val if x is not None and str(x) != "")
            built[key] = [joined] * adata.n_obs
        elif np.isscalar(val) or val is None:
            built[key] = [val] * adata.n_obs
        else:
            built[key] = val
    new_df = pd.DataFrame(built, index=adata.obs.index)
    overlap = [c for c in new_df.columns if c in adata.obs.columns]
    if overlap:
        adata.obs = adata.obs.drop(columns=overlap)
    adata.obs = pd.concat([adata.obs, new_df], axis=1)


# VDJ segment columns for full sequence assembly
VDJ_SEGMENT_COLS = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
VDJ_SEGMENT_NT_COLS = [c + "_nt" for c in VDJ_SEGMENT_COLS]


def _sample_prefix_from_contig_file(annotations_path: Path) -> str:
    """Extract sample filename prefix (e.g. ``'GSM1_'``) before ``filtered_contig`` / ``all_contig``."""
    lower = annotations_path.name.lower()
    for token in ("filtered_contig_annotations", "all_contig_annotations"):
        idx = lower.find(token)
        if idx > 0:
            return annotations_path.name[:idx]
    return ""


def _resolve_clonotypes_path(
    vdj_base_dir: Path,
    clonotypes_filename: str,
    *,
    sample_prefix: str = "",
) -> Path | None:
    """Locate an optional CellRanger ``clonotypes.csv(.gz)`` file in ``vdj_base_dir``."""
    if not vdj_base_dir.is_dir():
        return None
    if sample_prefix:
        for name in (
            f"{sample_prefix}{clonotypes_filename}",
            f"{sample_prefix}{clonotypes_filename}.gz",
        ):
            cand = vdj_base_dir / name
            if cand.is_file():
                return cand
        return None
    for name in (clonotypes_filename, f"{clonotypes_filename}.gz"):
        cand = vdj_base_dir / name
        if cand.is_file():
            return cand
    hits = sorted(
        p
        for p in (*vdj_base_dir.glob("*clonotypes.csv"), *vdj_base_dir.glob("*clonotypes.csv.gz"))
        if p.is_file()
    )
    return hits[0] if len(hits) == 1 else None


def _enrich_vdj_from_json(
    df: pd.DataFrame,
    vdj_base_dir: Path,
    *,
    sample_prefix: str = "",
    verbose: bool = True,
) -> pd.DataFrame:
    """Populate ``vdj_aa`` and ``vdj_nt`` from CellRanger ``*contig_annotations*.json(.gz)``
    when pre-CellRanger-6 CSV tables lack ``fwr1..fwr4`` segment columns."""
    if "vdj_aa" in df.columns or "contig_id" not in df.columns or not vdj_base_dir.is_dir():
        return df

    json_candidates: list[Path] = []
    if sample_prefix:
        patterns: tuple[str, ...] = (
            f"{sample_prefix}filtered_contig_annotations.json",
            f"{sample_prefix}filtered_contig_annotations.json.gz",
            f"{sample_prefix}all_contig_annotations.json",
            f"{sample_prefix}all_contig_annotations.json.gz",
        )
    else:
        patterns = (
            "filtered_contig_annotations.json",
            "filtered_contig_annotations.json.gz",
            "all_contig_annotations.json",
            "all_contig_annotations.json.gz",
            "*filtered_contig_annotations.json",
            "*filtered_contig_annotations.json.gz",
            "*all_contig_annotations.json",
            "*all_contig_annotations.json.gz",
        )
    for pat in patterns:
        found = sorted(p for p in vdj_base_dir.glob(pat) if p.is_file())
        if len(found) == 1:
            json_candidates = found
            break
    if not json_candidates:
        return df

    import gzip
    import json

    from .assemble import translate_dna
    from .leaders import GERMLINE_MIN_IDENTITY, germline_vgene_leaders
    from .species import detect_vdj_species

    json_path = json_candidates[0]
    needed_ids = set(df["contig_id"].dropna().astype(str))
    if not needed_ids:
        return df

    sp = detect_vdj_species(df) or "human"
    opener = gzip.open if json_path.name.lower().endswith(".gz") else open
    try:
        with opener(json_path, "rt", encoding="utf-8") as fh:
            records = json.load(fh)
    except Exception as exc:
        logger.warning("Could not parse VDJ JSON annotations %s: %s", json_path, exc)
        return df

    vdj_aa_map: dict[str, str] = {}
    vdj_nt_map: dict[str, str] = {}
    partial_vdj_aa_map: dict[str, str] = {}
    partial_vdj_nt_map: dict[str, str] = {}
    for rec in records:
        cid = rec.get("contig_name")
        if not cid or cid not in needed_ids:
            continue
        seq = str(rec.get("sequence") or "").upper()
        if not seq:
            continue
        anns = rec.get("annotations") or []
        j_ends: list[int] = []
        v_start_explicit = None
        v_cstart_vonly = None
        v_cstart_lv = None
        v_gene = None
        start_pos = rec.get("start_codon_pos")
        has_start = isinstance(start_pos, int) and start_pos >= 0
        for ann in anns:
            feat = ann.get("feature") or {}
            rtype = feat.get("region_type")
            if rtype == "J-REGION":
                cend = ann.get("contig_match_end")
                if isinstance(cend, int) and 0 < cend <= len(seq):
                    if not has_start or cend > start_pos:
                        j_ends.append(cend)
            elif rtype == "V-REGION":
                v_gene = v_gene or feat.get("gene_name")
                cstart = ann.get("contig_match_start")
                amatch = ann.get("annotation_match_start", 0)
                if isinstance(cstart, int) and 0 <= cstart < len(seq):
                    if v_cstart_vonly is None or cstart < v_cstart_vonly:
                        v_cstart_vonly = cstart
                if (
                    has_start
                    and v_start_explicit is None
                    and isinstance(cstart, int)
                    and cstart >= start_pos
                    and amatch == 0
                    and (cstart - start_pos) % 3 == 0
                ):
                    v_start_explicit = cstart
            elif rtype == "L-REGION+V-REGION":
                v_gene = v_gene or feat.get("gene_name")
                cstart = ann.get("contig_match_start")
                if isinstance(cstart, int) and 0 <= cstart < len(seq):
                    if v_cstart_lv is None or cstart < v_cstart_lv:
                        v_cstart_lv = cstart
        v_cstart_any = v_cstart_vonly if v_cstart_vonly is not None else v_cstart_lv
        if not j_ends:
            continue
        j_end = max(j_ends)
        cdr3 = str(rec.get("cdr3") or "")

        # Full-length path when leader start codon is present on the contig.
        if has_start and j_end > start_pos:
            raw_len = j_end - start_pos
            full_nt = seq[start_pos : j_end - (raw_len % 3)]
            if len(full_nt) >= 240:
                full_aa, _ = translate_dna(full_nt)
                if full_aa and "*" not in full_aa:
                    vdj_start_aa: int | None = None
                    if isinstance(v_start_explicit, int) and v_start_explicit >= start_pos:
                        vdj_start_aa = (v_start_explicit - start_pos) // 3
                    elif v_gene:
                        best_score = -1.0
                        for _allele, _func, g_leader in germline_vgene_leaders(v_gene, species=sp):
                            glen = len(g_leader)
                            max_k = min(15, max(0, len(full_aa) - glen - 60))
                            for k in range(max_k + 1):
                                if full_aa[k] != "M":
                                    continue
                                window = full_aa[k : k + glen]
                                if len(window) != glen:
                                    continue
                                ident = sum(a == b for a, b in zip(window, g_leader)) / glen
                                if ident >= GERMLINE_MIN_IDENTITY and ident > best_score:
                                    best_score = ident
                                    vdj_start_aa = k + glen

                    if vdj_start_aa is not None and vdj_start_aa < len(full_aa):
                        cand_aa = full_aa[vdj_start_aa:]
                        cand_nt = full_nt[3 * vdj_start_aa :]
                        if not cdr3 or cdr3 in cand_aa:
                            vdj_aa_map[str(cid)] = cand_aa
                            vdj_nt_map[str(cid)] = cand_nt.lower()
                            continue

        # 5'-partial contig path: extract in-frame 3' V-suffix + CDR3 + FR4 up to J-end.
        if cdr3 and isinstance(v_cstart_any, int) and j_end - v_cstart_any >= 3 * (len(cdr3) + 8):
            best_partial: tuple[tuple[int, int], str, str] | None = None
            for end_trim in (1, 0, 2):
                eff_end = j_end - end_trim
                if eff_end <= v_cstart_any:
                    continue
                for shift in (0, 1, 2):
                    c_start = v_cstart_any + shift
                    if eff_end <= c_start or (eff_end - c_start) % 3 != 0:
                        continue
                    p_nt = seq[c_start:eff_end]
                    p_aa, _ = translate_dna(p_nt)
                    if not p_aa or "*" in p_aa or "X" in p_aa or cdr3 not in p_aa:
                        continue
                    post_cdr3 = p_aa[p_aa.rindex(cdr3) + len(cdr3) :]
                    if not (8 <= len(post_cdr3) <= 12):
                        continue
                    if len(post_cdr3) >= 3 and not (
                        post_cdr3[0] in ("G", "A") and post_cdr3[2] == "G"
                    ):
                        continue
                    rank = (1 if end_trim == 1 else 0, len(p_aa))
                    if best_partial is None or rank > best_partial[0]:
                        best_partial = (rank, p_aa, p_nt.lower())
            if best_partial is not None:
                partial_vdj_aa_map[str(cid)] = best_partial[1]
                partial_vdj_nt_map[str(cid)] = best_partial[2]

    if vdj_aa_map or partial_vdj_aa_map:
        df = df.copy()
        cids = df["contig_id"].astype(str)
        if vdj_aa_map:
            df["vdj_aa"] = cids.map(vdj_aa_map)
            df["vdj_nt"] = cids.map(vdj_nt_map)
            if verbose:
                logger.info(
                    f"  Extracted full-length VDJ sequences for {len(vdj_aa_map):,} contigs from {json_path.name}"
                )
        if partial_vdj_aa_map:
            df["partial_vdj_aa"] = cids.map(partial_vdj_aa_map)
            df["partial_vdj_nt"] = cids.map(partial_vdj_nt_map)
            if verbose:
                logger.info(
                    f"  Extracted 5'-partial VDJ suffixes for {len(partial_vdj_aa_map):,} contigs from {json_path.name}"
                )
    return df


def load_cellranger_vdj(
    vdj_dir: str | Path,
    sample_name: str,
    annotations_filename: str = "filtered_contig_annotations.csv",
    clonotypes_filename: str = "clonotypes.csv",
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Load CellRanger VDJ output files.

    Parameters
    ----------
    vdj_dir : str or Path
        Path to CellRanger VDJ output directory (or contig CSV/.csv.gz file)
    sample_name : str
        Name to assign to this sample
    annotations_filename : str
        Name of the contig annotations CSV file
    clonotypes_filename : str
        Name of the clonotypes CSV file
    verbose : bool
        Print progress information

    Returns
    -------
    pd.DataFrame
        DataFrame with VDJ annotations for all cells
    """
    vdj_dir = Path(vdj_dir)

    try:
        annotations_path, unfiltered = resolve_cellranger_vdj_file(
            vdj_dir, annotations_filename=annotations_filename
        )
    except TCRsiftValidationError as exc:
        if "Multiple " in str(exc):
            raise
        raise TCRsiftValidationError(
            f"Invalid CellRanger VDJ directory for sample '{sample_name}': {vdj_dir}",
            hint="Make sure this is the 'outs' directory from 'cellranger vdj'. "
            "It should contain 'filtered_contig_annotations.csv(.gz)' or 'all_contig_annotations.csv(.gz)'.",
        ) from exc

    vdj_base_dir = annotations_path.parent
    sample_prefix = _sample_prefix_from_contig_file(annotations_path) if vdj_dir.is_file() else ""
    clonotypes_path = _resolve_clonotypes_path(
        vdj_base_dir, clonotypes_filename, sample_prefix=sample_prefix
    )
    if unfiltered and verbose:
        logger.info(f"  Using unfiltered contig table {annotations_path.name}")

    logger.info(f"Loading VDJ annotations from {annotations_path}")
    df = pd.read_csv(annotations_path, low_memory=False)

    # Validate the VDJ data
    if len(df) == 0:
        raise TCRsiftValidationError(
            f"VDJ annotations file is empty: {annotations_path}",
            hint="Check that CellRanger VDJ ran successfully. "
            "The file should contain contig annotations.",
        )

    required_cols = ["barcode", "chain"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise TCRsiftValidationError(
            f"VDJ annotations missing required columns: {missing_cols}",
            hint=f"Available columns: {list(df.columns)[:15]}. "
            "This doesn't look like a CellRanger VDJ output file.",
        )

    # all_contig_annotations.csv also lists contigs from barcodes that were not
    # called as cells and low-confidence contigs. Apply Cell Ranger's own
    # definition of the filtered table so the fallback means the same thing.
    from ._schema import clean_cdr3_series, truthy_flag_series

    if unfiltered:
        keep = pd.Series(True, index=df.index)
        for flag in ("is_cell", "high_confidence"):
            if flag in df.columns:
                keep &= truthy_flag_series(df[flag])
        if verbose and not keep.all():
            logger.info(f"  Kept {int(keep.sum()):,} of {len(df):,} contigs (called cells, high confidence)")
        df = df[keep]
    # A contig with no chain call (read as missing) cannot be assigned to a
    # receptor chain. Left in, it turned the per-chain entry numbers into
    # floats and produced columns such as ``TRA_1.0_cdr3``.
    no_chain = clean_cdr3_series(df, "chain") == ""
    if no_chain.any():
        if verbose:
            logger.info(f"  Dropped {int(no_chain.sum()):,} contigs with no chain call")
        df = df[~no_chain]
    # Cell Ranger's filtered table still lists non-productive contigs (frame
    # shifts, stop codons, no CDR3). Keep explicitly productive contigs, plus
    # 5'-partial contigs where Cell Ranger left ``productive="None"`` because
    # 5' coverage stopped downstream of the leader ATG despite a complete
    # V+J call and a stop-free canonical IMGT C...[FW] CDR3 (common in Cell
    # Ranger 2/3 and targeted TCR libraries such as GSE171037).
    if "productive" in df.columns:
        prod_str = df["productive"].astype(str).str.strip().str.lower()
        productive = truthy_flag_series(df["productive"])
        undetermined = df["productive"].isna() | prod_str.isin({"none", "", "nan", "<na>"})
        if undetermined.any() and {"v_gene", "j_gene", "cdr3"} <= set(df.columns):
            v_ok = clean_cdr3_series(df, "v_gene") != ""
            j_ok = clean_cdr3_series(df, "j_gene") != ""
            cdr3_ok = (
                df["cdr3"]
                .astype("string")
                .str.strip()
                .str.fullmatch(r"C[ACDEFGHIKLMNPQRSTVWY]{3,28}[FW]", na=False)
                .to_numpy(dtype=bool)
            )
            cell_ok = pd.Series(True, index=df.index)
            for flag in ("is_cell", "high_confidence"):
                if flag in df.columns:
                    cell_ok &= truthy_flag_series(df[flag])
            partial_rescue = undetermined & v_ok & j_ok & cdr3_ok & cell_ok
        else:
            partial_rescue = pd.Series(False, index=df.index)
        keep_prod = productive | partial_rescue
        if verbose and not keep_prod.all():
            logger.info(f"  Dropped {int((~keep_prod).sum()):,} non-productive contigs")
        if verbose and partial_rescue.any():
            logger.info(
                f"  Retained {int(partial_rescue.sum()):,} 5'-partial contigs "
                "(productive=None) with complete V+J genes and in-frame CDR3"
            )
        df = df[keep_prod].copy()
        df["_is_full_productive"] = productive.loc[df.index].to_numpy(dtype=bool)
    else:
        df = df.copy()
        df["_is_full_productive"] = True
    if len(df) == 0:
        raise TCRsiftValidationError(
            f"No usable contigs in: {annotations_path}",
            hint="Every contig was from a non-cell barcode, low confidence, non-productive, "
            "or had no chain call.",
        )

    # Log summary statistics
    n_contigs = len(df)
    n_cells = df["barcode"].nunique()
    if "productive" in df.columns:
        n_productive = "all productive"
    else:
        n_productive = "productivity not reported"
    if verbose:
        logger.info(
            f"  Loaded {n_contigs:,} contigs from {n_cells:,} cells ({n_productive})"
        )

    # Validate chain types
    valid_chains = {"TRA", "TRB", "TRD", "TRG", "IGH", "IGK", "IGL", "Multi"}
    invalid_chains = set(df["chain"].unique()) - valid_chains
    if invalid_chains:
        logger.warning(f"  Unexpected chain types found: {invalid_chains}")

    # Add sample information
    df["sample"] = sample_name
    df["vdj_dir"] = str(vdj_base_dir)

    # Load clonotypes if available
    if clonotypes_path is not None and clonotypes_path.exists():
        logger.info(f"Loading clonotypes from {clonotypes_path}")
        df_clonotypes = pd.read_csv(clonotypes_path)
        # Clonotypes contain MAIT/NKT evidence
        if (
            "clonotype_id" in df_clonotypes.columns
            and "raw_clonotype_id" in df.columns
            and (
                "mait_evidence" in df_clonotypes.columns
                or "inkt_evidence" in df_clonotypes.columns
            )
        ):
            # Merge clonotype info
            clonotype_cols = ["clonotype_id"]
            if "mait_evidence" in df_clonotypes.columns:
                clonotype_cols.append("mait_evidence")
            if "inkt_evidence" in df_clonotypes.columns:
                clonotype_cols.append("inkt_evidence")
            df_clonotypes_subset = df_clonotypes[clonotype_cols].copy()
            df_clonotypes_subset = df_clonotypes_subset.rename(
                columns={"clonotype_id": "raw_clonotype_id"}
            )
            # Drop duplicate clonotype ids before the left merge: a malformed /
            # concatenated clonotypes.csv with a repeated id would otherwise
            # fan-out (duplicate) every contig row for that clonotype.
            df_clonotypes_subset = df_clonotypes_subset.drop_duplicates(
                "raw_clonotype_id"
            )
            df = df.merge(df_clonotypes_subset, on="raw_clonotype_id", how="left")

    # Combine VDJ segments into full sequence ONLY on rows where all 7 segments
    # (fwr1..fwr4) are present and stop-free; 5'-partial contigs leave vdj_aa as None
    # rather than emitting an N-terminally truncated fragment.
    if all(col in df.columns for col in VDJ_SEGMENT_COLS):
        seg_complete = pd.Series(True, index=df.index)
        for col in VDJ_SEGMENT_COLS:
            s_val = df[col].astype("string").str.strip()
            seg_complete &= (
                df[col].notna()
                & ~s_val.str.lower().isin({"", "none", "nan", "<na>"})
                & ~s_val.str.contains(r"\*", regex=True, na=True)
            )
        joined_aa = df[VDJ_SEGMENT_COLS].fillna("").agg("".join, axis=1)
        df["vdj_aa"] = joined_aa.where(seg_complete, other=None)
        has_nt_segs = all(col in df.columns for col in VDJ_SEGMENT_NT_COLS)
        if has_nt_segs:
            nt_complete = seg_complete.copy()
            for col in VDJ_SEGMENT_NT_COLS:
                s_nt = df[col].astype("string").str.strip()
                nt_complete &= df[col].notna() & ~s_nt.str.lower().isin({"", "none", "nan", "<na>"})
            joined_nt = df[VDJ_SEGMENT_NT_COLS].fillna("").agg("".join, axis=1)
            df["vdj_nt"] = joined_nt.where(nt_complete, other=None)

        # Retain contiguous stop-free 3' segment suffixes (at least cdr3+fwr4) on
        # 5'-partial contigs so opt-in --germline-vdj-fallback can rescue only the
        # missing N-terminal framework prefix while keeping sequenced donor codons.
        if (~seg_complete).any():
            partial_aa_vals: list[str | None] = [None] * len(df)
            partial_nt_vals: list[str | None] = [None] * len(df)
            has_any_partial = False
            for row_pos, idx in enumerate(df.index):
                if bool(seg_complete.loc[idx]):
                    continue
                row_s = df.loc[idx]
                suffix_cols: list[str] = []
                for col in reversed(VDJ_SEGMENT_COLS):
                    val = row_s[col]
                    sval = str(val).strip() if pd.notna(val) else ""
                    if not sval or sval.lower() in {"none", "nan", "<na>"} or "*" in sval:
                        break
                    suffix_cols.append(col)
                if len(suffix_cols) >= 2:  # at least cdr3 + fwr4
                    suffix_cols.reverse()
                    p_aa = "".join(str(row_s[c]).strip() for c in suffix_cols)
                    partial_aa_vals[row_pos] = p_aa
                    has_any_partial = True
                    if has_nt_segs:
                        nt_cols = [c + "_nt" for c in suffix_cols]
                        nt_parts = [
                            str(row_s[c]).strip() if pd.notna(row_s[c]) else ""
                            for c in nt_cols
                        ]
                        if all(p and p.lower() not in {"none", "nan", "<na>"} for p in nt_parts):
                            partial_nt_vals[row_pos] = "".join(nt_parts).lower()
            if has_any_partial:
                df["partial_vdj_aa"] = partial_aa_vals
                if any(v is not None for v in partial_nt_vals):
                    df["partial_vdj_nt"] = partial_nt_vals
    if "vdj_aa" not in df.columns:
        df = _enrich_vdj_from_json(
            df, vdj_base_dir, sample_prefix=sample_prefix, verbose=verbose
        )

    return df


def load_cellranger_gex(
    gex_dir: str | Path,
    sample_name: str,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
    verbose: bool = True,
) -> ad.AnnData:
    """
    Load CellRanger gene expression output.

    Parameters
    ----------
    gex_dir : str or Path
        Path to CellRanger count output directory
    sample_name : str
        Name to assign to this sample
    min_genes : int
        Minimum genes detected per cell
    max_genes : int
        Maximum genes detected per cell
    min_counts : int
        Minimum UMI counts per cell
    max_counts : int
        Maximum UMI counts per cell
    max_mito_pct : float
        Maximum mitochondrial percentage
    min_mito_pct : float
        Minimum mitochondrial percentage (a FLOOR: cells below this are
        dropped, targeting near-zero-mito empty/ambient droplets). Default 2.0;
        the drop is logged as a warning so it isn't silent (#168). Set 0 to
        disable the floor.
    verbose : bool
        Print progress information

    Returns
    -------
    ad.AnnData
        AnnData object with gene expression data
    """
    # Validate numeric parameters
    validate_numeric_param(min_genes, "min_genes", min_value=0)
    validate_numeric_param(max_genes, "max_genes", min_value=1)
    validate_numeric_param(min_counts, "min_counts", min_value=0)
    validate_numeric_param(max_counts, "max_counts", min_value=1)
    validate_numeric_param(min_mito_pct, "min_mito_pct", min_value=0, max_value=100)
    validate_numeric_param(max_mito_pct, "max_mito_pct", min_value=0, max_value=100)

    if min_genes > max_genes:
        raise TCRsiftValidationError(
            f"min_genes ({min_genes}) cannot be greater than max_genes ({max_genes})",
            hint="Check your QC filter parameters.",
        )
    if min_counts > max_counts:
        raise TCRsiftValidationError(
            f"min_counts ({min_counts}) cannot be greater than max_counts ({max_counts})",
            hint="Check your QC filter parameters.",
        )
    if min_mito_pct > max_mito_pct:
        raise TCRsiftValidationError(
            f"min_mito_pct ({min_mito_pct}) cannot be greater than max_mito_pct ({max_mito_pct})",
            hint="Check your QC filter parameters.",
        )

    gex_dir = Path(gex_dir)

    try:
        source_kind, resolved_source = resolve_cellranger_gex_source(gex_dir)
    except TCRsiftValidationError as exc:
        if "Multiple " in str(exc):
            raise
        raise TCRsiftValidationError(
            f"Invalid CellRanger GEX directory for sample '{sample_name}': {gex_dir}",
            hint="Make sure this is the 'outs' directory from 'cellranger count' (or a 10x .h5 file). "
            "It should contain 'filtered_feature_bc_matrix' or 'filtered_feature_bc_matrix.h5'.",
        ) from exc

    if source_kind == "h5":
        logger.info(f"Loading GEX from h5 file: {resolved_source}")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Variable names are not unique",
                category=UserWarning,
            )
            adata = sc.read_10x_h5(str(resolved_source))
        # Index genes by Ensembl id, as the matrix-directory reader does.
        # Symbols repeat in 10x references, and a repeated name makes
        # every name-based gene lookup fail; an h5 sample combined with a
        # matrix-directory sample also needs the same gene axis.
        if "gene_ids" in adata.var.columns:
            adata.var["gene_symbols"] = adata.var_names.astype(str)
            adata.var_names = adata.var["gene_ids"].astype(str)
            adata.var = adata.var.drop(columns="gene_ids")
            adata.var.index.name = None
    else:
        logger.info(f"Loading GEX from matrix directory: {resolved_source}")
        adata = sc.read_10x_mtx(str(resolved_source), var_names="gene_ids")

    # Validate loaded data
    if adata.n_obs == 0:
        raise TCRsiftValidationError(
            f"Gene expression matrix contains no cells: {gex_dir}",
            hint="Check that CellRanger count ran successfully.",
        )

    if verbose:
        logger.info(f"  Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Add sample information
    adata.obs["sample"] = sample_name
    adata.obs["gex_dir"] = str(resolved_source.parent if gex_dir.is_file() else gex_dir)

    # Calculate QC metrics - detect mitochondrial genes across human (MT-*, ENSG*)
    # and mouse (mt-*, ENSMUSG*) inputs via the shared symbol resolver.
    adata.var["mt"] = adata_mito_mask(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Rename QC columns for consistency
    if "pct_counts_mt" in adata.obs.columns:
        adata.obs["percent_mt"] = adata.obs["pct_counts_mt"]
    if "n_genes_by_counts" in adata.obs.columns:
        adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    if "total_counts" in adata.obs.columns:
        adata.obs["n_counts"] = adata.obs["total_counts"]

    # Add QC filter flags
    adata.obs["filter:min_genes"] = adata.obs["n_genes"] >= min_genes
    adata.obs["filter:max_genes"] = adata.obs["n_genes"] <= max_genes
    adata.obs["filter:min_counts"] = adata.obs["n_counts"] >= min_counts
    adata.obs["filter:max_counts"] = adata.obs["n_counts"] <= max_counts
    adata.obs["filter:min_mito"] = adata.obs["percent_mt"] >= min_mito_pct
    adata.obs["filter:max_mito"] = adata.obs["percent_mt"] <= max_mito_pct

    # The min_mito_pct FLOOR discards *low*-mito cells (near-zero-mito
    # empty/ambient droplets). The default 2.0 floor trims a small tail on most
    # datasets — that's the intended behavior, so report it at info level rather
    # than warning on every run (#168). Only escalate to a warning when the floor
    # culls an unusually large fraction (>10%), which signals it's set too high.
    if min_mito_pct > 0:
        n_below_floor = int((~adata.obs["filter:min_mito"]).sum())
        if n_below_floor:
            frac = n_below_floor / adata.n_obs if adata.n_obs else 0.0
            msg = (
                f"  min_mito_pct={min_mito_pct} FLOOR dropped "
                f"{n_below_floor:,}/{adata.n_obs:,} ({frac * 100:.0f}%) low-mito cells"
            )
            if frac > 0.10:
                logger.warning(
                    msg + " — a >10% cull from the mito floor is unusual; "
                    "set min_mito_pct=0 to disable it."
                )
            elif verbose:
                logger.info(msg)

    adata.obs["filter:pass_qc"] = (
        adata.obs["filter:min_genes"]
        & adata.obs["filter:max_genes"]
        & adata.obs["filter:min_counts"]
        & adata.obs["filter:max_counts"]
        & adata.obs["filter:min_mito"]
        & adata.obs["filter:max_mito"]
    )

    # Apply the QC mask. Prior to issue #39 these flags were advisory only —
    # the column was computed but cells failing QC were not dropped, so the
    # user-facing --min-mito/--max-mito/--min-genes/etc. parameters were
    # silently no-ops.
    n_before = adata.n_obs
    adata = adata[adata.obs["filter:pass_qc"]].copy()
    if verbose:
        logger.info(
            f"  QC: {n_before:,} -> {adata.n_obs:,} cells pass "
            f"(dropped {n_before - adata.n_obs:,})"
        )
    # Catch "silently wrong" config mistakes (e.g. --min-mito 1.0 typo for
    # 10) where the filter eats most of the sample and downstream artifacts
    # just look weird rather than failing clearly (#41). 50% is a reasonable
    # line — peripheral blood routinely passes 90%+, tumor sits around 70%,
    # below half is almost always a config or dataset-mismatch issue.
    pass_rate = adata.n_obs / n_before if n_before else 0.0
    if pass_rate < 0.5:
        logger.warning(
            f"  QC dropped {(1 - pass_rate) * 100:.0f}% of cells "
            f"({n_before - adata.n_obs:,}/{n_before:,}) — "
            "check --min-genes/--max-genes/--min-counts/--max-counts/"
            "--min-mito/--max-mito if this is unexpected."
        )

    return adata


def _extract_tcell_markers(adata: ad.AnnData) -> pd.DataFrame:
    """
    Extract T cell marker gene expression from AnnData.

    Returns DataFrame with CD3, CD4, CD8 expression per cell.
    Uses version-robust ENSEMBL ID / symbol / mouse ortholog matching from genes.py.

    Resolves all present markers via a single column slice on `adata`
    rather than N per-marker slices — one anndata-view creation and one
    sparse-to-dense conversion instead of N (issue #35, finding #2).
    """
    var_names = list(adata.var_names)

    present_cols: list[str] = []
    present_symbols: list[str] = []
    unresolved_genes = []
    for gene in TCELL_MARKERS:
        col = find_column_for_gene(gene, var_names)
        if col is not None:
            present_cols.append(col)
            present_symbols.append(gene.symbol)
        else:
            unresolved_genes.append(gene)

    missing_symbols: list[str] = []
    if unresolved_genes:
        _, sym_to_var = resolve_adata_var_map(adata, [g.symbol for g in unresolved_genes])
        for gene in unresolved_genes:
            col = sym_to_var.get(gene.symbol)
            if col is not None:
                present_cols.append(col)
                present_symbols.append(gene.symbol)
            else:
                missing_symbols.append(gene.symbol)

    for sym in missing_symbols:
        logger.warning(f"T cell marker {sym} not found in gene expression data")

    if present_cols:
        block = adata[:, present_cols].X
        if hasattr(block, "toarray"):
            block = block.toarray()
        # asarray is defensive: covers the h5py-dataset path under backed
        # mode where .X is neither sparse nor an ndarray.
        block = np.asarray(block)
        marker_df = pd.DataFrame(
            block, index=adata.obs_names, columns=present_symbols
        )
        # Drop the local handle regardless of whether pandas shares or
        # copies its backing. If shared, harmless (marker_df still holds
        # the array); if copied (e.g. later re-blocked when we add
        # missing-marker columns below), this frees the source so peak
        # memory stays at ~n_obs*n_present, not 2x.
        del block
    else:
        marker_df = pd.DataFrame(index=adata.obs_names)
    for sym in missing_symbols:
        marker_df[sym] = 0

    return marker_df


def combine_gex_and_vdj(
    adata: ad.AnnData,
    vdj_df: pd.DataFrame,
) -> ad.AnnData:
    """
    Combine gene expression and VDJ data for a single sample.

    Parameters
    ----------
    adata : ad.AnnData
        Gene expression data
    vdj_df : pd.DataFrame
        VDJ annotations

    Returns
    -------
    ad.AnnData
        Combined AnnData with VDJ info in obs
    """
    # Extract T cell markers
    marker_df = _extract_tcell_markers(adata)
    for col in marker_df.columns:
        adata.obs[col] = marker_df[col].values

    # Calculate combined markers
    adata.obs["CD3"] = adata.obs["CD3D"] + adata.obs["CD3E"] + adata.obs["CD3G"]
    adata.obs["CD8"] = adata.obs["CD8A"] + adata.obs["CD8B"]

    if "sample" not in adata.obs.columns and "sample" in vdj_df.columns:
        uniq_samples = vdj_df["sample"].dropna().astype(str).unique()
        if len(uniq_samples) == 1:
            adata.obs["sample"] = uniq_samples[0]

    # Pivot VDJ data to get one row per barcode with chain info
    if len(vdj_df) > 0:
        vdj_pivoted = _pivot_vdj_by_barcode(vdj_df)

        # Build barcode mapping (handles gem group suffixes like -1, -2, -3)
        # Map each GEX barcode to matching VDJ barcode without modifying adata
        vdj_barcodes = set(vdj_pivoted.index)

        # Helper to strip gem group suffix (e.g., -1, -2, -3)
        def _strip_suffix(bc: str) -> str:
            return bc.rsplit("-", 1)[0] if "-" in bc else bc

        # Start with direct matches where possible
        barcode_to_vdj = {bc: bc for bc in adata.obs_names if bc in vdj_barcodes}

        # Suffix-insensitive matching is for a sample whose VDJ and GEX
        # libraries were given different GEM-group suffixes, so it applies
        # only when no barcode matches exactly. Applied per barcode it gave
        # the receptor of cell X-1 to an unrelated cell X-2 of another GEM
        # group. A base shared by several barcodes on either side is
        # ambiguous and is left unmatched.
        if not barcode_to_vdj:
            vdj_by_base: dict[str, list[str]] = {}
            for bc in sorted(vdj_barcodes):
                vdj_by_base.setdefault(_strip_suffix(bc), []).append(bc)
            gex_base_counts = pd.Series([_strip_suffix(bc) for bc in adata.obs_names]).value_counts()
            for bc in adata.obs_names:
                base = _strip_suffix(bc)
                candidates = vdj_by_base.get(base, [])
                if len(candidates) == 1 and gex_base_counts[base] == 1:
                    barcode_to_vdj[bc] = candidates[0]
        if not barcode_to_vdj and adata.n_obs:
            logger.warning(
                "  No GEX barcode matches any VDJ barcode "
                f"(GEX e.g. {adata.obs_names[0]!r}, VDJ e.g. {min(vdj_barcodes)!r}); "
                "every cell will have no receptor. Check that the two directories are the same sample."
            )

        # Add VDJ columns to adata.obs using the mapping. Build all columns,
        # then assign in one concat to avoid the fragmentation PerformanceWarning
        # from ~91 single-column inserts per donor.
        mapped_barcodes = [barcode_to_vdj.get(bc) for bc in adata.obs_names]
        new_cols: dict = {}
        for col in vdj_pivoted.columns:
            src = vdj_pivoted[col]
            reindexed = src.reindex(mapped_barcodes)
            values = reindexed.values
            # reindex upcasts bool/int to object when missing barcodes introduce NaN,
            # which then breaks h5ad serialization. Restore the original dtype with
            # the natural "absent" fill value.
            missing = pd.isna(values)
            if src.dtype == bool:
                values = np.where(missing, False, values).astype(bool)
            elif pd.api.types.is_integer_dtype(src.dtype):
                values = np.where(missing, 0, values).astype(src.dtype)
            new_cols[col] = values
        _add_obs_columns(adata, new_cols)

    return adata


def _pivot_vdj_by_barcode(vdj_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot VDJ data to get one row per barcode with TRA/TRB chain info.

    Handles doublets by keeping track of multiple chains per barcode.
    """
    # Sort explicitly productive full-length contigs ahead of 5'-partial rescued
    # contigs, then by UMI and read count.
    sort_cols = ["barcode", "chain"]
    sort_asc = [True, True]
    if "_is_full_productive" in vdj_df.columns:
        sort_cols.append("_is_full_productive")
        sort_asc.append(False)
    sort_cols.extend(["umis", "reads"])
    sort_asc.extend([False, False])
    vdj_df = vdj_df.sort_values(sort_cols, ascending=sort_asc)

    # Create entry ID for each chain per barcode (1, 2, etc.)
    vdj_df["entry_id"] = vdj_df.groupby(
        ["barcode", "chain"], observed=True,
    ).cumcount() + 1

    # Columns to pivot - include all CDR/FWR segments if available
    pivot_cols = ["cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "umis", "reads", "contig_id"]

    # Add VDJ segment columns (fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4)
    # Skip columns already in pivot_cols to avoid duplicates
    for seg_col in VDJ_SEGMENT_COLS:
        if seg_col in vdj_df.columns and seg_col not in pivot_cols:
            pivot_cols.append(seg_col)
    for seg_col in VDJ_SEGMENT_NT_COLS:
        if seg_col in vdj_df.columns and seg_col not in pivot_cols:
            pivot_cols.append(seg_col)

    # Add combined VDJ sequences (including 5'-partial suffixes) and clonotype innate-T evidence
    for extra_col in (
        "vdj_aa",
        "vdj_nt",
        "partial_vdj_aa",
        "partial_vdj_nt",
        "mait_evidence",
        "inkt_evidence",
    ):
        if extra_col in vdj_df.columns and extra_col not in pivot_cols:
            pivot_cols.append(extra_col)

    # Filter to only include first 2 entries per chain (to handle doublets)
    vdj_df = vdj_df[vdj_df["entry_id"] <= 2]

    # Pivot
    pivot_df = vdj_df.pivot_table(
        index="barcode",
        columns=["chain", "entry_id"],
        values=pivot_cols,
        aggfunc="first",
    )

    # Flatten column names
    pivot_df.columns = [f"{chain}_{entry}_{col}" for col, chain, entry in pivot_df.columns]

    # Add chain count and doublet flags
    for chain in ["TRA", "TRB"]:
        umi_cols = [
            c for c in pivot_df.columns if c.startswith(f"{chain}_") and c.endswith("_umis")
        ]
        pivot_df[f"{chain}_count"] = pivot_df[umi_cols].notna().sum(axis=1)
        pivot_df[f"has_{chain}"] = pivot_df[f"{chain}_count"] > 0
        pivot_df[f"multi_{chain}"] = pivot_df[f"{chain}_count"] > 1

    pivot_df["multi_chain"] = pivot_df["multi_TRA"] | pivot_df["multi_TRB"]
    pivot_df["has_both_chains"] = pivot_df["has_TRA"] & pivot_df["has_TRB"]

    # Create combined CDR3ab identifier
    pivot_df["CDR3_alpha"] = pivot_df.get("TRA_1_cdr3", pd.Series(index=pivot_df.index, dtype=object))
    pivot_df["CDR3_beta"] = pivot_df.get("TRB_1_cdr3", pd.Series(index=pivot_df.index, dtype=object))
    pivot_df["CDR3ab"] = pivot_df["CDR3_alpha"].fillna("") + "_" + pivot_df["CDR3_beta"].fillna("")

    return pivot_df


def load_sample(
    sample: Sample,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
    species: str | None = None,
) -> ad.AnnData | None:
    """
    Load all data for a single sample.

    Parameters
    ----------
    sample : Sample
        Sample object with paths and metadata
    min_genes, max_genes, min_counts, max_counts, max_mito_pct, min_mito_pct
        QC filter parameters for GEX data

    Returns
    -------
    ad.AnnData or None
        Combined AnnData with GEX and VDJ data, or None if sample has neither
        gex_dir nor vdj_dir.
    """
    adata = None
    vdj_df = None

    # Load GEX if available (from gex_dir or standalone gex_path)
    gex_source = sample.gex_dir or getattr(sample, "gex_path", None)
    if gex_source:
        adata = load_cellranger_gex(
            gex_source,
            sample.sample,
            min_genes=min_genes,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            max_mito_pct=max_mito_pct,
            min_mito_pct=min_mito_pct,
        )

    # Load VDJ if available
    if sample.vdj_dir:
        vdj_df = load_cellranger_vdj(sample.vdj_dir, sample.sample)

    # Combine or create from VDJ only
    if adata is not None and vdj_df is not None:
        adata = combine_gex_and_vdj(adata, vdj_df)
    elif vdj_df is not None:
        # Create minimal AnnData from VDJ data
        vdj_pivoted = _pivot_vdj_by_barcode(vdj_df)
        adata = ad.AnnData(obs=vdj_pivoted)
        adata.obs["sample"] = sample.sample

    # Add sample metadata. Skip fields the user didn't supply: an all-None
    # object column survives concat but breaks anndata.write_h5ad, and
    # downstream readers already guard with `if col in adata.obs.columns`.
    if adata is not None:
        metadata = [
            ("antigen_type", sample.antigen_type),
            ("antigen_description", sample.antigen_description),
            ("antigen_name", sample.antigen_name),
            ("antigen_sequence", sample.antigen_sequence),
            ("epitope_sequence", sample.epitope_sequence),
            ("mhc_allele", sample.mhc_allele),
            ("antigen_names", sample.antigen_names),
            ("antigen_sequences", sample.antigen_sequences),
            ("epitope_sequences", sample.epitope_sequences),
            ("source", sample.source),
            ("control_sample", getattr(sample, "control_sample", None)),
            ("patient_id", sample.patient_id),
            ("enrichment_method", sample.enrichment_method),
            ("timepoint", sample.timepoint),
            ("apc_type", sample.apc_type),
            ("tissue", sample.tissue),
            ("pre_sorted", getattr(sample, "pre_sorted", None)),
            ("expected_tcell_type", sample.get_expected_tcell_type()),
            ("species", getattr(sample, "species", None)),
        ]
        _add_obs_columns(
            adata, {col: val for col, val in metadata if val is not None}
        )
        effective_sp = (
            (species if species and species != "auto" else None)
            or getattr(sample, "species", None)
            or species
        )
        if effective_sp is not None:
            resolved_sp = resolve_species(
                effective_sp,
                adata=adata if adata.n_vars > 0 else None,
                df=vdj_df,
                context=f"sample '{sample.sample}'",
            )
            adata.uns["species"] = resolved_sp

    return adata


def _ensure_x(adata: ad.AnnData) -> None:
    """Normalize `adata.X` so the spilled h5ad merges via `concat_on_disk`.

    `concat_on_disk` requires every input store to have an X group, and
    along the obs axis only supports CSR sparse or dense. Two cases:

    1. VDJ-only AnnData has no expression matrix at all. Synthesize a
       zero-column CSR placeholder. The placeholder is stripped post-merge
       in `load_samples` when `combined.n_vars == 0`, preserving the
       contract that VDJ-only loads have `X is None`.

    2. Any non-CSR sparse format raises `NotImplementedError: Concat of
       following not supported: [...]` from `concat_on_disk`. Convert to
       CSR. The original report was CellRanger's CSC output (#37), but
       the same shape applies to COO, `csc_array` (scipy 1.11+ sparse
       arrays), etc.
    """
    if adata.X is None:
        adata.X = csr_matrix((adata.n_obs, 0))
    elif issparse(adata.X) and getattr(adata.X, "format", "csr") != "csr":
        adata.X = adata.X.tocsr()


def load_samples(
    sample_sheet_path: str | Path | SampleSheet,
    min_genes: int = 250,
    max_genes: int = 15000,
    min_counts: int = 500,
    max_counts: int = 100000,
    max_mito_pct: float = 8.0,
    min_mito_pct: float = 2.0,
    verbose: bool = True,
    show_progress: bool = True,
    tmpdir: str | Path | None = None,
    species: str | None = None,
) -> ad.AnnData:
    """
    Load all samples from a sample sheet into a single AnnData object.

    Parameters
    ----------
    sample_sheet_path : str or Path or SampleSheet
        Path to sample sheet (CSV or YAML), or a SampleSheet instance.
    min_genes, max_genes, min_counts, max_counts, max_mito_pct, min_mito_pct
        QC filter parameters.
    verbose : bool
        Print detailed progress information.
    show_progress : bool
        Show progress bar.
    tmpdir : str or Path or None
        Parent directory for the spill tempdir (see Notes). Defaults to the
        system temp location ($TMPDIR / /tmp). Pass an explicit disk-backed
        directory when the system temp is on tmpfs.

    Returns
    -------
    ad.AnnData
        Combined AnnData with all samples. `X` is None when every input was
        VDJ-only, otherwise it carries the outer-joined expression matrix.

    Notes
    -----
    Memory: each per-sample AnnData is spilled to a tempfile h5ad after
    load and the in-memory copy is freed, then merged with
    `anndata.experimental.concat_on_disk`, which streams inputs and output
    in bounded chunks (~400 MB default). The merged file is read back into
    memory once. Peak memory ≈ max(one sample, output), versus the prior
    in-memory `ad.concat` peak of ~2 × Σ(samples).

    Temp disk: spilled per-sample h5ads plus the merged output, so ~2× the
    total sparse dataset size. Redirect with `tmpdir=` when $TMPDIR is on
    tmpfs.

    VDJ-only samples: see `_ensure_x` — a zero-column placeholder is
    synthesized to keep the merge path uniform, then stripped after the
    merge so `combined.X` is None for VDJ-only loads.
    """
    # Load sample sheet (path or object)
    if isinstance(sample_sheet_path, SampleSheet):
        sample_sheet = sample_sheet_path
        sample_sheet_label = "<SampleSheet>"
    else:
        sample_sheet_path = validate_file_exists(sample_sheet_path, "sample sheet")
        sample_sheet = load_sample_sheet(sample_sheet_path)
        sample_sheet_label = str(sample_sheet_path)

    if len(sample_sheet) == 0:
        raise TCRsiftValidationError(
            f"Sample sheet is empty: {sample_sheet_label}",
            hint="Add sample entries to the sample sheet.",
        )

    logger.info(f"Loading {len(sample_sheet)} samples from {sample_sheet_label}")

    # Pre-validate all sample paths to fail fast
    validation_errors = []
    for sample in sample_sheet:
        if sample.vdj_dir and not Path(sample.vdj_dir).exists():
            validation_errors.append(
                f"Sample '{sample.sample}': VDJ directory not found: {sample.vdj_dir}"
            )
        if sample.gex_dir and not Path(sample.gex_dir).exists():
            validation_errors.append(
                f"Sample '{sample.sample}': GEX directory not found: {sample.gex_dir}"
            )
        if getattr(sample, "gex_path", None) and not Path(sample.gex_path).exists():
            validation_errors.append(
                f"Sample '{sample.sample}': GEX file not found: {sample.gex_path}"
            )

    if validation_errors:
        raise TCRsiftValidationError(
            f"Sample sheet validation failed with {len(validation_errors)} error(s):\n"
            + "\n".join(f"  - {e}" for e in validation_errors[:5]),
            hint="Check that all paths in the sample sheet are correct and accessible.",
        )

    sample_keys: list[str] = []
    gene_symbols: dict[str, str] = {}
    spill_paths: list[Path] = []
    total_cells = 0

    with tempfile.TemporaryDirectory(prefix="tcrsift_load_", dir=tmpdir) as spill_dir_str:
        spill_dir = Path(spill_dir_str)
        sample_iter = tqdm(
            sample_sheet,
            desc="Loading samples",
            unit="sample",
            disable=not show_progress,
        )

        for sample in sample_iter:
            sample_iter.set_postfix(sample=sample.sample[:20])

            if verbose:
                logger.info(f"Loading sample: {sample.sample}")

            try:
                adata = load_sample(
                    sample,
                    min_genes=min_genes,
                    max_genes=max_genes,
                    min_counts=min_counts,
                    max_counts=max_counts,
                    max_mito_pct=max_mito_pct,
                    min_mito_pct=min_mito_pct,
                    species=species or getattr(sample_sheet, "species", None),
                )
            except TCRsiftValidationError:
                raise
            except Exception as e:
                raise TCRsiftValidationError(
                    f"Failed to load sample '{sample.sample}': {e}",
                    hint="Check that the CellRanger output directories are valid and complete. "
                    f"VDJ: {sample.vdj_dir}, GEX: {sample.gex_dir or sample.gex_path}",
                ) from e

            if adata is None:
                continue
            if adata.n_obs == 0:
                logger.warning(
                    f"  Sample {sample.sample}: 0 cells remain after loading/QC; skipping"
                )
                continue

            _ensure_x(adata)
            rehydrate_obs(adata)
            sample_keys.append(str(sample.sample))
            if "gene_symbols" in adata.var.columns:
                gene_symbols.update(zip(adata.var_names, adata.var["gene_symbols"].astype(str)))
            total_cells += adata.n_obs

            # Number files by the load order of successful samples (skipped
            # samples don't leave gaps), so the spill dir is easy to inspect.
            path = spill_dir / f"sample_{len(spill_paths):04d}.h5ad"
            adata.write_h5ad(path)
            spill_paths.append(path)

            if verbose:
                logger.info(f"  Sample {sample.sample}: {adata.n_obs:,} cells")
            del adata  # release the in-memory copy now that it's on disk

        if not sample_keys:
            raise TCRsiftValidationError(
                "No samples loaded successfully",
                hint="Check that at least one sample has valid VDJ or GEX data.",
            )

        logger.info(f"Concatenating {len(sample_keys)} samples ({total_cells:,} total cells)")
        out_path = spill_dir / "combined.h5ad"
        with tqdm(
            total=1, desc="Concatenating samples", unit="step", disable=not show_progress
        ) as pbar:
            concat_on_disk(
                spill_paths, out_path,
                join="outer", label="sample", keys=sample_keys,
            )
            # Non-unique obs_names are EXPECTED and intentional here: 10x reuses
            # one barcode whitelist across samples, and we deliberately do not
            # uniquify (downstream code indexes positionally with a sample-column
            # mask — see qc.py `gex_vdj_overlap`). So anndata's "Observation names
            # are not unique" notice on read-back is noise, not a problem to fix.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", message="Observation names are not unique",
                    category=UserWarning,
                )
                combined = ad.read_h5ad(out_path)
            pbar.update(1)

        # The concatenation drops per-gene annotations. Matrix-directory inputs
        # index genes by Ensembl id and carry the reference's own symbols in
        # ``var["gene_symbols"]``; restore them so symbol lookups use the
        # sample's reference rather than a different Ensembl release (which
        # renames or lacks some genes and made multi-sample signature scoring
        # fail on symbols the matrices do contain).
        if gene_symbols:
            combined.var["gene_symbols"] = [gene_symbols.get(name, name) for name in combined.var_names]

    # Tempdir is cleaned up here; `combined` is fully in-memory.
    # Inverse of `_ensure_x`: an all-VDJ-only sheet produces n_vars == 0,
    # in which case the synthesized empty placeholders should look like the
    # original no-X state to downstream callers.
    if combined.n_vars == 0:
        combined.X = None

    # Store sample sheet as uns as a JSON string. The previous form,
    # `to_dataframe().to_dict()`, produced `{col: {row_idx: val}}` whose
    # int row-index keys break h5ad serialization. Other structured forms
    # (DataFrame directly, `to_dict(orient="list")`) hit anndata/h5py
    # mixed-None-object-column limitations. JSON is unambiguous about
    # None and stores as a single string. Read back with
    # `pd.read_json(io.StringIO(adata.uns["sample_sheet"]), orient="records")`.
    combined.uns["sample_sheet"] = sample_sheet.to_dataframe().to_json(
        orient="records"
    )
    sheet_species = species or getattr(sample_sheet, "species", None)
    combined.uns["species"] = resolve_species(
        sheet_species or "auto",
        adata=combined if combined.n_vars > 0 else None,
        df=combined.obs,
        warn_mismatch=bool(sheet_species and sheet_species != "auto"),
        context="load_samples",
    )

    logger.info(f"Successfully loaded {combined.n_obs:,} cells from {len(sample_keys)} samples")

    return combined
