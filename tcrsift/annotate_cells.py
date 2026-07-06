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

"""Generalized per-cell / per-cluster cell-type annotator (#312).

Supersedes the CD4/CD8-ratio-only :func:`tcrsift.phenotype.classify_tcell_type`
(which stays a thin back-compat shim): every scTCR project needs typing beyond
CD4/CD8 — myeloid, B/plasma, NK, DC, stroma, and T sub-states. Given a Leiden
clustering (from :mod:`tcrsift.embed`), each cluster is labeled by a
background-subtracted signature argmax over the shared
:data:`tcrsift.signatures.CELL_TYPE_SIGNATURES` registry, with biology-aware
gates so shared master-TF / activation genes can't win a specific label:

- **Treg requires FOXP3** (not IL2RA/CTLA4); **Th1/Th2/Th17** are gated on
  IFNG/IL13/IL17A and CD4-only; an unpolarized T cluster gets an honest
  ``activated`` / ``(unspecified)`` catch-all.
- **γδ T** is called on the δ-constant **TRDC only** (never TRGC1/2, which leak
  into αβ T cells); **MAIT** on **SLC4A10 + KLRB1**.
- **NK vs cytotoxic CD8** is resolved on NK-restricted receptors
  (KLRF1/NCR1/NCAM1) absent from T, gated on CD3.
- **Mast** requires granule proteases (already encoded in the registry, not KIT).

Tumor / cancer-testis-antigen typing is out of scope (oncoref, #310).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .signatures import (
    B_STATE_SIGNATURES,
    CELL_TYPE_SIGNATURES,
    LINEAGE_GENES,
    PROTECTED_SUBTYPES,
    T_STATE_SIGNATURES,
)

logger = logging.getLogger(__name__)

# Genes that must never become a cluster's distinguishing suffix: clonal TCR/Ig
# V/D/J/C segments (clonotype, not cell type).
_TOP_GENE_EXCLUDE = re.compile(r"^(TR[ABGD][VDJ]|IG[HKL][VDJ]|IGKC|IGLC|IGLL|JCHAIN)")


@dataclass
class AnnotationGates:
    """Biology-aware acceptance gates for :func:`annotate_clusters` (#312).

    All thresholds are cluster-mean **log1p CP10K** expression. Defaults match
    the validated reference; ``PhenotypeConfig`` fields override them.
    """

    # T-states that need their defining EFFECTOR gene above a floor (the master
    # TF / activation genes they share with any stimulated T cell are not
    # enough): name -> (gene, min_cluster_mean).
    state_defining_genes: dict = field(default_factory=lambda: {
        "Treg": ("FOXP3", 0.2), "Th1": ("IFNG", 0.2),
        "Th2": ("IL13", 0.2), "Th17": ("IL17A", 0.1),
    })
    helper_cd4_states: frozenset = frozenset({"Th1", "Th2", "Th17"})
    activation_genes: tuple = ("IL2RA", "TNFRSF9", "TNFRSF4", "CD69")
    activation_min: float = 0.5
    # γδ T is called ONLY on the δ-constant gene (never TRGC1/2).
    gd_tcr_genes: tuple = ("TRDC",)
    gd_tcr_min: float = 0.5
    mait_gene: str = "SLC4A10"
    mait_min: float = 0.3
    mait_klrb1_min: float = 1.0
    # Cell-type argmax gates.
    other_threshold: float = 0.04
    nk_specific_min: float = 0.3
    cd3_dominant: float = 1.0


DEFAULT_GATES = AnnotationGates()


def gates_from_phenotype_config(config) -> AnnotationGates:
    """Build :class:`AnnotationGates` from a :class:`tcrsift.config.PhenotypeConfig`
    (#312), so the annotator's floors are configured through the one phenotype
    section rather than a parallel config."""
    return AnnotationGates(
        state_defining_genes={
            "Treg": ("FOXP3", config.treg_foxp3_min),
            "Th1": ("IFNG", config.th1_ifng_min),
            "Th2": ("IL13", config.th2_il13_min),
            "Th17": ("IL17A", config.th17_il17a_min),
        },
        gd_tcr_min=config.gd_tcr_min,
        mait_min=config.mait_min,
        activation_min=config.activation_min,
    )


def _cluster_gene_mean(adata, gene: str, leiden_col: str) -> dict:
    """Per-cluster mean expression of one gene (empty dict if absent)."""
    if gene not in adata.var_names:
        return {}
    v = adata[:, gene].X
    v = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
    return (
        pd.Series(v, index=adata.obs.index)
        .groupby(adata.obs[leiden_col], observed=True)
        .mean()
        .to_dict()
    )


def _block_mean(adata, genes: list[str]) -> np.ndarray:
    """Per-cell mean expression over a gene block (zeros if none present)."""
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs)
    X = adata[:, present].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return X.mean(axis=1)


def score_reference(
    adata,
    *,
    cell_types: dict | None = None,
    t_states: dict | None = None,
    b_states: dict | None = None,
    lineage_genes: dict | None = None,
) -> None:
    """Score every cell-type / T-state / B-state signature and lineage mean into
    ``adata.obs`` (columns ``ct::`` / ``st::`` / ``bt::`` / ``lin::``).

    Signatures use scanpy ``score_genes`` (background-corrected, so a broadly
    expressed set doesn't always win the argmax); lineage means are plain
    log1p-CP10K gene-block means (for the absolute-threshold gates). Expects
    ``adata.X`` log-normalized (see :func:`annotate_cells`, which handles it).
    Gene sets come from :mod:`tcrsift.signatures` — defined once, shared with
    the per-clone/per-cell scorers.
    """
    import scanpy as sc

    cell_types = cell_types if cell_types is not None else CELL_TYPE_SIGNATURES
    t_states = t_states if t_states is not None else T_STATE_SIGNATURES
    b_states = b_states if b_states is not None else B_STATE_SIGNATURES
    lineage_genes = lineage_genes if lineage_genes is not None else LINEAGE_GENES

    present = set(adata.var_names)
    ctrl = min(50, max(1, adata.n_vars - 1))
    for prefix, registry in (("ct", cell_types), ("st", t_states), ("bt", b_states)):
        for name, genes in registry.items():
            col = f"{prefix}::{name}"
            g = [x for x in genes if x in present]
            if g:
                sc.tl.score_genes(adata, g, score_name=col, ctrl_size=ctrl)
            else:
                adata.obs[col] = np.nan
    for lin, genes in lineage_genes.items():
        adata.obs[f"lin::{lin}"] = _block_mean(adata, genes)


def annotate_clusters(
    adata,
    *,
    leiden_col: str = "leiden",
    reference: dict | None = None,
    allowed_types=None,
    gates: AnnotationGates = DEFAULT_GATES,
) -> dict:
    """Label each Leiden cluster: cell-type signature argmax → label; T clusters
    get CD8/CD4 + a gated dominant state; weakly-typed clusters → ``"Other"``.

    ``reference`` is the cell-type registry (default
    :data:`tcrsift.signatures.CELL_TYPE_SIGNATURES`); ``allowed_types`` restricts
    the argmax (e.g. :data:`tcrsift.signatures.PBMC_CULTURE_TYPES` for a culture,
    so a moDC/macrophage cluster can't win a stromal signature it shares an
    activation program with). Requires :func:`score_reference` to have run.
    Returns ``{cluster -> label}``.
    """
    reference = reference if reference is not None else CELL_TYPE_SIGNATURES
    types = list(reference) if allowed_types is None else [
        t for t in reference if t in allowed_types
    ]
    ct_cols = [f"ct::{t}" for t in types if f"ct::{t}" in adata.obs.columns]
    st_cols = [f"st::{s}" for s in T_STATE_SIGNATURES if f"st::{s}" in adata.obs.columns]
    bt_cols = [f"bt::{s}" for s in B_STATE_SIGNATURES if f"bt::{s}" in adata.obs.columns]
    grp = adata.obs.groupby(leiden_col, observed=True)
    ctm = grp[ct_cols].mean()
    stm = grp[st_cols].mean() if st_cols else pd.DataFrame(index=ctm.index)
    btm = grp[bt_cols].mean() if bt_cols else pd.DataFrame(index=ctm.index)
    lin_cols = [c for c in ("lin::CD8", "lin::CD4", "lin::CD3", "lin::NKspec")
                if c in adata.obs.columns]
    linm = grp[lin_cols].mean()

    gate_means = {
        gene: _cluster_gene_mean(adata, gene, leiden_col)
        for gene, _ in gates.state_defining_genes.values()
    }
    act_each = [_cluster_gene_mean(adata, g, leiden_col) for g in gates.activation_genes]
    act_each = [m for m in act_each if m]
    act_mean = (
        {cl: float(np.mean([m.get(cl, 0.0) for m in act_each])) for cl in ctm.index}
        if act_each else {}
    )
    gd_means = {g: _cluster_gene_mean(adata, g, leiden_col) for g in gates.gd_tcr_genes}
    mait_mean = _cluster_gene_mean(adata, gates.mait_gene, leiden_col)
    klrb1_mean = _cluster_gene_mean(adata, "KLRB1", leiden_col)
    has_nk = "ct::NK cell" in ctm.columns

    labels: dict = {}
    for cl in ctm.index:
        scores = ctm.loc[cl]
        best_v = scores.max()
        if pd.isna(best_v):  # every cell-type signature absent → untyped
            labels[cl] = "Other"
            continue
        best = scores.idxmax().split("::", 1)[1]

        # NK vs cytotoxic-CD8: resolve on the POSITIVE NK signal (NK-restricted
        # receptors CD8 lacks), gated so a strong-CD3 cluster stays T.
        if (has_nk and best in ("T cell", "NK cell")
                and ctm.loc[cl, "ct::NK cell"] > gates.other_threshold):
            is_nk = (
                linm.loc[cl, "lin::NKspec"] > gates.nk_specific_min
                and linm.loc[cl, "lin::CD3"] < gates.cd3_dominant
            )
            best = "NK cell" if is_nk else "T cell"

        # Innate-like T override (before Other/state): a TRDC+ (γδ) or
        # SLC4A10+/CD161+ (MAIT) lymphoid cluster is that T lineage even if the
        # NK gate grabbed it on shared GNLY. γδ is called ONLY on TRDC.
        if best in ("T cell", "NK cell"):
            if max((gd_means[g].get(cl, 0.0) for g in gates.gd_tcr_genes),
                   default=0.0) >= gates.gd_tcr_min:
                labels[cl] = "gd T"
                continue
            if (mait_mean.get(cl, 0.0) >= gates.mait_min
                    and klrb1_mean.get(cl, 0.0) >= gates.mait_klrb1_min):
                labels[cl] = "MAIT"
                continue

        if best_v < gates.other_threshold:
            labels[cl] = "Other"
        elif best == "T cell":
            lineage = ("CD8" if linm.loc[cl, "lin::CD8"] >= linm.loc[cl, "lin::CD4"]
                       else "CD4")
            ranked = stm.loc[cl].sort_values(ascending=False) if st_cols else pd.Series(dtype=float)
            state = None
            for scol in ranked.index:
                cand = scol.split("::", 1)[1]
                if cand in gates.helper_cd4_states and lineage != "CD4":
                    continue  # Th1/2/17 are CD4-only
                gate = gates.state_defining_genes.get(cand)
                if gate and gate_means.get(gate[0], {}).get(cl, 0.0) < gate[1]:
                    continue  # gated subset without its defining effector gene
                if ranked[scol] <= 0:
                    break  # no positively-supported state remains
                state = cand
                break
            if state is None:  # confounded / unpolarized → honest catch-all
                state = ("activated" if act_mean.get(cl, 0.0) >= gates.activation_min
                         else "(unspecified)")
            labels[cl] = f"{lineage} {state}"
        elif best == "B cell":
            bstate = btm.loc[cl] if bt_cols else pd.Series(dtype=float)
            labels[cl] = (bstate.idxmax().split("::", 1)[1]
                          if len(bstate) and bstate.max() > 0.0 else "B cell")
        else:
            labels[cl] = best
    return labels


def top_markers(adata, leiden_col: str = "leiden", *, n_genes: int = 30) -> dict:
    """Per-cluster top distinguishing gene (for the phenotype-label suffix).

    Ranks each cluster's DE genes (scanpy ``rank_genes_groups``, Wilcoxon) and
    returns the top one that is NOT a clonal TCR/Ig segment. Empty when a
    cluster can't be ranked.
    """
    import scanpy as sc

    out: dict = {}
    if adata.obs[leiden_col].nunique() < 2:
        return out
    try:
        sc.tl.rank_genes_groups(
            adata, leiden_col, method="wilcoxon", n_genes=n_genes,
        )
    except (ValueError, ZeroDivisionError):
        return out
    names = adata.uns["rank_genes_groups"]["names"]
    for cl in names.dtype.names:
        for gene in names[cl]:
            if gene and not _TOP_GENE_EXCLUDE.match(str(gene)):
                out[cl] = str(gene)
                break
    return out


def compose_phenotype_labels(
    base_by_cluster: dict,
    top_by_cluster: dict,
    size_by_cluster: dict,
    *,
    min_subtype_frac: float = 0.01,
    protect=PROTECTED_SUBTYPES,
) -> dict:
    """Compose ``"Type · GENE"`` labels, dropping the distinguishing-gene suffix
    for a cluster that is a negligible fraction (< ``min_subtype_frac``) of its
    cell type — a <1% same-type fragment shouldn't read as a distinct subtype.
    Such clusters collapse to the bare type label. Types in ``protect`` (e.g.
    pDC/cDC1) always keep their suffix; a type with a single cluster keeps it
    (it is 100% of the type).
    """
    type_total: dict = {}
    for cl, base in base_by_cluster.items():
        type_total[base] = type_total.get(base, 0) + size_by_cluster.get(cl, 0)
    labels: dict = {}
    for cl, base in base_by_cluster.items():
        gene = top_by_cluster.get(cl)
        total = type_total.get(base, 0) or 1
        frac = size_by_cluster.get(cl, 0) / total
        if gene and (base in protect or frac >= min_subtype_frac):
            labels[cl] = f"{base} · {gene}"
        else:
            labels[cl] = base
    return labels


def annotate_cells(
    adata,
    *,
    leiden_col: str = "leiden",
    reference: dict | None = None,
    allowed_types=None,
    gates: AnnotationGates = DEFAULT_GATES,
    normalize: bool = True,
    min_subtype_frac: float = 0.01,
    add_suffix: bool = True,
):
    """One-call annotator: score → per-cluster type/state → composed labels (#312).

    Writes ``adata.obs["phenotype_base"]`` (the per-cluster type/state, e.g.
    ``"CD8 effector/cytotoxic"``) and, when ``add_suffix``,
    ``adata.obs["phenotype"]`` (with the distinguishing-gene suffix for
    non-negligible sub-clusters). ``normalize=True`` scores on a log1p-CP10K
    copy (from the ``counts`` layer if present) so the embedded raw-count
    AnnData can be passed straight through. Returns ``adata``.
    """
    if leiden_col not in adata.obs.columns:
        raise ValueError(
            f"annotate_cells: {leiden_col!r} not in obs — run embed_cells first"
        )
    scored = adata
    if normalize:
        import scanpy as sc

        scored = adata.copy()
        if "counts" in scored.layers:
            scored.X = scored.layers["counts"].copy()
        sc.pp.normalize_total(scored, target_sum=1e4)
        sc.pp.log1p(scored)

    score_reference(scored)
    labels = annotate_clusters(
        scored, leiden_col=leiden_col, reference=reference,
        allowed_types=allowed_types, gates=gates,
    )
    sizes = scored.obs[leiden_col].value_counts().to_dict()
    adata.obs["phenotype_base"] = (
        adata.obs[leiden_col].map(labels).astype("object")
    )
    if add_suffix:
        tops = top_markers(scored, leiden_col)
        composed = compose_phenotype_labels(
            labels, tops, sizes, min_subtype_frac=min_subtype_frac,
        )
        adata.obs["phenotype"] = adata.obs[leiden_col].map(composed).astype("object")
    else:
        adata.obs["phenotype"] = adata.obs["phenotype_base"]
    logger.info(
        "annotate_cells: labeled %d clusters over %d cells",
        len(labels), adata.n_obs,
    )
    return adata
