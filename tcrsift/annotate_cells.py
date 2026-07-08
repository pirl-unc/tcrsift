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


@dataclass
class MarkerCountOverride:
    """Relabel a cluster on **positive marker-count evidence**, after the argmax (#325).

    The cell-type argmax can mis-call a cluster whose defining biology isn't in
    :data:`CELL_TYPE_SIGNATURES` — the motivating case is a solid-tumor cluster
    that shares collagen with fibroblasts and loses the argmax. Rather than teach
    the annotator tumor biology (that stays oncoref's job, #310), this is a
    generic mechanism: **a cluster is relabeled** ``label`` **when a fraction**
    ``min_cluster_frac`` **of its cells each express** ``>= min_distinct``
    **distinct genes from** ``gene_set``. The gene set is caller-supplied domain
    data (e.g. an oncoref cancer-testis-antigen panel, an ISG set, a cycling
    set) — tcrsift only provides the counting + override.

    ``rescue`` optionally saves a low-coverage cluster whose sparse markers drop
    below ``min_cluster_frac``: ``(score, min_mean, min_frac[, rescue_min_distinct])``
    fires the override when the cluster-mean of ``score`` is ``>= min_mean`` *and*
    at least ``min_frac`` of its cells carry ``>= rescue_min_distinct`` distinct
    markers. ``rescue_min_distinct`` defaults to ``min_distinct`` (so a 3-tuple
    behaves exactly as before), but a **lower** floor lets the rescue gate on a
    BROADER any-marker signal — ``rescue_min_distinct=1`` fires on "``score`` high
    AND ``>= 1`` marker in ``>= min_frac`` of cells", the natural low-coverage
    rescue when sparse dropout leaves few cells clearing the full ``>= min_distinct``
    bar (#339). ``score`` is either an ``adata.obs`` column name (a
    caller-precomputed lineage-TF score) or a gene list (its per-cell block-mean)
    — e.g. osteoblastic TFs rescuing a low-CTA-coverage osteosarcoma cluster.

    The per-cell "# distinct markers from ``gene_set``" is written to
    ``adata.obs[count_col]`` (default ``"n_markers_<label>"``) — reused for the
    marker-load UMAP and for QC (an immune cell carrying a tumor-level antigen
    load is a doublet).
    """

    label: str
    gene_set: tuple
    min_distinct: int = 2
    min_cluster_frac: float = 0.4
    rescue: tuple | None = None
    min_expr: float = 0.0
    count_col: str | None = None

    def resolved_count_col(self) -> str:
        if self.count_col:
            return self.count_col
        return "n_markers_" + re.sub(r"\W+", "_", self.label.strip().lower())


def tumor_override(
    gene_set=None,
    *,
    label: str = "Tumor",
    min_distinct: int = 2,
    min_cluster_frac: float = 0.4,
    lineage_tfs=None,
    tf_min: float = 1.0,
    rescue_frac: float = 0.1,
    rescue_min_distinct: int = 1,
    min_expr: float = 0.0,
    count_col: str | None = None,
) -> MarkerCountOverride:
    """Build a :class:`MarkerCountOverride` encoding the h37-style tumor call, so a
    new cancer type gets the rule in one line (#340 follow-up).

    The rule (h37's osteosarcoma relabel, generalized): a cluster is ``label``
    when **either**

    - **primary** — a fraction ``min_cluster_frac`` of its cells each express
      ``>= min_distinct`` distinct genes from ``gene_set``; **or**
    - **rescue** — the cluster-mean of ``lineage_tfs`` is ``>= tf_min`` *and* a
      fraction ``rescue_frac`` of its cells express ``>= rescue_min_distinct``
      (default **1**, the broad any-marker signal) — catching low-coverage tumor
      whose sparse markers dropped below the primary bar.

    ``gene_set`` defaults to oncoref's pan-cancer cancer-testis-antigen panel —
    :func:`oncoref.CTA_gene_names` (HPA-filtered + expressed) **unioned with**
    :func:`oncoref.cta_clinical_target_gene_names` (canonical clinically-expressed
    CTAs — NY-ESO-1/CTAG2, MAGEA11, … — that testis-restriction filtering excludes
    but are real tumor markers, #350). CTAs are shared across solid tumors, so the
    **panel is reused across cancers** and only ``lineage_tfs`` (the rescue) is
    tissue-specific. tcrsift owns the counting + relabel mechanism; the curated CTA
    gene sets live in oncoref (#310). Pass your own ``gene_set`` (e.g.
    ``oncoref.CTA_testis_restricted_gene_names()`` or a tuned
    ``oncoref.CTA_by_axes(...)``) to override. ``lineage_tfs`` is a gene list,
    single gene, or precomputed ``adata.obs`` score column:

    >>> # osteosarcoma: osteoblastic TFs rescue low-CTA-coverage tumor
    >>> tumor_override(lineage_tfs=["RUNX2", "SATB2"])
    >>> # melanoma
    >>> tumor_override(lineage_tfs=["MITF", "SOX10"])
    >>> # a carcinoma with no rescue signal — primary bar only, default panel
    >>> tumor_override()

    Pass the result to ``annotate_cells(overrides=[...])``. Omitting
    ``lineage_tfs`` disables the rescue (primary bar only).
    """
    if gene_set is None:
        import oncoref

        # oncoref's testis-restricted default (CTA_gene_names) drops several
        # canonical clinically-expressed CTAs (NY-ESO-1/CTAG2, MAGEA11, …) that are
        # bona-fide tumor markers; union in the clinical-target set so the default
        # doesn't undercount tumor (#350). Sorted for a reproducible override
        # (order is immaterial to the marker count).
        gene_set = sorted(
            set(oncoref.CTA_gene_names())
            | set(oncoref.cta_clinical_target_gene_names())
        )
    rescue = (
        (lineage_tfs, tf_min, rescue_frac, rescue_min_distinct)
        if lineage_tfs is not None else None
    )
    return MarkerCountOverride(
        label=label,
        gene_set=tuple(gene_set),
        min_distinct=min_distinct,
        min_cluster_frac=min_cluster_frac,
        rescue=rescue,
        min_expr=min_expr,
        count_col=count_col,
    )


def _distinct_marker_counts(adata, gene_set, min_expr: float = 0.0) -> np.ndarray:
    """Per-cell count of DISTINCT genes from ``gene_set`` expressed (> ``min_expr``).

    Genes absent from the matrix are skipped; detection (``> min_expr`` with the
    default ``0.0``) is orientation-invariant between raw counts and log1p.
    """
    present = [g for g in gene_set if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=int)
    X = adata[:, present].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return np.asarray((X > min_expr).sum(axis=1)).ravel().astype(int)


def _rescue_cluster_mean(adata, score, leiden_col: str) -> dict:
    """Per-cluster mean of a rescue ``score``: an obs column name, a single gene,
    or a gene list (block-mean). Empty dict if unresolvable."""
    if isinstance(score, str):
        if score in adata.obs.columns:
            v = pd.to_numeric(adata.obs[score], errors="coerce").to_numpy()
            return (
                pd.Series(v, index=adata.obs.index)
                .groupby(adata.obs[leiden_col], observed=True).mean().to_dict()
            )
        return _cluster_gene_mean(adata, score, leiden_col)
    v = _block_mean(adata, list(score))
    return (
        pd.Series(v, index=adata.obs.index)
        .groupby(adata.obs[leiden_col], observed=True).mean().to_dict()
    )


def _apply_marker_overrides(adata, labels: dict, leiden_col: str, overrides) -> dict:
    """Relabel clusters per :class:`MarkerCountOverride`, and write each override's
    per-cell distinct-marker count to ``adata.obs[count_col]``. First matching
    override wins; a cluster already carrying an override label is left alone."""
    leiden = adata.obs[leiden_col]
    overridden: set = set()
    for ov in overrides:
        counts = _distinct_marker_counts(adata, ov.gene_set, ov.min_expr)
        adata.obs[ov.resolved_count_col()] = counts
        hit = pd.Series(counts >= ov.min_distinct, index=adata.obs.index)
        frac = hit.groupby(leiden, observed=True).mean()
        rescue_mean: dict = {}
        rescue_frac = None
        if ov.rescue is not None:
            rescue_mean = _rescue_cluster_mean(adata, ov.rescue[0], leiden_col)
            # The rescue fraction may gate on a BROADER distinct-count floor than
            # the primary bar: rescue[3] (rescue_min_distinct) defaults to
            # ov.min_distinct — a 3-tuple keeps the old behavior — while e.g. 1
            # lets a low-coverage cluster fire on "score high AND >=1 marker in
            # >=min_frac of cells" (#339).
            rescue_min_distinct = ov.rescue[3] if len(ov.rescue) > 3 else ov.min_distinct
            rhit = pd.Series(counts >= rescue_min_distinct, index=adata.obs.index)
            rescue_frac = rhit.groupby(leiden, observed=True).mean()
        for cl in frac.index:
            if cl in overridden:
                continue
            fire = frac.loc[cl] >= ov.min_cluster_frac
            if not fire and ov.rescue is not None:
                rmin, rfrac = ov.rescue[1], ov.rescue[2]
                fire = (rescue_mean.get(cl, 0.0) >= rmin
                        and rescue_frac.loc[cl] >= rfrac)
            if fire:
                labels[cl] = ov.label
                overridden.add(cl)
    return labels


def _warn_dropped_states(adata, prefix: str, registry) -> None:
    """Warn when :func:`score_reference` wrote a ``<prefix>::`` column whose state
    the active registry won't read — the silent-drop case #337 is about. Fires
    only on genuine drift (a scored state absent from the registry), so the
    all-defaults path is quiet."""
    dropped = sorted(
        c for c in adata.obs.columns
        if c.startswith(prefix + "::") and c.split("::", 1)[1] not in registry
    )
    if dropped:
        lineage = "T" if prefix == "st" else "B"
        param = "t" if prefix == "st" else "b"
        logger.warning(
            "annotate_clusters: %d scored %s-state column(s) not in the active "
            "%s_state_reference and will be ignored: %s — pass the matching "
            "registry as %s_state_reference= to include them.",
            len(dropped), lineage, param, ", ".join(dropped), param,
        )


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
    t_state_reference: dict | None = None,
    b_state_reference: dict | None = None,
    allowed_types=None,
    gates: AnnotationGates = DEFAULT_GATES,
    overrides=None,
) -> dict:
    """Label each Leiden cluster: cell-type signature argmax → label; T clusters
    get CD8/CD4 + a gated dominant state; weakly-typed clusters → ``"Other"``.

    ``reference`` is the cell-type registry (default
    :data:`tcrsift.signatures.CELL_TYPE_SIGNATURES`); ``t_state_reference`` /
    ``b_state_reference`` are the T-state / B-state registries read for the
    per-cluster state call (defaults
    :data:`tcrsift.signatures.T_STATE_SIGNATURES` /
    :data:`~tcrsift.signatures.B_STATE_SIGNATURES`). Pass your own to have
    caller-defined states considered — the registries must match those given to
    :func:`score_reference`, or a scored ``st::``/``bt::`` state absent from the
    active registry is ignored (a WARNING flags this drift, #337).
    ``allowed_types`` restricts the argmax (e.g.
    :data:`tcrsift.signatures.PBMC_CULTURE_TYPES` for a culture, so a
    moDC/macrophage cluster can't win a stromal signature it shares an activation
    program with). ``overrides`` is a list of :class:`MarkerCountOverride` applied
    per-cluster **after** the argmax (before label composition), relabeling a
    cluster on positive marker-count evidence (e.g. CTA→Tumor) and writing each
    override's per-cell distinct-marker count to ``adata.obs``. Requires
    :func:`score_reference` to have run. Returns ``{cluster -> label}``.
    """
    reference = reference if reference is not None else CELL_TYPE_SIGNATURES
    t_state_reference = (t_state_reference if t_state_reference is not None
                         else T_STATE_SIGNATURES)
    b_state_reference = (b_state_reference if b_state_reference is not None
                         else B_STATE_SIGNATURES)
    types = list(reference) if allowed_types is None else [
        t for t in reference if t in allowed_types
    ]
    ct_cols = [f"ct::{t}" for t in types if f"ct::{t}" in adata.obs.columns]
    st_cols = [f"st::{s}" for s in t_state_reference if f"st::{s}" in adata.obs.columns]
    bt_cols = [f"bt::{s}" for s in b_state_reference if f"bt::{s}" in adata.obs.columns]
    _warn_dropped_states(adata, "st", t_state_reference)
    _warn_dropped_states(adata, "bt", b_state_reference)
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

    if overrides:
        labels = _apply_marker_overrides(adata, labels, leiden_col, overrides)
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
    t_state_reference: dict | None = None,
    b_state_reference: dict | None = None,
    allowed_types=None,
    gates: AnnotationGates = DEFAULT_GATES,
    overrides=None,
    normalize: bool = True,
    min_subtype_frac: float = 0.01,
    add_suffix: bool = True,
):
    """One-call annotator: score → per-cluster type/state → composed labels (#312).

    Writes ``adata.obs["phenotype_base"]`` (the per-cluster type/state, e.g.
    ``"CD8 effector/cytotoxic"``) and, when ``add_suffix``,
    ``adata.obs["phenotype"]`` (with the distinguishing-gene suffix for
    non-negligible sub-clusters). ``reference`` / ``t_state_reference`` /
    ``b_state_reference`` override the cell-type / T-state / B-state registries;
    they are threaded to **both** :func:`score_reference` and
    :func:`annotate_clusters`, so a caller-defined type or state is consistently
    scored *and* read (default registries are PBMC/blood-oriented — see #340).
    ``overrides`` (a list of :class:`MarkerCountOverride`) relabels clusters on
    positive marker-count evidence (e.g. CTA→Tumor, #325) and adds a per-cell
    distinct-marker-count column to ``adata.obs``. ``normalize=True`` scores on a
    log1p-CP10K copy (from the ``counts`` layer if present) so the embedded
    raw-count AnnData can be passed straight through. Returns ``adata``.
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

    score_reference(
        scored, cell_types=reference,
        t_states=t_state_reference, b_states=b_state_reference,
    )
    labels = annotate_clusters(
        scored, leiden_col=leiden_col, reference=reference,
        t_state_reference=t_state_reference, b_state_reference=b_state_reference,
        allowed_types=allowed_types, gates=gates, overrides=overrides,
    )
    # Propagate each override's per-cell distinct-marker count from the scored
    # copy back onto the returned adata (score_reference writes to `scored`).
    for ov in overrides or []:
        col = ov.resolved_count_col()
        if col in scored.obs.columns:
            adata.obs[col] = scored.obs[col].to_numpy()
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
