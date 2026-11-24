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

"""Cohort-level constant-allele audit and novel-allele detection."""

from __future__ import annotations

import logging

import pandas as pd

from .constants import (
    _MAX_DIVERGENCE_SCAN_POSITIONS,
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_DIVERGENT_CONTIG,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_SPARSE_CONTIG,
    CONSTANT_ALLELES_BY_SPECIES,
    CONSTANT_REGIONS_AA_BY_SPECIES,
    HUMAN_CONSTANT_ALLELES,
    HUMAN_CONSTANT_REGIONS_AA,
    _score_allele_against_contig,
)
from .qc import _infer_constant_species

logger = logging.getLogger(__name__)


def _parse_divergence_positions_string(s: str) -> list[tuple[int, str, str]]:
    """Inverse of :func:`_format_divergence_positions`."""
    if not isinstance(s, str) or not s:
        return []
    out: list[tuple[int, str, str]] = []
    for chunk in s.split(";"):
        if ":" not in chunk or "->" not in chunk:
            continue
        pos_str, change = chunk.split(":", 1)
        exp, obs = change.split("->", 1)
        try:
            out.append((int(pos_str), exp, obs))
        except ValueError:
            continue
    return out


_NOVEL_ALLELE_EMPTY_COLUMNS: tuple[str, ...] = (
    "chain", "gene", "reference_allele", "candidate_allele_label",
    "variant_description", "position", "expected_aa", "observed_aa",
    "n_clones", "n_observed_at_position", "pct_observed_at_position",
    "pct_chain_clones", "n_v_genes", "n_samples", "verdict",
)


def _expand_samples_unique(samples_series: pd.Series) -> int:
    """Count distinct semicolon-joined sample names in a Series."""
    sample_set: set[str] = set()
    for s in samples_series.dropna().astype(str):
        for piece in s.split(";"):
            piece = piece.strip()
            if piece:
                sample_set.add(piece)
    return len(sample_set)


def _candidate_constant_alleles_for_chain(
    chain: str, *, species: str = "human"
) -> list[tuple[str, str, str]]:
    """Return ``(gene, allele_label, aa)`` references compatible with chain."""
    pools = CONSTANT_ALLELES_BY_SPECIES.get(species, HUMAN_CONSTANT_ALLELES)
    if chain == "alpha":
        genes = ("TRAC",)
    elif chain == "beta":
        genes = ("TRBC1", "TRBC2")
    else:
        genes = tuple(pools)

    refs: list[tuple[str, str, str]] = []
    for gene in genes:
        for allele, aa in pools.get(gene, {}).items():
            refs.append((gene, f"{gene}*{allele}", aa))
    return refs


def _default_reference_for_gene(
    canonical_name: object, *, species: str = "human"
) -> tuple[str | None, str | None]:
    """Return the default packaged allele label + AA for a known gene."""
    gene = str(canonical_name)
    pools = CONSTANT_ALLELES_BY_SPECIES.get(species, HUMAN_CONSTANT_ALLELES)
    regions = CONSTANT_REGIONS_AA_BY_SPECIES.get(species, HUMAN_CONSTANT_REGIONS_AA)
    pool = pools.get(gene, {})
    if pool:
        allele, aa = next(iter(pool.items()))
        return f"{gene}*{allele}", aa
    aa = regions.get(gene)
    if aa:
        return gene, aa
    return None, None


def _closest_reference_for_observations(
    chain: str,
    observations: pd.Series,
    *,
    species: str = "human",
) -> tuple[str | None, str | None]:
    """Pick the nearest packaged constant allele for a novel gene label."""
    obs_values = [
        str(obs)
        for obs in observations.dropna()
        if isinstance(obs, str) and obs
    ]
    if not obs_values:
        return None, None

    best_label: str | None = None
    best_aa: str | None = None
    best_key: tuple[float, int, int] = (-1.0, -1, -1)
    for _, label, allele_aa in _candidate_constant_alleles_for_chain(
        chain, species=species
    ):
        n_agree_total = 0
        n_compared_total = 0
        for obs in obs_values:
            n_agree, n_compared = _score_allele_against_contig(obs, allele_aa)
            n_agree_total += n_agree
            n_compared_total += n_compared
        if n_compared_total == 0:
            continue
        score = n_agree_total / n_compared_total
        key = (score, n_agree_total, n_compared_total)
        if key > best_key:
            best_key = key
            best_label = label
            best_aa = allele_aa
    return best_label, best_aa


def _reference_for_gene_or_observations(
    chain: str,
    canonical_name: object,
    observations: pd.Series,
    *,
    species: str = "human",
) -> tuple[str | None, str | None]:
    """Resolve the reference allele used to describe divergences."""
    reference_label, reference_aa = _default_reference_for_gene(
        canonical_name, species=species
    )
    if reference_aa:
        return reference_label, reference_aa
    return _closest_reference_for_observations(
        chain, observations, species=species
    )


def _variant_description(position: int, expected_aa: str, observed_aa: str) -> str:
    """Compact protein-diff descriptor, e.g. ``p3N>K``."""
    return f"p{int(position)}{expected_aa}>{observed_aa}"


def _candidate_allele_label(
    gene: object,
    reference_allele: object,
    variant_description: object,
) -> str:
    """Human-readable candidate label anchored to a packaged reference."""
    gene_s = str(gene)
    ref_s = str(reference_allele)
    desc_s = str(variant_description)
    ref_gene = ref_s.split("*", 1)[0]
    if gene_s == ref_gene:
        return f"{ref_s}:{desc_s}"
    return f"{gene_s}~{ref_s}:{desc_s}"


def _count_observed_positions_for_chain(
    df: pd.DataFrame,
    chain: str,
) -> pd.DataFrame:
    """Count clones with contig coverage at each observed C-region position.

    #119's denominator is "all contigs observed at this position", not
    only no-call clones and not every clone in the cohort. Coverage can
    drop with depth into the constant region, so each position gets its
    own denominator.
    """
    obs_col = f"{chain}_observed_constant_aa_start"
    gene_col = f"{chain}_c_gene_canonical"
    columns = ["chain", "gene", "position", "n_observed_at_position"]
    if obs_col not in df.columns or gene_col not in df.columns:
        return pd.DataFrame(columns=columns)

    observed = df.loc[df[obs_col].notna(), [gene_col, obs_col]].copy()
    if observed.empty:
        return pd.DataFrame(columns=columns)

    observed[obs_col] = observed[obs_col].astype(str)
    observed = observed[observed[obs_col].str.len() > 0]
    if observed.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, object]] = []
    # ``observed=True`` skips unused categorical levels (AnnData/H5AD
    # round-trips leave them behind), so no empty groups reach the body.
    for canonical_name, gene_slice in observed.groupby(
        gene_col, dropna=True, observed=True,
    ):
        obs_series = gene_slice[obs_col]
        max_pos = min(int(obs_series.str.len().max()), _MAX_DIVERGENCE_SCAN_POSITIONS)
        for i in range(max_pos):
            n_observed = int(obs_series.str.get(i).notna().sum())
            if n_observed:
                rows.append({
                    "chain": chain,
                    "gene": str(canonical_name),
                    "position": i + 1,
                    "n_observed_at_position": n_observed,
                })
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns)


def _collect_divergences_for_chain(
    df: pd.DataFrame,
    chain: str,
) -> pd.DataFrame:
    """Return a long-form DataFrame of divergence rows for ``chain``.

    Two sources are unioned:

    * **Stored** — clones with ``_allele_called`` set AND non-empty
      ``_allele_divergence_positions``. Parsed via
      :func:`_parse_divergence_positions_string`. Captures
      ``auto_detected`` clones whose contig diverges from the called
      allele at tolerated positions — the heterozygous-donor case
      from #120.
    * **Recomputed** — clones with no-call reasons (no allele was
      committed, so the stored divergence column is ``None``).
      Compare ``_observed_constant_aa_start`` vs the gene's default
      canonical AA position-by-position. If the gene label is novel,
      compare against the nearest compatible packaged allele.
    """
    reason_col = f"{chain}_allele_called_reason"
    if reason_col not in df.columns:
        return pd.DataFrame()

    called_col = f"{chain}_allele_called"
    divs_col = f"{chain}_allele_divergence_positions"
    obs_col = f"{chain}_observed_constant_aa_start"
    gene_col = f"{chain}_c_gene_canonical"
    v_col = f"{chain}_v_gene"

    if gene_col not in df.columns:
        return pd.DataFrame()

    collected_frames: list[pd.DataFrame] = []

    # ---- Source 1: stored divergences from auto_detected clones ----
    if divs_col in df.columns and called_col in df.columns:
        keep_cols = [called_col, divs_col, gene_col]
        if v_col in df.columns:
            keep_cols.append(v_col)
        if "samples" in df.columns:
            keep_cols.append("samples")
        stored = df.loc[
            df[divs_col].notna() & df[called_col].notna(),
            keep_cols,
        ].copy()
        if "samples" not in stored.columns:
            stored["samples"] = ""
        if not stored.empty:
            # Parse each clone's divergence string into a list, then
            # explode so each (clone, position) becomes its own row.
            stored["__parsed"] = stored[divs_col].apply(
                _parse_divergence_positions_string
            )
            stored = stored[stored["__parsed"].str.len() > 0]
            if not stored.empty:
                stored = stored.explode("__parsed", ignore_index=False)
                stored = stored[stored["__parsed"].notna()]
                stored["position"] = stored["__parsed"].str[0].astype(int)
                stored["expected_aa"] = stored["__parsed"].str[1]
                stored["observed_aa"] = stored["__parsed"].str[2]
                stored = stored.drop(columns=["__parsed", divs_col])
                stored = stored.rename(columns={gene_col: "gene"})
                stored["gene"] = stored["gene"].astype(str)
                stored["reference_allele"] = stored[called_col].astype(str)
                if v_col in stored.columns:
                    stored = stored.rename(columns={v_col: "v_gene"})
                else:
                    stored["v_gene"] = None
                stored = stored.rename(columns={"samples": "sample"})
                stored["clone_idx"] = stored.index
                collected_frames.append(stored[[
                    "clone_idx", "gene", "reference_allele", "position",
                    "expected_aa", "observed_aa", "v_gene", "sample",
                ]])

    # ---- Source 2: recomputed for no-call clones ----
    # Compares against ``HUMAN_CONSTANT_REGIONS_AA`` (the *default*
    # canonical), not the user-overridden allele. That's safe here
    # because user-override clones get ``reason=overridden``, which
    # is not in ``no_call_reasons`` — so we never recompute for
    # them. The stored-source path above handles overridden clones
    # via their own ``_allele_divergence_positions`` column.
    no_call_reasons = {
        ALLELE_REASON_DIVERGENT_CONTIG,
        ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
        ALLELE_REASON_SPARSE_CONTIG,
    }
    if obs_col in df.columns:
        recompute_mask = (
            df[reason_col].isin(no_call_reasons)
            & df[obs_col].notna()
        )
        # Process per-gene so we have one canonical AA string to
        # compare against; vectorized within each gene group.
        # ``observed=True`` skips unused categorical levels, so no empty
        # groups reach the body (AnnData/H5AD round-trip robustness).
        for canonical_name, gene_slice in df.loc[
            recompute_mask
        ].groupby(gene_col, dropna=True, observed=True):
            const_col = f"{chain}_constant_aa"
            full_col = f"full_{chain}_aa"
            first_const = gene_slice[const_col].iloc[0] if const_col in gene_slice.columns else None
            first_full = gene_slice[full_col].iloc[0] if full_col in gene_slice.columns else None
            slice_sp = _infer_constant_species(
                first_const, first_full, str(canonical_name).upper(), default="human"
            )
            reference_label, reference_aa = _reference_for_gene_or_observations(
                chain, canonical_name, gene_slice[obs_col], species=slice_sp
            )
            if not reference_aa or not reference_label:
                continue
            obs_series = gene_slice[obs_col]
            v_series = gene_slice[v_col] if v_col in gene_slice.columns else pd.Series(
                [None] * len(gene_slice), index=gene_slice.index,
            )
            sample_series = gene_slice["samples"] if "samples" in gene_slice.columns else pd.Series(
                [None] * len(gene_slice), index=gene_slice.index,
            )
            max_pos = min(len(reference_aa), _MAX_DIVERGENCE_SCAN_POSITIONS)
            per_position: list[pd.DataFrame] = []
            for i in range(max_pos):
                canonical_residue = reference_aa[i]
                # str.get returns None when index is past the string length.
                obs_residue = obs_series.str.get(i)
                diverged_mask = obs_residue.notna() & (obs_residue != canonical_residue)
                if not diverged_mask.any():
                    continue
                per_position.append(pd.DataFrame({
                    "clone_idx": gene_slice.index[diverged_mask],
                    "gene": str(canonical_name),
                    "reference_allele": reference_label,
                    "position": i + 1,
                    "expected_aa": canonical_residue,
                    "observed_aa": obs_residue[diverged_mask].values,
                    "v_gene": v_series[diverged_mask].values,
                    "sample": sample_series[diverged_mask].values,
                }))
            if per_position:
                collected_frames.append(pd.concat(per_position, ignore_index=True))

    if not collected_frames:
        return pd.DataFrame()
    out = pd.concat(collected_frames, ignore_index=True)
    out.insert(0, "chain", chain)
    return out


def detect_novel_alleles(
    df: pd.DataFrame,
    *,
    min_pct: float = 0.05,
    min_v_spread: int = 3,
    min_samples: int = 2,
    min_cohort_size: int = 20,
) -> pd.DataFrame:
    """Aggregate per-clone divergence positions into cohort-level
    novel-allele candidates (#119).

    For each ``(chain, position, expected_aa, observed_aa)`` tuple
    across all clones whose contig diverges from canonical (either
    stored ``_allele_divergence_positions`` on auto-called clones or
    recomputed against canonical for no-call clones), count:

    * ``n_clones`` — number of clones carrying the variant.
    * ``pct_chain_clones`` — fraction of in-cohort clones for that
      chain showing the variant.
    * ``n_observed_at_position`` — number of clones with contig AA
      coverage at that chain/gene/position.
    * ``pct_observed_at_position`` — fraction of observed contigs at
      that position showing the variant. This is the candidate-calling
      denominator.
    * ``n_v_genes`` — number of distinct V genes the variant
      appears in.
    * ``n_samples`` — number of distinct CellRanger samples the
      variant appears in.
    * ``verdict`` — ``"novel_allele_candidate"`` when
      ``pct_observed_at_position >= min_pct`` AND
      ``n_v_genes >= min_v_spread`` AND ``n_samples >= min_samples``;
      ``"likely_artifact"`` otherwise.

    Returns a long-form DataFrame sorted by ``n_clones`` descending.
    Empty DataFrame when ``df`` has fewer than ``min_cohort_size``
    clones with the chain populated (the heuristic is unreliable
    below that floor) or when no divergences are present.

    Vectorized: ~100x faster than the original ``iterrows`` loop on
    cohorts of 10k+ clones (no per-row Python). For 100k clones the
    runtime is dominated by pandas groupby — still single-threaded
    but fully numpy under the hood.

    The heuristic encodes simple population genetics: a real
    polymorphism distributes across V-genes and samples (since V-gene
    usage is independent of constant-region allele); an artifact
    concentrates in one V or one sample (PCR / assembly noise tied to
    a specific V sequence). See #119 for the full rationale and pilot
    validation.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`assemble_full_sequences`.
    min_pct : float, default 0.05
        Minimum fraction of contigs observed at that position for the
        variant to be called a candidate. 5% is conservative for
        heterozygous donors (~50% expected) and dominates over
        single-V artifacts (~1%) but still surfaces minor population
        alleles.
    min_v_spread : int, default 3
        Minimum distinct V genes the variant must appear in. V-gene
        usage is independent of constant allele, so true polymorphisms
        spread across V's; V-concentrated variants are almost always
        contig-assembly artifacts.
    min_samples : int, default 2
        Minimum distinct samples. Single-sample variants may reflect
        a batch effect rather than a real allele.
    min_cohort_size : int, default 20
        Skip the aggregator entirely when neither chain has at least
        this many clones with ``_allele_called_reason`` populated —
        ``min_v_spread`` can't exceed the cohort size, so smaller
        cohorts produce noise verdicts. Set to 0 to force evaluation
        on any cohort.
    """
    empty = pd.DataFrame(columns=list(_NOVEL_ALLELE_EMPTY_COLUMNS))

    # Compute chain totals up-front and apply the min_cohort_size gate.
    chain_totals: dict[str, int] = {}
    for chain in ("alpha", "beta"):
        reason_col = f"{chain}_allele_called_reason"
        if reason_col in df.columns:
            chain_totals[chain] = int(df[reason_col].notna().sum())
    if (
        min_cohort_size > 0
        and max(chain_totals.values(), default=0) < min_cohort_size
    ):
        return empty

    per_chain: list[pd.DataFrame] = []
    observed_position_counts: list[pd.DataFrame] = []
    for chain in ("alpha", "beta"):
        counts = _count_observed_positions_for_chain(df, chain)
        if not counts.empty:
            observed_position_counts.append(counts)
        chain_df = _collect_divergences_for_chain(df, chain)
        if not chain_df.empty:
            per_chain.append(chain_df)
    if not per_chain:
        return empty

    long_df = pd.concat(per_chain, ignore_index=True)

    # Vectorized aggregation. ``unique`` -> nunique on clone_idx and
    # v_gene; samples need to be expanded via the helper.
    # NOTE: ``reference_allele`` is deliberately NOT a groupby key. The same
    # physical variant — identical (chain, gene, position, expected_aa,
    # observed_aa) — can arrive labelled against different reference alleles
    # (Source 1 tags it with the *called* allele, e.g. TRBC2*03; Source 2
    # recomputes against the gene *default*, TRBC2*01). ``expected_aa`` already
    # pins the reference residue at that position, so equal-expected rows are
    # the same variant and must aggregate together; keying on reference_allele
    # would split them and deflate n_clones/n_v_genes/n_samples below the
    # candidate thresholds. We pick a deterministic representative label below.
    grouped = long_df.groupby(
        [
            "chain", "gene", "position",
            "expected_aa", "observed_aa",
        ],
        dropna=False, sort=False, observed=True,
    )
    n_clones = grouped["clone_idx"].nunique().rename("n_clones")
    n_v_genes = grouped["v_gene"].nunique(dropna=True).rename("n_v_genes")
    n_samples = grouped["sample"].agg(_expand_samples_unique).rename("n_samples")
    reference_allele = grouped["reference_allele"].agg(
        lambda s: sorted(set(s.dropna().astype(str)))[0] if s.notna().any() else None
    ).rename("reference_allele")
    summary = pd.concat(
        [n_clones, n_v_genes, n_samples, reference_allele], axis=1
    ).reset_index()

    if observed_position_counts:
        position_counts = pd.concat(observed_position_counts, ignore_index=True)
        summary = summary.merge(
            position_counts,
            on=["chain", "gene", "position"],
            how="left",
        )
    else:
        summary["n_observed_at_position"] = pd.NA

    # Vectorized per-row pct: map each row's chain to its total, then
    # clip the denominator to ≥1 so zero-total chains (impossible
    # given the gate above but cheap to guard) don't produce inf.
    chain_total_series = summary["chain"].map(chain_totals).fillna(1).clip(lower=1)
    summary["pct_chain_clones"] = (
        summary["n_clones"] / chain_total_series
    ).round(4)
    # Prefer the position-specific observed-contig denominator. For
    # legacy CSVs that have divergence columns but not observed AA
    # columns, fall back to the old chain-level denominator. The
    # observed-position count and the variant count come from different
    # columns (``_observed_constant_aa_start`` vs the parsed divergence
    # string), so clamp the denominator to at least ``n_clones`` — a
    # variant clone is by definition observed at its own position, and
    # this guarantees ``pct_observed_at_position`` stays in [0, 1] even
    # if the two columns disagree.
    summary["n_observed_at_position"] = (
        summary["n_observed_at_position"]
        .fillna(chain_total_series)
        .astype(int)
    )
    summary["n_observed_at_position"] = summary[
        ["n_observed_at_position", "n_clones"]
    ].max(axis=1).clip(lower=1)
    summary["pct_observed_at_position"] = (
        summary["n_clones"] / summary["n_observed_at_position"]
    ).round(4)
    summary["variant_description"] = [
        _variant_description(pos, exp, obs)
        for pos, exp, obs in zip(
            summary["position"], summary["expected_aa"], summary["observed_aa"],
        )
    ]
    summary["candidate_allele_label"] = [
        _candidate_allele_label(gene, ref, desc)
        for gene, ref, desc in zip(
            summary["gene"],
            summary["reference_allele"],
            summary["variant_description"],
        )
    ]
    summary["verdict"] = (
        (summary["pct_observed_at_position"] >= min_pct)
        & (summary["n_v_genes"] >= min_v_spread)
        & (summary["n_samples"] >= min_samples)
    ).map(
        {True: "novel_allele_candidate", False: "likely_artifact"}
    )
    summary["position"] = summary["position"].astype(int)
    summary = summary[list(_NOVEL_ALLELE_EMPTY_COLUMNS)].sort_values(
        ["n_clones", "pct_observed_at_position"], ascending=[False, False],
    ).reset_index(drop=True)
    return summary


def allele_audit_report(
    df: pd.DataFrame,
    *,
    min_pct: float = 0.05,
    min_v_spread: int = 3,
    min_samples: int = 2,
    min_cohort_size: int = 20,
) -> str:
    """Human-readable cohort-level allele audit (#120 ask 3).

    Combines:
    * per-chain confident allele tally (calls vs. no-calls broken out
      by ``allele_called_reason``)
    * cohort-level novel-allele candidates from
      :func:`detect_novel_alleles`

    ``min_cohort_size`` (default 20) is forwarded to
    :func:`detect_novel_alleles` — smaller cohorts skip the
    aggregator entirely (it adds noticeable latency on every
    verbose assembly and produces noise verdicts below the
    ``min_v_spread`` floor).

    Output is suitable for printing at the end of a ``tcrsift run``
    pipeline or as the body of ``tcrsift audit-alleles``.
    """
    lines: list[str] = ["[allele audit]"]
    for chain in ("alpha", "beta"):
        called_col = f"{chain}_allele_called"
        reason_col = f"{chain}_allele_called_reason"
        if called_col not in df.columns or reason_col not in df.columns:
            continue
        chain_mask = df[reason_col].notna()
        n_total = int(chain_mask.sum())
        if n_total == 0:
            continue
        n_called = int(df[called_col].notna().sum())
        called_counts = df.loc[
            df[called_col].notna(), called_col
        ].value_counts()
        reason_counts = df[reason_col].value_counts(dropna=False)
        lines.append(
            f"  {chain} chain: {n_called}/{n_total} called "
            f"({100 * n_called / max(n_total, 1):.1f}%)"
        )
        for allele, n in called_counts.items():
            lines.append(f"    {allele}: {n}")
        no_calls = n_total - n_called
        if no_calls > 0:
            lines.append(f"    No-call breakdown ({no_calls} clones):")
            for reason, n in reason_counts.items():
                if reason in {ALLELE_REASON_AUTO_DETECTED, ALLELE_REASON_OVERRIDDEN}:
                    continue
                if pd.isna(reason):
                    continue
                lines.append(f"      {reason}: {n}")

    novel = detect_novel_alleles(
        df,
        min_pct=min_pct,
        min_v_spread=min_v_spread,
        min_samples=min_samples,
        min_cohort_size=min_cohort_size,
    )
    if not novel.empty:
        lines.append("")
        lines.append("  Observed polymorphisms (contig vs. canonical):")
        for _, row in novel.iterrows():
            verdict_tag = (
                "← novel_allele_candidate"
                if row["verdict"] == "novel_allele_candidate"
                else "→ likely_artifact"
            )
            lines.append(
                f"    {row['chain']} {row['candidate_allele_label']}  "
                f"{row['expected_aa']}→{row['observed_aa']}  "
                f"n={row['n_clones']}/{row['n_observed_at_position']} "
                f"observed ({100 * row['pct_observed_at_position']:.1f}%)  "
                f"V-genes={row['n_v_genes']}  samples={row['n_samples']}  "
                f"{verdict_tag}"
            )
        candidates = novel[novel["verdict"] == "novel_allele_candidate"]
        if candidates.empty:
            lines.append("  No novel-allele candidates detected.")
    else:
        lines.append("")
        lines.append("  No contig-vs-canonical divergences observed.")
    return "\n".join(lines)
