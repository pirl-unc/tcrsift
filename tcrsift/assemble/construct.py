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

"""Full-length TCR construct assembly, single-chain linking, and variant summaries."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from ..leaders import normalize_vgene
from ..validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_directory_exists,
)
from .alleles import allele_audit_report
from .codons import (
    DEFAULT_FORBIDDEN_MOTIFS,
    CodonOptimization,
    _add_single_chain,
    optimize_codons,
    recode_codons,
    stop_codons_nt,
)
from .constants import (
    DEFAULT_LEADERS,
    _apply_cohort_alpha_junctions,
    _count_divergence_positions,
    _describe_leader,
    _ensure_qc_warnings_column,
    _lookup_default_leader,
    _record_leader,
    _resolve_c_gene,
    _trbc_override_mask,
    _warn_constant_allele_divergence,
    get_constant_region_sequences,
)
from .contigs import _add_constant_from_row, _add_constant_regions, load_contigs
from .germline import (
    VDJ_SOURCE_5PRIME_RESCUED,
    VDJ_SOURCE_VJ_RECONSTRUCTED,
    reconstruct_germline_vdj,
)
from .leader_policy import (
    _FROM_CONTIG_SOURCES,
    _extract_leader_from_contigs_single,
    apply_leader_policy,
)
from .qc import _unverified_chains, add_synthesis_qc

logger = logging.getLogger(__name__)


def assemble_full_sequences(
    clonotypes: pd.DataFrame,
    contigs_dir: str | Path | None = None,
    alpha_leader: str | None = "CD28",
    beta_leader: str | None = "CD8A",
    include_constant: bool = True,
    constant_source: str = "ensembl",
    linker: str = "T2A",
    trac_allele: str = "auto",
    trbc1_allele: str = "auto",
    trbc2_allele: str = "auto",
    stop_codons: tuple[str, ...] = ("TAA", "TGA"),
    verbose: bool = True,
    show_progress: bool = True,
    sample_name_from: str = "parent",
    cellranger_dir: str | Path | None = None,
    sample_sheet: pd.DataFrame | None = None,
    leader_fallback: str = "germline",
    force_alpha_leader: str | None = None,
    force_beta_leader: str | None = None,
    secondary_alpha_leader: str | None = None,
    secondary_beta_leader: str | None = None,
    codon_optimization: str = "focal",
    max_codon_repeats: int = 1,
    avoid_enzymes: tuple[str, ...] = (),
    species: str = "human",
    constant_species: str | None = None,
    codon_species: str | None = None,
    germline_vdj_fallback: bool = False,
    data_dir: str | Path | None = None,
    auto_download_germline: bool = False,
) -> pd.DataFrame:
    """
    Assemble full-length TCR sequences.

    Parameters
    ----------
    clonotypes : pd.DataFrame
        Clonotype DataFrame with VDJ sequences (from fwr1/cdr1/fwr2/cdr2/fwr3/cdr3/fwr4)
    contigs_dir : str or Path, optional
        Directory with CellRanger contig FASTA files. Required if alpha_leader or
        beta_leader is set to "from_contig".
    alpha_leader : str or None
        Leader sequence for alpha chain. Options:
        - None: No leader sequence
        - "from_contig": Extract native leader from contig FASTA (requires contigs_dir)
        - Key from DEFAULT_LEADERS: "CD8A", "CD28", "IgK", "TRAC", "TRBC"
        Default is "CD28" to provide distinct sequences from beta chain.
    beta_leader : str or None
        Leader sequence for beta chain. Same options as alpha_leader.
        Default is "CD8A" to provide distinct sequences from alpha chain.
    leader_fallback : str or None
        Signal peptide to substitute when a ``from_contig`` leader is
        implausible — weak-Kozak over-capture, out-of-window length, missing
        Met/h-region (#263). Default ``"germline"`` uses the V gene's germline
        leader; curated keys CD8A/CD28/IgK/TRAC/TRBC are also accepted. ``None``
        currently resolves to the same germline fallback. Only affects
        ``from_contig`` leaders.
    include_constant : bool
        Include constant region sequences.
    constant_source : str
        Source for constant regions:

        - ``"canonical"``: splice in the packaged canonical
          human TRAC / TRBC1 / TRBC2 (:data:`HUMAN_CONSTANT_REGIONS_AA`).
          When CellRanger contigs are available, the canonical's first
          ~15 residues are verified against the observed contig start
          and a ``qc_warnings`` entry is emitted on mismatch.
        - ``"ensembl"`` (default): back-compat alias for ``"canonical"``. The
          earlier pyensembl-backed path was removed in #66 — it read
          the full mRNA at the wrong frame offset and silently
          truncated constants to 2–11 residues.
        - ``"from-data"``: read ``{chain}_constant_aa`` /
          ``{chain}_constant_nt`` directly from the input frame.
    linker : str
        Linker sequence for single-chain constructs: "T2A", "P2A", "E2A", "F2A"
    trac_allele, trbc1_allele, trbc2_allele : str
        Per-gene allele selection (#113). Default ``"auto"`` runs the
        allele picker over packaged alleles in
        :data:`HUMAN_CONSTANT_ALLELES` scoring each against the
        donor's contig translation; pass an explicit allele label
        (e.g. ``"01"`` or ``"03"``) to force a specific canonical.
        Invalid labels emit a QC warning and fall back to the
        default. ``constant_source="from-data"`` ignores these knobs
        (the constant comes from the input frame, not the FASTA);
        a warning is logged in that combination.
    stop_codons : tuple of str, default ``("TAA", "TGA")``
        Stop codons appended to the codon-optimized constant CDS.
        Default uses two non-redundant stops (recognized by
        different release factors → reduced read-through in
        synthesized constructs). Pass ``("TAA",)`` for single-stop
        (pre-2.4 behavior) or ``()`` to omit stops entirely. Each
        entry must be one of ``"TAA"``/``"TAG"``/``"TGA"``;
        invalid entries raise ``ValueError``. Stops are only
        appended to the ``_optimized`` columns and the final
        ``_constant_nt`` (alias) — not to ``_contig`` columns,
        which preserve the donor's CellRanger contig bytes verbatim
        (and are typically truncated past the contig coverage edge).
    verbose : bool
        Print progress information
    show_progress : bool
        Show progress bar

    Returns
    -------
    pd.DataFrame
        Clonotypes with full sequences added. Allele audit columns
        emitted per clone (#113):

        - ``{chain}_allele_called`` — e.g. ``"TRBC2*01"`` or ``None``
          when no contig + non-explicit override.
        - ``{chain}_allele_score`` — float ∈ [0, 1] (1.0 = perfect AA
          agreement; 1.0 also when the user forced an allele).
        - ``{chain}_allele_alternatives`` — semicolon-joined string
          of runner-up alleles + scores, sorted by score desc with
          FASTA order as tie-break. Example: ``"TRBC2*01:1.000;TRBC2*03:0.933"``.
          ``None`` when the picker no-decided or only one allele was
          packaged.

        Constant NT triad emitted per chain (#116):

        - ``{chain}_constant_nt_contig`` — pure CellRanger contig
          NT past the J→C junction (no canonical splicing, no
          codon optimization). Truncated where contig coverage
          ends. ``None`` when no contig is available for the clone.
        - ``{chain}_constant_nt_optimized`` — codon-optimized
          canonical CDS for the picked allele (motif-aware via
          :func:`optimize_codons`; avoids common restriction sites
          and homopolymer runs) plus the requested stop codons.
        - ``{chain}_constant_nt`` — the assembly-aware blend: uses
          donor-real bytes from the contig where they agree with
          the canonical AA, falls back to ``_optimized`` for the
          rest. This is the column most callers want; kept as the
          default for back-compat.

        The same triad is also exposed for ``full_{chain}_nt``
        (leader + VDJ + constant + stop).

    Examples
    --------
    >>> # Default: CD28 on alpha, CD8A on beta (distinct leaders)
    >>> assembled = assemble_full_sequences(clonotypes)

    >>> # No leader sequences
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader=None)

    >>> # Leader only on beta chain (first in 2A construct)
    >>> assembled = assemble_full_sequences(clonotypes, alpha_leader=None, beta_leader="CD8A")

    >>> # Extract native leaders from contig FASTAs
    >>> assembled = assemble_full_sequences(
    ...     clonotypes,
    ...     contigs_dir="/path/to/contigs",
    ...     alpha_leader="from_contig",
    ...     beta_leader="from_contig",
    ... )
    """
    from ..species import normalize_species, resolve_species

    # Validate inputs
    clonotypes = validate_clonotype_df(clonotypes, for_assembly=True)
    vdj_sp = resolve_species(species, df=clonotypes, warn_mismatch=False)
    const_sp = normalize_species(constant_species) if constant_species else vdj_sp
    cod_sp = normalize_species(codon_species) if codon_species else const_sp

    # Validate stop_codons early so users get the error before any
    # heavy lifting. Empty tuple is allowed (caller wants no stops).
    stops_nt = stop_codons_nt(stop_codons)
    recode = CodonOptimization(
        mode=codon_optimization,
        max_codon_repeats=max_codon_repeats,
        avoid_enzymes=tuple(avoid_enzymes or ()),
        codon_species=cod_sp,
    )

    valid_constant_sources = ["canonical", "ensembl", "from-data"]
    if constant_source not in valid_constant_sources:
        raise TCRsiftValidationError(
            f"Invalid constant_source: '{constant_source}'",
            hint=f"Valid options are: {valid_constant_sources}",
        )

    if constant_source == "from-data" and any(
        a != "auto" for a in (trac_allele, trbc1_allele, trbc2_allele)
    ):
        logger.warning(
            "constant_source='from-data' ignores trac_allele / "
            "trbc1_allele / trbc2_allele — the constant comes from "
            "the input frame, not the packaged canonical FASTA. "
            "Reset the allele kwargs to 'auto' or switch to "
            "constant_source='canonical' to use them."
        )

    # Validate and resolve leader options for each chain
    leader_config = {}
    for chain, leader_param in [("alpha", alpha_leader), ("beta", beta_leader)]:
        if leader_param is None:
            leader_config[chain] = None
        elif leader_param.lower() == "from_contig":
            if not contigs_dir and not cellranger_dir:
                raise TCRsiftValidationError(
                    f"{chain}_leader='from_contig' requires contigs_dir or "
                    "cellranger_dir to be specified",
                    hint="Provide contigs_dir/cellranger_dir with CellRanger FASTA files, or use a default leader like 'CD8A'",
                )
            leader_config[chain] = "from_contig"
        else:
            hit = _lookup_default_leader(leader_param)
            if hit is not None:
                leader_config[chain] = hit[1]
            else:
                raise TCRsiftValidationError(
                    f"Unknown {chain}_leader: '{leader_param}'",
                    hint=f"Valid options are: None, 'from_contig', or one of {list(DEFAULT_LEADERS.keys())}",
                )

    if verbose:
        alpha_desc = _describe_leader(alpha_leader, leader_config["alpha"])
        beta_desc = _describe_leader(beta_leader, leader_config["beta"])
        logger.info(
            f"Assembling full sequences for {len(clonotypes):,} clonotypes "
            f"(species={vdj_sp}, constant_species={const_sp}, codon_species={cod_sp})"
        )
        logger.info(f"  Alpha leader: {alpha_desc}")
        logger.info(f"  Beta leader: {beta_desc}")
        logger.info(f"  Constant regions: {include_constant} (source: {constant_source})")
        logger.info(f"  Linker: {linker}")

    df = clonotypes.copy()
    if germline_vdj_fallback:
        df = reconstruct_germline_vdj(
            df,
            species=vdj_sp,
            codon_species=cod_sp,
            data_dir=data_dir,
            auto_download=auto_download_germline,
        )
        if verbose:
            for ch in ("alpha", "beta"):
                sc = df.get(f"{ch}_vdj_source")
                if sc is not None:
                    n_5p = int((sc == VDJ_SOURCE_5PRIME_RESCUED).sum())
                    n_vj = int((sc == VDJ_SOURCE_VJ_RECONSTRUCTED).sum())
                    if n_5p or n_vj:
                        logger.info(
                            f"  Germline VDJ fallback ({ch}): "
                            f"{n_5p:,} 5'-partial rescued, {n_vj:,} V+J reconstructed"
                        )

    constant_seqs = {}
    if include_constant and constant_source in ("canonical", "ensembl"):
        if constant_source == "ensembl" and verbose:
            logger.info(
                "  constant_source='ensembl' is now an alias for 'canonical' "
                f"(#66) — splicing canonical {const_sp} TRAC / TRBC1 / TRBC2."
            )
        constant_seqs = get_constant_region_sequences(
            stop_codons=stop_codons, species=const_sp, codon_species=cod_sp
        )
        if verbose:
            logger.info(
                f"    Loaded {len(constant_seqs)} canonical constant region sequences"
            )
            if stop_codons:
                logger.info(
                    f"    Appending stop codons {list(stop_codons)} to "
                    "codon-optimized constants."
                )
            else:
                logger.info(
                    "    No stop codons configured (stop_codons=()); "
                    "constants will not be terminated."
                )

    # Warn if from-data constants requested but not present
    if include_constant and constant_source == "from-data":
        constant_cols = [
            "alpha_constant_aa",
            "alpha_constant_nt",
            "beta_constant_aa",
            "beta_constant_nt",
        ]
        if not any(col in df.columns for col in constant_cols):
            logger.warning(
                "  constant_source='from-data' but no constant region columns found in input. "
                "Constants will be omitted."
            )

    # Load contigs whenever a contigs_dir is given. Pre-1.3 we only
    # loaded when a leader was "from_contig"; that gate is too narrow
    # now that the C-region NT blend (#103) also consumes contigs to
    # splice donor-real bytes into the J→C boundary. Loading is cheap
    # relative to assembly and avoids a silent "contigs given but
    # blend silently noop'd" failure mode.
    sample_contigs: dict = {}
    # cellranger_dir is shorthand for contigs_dir + sample_name_from='grandparent'
    # (#124). Resolved here so the rest of the function only sees contigs_dir.
    if cellranger_dir is not None:
        if contigs_dir is not None:
            raise ValueError(
                "assemble_full_sequences: pass exactly one of contigs_dir or "
                "cellranger_dir, not both."
            )
        contigs_dir = cellranger_dir
        sample_name_from = "grandparent"
    if contigs_dir:
        contigs_dir = validate_directory_exists(Path(contigs_dir), "contigs directory")
        if verbose:
            logger.info(
                f"  Loading contigs from {contigs_dir} (sample_name_from={sample_name_from!r})..."
            )
        sample_contigs = load_contigs(
            contigs_dir,
            sample_name_from=sample_name_from,
            sample_sheet=sample_sheet,
        )
        if verbose:
            total_contigs = sum(len(c) for c in sample_contigs.values())
            logger.info(f"    Loaded {total_contigs:,} contigs from {len(sample_contigs)} samples")

    # Process each clonotype
    if verbose:
        logger.info("  Assembling sequences...")

    assembly_results = []

    # Build per-gene allele override dict for the picker (#113).
    allele_overrides = {
        "TRAC": trac_allele,
        "TRBC1": trbc1_allele,
        "TRBC2": trbc2_allele,
    }

    # Create iterator with optional progress bar
    row_iter = df.iterrows()
    if show_progress:
        row_iter = tqdm(
            list(df.iterrows()),
            desc="Assembling sequences",
            unit="clone",
        )

    for idx, row in row_iter:
        result = _assemble_clone(
            row,
            sample_contigs,
            constant_seqs,
            leader_config,
            include_constant,
            constant_source,
            allele_overrides=allele_overrides,
            stops_nt=stops_nt,
            recode=recode,
            vdj_species=vdj_sp,
            constant_species=const_sp,
            codon_species=cod_sp,
        )
        assembly_results.append(result)

    # Add assembly columns to dataframe
    result_df = pd.DataFrame(assembly_results)
    for col in result_df.columns:
        df[col] = result_df[col].values

    # Correct blanket-N α junctions using a per-J consensus from this cohort's
    # contig-verified clones BEFORE building single-chain constructs and leader
    # policy views so full_alpha_*, full_alpha_nt_optimized, and single_chain_*
    # remain byte-for-byte synchronized (#242).
    df = _apply_cohort_alpha_junctions(df, codon_species=cod_sp)

    # Add single-chain construct if requested
    if linker and "full_beta_aa" in df.columns and "full_alpha_aa" in df.columns:
        if verbose:
            logger.info(f"  Creating single-chain constructs with {linker} linker...")
        df = _add_single_chain(df, linker)

    # Leader policy (#270)
    df = apply_leader_policy(
        df, linker,
        leader_fallback=leader_fallback or "germline",
        force_alpha=force_alpha_leader, force_beta=force_beta_leader,
        secondary_alpha=secondary_alpha_leader, secondary_beta=secondary_beta_leader,
        verbose=verbose,
        species=vdj_sp,
        recode=recode,
    )

    # Summary
    if verbose:
        n_with_alpha = df["full_alpha_aa"].notna().sum() if "full_alpha_aa" in df.columns else 0
        n_with_beta = df["full_beta_aa"].notna().sum() if "full_beta_aa" in df.columns else 0
        n_single_chain = (
            df["single_chain_aa"].notna().sum() if "single_chain_aa" in df.columns else 0
        )
        logger.info("  Assembly complete:")
        logger.info(f"    With full alpha: {n_with_alpha:,}")
        logger.info(f"    With full beta: {n_with_beta:,}")
        logger.info(f"    Single-chain constructs: {n_single_chain:,}")

    # Per-row override audit column (#90). True where the canonical
    # picked by the J-family rule disagreed with CellRanger's raw
    # `beta_c_gene` call. Lets users filter / inspect post-hoc
    # without going through the log. Written even when no overrides
    # fired so the column is always present (defaults to False).
    if "beta_c_gene" in df.columns:
        df["beta_c_gene_overridden"] = _trbc_override_mask(df)

    # One aggregate line for J-family overrides of CellRanger TRBC
    # (#90). pick_canonical_constant is silent per-call to avoid a
    # log flood on cohorts where every other clone hits the case.
    # Gated by `verbose` so silent runs stay silent — the override
    # itself still happens, this is the audit trail. The per-row
    # `beta_c_gene_overridden` column above is the data-side audit
    # that survives the log.
    if verbose:
        n_overrides = (
            int(df["beta_c_gene_overridden"].sum())
            if "beta_c_gene_overridden" in df.columns
            else 0
        )
        if n_overrides:
            n_beta = df["beta_c_gene"].notna().sum() if "beta_c_gene" in df.columns else 0
            logger.warning(
                "  Overrode CellRanger TRBC call on %d / %d β clones to match "
                "J-family parity (TRBJ1→TRBC1, TRBJ2→TRBC2). CellRanger's "
                "TRBC1/2 discrimination is unreliable; the locus rule is "
                "authoritative.",
                n_overrides, n_beta,
            )

        # Cohort allele audit (#119, #120). Emits the per-chain
        # called/no-called breakdown + any novel-allele candidates
        # detected across the cohort. Silent on small cohorts where
        # the heuristic thresholds don't fire.
        if any(
            f"{c}_allele_called_reason" in df.columns for c in ("alpha", "beta")
        ):
            for line in allele_audit_report(df).splitlines():
                logger.info(line)

    # Constant-allele fidelity (#187). Derive a scalar divergence count per
    # chain from the per-position column the assembler recorded against the
    # donor's sequenced bases, and warn loudly (regardless of verbose) so an
    # allele-level fidelity loss isn't silent. Both are no-ops without contigs
    # (the divergence columns stay None / 0).
    for chain in ("alpha", "beta"):
        div_col = f"{chain}_allele_divergence_positions"
        if div_col in df.columns:
            df[f"{chain}_constant_allele_divergence"] = df[div_col].apply(
                _count_divergence_positions
            )
    _warn_constant_allele_divergence(df)

    # Synthesis-hazard QC (#206). Scores each assembled construct NT for
    # gene-synthesis pitfalls (GC window, homopolymers, long repeats, surviving
    # restriction sites) plus cross-construct hazards (duplicate CDS, α/β swap).
    # No-op when no construct-NT column exists. Warn (regardless of verbose) on
    # any flagged construct so a vendor-rejection risk isn't silent.
    df = add_synthesis_qc(df)

    # Surface dual-α (allelic-inclusion) clones explicitly (#237) — easy to miss
    # when the second α is packed into ``merged_alpha_partners``. ``alpha_count``
    # is the number of distinct α the clone carries; ``dual_alpha`` flags >1.
    if "merged_alpha_partners" in df.columns and len(df):
        def _alpha_partners(val: object) -> list:
            s = "" if val is None or (isinstance(val, float) and pd.isna(val)) else str(val)
            return [p.strip() for p in s.split(";") if p.strip() and p.strip().lower() != "nan"]
        partners = df["merged_alpha_partners"].map(_alpha_partners)
        df["alpha_count"] = partners.map(lambda p: len(p) if p else 1)
        df["dual_alpha"] = df["alpha_count"] > 1
        # Explicit second-α identity (#237): the partner α that isn't this row's
        # primary CDR3_alpha (else the 2nd listed), so the dual nature reads
        # without unpacking merged_alpha_partners. None for single-α clones.
        prim = df["CDR3_alpha"] if "CDR3_alpha" in df.columns else pd.Series([None] * len(df), index=df.index)
        df["CDR3_alpha_2"] = [
            next((a for a in ps if a != p), (ps[1] if len(ps) > 1 else None))
            for ps, p in zip(partners, prim)
        ]

    # Per-construct contig-verification provenance (#243/#244): a construct is
    # contig-verified only when every present chain's J→C junction residue was
    # read from the donor's contig (not the canonical fallback). Drives the
    # fail-closed fidelity gate.
    if len(df):
        df["construct_contig_verified"] = [
            not _unverified_chains(row) for _, row in df.iterrows()
        ]

    if "synth_gc_ok" in df.columns and len(df):
        n_hazard = int(
            (~df["synth_gc_ok"]).sum()
            + (df.get("synth_max_homopolymer", pd.Series(dtype=int)) >= 9).sum()
            + (df.get("synth_max_repeat", pd.Series(dtype=int)) >= 20).sum()
            + (df.get("synth_restriction_sites", pd.Series(dtype=str)) != "").sum()
        )
        if n_hazard:
            logger.warning(
                "  Synthesis-hazard QC flagged %d construct-level issue(s) "
                "(GC window / homopolymer / repeat / restriction site). See the "
                "synth_* columns and synthesis_qc_report() for details.",
                n_hazard,
            )

    # Framework (FR1–FR3) germline divergence (#286 follow-up). Compares each
    # clone's V-domain framework to the IMGT V-REGION germline and records the
    # divergence for the `framework` rows of collect_germline_variants. No-op
    # (adds no columns) unless the V-REGION reference is cached
    # (`tcrsift data download --db imgt_trav_vregion imgt_trbv_vregion`).
    df = annotate_vregion_divergence(df, species=vdj_sp, data_dir=data_dir)

    # Schema stability (#278). `qc_warnings` is stashed per-clone only when a
    # warning fires, so the column is absent on a clean cohort and present on a
    # noisy one — 144 vs 145 columns across donors, which breaks naive
    # `df["qc_warnings"]` consumers and cross-cohort concats. Always emit it,
    # with an empty list per row where no warning was raised (a list, not NaN,
    # so the surfacing loop's `isinstance(qcs, list)` and callers' `... or []`
    # both stay correct).
    df = _ensure_qc_warnings_column(df)

    return df


def annotate_vregion_divergence(
    df: pd.DataFrame, *, species: str = "human", data_dir=None
) -> pd.DataFrame:
    """Add ``{chain}_vregion_*`` framework-divergence columns (#286 follow-up)."""
    from ..species import resolve_species
    from ..vregion import germline_compare_vregion, reference_is_available

    sp = resolve_species(species, df=df, warn_mismatch=False)
    if not len(df) or not reference_is_available(data_dir, species=sp):
        return df
    for chain in ("alpha", "beta"):
        vdj_col, cdr_col, v_col = (
            f"VDJ_{chain}_aa", f"CDR3_{chain}", f"{chain}_v_gene"
        )
        if not {vdj_col, cdr_col, v_col} <= set(df.columns):
            continue
        alleles, germs, donors, idents, diffs = [], [], [], [], []
        for _, r in df.iterrows():
            res = germline_compare_vregion(
                r.get(vdj_col), r.get(cdr_col), r.get(v_col), species=sp, data_dir=data_dir
            )
            if res is None:
                for col in (alleles, germs, donors, idents, diffs):
                    col.append(None)
                continue
            allele, g_aa, identity, diff = res
            gene = normalize_vgene(r.get(v_col))
            alleles.append(f"{gene}*{allele}" if gene else allele)
            germs.append(g_aa)
            donor_fr = str(r.get(vdj_col))
            cdr = str(r.get(cdr_col))
            donors.append(donor_fr[: donor_fr.index(cdr)] if cdr in donor_fr else None)
            idents.append(identity)
            diffs.append(diff)
        df[f"{chain}_vregion_allele"] = alleles
        df[f"{chain}_vregion_germline_aa"] = germs
        df[f"{chain}_vregion_donor_aa"] = donors
        df[f"{chain}_vregion_identity"] = idents
        df[f"{chain}_vregion_vs_germline"] = diffs
    return df


def _assemble_clone(
    row: pd.Series,
    sample_contigs: dict,
    constant_seqs: dict,
    leader_config: dict,
    include_constant: bool,
    constant_source: str,
    allele_overrides: dict[str, str] | None = None,
    stops_nt: str = "TAATGA",
    recode: CodonOptimization | None = None,
    *,
    vdj_species: str = "human",
    constant_species: str = "human",
    codon_species: str = "human",
) -> dict:
    """Assemble full sequence for a single clone."""
    result: dict = {}
    if isinstance(row.get("qc_warnings"), list) and row["qc_warnings"]:
        result["qc_warnings"] = list(row["qc_warnings"])

    for chain in ["alpha", "beta"]:
        vdj_col = f"VDJ_{chain}_aa"
        vdj_nt_col = f"VDJ_{chain}_nt"

        if vdj_col in row and pd.notna(row.get(vdj_col)):
            result[f"vdj_{chain}_aa"] = str(row[vdj_col]).strip()
        if vdj_nt_col in row and pd.notna(row.get(vdj_nt_col)):
            vdj_nt = str(row[vdj_nt_col]).strip()
            vdj_aa = result.get(f"vdj_{chain}_aa")
            if isinstance(vdj_aa, str) and vdj_aa:
                vdj_nt = vdj_nt[: 3 * len(vdj_aa)]
            result[f"vdj_{chain}_nt"] = vdj_nt

        for gene in ("c_gene", "j_gene"):
            col = f"{chain}_{gene}"
            if col in row:
                result[col] = row[col]

    for chain in ["alpha", "beta"]:
        # Skip leader assembly when this chain has no VDJ sequence (unpaired clones).
        if not result.get(f"vdj_{chain}_aa"):
            continue
        chain_leader = leader_config.get(chain)
        if chain_leader is None:
            continue
        elif chain_leader == "from_contig":
            _extract_leader_from_contigs_single(
                row, sample_contigs, result, chain, species=vdj_species
            )
            if not result.get(f"{chain}_leader_aa"):
                result[f"{chain}_leader_source"] = "missing_from_contig"
        elif isinstance(chain_leader, dict):
            _record_leader(
                result, chain, chain_leader["aa"], chain_leader["dna"],
                "curated_default", v_gene=row.get(f"{chain}_v_gene"),
                species=vdj_species,
            )

    if include_constant:
        if constant_source == "from-data":
            _add_constant_from_row(
                row,
                result,
                species=constant_species,
                codon_species=codon_species,
                stops_nt=stops_nt,
            )
        else:
            _add_constant_regions(
                result, constant_seqs,
                row=row, sample_contigs=sample_contigs,
                allele_overrides=allele_overrides,
                stops_nt=stops_nt,
                vdj_species=vdj_species,
                constant_species=constant_species,
                codon_species=codon_species,
                recode=recode,
            )

    include_alpha_leader = leader_config.get("alpha") is not None
    include_beta_leader = leader_config.get("beta") is not None

    _build_full_sequences(
        result, include_alpha_leader, include_beta_leader, include_constant, recode=recode,
    )

    return result


def _build_full_sequences(
    result: dict,
    include_alpha_leader: bool,
    include_beta_leader: bool,
    include_constant: bool,
    recode: CodonOptimization | None = None,
):
    """Build complete sequences from parts.

    Emits the standard ``full_{chain}_aa`` / ``full_{chain}_nt`` plus
    (when constants are included) the NT triad introduced in #116:

    * ``full_{chain}_nt_contig`` — leader + VDJ + pure CellRanger
      contig bytes for the constant. **Truncated at contig coverage**
      (donor's contig typically reads ~30-150 nt into the C region,
      not the full ~430-530 nt). ``None`` when no contig is
      available — there is no canonical fallback for this view by
      design.
    * ``full_{chain}_nt_optimized`` — the synthesis-ready CDS. By default
      (``recode.mode == "focal"``) the leader+VDJ keep the donor's native
      codons with only minimal synonymous swaps to clear constraint violations
      (#293, :func:`recode_codons`), plus the codon-optimized canonical constant
      (+stops); ``mode == "full"`` recodes the whole ORF with
      :func:`optimize_codons`. The donor's untouched bytes stay in
      ``full_{chain}_nt_contig`` / ``full_{chain}_nt``.
    * ``full_{chain}_nt`` (legacy / blend) — leader + VDJ +
      donor bytes where they agree with canonical AA, ``_optimized``
      for the rest.
    * ``{chain}_vdj_recoded_positions`` / ``{chain}_leader_recoded_positions`` —
      which VDJ / leader codons the focal recode synonymously changed
      (``""`` when none / in full mode).
    """
    recode = recode or CodonOptimization()
    sites = recode.restriction_sites()
    include_leader_map = {"alpha": include_alpha_leader, "beta": include_beta_leader}

    for chain in ["alpha", "beta"]:
        parts_aa: list[str] = []
        parts_nt_blend: list[str] = []
        parts_nt_contig: list[str | None] = []
        # native leader+VDJ bytes that feed the focal-recode prefix
        opt_leader_nt = ""
        opt_vdj_nt = ""
        include_leader = include_leader_map[chain]
        result[f"{chain}_vdj_recoded_positions"] = ""
        result[f"{chain}_leader_recoded_positions"] = ""
        if not result.get(f"vdj_{chain}_aa"):
            continue

        if include_leader and f"{chain}_leader_aa" in result:
            parts_aa.append(result[f"{chain}_leader_aa"])
        if include_leader and f"{chain}_leader_nt" in result:
            leader_nt = result[f"{chain}_leader_nt"]
            parts_nt_blend.append(leader_nt)
            parts_nt_contig.append(leader_nt)
            opt_leader_nt = leader_nt or ""

        if f"vdj_{chain}_aa" in result:
            parts_aa.append(result[f"vdj_{chain}_aa"])
        if f"vdj_{chain}_nt" in result:
            vdj_nt = result[f"vdj_{chain}_nt"]
            parts_nt_blend.append(vdj_nt)
            parts_nt_contig.append(vdj_nt)
            opt_vdj_nt = vdj_nt or ""

        if include_constant and f"{chain}_constant_aa" in result:
            parts_aa.append(result[f"{chain}_constant_aa"])
        if include_constant and f"{chain}_constant_nt" in result:
            parts_nt_blend.append(result[f"{chain}_constant_nt"])
            # _contig is None when contig had no post-VDJ coverage —
            # propagate the None to full_{chain}_nt_contig so callers
            # can filter on it without parsing.
            parts_nt_contig.append(result.get(f"{chain}_constant_nt_contig"))

        if parts_aa:
            result[f"full_{chain}_aa"] = "".join(parts_aa)
        if parts_nt_blend:
            result[f"full_{chain}_nt"] = "".join(parts_nt_blend)
        if include_constant:
            # Contig view: any None constant means the donor's bytes
            # don't reach the C region; surface None instead of a
            # silently-truncated string.
            if all(p is not None for p in parts_nt_contig) and parts_nt_contig:
                result[f"full_{chain}_nt_contig"] = "".join(parts_nt_contig)
            else:
                result[f"full_{chain}_nt_contig"] = None
            if parts_aa:
                # The constant is the canonical-optimized bytes (+stops) in BOTH
                # modes — the same const_opt the standalone {chain}_constant_nt_
                # optimized column carries, so full_*_nt_optimized's C region
                # never diverges from it (review #294). Only the donor leader+VDJ
                # is (focal) recoded or (full) re-optimized.
                const_opt = result.get(f"{chain}_constant_nt_optimized") or ""
                const_aa = result.get(f"{chain}_constant_aa") or ""
                full_aa = result[f"full_{chain}_aa"]
                donor_aa = full_aa[: len(full_aa) - len(const_aa)] if const_aa else full_aa
                if recode.mode == "full":
                    full_forbidden = (
                        tuple(dict.fromkeys((*DEFAULT_FORBIDDEN_MOTIFS, *sites)))
                        if sites
                        else DEFAULT_FORBIDDEN_MOTIFS
                    )
                    donor_opt = optimize_codons(
                        donor_aa,
                        species=recode.codon_species,
                        forbidden_motifs=full_forbidden,
                    )
                    if sites or (recode.max_codon_repeats and recode.max_codon_repeats >= 1):
                        donor_opt, _ = recode_codons(
                            donor_opt,
                            max_codon_repeats=recode.max_codon_repeats,
                            restriction_sites=sites,
                            species=recode.codon_species,
                            right_context=const_opt,
                        )
                    result[f"full_{chain}_nt_optimized"] = donor_opt + const_opt
                else:
                    eff_leader_nt = opt_leader_nt or (
                        optimize_codons(
                            result[f"{chain}_leader_aa"], species=recode.codon_species
                        )
                        if include_leader and result.get(f"{chain}_leader_aa")
                        else ""
                    )
                    eff_vdj_nt = opt_vdj_nt or (
                        optimize_codons(
                            result[f"vdj_{chain}_aa"], species=recode.codon_species
                        )
                        if result.get(f"vdj_{chain}_aa")
                        else ""
                    )
                    native_prefix = (eff_leader_nt + eff_vdj_nt).upper()
                    recoded_prefix, changed = recode_codons(
                        native_prefix,
                        max_codon_repeats=recode.max_codon_repeats,
                        restriction_sites=sites,
                        species=recode.codon_species,
                        right_context=const_opt,
                    )
                    result[f"full_{chain}_nt_optimized"] = recoded_prefix + const_opt
                    # Record recodes by region — leader vs VDJ — so NONE are
                    # silently dropped (review #294). Offsets are construct-relative.
                    vdj_start = len(eff_leader_nt)

                    def _fmt(o, _np=native_prefix, _rp=recoded_prefix):
                        return f"{o}:{_np[o:o + 3]}>{_rp[o:o + 3]}"

                    result[f"{chain}_leader_recoded_positions"] = ";".join(
                        _fmt(o) for o in changed if o < vdj_start
                    )
                    result[f"{chain}_vdj_recoded_positions"] = ";".join(
                        _fmt(o) for o in changed if o >= vdj_start
                    )


_GERMLINE_VARIANT_COLS = [
    "region", "chain", "gene", "germline_allele", "donor_seq", "germline_seq",
    "variant", "identity", "source", "shipped",
]


def collect_germline_variants(
    df: pd.DataFrame, *, chains=("alpha", "beta")
) -> pd.DataFrame:
    """Collect donor-vs-germline divergences from assembled output, by region (#270/#276).

    Returns one row per DISTINCT variant with a ``region`` tag and ``n_constructs``
    (how many constructs carry it), sorted by ``n_constructs`` descending. Empty
    frame when nothing is germline-comparable. Regions:

    * ``leader`` — signal-peptide divergence from the IMGT V-leader. ``shipped``
      is ``kept`` (the donor variant IS the shipped leader) or ``switched`` (the
      policy replaced it with germline/curated; recovered from the
      ``{chain}_preswitch_*`` stash so it isn't lost, #276).
    * ``constant`` — donor contig constant diverging from the canonical TRAC/TRBC
      (the #187 ``{chain}_allele_divergence_positions``). ``shipped`` is
      ``canonical`` — tcrsift ships the canonical constant for TCR-T, so the donor
      variant is *detected and flagged but not shipped*.
    * ``framework`` — donor V-domain (FR1–FR3) diverging from the germline
      V-REGION, when ``{chain}_vregion_vs_germline`` is present (needs the bundled
      IMGT V-REGION reference). Forward-compatible: picked up automatically once
      that column exists.

    Columns: ``region, chain, gene, germline_allele, donor_seq, germline_seq,
    variant, identity, source, shipped, n_constructs``.
    """
    rows: list[dict] = []

    def _add(region, ch, *, gene, allele, donor, germ, variant, identity, source, shipped):
        rows.append({
            "region": region, "chain": ch, "gene": gene, "germline_allele": allele,
            "donor_seq": donor, "germline_seq": germ, "variant": variant,
            "identity": identity, "source": source, "shipped": shipped,
        })

    # Exclude secondary twin rows so donor variants are counted once per primary construct.
    base_df = (
        df[~df["leader_variant"].fillna("").astype(str).str.startswith("secondary:")]
        if "leader_variant" in df.columns
        else df
    )

    for ch in chains:
        # --- leader (kept donor contig leader only — exclude exogenous curated leaders) ---
        vs = f"{ch}_leader_vs_germline"
        if vs in base_df.columns:
            d = base_df[base_df[vs].notna() & ~base_df[vs].astype(str).isin(["identical"])]
            for _, r in d.iterrows():
                lsrc = r.get(f"{ch}_leader_source")
                if lsrc is not None and lsrc not in _FROM_CONTIG_SOURCES:
                    continue
                _add("leader", ch, gene=r.get(f"{ch}_v_gene"),
                     allele=r.get(f"{ch}_germline_allele"), donor=r.get(f"{ch}_leader_aa"),
                     germ=r.get(f"{ch}_germline_leader_aa"), variant=r.get(vs),
                     identity=r.get(f"{ch}_germline_identity"),
                     source=lsrc, shipped="kept")
        # --- leader (switched-away, preserved pre-switch, #276) ---
        pv = f"{ch}_preswitch_variant"
        if pv in base_df.columns:
            d = base_df[base_df[pv].notna() & ~base_df[pv].astype(str).isin(["identical"])]
            for _, r in d.iterrows():
                psrc = r.get(f"{ch}_preswitch_leader_source")
                if psrc is not None and psrc not in _FROM_CONTIG_SOURCES:
                    continue
                _add("leader", ch, gene=r.get(f"{ch}_v_gene"),
                     allele=r.get(f"{ch}_preswitch_germline_allele"),
                     donor=r.get(f"{ch}_preswitch_donor_aa"),
                     germ=r.get(f"{ch}_preswitch_germline_aa"), variant=r.get(pv),
                     identity=r.get(f"{ch}_preswitch_germline_identity"),
                     source=psrc, shipped="switched")
        # --- constant (detected donor-vs-canonical divergence, #187; canonical shipped) ---
        cpos = f"{ch}_allele_divergence_positions"
        if cpos in base_df.columns:
            d = base_df[base_df[cpos].notna() & ~base_df[cpos].astype(str).isin(["", "nan", "None"])]
            for _, r in d.iterrows():
                _add("constant", ch, gene=_resolve_c_gene(r, ch),
                     allele=_resolve_c_gene(r, ch), donor=None, germ=None,
                     variant=r.get(cpos), identity=None, source="contig",
                     shipped="canonical")
        # --- framework / V-REGION (forward-compatible; needs the V-REGION ref) ---
        fvs = f"{ch}_vregion_vs_germline"
        if fvs in base_df.columns:
            d = base_df[base_df[fvs].notna() & ~base_df[fvs].astype(str).isin(["identical"])]
            for _, r in d.iterrows():
                _add("framework", ch, gene=r.get(f"{ch}_v_gene"),
                     allele=r.get(f"{ch}_vregion_allele"),
                     donor=r.get(f"{ch}_vregion_donor_aa"),
                     germ=r.get(f"{ch}_vregion_germline_aa"), variant=r.get(fvs),
                     identity=r.get(f"{ch}_vregion_identity"),
                     # framework comes from the donor contig V(D)J and is shipped
                     # as-is (no keep/switch policy): source is the contig, and
                     # shipped="donor" (the donor variant IS in the construct).
                     source="contig", shipped="donor")

    if not rows:
        return pd.DataFrame(columns=[*_GERMLINE_VARIANT_COLS, "n_constructs"])
    out = pd.DataFrame(rows)
    grouped = (
        out.groupby(_GERMLINE_VARIANT_COLS, dropna=False).size()
        .reset_index(name="n_constructs")
        .sort_values(["region", "n_constructs"], ascending=[True, False], ignore_index=True)
    )
    return grouped

