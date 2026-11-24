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
Full-length TCR sequence assembly for TCRsift.

Builds complete TCR sequences including leader peptides and constant regions.
"""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from ._assemble_codons import (  # noqa: F401
    CODON_ALTERNATIVES_BY_SPECIES,
    CODON_TABLE,
    DEFAULT_FORBIDDEN_MOTIFS,
    DEFAULT_GC_TARGET_RANGE,
    DEFAULT_GC_TARGET_WINDOW,
    HUMAN_CODON_ALTERNATIVES,
    HUMAN_PREFERRED_CODONS,
    LINKERS,
    MOUSE_CODON_ALTERNATIVES,
    MOUSE_PREFERRED_CODONS,
    P2A_LINKER_AA,
    P2A_LINKER_DNA,
    PREFERRED_CODONS_BY_SPECIES,
    SAMPLE_NAME_FROM_CHOICES,
    SYNTHESIS_RESTRICTION_SITES,
    T2A_LINKER_AA,
    T2A_LINKER_DNA,
    CodonOptimization,
    _revcomp,
    back_translate,
    find_longest_orf,
    load_contigs,
    optimize_codons,
    parse_contig_json,
    parse_fasta,
    recode_codons,
    stop_codons_nt,
    translate_dna,
)
from ._assemble_constants import (  # noqa: F401
    _CANONICAL_ALLELES_FASTA,
    _CANONICAL_ALLELES_MOUSE_FASTA,
    _CANONICAL_CONSTANTS_FASTA,
    _CANONICAL_CONSTANTS_MOUSE_FASTA,
    _DEFAULT_LEADERS_CASEFOLD,
    _DEFAULT_MIN_PICKER_POSITIONS,
    _FALLBACK_JUNCTION_RESIDUE,
    _HYDROPHOBIC_AA,
    _MAX_DIVERGENCE_SCAN_POSITIONS,
    _VALID_AA_CHARS,
    _VALID_JUNCTION_RESIDUES,
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_DIVERGENT_CONTIG,
    ALLELE_REASON_FRAME_ERROR,
    ALLELE_REASON_INVALID_OVERRIDE,
    ALLELE_REASON_NO_ALLELE_POOL,
    ALLELE_REASON_NO_CONTIG,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_SPARSE_CONTIG,
    ALLELE_REASON_VALUES,
    BETA_JC_PARITY,
    CONSTANT_AA_FLOOR,
    CONSTANT_ALLELES_BY_SPECIES,
    CONSTANT_REGION_ENDINGS,
    CONSTANT_REGION_ENDINGS_BY_SPECIES,
    CONSTANT_REGIONS_AA_BY_SPECIES,
    DEFAULT_LEADERS,
    HUMAN_CONSTANT_ALLELES,
    HUMAN_CONSTANT_REGIONS_AA,
    HUMAN_TRAC_AA,
    HUMAN_TRBC1_AA,
    HUMAN_TRBC2_AA,
    KOZAK_ADEQUATE,
    MOUSE_CONSTANT_ALLELES,
    MOUSE_CONSTANT_REGION_ENDINGS,
    MOUSE_CONSTANT_REGIONS_AA,
    MOUSE_TRAC_AA,
    MOUSE_TRAJ_JUNCTION_RESIDUES,
    MOUSE_TRBC1_AA,
    MOUSE_TRBC2_AA,
    SP_LENGTH_HARD,
    SP_LENGTH_MAX,
    SP_LENGTH_MIN,
    TRAJ_JUNCTION_RESIDUES,
    _apply_cohort_alpha_junctions,
    _beta_c_base,
    _beta_j_family,
    _classify_leader,
    _count_divergence_positions,
    _count_trbc_overrides,
    _describe_leader,
    _divergence_positions,
    _ensure_qc_warnings_column,
    _fallback_junction_residue,
    _format_divergence_positions,
    _has_h_region,
    _kozak_score,
    _load_canonical_alleles_fasta,
    _load_canonical_constants_fasta,
    _lookup_default_leader,
    _pick_best_allele,
    _polymorphic_positions,
    _record_leader,
    _resolve_c_gene,
    _score_allele_against_contig,
    _traj_junction_residue,
    _trbc_override_mask,
    _verify_first_codon_frame,
    _warn_constant_allele_divergence,
    leader_qc,
    pick_canonical_constant,
    verify_canonical_constant_start,
)
from ._assemble_germline import (  # noqa: F401
    GERMLINE_IMPUTED_VDJ_SOURCES,
    VDJ_SOURCE_5PRIME_RESCUED,
    VDJ_SOURCE_CONTIG_FULL,
    VDJ_SOURCE_VJ_RECONSTRUCTED,
    reconstruct_germline_vdj,
    reconstruct_germline_vdj_chain,
)
from ._assemble_qc import (  # noqa: F401
    _NOVEL_ALLELE_EMPTY_COLUMNS,
    AssemblyQCCheck,
    AssemblyQCReport,
    ValidationMessage,
    _candidate_allele_label,
    _candidate_constant_alleles_for_chain,
    _closest_reference_for_observations,
    _collect_divergences_for_chain,
    _count_observed_positions_for_chain,
    _default_reference_for_gene,
    _expand_samples_unique,
    _expected_constant_start_from_full,
    _gc_fraction,
    _infer_constant_species,
    _longest_repeat,
    _max_homopolymer,
    _nt_contig_in_frame_no_premature_stop,
    _nt_translates_to,
    _parse_divergence_positions_string,
    _reference_for_gene_or_observations,
    _restriction_hits,
    _terminal_stops,
    _unverified_chains,
    _validate_nt_aa_roundtrip,
    _variant_description,
    add_synthesis_qc,
    allele_audit_report,
    assemble_fidelity_summary,
    assemble_qc_report,
    build_assembly_qc_report,
    detect_novel_alleles,
    enforce_contig_fidelity,
    export_fasta,
    fix_jc_parity,
    synthesis_qc_report,
    validate_sequences,
)
from .leaders import (
    characterize_divergence,
    germline_anchor_leader,
    germline_compare_leader,
    germline_leader_nt,
    germline_vgene_leaders,
    normalize_vgene,
)
from .sigpep import signal_peptide_features, sp_features_summary
from .validation import (
    TCRsiftValidationError,
    validate_clonotype_df,
    validate_directory_exists,
)

logger = logging.getLogger(__name__)


def get_constant_region_sequences(
    *,
    stop_codons: tuple[str, ...] = ("TAA", "TGA"),
    species: str = "human",
    codon_species: str | None = None,
) -> dict[str, str]:
    """Return TCR constant-region CDS (DNA) via codon-aware back-translation."""
    from .species import normalize_species

    const_sp = normalize_species(species)
    cod_sp = normalize_species(codon_species) if codon_species else const_sp
    stops_nt = stop_codons_nt(stop_codons)
    regions = CONSTANT_REGIONS_AA_BY_SPECIES[const_sp]
    out: dict[str, str] = {}
    for name, aa in regions.items():
        out[name] = optimize_codons(aa, species=cod_sp) + stops_nt
    return out


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
    from .species import normalize_species, resolve_species

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
    from .species import resolve_species
    from .vregion import germline_compare_vregion, reference_is_available

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
            result[f"vdj_{chain}_aa"] = row[vdj_col]
        if vdj_nt_col in row and pd.notna(row.get(vdj_nt_col)):
            vdj_nt = row[vdj_nt_col]
            vdj_aa = result.get(f"vdj_{chain}_aa")
            if isinstance(vdj_aa, str) and vdj_aa:
                vdj_nt = vdj_nt[: 3 * len(vdj_aa)]
            result[f"vdj_{chain}_nt"] = vdj_nt

        for gene in ("c_gene", "j_gene"):
            col = f"{chain}_{gene}"
            if col in row:
                result[col] = row[col]

    for chain in ["alpha", "beta"]:
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
            )

    include_alpha_leader = leader_config.get("alpha") is not None
    include_beta_leader = leader_config.get("beta") is not None

    _build_full_sequences(
        result, include_alpha_leader, include_beta_leader, include_constant, recode=recode,
    )

    return result


def _kozak_correct_leader(
    contig_seq: str, offset: int, orig_leader: str
) -> tuple[str, str, int, str]:
    """Choose the leader's start codon, auto-shortening ONLY an over-capture (#263).

    ``orig_leader`` is the ORF prefix before the VDJ (translated from the ATG at
    ``offset``). Policy:

    - **In-range leaders (len ≤ SP_LENGTH_MAX) are never altered.** The ORF
      start is kept. This avoids silently shortening a correct leader whose true
      start happens to have weak Kozak context (a downstream adequate Met must
      NOT win here). Weakness of the kept start is still surfaced as QC by the
      caller (``weak_kozak_start``) — leaky-scan informs, it doesn't act.
    - **Over-captures (len > SP_LENGTH_MAX)** trigger a leaky scan: the in-frame
      Mets are walked 5'→3' and initiation moves to the FIRST one with adequate
      Kozak context (``>= KOZAK_ADEQUATE``) whose resulting leader lands in the
      ``SP_LENGTH_MIN..MAX`` window — what the ribosome does, leaky-scanning past
      weak upstream AUGs. If no such start exists the over-capture is kept (the
      caller flags it ``too_long``); we never fabricate a worse leader.

    - **No adequate-Kozak start, but a structurally-valid nested SP exists**
      (over-capture only): fall back to a *best-available STRUCTURAL trim* —
      among in-range nested Mets whose leader carries a hydrophobic h-region
      (the SP's defining feature), take the best-available Kozak (even if weak),
      tie-broken toward the most DOWNSTREAM (shortest, "SP at the end"). This
      catches a real SP whose start is genuinely weak-Kozak — e.g. the TRAV1
      family leader ``MWG.FLLYVSMKMGGT.`` nested at the tail of a 28-34 aa
      over-capture. Flagged ``contig_hregion_trimmed`` so the weaker basis is
      visible. If no such nested SP exists (e.g. a genuinely long ~29 aa leader
      with no shorter Met, like TRBV13), the ORF start is kept and the caller
      bands it by length — ``long_leader`` (26-30 aa, accepted) or ``too_long``
      (> SP_LENGTH_HARD); we never fabricate one.

    Returns ``(leader_aa, leader_nt, kozak_score_of_chosen_start, source)`` where
    ``source`` is ``"contig"`` (ORF start kept), ``"contig_kozak_reselected"``
    (over-capture → in-range adequate-Kozak start), or
    ``"contig_hregion_trimmed"`` (over-capture → in-range h-region SP, weak Kozak).
    """
    orig_len = len(orig_leader)
    # In-frame Met starts (aa offsets into the ORF); nt position is offset+3*k.
    atgs = [k for k in range(orig_len) if orig_leader[k] == "M"]

    chosen_k, source = 0, "contig"
    if orig_len > SP_LENGTH_MAX:
        # 1. Leaky-scan for the first (most upstream) adequate Met whose leader
        #    lands in-range. atgs ascend in k, so the resulting length descends;
        #    once it drops below the floor, no later Met can qualify → stop.
        for k in atgs:
            if k == 0:
                continue  # can't shorten the over-capture to itself
            new_len = orig_len - k
            if new_len < SP_LENGTH_MIN:
                break
            if new_len <= SP_LENGTH_MAX and (
                _kozak_score(contig_seq, offset + 3 * k) >= KOZAK_ADEQUATE
            ):
                chosen_k, source = k, "contig_kozak_reselected"
                break
        # 2. No adequate-Kozak start → best-available STRUCTURAL trim: an
        #    in-range nested Met WITH a hydrophobic h-region (a real SP whose
        #    start is weak-Kozak). Best-available Kozak, tie-break most
        #    downstream (the SP sits at the tail of the over-capture).
        if source == "contig":
            h_cands = [
                k
                for k in atgs
                if k != 0
                and SP_LENGTH_MIN <= orig_len - k <= SP_LENGTH_MAX
                and _has_h_region(orig_leader[k:])
            ]
            if h_cands:
                chosen_k = max(
                    h_cands,
                    key=lambda k: (_kozak_score(contig_seq, offset + 3 * k), k),
                )
                source = "contig_hregion_trimmed"

    leader_aa = orig_leader[chosen_k:]
    leader_nt = contig_seq[offset + 3 * chosen_k : offset + orig_len * 3]
    score = _kozak_score(contig_seq, offset + 3 * chosen_k)
    return leader_aa, leader_nt, score, source


def _extract_leader_from_contigs_single(
    row: pd.Series,
    sample_contigs: dict,
    result: dict,
    chain: str,
    *,
    species: str = "human",
):
    """Extract + QC the signal peptide from contigs for a single chain (#263).

    Each supporting contig votes its (over-capture-corrected) leader; the
    consensus (most-common) leader is recorded with QC provenance:

    - ``{chain}_leader_qc`` — ``ok`` / ``weak_kozak_start`` (in-range, kept, but
      its start has weak Kozak context) / ``too_long`` / ``too_short`` /
      ``no_met`` / ``no_h_region`` / ``internal_stop`` / ``missing``.
    - ``{chain}_leader_source`` — ``contig`` (ORF start kept) or
      ``contig_kozak_reselected`` (over-capture auto-shortened to an in-range
      adequate-Kozak start).
    - ``{chain}_leader_len``, ``_leader_kozak_score`` (of the chosen start),
      ``_leader_support`` (n agreeing / n total).

    When no contig covers the VDJ the leader is simply left unset (no QC fields).
    """
    samples = str(row.get("samples", "")).split(";")

    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return

    contig_ids = str(row[contig_col]).split(";")
    vdj_aa = result.get(f"vdj_{chain}_aa", "")
    v_gene = row.get(f"{chain}_v_gene")

    leader_counter = Counter()
    # leader_aa -> (leader_nt, kozak_score, source). Metadata is consistent across
    # contigs voting the same leader_aa, so first sighting is fine. The germline
    # columns are added later by _record_leader's universal comparison.
    leader_meta: dict = {}
    n_total = 0

    for sample in samples:
        if sample not in sample_contigs:
            continue

        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue

            contig_seq = sample_contigs[sample][contig_id]
            translated, offset, ragged = find_longest_orf(contig_seq)

            if vdj_aa and vdj_aa in translated and offset is not None:
                orig_leader = translated.split(vdj_aa)[0]
                orig_len = len(orig_leader)
                n_total += 1
                # Layer 1: germline anchor (#267) decides the START — keep the
                # donor's native suffix whose length/start match the germline
                # V-leader. Falls back (None) when the gene/allele is absent or
                # the contig is too divergent → Kozak + h-region heuristic (#263).
                anchored = germline_anchor_leader(
                    orig_leader, v_gene, species=species
                )
                if anchored is not None:
                    leader_aa = anchored[0]
                    start_k = orig_len - len(leader_aa)
                    leader_nt = contig_seq[offset + 3 * start_k : offset + orig_len * 3]
                    score = _kozak_score(contig_seq, offset + 3 * start_k)
                    source = "contig_germline_anchored"
                else:
                    leader_aa, leader_nt, score, source = _kozak_correct_leader(
                        contig_seq, offset, orig_leader
                    )
                leader_counter[leader_aa] += 1
                leader_meta.setdefault(leader_aa, (leader_nt, score, source))

    if not leader_counter:
        return  # no contig covered the VDJ — leader simply absent

    leader_aa, support = leader_counter.most_common(1)[0]
    leader_nt, score, source = leader_meta[leader_aa]
    _record_leader(
        result, chain, leader_aa, leader_nt, source,
        v_gene=v_gene, kozak_score=score, support=f"{support}/{n_total}",
        species=species,
    )


def _add_constant_regions(
    result: dict,
    constant_seqs: dict,
    row: pd.Series | None = None,
    sample_contigs: dict | None = None,
    allele_overrides: dict[str, str] | None = None,
    stops_nt: str = "TAATGA",
    *,
    vdj_species: str = "human",
    constant_species: str = "human",
    codon_species: str = "human",
):
    """Add constant-region sequences and verify against the observed
    CellRanger contig start where possible (including chimeric constructs)."""
    allele_overrides = allele_overrides or {}
    is_chimeric = constant_species != vdj_species
    species_allele_pools = CONSTANT_ALLELES_BY_SPECIES.get(
        constant_species, HUMAN_CONSTANT_ALLELES
    )
    for chain in ["alpha", "beta"]:
        c_gene = result.get(f"{chain}_c_gene")
        j_gene = result.get(f"{chain}_j_gene")
        canonical_name, canonical_aa = pick_canonical_constant(
            chain, c_gene, j_gene, species=constant_species
        )

        canonical_nt_codon_opt = constant_seqs.get(
            canonical_name,
            optimize_codons(canonical_aa, species=codon_species) + stops_nt,
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name

        contig_nt_past_vdj = (
            _extract_c_region_nt_from_contig(
                row, sample_contigs, result.get(f"vdj_{chain}_nt", ""), chain
            )
            if (row is not None and sample_contigs)
            else None
        )

        allele_choice = allele_overrides.get(canonical_name, "auto")
        allele_pool = species_allele_pools.get(canonical_name, {})
        allele_called: str | None = None
        allele_score: float | None = None
        allele_alternatives: dict[str, float] = {}
        allele_called_reason: str | None = None
        divergence_ref_aa: str | None = None

        contig_post_junction_nt = (
            contig_nt_past_vdj[3:] if contig_nt_past_vdj else ""
        )
        contig_aa_past_junction, _ = (
            translate_dna(contig_post_junction_nt)
            if contig_post_junction_nt else ("", False)
        )

        if not allele_pool:
            allele_called_reason = ALLELE_REASON_NO_ALLELE_POOL
        elif allele_choice != "auto":
            if allele_choice in allele_pool:
                canonical_aa = allele_pool[allele_choice]
                canonical_nt_codon_opt = (
                    optimize_codons(canonical_aa, species=codon_species) + stops_nt
                )
                allele_called = allele_choice
                allele_score = 1.0
                allele_called_reason = ALLELE_REASON_OVERRIDDEN
                if not is_chimeric:
                    divergence_ref_aa = canonical_aa
            else:
                allele_called_reason = ALLELE_REASON_INVALID_OVERRIDE
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: invalid allele override "
                    f"{canonical_name}*{allele_choice!s} (available: "
                    f"{sorted(allele_pool.keys())}); falling back to "
                    f"default allele."
                )
        elif is_chimeric:
            # Heterologous constant region (e.g. human VDJ + mouse constant):
            # donor contigs encode the donor's native constant, not the chimeric
            # constant, so commit the canonical chimeric allele directly.
            default_allele = next(iter(allele_pool.keys()))
            canonical_aa = allele_pool[default_allele]
            canonical_nt_codon_opt = (
                optimize_codons(canonical_aa, species=codon_species) + stops_nt
            )
            allele_called = default_allele
            allele_score = 1.0
            allele_called_reason = ALLELE_REASON_OVERRIDDEN
        elif not contig_nt_past_vdj:
            allele_called_reason = ALLELE_REASON_NO_CONTIG
            if sample_contigs:
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — no contig coverage "
                    f"past VDJ for this clone."
                )
        else:
            best_allele, best_score, all_scores = _pick_best_allele(
                contig_aa_past_junction, canonical_name, allele_pool,
            )
            if best_allele is None:
                n_codons = len(contig_aa_past_junction) if contig_aa_past_junction else 0
                allele_called_reason = ALLELE_REASON_SPARSE_CONTIG
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — contig covers only "
                    f"{n_codons} codons past VDJ; need "
                    f"≥{_DEFAULT_MIN_PICKER_POSITIONS} for confident "
                    f"detection."
                )
            elif best_score <= 0:
                allele_called_reason = ALLELE_REASON_DIVERGENT_CONTIG
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: allele not called — contig disagrees with "
                    f"every packaged allele at every compared position."
                )
            else:
                candidate_aa = allele_pool[best_allele]
                if not _verify_first_codon_frame(
                    contig_post_junction_nt, candidate_aa,
                ):
                    allele_called_reason = ALLELE_REASON_FRAME_ERROR
                    result.setdefault("qc_warnings", []).append(
                        f"{chain}: allele picker selected "
                        f"{canonical_name}*{best_allele} (AA score "
                        f"{best_score:.3f}) but the contig's first "
                        "codon past the junction does not translate "
                        "to that allele's expected first residue — "
                        "possible frame error upstream. Falling "
                        "back to default."
                    )
                else:
                    divergence_ref_aa = candidate_aa
                    polymorphic = _polymorphic_positions(allele_pool)
                    disagreed_at_polymorphic = sorted(
                        i for i in polymorphic
                        if i < len(contig_aa_past_junction)
                        and i < len(candidate_aa)
                        and contig_aa_past_junction[i] != candidate_aa[i]
                    )
                    if disagreed_at_polymorphic:
                        allele_called_reason = (
                            ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION
                        )
                        positions_str = ", ".join(
                            f"{i + 1}({candidate_aa[i]}/"
                            f"{contig_aa_past_junction[i]})"
                            for i in disagreed_at_polymorphic
                        )
                        result.setdefault("qc_warnings", []).append(
                            f"{chain}: allele not called — contig "
                            f"disagrees with best-fit "
                            f"{canonical_name}*{best_allele} at "
                            f"allele-distinguishing position(s) "
                            f"{positions_str} (canonical AA / observed AA). "
                            f"Possible novel allele; falling back to default."
                        )
                        allele_alternatives = all_scores
                    else:
                        canonical_aa = candidate_aa
                        canonical_nt_codon_opt = (
                            optimize_codons(canonical_aa, species=codon_species)
                            + stops_nt
                        )
                        allele_called = best_allele
                        allele_score = best_score
                        allele_alternatives = all_scores
                        allele_called_reason = ALLELE_REASON_AUTO_DETECTED

        junction_residue: str | None = None
        junction_from_contig = False
        fallback_source = "canonical_fallback"
        fallback_reason: str | None = None
        if contig_nt_past_vdj and len(contig_nt_past_vdj) >= 3:
            donor_c_aa = CONSTANT_REGIONS_AA_BY_SPECIES.get(
                vdj_species, HUMAN_CONSTANT_REGIONS_AA
            ).get(canonical_name, "")
            if (
                is_chimeric
                and len(contig_post_junction_nt) >= 3
                and donor_c_aa
                and not _verify_first_codon_frame(contig_post_junction_nt, donor_c_aa)
            ):
                fallback_reason = (
                    "contig first post-junction codon does not match donor "
                    f"{vdj_species} {canonical_name} reading frame"
                )
                contig_nt_past_vdj = None
            else:
                junction_codon = contig_nt_past_vdj[:3]
                translated = CODON_TABLE.get(junction_codon, "X")
                if translated and translated not in {"*", "X"}:
                    junction_residue = translated
                    junction_from_contig = True
                else:
                    fallback_reason = (
                        f"contig junction codon {junction_codon!r} translates to "
                        f"{translated!r}"
                    )
                    contig_nt_past_vdj = None
        elif sample_contigs:
            fallback_reason = "no contig C-region coverage past VDJ"
        if junction_residue is None:
            junction_residue, fallback_source = _fallback_junction_residue(
                chain, result.get(f"{chain}_j_gene"), species=vdj_species
            )
            junction_residue = junction_residue or None
            if junction_residue is not None and fallback_reason is not None:
                if fallback_source == "j_inferred":
                    note = (f" — inferred from {result.get(f'{chain}_j_gene')!r} "
                            "(germline per-TRAJ residue)")
                elif chain == "alpha":
                    note = (" — α J gene unmapped; used most-common N "
                            "(low-confidence)")
                else:
                    note = " — β junction is invariantly E"
                result.setdefault("qc_warnings", []).append(
                    f"{chain}: J→C junction residue {junction_residue!r} from "
                    f"{fallback_source} ({fallback_reason}){note}."
                )
        if junction_residue:
            canonical_aa = junction_residue + canonical_aa
            canonical_nt_codon_opt = (
                back_translate(junction_residue, species=codon_species)
                + canonical_nt_codon_opt
            )
        result[f"{chain}_junction_residue"] = junction_residue
        result[f"{chain}_junction_residue_source"] = (
            "contig" if junction_from_contig else fallback_source
        )

        result[f"{chain}_allele_called"] = (
            f"{canonical_name}*{allele_called}" if allele_called else None
        )
        result[f"{chain}_allele_score"] = allele_score
        result[f"{chain}_allele_called_reason"] = allele_called_reason
        if allele_alternatives:
            result[f"{chain}_allele_alternatives"] = ";".join(
                f"{canonical_name}*{a}:{s:.3f}"
                for a, s in sorted(
                    allele_alternatives.items(), key=lambda kv: -kv[1]
                )
            )
        else:
            result[f"{chain}_allele_alternatives"] = None

        result[f"{chain}_observed_constant_aa_start"] = (
            contig_aa_past_junction[:15]
            if (contig_aa_past_junction and not is_chimeric)
            else None
        )

        if divergence_ref_aa is not None and contig_aa_past_junction and not is_chimeric:
            divs = _divergence_positions(divergence_ref_aa, contig_aa_past_junction)
            result[f"{chain}_allele_divergence_positions"] = (
                _format_divergence_positions(divs)
            )
        else:
            result[f"{chain}_allele_divergence_positions"] = None

        result[f"{chain}_constant_aa"] = canonical_aa
        contig_for_blend = (
            (contig_nt_past_vdj[:3] if (contig_nt_past_vdj and len(contig_nt_past_vdj) >= 3) else "")
            if is_chimeric
            else (contig_nt_past_vdj or "")
        )
        blended_nt, blend_debug = _blend_constant_nt_with_contig(
            contig_for_blend,
            canonical_aa,
            canonical_nt_codon_opt,
            codon_species=codon_species,
        )

        result[f"{chain}_constant_nt_contig"] = (
            contig_nt_past_vdj if (contig_nt_past_vdj and not is_chimeric) else None
        )
        result[f"{chain}_constant_nt_optimized"] = canonical_nt_codon_opt
        result[f"{chain}_constant_nt"] = blended_nt

        observed = (
            _extract_c_region_start_from_contig(
                row, sample_contigs, result.get(f"vdj_{chain}_aa", ""), chain
            )
            if (row is not None and sample_contigs and not is_chimeric)
            else None
        )
        if is_chimeric:
            result[f"{chain}_constant_source"] = (
                f"chimeric:{constant_species}:{canonical_name} [{blend_debug['source']}]"
            )
        elif observed is None:
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} [{blend_debug['source']}]"
            )
        elif verify_canonical_constant_start(observed, canonical_aa):
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (contig-verified) "
                f"[{blend_debug['source']}]"
            )
        else:
            result[f"{chain}_constant_source"] = (
                f"canonical:{canonical_name} (UNVERIFIED — start mismatch) "
                f"[{blend_debug['source']}]"
            )
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant start mismatch: observed "
                f"{observed!r} differs from canonical {canonical_name} "
                f"(expected start {canonical_aa[:len(observed)]!r}). "
                "Using canonical anyway; c_gene assignment may be wrong."
            )
        if not is_chimeric and blend_debug["aa_mismatch_at"] is not None:
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant: contig diverged from canonical at AA "
                f"position {blend_debug['aa_mismatch_at']}; switched to "
                f"codon-optimized canonical from that point."
            )
        if not is_chimeric and blend_debug.get("partial_codon_dropped"):
            result.setdefault("qc_warnings", []).append(
                f"{chain} constant: contig provided a partial codon at the "
                "C-region boundary, but its bytes were incompatible with "
                "the canonical residue; dropped donor fidelity at that "
                "codon and used canonical."
            )


def _extract_c_region_nt_from_contig(
    row: pd.Series,
    sample_contigs: dict,
    vdj_nt: str,
    chain: str,
) -> str | None:
    """Pull NT bytes observed immediately after the VDJ end in the CellRanger
    contig — the donor's real, unedited sequence at the J→C boundary.

    Returns the contig NT starting one nt after the last nt of ``vdj_nt`` in
    the contig. The returned NT is **frame-aligned** to the next codon: in
    the assembled construct, this is the first nt of the J→C junction codon
    (whose first 1-2 nt are biologically encoded by the J's terminal partial
    codon and whose remaining nt come from the C exon's start). The donor's
    actual junction codon and the first several downstream C-region codons
    live here, so :func:`_blend_constant_nt_with_contig` can splice them in
    before falling back to a codon-optimized canonical for the rest.

    Notes
    -----
    * ``vdj_nt`` is the frame-aligned (no +1 overshoot) lowercase
      ``vdj_{chain}_nt`` written by the assembler after #91. The original
      CellRanger contig (``filtered_contig.fasta`` parsed into ``sample_contigs``)
      retains the J→C boundary nt that were trimmed during clonotype-aggregation.
    * Returns ``None`` when no contig contains ``vdj_nt`` exactly — usually a
      sign that the canonical-from-contigs path isn't available for this
      clonotype (e.g., contigs_dir wasn't loaded, or VDJ_nt was re-derived
      from a different cell than the one whose contig we have).
    * When multiple contigs contain ``vdj_nt`` (multi-sample clonotype), the
      most-common post-VDJ NT is returned.
    * **Cross-source hazard mitigation (#94 archetype).** ``contig_ids`` is
      the semicolon-joined list of EVERY cell's contig — not just the
      representative cell whose ``vdj_nt`` was kept by the clonotype
      aggregator. We rely on ``contig_seq.find(vdj_nt)`` to filter: contigs
      whose own VDJ doesn't contain the representative ``vdj_nt`` byte-for-
      byte are skipped. Cells with cell-to-cell polymorphism inside the VDJ
      therefore contribute nothing to the post-VDJ consensus, and the
      blended NT stays sourced from a coherent biological signal. Don't
      remove the ``find`` filter without re-deriving the safety argument.
    """
    if not vdj_nt:
        return None
    samples = str(row.get("samples", "")).split(";")
    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return None

    contig_ids = str(row[contig_col]).split(";")
    # Match case-insensitively: the assembler's ``vdj_{chain}_nt`` is lowercase
    # (#91) while CellRanger ``filtered_contig.fasta`` is uppercase. When 5'
    # germline rescue prepended N-terminal FR1 codons onto a 5'-truncated contig,
    # also allow matching on ``partial_VDJ_{chain}_nt`` (the sequenced 3' suffix)
    # so the donor's real J→C junction and C-region codons are still recovered.
    search_nts = [vdj_nt.upper()]
    p_nt = row.get(f"partial_VDJ_{chain}_nt")
    p_aa = row.get(f"partial_VDJ_{chain}_aa")
    if row.get(f"{chain}_vdj_source") == VDJ_SOURCE_5PRIME_RESCUED:
        if isinstance(p_aa, str) and len(p_aa.strip()) >= 15:
            p_aa_clean = p_aa.strip()
            target_len = 3 * len(p_aa_clean)
            if len(vdj_nt) >= target_len:
                suffix_nt = vdj_nt[-target_len:].upper()
                if suffix_nt not in search_nts:
                    search_nts.append(suffix_nt)
            if isinstance(p_nt, str) and len(p_nt.strip()) >= target_len:
                p_nt_u = p_nt.strip().upper()
                for off in (0, 1, 2):
                    if off + target_len <= len(p_nt_u):
                        cand = p_nt_u[off : off + target_len]
                        tr_cand, _ = translate_dna(cand)
                        if tr_cand == p_aa_clean and cand not in search_nts:
                            search_nts.append(cand)
                            break
    candidates: Counter = Counter()
    for sample in samples:
        if sample not in sample_contigs:
            continue
        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue
            contig_seq = sample_contigs[sample][contig_id]
            contig_u = contig_seq.upper()
            for needle in search_nts:
                idx = contig_u.find(needle)
                if idx >= 0:
                    after = contig_seq[idx + len(needle):].upper()
                    if after:
                        candidates[after] += 1
                    break
    if candidates:
        return candidates.most_common(1)[0][0]
    return None


def _blend_constant_nt_with_contig(
    contig_nt_past_vdj: str,
    canonical_aa: str,
    canonical_nt_codon_opt: str,
    *,
    codon_species: str = "human",
) -> tuple[str, dict]:
    """Splice contig-derived NT for the J→C junction into a codon-optimized
    canonical C region NT.

    The result respects three principles, listed in order of precedence:

    1. **Junction-residue codon comes from the donor where possible.**
       The first nt(s) past the (trimmed) VDJ are biologically encoded by
       the J segment's terminal partial codon; the next nt(s) come from
       the C exon. Together they spell the J→C junction codon. Using the
       contig's actual bytes here keeps the assembled NT faithful to the
       donor's real sequence at this site.
    2. **Incomplete codons at the contig edge are completed using the
       canonical reference.** If the contig provides 1-2 nt of an
       otherwise-incomplete codon at its 3' end, we extend with the
       canonical codon's bytes (preferring the codon that matches both
       the contig prefix AND the canonical residue, so the protein output
       is preserved).
    3. **Everything past the contig coverage uses codon-optimized canonical
       NT.** Synthesis vendors prefer codon-optimized constructs, and the
       deep C region is invariant across donors — there's no donor-specific
       NT to preserve there.

    Stops at the first AA mismatch between the contig-translated codon and
    the canonical AA at that position — at that point the contig and
    canonical disagree about the protein, and we trust the canonical (the
    contig may be carrying a sequencing error or a non-canonical variant).
    The mismatch surfaces in the returned ``debug`` dict so callers can
    emit a QC warning.

    Parameters
    ----------
    contig_nt_past_vdj : str
        NT extracted from the contig starting one nt past the VDJ. Frame
        offset 0 — i.e., position 0 is the first nt of the J→C junction
        codon. May be empty (no contig coverage); in that case the entire
        ``canonical_nt_codon_opt`` is returned.
    canonical_aa : str
        Canonical AA sequence of the constant region (from the packaged
        FASTA, post-#100). Position 0 corresponds to the J→C junction
        residue (for β chains) or the first mature residue (for α).
    canonical_nt_codon_opt : str
        Codon-optimized back-translated NT for ``canonical_aa`` with
        stop codons already appended at the tail. Length
        ``3 * len(canonical_aa) + 3 * n_stops``. (Pre-#116 this was
        single-stop; from 2.4 onward it defaults to dual-stop.)

    Returns
    -------
    tuple[str, dict]
        ``(blended_nt, debug)`` where ``blended_nt`` is the assembled C
        region NT with stops preserved at the tail (caller does NOT
        re-append), and ``debug`` carries:

        - ``n_contig_codons``: number of full codons spliced from contig.
        - ``partial_codon_completed``: True when a partial codon at the
          contig 3' edge was extended with canonical bytes to spell the
          canonical residue (donor fidelity preserved in the partial nt).
        - ``partial_codon_dropped``: True when a partial codon WAS present
          but its bytes were incompatible with every canonical codon
          for the expected residue — partial bytes discarded and the
          full canonical codon used. Distinct from ``partial_codon_completed=False``
          when there simply was no partial codon to begin with.
        - ``aa_mismatch_at``: 1-indexed canonical AA position where the
          contig translation first disagreed with canonical (or ``None``
          if no mismatch was hit before the contig ran out).
        - ``source``: human-readable summary for the
          ``{chain}_constant_source`` column.

    Length invariant: ``len(blended_nt) == len(canonical_nt_codon_opt)``
    always, regardless of which path was taken. Equivalent to
    ``3 * len(canonical_aa) + 3 * n_stops`` since the input includes
    stops. Tested in ``TestBlendConstantNtWithContig``.
    """
    debug = {
        "n_contig_codons": 0,
        "partial_codon_completed": False,
        "partial_codon_dropped": False,
        "aa_mismatch_at": None,
        "source": "canonical-codon-opt",
    }
    if not contig_nt_past_vdj or not canonical_aa:
        return canonical_nt_codon_opt, debug

    n_full_codons_in_contig = len(contig_nt_past_vdj) // 3
    partial_nt = contig_nt_past_vdj[n_full_codons_in_contig * 3 :]

    # Find longest contig prefix whose per-codon translation matches the
    # canonical AA byte-for-byte. First mismatch ends the contig-derived
    # run; everything past that point uses canonical codon-optimized NT.
    n_matching = 0
    for i in range(min(n_full_codons_in_contig, len(canonical_aa))):
        codon = contig_nt_past_vdj[i * 3 : (i + 1) * 3]
        translated = CODON_TABLE.get(codon, "X")
        if translated == canonical_aa[i]:
            n_matching += 1
        else:
            debug["aa_mismatch_at"] = i + 1
            break

    contig_part = contig_nt_past_vdj[: n_matching * 3]
    debug["n_contig_codons"] = n_matching

    # Try to complete a partial codon at the contig 3' edge — only if
    # every full codon matched (otherwise we've already switched to
    # canonical for this position).
    completed_codon = ""
    canonical_offset = n_matching
    if (
        n_matching == n_full_codons_in_contig
        and partial_nt
        and canonical_offset < len(canonical_aa)
    ):
        expected_aa = canonical_aa[canonical_offset]
        # All codons that translate to expected_aa AND start with partial_nt.
        compatible = [
            cdn
            for cdn, aa in CODON_TABLE.items()
            if aa == expected_aa and cdn.startswith(partial_nt)
        ]
        if compatible:
            preferred_map = PREFERRED_CODONS_BY_SPECIES.get(
                codon_species, HUMAN_PREFERRED_CODONS
            )
            preferred = preferred_map.get(expected_aa)
            if preferred in compatible:
                completed_codon = preferred
            else:
                completed_codon = compatible[0]
            canonical_offset += 1
            debug["partial_codon_completed"] = True
        else:
            # Contig partial bytes incompatible with canonical residue —
            # discard the partial bytes and let canonical drive this codon.
            # Recorded separately from "no partial at all" so the audit
            # trail can flag clones where donor fidelity at the boundary
            # WAS available but had to be dropped.
            debug["partial_codon_dropped"] = True

    canonical_rest = canonical_nt_codon_opt[canonical_offset * 3 :]
    blended = contig_part + completed_codon + canonical_rest

    if n_matching > 0 and debug["aa_mismatch_at"] is None:
        suffix_parts = []
        if debug["partial_codon_completed"]:
            suffix_parts.append("partial completed")
        elif debug["partial_codon_dropped"]:
            suffix_parts.append("partial dropped — incompatible with canonical AA")
        suffix = (", " + "; ".join(suffix_parts)) if suffix_parts else ""
        debug["source"] = f"contig({n_matching} codons{suffix}) + canonical-codon-opt"
    elif n_matching > 0:
        debug["source"] = (
            f"contig({n_matching} codons) + canonical-codon-opt "
            f"(switched at AA mismatch pos {debug['aa_mismatch_at']})"
        )
    elif debug["partial_codon_completed"]:
        # No full codons, but a partial was completed — still donor-derived
        # at the J→C boundary even though nothing matched at the codon level.
        debug["source"] = "contig(partial completed) + canonical-codon-opt"
    elif debug["partial_codon_dropped"]:
        # No full codons AND partial dropped — purely canonical output, but
        # surface the dropped partial so the audit trail records that donor
        # fidelity was available and discarded.
        debug["source"] = (
            "canonical-codon-opt (partial dropped — incompatible with canonical AA)"
        )
    # else: no contig codons spliced and no partial — keep default
    #       "canonical-codon-opt"

    return blended, debug


def _extract_c_region_start_from_contig(
    row: pd.Series,
    sample_contigs: dict,
    vdj_aa: str,
    chain: str,
    n_aa: int = 15,
) -> str | None:
    """Pull the first ``n_aa`` residues observed after the VDJ in the
    CellRanger contig, if extractable.

    Mirrors the logic in :func:`_extract_leader_from_contigs_single` —
    longest-ORF translation, split on the VDJ AA, take what follows.
    Returns ``None`` when no contig contains the VDJ AA or the
    post-VDJ tail is shorter than ``n_aa``.
    """
    if not vdj_aa:
        return None
    samples = str(row.get("samples", "")).split(";")
    contig_col = f"{chain}_contig_ids"
    if contig_col not in row or pd.isna(row[contig_col]):
        return None

    contig_ids = str(row[contig_col]).split(";")
    search_aas = [vdj_aa]
    p_aa = row.get(f"partial_VDJ_{chain}_aa")
    if (
        row.get(f"{chain}_vdj_source") == VDJ_SOURCE_5PRIME_RESCUED
        and isinstance(p_aa, str)
        and len(p_aa.strip()) >= 15
        and p_aa.strip() not in search_aas
    ):
        search_aas.append(p_aa.strip())
    starts = Counter()
    for sample in samples:
        if sample not in sample_contigs:
            continue
        for contig_id in contig_ids:
            if contig_id not in sample_contigs[sample]:
                continue
            contig_seq = sample_contigs[sample][contig_id]
            translated, _, _ = find_longest_orf(contig_seq)
            for needle_aa in search_aas:
                if needle_aa in translated:
                    _, after = translated.split(needle_aa, 1)
                    if len(after) >= n_aa:
                        starts[after[:n_aa]] += 1
                    break
    if starts:
        return starts.most_common(1)[0][0]
    return None


def _add_constant_from_row(
    row: pd.Series,
    result: dict,
    *,
    species: str = "human",
    codon_species: str = "human",
    stops_nt: str = "TAATGA",
):
    """Add constant region sequences directly from row columns if present."""
    for chain in ["alpha", "beta"]:
        aa_col = f"{chain}_constant_aa"
        nt_col = f"{chain}_constant_nt"

        aa_val = row.get(aa_col) if aa_col in row else None
        nt_val = row.get(nt_col) if nt_col in row else None
        has_aa = pd.notna(aa_val) and str(aa_val).strip() != ""
        has_nt = pd.notna(nt_val) and str(nt_val).strip() != ""

        if has_aa:
            result[f"{chain}_constant_aa"] = str(aa_val).strip().rstrip("*")
        if has_nt:
            nt_str = str(nt_val).strip().upper()
            if len(nt_str) >= 3 and nt_str[-3:] not in {"TAA", "TAG", "TGA"} and stops_nt:
                nt_with_stop = nt_str + stops_nt
            else:
                nt_with_stop = nt_str
            result[f"{chain}_constant_nt"] = nt_with_stop
            result[f"{chain}_constant_nt_optimized"] = nt_with_stop
            if not has_aa:
                const_aa, _ = translate_dna(nt_str)
                const_aa = const_aa.rstrip("*")
                result[f"{chain}_constant_aa"] = const_aa
                has_aa = bool(const_aa)
        elif has_aa:
            opt_nt = (
                optimize_codons(result[f"{chain}_constant_aa"], species=codon_species)
                + stops_nt
            )
            result[f"{chain}_constant_nt"] = opt_nt
            result[f"{chain}_constant_nt_optimized"] = opt_nt

        if has_aa or has_nt:
            result[f"{chain}_constant_source"] = "from-data"
            result[f"{chain}_junction_residue_source"] = "from-data"

        canonical_name, _ = pick_canonical_constant(
            chain,
            c_gene=row.get(f"{chain}_c_gene"),
            j_gene=row.get(f"{chain}_j_gene"),
            species=species,
        )
        result[f"{chain}_c_gene_canonical"] = canonical_name


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
                    donor_opt = optimize_codons(donor_aa, species=recode.codon_species)
                    if sites or recode.max_codon_repeats != 1:
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


def _add_single_chain(df: pd.DataFrame, linker: str) -> pd.DataFrame:
    """Add single-chain construct (β-linker-α) including the NT triad.

    Emits four columns:

    * ``single_chain_aa``
    * ``single_chain_nt`` — built from the blend ``full_*_nt``. The
      column most callers reach for; donor-faithful at the J→C
      boundary, codon-optimized canonical otherwise.
    * ``single_chain_nt_optimized`` — built from
      ``full_*_nt_optimized``. The ready-to-order synthesis construct;
      same bytes for every donor with the same picked allele.

    No ``single_chain_nt_contig`` is emitted: the stitched contig construct is
    truncated at contig coverage (not a CDS, not synthesis-ready) and carries no
    QC signal beyond the per-chain ``{chain}_constant_nt_contig`` /
    ``full_{chain}_nt_contig`` provenance columns it was built from. Inspect
    those for the donor's raw bytes / contig backing.

    The β CDS in each variant has ALL trailing stops stripped before
    the linker — the 2A construct requires a continuous ORF across
    β → linker → α. The dual-stop default (``TAATGA``) made the
    one-codon-strip in 2.3 incorrect; the loop here is the fix.
    """
    if linker.upper() in LINKERS:
        linker_info = LINKERS[linker.upper()]
        linker_aa = linker_info["aa"]
        linker_nt = linker_info["dna"]
    else:
        # Custom linker sequence provided as amino acids; no NT.
        linker_aa = linker
        linker_nt = ""

    def strip_stop(seq):
        if not isinstance(seq, str):
            return None
        return seq.rstrip("*")

    beta_aa = df["full_beta_aa"].apply(strip_stop)
    alpha_aa = df["full_alpha_aa"].where(df["full_alpha_aa"].apply(lambda x: isinstance(x, str)))

    single_chain = beta_aa.fillna("") + linker_aa + alpha_aa.fillna("")
    missing_mask = df["full_beta_aa"].isna() | df["full_alpha_aa"].isna()
    df["single_chain_aa"] = single_chain.where(~missing_mask, pd.NA)

    def strip_stop_codon_dna(seq):
        if not isinstance(seq, str):
            return None
        # Strip ALL trailing stops. Single-stop strip pre-2.4 broke
        # frame across the linker once dual-stop became the default.
        while len(seq) >= 3 and seq[-3:] in {"TAA", "TAG", "TGA"}:
            seq = seq[:-3]
        return seq

    if linker_nt:
        # Build the triad. Each variant uses its matching beta/alpha
        # input column; if either input is missing or contains None
        # (contig view when contig didn't cover the C region) the
        # output is None — splicing a None into a 2A construct would
        # be silently wrong.
        # No ``_contig`` single-chain variant: the stitched contig construct is
        # truncated at contig coverage (not a CDS, not synthesis-ready) and adds
        # nothing over the per-chain ``{chain}_constant_nt_contig`` /
        # ``full_{chain}_nt_contig`` provenance columns it's built from.
        for suffix, out_col in (
            ("", "single_chain_nt"),
            ("_optimized", "single_chain_nt_optimized"),
        ):
            beta_col = f"full_beta_nt{suffix}"
            alpha_col = f"full_alpha_nt{suffix}"
            if beta_col not in df.columns or alpha_col not in df.columns:
                continue

            def _build_row(row, b=beta_col, a=alpha_col):
                bv = row.get(b)
                av = row.get(a)
                if not isinstance(bv, str) or not isinstance(av, str):
                    return None
                return strip_stop_codon_dna(bv) + linker_nt + av

            df[out_col] = df.apply(_build_row, axis=1)

    df["linker"] = linker_aa

    return df


# from_contig leader sources (the policy keeps/switches these; curated leaders
# the user chose explicitly are left alone except by an explicit force/secondary).
_FROM_CONTIG_SOURCES = frozenset(
    {"contig", "contig_kozak_reselected", "contig_hregion_trimmed", "contig_germline_anchored"}
)


def _resolve_named_leader(
    name: str, row, chain: str, *, species: str = "human"
) -> tuple[str, str | None, str] | None:
    """Resolve a leader directive → ``(aa, nt, source_tag)`` for ``chain``."""
    if name == "germline":
        g_aa = row.get(f"{chain}_germline_leader_aa")
        g_allele = row.get(f"{chain}_germline_allele")
        if isinstance(g_aa, str) and g_aa and isinstance(g_allele, str) and g_allele:
            gene, _, allele = g_allele.partition("*")
            return (
                g_aa,
                germline_leader_nt(gene, allele, species=species),
                "germline_reference_leader",
            )
        v_gene = row.get(f"{chain}_v_gene")
        cands = germline_vgene_leaders(v_gene, species=species)
        if cands:
            allele, _func, cand_aa = cands[0]
            gene_norm = normalize_vgene(v_gene) or str(v_gene)
            return (
                cand_aa,
                germline_leader_nt(gene_norm, allele, species=species),
                "germline_reference_leader",
            )
        return None
    hit = _lookup_default_leader(name)
    if hit is not None:
        canon_key, info = hit
        return info["aa"], info["dna"], f"curated:{canon_key}"
    return None


def _stash_preswitch(work: dict, chain: str) -> None:
    """Record a clone's donor leader variant BEFORE it's replaced (#276).

    When the policy switches (or a force overrides) a divergent contig leader,
    :func:`_substitute_chain_leader` refreshes ``{chain}_leader_vs_germline`` to
    describe the SHIPPED (germline/curated) leader — so the discarded donor
    variant would otherwise vanish from :func:`collect_germline_variants`. Stash
    it under ``{chain}_preswitch_*`` so the variant is still reported, tagged
    ``shipped="switched"``. No-op when the pre-switch leader already matched
    germline (nothing divergent to preserve).
    """
    diff = work.get(f"{chain}_leader_vs_germline")
    donor_aa = work.get(f"{chain}_leader_aa")
    germ_aa = work.get(f"{chain}_germline_leader_aa")
    # Prefer the recorded comparison; if it's absent (some paths don't stamp it)
    # but we have both sequences, derive it so the variant isn't lost on switch.
    if (not isinstance(diff, str) or not diff) and isinstance(donor_aa, str) and isinstance(germ_aa, str):
        diff = "identical" if donor_aa == germ_aa else characterize_divergence(donor_aa, germ_aa)
    if not isinstance(diff, str) or diff == "identical":
        return
    work[f"{chain}_preswitch_variant"] = diff
    work[f"{chain}_preswitch_donor_aa"] = work.get(f"{chain}_leader_aa")
    work[f"{chain}_preswitch_germline_aa"] = work.get(f"{chain}_germline_leader_aa")
    work[f"{chain}_preswitch_germline_allele"] = work.get(f"{chain}_germline_allele")
    work[f"{chain}_preswitch_germline_identity"] = work.get(f"{chain}_germline_identity")
    work[f"{chain}_preswitch_leader_source"] = work.get(f"{chain}_leader_source")


def _substitute_chain_leader(
    target: dict,
    orig,
    chain: str,
    new_aa: str,
    new_nt: str | None,
    source: str,
    *,
    species: str = "human",
    recode: CodonOptimization | None = None,
) -> None:
    """In ``target`` (a row dict), replace ``chain``'s leader with ``new_aa`` and
    re-stitch the affected full-sequence columns (including recoded ``_optimized``)."""
    eff_recode = recode or CodonOptimization(codon_species=species)
    donor_aa = orig.get(f"{chain}_leader_aa")
    donor_aa = donor_aa if isinstance(donor_aa, str) else ""
    donor_nt = orig.get(f"{chain}_leader_nt")
    donor_nt = donor_nt if isinstance(donor_nt, str) else ""
    eff_new_nt = (
        new_nt.upper()
        if isinstance(new_nt, str) and new_nt
        else optimize_codons(new_aa, species=eff_recode.codon_species)
    )
    target[f"{chain}_leader_aa"] = new_aa
    target[f"{chain}_leader_nt"] = new_nt
    target[f"{chain}_leader_len"] = len(new_aa)
    target[f"{chain}_leader_source"] = source
    target[f"{chain}_leader_qc"] = "substituted"
    feats = signal_peptide_features(new_aa)
    target[f"{chain}_sp_features_ok"] = feats["features_ok"]
    target[f"{chain}_sp_score"] = feats["score"]
    target[f"{chain}_sp_features"] = sp_features_summary(feats)
    cmp = germline_compare_leader(
        new_aa, orig.get(f"{chain}_v_gene"), species=species
    )
    if cmp is not None:
        allele, g_aa, identity, diff = cmp
        target[f"{chain}_germline_allele"] = f"{normalize_vgene(orig.get(f'{chain}_v_gene'))}*{allele}"
        target[f"{chain}_germline_leader_aa"] = g_aa
        target[f"{chain}_germline_identity"] = identity
        target[f"{chain}_leader_vs_germline"] = diff

    fa = orig.get(f"full_{chain}_aa")
    if isinstance(fa, str) and donor_aa and fa.startswith(donor_aa):
        target[f"full_{chain}_aa"] = new_aa + fa[len(donor_aa) :]
    elif isinstance(fa, str) and not donor_aa:
        target[f"full_{chain}_aa"] = new_aa + fa

    fn_blend = orig.get(f"full_{chain}_nt")
    if isinstance(fn_blend, str) and donor_nt and fn_blend.upper().startswith(donor_nt.upper()):
        target[f"full_{chain}_nt"] = (new_nt + fn_blend[len(donor_nt) :]) if new_nt else None
    elif isinstance(fn_blend, str) and not donor_nt:
        target[f"full_{chain}_nt"] = (new_nt + fn_blend) if new_nt else None

    fn_contig = orig.get(f"full_{chain}_nt_contig")
    if isinstance(fn_contig, str) and donor_nt and fn_contig.upper().startswith(donor_nt.upper()):
        target[f"full_{chain}_nt_contig"] = (new_nt + fn_contig[len(donor_nt) :]) if new_nt else None
    elif f"full_{chain}_nt_contig" in orig:
        target[f"full_{chain}_nt_contig"] = None

    # Rebuild full_{chain}_nt_optimized (and recoded position tags) from the new
    # leader + VDJ + constant_nt_optimized so translation always matches full_{chain}_aa.
    if f"full_{chain}_nt_optimized" in orig or orig.get(f"{chain}_constant_nt_optimized"):
        const_opt = target.get(f"{chain}_constant_nt_optimized") or orig.get(f"{chain}_constant_nt_optimized") or ""
        const_aa = target.get(f"{chain}_constant_aa") or orig.get(f"{chain}_constant_aa") or ""
        full_aa_now = target.get(f"full_{chain}_aa") or ""
        donor_aa_now = (
            full_aa_now[: len(full_aa_now) - len(const_aa)]
            if (const_aa and full_aa_now.endswith(const_aa))
            else full_aa_now
        )
        sites = eff_recode.restriction_sites()
        vdj_nt_raw = target.get(f"vdj_{chain}_nt") or orig.get(f"vdj_{chain}_nt") or orig.get(f"VDJ_{chain}_nt") or ""
        if isinstance(vdj_nt_raw, str):
            vdj_aa_raw = target.get(f"vdj_{chain}_aa") or orig.get(f"vdj_{chain}_aa") or orig.get(f"VDJ_{chain}_aa") or ""
            if isinstance(vdj_aa_raw, str) and vdj_aa_raw:
                vdj_nt_raw = vdj_nt_raw[: 3 * len(vdj_aa_raw)]
        else:
            vdj_nt_raw = ""
        if eff_recode.mode == "full" or not vdj_nt_raw:
            donor_opt = optimize_codons(donor_aa_now, species=eff_recode.codon_species)
            if sites or eff_recode.max_codon_repeats != 1:
                donor_opt, _ = recode_codons(
                    donor_opt,
                    max_codon_repeats=eff_recode.max_codon_repeats,
                    restriction_sites=sites,
                    species=eff_recode.codon_species,
                    right_context=const_opt,
                )
            target[f"full_{chain}_nt_optimized"] = (donor_opt + const_opt) if eff_new_nt else None
            target[f"{chain}_leader_recoded_positions"] = ""
            target[f"{chain}_vdj_recoded_positions"] = ""
        elif eff_new_nt:
            native_prefix = (eff_new_nt + vdj_nt_raw).upper()
            recoded_prefix, changed = recode_codons(
                native_prefix,
                max_codon_repeats=eff_recode.max_codon_repeats,
                restriction_sites=sites,
                species=eff_recode.codon_species,
                right_context=const_opt,
            )
            target[f"full_{chain}_nt_optimized"] = recoded_prefix + const_opt
            vdj_start = len(eff_new_nt)

            def _fmt(o, _np=native_prefix, _rp=recoded_prefix):
                return f"{o}:{_np[o : o + 3]}>{_rp[o : o + 3]}"

            target[f"{chain}_leader_recoded_positions"] = ";".join(
                _fmt(o) for o in changed if o < vdj_start
            )
            target[f"{chain}_vdj_recoded_positions"] = ";".join(
                _fmt(o) for o in changed if o >= vdj_start
            )
        else:
            target[f"full_{chain}_nt_optimized"] = None


def _unanimous_divergent_genes(df: pd.DataFrame) -> set[tuple[str, str]]:
    """``(chain, gene)`` whose from_contig leader is unanimous across the gene's
    clones AND diverges from germline — the consistent-divergent case.

    Divergence is ANY mismatch (substitution OR indel), not length only: a
    unanimous single-residue variant (e.g. a TRBV13*03-style donor SNP that the
    270-allele reference subset lacks) is the donor's likely real allele just as
    much as a unanimous indel is, so it must be protected from the keep/switch
    policy too. The SP-soundness gate in :func:`apply_leader_policy` still decides
    whether such a leader is actually kept."""
    out: set[tuple[str, str]] = set()
    for chain in ("alpha", "beta"):
        v_col, aa_col, g_col = (
            f"{chain}_v_gene", f"{chain}_leader_aa", f"{chain}_germline_leader_aa"
        )
        if not {v_col, aa_col, g_col} <= set(df.columns):
            continue
        sub = df[df[aa_col].notna() & df[v_col].notna()]
        for gene, grp in sub.groupby(v_col):
            leaders = set(grp[aa_col].astype(str))
            germs = set(grp[g_col].dropna().astype(str))
            if len(leaders) == 1 and len(germs) == 1:
                leader, germ = next(iter(leaders)), next(iter(germs))
                if leader != germ:  # any unanimous divergence (substitution or indel)
                    out.add((chain, str(gene)))
    return out


def _protected_divergent_leaders(
    df: pd.DataFrame, *, min_support: int = 2
) -> set[tuple[str, str, str]]:
    """``(chain, gene, leader_aa)`` for divergent leaders that look like a real
    donor allele rather than a one-off artifact — eligible to KEEP (the SP gate
    in :func:`apply_leader_policy` still decides the actual keep).

    A divergent leader qualifies when EITHER:

    * it is the gene's **unanimous** leader (whole-gene consistent — the original
      #270 case; covers single-clone genes), OR
    * it is carried by ``>= min_support`` clones of that gene (**het-aware**, #274):
      a heterozygous gene's minority allele has independent multi-clone support,
      so a TRBV6-2 G9A-style second allele is the donor's real allele even though
      the gene as a whole isn't unanimous. Only a singleton divergent leader in a
      non-unanimous gene (likely a sequencing artifact) is left unprotected.

    Divergence is judged per row (leader vs *that clone's* closest germline allele),
    so a het gene whose alleles resolve to different germline leaders is handled
    correctly. Returns ``leader_aa``-keyed tuples (not just gene) so only the
    supported allele is protected, not every leader the gene carries.
    """
    out: set[tuple[str, str, str]] = set()
    # Count support by distinct CLONE, not construct row: dual-α expansion
    # duplicates the β row (same β leader, two α variants), so a row-count would
    # promote a single dual-α clone's β leader to "multi-clone support" and keep
    # an artifact. `selected_clone` is shared across a clone's dual-α twins; fall
    # back to CDR3ab, then the row index, when it's absent.
    clone_col = next(
        (c for c in ("selected_clone", "CDR3ab") if c in df.columns), None
    )
    for chain in ("alpha", "beta"):
        v_col, aa_col, g_col = (
            f"{chain}_v_gene", f"{chain}_leader_aa", f"{chain}_germline_leader_aa"
        )
        if not {v_col, aa_col, g_col} <= set(df.columns):
            continue
        sub = df[df[aa_col].notna() & df[v_col].notna()]
        for gene, grp in sub.groupby(v_col):
            leaders = grp[aa_col].astype(str)
            keys = grp[clone_col].astype(str) if clone_col else grp.index.to_series()
            # distinct clones carrying each leader (de-duped over dual-α rows)
            support = keys.groupby(leaders).nunique()
            whole_gene_unanimous = leaders.nunique() == 1
            for _, r in grp.iterrows():
                leader = str(r[aa_col])
                germ = r[g_col]
                if not isinstance(germ, str) or leader == germ:
                    continue  # not germline-comparable, or matches germline (not divergent)
                if whole_gene_unanimous or support.get(leader, 0) >= min_support:
                    out.add((chain, str(gene), leader))
    return out


def _validate_leader_directive(val: str | None, param_name: str, *, allow_keep: bool = False) -> None:
    if val is None:
        return
    s = str(val).strip()
    if s == "germline" or (allow_keep and s == "keep"):
        return
    if _lookup_default_leader(s) is not None:
        return
    valid = ["germline", *(["keep"] if allow_keep else []), *DEFAULT_LEADERS.keys()]
    raise TCRsiftValidationError(
        f"Unknown {param_name}={val!r}; expected one of {valid}"
    )


def apply_leader_policy(
    df: pd.DataFrame,
    linker: str | None,
    *,
    leader_fallback: str = "germline",
    force_alpha: str | None = None,
    force_beta: str | None = None,
    secondary_alpha: str | None = None,
    secondary_beta: str | None = None,
    verbose: bool = False,
    species: str = "human",
    recode: CodonOptimization | None = None,
) -> pd.DataFrame:
    """Decide the shipped signal peptide for every clone (von Heijne policy, #270)."""
    _validate_leader_directive(leader_fallback, "leader_fallback", allow_keep=True)
    _validate_leader_directive(force_alpha, "force_alpha")
    _validate_leader_directive(force_beta, "force_beta")
    _validate_leader_directive(secondary_alpha, "secondary_alpha")
    _validate_leader_directive(secondary_beta, "secondary_beta")

    protected = _protected_divergent_leaders(df)
    force = {"alpha": force_alpha, "beta": force_beta}
    secondary = {"alpha": secondary_alpha, "beta": secondary_beta}
    df = df.copy()
    n_forced = n_switched = 0

    for idx, row in df.iterrows():
        work = row.to_dict()
        changed = False
        for chain in ("alpha", "beta"):
            aa = work.get(f"{chain}_leader_aa")
            has_leader = isinstance(aa, str) and bool(aa)
            if force[chain]:
                if not has_leader and not isinstance(work.get(f"full_{chain}_aa"), str):
                    continue
                res = _resolve_named_leader(
                    force[chain], work, chain, species=species
                )
                if res is not None:
                    new_aa, new_nt, src = res
                    if has_leader:
                        _stash_preswitch(work, chain)
                    _substitute_chain_leader(
                        work, work, chain, new_aa, new_nt, f"forced_{src}",
                        species=species, recode=recode,
                    )
                    changed = True
                    n_forced += 1
                continue
            if not has_leader:
                if (
                    work.get(f"{chain}_leader_source") == "missing_from_contig"
                    and isinstance(work.get(f"full_{chain}_aa"), str)
                    and work.get(f"full_{chain}_aa")
                    and leader_fallback != "keep"
                ):
                    res = _resolve_named_leader(
                        leader_fallback, work, chain, species=species
                    )
                    if res is None:
                        curated_fb = (
                            ("mCD28" if chain == "alpha" else "mCD8B")
                            if species == "mouse"
                            else ("CD28" if chain == "alpha" else "CD8A")
                        )
                        res = _resolve_named_leader(
                            curated_fb, work, chain, species=species
                        )
                    if res is not None:
                        new_aa, new_nt, src = res
                        _substitute_chain_leader(
                            work, work, chain, new_aa, new_nt, src,
                            species=species, recode=recode,
                        )
                        changed = True
                        n_switched += 1
                continue
            if work.get(f"{chain}_leader_source") not in _FROM_CONTIG_SOURCES:
                continue
            if leader_fallback == "keep":
                continue
            germ_aa = work.get(f"{chain}_germline_leader_aa")
            ident = work.get(f"{chain}_germline_identity")
            sp_ok = bool(work.get(f"{chain}_sp_features_ok"))
            matches_germline = (
                isinstance(germ_aa, str) and isinstance(ident, (int, float))
                and ident >= 0.999
            )
            divergent = isinstance(germ_aa, str) and not matches_germline
            gene = work.get(f"{chain}_v_gene")
            leader_now = str(work.get(f"{chain}_leader_aa"))
            is_consistent = (chain, str(gene), leader_now) in protected
            keep = matches_germline or (sp_ok and (not divergent or is_consistent))
            if keep:
                continue
            res = _resolve_named_leader(
                leader_fallback, work, chain, species=species
            )
            if res is None:
                continue
            new_aa, new_nt, src = res
            _stash_preswitch(work, chain)
            _substitute_chain_leader(
                work, work, chain, new_aa, new_nt, src, species=species, recode=recode
            )
            changed = True
            n_switched += 1
            if verbose:
                logger.warning(
                    "  %s leader rejected (sp_ok=%s, divergent=%s, consistent=%s) → "
                    "switched to %s", chain, sp_ok, divergent, is_consistent, src,
                )
        if changed:
            for k, v in work.items():
                if k not in df.columns:
                    df[k] = pd.Series([None] * len(df), index=df.index, dtype=object)
                df.at[idx, k] = v

    twins, orig_idx = [], []
    if any(secondary.values()):
        for idx, row in df.iterrows():
            sec = []
            for chain in ("alpha", "beta"):
                if not secondary[chain]:
                    continue
                cur_aa = row.get(f"{chain}_leader_aa")
                if not isinstance(cur_aa, str):
                    continue
                res = _resolve_named_leader(
                    secondary[chain], row, chain, species=species
                )
                if res is not None and res[0] != cur_aa:
                    sec.append((chain, secondary[chain], res))
            if not sec:
                continue
            twin = row.to_dict()
            for chain, name, (new_aa, new_nt, src) in sec:
                _substitute_chain_leader(
                    twin, row, chain, new_aa, new_nt, f"secondary_{src}",
                    species=species, recode=recode,
                )
            twin["leader_variant"] = "secondary:" + "+".join(s[1] for s in sec)
            twins.append(twin)
            orig_idx.append(idx)

    if twins:
        if "leader_variant" not in df.columns:
            df["leader_variant"] = pd.NA
        df.loc[orig_idx, "leader_variant"] = "primary"
        df = pd.concat([df, pd.DataFrame(twins)], ignore_index=True)

    # Rebuild single_chain (and refresh synthesis QC if already present) wherever leaders changed.
    if n_forced or n_switched or twins:
        if linker and {"full_beta_aa", "full_alpha_aa"} <= set(df.columns):
            df = _add_single_chain(df, linker)
        if "synth_gc_ok" in df.columns:
            df = add_synthesis_qc(df)
    if verbose and (n_forced or n_switched or twins):
        logger.info(
            "  Leader policy: %d forced, %d switched to fallback, %d secondary constructs",
            n_forced, n_switched, len(twins),
        )
    return df


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

