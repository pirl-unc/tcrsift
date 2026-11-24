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

"""Cell Ranger contig loading, constant-region extraction, and nucleotide blending."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import pandas as pd

from .codons import (
    CODON_TABLE,
    DEFAULT_FORBIDDEN_MOTIFS,
    HUMAN_PREFERRED_CODONS,
    PREFERRED_CODONS_BY_SPECIES,
    CodonOptimization,
    back_translate,
    find_longest_orf,
    optimize_codons,
    translate_dna,
)
from .constants import (
    _DEFAULT_MIN_PICKER_POSITIONS,
    ALLELE_REASON_AUTO_DETECTED,
    ALLELE_REASON_DIVERGENT_AT_POLYMORPHIC_POSITION,
    ALLELE_REASON_DIVERGENT_CONTIG,
    ALLELE_REASON_FRAME_ERROR,
    ALLELE_REASON_INVALID_OVERRIDE,
    ALLELE_REASON_NO_ALLELE_POOL,
    ALLELE_REASON_NO_CONTIG,
    ALLELE_REASON_OVERRIDDEN,
    ALLELE_REASON_SPARSE_CONTIG,
    CONSTANT_ALLELES_BY_SPECIES,
    CONSTANT_REGIONS_AA_BY_SPECIES,
    HUMAN_CONSTANT_ALLELES,
    HUMAN_CONSTANT_REGIONS_AA,
    _divergence_positions,
    _fallback_junction_residue,
    _format_divergence_positions,
    _pick_best_allele,
    _polymorphic_positions,
    _verify_first_codon_frame,
    pick_canonical_constant,
    verify_canonical_constant_start,
)
from .germline import VDJ_SOURCE_5PRIME_RESCUED

logger = logging.getLogger(__name__)


def parse_fasta(path: str | Path) -> dict[str, str]:
    """
    Parse a FASTA (or ``.fasta.gz``) file.

    Returns
    -------
    dict
        Mapping from sequence ID to sequence
    """
    import gzip

    path = Path(path)
    results = {}
    curr_id = None
    lines = []

    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_id and lines:
                    results[curr_id] = "".join(lines)
                    lines = []
                curr_id = line[1:].split()[0]  # Take first word after >
            else:
                lines.append(line)

        # Don't forget last entry
        if curr_id and lines:
            results[curr_id] = "".join(lines)

    return results


def parse_contig_json(path: str | Path) -> dict[str, str]:
    """Parse a CellRanger ``*contig_annotations*.json(.gz)`` file into ``{contig_id: sequence}``."""
    import gzip
    import json

    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        records = json.load(fh)
    results: dict[str, str] = {}
    if isinstance(records, list):
        for rec in records:
            if not isinstance(rec, dict):
                continue
            cid = rec.get("contig_name") or rec.get("contig_id")
            seq = rec.get("sequence")
            if cid and isinstance(seq, str) and seq:
                results[str(cid)] = seq.strip().upper()
    return results


SAMPLE_NAME_FROM_CHOICES = ("parent", "grandparent", "sheet")


def load_contigs(
    contig_dir: str | Path | None = None,
    *,
    sample_name_from: str = "parent",
    cellranger_dir: str | Path | None = None,
    sample_sheet: pd.DataFrame | None = None,
) -> dict[str, dict[str, str]]:
    """
    Load contig sequences from a CellRanger-style output tree (#124).

    Parameters
    ----------
    contig_dir : str or Path, optional
        Root directory to scan for ``*contig*.fasta``. Mutually exclusive
        with ``cellranger_dir``.
    sample_name_from : {'parent', 'grandparent', 'sheet'}, default 'parent'
        How to derive the sample name from each discovered FASTA's path:

        - ``'parent'`` (default, backward-compat) — sample = FASTA's
          immediate parent directory name. Matches the symlinked layout
          ``contig_dir/{sample}/filtered_contig.fasta``.
        - ``'grandparent'`` — sample = FASTA's grandparent directory name.
          Matches CellRanger ``multi`` output:
          ``per_sample_outs/{sample}/vdj_t/filtered_contig.fasta``.
        - ``'sheet'`` — match each FASTA against ``sample_sheet['vdj_dir']``
          and take the corresponding ``sample`` column. Most explicit;
          useful for non-standard layouts.
    cellranger_dir : str or Path, optional
        Shorthand for ``contig_dir=cellranger_dir, sample_name_from='grandparent'``.
        Pass a raw CellRanger ``per_sample_outs/`` directory and the sample
        names will resolve correctly without a symlink tree. Mutually
        exclusive with ``contig_dir``.
    sample_sheet : pandas.DataFrame, optional
        Required when ``sample_name_from='sheet'``. Must have columns
        ``sample`` and ``vdj_dir``; each FASTA's path must be inside one of
        the listed ``vdj_dir`` paths.

    Returns
    -------
    dict
        Nested dict: ``sample -> contig_id -> sequence``.
    """
    if cellranger_dir is not None:
        if contig_dir is not None:
            raise ValueError(
                "load_contigs: pass exactly one of contig_dir or cellranger_dir, not both."
            )
        contig_dir = cellranger_dir
        sample_name_from = "grandparent"
    if contig_dir is None:
        raise ValueError(
            "load_contigs: must pass either contig_dir or cellranger_dir."
        )
    if sample_name_from not in SAMPLE_NAME_FROM_CHOICES:
        raise ValueError(
            f"load_contigs: sample_name_from={sample_name_from!r} not in "
            f"{SAMPLE_NAME_FROM_CHOICES}"
        )
    if sample_name_from == "sheet" and sample_sheet is None:
        raise ValueError(
            "load_contigs: sample_name_from='sheet' requires sample_sheet=..."
        )

    contig_dir = Path(contig_dir)
    sample_contigs: dict[str, dict[str, str]] = {}

    sheet_dir_to_sample: dict[Path, str] = {}
    if sample_name_from == "sheet":
        missing = {"sample", "vdj_dir"} - set(sample_sheet.columns)
        if missing:
            raise ValueError(
                f"load_contigs: sample_sheet is missing required columns {sorted(missing)}"
            )
        sheet_dir_to_sample = {}
        for row in sample_sheet.itertuples(index=False):
            p = Path(row.vdj_dir).resolve()
            sheet_dir_to_sample[p.parent if p.is_file() else p] = row.sample

    def _sample_for(fasta_path: Path) -> str | None:
        if sample_name_from == "parent":
            return fasta_path.parent.name
        if sample_name_from == "grandparent":
            return fasta_path.parent.parent.name
        # sheet
        resolved = fasta_path.resolve()
        for sheet_path, sample in sheet_dir_to_sample.items():
            try:
                resolved.relative_to(sheet_path)
            except ValueError:
                continue
            return sample
        return None

    fasta_candidates = sorted(
        {
            *contig_dir.rglob("*contig*.fasta"),
            *contig_dir.rglob("*contig*.fasta.gz"),
        }
    )
    for fasta_path in fasta_candidates:
        if not fasta_path.is_file():
            continue
        sample_name = _sample_for(fasta_path)
        if sample_name is None:
            logger.warning(
                "load_contigs: %s did not match any sample_sheet vdj_dir; skipping",
                fasta_path,
            )
            continue
        sample_contigs.setdefault(sample_name, {}).update(parse_fasta(fasta_path))

    # Flat-layout fallback: bare ``*.fasta(.gz)`` files directly in contig_dir
    # that were not already matched above by ``*contig*.fasta(.gz)``.
    if sample_name_from == "parent":
        seen_fastas = {p.resolve() for p in fasta_candidates}
        for fasta_path in sorted({*contig_dir.glob("*.fasta"), *contig_dir.glob("*.fasta.gz")}):
            if not fasta_path.is_file() or fasta_path.resolve() in seen_fastas:
                continue
            stem = fasta_path.name[:-9] if fasta_path.name.lower().endswith(".fasta.gz") else fasta_path.stem
            sample_name = stem.split("_")[0]
            sample_contigs.setdefault(sample_name, {}).update(parse_fasta(fasta_path))

    # CellRanger JSON fallback: GEO deposits frequently ship
    # ``all_contig_annotations.json(.gz)`` instead of ``filtered_contig.fasta``.
    # Supplement any sample lacking FASTA coverage (or missing contig IDs) from JSON.
    json_candidates = sorted(
        {
            *contig_dir.rglob("*contig_annotations*.json"),
            *contig_dir.rglob("*contig_annotations*.json.gz"),
        }
    )
    for json_path in json_candidates:
        if not json_path.is_file():
            continue
        sample_name = _sample_for(json_path)
        if sample_name is None:
            continue
        try:
            parsed_json = parse_contig_json(json_path)
        except Exception as exc:
            logger.warning("load_contigs: could not parse %s: %s", json_path, exc)
            continue
        if parsed_json:
            target = sample_contigs.setdefault(sample_name, {})
            for cid, seq in parsed_json.items():
                target.setdefault(cid, seq)

    logger.info(f"Loaded contigs from {len(sample_contigs)} samples")
    return sample_contigs


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
    recode: CodonOptimization | None = None,
):
    """Add constant-region sequences and verify against the observed
    CellRanger contig start where possible (including chimeric constructs)."""
    allele_overrides = allele_overrides or {}
    is_chimeric = constant_species != vdj_species
    extra_sites = recode.restriction_sites() if recode is not None else ()
    const_forbidden = (
        tuple(dict.fromkeys((*DEFAULT_FORBIDDEN_MOTIFS, *extra_sites)))
        if extra_sites
        else DEFAULT_FORBIDDEN_MOTIFS
    )
    species_allele_pools = CONSTANT_ALLELES_BY_SPECIES.get(
        constant_species, HUMAN_CONSTANT_ALLELES
    )
    for chain in ["alpha", "beta"]:
        # Skip constant-region assembly when this chain has no VDJ sequence
        # (prevents spurious constant-allele warnings & no-call inflation on unpaired clones).
        if not result.get(f"vdj_{chain}_aa"):
            continue
        c_gene = result.get(f"{chain}_c_gene")
        j_gene = result.get(f"{chain}_j_gene")
        canonical_name, canonical_aa = pick_canonical_constant(
            chain, c_gene, j_gene, species=constant_species
        )

        if extra_sites:
            canonical_nt_codon_opt = (
                optimize_codons(
                    canonical_aa,
                    species=codon_species,
                    forbidden_motifs=const_forbidden,
                )
                + stops_nt
            )
        else:
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
                    contig_nt_past_vdj = None
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
                            optimize_codons(
                                canonical_aa,
                                species=codon_species,
                                forbidden_motifs=const_forbidden,
                            )
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

