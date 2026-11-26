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
Configuration management for TCRsift.

Provides a unified configuration system that works with both CLI and Python API.
Configuration can be loaded from YAML files, with CLI arguments taking precedence.
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class LoadConfig:
    """Configuration for the load step."""

    min_genes: int = 250
    max_genes: int = 15000
    min_counts: int = 500
    max_counts: int = 100000
    # Mitochondrial-content window. max_mito_pct removes high-mito dying cells
    # (standard scRNA-seq QC). min_mito_pct is a FLOOR — it discards cells
    # *below* the threshold, which targets near-zero-mito empty/ambient droplets.
    # The default 2.0 floor is retained; it is no longer silent — load_cellranger_gex
    # logs how many cells it drops so the cull is visible (#168). Set 0.0 to disable.
    min_mito_pct: float = 2.0
    max_mito_pct: float = 8.0


@dataclass
class PhenotypeConfig:
    """Configuration for the phenotype step."""

    cd4_cd8_ratio: float = 3.0
    # CD3 read gate, applied to cells entering clonotype aggregation in `run`.
    # Previously a dead parameter (set a column but never gated selection, #172);
    # now it actually drops low-CD3 (likely non-T / ambient) cells. Default 3 is
    # a gentle contamination floor validated on pilot data. Note a raw CD3 UMI
    # count scales with sequencing depth, so this gate is depth-dependent across
    # datasets — set 0 to disable, or raise it for deeper runs.
    min_cd3_reads: int = 3

    # Generalized per-cell annotator knobs (#312; used by
    # tcrsift.annotate_cells when the embedding is enabled). None cell-type
    # restriction = consider every reference type (whole-tissue TIL); pass
    # ["T cell", "NK cell", ...] (or the PBMC_CULTURE_TYPES names) for cultures.
    allowed_types: list[str] | None = None
    # Cluster-mean (log1p CP10K) floors for gated calls: a shared master-TF /
    # activation signature isn't enough — the defining gene must clear these.
    treg_foxp3_min: float = 0.2
    th1_ifng_min: float = 0.2
    th2_il13_min: float = 0.2
    th17_il17a_min: float = 0.1
    gd_tcr_min: float = 0.5      # δ-constant TRDC floor for a γδ T call
    mait_min: float = 0.3        # SLC4A10 floor for a MAIT call
    activation_min: float = 0.5  # activation-module floor for the "activated" catch-all
    # Sub-clusters below this fraction of their cell type lose the
    # distinguishing-gene suffix (a <1% fragment shouldn't read as a subtype).
    min_subtype_frac: float = 0.01


@dataclass
class ClonotypeConfig:
    """Configuration for the clonotype step."""

    group_by: str = "CDR3ab"
    handle_doublets: str = "flag"
    min_umi: int = 2
    # Warn (loudly, regardless of verbose) when the multi-chain/doublet rate
    # meets/exceeds this fraction (#165). 0 disables the warning.
    doublet_warn_rate: float = 0.1
    # Annotate every clonotype with data-driven Pgen/Ppost (per α/β chain) from
    # the shipped k-mer background, so ppost_alpha/ppost_beta are available for
    # PRISM selection without a separate step. Default on; no extra deps.
    add_pgen_ppost: bool = True


@dataclass
class FilterConfig:
    """Configuration for the filter step."""

    method: str = "threshold"
    tcell_type: str = "cd8"
    min_cells: int = 2
    min_frequency: float = 0.0
    require_complete: bool = True
    fdr_tiers: list[float] = field(default_factory=lambda: [0.15, 0.1, 0.01, 0.001, 0.0001])
    # Donor / method knobs (#15). Each is no-op when its source column isn't
    # on the clonotype table — preserves backwards compat for designs that
    # don't supply patient_id / enrichment_method.
    filter_mode: str = "fdr"
    min_donors: int = 0
    min_methods_per_donor: int = 0
    min_cells_per_method: int = 0
    min_frequency_per_method: float = 0.0
    # Timepoint / APC axis knobs (#9 chunk 3). No-op when their underlying
    # clonotype-table columns aren't present.
    min_timepoints: int = 0
    min_timepoints_per_donor: int = 0
    min_apcs: int = 0
    min_apcs_per_donor: int = 0
    min_til_cells_per_donor: int = 0
    # FDR scope (#26): {auto, global, per-donor, per-sample}.
    # 'auto' resolves to 'per-donor' when n_donors>1 and the sheet does
    # NOT set donors_share_antigen, else 'global'.
    fdr_scope: str = "auto"
    # Cohort-level signal that pooling across donors is biologically valid;
    # may be set explicitly via --donors-share-antigen or read from the YAML
    # sheet root.
    donors_share_antigen: bool = False
    # Paired stimulation-vs-control filtering thresholds (method="stimulation"):
    stim_fdr_threshold: float = 0.1
    stim_log2_fc_threshold: float = 1.0
    stim_delta_threshold: float = 0.5


@dataclass
class AttributionConfig:
    """Configuration for partial-information clonotype attribution (#176).

    Opt-in (``enabled=False`` by default → today's integer behavior is
    byte-identical). When on, each cell distributes a total weight of 1.0
    across complete paired alpha-beta clones: alpha-dropout cells are shared
    across clones with the same beta, droplet doublets are split, and recurrent
    dual-alpha cells (allelic inclusion) are merged into one clone. Weighted
    counts then replace integer ``cell_count``/``frequency`` everywhere, which
    also restores the per-sample frequency-sum invariant (subsumes #175).
    """

    enabled: bool = False
    # Share alpha-dropout cells across complete clones with the same beta.
    beta_sharing: bool = True
    # Split droplet doublets (dual-alpha/dual-beta singletons) across the
    # candidate complete clones rather than dropping or keeping primary only.
    doublet_split: bool = True
    # Merge recurrent dual-alpha cells (allelic inclusion) into one clone.
    dual_alpha_merge: bool = True
    # Minimum cells sharing a sorted (a1,a2,b) triple to treat it as one
    # biological clone rather than a per-cell droplet doublet.
    dual_alpha_min_cells: int = 2
    # Dominance gate (#198): only merge a dual-alpha triple when it is the beta's
    # *dominant* configuration — the merged pair must co-occur in at least this
    # fraction of the beta's paired cells. True allelic inclusion ~= every cell
    # carries both alphas (fraction ~1.0); a public/promiscuous beta's chance
    # recurrence is a tiny fraction, so 0.5 (majority) rejects those coincidences.
    dual_alpha_min_fraction: float = 0.5
    # Optional: also refuse to merge when the beta pairs with more than this many
    # distinct alphas (a promiscuous/public beta). 0 = disabled (rely on the
    # fraction gate). cf. the #148 pairing-promiscuity feature.
    dual_alpha_max_beta_partners: int = 0
    # How shared/split weight is apportioned: "proportional" (to each target
    # clone's complete-cell abundance) or "uniform".
    share_weighting: str = "proportional"
    # Reserved EM hook; v1 uses a single distribution pass (no iteration).
    em_iterations: int = 1
    # Orthogonal split controls (#181): split non-merged dual-alpha and dual-beta
    # independently. None -> inherit the deprecated `doublet_split` alias, so old
    # configs are unchanged. dual-alpha and dual-beta are biologically different
    # (TRA allelic inclusion is common/real; dual-beta is rarer, more doublet-like).
    dual_alpha_split: bool | None = None
    dual_beta_split: bool | None = None

    def __post_init__(self):
        if self.dual_alpha_split is None:
            self.dual_alpha_split = self.doublet_split
        if self.dual_beta_split is None:
            self.dual_beta_split = self.doublet_split


@dataclass
class AnnotateConfig:
    """Configuration for the annotate step."""

    vdjdb_path: str | None = None
    iedb_path: str | None = None
    cedar_path: str | None = None
    iedb_epitope_path: str | None = None
    match_by: str = "CDR3ab"
    match_strictness: str | None = None
    match_mode: str = "exact"
    max_distance: int | None = None
    exclude_viral: bool = False
    flag_only: bool = False
    host_species_filter: str = "all"
    auto_download: bool = False
    call_multimers: bool = False
    multimer_min_umi: float = 10.0
    multimer_min_nc_ratio: float = 5.0
    multimer_min_concordance: float = 0.5
    multimer_min_specific_cells: int = 2

    def __post_init__(self):
        if self.host_species_filter not in ("all", "same"):
            raise ValueError(
                f"Invalid host_species_filter: {self.host_species_filter!r}; "
                "expected 'all' or 'same'"
            )
        if self.match_strictness is not None and self.match_strictness not in (
            "strict_ab",
            "ab_with_partial",
            "b_only",
        ):
            raise ValueError(
                f"Invalid match_strictness: {self.match_strictness!r}; "
                "expected 'strict_ab', 'ab_with_partial', or 'b_only'"
            )
        if self.match_mode not in ("exact", "levenshtein", "tcrdist"):
            raise ValueError(
                f"Invalid match_mode: {self.match_mode!r}; "
                "expected 'exact', 'levenshtein', or 'tcrdist'"
            )
        if self.max_distance is not None:
            min_dist = 0 if self.match_mode == "tcrdist" else 1
            if (
                not isinstance(self.max_distance, int)
                or isinstance(self.max_distance, bool)
                or self.max_distance < min_dist
            ):
                raise ValueError(
                    f"Invalid max_distance: {self.max_distance!r}; "
                    f"expected an integer >= {min_dist}"
                )


@dataclass
class TILConfig:
    """Configuration for the TIL matching step."""

    match_by: str = "CDR3ab"
    min_til_cells: int = 1
    til_samples: list[str] = field(default_factory=list)


@dataclass
class SCTConfig:
    """Configuration for SCT (single-cell TCR) data loading."""

    min_snr: float = 2.0
    min_reads_per_chain: int = 10
    require_mutation_match: bool = True
    require_compact_match: bool = False


@dataclass
class GEXConfig:
    """Configuration for gene expression augmentation."""

    gene_list: list[str] | None = None  # None = use default T cell markers
    gene_groups: dict[str, list[str]] | None = None  # None = use default groups
    include_qc: bool = True
    aggregation_ops: list[str] = field(default_factory=lambda: ["sum", "mean"])


@dataclass
class EmbedConfig:
    """Configuration for the single-cell embedding step (#311).

    An OPTIONAL atlas step (``enabled=False`` by default) that runs before
    clonotype aggregation: analytic Pearson residuals on an informative-gene
    panel → PCA → optional Harmony batch integration → neighbors → UMAP →
    Leiden. Signature scoring/labeling stay on log1p CP10K (see
    :mod:`tcrsift.gex`); this step only computes structure.
    """

    enabled: bool = False
    # Clustering gene panel. None → pick highly-variable genes automatically
    # (``n_top_genes``). Provide a curated list so structure reflects known
    # biology rather than whatever is most variable (technical / ambient).
    informative_genes: list[str] | None = None
    n_top_genes: int = 2000
    n_pcs: int = 40
    # obs column to integrate over with Harmony (e.g. "sample"/"batch"). None
    # skips integration and embeds on raw PCA.
    batch_key: str | None = None
    leiden_resolution: float = 1.0
    n_neighbors: int = 15
    min_cells: int = 3
    # Drop TCR/Ig receptor transcripts from the clustering panel so structure
    # isn't driven by clonotype (anti-leakage, consistent with the rest of the
    # library).
    exclude_receptor_genes: bool = True
    seed: int = 0


@dataclass
class UnifyConfig:
    """Configuration for multi-experiment unification."""

    add_occurrence_flags: bool = True
    add_combined_stats: bool = True
    add_phenotype_confidence: bool = True
    phenotype_ratio_threshold: float = 10.0


@dataclass
class AssembleConfig:
    """Configuration for the assemble step."""

    # Whether `run` assembles full-length sequences for the whole filtered
    # repertoire. Set False (CLI: --no-assemble) for a select-then-assemble
    # flow — analysis over all clones, but synthesis sequences built only for
    # the selected basket via `select` + `report selected` (#246). Skips the
    # ~15–40× wasted assembly and keeps the fail-closed gate (#243) on the
    # shipped basket rather than aborting on non-selected clones.
    enabled: bool = True
    alpha_leader: str | None = "CD28"  # None, "from_contig", or key: CD8A, CD28, IgK, TRAC, TRBC
    beta_leader: str | None = "CD8A"  # None, "from_contig", or key: CD8A, CD28, IgK, TRAC, TRBC
    # Leader policy (#270). leader_fallback = what a REJECTED from_contig leader
    # (bad SP, or divergent-from-germline + inconsistent across the gene's clones)
    # switches to: "germline" (default — the gene's germline V-leader; if the gene
    # isn't in the reference, keep + flag) or a curated key (CD8A/CD28/IgK/TRAC/
    # TRBC). An SP-sound, consistent contig leader is kept.
    leader_fallback: str = "germline"
    # Force a specific leader (germline / CD8A / …) for every construct on a
    # chain, ignoring the contig. None = normal policy.
    force_alpha_leader: str | None = None
    force_beta_leader: str | None = None
    # Also emit a SECONDARY construct using this leader (germline / CD8A / …) per
    # chain, even for kept good calls. None = no secondary construct.
    secondary_alpha_leader: str | None = None
    secondary_beta_leader: str | None = None
    # Codon optimization of full_*_nt_optimized (#293). "focal" (default) keeps
    # the donor's native leader+VDJ codons and makes only the minimal synonymous
    # swaps that clear a constraint; "full" recodes the whole ORF. Constraints:
    # max_codon_repeats (default 1 = no immediate identical-codon repeat) and
    # avoid_enzymes (restriction enzyme names / ACGT sites; default none).
    codon_optimization: str = "focal"
    max_codon_repeats: int = 1
    avoid_enzymes: list[str] = field(default_factory=list)
    include_constant: bool = True
    constant_source: str = "ensembl"
    linker: str = "T2A"
    contigs_dir: str | None = None
    single_chain: bool = True
    # Contig sample-name policy (#124). "parent" (default) = immediate
    # dir name; "grandparent" = CellRanger per_sample_outs/{sample}/vdj_t
    # layout; "sheet" = match against the sample sheet's vdj_dir.
    sample_name_from: str = "parent"
    # Shorthand: a raw CellRanger per_sample_outs dir; implies
    # sample_name_from="grandparent".
    cellranger_dir: str | None = None
    # Optional species overrides for chimeric constructs (e.g. human VDJ +
    # mouse constant region via constant_species="mouse").
    constant_species: str | None = None
    codon_species: str | None = None
    trac_allele: str = "auto"
    trbc1_allele: str = "auto"
    trbc2_allele: str = "auto"
    stop_codons: list[str] = field(default_factory=lambda: ["TAA", "TGA"])
    germline_vdj_fallback: bool = True
    auto_download_germline: bool = False

    def __post_init__(self):
        from .species import normalize_species

        if isinstance(self.alpha_leader, str) and self.alpha_leader.strip().lower() == "none":
            self.alpha_leader = None
        if isinstance(self.beta_leader, str) and self.beta_leader.strip().lower() == "none":
            self.beta_leader = None
        if self.constant_species is not None:
            self.constant_species = normalize_species(self.constant_species)
        if self.codon_species is not None:
            self.codon_species = normalize_species(self.codon_species)
        if isinstance(self.stop_codons, tuple):
            self.stop_codons = list(self.stop_codons)


@dataclass
class OutputConfig:
    """Configuration for output options."""

    generate_plots: bool = True
    generate_report: bool = True
    report_format: str = "pdf"
    output_airr: bool = False
    output_fasta: bool = False
    # Long-format (clone, sample) CSV (#20 chunk 1). 'auto' = emit when the
    # run has >=2 samples; 'always' = always; 'never' = skip.
    emit_clone_sample_long: str = "auto"
    # Per-method ranked CSVs (#20 chunk 2). Top-N clones per (donor, method)
    # ranked by within-bucket frequency. Skipped when enrichment_method
    # axis isn't populated.
    per_method_top_n: int = 100
    # Method × method overlap matrix similarity metric (#27 chunk 3).
    # 'jaccard' / 'dice' / 'count'.
    method_overlap_similarity: str = "jaccard"
    # Figure output format (#169, #258): 'png' (default), 'pdf', 'svg', or
    # 'both'. A vector choice emits the vector file(s) alongside a PNG (the PDF
    # report embeds raster); 'both' emits PDF + SVG + PNG.
    plot_format: str = "png"
    # Baseline/exclusion markers pushed last in combined figure labels (#208):
    # e.g. ['CTY'] makes CTY⁻ read last in AIM⁺CTY⁻ and in overlap combinations.
    # None = use the format module default (("CTY",)). Set to [] to disable.
    baseline_markers: list[str] | None = None


_SECTION_CLASSES: dict[str, type] = {
    "load": LoadConfig,
    "phenotype": PhenotypeConfig,
    "clonotype": ClonotypeConfig,
    "attribution": AttributionConfig,
    "filter": FilterConfig,
    "annotate": AnnotateConfig,
    "til": TILConfig,
    "sct": SCTConfig,
    "gex": GEXConfig,
    "embed": EmbedConfig,
    "unify": UnifyConfig,
    "assemble": AssembleConfig,
    "output": OutputConfig,
}

_INVERTED_FLAT_ALIASES: frozenset[str] = frozenset(
    {"skip_plots", "no_pgen_ppost", "no_assemble"}
)


def _build_flat_to_nested() -> dict[str, tuple[str, str]]:
    """Auto-derive unambiguous flat -> (section, field) mappings from dataclasses,
    then overlay historical CLI/YAML aliases and shared-field disambiguations."""
    from collections import Counter
    from dataclasses import fields as dc_fields

    counts = Counter(
        f.name for sec_cls in _SECTION_CLASSES.values() for f in dc_fields(sec_cls)
    )
    mapping: dict[str, tuple[str, str]] = {
        f.name: (sec_name, f.name)
        for sec_name, sec_cls in _SECTION_CLASSES.items()
        for f in dc_fields(sec_cls)
        if counts[f.name] == 1
    }
    # Explicit aliases & disambiguations for shared field names (enabled, min_cells, match_by)
    mapping.update(
        {
            "min_mito": ("load", "min_mito_pct"),
            "max_mito": ("load", "max_mito_pct"),
            "no_pgen_ppost": ("clonotype", "add_pgen_ppost"),
            "attribution": ("attribution", "enabled"),
            "enabled": ("attribution", "enabled"),
            "min_cells": ("filter", "min_cells"),
            "vdjdb": ("annotate", "vdjdb_path"),
            "iedb": ("annotate", "iedb_path"),
            "cedar": ("annotate", "cedar_path"),
            "iedb_epitope": ("annotate", "iedb_epitope_path"),
            "match_by": ("annotate", "match_by"),
            "database_match": ("annotate", "match_strictness"),
            "til_match_by": ("til", "match_by"),
            "embed": ("embed", "enabled"),
            "no_assemble": ("assemble", "enabled"),
            "skip_plots": ("output", "generate_plots"),
        }
    )
    return mapping


_FLAT_TO_NESTED: dict[str, tuple[str, str]] = _build_flat_to_nested()


@dataclass
class TCRsiftConfig:
    """
    Unified configuration for TCRsift pipeline.

    All parameters have sensible defaults. Configuration can be:
    1. Created programmatically with keyword arguments
    2. Loaded from a YAML file
    3. Merged with CLI arguments (CLI takes precedence)

    Examples
    --------
    >>> # Default config
    >>> config = TCRsiftConfig()

    >>> # From YAML
    >>> config = TCRsiftConfig.from_yaml("config.yaml")

    >>> # Programmatic with overrides
    >>> config = TCRsiftConfig(
    ...     load=LoadConfig(min_genes=300),
    ...     filter=FilterConfig(tcell_type="cd4"),
    ... )
    """

    load: LoadConfig = field(default_factory=LoadConfig)
    phenotype: PhenotypeConfig = field(default_factory=PhenotypeConfig)
    clonotype: ClonotypeConfig = field(default_factory=ClonotypeConfig)
    attribution: AttributionConfig = field(default_factory=AttributionConfig)
    filter: FilterConfig = field(default_factory=FilterConfig)
    annotate: AnnotateConfig = field(default_factory=AnnotateConfig)
    til: TILConfig = field(default_factory=TILConfig)
    sct: SCTConfig = field(default_factory=SCTConfig)
    gex: GEXConfig = field(default_factory=GEXConfig)
    embed: EmbedConfig = field(default_factory=EmbedConfig)
    unify: UnifyConfig = field(default_factory=UnifyConfig)
    assemble: AssembleConfig = field(default_factory=AssembleConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    # Free-form selection-rule language (#122/#125). Nested by design
    # (rules + global_rank with arbitrary method/pair names), so it's a
    # raw dict rather than a typed section. Empty = no selection step.
    selection: dict[str, Any] = field(default_factory=dict)

    # Free-form in-silico filter layer (#149): an ordered list of per-clone
    # percentile predicates (on Ppost / GEX scores) that refine each assay
    # group into a named composite twin. Raw dict by design. Empty = off.
    insilico_filter: dict[str, Any] = field(default_factory=dict)

    # Global options
    species: str = "auto"
    verbose: bool = False

    @classmethod
    def from_yaml(cls, path: str | Path) -> TCRsiftConfig:
        """
        Load configuration from a YAML file.

        Parameters
        ----------
        path : str or Path
            Path to the YAML configuration file

        Returns
        -------
        TCRsiftConfig
            Configuration loaded from the file
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        import yaml

        with open(path) as f:
            data = yaml.safe_load(f) or {}

        return cls._from_dict(data)

    @classmethod
    def _from_dict(
        cls,
        data: dict[str, Any],
        *,
        _overrides: dict[str, Any] | None = None,
    ) -> TCRsiftConfig:
        """Create config from a dictionary (with optional flat CLI overrides)."""
        flat_to_nested = _FLAT_TO_NESTED

        # Initialize nested dictionaries
        nested: dict[str, dict[str, Any]] = {sec: {} for sec in _SECTION_CLASSES}
        global_opts: dict[str, Any] = {}

        def _apply_mapping(source: dict[str, Any], *, warn_unknown: bool = False) -> None:
            for key, value in source.items():
                # A key that names a section AND is given as a dict is a nested
                # block, even when a same-named flat scalar alias exists (e.g.
                # 'attribution', which is both a section and a --attribution bool).
                if key in nested and (isinstance(value, dict) or value is None):
                    if isinstance(value, dict):
                        nested[key].update(value)
                elif key in flat_to_nested:
                    section, param = flat_to_nested[key]
                    if key in _INVERTED_FLAT_ALIASES:
                        value = not value
                    nested[section][param] = value
                elif key == "verbose":
                    global_opts["verbose"] = value
                elif key == "species":
                    from .species import normalize_species

                    global_opts["species"] = normalize_species(value, allow_auto=True)
                elif key == "selection":
                    global_opts["selection"] = value if isinstance(value, dict) else {}
                elif key == "insilico_filter":
                    global_opts["insilico_filter"] = (
                        value if isinstance(value, dict) else {}
                    )
                elif warn_unknown:
                    logger.warning(
                        "Unknown configuration key %r ignored (check spelling or section nesting)",
                        key,
                    )

        _apply_mapping(data, warn_unknown=True)
        if _overrides:
            _apply_mapping(_overrides, warn_unknown=False)

        section_objs = {}
        for sec, sec_cls in _SECTION_CLASSES.items():
            try:
                section_objs[sec] = sec_cls(**nested[sec])
            except TypeError as exc:
                from .validation import TCRsiftValidationError

                valid = sorted(f.name for f in dataclasses.fields(sec_cls))
                unknown = sorted(set(nested[sec]) - set(valid))
                if unknown:
                    raise TCRsiftValidationError(
                        f"Unknown configuration key(s) in section '{sec}': {', '.join(unknown)}",
                        hint=f"Valid keys for '{sec}': {', '.join(valid)}",
                    ) from exc
                raise TCRsiftValidationError(
                    f"Invalid configuration in section '{sec}': {exc}",
                ) from exc
        return cls(**section_objs, **global_opts)

    def to_yaml(self, path: str | Path) -> None:
        """
        Save configuration to a YAML file.

        Parameters
        ----------
        path : str or Path
            Path to save the configuration
        """
        import yaml

        data = self.to_dict()
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration to a dictionary."""
        return {
            "load": dataclasses.asdict(self.load),
            "phenotype": dataclasses.asdict(self.phenotype),
            "clonotype": dataclasses.asdict(self.clonotype),
            "attribution": dataclasses.asdict(self.attribution),
            "filter": dataclasses.asdict(self.filter),
            "annotate": dataclasses.asdict(self.annotate),
            "til": dataclasses.asdict(self.til),
            "sct": dataclasses.asdict(self.sct),
            "gex": dataclasses.asdict(self.gex),
            "embed": dataclasses.asdict(self.embed),
            "unify": dataclasses.asdict(self.unify),
            "assemble": dataclasses.asdict(self.assemble),
            "output": dataclasses.asdict(self.output),
            "selection": self.selection,
            "insilico_filter": self.insilico_filter,
            "species": self.species,
            "verbose": self.verbose,
        }

    def merge_with_args(self, args: argparse.Namespace) -> TCRsiftConfig:
        """
        Merge configuration with CLI arguments.

        CLI arguments take precedence over config file values.
        Only non-None CLI arguments override config values.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed CLI arguments

        Returns
        -------
        TCRsiftConfig
            New config with CLI overrides applied
        """
        # Convert args to dict, filtering out None values
        args_dict = {k: v for k, v in vars(args).items() if v is not None}

        # Start with a deep copy of the nested config dictionary so shared
        # field names across sections (e.g. `enabled`, `min_cells`, `match_by`)
        # never collide and unmapped nested fields are never lost.
        merged_data = copy.deepcopy(self.to_dict())

        # CLI flags using `action="store_true"` without `default=None` evaluate
        # to `False` when omitted on the command line; ignore `False` for them
        # so they never clobber `True` values configured in YAML.
        store_true_default_false = {
            "verbose",
            "no_pgen_ppost",
            "no_assemble",
        }

        cli_overrides: dict[str, Any] = {}
        if args_dict.get("no_leaders"):
            cli_overrides["alpha_leader"] = None
            cli_overrides["beta_leader"] = None
        elif args_dict.get("leaders_from_contigs"):
            cli_overrides["alpha_leader"] = "from_contig"
            cli_overrides["beta_leader"] = "from_contig"

        for key, value in args_dict.items():
            if key in (
                "func",
                "command",
                "config",
                "sample_sheet",
                "input",
                "output",
                "output_dir",
                "no_leaders",
                "leaders_from_contigs",
                "til_sample",
                "no_sheet_contigs",
                "allow_canonical_fallback",
                "cache_dir",
                "report_name",
            ):
                continue
            if key in store_true_default_false and value is False:
                continue
            if key == "fdr_tiers" and isinstance(value, str):
                value = [float(x) for x in value.split(",")]
            if key == "til_samples" and isinstance(value, str):
                value = [x.strip() for x in value.split(",")]
            if key in ("alpha_leader", "beta_leader") and value == "none":
                value = None
            cli_overrides[key] = value

        # CLI override of one contig source clears the mutually exclusive
        # counterpart inherited from YAML config so assemble_full_sequences
        # doesn't fail with "pass exactly one of contigs_dir or cellranger_dir".
        if "cellranger_dir" in cli_overrides and "contigs_dir" not in cli_overrides:
            cli_overrides["contigs_dir"] = None
        elif "contigs_dir" in cli_overrides and "cellranger_dir" not in cli_overrides:
            cli_overrides["cellranger_dir"] = None

        return TCRsiftConfig._from_dict(merged_data, _overrides=cli_overrides)


def add_config_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the --config argument to a parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to add the argument to
    """
    parser.add_argument(
        "--config",
        "-c",
        help="YAML configuration file (CLI args override config values)",
        metavar="FILE",
    )


def load_config_with_args(args: argparse.Namespace) -> TCRsiftConfig:
    """
    Load configuration, applying CLI overrides.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments (may include --config)

    Returns
    -------
    TCRsiftConfig
        Configuration with CLI overrides applied
    """
    # Start with defaults or load from file
    if hasattr(args, "config") and args.config:
        config = TCRsiftConfig.from_yaml(args.config)
    else:
        config = TCRsiftConfig()

    # Apply CLI overrides
    return config.merge_with_args(args)


# Convenience function to generate example config
def generate_example_config(path: str | Path = "tcrsift_config.yaml") -> None:
    """
    Generate an example configuration file with all defaults.

    Parameters
    ----------
    path : str or Path
        Path to save the example config
    """
    config = TCRsiftConfig()
    config.to_yaml(path)
