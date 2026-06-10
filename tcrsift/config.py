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
import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


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
    match_by: str = "CDR3ab"
    exclude_viral: bool = False
    flag_only: bool = False


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

        with open(path) as f:
            data = yaml.safe_load(f) or {}

        return cls._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> TCRsiftConfig:
        """Create config from a dictionary."""
        # Map flat keys to nested structure for convenience
        flat_to_nested = {
            # Load
            "min_genes": ("load", "min_genes"),
            "max_genes": ("load", "max_genes"),
            "min_counts": ("load", "min_counts"),
            "max_counts": ("load", "max_counts"),
            "min_mito_pct": ("load", "min_mito_pct"),
            "max_mito_pct": ("load", "max_mito_pct"),
            "min_mito": ("load", "min_mito_pct"),  # alias
            "max_mito": ("load", "max_mito_pct"),  # alias
            # Phenotype
            "cd4_cd8_ratio": ("phenotype", "cd4_cd8_ratio"),
            "min_cd3_reads": ("phenotype", "min_cd3_reads"),
            # Clonotype
            "group_by": ("clonotype", "group_by"),
            "handle_doublets": ("clonotype", "handle_doublets"),
            "min_umi": ("clonotype", "min_umi"),
            "doublet_warn_rate": ("clonotype", "doublet_warn_rate"),
            # Attribution (#176). Field names are globally unique so they
            # survive merge_with_args' flatten/rebuild round-trip.
            "attribution": ("attribution", "enabled"),  # --attribution / flat bool
            "enabled": ("attribution", "enabled"),  # from flattened section
            "beta_sharing": ("attribution", "beta_sharing"),
            "doublet_split": ("attribution", "doublet_split"),
            "dual_alpha_merge": ("attribution", "dual_alpha_merge"),
            "dual_alpha_min_cells": ("attribution", "dual_alpha_min_cells"),
            "dual_alpha_min_fraction": ("attribution", "dual_alpha_min_fraction"),
            "dual_alpha_max_beta_partners": ("attribution", "dual_alpha_max_beta_partners"),
            "share_weighting": ("attribution", "share_weighting"),
            "em_iterations": ("attribution", "em_iterations"),
            "dual_alpha_split": ("attribution", "dual_alpha_split"),
            "dual_beta_split": ("attribution", "dual_beta_split"),
            # Filter
            "method": ("filter", "method"),
            "tcell_type": ("filter", "tcell_type"),
            "min_cells": ("filter", "min_cells"),
            "min_frequency": ("filter", "min_frequency"),
            "require_complete": ("filter", "require_complete"),
            "fdr_tiers": ("filter", "fdr_tiers"),
            "filter_mode": ("filter", "filter_mode"),
            "min_donors": ("filter", "min_donors"),
            "min_methods_per_donor": ("filter", "min_methods_per_donor"),
            "min_cells_per_method": ("filter", "min_cells_per_method"),
            "min_frequency_per_method": ("filter", "min_frequency_per_method"),
            "min_timepoints": ("filter", "min_timepoints"),
            "min_timepoints_per_donor": ("filter", "min_timepoints_per_donor"),
            "min_apcs": ("filter", "min_apcs"),
            "min_apcs_per_donor": ("filter", "min_apcs_per_donor"),
            "min_til_cells_per_donor": ("filter", "min_til_cells_per_donor"),
            "fdr_scope": ("filter", "fdr_scope"),
            "donors_share_antigen": ("filter", "donors_share_antigen"),
            # Annotate
            "vdjdb_path": ("annotate", "vdjdb_path"),
            "iedb_path": ("annotate", "iedb_path"),
            "cedar_path": ("annotate", "cedar_path"),
            "vdjdb": ("annotate", "vdjdb_path"),  # alias
            "iedb": ("annotate", "iedb_path"),  # alias
            "cedar": ("annotate", "cedar_path"),  # alias
            "match_by": ("annotate", "match_by"),
            "exclude_viral": ("annotate", "exclude_viral"),
            "flag_only": ("annotate", "flag_only"),
            # TIL
            "til_match_by": ("til", "match_by"),
            "min_til_cells": ("til", "min_til_cells"),
            "til_samples": ("til", "til_samples"),
            # SCT
            "min_snr": ("sct", "min_snr"),
            "min_reads_per_chain": ("sct", "min_reads_per_chain"),
            "require_mutation_match": ("sct", "require_mutation_match"),
            "require_compact_match": ("sct", "require_compact_match"),
            # GEX
            "gene_list": ("gex", "gene_list"),
            "gene_groups": ("gex", "gene_groups"),
            "include_qc": ("gex", "include_qc"),
            "aggregation_ops": ("gex", "aggregation_ops"),
            # Unify
            "add_occurrence_flags": ("unify", "add_occurrence_flags"),
            "add_combined_stats": ("unify", "add_combined_stats"),
            "add_phenotype_confidence": ("unify", "add_phenotype_confidence"),
            "phenotype_ratio_threshold": ("unify", "phenotype_ratio_threshold"),
            # Assemble
            "alpha_leader": ("assemble", "alpha_leader"),
            "beta_leader": ("assemble", "beta_leader"),
            "leader_fallback": ("assemble", "leader_fallback"),
            "force_alpha_leader": ("assemble", "force_alpha_leader"),
            "force_beta_leader": ("assemble", "force_beta_leader"),
            "secondary_alpha_leader": ("assemble", "secondary_alpha_leader"),
            "secondary_beta_leader": ("assemble", "secondary_beta_leader"),
            "include_constant": ("assemble", "include_constant"),
            "constant_source": ("assemble", "constant_source"),
            "linker": ("assemble", "linker"),
            "contigs_dir": ("assemble", "contigs_dir"),
            "sample_name_from": ("assemble", "sample_name_from"),
            "cellranger_dir": ("assemble", "cellranger_dir"),
            "single_chain": ("assemble", "single_chain"),
            # Output
            "generate_plots": ("output", "generate_plots"),
            "emit_clone_sample_long": ("output", "emit_clone_sample_long"),
            "per_method_top_n": ("output", "per_method_top_n"),
            "method_overlap_similarity": ("output", "method_overlap_similarity"),
            "plot_format": ("output", "plot_format"),
            "baseline_markers": ("output", "baseline_markers"),
            "generate_report": ("output", "generate_report"),
            "report_format": ("output", "report_format"),
            "output_airr": ("output", "output_airr"),
            "output_fasta": ("output", "output_fasta"),
            "skip_plots": ("output", "generate_plots"),  # inverted alias
        }

        # Initialize nested dictionaries
        nested: dict[str, dict[str, Any]] = {
            "load": {},
            "phenotype": {},
            "clonotype": {},
            "attribution": {},
            "filter": {},
            "annotate": {},
            "til": {},
            "sct": {},
            "gex": {},
            "unify": {},
            "assemble": {},
            "output": {},
        }
        global_opts: dict[str, Any] = {}

        for key, value in data.items():
            # A key that names a section AND is given as a dict is a nested
            # block, even when a same-named flat scalar alias exists (e.g.
            # 'attribution', which is both a section and a --attribution bool).
            if key in nested and isinstance(value, dict):
                nested[key].update(value)
            elif key in flat_to_nested:
                section, param = flat_to_nested[key]
                # Handle inverted aliases
                if key == "skip_plots":
                    value = not value
                nested[section][param] = value
            elif key in nested:
                # Nested section provided directly
                if isinstance(value, dict):
                    nested[key].update(value)
            elif key == "verbose":
                global_opts["verbose"] = value
            elif key == "selection":
                # Free-form nested selection-rule config (#122).
                global_opts["selection"] = value if isinstance(value, dict) else {}
            elif key == "insilico_filter":
                # Free-form in-silico filter layer (#149).
                global_opts["insilico_filter"] = (
                    value if isinstance(value, dict) else {}
                )
            # Ignore unknown keys

        return cls(
            load=LoadConfig(**nested["load"]),
            phenotype=PhenotypeConfig(**nested["phenotype"]),
            clonotype=ClonotypeConfig(**nested["clonotype"]),
            attribution=AttributionConfig(**nested["attribution"]),
            filter=FilterConfig(**nested["filter"]),
            annotate=AnnotateConfig(**nested["annotate"]),
            til=TILConfig(**nested["til"]),
            sct=SCTConfig(**nested["sct"]),
            gex=GEXConfig(**nested["gex"]),
            unify=UnifyConfig(**nested["unify"]),
            assemble=AssembleConfig(**nested["assemble"]),
            output=OutputConfig(**nested["output"]),
            **global_opts,
        )

    def to_yaml(self, path: str | Path) -> None:
        """
        Save configuration to a YAML file.

        Parameters
        ----------
        path : str or Path
            Path to save the configuration
        """
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
            "unify": dataclasses.asdict(self.unify),
            "assemble": dataclasses.asdict(self.assemble),
            "output": dataclasses.asdict(self.output),
            "selection": self.selection,
            "insilico_filter": self.insilico_filter,
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

        # Start with current config as dict
        config_dict = self.to_dict()

        # Flatten config for easier merging. The `selection` block is a
        # free-form nested dict (rules/global_rank), not a flat namespace
        # of scalar params, so keep it intact rather than spreading it —
        # otherwise its nested keys get hoisted and dropped on rebuild.
        flat_config = {}
        for section, params in config_dict.items():
            # Free-form nested dicts (rules / predicate lists), not flat
            # scalar-param namespaces — keep intact rather than spreading them,
            # else their nested keys get hoisted and dropped on rebuild.
            if section in ("selection", "insilico_filter"):
                flat_config[section] = params
            elif isinstance(params, dict):
                # `match_by` is the one field name shared across sections
                # (annotate + til). Flattening by bare name would let til's
                # value clobber annotate's (and til's would reset to default
                # on rebuild). Emit til.match_by under its unambiguous flat
                # alias so both round-trip independently.
                if section == "til" and "match_by" in params:
                    params = {
                        ("til_match_by" if k == "match_by" else k): v
                        for k, v in params.items()
                    }
                flat_config.update(params)
            else:
                flat_config[section] = params

        # Handle leader shortcut flags first
        if args_dict.get("no_leaders"):
            flat_config["alpha_leader"] = None
            flat_config["beta_leader"] = None
        elif args_dict.get("leaders_from_contigs"):
            flat_config["alpha_leader"] = "from_contig"
            flat_config["beta_leader"] = "from_contig"

        # Apply CLI overrides
        for key, value in args_dict.items():
            # Skip non-config args and shortcut flags
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
            ):
                continue
            # Handle special cases
            if key == "fdr_tiers" and isinstance(value, str):
                value = [float(x) for x in value.split(",")]
            if key == "til_samples" and isinstance(value, str):
                value = [x.strip() for x in value.split(",")]
            # Convert "none" string to None for leader options
            if key in ("alpha_leader", "beta_leader") and value == "none":
                value = None
            flat_config[key] = value

        return TCRsiftConfig._from_dict(flat_config)


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
