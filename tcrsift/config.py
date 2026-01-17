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
    min_mito_pct: float = 2.0
    max_mito_pct: float = 8.0


@dataclass
class PhenotypeConfig:
    """Configuration for the phenotype step."""

    cd4_cd8_ratio: float = 3.0
    min_cd3_reads: int = 10


@dataclass
class ClonotypeConfig:
    """Configuration for the clonotype step."""

    group_by: str = "CDR3ab"
    handle_doublets: str = "flag"
    min_umi: int = 2


@dataclass
class FilterConfig:
    """Configuration for the filter step."""

    method: str = "threshold"
    tcell_type: str = "cd8"
    min_cells: int = 2
    min_frequency: float = 0.0
    require_complete: bool = True
    fdr_tiers: list[float] = field(default_factory=lambda: [0.15, 0.1, 0.01, 0.001, 0.0001])


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

    alpha_leader: str | None = "CD28"  # None, "from_contig", or key: CD8A, CD28, IgK, TRAC, TRBC
    beta_leader: str | None = "CD8A"   # None, "from_contig", or key: CD8A, CD28, IgK, TRAC, TRBC
    include_constant: bool = True
    constant_source: str = "ensembl"
    linker: str = "T2A"
    contigs_dir: str | None = None
    single_chain: bool = True


@dataclass
class OutputConfig:
    """Configuration for output options."""

    generate_plots: bool = True
    generate_report: bool = True
    report_format: str = "pdf"
    output_airr: bool = False
    output_fasta: bool = False


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
    filter: FilterConfig = field(default_factory=FilterConfig)
    annotate: AnnotateConfig = field(default_factory=AnnotateConfig)
    til: TILConfig = field(default_factory=TILConfig)
    sct: SCTConfig = field(default_factory=SCTConfig)
    gex: GEXConfig = field(default_factory=GEXConfig)
    unify: UnifyConfig = field(default_factory=UnifyConfig)
    assemble: AssembleConfig = field(default_factory=AssembleConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

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
            # Filter
            "method": ("filter", "method"),
            "tcell_type": ("filter", "tcell_type"),
            "min_cells": ("filter", "min_cells"),
            "min_frequency": ("filter", "min_frequency"),
            "require_complete": ("filter", "require_complete"),
            "fdr_tiers": ("filter", "fdr_tiers"),
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
            # Amplify
            "min_snr": ("amplify", "min_snr"),
            "min_reads_per_chain": ("amplify", "min_reads_per_chain"),
            "require_mutation_match": ("amplify", "require_mutation_match"),
            "require_compact_match": ("amplify", "require_compact_match"),
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
            "include_constant": ("assemble", "include_constant"),
            "constant_source": ("assemble", "constant_source"),
            "linker": ("assemble", "linker"),
            "contigs_dir": ("assemble", "contigs_dir"),
            "single_chain": ("assemble", "single_chain"),
            # Output
            "generate_plots": ("output", "generate_plots"),
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
            if key in flat_to_nested:
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
            # Ignore unknown keys

        return cls(
            load=LoadConfig(**nested["load"]),
            phenotype=PhenotypeConfig(**nested["phenotype"]),
            clonotype=ClonotypeConfig(**nested["clonotype"]),
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
            "filter": dataclasses.asdict(self.filter),
            "annotate": dataclasses.asdict(self.annotate),
            "til": dataclasses.asdict(self.til),
            "sct": dataclasses.asdict(self.sct),
            "gex": dataclasses.asdict(self.gex),
            "unify": dataclasses.asdict(self.unify),
            "assemble": dataclasses.asdict(self.assemble),
            "output": dataclasses.asdict(self.output),
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

        # Flatten config for easier merging
        flat_config = {}
        for section, params in config_dict.items():
            if isinstance(params, dict):
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
            if key in ("func", "command", "config", "sample_sheet", "input", "output", "output_dir",
                      "no_leaders", "leaders_from_contigs"):
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
    print(f"Generated example config: {path}")
