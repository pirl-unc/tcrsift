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
TCRsift: TCR selection from antigen-specific culture and scRNA/VDJ sequencing data.

The public API is resolved lazily so light-weight imports such as
`tcrsift.sample_sheet` or `from tcrsift import load_sample_sheet` do not pull in
heavy optional runtime dependencies like Scanpy unless they are actually needed.
"""

from __future__ import annotations

from importlib import import_module

_MODULE_EXPORTS = {
    ".version": ("__version__",),
    ".config": (
        "TCRsiftConfig",
        "LoadConfig",
        "AssembleConfig",
        "SCTConfig",
        "GEXConfig",
        "UnifyConfig",
    ),
    ".sample_sheet": (
        "Sample",
        "SampleSheet",
        "load_sample_sheet",
        "validate_sample_sheet",
    ),
    ".loader": (
        "load_cellranger_vdj",
        "load_cellranger_gex",
        "load_sample",
        "load_samples",
    ),
    ".sct": (
        "load_sct",
        "aggregate_sct",
        "get_sct_specificities",
    ),
    ".gex": (
        "augment_with_gex",
        "aggregate_gex_by_clonotype",
        "compute_cd4_cd8_counts",
        "DEFAULT_GENE_LIST",
        "DEFAULT_GENE_GROUPS",
    ),
    ".phenotype": (
        "phenotype_cells",
        "classify_tcell_type",
        "filter_by_tcell_type",
        "get_phenotype_summary",
    ),
    ".clonotype": (
        "aggregate_clonotypes",
        "get_clonotype_summary",
        "export_clonotypes_airr",
    ),
    ".filter": (
        "filter_clonotypes",
        "filter_clonotypes_threshold",
        "assign_tiers_threshold",
        "split_by_tier",
        "get_filter_summary",
    ),
    ".annotate": (
        "load_vdjdb",
        "load_iedb",
        "load_cedar",
        "annotate_clonotypes",
        "get_annotation_summary",
    ),
    ".til": (
        "match_til",
        "get_til_summary",
        "identify_til_specific_clones",
        "load_til_specs",
        "summarize_til_clonotypes",
    ),
    ".til_select": (
        "load_from_consensus",
        "compute_marker_scores_for_timepoint",
        "build_harmonized_table",
        "run_selection_pipeline",
        "run_til_select",
    ),
    ".unify": (
        "merge_experiments",
        "add_phenotype_confidence",
        "compute_condition_statistics",
        "find_top_condition",
        "get_unify_summary",
    ),
    ".assemble": (
        "DEFAULT_LEADERS",
        "LINKERS",
        "assemble_full_sequences",
        "translate_dna",
        "validate_sequences",
        "export_fasta",
    ),
    ".plots": (
        "plot_funnel",
        "create_pipeline_funnel",
        "create_tcr_sequence_pdf",
    ),
    ".qc": (
        "QCReport",
        "QCResult",
        "find_repeated_kmers",
        "validate_sequence",
        "validate_clonotypes",
        "get_qc_summary",
    ),
    ".mnemonic": ("tcr_name",),
    ".validation": (
        "TCRsiftValidationError",
        "validate_cdr3_sequence",
        "validate_cdr3_dataframe",
        "safe_divide",
        "safe_percentage",
        "safe_mode",
    ),
}

_EXPORT_TO_MODULE = {
    export_name: module_name
    for module_name, export_names in _MODULE_EXPORTS.items()
    for export_name in export_names
}

__all__ = [
    "__version__",
    "TCRsiftConfig",
    "LoadConfig",
    "AssembleConfig",
    "SCTConfig",
    "GEXConfig",
    "UnifyConfig",
    "Sample",
    "SampleSheet",
    "load_sample_sheet",
    "validate_sample_sheet",
    "load_cellranger_vdj",
    "load_cellranger_gex",
    "load_sample",
    "load_samples",
    "load_sct",
    "aggregate_sct",
    "get_sct_specificities",
    "augment_with_gex",
    "aggregate_gex_by_clonotype",
    "compute_cd4_cd8_counts",
    "DEFAULT_GENE_LIST",
    "DEFAULT_GENE_GROUPS",
    "phenotype_cells",
    "classify_tcell_type",
    "filter_by_tcell_type",
    "get_phenotype_summary",
    "aggregate_clonotypes",
    "get_clonotype_summary",
    "export_clonotypes_airr",
    "filter_clonotypes",
    "filter_clonotypes_threshold",
    "assign_tiers_threshold",
    "split_by_tier",
    "get_filter_summary",
    "load_vdjdb",
    "load_iedb",
    "load_cedar",
    "annotate_clonotypes",
    "get_annotation_summary",
    "match_til",
    "get_til_summary",
    "identify_til_specific_clones",
    "load_til_specs",
    "summarize_til_clonotypes",
    "load_from_consensus",
    "compute_marker_scores_for_timepoint",
    "build_harmonized_table",
    "run_selection_pipeline",
    "run_til_select",
    "merge_experiments",
    "add_phenotype_confidence",
    "compute_condition_statistics",
    "find_top_condition",
    "get_unify_summary",
    "DEFAULT_LEADERS",
    "LINKERS",
    "assemble_full_sequences",
    "translate_dna",
    "validate_sequences",
    "export_fasta",
    "plot_funnel",
    "create_pipeline_funnel",
    "create_tcr_sequence_pdf",
    "QCReport",
    "QCResult",
    "find_repeated_kmers",
    "validate_sequence",
    "validate_clonotypes",
    "get_qc_summary",
    "tcr_name",
    "TCRsiftValidationError",
    "validate_cdr3_sequence",
    "validate_cdr3_dataframe",
    "safe_divide",
    "safe_percentage",
    "safe_mode",
]


def __getattr__(name: str):
    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module = import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
