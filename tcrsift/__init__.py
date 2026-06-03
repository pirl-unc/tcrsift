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
        "compute_signature_scores_per_clonotype",
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
        "build_clone_method_long",
        "build_clone_sample_long",
        "compute_sample_overlap_matrices",
        "get_clonotype_summary",
        "export_clonotypes_airr",
    ),
    ".candidate": (
        "compute_signature_picks_per_method",
        "select_candidates",
        "signature_picks_clone_to_methods",
    ),
    ".filter": (
        "DEFAULT_THRESHOLD_TIERS",
        "DEFAULT_FDR_TIERS",
        "clone_clears_tier",
        "filter_clonotypes",
        "filter_clonotypes_threshold",
        "assign_tiers_threshold",
        "split_by_tier",
        "strictest_tier_met",
        "per_sample_tier",
        "get_filter_summary",
    ),
    ".annotate": (
        "CATEGORY_BACTERIAL",
        "CATEGORY_CONTRADICTORY",
        "CATEGORY_OTHER",
        "CATEGORY_SELF",
        "CATEGORY_TUMOR_SELF",
        "CATEGORY_UNKNOWN",
        "CATEGORY_VIRAL",
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
        "AssemblyQCCheck",
        "AssemblyQCReport",
        "DEFAULT_LEADERS",
        "LINKERS",
        "assemble_full_sequences",
        "assemble_qc_report",
        "build_assembly_qc_report",
        "translate_dna",
        "validate_sequences",
        "export_fasta",
    ),
    ".plots": (
        "FUNNEL_LABEL_NICE",
        "normalize_funnel_label",
        "plot_funnel",
        "plot_assembly_qc",
        "plot_cells_per_sample_stacked",
        "plot_clone_freq_vs_signature_per_sample",
        "plot_pgen_distribution",
        "plot_publicness_vs_match_score",
        "create_pipeline_funnel",
        "create_tcr_sequence_pdf",
    ),
    ".format": (
        "pretty_method",
        "pretty_methods",
        "pretty_sample",
        "pretty_samples",
    ),
    ".pgen": (
        "annotate_publicness",
        "compute_pgen",
        "pgen_components",
        "pgen_single",
        "publicness_score",
    ),
    ".insilico_filter": (
        "FilterPredicate",
        "apply_insilico_filter",
        "average_percentile_rank",
        "expand_insilico_twins",
        "insilico_mask",
        "insilico_overlap_long",
        "percentile_rank",
        "predicates_from_config",
    ),
    ".olga_ppost": (
        "annotate_nearest_supported_allele",
        "compute_pgen_ppost",
        "flag_private_candidates",
        "load_chain_model",
        "nearest_supported_allele",
        "normalize_gene_name",
        "olga_sonia_available",
        "supported_alleles",
    ),
    ".seqprob": (
        "GeneAwareKmerModel",
        "KmerProbabilityModel",
        "SequenceProbabilityModel",
        "TCRpegProbabilityModel",
        "load_background_model",
        "score_log_pgen",
        "score_log_prob",
        "score_log_q",
    ),
    ".pgen_models": (
        "ensure_model",
        "fetch_repertoire",
        "train_model",
    ),
    ".annotate_tcrs": (
        "add_gex_signature_scores",
        "add_paired_ppost",
        "add_pgen_ppost",
        "annotate_tcrs",
        "naive_signature",
        "prism_score",
        "score_gex_signature_per_clone",
        "select_prism",
    ),
    ".signatures": (
        "ACTIVATION_GENES_HGNC",
        "EFFECTOR_GENES_HGNC",
        "NAIVE_STEM_GENES_HGNC",
        "ANTIGEN_RESPONSE_GENES_HGNC",
        "CYTOLYTIC_GENES_HGNC",
        "EXHAUSTION_GENES_HGNC",
        "TUMOR_REACTIVE_GENES_HGNC",
        "T_CELL_SIGNATURES",
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
        "most_common",
        "pick_representative_cell",
    ),
}

_EXPORT_TO_MODULE = {
    export_name: module_name
    for module_name, export_names in _MODULE_EXPORTS.items()
    for export_name in export_names
}

__all__ = list(_EXPORT_TO_MODULE)


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
