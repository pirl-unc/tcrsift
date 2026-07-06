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
        "EmbedConfig",
        "UnifyConfig",
    ),
    ".embed": (
        "embed_cells",
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
        "compute_signature_scores_per_cell",
        "DEFAULT_GENE_LIST",
        "DEFAULT_GENE_GROUPS",
    ),
    ".signature_methods": (
        "Signature",
        "SIGNATURES",
        "NEOANTIGEN_SIGNATURES",
        "score_signature",
        "score_by_name",
        "score_weighted_z",
    ),
    ".phenotype": (
        "phenotype_cells",
        "classify_tcell_type",
        "filter_by_tcell_type",
        "get_phenotype_summary",
    ),
    ".annotate_cells": (
        "annotate_cells",
        "annotate_clusters",
        "score_reference",
        "compose_phenotype_labels",
        "top_markers",
        "AnnotationGates",
        "DEFAULT_GATES",
    ),
    ".clonotype": (
        "aggregate_clonotypes",
        "build_clone_method_long",
        "build_clone_sample_long",
        "build_selection_sets",
        "compute_sample_overlap_matrices",
        "get_clonotype_summary",
        "set_overlap_table",
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
        "match_clonotypes",
        "viral_cdr3b_set",
    ),
    ".selection": (
        "build_selection_rules",
        "select_specificity_candidates",
        "select_by_dominant_specificity",
        "select_freq_prism_per_condition",
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
        "add_synthesis_qc",
        "assemble_full_sequences",
        "assemble_qc_report",
        "build_assembly_qc_report",
        "synthesis_qc_report",
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
        "plot_freq_prism_grid",
        "plot_set_overlap",
        "plot_pgen_distribution",
        "plot_publicness_vs_match_score",
        "plot_umap",
        "plot_provenance",
        "plot_signature_vs_background",
        "plot_raincloud",
        "save_figure",
        "set_plot_format",
        "set_polished_style",
        "create_pipeline_funnel",
        "create_tcr_sequence_pdf",
    ),
    ".format": (
        "BASELINE_MARKERS",
        "condition_sort_key",
        "order_conditions",
        "pdf_safe",
        "pretty_method",
        "pretty_methods",
        "pretty_sample",
        "pretty_samples",
        "pretty_antigen",
        "pretty_antigens",
        "set_antigen_labels",
        "set_baseline_markers",
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
        "AIM_GENES_HGNC",
        "CYTOLYTIC_GENES_HGNC",
        "EXHAUSTION_GENES_HGNC",
        "TUMOR_REACTIVE_GENES_HGNC",
        "EXPANSION_CORE_GENES_HGNC",
        "MARKER_PANEL_HGNC",
        "T_CELL_SIGNATURES",
        "MANASCORE_UP_HGNC",
        "MANASCORE_DOWN_HGNC",
        "MANASCORE_WEIGHTS_HGNC",
        "NEOTCR8_GENES_HGNC",
        "NEOTCR4_GENES_HGNC",
        "NEOTCRPBL_GENES_HGNC",
        "CELL_TYPE_SIGNATURES",
        "T_STATE_SIGNATURES",
        "B_STATE_SIGNATURES",
        "LINEAGE_GENES",
        "PBMC_CULTURE_TYPES",
    ),
    ".qc": (
        "QCReport",
        "QCResult",
        "cdr3_anchor_integrity",
        "clonal_expansion_metrics",
        "find_repeated_kmers",
        "gex_vdj_overlap",
        "inter_sample_barcode_jaccard",
        "read_cellranger_metrics",
        "sample_integrity_qc",
        "validate_sequence",
        "validate_clonotypes",
        "get_qc_summary",
        "cell_qc_funnel",
        "cross_lineage_doublets",
        "cd4_cd8_doublet_mask",
        "DEFAULT_LINEAGE_PROGRAMS",
    ),
    ".sort_qc": (
        "DEFAULT_SORT_SIGNATURE_MAP",
        "parse_sort_label",
        "sort_signature_consistency",
        "sort_signature_consistency_from_adata",
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
