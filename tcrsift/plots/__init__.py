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
"""Visualization subpackage for TCRsift.

Generates plots across repertoire QC, phenotyping, clonotyping, structured
annotations, rank/slope dynamics, single-cell atlases, selection/PRISM
comparison, assembly QC, and combined PDF/HTML reports.
"""

from __future__ import annotations

import seaborn as sns  # noqa: F401

from ..format import (  # noqa: F401
    pretty_method,
    pretty_methods,
    pretty_sample,
    pretty_samples,
)
from ..signatures import (  # noqa: F401
    ACTIVATION_GENES_HGNC,
    AIM_GENES_HGNC,
    ANTIGEN_RESPONSE_GENES_HGNC,
    CYTOLYTIC_GENES_HGNC,
    EXHAUSTION_GENES_HGNC,
    EXPANSION_CORE_GENES_HGNC,
    TUMOR_REACTIVE_GENES_HGNC,
)
from .annotations import (  # noqa: F401
    GRANULARITY_COLUMNS,
    _expand_sample_frequencies,
    _resolve_granularity,
    _top_n_clones,
    plot_annotations,
    plot_annotations_phase_b,
    plot_category_composition_bar,
    plot_match_strength_comparison,
    plot_matched_clone_heatmap,
)
from .atlas import (  # noqa: F401
    _categorical_colors,
    _coerce_row_mask,
    _resolve_cell_values,
    _resolve_subset_mask,
    _shared_umap_scale,
    _umap_coords,
    _violin_horizontal,
    plot_provenance,
    plot_raincloud,
    plot_signature_map,
    plot_signature_vs_background,
    plot_umap,
    plot_umap_facets,
)
from .ranks import (  # noqa: F401
    REFERENCE_FRACTIONS,
    _derive_n_methods,
    _per_sample_rank_table,
    draw_reference_fractions,
    make_method_panel_grid,
    plot_clone_freq_vs_signature_per_sample,
    plot_clone_tracking_slopes,
    plot_clones_by_n_methods,
    plot_cumulative_coverage,
    plot_rank_curves_per_sample,
    plot_rank_family,
    plot_top_n_cumulative_stacked,
    plot_top_n_labeled_bars_per_sample,
)
from .repertoire import (  # noqa: F401
    plot_assembly_qc,
    plot_cells_per_sample_stacked,
    plot_clonotypes,
    plot_filter,
    plot_pgen_distribution,
    plot_phenotype,
    plot_publicness_vs_match_score,
    plot_qc,
)
from .reports import (  # noqa: F401
    _TCR_SEQUENCE_COLUMN_CANDIDATES,
    FUNNEL_LABEL_NICE,
    _default_tcr_sequence_columns,
    _draw_funnel_side_labels,
    _expand_annotation_lines,
    _log_stage_widths,
    create_pipeline_funnel,
    create_tcr_sequence_pdf,
    generate_report,
    normalize_funnel_label,
    plot_funnel,
    plot_funnel_lollipop,
    plot_funnel_ribbon,
    plot_funnel_terrace,
    plot_funnel_variants,
)
from .selection import (  # noqa: F401
    _FREQ_CBAR_W,
    _FREQ_CELL_H,
    _FREQ_CELL_W,
    _FREQ_HEATMAP_PAGE_ROWS,
    _FREQ_LABEL_W,
    _FREQ_STRIP_W,
    _FREQ_XFONT,
    _FREQ_YFONT,
    _render_frequency_heatmap_page,
    plot_assembly,
    plot_cross_donor_venn,
    plot_freq_prism_grid,
    plot_freq_prism_scatter,
    plot_leader_summary,
    plot_method_overlap,
    plot_method_recovery,
    plot_prism_components,
    plot_selected_frequency_heatmap,
    plot_set_overlap,
    plot_til,
    plot_vgene_usage_by_method,
)
from .state import (  # noqa: F401
    DB_CATEGORY_PALETTE,
    VALID_PLOT_FORMATS,
    _category_color_strip,
    _grid_shape,
    _parse_sample_frequencies,
    _per_cell_signature,
    get_plot_format,
    save_figure,
    set_plot_format,
    set_polished_style,
)

__all__ = [
    "ACTIVATION_GENES_HGNC",
    "AIM_GENES_HGNC",
    "ANTIGEN_RESPONSE_GENES_HGNC",
    "CYTOLYTIC_GENES_HGNC",
    "DB_CATEGORY_PALETTE",
    "EXHAUSTION_GENES_HGNC",
    "EXPANSION_CORE_GENES_HGNC",
    "FUNNEL_LABEL_NICE",
    "GRANULARITY_COLUMNS",
    "REFERENCE_FRACTIONS",
    "TUMOR_REACTIVE_GENES_HGNC",
    "VALID_PLOT_FORMATS",
    "create_pipeline_funnel",
    "create_tcr_sequence_pdf",
    "draw_reference_fractions",
    "generate_report",
    "get_plot_format",
    "make_method_panel_grid",
    "normalize_funnel_label",
    "plot_annotations",
    "plot_annotations_phase_b",
    "plot_assembly",
    "plot_assembly_qc",
    "plot_category_composition_bar",
    "plot_cells_per_sample_stacked",
    "plot_clone_freq_vs_signature_per_sample",
    "plot_clone_tracking_slopes",
    "plot_clones_by_n_methods",
    "plot_clonotypes",
    "plot_cross_donor_venn",
    "plot_cumulative_coverage",
    "plot_filter",
    "plot_freq_prism_grid",
    "plot_freq_prism_scatter",
    "plot_funnel",
    "plot_funnel_lollipop",
    "plot_funnel_ribbon",
    "plot_funnel_terrace",
    "plot_funnel_variants",
    "plot_leader_summary",
    "plot_match_strength_comparison",
    "plot_matched_clone_heatmap",
    "plot_method_overlap",
    "plot_method_recovery",
    "plot_pgen_distribution",
    "plot_phenotype",
    "plot_prism_components",
    "plot_provenance",
    "plot_publicness_vs_match_score",
    "plot_qc",
    "plot_raincloud",
    "plot_rank_curves_per_sample",
    "plot_rank_family",
    "plot_selected_frequency_heatmap",
    "plot_set_overlap",
    "plot_signature_map",
    "plot_signature_vs_background",
    "plot_til",
    "plot_top_n_cumulative_stacked",
    "plot_top_n_labeled_bars_per_sample",
    "plot_umap",
    "plot_umap_facets",
    "plot_vgene_usage_by_method",
    "pretty_method",
    "pretty_methods",
    "pretty_sample",
    "pretty_samples",
    "save_figure",
    "set_plot_format",
    "set_polished_style",
]
