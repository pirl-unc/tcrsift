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

A tool for identifying antigen-specific T cell receptor clones from single-cell
sequencing data, with support for:

- Loading CellRanger VDJ and GEX outputs
- CD4/CD8 T cell phenotyping from gene expression
- Clonotype aggregation and frequency analysis
- Tiered filtering for antigen-specific clones
- Annotation with public TCR databases (VDJdb, IEDB, CEDAR)
- TIL (tumor-infiltrating lymphocyte) matching
- Full-length TCR sequence assembly

Example usage::

    # Run complete pipeline
    tcrsift run --sample-sheet samples.yaml --output-dir results/ --report

    # Or run individual steps
    tcrsift load --sample-sheet samples.yaml -o loaded.h5ad
    tcrsift phenotype -i loaded.h5ad -o phenotyped.h5ad
    tcrsift clonotype -i phenotyped.h5ad -o clonotypes.csv
    tcrsift filter -i clonotypes.csv -o filtered/
    tcrsift annotate -i filtered/tier4.csv -o annotated.csv --vdjdb /path/to/vdjdb
    tcrsift assemble -i annotated.csv -o full_sequences.csv --include-constant

"""

from .annotate import (
    annotate_clonotypes,
    get_annotation_summary,
    load_cedar,
    load_iedb,
    load_vdjdb,
)
from .assemble import (
    DEFAULT_LEADERS,
    LINKERS,
    assemble_full_sequences,
    export_fasta,
    translate_dna,
    validate_sequences,
)
from .clonotype import (
    aggregate_clonotypes,
    export_clonotypes_airr,
    get_clonotype_summary,
)
from .config import (
    AssembleConfig,
    GEXConfig,
    LoadConfig,
    SCTConfig,
    TCRsiftConfig,
    UnifyConfig,
)
from .filter import (
    assign_tiers_threshold,
    filter_clonotypes,
    filter_clonotypes_threshold,
    get_filter_summary,
    split_by_tier,
)
from .gex import (
    DEFAULT_GENE_GROUPS,
    DEFAULT_GENE_LIST,
    aggregate_gex_by_clonotype,
    augment_with_gex,
    compute_cd4_cd8_counts,
)
from .loader import (
    load_cellranger_gex,
    load_cellranger_vdj,
    load_sample,
    load_samples,
)
from .mnemonic import tcr_name
from .phenotype import (
    classify_tcell_type,
    filter_by_tcell_type,
    get_phenotype_summary,
    phenotype_cells,
)
from .plots import (
    create_pipeline_funnel,
    create_tcr_sequence_pdf,
    plot_funnel,
)
from .qc import (
    QCReport,
    QCResult,
    find_repeated_kmers,
    get_qc_summary,
    validate_clonotypes,
    validate_sequence,
)

# Core modules
from .sample_sheet import (
    Sample,
    SampleSheet,
    load_sample_sheet,
    validate_sample_sheet,
)
from .sct import (
    aggregate_sct,
    get_sct_specificities,
    load_sct,
)
from .til import (
    get_til_summary,
    identify_til_specific_clones,
    match_til,
)
from .unify import (
    add_phenotype_confidence,
    compute_condition_statistics,
    find_top_condition,
    get_unify_summary,
    merge_experiments,
)
from .validation import TCRsiftValidationError
from .version import __version__

__all__ = [
    # Version
    "__version__",
    # Configuration
    "TCRsiftConfig",
    "LoadConfig",
    "AssembleConfig",
    "SCTConfig",
    "GEXConfig",
    "UnifyConfig",
    # Sample sheet
    "Sample",
    "SampleSheet",
    "load_sample_sheet",
    "validate_sample_sheet",
    # Loading
    "load_cellranger_vdj",
    "load_cellranger_gex",
    "load_sample",
    "load_samples",
    # SCT (single-cell TCR platform)
    "load_sct",
    "aggregate_sct",
    "get_sct_specificities",
    # GEX
    "augment_with_gex",
    "aggregate_gex_by_clonotype",
    "compute_cd4_cd8_counts",
    "DEFAULT_GENE_LIST",
    "DEFAULT_GENE_GROUPS",
    # Phenotyping
    "phenotype_cells",
    "classify_tcell_type",
    "filter_by_tcell_type",
    "get_phenotype_summary",
    # Clonotyping
    "aggregate_clonotypes",
    "get_clonotype_summary",
    "export_clonotypes_airr",
    # Filtering
    "filter_clonotypes",
    "filter_clonotypes_threshold",
    "assign_tiers_threshold",
    "split_by_tier",
    "get_filter_summary",
    # Annotation
    "load_vdjdb",
    "load_iedb",
    "load_cedar",
    "annotate_clonotypes",
    "get_annotation_summary",
    # TIL
    "match_til",
    "get_til_summary",
    "identify_til_specific_clones",
    # Unify
    "merge_experiments",
    "add_phenotype_confidence",
    "compute_condition_statistics",
    "find_top_condition",
    "get_unify_summary",
    # Assembly
    "DEFAULT_LEADERS",
    "LINKERS",
    "assemble_full_sequences",
    "translate_dna",
    "validate_sequences",
    "export_fasta",
    # Plots
    "plot_funnel",
    "create_pipeline_funnel",
    "create_tcr_sequence_pdf",
    # QC
    "QCReport",
    "QCResult",
    "find_repeated_kmers",
    "validate_sequence",
    "validate_clonotypes",
    "get_qc_summary",
    # Utilities
    "tcr_name",
    # Exceptions
    "TCRsiftValidationError",
]
