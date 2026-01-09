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

from .version import __version__

# Core modules
from .sample_sheet import (
    Sample,
    SampleSheet,
    load_sample_sheet,
    validate_sample_sheet,
)

from .loader import (
    load_cellranger_vdj,
    load_cellranger_gex,
    load_sample,
    load_samples,
)

from .phenotype import (
    phenotype_cells,
    classify_tcell_type,
    filter_by_tcell_type,
    get_phenotype_summary,
)

from .clonotype import (
    aggregate_clonotypes,
    get_clonotype_summary,
    export_clonotypes_airr,
)

from .filter import (
    filter_clonotypes,
    filter_clonotypes_threshold,
    assign_tiers_threshold,
    split_by_tier,
    get_filter_summary,
)

from .annotate import (
    load_vdjdb,
    load_iedb,
    load_cedar,
    annotate_clonotypes,
    get_annotation_summary,
)

from .til import (
    match_til,
    get_til_summary,
    identify_til_specific_clones,
)

from .assemble import (
    assemble_full_sequences,
    translate_dna,
    validate_sequences,
    export_fasta,
)

from .mnemonic import tcr_name

__all__ = [
    # Version
    "__version__",
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
    # Assembly
    "assemble_full_sequences",
    "translate_dna",
    "validate_sequences",
    "export_fasta",
    # Utilities
    "tcr_name",
]
