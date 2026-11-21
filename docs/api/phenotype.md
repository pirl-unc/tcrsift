# Phenotyping API

Classify single cells into CD4+ and CD8+ T-cell categories from `CD3D/E/G`,
`CD4`, and `CD8A/B` expression counts.

::: tcrsift.phenotype
    options:
      members:
        - phenotype_cells
        - classify_tcell_type
        - filter_by_tcell_type
        - get_phenotype_summary
        - validate_phenotype_vs_expected
        - TCELL_TYPE_CATEGORIES
