# TIL Matching API

Load tumor-infiltrating lymphocyte (TIL) repertoires from CSV, `.h5ad`, or Cell
Ranger VDJ outputs and cross-reference them against culture clonotypes.

::: tcrsift.til
    options:
      members:
        - load_til_data
        - load_til_samples
        - load_til_specs
        - summarize_til_clonotypes
        - match_til
        - get_til_summary
        - identify_til_specific_clones
