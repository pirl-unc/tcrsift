# Clonotyping API

Aggregate paired single-cell VDJ calls into clonotypes and construct long-format
sample, method, overlap, and recovery tables.

::: tcrsift.clonotype
    options:
      members:
        - aggregate_clonotypes
        - build_clone_sample_long
        - build_clone_method_long
        - compute_sample_overlap_matrices
        - build_per_method_rankings
        - build_method_overlap_matrices
        - build_method_recovery_table
        - get_clonotype_summary
        - export_clonotypes_airr
        - calculate_clone_frequencies
