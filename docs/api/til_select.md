# Longitudinal TIL Selection API

Harmonize per-timepoint 10x VDJ and GEX files in the legacy v2 layout and apply
multi-branch TIL selection masks.

::: tcrsift.til_select
    options:
      members:
        - parse_sample_args
        - parse_config
        - default_timepoint_inputs
        - load_from_consensus
        - compute_marker_scores_for_timepoint
        - build_harmonized_table
        - run_selection_pipeline
        - run_til_select
