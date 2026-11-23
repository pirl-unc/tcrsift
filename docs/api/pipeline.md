# Pipeline & Prioritization API

Programmatic Python entry points backing `tcrsift run` (`tcrsift.pipeline`) and `tcrsift prioritize` (`tcrsift.prioritize` / `tcrsift.til_prioritize`).

::: tcrsift.pipeline
    options:
      members:
        - PipelineResult
        - run_pipeline

::: tcrsift.prioritize
    options:
      members:
        - PrioritizeConfig
        - resolve_signatures
        - signature_pass_mask
        - select_round_robin

::: tcrsift.til_prioritize
    options:
      members:
        - run_til_prioritize
