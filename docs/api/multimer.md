# Multimer Specificity API

Extract non-gene `Antibody Capture` / `Antigen Capture` feature barcodes from 10x
Genomics `.h5` or Matrix Market directories, call single-cell pMHC multimer
specificity against negative-control multimers, and compute clone-level
concordance metrics to separate genuine antigen-specific clones from sticky
multimer artifacts.

::: tcrsift.multimer
    options:
      members:
        - extract_feature_barcodes
        - classify_feature_columns
        - call_multimer_specificity
        - summarize_clone_multimer_specificity
