# Database Annotation API

Match clonotypes against VDJdb, IEDB, and CEDAR epitope references and flag
cataloged viral specificities.

::: tcrsift.annotate
    options:
      members:
        - load_vdjdb
        - load_iedb
        - load_cedar
        - load_databases
        - match_clonotypes
        - annotate_clonotypes
        - get_annotation_summary
        - VIRAL_SPECIES_PATTERNS
