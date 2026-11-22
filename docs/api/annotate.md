# Database Annotation API

Match clonotypes against VDJdb, IEDB, and CEDAR epitope references, separate
canonical protein/antigen labels (`db_protein_canonical`, e.g. `NY-ESO-1`,
`LAGE-1`, `MART-1`, `TRP-2`, `CEA`) from official HGNC/MGI encoding gene symbols
(`db_antigen_gene`, e.g. `CTAG1A`/`CTAG1B`, `CTAG2`, `MLANA`/`Mlana`,
`DCT`/`Dct`, `CEACAM5`; `None` for viral or exogenous proteins such as chicken
OVA), and track donor host species (`db_host_species`, `host_species_filter`).

::: tcrsift.annotate
    options:
      members:
        - load_vdjdb
        - load_iedb
        - load_cedar
        - load_databases
        - canonicalize_antigen
        - resolve_antigen_gene
        - match_clonotypes
        - annotate_clonotypes
        - get_annotation_summary
        - VIRAL_SPECIES_PATTERNS
