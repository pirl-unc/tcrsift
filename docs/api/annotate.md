# Database Annotation & TCRdist API

Match clonotypes against VDJdb, IEDB, and CEDAR epitope references using
exact CDR3 matching (`match_mode="exact"`), paired/β-only Levenshtein-1
(`match_mode="levenshtein"`), or native zero-dependency Dash/`tcrdist3`
V-gene loop + BLOSUM62 CDR3 similarity (`match_mode="tcrdist"`). Annotation
separates canonical protein/antigen labels (`db_protein_canonical`, e.g.
`NY-ESO-1`, `LAGE-1`, `MART-1`, `TRP-2`, `CEA`) from official HGNC/MGI
encoding gene symbols (`db_antigen_gene`, e.g. `CTAG1A`/`CTAG1B`, `CTAG2`,
`MLANA`/`Mlana`, `DCT`/`Dct`, `CEACAM5`; `None` for viral or exogenous
proteins such as chicken OVA), and tracks donor host species
(`db_host_species`, `host_species_filter`) alongside distance-aware match
confidence (`db_match_score`).

## Database loaders and matching

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

## Native TCRdist similarity (`tcrsift.tcrdist`)

Pure Python + NumPy implementation of the Dash (*Nature* 2017) / `tcrdist3`
(*eLife* 2021) metric for human and mouse αβ or β-only TCRs, backed by 448
bundled IMGT V-gene CDR1/CDR2/CDR2.5 loops (`tcrsift/refseqs/vgene_cdr_loops.csv.gz`):

::: tcrsift.tcrdist
    options:
      members:
        - trim_cdr3
        - cdr3_distance
        - v_gene_distance
        - tcr_distance
        - TCRDistIndex
