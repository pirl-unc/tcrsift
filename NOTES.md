# til-select dynamic notes

## v2 input contract observed
- Per-timepoint VDJ inputs:
  - `consensus_annotations.<TP>.csv`
  - `clonotypes.<TP>.csv`
- Per-timepoint marker/GEX inputs:
  - `filtered_contig_annotations.<TP>.csv`
  - `sample_filtered_feature_bc_matrix.<TP>.h5`

## v2 sample mapping format
- `LABEL=CONSENSUS_PATH,CLONOTYPES_PATH`
- Also supports YAML config with `timepoints:` mapping.

## Existing tcrsift overlap
- `til-clonotype` handles TIL-only clonotype aggregation, but not:
  - timepoint harmonization
  - marker CP10K/z scoring
  - multi-branch clone prioritization masks
- `annotate_clonotypes` now supports optional no-DB mode; reuse this.
- `load_databases` + `match_clonotypes` can be reused for DB annotation.

## Design decision
- Implement new high-level command: `til-select`.
- Keep `til-clonotype` low-level and unchanged.
- `til-select` should emit v2-style files (`subset_*.csv`, `selection_masks.csv`, top-k tables, etc.).

## Open checks while implementing
- Confirm count-column auto-detection behavior against v2 defaults.
- Confirm exact output names for top-k rank variants.
- Keep plotting/report minimal but compatible in filenames.
