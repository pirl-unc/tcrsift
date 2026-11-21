# Candidate Selection API

`tcrsift.candidate` builds culture shortlists that combine high-confidence
abundance tiers (`tier1` and `tier2`) with top-scoring signature picks from
lower abundance tiers (`tier3`–`tier5`):

```
Selected = tier1 ∪ tier2 ∪ top_n_per_signature(tier3 ∪ tier4 ∪ tier5)
```

## Example

```python
from tcrsift import (
    aggregate_clonotypes,
    build_clone_method_long,
    build_clone_sample_long,
    compute_signature_picks_per_method,
    compute_signature_scores_per_clonotype,
    filter_clonotypes,
    select_candidates,
)

# 1. Score signatures per clonotype (mean log1p expression on CD8 cells),
#    merge onto tiered clonotypes.
sig_scores = compute_signature_scores_per_clonotype(
    per_cell_df,
    group_col="CDR3ab",
    cd8_only=True,
)
clones = filter_clonotypes(aggregate_clonotypes(adata)).merge(sig_scores, on="CDR3ab")

# 2. Mark tier1 ∪ tier2 ∪ top-3 clones per signature from lower tiers.
clones = select_candidates(clones, top_n=3)
shortlist = clones[clones["is_selected"]]

# 3. Optional per-method picks for selection-route heatmaps.
tier3_plus_clones = clones.loc[
    clones["tier"].isin(["tier3", "tier4", "tier5"]), "CDR3ab"
]
clone_method = build_clone_method_long(build_clone_sample_long(adata))
clone_method = clone_method.merge(sig_scores, on="CDR3ab")
per_method_picks = compute_signature_picks_per_method(
    clone_method,
    top_n=1,
    pool_clones=set(tier3_plus_clones),
)
```

## API reference

::: tcrsift.candidate
    options:
      members:
        - select_candidates
        - compute_signature_picks_per_method
        - signature_picks_clone_to_methods
