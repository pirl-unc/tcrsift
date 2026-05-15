# Candidate Selection API

Signature-based shortlist generation (#75 / #84). After the tier
cascade has labeled clones at varying FDR thresholds, the typical
follow-up is a shortlist that combines strict-FDR picks with
phenotype-based rescues from the lower-tier pool:

```
Selected = tier1 ∪ tier2 ∪ top-N-per-signature(tier3 ∪ tier4 ∪ tier5)
```

## Workflow

```python
from tcrsift import (
    aggregate_clonotypes,
    build_clone_method_long,
    compute_signature_scores_per_clonotype,
    select_candidates,
    compute_signature_picks_per_method,
)

# 1) Clonotype-level signature scores (mean log1p, CD8-restricted).
sig_scores = compute_signature_scores_per_clonotype(per_cell_df)
clones = aggregate_clonotypes(adata).merge(sig_scores, on="CDR3_pair")

# 2) Tag tier1 ∪ tier2 ∪ top-N-by-signature.
clones = select_candidates(clones, top_n=3)
shortlist = clones[clones["is_selected"]]

# 3) Optional per-method picks for the selection-route heatmap.
clone_method = build_clone_method_long(build_clone_sample_long(adata))
clone_method = clone_method.merge(sig_scores, on="CDR3_pair")
per_method_picks = compute_signature_picks_per_method(
    clone_method, top_n=1, pool_clones=set(tier3_plus_clones),
)
```

## API

::: tcrsift.candidate
    options:
      members:
        - select_candidates
        - compute_signature_picks_per_method
        - signature_picks_clone_to_methods
