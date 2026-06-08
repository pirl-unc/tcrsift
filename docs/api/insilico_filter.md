# In-silico filter layer

A **stackable filter layer** (#149) that composes *on top of* the
assay/enrichment groups (the sorts) and the selection rules. Each layer is
a per-clone percentile predicate on an existing feature — Ppost
([#143](seqprob.md)) or a GEX signature score
([signatures](signatures.md)) — so a group can be progressively refined:

```
CTYneg
CTYneg · Ppost^low
CTYneg · Ppost^low · GEXnaive^low
CTYneg · Ppost^low · GEXnaive^low · GEXresponse^high
```

Applying the full stack yields one **named composite** subgroup (default
label `IS`, for *in-silico*), so each assay group gains a filtered twin —
`CTYneg` and `CTYneg+IS`. N assay groups → 2N groups, and the twins feed
the cross-group Jaccard-overlap analysis.

## Scoring conventions

- Each predicate is a **lower-is-better percentile rank** in `[0, 1]`
  (0 = best). `direction="high"` flips it so high scores win
  (e.g. `GEXresponse^high` → `−antigen_response`).
- Tie handling keeps `Ppost == 0` a valid **best** value and groups all
  equal scores together (`rank(method="average"/"min")`) — never dropped.
- A clone **passes** a predicate when its percentile rank ≤ the predicate's
  `percentile`. A missing score never passes.
- With a group column, ranks are computed **within each assay group**, so
  the filter refines each group against its own distribution.
- The composite's own ranking metric is the **average of the per-dimension
  percentile ranks** → top-k or a percentile cut.

## Config

```yaml
insilico_filter:
  label: IS
  predicates:
    - {score: log10_ppost,             direction: low,  percentile: 0.5}
    - {score: signature_naive_stem,    direction: low,  percentile: 0.5}
    - {score: signature_antigen_response, direction: high, percentile: 0.5}
```

## Usage

```python
from tcrsift.insilico_filter import (
    predicates_from_config, apply_insilico_filter, insilico_overlap_long,
)

preds = predicates_from_config(config.insilico_filter)

# Annotate each clone with pass-flag, avg rank, and composite group name:
clones = apply_insilico_filter(clones, preds, group_col="method", label="IS")

# Jaccard overlap across base groups AND their +IS twins:
mats = insilico_overlap_long(clones, preds, group_col="method")
jaccard = mats["jaccard"]   # e.g. does AIMpos+IS overlap tetpos+IS more?
```

::: tcrsift.insilico_filter
    options:
      members:
        - FilterPredicate
        - predicates_from_config
        - percentile_rank
        - insilico_mask
        - average_percentile_rank
        - apply_insilico_filter
        - expand_insilico_twins
        - insilico_overlap_long
