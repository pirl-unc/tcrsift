# In-Silico Filter Layer

`tcrsift.insilico_filter` applies composable percentile-rank predicates (such
as low `Ppost` in [sequence probability](seqprob.md) or high/low expression
signatures in [signatures](signatures.md)) within experimental sorting or
enrichment groups:

```
CTYneg
CTYneg + low Ppost
CTYneg + low Ppost + low naive_stem
CTYneg + low Ppost + low naive_stem + high antigen_response
```

Applying a predicate list creates a filtered counterpart subgroup (default
suffix `+IS`, for *in-silico*) alongside each original assay group (for
example, `CTYneg` and `CTYneg+IS`), so cross-method Jaccard overlap can be
compared before and after in-silico filtering.

## Ranking conventions

- Every predicate converts its feature to a **lower-is-better percentile rank**
  in `[0, 1]` (`0` = best). Setting `direction="high"` ranks values in
  descending order so high values rank near `0`.
- Equal values share a tied rank (`average` or `min`) rather than being
  dropped, preserving `Ppost == 0` as a valid minimum.
- A clone passes a predicate when its percentile rank is `<= percentile`; `NaN`
  values always fail.
- When `group_col` is supplied, percentile ranks are computed independently
  within each assay group.
- `average_percentile_rank` averages the per-predicate percentile ranks to
  provide a composite ranking score.

## Configuration

```yaml
insilico_filter:
  label: IS
  predicates:
    - {score: log10_ppost,                direction: low,  percentile: 0.5}
    - {score: signature_naive_stem,       direction: low,  percentile: 0.5}
    - {score: signature_antigen_response, direction: high, percentile: 0.5}
```

## Usage

```python
from tcrsift.insilico_filter import (
    apply_insilico_filter,
    insilico_overlap_long,
    predicates_from_config,
)

preds = predicates_from_config(config.insilico_filter)

# Add pass boolean, mean percentile rank, and composite subgroup label:
clones = apply_insilico_filter(clones, preds, group_col="method", label="IS")

# Compute pairwise overlap across base groups and their +IS counterparts:
mats = insilico_overlap_long(clones, preds, group_col="method")
jaccard = mats["jaccard"]
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
