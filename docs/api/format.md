# Display Formatting API

Helpers for rendering sample-sheet enrichment-method strings in
figures. The raw method names tcrsift accepts (`AIMpos`,
`CTYneg_tetpos`) are pipeline-friendly but ugly in plot axes.
`tcrsift.format` translates them to display-friendly forms:

| Raw | Pretty |
| --- | --- |
| `AIMpos` | `AIM⁺` |
| `CTYneg` | `CTY⁻` |
| `AIMpos_CTYneg` | `AIM⁺CTY⁻` |
| `CTYneg_tetpos` | `tet⁺CTY⁻` (CTY reorder) |
| `AIMpos_CTYneg-2` | `AIM⁺CTY⁻` (donor suffix dropped) |

Applied automatically by tcrsift's built-in plots (method × method
overlap heatmap, per-method recovery panels, per-sample scatter
titles). Call directly when adding labels to custom plots.

## Usage

```python
from tcrsift.format import pretty_method, pretty_samples

pretty_method("AIMpos_CTYneg")   # "AIM⁺CTY⁻"
pretty_method("CTYneg_tetpos")   # "tet⁺CTY⁻"

# Map over a list / Series of names.
pretty_samples(["AIMpos-1", "tetpos-1"])  # ["AIM⁺", "tet⁺"]
```

::: tcrsift.format
    options:
      members:
        - pretty_method
        - pretty_methods
        - pretty_sample
        - pretty_samples
