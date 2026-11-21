# Display Formatting API

`tcrsift.format` converts machine-friendly sample and `enrichment_method`
identifiers into compact plot labels with Unicode superscripts:

| Raw string | Formatted label | Rule applied |
| --- | --- | --- |
| `AIMpos` | `AIM⁺` | `pos` → `⁺` |
| `CTYneg` | `CTY⁻` | `neg` → `⁻` |
| `AIMpos_CTYneg` | `AIM⁺CTY⁻` | Underscore separator removed |
| `CTYneg_tetpos` | `tet⁺CTY⁻` | Baseline marker (`CTY`) moved last |
| `AIMpos_CTYneg-2` | `AIM⁺CTY⁻` | Trailing donor/replicate suffix stripped |

Built-in overlap heatmaps, method-recovery bar charts, and per-sample scatter
plots apply these helpers automatically; call them directly when labeling
custom figures.

## Usage

```python
from tcrsift.format import pretty_method, pretty_samples

pretty_method("AIMpos_CTYneg")            # "AIM⁺CTY⁻"
pretty_method("CTYneg_tetpos")            # "tet⁺CTY⁻"
pretty_samples(["AIMpos-1", "tetpos-1"])  # ["AIM⁺", "tet⁺"]
```

::: tcrsift.format
    options:
      members:
        - pretty_method
        - pretty_methods
        - pretty_sample
        - pretty_samples
