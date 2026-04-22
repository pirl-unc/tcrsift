# Filtering Strategies

TCRsift provides tiered filtering to prioritize antigen-specific clones while minimizing false positives.

## Overview

Antigen-specific T cell cultures typically contain:

- **True positives**: Clones expanded due to antigen recognition
- **Bystander clones**: Expanded due to cytokines, not antigen
- **Contaminants**: Low-frequency clones from sample handling

TCRsift uses multiple criteria to separate these populations.

## Threshold Method (Default)

The threshold method applies configurable criteria at each tier:

| Tier | Min Cells | Min Frequency | Max Conditions | Confidence |
|------|-----------|---------------|----------------|------------|
| Tier 1 | 10 | 1% | 2 | Highest |
| Tier 2 | 5 | 0.5% | 3 | High |
| Tier 3 | 3 | 0.1% | 5 | Medium |
| Tier 4 | 2 | 0.05% | 10 | Low |
| Tier 5 | 2 | 0% | unlimited | Lowest |

### Criteria Explained

**Min Cells**: Clones with more cells are more likely to be true positives.

**Min Frequency**: Higher frequency suggests active expansion.

**Max Conditions**: Antigen-specific clones should stay concentrated within a small number of antigen conditions. TCRsift uses antigen/condition counts when they are available and falls back to sample count otherwise.

### Using Threshold Filtering

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method threshold
```

### Custom Tier Definitions

```python
from tcrsift import filter_clonotypes

custom_tiers = {
    "tier1": {"min_cells": 20, "min_frequency": 0.02, "max_conditions": 1},
    "tier2": {"min_cells": 10, "min_frequency": 0.01, "max_conditions": 2},
    "tier3": {"min_cells": 5, "min_frequency": 0.005, "max_conditions": 3},
}

filtered = filter_clonotypes(clonotypes, tier_definitions=custom_tiers)
```

## Logistic Method

The logistic method uses logistic regression to adaptively set thresholds based on your data.

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --method logistic \
    --fdr-tiers 0.0001,0.001,0.01,0.1,0.15
```

### How It Works

1. Fits a logistic model predicting "high quality" based on frequency
2. Uses the model to estimate probability for each clone
3. Assigns tiers based on FDR thresholds

### FDR Tiers

| Tier | FDR | Interpretation |
|------|-----|----------------|
| Tier 1 | 0.0001 | < 0.01% false discovery rate |
| Tier 2 | 0.001 | < 0.1% false discovery rate |
| Tier 3 | 0.01 | < 1% false discovery rate |
| Tier 4 | 0.1 | < 10% false discovery rate |
| Tier 5 | 0.15 | < 15% false discovery rate |

### When to Use Logistic

The logistic method is useful when:

- You have many samples with varied characteristics
- You want data-driven threshold selection
- You need FDR-based confidence estimates

!!! warning "Note"
    The logistic method can be sensitive to noisy data. If you see unexpected results, try the threshold method instead.

## T Cell Type Filtering

Filter by phenotype classification:

```bash
# CD8+ only (default for short peptides)
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type cd8

# CD4+ only
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type cd4

# Both CD4+ and CD8+
tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type both
```

TCRsift uses the consensus T cell type from phenotyping. A clone is classified based on the majority of its cells.

## Viral Exclusion

Exclude clones matching known viral epitopes:

```bash
tcrsift filter -i clonotypes.csv -o filtered/ --exclude-viral
```

This requires prior annotation with `tcrsift annotate`. Clones matching CMV, EBV, HIV, influenza, and other common viruses are excluded.

## Combining with TIL Data

For tumor studies, the strongest evidence for tumor-specificity is:

1. Clone expanded in antigen culture
2. Clone also present in tumor (TIL)

```bash
# First filter culture data
tcrsift filter -i culture_clonotypes.csv -o filtered/ --method threshold

# Then match against TIL
tcrsift match-til -i filtered/tier1.csv --til-csv til_clonotypes.csv -o matched.csv
```

Clones in both culture and TIL are the highest confidence candidates.

## Recommended Workflow

1. **Start with default threshold filtering**
   ```bash
   tcrsift filter -i clonotypes.csv -o filtered/ --tcell-type cd8
   ```

2. **Review tier distributions**
   - Many clones in tier 1? Your cultures worked well.
   - Mostly tier 4-5? Consider more stringent culture conditions.

3. **Annotate and exclude viral**
   ```bash
   tcrsift annotate -i filtered/tier1.csv -o annotated.csv --vdjdb /path/to/vdjdb
   ```

4. **If available, validate with TIL**
   ```bash
tcrsift match-til -i annotated.csv --til-csv til.csv -o final.csv
   ```

5. **Prioritize clones for validation**
   - Tier 1 + TIL match = highest priority
   - Tier 1-2 + no viral match = high priority
   - Tier 3-4 = consider if functional validation capacity allows
