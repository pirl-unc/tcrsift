# Filtering API

Assign clonotypes to confidence tiers (`tier1`–`tier5`) using fixed abundance
thresholds or logistic FDR estimation.

::: tcrsift.filter
    options:
      members:
        - filter_clonotypes
        - filter_clonotypes_threshold
        - filter_clonotypes_logistic
        - assign_tiers_threshold
        - split_by_tier
        - clone_clears_tier
        - strictest_tier_met
        - get_filter_summary
        - DEFAULT_THRESHOLD_TIERS
        - DEFAULT_FDR_TIERS
