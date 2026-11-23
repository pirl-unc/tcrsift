# Filtering API

Assign clonotypes to confidence tiers (`tier1`–`tier5`) using fixed abundance
thresholds, logistic FDR estimation, or paired peptide-stimulation response.

::: tcrsift.filter
    options:
      members:
        - filter_clonotypes
        - filter_clonotypes_threshold
        - filter_clonotypes_logistic
        - filter_clonotypes_stimulation
        - assign_tiers_threshold
        - split_by_tier
        - clone_clears_tier
        - strictest_tier_met
        - get_filter_summary
        - score_stimulation_response
        - DEFAULT_THRESHOLD_TIERS
        - DEFAULT_FDR_TIERS
