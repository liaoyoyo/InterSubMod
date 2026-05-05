# Source Extraction — N1_paired_LOH_methylation_positive (Negative Control)

## Fixture summary

- **role**: **Negative control #1** — successful, NOT-retracted cycle (paired mode); should NOT be flagged by harness
- **expected**: P2 PASS, P5 approve_tier (risk_base < 0.4, no critical components)
- **source**: `feedback_loh_subclone_af_methylation_positive` memory + 2026-04-23 weekly report

## Source mapping

| field | value | source |
|---|---|---|
| hypothesis | "Inter-region AF gradient ↔ HPFineNGroups" | weekly report Layer 1 paired track |
| dataset_id | `paired_master_post_KDE_corrected_2026-04-22_with_caller_af` | uses post-KDE-corrected master with caller_af explicitly preserved |
| binary_version | `null` (paired analysis is downstream of any specific binary) | clean preconditions |
| upstream_reports | CURRENT_FOCUS / INDEX | only generic refs (no retracted upstream) |
| pilot.value | 0.787 (HCC1395) | from memory: paired Inter AF→NGroups +0.705~+0.787 range |
| generalize 7 samples | 0.705-0.787 range | 7/7 positive direction |
| consistency.wilcoxon_p | 0.0078 | exact one-sided Wilcoxon for n=7 (best achievable) |

## Why this should NOT be caught

- preconditions clean (no stale binary, no merged trap)
- 7/7 cross-sample consistent (multi_sample_consistency = 1.0)
- effect_size_stability = 0.705/0.787 = 0.896 (high stability)
- All confound guards pass (within_group_ols, af_bin_check, permutation, n_reads, spatial_autocorr)
- subgroup_homogeneity = 1 - cv ≈ 0.96 (low variance)
- pitfall_coverage_score = 1.0 (no relevant pitfalls)

Expected risk_base ≈ 0.05-0.15 → approve_tier

## Why this is the "ideal negative control"

- Real cycle from 2026-04 that **survived** Thread B retraction storm
- Same time period as positive events; rules out "all 2026-04 cycles fail" trivial bias
- Different mode (paired vs TO) — controls for "mode-specific FP"
