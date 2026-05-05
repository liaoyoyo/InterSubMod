# Source Extraction — N2_zone_aware_characterization (Negative Control)

## Fixture summary

- **role**: **Negative control #2** — characterization-only verdict, validated and retained at ⭐3
- **expected**: P2 PASS, P5 approve_tier
- **source**: `feedback_zone_aware_framework` memory + research_landscape/08_Zone_Aware.md

## Source mapping

| field | value | source |
|---|---|---|
| hypothesis | Zone-aware framework characterization | research_landscape/08_Zone_Aware.md |
| dataset_id | `zone_aware_master_post_KDE_corrected_caller_af_separate` | clean preconditions |
| binary_version | `null` | downstream characterization, not binary-dependent |
| pilot.value | 0.142 (HCC1395 zone TP rate diff) | typical effect size for characterization |
| generalize 7 samples | 0.117-0.151 range (CV ~10%) | tight cross-sample consistency |

## Why this should NOT be caught

- preconditions clean (post-KDE, caller_af separate)
- 7/7 cross-sample consistent
- effect_size_stability = 0.117/0.151 = 0.775 (high)
- All confound guards pass
- characterization-only verdict — does NOT claim QS adjustment or tier upgrade beyond ⭐3
- tier already at ⭐3 (current); not seeking upgrade to ⭐4-5 (so no /run-evaluator strict gate needed)

Expected risk_base ≈ 0.10-0.20 → approve_tier

## Difference from N1

- **N1**: paired-mode subclone marker (POSITIVE conclusion, biology-informed)
- **N2**: characterization-only (no biology-informed claim, just observation)
- Together they cover two distinct "good cycle" archetypes — biology-informed + observation-only
