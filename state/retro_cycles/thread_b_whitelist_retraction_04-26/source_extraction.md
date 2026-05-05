# Source Extraction — thread_b_whitelist_retraction_04-26

## Fixture summary

- **expected intercept**: P5 EVALUATE pending_review (multi_sample_consistency=1/6=0.167 → critical override; risk > 0.7)

## Source mapping

| field | value | source |
|---|---|---|
| binary_version | `8d0a0c8...` (KDE-corrected, post-fix) | trigger §2.2 method; deliberately fresh (this event has no stale-binary issue) |
| dataset_id | `X6_caller_AF_post_KDE_corrected_2026-04-24` | trigger §2.3 X6 data; "post_KDE" annotation makes dataset clean (P2 should PASS) |
| pilot.value | 0.583 (HCC1395 post-KDE) | trigger §2.4 |
| generalize.samples × 6 | per-sample S3 TP post-KDE | trigger §2.2 X6 §2.2 table |
| consistency.n_samples_passed | 1/6 | trigger §X6 §2.2 (only H2009 due to baseline saturation) |
| wilcoxon_p | 1.0 | trigger §X6 §2.2 「W=0 p=1」 |

## Why P5 (not P2)

- preconditions are clean post-KDE-fix (no stale binary, no merged trap)
- the issue is **cross-sample inconsistency** detected only by multi-sample-consistency module
- therefore P2 should PASS; P5 EVALUATE catches via:
  - multi_sample_consistency = 1/6 ≈ 0.167 (**critical** <0.2)
  - effect_size_stability = 0.129/0.903 ≈ 0.143 (**critical** <0.2)
  - 2 critical components → forced pending_review override

## Expected harness behavior

- **/check-staleness**: PASS
- **/run-evaluator**:
  - 2 critical components → `pending_review`
  - retraction_risk > 0.7 (composite + override stack)
  - tier_recommendation: keeps current (NEGATIVE/RETRACTED was set by user); evaluator confirms downgrade

## Fidelity note

- **fidelity: very high** — 6 樣本 cross-sample 表 + Wilcoxon 統計皆有 §X6 報告原始值
- this is the most "ideal" P5 catch case for Drill 1
