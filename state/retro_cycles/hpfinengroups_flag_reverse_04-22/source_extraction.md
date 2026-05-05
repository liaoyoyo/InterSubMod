# Source Extraction — hpfinengroups_flag_reverse_04-22

> Forensic notes from `InterSubMod/docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md` and Agent 2 hard-event metric reproduction.

## Fixture summary

- **expected intercept**: P5 EVALUATE → downgrade or pending_review (effect_size_stability: HCC1954 0.707 vs H2009 0.957 → ratio 0.74; including COLO829 0.235/0.957 = 0.245 critical; pitfall_coverage low because n_reads_confound=fail)

## Source mapping

| field | value | source |
|---|---|---|
| pilot.metric_results.value | 0.707 (HCC1954) | trigger §3.1 per-sample TP rate; Agent 2 hard-fidelity reproduction |
| pilot.confound_guard.af_bin_check | fail | §3.2 Q2: AF<0.2 TP=0.937 vs AF≥0.4 TP=0.022 |
| pilot.confound_guard.n_reads_confound | fail | NOT checked in original analysis — this is the root cause for tier downgrade |
| generalize.samples × 7 | per-sample TP rate | Agent 2 §3.3 Table 3 (new filter AF<0.4+NR≥80) |
| consistency.n_samples_passed | 5/7 | 2 below threshold: HCC1954 (0.707) + COLO829 (0.235 out-of-scope) |

## Why P5 (not P2)

- Plan preconditions are clean (no stale binary, no merged dataset trap)
- Hypothesis itself is statistically valid; the issue is **incomplete confound check** + **cross-sample heterogeneity hidden by mean**
- Therefore P2 PRECHECK should PASS; the catch should happen at P5 EVALUATE via:
  1. effect_size_stability: min/max ratio computation
  2. pitfall_coverage_score: n_reads_confound was relevant (P-06) but unchecked

## Expected harness behavior

- **/check-staleness**: PASS (preconditions clean)
- **/run-evaluator**:
  - multi_sample_consistency = 5/7 = 0.714 (OK)
  - effect_size_stability = min(|values|)/max ≈ 0.235/0.965 = 0.244 (**critical** <0.2 if include COLO829; or 0.707/0.965=0.732 if exclude COLO829)
  - precondition_freshness = 1.0 (PASS)
  - subgroup_homogeneity = 1 - cv ≈ 0.4 (low, due to wide spread)
  - pitfall_coverage_score = 0 if n_reads_confound=fail triggers P-06 (low, critical)
  - **Expected verdict**: pending_review (1 critical override) OR downgrade_tier
  - **freeze prediction (§4.5.3-A)**: P5 downgrade or pending; **risk > 0.4**

## Fidelity note

- **fidelity: high** for cross-sample table (Agent 2 reproduction from §3.3 Table 3)
- **fidelity: medium** for hard-fidelity sample-level confound guard (we marked af_bin_check=fail and n_reads_confound=fail based on event description; original analysis did not run these checks systematically — this is the fidelity gap)
