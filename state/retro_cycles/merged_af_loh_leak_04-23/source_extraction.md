# Source Extraction — merged_af_loh_leak_04-23

## Fixture summary

- **expected intercept**: P2 BLOCKED via merged dataset trap (dataset_id contains "merged" + missing "caller_af" annotation → `_probe_vcf_header` heuristic + dataset name probe both fire)
- **secondary**: even if P2 PASSed (override), P5 EVALUATE would catch via multi_sample_consistency=1/6=0.167 (critical)

## Source mapping

| field | value | source |
|---|---|---|
| dataset_id | `merged_master_2026-04-15_x6_caller_AF_S3S5` | trigger §2.3 X6 caller_af merge data; deliberately encodes "merged" trigger keyword for P-03 P-04 sweep |
| pilot.value | 0.955 (n=380) | trigger §2.4 「原 v2」 row |
| pilot.confound_guard | all skipped | original pilot did NOT cross-check KDE-corrected; this is the root cause |
| generalize.samples × 6 | per-sample S3 TP rate post-KDE | trigger §2.4 X6 cross-sample table |
| consistency.wilcoxon_p | 1.0 | trigger §2.2 X6 §2.2 "W=0 p=1" (one-sided S3 < baseline) |

## Why P2 should fire (primary intercept)

dataset_id contains "merged" → `check_dataset()` heuristic matches → schema_violations populated → P2 BLOCKED. Plan's preconditions binary_version is also stale (4dc2d73 pre-KDE-fix), so binary check should also flag stale.

## Expected harness behavior

- **/check-staleness**: BLOCKED (binary stale + dataset "merged" trap)
- **/run-evaluator** (if forced past P2): all components low → multi_sample_consistency=0.167 critical → pending_review

## Fidelity note

- **fidelity: high** — Agent 1 §2.2 提供完整 6 樣本 X6 cross-sample 表
- **caveat**: fidelity depends on user accepting Drill 1 as P2-primary intercept (which is consistent with §4.5.3-A freeze)
