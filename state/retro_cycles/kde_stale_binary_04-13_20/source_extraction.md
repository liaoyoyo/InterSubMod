# Source Extraction — kde_stale_binary_04-13_20

> Forensic notes from `InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md` and CovM baseline report.

## Fixture summary

- **event_id**: `kde_stale_binary_04-13_20`
- **fixture mode**: full 5-artifact set (state + plan + pilot + generalize + this source_extraction)
- **expected intercept**: P2 BLOCKED via binary stale (HEAD has KDE fix, plan-stated SHA predates fix → stale_distance ≥ 1)

## Source mapping

| field | value | source |
|---|---|---|
| binary_version | `4dc2d732c7230d71f7147692de3fd1dafcf88a7a` | git log: pre-KDE-fix commit (HEAD when fixture written = bc454ab; KDE fix at 8d0a0c8 in plan reference) |
| dataset_id | `master_2026-03-30_stale_binary_pre_KDE_fix` | Agent 1 forensic §2.1; deliberately encodes "stale" keyword for debugging clarity (P2 will detect via binary version mismatch) |
| pilot.metric_results.value | -0.019 (−1.9% bias) | trigger §1.2 / §2.2 row |
| pilot.metric_results.n_samples | 28495 regions | trigger §1.1 |
| pilot.confound_guard | uniform rescaling 0.880 → 1.245 (×1.415) | trigger §4.1 |
| generalize.samples × 7 | bias values | Agent 1 forensic §5.0 (post-fix per-sample bias) |
| consistency.direction | "7/7 bias < 3%" | trigger §5.0 verdict; pre-fix range −17.6% to +158.6% (COLO829 worst) |

## Why P2-only intercept (not P5)

- Pilot/generalize artifacts post-KDE-fix show **good results** (7/7 bias <3%)
- The retraction is **not** about the conclusion — the conclusion (KDE auto_kde works) is correct
- The retraction is about **all upstream cycles built on stale 75.0 binary** that should have been blocked
- Therefore P2 PRECHECK should detect: plan-stated binary_version (pre-fix) vs HEAD (post-fix) → STALE → all derived analyses should have been blocked or re-run

## Expected harness behavior (freeze table §4.5.3-A)

- **/check-staleness**: BLOCKED (`binary` status=stale, stale_distance >= 1, expected != HEAD)
- **/run-evaluator**: would run (full 5 artifacts present); expected `precondition_freshness=0.3` (BLOCKED) → composite drops; risk_base likely > 0.4 (downgrade tier)

## Fidelity note

- **fidelity: high** — 7/7 sample data has explicit Agent 1 forensic
- **caveat**: SHA `4dc2d732...` is a specific git commit; if HEAD moves further, stale_distance grows but P2 verdict still BLOCKED (correct behavior)
