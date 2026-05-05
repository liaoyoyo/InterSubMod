# Source Extraction — vcf_source_error_04-04

> Forensic notes mapping fixture values to original 2026-04-04 reports.
> Used by Day 6 Drill 1 to verify retro_cycle artifacts faithfully reflect the actual incident.

## Fixture summary

- **event_id**: `vcf_source_error_04-04`
- **fixture mode**: **P2-only** (no pilot.json / generalize.json — original event was diagnostic, not executed)
- **expected intercept**: P2 BLOCKED (dataset_id contains "pileup symlink" + caller mismatch keywords → triggers `_probe_vcf_header` heuristic)

## Source reports

| field | value | source |
|---|---|---|
| hypothesis | "pileup symlink → ClairS paired for TO analysis" | `InterSubMod/docs/experiments/in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md` §1.1-1.3 |
| dataset_id | `pileup_symlink_ClairS_paired_for_TO_analysis` | trigger §1.1 (symlink resolved to paired caller); deliberately encodes the mismatch keyword to trigger P-04 |
| upstream_reports | 3 affected | trigger §1.1 lists 4 affected: HPFineP_QS_research / HPFineP_QS_integrated / LOH_Strong_Weak / LOH_AlleleDelta (we use first 3) |
| samples | `["HCC1395"]` | trigger §3.1 — single-sample diagnostic |
| binary_version | `null` | event was pure diagnostic; no binary commit triggered |

## Why no pilot.json / generalize.json

Original event was a **planned correction action** (re-run ISM with correct VCF), but the corrective re-run was **not executed within this event**. The event ends at "diagnose + recommendation" stage. Therefore:

- `pilot.json`: N/A (no actual run)
- `generalize.json`: N/A (no cross-sample data — diagnostic was on single-sample HCC1395)
- Only `plan.json` carries the preconditions; expected to BLOCK at P2 PRECHECK before any pilot runs

## Expected harness behavior (freeze table §4.5.3-A)

- **/check-staleness** verdict: **BLOCKED**
  - Binary check: status=`skipped` (binary_version=null)
  - Dataset check: should detect `pileup_symlink` keyword and flag VCF source mismatch (P-04 pitfall)
  - Upstream reports: 3/3 fresh (we are not testing report freshness here)
- **/run-evaluator**: **NOT RUN** (no pilot/generalize — fixture is P2-only)

## Caveats / fidelity note

- fidelity: **medium-high** for plan.json (95% per Agent 1 forensic)
- fidelity: **N/A** for pilot/generalize (event by nature has no run)
- **Note**: this is the only event in 6 where evaluator does not run; comparison against negative controls in confusion matrix will be sensitivity-only for this case
