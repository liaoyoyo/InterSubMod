---
name: infra-ops
description: Runs disk and infrastructure preflight when user starts long pipeline OR after disk_guard alert OR diagnosing OOM/disk-full incidents. USE WHEN pipeline preflight, post-OOM, server migration, /tmp suspicions, output dir cleanup proposals. SKIP WHEN pure code edits, doc writing, single-file analysis.
type: skill
allowed-tools: Read, Bash, Grep, AskUserQuestion
---

# Infra-Ops Skill

Defensive infrastructure operations skill — pipeline preflight, disk-full / OOM diagnosis, and state-machine-guarded cleanup of large output dirs. Built around the 2026-05-08 /tmp 800 GB disaster (see `references/tmp_disaster_2026_05_08.md`).

## Phase & Chain Position

Defensive layer; fires **before** long-running pipelines, **after** `disk_guard.sh` alerts, or **during** OOM / disk-full diagnosis. Sits parallel to `confirmation-protocol` Hard Gate for any filesystem-mutating action.

## Dependencies

- **Uses**: existing `scripts/infra/disk_guard.sh` (monitor) + `scripts/hooks/pipeline_block_check.sh` (pre-pipeline gate)
- **Used by**: any long-pipeline workflow (`run_vcf_all_snv.sh`, `run_batch_vcf_analysis.sh`), `confirmation-protocol` Hard Gate for filesystem changes
- **Reads**: `df -h`, `mount`, `.claude/state/disk_guard.alert` (if exists)
- **Writes**: `.claude/state/infra_operations.jsonl` (append-only audit trail; NOT git-tracked)

## Failure Mode & Diagnostics

| Symptom | Likely cause | Quick check |
|---------|-------------|-------------|
| Bash commands silently fail with exit 1 | /tmp full | `df -h /` → if 100% confirmed |
| Pipeline OOM-killed mid-run | TMPDIR not set, output to root volume | `mount \| grep tmp` + check pipeline log for /tmp paths |
| `make` fails after long absence | CMakeCache stale | `ls -lt build/CMakeCache.txt` vs source mtimes |

## 6 Scenario Index → references/known_pitfalls.md

1. /tmp cleanup (root >85% AND /tmp <1GB free)
2. CMake cache reset (CMakeCache stale)
3. Output archive (`output/synthesis/` >500GB)
4. Conda env reset (broken activation)
5. Docker dangling images (>50GB)
6. Log rotation (`/var/log` >50GB)

## State Machine for Destructive Ops

Never auto-delete. Always:

1. Dry-run list (read-only) → `tools/disk_audit.sh`
2. `AskUserQuestion(["confirm cleanup", "archive instead", "cancel"])`
3. Execute only after confirm
4. Append to `.claude/state/infra_operations.jsonl`

## Quick Reference

- Pipeline preflight → `tools/tmp_check.sh`
- Project disk audit → `tools/disk_audit.sh`
- Step-by-step recovery → `playbook.md` (8-step workflow)
- Case study → `references/tmp_disaster_2026_05_08.md`
- Per-scenario fix → `references/known_pitfalls.md`
