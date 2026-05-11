# Infra-Ops Playbook

## When to invoke

Per SKILL.md USE WHEN. Common triggers:

- User says "我要跑全量 pipeline"
- `disk_guard.sh` alert fired (`.claude/state/disk_guard.alert` exists)
- Bash exit-1 with no stdout (silent failure → likely disk full)
- "為什麼 build 一直壞"

## 8-step workflow

1. **Preflight** → run `tools/tmp_check.sh` + `tools/disk_audit.sh` (both read-only)
2. **Identify** → which of 6 scenarios applies? (see `references/known_pitfalls.md`)
3. **Dry-run** → for cleanup: list what would be removed, no action
4. **Confirm** → `AskUserQuestion` 3 options (confirm / archive instead / cancel)
5. **Execute** → only if confirmed; use atomic ops (`mv` first then `rm`)
6. **Audit log** → append `.claude/state/infra_operations.jsonl` with action + before/after metric
7. **Verify** → re-run preflight, confirm metric improved (e.g., `df -h` shows recovered space)
8. **Done** → return summary (X GB freed, Y dirs touched, audit log path)

## Anti-patterns (DO NOT)

- DO NOT `rm -rf` without dry-run + `AskUserQuestion` confirm
- DO NOT delete from `/tmp` blindly (might kill running processes)
- DO NOT proceed if disk_guard alert is active without addressing root cause
- DO NOT skip the audit log (forensic trail mandatory for destructive ops)

## Audit log format (jsonl)

Each line is one JSON object:

```json
{"ts": "2026-05-11T12:34:56Z", "scenario": 3, "action": "archive", "target": "output/synthesis/round_03", "size_gb_before": 87, "size_gb_after": 0, "user_confirmed": true, "session_id": "..."}
```

Append-only. Never edit existing lines. File lives at `.claude/state/infra_operations.jsonl` (gitignored).
