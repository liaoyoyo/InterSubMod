# `state/` — Resilient Waterfall Cycle State

> **Purpose**: Single Source of Truth for active research cycles in the 7-phase Resilient Waterfall harness.
>
> **Plan reference**: `~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md`
> **Created**: 2026-05-03 (Path A · Week 1 · Day 1-2)

---

## Directory Layout

```
state/
├── README.md                         # This file
├── active.json                       # Index of active cycles (≤5)
├── schemas/                          # JSON Schema definitions (draft 2020-12)
│   ├── active.schema.json
│   ├── state.schema.json             # Per-cycle main state
│   ├── plan.schema.json              # P1 PLAN artifact
│   ├── precheck.schema.json          # P2 PRECHECK artifact
│   ├── pilot.schema.json             # P3 PILOT artifact
│   ├── generalize.schema.json        # P4 GENERALIZE artifact
│   ├── evaluation.schema.json        # P5 EVALUATE artifact
│   └── reflection.schema.json        # reflection.log entry
├── cycles/                           # Active cycles
│   └── {cycle_id}/                   # cycle_id = YYYYMMDD-HHMM-{slug}
│       ├── state.json                # Main state (matches state.schema.json)
│       ├── plan.json
│       ├── precheck.json
│       ├── pilot.json
│       ├── generalize.json
│       ├── evaluation.json
│       └── reflection.log            # JSON Lines, each line per reflection.schema.json
├── cycles_archived/                  # Cycles after P6 COMMIT
└── invalidation/
    ├── stale_marks.jsonl             # Reports / cycles that need rerun
    └── binary_versions.jsonl         # C++ binary commit hash + timestamp index
```

## Conventions

- `cycle_id` format: `YYYYMMDD-HHMM-{slug}` (e.g. `20260503-2330-loh-kde-quantify`)
- All timestamps in ISO 8601 UTC (`2026-05-03T23:30:00Z`)
- All JSON files validated against their schema in `schemas/`
- `*.jsonl` = JSON Lines (one JSON object per line, append-only)
- Schemas use `$schema: https://json-schema.org/draft/2020-12/schema`

## Source-of-Truth Rule

When state and `evidence_ledger.jsonl` disagree:

1. **`evidence_ledger.jsonl` is the historical SoT** (append-only audit log)
2. **`state/cycles/{id}/state.json` is the current snapshot** (mutable)
3. On crash recovery, rebuild state from ledger by replaying entries

## File Lifecycle

| Phase | Action |
|---|---|
| `/cycle-init` | Create `cycles/{id}/state.json`, append entry to `active.json` |
| Phase advance | Update `state.json.phase` + write phase artifact, append `reflection.log` |
| `/check-staleness` | Read from `invalidation/*.jsonl`, write `precheck.json` |
| `/run-evaluator` | Read all artifacts, write `evaluation.json` |
| `/cycle-state` | Read-only dashboard over `active.json` + each `state.json` |
| P6 COMMIT | Move `cycles/{id}/` → `cycles_archived/{id}/`, remove from `active.json` |

## Hooks That Touch This Directory

| Hook | Trigger | Action |
|---|---|---|
| `post_cpp_commit_invalidate.sh` | PostToolUse on git commit of `src/*.cpp` | Append to `invalidation/binary_versions.jsonl` + scan reports → `stale_marks.jsonl` |
| `pre_tier_upgrade_check.sh` | PreToolUse on edit of INDEX.md / state.json that changes tier to ⭐4-5 | Block if `evaluation.json` missing |
| (SessionStart) | Claude Code session start | Read `active.json` + `stale_marks.jsonl` summary |

## Git Tracking

- ✅ **Tracked**: `schemas/`, `README.md`, `cycles/{id}/*` (final artifacts), `cycles_archived/`, `invalidation/binary_versions.jsonl`
- ❓ **Optional gitignore**: `active.json` (high-frequency churn) — decision deferred to Day 3+
- ❌ **Never tracked**: nothing yet

## Schema Versioning

All schemas use `schema_version` field. Current: **`1.0`**.

**Rule**: never edit a schema in place once a cycle is using it. To extend:
1. Bump `schema_version` to `1.1`
2. Make new fields optional with sensible defaults
3. If breaking, write a migration script `state/schemas/migrate_v1_0_to_v1_1.py`
