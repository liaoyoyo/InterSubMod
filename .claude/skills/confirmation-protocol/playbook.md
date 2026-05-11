# Confirmation Protocol Playbook

## Decision flow

```
1. Is this an irreversible action (Hard Gate)?
   ├─ YES → Look up references/askuserquestion_patterns.md §1-§5
   │        → Call AskUserQuestion with structured options
   │        → Branch on user's choice (see § 後續邏輯)
   └─ NO → continue

2. Mode check:
   ├─ 互動模式 → Use prose pause for Gate/Review (existing behavior)
   └─ 全自動模式 → Skip Gate/Review, only Hard Gate stops

3. Action:
   - 🔴 Hard Gate → AskUserQuestion (always)
   - 🟠 Gate (互動) → Prose pause "請確認...?"
   - 🟡 Review (互動) → Brief node summary + 「繼續嗎？」
   - 🟢 FYI → One-line decision trace, no pause
```

## Quick reference

| Decision | Tool |
|----------|------|
| Hard Gate (5 scenarios) | AskUserQuestion (structured) |
| Gate / Review pause | Prose question |
| FYI / 一行告知 | Plain text output |
| Plan approval | ExitPlanMode (not AskUserQuestion) |

## Anti-patterns

- DO NOT use AskUserQuestion for non-Hard-Gate confirmations (over-formalizes simple flows)
- DO NOT use prose 「請確認」 for the 5 Hard Gates (regression to old behavior)
- DO NOT reference plan content in AskUserQuestion question text (user can't see plan there — use ExitPlanMode for plan)
