# AI-Session-Companion Skeleton (= 既有 report)

> Timeline + Key decisions + Provenance — AI 對話紀錄

```markdown
---
framework: AI-Session-Companion
when: <AI 對話過程決策紀錄 / 非實驗 / 非工程結論>
對應既有 skill: /report (thin wrapper)
---

# AI Session Report: <Title>

- **Session ID**: <UUID>
- **Date**: YYYY-MM-DD HH:MM-HH:MM
- **Participant**: <AI model + user>
- **Topic**: <一句主題>

## TL;DR

<3 行 — 用戶問題 + AI 主要協助 + 結論>

## Timeline

| Time | Event | Outcome |
|------|-------|---------|
| 0:00 | <user prompt 1> | <AI response 摘要> |
| 0:15 | <decision moment> | <chose X over Y> |
| 0:30 | <pivot> | <reason> |
| 1:00 | <conclusion> | <verdict> |

## Key Decisions

### Decision 1: <name>
- Context: ...
- Considered: A / B / C
- Chosen: X
- Rationale: ...

### Decision 2: <name>
...

## Code / Doc Modifications

| File | Change | Why |
|------|--------|-----|
| `InterSubMod/.../X.cpp` | function Y refactor | ADR-42 |
| `InterSubMod/docs/.../Z.md` | section update | cycle 5 transition |

## Follow-up Actions

- [ ] <action 1 + owner + deadline>
- [ ] <action 2>

## Provenance Footer

- **Commit hash**: <git rev-parse HEAD>
- **AI model**: claude-opus-4-7
- **Skills invoked**: <list>
- **Session length**: <Y> min
- **Token usage**: <if relevant>

---

Framework: AI-Session-Companion (InterSubMod report)
對應既有 skill: /report (thin wrapper 預設套此 framework)
```
