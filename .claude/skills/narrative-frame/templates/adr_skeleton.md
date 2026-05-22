# ADR Skeleton — Architecture Decision Record

> Michael Nygard 2011; AWS Prescriptive Guidance

```markdown
---
framework: ADR (Architecture Decision Record)
when: <工程架構決策紀錄 / 技術選型 / 長期可追溯>
---

# ADR-<NN>: <Decision Title>

- **Date**: YYYY-MM-DD
- **Status**: Proposed / Accepted / Deprecated / Superseded by ADR-XX

## Context

<Why is this decision needed? 背景 + 約束 + stakeholder + 目標>

## Decision Drivers

- Driver 1: <e.g., performance, cost, maintainability>
- Driver 2: ...
- Driver 3: ...

## Considered Options

### Option A: <name>
- Pros: ...
- Cons: ...

### Option B: <name>
- Pros: ...
- Cons: ...

### Option C: <name>
- Pros: ...
- Cons: ...

## Decision

**Chosen**: Option <X>

**Rationale**: <為什麼選這個；對應哪些 Decision Drivers>

## Consequences

- ✅ Positive: <e.g., what we gain>
- ⚠ Negative: <what we lose; trade-off>
- 🔄 Follow-up: <what to revisit>

## Links

- Related ADR: ADR-XX
- PR / commit: <hash>
- Discussion: <issue link>

---

範例:

# ADR-042: 採用 /narrative-frame 取代 7 個固定範本 skill

- Date: 2026-05-22
- Status: Accepted

## Context
既有 7 報告類 skill 用固定範本（13/17/5 段）；50+ 業界 framework 散在各處無 SoT；用戶要求動態挑 framework 減少理解負擔。

## Drivers
- 1. 動態框架支援多場景
- 2. 不破壞既有 INDEX 引用
- 3. 對話層級啟用

## Options
- A: 全替換（刪 7 skill）— 高 migration cost
- B: thin wrapper（保 skill 名 + 內部 forward）— moderate cost
- C: 並存（新 skill + 舊 skill 共存）— 用戶混淆

## Decision
**Chosen**: B (thin wrapper)

## Consequences
- ✅ INDEX / 引用全保留
- ⚠ /narrative-frame 是 single point of failure；catalog drift 風險
- 🔄 V6 audit 每月跑 catalog drift check

---

Framework: ADR (Nygard 2011; AWS Prescriptive Guidance)
```
