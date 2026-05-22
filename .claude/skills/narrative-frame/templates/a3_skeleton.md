# A3 Report Skeleton

> Toyota internal practice; John Shook《Managing to Learn》(2008) — 一張 A3 紙講完問題 + 提案

```markdown
---
framework: A3 Report
when: <Kaizen / bug fix postmortem / 工程改善敘事>
---

# A3 — <Title>

## 1. Background

<為什麼這個問題重要 — 戰略 context>

## 2. Current State

<量化現況 — 現在發生什麼 + 影響>

(source: ...)

## 3. Goal

<量化目標 — 我們想要達到什麼 metric>

## 4. Analysis（根因）

### 4.1 5 Whys

- Why 1: <表象問題> → because <Why 2>
- Why 2: ... → because <Why 3>
- Why 3: ... → because <Why 4>
- Why 4: ... → because <Why 5>
- **Root cause**: <Why 5 的答案>

### 4.2 Fishbone（6M）

- Man: <人因>
- Machine: <機器 / tool>
- Method: <方法>
- Material: <資料 / input>
- Measurement: <測量>
- Milieu: <環境>

## 5. Proposal

<候選方案 + 採用理由 (ADR style)>

- Option A: <pros / cons>
- Option B: <pros / cons>
- **Chosen**: A (因為 ...)

## 6. Plan

| Step | Who | When | Verify |
|------|-----|------|--------|
| 1 | ... | ... | ... |
| 2 | ... | ... | ... |

## 7. Follow-up

<驗證計畫 + 後續審查 cadence>

---

Framework: A3 Report (Toyota TPS; Shook《Managing to Learn》2008)
```
