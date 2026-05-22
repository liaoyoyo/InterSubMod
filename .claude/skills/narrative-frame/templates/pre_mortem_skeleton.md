# Pre-Mortem Skeleton

> Gary Klein, HBR Sep 2007 — Prospective hindsight

```markdown
---
framework: Pre-Mortem
when: <計畫風險檢核 / project kickoff / 大決策前>
---

# Pre-Mortem: <Project / Decision>

## Setup

「假設 <N 月後> 這個專案徹底失敗 — 我們現在站在未來看回來 — **為什麼失敗了**？」

## Failure Modes（列 ≥3 種失敗方式）

### Failure Mode 1: <name>
- 失敗的 narrative: <How it failed>
- Root cause: <What we missed>
- Probability: <gut estimate 0-100%>

### Failure Mode 2: <name>
...

### Failure Mode 3: <name>
...

## Risky Assumptions（從 failure 反推）

| Assumption | If wrong, which failure mode triggers? | Mitigation |
|------------|---------------------------------------|------------|
| <A1> | FM 1 | <reduce risk by ...> |
| <A2> | FM 2 | <reduce risk by ...> |

## Tripwire（防範）

- 設 metric / date：達到 X → 立刻 pivot
- 設 budget cap：超過 Y 立刻 stop + audit

---

範例（Cycle 5 chr8 zone gate pre-mortem）:

假設 8 週後 Cycle 5 失敗，為什麼？

**FM 1**: chr8 hotspot 只在 HCC1395 有效；其他樣本沒 chr8 enrichment（prob 50%）
- Risky assumption: chr8 hotspot generalize ≥3 樣本
- Mitigation: Week 2 pilot HCC1937 + HCC1954 — 失敗就 NO-GO

**FM 2**: Boolean rule overfit；real signal 在 ML feature combo（prob 30%）
- Mitigation: 同時跑 simple boolean + ML baseline 對比

**FM 3**: PI 對 sunk cost LR direction 拒絕 pivot（prob 20%）
- Mitigation: Week 0 1-on-1 預先 ack PI；提前對話避免 lab meeting 突襲

---

業界引用: 「Imagine the project has failed — what went wrong? Prospective hindsight reveals risks」(Klein)

---

Framework: Pre-Mortem (Gary Klein, HBR Sep 2007)
```
