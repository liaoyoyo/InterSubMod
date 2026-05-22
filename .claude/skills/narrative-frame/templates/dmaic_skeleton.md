# DMAIC Skeleton

> Six Sigma — Bill Smith / Motorola → Jack Welch GE

```markdown
---
framework: DMAIC
when: <流程改善 / quality control / operational excellence>
---

# <Process Improvement Title>

## 1. Define（定義）

- **Problem statement**: <一句問題>
- **Goal**: <量化目標 — X → Y>
- **Scope**: <process 範圍>
- **Stakeholders**: <涉事者>

## 2. Measure（測量）

- **Current state metric**: <baseline 數據>
- **Data collection plan**: <怎麼測 + 多久>
- **Capability index**: <Cp / Cpk if applicable>

## 3. Analyze（分析）

- **Root cause analysis**: 5 Whys / Fishbone
- **Statistical analysis**: <hypothesis test / regression / control chart>
- **Key drivers**: <影響最大的 1-3 個因子>

## 4. Improve（改善）

- **Proposed solution**: <方案>
- **Pilot test**: <小規模驗證>
- **Implementation plan**: <推全 plan>

## 5. Control（控制）

- **Monitoring metric**: <持續追蹤 metric>
- **Control chart**: <SPC / dashboard>
- **Handoff**: <移交 operations + documentation>

---

範例（pipeline reliability 改善）:

1. Define: paired pipeline run 失敗率 12%；目標 < 2%
2. Measure: 過去 30 天 200 runs 中 24 失敗；fail mode classification
3. Analyze: 80% failures = /tmp 寫滿（infrastructure root cause）
4. Improve: 強制 export TMPDIR=/big7_disk; pipeline preflight check
5. Control: disk_guard hook + 月度 audit

---

Framework: DMAIC (Six Sigma; George《Lean Six Sigma》2002)
```
