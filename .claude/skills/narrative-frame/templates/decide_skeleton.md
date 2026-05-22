# DECIDE Skeleton

> Kristina Guo (2008) — 醫療緊急決策

```markdown
---
framework: DECIDE
when: <醫療緊急 / 急救 / 生死決策 / structured decision>
---

【D】Define problem: <清楚 problem statement>

【E】Establish criteria: <評估維度 — must / should / nice-to-have>

【C】Consider alternatives: <列 ≥3 個 options>

【I】Identify best alternative: <評分 + 選最佳>

【D】Develop and implement plan: <action + timeline + responsibility>

【E】Evaluate and monitor: <outcome metric + 回饋 loop>

---

範例（cycle 跑到一半要不要 abort）:

【D】Cycle 5 chr8 zone gate week 2 pilot HCC1937 = NEGATIVE — 要不要 abort 還是繼續？
【E】
- Must: 不超出 4 週 budget；不浪費 PI tolerance
- Should: 保留 chr8 hypothesis 可後續驗證
- Nice: 同 cycle 內 fall back to backup direction
【C】
- A: 立刻 abort + NO-GO（純損失）
- B: 多 1 樣本 HCC1954 確認 (1 week)
- C: 直接 pivot Path B low-F1 panel（同 cycle 內 continue）
【I】Best = B（最低風險 + 保留可選）
【D】Plan: 跑 HCC1954 chr8 zone 1 week；若 NEGATIVE → C 走 Path B
【E】Monitor: week 3 結束時 F1 metric；HCC1954 結果

---

Framework: DECIDE (Guo 2008《DECIDE: A Decision-Making Model for Health Care Managers》)
```
