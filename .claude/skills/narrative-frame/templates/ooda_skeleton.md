# OODA Loop Skeleton

> John Boyd, USAF — Observe / Orient / Decide / Act

```markdown
---
framework: OODA Loop
when: <高速決策 / 紅藍對抗 / incident response / 軍事戰術>
---

# OODA Cycle <N>

## 1. Observe（觀察）

<目前狀況 — sensor data / metric / 對手動作>

## 2. Orient（定向）

- **Pattern match**: <對應已知 pattern？>
- **Mental model update**: <基於新觀察修正世界觀>
- **Cultural / heritage / new info / analysis / synthesis**: <Boyd 5 inputs>

## 3. Decide（決策）

<從可選 actions 中選 1 個>

## 4. Act（行動）

<執行 — fast；不完美也優於 hesitation>

→ 立刻回到 Observe，跑下一個 OODA cycle

---

業界引用: 「The side that completes OODA faster wins」(Boyd)

範例（incident response）:

- Observe: pipeline 跑 30min 後 OOM；磁碟剩 0%
- Orient: pattern match → 上次也是 /tmp 寫滿；mental model: 此 pipeline 預設 /tmp out
- Decide: 改 export TMPDIR=/big7_disk 重跑 (不等 root cause 完整分析)
- Act: 立刻 kill + restart with TMPDIR set
- → 新 cycle 觀察：是否 30min 後再 OOM？

---

Framework: OODA Loop (John Boyd, USAF 1976-1995)
```
