# Example: 説服 PI 換方向 — Golden Circle + Pre-Mortem + WRAP

> 場景: lab meeting 説服 PI pivot, Who=PI + 同儕, Why=説服, What=方向選擇

```markdown
---
framework: Golden Circle (主) + Pre-Mortem (備) + WRAP (決策評估)
audience: PI + lab 同儕 (8 人)
when: 15min lab meeting
applied_by: narrative-frame v0.1
---

# 為什麼我們應該 pivot Cycle 5 至 chr8-specific zone gate（不再做 LR direction）

## WHY — 核心信念

我相信 **non-LR framework 比 LR threshold 更適合 cross-sample 推廣**，因為過去 3 個 cycle 顯示 caller_af direction-inconsistency 是 LOSO 災難主因。

> 業界引用「People don't buy WHAT you do; they buy WHY you do it」(Sinek)

## HOW — 我們怎麼做

- **Path A**（推薦）: chr8-specific zone gate
  - 不依賴 LR threshold — 用 boolean rule
  - chr8 hotspot 已有 7.4× FP enrichment evidence (HCC1395)
  - 與 caller_af direction 解耦

- **Path B**（備用）: low-F1 panel 驗 HCC1395 +0.00699 generalize
  - 仍 caller_af-bounded
  - 風險較低但 ROI 不如 Path A

## WHAT — 具體交付

- 1 個新 cycle (cycle_20260524-XXXX) registered
- 1 個 chr8 zone rule .py script
- 1 份 cross-sample evaluation report

---

## ⚠ Pre-Mortem（反推風險）

假設 Cycle 5 結束後失敗，我會看到什麼？

1. **chr8 zone rule 過 narrow** — 只在 HCC1395 有效，其他樣本沒 chr8 hotspot
2. **Boolean rule overfitting** — chr8 HP ratio 在 LOH 區隨機 noise
3. **PI 拒絕 pivot** — 因為 Cycle 1-3 已投入 LR direction，沉沒成本心理

**Reality-test 三反例**:
- (1) 跑 HCC1937 / HCC1954 — 看 chr8 hotspot 是否存在
- (2) shuffle chr8 region — rule 效果是否消失
- (3) 對 PI 直接 ack 沉沒成本，提早 1-on-1 對話而非 lab meeting 突襲

---

## WRAP 決策評估

### W — Widen

- Path A: chr8 zone gate（推薦）
- Path B: low-F1 panel
- Path C: 停 Cycle 5，做 Cycle 4 v2 重跑（vanishing option）
- Path D: 引入 outside data 重訓 LR（新方向）

### R — Reality-test

每 path 必須成立的假設：
- A: chr8 hotspot generalize ≥ 3 樣本 — L3 ⭐⭐⭐ partial evidence (HCC1395 only)
- B: low-F1 panel 仍 caller_af-bounded — L4 strong negative prior
- C: vanishing option — 對研究進度負影響
- D: outside data 引入 — 6 月後才到位

### A — Attain distance

10 分鐘後: 在意 Path A vs B
10 個月後: 在意 Cycle 5 是否解 cross-sample
10 年後: 在意 ISM framework 是否成型（不在意 single cycle）
→ 偏向 Path A（high upside if right）

### P — Prepare to be wrong

- **Tripwire**: 跑 2 樣本（HCC1937 + HCC1954），若 chr8 zone rule 在 ≥1 樣本 ΔF1 ≤ 0 → 立刻 pivot 至 Path B
- 設 2 週時間 budget；過了就 fall back

---

## Call to Action

請 PI ack:
1. Cycle 5 走 Path A（chr8 zone gate）
2. 2 週 tripwire 進入後若失敗 → pivot Path B
3. 我下週開始建 chr8 rule .py + 跑 HCC1937 / HCC1954

---

Framework: Golden Circle (Sinek《Start with Why》2009) + Pre-Mortem (Klein HBR 2007) + WRAP (Heath《Decisive》2013)
```

---

## 為什麼選此 hybrid（N2 推薦理由）

**5W 識別**:
- Who: PI + 同儕（混合 — 主要 PI）
- Why: **説服 pivot**（核心）
- What: 方向選擇 + 風險評估
- When: 15min lab meeting
- How: slide 或 .md companion

**主框架 Golden Circle** — Why-How-What 結構 適合 説服 + 起始於信念而非結論
**Sub: Pre-Mortem** — 反推失敗風險增加 credibility（不只樂觀）
**Sub: WRAP** — 決策評估 不只 1 path 而是多 path 對比 + tripwire 設計

**不選擇**:
- SCQA: 太結論先行，説服力較弱（不如 Golden Circle 的 emotional why）
- AIDA: 行銷框架對 PI 過於 sales-y
- Pure Pixar Spine: 缺決策評估 + tripwire
