# Data-Showcase Skeleton (= 既有 results-report)

> Hypothesis → Data → Stats → Caveats — 單實驗結果分析

```markdown
---
framework: Data-Showcase
when: <跑完單一實驗 / metric 解讀>
對應既有 skill: /results-report (thin wrapper)
---

# <Title — Assertion-Evidence>

## 1. Hypothesis

<Pre-registered hypothesis — what we expected to see>

- H0: ...
- H1: ...
- Pre-reg prior probability: <%>

## 2. Data

- Sample: <N + 樣本 list + selection criteria>
- Metric: <primary + secondary>
- Procedure: <method 摘要 + reference 詳述 doc>

## 3. Statistical Analysis

- Test: <e.g., Mann-Whitney U / Cohen d / LR cross-val>
- Effect size: <Cohen d / R² / OR>
- CI: <95% CI>
- p-value: <if applicable>
- Multiple-testing correction: <Bonferroni / FDR>

## 4. Result

<Verdict: POSITIVE / NEGATIVE / PARTIAL / INCONCLUSIVE>

對照 Pre-reg:
- ✓ Match: <observation matched prior expectation>
- ⚠ Surprise: <unexpected observation>

## 5. Caveats（誠實揭露）

- Sample limitation: <e.g., single-sample / cell-line bias>
- Confound: <未控制的變數>
- Generalization scope: <can / cannot extend to>
- Reproducibility: <other groups can replicate?>

## 6. Decision (對齊 /scientific-rigor §2-§7)

- Tier upgrade: <⭐X 不變 / ⭐X→Y>
- Next cycle: <follow-up plan>
- HARKing avoidance: 結論不超 pre-reg 範圍

---

Framework: Data-Showcase (InterSubMod results-report)
對應既有 skill: /results-report (thin wrapper 預設套此 framework)
```
