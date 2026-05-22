# CPP Skeleton — Concept / Procedure / Principle

> Instructional Design Theory

```markdown
---
framework: CPP (Concept-Procedure-Principle)
when: <教學講義 / tutorial / 手冊 / SOP>
---

# <Topic>

## Concept（是什麼）

<一句定義 — what is X>

<2-3 句展開：scope / boundary / 為什麼存在>

## Procedure（怎麼做）

Step-by-step：
1. <Step 1>
2. <Step 2>
3. <Step 3>

✓ Verification: <怎麼驗證做對了>

## Principle（為什麼）

<深層機制 / 底層原理 — 解釋 Procedure 為什麼 work>

當 X 發生時，Y 必跟著發生，因為 <causal chain>。

---

範例（ISM filter 教學）:

**Concept**: ISM filter 用 epigenetic context 強化 cancer caller — 從 read-level methylation pattern 判斷 variant 真假。

**Procedure**:
1. Load ONT BAM + 5mC tags
2. 對每 variant 抓周圍 ±100bp CpG methylation
3. 跨 read 計算 methylation pattern 一致性 score
4. Score < threshold → 過濾為 FP

✓ Verification: TP recall ≥ 0.95 + FP removal ≥ 30%

**Principle**: 真 somatic mutation 起源於 1 個 founder cell；其 methylation pattern 在 division 中遺傳 → 所有 derived reads 帶相同 methyl signature。Sequencing error 隨機分佈，沒有 coherent signature → 一致性 score 低。

---

Framework: CPP (Instructional Design theory)
```
