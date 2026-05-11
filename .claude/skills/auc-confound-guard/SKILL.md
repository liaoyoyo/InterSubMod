---
name: auc-confound-guard
description: AUC 結論三關 confound 驗證協議 — within-group OLS + AF-bin 交叉 + permutation test。宣告任何 AUC>0.58 的特徵顯著性前必經此流程。觸發：「feature AUC」「discriminative」「predictive power」「TP vs FP 區分」「residualize」「confound」「特徵能否區分 TP/FP」。SKIP WHEN AUC<0.58 不需宣告顯著性、純描述性統計（無 TP/FP 分類目標）、未提及 AUC 的特徵觀察、bug fix / build / commit、純 docs 編輯。
allowed-tools: Read, Write, Bash, Grep
user-invocable: true
---

# AUC Confound Guard（AUC 結論三關驗證協議）

本 skill 定義 InterSubMod 宣告「某特徵能區分 TP/FP」之前必經的三關驗證。違反任一關 → 結論降級為 `annotation_only` 或 `confound_suspected`，不得寫入 MEMORY.md Concluded。

## 為何需要

歷史上 11 個 AUC 方向被 confound 翻盤（Agent 3 統計）：
- O11 epipolymorphism → 全因 n_reads confound
- O12 AlleleDelta → 與 AF 高度共線（L2 collider bias）
- O13 跨區域甲基化 → shared read count confound
- G1-G7 germline FP → AUC<0.64 但 FP removal=0%

**共同失敗模式**：pooled 統計下 AUC>0.58，分層/殘差化後歸零或消失。

## 三關流程

```
原始 AUC > 0.58?
  │
  ├─ 否 → 結論：not discriminative（tier 1，allowed）
  │
  └─ 是 → Gate 1: Within-Group OLS Residualization
            │
            ├─ 殘差 AUC > 0.55 → Gate 2: AF-Bin Stratification
            │                      │
            │                      ├─ AUC stable across AF bins → Gate 3: Permutation Test
            │                      │                              │
            │                      │                              ├─ p<0.05 → tier 3+
            │                      │                              └─ p≥0.05 → tier 1 (chance)
            │                      │
            │                      └─ AUC varies > 0.10 → CONFOUND（AF driven）
            │
            └─ 殘差 AUC ≤ 0.55 → CONFOUND（confound driven）
```

## Gate 1: Within-Group OLS Residualization

**目的**：移除 confound 變數（n_reads, AF, purity）對特徵的線性影響。

**必做**：
```python
# WRONG: pooled OLS (保留組間差)
residuals = feature - OLS(feature ~ n_reads)

# RIGHT: within-group OLS (每組各 fit)
residuals_tp = feature[tp_mask] - OLS(feature[tp_mask] ~ n_reads[tp_mask])
residuals_fp = feature[fp_mask] - OLS(feature[fp_mask] ~ n_reads[fp_mask])
residuals = concat(residuals_tp, residuals_fp)
```

**pitfall**：見 `/known-pitfalls` P-02。pooled OLS 會產生假信號。

**必 residualize 的 confound**：
- `n_reads`（read count 影響統計功效）
- `AF`（allele frequency 與 quality 共線）
- `coverage`（影響所有 methylation 量化精度）

**通過準則**：殘差 AUC > 0.55（相較原始 AUC 降幅 < 0.05）

## Gate 2: AF-Bin Stratification

**目的**：確認效應不是被 AF 分佈差異驅動。

**必做**：
```python
AF_BINS = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.0)]
for bin_low, bin_high in AF_BINS:
    mask = (AF >= bin_low) & (AF < bin_high)
    if sum(mask) < 50: continue  # 樣本量不足
    auc_bin = compute_auc(feature[mask], label[mask])
    print(f"AF [{bin_low},{bin_high}): AUC={auc_bin:.3f} n={sum(mask)}")
```

**通過準則**：
- 所有 bin 的 AUC 方向一致（都 >0.55 或都 <0.45）
- Bin 間 AUC 極差 < 0.10

**失敗情境**：AUC 在 low-AF bin 為 0.70、high-AF bin 為 0.50 → 信號來自 AF 分佈差異，不是特徵本身。

## Gate 3: Permutation Test

**目的**：確認 AUC 顯著性不是偶然。

**必做**：
```python
from sklearn.metrics import roc_auc_score
import numpy as np

observed_auc = roc_auc_score(label, feature)

null_aucs = []
for _ in range(1000):
    shuffled_label = np.random.permutation(label)
    null_aucs.append(roc_auc_score(shuffled_label, feature))

p_value = np.mean(np.array(null_aucs) >= observed_auc)
```

**通過準則**：`p_value < 0.05`

**採樣數**：最少 1000 次；若樣本量 <500 時用 5000 次。

## 輸出規範

任何通過三關的 AUC 結論，**必須在 evidence_ledger.jsonl 同時填寫**：

```json
{
  "tier": 3,
  "tier_flags": ["within_group_ols", "af_bin_stratified", "permutation_tested"],
  "confidence_cap": 3,
  "auc_raw": 0.XX,
  "auc_residualized": 0.XX,
  "auc_by_af_bin": {"0.0-0.1": 0.XX, "0.1-0.2": 0.XX, ...},
  "permutation_p": 0.0XX,
  "confounds_controlled": ["n_reads", "AF"]
}
```

## 快速檢核（AI 自我審查用）

提到 AUC 結論前先問自己：

```
☐ Gate 1: residualize 用 within-group OLS 嗎？（不是 pooled OLS）
☐ Gate 1: residualize 了 n_reads、AF、coverage 嗎？
☐ Gate 2: AF-binned AUC 方向一致嗎？
☐ Gate 3: permutation test 跑了嗎？p<0.05?
☐ 三關都過才能宣告 tier ≥ 3
☐ 任一關失敗 → verdict=annotation_only（非 positive）
```

## 回報範本

結論寫入報告時引用：

> 本結論通過 auc-confound-guard 三關驗證（within-group OLS 殘差 AUC={X}；AF-bin 極差={Y}；permutation p={Z}）。tier=3。

## 相關資源

- `/known-pitfalls` P-01, P-02, P-06（collider bias、pooled OLS、n_reads confound）
- `docs/standards/evidence_tier_rubric.md`（tier 分級規範）
- MEMORY: `feedback_L2_collider_bias.md`、`feedback_pooled_ols_residualization_trap.md`
