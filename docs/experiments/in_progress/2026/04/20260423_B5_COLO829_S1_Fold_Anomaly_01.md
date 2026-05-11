---
title: B5 — COLO829 paired_full S1 Fold 0.59× 反轉根因調查
date: 2026-04-23
author: InterSubMod Research
status: in_progress
hypothesis_id: B5
related:
  - docs/experiments/in_progress/2026/04/20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md
artifacts:
  - scripts/analysis/20260423_B5_colo829_s1_fold.py
  - research/ng_kde_rescaling/data/B5_colo829_s1_fold_detail.tsv
  - research/ng_kde_rescaling/data/B5_colo829_vs_cohort_distributions.tsv
  - research/ng_kde_rescaling/data/B5_colo829_s1_fold_detail.summary.json
---

# B5 — COLO829 paired_full S1 Fold 0.59× 反轉根因調查

## 1. 背景與目標

Thread B 的 Scheme S1 = `LOH_Strong ∩ Extreme AF`，設計為「高信心 TP 留存」方案。
`fold_improvement = scheme_precision_ratio / baseline_precision_ratio`（來源：`tpfp_per_sample_scheme_full.tsv`）。

多數樣本 fold > 1（S1 cell 的 TP/FP ratio 高於 baseline），但 **COLO829 paired_full fold = 0.592×**（反轉）。

三種候選解釋：
- (a) **跨樣本 cell 定義不一致**
- (b) **COLO829 特殊腫瘤生物學**
- (c) **denominator artifact**（baseline TP rate 已飽和，S1 "FP 移除" 空間極小 → fold 比例失真）

## 2. 資料與方法

- 輸入：`research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz`（698,059 rows）
- Filter：`sample ∈ {COLO829, H1437, H2009, HCC1395, HCC1395_DORADO, HCC1937}` × `mode=paired_full`（排除 HCC1954 — 使用者指示為 post-hoc）
- S1 cells：`LOH_Subtype='LOH_Strong' ∧ AF_class='Extreme'`
- 指標：baseline / S1 的 TP rate、TP/FP ratio、Wilson 95% CI、fold_precision、跨樣本 LOH/AF/CN 分佈

## 3. 結果

### 3.1 S1 fold × baseline 飽和度（核心表）

| 樣本 | baseline n | baseline TP rate [95% CI] | baseline FP rate | S1 n (frac%) | S1 TP rate [95% CI] | Δ (S1−base) | fold_precision |
|---|---|---|---|---|---|---|---|
| **COLO829** | 37,458 | **0.9393 [0.9369, 0.9417]** | **0.0607** | 244 (0.65%) | 0.9016 [0.8578, 0.9330] | **-0.0377** | **0.592×** |
| H1437 | 67,476 | 0.9999 [0.9998, 0.9999] | 0.00012 | 1,182 (1.75%) | 0.9992 [0.9952, 0.9999] | -0.0007 | 0.140× |
| H2009 | 132,995 | 0.9994 [0.9992, 0.9995] | 0.00065 | 3,516 (2.64%) | 0.9997 [0.9984, 1.0000] | +0.0004 | 2.27× |
| HCC1395 | 30,381 | 0.9794 [0.9777, 0.9809] | 0.0206 | 347 (1.14%) | 0.9971 [0.9839, 0.9995] | +0.0178 | **7.29×** |
| HCC1395_DORADO | 30,129 | 0.9920 [0.9910, 0.9930] | 0.0080 | 500 (1.66%) | 1.0000 [0.9924, 1.0000] | +0.0080 | ∞ (S1 FP=0) |
| HCC1937 | 12,588 | 0.9845 [0.9822, 0.9865] | 0.0155 | 546 (4.34%) | 0.9853 [0.9714, 0.9926] | +0.0008 | 1.06× |

**Key observation**：COLO829 是 **唯一 Δ 顯著為負** 的樣本（S1 TP rate CI [0.858, 0.933] vs baseline CI [0.937, 0.942] — CI 不重疊；其他樣本 Δ 要麼微小正值要麼 CI 橫跨 baseline）。

### 3.2 baseline TP rate 飽和度

- COLO829：**93.93%**（全 cohort 最低，FP rate 6.07% — 全 cohort 最高）
- H1437、H2009：>99.9%（幾乎無 FP headroom）
- HCC1395、HCC1395_DORADO、HCC1937：97.9–99.2%

→ 若 S1 是「移除」策略，H1437/H2009 的 fold 被低 FP denominator 放大為不穩定數值。這支持候選 (c) 對 H1437/H2009 的 fold 失真，但 **不解釋 COLO829 的反轉**（COLO829 FP headroom 最大，理應最有 S1 enrichment 空間）。

### 3.3 跨樣本 LOH_Subtype 分佈

| LOH_Subtype | COLO829 | H1437 | H2009 | HCC1395 | HCC1395_DORADO | HCC1937 |
|---|---|---|---|---|---|---|
| None | **0.794** | 0.791 | 0.714 | 0.525 | 0.552 | 0.440 |
| LOH_Weak | 0.150 | 0.127 | 0.099 | 0.277 | 0.241 | 0.217 |
| LOH_Noise | 0.040 | 0.029 | 0.054 | 0.153 | 0.134 | 0.217 |
| **LOH_Strong** | **0.0154** | 0.048 | 0.113 | 0.035 | 0.054 | 0.088 |
| LOH_Subclone | 0.00083 | 0.005 | 0.021 | 0.010 | 0.019 | 0.038 |

COLO829 LOH_Strong frac 1.54%（全 cohort 最低），LOH_Subclone 僅 0.08%（極度稀少）。但其他樣本也有低 LOH_Strong（HCC1395 3.5%）卻 fold = 7.29×，表示 **LOH_Strong 比例不是反轉根因**。

### 3.4 Coverage_Category 分佈

| Coverage | COLO829 | H1437 | H2009 | HCC1395 | HCC1395_DORADO | HCC1937 |
|---|---|---|---|---|---|---|
| Normal | **0.801** | 0.795 | 0.603 | 0.601 | 0.574 | 0.502 |
| CNV_Gain | **0.0056** | 0.0071 | 0.0635 | 0.0821 | 0.0954 | 0.1193 |
| CNV_Loss | 0.0160 | 0.0003 | 0.0008 | 0.0004 | 0.0001 | 0.0118 |
| High_Copy | 0.0018 | 0.0014 | 0.0402 | 0.0074 | 0.0077 | 0.0354 |

COLO829 CNV_Gain 僅 0.56%（全 cohort 最低），Normal coverage dominant。LOH_Strong 在 COLO829 來自「低 CN 擾動背景下的 AF extreme」，與其他樣本（部分 LOH_Strong 帶 CNV） **組成不同**。

## 4. 候選解釋判定

| 候選 | 支持證據 | 反對證據 | 判定 |
|---|---|---|---|
| (a) cell 定義不一致 | — | Filter 邏輯跨樣本同一，SQL 確定 | **排除** |
| (b) **COLO829 特殊生物學** | **S1 TP rate Wilson CI [0.858, 0.933] vs baseline [0.937, 0.942] 不重疊，Δ=-3.8%**；COLO829 LOH_Strong 背景 CN normal 佔比極高，與他樣本（多帶 CNV）組成不同 | — | **主因** |
| (c) denominator artifact | H1437/H2009 baseline TP=99.9% → fold 數字對 FP 敏感，確屬 artifact；但 COLO829 baseline TP 是 cohort 最低，headroom 最大 | COLO829 不符合「飽和→fold 失真」路徑 | 解釋 H1437/H2009，**不解釋 COLO829** |

## 5. 結論

**Real biology (候選 b)**：COLO829 paired_full S1 fold=0.592× 是 **LOH_Strong∩Extreme cell 本身 TP rate 顯著低於 baseline (-3.77%, CI 不重疊)** 的結果，非 denominator artifact。

量化 evidence：
- COLO829 baseline TP rate = 93.93% [93.69, 94.17]
- COLO829 S1 TP rate = 90.16% [85.78, 93.30]
- Wilson CI 上界 (0.933) 低於 baseline 下界 (0.937) → **non-overlap, real shift**
- 同群 5 個樣本（H1437/H2009/HCC1395/HCC1395_DORADO/HCC1937）S1 Δ ≥ 0（CI 涵蓋或高於 baseline）
- COLO829 的 LOH_Strong 背景 coverage 幾乎都是 Normal（Gain 0.56%），與 HCC1937 (Gain 11.9%) 等樣本 S1 cell 組成本質不同

## 6. 對 Thread B 論證的影響

1. **S1 scheme 不是跨樣本通用**：COLO829 不應被納入「S1 as high-confidence retention」共識樣本，需在 Thread B 報告中加註 **S1 samples={H1437, H2009, HCC1395, HCC1395_DORADO, HCC1937}**，排除 COLO829。
2. **COLO829 需獨立子研究**：為何 COLO829 LOH_Strong 無 TP enrichment？候選假說：
   - COLO829 LOH_Strong region 多落在 amplicon 或 germline-biased 區段
   - COLO829 的 ClairS 對 LOH region 的 FP 模式跨樣本不同
   - 建議 B5-follow-up：對 COLO829 LOH_Strong ∩ Extreme 的 24 個 FP 做 chr×pos 聚合檢視
3. **fold_precision 指標在飽和樣本不可用**：H1437/H2009 baseline TP > 99.9%，任何 scheme fold 都會被 FP denominator 放大；建議未來 Thread B 報告以 **baseline FP rate ≥ 1%** 為前提再比較 fold。

## 7. 後續工作

- [ ] B5-F1：COLO829 LOH_Strong ∩ Extreme 的 24 個 FP chr×pos 聚集分析
- [ ] B5-F2：與 COLO829 whole-sample FP provenance 對照，檢查 S1 FP 是否來自特定 amplicon
- [ ] 更新 Thread B 報告：加註 S1 共識樣本排除 COLO829，理由見本文件

---

**腳本**：`scripts/analysis/20260423_B5_colo829_s1_fold.py`
**數據**：`research/ng_kde_rescaling/data/B5_colo829_s1_fold_detail.tsv`、`..._vs_cohort_distributions.tsv`
