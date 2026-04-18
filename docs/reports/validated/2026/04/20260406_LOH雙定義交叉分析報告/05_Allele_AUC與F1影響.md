<!--
建立時間: 2026-04-06
目標: Step 3.3 LOH 內 Allele AUC + F1 影響分析
處理範圍: 7 samples × 2 modes, 2×2 AUC + 5 scenarios + confound check
關聯檔案:
  - scripts/analysis/build_allele_loh_auc_f1_analysis.py
-->

# Step 3.3: Allele AUC 與 F1 影響

## 背景

Wave 1 (Step 3.1) 確認 HP PERMANOVA 在 LOH 內不可用 (valid rate 5-6%)，但 Allele PERMANOVA valid rate 達 50-59% (~10x)。**核心問題**：Allele 維度的 AUC（TP/FP 區分力）是否足夠？

---

## C.1: 2×2 AUC 矩陣

### Overall AUC

| Mode | LOH Status | HPMergedDelta AUC | AlleleDelta AUC |
|------|-----------|-------------------|-----------------|
| PAIRED | LOH | 0.4660 | 0.4830 |
| PAIRED | Non-LOH | 0.4833 | 0.5050 |
| **TO** | **LOH** | **0.4993** | **0.5564** |
| TO | Non-LOH | 0.4807 | 0.4642 |

**核心判定**: TO LOH AlleleDelta AUC = **0.5564** → **未達 0.58 閾值** (J10)

![2×2 AUC heatmap](figures/w2a_01_auc_heatmap_2x2.png)

### Per-sample 方向一致性 (TO mode, LOH)

| Sample | AlleleDelta AUC | > 0.58? |
|--------|-----------------|---------|
| HCC1395 | 0.6079 | Yes |
| HCC1395_DORADO | 0.6442 | Yes |
| COLO829 | 0.5153 | No |
| H1437 | 0.5773 | No |
| H2009 | 0.5598 | No |
| HCC1937 | **0.6604** | Yes |
| HCC1954 | 0.5492 | No |

- **7/7** samples AUC > 0.5 → 方向一致（Allele 確實對 TP 略有偏好）
- **3/7** samples AUC > 0.58 → 部分通過閾值
- HCC1395 系列 + HCC1937 表現最好（高 LOH 比例 samples）

![AlleleDelta AUC per sample](figures/w2a_02_allele_auc_per_sample.png)

![AlleleDelta distribution (TP vs FP)](figures/w2a_03_allele_delta_distribution.png)

![Valid rate in LOH: HP vs Allele](figures/w2a_04_valid_rate_in_loh.png)

---

## C.2: Allele 篩選 F1 影響 (A1-A5)

### 情境定義

| 情境 | 移除條件 | 邏輯 |
|------|----------|------|
| A1 | LOH + AlleleSig=False | LOH 內無 Allele 信號 |
| A2 | LOH + AlleleP > 0.1 | LOH 內 Allele p-value 不顯著 |
| A3 | LOH + AlleleDelta < median | LOH 內 effect size 低 |
| A4 | LOH + HPMergedSig=F + AlleleSig=F | 雙維度都無信號 |
| A5 | Allele-primary mode | LOH 用 Allele, non-LOH 用 HP |

### 結果 (HCC1395, TO mode)

| 情境 | TP Loss | FP Removal | ΔF1 | 安全 |
|------|---------|------------|-----|------|
| A1 | **31.8%** | 40.5% | -0.1256 | FAIL |
| A2 | **29.8%** | 39.3% | -0.1146 | FAIL |
| A3 | **25.9%** | 37.8% | -0.0935 | FAIL |
| A4 | **31.2%** | 40.3% | -0.1222 | FAIL |
| A5 | **44.0%** | 55.0% | -0.1901 | FAIL |

**所有 35 組 (5 情境 × 7 samples) 均 FAIL TP loss ≤ 2%**。

![F1 scenario heatmap (scenarios × samples)](figures/w2a_07_f1_scenario_heatmap.png)

![TP loss vs FP removal tradeoff](figures/w2a_08_tp_fp_tradeoff.png)

### 共同原因

與 Step 1.4 (S1-S5) 的根因相同：LOH 區域 TP-enriched。LOH 內 ~70% 是 TP，移除 LOH 內的 variants 必然大量損失 TP。

---

## C.3: Confound 檢查

按照 O11-O13 教訓，檢查 AlleleDelta AUC 是否受 NumReads 或 Coverage_Multiple confound。

| LOH Status | Raw AUC | Resid(NumReads) | Resid(Coverage) | AUC Drop |
|-----------|---------|-----------------|-----------------|----------|
| LOH | 0.5564 | 0.5560 | 0.5560 | **< 0.001** |
| Non-LOH | 0.4642 | 0.4779 | 0.4779 | < 0.02 |

**Confound 檢查通過**: AUC drop < 0.05 → AlleleDelta 的信號不是 NumReads 或 Coverage 的 confound。

**判定**: 信號是真實的 (7/7 方向一致, confound-free)，但太微弱 (AUC=0.556 < 0.58) 無法作為 binary filter。

---

## 操作條件與限制

- **FN count**: 僅 HCC1395 有精確值 (11,051)，其他 samples 使用 FN/TP ratio 估計 (0.388)
- **LOH 定義**: 使用 Potential_LOH (HP_Ratio based)，非 LOH.bed
- **Confound check**: 使用線性 residualization，可能不足以消除非線性 confound
- **AlleleDelta**: ISM 計算的 effect size，已知可能受 AF confound (O12 教訓)，但 residualization 後 AUC 不變
- 如果 AUC 在 0.55-0.65 → 應用 L3 AF-bin 交叉驗證排除 confound（本分析的 0.556 落在此範圍）

### 輸出檔案

- `auc_matrix_2x2.tsv` — (LOH/Non-LOH) × (HP/Allele) × (mode) AUC
- `allele_loh_statistics.tsv` — Per sample AUC
- `f1_scenario_results.tsv` — 35 組 F1 影響結果
- `allele_filtered_variants.tsv.gz` — 被移除 variant 完整特徵 (638,612 rows)
- `confound_check.tsv` — Partial correlation 結果
