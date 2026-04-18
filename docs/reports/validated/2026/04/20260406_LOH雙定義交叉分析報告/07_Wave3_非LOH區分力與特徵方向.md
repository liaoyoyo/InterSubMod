<!--
建立時間: 2026-04-06
目標: Wave 3 M2+M3 結果 — 非 LOH 區域 TP/FP 區分力量化 + 特徵方向全景圖
處理範圍: 31 特徵 × (LOH/Non-LOH) × (Paired/TO) AUC 矩陣、mode-flip 特徵、修正模擬
關聯檔案:
  - scripts/analysis/build_non_loh_discrimination.py
  - scripts/analysis/build_feature_direction_map.py
  - big7_disk_output/synthesis/observation_workspaces/20260406_non_loh_discrimination/
  - big7_disk_output/synthesis/observation_workspaces/20260406_feature_direction_map/
-->

# Wave 3: 非 LOH 區分力與特徵方向全景圖

## M2: 非 LOH 區域 TP/FP 區分力量化

### 背景

Wave 1+2 聚焦 LOH 區域，但 **Q4 (Neither LOH) 佔 58.1%**。如果 non-LOH 已高度可解，就可以集中資源在 LOH 殘餘問題。

### 結果

#### TO mode 特徵 AUC 對比

| 排名 | 特徵 | Non-LOH AUC | LOH AUC | 差異 | 註記 |
|------|------|-------------|---------|------|------|
| 1 | HPFineNGroups | **0.6428** | 0.5997 | +0.043 | read count proxy |
| 2 | Coverage_Multiple | **0.6046** | 0.5500 | +0.055 | confound 嫌疑 |
| 3 | NumReads | **0.6046** | 0.5500 | +0.055 | confound 嫌疑 |
| 4 | HP1FamilyN | 0.5839 | 0.5507 | +0.033 | HP-derived |
| 5 | HP2FamilyN | 0.5603 | 0.4753 | +0.085 | |
| 6 | Quality_Score | 0.5487 | 0.5538 | -0.005 | |
| 7 | PairwiseMeanDist | 0.5474 | 0.5418 | +0.006 | |
| 8 | PairwiseMedianDist | 0.5404 | 0.5332 | +0.007 | |

#### 關鍵判定

| 判定 | 結論 |
|------|------|
| **J11: Non-LOH 區分力同樣有限** | 最高 AUC = 0.6428 (HPFineNGroups)，但為 read count proxy；排除 confound 後無特徵 > 0.58 |
| **J12: LOH/Non-LOH 差異不大** | 差異 < 0.06 個 AUC 單位，問題是全域性的而非 LOH 特異的 |

### 圖表

![LOH vs Non-LOH AUC 對比 heatmap](figures/w3_m2_01_auc_loh_vs_nonloh_heatmap.png)

![Non-LOH top-10 features ranked bar](figures/w3_m2_02_nonloh_top10_features.png)

![Per-sample Non-LOH AUC heatmap](figures/w3_m2_03_per_sample_nonloh_auc.png)

![FP rate 對比: LOH vs Non-LOH](figures/w3_m2_04_fp_rate_loh_vs_nonloh.png)

![AUC 差異分佈 (LOH minus Non-LOH)](figures/w3_m2_05_auc_delta_loh_minus_nonloh.png)

![TP/FP stacked bar by LOH status](figures/w3_m2_06_tp_fp_stacked_by_loh.png)

---

## M3: 特徵方向全景圖 + 可修正性評估

### Mode-Flip 特徵

在 LOH 區域中，以下特徵在 Paired vs TO 方向翻轉：

| 特徵 | Paired LOH AUC | TO LOH AUC | 方向 | 根因 |
|------|---------------|------------|------|------|
| **caller_af** | 0.535 (TP-favored) | **0.364** (FP-favored) | 翻轉 | Self-phasing: LOH 高 AF → FP 也高 AF |
| PairwiseMedianDist | 0.432 (FP-favored) | 0.533 (TP-favored) | 翻轉 | mode 間 label 比例差異 |
| PairwiseMeanDist | 0.436 (FP-favored) | 0.542 (TP-favored) | 翻轉 | 同上 |

### TO LOH 反向特徵 (AUC < 0.48)

| 特徵 | TO LOH AUC | 對 QS 的影響 |
|------|-----------|-------------|
| caller_af | 0.364 | QS 使用 AF 會反向加分 |
| NumCpGs | 0.443 | CpG 數量→FP 反而高 |
| GlobalP | 0.443 | 全域 p-value 反向 |
| HP2FamilyN | 0.475 | HP 第二族群大小 |

### 修正模擬結果

| 修正策略 | 原始 AUC | 修正後 AUC | 提升 |
|----------|----------|-----------|------|
| HPMergedDelta → AlleleDelta (LOH 區域) | 0.465 | 0.501 | +0.036 |
| caller_af → -abs(af-0.5) | 0.418 | 0.526 | +0.108 |
| Quality_Score LOH+AlleleDelta×20 | 0.497 | 0.501 | +0.004 |

**結論**: 修正後最高 AUC 僅 0.526，仍遠低於 0.58 可用門檻。特徵方向修正無法拯救 TO LOH 維度。

### AF-bin 交叉驗證 (L3)

已確認 caller_af 反向在所有 AF bin 中一致（非 bin-specific artifact），進一步驗證 self-phasing 根因。

### 圖表

![4-way AUC heatmap: features × (Paired/TO) × (LOH/Non-LOH)](figures/w3_m3_01_four_way_auc_heatmap.png)

![Mode-flip scatter: Paired AUC vs TO AUC](figures/w3_m3_02_mode_flip_scatter.png)

![AF-bin AUC 折線圖 — 反向特徵在各 AF bin 的表現](figures/w3_m3_03_af_bin_auc_lines.png)

![修正前後 AUC 對比](figures/w3_m3_04_correction_comparison.png)

![Per-sample 方向一致性 heatmap](figures/w3_m3_05_per_sample_consistency.png)
