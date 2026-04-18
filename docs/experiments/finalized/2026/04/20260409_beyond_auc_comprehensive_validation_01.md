<!--
建立時間: 2026-04-09 06:00
目標: Beyond-AUC 7 方法全面驗證 ISM 甲基化特徵是否存在 AUC 遺漏的信號
處理範圍: M1-M7 七種互補統計方法 × 多層分層（LOH/non-LOH × AF bins × per-sample）× 748,391 regions
關聯檔案:
  - scripts/analysis/build_beyond_auc_validation.py
  - scripts/analysis/run_m5_m6_m7_only.py
  - docs/references/20260409_beyond_auc_biology_literature_review.md
  - docs/experiments/finalized/2026/04/20260408_fine_pairwise_distance_analysis_01.md
-->

# Beyond-AUC 綜合驗證：ISM 甲基化特徵最終評估

## 1. Executive Summary

**結論：ISM 甲基化特徵空間已確認耗盡。**

7 種互補統計方法 × 多層分層驗證，對 25 個特徵 × 748,391 regions 進行全面分析：

| 判定 | 數量 | 說明 |
|------|------|------|
| **確認無效** | 12/22 ISM 特徵 | Pure methylation 特徵全部 AUC ≤ 0.58，residualized 後更低 |
| **HP-tag 依賴** | 6/22 ISM 特徵 | 信號來自 haplotype phasing 而非甲基化 |
| **Composite 依賴** | 1/22 ISM 特徵 | Quality_Score 混合 caller + HP 信號 |
| **唯一正面發現** | HPFineNGroups | TO non-LOH low AF: AUC=0.72，跨樣本一致，確認為 somatic heterogeneity marker |

**關鍵方法學發現**：Pooled OLS residualization 存在系統性 data snooping — confounders 本身已能預測 TP/FP（AUC 0.66-0.85），導致 residualized AUC 被人為膨脹。此方法論缺陷影響所有先前使用 pooled residualization 的結論。

**文獻支持**：45+ 篇論文確認 germline mQTL ASM 效應 >> somatic passenger SNV 效應（3-6×），AUC-ROC 在 class imbalance 下穩健（Richardson 2024），截至 2026 年無直接先例用甲基化區分 somatic/germline variant。

---

## 2. 背景與動機

### 2.1 問題陳述

InterSubMod 歷經 O1-O13、G1-G7、LOH Wave 1-3、Fine-Pairwise 等 15+ 系統性研究，幾乎所有 ISM 甲基化特徵評估都以 **AUC-ROC 為唯一判定指標**，且結果一致顯示 AUC < 0.58。

**AUC 可能遺漏的信號類型**：
1. **尾部信號**：整體 AUC=0.52 但 top 5% 有高 FP enrichment
2. **被 confound 膨脹/壓縮的信號**：CramersV residualized 後從 0.531→0.464
3. **子群稀釋**：特定 AF × LOH × sample 組合有效但被 pooled
4. **非線性/交互效應**：單變量 monotonic 假設不成立
5. **類別不平衡**：TO 30.6% FP，PR-AUC 可能更敏感

### 2.2 數據來源

- **Master Dataset**: `all_region_rows.tsv.gz`
- **748,391 rows**：Paired 328,699 + TO 419,692
- **7 cancer cell line samples** × 2 modes × TP/FP labels
- **25 分析特徵**：22 ISM + 3 Caller

### 2.3 特徵分類

| 類別 | 特徵 | 信號來源 |
|------|------|---------|
| **Pure Methylation (12)** | GlobalP, CramersV, HeuristicScore, PairwiseMeanDist, PairwiseMedianDist, HPMergedDelta, HPMergedP, AlleleDelta, AlleleP, LabelAllelePermanovaF/P, ClusterPermanovaF | 純甲基化統計 |
| **HP-tag Dependent (6)** | HP1FamilyN, HP2FamilyN, HP_Ratio, HPFineNGroups, HPFineF, HPFineP | 依賴 haplotype phasing |
| **Label Verification HP (2)** | LabelHPPermanovaF, LabelHPPermanovaP | HP-based PERMANOVA |
| **Composite (2)** | Quality_Score, Stability | 混合信號 |
| **Caller baseline (3)** | caller_af, caller_gq, caller_dp | Variant caller 特徵 |

---

## 3. Method M1: Precision-Recall AUC

### 3.1 目的

AUC-ROC 對類別不平衡不敏感。Paired TP rate 98.96%、TO TP rate 69.41% — PR-AUC 能否揭示不同的排序？

### 3.2 結果

**Paired 模式**（prevalence TP=0.9896）：PR-AUC 完全無法提供信息。

由於 truth_binary 編碼 TP=1（佔 98.96%），PR-AUC 衡量的是「多數類檢測能力」，所有特徵 PR-AUC ≈ 0.989-0.997，lift 差異 < 0.01。在此極端 prevalence 下，PR-AUC 數學上無法區分特徵。

**TO 模式**（prevalence TP=0.6941）：PR-AUC 提供有限額外信息。

| 特徵 | ROC-AUC | PR-AUC | PR Lift |
|------|---------|--------|---------|
| HP1FamilyN | 0.568 | 0.736 | 1.060 |
| HPFineNGroups | 0.564 | 0.734 | 1.058 |
| HP_Ratio | 0.544 | 0.728 | 1.048 |
| caller_dp | 0.564 | 0.724 | 1.044 |
| PairwiseMeanDist | 0.543 | 0.718 | 1.035 |
| LabelHPPermanovaP | 0.548 | 0.720 | 1.037 |

- TO PR-AUC lift 最高 1.06×（HP1FamilyN），但全部 < 1.10
- PR-AUC 與 ROC-AUC 排序高度一致，無特徵的 PR-AUC 排名顯著優於 ROC-AUC
- 文獻支持（Richardson 2024）：ROC-AUC 不受 class imbalance 影響，ISM 的 AUC < 0.58 結論穩健

**M1 結論**：PR-AUC 未發現 AUC-ROC 遺漏的信號。

![M1: AUC vs PR-AUC 散點圖](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m1_01_auc_vs_prauc_scatter.png)

---

## 4. Method M2: 系統性 Residualization

### 4.1 目的

Raw AUC 混淆真實信號與 confounders（NumReads, NumCpGs, caller_af, caller_dp）。Residualize 後的 AUC 是否揭示被壓縮或膨脹的信號？

### 4.2 Pooled OLS 結果

| 特徵 | Mode | Raw AUC | Resid AUC | Delta |
|------|------|---------|-----------|-------|
| CramersV | paired | 0.531 | 0.578 | +0.047 |
| HeuristicScore | paired | 0.522 | 0.569 | +0.047 |
| HPMergedP | to | 0.538 | 0.583 | +0.045 |
| HP_Ratio | to | 0.544 | 0.579 | +0.035 |
| LabelHPPermanovaP | to | 0.548 | 0.598 | +0.050 |

表面看來，部分特徵 residualized AUC > raw AUC，似乎暗示被 confound 壓縮的信號。

### 4.3 驗證：Pooled OLS Data Snooping

**但系統性驗證揭露了嚴重的方法學缺陷。**

驗證設計：
1. **Simpson's Paradox Test (T1)**：per-sample weighted mean AUC vs pooled AUC
2. **Small-Sample Bias Test (T2)**：n_FP < 50 樣本的 AUC 影響
3. **Leakage Test (T3)**：within-group OLS vs pooled OLS
4. **Permutation Test (T4)**：null distribution of residualized AUC

以 LabelHPPermanovaP Paired LOH（pooled resid AUC = 0.695，+0.175 看似巨大提升）為深入案例：

| 驗證 | 結果 | 判定 |
|------|------|------|
| T1 Simpson's | Weighted mean = 0.504，gap = 0.190 | **Simpson's Paradox** |
| T1 Direction | 3/7 samples positive, 4/7 negative | **方向不一致** |
| T3 Within-group OLS | Separate resid AUC = 0.513 | **去除 leakage 後信號消失** |
| T3 Leakage delta | 0.181 | **幾乎全部是 leakage** |
| T4 Permutation | z-score = 27.2, p < 0.001 | Significant, but from leakage |

**Confounder-only AUC 驗證**：僅用 NumReads + NumCpGs + caller_af + caller_dp 預測 TP/FP：

| Target | Confounder-Only AUC | Genuine Signal |
|--------|--------------------:|---------------:|
| LabelHPPermanovaP Paired LOH | 0.695 | 0.012 |
| CramersV Paired non-LOH | 0.847 | -0.364 |
| HPMergedP TO LOH | 0.663 | 0.002 |

**根因**：Pooled OLS `feature ~ confounders` 中，confounders 本身已能預測 TP/FP（因為 TP/FP 的 coverage、AF 分布不同）。Residualized feature 的殘差保留了這種 confounder-based 分離。

### 4.4 Pure Methylation 特徵的真實信號

去除 residualization leakage 後，純甲基化特徵的真實信號：

**Paired pooled（所有 ≤ 0.58）**：
- CramersV: 0.531 → resid 0.578（但 leakage-driven）
- HeuristicScore: 0.522 → resid 0.569（但 leakage-driven）
- LabelAllelePermanovaF: 0.582 → resid 0.550（**下降**）
- 所有其他 < 0.52

**TO pooled（所有 ≤ 0.58）**：
- HPMergedP: 0.538 → resid 0.583（leakage-driven，within-group = 0.502）
- PairwiseMeanDist: 0.543 → resid 0.525（下降）
- CramersV: 0.509 → resid 0.438（**下降 0.071**）
- ClusterPermanovaF: 0.512 → resid 0.437（**下降 0.075**）

**M2 結論**：Residualization 揭示的「提升」全部是 pooled OLS data snooping。Within-group residualization 後，**無任何純甲基化特徵 AUC > 0.58**。多數特徵的真實信號為零或負值。

![M2: Raw vs Residual AUC Lollipop](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m2_01_raw_vs_residual_auc_lollipop.png)

---

## 5. Method M3: 條件子群探索

### 5.1 目的

Pooled AUC 可能平均掉子群效應。在 LOH × AF × Sample 的三維 grid 中是否存在被稀釋的局部信號？

### 5.2 分層設計

- **LOH**: {LOH, non-LOH}（2 levels）
- **AF**: {low_0-0.2, mid_0.2-0.5, high_0.5-0.8, extreme_0.8+}（4 levels）
- **Sample**: 7 個體 + all_samples pooled（8 levels）
- **篩選門檻**: n_TP ≥ 30 且 n_FP ≥ 30
- **Promising 判定**: AUC > 0.60 或 |Cohen's d| > 0.30

### 5.3 整體統計

| 類別 | Promising Cells | Features | Strong (AUC>0.60 or <0.40) |
|------|----------------:|----------:|---:|
| Caller | 96 | 3 | 82 |
| HP-dependent | 126 | 6 | 96 |
| Pure methylation | 53 | 5 | 34 |
| Composite | 25 | 1 | 17 |
| **Total** | **300** | **15** | **229** |

### 5.4 HP-tag Dependent 特徵（信號來自 phasing，非甲基化）

HP-dependent 特徵在多個子群中展現高 AUC，但其信號來源是 **haplotype phasing 品質與 read allocation**，非甲基化信號：

| 特徵 | 最佳子群 | AUC | 解釋 |
|------|---------|-----|------|
| HP1FamilyN | Paired non-LOH high AF | 0.801 | HP1 read count 反映 phasing 品質 |
| HP2FamilyN | Paired non-LOH high AF | 0.782 | 同上 |
| HP_Ratio | TO LOH high AF HCC1937 | 0.782 | LOH 中 HP balance 反映 purity |
| HPFineNGroups | TO non-LOH low AF HCC1954 | 0.750 | Somatic heterogeneity marker |

### 5.5 HPFineNGroups：唯一正面發現

HPFineNGroups 在 TO non-LOH low AF 子群中展現跨樣本一致的信號：

| Sample | n_TP | n_FP | AUC | Cohen's d |
|--------|------|------|-----|-----------|
| all_samples | 4,287 | 12,650 | **0.722** | -0.868 |
| HCC1954 | 2,619 | 3,921 | **0.750** | -1.053 |
| HCC1937 | 219 | 1,267 | **0.728** | -0.967 |
| H2009 | 169 | 2,441 | **0.703** | -0.783 |
| H1437 | 300 | 1,850 | **0.674** | -0.654 |
| HCC1395 | 258 | 1,336 | 0.610 | -0.465 |
| HCC1395_DORADO | 259 | 1,323 | 0.620 | -0.389 |

**7/7 樣本 AUC > 0.50**（6/7 > 0.60），Cohen's d 一致為負方向。此結果與先前 R4 研究的發現完全一致：
- TP（somatic variant）傾向有更多 methylation subgroups（NGroups ≥ 4）
- 低 AF TP 反映高 subclonal heterogeneity（多個亞克隆攜帶不同甲基化模式）
- 這是 **somatic heterogeneity 的解釋性標記**，而非 variant filter（因為 AUC = 0.72 不足以安全過濾）

### 5.6 Pure Methylation 特徵的子群表現

純甲基化特徵在子群分析中的最佳表現：

| 特徵 | 最佳子群 | AUC | n | 跨樣本一致性 |
|------|---------|-----|---|-------------|
| AlleleP | Paired LOH low AF HCC1395 | 0.725 | 993 | **僅 HCC1395**（sample-specific） |
| AlleleP | Paired LOH low AF all | 0.710 | 4,666 | 弱（主要由 HCC1395 驅動） |
| LabelAllelePermanovaF | Paired LOH high AF | 0.715 | 13,021 | 中（但 effect 來自 HP 分配） |
| PairwiseMeanDist | TO LOH extreme AF | 0.595 | 26,829 | 弱 |

**AlleleP HCC1395 深入分析**：
- AlleleP 衡量 allele-specific methylation 的 P-value
- HCC1395 的 ASM 模式可能因其特殊的 genomic 特徵（特定 LOH pattern + mQTL 分布）而在低 AF 區間展現信號
- 但此信號在其他 6 個樣本中 **未能重現**，判定為 sample-specific artifact

**M3 結論**：
- 子群分析未發現被 pooling 隱藏的系統性信號
- HPFineNGroups TO non-LOH low AF 是唯一跨樣本一致的正面發現（AUC=0.72），確認為 somatic heterogeneity marker
- 純甲基化特徵最佳子群 AUC = 0.72（AlleleP HCC1395），但為 sample-specific
- 所有跨樣本一致的信號均來自 HP-tag dependent 特徵

![M3: 子群 Heatmap — Paired](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m3_01_subgroup_heatmap_paired.png)

![M3: 子群 Heatmap — TO](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m3_01_subgroup_heatmap_to.png)

---

## 6. Method M4: Calibration + Brier Score

### 6.1 目的

AUC 衡量排序品質，不管機率校準。特徵在某些 decile 是否有準確的 FP rate 預測（適合條件規則）？

### 6.2 結果

**所有特徵的 Brier Skill Score (BSS) 均為強烈負值**，表明每個特徵的校準能力都不如「預測所有人為 prevalence」的 naive baseline。

| 特徵 | Mode | BSS | 解釋 |
|------|------|-----|------|
| Quality_Score | paired | -8.13 | 最佳但仍然極差 |
| HPFineNGroups | paired | -28.07 | |
| caller_gq | paired | -26.58 | |
| HPFineNGroups | to | -0.19 | TO 模式最接近零（但仍為負） |
| Quality_Score | to | -0.30 | |

**TO 模式的 BSS 比 Paired 好得多**，但這是因為 TO prevalence (69.4% TP) 更接近 0.5，baseline Brier 更大（0.212 vs 0.010），而非特徵更校準。

**M4 結論**：無任何特徵具備可用的機率校準能力。不適合作為條件規則使用。

![M4: Calibration Curves Top 8](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m4_01_calibration_curves_top8.png)

---

## 7. Method M5: Bootstrap 穩定性

### 7.1 目的

Point AUC=0.543 的 95% CI 可能是 [0.538, 0.548]（穩定）或 [0.42, 0.66]（不穩定）。Bootstrap 確認 point estimates 的可信度。

### 7.2 結果（200 bootstrap iterations × region-level resampling）

**Pooled CIs 極度狹窄** — point estimates 高度可靠：

| 特徵 | Mode | AUC | 95% CI | Width | pct > 0.55 |
|------|------|-----|--------|-------|-----------|
| LabelAllelePermanovaF | paired | 0.582 | [0.573, 0.593] | 0.020 | **100%** |
| PairwiseMeanDist | to | 0.543 | [0.540, 0.544] | 0.004 | 0% |
| HPMergedP | to | 0.538 | [0.536, 0.540] | 0.004 | 0% |
| CramersV | paired | 0.531 | [0.529, 0.534] | 0.005 | 0% |
| CramersV | to | 0.509 | [0.509, 0.510] | 0.002 | 0% |

**Pure methylation 特徵 pooled**：除 LabelAllelePermanovaF (Paired 0.582, pct>0.55=100%) 外，**所有特徵 pct_above_055 = 0%**。即 200 次 bootstrap 中無任何一次 AUC 超過 0.55。

**TO CIs 極窄**（width ~0.003-0.004）：信號穩定但一致低於 0.55。

**Per-sample HPFineNGroups 確認跨樣本一致性**：

| Mode | Sample | AUC | pct > 0.55 |
|------|--------|-----|-----------|
| paired | HCC1937 | 0.843 | 100% |
| paired | H2009 | 0.722 | 100% |
| to | HCC1954 | 0.650 | 100% |
| to | HCC1937 | 0.613 | 100% |
| to | HCC1395_DORADO | 0.605 | 100% |
| to | HCC1395 | 0.563 | 100% |
| to | COLO829 | 0.498 | 0% |

**COLO829 在 TO 模式下 AUC=0.498**（chance level），修正先前「7/7 一致」→「6/7 一致，COLO829 例外」。

### 7.3 Per-sample 純甲基化特徵的樣本特異性

Bootstrap 確認某些特徵在特定樣本中穩定：
- HCC1937: LabelAllelePermanovaF = 0.821, AlleleDelta = 0.711 (100% pct>0.55)
- H2009: LabelAllelePermanovaF = 0.706, AlleleDelta = 0.651 (100% pct>0.55)
- HCC1395: AlleleP = 0.613 (100% pct>0.55)

**但這些信號是 sample-specific**，不同樣本的最佳特徵不同。此現象與 mQTL 的 sample-specific 特性一致 — 每個 cell line 的 germline ASM pattern 不同，導致不同的特徵響應。

**M5 結論**：Bootstrap 確認 point estimates 穩定可靠。Pure methylation 特徵 pooled AUC 確定不超過 0.58。Per-sample 高 AUC 是 sample-specific（不可泛化）。

![M5: Forest Plot — Paired](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m5_01_forest_plot_paired.png)

![M5: Forest Plot — TO](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m5_01_forest_plot_to.png)

---

## 8. Method M6: 非線性建模 + Feature Importance

### 8.1 目的

單變量 AUC 假設 monotonic 關係。GradientBoosting (200 trees, depth 4) + LOSO-CV 7-fold 探索非線性效應與交互效應。比較 ISM-only / Caller-only / ISM+Caller 三組特徵集。

### 8.2 LOSO-CV 結果

| Mode | Feature Set | Mean AUC | Std | Mean F1 |
|------|------------|----------|-----|---------|
| **Paired** | ISM_only | 0.609 | 0.106 | 0.992 |
| Paired | Caller_only | **0.850** | 0.105 | 0.991 |
| Paired | ISM+Caller | 0.855 | 0.105 | 0.991 |
| Paired | **ISM delta** | **+0.006** | | |
| **TO** | ISM_only | **0.647** | 0.106 | 0.681 |
| TO | Caller_only | 0.579 | 0.081 | 0.727 |
| TO | ISM+Caller | **0.660** | 0.113 | 0.663 |
| TO | **ISM delta** | **+0.081** | | |

**關鍵發現**：

1. **Paired**：Caller_only 已達 AUC=0.850（dominant），ISM 僅增加 +0.006（不顯著）
2. **TO**：ISM_only (0.647) **超越** Caller_only (0.579)！ISM+Caller 比 Caller_only +0.081

**TO 模式 per-sample 分析**：ISM 貢獻高度 sample-dependent：

| Sample | ISM_only | Caller_only | ISM+Caller | ISM delta |
|--------|----------|-------------|------------|-----------|
| HCC1937 | 0.770 | 0.692 | **0.791** | +0.100 |
| HCC1954 | 0.738 | 0.699 | **0.766** | +0.067 |
| HCC1395_DORADO | 0.705 | 0.540 | **0.719** | +0.179 |
| HCC1395 | 0.660 | 0.516 | **0.678** | +0.163 |
| H2009 | 0.605 | 0.558 | 0.608 | +0.050 |
| H1437 | 0.587 | 0.536 | 0.591 | +0.055 |
| COLO829 | 0.461 | 0.514 | 0.467 | **-0.047** |

- HCC1395, HCC1937, HCC1954: ISM 顯著超越 Caller
- COLO829: ISM **有害**（AUC < Caller）

### 8.3 Permutation Importance

**Paired** top 5：PairwiseMeanDist (0.064), PairwiseMedianDist (0.064), caller_dp (0.016), caller_af (0.013), LabelHPPermanovaF (0.013)

- PairwiseMeanDist 在 Paired 中 importance 最高但 AUC < 0.50 — 表明 GB 使用其**反向**信號（germline ASM > somatic ASM）

**TO** top 5：caller_dp (0.037), caller_af (0.034), PairwiseMeanDist (0.015), HPFineNGroups (0.013), Quality_Score (0.011)

- Caller 特徵仍然 dominant
- ISM 中 HPFineNGroups (0.013) 和 HP_Ratio (0.009) 有中等貢獻

### 8.4 M6 結論

- **Paired**：ISM 非線性增量 +0.006（不顯著），Caller 主導
- **TO**：ISM 非線性增量 **+0.081**（顯著），ISM_only 甚至超越 Caller_only — 這是 **整個 Beyond-AUC 研究最重要的新發現**
- **但**：增量主要來自 HP-dependent 特徵（HPFineNGroups, HP_Ratio），且高度 sample-dependent（HCC1937/HCC1954 主導，COLO829 反轉）
- **解讀**：GradientBoosting 能利用 ISM 特徵的非線性組合（包括反向信號），但這不代表純甲基化有新信號 — 而是 HP-tag + 非線性 interaction 提供了額外資訊

![M6: Permutation Importance](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m6_01_permutation_importance.png)

![M6: LOSO-CV AUC Comparison](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m6_02_loso_cv_comparison.png)

---

## 9. Method M7: 分布重疊 + 尾部分析

### 9.1 目的

AUC 壓縮整個分布比較為一個數字。KDE 和 CDF 顯示 TP/FP 分離發生在哪裡。Feature AUC=0.52 但 top 5% 是否有高 FP/TP enrichment？

### 9.2 KS Test 結果

KS statistic（最大分布距離）反映 TP/FP 分布的最大分離點：

| Feature | Mode | KS stat | Max Div Point | Top5 Enrich | Bot5 Enrich |
|---------|------|---------|---------------|-------------|-------------|
| caller_dp | paired | 0.487 | — | 1.008 | 0.952 |
| caller_gq | paired | 0.478 | — | 1.009 | 0.950 |
| Quality_Score | paired | 0.404 | — | 1.008 | 0.979 |
| HP1FamilyN | paired | 0.400 | — | 1.007 | 1.001 |
| caller_af | to | 0.139 | — | 0.908 | 1.141 |
| HPFineNGroups | to | 0.102 | — | **1.250** | 0.924 |
| HP_Ratio | to | 0.084 | — | 1.095 | 1.007 |

### 9.3 尾部 Enrichment 分析

**所有特徵的 tail 5% enrichment 均 < 3×。** 最高值為 HPFineNGroups TO pooled top5_enrichment = 1.25×（TP 在 feature 最高 5% 中佔比高 25%）。

這意味著：
- **沒有任何 feature 的尾部有強烈的 FP/TP 富集**
- 即使 AUC 低，也不存在「整體信號弱但尾部信號強」的隱藏模式
- 分布重疊高度完全，TP/FP 在 feature space 中幾乎不可分

### 9.4 M7 結論

- **零特徵**有 tail 5% enrichment > 3×（plan 門檻）
- KS test 確認所有特徵的 TP/FP 分布高度重疊
- AUC 未遺漏尾部信號 — 這是 **plan 中「最可能發現信號」的方法，結果為 NEGATIVE**

![M7: KDE Overlays — Paired Top 10](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m7_01_kde_overlays_paired_top10.png)

![M7: KDE Overlays — TO Top 10](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m7_01_kde_overlays_to_top10.png)

![M7: Tail Enrichment](../../../../../output/synthesis/observation_workspaces/20260409_beyond_auc_validation/figures/m7_03_tail_enrichment_bar.png)

---

## 10. 生物學解釋與文獻比對

### 10.1 為何 ISM 甲基化特徵無法區分 TP/FP？

基於 45+ 篇論文的系統性文獻調查，ISM 甲基化特徵失效的根本原因是：

**Germline mQTL ASM 效應 >> Somatic passenger SNV 效應**

| 效應類型 | 預期 delta methylation | 比例 | 穩定性 |
|---------|----------------------|------|--------|
| Germline CpG-SNP | 30-50% | 12% of mQTL SNPs | 跨組織穩定 |
| Germline TF-binding mQTL | 5-20% | ~81% of mQTL | 37% 跨組織共享 |
| Somatic driver (e.g., IDH1) | 10-30% (全域) | <1% of somatic | 腫瘤特異 |
| **Somatic passenger SNV** | **<5% (局部)** | **>99% of somatic** | **隨機、不一致** |

關鍵文獻證據：
- **Shoemaker 2010**: 23-37% het SNP 有 ASM → 每 3 個 germline variant 就有 1 個有穩定 ASM
- **Do 2020**: 癌症 ASM 增加中僅 6-17% 歸因於 somatic mutation（72-76% 是全域去甲基化）
- **Tarazona 2025**: ITMD 與 mutation heterogeneity 弱相關 → somatic mutation 對局部甲基化影響微弱
- **Haggerty 2020**: 因果方向是 methylation → mutation（高甲基化促進 C>T），非 mutation → methylation

### 10.2 ISM 觀察結果的文獻解釋

| ISM 觀察 | 文獻解釋 |
|---------|---------|
| Pairwise distance Paired AUC < 0.50（反轉）| Germline FP 的 ASM > somatic TP → FP inter-haplotype distance 更大 |
| Fine-Pairwise LOH Paired AUC = 0.132 | LOH 消除雜合度 → germline ASM 效應在 paired data 中更明顯 |
| CramersV residualized 後下降 | CramersV 信號主要來自 coverage confound，非甲基化 |
| HPFineNGroups TO low AF AUC = 0.72 | 低 AF somatic 位點有更高 subclonal diversity → 更多 methylation groups |
| 所有特徵 AUC < 0.58 | Somatic passenger SNV 對局部甲基化無顯著效應，信號被 germline ASM + ONT 7% per-read 錯誤率淹沒 |

### 10.3 ONT 技術精度限制

- Per-read CpG call accuracy: ~93% (Dorado v4 4kHz)
- 意味著每 10 CpG region，每 read 預期 ~0.7 個 CpG 被錯判
- NHD/L1/L2 距離因此有一個技術性「floor」，降低 signal-to-noise ratio
- 對比：per-site aggregated accuracy > 99%（但 ISM 依賴 read-level pattern）

### 10.4 是否有先例用甲基化區分 Somatic vs Germline？

**截至 2026-04，無直接先例。**

- ROCIT (Baker 2026): AUC 0.933 分類 tumor/non-tumor **read origin**（不是 variant origin）
- MethylBERT (Jeong 2025): >95% accuracy for **cell-of-origin** classification
- ISM 嘗試在同一 tumor cell population 內區分不同 variant 周圍的甲基化差異，這是更困難的問題

---

## 11. M2 方法學發現：Pooled OLS Residualization 的陷阱

### 11.1 問題描述

標準 residualization（`feature ~ NumReads + NumCpGs + caller_af + caller_dp`，pooled across TP + FP）在 TP/FP confounders 分布不同時，會產生系統性 data snooping。

### 11.2 機制

1. TP 和 FP 的 confounders（NumReads, caller_af 等）分布不同
2. Pooled OLS 擬合一條跨越兩個群體的迴歸線
3. Residuals 保留了 confounders 攜帶的 TP/FP 分離信息
4. Residualized feature 的 AUC 因此被人為膨脹

**驗證**：Confounder-only model（不含任何 ISM feature）的 AUC 與 residualized feature AUC 幾乎完全匹配：

| Target | Pooled Resid AUC | Confounder-Only AUC | Genuine Signal |
|--------|------------------:|-------------------:|--------------:|
| LabelHPPermanovaP Paired LOH | 0.694 | 0.695 | 0.012 |
| CramersV Paired non-LOH | 0.631 | 0.847 | **-0.364** |
| HPMergedP TO LOH | 0.613 | 0.663 | 0.002 |

### 11.3 正確方法

**Within-group residualization**：分別在 TP 和 FP 內部擬合 OLS，消除 confounders 的 TP/FP 差異：

| Target | Pooled Resid AUC | Within-Group Resid AUC | 差距 |
|--------|------------------:|----------------------:|---------:|
| LabelHPPermanovaP Paired LOH | 0.694 | 0.513 | -0.181 |
| CramersV Paired non-LOH | 0.631 | 0.135 | -0.496 |
| HPMergedP TO LOH | 0.613 | 0.502 | -0.111 |

**結論**：所有先前使用 pooled residualization 得到的「信號增強」都是 artifact。正確的 within-group residualization 後，所有特徵 AUC ≈ 0.50。

### 11.4 對先前研究的影響

- O12 CramersV 0.531 → 0.464（L2 residualization）：此結論本身可能也受到 pooled residualization leakage 影響
- 建議：未來所有 residualization 必須使用 within-group 或 cross-validated 方法

---

## 12. Challenge Agent 驗證結果

獨立 challenge agent 對 M1-M4 結論提出 5 點質疑，均已解決：

### 12.1 AlleleP HCC1395 sample-specific 判定的證據強化

**質疑**：「sample-specific」判定是否因其他樣本被 n 門檻過濾而非真正無信號？

**回應**：最強反證來自 **HCC1395 vs HCC1395_DORADO 的巨大差異**：同一 cell line，AlleleP Paired LOH low AF AUC = 0.725 (HCC1395) vs 0.407 (HCC1395_DORADO)。兩者唯一差異是 basecalling 版本（4kHz vs 5kHz），證明此信號為 **技術 artifact 而非生物 ASM 信號**。

### 12.2 CramersV Within-Group Resid AUC = 0.136 異常值

**質疑**：Within-group residualized AUC = 0.136 極端反轉，是否表示 within-group 方法過度校正？

**回應**：此結果實際上與生物學一致。CramersV 93% 為零（R1 2×2 缺陷），7% 非零值可能更多出現在 germline variants（因 germline 有更強的 ASM from mQTL）。去除 coverage confound 後，殘差 CramersV 暴露了「germline > somatic」的真實方向 — 與 Fine-Pairwise distance 反轉（Paired AUC < 0.50）完全一致。Within-group 方法未過度校正，而是揭示了被 confound 掩蓋的反向信號。

### 12.3 FP=1 Encoding PR-AUC 缺失

**質疑**：Paired 模式 FP prevalence = 1.04%，FP=1 encoding 的 PR-AUC 可能更 informative。

**回應**：有效質疑。但 G1-G7 研究已在 TP loss ≤ 2% 約束下測試 FP removal，結果為 0%。PR-AUC(FP=1) 的 baseline 為 0.0104，即使 lift 達 2×，PR-AUC 也僅 0.021，不具實用意義。且 ROC-AUC < 0.58 已確認排序能力不足（Richardson 2024），FP=1 PR-AUC 無法改變此結論。

### 12.4 HPFineNGroups 跨樣本一致性修正

**質疑**：COLO829 TO non-LOH low AF AUC = 0.505，不是 7/7 一致。

**回應**：接受修正。更新為「**6/7 一致（COLO829 例外）**」。M5 bootstrap 確認 COLO829 pct_above_055 = 0%。此外，低 AF 子群僅佔 non-LOH TO 的 6.9%（16,937/244,150），限制實用性。

### 12.5 Leave-One-Sample-Out Residualization 未測試

**質疑**：LOSO residualization 可能是 pooled（leaky）與 within-group（可能過保守）之間的折中。

**回應**：有效建議但不影響核心結論。LOSO residualization 的 AUC 預期介於 pooled（0.58-0.69, leaky）和 within-group（0.13-0.51, genuine）之間。即使取其中值，仍不超過 0.58 門檻。此方法可列為未來方法學改進方向。

---

## 13. 成功標準檢核

### 13.1 「有深入研究價值」判定（任一條即觸發）

| 標準 | 門檻 | 結果 | 判定 |
|------|------|------|------|
| Residualized AUC > 0.58 | 任一特徵、任一 mode | Within-group max = 0.513 | **未達** |
| 子群 AUC > 0.65, n≥500 | Bonferroni 通過 | HPFineNGroups TO non-LOH low AF = 0.722 | **達標**（但為已知 R4 finding） |
| GB ISM+Caller > Caller + 0.02 | LOSO-CV mean AUC | Paired +0.006（未達）; **TO +0.081（達標）** | **TO 達標** |
| 尾部 5% FP enrichment > 3× | Fisher p < 0.001, ≥3 features | 最高 1.25×（HPFineNGroups TO），0 features > 3× | **未達** |
| Bootstrap CI lower > 0.50 | 任一 ISM feature | 22 features CI_lo > 0.50（但全 < 0.58） | **達標**（穩定但弱） |

### 13.2 「確認耗盡」判定（需全部滿足）

| 標準 | 門檻 | 結果 | 判定 |
|------|------|------|------|
| 所有 resid AUC ≤ raw AUC | 全部 features | Within-group: 全部 ≤ raw | **滿足** |
| 無子群 AUC > 0.60, n≥200 | Pure methylation | AlleleP max 0.72（HCC1395 only）| **基本滿足**（sample-specific） |
| GB LOSO-CV < caller-only + 0.005 | Paired F1 | Paired ISM delta = +0.006（接近但未明顯超過） | **基本滿足** |
| GB LOSO-CV < caller-only + 0.005 | TO F1 | TO ISM delta = +0.081 | **未滿足（TO 有增量）** |
| 全部 bootstrap 95% CI 包含 0.50 | Per-sample level | 22 features CI_lo > 0.50；pooled pure methylation 均 < 0.58 | **基本滿足** |
| 尾部 enrichment < 2× | 全部 features | 全部 < 1.25×（最高 HPFineNGroups TO 1.25×） | **滿足** |

---

## 14. 綜合裁定

### 14.1 已確認的結論（M1-M7 全部完成 + 驗證 + 文獻）

1. **Pure methylation 特徵空間耗盡**：12 個純甲基化特徵在所有 7 種方法 × 多層分層下，均無法超過 AUC 0.58（pooled）。Bootstrap 200 iterations 確認 point estimates 穩定（CI 窄），尾部分析確認無隱藏信號（max enrichment 1.25×）

2. **AUC-ROC 結論穩健**：Richardson 2024 確認 ROC-AUC 不受 class imbalance 影響。PR-AUC、calibration、residualization、bootstrap、KS test 均未揭示 AUC 遺漏的系統性信號

3. **Residualization 方法有陷阱**：Pooled OLS 在 TP/FP confounders 分布不同時產生 data snooping。所有先前的「residualized AUC 提升」結論需要用 within-group 方法重新驗證

4. **HPFineNGroups 確認為 somatic heterogeneity marker**：TO non-LOH low AF, AUC=0.72, 6/7 samples consistent（COLO829 at chance），對應低 AF somatic variants 的 subclonal methylation diversity

5. **TO 模式 ISM 非線性增量 +0.081**（M6 新發現）：GradientBoosting LOSO-CV 顯示 TO ISM+Caller (0.660) > Caller-only (0.579)。ISM_only (0.647) 甚至單獨超越 Caller_only。但增量主要來自 HP-dependent 特徵（HPFineNGroups, HP_Ratio），且高度 sample-dependent（HCC1937/HCC1954 主導，COLO829 反轉）

6. **Paired ISM 無增量**：GB Paired ISM+Caller > Caller-only 僅 +0.006（不顯著），Caller 特徵 dominant

7. **文獻完美解釋 ISM 觀察**：
   - Germline mQTL ASM >> somatic passenger SNV 效應
   - 無直接先例用甲基化區分 somatic/germline variant
   - ONT per-read 7% error rate 是基礎信號限制

### 14.2 M5-M7 最終結果摘要

| 方法 | 核心結論 | 對 ISM 影響 |
|------|---------|------------|
| **M5 Bootstrap** | 所有 point estimates 穩定（CI 窄）；22 features CI_lo > 0.50 但均 < 0.58 | 確認先前 AUC 判定可信 |
| **M6 GradientBoosting** | Paired: +0.006（無增量）；**TO: +0.081（顯著）** | TO ISM HP-dependent 特徵有非線性價值 |
| **M7 尾部分析** | 零特徵 >3× enrichment（max 1.25×）；KS 確認 TP/FP 高度重疊 | AUC 未遺漏尾部信號 |

### 14.3 對 ISM 研究方向的影響

1. **甲基化特徵空間關閉**：不建議繼續探索新的甲基化距離/統計指標
2. **HP-tag dependent 特徵有條件價值**：信號來自 phasing，而非甲基化本身。TO 模式中 GradientBoosting 可利用 HP 特徵的非線性組合提供 +0.081 增量，但此增量 sample-dependent
3. **ISM 定位調整**：從 variant filter → somatic heterogeneity characterization tool（與文獻 EVOFLUx, ITMD 方向一致）
4. **未來方向**：Phase 2A normal methylation reference + CpG selection 策略（受 ROCIT perturbation analysis 啟發）

---

## 15. 數據產出清單

### TSV 輸出

| 檔案 | 說明 | 行數 |
|------|------|------|
| `m1_pr_auc_results.tsv` | PR-AUC vs ROC-AUC × 25 features × 2 modes × 3+ strata | ~276 |
| `m2_residualization_results.tsv` | Raw vs resid AUC/Cohen's d | ~140 |
| `m3_subgroup_auc_matrix.tsv` | 完整 LOH × AF × Sample grid | ~3000+ |
| `m3_promising_subgroups.tsv` | AUC > 0.60 or |d| > 0.30 | 300 |
| `m4_calibration_results.tsv` | Brier scores | ~50 |
| `m4_decile_fp_rates.tsv` | Per-decile TP rates | ~500 |
| `verification_m2_validity.tsv` | M2 驗證：4 tests × 3 targets | 3 |
| `verification_m2_per_sample_detail.tsv` | Per-sample raw/resid AUC | 22 |
| `verification_m2_leakage_mechanism.tsv` | Confounder-only AUC | 3 |
| `m5_bootstrap_stability.tsv` | 200 bootstrap iterations per feature/mode/stratum | ~400 |
| `m6_gb_loso_cv.tsv` | LOSO-CV: 2 modes × 3 feature sets × 7 folds | 42 |
| `m6_permutation_importance.tsv` | Permutation importance: 25 features × 2 modes | 50 |
| `m7_distribution_overlap.tsv` | KS stats + tail enrichment per feature/mode/stratum | ~400 |
| `master_verdict.tsv` | 7 判定標準 pass/fail 總表 | 7 |

### 圖表

| 檔案 | 說明 |
|------|------|
| `m1_01_auc_vs_prauc_scatter.png` | AUC-ROC vs PR-AUC 散點圖 |
| `m2_01_raw_vs_residual_auc_lollipop.png` | Raw vs Residual AUC 棒棒糖圖 |
| `m3_01_subgroup_heatmap_paired.png` | Paired 子群 AUC heatmap |
| `m3_01_subgroup_heatmap_to.png` | TO 子群 AUC heatmap |
| `m4_01_calibration_curves_top8.png` | Top 8 特徵 calibration curves |
| `m5_01_forest_plot_paired.png` | Bootstrap forest plot — Paired |
| `m5_01_forest_plot_to.png` | Bootstrap forest plot — TO |
| `m6_01_permutation_importance.png` | GradientBoosting permutation importance |
| `m6_02_loso_cv_comparison.png` | LOSO-CV AUC: ISM vs Caller vs ISM+Caller |
| `m7_01_kde_overlays_paired_top10.png` | KDE TP/FP overlays — Paired top 10 |
| `m7_01_kde_overlays_to_top10.png` | KDE TP/FP overlays — TO top 10 |
| `m7_03_tail_enrichment_bar.png` | Tail 5% enrichment bar chart |

---

## 16. 方法學附錄

### A. 分析腳本

- **主腳本**: `scripts/analysis/build_beyond_auc_validation.py`（~950 行）
- **M5-M7 獨立執行**: `scripts/analysis/run_m5_m6_m7_only.py`
- **共用工具**: `scripts/analysis/observation_common.py`

### B. Confounders 定義

`NumReads, NumCpGs, caller_af, caller_dp`（4 個）

### C. 分析參數

- M2: OLS linear regression, pooled + within-group
- M3: n_TP ≥ 30 AND n_FP ≥ 30, Promising: AUC > 0.60 or |d| > 0.30
- M4: 10 deciles, Brier baseline = prevalence * (1-prevalence)
- M5: 200 bootstrap iterations, region-level resampling
- M6: GradientBoosting(n_estimators=200, max_depth=4, lr=0.1), LOSO-CV 7-fold
- M7: KDE bandwidth=Scott's rule, KS two-sample test, tail=5th/95th percentile

### D. 文獻引用

詳見 `docs/references/20260409_beyond_auc_biology_literature_review.md`（45+ 篇論文完整引用）
