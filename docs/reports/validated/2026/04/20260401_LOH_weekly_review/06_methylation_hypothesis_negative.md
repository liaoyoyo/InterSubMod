<!--
建立時間: 2026-04-01 10:00
目標: 記錄 O11/O12 甲基化假說的完整否決過程，包含方法學教訓（L2 Collider Bias）
處理範圍: O11 Within-Group Methylation Heterogeneity + O12 TO LOH Methylation Scenario Discrimination
關聯檔案:
  - observation_workspaces/20260331_O11_methylation_heterogeneity/
  - observation_workspaces/20260331_O12_loh_methylation_scenarios/
  - docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md
-->

# 06 — 甲基化假說否決報告：O11 Heterogeneity + O12 LOH Scenario Discrimination

> **結論前置**：兩個方向均為 **NEGATIVE**。TO 模式下甲基化特徵正式排除用於 germline/somatic 區分。唯一剩餘路徑 = Phase 2 Normal Methylation Reference。

**為什麼 negative result 也是重要的？** 在系統性研究中，明確排除不可行的方向與發現可行方向具有同等價值。O11/O12 的否決：(1) 防止未來研究者重複投入已證偽的方向；(2) 揭示了 n_reads confound 和 L2 collider bias 兩個方法學陷阱，這些教訓適用於所有涉及 residualization 的生物資訊分析；(3) 將研究資源重新聚焦到唯一剩餘的可行路徑上。

---

## 一、O11: Within-Group Methylation Heterogeneity

### 1.1 假說

| 項目 | 內容 |
|------|------|
| **原始假說** | Germline ASM（Allele-Specific Methylation，由 mQTL 驅動）在同一 HP（haplotype）組內的 reads 比 cancer ASM 更一致（更低 heterogeneity） |
| **預期** | TP（somatic variant with cancer methylation）有更高 within-group variance，heterogeneity metrics 可區分 TP/FP |
| **理論基礎** | 如果 germline 的甲基化是穩定的 ASM pattern（受 cis-regulatory element 控制），同一 haplotype 內的 reads 應有高度一致的甲基化狀態；但 cancer 甲基化是隨機/不穩定的（stochastic），因此同組 reads 之間的 heterogeneity 應較高 |
| **數據規模** | 561 regions（TP=337, FP=224），來自 paired multibio 637 manifest，7 個樣本 |
| **來源** | Paired-mode data，min_group_reads=5，methyl_threshold=0.5 |

### 1.2 方法

#### 六個 Heterogeneity 特徵

| # | 特徵名稱 | 計算方式 | 物理意義 |
|---|---------|---------|---------|
| 1 | **epipolymorphism** | 一個 region 內所有 reads 的甲基化 pattern（binary CpG vector）的 Shannon entropy；pattern 越多樣，entropy 越高 | 衡量一個位點的表觀遺傳多樣性——多少種不同的甲基化「型態」共存 |
| 2 | **per_cpg_concordance** | 每個 CpG 位點在同組 reads 中一致性的平均值；計算每個 CpG 的 majority allele fraction（最多數狀態的 read 比例），取所有 CpG 的 mean | 衡量同 HP 組內各 CpG 的方向一致程度；高值 = 一致，低值 = 混亂 |
| 3 | **per_cpg_entropy** | 每個 CpG 位點在同組 reads 中二元甲基化狀態的 Shannon entropy，取所有 CpG 的 mean | 與 concordance 互補，衡量每個 CpG 的不確定性 |
| 4 | **within_dist_mean** | 同 HP 組內所有 read pairs 的 pairwise distance（Hamming distance on methylation vectors）的 mean | 直接衡量同組 reads 之間的整體差異程度 |
| 5 | **between_within_ratio** | (between-group mean distance) / (total mean distance)，類似 PERMANOVA 的 R-squared | 衡量 HP 之間的差異相對於總體差異的比例；高值 = HP 間差異顯著 |
| 6 | **within_dist_var_mean** | 各 HP 組內 pairwise distance 的 variance 的 mean | 衡量同組 reads 之間差異的「穩定性」；高 variance = 組內 reads 的一致程度不穩定 |

#### Residualization（殘差化）流程

**為什麼需要**：n_reads（read 數量）是嚴重的 confound。更多 reads 機械性地產生更多 unique methylation patterns，導致 epipolymorphism 自然偏高。在本資料集中，TP regions 的 median n_reads=157（FP=92），n_reads AUC(TP)=0.926，即 TP 和 FP 之間的 read 數差異極大。

**方法**：
1. 以 OLS regression 建模 `feature ~ n_reads + num_cpgs`
2. 取 residual（殘差）作為 confound-corrected feature
3. 對 residualized feature 重新計算 AUC
4. 額外以 read-count-matched bin（n_reads 在 [81, 120] 範圍內，TP=42, FP=117）驗證

### 1.3 關鍵數據

#### Raw vs Residualized AUC 對比

| 特徵 | Raw AUC | Residualized AUC | Cohen's d | 判定 |
|------|---------|------------------|-----------|------|
| epipolymorphism | **0.845** | **0.530** | ~0 | 完全被 confound 驅動 |
| per_cpg_concordance | 0.642 | 0.523 | ~0 | confounded |
| per_cpg_entropy | 0.633 | 0.513 | ~0 | confounded |
| within_dist_mean | 0.623 | 0.509 | ~0 | confounded |
| between_within_ratio | 0.601 | 0.536 | ~0 | confounded |
| within_dist_var_mean | 0.566 | **0.578** | ~0 | 最弱 confound，但仍 < 0.60 |

**Confound 證據鏈**：
- n_reads AUC(TP) = 0.926（TP regions 有 1.87x 更多 reads）
- n_hp1 AUC(TP) = 0.830（TP 的 HP1 reads 是 FP 的 2.9x）
- Epipolymorphism 和 n_reads 的 Spearman r = **0.79**（機械性 artifact）
- Read-count-matched bin [81-120] 中 epipolymorphism AUC = **0.560**（確認殘差化結論）

**為什麼 TP 有更多 reads？** 這是 paired manifest 的 **sample selection artifact**：TP regions 以 true somatic variants 為中心，這些位點有高信心的 variant calling 結果，因此傾向於高覆蓋率區域；FP regions 包含低品質的 germline variants 誤判為 somatic，覆蓋率較不一致。此覆蓋率差異在 TO mode 中不存在（TP 和 FP 來自同一 caller）。

#### 特徵間相關性

epipolymorphism 與 within_dist_mean/per_cpg_entropy 等特徵的高相關性（per_cpg_concordance 和 per_cpg_entropy Spearman r = -0.99）表明這些特徵捕捉的是同一個信號維度，因此其 residualized AUC 的集體崩潰是一致的。

### 1.4 圖片清單

以下圖片位於 `observation_workspaces/20260331_O11_methylation_heterogeneity/`：

| 檔案 | 內容 |
|------|------|
| `fig01_within_dist_variance.png` | Within-group distance variance：TP vs FP 分佈 |
| `fig02_per_cpg_concordance.png` | Per-CpG concordance：TP vs FP 分佈 |
| **`fig03_epipolymorphism.png`** | **[關鍵圖] Epipolymorphism raw vs residualized 對比，展示 n_reads confound 的效果** |
| **`fig04_feature_auc_ranking.png`** | **[關鍵圖] 6 個特徵的 Raw AUC vs Residualized AUC ranking，一目了然的信號崩潰** |
| `fig05_cpg_density_stratified.png` | CpG density 分層後的特徵表現 |
| `fig06_feature_correlation.png` | 特徵間相關矩陣 |
| **`fig07_cross_feature_scatter.png`** | **[關鍵圖] Epipolymorphism vs n_reads 的 scatter plot，展示 r=0.79 的線性關係** |
| **`fig08_combined_roc.png`** | **[關鍵圖] 所有特徵的 Combined ROC 曲線（residualized），展示所有曲線趨近 diagonal** |

### 1.5 O11 結論

**Within-group methylation heterogeneity（concordance、entropy、epipolymorphism、within-group distance variance）無法區分 germline 和 somatic variants。**

mQTL 驅動的「穩定 germline ASM」vs「隨機 cancer ASM」假說在 ISM 的解析度下（10kb window、~50-200 reads）不被數據支持。Germline 和 somatic variants 在控制 read 數量後，展示出相似的 within-group 甲基化模式。

**在當前線性特徵框架與 ISM 解析度下，此方向正式排除。**

---

## 二、O12: TO LOH Methylation Scenario Discrimination

### 2.1 假說

TO（Tumor-Only）LOH（Loss of Heterozygosity）位點存在三種可能的生物學場景：

| 場景 | 事件順序 | 甲基化預測 | Truth Label |
|------|---------|-----------|-------------|
| **S1: LOH → sSNV** | 先發生 LOH，然後在 LOH 區域上發生體細胞變異 | Uniform（僅單一 allele 的甲基化模式殘留） | TP |
| **S2: sSNV → LOH** | 先發生體細胞變異，然後 LOH 消除了 wild-type allele | 可能 heterogeneous（如果 LOH 發生前有 biallelic methylation） | TP |
| **S3: Germline LOH** | 看起來像 LOH 但實際是 germline variant 被 TO caller 誤判為 somatic | Uniform（mQTL 穩定） | FP |

**核心挑戰**：S1 和 S3 預測相同的甲基化模式（均為 uniform），理論上僅 S2 可能有不同的甲基化特徵。

**預期**：如果三場景的甲基化模式有系統性差異，ISM 特徵可用於在 TO mode 中區分 LOH 區域的 TP 和 FP。

**數據規模**：175,542 TO LOH rows（TP=129,589, FP=45,953） x 22 features x 7 samples x 4-level confound control

### 2.2 方法：四層 Confound 控制

| 層級 | 名稱 | 方法 | 目的 |
|------|------|------|------|
| **L0** | Raw AUC | 直接在所有資料上計算 feature AUC | 基準線——可能包含 confound |
| **L1** | 控制 n_reads | Residualize feature on n_reads | 移除 read 數量帶來的機械性 artifact |
| **L2** | 控制 caller_af | OLS regression `feature ~ caller_af`，取 residual 計算 AUC | 移除 AF（allele frequency）confound |
| **L3** | AF-bin 內 AUC | 將 AF 分成 4 個 bins（[0, 0.3], [0.3, 0.6], [0.6, 0.9], [>0.95]），在每個 bin 內直接計算 AUC | 最可靠的 confound 控制——完全隔離 AF 影響 |

**L3 比 L2 更可靠的原因**：L2（OLS residualization）假設 feature 和 confound 之間是線性關係，且 feature 本身有足夠的 variance。當 feature 近乎常數時，residualization 會引入 artifact（見下節 collider bias）。L3 直接在 AF-matched subsets 中計算，不受此限制。

### 2.3 Step 1: 既有特徵的 Per-Sample AUC

| 樣本 | n_TP | n_FP | FP rate | Best Feature (L0) | AUC L0 | max L3 | 判定 |
|------|------|------|---------|-------------------|--------|--------|------|
| HCC1395 | 17,269 | 6,303 | 26.7% | AlleleDelta | 0.608 | 0.557 | Borderline |
| HCC1395_DORADO | 17,764 | 6,477 | 26.7% | AlleleDelta | 0.644 | 0.568 | Borderline |
| COLO829 | 11,697 | 5,922 | 33.6% | PairwiseMeanDist | 0.569 | 0.579 | Borderline |
| H1437 | 18,988 | 5,160 | 21.4% | AlleleDelta | 0.577 | 0.594 | Signal |
| H2009 | 51,429 | 4,528 | 8.1% | PairwiseMeanDist | 0.569 | 0.597 | Signal |
| HCC1937 | 8,181 | 6,880 | 45.7% | AlleleDelta | 0.660 | 0.586 | Signal |
| HCC1954 | 4,261 | 10,683 | 71.5% | PairwiseMeanDist | 0.551 | 0.553 | Borderline |

**AF Confound 分析**：
- AlleleDelta 和 caller_af 的 Spearman r 範圍：-0.23 ~ -0.57（各樣本不同）
- Residualize on AF 後，所有 per-sample AlleleDelta AUC 降至 < 0.55
- AF-bin 內（完全控制 AF）：AlleleDelta AUC 全部 < 0.59

**結論**：AlleleDelta 的表面信號**完全由 AF confound 驅動**。在 AF 被控制的條件下，殘餘的區分能力不足 0.59。

### 2.4 L2 Collider Bias 發現（重要方法學發現）

在 L2 分析中，出現了一個令人驚訝的現象：部分特徵在 residualize on AF 後，AUC 不降反升，且飆升幅度極大。

#### 被偵測到的 Collider Bias 案例（部分）

| 特徵 | 樣本 | AUC L2 | max L3 | 膨脹量 |
|------|------|--------|--------|--------|
| CramersV | HCC1937 | **0.801** | 0.506 | **+0.295** |
| HPMergedDelta | HCC1937 | **0.777** | 0.520 | **+0.257** |
| CramersV | HCC1954 | **0.762** | 0.509 | **+0.253** |
| HPMergedDelta | HCC1954 | **0.738** | 0.510 | **+0.228** |
| HPMergedSig | HCC1954 | **0.735** | 0.502 | **+0.232** |
| HeuristicScore | HCC1937 | **0.734** | 0.531 | **+0.203** |
| CramersV | HCC1395 | **0.727** | 0.515 | **+0.212** |
| PassedGating | HCC1937 | **0.715** | 0.517 | **+0.198** |

#### 為什麼會發生 Collider Bias？

**前提條件**：
1. Feature Y 在 TO 模式下**近乎常數**（例如 CramersV 在 TO 中無 HP 分化，值幾乎相同）
2. Confound X（caller_af）與 outcome（TP/FP）之間有相關性

**數學直覺**：
假設 CramersV 在某 sample 中的值分佈為 Y ~ 0.05 +/- 0.01（近乎常數），而 caller_af 的值分佈為 X ~ 0.50 +/- 0.20。

當我們做 OLS regression `Y ~ X` 時：
- 因為 Y 幾乎沒有 variance，OLS fit 會找到一條幾乎水平的線，但殘差的微小 pattern 會 **inherit X 的信號方向**
- Residual = Y - beta * X。由於 Y ≈ constant，residual ≈ constant - beta * X = -beta * X + noise
- 如果 X（AF）與 outcome 有相關，那麼 residual 也自動與 outcome 相關——**即使 Y 本身與 outcome 完全無關**

**具體數字說明**：
- CramersV (HCC1937)：L2 AUC = 0.801，看起來是非常強的信號
- 但 AF-bin = [0.3, 0.6]（控制 AF 後）：AUC = 0.501（完全隨機）
- AF-bin = [0.6, 0.9]：AUC = 0.505（完全隨機）
- L3 max AUC = 0.506（全部 AF bins 中最高）
- **膨脹量 = 0.801 - 0.506 = 0.295**——完全是 artifact

**受影響的特徵**：CramersV, HeuristicScore, HPMergedDelta, PassedGating, HPMergedSig——這 5 個特徵在 TO mode 中都是近常數或 binary（因為 TO 無 HP 資訊，與 haplotype 相關的統計量無法計算，回傳預設值）。

#### L2 Collider Bias 的教訓

| 面向 | 內容 |
|------|------|
| **何時出現** | 對近乎常數（low variance）的特徵做 OLS residualization，且 confound 與 outcome 有相關時 |
| **如何偵測** | 計算 `inflation = AUC_L2 - max(AUC_L3)`；若 > 0.10，flag 為 collider bias |
| **如何預防** | 所有 L2 結果必須用 L3（AF-bin stratification）交叉驗證 |
| **保守估計** | `auc_decision = min(auc_L2, max_L3 + 0.03)` |
| **適用範圍** | 不限於本分析——任何使用 OLS residualization 的生物資訊分析都可能遇到此問題 |

### 2.5 Step 2: Novel Feature AUC

除了既有的 9 個 ISM 特徵外，O12 額外計算了 13 個新特徵，涵蓋 region-level 甲基化統計和空間特徵。

| 類別 | 特徵 | 定義 |
|------|------|------|
| Region 統計 | region_methyl_mean | Region 內所有 CpG 的平均甲基化水平 |
| Region 統計 | region_methyl_std | Region 內甲基化水平的標準差 |
| Region 統計 | region_low_methyl_fraction | 低甲基化（<0.3）reads 的比例 |
| Region 統計 | region_high_methyl_fraction | 高甲基化（>0.7）reads 的比例 |
| Read 統計 | read_methyl_uniformity | 每個 read 內 CpG 甲基化狀態的一致性 |
| 全域距離 | global_mean_dist | 所有 read pairs（不分 HP）的平均 distance |
| 全域距離 | global_dist_var | 所有 read pairs distance 的 variance |
| 全域統計 | cpg_concordance_global | 不分 HP 的 per-CpG concordance |
| 全域統計 | cpg_entropy_global | 不分 HP 的 per-CpG entropy |
| 空間特徵 | transition_count | 相鄰 CpG 之間甲基化狀態翻轉的次數 |
| 空間特徵 | cpg_autocorrelation | 相鄰 CpG 甲基化狀態的自相關係數 |
| 空間特徵 | block_length_mean | 連續相同甲基化狀態的 CpG block 平均長度 |
| 控制特徵 | strand_delta | 正負股 reads 甲基化差異 |

#### Per-Sample 最佳 Novel Feature（L2 AUC）

| 樣本 | Best Novel Feature | AUC L2 |
|------|-------------------|--------|
| HCC1395 | cpg_entropy_global | 0.512 |
| HCC1395_DORADO | global_mean_dist | 0.526 |
| COLO829 | region_high_methyl_fraction | **0.637** |
| H1437 | region_low_methyl_fraction | **0.629** |
| H2009 | region_low_methyl_fraction | 0.557 |
| HCC1937 | block_length_mean | 0.539 |
| HCC1954 | region_methyl_mean | 0.575 |

**關鍵觀察**：
- **region_methyl_mean / region_low_methyl_fraction / region_high_methyl_fraction** 在 COLO829, H1437, HCC1954 三個樣本中有 L2 AUC > 0.57，但在其他 4 個樣本中 < 0.56
- **穩定性判定**：需 >=5/7 samples 一致才算「cross-sample stable」。region_methyl_mean 僅 4/7 樣本 AUC > 0.55，2/7 > 0.58 → **不穩定**
- **空間特徵完全失敗**：transition_count（max L2 = 0.583, COLO829 only）、cpg_autocorrelation（max L2 = 0.533）、block_length_mean（max L2 = 0.575）——均 < 0.55 in majority of samples

### 2.6 Step 3: Integration（Collider-Corrected）

在排除 collider-biased 特徵後的跨特徵整合分析：

#### 特徵穩定性排名（corrected AUC, >=5/7 samples above threshold）

| 特徵 | >= 0.55 的 samples 數 | >= 0.58 的 samples 數 | Median AUC | Collider Bias | 判定 |
|------|---------------------|-----------------------|------------|---------------|------|
| region_methyl_mean | 4/7 | 2/7 | 0.556 | No | 不穩定 |
| region_low_methyl_fraction | 4/7 | 2/7 | 0.557 | No | 不穩定 |
| region_high_methyl_fraction | 3/7 | 2/7 | 0.547 | No | 不穩定 |
| HeuristicScore | 4/7 | 2/7 | 0.561 | **Yes** | Artifact |
| PassedGating | 4/7 | 1/7 | 0.557 | **Yes** | Artifact |
| AlleleDelta | 0/7 | 0/7 | 0.526 | No | 無信號 |
| PairwiseMeanDist | 2/7 | 0/7 | 0.504 | No | 無信號 |
| cpg_autocorrelation | 0/7 | 0/7 | 0.518 | No | 無信號 |
| transition_count | 1/7 | 1/7 | 0.521 | No | 無信號 |

**無任何特徵達到「>=5/7 samples, AUC > 0.58」的穩定信號標準。**

### 2.7 Step 4: Cross-Mode Validation

| 指標 | 數值 |
|------|------|
| Total same-locus variants | 459,782 |
| TO LOH same-locus | 175,542 |
| TO-only FP | 45,328 |
| Both-mode FP | 625 |
| Both-mode TP | 127,462 |

Cross-mode 相關性驗證：
- allele_delta Spearman r = 1.00（same locus = same methylation）
- pairwise_median_dist Spearman r = 1.00
- quality_score Spearman r = 0.794（QS 在 TO/paired 中計算方式不同）

### 2.8 圖片清單

以下圖片位於 `observation_workspaces/20260331_O12_loh_methylation_scenarios/`：

| 檔案 | 內容 |
|------|------|
| **`fig01_per_sample_auc_heatmap.png`** | **[關鍵圖] 7 samples x 22 features 的 AUC heatmap（L0/L2/L3 對比），展示特徵信號在跨樣本中的不穩定性** |
| **`fig02_af_confound_demonstration.png`** | **[關鍵圖] AlleleDelta 與 caller_af 的 confound 展示：scatter plot + AF-bin AUC 對比** |
| `fig03_af_distribution.png` | TP vs FP 的 caller_af 分佈 |
| **`fig04_af_bin_auc.png`** | **[關鍵圖] AF-bin 分層後的 AUC 圖表，展示控制 AF 後信號消失** |
| `fig05_novel_feature_ranking.png` | 13 個 novel features 的 AUC ranking |
| `fig08_integration_heatmap.png` | Collider-corrected 整合分析 heatmap |
| `fig09_stability_bar.png` | 特徵跨 sample 穩定性 bar chart |
| **`fig14_decision_summary.png`** | **[關鍵圖] 最終決策圖：所有特徵 corrected AUC < 0.58 的結論展示** |
| `fig_sample_hcc1395.png` | HCC1395 per-sample 詳細圖 |
| `fig_sample_hcc1395_dorado.png` | HCC1395_DORADO per-sample 詳細圖 |
| `fig_sample_colo829.png` | COLO829 per-sample 詳細圖 |
| `fig_sample_h1437.png` | H1437 per-sample 詳細圖 |
| `fig_sample_h2009.png` | H2009 per-sample 詳細圖 |
| `fig_sample_hcc1937.png` | HCC1937 per-sample 詳細圖 |
| `fig_sample_hcc1954.png` | HCC1954 per-sample 詳細圖 |

### 2.9 O12 結論

**甲基化無法區分 TO LOH 區域中的 germline 和 somatic variants。**

三種生物學場景（LOH->sSNV, sSNV->LOH, germline LOH）在 ISM 解析度下，其甲基化模式不可區分。所有既有特徵（9 個）和新特徵（13 個），在 4 層 confound 控制後，corrected AUC 全部 < 0.58。

**在當前線性特徵框架與 ISM 解析度下，此方向正式排除。**

---

## 三、總結論

### 正式排除的方向

| 方向 | 假說 | 結果 | 原因 |
|------|------|------|------|
| O11 | Within-group heterogeneity 可區分 germline/somatic | NEGATIVE | 完全是 n_reads confound（AUC 0.845 → 0.530） |
| O12 | TO LOH 三場景可透過甲基化區分 | NEGATIVE | AF confound + collider bias；corrected AUC 全 < 0.58 |

### 方法學收穫

1. **n_reads confound**：任何 region-level TP/FP 分析必須控制 read count
2. **L2 Collider Bias**：OLS residualization 對近常數特徵會產生虛假信號；所有 L2 結果必須以 L3（AF-bin stratification）交叉驗證
3. **AF-bin stratification (L3)** 是最可靠的 confound 控制方法

### 剩餘路徑

TO 模式下甲基化特徵正式排除用於 germline/somatic 區分。**唯一剩餘路徑 = Phase 2 Normal Methylation Reference**（需要 matched normal BAM 作為 baseline）。

### 一致性確認

- O11 結論與 O1-O10 一致：no single ISM feature AUC > 0.58 (paired, after confound correction)
- O12 結論與 O5 一致：methylation features 在 TO 和 paired 中行為不同
- O12 結論與 O1 一致：no single TO feature AUC > 0.58 (after confound correction)
- L2 collider bias 發現為新的方法學貢獻

---

## 待驗證問題（已驗證 / 已更新）

### 部分解決

2. **region_methyl_mean 在 COLO829/H1437 的弱信號**（L2 AUC ~0.63） — COLO829 部分已解釋：C11 確認 COLO829 為低測序深度 (~30x) outlier，TP/FP AF 分佈高度重疊（KS=0.122 最低），所有特徵在此樣本上判別力減弱是深度效應的延伸。H1437 部分仍需單獨確認 purity。

### 尚未解決（全為 P2 研究方向）

1. **Phase 2 Normal Methylation Reference 的可行性**：需評估 matched normal BAM 數據可用性和 CpG 覆蓋率。
3. **AF-bin [0.3-0.6] 內的殘餘微弱信號**：corrected AUC 全 < 0.58，此方向優先級很低。
4. **O12 的 TP/FP label 品質**：TO truth set 是否有 systematic labeling bias。
5. **L2 Collider Bias 的理論邊界**：什麼樣的 variance 水平才是安全的？需 formal diagnostic。

---

## 認知門檻補充建議

### 給未來研究者的建議

1. **不要被 raw AUC 欺騙**：O11 的 epipolymorphism AUC=0.845 和 O12 的 CramersV L2 AUC=0.801 都看起來很「好」，但都是 artifact。任何未經 confound 控制的 AUC 都不應被信任。

2. **Residualization 不是萬能的**：OLS residualization（L2）在大多數情況下有效，但對近常數特徵會失效。永遠需要一個獨立的驗證方法（如 L3 stratification）。

3. **Negative results 需要嚴格的 confound 控制才有說服力**：如果只報告 raw AUC < 0.60 就宣告 negative，可能遺漏被 confound 壓制的真信號。O11/O12 的價值在於它們系統性地控制了所有已知 confound 後仍然 negative。

4. **生物學直覺 vs 數據現實**：S1/S3 場景預測相同的甲基化模式（uniform）這一理論觀察，在實驗前就暗示了此方向的困難。未來設計實驗時，應先評估理論上的可區分性。

5. **跨 sample 穩定性是關鍵標準**：一個特徵在 2-3 個 sample 中有信號不代表它有效。至少需 5/7（~70%）的 samples 一致顯示 AUC > 0.58 才值得追蹤。

### 推薦閱讀順序

對於首次接觸此報告的讀者，建議依序閱讀：

1. 本報告 Section 1.3（O11 核心數據）→ 理解 confound 的威力
2. 本報告 Section 2.4（L2 Collider Bias）→ 理解方法學陷阱
3. `docs/experiments/INDEX.md` → 了解 O11/O12 在整體研究框架中的位置
4. `docs/references/20260331_甲基化區分germline_somatic_variant文獻調查_01.md` → 了解文獻中相關方法的局限性

---

*本報告基於 O11/O12 observation workspaces 數據撰寫，原始數據與圖表位於：*
- *`observation_workspaces/20260331_O11_methylation_heterogeneity/`*
- *`observation_workspaces/20260331_O12_loh_methylation_scenarios/`*
- *分析腳本：`scripts/analysis/build_observation_O11_methylation_heterogeneity.py`、`scripts/analysis/build_observation_O12_loh_methylation_scenarios.py`*
