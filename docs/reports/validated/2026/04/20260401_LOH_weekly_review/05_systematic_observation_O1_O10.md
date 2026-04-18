<!--
建立時間: 2026-04-01 04:00
目標: 系統性觀察 O1-O10 研究審閱（LOH Weekly Review 第五章節）
處理範圍: 748,391 rows × 116 columns, 9/10 觀察主題, 82 張圖表
關聯檔案:
  - docs/reports/validated/2026/04/20260401_systematic_observation_O1_O10_summary_01.md
  - docs/reports/validated/2026/04/20260401_comprehensive_observation_O1_O10_report_01.md
  - big7_disk_output/synthesis/observation_workspaces/OBSERVATION_INDEX.md
-->

# 05. 系統性觀察 O1-O10 審閱

---

## 觀察方法論

### 為什麼需要系統性觀察

InterSubMod 的 master dataset 包含 **748,391 rows x 116 columns**，涵蓋 7 個癌症樣本、2 種分析模式（paired/TO）、14 個 dataset。如此龐大的數據量無法僅靠單一全域統計理解——特徵間的交互作用、模式間的方向翻轉、樣本間的異質性，都需要逐主題拆解才能浮現。

系統性觀察的核心策略是：**將 116 個欄位依主題分為 9 個觀察單元**（O1-O8, O10），每個觀察聚焦一個維度（全域分佈、QS 組件、LOH 特徵、caller 特徵、甲基化特徵、VerificationClass、TO phasing 品質、樣本異質性、read-level 甲基化），產出獨立的圖表與結論後，再以交叉驗證整合。

### 觀察等級定義

| 等級 | 定義 | 條件 |
|------|------|------|
| **Level A** | 跨觀察交叉驗證通過的結論 | 至少 2 個獨立觀察確認，且零矛盾 |
| **Level B** | 單觀察確認的結論 | 僅 1 個觀察支持，尚需後續驗證 |

### 交叉驗證方法

採用**零矛盾準則**：針對 9 個觀察之間可能存在矛盾的 6 組結論進行逐一檢查。判定標準為數字一致（AUC 四捨五入內吻合）且結論方向一致。最終結果：**6 組全部判定一致，零矛盾**。數字一致性 8/9 完全吻合。

> 交叉驗證詳細報告見 `20260401_cross_validation_report.md`

---

## O1: Master Distribution — 全域數據概覽

**目的**：建立 748K rows x 116 columns 的全域統計基線，識別 paired/TO 兩種模式的根本差異。

**方法**：對 20 個核心數值欄位繪製分佈圖、相關矩陣、缺失值模式，計算每個特徵在 paired/TO 模式下的 AUC，並比較 TP/FP 組成。

**關鍵發現**：

- 資料涵蓋 **748,391 rows, 14 datasets, 7 samples x 2 modes**
- Paired FP rate 僅 **1.04%**（3,429/328,699），TO FP rate 高達 **30.6%**（128,382/419,692）——相差 29 倍
- Paired 最強特徵 caller_gq AUC=0.811，TO 所有特徵 AUC < 0.58
- Quality_Score 在 TO 下 AUC=0.497（等同隨機猜測）
- H2009 佔全域 36.2% rows，樣本間數量差異達 4-10 倍

**代表圖表**：

![fig01 — 20 個核心數值特徵分佈](../figures/20260401_observation_figures/O01/fig01_numeric_feature_distributions.png)

*NumReads 中位數 74，CramersV 零膨脹率 >90%，HP_Ratio 呈多峰結構。*

![fig10 — Paired/TO 各特徵 AUC 排名](../figures/20260401_observation_figures/O01/fig10_top_discriminative_features.png)

*Paired 三強: caller_gq (0.811), NumReads (0.784), Quality_Score (0.754)。TO 全面崩潰，最高僅 0.572。*

![fig11 — Paired vs TO 特徵分佈偏移](../figures/20260401_observation_figures/O01/fig11_paired_vs_to_feature_shifts.png)

*效應量普遍很小 (|r| < 0.12)，但微小偏移足以讓分類器跨模式失效。*

---

## O2: Quality Score Decompose — QS 組件拆解

**目的**：拆解 Quality Score 各組件（LOH penalty, CNV Loss, Low Reads, Dual Sig bonus 等），量化每個組件對 TP/FP 區分的正面或負面貢獻。

**方法**：完美重建 QS（corr=1.0000, MAE=0.0000），逐一打開/關閉各組件，測量 AUC 變化（Leave-One-Component-Out sensitivity analysis）。

> **AUC 簡述**：AUC (Area Under the ROC Curve) 衡量分類器區分正負樣本的能力。0.5 = 隨機猜測（無區分力），1.0 = 完美區分，<0.5 表示方向反轉（預測錯誤方向）。

**關鍵發現**：

- LOH penalty 觸發率：TP=44.5% vs FP=35.8%——**反向懲罰 TP 多於 FP**
- 移除 LOH penalty 後 TO AUC 從 0.497 升至 **0.542**（+0.045，僅 LOH penalty）；實際 commit 同時移除 verify bonus，預估 ~0.546（+0.049，待 benchmark 驗證）
- 即使完全優化，TO QS ceiling 僅 ~0.55，遠低於實用門檻 0.70+
- Paired QS AUC=0.754（Cohen's d=1.018），驅動力為 CNV Loss 和 Low Reads

**代表圖表**：

![fig03 — LOH Penalty 觸發率（TP vs FP）](../figures/20260401_observation_figures/O02/fig03_qs_loh_penalty_trigger_rate.png)

*TO 模式 LOH penalty 不成比例地懲罰 TP，設計假設在 TO 下不成立。*

![fig04 — 移除 LOH Penalty 前後 ROC 比較](../figures/20260401_observation_figures/O02/fig04_qs_without_loh_penalty_roc.png)

*移除後 TO AUC +0.045，但仍遠低於實用水準。*

![fig05 — 各組件個別 AUC 排名](../figures/20260401_observation_figures/O02/fig05_qs_component_auc_ranking.png)

*TO 中 LOH 是最差組件 (AUC=0.457, 反向)，CNV Loss 是唯一正向貢獻者。*

---

## O3: LOH Features Post-Fix — HP fix 後 LOH 特徵全面評估

**目的**：在 HP integer tag fix 完成後，全面重新評估所有 LOH 相關特徵的分佈、enrichment 方向、以及 TO/paired concordance。

**方法**：計算 LOH rate by truth label/mode/tier/sample，繪製 enrichment heatmap，建立 TO vs paired LOH concordance matrix（288,609 匹配位點）。

**關鍵發現**：

- 整體 LOH 率 36.4%（272,147/748,391）
- **TO 過度判定 LOH**：concordance 85.5%，不一致的 41,852 位點中 **95.5% 是 TO=LOH where paired=nonLOH**
- LOH enrichment 方向完全翻轉：paired LOH = FP-enriched (1.02-3.18x)，TO LOH = TP-enriched (0.85-0.96x)
- effective_hp_reads 是唯一有用的 LOH 特徵（paired AUC=0.727, TO AUC=0.564）
- TO LOH 過判三大原因：partial genotype (51.5%)、phase block 碎片化 (29.3% singleton)、PS 缺失率偏高 (11.1%)

**代表圖表**：

![fig04 — LOH Enrichment Heatmap（Sample x Mode）](../figures/20260401_observation_figures/O03/fig04_loh_enrichment_heatmap.png)

*Paired 所有樣本 enrichment >1.0（LOH=FP-enriched），TO 全部 <1.0（LOH=TP-enriched）。方向完全翻轉。*

![fig06 — TO vs Paired LOH Concordance Matrix](../figures/20260401_observation_figures/O03/fig06_to_vs_paired_loh_concordance.png)

*95.5% 不一致來自 TO 過度判定 LOH，僅 4.5% 為 paired=LOH 而 TO=nonLOH。*

---

## O4: Caller Features — Caller 特徵區分能力

**目的**：評估 variant caller 提供的 metadata（AF, GQ, DP, AD）在 paired/TO 兩種模式下區分 TP/FP 的能力。

**方法**：計算每個 caller 特徵的 per-mode AUC 與 Cohen's d，繪製 AF bin FP rate 分析，進行跨樣本穩定性檢驗。

**關鍵發現**：

- **GQ 是 paired 最強單一特徵**：AUC=0.811, Cohen's d=1.314，跨樣本穩定 0.755-0.947
- TO 無 caller feature AUC > 0.60，最佳為 AD_ref (0.597)
- **AF 在 TO 方向反轉**：AUC=0.418（高 AF = 更多 FP），FP rate 在 AF 0.8-0.9 區間高達 55.2%
- Paired 可用 GQ+DP 達 AUC~0.85，TO 完全不可行

**代表圖表**：

![fig04 — Caller Feature AUC by Mode](../figures/20260401_observation_figures/O04/fig04_caller_feature_auc_by_mode.png)

*Paired GQ/DP/AD_alt 均強（>0.78），TO 全面低於 0.60。*

![fig06 — AF Bin TP/FP Rate](../figures/20260401_observation_figures/O04/fig06_af_bin_tp_fp_rate.png)

*TO FP rate 隨 AF 單調遞增：AF 0.0-0.1 僅 14%，AF 0.8-0.9 高達 55.2%。AF 硬閾值在 TO 有害。*

---

## O5: Methylation Features — ISM 甲基化特徵效果

**目的**：評估 InterSubMod 產出的 9 個甲基化統計特徵（PairwiseMeanDist, CramersV, HeuristicScore, AlleleDelta, HPMergedDelta, PassedGating 等）對 TP/FP 的區分能力，以及與 caller 特徵的互補性。

**方法**：計算每個甲基化特徵的 AUC 排名、Cohen's d 效應量、paired/TO 方向比較，以及與 caller 特徵的 Spearman 相關矩陣。

**關鍵發現**：

- 所有 9 個甲基化特徵 AUC < 0.55，最佳為 PairwiseMeanDist TO=0.543
- **5/9 特徵在 paired/TO 間方向反轉**——無特徵可跨模式使用
- 甲基化與 caller 特徵近乎正交 (|r| < 0.26)，可提供獨立增量，但個別信號太弱
- PassedGating AUC=0.503/0.512（無效的二元判別器）

**代表圖表**：

![fig02 — 甲基化特徵 AUC 排名](../figures/20260401_observation_figures/O05/fig02_methyl_feature_auc_ranking.png)

*Paired 0.42-0.53，TO 0.47-0.54。5 個特徵方向翻轉。最強效應僅為 GQ 的 1/4.6。*

---

## O6: VerificationClass — 標籤有效性分析

**目的**：檢驗 ISM VerificationClass 四分類標籤（Noise/Weak/Strong/Subclone）是否能有效區分 TP/FP。

**方法**：計算 VerificationClass 與 truth_label 的 Cramer's V 效應量，分析 class 組成穩定性、paired-TO transition matrix、以及驅動 VC 的底層特徵。

**關鍵發現**：

- **Cramer's V=0.023**（negligible effect）——TP/FP 區分力近乎零
- Strong TP rate=98.98% vs Noise=98.73%（paired），差距 <1%
- VC 本質上是 AlleleDelta/HPMergedDelta 的離散化（eta-sq=0.604），但這些連續特徵本身對 truth 無區分力
- Paired-TO transition kappa=0.854（高度一致），但下行遷移多於上行

**代表圖表**：

![fig01 — VerificationClass 組成（跨 Dataset）](../figures/20260401_observation_figures/O06/fig01_verification_class_composition.png)

*四類比例跨 14 個 dataset 極穩定：Noise ~48-50%, Weak ~31-33%, Strong ~16-17%, Subclone ~2.1-2.4%。*

![fig02 — VerificationClass TP/FP Precision](../figures/20260401_observation_figures/O06/fig02_verification_class_tp_fp_rate.png)

*所有 class 的 TP rate 幾乎相同（paired 差距 <1%, TO 差距 <4%）。VC 對 TP/FP 判別無用。*

---

## O7: TO Phasing Quality — TO 模式 phasing 品質

**目的**：評估 tumor-only 模式下 haplotype phasing 的品質，以及 HP-based 指標在 TO 下的可靠性。

**方法**：分析 HP tag 分配率、phase block 分佈、TO/paired HP_Ratio concordance（288,609 匹配位點），以及 HP imbalance 與 LOH 的關係。

**關鍵發現**：

- TO hp_assign_rate 均值 0.924（paired 0.853），但 TP/FP 差異為零 (0.9241 vs 0.9245)
- **TO/paired HP_Ratio 連續值不相關但 LOH binary 分類一致**：Pearson r=0.001（Spearman ρ=0.006），但 LOH binary concordance=85.5%。兩者不矛盾——散點圖呈 **"X" 形交叉**（見 O07 fig05），原因是 **haplotype swap**：TO 和 Paired 使用不同 phasing anchor（germline vs somatic het SNV），對同一 LOH 區域可能分配相反的 haplotype 方向（Paired HP_Ratio=0.05 → TO HP_Ratio=0.95 或反之）。LOH binary 分類（HP_Ratio < 0.1 OR > 0.9）不受 swap 影響，因此 concordance 高；但 continuous HP_Ratio 的正/負相關互相抵消，Pearson r ≈ 0
- HP fix 後 88% regions 在 Tier A/A+
- 29.3% phase block 為 singleton（碎片化嚴重），LOH 位點 PS 缺失率 11.1%（nonLOH 僅 1.2%）
- HP imbalance 是 LOH 的結果而非獨立特徵（d=4.35 但 TP/FP AUC 僅 0.531）

**代表圖表**：

![fig06 — TO Quality Tier 分佈（Post HP-Fix）](../figures/20260401_observation_figures/O07/fig06_to_tier_distribution_post_fix.png)

*Paired: High tier TP/FP 差距 35.8pp。TO: 差距僅 1.3pp——Quality_Tier 在 TO 失去區分力。*

![fig08 — HP Imbalance vs LOH](../figures/20260401_observation_figures/O07/fig08_to_hp_imbalance_vs_loh.png)

*LOH vs nonLOH HP imbalance d=4.35（極大），但作為 TP/FP 分類器 AUC 僅 0.531。*

---

## O8: Sample Heterogeneity — 跨樣本差異評估

**目的**：量化不同癌症樣本之間的特徵分佈差異、FP rate 差異，評估哪些特徵在跨樣本場景下穩定。

**方法**：計算 per-sample FP rate、feature CV、PCA 聚類、per-sample per-feature AUC stability heatmap、平台效應（5kHz vs DORADO）。

**關鍵發現**：

- **TO FP rate 跨樣本差異 8.6 倍**：8.7% (H2009) 至 74.6% (HCC1954)
- H2009 佔全域 36.2%，depth (89x) 是其 TO 表現最佳的主因
- Depth 是 PCA 第一驅動因素（PC1=34.6%），超過癌症類型影響
- caller_gq 是唯一跨樣本穩定的 high-AUC 特徵（paired only, 0.755-0.947）
- **Basecaller 對甲基化距離有大效應**（PairwiseMeanDist r=-0.499），跨平台需 normalization
- HCC1395 chr8：LOH rate 90.8%，paired FP enrichment 23.3x（sample-specific hotspot）

**代表圖表**：

![fig02 — Sample-Specific TP/FP Rate](../figures/20260401_observation_figures/O08/fig02_sample_specific_tp_fp_rate.png)

*TO FP rate 從 8.7% 到 74.6%，65.9pp 跨度，Cohen's h=1.47（extreme effect）。*

![fig07 — Feature AUC Stability Heatmap](../figures/20260401_observation_figures/O08/fig07_feature_stability_heatmap.png)

*Paired caller_gq 跨 sample 穩定 (0.755-0.947)。TO 無任何特徵跨 sample 可靠。*

![fig09 — HCC1395 chr8 Hotspot](../figures/20260401_observation_figures/O08/fig09_hcc1395_chr8_hotspot.png)

*chr8 LOH rate=90.8%，paired FP enrichment=23.3x，但 TO chr8 FP rate 反而偏低 (0.61x)。*

---

## O10: Read-Level Methylation — Read 層級甲基化特徵

**目的**：評估 read-level 甲基化特徵（methyl_mean, methyl_std, methyl_low/high_fraction 等 12 個特徵）能否區分 ALT/REF support reads 或 TP/FP regions。

**方法**：以 86,521 reads（paired-pure only）為分析對象，計算 12 個 read-level 特徵對 ALT/REF 分類與 TP/FP region 分類的 AUC。

**關鍵發現**：

- **ALT/REF 分類完全無效**：12 個特徵 AUC 0.500-0.547，最佳 methyl_bimodal_score (0.547)
- TP/FP region 分類中等但受 clustering 膨脹：methyl_low_fraction AUC=0.737, methyl_mean AUC=0.733
- **FP variants 傾向高甲基化（非活性）區域**：FP mean methyl=0.679 vs TP=0.463
- ALT-REF 甲基化差異方向跨樣本不一致——非普遍規律

**代表圖表**：

![fig08 — Read-Level Feature AUC 排名](../figures/20260401_observation_figures/O10/fig08_read_feature_importance_for_alt_support.png)

*ALT/REF 所有 AUC 在 0.50-0.55（無效）。TP/FP region 中 methyl_low_fraction 0.737（但需注意 clustering 膨脹）。*

---

## Top 10 Cross-Observation Level A 發現

| # | 發現 | 來源 | 關鍵數字 | 影響 |
|---|------|------|---------|------|
| 1 | TO 模式下無任何單一特徵有效區分 TP/FP | O1, O4, O5, O8 | 所有 AUC < 0.58 | TO 必須依賴多特徵組合或全新特徵工程 |
| 2 | LOH Penalty 是 TO QS 失效的根本原因 | O2, O3 | TP trigger 44.5% vs FP 35.8%; 移除後 +0.045 | LOH penalty 必須在 TO 移除或反轉 |
| 3 | Paired 與 TO 是根本不同的問題空間 | O1, O5, O7 | FP rate 1% vs 31%; HP r=0.001; 5/9 方向反轉 | 必須分離 paired/TO 模型 |
| 4 | GQ 是 Paired 最強特徵，超越 QS Composite | O4, O8 | AUC=0.811; d=1.314; 跨樣本 0.755-0.947 | Paired ML baseline 以 GQ 為首要特徵 |
| 5 | 樣本間異質性極大 | O8 | TO FP rate 8.7%-74.6% (8.6x); H2009 佔 36.2% | 需要 sample-aware 策略或 LOSO-CV |
| 6 | ISM 甲基化特徵 region-level 鑑別力極弱 | O5, O10 | 最佳 AUC=0.543; 但與 caller 正交 | 需要新特徵工程方向 |
| 7 | VerificationClass 無法作為 TP/FP 過濾器 | O6 | Cramer's V=0.023; 各 class TP rate 差距 <4% | 應以連續 AlleleDelta 取代 |
| 8 | AF 在 TO 反向（高 AF = 更多 FP） | O4, O8 | AUC=0.418; FP rate @ AF 0.8-0.9 = 55.2% | AF 硬閾值在 TO 有害 |
| 9 | TO 過度判定 LOH（85.5% concordance） | O3, O7 | 95.5% discordant = TO=LOH where paired=nonLOH | TO LOH 判定不可信 |
| 10 | Read-level 甲基化對 ALT/REF 分類無用 | O10 | 12 特徵 AUC 0.500-0.547; 跨樣本 delta 不一致 | 不應投入 read-level ALT/REF 分類器 |

---

## 行動建議

| 優先級 | 行動 | 依據 | 預期效果 |
|--------|------|------|---------|
| **P0** | 移除 TO LOH penalty | O2, O3 | QS AUC +0.045（已完成） |
| **P0** | 建立 Paired/TO 分離模型策略 | O1, O5, O7 | 避免方向反轉特徵汙染模型 |
| **P1** | Phase 1A ML 特徵集: GQ + DP + 5 甲基化 + effective_hp_reads | O4, O5, O3 | Paired baseline AUC ~0.85+ |
| **P1** | 移除 VerificationClass 從 QS 決策 | O6 | 降低離散化噪音 |
| **P2** | Sample-aware calibration 或 LOSO-CV 評估 | O8 | 處理 8.6x TO FP 率差異 |
| **P2** | 執行 O9 FN 觀察（需 FN ISM 數據） | 待 FN ISM | 量化 LOH rescue 潛力 |
| **P2** | 探索 genomic context（高/低甲基化區域）作為新特徵 | O10 | FP=0.679 vs TP=0.463 暗示可用性 |

---

## 待驗證問題（已驗證 / 已更新）

### 已解決

1. **~~O10 TP/FP AUC 膨脹程度~~** ✅ **膨脹極小（-0.009），AUC 有效。** Region-level bootstrap 驗證結果：
   - Read-level AUC = 0.737 → Region-level AUC = **0.728** [95% CI: 0.686, 0.771]（N=620 regions, 2000 bootstrap）
   - 膨脹量僅 0.009（1.2%），遠小於 CI 寬度 0.085
   - ICC = 0.805（cluster 內高相關），但 AUC 未被顯著膨脹，因 region-level 聚合後信號保持
   - Per-sample 差異大：COLO829 0.600（最差，C11 低深度效應延伸）、H2009 0.840（最佳）、HCC1395_DORADO 0.823
   - **結論**：methyl_low_fraction 在 paired mode 的 region-level 鑑別力確認有效（AUC ~0.73），但遠低於 GQ（0.811），且跨樣本不穩定

### 尚未解決

2. **多特徵組合在 TO 的上限**：單一特徵全部 <0.58，組合上限待 Phase 1A ML 評估。定位 P1。
3. **Basecaller normalization 可行性**：PairwiseMeanDist 跨平台 r=-0.499。定位 P2。
4. **HCC1395 chr8 ASM+LOH block 機制**：定位 P2 待資源。
5. **O9 FN 特徵觀察**：需 FN ISM 數據。定位 P2。
6. **未覆蓋欄位補充分析**：~26 個待覆蓋。定位 P2。

---

## 認知門檻補充建議

### 給初次閱讀者的核心概念

- **AUC (Area Under ROC Curve)**：衡量分類能力的指標。0.5 = 隨機猜測，0.7 = 可接受，0.8 = 良好，0.9+ = 優秀。低於 0.5 表示預測方向反轉（用來排除反而更好）。
- **Cohen's d**：效應量指標。|d| < 0.2 = 微小，0.2-0.5 = 小，0.5-0.8 = 中等，>0.8 = 大。本報告中 GQ d=1.314（非常大）vs 甲基化最強 d=0.285（小）。
- **Cramer's V**：類別變數關聯強度。0 = 完全無關，1 = 完全關聯。V < 0.1 通常視為 negligible。
- **Paired vs TO 模式**：Paired = 有配對正常樣本作為參照；TO (Tumor-Only) = 僅有腫瘤樣本。TO 因缺乏正常參照，FP rate 從 1% 暴增至 31%。
- **LOH (Loss of Heterozygosity)**：某基因座的兩個等位基因中一個遺失，導致只剩一個 haplotype。在 paired 模式下 LOH 區域容易產生 FP，但在 TO 模式下 LOH 的判定本身不可靠（concordance 85.5%）。
- **方向反轉**：同一特徵在 paired 模式下「高值=TP」但在 TO 模式下變成「高值=FP」。這意味著該特徵不可跨模式共用，必須分離建模。

### 閱讀順序建議

1. 先讀 Top 10 發現（本文上方表格）掌握全局
2. 若關注 QS 修復：O2 -> O3
3. 若關注 ML 特徵選擇：O4 -> O5 -> O8 fig07
4. 若關注 TO 獨立主線：O1 -> O7 -> O8
5. 完整報告見 `20260401_comprehensive_observation_O1_O10_report_01.md`（含全部 82 張圖表）
