# 20251223_Significance_Analysis_Strategy

## 1. 概述與目標 (Overview & Goals)

### 1.1 背景與核心假設修正
InterSubMod 專案旨在探討 Somatic 變異與甲基化結構的關聯。我們採用穩健的假設：
> **假設**：在部分高信心的 Somatic 變異區域，其真實變異可能伴隨著可重現的 Read-level 甲基化結構差異。
> **限制**：此現象受 PMD、CNV 或低覆蓋度影響，本模組目標是**篩選出具備此特徵的高信心位點**。

### 1.2 目標功能
本模組將執行以下層次的統計驗證：
1.  **全域關聯性**: 檢定甲基化分群與 **ALT/REF**、**HP**、**Tumor/Normal** 標籤的關聯。
2.  **局部解釋性**: 識別驅動關聯的特定分群與其效應量。
3.  **結構真實性**: 排除隨機噪音與群內分散度 (Dispersion) 對分群的假性影響。
4.  **結果校準**: 轉化統計指標為經 SEQC2 真值校準的「信心分數」。

---

## 2. 統計方法論 (Statistical Methodology)

### Phase 1: 全域關聯性檢定 (Global Association Test)

**優先順序**: **ALT/REF** > **HP** > **Tumor/Normal**。

**方法**: 
*   **RxC 列聯表分析**: 建立 $K$ (Clusters) $\times M$ (Labels) 的列聯表。
    *   **Fisher-Freeman-Halton Test (Monte Carlo)** (主力): 
        *   **Stopping Rule**: 採用置信區間 (Confidence Interval) 判定。當 $\hat{p}$ 的 99% CI 完全高於或低於顯著門檻 (0.05) 時才提早停止，避免估壓不準。
        *   **預設次數**: 2,000 ~ 10,000 次。
    *   **Chi-square Test** (輔助): 僅在期望值 $\ge 5$ 的格子比例 > 80% 時報告，否則標記為不可靠。
*   **效應量 (Effect Size)**:
    *   **Cramér’s V**: 若 Chi-square 不適用 (稀疏)，標記 `v_reliable = false`，僅供參考。
*   **直接距離與離散度檢定 (Distance & Dispersion)**:
    *   **Data Preparation**: **必須先過濾 Reads**。僅保留與其他 Reads 有足夠 "Effective Overlap" 的 Reads，確保距離矩陣結構完整，**嚴禁**對無效距離進行插補 (Imputation)。
    *   **PERMANOVA (Pseudo-F)**: 檢驗群間中心差異。
    *   **Dispersion Check (關鍵)**: 計算每群到其幾何中心的**平均距離 (Mean Distance to Centroid)**。若各群平均距離差異過大 (例如 ANOVA P < 0.05)，標記 Dispersion 警告，避免將「分散度差異」誤判為「中心差異」。

### Phase 2: 局部解釋與效應量 (Local Interpretation & Effect Size)

**前提**: 僅對通過 Global 檢定 (e.g., P < 0.05) 的 Region 執行。

**方法**:
*   **Post-hoc Fisher's Exact Test**: 採用 **One-vs-Rest** 策略。
    *   定義: Rows=[Target Label, Other Labels], Cols=[Target Cluster, Other Clusters]。
    *   **Two-sided 定義**: 嚴格比照 R `fisher.test` (sum of probabilities $\le P_{obs}$)，確保與標準統計軟體一致。
*   **效應量**:
    *   **Odds Ratio (OR)**: 使用 **Haldane-Anscombe correction** (+0.5) 必免除零。
    *   **Delta Proportion ($\Delta P$)**: $P(L|C_k) - P(L|\text{Rest})$。
*   **多重校正 (Hierarchy)**:
    *   **C++ 端**: 輸出原始 P-value。
    *   **Python 端**: 執行 Global BH 校正 -> 產出 $q_{global}$ -> 對通過者做 Local BH 校正。

### Phase 3: 結構與真實性驗證 (Structural & Realism Validation)

**Gating**: 僅對 Top 候選位點或 Global 顯著位點執行 (避免計算爆炸)。

*   **Association Realism (Permutation)**:
    *   打散 Labels (Shuffle) 999 次，計算 Cramér’s V 或 Pseudo-F 的 Empirical P-value。
    *   **Null Check**: 確認 P-value 分布在隨機資料上近似 Uniform(0,1)。
*   **Clustering Realism (Bootstrap)**:
    *   **Bootstrap**: 對 **CpG sites** 進行重抽樣 (N=200)，設 Early Stop (若 ARI 已穩定低則停)。
    *   **指標**: **ARI (Adjusted Rand Index)**。ARI < 0.2 視為隨機結構。

### Phase 4: 評分與校準 (Scoring & Calibration)

**方法**:
*   **特徵工程 (Feature Schema)**:
    *   **品質特徵**: `n_reads_valid`, `n_cpg`, `effective_overlap_median`, `mapq_mean`, `hp0_rate` (避免模型學習到覆蓋度偏差)。
    *   **統計特徵**: `-log10(p_global)`, `cramers_v`, `max_log_or`, `silhouette_mean`.
*   **校準模型**:
    *   **Truth Alignment**: 
        *   Anchor SNV 在 SEQC2 Truth Set (TP) -> Region TP。
        *   Anchor SNV 在 Callable Region 但無變異 -> Region FP。
        *   其他 (Unknown) -> 排除。
    *   **Metric**: 使用 **AUPRC** (相對於 Baseline) 與 **Calibration Curve** (Brier Score)，不單看 AUROC。

---

## 3. 實作架構 (Implementation Architecture)

### 3.1 核心模組 (C++)

*   **`Statistics` Module**:
    *   `FisherExact`: 實作 R-compatible two-sided logic, MC with CI-stop.
    *   `Metrics`: Cramér’s V (w/ sparsity check), Haldane OR.
*   **`StructureTest` Module**:
    *   `PERMANOVA`: 實作 F-test，需處理 Distance Matrix filtering。
    *   `Dispersion`: 實作 Mean-to-centroid distance 比較。

### 3.2 輸出格式 (CSV)

*   **`significance_report.csv`**:
    *   **Metadata**: `region_id`, `valid_flag` (是否可算), `invalid_reason` (e.g. low_overlap).
    *   **Quality**: `n_reads`, `n_cpg`, `overlap_score`, `dispersion_val`.
    *   **Global**: `p_alt`, `v_alt`, `p_perm_alt`.
    *   **Local (Best)**: `best_cluster_id`, `p_local`, `log_or`.

---

## 4. 執行步驟 (Roadmap Summary)

1.  **Stats Library**: 實作高精度 Fisher 與基礎統計。
2.  **Global/Local Tests**: 實作核心檢定與 Gating 邏輯 (Global Pass -> Local)。
3.  **Structural Validation**: 實作 PERMANOVA, Dispersion, Bootstrap (注意效能管控)。
4.  **Calibration**: Python 端訓練與驗證。

## 5. 結論
本策略已整合三次評估報告的建議，建立了多層次 (Global/Local/Structural) 的嚴謹驗證體系。特別強調了 **Dispersion Check**、**停止規則的統計有效性** 與 **完整的特徵工程**，以確保最終輸出的信心分數具有高度的生物統計公信力。
