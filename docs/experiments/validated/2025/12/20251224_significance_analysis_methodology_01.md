<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 甲基化距離矩陣之群集顯著性分析方法論 (Methodology of Significance Analysis from Methylation Distance Matrices)

**日期**: 2025-12-24  
**作者**: InterSubMod 開發團隊  
**版本**: 1.0  
**狀態**: 已驗證 (Verified)

## 1. 概述 (Overview)

本文件詳細說明 InterSubMod 系統如何從讀取序列 (Read) 間的「甲基化距離矩陣」出發，經過分群、標籤關聯性檢定，最終判斷該基因位點是否具有統計顯著的「等位基因特異性甲基化 (ASM)」或「單倍型特異性甲基化 (HSM)」。

目標是解決以下問題：
> 當我們觀察到甲基化模式分為兩群時，這種分群是否與生物學標籤（如 Allele 或 Haplotype）有真實的統計關聯，亦或是隨機雜訊？

---

## 2. 分析流程圖 (Analysis Workflow)

如果不透過圖形，文字流程如下：

1.  **輸入 (Input)**: $N \times N$ 甲基化距離矩陣 (Distance Matrix)。
2.  **分群 (Clustering)**: 階層式分群法 (Hierarchical Clustering) + 動態剪枝 (Dynamic Tree Cutting)。
3.  **標籤映射 (Label Mapping)**: 將每條 Read 的生物資訊（Allele, HP, Sample）映射到分群結果。
4.  **統計檢定 (Statistical Testing)**:
    *   **Phase 1: 全域檢定 (Global Test)** - 檢驗標籤分佈與分群是否獨立。
    *   **Phase 2: 局部檢定 (Local Test)** - 找出最具特異性的群集。
    *   **Phase 3: 結構檢定 (Structure Test)** - 驗證分群的幾何合理性 (PERMANOVA)。
5.  **綜合評分 (Scoring & Decision)**: 計算啟發式分數並進行二元分類 (顯著/不顯著)。

---

## 3. 詳細步驟說明 (Detailed Steps)

### 3.1 步驟一：從距離矩陣到分群 (Clustering)

**目標**: 將具有相似甲基化模式的 Reads 聚集成群。

*   **演算法**: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)。
    *   選擇原因：相比於 Ward 或 Single-linkage，UPGMA 在處理甲基化距離時通常能產生較平衡且具生物意義的群集。
*   **最佳 $k$ 值選擇 (Automatic $k$ Selection)**:
    *   **預設**: 計算 $k=2$ 到 $k=6$ 的輪廓係數 (Silhouette Score)，選擇分數最高者。
    *   **異常值處理 (Outlier Heuristic)**: 若 $k=2$ 時產生極度不平衡的分群（如 1 vs 99），系統會自動嘗試 $k=3$。若 $k=3$ 能成功將大群再次分割且結構更合理，則採用 $k=3$。這是為了避免單一離群值掩蓋了主要的生物結構。

### 3.2 步驟二：全域關聯性檢定 (Global Association Test)

**目標**: 詢問「分群結果 (Cluster ID)」與「生物標籤 (Label)」是否統計相關？

我們針對三個維度分別進行檢定：
1.  **Allele (ALT/REF)**: 檢測 ASM (Allele-Specific Methylation)。
2.  **Haplotype (HP1/HP2)**: 檢測 HSM (Haplotype-Specific Methylation)。
3.  **Sample (Tumor/Normal)**: 檢測腫瘤特異性變化。

**使用的統計方法**:
*   **Fisher-Freeman-Halton Test (Monte Carlo)**:
    *   適用於 $R \times C$ 列聯表（例如 3個群集 $\times$ 2種單倍型）。
    *   **算法**: 使用蒙地卡羅模擬 (Monte Carlo Simulation, $10^5 \sim 10^7$ 次採樣) 來估算精確 P-value。
    *   優點：即使在樣本數少或表格稀疏時依然準確，優於卡方檢定。
*   **效應值 (Effect Size - Cramér's V)**:
    *   衡量關聯強度的指標 (0~1)。
    *   公式: $V = \sqrt{\frac{\chi^2}{n \cdot \min(k-1, r-1)}}$
    *   判定: $V > 0.3$ 表示中度關聯， $V > 0.5$ 表示高度關聯。

### 3.3 步驟三：局部特異性檢定 (Local Specificity Test)

**目標**: 若全域顯著，進一步找出「哪一個群集」是特定的。

*   **方法**: One-vs-Rest Fisher Exact Test。
*   **邏輯**: 對於每個群集 $C_i$ 與每個標籤 $L_j$ (如 HP1)，建立 $2 \times 2$ 表格：
    *   (In $C_i$, Is $L_j$)
    *   (In $C_i$, Not $L_j$)
    *   (Not $C_i$, Is $L_j$)
    *   (Not $C_i$, Not $L_j$)
*   **輸出**: 算出每個群集的 "Purity" (純度) 與 "Log Odds Ratio"。

### 3.4 步驟四：結構與離散度檢定 (Structure & Dispersion)

**目標**: 排除「雖然標籤有分群，但甲基化模式本身分不開」的偽陽性，或「群內變異過大」的情況。

1.  **PERMANOVA (Permutational MANOVA)**:
    *   利用距離矩陣檢驗群集間的幾何中心是否顯著不同。
    *   如果不顯著 ($P > 0.05$)，表示分群在幾何空間中重疊嚴重，可能是過度切割。
2.  **Homogeneity of Multivariate Dispersions (PERMDISP 概念)**:
    *   檢查群集內的變異程度 (Dispersion) 是否一致。若某群極度發散，可能會導致統計偽陽性。

---

## 4. 顯著性判斷邏輯 (The Decision Logic)

系統最終如何決定一個位點是「顯著 (Significant)」的？

我們採用 **"Gating + Scoring"** 雙重機制：

### 4.1 Gating (門檻過濾)
必須滿足以下**任一**條件才算通過初步篩選：
*   **Global P-value (Allele)** $\le$ 0.1
*   **Global P-value (Haplotype)** $\le$ 0.1

### 4.2 Heuristic Score (啟發式評分)
將多個指標綜合成單一分數 $S$，用於排序重要性。

$$ S = -\log_{10}(P_{best}) + (2.0 \times V_{best}) + \text{Penalties} $$

其中：
*   $P_{best} = \min(P_{allele}, P_{hp})$
    *   **修正**: 若 $P=0$ (蒙地卡羅模擬無更極端值)，則 $S$ 直接設為最大值 (20.0)。
*   $V_{best} = \max(V_{allele}, V_{hp})$
*   **Penalties**:
    *   若 PERMANOVA 不顯著，分數 $\times 0.7$。
    *   若群內離散度 (Dispersion) 異常，分數 $\times 0.5$。

### 4.3 最終判定 (Categorization)
*   **Significant**: Passed Gating **AND** ($P_{best} \le 0.05$).
*   **Non-Significant**: 其他。

---

## 5. 延伸驗證建議 (Verification & Extension)

為了驗證此流程的可靠性，建議進行以下實驗：

1.  **混洗測試 (Shuffle Test)**:
    *   隨機打亂 Reads 的標籤 (Permutation)，重新計算 P-value。
    *   預期：打亂後的 P-value 分佈應為均勻分佈 (Uniform)，且高分比例極低。
2.  **降採樣測試 (Subsampling)**:
    *   隨機移除 50% Reads，檢查顯著性是否保持穩定。
3.  **合成數據 (Synthetic Data)**:
    *   人工生成完全分離 (Perfect separation) 與完全混合 (Perfect mixture) 的距離矩陣，驗證 $P$ 值是否分別趨近於 0 與 1。

---

## 6. 結論 (Conclusion)

本方法論結合了非參數統計 (Monte CarloFisher) 與幾何距離分析 (UPGMA/PERMANOVA)，能夠穩健地同時檢測 ASM 與 HSM。

*   **優勢**:
    *   **雙重檢定**: 同時考慮 Allele 與 Haplotype，避免漏掉僅在單倍型層次顯著的位點 (如 `chr19-verification` 案例)。
    *   **抗噪聲**: 透過 Outlier Heuristic 與 Structure Test 排除偽陽性。
    *   **適應性**: 適用於不同覆蓋度與甲基化水平的區域。

---
*文件生成於: /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/20251224_significance_analysis_methodology.md*
