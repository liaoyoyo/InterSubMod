<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 2026-01-06 TP/FP 批量分析計畫

本計畫旨在透過 InterSubMod 的批量分析流程，針對 **真陽性 (TP, True Positive)** 與 **假陽性 (FP, False Positive)** 的 VCF 資料進行深入的比較分析，以評估模型在不同類型位點上的鑑別能力與特徵分佈差異。

## 1. 輸入資料與格式

分析將基於 `run_batch_vcf_analysis.sh` 產出的標準輸出目錄結構，主要依賴各分析目錄下的 `significance_summary.csv` 匯總檔案。

### 資料來源
- **TP (True Positive)**: `/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/filtered_snv_tp/significance_summary.csv`
- **FP (False Positive)**: `/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/filtered_snv_fp/significance_summary.csv`

### 關鍵欄位 (Metrics)
| 欄位名稱 | 說明 | 用途 |
| :--- | :--- | :--- |
| **NumReads** | 定序深度 (Reads count) | 分析深度對顯著性的影響 |
| **GlobalP** | 全域 FISHER 檢定 P-value | 評估顯著性的主要指標 |
| **CramersV** | 關聯強度 (Effect Size) | 評估甲基化與變異關聯的強度 (0~1) |
| **HeuristicScore** | 綜合啟發式分數 | 結合 P-value 與 V 值的綜合評分 |
| **LabelDelta** | 標籤差異 (Haplotype/Allele) | 評估單倍體/等位基因間的甲基化差異 |
| **Significant** | 是否顯著 (Boolean) | 最終判斷結果 (通常 P<=0.05 且通過過濾) |

## 2. 分析目標與策略

### 2.1 顯著性分佈與差異分析 (Significance Distribution)
**目標**: 比較 TP 與 FP 在顯著性指標上的分佈差異，驗證模型對 TP 的敏感度與對 FP 的特異度。

*   **分析方法**:
    *   計算 TP 與 FP 群組中 `Significant` (顯著位點) 的比例。
    *   繪製 `CramersV` 與 `HeuristicScore` 的直方圖 (Histogram/KDE)，觀察兩者的分佈重疊情況。
    *   **預期**: TP 的 V 值與分數應顯著高於 FP。

### 2.2 特徵關聯分析 (Feature Correlation)
**目標**: 探討定序深度 (Read Depth) 等特徵是否會造成 FP 或影響 TP 的檢出。

*   **分析方法**:
    *   **Scatter Plot**: X軸為 `NumReads`，Y軸為 `CramersV` 或 `HeuristicScore`。
    *   **Box Plot**: 比較不同深度區間 (e.g., <20, 20-50, >50) 的顯著性差異。
    *   **觀察**: 是否在低深度下容易產生高 V 值的 FP？

### 2.3 鑑別力評估 (Discrimination Power)
**目標**: 將 TP 視為正樣本，FP 視為負樣本，評估各指標的分類能力。

*   **分析方法**:
    *   **ROC Curve**: 以 `CramersV`、`GlobalP` (取 -log10)、`HeuristicScore` 作為閾值變數，繪製 ROC 曲線並計算 AUC。
    *   **Precision-Recall Curve**: 針對不平衡資料 (TP 數量通常遠大於 FP) 進行評估。

## 3. 視覺化圖表規劃

分析工具將產生以下關鍵圖表，並存放於 `analysis/plots/`：

1.  **`dist_cramers_v.png`**: TP/FP Cramer's V 分佈比較圖 (KDE Plot)。
2.  **`dist_score.png`**: TP/FP Heuristic Score 分佈比較圖。
3.  **`scatter_reads_vs_v.png`**: 深度與 V 值散佈圖 (依 TP/FP 著色)。
4.  **`roc_curves.png`**: 包含 V 值與 Score 的 ROC 曲線比較。
5.  **`bar_significance_rate.png`**: TP/FP 顯著比例長條圖。

## 4. 工具修改計畫 (`compare_vcf_results.py`)

修改現有的 Python 比較工具以支援上述分析：

1.  **資料載入優化**: 優先讀取 `significance_summary.csv`，若無則降級使用舊版目錄掃描。
2.  **欄位標準化**: 確保不同版本的 CSV 欄位名稱對齊。
3.  **新增繪圖模組**: 引入 `scikit-learn` (計算 ROC/AUC) 與 `seaborn` (進階繪圖)。
4.  **輸出結構**:
    ```text
    analysis/
    ├── tables/
    │   ├── summary_stats.csv      # 彙總統計
    │   └── discrimination_metrics.csv # AUC 等指標
    └── plots/
        ├── dist_*.png
        ├── scatter_*.png
        └── roc_*.png
    ```

## 5. 預期產出與報告

最終將於 `/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/tests/20260106_TP_FP_分析` 產出完整報告，包含：
- 數據彙總表。
- 關鍵圖表解讀。
- TP/FP 鑑別的閾值建議 (基於 ROC 分析)。
