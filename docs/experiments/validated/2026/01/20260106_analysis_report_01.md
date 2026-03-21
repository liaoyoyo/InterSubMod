<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 2026-01-07 TP/FP 完整分析報告

> **分析版本**: v2.0 (增強版)  
> **執行日期**: 2026-01-07  
> **資料來源**: `/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/`

---

## 1. 分析概況

| 項目 | TP (True Positive) | FP (False Positive) |
|:---|:---:|:---:|
| **分析位點數** | 30,476 | 4,822 |
| **顯著位點數** | 8,710 | 1,968 |
| **顯著比例** | 28.6% | **40.8%** ⚠️ |
| **平均深度** | 71.4 reads | 59.8 reads |
| **平均 CpGs** | 97.4 | 98.0 |
| **平均 Cramér's V** | 0.052 | 0.017 |
| **平均 LabelDelta** | 0.039 | 0.077 |

> [!WARNING]
> **異常發現**：FP 的顯著比例 (40.8%) 反而高於 TP (28.6%)。這顯示目前的顯著性閾值 (`PassedGating` + `GlobalP ≤ 0.05`) **無法有效區分生物學上的真假陽性位點**。

---

## 2. VerificationClass 分布分析

| 類別 | FP | TP | 差異 |
|:---|:---:|:---:|:---|
| **Strong** | 37.1% | 24.5% | FP 高 12.6% |
| **Weak** | 31.6% | 43.5% | TP 高 11.9% |
| **Subclone** | 3.7% | 4.0% | 相近 |
| **Noise** | 27.6% | 27.9% | 相近 |

### 解讀
- **FP 在 `Strong` 類別的比例反而較高**，這是反直覺的。可能原因：
  1. FP 位點位於高甲基化變異區域（如 CpG islands 邊界）
  2. 假陽性 SNV 位點可能與真實的 ASM 區域重疊
  3. 評分機制對某些 artifact 過於敏感

- **TP 更傾向於 `Weak` 類別**，顯示許多真實體細胞變異位點周圍的甲基化差異信號較弱

---

## 3. DominantLabel 分布分析

| 標籤類型 | FP | TP | 差異 |
|:---|:---:|:---:|:---|
| **hp** (單倍體) | 69.3% | 68.5% | 相近 |
| **allele** (等位基因) | 26.1% | 17.0% | FP 高 9.1% |
| **none** | 4.6% | 14.5% | TP 高 9.9% |

### 解讀
- **FP 更常顯示 `allele` 差異**，這可能是因為：
  1. 假陽性 SNV call 本身的 allele 資訊不正確，導致基於 allele 的甲基化分群產生假差異
  2. FP 位點可能位於重複序列區域，mapping 錯誤造成 allele 混亂

- **TP 更常為 `none`**，顯示真實變異位點周圍可能沒有明顯的等位基因特異性甲基化

---

## 4. ROC 鑑別力評估

| 指標 | AUC | 最佳閾值 | Sensitivity | Specificity |
|:---|:---:|:---:|:---:|:---:|
| Cramér's V | 0.519 | 0.257 | 6.1% | 98.0% |
| Heuristic Score | 0.443 | 20.47 | 5.9% | 98.0% |
| -log10(P-value) | 0.438 | ∞ | 0.0% | 100% |
| Label Delta | 0.409 | -0.008 | 96.6% | 6.1% |

### 解讀
- **所有 AUC 值接近或低於 0.5**，表示這些指標作為 TP/FP 分類器的效果等同或不如隨機猜測
- **Heuristic Score 和 P-value 的 AUC < 0.5**，意味著這些指標實際上**反向預測**（FP 的值反而更極端）
- **Label Delta AUC = 0.409** 也顯示反向趨勢，這與上述發現一致：FP 的 LabelDelta 平均值 (0.077) 高於 TP (0.039)

---

## 5. 圖表產出

分析工具產生了以下圖表，存放於 `analysis_enhanced/plots/`：

### 分布圖 (`distribution/`)
| 圖表 | 說明 |
|:---|:---|
| `dist_cramers_v.png` | TP/FP 的 Cramér's V 分布 (KDE) |
| `dist_score.png` | Heuristic Score 分布 |
| `dist_p_value.png` | -log10(P-value) 分布 |
| `dist_reads.png` | 定序深度 Box Plot |
| `dist_label_delta.png` | LabelDelta 分布 |
| `dist_verification_class.png` | VerificationClass 堆疊長條圖 |
| `dist_dominant_label.png` | DominantLabel 群組長條圖 |

### 相關性圖 (`correlation/`)
| 圖表 | 說明 |
|:---|:---|
| `scatter_reads_vs_v.png` | 深度 vs V 值散佈圖 |
| `scatter_cpgs_vs_v.png` | CpG 數量 vs V 值 |
| `scatter_delta_vs_v.png` | LabelDelta vs V 值 |
| `scatter_reads_vs_score.png` | 深度 vs Heuristic Score |

### 評估圖 (`evaluation/`)
| 圖表 | 說明 |
|:---|:---|
| `roc_curves.png` | ROC 曲線（多指標比較） |
| `precision_recall.png` | Precision-Recall 曲線 |

### 熱圖 (`heatmaps/`)
| 圖表 | 說明 |
|:---|:---|
| `heatmap_class_vs_label_TP.png` | TP 的 VerificationClass × DominantLabel 交叉熱圖 |
| `heatmap_class_vs_label_FP.png` | FP 的交叉熱圖 |

---

## 6. 結論與建議

### 6.1 核心問題

目前的甲基化異質性分析方法 **無法區分 TP 與 FP 體細胞變異**。這表示：
1. 甲基化信號本身不足以作為區分真假變異的唯一依據
2. 現有的評分機制（P-value, Cramér's V, Heuristic Score）需要重新設計
3. FP 位點可能並非「沒有甲基化信號」，而是「有不同來源的甲基化信號」

### 6.2 改進方向

1. **重新定義 FP 集合**  
   目前的 FP 可能包含真實的甲基化異質性區域。建議檢視 V > 0.3 的 FP 位點，確認其是否為：
   - 定序 artifact (可用 IGV 檢視)
   - 潛在的 germline 變異
   - 或確實位於 ASM 區域

2. **加入深度懲罰**  
   低深度 (< 30 reads) 可能更容易產生假顯著。建議在評分中加入深度權重。

3. **結合多指標過濾**  
   單一指標無法區分 TP/FP，但可考慮組合過濾：
   ```
   Significant AND VerificationClass = Strong AND DominantLabel = hp
   ```

4. **改用局部最佳聚類指標**  
   `LocalBestP` 尚未分析，可作為下一步探索方向。

---

## 7. 檔案位置

| 檔案 | 路徑 |
|:---|:---|
| 輸出格式說明 | `docs/development/plans/2026_01_06/20260107_output_data_format_analysis.md` |
| 增強版分析工具 | `tools/compare_vcf_results.py` |
| 自動生成報告 | `output/.../analysis_enhanced/report.md` |
| 本分析報告 | `docs/reports/tests/20260106_TP_FP_分析/analysis_report.md` |
