# InterSubMod 批量 VCF 分析輸出資料格式與 TP/FP 分析指南

> **版本**: v1.0  
> **日期**: 2026-01-07  
> **目的**: 說明 `run_batch_vcf_analysis.sh` 產出的資料格式，並提供透過 Python 進行 TP/FP 比較分析的完整指引。

---

## 1. 輸出資料結構說明

### 1.1 目錄結構

執行 `scripts/run_batch_vcf_analysis.sh` 後，會在指定的輸出基礎目錄下建立如下結構：

```
{BASE_OUTPUT_DIR}/{YYYYMMDD}_{MODE}_{SEQ}/
├── filtered_snv_tp/                          # TP 分析結果
│   ├── filtered_snv_tp/                      # 以 VCF 名稱為子目錄（內含各染色體資料）
│   │   └── chr{N}/{position}/{start_end}/    # 各位點詳細資料
│   ├── debug/                                # 除錯資訊
│   ├── significance_summary.csv              # ★ 關鍵：彙總統計 CSV
│   ├── significance_statistics.txt           # 統計摘要文字檔
│   ├── full_execution_analysis.log           # 完整執行日誌
│   └── heatmap_generation.log                # 圖表產生日誌
│
├── filtered_snv_fp/                          # FP 分析結果（結構同上）
│   ├── significance_summary.csv
│   └── ...
│
├── analysis/                                 # 比較分析結果
│   ├── tables/                               # 統計表格
│   │   ├── summary_stats.csv
│   │   └── discrimination_metrics.csv
│   └── plots/                                # 圖表
│       ├── dist_*.png
│       ├── scatter_*.png
│       └── roc_curves.png
│
├── filtered_snv_tp.log                       # TP 執行日誌
└── filtered_snv_fp.log                       # FP 執行日誌
```

---

## 2. significance_summary.csv 欄位詳解

此為核心分析檔案，每行代表一個 SNV 位點的甲基化異質性分析結果。

### 2.1 欄位一覽表

| # | 欄位名稱 | 資料類型 | 說明 | 範例值 |
|:---:|:---|:---|:---|:---|
| 1 | **RegionID** | int | 位點流水號（0-indexed） | `0`, `1`, `2` |
| 2 | **Chr** | string | 染色體名稱 | `chr1`, `chr19` |
| 3 | **Pos** | int | 染色體位置（1-based） | `877772` |
| 4 | **Ref** | char | 參考鹼基 | `G`, `A`, `C`, `T` |
| 5 | **Alt** | char | 變異鹼基 | `C`, `T`, `A`, `G` |
| 6 | **NumReads** | int | 定序深度（涵蓋該位點的 reads 數量） | `46`, `71` |
| 7 | **NumCpGs** | int | 分析範圍內的 CpG 位點數量 | `125`, `506` |
| 8 | **GlobalP** | float | 全域 Fisher 檢定 P-value | `0.0`, `1.0`, `1.845e-01` |
| 9 | **CramersV** | float | Cramér's V 效應量 (0~1) | `0.0`, `0.257` |
| 10 | **HeuristicScore** | float | 綜合啟發式評分 | `20.0`, `0.5138` |
| 11 | **PassedGating** | bool | 是否通過前置過濾條件 | `true`, `false` |
| 12 | **LabelDelta** | float | 標籤間甲基化差異程度 | `0.2897`, `0.0287` |
| 13 | **LabelP** | float | 標籤檢定 P-value | `1e-03`, `0.469` |
| 14 | **LabelSig** | bool | 標籤檢定是否顯著 | `true`, `false` |
| 15 | **DominantLabel** | string | 主要差異來源標籤類型 | `hp`, `allele`, `none` |
| 16 | **Stability** | float | 結構穩定性指標 | `0.0` |
| 17 | **VerificationClass** | string | 驗證分類 | `Strong`, `Subclone`, `Weak`, `Noise` |
| 18 | **LocalBestCluster** | int | 最佳局部聚類編號 (-1 表示無) | `0`, `1`, `-1` |
| 19 | **LocalBestP** | float | 最佳局部聚類 P-value | `5e-10`, `1.0` |
| 20 | **Significant** | bool | **最終顯著性判定** | `true`, `false` |

---

### 2.2 關鍵欄位詳細說明

#### 2.2.1 顯著性相關

| 欄位 | 說明 | 使用時機 |
|:---|:---|:---|
| `GlobalP` | 全域 Fisher 精確檢定的 P-value。測試所有 reads 的甲基化模式是否顯著不同於隨機分布。 | 基礎顯著性篩選 (< 0.05) |
| `CramersV` | 效應量指標，衡量甲基化狀態與 reads 分組間的關聯強度。**0** 表示無關聯，**1** 表示完全關聯。 | 過濾弱信號 (> 0.2 為中等效應) |
| `HeuristicScore` | 結合 P-value 與效應量的綜合評分。計算公式：`-log10(P) * (1 + V)`。 | 排序優先處理的位點 |
| `Significant` | 最終判定：`PassedGating = true` 且 `GlobalP <= 0.05` | 直接篩選顯著位點 |

#### 2.2.2 驗證分類 (VerificationClass)

此欄位根據多個指標綜合判定位點的可信度：

| 類別 | 意義 | 預期 TP 比例 |
|:---|:---|:---|
| `Strong` | 強信號：高顯著性 + 高效應量 + 一致的標籤差異 | **高** |
| `Subclone` | 亞克隆信號：顯著但可能來自腫瘤異質性（非單倍體） | 中 |
| `Weak` | 弱信號：P-value 邊界顯著，效應量低 | 低 |
| `Noise` | 噪音：不顯著或未通過 Gating | **低** |

#### 2.2.3 主導標籤 (DominantLabel)

表示甲基化差異的主要來源：

| 類別 | 意義 | 生物解釋 |
|:---|:---|:---|
| `hp` | Haplotype (單倍體標籤) 差異最大 | 等位基因特異性甲基化 (ASM) |
| `allele` | Allele (等位基因) 差異最大 | 突變相關甲基化改變 |
| `none` | 無顯著差異 | 無明確甲基化異質性 |

---

## 3. TP/FP 比較分析方法

### 3.1 可分析的差異特徵

以下列出可透過 Python 對 TP/FP 群組進行比較分析的特徵維度：

#### 3.1.1 顯著性指標分布

| 分析目標 | 分析方法 | 預期行為 |
|:---|:---|:---|
| 顯著率比較 | 計算 `Significant=true` 的比例 | TP 應有較高顯著率（但目前數據反向） |
| P-value 分布 | KDE/直方圖比較 `-log10(GlobalP)` | TP 應偏向低 P-value |
| V 值分布 | KDE 比較 `CramersV` | TP 應有較高 V 值分布 |

#### 3.1.2 驗證分類交叉分析

```python
# 建立 VerificationClass 與 TP/FP 的交叉表
pd.crosstab(df['Label'], df['VerificationClass'], normalize='index')
```

**預期**：TP 在 `Strong` 類別的比例應高於 FP。

#### 3.1.3 標籤來源分析

```python
# DominantLabel 在 TP/FP 間的分布
df.groupby(['Label', 'DominantLabel']).size().unstack()
```

**預期**：真正的體細胞變異 (TP) 應更常顯示 `hp` 或 `allele` 差異。

#### 3.1.4 深度與 CpG 相關性

| 特徵 | 分析問題 |
|:---|:---|
| `NumReads` | 低深度是否更容易產生假顯著？ |
| `NumCpGs` | CpG 密度是否影響檢測能力？ |
| `LabelDelta` | 效應量是否比 P-value 更能區分 TP/FP？ |

---

### 3.2 ROC/AUC 鑑別力評估

將 TP 視為正樣本，FP 視為負樣本，評估各指標作為分類器的能力：

```python
from sklearn.metrics import roc_curve, auc

# 以 CramersV 為例
fpr, tpr, thresholds = roc_curve(y_true, df['CramersV'])
roc_auc = auc(fpr, tpr)
```

**可測試指標**：
- `CramersV`
- `HeuristicScore`
- `-log10(GlobalP)`
- `LabelDelta`
- `LocalBestP`

---

## 4. 可視化圖表規劃

### 4.1 圖表類型與用途

| 圖表類型 | 檔案名稱 | 用途 |
|:---|:---|:---|
| **分布比較圖** | `dist_cramers_v.png` | 比較 TP/FP 的 V 值分布 (KDE) |
| | `dist_score.png` | 比較 TP/FP 的 Heuristic Score 分布 |
| | `dist_p_value.png` | 比較 -log10(P) 分布 |
| | `dist_reads.png` | 比較定序深度分布 (Box Plot) |
| **類別分布圖** | `dist_verification_class.png` | VerificationClass 堆疊長條圖 |
| | `dist_dominant_label.png` | DominantLabel 群組長條圖 |
| **相關性散佈圖** | `scatter_reads_vs_v.png` | 深度 vs V 值（依 TP/FP 著色） |
| | `scatter_cpgs_vs_v.png` | CpG 數量 vs V 值 |
| | `scatter_delta_vs_v.png` | LabelDelta vs CramersV |
| **評估圖表** | `roc_curves.png` | 多指標 ROC 曲線比較 |
| | `precision_recall.png` | PR 曲線（處理不平衡資料） |
| **交叉分析熱圖** | `heatmap_class_vs_label.png` | VerificationClass × DominantLabel 熱圖 |

---

### 4.2 輸出目錄分類結構

```
analysis/
├── tables/                        # 統計表格
│   ├── summary_stats.csv          # 基本統計彙總
│   ├── verification_class.csv     # VerificationClass 分布
│   ├── dominant_label.csv         # DominantLabel 分布
│   └── discrimination_metrics.csv # ROC AUC 指標
│
├── plots/
│   ├── distribution/              # 分布類圖表
│   │   ├── dist_cramers_v.png
│   │   ├── dist_score.png
│   │   ├── dist_p_value.png
│   │   ├── dist_reads.png
│   │   ├── dist_verification_class.png
│   │   └── dist_dominant_label.png
│   │
│   ├── correlation/               # 相關性圖表
│   │   ├── scatter_reads_vs_v.png
│   │   ├── scatter_cpgs_vs_v.png
│   │   └── scatter_delta_vs_v.png
│   │
│   ├── evaluation/                # 評估類圖表
│   │   ├── roc_curves.png
│   │   └── precision_recall.png
│   │
│   └── heatmaps/                  # 熱圖
│       └── heatmap_class_vs_label.png
│
└── report.md                      # 自動生成的分析報告
```

---

## 5. Python 分析執行方式

### 5.1 基本使用

```bash
python3 tools/compare_vcf_results.py \
    --output-dir <OUTPUT_ANALYSIS_DIR> \
    --labels TP FP \
    --paths <TP_RESULT_DIR> <FP_RESULT_DIR>
```

### 5.2 完整範例

```bash
python3 /big8_disk/liaoyoyo2001/InterSubMod/tools/compare_vcf_results.py \
    --output-dir /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/analysis \
    --labels TP FP \
    --paths \
        /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/filtered_snv_tp \
        /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260106_all-with-w5000_1/filtered_snv_fp
```

---

## 6. 分析結論撰寫指引

完成分析後，應將結論撰寫至：
```
docs/reports/tests/20260106_TP_FP_分析/analysis_report.md
```

報告應包含：
1. **統計彙總表**：TP/FP 的位點數、顯著率、平均 V 值等
2. **關鍵發現**：兩組間的顯著差異或異常現象
3. **ROC 分析結果**：各指標的 AUC 與最佳閾值建議
4. **圖表解讀**：每張關鍵圖表的說明與生物學意義
5. **結論與建議**：是否需要調整評分機制或閾值
