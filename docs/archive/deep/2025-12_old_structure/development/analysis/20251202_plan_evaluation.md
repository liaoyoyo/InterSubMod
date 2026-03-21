# 1202_研究計劃評估與確認事項

## 1. 計劃整體評估

### 1.1 優點與合理性

✅ **架構清晰**：五個階段劃分合理，從距離計算 → 聚類 → 視覺化 → 統計分析 → 報告輸出，邏輯流程完整。

✅ **方法選擇適當**：

- 距離度量選擇多樣且適用於甲基化資料（NHD, Euclidean, Pearson, Jaccard）
- 聚類方法（UPGMA, NJ, Ward）適合演化樹構建
- Bootstrap 驗證是標準的穩定性評估方法
- 統計方法（PERMANOVA, Mantel Test, Phylogenetic Signal）對應研究問題

✅ **整合現有進度**：計劃與現有實作相符

- 已完成的距離矩陣計算功能（參見 `DistanceMatrix.cpp`）
- 現有配置系統（`Config.hpp`）已支援多種距離度量
- 已有 strand-aware 設計

✅ **技術棧選擇合理**：C++ 處理核心計算，Python 處理視覺化與統計分析

---

## 2. 關鍵問題與需要確認的部分

### 2.1 ⚠️ **階段一：距離矩陣計算**

#### 已完成的部分

根據 `docs/reports/2025_12_01_distance_matrix_implementation.md`，距離矩陣功能已經實作完成：

- ✅ 支援 5 種距離度量（NHD, L1, L2, CORR, JACCARD）
- ✅ Strand-aware 分離計算
- ✅ 缺失值處理（`min_common_coverage`）
- ✅ CSV 輸出

#### 需要確認

❓ **Q1: 是否需要同時支援多種距離度量的批次計算？**

- 當前 `Config.hpp` 支援 `std::vector<DistanceMetricType> distance_metrics`
- 但實作可能只處理單一 metric
- **建議**：確認是否需要一次性輸出所有度量的距離矩陣，還是分別執行

❓ **Q2: 距離矩陣輸出格式是否足夠？**

- 當前輸出 CSV 格式
- Python 聚類分析可能需要額外的 metadata（Read tags, 位點資訊）
- **建議**：設計一個統一的輸出格式，包含：
  - 距離矩陣（CSV 或 NPY）
  - Read metadata（TSV：read_id, tumor/normal, HP, strand, ALT/REF）
  - CpG 位點資訊（TSV：位置, 覆蓋度）

---

### 2.2 ⚠️ **階段二：聚類與 Bootstrap**

#### 主要問題

❓❓ **Q3: C++ 聚類庫的選擇**（高優先級）
當前計劃提到「使用現有庫如 FastTree 或自行實作」，需要明確：

**選項 A：使用 C++ 聚類庫**

- **優點**：效能高，整合簡單
- **缺點**：
  - FastTree 是系統發育樹工具，需要序列資料，不直接支援距離矩陣輸入
  - 其他 C++ 聚類庫（如 `mlpack`）功能有限，不一定支援 Bootstrap
- **可行性**：★★☆☆☆

**選項 B：呼叫 Python 聚類套件**（推薦）

- **優點**：
  - `scipy.cluster.hierarchy` 成熟穩定，支援多種 linkage 方法
  - `scikit-learn` 提供豐富的聚類演算法
  - 易於實作 Bootstrap（重抽樣後重新計算）
- **缺點**：需要跨語言呼叫
- **可行性**：★★★★★

**選項 C：自行實作 UPGMA/NJ**

- **優點**：完全控制
- **缺點**：開發時間長，容易有 bug
- **可行性**：★★☆☆☆

**🔴 建議決策**：採用選項 B（Python），工作流程為：

1. C++ 輸出距離矩陣與 metadata
2. Python 腳本讀取並執行聚類
3. Python 輸出演化樹（Newick 格式）與圖片

---

❓❓ **Q4: Bootstrap 的實作細節**（高優先級）

計劃提到「對 CpG 位點進行有放回重抽樣」，需要確認：

**方案 A：在 C++ 層級實作 Bootstrap**

- 重複 N 次：重抽樣 CpG → 重新計算距離矩陣 → 輸出給 Python
- **問題**：每次都要重新執行 C++ 主程式，效率低

**方案 B：在 Python 層級實作 Bootstrap**（推薦）

- C++ 輸出「Read-Methylation 原始矩陣」（非距離矩陣）
- Python 執行：
  1. 讀取原始矩陣
  2. 重抽樣 CpG columns N 次
  3. 每次重新計算距離矩陣 → 聚類
  4. 比較拓撲結構，計算 support values
- **優點**：靈活、快速測試不同參數

**🔴 建議決策**：採用方案 B，需要：

- C++ 額外輸出 `methylation_matrix.csv`（Reads × CpGs，值為 0/1/-1）
- Python 腳本實作 Bootstrap 邏輯

---

❓ **Q5: 演化樹格式與儲存**

- **Newick 格式**：標準的樹結構文字格式，適合儲存拓撲與分支長度
- **需要包含**：
  - Node labels（Read IDs）
  - Branch lengths（對應距離）
  - Bootstrap support values（如 `(A:0.1,B:0.2)95:0.3`，95 為 support）
- **🟢 無問題**，scipy 可直接輸出

---

### 2.3 ⚠️ **階段三：視覺化**

#### 工具確認

✅ **Python 視覺化方案無問題**：

- `seaborn.clustermap` 是最佳選擇
- 可自動處理聚類 + heatmap + dendrogram
- 支援 annotation bars（標示 tags）

#### 需要確認

❓ **Q6: Annotation Bars 的設計**
需要在 Heatmap 旁顯示多個標籤維度：

- Tumor/Normal（2 色）
- ALT/REF（2 色）
- HP Tag（可能 6-7 種：HP1, HP2, HP1-1, HP2-1, HP3, Unphased, None）
- Strand（2 色）

**建議**：

- 使用 `seaborn.clustermap` 的 `row_colors` 參數
- 每個維度一個顏色條（共 4 條）
- 需要設計配色方案（colorblind-friendly）

❓ **Q7: 圖片輸出規格**

- **解析度**：建議 300 DPI（出版品質）
- **格式**：PNG（展示）、PDF（可編輯）、SVG（向量圖）
- **尺寸**：依據 Read 數量動態調整（避免過度壓縮）

**🟢 補充建議**：在計劃中明確列出圖片輸出參數

---

### 2.4 ⚠️ **階段四：統計分析**

#### 方法可行性評估

✅ **PERMANOVA**（Permutational MANOVA）

- **套件**：`scikit-bio.diversity.beta_diversity` 的 `permanova`
- **輸入**：距離矩陣 + 分組標籤（如 Tumor vs Normal）
- **輸出**：F-statistic, p-value
- **適用性**：★★★★★，直接適用於距離矩陣

✅ **Mantel Test**

- **套件**：`scikit-bio.stats.distance.mantel`
- **輸入**：兩個距離矩陣（甲基化距離 vs 標籤距離）
- **輸出**：相關係數, p-value
- **適用性**：★★★★☆，需要定義「標籤距離」

⚠️ **Phylogenetic Signal (Lambda 值)**

- **問題**：Lambda 是用於「連續性狀」在演化樹上的分佈，不一定適合離散標籤（如 Tumor/Normal）
- **替代方案**：
  - **Parsimony Score**：計算演化樹上標籤變化的最小次數
  - **Clade Enrichment**：Fisher's Exact Test（計劃已提到）
- **🔴 需要釐清**：具體要用哪種方法測量「Tag 與演化樹拓撲的關聯」

#### 需要確認

❓❓ **Q8: 多維度標籤的統計分析策略**（高優先級）

現有 4 個維度的標籤：

1. Tumor/Normal
2. ALT/REF
3. HP Tag（多類別）
4. Strand（Forward/Reverse）

**問題**：

- 是否每個維度獨立分析？
- 是否考慮組合效應（如 Tumor+ALT vs Normal+REF）？
- HP Tag 是多類別（6-7 種），如何處理？

**建議策略**：

1. **單維度分析**：每個 Tag 獨立執行 PERMANOVA
2. **多維度分析**：使用 PERMANOVA 的多因子設計（類似 Two-way ANOVA）
3. **分層分析**：先按 Tumor/Normal 分組，再分析 HP 的效應

**🔴 需要決策**：在計劃中明確分析層次與優先級

---

❓ **Q9: Strand 的處理**
計劃提到分析「正反股」標籤，但距離矩陣已經是 strand-specific（分離計算）：

- `distance_forward.csv`：只包含 Forward reads
- `distance_reverse.csv`：只包含 Reverse reads

**問題**：如果已經分離，那「Strand」作為統計分析的維度意義何在？

**可能情境**：

1. **在合併的距離矩陣中分析**：使用 `distance_matrix.csv`（包含所有 reads），測試 Forward vs Reverse 是否聚類分離
2. **分別分析兩個矩陣**：Forward 與 Reverse 各自進行聚類與統計

**🔴 需要釐清**：Strand 的分析目的是什麼？

---

### 2.5 ⚠️ **階段五：報告輸出**

#### 架構合理性

✅ **單一位點分析 + 多位點綜合分析** 的設計合理

#### 需要確認

❓ **Q10: 輸出檔案組織結構**

建議結構（需確認）：

```
output/
├── single_locus/
│   ├── chr19_29283968/
│   │   ├── methylation_matrix.csv
│   │   ├── distance_matrix.csv
│   │   ├── distance_forward.csv
│   │   ├── distance_reverse.csv
│   │   ├── tree.newick
│   │   ├── heatmap.png
│   │   ├── statistics.txt
│   │   └── bootstrap_support.txt
│   └── chr19_XXXXXXX/
│       └── ...
└── aggregate/
    ├── summary_statistics.tsv
    ├── all_permanova_results.tsv
    ├── tag_enrichment_summary.tsv
    └── figures/
        ├── overall_tag_distribution.png
        └── p_value_histogram.png
```

**🔴 需要確認**：輸出結構是否符合後續分析需求

---

## 3. 資料流程與介面設計

### 3.1 C++ → Python 資料介面

**當前狀態**：

- C++ 已輸出 `distance_matrix.csv`, `reads.tsv`

**需要新增**（為了支援 Bootstrap 與統計分析）：

| 檔案 | 格式 | 內容 | 用途 |
|------|------|------|------|
| `methylation_matrix.csv` | CSV | Reads × CpGs (0/1/-1) | Bootstrap 重抽樣 |
| `cpg_sites.tsv` | TSV | chr, pos, coverage | 位點資訊 |
| `read_metadata.tsv` | TSV | read_id, tumor_normal, alt_ref, hp, strand | 統計分析分組 |

**🔴 需要確認**：這些檔案是否已經存在？若否，需要在 C++ 新增輸出功能

---

### 3.2 Python 分析流程

建議模組化設計：

```python
scripts/
├── clustering/
│   ├── cluster_analysis.py       # 主要聚類腳本
│   ├── bootstrap.py               # Bootstrap 驗證
│   └── tree_utils.py              # 演化樹處理
├── visualization/
│   ├── plot_heatmap.py            # Heatmap 繪製
│   └── plot_tree.py               # 演化樹繪製
├── statistics/
│   ├── permanova_analysis.py      # PERMANOVA
│   ├── mantel_test.py             # Mantel Test
│   └── enrichment_test.py         # Clade enrichment
└── reporting/
    ├── single_locus_report.py     # 單一位點報告
    └── aggregate_report.py        # 綜合報告
```

**🔴 需要確認**：這樣的模組劃分是否合理

---

## 4. 技術挑戰與風險

### 4.1 效能與記憶體

⚠️ **Bootstrap 計算量**

- 假設 N=100 bootstrap iterations
- 每個 iteration 需要重新計算距離矩陣與聚類
- **估計時間**：若單次聚類 < 1秒，總計 100 秒/位點
- **多位點處理**：30,000 位點 × 100 秒 = 833 小時（34 天）

**🔴 風險緩解**：

1. 僅對「有足夠 reads」的位點執行 Bootstrap（設定閾值，如 N_reads > 20）
2. 使用多進程平行化（Python `multiprocessing`）
3. 降低 Bootstrap 次數（如 N=50）

---

### 4.2 統計檢定的多重檢驗校正

⚠️ **多重比較問題**

- 針對 30,000 個位點，每個執行 PERMANOVA
- 需要 p-value 校正（Bonferroni, FDR）

**🔴 建議**：在計劃中加入「多重檢驗校正」段落

---

### 4.3 演化樹與甲基化資料的適配性

⚠️ **理論考量**

- 傳統演化樹假設「垂直遺傳」（親子關係）
- 甲基化資料反映的是「表觀遺傳狀態」，不一定有演化關係
- Reads 來自同一細胞群體，更像是「狀態聚類」而非「演化樹」

**問題**：使用 UPGMA/NJ 這類演化樹方法是否合適？

**替代方案**：

- **一般聚類方法**：Ward's hierarchical clustering（不假設演化）
- **Network-based methods**：Minimum Spanning Tree

**🔴 需要釐清**：研究的理論框架是什麼？是否真的需要「演化樹」的概念？

---

## 5. 建議補充的內容

### 5.1 明確的實作優先級

建議在計劃中加入：

**Phase 1（核心功能）：**

1. 確認距離矩陣輸出格式
2. 開發基礎聚類腳本（無 Bootstrap）
3. 開發基礎 Heatmap 繪圖
4. 單一位點的完整流程測試

**Phase 2（統計驗證）：**

1. 實作 PERMANOVA 分析
2. 實作 Bootstrap 驗證
3. 多位點批次處理

**Phase 3（進階分析）：**

1. Mantel Test
2. Enrichment Analysis
3. 綜合報告生成

---

### 5.2 測試資料與驗證策略

建議加入：

- **單元測試**：每個 Python 模組的測試
- **整合測試**：使用少量位點（如 10 個）測試完整流程
- **邊界案例**：
  - 只有 1-2 個 reads
  - 所有 reads 都是同一 Tag
  - 缺失值極多

---

### 5.3 參數調優計劃

需要調整的參數：

- 距離度量選擇（5 種）
- 聚類方法選擇（UPGMA, Ward's, etc.）
- Bootstrap 次數（50? 100? 1000?）
- 最小 Read 數閾值
- 最小共同 CpG 覆蓋數

**建議**：設計參數掃描實驗（先用小規模資料測試）

---

## 6. 具體行動建議

### 立即需要確認的問題（按優先級排序）

1. **Q3**: 聚類實作方案（C++ vs Python）→ **建議 Python**
2. **Q4**: Bootstrap 實作層級（C++ vs Python）→ **建議 Python**
3. **Q8**: 多維度標籤的統計策略 → **需要使用者決策**
4. **Q9**: Strand 分析的目的 → **需要使用者釐清**
5. **Q1**: 是否需要批次支援多種距離度量
6. **Q2**: 確認 C++ 輸出格式是否完整

### 建議修訂計劃的部分

1. **階段二**：明確採用 Python 聚類方案，列出具體套件與函式
2. **階段四**：細化統計分析的層次與方法選擇
3. **新增階段 0**：資料準備與格式驗證
4. **新增附錄**：參數調優計劃、測試策略

---

## 7. 總結

### 計劃的整體評價

- **清晰度**：★★★★☆（4/5）- 架構清楚但細節需補充
- **完整性**：★★★☆☆（3/5）- 缺少資料介面設計與實作優先級
- **可行性**：★★★★☆（4/5）- 技術方案可行，但需明確工具選擇
- **邏輯性**：★★★★★（5/5）- 流程邏輯嚴謹

### 關鍵決策點（需要使用者確認）

1. ✅ 採用 Python 進行聚類與統計分析（強烈建議）
2. ❓ 多維度標籤的分析策略（單獨 vs 組合）
3. ❓ Strand 的分析定位（已分離計算，是否還需作為統計維度）
4. ❓ 演化樹 vs 一般聚類的理論定位
5. ❓ Bootstrap 的計算範圍（所有位點 vs 篩選位點）

### 建議下一步

1. 根據本評估文件，與使用者討論上述關鍵問題
2. 修訂開發計劃，補充技術細節與實作優先級
3. 設計 C++ → Python 的資料介面規格
4. 開發 Proof-of-Concept：單一位點的完整分析流程（不含 Bootstrap）
