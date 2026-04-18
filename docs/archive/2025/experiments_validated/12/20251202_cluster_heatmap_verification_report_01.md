<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Cluster Heatmap 實作驗證與修正報告

**日期**：2025-12-03  
**驗證範圍**：C++ 輸出資料、Python 繪圖程式、Shell 腳本整合  
**狀態**：🔴 發現多個關鍵問題需修正

---

## 1. 執行摘要

### 1.1 驗證目標

檢查 InterSubMod 系統是否能夠：

1. ✅ C++ 輸出完整的聚類分析資料
2. ❌ 生成演化樹檔案 (.nwk格式)
3. ❌ Python 正確繪製 read-read 距離的 Cluster Heatmap
4. ✅ Shell 腳本自動整合兩階段流程

### 1.2 主要發現

| 項目 | 狀態 | 問題描述 |
|------|------|---------|
| **C++ 距離矩陣輸出** | ✅ 正常 | 完整輸出 `distance/NHD/matrix.csv` 等檔案 |
| **C++ 演化樹建構** | 🔴 **未啟用** | `HierarchicalClustering.cpp` 已實作但未在 `RegionProcessor` 中調用 |
| **Newick 樹檔案** | 🔴 **缺失** | 無 `.nwk` 檔案輸出 |
| **Python 繪圖邏輯** | 🔴 **錯誤** | 繪製的是 reads × CpGs，而非 reads × reads 距離矩陣 |
| **Dendrogram 顯示** | 🔴 **缺失** | 無演化樹視覺化 |
| **Reads 聚類排序** | 🔴 **未執行** | Reads 按原始順序排列，未依聚類結果重排 |

---

## 2. 詳細問題分析

### 2.1 C++ 端問題

#### 問題 1：聚類功能未整合到主流程

**現況**：

- ✅ `/big8_disk/liaoyoyo2001/InterSubMod/src/core/HierarchicalClustering.cpp` **已完整實作**
  - `build_upgma()`, `build_ward()`, `build_single()`, `build_complete()`
  - 支援多種連結方法
- ✅ `/big8_disk/liaoyoyo2001/InterSubMod/src/io/TreeWriter.cpp` **已實作**
  - `write_newick()` 方法可輸出 Newick 格式
- 🔴 `/big8_disk/liaoyoyo2001/InterSubMod/src/core/RegionProcessor.cpp` **未調用聚類**
  - `process_single_region()` 只計算距離矩陣（L398-496）
  - **沒有呼叫** `HierarchicalClustering::build_tree()`
  - **沒有呼叫** `TreeWriter::write_newick()`

**證據**：

```bash
# 檢查輸出目錄
$ find /big8_disk/liaoyoyo2001/InterSubMod/output -name "*.nwk" | wc -l
0  # 沒有任何 Newick 樹檔案

# 檢查 RegionProcessor.cpp
$ grep -n "build_tree" src/core/RegionProcessor.cpp
# 無結果
```

**需要修改的程式碼位置**：

- `src/core/RegionProcessor.cpp`：L496 之後（距離矩陣計算完成後）
- 需要加入聚類建樹與輸出邏輯

---

#### 問題 2：缺少 Config 參數控制聚類功能

**現況**：

- `include/core/Config.hpp` 已有部分聚類相關參數（如 `linkage_method`）
- 但缺少以下關鍵參數：
  - `bool compute_clustering`：是否執行聚類
  - `bool output_tree_files`：是否輸出 Newick 樹
  - `int clustering_min_reads`：最小 Read 數閾值

**需要修改的檔案**：

- `include/core/Config.hpp`
- `src/utils/ArgParser.cpp`（新增命令列參數）

---

### 2.2 Python 端問題

#### 問題 3：繪製的不是距離矩陣 Heatmap

**當前 `plot_cluster_heatmap.py` 的行為**：

```python
# L276-348: plot_cluster_heatmap() 函式
def plot_cluster_heatmap(
    meth_matrix: pd.DataFrame,  # ← 甲基化矩陣 (Reads × CpGs)
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    ...
):
    # L310: 按聚類順序重排
    meth_ordered = meth_matrix.iloc[order]
    
    # L335: 繪製 seaborn clustermap
    g = sns.clustermap(
        data_filled,  # ← 這是 Reads × CpGs 的甲基化矩陣
        row_cluster=False,  # 已排序，不重複聚類
        col_cluster=False,  # CpG 位點不聚類
        ...
    )
```

**問題所在**：

1. **輸入數據錯誤**：使用的是 `methylation_matrix.csv`（Reads × CpGs）
   - 正確應該使用：`distance/NHD/matrix.csv`（Reads × Reads）
2. **缺少 Dendrogram**：
   - `row_cluster=False` 關閉了行聚類
   - 即使提供了 `linkage_matrix`，也沒有用它來繪製 dendrogram
3. **Y 軸是 Reads，X 軸是 CpG 位點**：
   - 正確的 Cluster Heatmap 應該是：Y 軸 = Reads，X 軸 = Reads
   - 顏色 = 距離值

**預期的正確行為**：

應該繪製的是類似這樣的圖：

```
        Read1  Read2  Read3  Read4  ...
Read1    0.0    0.2    0.8    0.9
Read2    0.2    0.0    0.7    0.85  ← 距離矩陣
Read3    0.8    0.7    0.0    0.3
Read4    0.9    0.85   0.3    0.0
...
```

並在左側/上方顯示 Dendrogram（演化樹）

---

#### 問題 4：與文檔說明不符

**文檔中的描述** (`clustering_evolution_implementation_guide.md` L69-73)：

```python
**實作細節**：
- **Heatmap**：顯示 Read (Y軸) x CpG Site (X軸) 的甲基化狀態。  # ← 錯誤
- **Dendrogram**：顯示在 Heatmap 側邊，呈現聚類結構。
- **Annotation Bars**：在 Heatmap 旁添加顏色條，標示 Read 的 Tag 資訊...
```

**實際應該是**：

- **❌ Methylation Heatmap**（甲基化熱圖）：Reads × CpGs
  - 這是原始資料視覺化，不是聚類分析結果
- **✅ Cluster Hea

tmap**（聚類熱圖）：Reads × Reads（距離矩陣）

- 這才是聚類結果的視覺化

**結論**：文檔描述混淆了兩種不同的 Heatmap 類型

---

### 2.3 整合問題

#### 問題 5：Shell 腳本假設資料已存在

**`run_vcf_all_snv.sh` L281-289**：

```bash
python3 "${PLOT_SCRIPT}" \
    --output-dir "${OUTPUT_DIR}" \
    --threads "${PLOT_THREADS}" \
    --metric "${FIRST_METRIC}" \
    --linkage average \  # ← Python 重新計算聚類
    ...
```

**問題**：

1. Python 從距離矩陣重新計算聚類（`compute_linkage()` L184-209）
2. C++ 已有聚類實作，但未使用
3. 重複計算浪費資源，且可能產生不一致的結果

**理想流程**：

1. C++ 計算距離矩陣 **並執行聚類**
2. C++ 輸出演化樹 (.nwk) 與排序後的 Read 順序
3. Python 讀取樹檔案，直接繪圖（不重新聚類）

---

## 3. 資料檔案完整性檢查

### 3.1 C++ 當前輸出檔案清單

以 `/output/20251202_vcf_all_w1000/chr1_877772/chr1_876772_878772/` 為例：

```
✅ metadata.txt                           # 區域元資訊
✅ reads/reads.tsv                        # Read 標籤 (HP, Strand, etc.)
✅ methylation/methylation.csv            # 甲基化矩陣
✅ methylation/methylation_forward.csv    # Forward strand
✅ methylation/methylation_reverse.csv    # Reverse strand
✅ methylation/cpg_sites.tsv              # CpG 位點列表
✅ distance/NHD/matrix.csv                # 完整距離矩陣
✅ distance/NHD/matrix_forward.csv        # Forward strand 距離矩陣
✅ distance/NHD/matrix_reverse.csv        # Reverse strand 距離矩陣
✅ distance/NHD/stats.txt                 # 距離統計資訊
❌ tree.nwk                                # 演化樹 (缺失！)
❌ tree_forward.nwk                        # Forward 演化樹 (缺失！)
❌ tree_reverse.nwk                        # Reverse 演化樹 (缺失！)
❌ clustered_read_order.txt                # 聚類後的 Read 順序 (缺失！)
```

### 3.2 需要新增的資料輸出

| 檔案名稱 | 格式 | 用途 | 優先級 |
|---------|------|------|--------|
| `tree.nwk` | Newick | 完整演化樹（所有 reads） | ⭐⭐⭐⭐⭐ |
| `tree_forward.nwk` | Newick | Forward strand 演化樹 | ⭐⭐⭐⭐ |
| `tree_reverse.nwk` | Newick | Reverse strand 演化樹 | ⭐⭐⭐⭐ |
| `clustering/read_order.txt` | Text | 聚類後的 Read ID 順序 | ⭐⭐⭐ |
| `clustering/clusters.tsv` | TSV | 切割樹後的 Cluster 標籤 | ⭐⭐ |

**Newick 格式範例**：

```
((Read_1:0.05,Read_2:0.05)95:0.03,(Read_3:0.08,Read_4:0.08)80:0.02)100:0.0;
```

- `Read_X` = 葉節點（實際 Read ID）
- `:0.05` = 分支長度
- `)95` = Bootstrap 支持度（若有執行 Bootstrap）

---

## 4. Python 腳本問題詳細分析

### 4.1 當前執行結果

**已生成的圖片**：

```bash
$ find /big8_disk/liaoyoyo2001/InterSubMod/output -name "cluster_heatmap.png" | wc -l
# 已生成數千張圖片
```

**問題**：這些圖片顯示的是什麼？

讓我們追蹤程式碼執行流程：

1. **L426**: 載入 `methylation_matrix.csv`

   ```python
   meth_df, cpg_positions = load_methylation_matrix(region_dir, strand)
   ```

   - `meth_df` = Reads × CpGs 矩陣（例如：80 reads × 30 CpGs）

2. **L443**: 載入距離矩陣

   ```python
   dist_matrix = load_distance_matrix(region_dir, distance_metric, strand)
   ```

   - `dist_matrix` = Reads × Reads 距離矩陣（例如：80 × 80）

3. **L452**: 從距離矩陣計算聚類

   ```python
   Z = compute_linkage(dist_matrix, method=linkage_method)
   ```

   - `Z` = Linkage matrix（用於確定樹的拓撲）

4. **L477**: 繪製 Heatmap

   ```python
   success = plot_cluster_heatmap(
       meth_df,  # ← 傳入的是甲基化矩陣！
       reads_df, Z, output_path, ...
   )
   ```

5. **L335 (plot_cluster_heatmap內部)**: Seaborn clustermap

   ```python
   g = sns.clustermap(
       data_filled,  # ← data_filled 來自 meth_ordered (甲基化矩陣)
       row_cluster=False,
       col_cluster=False,
       ...
   )
   ```

**結論**：

- ✅ **聚類計算是正確的**（使用距離矩陣）
- ❌ **視覺化是錯誤的**（繪製甲基化矩陣）
- ❌ **缺少 Dendrogram**（`row_cluster=False`）

### 4.2 正確的繪圖邏輯應該是什麼？

**選項 A：繪製距離矩陣 Heatmap（標準做法）**

```python
def plot_distance_heatmap(
    dist_matrix: np.ndarray,
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    ...
):
    # 1. 按聚類順序重排距離矩陣
    order = leaves_list(linkage_matrix)
    dist_ordered = dist_matrix[order, :][:, order]
    
    # 2. 建立 DataFrame
    read_ids_ordered = [read_ids[i] for i in order]
    dist_df = pd.DataFrame(
        dist_ordered,
        index=read_ids_ordered,
        columns=read_ids_ordered
    )
    
    # 3. 繪製 Clustermap
    g = sns.clustermap(
        dist_df,
        row_linkage=linkage_matrix,  # ← 使用預先計算的 linkage
        col_linkage=linkage_matrix,  # ← 同樣的 linkage（對稱矩陣）
        row_colors=annotation_colors,
        cmap='viridis',
        vmin=0, vmax=1,
        ...
    )
```

**輸出效果**：

- X 軸、Y 軸都是 Reads
- 左側與上方顯示 Dendrogram
- Heatmap 顏色表示距離（0 = 相似，1 = 不同）

---

**選項 B：繪製甲基化矩陣 + Dendrogram（也可接受）**

```python
def plot_methylation_with_dendrogram(
    meth_matrix: pd.DataFrame,
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    ...
):
    # 按聚類順序重排
    order = leaves_list(linkage_matrix)
    meth_ordered = meth_matrix.iloc[order]
    
    # 繪製 Clustermap
    g = sns.clustermap(
        meth_ordered,
        row_linkage=linkage_matrix,  # ← 提供預先計算的 linkage
        row_cluster=False,  # ← 已排序，不重複聚類（但仍顯示 dendrogram）
        col_cluster=False,  # CpG 位點不聚類
        dendrogram_ratio=(0.15, 0),  # ← 顯示行 dendrogram
        ...
    )
```

**問題**：`seaborn.clustermap` 在 `row_cluster=False` 時，**不會顯示 dendrogram**

**解決方案**：需要手動繪製 dendrogram

---

## 5. 文檔錯誤檢查

### 5.1 `clustering_evolution_implementation_guide.md`

#### 錯誤 1：Heatmap 定義混淆（L69-73）

**文檔原文**：

```markdown
**實作細節**：
- **Heatmap**：顯示 Read (Y軸) x CpG Site (X軸) 的甲基化狀態。
- **Dendrogram**：顯示在 Heatmap 側邊，呈現聚類結構。
- **Annotation Bars**：在 Heatmap 旁添加顏色條...
```

**問題**：

- 這描述的是 **Methylation Heatmap**，不是 **Cluster Heatmap**
- Cluster Heatmap 應該顯示 Read × Read 的距離或相似性

**建議修正**：

```markdown
**兩種視覺化類型**：

1. **Distance Heatmap（距離熱圖）**：
   - X 軸、Y 軸：Reads
   - 顏色：Read 間的距離值
   - Dendrogram：顯示在左側與上方
   - 用途：直接視覺化聚類結果

2. **Methylation Heatmap（甲基化熱圖）+ Dendrogram**：
   - Y 軸：Reads（按聚類順序排列）
   - X 軸：CpG Sites
   - 顏色：甲基化程度（0-1）
   - Dendrogram：僅顯示在左側
   - 用途：展示聚類後的甲基化模式
```

---

#### 錯誤 2：Python 程式碼範例（L415 onwards）

**文檔提供的範例**：

```python
def plot_clustermap(...):
    # 讀取資料
    meth_df = pd.read_csv(methylation_matrix_path, index_col=0)
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    
    # ... 建立 annotation colors ...
    
    # 繪製 Clustermap（使用預先計算的樹）
    g = sns.clustermap(
        meth_df,  # ← 甲基化矩陣
        row_colors=row_colors,
        method='average',  # UPGMA
        cmap='RdBu_r',
        ...
    )
```

**問題**：

1. 沒有使用預先計算的距離矩陣
2. 沒有使用預先計算的演化樹 (.nwk)
3. `method='average'` 會重新計算聚類

**應該修正為**：

```python
def plot_distance_heatmap(...):
    # 讀取距離矩陣
    dist_df = pd.read_csv(distance_matrix_path, index_col=0)
    
    # 從距離矩陣計算 linkage
    Z = linkage(squareform(dist_df.values), method='average')
    
    # 繪製
    g = sns.clustermap(
        dist_df,
        row_linkage=Z,
        col_linkage=Z,
        ...
    )
```

---

### 5.2 `clustering_evolution_research_overview.md`

#### 潛在混淆：未明確區分兩種視覺化

文檔中提到「Cluster Heatmap」但未明確說明是距離矩陣還是甲基化矩陣。

**建議**：在第 4.1 節明確說明。

---

### 5.3 `heatmap_generation_strategy_analysis.md`

#### 無明顯錯誤

該文檔主要討論生成策略（即時 vs 批次），邏輯正確。

---

## 6. 修正方案

### 6.1 C++ 端修正（優先級：⭐⭐⭐⭐⭐）

#### 修正 1：RegionProcessor 加入聚類建樹

**檔案**：`src/core/RegionProcessor.cpp`

**位置**：`process_single_region()` 方法，L496 之後（距離矩陣輸出完成後）

**新增程式碼**：

```cpp
// 已存在：L473-483 輸出距離矩陣

// ==== 新增：執行聚類建樹 ====
if (config_.compute_clustering && result.num_reads >= config_.clustering_min_reads) {
    HierarchicalClustering::ClusteringConfig cluster_config;
    cluster_config.method = config_.linkage_method;  // UPGMA, WARD, etc.
    
    HierarchicalClustering clusterer(cluster_config);
    
    // 準備 Read IDs
    std::vector<std::string> read_ids;
    for (const auto& r : read_list) {
        read_ids.push_back(r.read_id);
    }
    
    // 建構演化樹（所有 reads）
    Tree tree_all = clusterer.build_tree(all_dist, read_ids);
    
    // 輸出 Newick 格式
    TreeWriter tree_writer;
    tree_writer.write_newick(tree_all, region_dir + "/tree.nwk");
    
    // Strand-specific 樹（若啟用）
    if (output_strand_distance_matrices_) {
        // Forward strand
        if (forward_dist.size() > 0) {
            std::vector<std::string> fwd_ids;
            for (const auto& r : read_list) {
                if (r.strand == Strand::FORWARD) {
                    fwd_ids.push_back(r.read_id);
                }
            }
            if (fwd_ids.size() >= 2) {
                Tree tree_fwd = clusterer.build_tree(forward_dist, fwd_ids);
                tree_writer.write_newick(tree_fwd, region_dir + "/tree_forward.nwk");
            }
        }
        
        // Reverse strand
        if (reverse_dist.size() > 0) {
            std::vector<std::string> rev_ids;
            for (const auto& r : read_list) {
                if (r.strand == Strand::REVERSE) {
                    rev_ids.push_back(r.read_id);
                }
            }
            if (rev_ids.size() >= 2) {
                Tree tree_rev = clusterer.build_tree(reverse_dist, rev_ids);
                tree_writer.write_newick(tree_rev, region_dir + "/tree_reverse.nwk");
            }
        }
    }
    
    if (log_level_ >= LogLevel::LOG_DEBUG) {
        #pragma omp critical
        {
            std::cout << "  Clustering tree built with " 
                      << tree_all.num_leaves() << " leaves" << std::endl;
        }
    }
}
```

---

#### 修正 2：Config 新增參數

**檔案**：`include/core/Config.hpp`

**新增成員變數**：

```cpp
struct Config {
    // ... 現有參數 ...
    
    // === Clustering Configuration ===
    bool compute_clustering = true;               ///< 是否執行聚類
    bool output_tree_files = true;                ///< 是否輸出 Newick 樹檔案
    LinkageMethod linkage_method = LinkageMethod::UPGMA;  ///< 連結方法
    int clustering_min_reads = 10;                ///< 最小 Read 數閾值
};
```

**檔案**：`src/utils/ArgParser.cpp`（需確認是否存在此檔案）

**新增命令列參數**：

```cpp
app.add_flag("--compute-clustering", config.compute_clustering,
    "Enable hierarchical clustering analysis");

app.add_option("--linkage-method", config.linkage_method,
    "Linkage method: UPGMA, WARD, SINGLE, COMPLETE")
    ->transform(CLI::CheckedTransformer(linkage_map, CLI::ignore_case));

app.add_option("--clustering-min-reads", config.clustering_min_reads,
    "Minimum reads required for clustering");
```

---

### 6.2 Python 端修正（優先級：⭐⭐⭐⭐⭐）

#### 修正 3：修復 plot_cluster_heatmap.py

**問題檔案**：`/big8_disk/liaoyoyo2001/InterSubMod/tools/plot_cluster_heatmap.py`

**修正方案**：提供兩種繪圖模式

**新增函式 1：繪製距離矩陣 Heatmap**

```python
def plot_distance_heatmap(
    dist_matrix: np.ndarray,
    read_ids: List[str],
    reads_df: pd.DataFrame,
    linkage_matrix: np.ndarray,
    output_path: str,
    region_info: Dict,
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 150
) -\u003e bool:
    """
    繪製 Read-Read 距離矩陣 Heatmap
    
    這是真正的 Cluster Heatmap：
    - X 軸、Y 軸都是 Reads
    - 顏色表示距離
    - 左側與上方顯示 Dendrogram
    """
    try:
        # 1. 按聚類順序重排
        order = get_cluster_order(linkage_matrix)
        dist_ordered = dist_matrix[order, :][:, order]
        read_ids_ordered = [read_ids[i] for i in order]
        reads_ordered = reads_df.loc[read_ids_ordered]
        
        # 2. 建立 DataFrame
        dist_df = pd.DataFrame(
            dist_ordered,
            index=read_ids_ordered,
            columns=read_ids_ordered
        )
        
        # 3. 建立 annotation
        annotations, color_dict = create_annotation_colors(reads_ordered)
        row_colors = None
        if not annotations.empty:
            row_colors_list = []
            for col in annotations.columns:
                row_colors_list.append(annotations[col].map(color_dict[col]))
            row_colors = pd.concat(row_colors_list, axis=1)
            row_colors.columns = annotations.columns
        
        # 4. 繪製 Clustermap
        g = sns.clustermap(
            dist_df,
            row_linkage=linkage_matrix,
            col_linkage=linkage_matrix,
            row_colors=row_colors,
            col_colors=row_colors,
            cmap='viridis',  # 距離用不同配色
            vmin=0, vmax=1,
            figsize=figsize,
            cbar_kws={'label': 'Distance'},
            xticklabels=False,  # 太多 reads，不顯示標籤
            yticklabels=False,
        )
        
        # 5. 標題
        title = f\"Read-Read Distance Heatmap: {region_info.get('snv', 'Unknown')}\\n\"
        title += f\"Reads: {len(read_ids)}, Linkage: UPGMA\"
        g.fig.suptitle(title, y=0.98, fontsize=12, fontweight='bold')
        
        # 6. 儲存
        g.savefig(output_path, dpi=dpi, bbox_inches='tight')
        plt.close(g.fig)
        
        return True
    except Exception as e:
        print(f\"Error creating distance heatmap: {e}\")
        return False
```

**新增函式 2：讀取 Newick 樹**

```python
def load_newick_tree(region_dir: str) -\u003e Optional[str]:
    """
    載入 Newick 格式演化樹檔案
    """
    tree_file = os.path.join(region_dir, \"tree.nwk\")
    if not os.path.exists(tree_file):
        return None
    
    try:
        with open(tree_file, 'r') as f:
            return f.read().strip()
    except:
        return None
```

**修改主函式 process_single_region**：

```python
def process_single_region(..., plot_type: str = \"distance\"):
    # ... 現有載入邏輯 ...
    
    # 載入距離矩陣
    dist_matrix = load_distance_matrix(region_dir, distance_metric, strand)
    
    # 嘗試載入 Newick 樹
    newick_tree = load_newick_tree(region_dir)
    
    # 計算 linkage（若無預先計算的樹）
    if newick_tree is None:
        Z = compute_linkage(dist_matrix, method=linkage_method)
    else:
        # TODO: 從 Newick 解析 linkage matrix
        # 目前仍需重新計算
        Z = compute_linkage(dist_matrix, method=linkage_method)
    
    # 選擇繪圖類型
    if plot_type == \"distance\":
        # 繪製距離矩陣 Heatmap
        read_ids = list(meth_df.index)
        success = plot_distance_heatmap(
            dist_matrix, read_ids, reads_df, Z,
            output_path, region_info, figsize=figsize, dpi=dpi
        )
    else:
        # 繪製甲基化矩陣 Heatmap（現有功能）
        success = plot_cluster_heatmap(
            meth_df, reads_df, Z, output_path,
            region_info, figsize=figsize, dpi=dpi
        )
    
    return region_dir, success, output_path if success else \"Failed\"
```

---

### 6.3 Shell 腳本修正（優先級：⭐⭐⭐）

**檔案**：`scripts/run_vcf_all_snv.sh`

**修改點 1**：加入繪圖類型參數（L288）

```bash
python3 "${PLOT_SCRIPT}" \
    --output-dir "${OUTPUT_DIR}" \
    --threads "${PLOT_THREADS}" \
    --metric "${FIRST_METRIC}" \
    --linkage average \
    --min-reads 10 \
    --min-cpgs 3 \
    --format png \
    --dpi 150 \
    --plot-type distance  # ← 新增：指定繪製距離矩陣
```

**修改點 2**：更新幫助訊息（L90-92）

```bash
echo \"Cluster Heatmap:\"
echo \"  By default, distance-based cluster heatmaps are generated.\"
echo \"  Use --no-plots to skip this step.\"
```

---

## 7. 實作優先級與時程

### 7.1 立即執行（本週內）

| 任務 | 檔案 | 預估時間 | 優先級 |
|------|------|---------|--------|
| **1. C++ 整合聚類** | `RegionProcessor.cpp` | 2 小時 | ⭐⭐⭐⭐⭐ |
| **2. Config 參數** | `Config.hpp` | 30 分鐘 | ⭐⭐⭐⭐⭐ |
| **3. Python 距離 Heatmap** | `plot_cluster_heatmap.py` | 3 小時 | ⭐⭐⭐⭐⭐ |
| **4. 測試驗證** | 執行小規模測試 | 1 小時 | ⭐⭐⭐⭐⭐ |

---

### 7.2 次要優化（下週）

| 任務 | 說明 | 預估時間 | 優先級 |
|------|------|---------|--------|
| **5. Bootstrap 整合** | C++ Bootstrap 分析 | 4 小時 | ⭐⭐⭐⭐ |
| **6. Newick 解析** | Python 讀取樹檔案 | 2 小時 | ⭐⭐⭐ |
| **7. 文檔修正** | 更新所有文檔 | 2 小時 | ⭐⭐⭐ |

---

## 8. 測試驗證計劃

### 8.1 單元測試

**C++ 測試**：

```cpp
// tests/test_clustering_integration.cpp
TEST(RegionProcessor, ClusteringWithTreeOutput) {
    // 使用測試資料
    Config config;
    config.compute_clustering = true;
    config.output_tree_files = true;
    
    RegionProcessor processor(config);
    // ... 處理一個測試位點 ...
    
    // 驗證樹檔案存在
    EXPECT_TRUE(std::filesystem::exists(output_dir + \"/tree.nwk\"));
    
    // 驗證 Newick 格式正確
    std::ifstream tree_file(output_dir + \"/tree.nwk\");
    std::string newick;
    std::getline(tree_file, newick);
    EXPECT_TRUE(newick.back() == ';');  // Newick 以分號結束
}
```

**Python 測試**：

```python
# tests/test_distance_heatmap.py
def test_distance_heatmap_generation():
    # 建立測試資料
    dist_matrix = np.array([...])
    read_ids = [\"Read1\", \"Read2\", \"Read3\"]
    
    # 繪製
    success = plot_distance_heatmap(...)
    
    assert success
    assert os.path.exists(output_path)
    
    # 檢查圖片尺寸
    img = Image.open(output_path)
    assert img.size[0] > 0
```

---

### 8.2 整合測試

**小規模測試**（100 位點）：

```bash
# 1. 執行 C++ 處理
./build/bin/inter_sub_mod \
    --vcf test_data/100_snvs.vcf \
    --output-dir test_output \
    --threads 8 \
    --compute-clustering

# 2. 檢查輸出
find test_output -name \"*.nwk\" | wc -l  # 應有 100 個樹檔案

# 3. 執行 Python 繪圖
python3 tools/plot_cluster_heatmap.py \
    --output-dir test_output \
    --threads 4 \
    --plot-type distance

# 4. 檢查圖片
find test_output -name \"cluster_heatmap.png\" | wc -l  # 應有 ~100 張圖
```

---

### 8.3 視覺驗證

**檢查清單**：

1. ✅ Heatmap 是否為方形（N × N）？
2. ✅ 對角線是否為深色（距離 = 0）？
3. ✅ 左側與上方是否有 Dendrogram？
4. ✅ Annotation bars 是否正確顯示？
5. ✅ 顏色映射是否合理（近似 reads 顏色相近）？

---

## 9. 預期效能

### 9.1 C++ 聚類增加的時間

**單一位點**（100 reads）：

- 距離矩陣計算：1 ms（已有）
- UPGMA 聚類：**~2 ms**（新增）
- Newick 輸出：**< 0.5 ms**（新增）
- **總增加**：~2.5 ms

**40,000 位點**：

- 總增加時間：40,000 × 2.5 ms = **100 秒**（約 1.7 分鐘）
- 相對原本 10 分鐘的執行時間，增幅約 **+17%**

**結論**：影響可接受

---

### 9.2 Python 繪圖時間

**距離 Heatmap** vs **甲基化 Heatmap**：

| 指標 | 甲基化 Heatmap | 距離 Heatmap |
|------|---------------|-------------|
| 矩陣大小 | 100 × 50 | 100 × 100 |
| 繪圖時間 | ~300 ms | ~400 ms |

**預估**（10,000 位點，16 執行緒）：

- 總時間：10,000 / 16 × 0.4s = **250 秒**（約 4 分鐘）

---

## 10. 結論與建議

### 10.1 關鍵問題匯總

| 問題 | 嚴重性 | 狀態 | 修復時間 |
|------|--------|------|---------|
| C++ 未呼叫聚類功能 | 🔴 高 | 待修復 | 2 小時 |
| 無 Newick 樹檔案輸出 | 🔴 高 | 待修復 | 30 分鐘 |
| Python 繪製錯誤的 Heatmap | 🔴 高 | 待修復 | 3 小時 |
| 文檔描述不準確 | 🟡 中 | 待更新 | 2 小時 |
| 缺少 Dendrogram 顯示 | 🟡 中 | 待修復 | 1 小時 |

---

### 10.2 立即行動項

**第一優先**（今日完成）：

1. 修改 `RegionProcessor.cpp`，加入聚類建樹
2. 修改 `Config.hpp`，加入 clustering 參數
3. 編譯測試，確認無錯誤

**第二優先**（明日完成）：

1. 修改 `plot_cluster_heatmap.py`，新增距離 Heatmap 功能
2. 執行小規模測試（100 位點）
3. 視覺驗證圖片正確性

**第三優先**（本週內）：

1. 更新文檔，修正錯誤描述
2. 執行大規模測試（40,000 位點）
3. 撰寫完整測試報告

---

### 10.3 長期改進

1. **Bootstrap 支持度**：在演化樹上標註 Bootstrap 值
2. **互動式視覺化**：使用 Plotly 生成可縮放的 HTML 圖片
3. **自動化報告**：整合統計分析結果到圖片中

---

**文件狀態**：✅ 驗證完成  
**下一步**：開始實作修正方案  
**預計完成時間**：2025-12-05
