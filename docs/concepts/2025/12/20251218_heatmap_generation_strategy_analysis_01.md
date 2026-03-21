<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Cluster Heatmap 生成策略分析與建議

**日期**：2025-12-02  
**目的**：評估三種 Cluster Heatmap 生成方案的優劣

---

## 1. 現有程式碼架構分析

### 1.1 目前已實作的功能

根據 `/big8_disk/liaoyoyo2001/InterSubMod/src` 的程式碼：

✅ **已完成**：

- `HierarchicalClustering.cpp`：C++ 層次聚類實作（UPGMA, Ward, Single, Complete）
- `TreeWriter.cpp`：Newick 格式樹輸出
- `DistanceMatrix.cpp`：距離矩陣計算
- `RegionWriter.cpp`：各類資料輸出（reads, CpG sites, distance matrices）
- `RegionProcessor.cpp`：主要處理流程，平行化處理多個位點

### 1.2 當前資料流程

```
RegionProcessor::process_single_region() (每個位點獨立處理)
  ↓
[1] 讀取 Reads、解析甲基化
  ↓
[2] MatrixBuilder 建構甲基化矩陣
  ↓
[3] RegionWriter::write_region()  ← 輸出 methylation_matrix.csv, reads.tsv
  ↓
[4] DistanceCalculator::compute()  ← 計算距離矩陣
  ↓
[5] RegionWriter::write_distance_matrices()  ← 輸出 distance_matrix.csv
  ↓
[目前未實作] 聚類 + Heatmap 生成
```

**關鍵觀察**：

- 採用 **OpenMP 平行化**（L208-261），每個執行緒處理一個位點
- 每個執行緒有獨立的 `BamReader` 和 `FastaReader`（L214-220）
- 處理完單一位點立即輸出結果（L378-396）

---

## 2. 三種方案評估

### 方案 A：每個位點處理時立即呼叫 Python 繪製

**實作方式**：

```cpp
// 在 RegionProcessor::process_single_region() 的最後
if (config_.generate_heatmaps && result.num_reads >= 20) {
    // 呼叫 Python 腳本
    std::ostringstream cmd;
    cmd << config_.python_executable << " scripts/visualization/plot_clustermap.py"
        << " --methylation " << region_dir << "/methylation_matrix.csv"
        << " --metadata " << region_dir << "/read_metadata.tsv"
        << " --tree " << region_dir << "/tree.nwk"
        << " --output " << region_dir << "/heatmap.png";
    
    int ret = std::system(cmd.str().c_str());
}
```

#### 優點 ✅

1. **即時視覺化**：處理完立即有圖片，方便除錯與監控
2. **記憶體友善**：每個位點處理完即釋放記憶體，不累積
3. **失敗隔離**：某個位點的繪圖失敗不影響其他位點
4. **易於實作**：直接在現有流程末端加入一行呼叫

#### 缺點 ❌

1. **效率低**：
   - 每次呼叫 Python 都需要啟動直譯器（~100-500ms 開銷）
   - 40,000 位點 × 500ms = **5.5 小時**額外開銷
   - 平行化困難：`system()` 會阻塞當前執行緒
2. **資源競爭**：
   - 64 個執行緒同時呼叫 Python → 64 個 Python 行程
   - 每個 Python 行程可能使用 500MB-1GB 記憶體
   - 總記憶體 > 32GB，可能 OOM
3. **I/O 瓶頸**：
   - 大量程序同時讀寫檔案 → 磁碟 I/O 飽和
   - 影響主程式的資料輸出效能

#### 適用情境

- **小規模測試**（< 100 位點）
- **單執行緒模式**（方便除錯）
- **需要即時監控**的情境

---

### 方案 B：C++ 先輸出所有資料，完成後批次執行 Python

**實作方式**：

```bash
# Step 1: C++ 處理所有位點，輸出資料
./inter_sub_mod ... --output-dir output/

# Step 2: 使用 Python 批次腳本
python3 scripts/visualization/batch_plot_heatmaps.py \
    --input-pattern "output/single_locus/*" \
    --n-jobs 16
```

Python 批次腳本：

```python
from multiprocessing import Pool
import glob

def plot_single_locus(locus_dir):
    # 讀取資料並繪製
    plot_clustermap(...)

locus_dirs = glob.glob("output/single_locus/*")

with Pool(16) as pool:
    pool.map(plot_single_locus, locus_dirs)
```

#### 優點 ✅

1. **高效率**：
   - C++ 全速執行，不被 Python 拖慢
   - Python 批次處理可使用 `multiprocessing`，充分利用多核
   - 預估時間：C++ 10 分鐘 + Python 1 小時 = **總計 70 分鐘**
2. **資源控制**：
   - 可精確控制平行數（如 16 個 Python worker）
   - 避免記憶體爆炸
3. **靈活調整**：
   - 可重複繪圖而無需重新執行 C++
   - 可調整繪圖參數、配色、篩選條件

#### 缺點❌

1. **兩階段流程**：需要手動執行兩次（或用腳本自動化）
2. **延遲反饋**：需等待 C++ 完成才能看到圖片
3. **重複 I/O**：Python 需重新讀取 C++ 已輸出的檔案

#### 適用情境 ⭐⭐⭐⭐⭐

- **生產環境**（40,000 位點）
- **需要調整繪圖參數**的實驗性分析
- **記憶體受限**的環境

---

### 方案 C：直接輸出各自的熱圖與演化樹（PNG圖片）

**實作方式**：C++ 使用圖形函式庫（如 matplotlib-cpp 或 SVG 生成）直接繪製圖片

```cpp
#include "matplotlibcpp.h"

void RegionProcessor::plot_heatmap(...) {
    namespace plt = matplotlibcpp;
    
    // 繪製 Heatmap
    plt::imshow(...);
    plt::savefig(output_path);
}
```

#### 優點 ✅

1. **單一流程**：C++ 一次性完成所有工作
2. **無 Python 依賴**：減少外部依賴

#### 缺點 ❌❌❌

1. **開發成本極高**：
   - `matplotlib-cpp` 仍需 Python 環境（只是包裝）
   - 純 C++ 繪製 dendrogram + heatmap + annotations 需數千行程式碼
2. **美觀度差**：
   - C++ 圖形庫遠不如 Seaborn 精緻
   - 顏色映射、字體、佈局調整極其困難
3. **維護困難**：
   - 調整圖片樣式需修改 C++ 並重新編譯
   - 無法快速測試不同視覺化參數

#### 適用情境

- **非常罕見**：除非完全無 Python 環境且不允許安裝

---

## 3. 綜合建議與最佳實踐

### 🏆 推薦方案：**方案 B（兩階段批次處理）**

**理由**：

1. **效能最優**：C++ 與 Python 各司其職，無相互阻塞
2. **資源可控**：避免記憶體與 I/O 競爭
3. **靈活性高**：可重複繪圖、調整參數

---

### 📋 具體實作方案

#### 階段一：C++ 輸出完整資料

**修改 `RegionProcessor::process_single_region()`**：

```cpp
// 在距離矩陣輸出後（L483 之後）

// 執行聚類建樹（若啟用）
if (config_.compute_clustering && result.num_reads >= config_.clustering_min_reads) {
    HierarchicalClustering clusterer(config_.linkage_method);
    
    // 建構演化樹
    std::vector<std::string> read_ids;
    for (const auto& r : read_list) {
        read_ids.push_back(r.read_id);
    }
    
    Tree tree = clusterer.build_tree(all_dist, read_ids);
    
    // 輸出 Newick 樹
    TreeWriter tree_writer;
    tree_writer.write_newick(tree, region_dir + "/tree.nwk");
    
    // 若有 strand-specific 矩陣，也建樹
    if (output_strand_distance_matrices_ && forward_dist.size() > 0) {
        // 提取 forward reads
        std::vector<std::string> fwd_ids;
        for (const auto& r : read_list) {
            if (r.strand == Strand::FORWARD) {
                fwd_ids.push_back(r.read_id);
            }
        }
        if (fwd_ids.size() >= 2) {
            Tree fwd_tree = clusterer.build_tree(forward_dist, fwd_ids);
            tree_writer.write_newick(fwd_tree, region_dir + "/tree_forward.nwk");
        }
        
        // 同理處理 reverse
        // ...
    }
}
```

**輸出資料清單**：

- ✅ `methylation_matrix.csv`（已有）
- ✅ `distance_matrix.csv`（已有）
- ✅ `read_metadata.tsv`（已有：包含 read_id, hp_tag, tumor_normal, strand, alt_ref）
- ✅ `cpg_sites.tsv`（已有）
- ✅ `tree.nwk`（新增）
- ✅ `tree_forward.nwk`, `tree_reverse.nwk`（新增）

---

#### 階段二：Python 批次繪圖與統計

**主控腳本**：`scripts/run_analysis_pipeline.py`

```python
#!/usr/bin/env python3
"""
完整分析流程：批次繪圖 + 統計分析
"""

import argparse
import glob
from multiprocessing import Pool
import os

def plot_locus(locus_dir):
    """繪製單一位點的 Heatmap"""
    from visualization.plot_clustermap import plot_clustermap
    
    methylation_csv = os.path.join(locus_dir, 'methylation_matrix.csv')
    metadata_tsv = os.path.join(locus_dir, 'read_metadata.tsv')
    tree_nwk = os.path.join(locus_dir, 'tree.nwk')
    output_png = os.path.join(locus_dir, 'heatmap.png')
    
    if not os.path.exists(methylation_csv):
        return None
    
    try:
        plot_clustermap(methylation_csv, metadata_tsv, output_png, tree_nwk)
        return locus_dir
    except Exception as e:
        print(f"Failed to plot {locus_dir}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', required=True, help='Output dir from C++')
    parser.add_argument('--n-jobs', type=int, default=16, help='Number of parallel jobs')
    parser.add_argument('--plot-only-significant', action='store_true')
    args = parser.parse_args()
    
    # Step 1: 收集所有位點目錄
    pattern = os.path.join(args.input_dir, 'single_locus', '*')
    locus_dirs = glob.glob(pattern)
    print(f"Found {len(locus_dirs)} loci")
    
    # Step 2: 篩選（若需要）
    if args.plot_only_significant:
        # TODO: 先執行統計分析，僅繪製顯著位點
        pass
    
    # Step 3: 平行繪製 Heatmap
    print(f"Plotting heatmaps with {args.n_jobs} workers...")
    with Pool(args.n_jobs) as pool:
        results = pool.map(plot_locus, locus_dirs)
    
    success_count = sum(1 for r in results if r is not None)
    print(f"Successfully plotted {success_count}/{len(locus_dirs)} heatmaps")
    
    # Step 4: 執行統計分析
    print("Running statistical analysis...")
    from statistics.run_permanova import batch_permanova
    batch_permanova(
        locus_dirs,
        os.path.join(args.input_dir, 'aggregate', 'permanova_results.tsv')
    )
    
    print("Analysis pipeline completed!")

if __name__ == "__main__":
    main()
```

**使用方式**：

```bash
# C++ 處理
./inter_sub_mod --vcf data.vcf ... --output-dir output/

# Python 分析（自動批次繪圖 + 統計）
python3 scripts/run_analysis_pipeline.py \
    --input-dir output/ \
    --n-jobs 16 \
    --plot-only-significant  # 可選：僅繪製顯著位點
```

---

### 🔧 優化策略

#### 1. 篩選機制（減少無效繪圖）

**在 C++ 側預先標記**：

```cpp
// 在 RegionProcessor 中加入簡單篩選
bool RegionProcessor::should_generate_heatmap(const RegionResult& result) {
    // 條件 1: Read 數量足夠
    if (result.num_reads < 20) return false;
    
    // 條件 2: HP 標籤多樣性
    // (需從 read_list 統計)
    std::set<std::string> unique_hp;
    for (const auto& r : read_list) {
        unique_hp.insert(r.hp_tag);
    }
    if (unique_hp.size() < 2) return false;
    
    // 條件 3: 有效距離對比例
    if (result.num_valid_pairs < result.num_reads * 0.5) return false;
    
    return true;
}

// 輸出標記檔案
if (should_generate_heatmap(result)) {
    std::ofstream flag(region_dir + "/.should_plot");
    flag.close();
}
```

**Python 側檢查**：

```python
def should_plot(locus_dir):
    return os.path.exists(os.path.join(locus_dir, '.should_plot'))

# 篩選
locus_dirs = [d for d in locus_dirs if should_plot(d)]
```

#### 2. 分級繪圖策略

```python
# 快速預覽（低解析度）
plot_clustermap(..., dpi=100, figsize=(8, 6))

# 顯著位點（高解析度）
if is_significant(locus):
    plot_clustermap(..., dpi=300, figsize=(14, 10))
```

---

## 4. 時間與資源估算

### 方案對比

| 指標 | 方案 A（立即呼叫）| 方案 B（批次處理）⭐ | 方案 C（C++ 繪圖）|
|------|-----------------|---------------------|------------------|
| **C++ 執行時間** | 10 min + Python等待 5.5 hr | 10 min | 30 min（含繪圖）|
| **Python 繪圖時間** | - | 60 min（16核平行）| - |
| **總時間** | **5.7 hr** | **70 min** | 30 min |
| **峰值記憶體** | 32+ GB（64個Python）| 12 GB（16個Python）| 8 GB |
| **開發成本** | 低（1小時）| 中（4小時）| **極高（2週+）** |
| **維護成本** | 低 | 低 | **極高** |
| **美觀度** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ |

### 方案 B 細節估算（40,000 位點，篩選後 10,000 個）

1. **C++ 處理**：10 分鐘（已知）
2. **Python 繪圖**（10,000 位點）：
   - 單張圖：500ms
   - 16 個 worker：10,000 / 16 × 0.5s = **312 秒**（5.2 分鐘）
3. **統計分析**：
   - PERMANOVA（40,000 位點）：~10 分鐘
4. **總計**：10 + 5 + 10 = **25 分鐘**

---

## 5. 最終建議

### ✅ 採用方案 B（兩階段批次處理）

**實作步驟**：

1. **Phase 1**（本週）：
   - C++ 新增聚類建樹功能（已有 `HierarchicalClustering`，僅需整合）
   - 輸出 Newick 格式（已有 `TreeWriter`）

2. **Phase 2**（下週）：
   - 開發 Python 批次繪圖腳本
   - 實作篩選機制與分級繪圖

3. **Phase 3**（測試）：
   - 小規模測試（100 位點）
   - 大規模驗證（40,000 位點）

### 🎯 關鍵優勢

1. **效能最佳**：總時間 < 30 分鐘（含篩選）
2. **資源可控**：記憶體 < 16 GB
3. **靈活性高**：可隨時調整繪圖參數而無需重新執行 C++
4. **易於除錯**：C++ 與 Python 失敗互不影響

### ⚠️ 注意事項

1. **自動化腳本**：建議寫一個 Shell 腳本自動化兩階段流程
2. **錯誤處理**：Python 繪圖失敗不應中斷整體流程
3. **進度監控**：使用 `tqdm` 顯示進度條

---

**結論**：方案 B 在效能、資源使用、開發成本上達到最佳平衡，強烈建議採用。
