# InterSubMod - 快速開始 (Quick Start)

本文檔將協助您快速配置環境並執行 InterSubMod 分析流程。

---

## 1. 環境準備與編譯 (Setup & Build)

### 系統需求
*   **OS**: Linux (Ubuntu 18.04+ Recommended)
*   **Compiler**: GCC 7+ 或 Clang 5+ (支援 C++17)
*   **Dependencies**: HTSlib 1.10+, CMake 3.14+, Python 3 (Matplotlib, Seaborn, Pandas)

### 編譯步驟

```bash
# 1. Clone 專案
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod

# 2. 建立並編譯 (使用最大並行數加速)
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

編譯完成後，主程式位於 `build/bin/inter_sub_mod`。

---

## 2. 執行分析 (Run Analysis)

我們提供了一個整合腳本，能夠自動完成核心運算與圖表繪製。

### 預設快速執行

請使用以下指令執行完整的 VCF 測試流程：

```bash
# 回到專案根目錄
cd /big8_disk/liaoyoyo2001/InterSubMod

# 執行全流程測試 (包含 reads 過濾、矩陣建構、視覺化)
./scripts/run_full_vcf_test.sh --mode all-with-w1000
```

### 指令說明

*   `--mode all-with-w1000`: 啟用標準過濾器，並設定分析窗口為 SNV 前後 ±1000bp。這是最推薦的標準分析模式。
*   此腳本會自動執行以下步驟：
    1.  呼叫 C++ 核心程式 `inter_sub_mod` 處理數據。
    2.  產出甲基化矩陣與 reads 資訊。
    3.  呼叫 Python 腳本生成距離熱圖與分群熱圖。

### 進階選項

您也可以自定義參數：

```bash
# 指定輸出目錄
./scripts/run_full_vcf_test.sh --mode all-with-w1000 -o output/my_custom_run

# 僅生成距離熱圖 (Skip Cluster Heatmap)
./scripts/run_full_vcf_test.sh --mode all-with-w1000 --plot-type distance

# 使用更多執行緒 (C++ Core: 64, Python Plotting: 64)
./scripts/run_full_vcf_test.sh --threads 64 --plot-threads 64
```

---

## 3. 輸出結果檢視 (Output)

執行完成後，請前往輸出目錄（預設為 `output/yyyymmdd_vcf_all_w1000_t*`）。
目錄結構如下：

```text
output/
├── full_execution_analysis.log  # 完整執行日誌
├── region_0000/                 # 第一個 SNV 區域
│   ├── metadata.txt             # 統計資訊
│   ├── reads.tsv                # Reads 詳細列表
│   ├── methylation.csv          # 原始甲基化矩陣
│   ├── distance_matrix_L1.csv   # L1 距離矩陣
│   ├── distance_heatmap_L1.png  # [圖表] Read-Read 距離熱圖
│   └── cluster_heatmap_L1.png   # [圖表] 甲基化分群熱圖
├── region_0001/
└── ...
```

---

## 4. 進階：手動執行 C++ 核心 (Manual Execution)

若需直接控制底層參數，可直接呼叫執行檔：

```bash
./build/bin/inter_sub_mod \
    --tumor-bam data/tumor.bam \
    --reference data/ref.fa \
    --vcf data/somatic.vcf \
    --output-dir results \
    --threads 32 \
    --window-size 1000 \
    --log-level debug
```

> [!NOTE]
> 手動執行僅產出數據檔案 (TSV/CSV)，不會自動生成熱圖。如需圖表，請接續執行 `tools/` 下的 Python 腳本。
