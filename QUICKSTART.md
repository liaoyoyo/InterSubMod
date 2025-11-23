# InterSubMod - Quick Start Guide

## 快速開始

### 1. 編譯專案
```bash
cd /big8_disk/liaoyoyo2001/InterSubMod
mkdir -p build && cd build
cmake ..
make -j4
```

### 2. 準備輸入檔案

**必需檔案**:
- `tumor.bam` (+ `.bai` index) - 包含 MM/ML tags 的 BAM
- `hg38.fa` (+ `.fai` index) - 參考基因組
- `snvs.tsv` - SNV 列表 (格式: chr pos ref alt qual)

**SNV Table 格式範例**:
```tsv
chr	pos	ref	alt	qual
chr17	7578000	C	T	100.0
chr17	7580000	G	A	95.5
```

### 3. 執行分析

#### 方法 A: 使用測試程式 (推薦入門)
```bash
cd /big8_disk/liaoyoyo2001/InterSubMod

# 測試 Phase 1 & 2: BAM 讀取與甲基化解析
./build/bin/test_phase1_2

# 測試 Phase 3: 矩陣建構與輸出
./build/bin/test_phase3

# 測試 Phase 4 & 5: 平行化處理 32 SNVs
./build/bin/test_phase4_5
```

#### 方法 B: 使用 C++ API
```cpp
#include "core/RegionProcessor.hpp"

int main() {
    // 初始化
    RegionProcessor processor(
        "tumor.bam",
        "normal.bam",
        "hg38.fa",
        "output/",
        4,      // threads
        2000    // window size (±2000 bp)
    );
    
    // 載入 SNVs
    processor.load_snvs("snvs.tsv");
    
    // 平行處理所有 regions
    auto results = processor.process_all_regions();
    
    // 輸出摘要
    processor.print_summary(results);
    
    return 0;
}
```

### 4. 輸出結構

```
output/
├── region_0000/
│   ├── metadata.txt      # Region 資訊 & 統計
│   ├── reads.tsv         # Read 列表
│   ├── cpg_sites.tsv     # CpG 位點座標
│   └── methylation.csv   # Read × CpG 矩陣
├── region_0001/
└── ...
```

### 5. 驗證輸出

```bash
./scripts/verify_output.sh
```

---

## 輸出檔案說明

### metadata.txt
```
Region ID: 0
Region: chr17:7576000-7580000
Region Size: 4001 bp

SNV ID: 0
SNV Position: chr17:7578000
SNV: C -> T
SNV Quality: 100

Num Reads: 154
Num CpG Sites: 59
Matrix Dimensions: 154 × 59

Processing Time: 237.98 ms
Peak Memory: 0.00 MB
```

### reads.tsv
```
read_id	read_name	chr_id	start	end	mapq	hp	alt_support	is_tumor
0	read_1	17	7576100	7579800	60	1	UNKNOWN	1
1	read_2	17	7576200	7579500	60	2	UNKNOWN	1
...
```

### cpg_sites.tsv
```
cpg_id	chr_id	position
0	17	7576001
1	17	7576123
2	17	7576456
...
```

### methylation.csv
```
read_id,7576001,7576123,7576456,...
0,NA,0.9234,0.8765,...
1,0.6543,NA,0.9012,...
...
```
- **NA**: 該 read 未覆蓋此 CpG 位點
- **數值**: 甲基化機率 [0.0, 1.0]

---

## 常見問題

### Q: 如何調整平行化執行緒數？
A: 修改 `RegionProcessor` 建構函式的 `num_threads` 參數
```cpp
RegionProcessor processor(..., 8, ...);  // 使用 8 threads
```

### Q: 如何修改 Region 窗口大小？
A: 修改 `window_size` 參數
```cpp
RegionProcessor processor(..., 4, 5000);  // ±5000 bp
```

### Q: BAM 沒有 MM/ML tags 怎麼辦？
A: 使用 ONT basecaller (如 Guppy/Dorado) 重新 basecall 並加上 `--modified-bases 5mC` 選項

### Q: 如何只處理特定數量的 SNVs？
A: 使用 `process_all_regions()` 的參數
```cpp
auto results = processor.process_all_regions(10);  // 只處理前 10 個
```

---

## 效能優化建議

1. **增加執行緒數** (適用於多核心 CPU)
   ```cpp
   RegionProcessor processor(..., 16, ...);  // 使用更多 threads
   ```

2. **調整 BAM reader threads**
   ```cpp
   bam_reader.set_hts_threads(4);  // HTSlib 解壓縮 threads
   ```

3. **批次處理** (一次處理多個 SNVs)
   ```cpp
   // 將大量 SNVs 分批處理
   for (int batch = 0; batch < num_batches; batch++) {
       auto results = processor.process_batch(batch);
   }
   ```

---

## 支援與文件

- **完整文件**: `docs/`
  - 設計文件: `docs/design/`
  - 開發計畫: `docs/dev/`
  - 除錯紀錄: `docs/debug/`
  
- **測試程式**: `src/test_*.cpp`
  
- **驗證腳本**: `scripts/verify_output.sh`

---

## 系統需求

- **CPU**: 多核心 (建議 ≥4 cores)
- **RAM**: ≥8 GB
- **Storage**: 視 BAM 大小而定
- **OS**: Linux (測試於 Ubuntu)
- **Compiler**: GCC ≥7.0 或 Clang ≥5.0 (支援 C++17 & OpenMP)

