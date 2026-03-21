# 實作細節與技術規格 (Implementation Details)

本文件詳細描述 `InterSubMod` 系統的技術選型、建置系統、演算法實作細節與效能優化策略。

---

## 1. 技術選型與依賴庫 (Libraries)

| Library | 用途 | 選擇理由 |
| :--- | :--- | :--- |
| **htslib** | BAM/CRAM/VCF I/O | 生物資訊領域標準 C API，效能最佳，支援 CRAM 壓縮與索引隨機存取。 |
| **ArgParser** | Command Line Parser | 自製封裝 `getopt_long` (C Standard)，輕量且無額外依賴。 |
| **Eigen3** | 矩陣運算與儲存 | 高效能 C++ 模板庫，支援 SIMD 向量化，方便處理 `MatrixXd` 等結構。 |
| **GoogleTest** | 單元測試 | 業界標準 C++ 測試框架，透過 CMake `FetchContent` 自動整合。 |
| **OpenMP** | 平行運算 | 編譯器內建，適合 Loop-level 平行化 (Per-Region)。 |
| **nlohmann/json** | JSON 序列化 | (Planned) 現代 C++ 標配，語法直觀，用於輸出 Per-read 詳細結果。 |
| **spdlog** | Logging | (Planned) 極快且 thread-safe 的 C++ logging 庫。 |

---

## 2. 建置系統 (CMake)

專案採用 CMake 進行管理，已配置自動下載 GoogleTest 並透過 PkgConfig 尋找系統庫。

### CMakeLists.txt 核心架構

```cmake
cmake_minimum_required(VERSION 3.14)
project(InterSubMod CXX)

# ... Compiler Options ...

# --- Dependencies ---
find_package(OpenMP REQUIRED)

# HTSlib via PkgConfig
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSlib REQUIRED htslib)

# Eigen3
find_package(Eigen3 3.3 QUIET)
# ... fallback logic ...

# --- Tests (GoogleTest) ---
include(FetchContent)
FetchContent_Declare(googletest ...)
FetchContent_MakeAvailable(googletest)
```

---

## 3. 核心演算法實作細節

### 3.1 BAM CIGAR 與 SNV 座標對應

在 `S2` 步驟中，必須精確判斷 Read 上的哪個 base 對應到 Genome 上的 SNV 位置。

* **挑戰**：INDELs 與 Soft-clipping 會導致 Read 座標與 Reference 座標非線性對應。
* **實作邏輯**：
    1. 初始化 `ref_pos = align_start`, `read_pos = 0`。
    2. 遍歷 CIGAR ops：
        * `M/=/X` (Match/Mismatch): `ref_pos` 與 `read_pos` 同步增加。若 SNV 在此區間，則 Read Base = `seq[read_pos + (snv_pos - ref_pos)]`。
        * `I` (Insertion): `read_pos` 增加，`ref_pos` 不變。
        * `D` (Deletion): `ref_pos` 增加，`read_pos` 不變。若 SNV 在此區間，標記為 **REF** (若視為無變異) 或 **UNKNOWN** (視 policy 而定)。
        * `S` (Soft Clip): `read_pos` 增加，`ref_pos` 不變。
    3. **Base Quality 檢查**：取得 Base 後，必須檢查對應的 Quality Score。若 `< min_base_quality`，則標記為 `AltSupport::UNKNOWN`。

### 3.2 MM/ML 甲基化標籤解析

在 `S3` 步驟中，需將 `MM` (Modified Base) tags 轉換為 Genomic 座標。

* **MM 格式**：`C+m?,5,0,3;` (Delta encoding)。
* **解析流程**：
    1. **確認方向**：若 `bam_is_rev(b)`，需確認 MM tag 是否已由 Caller 轉為 5'->3' (SAM Spec v1.6+ 建議一致性)。通常 Dorado/Guppy 輸出的 MM 對應 Read Sequence (儲存方向)。
    2. **定位 C 位點**：遍歷 Read Sequence 尋找 `C` (或 `G` if mapped to rev and logic requires)。
    3. **跳過計數**：依據 MM 的數字跳過非甲基化的 C，鎖定有紀錄的 C。
    4. **機率查找**：使用計數器 index 去 `ML` array 查機率值。
    5. **座標轉換**：利用 3.1 的 CIGAR 邏輯，將 Read 上的 C index 映射回 Genome 座標 `(chr, pos)`。

### 3.3 相位區塊 (Phase Set) 處理

* **問題**：Region Window (±1000bp) 可能跨越 Phase Block 邊界，導致 HP=1 在左側與右側代表不同單倍體。
* **實作檢查**：
  * 讀取每條 Read 的 `PS` tag。
  * 統計 Window 內的主流 `PS`。
  * 若發現多個 `PS` 且佔比顯著，標記該 Region 為 `PHASE_AMBIGUOUS`，或僅保留主流 `PS` 的 Reads。

---

## 4. 效能與記憶體優化

### 4.1 記憶體佈局 (Memory Layout)

* **矩陣儲存**：使用 **Row-Major** (`Eigen::Matrix<..., RowMajor>`)。
  * **原因**：距離計算是比較兩個 Reads (Rows)。Row-Major 讓同一條 Read 的 CpG 數據在記憶體中連續，大幅提升 CPU Cache Hit rate。
* **Binary Matrix**：不使用 `int` 存 0/1，改用 `std::vector<uint64_t>` 或 `boost::dynamic_bitset` 實作 Bitset。
  * **NHD 計算**：使用 `XOR` + `popcount` 指令，可一次比較 64 個 CpG sites，速度提升數十倍。

### 4.2 平行化策略 (Parallelism)

* **Region-Level (粗粒度)**：
  * 使用 OpenMP `#pragma omp parallel for` 在最外層的 Region 迴圈。
  * 這是最無痛且效率最高的平行化方式，因為每個 Region 計算完全獨立。
  * **記憶體注意**：每個 Thread 會擁有自己的一份 Matrix Buffer，需確保總記憶體不超過機器限制 (`jemalloc` 在此有助於減少 fragmentation)。

### 4.3 I/O 優化

* **HTSlib Thread Pool**：
  * 啟用 `hts_set_threads(infile, n_threads)`。
  * 讓 BGZF 解壓縮在背景執行，主執行緒專注於 Parsing 與 Data Structuring。
