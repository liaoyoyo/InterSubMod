# ISM Config / Threshold Extract — 二值化閾值 + 讀取過濾

> aspect: `config` | 稽核日期 2026-06-13 | 來源：實際讀取 repo 源碼，逐條附 檔:行
> 範圍：`include/core/Config.hpp` + `src/core/Config.cpp` + CLI (`include/utils/ArgParser.hpp`) + 下游使用點 (`RegionProcessor.cpp` / `DistanceMatrix.cpp` / `ReadParser.cpp`)

---

## L0 一眼結論

- 二值化閾值預設：**high=0.8 / low=0.2**（中間 0.2~0.8 = ambiguous 標 -1）。`include/core/Config.hpp:33-34`
- 二值化只在 **building methylation matrix 的 `binary_matrix`** 用到（`RegionProcessor.cpp:1417-1422`）。
- **距離計算「依 metric 而定」用連續或二值**：NHD / JACCARD 用 `binary_matrix`（二值）；L1 / L2 / CORR / BERNOULLI 用 `raw_matrix`（連續）。`DistanceMatrix.cpp:309-345`
- 讀取過濾預設：**min_mapq=20 / min_read_length=1000 / min_base_quality=20**（後者只作用於 SNV 位點的鹼基品質，**不**過濾 CpG/甲基讀取）。`Config.hpp:29-31` + `ReadParser.cpp:45,51,245`

---

## L1 二值化閾值（methyl-high / methyl-low）

### 預設值與語意
- `double binary_methyl_high = 0.8;` // Threshold for methylated (1) call —— `include/core/Config.hpp:33`
- `double binary_methyl_low = 0.2;` // Threshold for unmethylated (0) call —— `include/core/Config.hpp:34`

### 二值化規則（哪裡用到？）
唯一二值化發生在 `RegionProcessor::build_methylation_matrix`，把 `raw_matrix[i][j]`（連續甲基化值）映成 `binary_matrix`：

```cpp
// src/core/RegionProcessor.cpp:1410-1423
double val = raw_matrix[i][j];
if (val < 0) {  // -1.0 indicates no coverage
    meth_mat.raw_matrix(i, j) = NAN;
    meth_mat.binary_matrix(i, j) = -1;
} else {
    meth_mat.raw_matrix(i, j) = val;
    // Binary threshold (configurable via --methyl-high/--methyl-low)
    if (val >= binary_methyl_high_) {
        meth_mat.binary_matrix(i, j) = 1;
    } else if (val <= binary_methyl_low_) {
        meth_mat.binary_matrix(i, j) = 0;
    } else {
        meth_mat.binary_matrix(i, j) = -1;  // Ambiguous
    }
}
```

- `val >= high(0.8)` → 1（methylated）
- `val <= low(0.2)` → 0（unmethylated）
- `low < val < high`（0.2~0.8）→ **-1（Ambiguous，等同 no-coverage 被排除）** `RegionProcessor.cpp:1422`
- `val < 0`（-1.0 = no coverage）→ raw=NaN, binary=-1 `RegionProcessor.cpp:1411-1413`
- 注意：**邊界值含端點**（`>=` high 與 `<= low`），val=0.8 算 1、val=0.2 算 0。

成員初始化（閾值從 Config 帶入 RegionProcessor）：
```cpp
// src/core/RegionProcessor.cpp:390-391
binary_methyl_high_(config.binary_methyl_high),
binary_methyl_low_(config.binary_methyl_low) {
```
宣告：`include/core/RegionProcessor.hpp:509-510`

### 驗證/約束
- `validate()`：`binary_methyl_high <= binary_methyl_low` → 報錯（high 必須嚴格大於 low）。`Config.cpp:108-111`
- `validate()`：兩閾值必須落在 [0.0, 1.0]。`Config.cpp:113-116`
- `print()` 輸出 `Methylation Thresholds: Low=..., High=...`。`Config.cpp:131`

---

## L1 距離計算：連續 vs 二值（依 metric）

`Config` 有 `bool distance_use_binary = true;`（`Config.hpp:50`）對應 `DistanceConfig::use_binary_matrix = true`（`DistanceMatrix.hpp:26`），但**實際 dispatcher 是按 metric 硬選矩陣，不是讀此旗標**：

```
// src/core/DistanceMatrix.cpp:309-345（calculate_distance_impl dispatcher）
NHD       → mat.binary_matrix.row(...)   // 二值        (309-313)
L1        → mat.raw_matrix.row(...)      // 連續        (316-319)
L2        → mat.raw_matrix.row(...)      // 連續        (322-325)
CORR      → mat.raw_matrix.row(...)      // 連續        (328-331)
JACCARD   → mat.binary_matrix.row(...)   // 二值        (334-338)
BERNOULLI → mat.raw_matrix.row(...)      // 連續        (341-344)
```

→ **結論**：「距離計算用連續還是二值」**取決於選用的 metric**。NHD/JACCARD 吃二值矩陣（因此受 methyl-high/low 閾值影響，含 ambiguous→-1 的排除）；L1/L2/CORR/BERNOULLI 吃連續 raw 值（不受二值化閾值影響）。
（**uncertain**：`distance_use_binary` / `use_binary_matrix` 旗標在 dispatcher 中未見被讀取，疑為未生效或舊路徑遺留；見 uncertain 區。）

### 預設 metric 不一致（潛在 bug/陷阱）
- `Config.hpp:40` 預設 `distance_metrics = {DistanceMetricType::BERNOULLI}`（→ 連續 raw）
- `ArgParser.hpp:86` CLI 預設字串 `distance_metric_strs = {"NHD"}`（→ 二值）

兩處預設不同（**uncertain**：實際生效者取決於 ArgParser 是否在未給 `--distance-metric` 時用 "NHD" 覆寫 Config 的 BERNOULLI，需追 ArgParser 後段轉換邏輯確認）。

---

## L1 讀取過濾閾值（min-mapq / min-base-quality / min-read-length）

### 預設值
- `int min_mapq = 20;` —— `Config.hpp:29`
- `int min_read_length = 1000;` // bp —— `Config.hpp:30`
- `int min_base_quality = 20;` // for SNV/CpG sites（comment 寫 SNV/CpG，但實作只用於 SNV，見下）—— `Config.hpp:31`
- （無 `--min-site-coverage` 對應；`min_site_coverage = 5` 僅宣告見 uncertain 區）

### CLI 範圍檢查（ArgParser）
- `--methyl-high` `--methyl-low`：`CLI::Range(0.0, 1.0)` —— `ArgParser.hpp:60-64`
- `--min-mapq`：`CLI::Range(0, 60)`，help 寫 Default: 20 —— `ArgParser.hpp:67-68`
- `--min-read-length`：`CLI::PositiveNumber`，help 寫 Default: 1000 —— `ArgParser.hpp:70-71`
- `--min-base-quality`：`CLI::Range(0, 93)`，help 寫 "at SNV"，Default: 20 —— `ArgParser.hpp:73-74`

### 如何影響甲基讀取 / 矩陣（實作點）
從 Config 帶到 ReadParser 的 filter_config：
```cpp
// src/core/RegionProcessor.cpp:411-413
filter_config_.min_mapq = config.min_mapq;
filter_config_.min_read_length = config.min_read_length;
filter_config_.min_base_quality = config.min_base_quality;
```
ReadParser 預設同步：`include/core/ReadParser.hpp:16-18`（min_mapq=20 / min_read_length=1000 / min_base_quality=20）。

實際過濾：
- **MAPQ**：`if (b->core.qual < config_.min_mapq) reasons |= LOW_MAPQ;` —— `ReadParser.cpp:45`（整條 read 層級過濾 → 被濾掉的 read 不進甲基矩陣任何列）
- **read length**：`if (read_len < config_.min_read_length) reasons |= SHORT_READ;`（`read_len = bam_cigar2qlen(...)`）—— `ReadParser.cpp:50-52`（read 層級過濾）
- **base quality**：`if (qual[read_offset] < config_.min_base_quality) return ...LOW_BASE_QUALITY;` —— `ReadParser.cpp:245`。**此處 `read_offset` 是 SNV 位點 offset**（緊接 `seq[read_offset]` 比對 REF/ALT，`ReadParser.cpp:249-259`）→ 僅 gating「該 read 對 SNV 的 ALT/REF 支持判定」，**不直接過濾 CpG 甲基讀值**。

→ **影響鏈總結**：
1. min_mapq / min_read_length = **read 級**過濾，決定哪些 read 進得了甲基化矩陣（影響矩陣列數 / 覆蓋深度 → 間接影響每個 CpG site 的可用覆蓋與距離可計算性）。
2. min_base_quality = **SNV 位點級**過濾，只影響 ALT/REF 支持分類（haplotype/somatic 判定），**不影響 CpG 甲基化值二值化**（comment 「for SNV/CpG sites」與實作不符，見 uncertain）。

---

## L1 其他相關 coverage 閾值
- `int min_site_coverage = 5;` // Minimum reads covering a CpG site to keep it —— `Config.hpp:36`（**uncertain**：grep `src include` 無任何使用點，疑未接線/dead config）
- `int min_common_coverage = 3;` // C_min，計算距離所需最少共同 CpG —— `Config.hpp:37`；CLI `--min-common-coverage` Default: 3，`ArgParser.hpp:92-94`；dispatcher 各 metric 傳入 `config.min_common_coverage`（`DistanceMatrix.cpp:313,319,...`）→ **此為實際生效的「最少共同覆蓋」門檻**。
