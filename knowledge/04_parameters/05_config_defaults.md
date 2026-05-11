---
id: ism-kb-04-parameters-config-defaults
name: "Config Internal Defaults"
description: "無法由 CLI 調整的內部 Config 常數：UPGMA linkage、Cramér's V threshold、最小分群 reads、kNumPermutations 等。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "internal constants against include/core/Config.hpp:17-109"
related_ids:
  - ism-kb-04-parameters-index
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-04-parameters-statistical-methods
tags: [parameters, config, defaults, internal, constants]
canonical_paths: [04_parameters/05_config_defaults.md]
alias_paths: []
---

# Config Internal Defaults

- 一句結論：10+ 個內部常數（linkage、閾值、permutations）寫死於 Config.hpp，**無法由 CLI 改**；需改須修 source
- 適用對象：修改 Config 前、理解 ISM 內部行為
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -nE 'constexpr|DEFAULT_|k[A-Z]' /big7_disk/liaoyoyo2001/InterSubMod/include/core/Config.hpp | head -30
  ```

---

## 內部常數清單

| 常數 | 值 | 用途 | 位置 |
|------|-----|------|------|
| 預設距離度量（fallback） | NHD | ArgParser 若解析失敗的 fallback | ArgParser.hpp:173-176 |
| Linkage method | UPGMA | 階層式分群 | Config.hpp:66 |
| Min reads for clustering | 10 | 少於此數跳過分群 | Config.hpp:67 |
| Min reads per group (fine-grained) | 3 | Fine-grained PERMANOVA 每群最小 reads | RegionProcessor.cpp |
| Cramér's V threshold | 0.1 | 判定有意義效應 | **GlobalTest.hpp:32** (`min_cramers_v`；**非** Config.hpp) |
| PERMANOVA permutations | 999 struct default / **99 runtime** | Permutation 次數（RegionProcessor.cpp:1573 override 為 99） | StructureTest.hpp:26 / RegionProcessor.cpp:1573 |
| Heuristic score adjustments | ×0.5, ×0.7 | Dispersion warn / PERMANOVA ns 調整 | RegionProcessor.cpp |
| MM/ML parsing strand handling | 依 BAM spec | MM/ML tag 解析方式 | MethylationParser |
| Debug filtered reads auto-enable | ≥ LOG_DEBUG | 自動啟用 --output-filtered-reads | ArgParser.hpp:190-192 |

---

## 為何這些不放 CLI？

| 常數 | 為何不開放 |
|------|-----------|
| UPGMA linkage | 換 linkage（single/complete/ward）結果差異大，需充分研究後才能開放 |
| Cramér's V threshold 0.1 | 統計共識值；改了難以跨樣本比較 |
| 999 permutations | 標準；改小會影響 p-value 精度 |
| Min reads 10 / 3 | 統計有效性底線 |

---

## 修改流程（若真需改）

1. 讀 `include/core/Config.hpp` 找對應常數
2. 修改 default value
3. 重新編譯（`cd build && make -j`）
4. **關鍵**：加 unit test 確保行為正確
5. 在 `09_conclusions/` 或 `docs/solutions/` 留文件說明改動理由
6. 重跑 canonical benchmark 對比前後 F1

**警告**：改內部常數可能導致 canonical run 結果無法對比；慎重！

---

## Config 結構（概覽）

```cpp
// include/core/Config.hpp（節錄）
struct Config {
    // CLI-settable
    std::string tumor_bam_path;
    std::string normal_bam_path;
    std::string reference_fasta_path;
    std::string somatic_vcf_path;
    std::string output_dir = "output";
    int window_size_bp = 1000;
    int threads = 1;
    double binary_methyl_high = 0.8;
    double binary_methyl_low = 0.2;
    int min_mapq = 20;
    int min_read_length = 1000;
    int min_base_quality = 20;
    std::vector<DistanceMetricType> distance_metrics;
    int min_common_coverage = 3;
    NanDistanceStrategy nan_distance_strategy;
    double max_distance_value = 1.0;
    double expected_coverage = 0.0;
    bool germline_hp_only = false;
    bool use_full_read_span = false;
    LogLevel log_level = LogLevel::LOG_INFO;
    // ...

    // Non-CLI (internal)
    std::string linkage_method = "UPGMA";    // :66
    int min_reads_for_clustering = 10;        // :67
    // ...
};
```

**完整 Config**：直接 `cat include/core/Config.hpp`

---

## 相關

- CLI 參數：[01_cli_arguments.md](01_cli_arguments.md)
- 統計方法：[03_statistical_methods.md](03_statistical_methods.md)
- 原始碼：`include/core/Config.hpp:17-109`
