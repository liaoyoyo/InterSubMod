---
id: ism-kb-04-parameters-index
name: "Parameters 目錄索引"
description: "ISM C++ CLI 參數、距離度量、統計方法、過濾規則、內部 Config 常數的完整索引。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "parameters directory structure"
related_ids:
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-04-parameters-distance-metrics
  - ism-kb-04-parameters-statistical-methods
  - ism-kb-04-parameters-filtering-rules
  - ism-kb-04-parameters-config-defaults
tags: [parameters, index, cli, config]
canonical_paths: [04_parameters/00_index.md]
alias_paths: []
---

# Parameters 目錄索引

- 一句結論：查 CLI 參數用 01，距離度量 02，統計檢驗 03，過濾閾值 04，內部常數 05
- 適用對象：使用 ISM binary 時、調整參數時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  /big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod --help 2>&1 | head -60
  ```

---

## 文件列表

| 檔案 | 主題 | 對應原始碼 |
|------|------|----------|
| [01_cli_arguments.md](01_cli_arguments.md) | 全 CLI 參數表（input/output、視窗、過濾、距離、日誌） | `include/utils/ArgParser.hpp:29-195` |
| [02_distance_metrics.md](02_distance_metrics.md) | NHD/L1/L2/CORR/JACCARD/BERNOULLI 數學定義 | `src/core/DistanceMatrix.cpp:1-262` |
| [03_statistical_methods.md](03_statistical_methods.md) | FFH/PERMANOVA/PERMDISP2/Cramér's V | `include/core/{GlobalTest,StructureTest,LabelTest}.hpp` |
| [04_filtering_rules.md](04_filtering_rules.md) | mapq/read_length/base_quality/methyl 閾值 | `include/utils/ArgParser.hpp:60-74` |
| [05_config_defaults.md](05_config_defaults.md) | Config 內部常數（不可由 CLI 設定） | `include/core/Config.hpp:17-109` |

---

## 快速速查：最常問的參數

| 問題 | 答案 | 詳見 |
|------|------|------|
| 預設視窗大小 | `-w 1000`（Canonical 用 5000） | 01_cli_arguments.md |
| 預設距離度量 | `--distance-metric NHD`（Canonical 用 BERNOULLI） | 02_distance_metrics.md |
| 預設 threads | `-j 1`（批次腳本用 120） | 01_cli_arguments.md |
| 預設 min-mapq | 20 | 04_filtering_rules.md |
| methylation 閾值 | 0.8 / 0.2 | 04_filtering_rules.md |
| PERMANOVA permutations | 999 | 03_statistical_methods.md |
| Self-phasing 降級 | `--germline-hp-only`（預設 off） | 01_cli_arguments.md |

---

## 相關

- Pipelines：[../03_pipelines/00_index.md](../03_pipelines/00_index.md)
- 執行流程：[../06_workflows/00_index.md](../06_workflows/00_index.md)
