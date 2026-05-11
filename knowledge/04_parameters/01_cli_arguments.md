---
id: ism-kb-04-parameters-cli-arguments
name: "CLI Arguments 完整表"
description: "ISM binary (`build/bin/inter_sub_mod`) 所有命令列參數：輸入/輸出、視窗、過濾、距離度量、日誌、self-phasing 降級；含預設值、取值範圍、原始碼 file:line。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "CLI args against include/utils/ArgParser.hpp:29-195 on HEAD"
related_ids:
  - ism-kb-04-parameters-index
  - ism-kb-04-parameters-distance-metrics
  - ism-kb-04-parameters-filtering-rules
  - ism-kb-04-parameters-config-defaults
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-paired-pileup
  - ism-kb-06-workflows-debug-and-logging
tags: [parameters, cli, arguments, ArgParser]
canonical_paths: [04_parameters/01_cli_arguments.md]
alias_paths: []
---

# CLI Arguments 完整表

- 一句結論：28+ 個 CLI 參數分為輸入/輸出、視窗、過濾、距離、統計、日誌、self-phasing 7 類；定義於 `include/utils/ArgParser.hpp:29-195`
- 適用對象：執行 `./build/bin/inter_sub_mod` 者、寫 wrapper script 者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  /big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod --help
  ```

---

## 1. 輸入／輸出路徑

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--tumor-bam` | `-t` | string | **必填** | 有效 BAM | Tumor BAM | ArgParser.hpp:33-35 |
| `--normal-bam` | `-n` | string | (空) | 有效 BAM 或空 | Normal BAM（選填） | ArgParser.hpp:37-38 |
| `--reference` | `-r` | string | **必填** | 有效 FASTA | 參考基因組 | ArgParser.hpp:43-45 |
| `--vcf` | `-v` | string | **必填** | 有效 VCF.gz | 體細胞變異清單 | ArgParser.hpp:47-49 |
| `--output-dir` | `-o` | string | `"output"` | 任意目錄 | 輸出根目錄 | ArgParser.hpp:51 |
| `--loh-bed` | — | string | (空) | 有效 BED | LongPhase LOH 註記（選填） | ArgParser.hpp:40-41 |
| `--debug-output-dir` | — | string | `{output}/debug` | 任意目錄 | Debug 輸出位置 | ArgParser.hpp:109-110 |

## 2. 分析視窗

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--window-size` | `-w` | int | 1000 | > 0 | SNV 前後 ±bp 視窗（Canonical 用 5000） | ArgParser.hpp:54-55 |
| `--full-read` | — | flag | false | — | 啟用動態視窗涵蓋全讀取跨度 | ArgParser.hpp:119 |

## 3. 讀取過濾

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--min-mapq` | — | int | 20 | [0, 60] | 最小映射品質 | ArgParser.hpp:67-68 |
| `--min-read-length` | — | int | 1000 | > 0 | 最小讀取長度 (bp) | ArgParser.hpp:70-71 |
| `--min-base-quality` | — | int | 20 | [0, 93] | SNV 位置鹼基品質 | ArgParser.hpp:73-74 |
| `--no-filter` | — | flag | false | — | 不過濾輸出所有 reads（驗證用） | ArgParser.hpp:115-116 |

## 4. 甲基化閾值

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--methyl-high` | — | double | 0.8 | [0.0, 1.0] | 二元甲基化「甲基化」判定閾值 | ArgParser.hpp:60-61 |
| `--methyl-low` | — | double | 0.2 | [0.0, 1.0] | 二元甲基化「未甲基化」判定閾值 | ArgParser.hpp:63-64 |

## 5. 距離度量

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--distance-metric` | — | list[string] | `["NHD"]` | NHD / L1 / L2 / CORR / JACCARD / BERNOULLI（大小寫不敏） | 計算距離度量列表 | ArgParser.hpp:87-90 |
| `--min-common-coverage` | — | int | 3 | > 0 | 計算距離所需最小 CpG 重疊（C_min） | ArgParser.hpp:92-94 |
| `--nan-distance-strategy` | — | string | `"MAX_DIST"` | MAX_DIST / SKIP | 重疊不足的距離對策 | ArgParser.hpp:97-99 |
| `--max-distance-value` | — | double | 1.0 | [0.0, 1000.0] | MAX_DIST 策略的距離值 | ArgParser.hpp:101-102 |

## 6. 距離矩陣輸出

| 參數 | 類型 | 預設 | 用途 | 原始碼 |
|------|------|------|------|--------|
| `--compute-distance-matrix` / `--no-distance-matrix` | flag | true | 計算距離矩陣 | ArgParser.hpp:77-78 |
| `--output-distance-matrix` / `--no-output-distance-matrix` | flag | true | 輸出 CSV 距離矩陣 | ArgParser.hpp:80-81 |
| `--output-strand-distance-matrices` | flag | true | 輸出正反股分開的距離矩陣 | ArgParser.hpp:83-84 |

## 7. Coverage 正規化

| 參數 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|--------|
| `--expected-coverage` | double | 0.0 | [0.0, 1000.0] | 期望雙倍體覆蓋度；0 = 自動估計（KDE mode） | ArgParser.hpp:122-124 |

## 8. Self-phasing 降級

| 參數 | 類型 | 預設 | 用途 | 原始碼 |
|------|------|------|------|--------|
| `--germline-hp-only` | flag | false | 忽略 somatic HP tags (HP:i:11/21/33)，僅用 germline HP:i:1/2 | ArgParser.hpp:127-128 |

**注意**：flag=on 會影響 HPFineNGroups 結論；見 [03_pipelines/03_tumor_only.md](../03_pipelines/03_tumor_only.md)。

## 9. 平行化與日誌

| 參數 | 別名 | 類型 | 預設 | 取值 | 用途 | 原始碼 |
|------|------|------|------|------|------|--------|
| `--threads` | `-j` | int | **16** | > 0 | OpenMP 執行緒數（Config.hpp:43 struct 預設 16；ArgParser help string 誤植 "Default: 1"；批次腳本通常覆蓋為 120） | ArgParser.hpp:57 / Config.hpp:43 |
| `--log-level` | — | string | `"info"` | fatal / error / warn / info / debug / trace | 日誌級別 | ArgParser.hpp:106-107 |
| `--output-filtered-reads` | — | flag | auto(≥debug) | — | 輸出過濾讀取詳細日誌 | ArgParser.hpp:112-113 |

---

## 完整範例命令

**Canonical paired_full 執行**（含所有主要參數）：
```bash
./build/bin/inter_sub_mod \
  -t tumor_haplotagged.bam \
  -n normal.bam \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -v filtered_snv_tp.vcf.gz \
  --loh-bed phased.bed \
  -w 5000 \
  -j 120 \
  --distance-metric BERNOULLI \
  --min-common-coverage 3 \
  --min-mapq 20 \
  --methyl-high 0.8 \
  --methyl-low 0.2 \
  --log-level info \
  -o output/paired_full_run
```

**Debug + TO + self-phasing 降級**：
```bash
./build/bin/inter_sub_mod \
  -t tumor_tagged.bam \
  -r ref.fa \
  -v clairsto_tp.vcf.gz \
  --germline-hp-only \
  --log-level debug \
  --output-filtered-reads \
  -o output/to_debug_run
```

---

## 與 Config defaults 的差異

CLI 只能設定**部分**參數。其他內部參數（如 UPGMA linkage、Cramér's V threshold=0.1）在 Config 寫死，見 [05_config_defaults.md](05_config_defaults.md)。

---

## 相關

- 距離度量數學：[02_distance_metrics.md](02_distance_metrics.md)
- 過濾規則：[04_filtering_rules.md](04_filtering_rules.md)
- Config 常數：[05_config_defaults.md](05_config_defaults.md)
- Paired Full 使用：[../03_pipelines/01_paired_full.md](../03_pipelines/01_paired_full.md)
- 原始碼：`include/utils/ArgParser.hpp`
