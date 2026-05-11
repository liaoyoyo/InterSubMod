---
id: ism-kb-06-workflows-debug-and-logging
name: "Debug & Logging Workflow"
description: "`--log-level debug` + `--output-filtered-reads` 排查 read 過濾問題；trace 級深入；debug_output_dir 配置。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: howto
verified_scope: "debug flags against include/utils/ArgParser.hpp:105-116"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-04-parameters-filtering-rules
  - ism-kb-04-parameters-cli-arguments
tags: [workflow, debug, logging, filtered-reads, trace]
canonical_paths: [06_workflows/05_debug_and_logging.md]
alias_paths: []
---

# Debug & Logging Workflow

- 一句結論：`--log-level debug` 自動啟用 `--output-filtered-reads`；debug/filtered_reads.tsv 記錄過濾原因
- 適用對象：排查 read 數異常、過濾問題
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ./build/bin/inter_sub_mod --log-level debug --output-filtered-reads \
    -t tumor.bam -r ref.fa -v small.vcf.gz -o debug_run
  ls debug_run/debug/
  ```

---

## Log Level 六級

| 級別 | 用途 | 輸出量 |
|------|------|--------|
| `fatal` | 僅致命錯誤 | 極少 |
| `error` | 錯誤 | 少 |
| `warn` | 警告 | 少 |
| `info` | 預設；正常流程訊息 | 中 |
| `debug` | 詳細（自動啟用 filtered_reads 輸出） | 多 |
| `trace` | 極詳細（每 read 處理） | 極多 |

**設定**：`--log-level debug`

---

## Debug 輸出目錄

**參數**：`--debug-output-dir`（預設 `{output_dir}/debug/`）

**內容**：
```
<debug_output_dir>/
├── filtered_reads.tsv           # 被過濾的 reads 清單
├── ism_run_<timestamp>.log       # 完整 log
└── (其他 debug artifacts)
```

### filtered_reads.tsv 欄位
| 欄位 | 意義 |
|------|------|
| read_id | 內部編號 |
| read_name | BAM QNAME |
| filter_reason | 被過濾原因（見下表） |
| chr, pos | SNV 位置 |
| mapq | 該 read MAPQ |
| read_length | 該 read 長度 |

### filter_reason 取值
| 值 | 說明 |
|----|------|
| `LOW_MAPQ` | MAPQ < `--min-mapq` |
| `SHORT_READ` | 長度 < `--min-read-length` |
| `LOW_BASE_QUAL_AT_SNV` | SNV 位置鹼基品質 < `--min-base-quality` |
| `MISSING_MM_ML_TAG` | 缺 MM/ML tag（無法解甲基化） |
| `CIGAR_ABNORMAL` | CIGAR 異常 |
| 其他 | 專案定義的過濾原因 |

---

## 典型 Debug 流程

### 情境 1：某 region read 數異常低
```bash
# 1. 關閉過濾重跑
./build/bin/inter_sub_mod --no-filter --log-level debug \
  -t tumor.bam -r ref.fa -v region.vcf.gz -o diag_nofilter

# 2. 開啟過濾但記錄
./build/bin/inter_sub_mod --log-level debug --output-filtered-reads \
  -t tumor.bam -r ref.fa -v region.vcf.gz -o diag_withfilter

# 3. 對比 read 數
wc -l diag_nofilter/TP/region_0/reads.tsv diag_withfilter/TP/region_0/reads.tsv

# 4. 看過濾原因分布
awk -F'\t' 'NR>1{print $3}' diag_withfilter/debug/filtered_reads.tsv | sort | uniq -c
```

### 情境 2：甲基化 NaN 太多
```bash
# 1. 確認 MM/ML tag 存在
samtools view tumor.bam | head | grep -oE 'MM:Z:[^\t]*|ML:B:[^\t]*'

# 2. trace 級別看 MM/ML 解析
./build/bin/inter_sub_mod --log-level trace \
  -t tumor.bam -r ref.fa -v region.vcf.gz -o trace_run 2>&1 | \
  grep -i "methyl\|MM\|ML"
```

### 情境 3：PERMANOVA `Valid=false`
- 可能原因：reads 不足 / groups 不足（見 [../04_parameters/03_statistical_methods.md](../04_parameters/03_statistical_methods.md)）
- 檢查：`ClusterPermanovaValid` 欄位 + reads.tsv 行數

---

## Log 輸出到檔案

```bash
./build/bin/inter_sub_mod --log-level debug ... -o run \
  2>&1 | tee run/debug/ism_verbose.log
```

---

## 與 /test-quick 結合

```bash
# Debug 模式快速驗證
./scripts/run_vcf_all_snv.sh --mode chr19-verification -- \
  --log-level debug --output-filtered-reads
```

---

## 效能注意

| Log Level | 效能影響 |
|-----------|---------|
| info（預設） | 基準 |
| debug | ~1.2x slower |
| trace | ~2-5x slower，生成大 log 檔 |

**建議**：trace 只用在單一 region 測試，不要跑全量

---

## 相關

- 過濾規則：[../04_parameters/04_filtering_rules.md](../04_parameters/04_filtering_rules.md)
- CLI 參數：[../04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)
- 原始碼：`include/utils/ArgParser.hpp:105-116`
