---
id: ism-kb-04-parameters-filtering-rules
name: "Filtering Rules"
description: "Read-level 過濾規則：mapq ≥ 20、read length ≥ 1000、base quality ≥ 20、methylation 閾值 0.8/0.2；--no-filter 驗證用途。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "filter rules against include/utils/ArgParser.hpp:60-74"
related_ids:
  - ism-kb-04-parameters-index
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-06-workflows-debug-and-logging
tags: [parameters, filter, mapq, base-quality, methylation, threshold]
canonical_paths: [04_parameters/04_filtering_rules.md]
alias_paths: []
---

# Filtering Rules

- 一句結論：4 個讀取過濾閾值（mapq, read_length, base_quality, methyl high/low）；另有 `--no-filter` 驗證用途
- 適用對象：調整過濾嚴謹度、Debug 過濾造成的 read 流失
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 輸出過濾詳細日誌
  ./build/bin/inter_sub_mod --log-level debug --output-filtered-reads ... -o debug_run
  head debug_run/debug/filtered_reads.tsv
  ```

---

## 過濾閾值總表

| 參數 | 預設 | 取值 | 意義 | 原始碼 |
|------|------|------|------|--------|
| `--min-mapq` | 20 | [0, 60] | Minimum MAPQ（reads 整體映射品質） | ArgParser.hpp:67-68 |
| `--min-read-length` | 1000 | > 0 | Minimum read length (bp) | ArgParser.hpp:70-71 |
| `--min-base-quality` | 20 | [0, 93] | SNV 位置鹼基品質 | ArgParser.hpp:73-74 |
| `--methyl-high` | 0.8 | [0.0, 1.0] | 二元化「甲基化」閾值（p ≥ 0.8） | ArgParser.hpp:60-61 |
| `--methyl-low` | 0.2 | [0.0, 1.0] | 二元化「未甲基化」閾值（p ≤ 0.2） | ArgParser.hpp:63-64 |

---

## 過濾邏輯

### Read-level 過濾順序
```
1. MAPQ < --min-mapq                        → 過濾
2. Read length < --min-read-length           → 過濾
3. SNV 位置 base quality < --min-base-quality → 過濾（該 SNV 不計入此 read）
4. 缺 MM/ML tag                              → 過濾
5. (可選) soft-clip 比例 / CIGAR 異常等      → 過濾
```

### CpG-level 二元化
```
p_methylation (從 ML tag 解碼)：
  ├── p >= --methyl-high  → "methylated" (1)
  ├── p <= --methyl-low   → "unmethylated" (0)
  └── else (0.2 < p < 0.8) → "ambiguous" (NaN, 影響 NHD / JACCARD)
```

**注意**：L1/L2/BERNOULLI 不二元化，直接用原始 p。

---

## `--no-filter` 旗標

**用途**：跳過所有 read-level 過濾，輸出所有 reads

**場景**：
- 驗證過濾導致的 read 數異常
- Debug 為何某 region 無 reads

**注意**：統計結果**不可信**，僅作診斷

---

## `--output-filtered-reads` 旗標

**用途**：輸出被過濾的 reads 清單含原因

**自動啟用條件**：`--log-level debug` 或 `trace` 時自動開啟（ArgParser.hpp:190-192）

**輸出位置**：`<debug_output_dir>/filtered_reads.tsv`

**欄位**：read_id, filter_reason, chr, pos, mapq, read_length 等

---

## 典型診斷流程

### 問題：某 region read 數太少
```bash
# 1. 關閉過濾重跑同 region
./build/bin/inter_sub_mod --no-filter --log-level debug \
  -t tumor.bam -r ref.fa -v target_region.vcf.gz -o diag_run

# 2. 對比過濾前後 read 數
diff <(wc -l diag_run/TP/region_0/reads.tsv) \
     <(wc -l prod_run/TP/region_0/reads.tsv)

# 3. 檢查過濾原因分布
awk -F'\t' '{print $2}' diag_run/debug/filtered_reads.tsv | sort | uniq -c
```

### 問題：甲基化狀態大量 NaN
- 可能原因：0.2 < p < 0.8 的位點太多（ambiguous）
- 解法：降低 `--methyl-high` 或提高 `--methyl-low`（收窄 ambiguous 區間）
- **但**：鬆綁閾值會放大雜訊，通常不建議

---

## 閾值對 F1 的影響

| 閾值調整 | 預期影響 |
|---------|---------|
| `--min-mapq 30`（嚴） | Read 少；F1 可能略升（雜訊降） |
| `--min-mapq 10`（鬆） | Read 多；F1 可能略降 |
| `--methyl-high 0.7` / `--methyl-low 0.3`（鬆） | ambiguous 變少；NHD 結果改變 |
| `--no-filter` | Read 大量增加；統計結構破壞 |

**建議**：**不要改預設**，除非有明確理論理由；Canonical benchmark 一律用預設。

---

## 與 VCF 過濾的區別

- **Read filter**（本文件）：BAM 層過濾 reads
- **VCF filter**（`scripts/pipeline/steps/03_filter_analysis.py`）：在 VCF 層分類 TP/FP

**VCF filter 需求**：truth VCF + caller VCF → 產生 `filtered_snv_tp.vcf.gz` / `filtered_snv_fp.vcf.gz`

詳見：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)

---

## 相關

- CLI 參數：[01_cli_arguments.md](01_cli_arguments.md)
- Debug 日誌：[../06_workflows/05_debug_and_logging.md](../06_workflows/05_debug_and_logging.md)
- 原始碼：`include/utils/ArgParser.hpp:60-74`
