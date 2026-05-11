---
id: ism-kb-05-data-formats-index
name: "Data Formats 目錄索引"
description: "ISM 輸出檔案格式索引：significance_summary.csv (59 欄)、master_dataset (116 欄)、per-region 檔案、VCF 來源對照、輸出目錄 layout。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "data formats directory structure"
related_ids:
  - ism-kb-05-data-formats-significance-summary-schema
  - ism-kb-05-data-formats-master-dataset-schema
  - ism-kb-05-data-formats-per-region-files
  - ism-kb-05-data-formats-vcf-sources
  - ism-kb-05-data-formats-output-directory-layout
tags: [data-format, index, schema, output]
canonical_paths: [05_data_formats/00_index.md]
alias_paths: []
---

# Data Formats 目錄索引

- 一句結論：主要產出 significance_summary.csv (59 欄) + per-region 目錄；Python 聚合後成 master_dataset (116 欄)
- 適用對象：分析 ISM 輸出、寫後處理 script、改輸出格式
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  head -1 /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/*/significance_summary.csv | tr ',' '\n' | head
  ```

---

## 文件列表

| 檔案 | 主題 | 欄位數 |
|------|------|--------|
| [01_significance_summary_schema.md](01_significance_summary_schema.md) | C++ 主要輸出 | 59 |
| [02_master_dataset_schema.md](02_master_dataset_schema.md) | Python 聚合 | 116 |
| [03_per_region_files.md](03_per_region_files.md) | region 子目錄內容 | — |
| [04_vcf_sources.md](04_vcf_sources.md) | 7 樣本 VCF 符號連結對照 | — |
| [05_output_directory_layout.md](05_output_directory_layout.md) | canonical / synthesis / archive 頂層 | — |
| [06_merged_dataset_pitfalls.md](06_merged_dataset_pitfalls.md) ⚠ | **Merged 合成檔已知陷阱（AF 非 caller_af + HCC1395 phase1_new LOH 殘缺）** | 2 pitfalls |

---

## 檔案關係

```
ISM C++ binary
      │
      ├── significance_summary.csv (59 cols, per run)
      ├── TP/region_<N>/ {reads.tsv, methylation/, distance_<METRIC>/}
      └── FP/region_<N>/ {同上}

Python 後處理（scripts/analysis/build_loh_round1_cross_sample_audit.py）
      │
      └── all_region_rows.tsv.gz (116 cols, cross-run aggregated)
```

---

## 權威來源

本 KB 提供概覽；詳細欄位字典在 `docs/data_specs/`：
- `docs/data_specs/20260411_significance_summary欄位字典_01.md`（59 欄完整）
- `docs/data_specs/20260411_master_dataset欄位字典_01.md`（116 欄完整）

---

## 相關

- Pipelines 輸出：[../03_pipelines/](../03_pipelines/)
- 距離度量產出：[../04_parameters/02_distance_metrics.md](../04_parameters/02_distance_metrics.md)
