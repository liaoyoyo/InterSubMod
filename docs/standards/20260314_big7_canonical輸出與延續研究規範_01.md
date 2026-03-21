<!--
建立時間: 2026-03-14
目標: 定義 big7 canonical output、archive/synthesis 分層、命名與完整性判定
-->

# Big7 Canonical 輸出與延續研究規範

## 1. 根目錄

InterSubMod 在本機器上的主輸出根目錄固定為：

- canonical：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical`
- synthesis：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis`
- archive physical roots：
  - `/big7_disk/liaoyoyo2001/big7_disk_output/big8_output_archive`
  - `/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive`

## 2. Canonical 結構

正式 sample/mode/run bundle 一律放在：

`canonical/{sample}/{canonical_mode}/{run_id}/`

其中 `canonical_mode` 固定只允許：

- `paired_full`
- `paired_pileup`
- `to_full`
- `to_pileup`

每個 run bundle 至少包含：

- `manifest/run_context.json`
- `manifest/file_manifest.tsv`
- `manifest/upstream_dependency.tsv`
- `manifest/migration_provenance.tsv`
- `metrics/metrics.json`
- `metrics/benchmark_comparison.tsv`
- `metrics/completeness.tsv`
- `reports/README.md`

## 3. 命名規範

`run_id` 固定格式：

`{date}_{sample}_{canonical_mode}_{caller_model}`

範例：

- `20260307_HCC1395_paired_full_full`
- `20260212_H1437_paired_pileup_pileup`
- `20260307_hcc1395_to_pilot_1_HCC1395_to_pileup_pileup`

## 4. 不可做的事

1. 不把 symlink 指向的來源檔再複製成新的大型輸出。
2. 不把 archive/scratch/debug 目錄直接當 canonical baseline run。
3. 不只在報告文字中寫 F1，而不輸出 `TP/FP/FN/precision/recall/F1` 表格。

## 5. 完整性判定

`completeness.tsv` 為唯一 canonical 完整性判定入口，最少欄位：

- `tagged_bam_ready`
- `truth_ready`
- `baseline_ready`
- `intersubmod_ready`
- `metrics_ready`
- `decision_metrics_ready`
- `completeness_state`
- `blocking_reason`

現階段 `tagged_bam_ready=false` 的舊 run 不等於結果失效，只表示 archive 中未保留實體 tagged BAM，需要後續正式重跑才能升級為 completed。
