<!--
建立時間: 2026-04-14
目標: v1 錯誤數據暫存區，待清理刪除
-->

# to_pipeline_staging v1 — 已棄用（待刪除）

這些檔案使用了錯誤的數據來源，已於 2026-04-14 從 `research/to_pipeline_staging/` 移至此處。

## 錯誤原因

v1 混用三個不一致的數據來源：

| 分析層 | v1 錯誤來源 | 正確來源 |
|--------|-----------|---------|
| VCF | Clair 團隊 run `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz` (ONT BAM 無 MM/ML, 47,798 PASS, Sep 2025) | canonical TO pipeline `step01_clairs_to/snv.vcf.gz` (ONT_5kHz BAM 有 5mCG+5hmCG, 48,085 PASS, Mar 2026) |
| ISM | `three_way_comparison/tumor_only_full/` (112,387 rows, 未分 TP/FP) | canonical `step05_intersubmod/` (TP 28,383 + FP 11,830) |
| BED | 未限制 truth BED → 19,689 FP | 限制 truth BED → 11,843 FP |

## 錯誤檔案清單

```
scripts/
├── 01_multi_stage_characterization.py   # 使用 zhenyu112 VCF + three_way ISM
├── 02_generate_plots.py                 # 讀取 v1 數據 (hcc1395_pass_*)
└── 03_updated_waterfall_plot.py         # 讀取 v1 數據

data/
├── hcc1395_pass_multimodal.csv          # 基於錯誤 VCF + ISM 的合併表
├── hcc1395_feature_auc_by_stage.csv     # 基於錯誤數據的 AUC
└── hcc1395_stage_metrics.json           # 錯誤的 F1 (0.649 而非 0.7127)

figures/
├── 01_stage_pipeline_f1.png
├── 02_feature_auc_by_stage.png
├── 03_cnv_loh_stratification.png
├── 04_vcf_feature_distributions.png
├── 05_ism_feature_distributions.png
├── 06_multimodal_summary.png
├── 07_af_cnv_scatter.png
├── 08_corrected_waterfall_and_middle_block.png
└── 09_middle_block_ism_features.png

reports/
└── 20260412_TO_pipeline_multi_stage_characterization_01.md  # v1 報告
```

## 正確版本

校正後的 v2 檔案位於：`research/to_pipeline_staging/`（scripts 04-07, C*/M* figures, canonical data）
