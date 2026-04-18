<!--
建立時間: 2026-04-13 19:00
更新時間: 2026-04-14 
目標: TO Pipeline 多階段 TP/FP 特徵分析研究工作區（v2 校正版）
-->

# TO Pipeline Multi-Stage Characterization (v2)

分析 ClairS-TO → LongPhase-TO → ISM 三階段中 TP/FP 的多模態特徵差異。

> **v1 已完全棄用並移至 `docs/trash/to_pipeline_staging_v1/`。**
> v1 錯誤根因：使用 zhenyu112 VCF + three_way_comparison ISM + 未限制 truth BED，
> 導致三層數據來源不一致（詳見下方「v1 → v2 校正說明」）。

---

## 數據產生流程

### 原始數據來源（全部來自 canonical pipeline）

```
/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/research_rounds/20260307_hcc1395_to_pilot_1/
├── step01_clairs_to/
│   └── snv.vcf.gz                          # ClairS-TO VCF (liaoyoyo2001, ONT_5kHz, Mar 2026)
│                                            # BAM: /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam
│                                            # 3,203,048 total → 48,085 PASS
├── step02_benchmark_clairs_to/
│   ├── tp.vcf                              # 28,396 TP (bcftools isec with truth BED)
│   ├── fp.vcf                              # 11,843 FP (bcftools isec with truth BED)
│   └── variant_counts.txt                  # 計數確認
├── step03_longphase_to/
│   ├── tumor_phased.vcf.gz                 # LongPhase-TO 輸出（添加 PS/GT2/GT3 FORMAT）
│   └── tumor_phased_LOH.bed               # 1,100 LOH regions
└── step05_intersubmod/
    ├── intersubmod_tp/significance_summary.csv  # 28,383 TP ISM 結果
    └── intersubmod_fp/significance_summary.csv  # 11,830 FP ISM 結果
```

### 其他 6 樣本 ClairS-TO VCF（由 zhenyu112 產生，無 ISM）

```
/big8_disk/data/{SAMPLE}/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
  - HCC1395_DORADO: /big8_disk/data/HCC1395/ONT_Dorado/ClairS_TO_v0_3_0/snv.vcf.gz
  - COLO829:        /big8_disk/data/COLO829/ONT_PAO/ClairS_TO_v0_3_0/snv.vcf.gz
  - H1437:          /big8_disk/data/H1437/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
  - H2009:          /big8_disk/data/H2009/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
  - HCC1937:        /big8_disk/data/HCC1937/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
  - HCC1954:        /big8_disk/data/HCC1954/ONT/ClairS_TO_v0_3_0/snv.vcf.gz
```

> 注意：這 6 樣本僅有 ClairS-TO VCF，無 LongPhase-TO 或 ISM 下游分析。
> 多樣本比較僅使用 VCF-level 特徵（AF, QUAL, GQ, DP, SB, haplotype_flag）。

---

## 數據產生流程圖

```
[Script 04] 04_canonical_multi_stage_analysis.py
  輸入: step01 VCF + step02 TP/FP + step03 phased VCF + step05 ISM
  輸出: hcc1395_canonical_multimodal.csv (40,239 rows)
        hcc1395_canonical_feature_auc.csv (38 features)
        hcc1395_canonical_stage_metrics.json

[Script 05] 05_multi_sample_to_benchmark.py
  輸入: HCC1395 canonical VCF + 6 樣本 zhenyu112 VCFs + 各樣本 truth sets
  輸出: multi_sample_to_summary.json
        multi_sample_to_features.csv (420,254 rows)

[Script 06] 06_canonical_plots.py
  輸入: hcc1395_canonical_*.csv/json
  輸出: C01~C06 圖表

[Script 07] 07_multi_sample_plots.py
  輸入: multi_sample_to_*.csv/json
  輸出: M01~M04 圖表
```

---

## 核心結論

1. **ISM 隱性過濾在 canonical 數據中可忽略** — 僅移除 26 variants (13 TP + 13 FP)
2. **ISM SuggestFilter 反效果** — 移除 124 TP + 95 FP，F1 下降 0.0012
3. **LongPhase-TO 不改變 PASS/FILTER** — 只添加 phasing info
4. **AF 是最有效的 VCF 特徵** — AUC=0.666 (inverted)，跨 6/7 樣本一致
5. **H-flag 最穩定** — 7/7 樣本 AUC > 0.74
6. **TO F1 受樣本特性影響極大** — 0.378 (HCC1954) 到 0.899 (H2009)

---

## v1 → v2 校正說明

| 項目 | v1（已棄用） | v2（校正） | 錯誤原因 |
|------|-------------|-----------|---------|
| ClairS-TO VCF | zhenyu112 (47,798 PASS) | liaoyoyo2001, ONT_5kHz canonical (48,085 PASS) | 不同 BAM、不同 run |
| ISM 來源 | three_way_comparison (112K rows, 未分 TP/FP) | canonical step05 (TP 28,383 + FP 11,830) | 管線來源不對應 |
| TP/FP 分類 | 自行 isec 未限 BED → 19,689 FP | step02 pre-benchmarked 限 BED → 11,843 FP | FP 膨脹 66% |
| ISM FP 移除率 | 39.9% | 0.11% | 三層錯誤疊加 |
| Neutral TP rate | 16.9% | 63.6% | 三層錯誤疊加 |
| 基線 F1 | 0.649 | 0.7127 | 三層錯誤疊加 |

v1 檔案已移至：`docs/trash/to_pipeline_staging_v1/`

---

## 目錄結構

```
scripts/
├── 04_canonical_multi_stage_analysis.py    # HCC1395 canonical 三階段分析
├── 05_multi_sample_to_benchmark.py         # 7 樣本 ClairS-TO benchmark
├── 06_canonical_plots.py                   # HCC1395 canonical 圖表 (C01-C06)
└── 07_multi_sample_plots.py                # 多樣本比較圖表 (M01-M04)

data/
├── hcc1395_canonical_multimodal.csv        # 40,239 variants 全特徵
├── hcc1395_canonical_feature_auc.csv       # 38 特徵 AUC
├── hcc1395_canonical_stage_metrics.json    # 三階段 F1 指標
├── multi_sample_to_summary.json            # 7 樣本 benchmark 結果
└── multi_sample_to_features.csv            # 7 樣本 VCF 特徵 (420K rows)

figures/
├── C01_canonical_pipeline_waterfall.png    # Pipeline waterfall
├── C02_canonical_feature_auc_ranking.png   # 全特徵 AUC 排名
├── C03_canonical_vcf_features.png          # VCF 特徵分布
├── C04_canonical_cnv_stratification.png    # CNV/LOH 分層
├── C05_canonical_ism_features.png          # ISM 特徵分布
├── C06_canonical_ism_filter_detail.png     # ISM filter 細節
├── M01_multi_sample_f1_comparison.png      # 7 樣本 F1 比較
├── M02_multi_sample_af_distribution.png    # AF 分布比較
├── M03_multi_sample_feature_auc_heatmap.png # 特徵 AUC 熱圖
└── M04_to_vs_paired_f1.png                 # TO vs Paired 比較

reports/
└── 20260413_TO_pipeline_canonical_analysis_01.md  # v2 校正版報告
```
