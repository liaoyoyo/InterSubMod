---
id: ism-kb-03-pipelines-paired-full
name: "Paired Full Pipeline（Canonical Benchmark）"
description: "ISM canonical benchmark pipeline：ClairS paired full VCF + LongPhase-S haplotag + ISM 分析。Phase 1A delta F1=+0.0112 locked。含輸入鏈路、完整命令、輸出結構、驗證方法。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "pipeline command against scripts/run_vcf_all_snv.sh; file:line refs against HEAD"
related_ids:
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-pipeline-comparison
  - ism-kb-04-parameters-cli-arguments
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-05-data-formats-output-directory-layout
  - ism-kb-02-samples-hcc1395
  - ism-kb-03-pipelines-paired-pileup
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [pipeline, paired, full, canonical, benchmark, clairs, longphase-s]
canonical_paths: [03_pipelines/01_paired_full.md]
alias_paths: []
---

# Paired Full Pipeline（Canonical Benchmark）

- 一句結論：ClairS paired full VCF + LongPhase-S haplotag + ISM，本專案 F1 主表的唯一來源
- 適用對象：Benchmark 執行者、論文數據重現者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cd /big7_disk/liaoyoyo2001/InterSubMod
  ./scripts/run_vcf_all_snv.sh --caller-mode paired --mode chr19-verification
  ```

---

## 輸入鏈路

```
┌─────────────────┐   ┌─────────────────┐   ┌──────────────┐
│ Tumor BAM       │   │ Normal BAM       │   │ Reference     │
│ (HP tagged by    │   │ (unphased)       │   │ hg38.fa       │
│  LongPhase-S)    │   │                  │   │               │
└────────┬────────┘   └────────┬─────────┘   └──────┬───────┘
         │                      │                    │
         ▼                      ▼                    ▼
┌──────────────────────────────────────────────────────────┐
│  ClairS paired full                                       │
│  → filtered_snv_tp.vcf.gz / filtered_snv_fp.vcf.gz       │
└──────────────────────────┬───────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│  LongPhase-S                                              │
│  → tumor.haplotagged.bam + phased.bed (LOH)               │
└──────────────────────────┬───────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────┐
│  InterSubMod (build/bin/inter_sub_mod)                    │
│  → significance_summary.csv (59 cols)                     │
│  → per-region: reads.tsv, methylation.csv, distance/      │
└──────────────────────────────────────────────────────────┘
```

---

## 完整命令（直接執行 binary）

```bash
./build/bin/inter_sub_mod \
  -t <tumor_haplotagged.bam> \
  -n <normal.bam> \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -v <filtered_snv_tp.vcf.gz> \
  --loh-bed <phased.bed> \
  -w 5000 \
  -j 120 \
  --distance-metric BERNOULLI \
  --min-common-coverage 3 \
  -o output/canonical/<sample>/paired_full/<timestamp>_TP
```

**關鍵參數**：
- `-w 5000`：±5kb 視窗（Phase 1A canonical）
- `-j 120`：120 threads（OpenMP）
- `--distance-metric BERNOULLI`：預設距離度量（confidence-weighted）
- `--min-common-coverage 3`：最小 CpG 重疊

完整參數見 [../04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)。

---

## 包裝腳本命令（推薦）

```bash
./scripts/run_vcf_all_snv.sh \
  --caller-mode paired \
  --mode all-with-w5000 \
  --metrics BERNOULLI \
  --threads 120 \
  -o output/paired_full_run \
  --plot-type all
```

**時間**：~5 分鐘（全量 HCC1395 TP+FP）

---

## 輸出結構

```
output/canonical/<sample>/paired_full/<timestamp>/
├── significance_summary.csv         # 59 欄，所有 regions
├── TP/
│   ├── region_0/
│   │   ├── reads.tsv
│   │   ├── methylation/
│   │   │   ├── methylation.csv
│   │   │   ├── methylation_forward.csv
│   │   │   └── methylation_reverse.csv
│   │   └── distance_BERNOULLI/
│   │       ├── matrix.csv
│   │       ├── matrix_forward.csv
│   │       └── matrix_reverse.csv
│   └── region_<N>/ ...
└── FP/
    └── <same structure>
```

欄位字典：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)

---

## 輸入檔案權威來源

| 檔案 | 典型路徑 | 生成工具 |
|------|---------|---------|
| Tumor BAM (haplotagged) | `/big8_disk/liaoyoyo2001/longphase-s/output/<sample>/paired_v5/tumor_tagged.bam` | LongPhase-S |
| Normal BAM | `/big8_disk/data/<SAMPLE>/bam/normal.bam` | 原始樣本 |
| Filtered VCF (TP) | LongPhase-S 對 truth VCF 比對後分類 | `scripts/pipeline/steps/03_filter_analysis.py` |
| LOH BED | `<longphase_s_output>/phased.bed` | LongPhase-S |
| Reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta` | — |

---

## 可信度與驗證

**資料信任度**：🟢 Canonical（論文標準）

**Phase 1A locked 結論**（2026-04-17 凍結）：
- Paired-pure ΔF1 = **+0.0112**（95% CI [+0.0044, +0.0188]）
- **完整 provenance（方法、樣本、驗證、限制）** → [05_f1_baseline_canonical.md#cl-002](05_f1_baseline_canonical.md)（SoT）

**驗證方式**：
```bash
# 1. 跑短路徑測試（<30秒）
./scripts/run_vcf_all_snv.sh --mode chr19-verification --caller-mode paired

# 2. 跑 F1 benchmark
python3 scripts/analysis/compare_benchmark_f1.py \
  --run-dir output/canonical/HCC1395/paired_full/<timestamp>/ \
  --output-tsv benchmark_comparison.tsv

# 3. 比對上一次 canonical：F1 差異應 < 0.01
```

---

## 已知 canonical run 位置（7 樣本）

| 樣本 | Canonical Run 路徑 |
|------|-------------------|
| HCC1395 | `output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/` |
| HCC1395_DORADO | `output/canonical/HCC1395_DORADO/paired_full/20260315_*` |
| COLO829 | `output/canonical/COLO829/paired_full/` |
| H1437 | `output/canonical/H1437/paired_full/` |
| H2009 | `output/canonical/H2009/paired_full/`（負向案例：caller 已完美） |
| HCC1937 | `output/canonical/HCC1937/paired_full/` |
| HCC1954 | `output/canonical/HCC1954/paired_full/` |

---

## 與其他 pipeline 對比

| 面向 | paired_full | paired_pileup | TO |
|------|-------------|---------------|-----|
| VCF 完整性 | Full (pileup + full) | 僅 pileup | TO |
| Haplotag bias | 低 | 低 | **高（94.6% HP1）** |
| ISM ΔF1 | **+0.0112** | 類似 paired_full | **-0.0206** |
| 適用場景 | Canonical | 模型評估 | TO 場景探索 |

詳見 [04_pipeline_comparison.md](04_pipeline_comparison.md)

---

## 相關

- CLI 參數：[../04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)
- 距離度量：[../04_parameters/02_distance_metrics.md](../04_parameters/02_distance_metrics.md)
- 執行流程：[../06_workflows/02_full_vcf_analysis.md](../06_workflows/02_full_vcf_analysis.md)
- F1 計算：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
- 輸出欄位：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)
