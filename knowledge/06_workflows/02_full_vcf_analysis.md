---
id: ism-kb-06-workflows-full-vcf-analysis
name: "Full VCF Analysis Workflow"
description: "`run_vcf_all_snv.sh` 的三個 mode：all-with-w5000（完整 5kb 視窗）、chr19-verification（快速驗證）、自定義 mode；含 caller-mode 切換。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: howto
verified_scope: "scripts/run_vcf_all_snv.sh mode options"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-06-workflows-batch-tp-fp-analysis
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-06-workflows-build-and-test
tags: [workflow, vcf-analysis, run_vcf_all_snv, mode]
canonical_paths: [06_workflows/02_full_vcf_analysis.md]
alias_paths: []
---

# Full VCF Analysis Workflow

- 一句結論：`./scripts/run_vcf_all_snv.sh --mode <mode> --caller-mode <paired|to>` 為主要執行入口
- 適用對象：執行 ISM 分析、Debug pipeline
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  /big7_disk/liaoyoyo2001/InterSubMod/scripts/run_vcf_all_snv.sh --help 2>&1 | head -30
  ```

---

## 腳本概要

**位置**：`scripts/run_vcf_all_snv.sh`
**功能**：單一樣本、單 VCF、執行 ISM 分析
**常用模式**：`all-with-w5000`（canonical）、`chr19-verification`（快速）

---

## 三個主要 Mode

### 1. `all-with-w5000`（Canonical）
```bash
./scripts/run_vcf_all_snv.sh \
  --caller-mode paired \
  --mode all-with-w5000 \
  --metrics BERNOULLI \
  --threads 120 \
  -o output/paired_full_run \
  --plot-type all
```
**用途**：完整 benchmark run（論文標準 -w 5000）
**時間**：~5 分鐘（HCC1395 TP+FP）

### 2. `chr19-verification`（Quick）
```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```
**用途**：快速驗證主流程（chr19 僅 10-20 個 SNV）
**時間**：<30 秒

### 3. 自定義 mode
```bash
# 用 --help 看完整選項
./scripts/run_vcf_all_snv.sh --help
```

---

## `--caller-mode` 切換

| 值 | 使用 VCF | 對應 pipeline |
|----|---------|--------------|
| `paired` | ClairS paired | paired_full / paired_pileup |
| `to` | ClairS-TO | tumor_only |

**範例**：
```bash
# paired
./scripts/run_vcf_all_snv.sh --caller-mode paired --mode all-with-w5000

# tumor-only
./scripts/run_vcf_all_snv.sh --caller-mode to --mode all-with-w5000 \
  -v /big8_disk/liaoyoyo2001/longphase-to-mod/output/clairsto_v3fixed_work/clairsto_tp.vcf.gz \
  -t /big8_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam
```

---

## 輸出結構

```
<output_dir>/
├── significance_summary.csv      # 59 欄（見 05_data_formats/01）
├── TP/                           # TP VCF 分析
│   └── region_<N>/ ...
├── FP/                           # FP VCF 分析
│   └── region_<N>/ ...
├── logs/
│   └── ism_run_<timestamp>.log
└── (可選) figures/               # 若 --plot-type all
```

詳見 [../05_data_formats/05_output_directory_layout.md](../05_data_formats/05_output_directory_layout.md)

---

## 常用組合

### 單次 BERNOULLI benchmark
```bash
./scripts/run_vcf_all_snv.sh \
  --mode all-with-w5000 \
  --metrics BERNOULLI \
  -j 120 -o out/bernoulli
```

### 多度量並行（一次跑完）
```bash
./scripts/run_vcf_all_snv.sh \
  --mode all-with-w5000 \
  --metrics "NHD,L1,BERNOULLI" \
  -j 120 -o out/multi_metric
```

### 加上圖表
```bash
./scripts/run_vcf_all_snv.sh \
  --mode all-with-w5000 \
  --plot-type all \       # 距離熱圖 + 分群熱圖
  -o out/with_plots
```

---

## 相關腳本

- **Batch（TP+FP 並行）**：[03_batch_tp_fp_analysis.md](03_batch_tp_fp_analysis.md)
- **F1 benchmark**：[04_f1_benchmark.md](04_f1_benchmark.md)
- **原始腳本**：`scripts/run_vcf_all_snv.sh`

---

## Debug 模式

```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification \
  -- \
  --log-level debug \
  --output-filtered-reads
```

詳見 [05_debug_and_logging.md](05_debug_and_logging.md)

---

## 常見問題

| 問題 | 解法 |
|------|------|
| `binary not found` | 先 `cd build && make -j` |
| `VCF not indexed` | `tabix -p vcf <vcf.gz>` |
| Thread 太高 OOM | 降 `-j`（ISM per-thread 記憶體高） |
| Output dir 已存在 | 腳本會拒絕覆蓋；手動刪或用新目錄 |

---

## 相關

- Batch 腳本：[03_batch_tp_fp_analysis.md](03_batch_tp_fp_analysis.md)
- F1 計算：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
- 原始碼：`scripts/run_vcf_all_snv.sh`
