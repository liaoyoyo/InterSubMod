---
id: ism-kb-08-truth-and-benchmark-benchmark-protocols
name: "Benchmark Protocols"
description: "som.py (hap.py 框架) vs 本專案內建 F1 計算；何時用何者、輸出差異、論文標準方法。"
status: active
last_verified: 2026-04-22
content_nature: reference
doc_type: reference
verified_scope: "benchmark protocols comparison against docs/workflows/benchmark-workflow.md"
related_ids:
  - ism-kb-08-truth-and-benchmark-index
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-06-workflows-f1-benchmark
tags: [benchmark, som.py, hap.py, protocol, canonical]
canonical_paths: [08_truth_and_benchmark/03_benchmark_protocols.md]
alias_paths: []
---

# Benchmark Protocols

- 一句結論：som.py（hap.py 框架）是跨工具公平比較的 gold standard；本專案內建 F1 用於 pipeline 內部
- 適用對象：論文 benchmark、跨 caller 比較
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  som.py --version  # 應為 hap.py v0.3.15 或以上
  ```

---

## 兩種 Benchmark 協議對照

| 面向 | som.py (hap.py) | 本專案內建 F1 |
|------|----------------|---------------|
| **來源** | Illumina hap.py repo | `scripts/pipeline/steps/03_filter_analysis.py` |
| **用途** | 跨工具比較、論文 | pipeline 內部 F1 計算 |
| **GT 匹配** | 嚴格（genotype-aware） | 相對寬鬆（變異位置匹配） |
| **Normalization** | bcftools norm + 左對齊 | 視 pipeline 預處理 |
| **Region 篩選** | 用 `-R <BED>` | 用 HC BED 過濾 |
| **輸出格式** | `.stats.csv` + `.vcf.gz` | `benchmark_comparison.tsv` |
| **時間** | 10-20 分鐘 | 1-2 分鐘 |
| **推薦場景** | 論文表格、正式 benchmark | Iterative 驗證、Quick check |

---

## som.py（外部 gold standard）

### 安裝
```bash
# hap.py 通常透過 docker
docker pull pkrusche/hap.py

# 或 conda
conda install -c bioconda hap.py
```

### 標準命令
```bash
som.py \
  <truth.vcf.gz> \
  <caller.vcf.gz> \
  -R <HC_regions.bed> \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -o out/som_py_result
```

### 輸出
```
out/som_py_result/
├── .stats.csv          # F1、precision、recall 主表
├── .counts.csv         # TP/FP/FN 原始計數
├── .vcf.gz             # 分類後的 VCF（含 TP/FP/FN tag）
└── .summary.csv        # 簡要摘要
```

### 論文用途
- Methods 章節：
  > "Somatic SNV benchmarking was performed using som.py (v0.3.15) with default settings, using GIAB v2.5 / SEQC2 v1.2.1 as truth sets, restricted to high-confidence regions provided by each consortium."

---

## 本專案內建 F1

### 調用
```bash
./scripts/pipeline/run_benchmark.sh \
  --mode s-pure \
  --sample HCC1395
```

### 公式與原理
見 [02_f1_calculation.md](02_f1_calculation.md)

### 何時用
- Pipeline 內部迭代驗證（快）
- ISM filter 前後 Δ F1 計算
- Iterative research（Phase 1A 迭代）

### 何時不用
- 論文 benchmark 表格 → 用 som.py
- 跨 caller 公平比較 → 用 som.py（因 normalization 一致）

---

## Pipeline 整合方式

### 1. ClairS 原始 VCF 分類
```
Caller VCF
    │
    ├─ vs truth (內建 03_filter_analysis.py)
    │      └─ filtered_snv_tp.vcf.gz, filtered_snv_fp.vcf.gz
    │
    └─ vs truth (som.py)
           └─ som_py_result.stats.csv
```

兩種協議**並行執行**，內建用於 iterative，som.py 用於最終報告。

### 2. ISM filter 後的 F1
ISM filter 輸出篩選後的 VCF → 同樣用內建 03_filter_analysis.py + som.py 兩協議計算

---

## Δ F1 的兩種算法

### 方法 A：內建 F1
```
Δ F1 = F1_ISM_filtered - F1_baseline
     = (from pipeline internal) - (from pipeline internal)
```

### 方法 B：som.py F1
```
Δ F1 = F1_ISM_filtered (som.py) - F1_baseline (som.py)
```

**兩者差異**：
- 絕對 F1 可能略不同（~0.01 範圍）
- Δ F1 方向一致
- 論文應用 som.py 報告

---

## 7 樣本 benchmark 批次

```bash
for sample in HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954; do
  # 1. 跑 pipeline + ISM
  ./scripts/pipeline/run_benchmark.sh --sample $sample --mode s-pure

  # 2. 跑 som.py
  som.py \
    <truth.vcf.gz> \
    output/canonical/$sample/paired_full/*/caller.vcf.gz \
    -R <hc.bed> -r <ref.fa> \
    -o benchmark/${sample}_sompy
done
```

**時間**：~2-3 小時（7 樣本連跑 som.py）

---

## Pipeline 資源需求

| 協議 | CPU | RAM | Disk I/O |
|------|-----|-----|---------|
| 內建 F1 | 低 | 低 | 中 |
| som.py | 中 | 中 | 低 |

---

## 推薦做法

1. **開發階段**：用內建 F1 快速迭代
2. **正式驗證**：跑 som.py 二次確認
3. **論文表格**：som.py 為主，內建 F1 為輔（提 supplementary）
4. **Δ F1 報告**：兩種協議各給一份，確認方向一致

---

## 常見誤區

- ❌ **混合兩協議的絕對 F1**：同一表格不可混用
- ❌ **省略 HC BED 篩選**：som.py 的 `-R` 是標準做法
- ❌ **忘記 normalization**：som.py 內建，但內建 F1 依賴 pipeline 預處理

---

## 權威文件

👉 **`docs/workflows/benchmark-workflow.md`**（若存在）
👉 **`docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`**

---

## 相關

- Truth set：[01_truth_set_registry.md](01_truth_set_registry.md)
- F1 公式：[02_f1_calculation.md](02_f1_calculation.md)
- F1 Workflow：[../06_workflows/04_f1_benchmark.md](../06_workflows/04_f1_benchmark.md)
