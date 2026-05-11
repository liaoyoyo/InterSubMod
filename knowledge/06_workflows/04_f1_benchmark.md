---
id: ism-kb-06-workflows-f1-benchmark
name: "F1 Benchmark Workflow"
description: "F1 比對協議：som.py（hap.py 框架）+ pipeline/run_benchmark.sh；公式 2*P*R/(P+R)；典型 command。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: howto
verified_scope: "F1 formula against scripts/pipeline/steps/03_filter_analysis.py:229-234"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-08-truth-and-benchmark-truth-set-registry
  - ism-kb-08-truth-and-benchmark-benchmark-protocols
  - ism-kb-06-workflows-batch-tp-fp-analysis
tags: [workflow, f1, benchmark, som.py, hap.py]
canonical_paths: [06_workflows/04_f1_benchmark.md]
alias_paths: []
---

# F1 Benchmark Workflow

- 一句結論：本專案 F1 由 `scripts/pipeline/steps/03_filter_analysis.py:229-234` 計算；somatic benchmark 用 som.py
- 適用對象：執行 F1 主表、論文數據重現
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 本專案 F1 計算函式
  grep -n "def compute_metrics" /big7_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/steps/03_filter_analysis.py
  ```

---

## F1 計算公式

**位置**：`scripts/pipeline/steps/03_filter_analysis.py:229-234`

```python
def compute_metrics(tp, fp, truth_total):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / truth_total
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"tp": tp, "fp": fp, "precision": precision, "recall": recall, "f1": f1}
```

**參數**：
- `tp`：LongPhase-S 產出的 `filtered_snv_tp.vcf.gz` 變異數
- `fp`：`filtered_snv_fp.vcf.gz` 變異數
- `truth_total`：整個 truth VCF 的變異數（**不扣 HC BED 之外**的變異）

---

## Pipeline 執行

### 內建：LongPhase-S `--truth-vcf` 旗標
```bash
./scripts/pipeline/run_benchmark.sh \
  --mode s-pure \
  --sample HCC1395
```
**內部**：haplotag 時同時比對 truth VCF 分類 TP/FP → 產生 filtered_snv_{tp,fp}.vcf.gz

### 外部：som.py（hap.py 框架，論文 gold standard）

```bash
# 完整命令見 docs/workflows/benchmark-workflow.md
som.py \
  <truth.vcf.gz> \
  <caller.vcf.gz> \
  -R <HC_regions.bed> \
  -r <reference.fa> \
  -o out/som_py_result
```

**用途**：跨工具公平比較、論文表格
**輸出**：`out/som_py_result.stats.csv`（F1、precision、recall）

---

## 跨 run 聚合

```bash
python3 scripts/analysis/compare_benchmark_f1.py \
  --run-dir output/canonical/HCC1395/paired_full/20260314_*/ \
  --output-tsv benchmark_comparison.tsv
```

**輸出**：`benchmark_comparison.tsv`（dataset_id, f1_before, f1_after, delta_f1）

---

## 三種 Benchmark 情境

### 情境 1：跑一次 canonical paired_full
```bash
./scripts/pipeline/run_benchmark.sh --mode s-pure --sample HCC1395
```

### 情境 2：paired_pileup 對照
```bash
./scripts/pipeline/run_benchmark.sh --mode s-pure-pileup --sample HCC1395
```

### 情境 3：TO 對照
```bash
./scripts/pipeline/run_benchmark.sh --mode to-pure --sample HCC1395 \
  --vcf-source pileup
```

---

## TP / FP / FN 定義

- **TP**：caller 輸出 ∩ truth VCF（適用 HC BED 篩選後）
- **FP**：caller 輸出 ∖ truth VCF
- **FN**：truth VCF 有 + caller 未輸出（需 `bcftools isec` 推得）

```bash
# 手動計算 FN
bcftools isec -p isec_output <truth.vcf.gz> <caller.vcf.gz>
# 0000.vcf = truth only = FN candidates
```

---

## 7 樣本 Truth Set

完整路徑表：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)

**摘要**：
| 樣本 | Truth Set | TRUTH_TOTAL |
|------|----------|-------------|
| HCC1395 | SEQC2 high-confidence | 39,447 |
| COLO829 | NYGC | 41,427 |
| H1437 | Orthogonal tools | 90,129 |
| H2009 | Orthogonal tools | 168,529 |
| HCC1937 | Orthogonal tools | 60,691 |
| HCC1954 | Orthogonal tools | 26,567 |

---

## Phase 1A Locked Numbers

- **paired_full ΔF1**：+0.0112（locked）
- **TO ΔF1**：-0.0206（NEGATIVE locked）

👉 **完整 provenance（CI、方法、60+ features × 748K regions、運行條件、限制）** → [../03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)（SoT）

---

## 常見陷阱

- ❌ **用 som.py 結果 vs 內建 F1 直接比較**：兩者細節不同（GT 匹配、normalization 等）
- ❌ **不用 HC BED 過濾就算 F1**：FP 會暴增，F1 假低
- ❌ **truth_total 用錯**：Paired 與 TO 共用同一 truth_total，**不是**各自的 TP+FN

---

## 相關

- Truth set 全表：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)
- F1 計算細節：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
- Benchmark protocols：[../08_truth_and_benchmark/03_benchmark_protocols.md](../08_truth_and_benchmark/03_benchmark_protocols.md)
- 原始碼：`scripts/pipeline/steps/03_filter_analysis.py`
