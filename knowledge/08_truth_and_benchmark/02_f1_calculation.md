---
id: ism-kb-08-truth-and-benchmark-f1-calculation
name: "F1 Calculation"
description: "F1 公式、計算腳本位置 (`scripts/pipeline/steps/03_filter_analysis.py:229-234`)、TP/FP/FN 定義、Δ F1 計算方式。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: reference
verified_scope: "F1 formula against scripts/pipeline/steps/03_filter_analysis.py:229-234"
related_ids:
  - ism-kb-08-truth-and-benchmark-index
  - ism-kb-08-truth-and-benchmark-truth-set-registry
  - ism-kb-08-truth-and-benchmark-benchmark-protocols
  - ism-kb-06-workflows-f1-benchmark
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-f1-baseline-canonical
tags: [f1, formula, calculation, precision, recall, delta-f1]
canonical_paths: [08_truth_and_benchmark/02_f1_calculation.md]
alias_paths: []
---

# F1 Calculation

- 一句結論：F1 = 2 × precision × recall / (precision + recall)；定義於 `scripts/pipeline/steps/03_filter_analysis.py:229-234`
- 適用對象：計算 ISM 效能、重現論文數據
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -n "def compute_metrics" /big7_disk/liaoyoyo2001/InterSubMod/scripts/pipeline/steps/03_filter_analysis.py
  ```

---

## F1 公式（Python 實作）

**位置**：`scripts/pipeline/steps/03_filter_analysis.py:229-234`

```python
def compute_metrics(tp, fp, truth_total):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / truth_total
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"tp": tp, "fp": fp, "precision": precision, "recall": recall, "f1": f1}
```

---

## 關鍵量定義

| 量 | 定義 | 取得方式 |
|----|------|---------|
| **TP** (True Positive) | caller 輸出 ∩ truth VCF | `filtered_snv_tp.vcf.gz` 變異數 |
| **FP** (False Positive) | caller 輸出 ∖ truth VCF | `filtered_snv_fp.vcf.gz` 變異數 |
| **FN** (False Negative) | truth VCF ∖ caller 輸出 | 需 `bcftools isec` 計算 |
| **TRUTH_TOTAL** | truth VCF 的變異總數 | 見 [01_truth_set_registry.md](01_truth_set_registry.md) |

---

## 三個指標

### Precision
```
Precision = TP / (TP + FP)
```
**意義**：caller 輸出中有多少比例是真正變異

### Recall (Sensitivity)
```
Recall = TP / TRUTH_TOTAL
```
**意義**：truth set 中有多少比例被 caller 抓到

**注意**：分母是 `TRUTH_TOTAL`，**不是** `TP + FN`（若 caller 完全失敗，FN 無從得知）

### F1 Score
```
F1 = 2 × Precision × Recall / (Precision + Recall)
```
**意義**：Precision 與 Recall 的調和平均

---

## Δ F1 (ISM-filtered F1 - baseline F1)

```
Δ F1 = F1_after_ISM - F1_before_ISM
```

**使用**：
- `F1_before_ISM`：ClairS 原始輸出的 F1
- `F1_after_ISM`：ISM filter 後的 F1

**Phase 1A Locked**：
- paired_full ΔF1 = **+0.0112**（locked）
- TO ΔF1 = **-0.0206**（NEGATIVE locked）
- **完整 provenance（CI、方法、樣本、限制）** → [../03_pipelines/05_f1_baseline_canonical.md](../03_pipelines/05_f1_baseline_canonical.md)（SoT）

---

## TP/FP VCF 生成流程

```bash
# scripts/pipeline/steps/02_filter_analysis.py 類
# 對 caller VCF 分類為 TP / FP

# 簡化流程（實際用 bcftools isec）
bcftools isec -p isec/ <truth.vcf.gz> <caller.vcf.gz>
# 0000.vcf.gz = truth ∖ caller   (FN candidates)
# 0001.vcf.gz = caller ∖ truth   (FP candidates)
# 0002.vcf.gz = intersection      (TP candidates)
```

---

## HC BED 過濾

若 truth set 有 HC BED，**必須**先用 HC BED 過濾 caller 輸出：

```bash
bcftools view -R <hc.bed> <caller.vcf.gz> > <caller_hc.vcf.gz>
# 然後再分 TP/FP
```

**為何**：HC BED 外的 caller 輸出不屬於 truth set 定義範圍，不應計入 FP

---

## 計算範例

### 範例 1：HCC1395 paired_full
```
TP = 24,500
FP = 3,200
TRUTH_TOTAL = 39,447

Precision = 24500 / (24500 + 3200) = 0.8845
Recall    = 24500 / 39447 = 0.6211
F1        = 2 × 0.8845 × 0.6211 / (0.8845 + 0.6211) = 0.7297
```

### 範例 2：ISM filter 後
```
ISM filter 去除 500 FP，誤刪 100 TP
TP_new = 24500 - 100 = 24400
FP_new = 3200 - 500 = 2700

Precision_new = 24400 / (24400 + 2700) = 0.9004
Recall_new    = 24400 / 39447 = 0.6186
F1_new        = 0.7334

Δ F1 = 0.7334 - 0.7297 = +0.0037
```

---

## 聚合跨 sample 的 F1

```bash
python3 scripts/analysis/compare_benchmark_f1.py \
  --run-dir output/canonical/ \
  --output-tsv all_samples_f1.tsv
```

**輸出欄位**：
```
dataset_id, truth_label, precision, recall, f1, tp, fp, truth_total
```

---

## 常見陷阱

- ❌ **TRUTH_TOTAL 用錯**：應為 truth VCF 變異數，不是 TP+FN
- ❌ **不用 HC BED 就算 F1**：結果 FP 膨脹
- ❌ **混用不同版本 truth VCF**：F1 不同
- ❌ **TP+FP 不 match**：需 TP、FP、truth 三者使用**同一 truth VCF 定義**

---

## 97.77% PON 純度（7 樣本交叉驗證）

- 2026-04-11 獨立分析結果
- PON 去除後 germline FP 純度 97.77% ± 1.12%（7 樣本）
- 穩定度 3→4（提升）

（與 F1 無直接關係，屬 PON 獨立驗證）

---

## 相關

- Truth set：[01_truth_set_registry.md](01_truth_set_registry.md)
- Benchmark 協議：[03_benchmark_protocols.md](03_benchmark_protocols.md)
- F1 Workflow：[../06_workflows/04_f1_benchmark.md](../06_workflows/04_f1_benchmark.md)
- 原始碼：`scripts/pipeline/steps/03_filter_analysis.py:229-234`
