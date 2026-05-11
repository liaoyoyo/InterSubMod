---
id: ism-kb-06-workflows-batch-tp-fp-analysis
name: "Batch TP/FP Analysis"
description: "`run_batch_vcf_analysis.sh` 並行執行 TP 與 FP VCF，生成對比報表；Full test 標準流程。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: howto
verified_scope: "scripts/run_batch_vcf_analysis.sh flow"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-06-workflows-f1-benchmark
tags: [workflow, batch, tp-fp, parallel]
canonical_paths: [06_workflows/03_batch_tp_fp_analysis.md]
alias_paths: []
---

# Batch TP/FP Analysis

- 一句結論：`./scripts/run_batch_vcf_analysis.sh` 同時跑 TP 與 FP VCF，輸出對比；~5 min（HCC1395）
- 適用對象：Full test、F1 benchmark 主流程
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  /big7_disk/liaoyoyo2001/InterSubMod/scripts/run_batch_vcf_analysis.sh --help 2>&1 | head -20
  ```

---

## 腳本概要

- **位置**：`scripts/run_batch_vcf_analysis.sh`
- **功能**：內部呼叫 `run_vcf_all_snv.sh` 兩次（TP + FP），並行化加速
- **典型時間**：~5 分鐘（HCC1395 paired_full）

---

## 典型命令

```bash
./scripts/run_batch_vcf_analysis.sh \
  --caller-mode paired \
  --mode all-with-w5000 \
  --threads 120 \
  --plot-type all
```

**等價於**：
```bash
# 並行執行
./scripts/run_vcf_all_snv.sh ... -v filtered_snv_tp.vcf.gz -o out/TP &
./scripts/run_vcf_all_snv.sh ... -v filtered_snv_fp.vcf.gz -o out/FP &
wait
# 然後生成對比報表
```

---

## 輸出

```
<output>/
├── TP/
│   └── <ism output>
├── FP/
│   └── <ism output>
└── benchmark_comparison.tsv        # ← 主要對比表
```

### benchmark_comparison.tsv 欄位
| 欄位 | 意義 |
|------|------|
| dataset_id | {sample}_{mode}_{timestamp} |
| truth_label | TP / FP |
| num_regions | region 總數 |
| num_reads_mean | 平均每 region read 數 |
| f1_before_ism | ClairS 原始 F1 |
| f1_after_ism | ISM 過濾後 F1 |
| delta_f1 | f1_after_ism - f1_before_ism |

---

## /test-full

此 workflow 對應 `/test-full` slash command：
```bash
/test-full  # = ./scripts/run_batch_vcf_analysis.sh
```

**驗收**：Δ F1 與上次 canonical 差異 < 0.01

---

## 7 樣本批次（手動迴圈）

若要跑全 7 樣本：
```bash
for sample in HCC1395 HCC1395_DORADO COLO829 H1437 H2009 HCC1937 HCC1954; do
  ./scripts/run_batch_vcf_analysis.sh \
    --caller-mode paired --mode all-with-w5000 \
    --sample $sample \
    -o output/batch_${sample}
done
```

**時間**：~35-40 分鐘（7 樣本連跑）

---

## 平行化 vs 記憶體

| Threads | 預期加速 | 記憶體峰值 |
|---------|---------|-----------|
| -j 60 | ~6x | ~60 GB |
| -j 120 | ~10-12x | ~120 GB |
| -j 240 | ~12-15x | OOM risk |

**建議**：server 128 GB RAM → `-j 120`

---

## 常見誤區

- ❌ **跑 TP 未跑 FP 就計算 F1**：需要兩者
- ❌ **TP 與 FP 用不同參數**：必須完全一致的 ISM 參數
- ❌ **忽略 benchmark_comparison.tsv 報表**：那是主要 F1 對比來源

---

## 相關

- Full VCF 分析：[02_full_vcf_analysis.md](02_full_vcf_analysis.md)
- F1 benchmark 進階：[04_f1_benchmark.md](04_f1_benchmark.md)
- 原始腳本：`scripts/run_batch_vcf_analysis.sh`
