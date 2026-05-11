---
id: ism-kb-08-truth-and-benchmark-index
name: "Truth & Benchmark 目錄索引"
description: "7 樣本 truth set 清單、F1 計算、som.py / hap.py benchmark 協議。論文數據重現的必備查詢入口。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "truth and benchmark directory structure"
related_ids:
  - ism-kb-08-truth-and-benchmark-truth-set-registry
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-08-truth-and-benchmark-benchmark-protocols
tags: [truth, benchmark, index, f1, som.py]
canonical_paths: [08_truth_and_benchmark/00_index.md]
alias_paths: []
---

# Truth & Benchmark 目錄索引

- 一句結論：查 truth set 路徑用 01，算 F1 看 02，用 som.py / pipeline 跑 benchmark 用 03
- 適用對象：Benchmark 執行、論文數據
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cat /big7_disk/liaoyoyo2001/InterSubMod/docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md | head -40
  ```

---

## 文件列表

| 檔案 | 主題 |
|------|------|
| [01_truth_set_registry.md](01_truth_set_registry.md) | 7 樣本 truth VCF / HC BED / TRUTH_TOTAL 完整表 |
| [02_f1_calculation.md](02_f1_calculation.md) | F1 公式 + 計算邏輯 + pipeline 整合 |
| [03_benchmark_protocols.md](03_benchmark_protocols.md) | som.py (hap.py) vs 內建 F1；何時用何者 |

---

## 權威文件

👉 **`docs/data_specs/20260422_truth_sets_and_f1_protocol_01.md`**

此 KB 提供結構化入口；全部細節在該權威文件。

---

## 相關

- Pipelines：[../03_pipelines/](../03_pipelines/)
- F1 Workflow：[../06_workflows/04_f1_benchmark.md](../06_workflows/04_f1_benchmark.md)
