---
id: ism-kb-02-samples-index
name: "Samples 目錄索引"
description: "7 主樣本快速對照：HCC1395/HCC1395_DORADO/COLO829/H1437/H2009/HCC1937/HCC1954；各自 truth set、canonical run、特殊 caveat。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "7 samples directory index"
related_ids:
  - ism-kb-02-samples-hcc1395
  - ism-kb-08-truth-and-benchmark-truth-set-registry
  - ism-kb-02-samples-hcc1395-dorado
  - ism-kb-02-samples-colo829
  - ism-kb-02-samples-h1437
  - ism-kb-02-samples-h2009
  - ism-kb-02-samples-hcc1937
  - ism-kb-02-samples-hcc1954
  - ism-kb-06-workflows-adding-new-sample
tags: [samples, index, overview]
canonical_paths: [02_samples/00_index.md]
alias_paths: []
---

# Samples 目錄索引

- 一句結論：7 主樣本中 HCC1395 為 benchmark 主 canonical；每樣本有 paired_full + paired_pileup + TO 三 pipeline
- 適用對象：選擇測試樣本、理解 sample-specific caveats
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/big7_disk_output/canonical/
  ```

---

## 7 樣本速查表

| 樣本 | Platform | Truth Source | TRUTH_TOTAL | 特殊情況 | 文件 |
|------|---------|-------------|-------------|---------|------|
| **HCC1395** | ONT 5kHz | SEQC2 | 39,447 | 主 benchmark；最完整 | [01_hcc1395.md](01_hcc1395.md) |
| **HCC1395_DORADO** | ONT Dorado | SEQC2（共用） | 39,447 | Basecall 變體 | [02_hcc1395_dorado.md](02_hcc1395_dorado.md) |
| **COLO829** | ONT PAO/R10 | NYGC | 41,427 | 無 HC BED；TO 無 methylation | [03_colo829.md](03_colo829.md) |
| **H1437** | ONT | Orthogonal tools | 90,129 | 標準 | [04_h1437.md](04_h1437.md) |
| **H2009** | ONT | Orthogonal tools | 168,529 | 負向案例（caller 已完美） | [05_h2009.md](05_h2009.md) |
| **HCC1937** | ONT | Orthogonal tools | 60,691 | 標準 | [06_hcc1937.md](06_hcc1937.md) |
| **HCC1954** | ONT | Orthogonal tools | 26,567 | Amplicon artifact | [07_hcc1954.md](07_hcc1954.md) |

---

## 選擇樣本建議

| 目的 | 推薦樣本 |
|------|---------|
| 新功能開發 | HCC1395（資料最全） |
| Self-phasing 研究 | HCC1395, COLO829 |
| LOH × AF 研究 | 全 7 樣本（已驗證 7/7 一致） |
| Amplicon 處理 | HCC1954（已知 artifact） |
| 負向控制 | H2009（caller 已完美，ISM 難提升） |

---

## 相關

- Truth set 全表：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)
- 新增樣本流程：[../06_workflows/06_adding_new_sample.md](../06_workflows/06_adding_new_sample.md)
