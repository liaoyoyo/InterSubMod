---
id: ism-kb-07-derived-features-index
name: "Derived Features 索引"
description: "衍生特徵索引：HPFineNGroups、Coverage_Multiple、Zone-Aware、LOH×AF×Methylation、Fisher_Frac_Sig。各特徵的 positive/characterization/NEGATIVE 狀態對照表。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "derived features status table"
related_ids:
  - ism-kb-07-derived-features-hpfinengroups
  - ism-kb-07-derived-features-coverage-multiple
  - ism-kb-07-derived-features-zone-aware-framework
  - ism-kb-07-derived-features-loh-af-methylation
  - ism-kb-09-conclusions-index
  - ism-kb-07-derived-features-fisher-frac-sig
  - ism-kb-09-conclusions-positive-findings
  - ism-kb-06-workflows-analysis-scripts-index
tags: [features, index, derived, status]
canonical_paths: [07_derived_features/00_index.md]
alias_paths: []
---

# Derived Features 索引

- 一句結論：5 個主要衍生特徵，狀態分為 positive（可用）、characterization（能描述但不可做 filter）、NEGATIVE（已證偽）
- 適用對象：選特徵做分析、理解結論狀態
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  ls /big7_disk/liaoyoyo2001/InterSubMod/knowledge/07_derived_features/
  ```

---

## 特徵狀態對照表

| 特徵 | 狀態 | Δ F1 / 效應 | 文件 |
|------|------|------------|------|
| **HPFineNGroups** | 🟡 characterization | N≥4+NR≥80 TP rate 92.81%（**不可 filter**） | [01_hpfinengroups.md](01_hpfinengroups.md) |
| **Coverage_Multiple** | 🟢 CN proxy 可用 | r=0.83 vs CN | [02_coverage_multiple.md](02_coverage_multiple.md) |
| **Zone-Aware Framework** | 🟡 characterization | Zone TP rate 差異真實但 QS 調整 NEGATIVE | [03_zone_aware_framework.md](03_zone_aware_framework.md) |
| **LOH × AF × Methylation** | 🟢 positive | Inter AF→NGroups +0.705 (7/7 p<10^-39) | [04_loh_af_methylation.md](04_loh_af_methylation.md) |
| **Fisher_Frac_Sig** | 🔴 NEGATIVE (F1-filter) | CI 跨隨機；F pilot TP 99.5% 飽和 | [05_fisher_frac_sig.md](05_fisher_frac_sig.md) |

---

## 狀態圖示

- 🟢 **Positive**：可作為分析特徵使用
- 🟡 **Characterization only**：能描述真實現象，**不可做 variant filter**
- 🔴 **NEGATIVE / NO-GO**：已證偽，別重新調查

---

## 選擇建議

| 目的 | 推薦特徵 |
|------|---------|
| 做 variant filter | **❌ 無**（所有 filter 方向皆 NEGATIVE） |
| 做 characterization | HPFineNGroups, Zone-Aware |
| CN / copy number 分析 | Coverage_Multiple |
| LOH × epigenetic 研究 | LOH × AF × Methylation |

---

## 與結論的交叉關係

此目錄的每個特徵文件**只解釋特徵本身**；研究結論（filter 可行性、Δ F1）在 `09_conclusions/`：
- Positive findings：[../09_conclusions/01_positive_findings.md](../09_conclusions/01_positive_findings.md)
- Characterization only：[../09_conclusions/02_characterization_only.md](../09_conclusions/02_characterization_only.md)
- NEGATIVE 目錄：[../09_conclusions/03_concluded_negative.md](../09_conclusions/03_concluded_negative.md)

---

## 相關

- C++ 輸出欄位：[../05_data_formats/01_significance_summary_schema.md](../05_data_formats/01_significance_summary_schema.md)
- Master dataset：[../05_data_formats/02_master_dataset_schema.md](../05_data_formats/02_master_dataset_schema.md)
- Research landscape：[../../docs/reports/research_landscape/00_INDEX.md](../../docs/reports/research_landscape/00_INDEX.md)
