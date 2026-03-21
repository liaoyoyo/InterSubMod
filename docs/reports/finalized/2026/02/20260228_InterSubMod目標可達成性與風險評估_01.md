<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# InterSubMod 目標可達成性與風險評估

- 文件日期：2026-02-28
- 目的：判斷「是否能達成專案最初目標」、缺什麼、風險在哪裡、決策門檻是什麼

## 1. 專案最初目標（對齊來源）

1. 根據 `README_PROJECT_SUMMARY.md`
- 整合 long-read + somatic SNV + methylation，解析腫瘤亞克隆與表觀遺傳異質性。

2. 根據 `Knowledge/05_tools/InterSubMod.md`
- 提供 read-level methylation pattern 分析、距離矩陣、分群與 subclonal resolution。

3. 根據 `Knowledge/06_workflows/benchmark_workflow.md`
- 以可重現 benchmark（TP/FP/FN、HC regions、固定口徑）建立可比較成效。

## 2. 目前達成度評估

| 目標 | 當前狀態 | 達成度 | 主要證據 |
|---|---|---:|---|
| G1：資料讀取與甲基化解析穩定 | 已達成（工程層） | 85% | `test_phase1_2` 可穩定解析 MM/ML；`run_tests` 100 pass |
| G2：以甲基化訊號提升 TP/FP 區分 | 部分達成 | 55% | `HPMergedDelta` 類規則可改善；現行 `Significant` F1 僅 0.0904 |
| G3：亞克隆辨識（Subclone）有效 | 未達成 | 30% | `Subclone` 顯著率 1.48%，召回不足 |
| G4：跨樣本/跨純度可泛化 | 部分達成（文件層） | 35% | 已有方向與部分分析報告，缺統一重跑與同口徑證據鏈 |
| G5：文件與流程可持續維護 | 部分達成 | 60% | 文檔架構完整；但主 README 有舊命令/舊輸出結構殘留 |

## 3. 是否可達成原始目標？

### 結論（2026-02-28）

**可以達成，但有條件。**

必須先完成三個必要條件：
1. 以同口徑重跑跨樣本/跨純度評估，形成可比較基準。
2. 把目前 rule-based 策略升級為「多特徵融合」策略（顯著性只做 feature，不再單獨硬切）。
3. 補齊 purity-aware 且含 MM/ML 的可驗證資料，否則核心假設無法閉環。

若上述三點未完成，專案可達到「工程可用」，但難以達到「研究結論可信且可外推」。

## 4. 風險矩陣

| 風險ID | 風險描述 | 機率 | 影響 | 目前訊號 | 建議處置 |
|---|---|---|---|---|---|
| R1 | 顯著性規則過嚴，長期低召回 | 高 | 高 | Current F1=0.0904 | 改為多特徵融合 + 閾值校正 |
| R2 | `Strong` 類別混入高比例非目標樣本 | 中高 | 高 | FP Strong 量仍高 | 加入 region/CNV/coverage 註解再分類 |
| R3 | purity 分析資料不含 MM/ML | 高 | 高 | Knowledge 已註明限制 | 建立新 purity methylation dataset |
| R4 | 跨樣本泛化不足 | 中高 | 高 | 多數結果集中 HCC1395 | 制定固定 benchmark protocol 並重跑 |
| R5 | 文件與實作脫節導致重現成本高 | 中 | 中 | README 舊命令 | 文件同步與自動檢查 |
| R6 | 測試有長期 skip 導致隱性退化 | 中 | 中 | BamReader 5 tests skip | 補 `data/bam/test.bam` 或替代 fixture |

## 5. 決策 Gate（是否進入下一階段）

### Gate A：研究可驗證性
1. 至少 3 個樣本、同口徑評估（固定 BED、PASS、SNV-only）
2. 每樣本都輸出 TP/FP/FN + methylation 特徵表

### Gate B：方法有效性
1. 新策略在至少 2/3 樣本中優於目前正式基準
2. `Subclone` 召回有統計上可見提升（且 precision 不崩）

### Gate C：可維運性
1. README/Quickstart/主腳本一致
2. CI 或固定測試腳本不再有長期未知 skip

## 6. 建議的目標重定義（避免目標過大過散）

### 6.1 近期（4-6 週）
- 目標：把「甲基化訊號」從探索性結果變成穩定可重現的特徵增益。

### 6.2 中期（6-12 週）
- 目標：完成跨樣本與 purity-aware 的泛化證據，確認方法是否具可轉移性。

### 6.3 長期
- 目標：形成可發表/可共享的 benchmark protocol 與可重現資料包。

## 7. 本文件依據

### 內部證據
- `output/20260118_vcf_all_w5000_t120/significance_summary.csv`
- `output/f1_evaluation_20260119/criteria_comparison.csv`
- `docs/reports/validated/2026/02/20260216_跨樣本跨純度TP_FP_F1綜合分析報告_01.md`
- `docs/reports/validated/2026/02/20260217_purity_aware_analysis_report_01.md`

### 知識庫（根據）
- `Knowledge/05_tools/InterSubMod.md`
- `Knowledge/06_workflows/methylation_analysis.md`
- `Knowledge/04_databases/seqc2_truth_set.md`
- `Knowledge/06_workflows/benchmark_workflow.md`

### 外部對照
- Nature Communications (ClairS-TO): https://www.nature.com/articles/s41467-025-64547-z
- Nature Biotechnology (DeepSomatic): https://www.nature.com/articles/s41587-025-02839-x
- bioRxiv (LongPhase-S): https://www.biorxiv.org/content/10.1101/2025.11.20.689492v1
- BMC Genomics (long-read benchmark): https://link.springer.com/article/10.1186/s12864-025-12259-5
