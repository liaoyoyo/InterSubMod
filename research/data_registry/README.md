---
title: KDE-Corrected Data Provenance Registry
created: 2026-04-26
updated: 2026-04-26
status: active
owner: InterSubMod Research
---

# KDE-Corrected Data Provenance Registry

## 目的

InterSubMod 過去（2026-03-30 之前）的 master dataset 含 **stale-binary KDE artifact**：`Diploid_Coverage_Used = 75.0×` 為硬編碼預設值（並非每樣本實測 KDE 估計）。新 KDE binary（commit `374fad4` / `12d9b3e`）已修正 KDE fallback WARN 與 audit column 機制。

本 registry 標明哪些 dataset 為 KDE-corrected、哪些仍為 stale，以避免下游分析誤用污染資料。

## 檔案

- `kde_corrected_provenance_20260426.tsv` — 主 registry（tab-separated，9 欄）

## 欄位定義

| 欄位 | 說明 |
|------|------|
| `dataset_path` | dataset 絕對路徑（路徑為目錄或 .tsv.gz 檔） |
| `sample` | 樣本名（HCC1395 / HCC1395_DORADO / HCC1937 / HCC1954 / H2009 / H1437 / COLO829 / MIXED_*） |
| `mode` | 分析模式（paired_full / paired_pileup / to_pileup / legacy_master_*） |
| `build_commit` | 產生此 dataset 的 binary git commit（hash 短碼） |
| `Diploid_Coverage_Used` | KDE 估計的 diploid coverage 中位數（單位：×）。`75.0` = stale-binary fallback |
| `is_stale` | `True` = 含 stale-KDE artifact 或未通過 KDE-fix 重跑驗證；`False` = 已驗證 KDE-corrected |
| `last_verified` | 最近一次驗證日期（YYYY-MM-DD） |
| `do_not_use` | `True` = 嚴禁用於新分析；`False` = 可使用；`with_caveat` = 可參考但需在報告註明 KDE 狀態尚未獨立 audit |
| `downstream_dependencies` | 受此 dataset 影響的下游分析（分號分隔識別字） |

## KDE Fix 提交歷史

| Commit | 內容 |
|--------|------|
| `ec0608b` | snapshot baseline before KDE expected_coverage logging+audit change（**stale baseline**） |
| `374fad4` | feat(covm): KDE fallback WARN + diploid_coverage_used audit column（**KDE-corrected baseline**） |
| `5abc659` | chore: evidence_ledger record H_KDE_001 cpp_improvement KDE logging+audit |
| `485075f` | fix(scripts): point run_vcf_all_snv.sh EXECUTABLE to big7 KDE-fixed binary |

## 識別 stale 資料的方式

1. **Diploid_Coverage_Used 欄位空值** → significance_summary.csv 來自 `ec0608b` 之前；無 KDE audit
2. **Diploid_Coverage_Used = 75.0** → stale-binary fallback；資料未經實測 KDE
3. **Diploid_Coverage_Used 為樣本特異實值**（例 HCC1954=109.0, COLO829=99.0, HCC1395 TO=115.0）→ KDE-corrected

## 如何更新 registry

當 archive rerun 完成後：

1. 將該 row 的 `is_stale` 由 `True` 改 `False`
2. 將 `do_not_use` 由 `with_caveat` 改 `False`
3. 更新 `last_verified` 為當日日期
4. 在 `evidence_ledger.jsonl` 寫入對應 audit 紀錄（type=`kde_provenance_audit`）
5. 同步 `research/feature_layered_observation/data/_manifest.md` 的 KDE 欄位

## 與 evidence_ledger.jsonl 的關係

- **registry**（本檔）= **dataset-level** provenance：哪個目錄哪個 commit 哪個 KDE 值
- **evidence_ledger.jsonl** = **evidence-level** trace：每個假說/實驗的 input dataset 對應到此 registry 的 path

evidence_ledger 引用本 registry：寫 evidence 時 `input_data: <dataset_path from registry>`，可透過 path 反查 registry 的 `is_stale` 與 `do_not_use`。

## Thread D 主軸報告連結

- Thread D 主軸 Validated 報告 §10 證據檔案路徑：**TBD**（待 P0-1 撰寫完成後補連結至 `docs/reports/thread_d/10_evidence_files.md`）
- HCC1954 Standalone Case Panel：`docs/reports/case_studies/HCC1954/`（P0-3）
- KDE Fix Acceptance Validation：`docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md`

## 已知限制

1. **B 段 5 row（HCC1937, HCC1954, H2009, H1437, HCC1395_DORADO TO）的 `is_stale=True`** 是因 archive rerun（位於 `kde_smoke_test/x1_archive_to_rerun/`，2026-04-24 完成）尚未通過獨立 audit；資料本身可能已 KDE-corrected 但未進 master 流程。
2. **C 段 master TSV** 只列為單一 row 代表 legacy；實際 stale 區塊為 7 樣本 × 多模式 mixed，無法以單一 sample/mode 標識。
3. **paired_pileup 模式未列入** —— 該模式為次要驗證軌（per docs/CURRENT_FOCUS.md），暫不需 KDE provenance 標記；如需請於下一輪擴展。
