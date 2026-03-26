<!--
建立時間: 2026-03-26
目標: Phase 1 重構後驗證報告（安全性與正確性修正）
處理範圍: Tasks 1.1–1.5 完成後的驗證
關聯檔案:
  - docs/refactor_baseline/BASELINE_SUMMARY.md
  - docs/refactor_baseline/post_phase1/test_results.txt
  - docs/refactor_baseline/post_phase1/code_metrics.txt
  - docs/refactor_baseline/post_phase1/build_time.txt
-->

# Phase 1 重構後驗證報告

**分支**：`refactor/phase1-safety`
**完成時間**：2026-03-26
**Commits（本 branch）**：
- `a657462` — chore: establish refactoring baseline checkpoints
- `4f57ba1` — feat(BamReader): add RAII fetch_reads_safe() to prevent HTSlib memory leaks
- `291ed2a` — fix(Config): eliminate HTSlib resource leaks in validate()
- `5b7efc8` — fix(FisherExact): use log-sum-exp for numerically stable p-value accumulation
- `f9fee1d` — fix(ReadParser): guard against uint32_t→int32_t overflow in SNV position cast

---

## 測試結果

| 指標 | Baseline | Phase 1 後 | 差異 |
|------|---------|-----------|------|
| 測試案例數 | 107 | **121** | +14 新增 |
| 通過率 | 107/107 (100%) | **121/121 (100%)** | 持平 |
| 測試執行時間 | 2.83s | 3.62s | +0.79s（新增測試所致）|

**新增測試明細**：
- `BamRecordPtrTest::TypeExists` — RAII 型別存在確認
- `BamRecordPtrTest::DeleterReleasesMemory` — BamRecordDeleter 正確呼叫 bam_destroy1
- `BamReaderTest::FetchReadsSafe_ReturnsUniquePtr` — fetch_reads_safe() 回傳 BamRecordPtr
- `BamReaderTest::FetchReadsSafe_InvalidChromosome` — 無效染色體回傳空向量
- `ConfigTest::ValidationRepeatedCallsDoNotLeak` — 50 次連呼 validate() 不造成 fd 耗盡
- `MathUtilsTest::LogSumExp_SingleValue` — 單值 log-sum-exp
- `MathUtilsTest::LogSumExp_TwoEqualValues` — 兩個相等值的 log-sum-exp
- `MathUtilsTest::LogSumExp_NumericallyStable_LargeValues` — 大數值穩定性
- `MathUtilsTest::LogSumExp_NumericallyStable_VeryNegativeValues` — 極負值穩定性
- `MathUtilsTest::LogSumExp_Empty` — 空向量邊界
- `FisherExactTest::Fisher2x2_ExtremeSeparation_PValueNotZero` — 極端分離表格 p > 0
- `FisherExactTest::Fisher2x2_StrongAssociation_PValuePrecision` — 強關聯精度
- `ReadParserPositionSafetyTest::ZeroPosition_ReturnsUnknown` — pos=0 邊界防護
- `ReadParserPositionSafetyTest::OverflowPosition_ReturnsUnknown` — pos>INT32_MAX 防護

---

## 編譯時間

| 指標 | Baseline | Phase 1 後 | 差異 |
|------|---------|-----------|------|
| 增量編譯 | ~60s | ~61s | 持平（新增程式碼微幅增加）|

---

## 程式碼規模

| 模組 | Baseline (LoC) | Phase 1 後 (LoC) | 差異 |
|------|---------------|-----------------|------|
| src/core 合計 | 6,693 | ~6,740 | +47（防護/RAII 邏輯）|
| BamReader.cpp | 110 | 126 | +16（fetch_reads_safe）|
| Config.cpp | 120 | 142 | +22（RAII lambda）|
| FisherExact.cpp | 395 | 418 | +23（log-sum-exp）|
| ReadParser.cpp | 234 | 240 | +6（overflow guard）|
| MathUtils.cpp | 240 | 261 | +21（log_sum_exp 實作）|
| tests 合計 | — | 2,171 | 新增 14 測試 |

---

## 修正內容摘要

### Task 1.1：BamReader RAII（commit 4f57ba1）

**問題**：`fetch_reads()` 返回原始 `bam1_t*`，例外路徑中靜默洩漏。

**修正**：
- 新增 `BamRecordDeleter`（自訂刪除器）與 `BamRecordPtr`（unique_ptr 別名）
- 新增 `fetch_reads_safe()` — 返回 `vector<BamRecordPtr>`，舊介面保留向後相容

**驗證**：編譯時靜態檢查（static_assert）+ 執行時刪除器呼叫確認

---

### Task 1.2：Config::validate() 資源洩漏（commit 291ed2a）

**問題**：HTSlib 檔案控制碼（samFile、hts_idx_t、bcf_hdr_t）在錯誤路徑可能未釋放。

**修正**：全部改為 lambda 自訂刪除器 + `unique_ptr` RAII 管理，所有退出路徑自動釋放。

**驗證**：50 次連呼 validate()（無效路徑），無崩潰或 fd 耗盡。

---

### Task 1.3：FisherExact log-sum-exp 數值穩定性（commit 5b7efc8）

**問題**：直接累加 `exp(log_p_k)` 在極端 2×2 表格下精度損失。

**修正**：
- `MathUtils::log_sum_exp()` — max-shift trick，避免 overflow/underflow
- `FisherExact::test_2x2()` — 收集尾部 log 機率後一次呼叫 log-sum-exp

**驗證**：對比 R `fisher.test()` 的極端表格結果（例如 [[100,2],[3,95]]，期望 p < 1e-10）

---

### Task 1.4：PERMANOVA 災難性消去（暫緩）

> 等 Phase 1 ML read classification 研究結論確定後再執行，避免 p-value 改變影響當前研究基準。

---

### Task 1.5：ReadParser 型別安全（commit f9fee1d）

**問題**：`static_cast<int32_t>(snv.pos) - 1` 在 `snv.pos > INT32_MAX` 時為未定義行為（有號整數溢出）。

**修正**：在轉型前加入明確邊界檢查（pos==0 與 pos>INT32_MAX），兩者皆返回 `UNKNOWN/SNV_NOT_COVERED`。

**驗證**：2 個合成 bam1_t 測試，不需要真實 BAM 資料。

---

## Phase 1 完成結論

✅ **所有 5 個修正目標中 4 個已完成**（Task 1.4 暫緩屬計劃內）
✅ **全部 121 測試通過（100%）**
✅ **無功能退化**（新增測試覆蓋所有修改路徑）
✅ **資源安全性提升**：HTSlib 洩漏已消除（BamReader + Config）
✅ **數值精度提升**：FisherExact p-value 對極端表格更精確
✅ **型別安全**：ReadParser 位置轉型 UB 已修正

**後續**：待使用者確認，可進行 Phase 2 架構改善（分解 God Class、消除重複程式碼）。
