---
title: ISM 參數實驗框架 — 位點排除記錄 + 參數改動→結果差異追溯
date: 2026-06-12
type: governance / framework-design
status: draft（與用戶確認中）
branch: chore/ism-review-param-governance-202606
related:
  - 實例 #1: docs/experiments/in_progress/2026/06/20260612_2201_sampling_threshold_audit/manifest.yaml
  - 母設計: docs/plans/20260612_ism_review_param_governance_design_01.md
---

# ISM 參數實驗框架（用戶 #4 + #5）

> 目標：讓**未來任何參數想法/自動測試**都能 (a) 知道每個位點被保留/排除的原因與數值；(b) 量化「某個改動讓多少位點的結果改變、因為哪個門檻」。

## L0 兩個框架

| 框架 | 解決的痛點 | 需要 |
|---|---|---|
| **A 位點排除原因記錄** | 現在被前置 gate 擋的位點（如 **4215 個 FP, 87%**）`significance_computed=false` 就跳過寫出 → 完全不可見、無原因、無數值 | C++ 改動（小）|
| **B 參數改動→結果差異追溯** | 改一個參數重跑後，無法系統化回答「多少位點歸類變了、因為哪個門檻」 | manifest 規範 + diff 工具 |

---

## A. 位點排除原因記錄框架

### A.1 現狀（root cause）
`RegionProcessor.cpp:1134` — `write_significance_summary` 跳過 `!significance_computed` 的位點。被位點前置 gate（`reads < clustering_min_reads(10)` 或 `cpgs < 1` 或 SNV 無覆蓋）擋掉的位點**根本不寫任何檔**。
→ summary 只含「算過的」；被擋的（TP 736 / **FP 4215** / FN ?）連 RegionID 都沒有。

### A.2 設計：全位點 accounting 表（新增輸出）
讓**每個進入 pipeline 的位點**（含被擋的）都寫一筆到 `region_accounting.csv`：

| 欄 | 內容 |
|---|---|
| RegionID, Chr, Pos, Ref, Alt | 位點識別 |
| Class | TP / FP / FN（由輸入 VCF 決定）|
| NumReads, NumCpGs | 進入時的數值（即使被擋也記）|
| Stage | `PRE_GATE` / `ANALYZED` |
| Outcome | `RETAINED` / 各 `REJECT_*` |
| RejectReason | 枚舉（見下）|

**RejectReason 枚舉**（對齊已盤點的門檻層）：
```
RETAINED                       # 通過全部 gate（= 現在的 Significant=true）
REJECT_PREGATE_READS           # reads < clustering_min_reads(10)
REJECT_PREGATE_CPGS            # num_cpgs < 1
REJECT_PREGATE_NO_SNV_COVERAGE # tumor read 不覆蓋 SNV
REJECT_GATING                  # 未過 Fisher×CramersV gating
REJECT_GLOBALP                 # GlobalP > 0.05
REJECT_CRAMERSV                # CramersV(gated) < 0.1
REJECT_READS20                 # NumReads < 20（保留 depth gate）
```

### A.3 落地（走 methodology-audit → cpp-change）
- 改動點：`RegionProcessor.cpp` 在前置 gate 判斷處（:833 / :857）+ summary 寫出處（:1134），無論是否 computed 都 append 一筆 accounting。
- 影響：純新增輸出欄/檔，不改任何判定邏輯 → 低風險。
- 收益：`groupby(RejectReason)` 即得各層排除統計；**未來改任一門檻重跑，accounting 的 RejectReason 分佈 diff = 該門檻影響的位點數**。

> ⚠ 本輪（這次全量觀察）尚無此欄位 → 各層排除數用「VCF 總 − summary computed」推算總數（TP 736/FP 4215），**精確分因（reads vs cpgs vs no-coverage）需 A.2 落地後才有**。

---

## B. 參數改動→結果差異追溯框架

### B.1 一個實驗 = 一份 manifest（已 dogfood 實例 #1）
每次「為某參數跑 ISM」建一個帶時間戳目錄 + `manifest.yaml`（schema 見實例 #1）：
- `parameters_under_test`：這次的參數值（baseline）
- `provenance`：git commit + binary commit（鎖版本，對齊 web 最佳實踐）
- `exclusion_accounting` / `results`：結果數字（含 A 的 RejectReason 分佈）
- `artifacts`：產物分類（baseline 留 / ephemeral 刪）
- `cost`：時間/磁碟（可預估）

### B.2 改參數 → 對比
```
實驗A (C_min=3) ──┐
                  ├─ scripts/param_experiment_diff.py manifestA manifestB
實驗B (C_min=5) ──┘   → 輸出：哪些位點 RejectReason 變了、各門檻 ±幾個位點
```
diff 工具比較兩 manifest 的 `exclusion_accounting` + `results` + （若有）逐位點 accounting → 量化改動影響。

### B.3 產物分類規則（你的 #3，寫成規範）
| 類別 | 範例 | 處置 |
|---|---|---|
| **baseline**（重複比對主要）| 彙總統計表、golden 輸出 | 入 git，長期保留 |
| **observation**（重複觀察）| 圖表 | 入 git |
| **record** | manifest / README | 入 git |
| **ephemeral**（知道結果即移除）| per-region 矩陣（~12GB）| 跑完分析即刪，**不入 git、不長期佔磁碟** |

### B.4 成本記錄（可預估）
manifest `cost` 欄記 `expected_runtime` + `actual_runtime` + `peak_disk` → 累積後可預估「全量 ISM ≈ N min / M GB」。本實例基準：全量 ISM ≈ 3-5min（純彙總）/ 開矩陣 ~12GB 臨時。

---

## C. 與階段 A 機械守衛（母設計 §4.A）的關係
- A 排除記錄 + B 追溯 = 階段 A「參數治理」的**資料層**。
- 階段 A 的**機械守衛層**（hook）在此之上：param 改動納入 invalidation ledger、結果按 param 隔離、binary 內嵌 commit hash。
- 落地順序：先這份框架（資料層，本 session）→ 再機械守衛（hook，後續）。
