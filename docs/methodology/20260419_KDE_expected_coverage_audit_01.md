---
title: KDE expected_coverage 方法學審查
date: 2026-04-19
type: methodology_audit
status: pending_user_decision
related_commits:
  - 8d0a0c8 (2026-04-13): feat — Normal BAM + Phase A-D 引入 KDE 二次 pass
related_reports:
  - docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md
  - research/coverage_multiple_validation/reports/20260419_covm_vs_seqc2_truth_step1_2_01.md
---

# KDE expected_coverage 方法學審查

## 1. 審查緣由

CovM validation 報告（2026-04-19）發現 7 樣本 × 2 模式的 `expected_coverage` 皆等於 **硬編碼值 75.0**，per-sample bias −28% 到 +150%。結論：CovM 跨樣本不可比，Z3 amplicon blacklist pilot 無法泛化。

本審查目的：定位 **「75.0 hardcoded」是 C++ 缺陷還是操作層 artifact**，並決定是否需要修改 C++。

---

## 2. 關鍵時間線

| 日期時間 | 事件 | 證據 |
|---------|------|------|
| 2026-03-30 17:43 | Master dataset `all_region_rows.tsv.gz` 產生 | `ls -la` 時間戳 |
| **2026-04-13 03:21** | **commit 8d0a0c8 首次引入 KDE 二次 pass** | `git log src/core/RegionProcessor.cpp` |
| 2026-04-19 | CovM validation 檢出 75.0 hardcoded | Step 2 報告 |

**結論**：Master dataset 比 KDE 實作 **早 14 天**。報告所用數據並非由現行 C++ 產生，而是由 **不含 KDE 的舊 binary** 產生 → 當時所有 region 均走 `compute_coverage_multiple(num_reads, 75.0)` 的 default 路徑。

---

## 3. 現行 C++ 實作（無缺陷）

### 3.1 Two-pass 架構

- **Pass 1（parallel）** `src/core/RegionProcessor.cpp:1654` — 每個 region 以 `expected_coverage=75.0` default 預估 CovM
- **Pass 2（serial post-processing）** `RegionProcessor.cpp:662-684`（新增於 8d0a0c8）：
  ```cpp
  double diploid_coverage;
  if (expected_coverage_ > 0.0) {
      diploid_coverage = expected_coverage_;            // 使用者指定
  } else {
      diploid_coverage = estimate_diploid_coverage(results);  // auto-KDE mode
  }
  for (auto& r : results) {
      if (r.success) {
          r.coverage_multiple = compute_coverage_multiple(r.num_reads, diploid_coverage);
          r.coverage_category = determine_coverage_category(r.coverage_multiple);
          ...
      }
  }
  ```

### 3.2 KDE 估計方法（`RegionProcessor.cpp:74-149`）

1. 收集 `num_reads > 5` 的全部 region
2. 以 p80 上限排除極端 gain 影響
3. 直方圖（bin=2）+ Gaussian 平滑（σ=3 bins）
4. 取峰值為 diploid mode
5. Fallback 回 75.0 的三種情況：
   - 有效 region < 100（資料不足）
   - n_bins < 10（範圍太窄）
   - estimate 落在 10-300 區間之外（數值異常）

### 3.3 Config default

`include/core/Config.hpp:56` — `expected_coverage = 0.0` 表示走 auto-KDE 路徑。Scripts 亦無傳入 `--expected-coverage`，驗證 auto-KDE 為現行預設行為。

### 3.4 審查判定

**C++ 實作無缺陷**。報告觀察到的 75.0 hardcoded **純由操作層原因（舊 binary）產生**。

---

## 4. 次要觀察（KDE 健壯性風險）

即使 C++ 無 bug，KDE 實作仍有兩個值得改善的點：

### 4.1 Fallback 靜默

三種 fallback 回 75.0 的情況僅以 `LOG_INFO` 記錄，在大量日誌中難以察覺。若未來某樣本資料量恰低於 100 regions 或 num_reads 分布極端，使用者**不會被明確提醒** CovM 不可信。

### 4.2 高 ploidy 樣本的 mode 上抬

KDE mode 假設二倍體區域在 num_reads 分布上形成最高峰。當樣本含大量 gain 區域（如 HCC1954 HER2 amplicon），p80 切截仍可能把 mode 拉到超二倍體區。現行 10-300 sanity 過寬，無法捕捉此類 bias。

### 4.3 下游審計不便

TSV 輸出**不含** `diploid_coverage_used` 欄位，下游分析無法判斷該 region 的 CovM 是以哪個 baseline 計算。若要 audit 跨樣本 CovM scale 是否可比，必須重新 parse log。

---

## 5. 影響範圍（如不修）

| 受影響分析 | 是否需重跑 Master | 原因 |
|-----------|------------------|------|
| CovM validation (H1/H2) | **必須** | 結論核心基於 per-sample expected_coverage 估計 |
| Z3 amplicon blacklist pilot | 必須（若 CovM 分層為關鍵變數） | CovM 數值跨樣本不可比 |
| HPFineNGroups canonical filter | 建議（CovM 未作為直接分層，影響較低） | HPFineNGroups 為主要變數 |
| Round-level LOH enrichment | 建議（不直接用 CovM） | CovM 僅為 QC 輔助 |

**重點**：Master dataset 需以現行 binary 重跑一次，後續所有 CovM 相關結論才能以 per-sample KDE baseline 為基準。

---

## 6. 方案選項

### A) 不做任何事 ❌

已被 cpp-change skill 前提排除，此處僅陳述原因：**若不更新 master dataset 並修正結論，CovM validation H1 negative 的核心論據不成立**（75.0 hardcoded 本身不是 C++ bug，是 stale artifact）。

### B) 最小干預（僅結論修正 + 重跑） ✅ 可選

- **不改 C++**，只改文件
- **Deliverables**：
  1. 更新 `20260420_CovM_Baseline_Accuracy_Validation_01.md`：
     - H-CN1 由 POSITIVE 改為 **CONDITIONAL**（需以新 master dataset 重測）
     - 標註此報告依據 stale master（2026-03-30）
  2. 更新 `research/coverage_multiple_validation/00_PLAN.md` 與 manifest：
     - 明記數據源時間戳問題，Step 1/2 需在新 master 上重做
  3. 於 `docs/CURRENT_FOCUS.md` 新增一項 blocker：master dataset 重跑
  4. 新增 commit `chore:` 記錄 KDE 驗證結論
- **優點**：零 C++ 風險、可立即完成
- **缺點**：未預防 KDE 下次 silent fallback；未增加下游審計便利性
- **預估工時**：30-60 min（純文件）

### C) 最小 C++ 改善 + 結論修正 + 重跑 ✅ 可選

- **B 的全部** + 下列 C++ 改動：
  1. **C1（Logging 提級）** `RegionProcessor.cpp:74-149`
     - 三個 fallback return 75.0 處改為 `LOG_WARN` 並印出原因（資料不足／範圍太窄／estimate 異常）
     - Step 2 pass 開始時印出 `[KDE] Using diploid_coverage=<value> (source=auto_kde|user_specified|fallback_default)`
  2. **C2（Audit column）** `include/core/RegionResult.hpp` + `src/io/OutputWriter.cpp`
     - 新增 `double diploid_coverage_used` 欄位
     - TSV 輸出多一欄 `diploid_coverage_used`（float, 2 decimals）
  3.（選做）**C3（KDE 健壯性）** 可**暫緩**，因現行 mode 估計對 bulk tumor 通常可接受；若 C2 上線後觀察到高 ploidy 樣本 mode 偏高，再專開一輪 audit
- **優點**：下次 rerun 可直接從 TSV audit；未來 stale/配置問題更易偵測
- **缺點**：需 1.5-2h 實作 + 測試 + rerun 配合
- **預估工時**：
  - C1 + C2 程式碼：~30-45 min
  - Unit test：~30 min
  - Full rebuild + quick test：~30 min
  - 文件 + commit：~30 min
  - Total：2-2.5h
- **風險**：C2 變更 output schema，下游 pandas 腳本若指定欄位清單需調整

### D) 進一步改 KDE 演算法 ❌ 暫不建議

- 改用 trimmed median / 2nd-derivative peak / per-chromosome weighted mode 等方法
- 風險：目前無證據顯示 auto-KDE 在常規樣本失準；改演算法前應先讓 C2 上線、取得 audit 資料再決定
- **建議歸入獨立後續研究，而非本次審查範圍**

---

## 7. 建議與決策請求

**個人建議方案 C**：
- C1/C2 改動範圍小、風險低
- 修完後的 rerun 可一次完成「驗證 KDE 正確 + 產生可 audit 的新 master + 修正結論」
- 比 B 多花 1-1.5h，但可一勞永逸建立 KDE 透明度

**若時間緊迫則選 B**：
- 接受當前 audit 不便
- 僅修文件 + rerun，結論層面影響與 C 相同

### ⚠ Hard Gate — 需用戶選擇

請選擇：
- **B**：純文件修正 + 重跑（不動 C++）
- **C**：C++ 小幅改善（C1+C2）+ 文件修正 + 重跑
- **其他**（D 或其他方向，請說明）

選定後進入 `/cpp-change` Step 1（baseline commit）。
