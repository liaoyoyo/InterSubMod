<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# 20251223_Significance_Implementation_Roadmap

本文件詳述 **Significance Analysis Strategy** 的實作開發計畫。

**目標**: 在 C++ Core 與 Python Tools 中實作全套統計驗證流程，產出具備生物學信心的 Somatic 支持度分數。

---

## 總覽 (Overview)

開發流程分為五個階段，建議依序執行：
1.  **Phase 1: 基礎統計架構 (C++ Core Stats)** - 重點在數值運算基礎與 R 兼容性。
2.  **Phase 2: 全域與局部檢定 (Global & Local Tests)** - 實作核心業務邏輯與 Gating 規則。
3.  **Phase 3: 進階結構驗證 (Structural Validation)** - PERMANOVA 與 Bootstrap。
4.  **Phase 4: 系統整合與輸出 (Integration)** - 串接 RegionProcessor 與 CSV 輸出。
5.  **Phase 5: 文檔與校準 (Doc & Calibration)** - Python 訓練與最終報告。

---

## Phase 1: 基礎統計函式庫 (Basic Statistics Library)

**目標**: 建立穩健的數學底層，確保與 R `fisher.test` 結果一致。

### 1.1 `MathUtils` 擴充
- [ ] **Log-Factorial 實作**: 實作 `log_factorial(n)`。
- [ ] **Haldane Correction**: 實作 `odds_ratio(a, b, c, d)`，內建 `+0.5` 平滑。
- [ ] **Cramér's V**: 實作 `calculate_cramers_v`，增加 `is_reliable` 檢查 (若太多格期望值 < 5 則標記不可靠)。

### 1.2 `FisherExact` 模組
- [ ] **2x2 Fisher Test (Two-sided)**: 
    *   **關鍵規格**: 採用 R 定義 (Sum of probabilities $\le P_{observed}$)。
    *   回傳 p-value 與 log-odds。
- [ ] **RxC Fisher Test (MC)**: 實作 Fisher-Freeman-Halton 演算法：
    *   **Stopping Rule (NEW)**: 計算 $\hat{p}$ 的 99% Binomial CI。若 CI 上界 < 0.05 或 CI 下界 > 0.05 則提早停止。
    *   預設抽樣 2,000 ~ 10,000 次。

### 1.3 驗證方式 (Verification)
*   **Unit Test**: 建立 `test_fisher_concordance.cpp`。
    *   構造 Borderline 2x2 表格，比較 C++ 輸出與 R `fisher.test` 輸出，誤差需 < 1e-6。

---

## Phase 2: 全域與局部檢定實作 (Association Tests)

**目標**: 實作業務邏輯，並嚴格執行 Gating 規則。

### 2.1 全域檢定 (`GlobalTest.cpp`)
- [ ] 實作 `compute_global_association`。
- [ ] **Gating Logic**: 
    *   若 Global P > 0.1 (寬鬆門檻)，標記 `passed_gate = false`。後續不執行 Bootstrap/Local Tests 以節省資源。

### 2.2 局部檢定 (`LocalTest.cpp`)
- [ ] 實作 **One-vs-Rest** 邏輯 (Two-sided Fisher)。
- [ ] 計算 `log_odds_ratio` (Haldane) 與 `delta_proportion`。

### 2.3 驗證方式 (Verification)
*   驗證 Gating 邏輯：確認被 Gating 掉的 Region 沒有執行後續計算。

---

## Phase 3: 進階結構驗證 (Structural Validation)

**目標**: 實作 Phase 3，排除 Dispersion 造成的誤判。

### 3.1 距離檢定 (`StructureTest.cpp`)
- [ ] **Data Filtering (關鍵)**: 在計算前，剔除 Overlap 不足的 Reads (invalid distance)。禁止插補。
- [ ] **Dispersion Check**: 計算 `mean_dist_to_centroid`。若群間差異顯著 (ANOVA F-test P < 0.05)，標記 `dispersion_warning`。
- [ ] **PERMANOVA**: 僅對 `passed_gate == true` 的 Region 執行。

### 3.2 Cluster Stability (Bootstrap)
- [ ] **Bootstrap**: 對 CpG sites 重抽樣 (N=200)。
- [ ] **Early Stop**: 若 ARI 連續 10 次 < 0.2 或 > 0.9，提早停止。

---

## Phase 4: 系統整合與輸出 (Integration)

**目標**: 整合至 RegionProcessor，輸出完整 CSV。

### 4.1 Output CSV Definition
- [ ] **`significance_report.csv` (One-line-per-region)**:
    *   `valid_flag`: 是否可計算 (排除 low coverage)。
    *   `passed_gate`: 是否通過 Global 篩選。
    *   `p_global_alt`, `cramers_v_alt`, `dispersion_warning`.
    *   `best_cluster_p`, `best_cluster_or`.
- [ ] **`significance_features.csv` (For ML)**:
    *   `n_reads`, `n_cpg`, `effective_overlap`, `mapq_mean`, `hp0_rate` (重要資訊量特徵)。
    *   所有統計 P-value (Log scale)。

---

## Phase 5: Python 校準模型 (Calibration)

**目標**: 訓練模型並產出 AUPRC 報告。

### 5.1 Training Pipeline
- [ ] **Truth Alignment**: 實作 Variant-level 對齊規則 (Anchor TP -> Region TP)。
- [ ] **FDR Control**: 實作 Global BH 校正 -> 計算 $q_{global}$。
- [ ] **Evaluation**: 計算 AUPRC (相對於 Baseline) 與 Brier Score。

---

## 注意事項 (Critical Notes)

1.  **Stop Rule**: 務必嚴格測試 MC Fisher 的 CI stopping，避免在 P=0.05 附近過早停止導致誤判。
2.  **NaN Handling**: 距離矩陣運算前必須 Filter Reads，不可讓 NaN 進入 SS 計算。
3.  **Defaults**: MC=2000, Boot=200, Perm=999。

---

開發順序：Phase 1 -> Phase 2 -> Phase 3 -> Phase 4 -> Phase 5。
