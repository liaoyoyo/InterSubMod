# 20241224_Phase2_GlobalTest_LocalTest_Test

## 測試日期
2025-12-24

## 測試範圍
Phase 2: 全域與局部檢定
- `GlobalTest`: Fisher-Freeman-Halton RxC 檢定、Gating 邏輯、Cramér's V
- `LocalTest`: One-vs-Rest 2x2 Fisher 檢定、delta proportion、最佳分群識別

## 測試結果

**全部通過：10/10 測試**

### GlobalTestTest (6 tests)

| 測試名稱 | 結果 | 驗證內容 |
|---------|------|---------|
| AlleleTest_StrongAssociation | ✅ PASS | 強關聯 p < 0.01，Gating 通過 |
| AlleleTest_NoAssociation | ✅ PASS | 無關聯 p > 0.1，Gating 失敗 |
| HPTest_Works | ✅ PASS | HP1 vs HP2 檢定運作正常 |
| GatingLogic | ✅ PASS | 嚴格門檻 0.05 正確過濾 |
| EmptyInput | ✅ PASS | 空輸入標記為 invalid |
| TestAll | ✅ PASS | 同時執行三維度檢定 |

### LocalTestTest (4 tests)

| 測試名稱 | 結果 | 驗證內容 |
|---------|------|---------|
| BasicLocalTest | ✅ PASS | One-vs-Rest 檢定運作，best_p_value < 0.001 |
| ClusterCounts | ✅ PASS | 分群計數正確（count_alt, count_hp1, count_tumor） |
| DeltaProportion | ✅ PASS | delta_proportion_alt = 1.0 正確 |
| BestDimensionIdentified | ✅ PASS | 正確識別 HP 為最佳維度 |

## 功能驗證

### Gating 邏輯
- 預設門檻：p > 0.1 標記為 `passed_gate = false`
- 通過 Gating 的 Region 才會執行後續 Local Tests / Structure Tests

### One-vs-Rest 策略
- 對每個 Cluster 建立 2x2 表格
- 計算 delta proportion：P(target|cluster) - P(target|rest)
- 追蹤最佳分群（最低 p-value）與對應維度

## 累計測試統計

| Phase | 測試數 | 通過數 | 通過率 |
|-------|-------|-------|--------|
| Phase 1 | 25 | 25 | 100% |
| Phase 2 | 10 | 10 | 100% |
| **總計** | **35** | **35** | **100%** |

## 結論

Phase 2 全域與局部檢定模組實作完成：
1. **GlobalTest**: Fisher-Freeman-Halton MC 檢定運作正常，Gating 邏輯正確
2. **LocalTest**: One-vs-Rest 策略正確識別富集分群

可進入 Phase 3 開發（PERMANOVA、Dispersion、Bootstrap）。
