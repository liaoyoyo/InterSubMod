<!--
建立時間: 2026-03-23
問題類型: C++ 方法 | 特徵設計
影響 track: TO | 兩者
狀態: pending_decision
-->

# Quality Score 懲罰偏重分析 審查文件

## 問題描述

ISM 計算 `heuristic_score`（啟發式得分）用於排序和輔助判斷。公式位於：

**程式碼位置**：`src/core/SignificanceAnalyzer.cpp:375-407`

```cpp
double SignificanceAnalyzer::compute_heuristic_score(const SignificanceResult& result) {
    double score = 0.0;
    double best_p = std::min(p_val_alt, p_val_hp);

    // Primary: -log10(best_p)
    if (best_p > 0) score = -std::log10(best_p);
    else score = 20.0;  // p=0 上限

    // Bonus: Cramér's V × 2.0
    score += best_v * 2.0;

    // Penalty: dispersion warning → × 0.5
    if (result.dispersion_warning) score *= 0.5;

    // Penalty: PERMANOVA p > 0.05 → × 0.7
    if (result.permanova.valid && result.permanova.p_value > 0.05) score *= 0.7;

    return score;
}
```

**注意**：以上是 `heuristic_score`（在 significance.json 輸出）。

另外，label_first_metrics.tsv 中有一個不同的 `QualityScore` 欄位，其計算邏輯不同（位於 RegionProcessor.cpp）。需確認兩者的關係。

---

## 問題 1：PERMANOVA 懲罰無效

```cpp
if (result.permanova.valid && result.permanova.p_value > 0.05) score *= 0.7;
```

**問題**：`enable_permanova = false`（RegionProcessor.cpp:1146），所以 `result.permanova.valid` 永遠是 false。

**結果**：這個 PERMANOVA 懲罰條件**永遠不會被觸發**，是死代碼（dead code）。

---

## 問題 2：QualityScore（label_first_metrics.tsv）的懲罰設計

**待確認**：`QualityScore` 欄位的計算邏輯位置，以及 -20/-15/-25/+10/+15 懲罰分數的來源。

可能位置：`src/core/RegionProcessor.cpp`（需搜尋 quality_score）

---

## 量化影響

| 問題 | 類型 | 影響 |
|------|------|------|
| PERMANOVA 懲罰死代碼 | Bug（死代碼） | 無實際影響（已停用） |
| heuristic_score 未連接最終過濾 | 架構問題 | 無直接 F1 影響 |
| QualityScore 設計待確認 | 待審查 | 未知 |

---

## 修改選項

### 方案 A：不修改
- heuristic_score 和 QualityScore 均非最終過濾依據
- **F1 影響**：0

### 方案 B：移除死代碼（PERMANOVA 懲罰條件）
- 修改位置：`src/core/SignificanceAnalyzer.cpp:402-405`
- 移除或加上 `// TODO: re-enable when permanova is active` 注解
- **F1 影響**：0（無行為改變）
- **成本**：低（refactor 類型 commit）

### 方案 C：重新設計 heuristic_score 使其有效
- 連接到最終判決或排序
- **成本**：高（需整體架構設計）

---

## 用戶決策

**選擇**：[ ] A / [ ] B（清理死代碼）/ [ ] C
**日期**：
**理由**：
