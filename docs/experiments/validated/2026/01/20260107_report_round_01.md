<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Round 1 分析報告

> **日期**: 2026-01-07  
> **輸出目錄**: `output/bip8_disk_output/20260107_all-with-w5000_1/`

---

## 1. 此次程式碼修改

### 修改檔案
| 檔案 | 修改內容 |
|:---|:---|
| `GlobalTest.hpp` | 新增 `min_cramers_v = 0.1` 配置項 |
| `GlobalTest.cpp` | 修改 `apply_gating()`: `passed_gate = p_pass && v_pass` |

### 修改邏輯
```cpp
// 修改前
result.passed_gate = (result.fisher_ffh.p_value <= config_.gating_p_threshold);

// 修改後
bool p_pass = (result.fisher_ffh.p_value <= config_.gating_p_threshold);
bool v_pass = (result.cramers_v >= config_.min_cramers_v);
result.passed_gate = p_pass && v_pass;
```

---

## 2. 分析重點與結果

### 主要指標比較

| 指標 | 修改前 (20260106) | 修改後 (Round 1) | 變化 |
|:---|:---:|:---:|:---|
| TP 顯著率 | 28.58% | 28.49% | ↓0.09% |
| FP 顯著率 | 40.81% | 40.92% | ↑0.11% |
| TP Count | 8,710 | 8,683 | -27 |
| FP Count | 1,968 | 1,973 | +5 |
| ROC AUC (V) | 0.519 | 0.519 | 無變化 |

### VerificationClass 分布

| 類別 | FP (修改後) | TP (修改後) | 差異 |
|:---|:---:|:---:|:---|
| Strong | 37.12% | 24.50% | FP 高 12.6% |
| Weak | 31.34% | 43.54% | TP 高 12.2% |

---

## 3. 結果觀察總結

### 問題發現
1. **修改幾乎無效**：FP 顯著率未下降，反而微幅上升
2. **原因分析**：
   - `apply_gating()` 只影響 `passed_gate` 欄位
   - 最終 `Significant` 判定在 `write_significance_summary()` 中為：
     ```cpp
     bool is_significant = r.passed_gating && (r.global_p_value <= 0.05);
     ```
   - **這個判定沒有直接使用 Cramér's V 閾值**
3. **Gating 效果**：雖然 `passed_gating` 可能減少，但許多 FP 仍然有 p ≤ 0.05，所以仍被標記為顯著

### 核心問題
Gating 機制主要控制是否進行後續分析（LocalTest, PERMANOVA），但不直接影響最終 `Significant` 判定。FP 的 GlobalP 通常較低（平均 0.429），即使 V 值低也能通過 p ≤ 0.05。

---

## 4. 下次修改重點

### Round 2 計畫

1. **直接在 `Significant` 判定加入 V 值閾值**：
   ```cpp
   // RegionProcessor.cpp write_significance_summary()
   bool is_significant = r.passed_gating && 
                         (r.global_p_value <= 0.05) && 
                         (r.cramers_v >= 0.1);  // 新增
   ```

2. **考慮更嚴格的閾值**：
   - V ≥ 0.15 或 V ≥ 0.2
   - 深度 ≥ 20 reads

### 修改檔案
- `src/core/RegionProcessor.cpp` (line 907 附近)
