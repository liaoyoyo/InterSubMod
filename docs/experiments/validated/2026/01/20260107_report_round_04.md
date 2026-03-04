<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Round 4 分析報告

> **日期**: 2026-01-07  
> **輸出目錄**: `output/bip8_disk_output/20260107_all-with-w5000_4/`

---

## 1. 此次程式碼修改

### 修改內容
改用 VerificationClass 篩選，取代 V 值閾值：

```cpp
// Round 3 (V 值篩選)
bool is_significant = r.passed_gating && 
                      (r.global_p_value <= 0.05) && 
                      (r.cramers_v >= 0.05);

// Round 4 (VerificationClass 篩選)
bool is_significant = r.passed_gating && 
                      (r.global_p_value <= 0.05) && 
                      (r.verification_class == "Strong");
```

---

## 2. 分析重點與結果

### ⚠️ 主要指標 — 結果惡化！

| 指標 | Round 2/3 (V≥0.1) | **Round 4 (Strong)** | 變化 |
|:---|:---:|:---:|:---|
| TP 顯著率 | 6.10% | 24.55% | ↑18.4% |
| FP 顯著率 | 2.09% | **37.08%** | ↑35.0% |
| TP/FP 比 | 2.92 | **0.66** | ↓4.4x |

### VerificationClass 原始分布
| 類別 | FP | TP |
|:---|:---:|:---:|
| **Strong** | **37.08%** | 24.55% |

---

## 3. 結果觀察總結

### ❌ 失敗原因
1. **VerificationClass 判定本身有問題**：FP 在 Strong 類別的比例 (37%) 反而高於 TP (24.5%)
2. **這與直覺完全相反**：Strong 原本應代表「最可信」的位點
3. **根本原因**：VerificationClass 的 "Strong" 定義是 `cluster_significant && label_significant`，但 FP 的聚類特性可能更容易滿足這個條件

### 關鍵發現
**VerificationClass 不能用於區分 TP/FP**，它只描述單倍體與聚類的一致性，不代表該位點是真正的生物信號

---

## 4. 下次修改重點

### Round 5 計畫：回到 V ≥ 0.1 + 深度過濾

最佳方案仍是 Round 2，但可嘗試加入深度過濾：

```cpp
bool is_significant = r.passed_gating && 
                      (r.global_p_value <= 0.05) && 
                      (r.cramers_v >= 0.1) &&
                      (r.num_reads >= 30);  // 新增深度閾值
```

**預期**：
- 進一步減少低深度假顯著
- 維持 TP/FP 比 > 2
