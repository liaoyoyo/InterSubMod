<!--
建立時間: 2026-01-13 23:20
目標: RegionProcessor.cpp 重構報告
處理範圍: process_single_region 方法拆分與簡化
關聯檔案:
  - src/core/RegionProcessor.cpp
  - include/core/RegionProcessor.hpp
-->

# RegionProcessor 重構報告

## 重構目標

將 `process_single_region` 方法從約 570 行簡化到約 120 行，提升程式碼可讀性與維護性。

## 重構策略

採用 **Extract Method** 重構模式，將大型方法拆分為多個職責單一的輔助方法。

## 新增輔助方法

| 方法名稱 | 職責 | 行數 |
|---------|------|------|
| `process_reads()` | 處理 BAM reads、過濾、甲基化解析 | ~55 |
| `build_methylation_matrix()` | 建構 Eigen 格式甲基化矩陣 | ~50 |
| `compute_and_write_distance_matrix()` | 計算並輸出距離矩陣 | ~40 |
| `perform_clustering_and_significance()` | 階層聚類與顯著性分析 | ~130 |
| `write_strand_specific_trees()` | 輸出鏈特異性樹檔案 | ~40 |

## 程式碼變更摘要

### RegionProcessor.hpp
- 修改 `process_reads` 簽名：接受已 fetch 的 reads 向量（避免重複 fetch）
- 新增 4 個輔助方法宣告

### RegionProcessor.cpp
- `process_single_region`: 570 行 → 117 行（減少 79%）
- 總檔案行數: ~1460 行 → 1005 行（減少 31%）

## 測試驗證

### 編譯測試
```
✓ make -j$(nproc) 成功
```

### 功能驗證
| 測試項目 | 基準值 | 重構後 | 狀態 |
|---------|--------|--------|------|
| TP Records | 30476 | 30476 | ✓ |
| FP Records | 4822 | 4822 | ✓ |
| TP Sig Rate | 6.103163% | 6.103163% | ✓ |
| FP Sig Rate | 2.094567% | 2.094567% | ✓ |
| Cramér's V AUC | 0.519386 | 0.519386 | ✓ |

## 改進效果

1. **可讀性提升**: 主方法流程清晰，一目了然
2. **可測試性提升**: 輔助方法可獨立單元測試
3. **可維護性提升**: 修改特定功能無需觸及整體邏輯
4. **程式碼復用**: 輔助方法可被其他方法調用

## 後續建議

1. **新增單元測試**: 針對各輔助方法撰寫邊界條件測試
2. **效能剖析**: 確認重構未引入額外開銷
3. **文檔更新**: 更新 API 文檔說明新方法

## 修改檔案清單

- `include/core/RegionProcessor.hpp` (修改)
- `src/core/RegionProcessor.cpp` (修改)
