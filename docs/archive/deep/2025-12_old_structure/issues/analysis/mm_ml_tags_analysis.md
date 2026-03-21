# MM/ML Tags 分析除錯紀錄

## 日期: 2025-11-22

## 問題描述
在測試 Phase 1 & 2 功能時，發現 HCC1395 tumor.bam 沒有檢測到 MM/ML tags。

## 調查結果

### BAM 檔案差異

1. **tumor.bam** (`HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`)
   - 包含標準 ONT tags (NM, ms, AS, tp, cm, etc.)
   - **不包含** MM/ML methylation tags
   - 原因：可能是經過 variant calling 處理後移除了甲基化標籤

2. **normal.bam** (`HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam`)
   - **包含** MM/ML tags
   - **包含** HP/PS phasing tags
   - MM tag 格式：`C+h?,...;C+m?,...;` (同時包含 5hmC 和 5mC)

### MM Tag 格式範例

```
MM:Z:C+h?,16,0,1,0,20,...;C+m?,16,0,1,0,20,...;
ML:B:C,2,4,117,158,4,1,28,15,...  (非常長的 uint8 陣列)
HP:i:1
PS:i:7272060
```

### 關鍵發現

1. **多種修飾類型**：MM tag 包含兩種修飾類型 (5hmC 和 5mC)，delta arrays 相同
2. **Delta encoding**：使用逗號分隔的整數表示 skip counts
3. **ML array 長度**：對於有兩種修飾類型的情況，ML array 長度 = 2 × delta array 長度

## 程式碼調整

### MethylationParser 設計確認

當前設計已能處理這種情況：
- `parse_mm_tag()` 函式專門搜尋 "C+m?" 段落
- 只處理 5mC，忽略 5hmC
- 正確處理 ML array 與 delta array 長度不匹配的情況（當有多種修飾時）

### 測試程式調整

需要使用 normal.bam 進行測試，因為它包含完整的甲基化標籤。

## 解決方案

1. 使用 normal.bam 作為測試數據
2. 確認 MethylationParser 能正確解析多修飾類型的 MM tags
3. 後續若需要使用 tumor BAM 的甲基化資訊，需找到未經處理的原始 BAM 檔案

## 後續測試計畫

- [ ] 使用 normal.bam 測試完整的 methylation parsing
- [ ] 驗證 CpG 座標與參考基因組的一致性
- [ ] 測試 HP/PS tags 的正確解析
- [ ] 確認多種修飾類型時 ML array 索引的正確性

