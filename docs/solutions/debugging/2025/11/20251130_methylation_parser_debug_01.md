# MethylationParser Debug - Why No Calls?

## Issue
tumor.bam 確實包含 MM/ML tags，經 samtools 確認第一個 read 有完整的標籤：
```
MM:Z:C+h?,...;C+m?,...;
ML:B:C,17,3,10,102,108,28,166,42,3,5,3,26,4,2,16,...
```

但 `MethylationParser::parse_read()` 回傳空的 vector。

## 可能原因

1. **MM tag 解析錯誤**
   - `C+m?` 片段的分隔符號處理不正確
   - delta encoding 計算錯誤

2. **CpG 座標映射錯誤**
   - CIGAR 遍歷邏輯問題
   - ref_pos 與 ref_seq 索引不一致

3. **CpG context 驗證失敗**
   - ref_seq 範圍不正確
   - 大小寫轉換問題

4. **ML array 長度不匹配**
   - C+h? 和 C+m? 共用同一個 ML array
   - 需要正確計算 offset

## 解決策略

1. 加入 debug 輸出到 `parse_read` 檢視：
   - MM string 內容
   - delta array 長度
   - ML array 長度
   - 找到的 CpG 數量

2. 檢查 MM parsing 是否正確提取 `C+m?` 片段

3. 驗證 CIGAR 映射邏輯

## Next Steps
創建一個詳細的 debug 版本 來找出問題

