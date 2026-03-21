# MM/ML Tags 檢測問題調查

## 問題
儘管 samtools 能看到 normal.bam 中的 MM/ML tags，但程式仍然回報 "Reads with MM/ML tags: 0"。

## 可能原因

### 1. 過濾邏輯問題
`ReadParser::should_keep()` 檢查 `require_mm_ml = true`，可能在檢測 tags 時有問題。

### 2. 測試方式
需要檢查：
- `bam_aux_get(b, "MM")` 是否正確返回指標
- reads 可能因為其他原因（MAPQ, length）被過濾掉了
- 重複的 reads？（看到很多相同的 read name）

### 3. 發現：Reads 重複
從測試輸出看到：
- Read 2 和 3: 同一個 read `0c25d1ae-b184-4e79-8c14-367a69bedae8`
- Read 4 和 5: 同一個 read `6d164505-0425-431e-bf15-80e7d8088da9`

這可能意味著：
- 有 split alignments（supplementary/secondary）
- 或者是真的重複

但我們的過濾器應該已經過濾掉 secondary/supplementary...

## 解決方案

需要修改測試程式來詳細輸出：
1. 檢查每個 read 是否有 MM/ML tags
2. 在 should_keep() 之前檢查
3. 顯示過濾原因

