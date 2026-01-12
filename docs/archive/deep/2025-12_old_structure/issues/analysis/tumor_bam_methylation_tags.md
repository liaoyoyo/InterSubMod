# Tumor BAM 甲基化標籤問題

## 現狀

經過檢查多個 BAM 檔案，發現：

### 有 MM/ML Tags 的檔案
1. `HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (Normal樣本)
   - 檔案大小：136G
   - 包含：MM, ML, HP, PS tags
   - MM 格式：`C+h?,...;C+m?,...;` (5hmC + 5mC)

### 沒有 MM/ML Tags 的檔案
1. `HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` (Tumor, 經ClairS處理)
   - 檔案大小：264G
   - 只有標準 ONT tags
   - **問題**：經過 variant calling 後可能移除了甲基化標籤

2. `HCC1395_Tumor_ONT.GRCh38.sorted.bam` (Google somatic data)
   - 沒有甲基化標籤

## 可能的解決方案

### 選項 1：使用 Normal BAM 進行開發測試
- **優點**：有完整的 MM/ML tags，可以測試甲基化解析功能
- **缺點**：不是 tumor 樣本

### 選項 2：尋找原始的 Tumor BAM
- 需要找到未經 ClairS 處理的原始 ONT Tumor BAM
- 檔案名稱可能類似：`HCC1395_*_tumor_*_tagged.bam`

### 選項 3：修改測試策略
- Phase 1-2：使用 normal.bam 測試甲基化解析功能
- 後續階段：結合 tumor + normal 的實際分析

## 當前決策

**採用選項 3**：
1. 先用 normal.bam 完成 Phase 1-2 的功能開發與測試
2. 確保所有模組（BamReader, FastaReader, MethylationParser）都正常工作
3. 在後續的整合階段，若需要 tumor 甲基化資訊，再尋找原始檔案

## 測試程式調整

將 `ReadFilterConfig.require_mm_ml` 設為 `false`，這樣：
- 可以測試 BAM reading、CIGAR parsing、AltSupport 判斷
- 對有 MM/ML tags 的 reads 仍會嘗試解析甲基化
- 不會因為缺少 tags 而過濾掉所有 reads

