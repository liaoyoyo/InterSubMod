# Phase 3 完成報告

## 完成日期
2025-11-23

## 已完成項目

### MatrixBuilder 類別 ✅
- **功能**：建構 Read × CpG 甲基化矩陣
- **數據結構**：使用 `std::vector<std::vector<double>>` 代替 Eigen（避免依賴）
- **特點**：
  - 使用 `-1.0` 表示無覆蓋（而非 NaN）
  - 自動收集唯一 CpG 位點並排序
  - 稀疏矩陣設計（平均 56.6% sparsity）
  - 支援 `clear()` 重複使用

### RegionWriter 類別 ✅
- **功能**：將 region 資料輸出到結構化目錄
- **輸出結構**：
  ```
  output/region_0000/
    ├── metadata.txt      # Region 與 SNV 資訊
    ├── reads.tsv         # Read 列表（9欄）
    ├── cpg_sites.tsv     # CpG 位點列表
    └── methylation.csv   # 甲基化矩陣（CSV）
  ```
- **格式**：
  - TSV 用於表格數據
  - CSV 用於矩陣（使用 "NA" 表示無覆蓋）
  - 固定精度（4位小數）

### 測試結果 ✅
- ✅ 154 reads 成功處理
- ✅ 147 reads 包含甲基化數據
- ✅ 59 個唯一 CpG 位點
- ✅ 矩陣維度：154 × 59 (0.07 MB)
- ✅ 稀疏度：56.6%
- ✅ 處理時間：238 ms
- ✅ 所有檔案成功輸出

## 技術決策

### 1. 使用 std::vector 代替 Eigen
**原因**：
- Eigen3 未安裝於系統
- `std::vector<std::vector<double>>` 足夠簡單且高效
- 避免額外依賴，提升可移植性

### 2. 使用 -1.0 代替 NaN
**原因**：
- 簡化比較邏輯（`val < 0.0` vs `std::isnan(val)`）
- 更好的 CSV 輸出兼容性
- 避免浮點數精度問題

### 3. 輸出格式設計
- **TSV**：便於 pandas/R 讀取
- **CSV**：便於矩陣可視化與分析
- **Metadata**：記錄處理資訊與性能指標

## 性能指標

### 矩陣統計
- Total cells: 9,086
- NA (no coverage): 5,143 (56.6%)
- Zero methylation: 7 (0.08%)
- Non-zero methylation: 3,936 (43.3%)

這表示：
- 每個 read 平均覆蓋 ~26 CpG 位點（43.3% of 59）
- 大部分甲基化值 > 0（高度甲基化區域）

### 記憶體使用
- Matrix: 0.07 MB
- 平均每個 cell: ~8 bytes (sizeof(double))

## 檔案輸出驗證

```bash
$ ls -lh output/region_0000/
-rw-rw-r-- cpg_sites.tsv     (839 bytes)   # 59 CpG × 3 columns
-rw-rw-r-- metadata.txt      (244 bytes)   # 元數據
-rw-rw-r-- methylation.csv   (43K)         # 154 × 59 matrix
-rw-rw-r-- reads.tsv         (12K)         # 154 reads × 9 columns
```

### metadata.txt 內容
```
Region ID: 0
Region: chr17:7576000-7580000
Region Size: 4001 bp

SNV ID: 0
SNV Position: chr17:7578000
SNV: C -> T
SNV Quality: 100

Num Reads: 154
Num CpG Sites: 59
Matrix Dimensions: 154 × 59

Processing Time: 237.98 ms
Peak Memory: 0.00 MB  # 尚未實作監控
```

## 下一步驟

### Phase 4: 平行化處理
- [ ] 實作 RegionProcessor class
- [ ] 支援多 region 並行處理
- [ ] Thread-local BAM readers
- [ ] 整合所有模組

### Phase 5: 完整測試
- [ ] 載入真實 SNV table
- [ ] 測試 32 SNVs × 4 threads
- [ ] 資源監控（時間、記憶體）
- [ ] 驗證所有輸出正確性

## 總結

Phase 3 成功實作矩陣建構與結構化輸出功能。所有數據正確存儲並輸出到指定格式。程式碼模組化設計良好，為後續的平行化處理奠定了堅實基礎。

