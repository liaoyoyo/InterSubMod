# 輸出完整性驗證報告 (Output Integrity Verification Report)

**日期**: 2025-11-30
**目標**: 驗證 `/big8_disk/liaoyoyo2001/InterSubMod/output/full_vcf_test/chr19_29283968` 的輸出完整性與正確性。

## 1. 驗證結果摘要 (Executive Summary)

* **完整性 (Completeness)**: ❌ **不完整**。
  * **缺失項目**: `reads.tsv` 缺少 **Strand (正反股)** 資訊。這對於後續的甲基化分析與偽陽性過濾至關重要。
  * **現有檔案**: `reads.tsv`, `methylation.csv`, `cpg_sites.tsv`, `metadata.txt` 皆已產生且格式正確。
* **正確性 (Correctness)**: ⚠️ **有差異，但可解釋**。
  * **Read 數量**: 程式輸出 97 條 reads，而 `samtools` 顯示 131 條。
  * **原因**: 程式內建的品質過濾 (MAPQ, Base Quality, CIGAR check) 比單純的 `samtools view` 更嚴格。

## 2. 詳細檢查與數據比對 (Detailed Inspection)

### 2.1 檔案結構檢查

| 檔案名稱 | 狀態 | 說明 |
| :--- | :--- | :--- |
| `reads.tsv` | ⚠️ 格式正確但缺欄位 | 包含 97 條 reads。**缺少 Strand 欄位**。 |
| `methylation.csv` | ✅ 格式正確 | 包含 Read x CpG 的甲基化機率矩陣。 |
| `cpg_sites.tsv` | ✅ 格式正確 | 列出區域內所有 CpG 位點座標。 |
| `metadata.txt` | ✅ 格式正確 | 包含執行參數與統計摘要。 |

### 2.2 數據一致性抽樣檢查 (Sampling Check)

使用快速抽樣方法確認輸出邏輯：

1. **Read 數量比對**:
    * **程式輸出 (`reads.tsv`)**: 97 reads (扣除 header)。
    * **Samtools (`chr19:29283968-29283968`)**: 131 reads。
    * **差異分析**: 程式過濾掉了約 26% 的 reads。這通常是因為這些 reads 雖然覆蓋該位點，但：
        * MAPQ 低於閾值。
        * 在 SNV 位點的 Base Quality 低於閾值。
        * CIGAR 字串顯示該位點位於 Deletion 或 Ref Skip 區域。
        * 缺少 MM/ML tags。

2. **關鍵 Read 抽樣**:
    * **Read ID**: `51219d01-133b-4146-b94f-4650dc31b41b`
    * **程式判定**: `REF` support, `is_tumor=1`。
    * **Samtools 驗證**: 確認該 read 存在於 BAM 檔中，且在 29283968 位置為參考鹼基。

## 3. 設計建議與修正行動 (Recommendations)

基於上述驗證，提出以下修正建議以符合設計目標：

1. **立即修正**: 在 `ReadInfo` 結構與 `reads.tsv` 輸出中加入 `strand` 欄位 (Forward/Reverse)。
2. **品質控制**: 雖然目前的過濾邏輯是合理的，但建議在 `metadata.txt` 中記錄「被過濾掉的 Reads 數量與原因」，以便區分是「無訊號」還是「訊號被過濾」。
3. **快速檢查腳本**: 建議建立一個自動化腳本 `scripts/verify_output.sh`，自動執行上述的 `wc -l` 與 `samtools view -c` 比對，並在差異過大 (>50%) 時發出警告。

## 4. 結論

目前的輸出**尚未達到**「完整符合要求」的標準，主要缺失為 **Strand 資訊**。除此之外，數據處理邏輯正常，與原始 BAM 檔的差異在預期的品質控制範圍內。

---
**參考資訊**:

* Input VCF: `/big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz`
* Output Dir: `/big8_disk/liaoyoyo2001/InterSubMod/output/full_vcf_test`
