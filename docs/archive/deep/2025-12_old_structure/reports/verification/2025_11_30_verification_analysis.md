# 流程、輸出與程式碼驗證報告

本報告針對 `chr19_29283968` 位點的輸出結果進行詳細驗證，並結合程式碼分析回答關於設計目標、篩選標準及數據差異的問題。

## 1. 輸出完整性驗證 (Output Verification)

檢查輸出資料夾 `/big8_disk/liaoyoyo2001/InterSubMod/output/full_vcf_test/chr19_29283968/chr19_29282968_29284968/` 中的檔案：

* **`reads.tsv`**:
  * **現狀**: 包含 `read_id`, `read_name`, `chr`, `start`, `end`, `mapq`, `hp`, `alt_support`, `is_tumor`。
  * **缺失**: **缺少 Strand (正反股) 資訊**。
  * **分析**: 在 `ReadParser.cpp` 中，雖然有解析 FLAG，但未將 `BAM_FREVERSE` (0x10) 的狀態儲存到 `ReadInfo` 結構或寫入 `reads.tsv`。
  * **結論**: 目前輸出**不符合**「詳細確認...甲基矩陣資訊有符合要求實際完整輸出」的要求，因為缺乏正反股資訊，無法支援後續針對正反股的差異分析。

* **`methylation.csv`**:
  * **現狀**: 包含 Read x CpG 的甲基化機率值。
  * **完整性**: 格式正確，但同樣無法區分該 Read 來自哪一股。

## 2. 正反股分開處理的設計分析 (Design Analysis: Strand Separation)

針對 `/big8_disk/liaoyoyo2001/InterSubMod/docs/design` 的設計目標，分析是否需要將正反股數據分開處理：

* **建議**: **強烈建議區分正反股 (Strand-aware processing)**。
* **理由**:
    1. **生物學意義**: 雖然 CpG 甲基化通常是對稱的，但在複製過程中或特定調控下可能出現半甲基化 (Hemimethylation)。
    2. **技術誤差 (Artifacts)**: 定序誤差 (Sequencing bias) 往往具有股特異性 (Strand bias)。如果變異或甲基化模式只出現在某一條股上，很可能是偽陽性 (False Positive)。
    3. **驗證需求**: 您提到的「單股 3 條以上，相鄰雙股各有 4 條」是一個非常好的品質控制標準。這能確保訊號不是來自單一股的隨機誤差。
* **執行策略**:
  * **修改資料結構**: 在 `ReadInfo` 中增加 `strand` 欄位 (Forward/Reverse)。
  * **輸出更新**: `reads.tsv` 增加 `strand` 欄。
  * **分析層面**: 在聚類或統計時，可以選擇將正反股視為獨立的證據來源，或在過濾階段要求雙股皆有支持。詳見 `docs/architecture/strand_aware_processing_design.md`。

## 3. 篩選標準建立與合理性分析 (Filtering Criteria Analysis)

檢視 `ReadParser.cpp` 中的 `should_keep` 與 `determine_alt_support` 函式：

* **現有篩選條件**:
    1. **FLAG 過濾**: 排除 Secondary, Supplementary, Duplicate, Unmapped。 (標準且合理)
    2. **MAPQ**: `config_.min_mapq` (預設值需確認，通常建議 >= 20)。
    3. **Read Length**: `config_.min_read_length`。
    4. **MM/ML Tags**: 必須存在。 (必要條件)
    5. **Alt Support 判定**:
        * Read 必須覆蓋 SNV 位點。
        * 該位點的 Base Quality 必須 >= `min_base_quality`。
        * CIGAR 操作必須是 Match/Mismatch，不能是 Deletion 或 Ref Skip。

* **潛在過嚴的標準與影響**:
  * **Base Quality**: 如果 `min_base_quality` 設得太高 (例如 > 30)，可能會丟失很多 Nanopore reads (雖然近年品質提升，但部分區域仍較低)。這會導致 `reads.tsv` 中的 read 數量遠少於 `samtools view` 看到的數量。
  * **CIGAR 處理**: 如果 SNV 剛好位於 Indel 附近，複雜的 CIGAR 可能導致判定為 `UNKNOWN` 而被忽略。

## 4. 與 IGV/Samtools 數據差異之原因與程式碼對照

您觀察到 `reads.tsv` 的數量 (約 19 條) 與 IGV/Samtools (46 或 131 條) 有顯著差異。

**主要原因分析**:

1. **範圍定義不同**:
    * **程式**: 嚴格限制在 SNV 位點 ± Window Size (例如 2000bp) **且** 必須覆蓋 SNV 本身。
    * **IGV**: 顯示視窗內所有 Reads，包含未覆蓋到 SNV (例如斷在中間) 的 Reads。

2. **Alt Support 的嚴格判定 (`ReadParser.cpp:88`)**:
    * 程式碼邏輯：

        ```cpp
        // ReadParser.cpp
        if (snv_pos_0based < read_start || snv_pos_0based >= read_end) {
            return AltSupport::UNKNOWN; // 未覆蓋 SNV 的 reads 會被標記為 UNKNOWN
        }
        ```

    * **影響**: 許多 Reads 雖然在該區域，但如果沒有覆蓋到 `29283968` 這個特定點，或者在該點是 Deletion，或者 Base Quality 太低，它們可能被過濾掉或歸類為 `UNKNOWN`。如果後續步驟只選取 `REF` 或 `ALT` 的 reads，這些 reads 就會消失。

3. **Samtools 驗證**:
    * 指令 `samtools view ... chr19:29283968-29283968` 會找出**覆蓋該點**的 reads。
    * 指令 `samtools view ... chr19:29282968-29284968` 會找出**區域內**的 reads (數量會多很多)。
    * **程式碼行為**: 程式是先抓區域 (`querys`)，再過濾 (`should_keep`)，再判定 SNV (`determine_alt_support`)。

**程式碼相關部分**:

* **`src/core/ReadParser.cpp`**:
  * **Line 10-39 (`should_keep`)**: 負責初步過濾 (MAPQ, Length, Tags)。如果這裡太嚴，Read 直接丟棄。
  * **Line 88-181 (`determine_alt_support`)**: 負責判定 SNV 狀態。如果 CIGAR 解析失敗或品質低，回傳 `UNKNOWN`。

**總結與建議**:

1. **修正程式碼**: 立即在 `ReadInfo` 與 `ReadParser` 中加入 `strand` 處理。
2. **調整參數**: 檢查 `min_mapq` 與 `min_base_quality` 設定，確認是否過濾了太多有效 reads。
3. **驗證方法**: 使用 `samtools view -q [MIN_MAPQ] -F 3844` (排除次要/重複/未比對) 來模擬程式的過濾邏輯，比較數量是否接近。
