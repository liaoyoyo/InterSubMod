# ISM MM/ML 甲基化解析源碼稽核 — `MethylationParser`

**稽核日期**: 2026-06-13
**源檔**:
- `src/core/MethylationParser.cpp` (291 行, mtime Mar 27 02:21)
- `include/core/MethylationParser.hpp` (103 行, mtime Jan 13 18:33)
**caller**: `src/core/ReadAggregator.cpp:62` → `methyl_parser.parse_read(b, ref_seq, region_start)`

> 紀律：以下每條結論皆引行號 + 原始碼片段。只讀實際檔案，未跑、未臆測。讀不確定處標 **[uncertain]**。

---

## L0 一眼結論

ISM 的 MM/ML 解析**只取 5mC (`C+m?`) 一個通道**，**完全忽略 5hmC (`C+h?`)** — 不是 max-collapse、不是雙列，而是「跳過 5hmC 的機率、只輸出 5mC 機率」。ML byte→機率為 **線性 `ML[i] / 255.0`**（無 0.5 偏移、無分箱）。每個 call 只帶 `{ref_pos (1-based), probability}`，**不帶 modified-base 種類、不帶 strand**。CpG 對位用 reference 序列驗證（forward 查 `C`+`G`，reverse 查 `G` 且前一鹼基 `C`）。**未覆蓋的 CpG 不會被標 NA — 根本不會被產出**（parse 只列舉「read 上被標記為修飾的 C」，非「所有 CpG」）。

---

## L1 關鍵發現（每條附 檔:行）

### F1 — MM/ML tag 取得與 ML 格式硬性檢核
`parse_read` 用 htslib `bam_aux_get` 取 `"MM"` 與 `"ML"`；缺任一即回空。ML 只接受 `B`(array) + `C`(uint8) 型別，否則回空。
- `MethylationParser.cpp:20-25`：
  ```cpp
  uint8_t* mm_aux = bam_aux_get(b, "MM");
  uint8_t* ml_aux = bam_aux_get(b, "ML");
  if (!mm_aux || !ml_aux) { return calls; }
  ```
- `MethylationParser.cpp:31-33`：
  ```cpp
  if (ml_aux[0] != 'B' || ml_aux[1] != 'C') { return calls; }  // Invalid ML format
  ```
- ML 長度為 little-endian uint32（offset +2 起 4 bytes），資料指標在 offset +6：
  - `MethylationParser.cpp:35-36`：`uint32_t ml_len = le_to_u32(ml_aux + 2); const uint8_t* ml_data = ml_aux + 6;`
  - `le_to_u32` 定義：`MethylationParser.cpp:10-13`

### F2 — 只解析 5mC `C+m?`，5hmC `C+h?` 被「跳過機率、不輸出」
程式以指標掃描 MM tag 的 `;`-分隔 block，唯一目標是 `"C+m?"`（hard-coded，`tlen=4`）。非目標 block（含 `C+h?`）**不解析其 delta，只數其 delta 個數來推進 `ml_offset`**，藉此把 ML 陣列中對應 5hmC 的那段機率「跳過」。
- target 常數：`MethylationParser.cpp:60-63`：
  ```cpp
  const char* target = "C+m?";
  const size_t tlen  = 4;
  bool is_target = (... >= tlen) && std::memcmp(block_start, target, tlen) == 0;
  ```
- 命中 target 才呼叫 `parse_mm_tag(mm_str, "C+m?")` 並 break：`MethylationParser.cpp:65-69`
- 非 target block 數 delta 個數推進 `ml_offset`：`MethylationParser.cpp:71-82`（`delta_count++` 補最後一個無逗號值，`ml_offset += delta_count`）
- 註解明確：ML 陣列順序為 `[C+h? probs..., C+m? probs...]`（`MethylationParser.cpp:41`），故需 `ml_offset` 跳過前段。
- `parse_mm_tag` 預設 `mod_code = "C+m?"`：`MethylationParser.hpp:76`
- **結論**：5hmC 機率**不參與任何運算**（非 max-collapse、非 sum、非雙列），純粹被 offset 略過。每個輸出 CpG 只有「5mC 機率」一個值。

### F3 — ML byte→機率為線性 `ML[i] / 255.0f`（無偏移、無分箱、無 clip）
forward 與 reverse 兩路皆同公式：
- forward：`MethylationParser.cpp:180`：`float prob = ml_data[ml_offset + delta_idx] / 255.0f;`
- reverse：`MethylationParser.cpp:140`：`float prob = ml_data[ml_offset + delta_idx] / 255.0f;`
- header 文件一致：`MethylationParser.hpp:39`：`Probability = ML[i] / 255.0`
- 範圍：byte 0→0.0，byte 255→1.0；`MethylCall.probability` 文件標 `[0.0, 1.0]`（`MethylationParser.hpp:16`）。
- **無 0.5/255 居中偏移、無 thresholding、無 binarization** — 直接除 255。

### F4 — strand 處理：reverse 反向迭代 + 目標鹼基切換為 'G'
MM tag delta 以 **原始 read 5'→3'** 排序；BAM SEQ 對 reverse strand 是 RevComp（3'→5'），故 reverse 必須**從 `seq_len-1` 倒著數**並以 `'G'` 為目標鹼基（原 read 的 C 在 RevComp SEQ 變 G）。
- 設計註解：`MethylationParser.cpp:103-110`
- forward：正向迭代、target `'C'`：`MethylationParser.cpp:162-164`
- reverse：反向迭代、target `'G'`：`MethylationParser.cpp:124-126`
  ```cpp
  char target_base = 'G';
  for (int seq_idx = seq_len - 1; seq_idx >= 0; seq_idx--) { ... }
  ```
- strand 由 `bool is_reverse = bam_is_rev(b);` 決定：`MethylationParser.cpp:114`

### F5 — CpG 對位：以 reference 序列驗證（forward C-then-G / reverse G-preceded-by-C）
delta 命中後，用 `seq_to_ref` 把 read index 映回 ref 0-based 座標，再用 `ref_offset = ref_pos_0based - ref_start_pos` 在傳入的 `ref_seq` 上**驗證 CpG dinucleotide 才收錄**。
- forward 驗證：`MethylationParser.cpp:176-178`：
  ```cpp
  if (ref_offset >= 0 && (size_t)(ref_offset + 1) < ref_seq.size() &&
      ref_seq[ref_offset] == 'C' && ref_seq[ref_offset + 1] == 'G') {  // Simple is_cpg_site inline
  ```
- reverse 驗證（C 在前、G 在後，回報 C 位）：`MethylationParser.cpp:138-139`：
  ```cpp
  if (ref_offset > 0 && (size_t)(ref_offset) < ref_seq.size() &&
      ref_seq[ref_offset] == 'G' && ref_seq[ref_offset - 1] == 'C') {
  ```
- 等價的 helper（未在主路徑被呼叫，主路徑用 inline 版）：`is_cpg_site` `MethylationParser.cpp:283-288`（`ref_seq[offset]=='C' && ref_seq[offset+1]=='G'`）。
- **CpG 驗證使用 reference 序列**（非 read SEQ），故 read 上的 mismatch/變異不會被誤判為 CpG。

### F6 — 座標系：回報 1-based ref_pos，forward 與 reverse 回報位置不同
- forward 收 `ref_pos_0based + 1`（C 位轉 1-based）：`MethylationParser.cpp:181`：`calls.emplace_back(ref_pos_0based + 1, prob);`
- reverse 收 `ref_pos_0based`（其註解：`ref_pos_0based - 1 + 1 (1-based) = ref_pos_0based`，即 reverse 的修飾 G 位減 1 得 C 位、再轉 1-based）：`MethylationParser.cpp:142-143`
- `MethylCall.ref_pos` 文件標 "1-based genomic position (hg38)"：`MethylationParser.hpp:15`

### F7 — seq→ref 映射用 CIGAR；insertion/soft-clip 映 -1（不在 reference）
`build_seq_to_ref_map` 走訪 CIGAR：M/=/X 同步推進 seq 與 ref 並記錄；I/S 只推進 seq（保持 -1）；D/N 只推進 ref；H 不動。
- `MethylationParser.cpp:232-281`
- 預設全 -1：`MethylationParser.cpp:233`：`std::vector<int32_t> seq_to_ref(b->core.l_qseq, -1);`
- M/=/X：`:245-256`；I/S：`:258-262`；D/N：`:264-268`；H：`:270-272`
- 主路徑在 `ref_pos_0based >= 0` 才處理 → **落在 insertion/soft-clip 的修飾 C 被丟棄**（forward `:172`、reverse `:134`）。

### F8 — 未覆蓋 CpG 的「NA」語意：不會被產出（非顯式 NA 標記）
此 parser **不列舉所有參考 CpG**，只列舉「read MM tag 標記為修飾的 C」且通過 ref-CpG 驗證者。因此：
- read 未覆蓋的 CpG → 不在輸出中（**parser 層無 NA sentinel**）。
- 覆蓋但落在 insertion → 丟棄（`ref_pos_0based < 0`，`:172`/`:134`）。
- 覆蓋但 ref 非 CpG context → 丟棄（F5 驗證不過）。
- **NA / 缺值的填補語意由下游決定**（`ReadAggregator` / `MatrixBuilder`），不在本檔；本檔僅產 sparse 的 `vector<MethylCall>`。**[uncertain]** matrix 層如何填 NA 不在本稽核範圍。

### F9 — MM/ML 長度一致性與越界保護
- 要求 `ml_offset + deltas.size() <= ml_len`，否則回空：`MethylationParser.cpp:94-96`：
  ```cpp
  if (static_cast<size_t>(ml_offset + deltas.size()) > ml_len) { return calls; }
  ```
- delta 為空（無 5mC 修飾）回空：`MethylationParser.cpp:89-91`
- `parse_mm_tag` 用 C++17 `std::from_chars` 解析逗號分隔整數，遇 `;` 或 parse error 停止：`MethylationParser.cpp:220-227`；要求 mod_code 後緊接 `,`：`:211-213`。
- delta 解碼：`next_target += deltas[delta_idx] + 1`（delta = 兩個修飾間「未修飾目標鹼基」數）：forward `:187-191`、reverse `:149-153`；header 註解 `MethylationParser.hpp:73-74`。

---

## L2 關鍵數字 / 常數 / 閾值（附來源）

| 項目 | 值 | 來源 |
|------|----|------|
| ML byte→機率除數 | `255.0f`（線性） | `MethylationParser.cpp:140, 180`；`.hpp:39` |
| 目標 modification code | `"C+m?"`（5mC，hard-coded） | `MethylationParser.cpp:60`；`parse_mm_tag` 預設 `.hpp:76` |
| target code 長度 | `tlen = 4` | `MethylationParser.cpp:61` |
| ML 型別檢核 | `'B'` then `'C'`（array of uint8） | `MethylationParser.cpp:31` |
| ML 長度欄位 offset | `+2`（4-byte LE uint32） | `MethylationParser.cpp:35` |
| ML 資料 offset | `+6` | `MethylationParser.cpp:36` |
| forward 回報座標 | `ref_pos_0based + 1`（1-based C 位） | `MethylationParser.cpp:181` |
| reverse 回報座標 | `ref_pos_0based`（C 位 1-based） | `MethylationParser.cpp:142-143` |
| delta reserve hint | `64` | `MethylationParser.cpp:217` |
| 機率輸出範圍 | `[0.0, 1.0]` | `MethylationParser.hpp:16` |
| 無顯式機率閾值/二值化 | （無） | 全檔無 threshold 常數 |

---

## L3 讀不確定 / 邊界註記（uncertain）

1. **[uncertain]** 主路徑命中 target block 後呼叫 `parse_mm_tag(mm_str, "C+m?")`（`:66`），該 helper 用 `strstr` 從 **MM 字串開頭**再找一次 `"C+m?"`（`:207`）。若 MM 中存在多個 `C+m?` 或 `C+m?` 為其他 code 子字串，可能有歧義；但標準 ONT MM tag 不會出現此情形，正常輸入無影響。屬重複掃描的微冗餘，非正確性 bug（依標準輸入）。
2. **[uncertain]** reverse 分支 CpG 驗證用 `ref_offset > 0`（嚴格大於 0），forward 用 `ref_offset >= 0`。reverse 需檢查 `ref_offset - 1`，故 `> 0` 合理（避免越界），但意味 ref_seq 最左端的 reverse CpG 邊界鹼基可能漏收 — 屬區域邊界 1-bp 效應，需下游 `ref_start_pos` 是否有 padding 才能判定影響。
3. **[uncertain]** 5hmC 被「offset 跳過」依賴 MM block 順序為 `C+h?` 在前、`C+m?` 在後（`:41` 註解假設），且非 target block 必含 `,` 才推進 offset（`:75`）。若某 mod block 無 delta（無逗號），`ml_offset` 不推進 — 對 5hmC 通常有 delta，正常無影響，但極端空 block 情形未測。
4. **[uncertain]** 「未覆蓋 CpG → NA」的最終呈現由下游 matrix builder 決定，不在本兩檔內，無法在此確認 NA sentinel 值。
5. 主路徑用 inline CpG 檢查（`:176-178`/`:138-139`），**未呼叫** `is_cpg_site()`（`:283`）；後者僅 forward 語意（`C` then `G`），與 reverse inline 邏輯不同 — 為 dead-ish helper，不影響輸出。

---

## 對「論文方法描述」的直接可用結論

- ISM 甲基化定量 = **ONT 5mC 機率（`ML/255`，linear）**，**不含 5hmC**（5hmC 在 ML 陣列中被 offset 略過，不 collapse、不雙列、不取 max）。
- 每條 read 產出 sparse 的 per-CpG 5mC 機率，CpG 以 **reference 序列**驗證、座標回報 1-based、strand-aware（reverse 反向迭代 + G 目標）。
- 未覆蓋 / insertion / 非 CpG-context 的位點**不進輸出**，NA 語意在下游 matrix 層補（本層無 NA sentinel）。
