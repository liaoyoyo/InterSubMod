# ISM 甲基數據真值稽核 — methylation.csv / cpg_sites.tsv 格式驗證

> 稽核日期 2026-06-13。只讀實際檔案 + C++ 原始碼，引行號佐證。所有數字附來源檔:行。

## L0 一眼結論

`methylation.csv` = **read × CpG 連續甲基化機率矩陣**（值域 [0,1]，granularity 1/255 ≈ 0.0039），
**非二元 call**。`NA` = 該 read 在該 CpG 無覆蓋（C++ 內部 sentinel `-1.0`）。
矩陣維度與 metadata.txt 宣告完全一致（flagship 57×203 / 比較檔 32×66）。
值來源 = BAM `ML` tag 的 8-bit 機率 `ml_data[...] / 255.0f`（`MethylationParser.cpp:140`）。

---

## 1. 檔案盤點（兩個 run）

| 檔案 | flagship_chr2_18086020 | 比較檔 V6_on_fp chr2_18096269 |
|---|---|---|
| methylation.csv | 57 reads × 203 CpG（204 欄含 read_id；58 行含 header） | 32 reads × 66 CpG（67 欄；33 行含 header） |
| cpg_sites.tsv | 204 行（header + 203 CpG，cpg_id 0–202） | 67 行（header + 66 CpG） |
| 同層伴隨檔 | methylation_forward.csv / methylation_reverse.csv（±strand 拆分） | 同左 |

flagship 路徑：`research/flagship_chr2_18086020_20260612/ism_out/anchor_18086020/chr2/chr2_18086020/chr2_18071020_18101020/methylation/`
比較路徑：`research/paired_priority_bug_audit/phaseC_genome_three_way_with_significance/V6_on_fp/filtered_snv_fp/chr2/chr2_18096269/chr2_18091269_18101269/methylation/`

---

## 2. methylation.csv 格式（header + 值語意）

### Header 結構
- 第一欄 = `read_id`（整數 0..N-1，**非 read name**；read name 在 reads/ 子目錄的 metadata，非此檔）。
  - flagship: `read_id` 範圍 0..56（57 reads）。比較檔 0..31。
- 其餘欄 = **genomic CpG 座標**（hg38 絕對位置，非 cpg_id），逐欄遞增。
  - flagship header[1]=`18071239`, header[-1]=`18100977`（與 cpg_sites.tsv position 欄一字不差對齊）。
  - 比較檔 header[1]=`18091479`, header[-1]=`18100977`。
- 源碼佐證：`src/io/RegionWriter.cpp:222-225` —
  `ofs << "read_id"; for (auto pos : cpg_positions) ofs << "," << pos;` → header 是 read_id 後接 **genomic position**（不是 cpg_id）。

### 值範圍（實測，全 numeric cell）
- flagship: `min=0.0`, `max=1.0`, **0 個值落在 [0,1] 之外**。
- 比較檔: `min=0.0039`, `max=1.0`, 同樣全在 [0,1]。
- granularity = **1/255**：實測最小非零相鄰步進 = 0.0039（= 1/255 = 0.00392…），distinct 值集合 `[0.0, 0.0039, 0.0078, 0.0118, …, 0.9961, 1.0]`（flagship 240 個 distinct 值）→ 8-bit 量化機率。
- 印出精度 = `setprecision(4)` 定點 4 位小數（`RegionWriter.cpp:236`）。

### NA 語意
- `NA` = 該 read 在該 CpG **無覆蓋 / 無 ML call**，非 0。
- 源碼佐證：`RegionWriter.cpp:233-234` — `if (matrix[r][c] < 0.0) ofs << "NA";`。矩陣初值 `-1.0`（`MatrixBuilder.cpp:79-82`：`// Initialize with -1.0 (NaN/No Coverage)` … `matrix_[r].resize(num_cols, -1.0);`）。
- flagship NA = 5807 / 11571 cells（50.2%），numeric = 5764；屬正常稀疏（read 只覆蓋自身對齊跨度內的 CpG，矩陣右上/左下大量 NA）。
- 比較檔 NA = 478 / 2112 cells，numeric = 1634。

---

## 3. 「連續機率」語意 — 3 個 read 實值佐證（flagship）

非二元的最直接證據：大量中間值（既非 ~0 也非 ~1）。實測 flagship 全矩陣分布：
- frac ≤0.05（unmethylated-like）= 0.254
- frac ≥0.95（methylated-like）= 0.494
- **frac 嚴格 0.05<v<0.95（中間機率）= 0.252** ← 1/4 cell 是中間機率，二元 call 不可能出現

逐 read 中間機率實例（methylation.csv 行；flagship）：
- **read_id=0**（檔行 2）：pos `18074603`=**0.5882**、pos `18081566`=**0.6627**、pos `18076408`=**0.7137**、pos `18072863`=**0.8941**、pos `18074539`=0.1059。
- **read_id=1**（檔行 3）：pos `18074539`=**0.5882**、pos `18072863`=**0.6784**、pos `18074199`=**0.7137**、pos `18074321`=**0.8824**、pos `18074603`=0.1373。
- **read_id=2**（檔行 4）：pos `18074321`=**0.3569**、pos `18074478`=0.1490、pos `18073987`=**0.7451**、pos `18074144`=**0.7725**、pos `18071954`=0.1098。

這些 0.36 / 0.59 / 0.66 / 0.71 等值 = per-read per-CpG 的**機率**（modkit/ML 的 likelihood），不是 0/1 hard call。

---

## 4. cpg_sites.tsv 欄位

- 分隔 = tab；3 欄：`cpg_id\tchr\tposition`（`RegionWriter.cpp:207`）。
- 每列 = `i \t chr_name \t cpg_positions[i]`（`RegionWriter.cpp:210-211`），cpg_id 從 0 連續遞增。
- flagship: cpg_id 0..202（203 CpG），position 18071239..18100977；末 3 列 `200 chr2 18100682 / 201 chr2 18100819 / 202 chr2 18100977`。
- **cpg_id ↔ methylation.csv 欄對應**：methylation.csv 第 (k+1) 個資料欄的 genomic position = cpg_sites.tsv `cpg_id=k` 的 position（兩檔 position 序列逐一相等，已實測對齊）。即 cpg_sites.tsv 是 methylation.csv 欄頭的「id↔座標」對照表。

---

## 5. metadata.txt 對齊（flagship 區域目錄 metadata.txt）

```
Num Reads: 57   (Forward + : 32, Reverse - : 25)
Num CpG Sites: 203
Matrix Dimensions: 57 × 203
SNV: G -> A @ chr2:18086020   Region: chr2:18071020-18101020 (30001 bp)
```
→ 與 methylation.csv 實測 57 reads × 203 CpG **完全一致**；32+25=57 對齊 forward/reverse 拆分檔。

---

## 6. 數據真值來源鏈（C++ 原始碼）

1. **ML tag → 機率**：`src/core/MethylationParser.cpp:140`（亦 :180）
   `float prob = ml_data[ml_offset + delta_idx] / 255.0f;`
   → BAM `ML` tag 為 0–255 的 8-bit 機率，除以 255 得 [0,1] 連續機率（解釋了實測 1/255 granularity）。
2. **MethylCall 結構**：`include/core/MethylationParser.hpp:14-22`
   `int32_t ref_pos;` + `float probability;  ///< Methylation probability [0.0, 1.0]`。
3. **存入矩陣**：`src/core/MatrixBuilder.cpp:31-35`（`calls.emplace_back(call.ref_pos, call.probability)`）、:79-95（矩陣初值 -1.0，填入 `prob`）。
4. **寫出 CSV**：`src/io/RegionWriter.cpp:217-243`（header read_id+positions；`<0.0`→`NA`；else `setprecision(4)`）。
5. **CpG 判定**：`MethylationParser.cpp:137-139` 要求 ref 為 `C`(ref_offset-1) + `G`(ref_offset) 才算 CpG（context-validated）。
6. **MM/ML 必要性**：`ReadParser.cpp:55-63` — `require_mm_ml` 時無 MM/ML tag 的 read 被 filter（`MISSING_MM_TAG` / `MISSING_ML_TAG`）。

---

## 7. 不確定 / 邊界

- header 的 CpG 欄頭 = genomic position（不是 cpg_id）；若下游程式誤把欄頭當 0-based index 會錯位（兩者不同：position 是絕對座標）。
- `MethylCall.ref_pos` 註解寫 "1-based genomic position"（hpp:15），但 `MatrixBuilder` 直接用 `ref_pos_0based`（MethylationParser.cpp:132,142 變數名 `ref_pos_0based` 傳入 emplace）→ **1-based vs 0-based 命名在 struct 註解與 parser 變數名之間不一致**，未實測下游座標是否確實偏移 1（標 uncertain；不影響本次「值語意/維度」結論，但若做 CpG 位置精確比對需查證）。
- 比較檔 region 跨度 18091269-18101269 與 flagship 18071020-18101020 部分重疊，故末欄 position 相同（18100977）；屬不同 run/不同 SNV anchor，非同一矩陣子集。
- forward/reverse 拆分檔（methylation_forward/reverse.csv）多一欄 `original_read_id`（`RegionWriter.cpp:265,294`），與主檔 schema 略不同；本次未逐值核對拆分檔。
