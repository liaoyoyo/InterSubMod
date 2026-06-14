<!--
建立時間: 2026-06-13
狀態: verified_mechanism_doc
data_sources:
  # 源碼檔（C++ binary）
  - src/core/MethylationParser.cpp
  - include/core/MethylationParser.hpp
  - src/core/MatrixBuilder.cpp
  - include/core/MatrixBuilder.hpp
  - include/core/MethylationMatrix.hpp
  - src/core/RegionProcessor.cpp
  - include/core/RegionProcessor.hpp
  - src/io/RegionWriter.cpp
  - src/core/PerCpgAsm.cpp
  - include/core/PerCpgAsm.hpp
  - src/core/DistanceMatrix.cpp
  - src/core/ReadParser.cpp
  - include/core/ReadParser.hpp
  - include/core/Config.hpp
  - src/core/Config.cpp
  - include/utils/ArgParser.hpp
  - src/core/ReadAggregator.cpp
  - tests/test_methylation_parser.cpp
  - scripts/run_vcf_all_snv.sh
  - scripts/fill_report.py
  # 數據檔（實測真值）
  - research/flagship_chr2_18086020_20260612/ism_out/anchor_18086020/chr2/chr2_18086020/chr2_18071020_18101020/methylation/methylation.csv
  - research/flagship_chr2_18086020_20260612/.../methylation/cpg_sites.tsv
  - research/flagship_chr2_18086020_20260612/.../metadata.txt
  # 稽核 doc
  - docs/methodology/20260610_code_methodology_detail_audit_01.md
provenance_note: >
  本檔由 ism-methylation-reading-verify workflow 產出。每項細節皆附「檔:行」溯源，
  並經獨立對抗驗證（extract 階段讀源碼 → verify 階段預設懷疑、回原碼複核）。
  5 條核心 claim 各有對應 verify/*.md 裁決（SUPPORTED / PARTIAL / REFUTED / OUTDATED）。
  數據真值來自本日 Python 全矩陣掃描 flagship methylation.csv（非記憶/預期）。
  中間檔: docs/paper_focus/_ism_methyl_verify_20260613/{extract,verify}/。
-->

# ISM 甲基讀取 — 機制細節與數字（verified 校正版）

> 紀律：本檔只整理 extract/verify 的源碼+數據佐證結果，不新增未經佐證的數字；有口徑衝突誠實並列。
> 標記：🔴 最該校正/高重要 ・ ⭐ 口徑 caveat ・ ◽ 細節 ・ [uncertain] 讀不確定。

---

## L0 一眼

ISM 甲基讀取 = **直接讀 BAM `MM/ML` tag → 解碼成 read×CpG 的「連續甲基化機率」矩陣（值 ∈ [0,1]，granularity 1/255）→ 再餵距離計算 / 分群 / per-CpG ASM**。每個 cell 是一條 read 在一個 CpG 位點的 ONT 機率，不是 0/1 hard call；未覆蓋的格子 = NA（內部 sentinel `-1.0`）。

**最該校正的 1-2 點**：

1. 🔴 **ISM C++ 是 5mC-only，不是「5mC+5hmC max-collapse / 雙列砍半」**。C++ 解析只取 `C+m?`（5mC），5hmC（`C+h?`）的機率被 `ml_offset` 純粹**跳過、從不產生 call**。所謂「max-collapse / 雙列砍半 dup-bug」是**外部 Python/MSA（MethylSomaticAnalysis）抽取工具**的 Level-1 行為，**不屬本 C++ binary**（`MethylationParser.cpp:60-82`；verify `5mc_only`=SUPPORTED、`dup_bug`=REFUTED）。
2. 🔴 **methylation.csv 存的是連續機率 P(modified)，不是 P(5mC) 也不是二值**。把它標成「P(5mC)」是 over-specification（源碼層只能斷定是「一個甲基化機率」）；二值化（閾值 0.8/0.2）發生在**下游 `binary_matrix`**，**不寫進 methylation.csv**（verify `matrix_semantics`=PARTIAL、`threshold_use`=PARTIAL）。

---

## L1 機制流程（verified，逐步附 檔:行）

資料流主鏈（caller：`src/core/ReadAggregator.cpp:62` → `methyl_parser.parse_read(b, ref_seq, region_start)`）：

```
BAM MM/ML tag
   │  MethylationParser::parse_read           [MethylationParser.cpp:15-199]
   ▼
vector<MethylCall>{ref_pos(1-based), probability∈[0,1]}   ← sparse，只含「被標記的 CpG」
   │  MatrixBuilder::add_read / finalize      [MatrixBuilder.cpp:17-104]
   ▼
matrix_ : vector<vector<double>>  連續機率，未覆蓋=-1.0
   ├─► RegionWriter::write_matrix_csv  → methylation.csv（連續，<0→"NA"）  [RegionWriter.cpp:217-243]
   └─► RegionProcessor::build_methylation_matrix → MethylationMatrix{raw(NaN), binary(1/0/-1)}  ← 二值化在此  [RegionProcessor.cpp:1381-1429]
                                                        │
                          DistanceMatrix（依 metric 吃 raw 或 binary） / PerCpgAsm（Fisher 吃 raw>0.5）
```

### ① MM/ML 解析（含 5mC / 5hmC 處理）

- 用 htslib `bam_aux_get("MM")` / `bam_aux_get("ML")` 取 tag，缺任一即回空（`MethylationParser.cpp:20-25`）。ML 只接受型別 `B`（array）+ `C`（uint8），否則回空（`MethylationParser.cpp:31-33`）；ML 長度為 offset `+2` 起的 4-byte LE uint32、資料指標在 offset `+6`（`MethylationParser.cpp:35-36`）。
- 掃描 MM tag 的 `;`-分隔 block，**唯一硬編碼 target = `"C+m?"`（5mC，`tlen=4`）**（`MethylationParser.cpp:60-63`）；命中才呼叫 `parse_mm_tag(mm_str, "C+m?")` 並 break（`:65-69`），找不到 `C+m?` 整條 read 回空（`:85-87`）。
- 🔴 **5hmC（`C+h?`）只用來推進 `ml_offset`、從不 emit call**：非 target block 數其 delta 個數（`delta_count++; ml_offset += delta_count;`，`:71-82`），藉此把 ML 陣列中 5hmC 那段機率「跳過」。註解明列 ML 順序為 `[C+h? probs..., C+m? probs...]`（`:41`）。`C+h?`/`5hmC` 在全 parser+matrix 程式碼**僅出現於註解**（`MethylationParser.cpp:39,41,104`、`.hpp:33`），無任何 active 解析路徑（verify `dup_bug` E1 grep 佐證）。

### ② ML 機率轉換

- forward 與 reverse 兩路皆同公式：`float prob = ml_data[ml_offset + delta_idx] / 255.0f;`（forward `:180`、reverse `:140`；header `.hpp:39` `Probability = ML[i] / 255.0`）。
- **線性換算，無 0.5/255 居中偏移、無 thresholding、無 binarization、無 clip**（byte 0→0.0，byte 255→1.0）；輸出範圍文件標 `[0.0, 1.0]`（`.hpp:16`）。

### ③ CpG 對位 + strand

- delta 命中後用 `build_seq_to_ref_map`（走 CIGAR）把 read index 映回 ref 0-based 座標（`MethylationParser.cpp:232-281`；I/S 映 -1 故落 insertion/soft-clip 的修飾 C 被丟棄，`:172`/`:134` 的 `>=0` 守衛）。
- **用 reference 序列驗證 CpG dinucleotide 才收錄**（非 read SEQ，故 read 上 mismatch/變異不會被誤判）：forward `ref_seq[ref_offset]=='C' && [ref_offset+1]=='G'`（`:176-178`）；reverse `ref_seq[ref_offset]=='G' && [ref_offset-1]=='C'`（`:138-139`）。
- **strand-aware**：`bool is_reverse = bam_is_rev(b);`（`:114`）。MM delta 以原始 read 5'→3' 排序，BAM SEQ 對 reverse 是 RevComp，故 reverse **從 `seq_len-1` 倒著掃**、target 鹼基切為 `'G'`（`:124-126`，碼內自標 CRITICAL `:103-110`）；forward 正向掃、target `'C'`（`:162-164`）。
- 座標回報 1-based：forward `ref_pos_0based + 1`（`:181`）、reverse `ref_pos_0based`（C 位 1-based，`:142-143`）。reverse 與 forward 的同一 CpG 會聚到**同一座標欄**（two-strand collapse，test TC2 已驗）。

### ④ 矩陣值語意 + NA

- `MatrixBuilder` 的 `matrix_`（`vector<vector<double>>`）**直接存連續機率、無二值化**：`add_read` 存 `(ref_pos, probability)` pair（`MatrixBuilder.cpp:34-36`），`finalize` 填 `row[col] = prob`（`:89-96`）。
- **Row = read（讀入順序）；Col = unique CpG genomic position（sort+unique 去重、座標排序）**（`MatrixBuilder.hpp:22-23`；`MatrixBuilder.cpp:62-66, 75-83`）。column 軸 = 該 region 內被任一 read 覆蓋過的 CpG 聯集，非固定參考 CpG list。
- **NA = 未覆蓋**：矩陣先全填 `-1.0`（`MatrixBuilder.cpp:82`），只有「該 read 在該 pos 有 MethylCall」才覆寫成機率；沒 call 的格保持 `-1.0`。CSV 以 `matrix[r][c] < 0.0` 寫成字串 `"NA"`、否則 4 位小數機率（`RegionWriter.cpp:233-237`）。

### ⑤ 閾值 / 二值化何時用

- **唯一二值化點在下游** `RegionProcessor::build_methylation_matrix`（`RegionProcessor.cpp:1410-1423`），由 `--methyl-high`/`--methyl-low`（預設 0.8 / 0.2）切三態：`val ≥ high → 1`、`val ≤ low → 0`、`low < val < high → -1`（Ambiguous）、`val < 0 → raw=NaN, binary=-1`。
- ⭐ **二值化不發生在 methylation.csv**：CSV 走連續 raw（verify `threshold_use` 確認）。閾值的下游消費者只有 NHD/JACCARD 距離 + binary_matrix-based 度量。
- ⭐ **另一道二值化口徑不一致**：`PerCpgAsm` 的 per-CpG Fisher 用 `raw > 0.5` 硬切 meth/unmeth（`PerCpgAsm.cpp:282,283,289,290`），與 binary_matrix 的 0.2/0.8 灰帶**口徑不同**（稽核 ASM-4）。

### ⑥ 距離吃連續 or 二值

- **依 metric 而定，非無條件連續**（`DistanceMatrix.cpp:309-345`）：`NHD`/`JACCARD` 吃 `binary_matrix`（受 0.8/0.2 閾值影響）；`L1`/`L2`/`CORR`/`BERNOULLI` 吃 `raw_matrix`（連續，不受閾值影響）。
- 🔴 **生產實況用 `BERNOULLI`（連續）** — 生產腳本顯式傳 `--distance-metric BERNOULLI`（`scripts/run_vcf_all_snv.sh:65`；`run_pure_research_round.sh:17`）。但 **ISM binary 的 CLI 預設（不傳 `--distance-metric`）卻是 NHD（二值）**（`ArgParser.hpp:86,163-176` 無條件清空 Config.hpp:40 的 BERNOULLI 預設）。故「ISM 用連續距離」是**生產配置選的，不是 binary 內建保證**。

---

## L2 關鍵數字 / 常數表（每個附來源）

| 項目 | 值 | 來源 檔:行 |
|------|----|-----------|
| ML byte→機率換算 | `ML[i] / 255.0f`（線性，無偏移/分箱/clip） | `MethylationParser.cpp:140, 180`；`.hpp:39` |
| 機率量化粒度（granularity） | 1/255 ≈ 0.0039（8-bit） | 數據實測：flagship min positive step 0.003900 |
| 矩陣值範圍 | `[0.0, 1.0]`；未覆蓋 = `-1.0`（連續層 sentinel） | `MatrixBuilder.hpp:25-28`；`MatrixBuilder.cpp:82` |
| methylation.csv 印出精度 | 定點 4 位小數（`setprecision(4)`） | `RegionWriter.cpp:236` |
| 目標 modification code | `"C+m?"`（5mC，hard-coded；`tlen=4`） | `MethylationParser.cpp:60-61`；`parse_mm_tag` 預設 `.hpp:76` |
| methyl-high（二值 1 閾值） | 預設 **0.8**（`>=` 含端點） | `Config.hpp:33`；`RegionProcessor.cpp:1417` |
| methyl-low（二值 0 閾值） | 預設 **0.2**（`<=` 含端點） | `Config.hpp:34`；`RegionProcessor.cpp:1419` |
| 二值中間態（Ambiguous） | `low < val < high → -1`（與未覆蓋共用 -1 編碼） | `RegionProcessor.cpp:1422` |
| 閾值約束 | `high > low`，兩者 ∈ [0,1] | `Config.cpp:108-116` |
| Fisher per-CpG 二值化閾值 | `raw > 0.5`（與 0.2/0.8 口徑不同） | `PerCpgAsm.cpp:282,283,289,290` |
| Fisher 顯著閾值（BH-FDR） | FDR `< 0.05` | `PerCpgAsm.cpp:326` |
| min_reads_per_group（ASM 分組門檻） | 預設 **5** | `PerCpgAsm.hpp:66` |
| min-mapq（read 級過濾） | 預設 **20**（range 0-60） | `Config.hpp:29`；`ReadParser.cpp:45`；`ArgParser.hpp:67-68` |
| min-read-length（read 級過濾） | 預設 **1000** bp | `Config.hpp:30`；`ReadParser.cpp:50-52` |
| min-base-quality | 預設 **20**（range 0-93）— **只 gate SNV 位點，不過濾 CpG/甲基讀值** | `Config.hpp:31`；`ReadParser.cpp:245` |
| min_common_coverage（C_min，距離最少共同 CpG） | 預設 **3**（實際生效） | `Config.hpp:37`；`DistanceMatrix.cpp:313` |
| 維度規則 | Row=read（讀入序）；Col=unique CpG genomic pos（sort+unique） | `MatrixBuilder.hpp:22-23`；`MatrixBuilder.cpp:62-66` |
| 生產距離 metric | `BERNOULLI`（連續 raw） | `run_vcf_all_snv.sh:65` |
| CLI 預設距離 metric（不傳時） | `NHD`（二值）— 覆蓋 Config.hpp:40 的 BERNOULLI | `ArgParser.hpp:86,163-176` |
| 無顯式機率閾值/二值化於 parser | （parser 全檔無 threshold 常數） | `MethylationParser.cpp` 全檔 |

**實測數據真值（flagship methylation.csv，2026-06-13 Python 全矩陣掃描）**：

| 量 | 實測值 | 來源 |
|----|--------|------|
| 維度 | **57 reads × 203 CpG** | methylation.csv + metadata.txt `Matrix Dimensions: 57 × 203` |
| 值範圍 | min=0.0, max=1.0，out-of-[0,1] = 0 | flagship 全矩陣掃描 |
| distinct 值數 | **240**；min positive step 0.0039（≈1/255） | flagship 掃描 |
| frac 嚴格中間 (0.05<v<0.95) | **0.252**（1/4 cell 是中間機率 → 二值 call 不可能出現） | flagship 掃描 |
| NA 比例 | **50.2%**（5807 / 11571 cells，正常稀疏） | flagship 掃描 |

---

## L2 校正既有理解（對抗驗證結果）

下表為 5 條核心 claim 的 verify 裁決（每條對應 `verify/*.md`），外加 extract 額外發現的口徑差異。

| # | 舊理解 | verdict | 校正後 | 原碼證據 |
|---|--------|---------|--------|----------|
| 1 | ISM C++ 用「5mC+5hmC max-collapse / 雙列砍半」 | **SUPPORTED（核心）+ 措辭校正** | ISM C++ **5mC-only**：只解析 `C+m?`，5hmC 機率被 `ml_offset` 純粹**跳過、不 emit call**。不是 collapse、不是 merge、不是雙列。max-collapse/雙列是外部 MSA（Python）行為。 | `MethylationParser.cpp:60-69`（target=`C+m?`）、`:71-82`（5hmC 只推 offset）、`:85-87`（找不到即回空）；verify `5mc_only` |
| 2 | methylation.csv 值 = 「P(5mC)」、可能是二值 | **PARTIAL** | 值 = **連續甲基化機率 P(modified) ∈ [0,1]**（非二值，實測 240 distinct / 25.2% 中間值）。標成「P(**5mC**)」是 over-specification — 源碼層只能斷定是「一個甲基化機率」。per-read×per-CpG、NA=未覆蓋皆成立。 | `MethylationParser.cpp:140,180`（ML/255）；`MatrixBuilder.cpp:34-35,91-95`（無二值化）；`RegionWriter.cpp:229-239`；verify `matrix_semantics` |
| 3 | MM/ML 解析「零測試」（舊稽核説法） | **OUTDATED** | 06-10 早稽核當下為真（HEAD `bbcc7a9` 測試檔 ABSENT）；同日稍晚 commit `6593f96`（13:51）補上 `tests/test_methylation_parser.cpp`（3 golden TC，已 wire `CMakeLists.txt:169`）。正確表述 = 「**已修正的歷史缺口（happy-path only）**」。 | `tests/test_methylation_parser.cpp:66,97,139`；`git cat-file -e bbcc7a9:...` = ABSENT；verify `mmml_tested` |
| 4 | 5mC+5hmC「雙列砍半 dup-bug」在 ISM C++ | **REFUTED** | dup-bug **不存在於本 ISM C++**：parser 5mC-only 單通道、每 CpG 單機率，`MatrixBuilder` 以 position 為唯一鍵去重（1 pos = 1 col）、填值是覆寫非平均/砍半 → 結構上無從雙列。dup-bug 屬外部 MSA。 | `MatrixBuilder.cpp:55-72`（去重）、`:89-96`（覆寫）；`ReadAggregator.cpp:62`（唯一生產者）；verify `dup_bug` |
| 5 | 距離矩陣「無條件用連續、閾值只在某些輸出二值化」 | **PARTIAL** | 距離**依 metric 而定**：NHD/JACCARD 吃二值、L1/L2/CORR/BERNOULLI 吃連續。生產用 BERNOULLI（連續）→ 結論成立；但 **CLI 不傳 `--distance-metric` 時預設落 NHD（二值）**，故「連續」是生產配置選的非 binary 內建。閾值唯一寫進 binary_matrix（後半成立）。 | `DistanceMatrix.cpp:309-345`；`run_vcf_all_snv.sh:65`；`ArgParser.hpp:86,163-176`；verify `threshold_use` |

**extract 額外發現的口徑差異 / 細節**：

- ⭐ **binary_matrix 的 `-1` 有歧義**：同時表「未覆蓋」與「中間模糊 (low<val<high)」兩種情形（`RegionProcessor.cpp:1413` vs `1422`）。下游若要區分「真未覆蓋」vs「覆蓋但模糊」需回 `raw_matrix`（NaN vs 非NaN）判斷。header 註解只寫 "-1 (missing)" 未提 Ambiguous 也併入。
- ⭐ **per-CpG Fisher 二值化 `raw>0.5` vs binary_matrix 0.2/0.8 口徑不一**（稽核 ASM-4，◽ caveat，非 bug）。
- ◽ **1-based vs 0-based 命名不一致**：`MethylCall.ref_pos` 註解寫 "1-based"（`.hpp:15`），forward `emplace_back(ref_pos_0based + 1, ...)` 對齊 1-based，reverse `emplace_back(ref_pos_0based, ...)`（註解稱 `-1+1`）。不影響值語意/維度/NA 結論；CpG 位置精確比對才需查證 [uncertain]。
- ◽ **min-base-quality comment 與實作不符**：comment 寫 "for SNV/CpG sites"，實作只 gate SNV 位點 ALT/REF 支持，**不過濾 CpG 甲基讀值**（`ReadParser.cpp:245,249-259`）。
- ◽ **methylation.csv 第一欄 = row index（整數 0..N-1），非 read name**（`RegionWriter.cpp:230`）；header 其餘欄 = genomic position（非 cpg_id），`cpg_sites.tsv` 為 id↔座標對照（`RegionWriter.cpp:222-225`；data extract §4）。

---

## L2 既有稽核項狀態（2026-06-10 甲基 finding 逐條現狀）

對象稽核：`docs/methodology/20260610_code_methodology_detail_audit_01.md`（12-reviewer workflow，稽核 HEAD `bbcc7a9`）。下表為 6 條甲基相關 finding + 1 引用項的現行 repo 狀態。

| Finding | 稽核分級 | 現狀 | 證據 |
|---------|---------|------|------|
| **TESTS-1** MM/ML 解析最易錯碼零測試 | 🔴 必修 CONFIRMED HIGH | **已修（POST-AUDIT）** | commit `6593f96`（06-10 13:51，稽核 04:50 之後）新增 `test_methylation_parser.cpp` 3 golden TC；wired `CMakeLists.txt:169`。⚠ 僅 happy-path（見 L3）。 |
| **PYTOOL-1** fill_report.py 對 null/NaN 失效（反捏造洞） | 🔴 必修 CONFIRMED HIGH | **已修（POST-AUDIT）** | commit `eed4300`（06-10 13:51）新增 `is_nullish()`，`fill_report.py:94` `if not found or is_nullish(val): missing.append(key)`。 |
| **引用 De Waele→Orjuela** | ✅ 已修（附錄 A） | **已修（PRE-AUDIT）** | commit `891e04b`（06-10 01:11，稽核前）；`PerCpgAsm.cpp:6` 現 = `DAMEfinder (Orjuela 2020)`；`grep "De Waele" src/ include/` = 0 命中。舊 memory 標的「PerCpgAsm.cpp:6 誤引」屬此修之前的舊紀錄、現已過時。 |
| **FISHER-1** per-CpG Fisher 忽略 read over-dispersion | ⭐ CONFIRMED→MEDIUM | **仍在（未改）** | `PerCpgAsm.cpp:279-305` 每 read×CpG binary call 當獨立 Bernoulli 餵 Fisher，無 beta-binomial/quasi/GEE。屬 Hard Gate C++ 改動（需 /methodology-audit→/cpp-change）。⚠ 稽核的 53–68% FP 膨脹數字**自標 ILLUSTRATIVE（合成示意非真資料）**，真量級 needs_rerun — 引用必帶此 caveat。 |
| **ASM-2 / REGION-5** ASM 軸別 = germline/ALLELE 軸（1-1 折進 family 1） | ⭐ UNCERTAIN→MEDIUM | **仍在（口徑治理，未改）** | `PerCpgAsm.cpp:244-250` 把 `{1,HP1,1-1}→0`、`{2,HP2,2-1}→1` 合併 → 量的是 HP1-family vs HP2-family（germline/ALLELE 軸，會被 baseline allelic methylation confound），非 somatic-controlled HP1-vs-HP1-1。輸出欄名 `tumor_hp_signed_delta` 未明示軸別 → 下游易誤用。與 memory `feedback_asm_allele_axis_baseline_confound` 一致。 |
| **ASM-3** C++ 為 5mC-only | ◽ caveat | **仍在（非 bug）** | `MethylationParser.cpp:60`（target=`C+m?`）、`:71-82`（5hmC 只推 offset）。稽核明確：C++ 層無「雙列砍半」bug，真議題是 C++ **忽略 5hmC**（與論文要的 any-mod max-collapse 有口徑差）。 |
| **ASM-4** Fisher `raw>0.5` 與 binary_matrix 0.2/0.8 口徑不一 | ◽ caveat | **仍在（非 bug）** | `PerCpgAsm.cpp:282` `if (val > 0.5)`，無 0.2~0.8 灰帶（NaN 在二值化前 skip，`:281`）。 |

> 稽核總結（附錄）：兩個 🔴 必修均在稽核**當日稍晚**落地修復；引用 bug 在稽核**前**已修。仍未動的甲基段議題（FISHER-1 / ASM-2 / ASM-3 / ASM-4）皆口徑治理 / 統計校準，稽核明示**不翻轉任何已 concluded 的 NEGATIVE 結論方向**。

---

## L3 邊界 / 未測（誠實列出）

**測試覆蓋的 gap（補測 commit `6593f96` 僅 happy-path，以下仍零測試）**：

1. 🔴 **5hmC `C+h?` 在前的 `ml_offset` 前進**：3 TC 的 MM 全是單一 `C+m?` block（`test:75,104,114`），`grep "C+h" tests/` = 0 命中 → 源碼 `MethylationParser.cpp:71-82`（5hmC offset 跳過）**未被測**。此即「5mC/5hmC max-collapse」相關路徑。
2. **非零 delta（skip-count>0）**：test 所有 MM delta 皆 `0`，未測 delta≥1（C 間有 skip）的解碼正確性（源碼 `:201-230` / `next_target += deltas[delta_idx]+1` `:150,188`）。
3. **非全-match CIGAR**：`make_meth_bam` 只造 `BAM_CMATCH`；`build_seq_to_ref_map`（`:232-281`）的 I/D/S/N/H 分支（含 soft-clip 偏移、insertion 不推 ref）**完全未測**。
4. **格式/邊界防呆 early-return**：`ml_aux[0]!='B'||[1]!='C'`（`:31-33`）、`ml_offset+deltas.size()>ml_len`（`:94`）、`deltas.empty()`（`:89`）、`found_target=false`（`:85`）皆未測。
5. **`ref_pos_0based<0` 守衛**（`:134,172`）未測。
6. **未實際編譯/執行測試**：verify 為 read-only，未 `make`+`ctest`；「3 TC 全 PASS」**未經 runtime 證實**（只證 source 邏輯與 spec 一致）。要 SUPPORTED-with-green 需跑 `ctest -R MethylationParser`。

**uncertain 路徑（讀不確定，需另查）**：

7. **5mC vs 5hmC 在矩陣值的最終口徑**：methylation.csv 值在源碼層只能斷定是「一個甲基化機率」；要主張純 P(5mC) 需另查 MM/ML 多修飾解析口徑（不在本 claim 已證範圍）。本 C++ binary = 5mC-only 已定案，但若數據經 MSA Python 路徑則口徑不同。
8. **5hmC offset 跳過依賴 MM block 順序**（`C+h?` 在前、`C+m?` 在後，`:41` 註解假設）+ 非 target block 含 `,` 才推進。順序若反 / 空 block，offset 對位邏輯未測（對標準 ONT dorado MM tag 正常無影響）。
9. **`distance_use_binary`(`Config.hpp:50`) / `use_binary_matrix`(`DistanceMatrix.hpp`) 旗標**在 dispatcher 中未見被讀取（按 enum metric 選），疑為遺留 dead flag。
10. **`min_site_coverage = 5`**（`Config.hpp:36`）grep `src include` 無使用點，疑未接線 / dead config。
11. **`FisherExact::test_2x2` 雙尾/單尾 + 是否含 mid-p**：未讀 `FisherExact.cpp`；`seed=42`（`PerCpgAsm.cpp:267`）暗示某種隨機/Monte-Carlo 大表路徑，未驗證。
12. **「全碼無 beta-binomial/quasi/GEE」** 為稽核宣稱；本萃取只獨立確認 `PerCpgAsm.cpp` 該段是獨立 Bernoulli Fisher，未重 grep 全 src/core/*.cpp。
13. **不傳 `--distance-metric` 的 handoff/HKU pipeline 路徑**未逐一審；若存在會落回 NHD（二值）。

---

## Provenance footer

- **build_commit**: `069cadb`（session anchor；本檔撰寫時 working HEAD = `4ef557b`）。
- **稽核相關 commit**: `891e04b`（pre-audit 引用修）/ `bbcc7a9`（稽核 HEAD）/ `6593f96` + `eed4300`（post-audit 兩 🔴 修）。
- **中間檔**: `InterSubMod/docs/paper_focus/_ism_methyl_verify_20260613/{extract,verify}/`（extract 7 檔 + verify 5 檔；每條結論回源碼/數據對證）。
- **對抗驗證裁決**: 5mc_only=SUPPORTED ・ matrix_semantics=PARTIAL ・ mmml_tested=OUTDATED ・ dup_bug=REFUTED ・ threshold_use=PARTIAL。
- **數據真值**: flagship methylation.csv 本日 Python 全矩陣掃描（57×203 / 240 distinct / 50.2% NA），與 metadata.txt 完全一致。
