# ISM 甲基測試覆蓋稽核 — `test_methylation_parser.cpp` 及鄰近測試

**稽核日期**: 2026-06-13
**範圍**: `tests/test_methylation_parser.cpp` + `tests/test_per_cpg_asm.cpp` + `tests/test_distance_matrix.cpp`
**紀律**: 只引實檔行號 + 原始碼片段；不臆測；數字附 `檔:行`

---

## L0 一眼結論

舊稽核説法「**MM/ML 解析零測試**」**已被反駁 / 過時** —— `tests/test_methylation_parser.cpp` 存在且含 **3 個 gtest test case**，直接覆蓋 MM delta 解碼、ML→機率、CpG-context filter、forward/reverse strand collapse、no-MM-tag 三條路徑。

但須**限定**：這 3 個測試是 **2026-06-10 13:51（與 06-10 全碼稽核同日）**新增的**補救 commit `6593f96`（"audit TESTS-1"）**。檔頭 docstring 明文寫「It was previously UNTESTED (audit 2026-06-10, finding TESTS-1)」（`test_methylation_parser.cpp:7`）。亦即：

- 舊稽核「MM/ML 解析零測試」在**該稽核當下（06-10 早）為真**；
- 該稽核**本身就是觸發補測的 finding TESTS-1**；
- 補測落地後（06-10 晚 commit `6593f96`），説法**從此不再成立**。

故正確表述是「**已修正的歷史缺口**」，而非「現存零測試」。

---

## L1 重點邏輯

| 維度 | 結論 | 來源 |
|------|------|------|
| 測試檔是否存在 | **是**，152 行 | `tests/test_methylation_parser.cpp`（`wc`/Read） |
| test case 數 | **3** | `test_methylation_parser.cpp:66,97,139`（grep `^TEST`=3） |
| 是否 wired 進 build | **是** | `CMakeLists.txt:169` 列入 |
| 新增 commit | `6593f96` `test(methylation-parser): add MM/ML golden characterization tests (audit TESTS-1)` | `git log --follow` |
| commit 時間 | **2026-06-10 13:51:39 +0800** | `git show -s 6593f96` |
| 被測單元 | `MethylationParser::parse_read` | `src/core/MethylationParser.cpp:15`；`test_*.cpp:79,107,116,148` |
| 測試性質 | golden / characterization（由 SAM MM/ML spec 獨立推導真值，非抄實作） | `test_methylation_parser.cpp:6-9` docstring |

---

## L2 細節：3 個 test case 各測什麼

### TC1 `ForwardStrand_DecodesDeltas_FiltersNonCpg_ProbFromML`（`test_methylation_parser.cpp:66-90`）
ref window `ACGAACGTTC` @ 0-based 1000（`:71`），C 在 seq idx 1/5/9，MM `"C+m?,0,0,0;"` 標三個 C 全 modified，ML `{230,130,255}`（`:74-75`）。
**斷言**:
- `calls.size()==2u`（`:82`）—— offset 1、5 是 CpG（C-then-G）保留；**offset 9 的 C 在 window 邊緣（無 offset 10）→ 必須 drop**（`:70,81`）。 → **測 CpG-context filter + 邊界 drop**。
- `calls[0].ref_pos==1002`（`:84`）、`calls[1].ref_pos==1006`（`:86`）→ **1-based 座標 = ref_start+offset+1**（對應源碼 `:181` `ref_pos_0based + 1`）。
- `probability ≈ 230/255`（`:85`）、`130/255`（`:87`）→ **ML[i]/255 機率轉換**（對應源碼 `:180`）。
- **測 MM delta 解碼**：`"C+m?,0,0,0"` 三個 delta=0（連續 C，零 skip）。

### TC2 `ReverseStrand_CollapsesToForwardCpgCoordinate`（`test_methylation_parser.cpp:97-133`）
ref `AACGTTACGT`，CpG 在 offset 2、7（`:101-102`）。先跑 forward read（`:104-111`）得 ref_pos `1003`、`1008`；再跑 **reverse read**（`flag=BAM_FREVERSE`，BAM SEQ 相同，ML `{200,100}`，`:114`）。
**斷言**:
- `rcalls.size()==2u`（`:117`）。
- reverse 落在 **與 forward 同一座標** `{1003,1008}`（`:124-125`）→ **two-strand CpG collapse**（雙股聚到 read×CpG 矩陣同一欄）。docstring 標此為「most error-prone branch（backward iteration over BAM SEQ）」（`test_methylation_parser.cpp:16-17`）。
- 機率反序映射：`rev_by_pos[1008]≈200/255`、`rev_by_pos[1003]≈100/255`（`:128-129`）→ **測 reverse strand backward iteration（MM 5'→3' 對 BAM SEQ 反向掃）**，對應源碼 `:120-157`（`is_reverse` 分支、`seq_idx` 由 `seq_len-1→0`、target_base `'G'`、CpG 驗證 `ref_seq[ref_offset]=='G' && ref_seq[ref_offset-1]=='C'` `:139`）。

### TC3 `NoMmTag_ReturnsEmpty`（`test_methylation_parser.cpp:139-151`）
建一條無 MM/ML tag 的 BAM record（`:141-146`，未 `bam_aux_append`）。
**斷言**: `calls.empty()`（`:149`）→ **測 no-MM-tag → 空 calls**，對應源碼 `:23-25`（`if (!mm_aux || !ml_aux) return calls;`）。

---

## L2b 已測 vs 未測（誠實對照）

### ✅ 已測路徑（3 TC 覆蓋）
1. **MM `C+m?` delta 解碼**（delta=0 連續情形）— TC1/TC2。
2. **ML[i]/255 → 機率**（forward 與 reverse 兩序方向）— TC1（`:85,87`）、TC2（`:128-129`）。
3. **CpG-context filter**（C-then-G 保留；非 CpG / window 邊緣 C drop）— TC1（`:82`）。
4. **forward strand 正向掃 + 1-based 座標報告**（`+1`）— TC1（`:84,86`）。
5. **reverse strand 反向掃 + two-strand collapse 到 forward-CpG 座標**（最易錯分支）— TC2。
6. **no-MM/ML-tag → 空**— TC3。

### ❌ 未測 / 未覆蓋路徑（仍是 gap，誠實列出）
1. **多 mod-code MM 字串的 `ml_offset` 前進**：源碼 `:43-83` 處理「先 `C+h?`（5hmC）後 `C+m?`」需跳過前一段 ML 偏移；3 個 TC 的 MM 全是**單一 `C+m?` 區塊**（`:75,104,114`），**從未測 `C+h?` 在前時 ml_offset 是否正確跳過**。→ 此為「5mC/5hmC max-collapse」相關路徑，**零測試**。 🔴
2. **`parse_mm_tag` 非零 delta（skip-count>0）**：所有 MM delta 皆為 `0`（`:75,104,114`），**未測** delta=1/2…（C 之間有 skip）的解碼正確性。源碼 `:201-230` / `next_target += deltas[delta_idx]+1`（`:150,188`）未被非零值驗證。
3. **CIGAR 非全-match 的 seq→ref 映射**：`make_meth_bam` 只造 `LxM`（全 BAM_CMATCH，`:45`）。`build_seq_to_ref_map`（源碼 `:232-281`）的 **I/D/S/N/H** 分支（`:258-272`）**完全未測**（含 soft-clip 起點偏移、insertion 不推進 ref 等真實 long-read 常見情形）。
4. **ML 格式 / 邊界防呆**：`ml_aux[0]!='B' || ml_aux[1]!='C'`（源碼 `:31-33`）、`ml_offset+deltas.size()>ml_len` 不足（`:94`）、`deltas.empty()`（`:89`）、`found_target=false`（`:85`，有 MM 但無 `C+m?` 區塊）→ 這些 early-return 防呆分支**皆未測**。
5. **`ref_pos_0based < 0`（seq 落在 insertion/clip → map=-1）跳過**：源碼 `:134,172` 的 `>=0` 守衛**未測**。
6. **二進位化閾值 / NA**：本檔**不測** raw→binary 或 NA(-1) 語意。註：MethylCall 只帶 `(ref_pos, probability)`（源碼 `:142,181`），**binary/NA 轉換不在 `parse_read` 內** → 不屬此 parser 測試職責（NA(-1) 由下游 `MethylationMatrix`/`DistanceMatrix` 處理，見 `test_distance_matrix.cpp` 的 `-1`/`NAN` 路徑）。

---

## L2c 鄰近檔覆蓋（佐證甲基測試非孤島，但與「MM/ML 解析」職責不同）

- **`test_per_cpg_asm.cpp`**：**20 個 TEST**（grep `^TEST`=20）。測 `bh_fdr_correction`（`:55-93`）、`compute_nme`（`:99-153`）、`compute_epipolymorphism`（`:159-195`）、`compute_per_cpg_asm` 全管線（`:201-369`）。含 NaN/missing 處理（`:142-153,301-325`）、HP1-1/HP2-1 somatic tag 映射（`:281-299`）、clear-ASM/no-ASM（`:216-279`）。**這是「下游 per-CpG ASM 統計」測試，非 MM/ML 解析**。
- **`test_distance_matrix.cpp`**：**23 個 TEST**（grep `^TEST`=23；含 `TEST_F`）。測 NHD/L1/L2/CORR/JACCARD/BERNOULLI 6 metric、min_common_coverage、strand-specific 拆分（`:293-316`）、NA(-1)/NAN missing（fixture `:36-101`）、對稱性/對角線/非負三大數學性質（`:542-601`）、CSV 輸出（`:446-464`）。BERNOULLI 為新增（"Task 2.1 previously untested metric"，`:467-468`），含 python-oracle golden 值（`:481-488`）。**這是「read-read 距離」測試，非 MM/ML 解析**。

> 三檔合計形成「**解析（parser）→ 距離（distance）→ 統計（per-CpG ASM）**」三段甲基管線的測試鏈；MM/ML 解析這一最上游段在 06-10 前確實是斷點，現已補上但僅 happy-path。

---

## L3 溯源

- commit: `6593f96`（`test_methylation_parser.cpp` 唯一 commit，`git log --follow`）；時間 `2026-06-10 13:51:39 +0800`（`git show -s`）。
- 被測源碼: `src/core/MethylationParser.cpp:15-199`（`parse_read`）、`:201-230`（`parse_mm_tag`）、`:232-281`（`build_seq_to_ref_map`）、`:283-288`（`is_cpg_site`）。
- header: `include/core/MethylationParser.hpp`（存在，3694 bytes）。
- build wiring: `CMakeLists.txt:169`。
