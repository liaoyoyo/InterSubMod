<!--
建立時間: 2026-06-13
任務: 對抗審查 claim 'mmml_tested' — 舊稽核「MM/ML 解析零測試」是否仍成立
紀律: 只讀實際檔案、引行號、引原始碼片段；不憑記憶；數字附 檔:行；有疑標 uncertain
驗證者: 獨立對抗審查員（預設懷疑，回原碼複核 extract claim）
-->

# 對抗驗證 — claim `mmml_tested`

## verdict: OUTDATED（舊稽核「MM/ML 解析零測試」在稽核當下為真，補測落地後已不成立）

舊稽核説法「**MM/ML 解析零測試**」**在 2026-06-10 早稽核時點（HEAD `bbcc7a9` @ 04:50）為真**，但**已被同日稍晚的補救 commit `6593f96`（13:51）反駁**。現行 repo（HEAD `4ef557b`）`tests/test_methylation_parser.cpp` 存在、含 3 個 gtest case、已 wire 進 build。故正確表述 = 「**已修正的歷史缺口**（happy-path only）」而非「現存零測試」。

我**獨立回原碼複核**了 extract（`../extract/tests.md` + `../extract/audit.md`）的每一條斷言，**全部對證成立**，未發現需要翻轉或校正的錯誤。

---

## L1 對抗複核表（extract claim → 我獨立對證的結果）

| extract 斷言 | 獨立對證 | 來源 檔:行 |
|------|:---:|------|
| 測試檔存在，含 3 TEST | ✅ 證實 | `tests/test_methylation_parser.cpp:66,97,139`（grep `^TEST`=3；wc=151 行，extract 寫 152 含 EOF 差 1 行，不影響） |
| wired 進 build | ✅ 證實 | `CMakeLists.txt:169` `tests/test_methylation_parser.cpp` |
| 新增 commit `6593f96` @ 06-10 13:51 | ✅ 證實 | `git log --follow` → `6593f96 2026-06-10 13:51:39 +0800` |
| 稽核 HEAD `bbcc7a9` 時測試 **ABSENT** | ✅ 證實 | `git cat-file -e bbcc7a9:tests/test_methylation_parser.cpp` → `fatal: Not a valid object name`（exit 128 = ABSENT） |
| commit 為 **post-audit** | ✅ 證實 | `git merge-base --is-ancestor bbcc7a9 6593f96` → exit 0（bbcc7a9 是 6593f96 祖先）；稽核 HEAD 時間 04:50 < commit 13:51 |
| 被測單元 = `parse_read` | ✅ 證實 | test `:79,107,116,148` 全呼 `parser.parse_read(...)`；源碼 `src/core/MethylationParser.cpp:15` |
| 性質 = golden/characterization（由 spec 推導非抄實作） | ✅ 證實 | docstring `test_methylation_parser.cpp:6-9,19-20`「derived independently from the SAM MM/ML specification, NOT from the implementation」 |

> 結論：extract 的 L0「OUTDATED / 已修正的歷史缺口」**站得住**。我未發現 extract 有過度宣稱或捏造行號。

---

## L2 — 3 個 test case 各測什麼（逐行對證源碼）

### TC1 `ForwardStrand_DecodesDeltas_FiltersNonCpg_ProbFromML`（`test_methylation_parser.cpp:66-90`）
- ref `ACGAACGTTC` @ 0-based 1000，MM `"C+m?,0,0,0;"` ML `{230,130,255}`（`:71,74-75`）。
- 斷言 `calls.size()==2u`（`:82`）—— offset 1/5 是 CpG 保留，offset 9 邊緣 C drop。**對證源碼**：forward 分支 CpG 驗證 `ref_seq[ref_offset]=='C' && ref_seq[ref_offset+1]=='G'`（`MethylationParser.cpp:177-178`），邊界守衛 `ref_offset+1 < ref_seq.size()`（`:176`）→ offset 9 無 offset 10 → drop。**測 CpG-context filter + 邊界 drop ✅**。
- 斷言 `calls[0].ref_pos==1002`、`calls[1].ref_pos==1006`（`:84,86`）。**對證源碼** `:181` `calls.emplace_back(ref_pos_0based + 1, prob)` → 1-based = ref_start+offset+1。**測座標報告 ✅**。
- 斷言 `probability≈230/255`、`130/255`（`:85,87`）。**對證源碼** `:180` `ml_data[ml_offset + delta_idx] / 255.0f`。**測 ML→機率 ✅**。
- MM `0,0,0` 三 delta=0（連續 C 零 skip）。**對證 `parse_mm_tag`**（`:201-230`，`from_chars` 逐 delta）。**測 MM delta 解碼（delta=0 情形）✅**。

### TC2 `ReverseStrand_CollapsesToForwardCpgCoordinate`（`test_methylation_parser.cpp:97-133`）
- ref `AACGTTACGT`，CpG offset 2/7。先 forward（`:104-111`）得 ref_pos `{1003,1008}`，再 **reverse**（`flag=BAM_FREVERSE`，BAM SEQ 同，ML `{200,100}`，`:114`）。
- 斷言 reverse 落同座標 `{1003,1008}`（`:124-125`）。**對證源碼**：reverse 分支倒序掃 `for (int seq_idx = seq_len - 1; seq_idx >= 0; seq_idx--)`（`MethylationParser.cpp:126`），target `'G'`（`:124`），CpG 驗證 `ref_seq[ref_offset]=='G' && ref_seq[ref_offset-1]=='C'`（`:139`），report at `ref_pos_0based`（`:142`，等同 C 位 1-based）。此為碼內自標 `CRITICAL`（`:107-110`）的最易錯分支。**測 two-strand collapse + 反股倒序迭代 ✅**。
- 斷言反序機率映射 `rev_by_pos[1008]≈200/255`、`rev_by_pos[1003]≈100/255`（`:128-129`）。**對證源碼** `:140` `ml_data[ml_offset + delta_idx]`，倒序掃使 ML[0]→較大座標。**測 reverse ML 反序映射 ✅**。

### TC3 `NoMmTag_ReturnsEmpty`（`test_methylation_parser.cpp:139-151`）
- 建無 MM/ML tag 的 record（`:141-146`，未 append）。斷言 `calls.empty()`（`:149`）。**對證源碼** `:23-25` `if (!mm_aux || !ml_aux) return calls;`。**測 no-tag→空 ✅**。

---

## L2b 已測 vs 未測（我獨立 grep 複核 extract 的 gap 清單 — 全部成立）

### ✅ 已測（3 TC 覆蓋）
1. MM `C+m?` delta 解碼（delta=0 連續）— TC1/TC2。
2. ML[i]/255→機率（forward + reverse 兩序）— TC1/TC2。
3. CpG-context filter（保留 C-then-G；drop 非 CpG / 邊緣 C）— TC1。
4. forward 正向掃 + 1-based 座標 — TC1。
5. reverse 倒序掃 + two-strand collapse（最易錯分支，碼內標 CRITICAL）— TC2。
6. no-MM/ML-tag→空 — TC3。

### ❌ 未測（仍是 gap — 我獨立 grep 證實，非僅信任 extract）
- **5hmC `C+h?` 在前的 `ml_offset` 前進**：`grep "C+h" tests/...` = **0 命中**；3 TC 全是單一 `C+m?` block（`test:75,104,114`）。源碼 `MethylationParser.cpp:71-82`（非 target block 數 delta 推進 `ml_offset += delta_count`）**未被測**。🔴（此即「5mC/5hmC max-collapse」相關路徑）
- **非零 delta（skip>0）**：test 所有 MM delta 皆 `0`（`:75,104,114`），未測 delta≥1（C 間有 skip）。源碼 `parse_mm_tag`（`:201-230`）+ `next_target += deltas[delta_idx]+1`（`:150,188`）的非零值未驗。
- **非全-match CIGAR**：`grep "BAM_C(INS|DEL|SOFT|REF_SKIP|HARD|EQUAL|DIFF)"` = **0 命中**；`make_meth_bam` 只造 `BAM_CMATCH`（`test:45,142`）。`build_seq_to_ref_map`（源碼 `:232-281`）的 I/D/S/N/H 分支（`:258-272`）**完全未測**（含 soft-clip 偏移、insertion 不推進 ref）。
- **格式/邊界防呆**：`ml_aux[0]!='B'||[1]!='C'`（`:31-33`）、`ml_offset+deltas.size()>ml_len`（`:94`）、`deltas.empty()`（`:89`）、`found_target=false`（`:85`）— early-return 分支皆未測。
- **`ref_pos_0based<0` 守衛**（`:134,172`）未測。
- **binary/NA 二值化**：不在 `parse_read` 內（MethylCall 只帶 `(ref_pos, probability)`，`include/core/MethylationParser.hpp:14-20`）→ 非此 parser 測試職責。

> ⇒ 補測為 **happy-path / 主要分支 only**：覆蓋了「正股 + 反股 + 多 CpG + window 邊界 + no-tag」，但 5hmC 多 block、非零 delta、indel/clip CIGAR、early-return 防呆 **仍零測試**。

---

## L3 溯源（時間軸）

| commit | 日期 | 性質 | 對 TESTS-1 |
|---|---|---|---|
| `bbcc7a9` | 2026-06-10 04:50 | **= 稽核 HEAD** | 測試檔 ABSENT（`git cat-file -e` exit 128）→ 「零測試」當下為真 |
| `6593f96` | 2026-06-10 13:51 | **post-audit**（`merge-base --is-ancestor bbcc7a9 6593f96`=YES）| 新增 `test_methylation_parser.cpp`（3 golden TC）|
| `4ef557b` | （現行 HEAD）| — | 測試在線、wired（`CMakeLists.txt:169`）|

- 被測源碼：`src/core/MethylationParser.cpp:15-199`（`parse_read`）/`:201-230`（`parse_mm_tag`）/`:232-281`（`build_seq_to_ref_map`）/`:283-288`（`is_cpg_site`）。
- header：`include/core/MethylationParser.hpp:14-20`（`MethylCall{ref_pos, probability}`）。

---

## uncertain / 未做之事（誠實標註）
1. **未實際編譯/執行測試**：`build/tests/test_methylation_parser` 無 prebuilt binary（本稽核 read-only，未 `make`+`ctest`）→ 斷言「3 TC 全 PASS」**未經 runtime 證實**，只證實「source 邏輯與源碼一致、斷言值與 spec 一致」。標 uncertain（若要 SUPPORTED-with-green，需跑 `ctest -R MethylationParser`）。
2. 行數小差：extract 寫 152 行，實測 `wc -l`=151（EOF 換行差 1），不影響結論。
3. 反股 `bam_set1` 的 BAM SEQ 在 reverse case 是否真等同 ref-forward（test 註解 `:41` 宣稱「BAM SEQ always stored reference-forward」）—— htslib 行為層面未獨立驗，信任 test 設計者；但 TC2 forward/reverse 用同一 SEQ 得同座標，邏輯自洽。
