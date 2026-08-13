<!--
建立時間: 2026-08-13
狀態: validated
目標: 判定 frozen binary exact_ps_topology_af 的真實源碼，並更正 20260801 登錄表的錯誤
處理範圍: InterSubMod research CLI + LongLineage solver 全部 5 個模組
build_branch: chore/handoff-20260813
build_commit: 73afaeac
worktree: dirty
data_sources:
  - /bip7_disk/liaoyoyo2001/build-exact-af-20260724/（frozen build tree，含 41 個 .o 與 compile_commands.json）
  - /big7_disk/liaoyoyo2001/LongLineage @ b979760 與 @ HEAD
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/20260801_權威資料與過期來源登錄_01.md
驗證方式: 用 compile_commands.json 的原始旗標逐字重編兩個候選，三層比對（cmp / objdump / nm）
-->

# Frozen binary 源碼判定：`exact_ps_topology_af`

## Verdict：`RESOLVED_CANDIDATE_A`

**論文核心數字的重建鏈完整，且已用 byte-identical 證明。**

先前 `20260805/handoff_manifest.json` 標記的 `not_proven_byte_reproducible` 現已可解除；
`20260801_權威資料與過期來源登錄_01.md` §3 登錄的 LongLineage 側 SHA **有誤**，更正見 §4。

---

## 1. 判定對象

| 項目 | 值 |
|---|---|
| binary | `/bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af` |
| size / mtime | 280,176 B ／ 2026-07-24 16:50 |
| sha256 | `ba13ccf23d091854c191f81dd97fa891368d11179df9a69e915f0340b7233b2e` |
| 產出 | 論文的 71,955 read-AF ranked units、exact topology signature census |

這支 binary 的源碼**橫跨兩個 repo**，這是先前判定困難的根源。

---

## 2. 判定方法

`build-exact-af-20260724/` 整個 build tree 完整保留（含 41 個 `.o` 與 `compile_commands.json`），
使得判定不必靠推測 —— **可以逐字重跑當初的編譯指令，再比對中間產物 `.o`**。

原始編譯旗標（取自 `compile_commands.json`）：

```
/usr/bin/c++ -DBOOST_ALL_NO_LIB -I<LL>/include -O3 -DNDEBUG \
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion \
  -Wformat=2 -Wundef -Werror=return-type -fno-fast-math -ffp-contract=off \
  -Werror -std=c++17 -o <out>.o -c <LL>/src/solver/<module>.cpp
```

兩個候選（在 `/bip7` 的臨時 tree 編譯，**未觸動 LongLineage working tree**）：

| 候選 | 來源 | `obligation_bnb.cpp` sha256 | 時序 |
|---|---|---|---|
| **A** | LL commit `b979760`（2026-07-22 22:55） | `ba4924bf3b7d…` | binary 建於 07-24，其間無 commit |
| **B** | LL worktree（2026-07-26 04:12） | `c499954e1083…` | **晚於 binary 兩天** |

註：候選 B 的 sha256 正是 20260801 登錄表所記載的值。

---

## 3. 判定結果

### 3.1 `obligation_bnb.cpp` 三層比對

| 層級 | Candidate A | Candidate B |
|---|---|---|
| **L1 byte-identical** (`cmp`) | ✅ **完全相同** | ✗ 28,665 bytes 相異 |
| **L2 指令序列** (`objdump -d`) | ✅ 相同（4 行差異僅為 objdump 輸出的檔名標頭） | ✗ 3,868 行相異 |
| **L3 符號表** (`nm --defined-only`) | ✅ 相同 | ✗ 6 行相異 |
| `.o` 大小 | 57,392 B（＝frozen） | 57,416 B |
| `.o` sha256 | `37e203682c74c678cce28d08e7024c4e8e98fdda73cd6d71cee7886723c82e26` | `024febbc7c024d58…` |

frozen `.o` 的 sha256 = `37e203682c74c678cce28d08e7024c4e8e98fdda73cd6d71cee7886723c82e26` — **與 Candidate A 完全一致**。

### 3.2 全部 5 個 solver 模組（本次新發現）

先前文件只提到 `obligation_bnb` 與 `parent_mapping` 兩個模組，但 frozen build tree 內實際有 **5 個** solver `.o`。
以 Candidate A 逐一重編並比對：

| 模組 | 結果 |
|---|---|
| `obligation_bnb` | ✅ BYTE-IDENTICAL |
| `parent_mapping` | ✅ BYTE-IDENTICAL |
| `evidence_builder` | ✅ BYTE-IDENTICAL |
| `small_q_oracle` | ✅ BYTE-IDENTICAL |
| `terminal_subset_dp` | ✅ BYTE-IDENTICAL |

**5/5 全數 byte-identical。**

### 3.3 為什麼 B 不可能

Candidate B 相對 A 的差異是新增 `options.seed_incumbent` —— 用 certified optimum（objective subset DP）
預先播種 B&B 上界以剪枝。這個功能於 2026-07-26 加入，**晚於 binary 建置時間 2 天**，
因此那顆 binary 不含此功能。編譯結果的 24 bytes 差距與 3,868 行指令差異印證了這一點。

該功能後續已 commit 為 LongLineage `9ad976b`（2026-08-13）。

---

## 4. 對 20260801 登錄表的更正

依 evidence_ledger 的 append-only 慣例，**不修改 20260801 的凍結文件**，於此登錄更正：

| 欄位 | 20260801 登錄值 | 更正後 |
|---|---|---|
| `src/solver/obligation_bnb.cpp` | `c499954e10836b7d0fea0888c5b68779456bec03dadb84d481a89bb10985d97a` | `ba4924bf3b7d…`（LL `b979760`） |
| `include/…/obligation_bnb.hpp` | 未登錄 | `1d19dfb9a107…`（LL `b979760`） |
| `src/solver/parent_mapping.cpp` | `95aa509b9bca1e3eef67be25f56b48b42373c5e127decb8a31c7485414163862` | **不變**（兩候選一致，登錄正確） |
| 未登錄的三個模組 | — | `evidence_builder` / `small_q_oracle` / `terminal_subset_dp`，皆屬 `b979760` |

**錯誤成因**：登錄當時（08-01）取的是 LongLineage working tree 的當下狀態，
而該 working tree 已在 07-26 被修改，晚於 binary 的 07-24。這是「用當下 worktree 當作歷史產物的來源」的典型陷阱。

---

## 5. 完整重建指令

```bash
# ① InterSubMod 側 CLI（已 tracked，sha256 629fc309…，mtime 07-24 16:48）
ISM=/big7_disk/liaoyoyo2001/InterSubMod
CLI=$ISM/research/20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp

# ② LongLineage 側：checkout 到 b979760
git -C /big7_disk/liaoyoyo2001/LongLineage archive b979760 | tar -x -C <TMP>

# ③ 用原始旗標編譯（見 §2）
```

時序一致性佐證：CLI 源碼 mtime 16:48 → binary mtime 16:50，相隔 2 分鐘，合理。

**永久凍結副本**：`/bip7_disk/liaoyoyo2001/_frozen/exact_ps_topology_af_20260724/`（見 P1.3 產出的 bundle receipt）

---

## 6. 引用口徑：「只引 InterSubMod」與「binary 連結了 LongLineage」如何並存

`20260805/handoff_manifest.json:59` 的 `citation_rule` 規定「論文軟體引用只引 InterSubMod + commit hash，
絕不引 LongLineage」。而本判定確認 frozen binary 確實靜態連結了 LongLineage 的 5 個 solver 模組。
兩者並不矛盾，但**先前沒有任何文件處理這個張力**。

### 定案說法

> 引用鐵則約束的是**研究結論的出處歸屬**（哪個系統的分析產生了這個數字），
> 不是**編譯期的程式碼複用**。
>
> frozen binary `exact_ps_topology_af` 是 InterSubMod 的 research CLI —— 其分析語意、輸入契約、
> 輸出 schema 與 claim boundary 全部定義在 InterSubMod 側。它在編譯期靜態連結了 LongLineage 的
> 5 個純組合最佳化模組。

### 支持依據（本次實測）

對 Candidate A 的 `src/solver/` 8 個 `.cpp` 掃描甲基化相關語意：

| 關鍵字 | `methyl` | `Methyl` | `cpg` | `CpG` | `beta_value` | `m1_` | `label` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 命中數 | **0** | **0** | **0** | **0** | **0** | **0** | **0** |

solver 的輸入型別只有 `pattern_raox`、`base_qualities`、`multiplicity`
（`struct TopologyProblem` / `struct ExactTopologyProblem`）—— 不引入任何甲基化語意，
因此不構成「引用 LongLineage 的研究結論」。

### 論文 Methods 章節的正確寫法

- **演算法實作說明處**：註明 solver 的來源與版本（LongLineage `b979760`，5 個模組，byte-identical 已驗證）
- **結果歸屬處**：只引 InterSubMod + commit hash

### 同時必須誠實記載

LongLineage 自身狀態：`release_attestation = NOT_READY`、P0–P8 九個 gate 全開、0 個 VERIFIED、
formal topology 端點在全 autosome 上輸出 **0 units**（見
`InterSubMod/docs/handoff/20260806_LongLineage充分性稽核與路線裁決_01.md`）。

這不影響上述論述 —— frozen binary 用的是 solver **模組本身**，不是 LongLineage 的 pipeline 端點。

---

## 7. 對交接的意義

| 先前狀態 | 現在 |
|---|---|
| `not_proven_byte_reproducible` | **5/5 模組 byte-identical，已證明** |
| solver 源碼未進版本控制 | 已 commit（LL `9ad976b`／`6ce62b2`） |
| 登錄的 SHA 指向錯誤版本 | 已更正並記錄成因 |
| 「引用鐵則 vs 事實」無人處理 | 已定案，附實測依據 |

**論文核心數字現在可以從版本控制完整重建。** 這是本次交接整理中風險最高、也最有價值的一項。
