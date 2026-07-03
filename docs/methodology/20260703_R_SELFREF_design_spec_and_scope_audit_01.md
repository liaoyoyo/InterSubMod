<!--
建立時間: 2026-07-03
問題類型: C++ 方法 / 負控實驗設計 (HD-1 phasing by-construction 循環)
影響 track: paired (sSNV backbone + HP discriminator)
狀態: pending_decision
data_sources: src/core/ReadParser.cpp, include/utils/ArgParser.hpp, include/core/Config.hpp, memory/project_self_phasing_causal_chain_confirmed.md, docs/references/20260701_thesis_research_handoff_onboarding_01.md
-->

# R-SELFREF flag-on circularity 負控 — 設計 spec + scope 稽核

> 框架：Verdict-Pyramid。methodology-audit Step 1-4。目的：解 HD-1（phasing by-construction 循環）前，先釐清 R-SELFREF 要測什麼、規模、與現主軸框架的關係。**本輪 Step 1 兩個發現改變了原「~25-50hr from-scratch C++」的規模假設。**

## §1 TL;DR（裁決）

| 問題 | 裁決 | 證據 |
|---|---|---|
| R-SELFREF 要 from-scratch 寫 C++? | ❌ **否**——核心 flag `germline_hp_only` 已完整 wired | `ArgParser.hpp:127`（CLI `--germline-hp-only`）+ `Config.hpp:61` + `ReadParser.cpp:158` + `RegionProcessor.cpp:418` |
| R-SELFREF 實際工作量? | **~25-50hr「flag on/off × 7 樣本」對照跑 + 分析**（+ 可選 permutation null 小 C++/script）| flag 已在 → 主成本是 run 非 implement |
| R-SELFREF 防的 claim 現在承重嗎? | ⚠ **不承重**——現主軸 positive = sSNV genetic 骨幹（非循環），phasing/HP 是鑑別器非 headline | handoff §0；backbone audit §2「樹由 genotype 軸建、不靠 HP、不靠甲基」；grep「positive-spine 依賴 phasing」現框架無命中 |

**一句話**：R-SELFREF 不是數週 C++ 建構（flag 已在），是 ~25-50hr 對照跑；但它測的是**舊 G6「phasing 當 positive headline」框架的循環**，現主軸已 reframe 到 sSNV 骨幹 positive，所以 R-SELFREF 現在的角色從「解鎖 positive headline 的必需 gate」降為「HP 鑑別器的**經驗性非循環 robustness 檢查**」。

## §2 問題描述 — S1 自我參照循環機制

**循環（memory `project_self_phasing_causal_chain_confirmed`，2026-04-02 CONFIRMED）**：LongPhase somatic HP tag（`1-1`/`2-1`/`3`，整數 11/21/33）標記 somatic-only phase block —— **somatic 變異同時當 phasing anchor 又當被評估對象** → ALT reads 系統性偏向一個 HP → 虛假 LOH / 虛假 positive 結構。量化：**62.0% TO TP LOH 移除 self-phasing 後消失**（AF 0.1-0.8 近 100%）；31.2% TO LOH 是 self-phasing artifact；同位點 HP_Ratio 跨模式 r≈0。

**現有 C++ 緩解（已實作）**：`ReadParser.cpp:153-160` — `config_.germline_hp_only=true` 時，把 somatic HP tag（1-1/2-1/3）降級為 "0"，下游 HP 特徵只信 germline phasing（1/2）。CLI `--germline-hp-only`。

## §3 R-SELFREF 設計 — 要測什麼

**核心對照**：同一 pipeline 跑兩態 ×7 樣本：
- **flag-OFF**（預設）：信 somatic HP tag（含循環的 1-1/2-1/3）
- **flag-ON**（`--germline-hp-only`）：只信 germline phasing（1/2），循環 somatic tag 降級

**判準**：circularity-sensitive 指標（HP 鑑別 57% allelic 移除、same-HP 結構、任何 phasing-dependent claim、backbone determinacy）在 flag-ON 下**是否存活**。
- 存活（flag-ON 仍成立）→ **非 by-construction 循環**（防彈）
- 崩塌（flag-ON 消失）→ by-construction 循環 artifact，該 claim 須降級

**可選加強（permutation null）**：置換 somatic-flag 指派 ×N 重跑 → 真值 vs null 分佈，證觀測超越 null。

## §4 修改選項 + SWOT

### 方案 A：不做（承認現主軸已非循環，R-SELFREF 列 Future Work）
現主軸 positive = sSNV genetic 骨幹（backbone audit §2 已證非循環：TP/FP 只觀察、樹由 genotype 軸建、HP 只鑑別器、甲基 off-ladder）。phasing 非 positive headline → R-SELFREF 防的循環不承重。

|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** 零成本、資源投定稿 ⭐3 | **W** HP 鑑別器的非循環仍只「宣稱」非「經驗證」|
| External | **O** 對齊紅隊策略（先定稿⭐3、R-SELFREF Future Work）| **T** reviewer 可能仍問「HP 鑑別是否循環」，無經驗數據回 |

### 方案 B：跑既有 flag on/off × 7 樣本 + 分析（經驗非循環檢查）— **無 C++ 建構**
用既有 `--germline-hp-only` 跑兩態 ×7，比 HP 鑑別 / backbone determinacy 存活率。~25-50hr run。

|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** flag 已在、無 C++ 風險；直接經驗回答「非循環」；driver+分析 script 即可 | **W** ~25-50hr compute（7 樣本大 BAM ISM 全跑）；/big8 唯讀需 output 導 /big7 |
| External | **O** 把 handoff「first-digit germline anchor 非循環」從宣稱升經驗證據；防彈 HP 鑑別器段 | **T** 若 flag-ON 下 HP 鑑別大幅崩塌 → 反證循環，需重寫 §2（但這正是誠實科學該知道的）|

### 方案 C：B + permutation null（小 C++/script）— 最嚴謹
B 加置換 somatic-flag null 分佈。需小改（置換模式）+ ×N 重跑。

|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** 形式化 null、最強論證 | **W** +C++/script 小改（PDD + Hard Gate）+ ×N compute 倍增（數天）|
| External | **O** 方法學最防彈、可寫成 methods 亮點 | **T** 投入 vs 現框架不承重 → 報酬率低（紅隊：R-SELFREF 不解鎖 tier）|

## §5 驗收標準（B/C）
- [ ] flag-ON/OFF ×7 樣本全跑完、output 落 /big7（非 /big8 唯讀）
- [ ] 比對表：HP 鑑別 57% / backbone determinacy / same-HP 結構 在兩態的值 + 存活率
- [ ] 裁決：非循環（存活）or 循環 artifact（崩塌）+ 各 claim 對應處置
- [ ] （C）permutation null 分佈 + 真值百分位

## §6 成本 vs 報酬（誠實）
- **成本**：B = ~25-50hr run（無 C++）；C = +數天（C++ 小改 + ×N）。
- **報酬**：把 HP 鑑別器非循環從「宣稱」→「經驗證據」（robustness 加值）。**但不解鎖 tier**（紅隊 C4：真天花板是 single-cell/正交真值，非 phasing 循環）。
- **關鍵**：R-SELFREF 是舊 G6「phasing positive headline」的 gate；現主軸 reframe 後，它是 nice-to-have robustness 非 load-bearing。

## §7 用戶決策
**選擇**：[x] **C（B + permutation null，小 C++）**
**日期**：2026-07-03
**理由**：知道 flag 已 wired（B 無 C++）+ 不承重 + 真規模 ~25-50hr run 後，用戶仍選最嚴謹 C —— 把 HP 鑑別器非循環升為方法學防彈的經驗 + null 證據。

> ⚠ **重新告知已執行**：用戶原選「全套做（數週 C++）」基於錯誤前提；Step 1 揭 flag 已 wired + 現框架不承重後 re-decision → 確認走 C。

## §8 實作 spec（供 /cpp-change PDD）

**新增 C++ 能力 = permutation-null 模式**（核心 germline_hp_only flag 已在，不動）：
- **CLI**：`--permute-hp-seed <N>`（0=關；>0=以 seed N 做置換 null）。
- **位置**：`ReadParser.cpp` HP tag 解析後 / `RegionProcessor` 分群前——**per-region 內把 HP tag 值在 reads 間洗牌**（保 per-region HP 組成 marginal 不變，打斷 read↔somatic-tag 自我參照連結）。determinism：seed N + region key 派生 per-region RNG（不可用 `Math.random`/全域 rand，須 reproducible）。
- **不碰**：sSNV genotype 軸（backbone 非 HP-based，permute HP 不影響 determinacy 本體）；只影響 HP-鑑別 / same-HP enrichment 指標。
- **輸出**：每 region + 全域的 same-HP mutual_excl enrichment（觀測 0.86×）在 {觀測, flag-ON, permute×N} 三態的值 → null 分佈 + 觀測百分位。

**Step→Verify（cpp-change 6 步預告）**：
1. spec + 測試設計（本 §8）→ 2. 寫 `--permute-hp-seed` + per-region shuffle（reproducible RNG）→ 3. 單元測試（同 seed byte-identical、marginal 保持、seed=0 等同現行）→ 4. 編譯（Hard Gate）→ 5. 小規模 chr19 驗證（觀測 vs null 位移合理）→ 6. commit。
**驗收**：seed=0 與現行 byte-identical；permute 保 per-region HP marginal；chr19 null 分佈產出；full 7-樣本 ×{OFF, ON, permute×N} run 落 /big7。

## §9 2026-07-03 HD-1 收斂更新（S-vs-TO 澄清 → R-SELFREF 降 optional）

> 用戶 catch + 三方源碼驗證（`01_longphase_s.sh:156-160` + memory `reference_hp_tag_definition_and_subclone_caveat` + KB longphase-s/to）：**「phasing by-construction 循環」大半是 LongPhase-TO(tumor-only 真 self-phasing) 與 LongPhase-S(現 canonical, germline-anchored) 混淆**。詳 memory `project_hd1_circularity_longphase_s_vs_to_resolution`。

**精確裁決（非「循環不存在」，是「非承重」）——兩層拆開**：
1. **重建骨幹非循環 by design**：LongPhase-S 的 HP1/HP2 由 germline het SNP 定（`-s germline_phased_vcf -b normal_bam`），somatic 從獨立輸入（`--tumor-snv-file`）進來只做上層標記；重建骨幹＝germline-HP 上 **sSNV read-level 共現**（genotype 軸，非 HP tag）→ **非循環**（源碼證，強於 empirical）。
2. **循環的 phasing-spine 已 reframe 降階**：舊 G6 positive headline（NG=2 Inner same-HP1>Outer）用 somatic 子標記 `1-1/2-1` 的相內子結構，**那層當 positive headline 確有 by-construction 循環**；現主軸已把它降為 **§2.6 支撐 observation**（非 positive headline）→ 論文 positive 主張不依賴此循環元素。

→ **HD-1 從「唯一決定性卡關（blocking gate）」收斂為「非承重（resolved-as-non-blocking）」**：骨幹設計證非循環 + 循環 spine 已降階。**R-SELFREF（本 doc 工具，commit 95e6f62）降為 optional 經驗附錄**（可量化 spine 循環幅度，但不承重、不解鎖 tier）。論文 methods 引 `01_longphase_s.sh` 命令證骨幹非循環（design-level > empirical）。**RUN 無急迫性**，tooling 就緒待跑。
