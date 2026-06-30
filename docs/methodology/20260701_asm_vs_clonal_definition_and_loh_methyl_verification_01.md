<!--
建立時間: 2026-07-01
類型: methodology 定稿 — 前提定義(ASM vs clonal 甲基差異)+ LOH genotype-向量甲基直接驗證 + bounded-auxiliary 定量收斂
狀態: concluded(定義已釐清 + 3 個 compute 實證完成;甲基維持 bounded-auxiliary)
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/founder_imputation_loocv.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/methyl_aux_flag.json,docs/methodology/_assets/20260627_subclone_4axis_teaching/data/loh_genotype_methyl_verify.json,docs/methodology/20260628_cis_control_scope_pilot_verdict_01.md
provenance: HCC1395 canonical tagged BAM(big7 canonical/HCC1395/paired_full/*_complete_matrix/longphase_s)。腳本已 commit(founder_imputation_loocv.py / methyl_aux_flag.py / loh_genotype_methyl_verify.py);data JSON 為 gitignore build artifact,可由腳本重跑重現(seed=20260630)。每數字標來源檔(§13-C)。
-->

# ASM vs Clonal 甲基差異 — 前提定義釐清 + LOH genotype-向量直接驗證

> 框架：Verdict-Pyramid（裁決先行 → 定義 → 證據 → 邊界）。HCC1395 ⭐3 單樣本。每個數字標來源檔（§13-A）。
> 起因：2026-06-30~07-01 用戶連續 catch 兩個定義層問題 —— ①「同 family 內(HP1-1 vs HP1-2)甲基差是 clonal 不是 ASM」②「LOH 下 germline-ASM 抵消,genotype-向量甲基差是否 = somatic 支持」。本文件釐清前提定義 + 用 3 個 compute 直接驗證。

---

## §1 TL;DR（裁決）

| 問題 | 裁決 | 證據級 |
|---|---|---|
| ASM = HP1 family vs HP2 family(跨 germline allele)的定義對嗎？ | ✅ **正確**(文獻 axis-relative + 源碼 HP-axis 操作一致) | L2 |
| HP1-1 vs HP1-2(同 family、不同 somatic 子相)= clonal 非 ASM？ | ✅ **正確**(同 germline allele → germline-ASM 抵消、殘差為 somatic) | L2 |
| 我前述把「germline vs mutated」Δβ 叫「cis-ASM」對嗎？ | ❌ **用詞錯**(該分組是 ASM+clonal 混合 hybrid,已修正) | — |
| LOH 下 RR vs AR / AR vs RA 有甲基差異嗎？ | ✅ **真實存在**(LOH 276 配對 near Δβ 中位 0.0372) | L2 |
| 該差異是 germline-ASM 嗎？ | ❌ **不是**(跨 CN 分層 near Δβ 幾乎一致 → 普遍 somatic-cis,非 germline-ASM) | L2 |
| 該差異 = subclone lineage marker 嗎？ | ⚠ **主要是突變局部 cis 足跡**(near/distal=1.46;僅 6-11% 遠端顯著) | L3 |
| 甲基整體角色 | **bounded-auxiliary / characterization-only**(負篩 ✅ + L3 軟旗標,非確認器) | L2 |

**一句話**：用戶的 ASM(跨 family)vs clonal(family 內)定義**完全正確,無認知錯誤**;LOH 下 genotype-向量甲基差**真實且是 somatic(非 germline-ASM)**——但**主要是突變局部 cis 足跡,lineage 殘餘訊號弱**,只 6-11% 區(LOH 21 區)有超出足跡的遠端訊號。甲基維持 bounded-auxiliary。

---

## §2 前提定義釐清 — ASM vs Clonal 甲基差異

### §2.1 源碼確認的 HP tag 語意（`src/core/ReadParser.cpp:121-149`）

| tag(longphase-s 字串)| 語意 | 軸 |
|---|---|---|
| `HP:Z:1` / `HP:Z:2` | germline 親代單體型(HP1 family / HP2 family) | **germline allele 軸** |
| `HP:Z:1-1` / `HP:Z:2-1` / `HP:Z:3` | somatic-integrated phase(演算法指派的 somatic 子相) | somatic/clonal 軸 |

> 規則：**第一個 dash 前的數字 = germline family**。`ReadParser.cpp:145-149` 在 self-phasing fallback 時把 somatic tag demote 成 "0"，理由原文「these reflect [self-phased somatic blocks], not [germline] haplotypes」——程式碼本身就把兩軸分開對待。

### §2.2 用戶定義（✅ 正確）

| 定義 | 對否 | 依據 |
|---|:--:|---|
| **ASM** = HP1 family vs HP2 family(不同 germline 親代 allele 間的甲基差) | ✅ | 文獻 ASM 永遠 allele/haplotype 軸(Shoemaker 2010 PMID20418490 Cat I/II/III, axis-relative)；本專案字典 line 52「ASM 哪兩組：HP-axis(HP1/HP2)」 |
| HP1-1 vs HP1-2（同 HP1 family、不同 somatic 子相）= **clonal 甲基差異**，非 ASM | ✅ | 同一 germline allele → germline-ASM **抵消**；殘差只能是 somatic/clonal。字典 line 65「-2 = 第二 somatic 子相(multi-subclone marker)」 |
| HP1-1 vs HP1-1-1（巢狀更深子相）= clonal | ✅(概念) | 原理同上；惟 longphase-s 實際 tag 是扁平(`1-1`/`2-1`/`3`)，無 `1-1-1` 這層；`HP3`(無 dash)無明確 germline family，不適用此規則 |

### §2.3 修正：我前述「cis-ASM」用詞錯誤

founder-imputation / aux_flag 的 `germline(全R) vs mutated(≥1 ALT)` Δβ，我前一輪叫它「cis-ASM」**不精確**。該分組是 **hybrid**：
```
mutated 群  (≥1 ALT)  = HP1-1                ← somatic 子相
germline 群 (全 R)    = HP1-ancestral + HP2  ← 混進了 HP2
```
→ Δβ = **ASM 成分**(HP2 漏進 germline 群)+ **clonal 成分**(HP1-1 vs HP1-ancestral)混在一起。正確叫法 = **「genotype-軸甲基分離(ASM + clonal 混合)」**，非純 cis-ASM。

---

## §3 LOH 改變 confound 結構（加法模型）

已驗證加法模型（`20260628_cis_control_scope_pilot_verdict_01.md` §10.5，4 角度對抗一致）：
```
β = μ + g(germline-HP) + c(somatic-cis 突變足跡) + t(clonal-lineage) + κ·CN + ε
```

| 比較 | g(HP) germline-ASM | c(somatic-cis) | t(lineage) |
|---|---|---|---|
| CN-neutral het：RR vs AR | ❌ 混進來(RR 含 HP2) | ✅ 在 | ✅ 在 |
| **LOH：RR vs AR** | ✅ **抵消**(只剩一條 HP) | ✅ 在 | ✅ 在 |
| **LOH：AR vs RA（遠端 CpG）** | ✅ 抵消 | ✅ 抵消(差異位點不同、遠端足跡衰減) | ✅ **殘餘 = 最乾淨 lineage** |

**用戶 LOH 直覺對的部分**：LOH 只剩一條 germline allele → germline-ASM 自動抵消 → genotype-向量甲基差**不是 ASM**(是 somatic)。
**要補的 caveat**：LOH **沒移除** `c(somatic-cis)`——突變本身在**局部**改變甲基。要從 lineage 拆出，靠 **near/distal**：遠端 CpG 無突變足跡 → 遠端仍有 Δβ 才算 lineage。

---

## §4 實證結果（3 個 compute，HCC1395 canonical full long-read all-sites）

### §4.1 Founder-imputation LOOCV（甲基能否 impute uncovered read 的 genotype）

半監督非循環：群由 sSNV 定義(germline=全R vs mutated=≥1 ALT)，甲基當分類器，LOOCV 驗準確度。`founder_imputation_loocv.json`（1598 可測區）：

| 指標 | 值 |
|---|---|
| balanced_acc 中位 / 平均 | **0.566 / 0.560**（0.5=隨機）→ 整體僅略勝隨機 |
| **corr(Δβ_centroid, balanced_acc)** | **0.617**（強正相關 → 預測力由 ASM 強度驅動）|
| 弱 ASM <0.1（n=**1513**, 95%）| balanced_acc **0.560 ≈ 隨機** |
| 中 ASM 0.1-0.2（n=68, 4%）| **0.778** |
| 強 ASM ≥0.2（n=17, 1%）| **0.937** |

**裁決**：甲基 impute genotype **整體 ≈ 隨機**；僅在強 ASM 區(~5%)有效。能力 = cis-ASM 強度的函數（非循環但 ASM-conditional）。

### §4.2 aux_flag — 差異分群顯著性（permutation-FDR）

每區 germline-vs-mutated 群間甲基 |Δβ| + 2000× 標籤置換 + BH-FDR。`methyl_aux_flag.json`（1598 區）：

| 指標 | 值 |
|---|---|
| perm-p FDR<0.05（顯著差異分群）| **178 (11.1%)** |
| **FLAG=True**（FDR<0.05 AND Δβ≥0.1，真輔助標記可用）| **71 (4.4%)** |
| flagged 區 LOOCV balanced_acc 中位 | **0.817** |
| 非 flagged 區 balanced_acc 中位 | **0.561** |
| FDR<0.05 但 Δβ<0.1（顯著但效應微小，不入 flag）| 107 |

**裁決**：71 區(4.4%)有「真差異分群」(顯著 + 有效應)。raw balanced_acc 會把 126 個小樣本假高誤計，FDR + 效應量門檻濾掉(僅 71 入 flag)。可當 **L3 輔助標記欄**(cis-ASM 一致性旗標)，非 subclone 確認器。

### §4.3 LOH genotype-向量甲基差異 — near/distal 直接驗證

LOH 全集，BAM 自算 genotype，RR_vs_1ALT 與 sib_1ALT(AR vs RA)配對，near(±1kb)/distal Δβ + 置換 + BH-FDR。`loh_genotype_methyl_verify.json`（2048 配對）：

| CN ＼ 配對 | n | Δβ_near 中位 | Δβ_distal 中位 | distal 顯著(FDR<.05) |
|---|--:|--:|--:|--:|
| **loh** RR_vs_1ALT | 197 | 0.040 | 0.0254 | **6.1%** (12) |
| **loh** sib(AR vs RA) | 79 | 0.031 | 0.0261 | **11.4%** (9) |
| neutral RR_vs_1ALT | 49 | 0.040 | 0.0227 | 10.2% (5) |
| neutral sib | 32 | 0.030 | 0.0193 | 15.6% (5) |
| gain RR_vs_1ALT | 1135 | 0.041 | 0.0247 | 13.0% (148) |
| gain sib | 550 | 0.036 | 0.0280 | 16.4% (90) |

**LOH 全部：Δβ_near 0.0372 vs Δβ_distal 0.0255，near/distal = 1.46**

**四個關鍵讀法**：
1. **現象存在 ✅** — LOH 區 RR vs AR、AR vs RA 確實有甲基差異。
2. **差異小** — 中位 Δβ_near ~0.04、distal ~0.025，遠低於 ASM 門檻 0.2。
3. **🔑 是 somatic 非 germline-ASM** — near Δβ 跨 loh(0.037)/neutral(0.040)/gain(0.041)**幾乎一致** → 普遍 somatic-cis 效應，**非 germline-ASM(否則 LOH 應顯著較小)、非 CN-multiplicity**。用戶 LOH 區分獲數據支持。
4. **主要是突變局部 cis 足跡** — near/distal=1.46；distal 中位 0.0255 但 ~90% 不顯著(只 6.1%/11.4% FDR<.05) → 多數遠端差 ≈ 雜訊地板。僅 **LOH 21 區**(12 RR_vs_1ALT + 9 sib)有顯著遠端 Δβ = lineage-candidate(仍小、需 orthogonal 驗證)。

---

## §5 結論

1. **前提定義（§2）**：用戶 ASM(跨 family) vs clonal(family 內) 定義正確、與文獻 + 源碼一致，**無認知錯誤需修正**；修正的是我「cis-ASM」用詞(該為 genotype-軸 hybrid)。
2. **LOH（§3-4.3）**：LOH 下 genotype-向量甲基差**真實且 somatic(非 germline-ASM)**——跨 CN 一致性是直接證據。但**主要是突變局部 cis 足跡**(near/distal 1.46)，lineage 殘餘訊號弱(6-11% 遠端顯著)。
3. **甲基整體**：維持 **bounded-auxiliary**——
   - ✅ **負篩**(95.8% 群 unimodal、95% 區弱 ASM 無分群)可信排除；
   - ✅ **L3 輔助標記**(aux_flag 71 區 cis-ASM 一致性旗標 / LOH 21 區 lineage-candidate)；
   - ❌ **不可**當 subclone 確認器、impute genotype(整體隨機)、lineage 排序(弱未顯著)。
4. **論文 A-framing 獲強化**：reconstruction 歸 sSNV 共現骨幹；甲基 characterize 有界 + somatic-cis 真實但單 bulk 不可識別 footprint vs lineage（需 single-cell / multi-region，Tarabichi LEARN）。

---

## §6 數字溯源表（§13-C）

| 數字 | 值 | 來源檔(腳本可重跑) |
|---|---|---|
| founder N / ba 中位 / corr | 1598 / 0.566 / 0.617 | `founder_imputation_loocv.json`（`founder_imputation_loocv.py`）|
| founder ASM 分層 ba(弱/中/強) | 0.560(n1513) / 0.778(n68) / 0.937(n17) | 同上 |
| aux_flag FDR<.05 / FLAG | 178(11.1%) / 71(4.4%) | `methyl_aux_flag.json`（`methyl_aux_flag.py`）|
| aux flagged / 非flagged ba | 0.817 / 0.561 | 同上 |
| LOH 總配對 / loh 配對 | 2048 / 276 | `loh_genotype_methyl_verify.json`（`loh_genotype_methyl_verify.py`）|
| LOH near / distal 中位 | 0.0372 / 0.0255（比 1.46）| 同上 |
| LOH distal 顯著率(RR_vs_1ALT / sib) | 6.1% / 11.4% | 同上 |
| cis-control LOH cross-tab(CROSS/SAME/MIXED) | 10 / 467 / 27 | `hp_alignment_fullscan.json`（見 20260628 verdict）|

> ⚠ data JSON 為 gitignore build artifact；可由已 commit 腳本 + canonical tagged BAM 重跑重現(seed=20260630)。

---

## §7 REFLECTION（給下次同方向 agent）

**警示指標**：若再問「甲基能否確認/排序/impute subclone」→ 本輪已 3 路實證封頂：impute 整體隨機(0.566)、ordering 弱未顯著(ρ0.18 p0.06)、LOH lineage 殘餘僅 6-11%。甲基 = bounded-auxiliary，勿再開「甲基驅動重建」。

**根因（double-loop）**：甲基沿 genotype 軸的差異**真實**，但其成分(germline-ASM / somatic-cis / lineage)在單一 bulk **結構性不可識別**；LOH 抵消 germline-ASM 後仍剩 somatic-cis ⊥ lineage 混淆。前提（單 bulk 可分離成分）錯，非方法錯。

**Reopen 條件**（scientific-rigor §8.3.1）：C1 COLO829/第二樣本 multi-region；C3 single-cell 甲基真值。

**Spaced recall**：2026-07-31（30d）。
