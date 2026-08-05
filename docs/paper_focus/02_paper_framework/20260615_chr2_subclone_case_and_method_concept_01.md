<!--
建立時間: 2026-06-15
狀態: 整理 (chr2:18M subclone 旗艦案例最終 verified 敘述 + 方法實作構想基礎; 供 Ch4 4.6 + Ch3 Methods C3 + handoff B2)
報告類型: subclone_case_consolidation + method_implementation_concept
受眾: 廖子游 · PI · 碩論 Ch3/Ch4 撰寫 · 系統化 session (handoff B2)
framework: Verdict-Pyramid + Assertion-Evidence (Tier 4 paper-scope)
data_sources:
  - InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_HCC1395_01.md (主報告)
  - InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_independent_verdict_02.md (第二次獨立 audit 最終裁決)
  - .../20260615_chr2_18M_subclone_verification_assets/data/independent_audit.{json,md} (機器輸出, byte-identical 可重現)
  - memory: project_chr2_18m_subclone_locus_verification (2026-06-15 cross-AI verified)
provenance_note: 所有數字引自上列 verified 報告 + 機器 audit JSON(獨立 session 產, 本 session memory §47-53 已 cross-AI 重算核對 byte-identical/逐位吻合)。本檔=整理/consolidation, 不產新數字; §13.0 撰寫與分析分離。最終敘述採 verdict_02 定案用語。
-->
<!-- provenance-verified: linkage 00/10/01/11 / 甲基 q-value / ASM-control / homopolymer context 全引自 independent_verdict_02 + audit JSON; 採 verdict_02 最終裁決。 -->

# chr2:18M Subclone 旗艦案例 — 最終敘述 + 方法實作構想基礎

> **這份是什麼**：把 chr2:18,068,480–18,108,828（HCC1395）subclone 重建案例的**已驗證完成資訊**整理成兩部分 —— **Part A 案例**（最終可防守敘述 + 觀察資料）與 **Part B 方法實作構想基礎**（從案例萃取的 ISM subclone 分析 pipeline 概念，供系統化）。
> **驗證狀態**：L2（單位點 × 單樣本 × 單 pipeline；HKU_MOD baseline + DORADO cross-basecaller 為主要 uplift）。經 4-agent 對抗複核 + 第二次獨立 audit（機器輸出 byte-identical 可重現）。
> **不改原件**：原驗證報告/audit/教學 HTML 在 `docs/experiments/.../20260615_chr2_18M_*` + `docs/explain/04,05_*`；本檔連結整理，不修改。

---

## L0 — 最終可防守敘述（一句話，採 verdict_02 定案）

> **chr2:18M = 在 SEQC2-confirmed LOH 區中，由多個 somatic allele linkage 定義、帶一致局部 5mC 結構的 regional subclonal structure；可防守用語 = 「a regional, LOH-constrained, somatic-haplotag-conditioned subclonal structure with cross-basecaller methylation coherence」。**
>
> ❌ **不可寫**：「confirmed five-subclone evolutionary reconstruction」/「甲基造成突變」/「甲基獨立重建完整 tree」/「五個 subclone 已證實」/「完整演化順序已確認」/「DORADO = 獨立病人/生物 replicate」（它是同細胞株技術重現）。

---

# Part A — Subclone 案例（verified 觀察 + 最終敘述）

## A1 區域與資料
- **區域**：chr2:18,068,480–18,108,828（~40 kb），HCC1395 乳癌細胞株。
- **資料**：HKU_MOD（5mCG_5hmCG, haplotagged）baseline + DORADO（raw, 5mC）cross-check；SEQC2 truth VCF + HC BED + CNV/LOH BED；GRCh38 ref。
- **方法**：pysam 萃取每 read 在 6 sSNV 的 base（含 DEL）+ HP/PS + MM/ML 解碼 5mC → tumor+normal；數字先落 JSON→讀回→才撰寫。

## A2 六個 sSNV（5 TP + 1 truth-gap，ref/alt 全符 SEQC2）
| # | 座標 | 變異 | SEQC2 狀態 |
|:--:|---|:--:|---|
| (1) | 18,068,480 | C>G | TP HighConf |
| (2) | 18,072,546 | G>C | TP HighConf |
| **(3)** | **18,086,020** | **G>A** | **落 HC BED 空隙（18,085,984–18,086,058）= truth-unevaluable**，非可由 SEQC2 判定的 FP |
| (4) | 18,096,269 | C>G | TP MedConf（ClairS 抓 T 標 FP；後接 20bp poly-T）|
| (5) | 18,099,697 | G>C | TP（原 LikelyFP）|
| (6) | 18,108,828 | C>G | TP HighConf |

🔴 **(3) 校正**：不是 SEQC2 FP，是 **out-of-HC 候選**。兩 basecaller tumor 皆 G>A PASS、matched normal 0% A（HKU tumor 29A/30G normal 0A/18G；DORADO tumor 21A/27G normal 0A/27G）、context `GAGACGG` 非 homopolymer → **強 tumor-specific somatic candidate，僅落 truth-gap 不可升 truth-confirmed**。

## A3 LOH 背景（✅ confirmed within benchmark）
- SEQC2 CNV BED：`chr2:16,146,119–22,100,000 = loh`，完整涵蓋本區。
- HP imbalance：HKU HP2-parent 232 / HP1-parent 4（**HP2 ~98.3%**）；DORADO HP2 279 / HP1 1（**~99.6%**）。
- **normal 6 點全 REF**（pos3 normal 100% G，兩 basecaller）。
- 🔴 「HP2 突變到不見/沒抽到 HP2 read」**否定**：HP2 是主體；HKU 17 條 reads 在 ≥2 已覆蓋 SNV 全 REF、≥3 點全 REF 有 4 條；缺跨 ≥4 點 all-REF read 是 **coverage/read-length 限制**，非 ancestral HP2 消失。

## A4 read linkage → ≥3 regional molecular states（互斥 0 違反）
**linkage 00/10/01/11**（unique primary, MAPQ≥20, BQ≥10；`11`=兩事件同 read）：

| 事件組合 | HKU | DORADO | 判讀 |
|---|---|---|---|
| (1)G vs (2)C | 18/0/0/**13** | 11/1/0/**6** | 同支強支持 |
| (1)G vs (3)A | 3/3/7/**0** | 跨距不足 | 互斥 |
| (2)C vs (3)A | 3/5/9/**0** | 0/1/0/**0** | 互斥 |
| (3)A vs (5)C | 8/1/0/**4** | 1/2/0/**2** | **(5)C 巢狀於 (3)A** |
| (3)A vs (4)non-ref | 0/9/13/**0** | 0/7/5/**1** | 互斥（DORADO 1 discordant）|
| (3)A vs (6)G | 0/3/3/**0** | 跨距不足 | 互斥 |
| (4)non-ref vs (6)G | 9/0/0/**10** | 1/0/0/**1** | 同 state（順序不可分）|
| (5)C vs (6)G | 3/7/11/**0** | 1/7/0/**0** | α-1 與 beta-like 互斥 |

→ **可支持 ≥3 regional molecular states**：
- **α**：(3)=A, (5)=G
- **α-1**：(3)=A, (5)=C（nested in α）
- **beta-like**：(3)=G, (4)=non-ref, (6)=G
- **+ beta-left block**：(1)=G, (2)=C（已觀察，但與 beta-right **無直接 genetic-link read**，僅以互斥 + 甲基 coherence 暫接）

## A5 甲基判別（cross-basecaller 複製 + 🔴 normal-HP ASM-control 拆解）
allele-CpG association（Mann-Whitney BH-FDR，10 候選 CpG）—— **這是本案例最關鍵的誠實拆解**：

| 類別 | CpG | 證據 | 結論 |
|---|---|---|---|
| ✅ **乾淨 tumor-associated**（normal HP **無** ASM + 兩 basecaller 複製 + 過 FDR）| **2.1, 2.2, 3.1, 3.2, 5.1, 5.2** | 例：3.2 HKU ALT 0.97/REF 0.14 q4.1e-5、DORADO 1.00/0.006 q2.4e-5；2.1 HKU q7.2e-4/DORADO 1.5e-5 | **真 tumor-acquired subclone-甲基訊號**（最強 3.1/3.2）|
| 🔴 **germline-ASM confounded**（normal HP1-vs-HP2 已顯著）| **3.3, 3.4, 3.5, corrected-4.1** | normal HP1/HP2：3.3=0.814/0.046、3.4=0.033/0.850、3.5=0.583/0.012、4.1=0.887/0.078（MW FDR≤0.05）；且 DORADO 3.4(p=0.060)/3.5(p=0.084) **未過 FDR** | **部分受既存 haplotype ASM 影響，不可當 subclone 形成後新 epimutation** |

🔴 **座標校正**：(4.1) 圖標 18,096,341（CpG 的 G）→ 正確 CpG C = **18,096,340**；文字 18,096,041 在 GRCh38 是 A 非 CpG（原腳本因此漏算，已重算）。
🟢 **tumor/normal differential 也跨 basecaller 複製**（2.1/3.1/5.1 等同方向過 FDR）—— 支持「tumor/normal 差異甲基位點」輸出目標，但 all-read 差異仍混 composition/LOH/ASM，須 allele-conditioned + normal-HP-conditioned 拆解。

> **核心訊息**：此區同時含 **normal haplotype ASM** 與 **tumor-associated allele-methylation structure**；**ISM 的價值是把兩者拆開**（normal-anchored control），而非把所有甲基差異都歸 subclone。subclone-甲基**真實存在但僅限無-germline-ASM 的 CpG 子集**。

## A6 可防守的候選樹（parsimony；註記嚴格）
```
all-REF / unsampled common ancestor   ← root = parsimony 假說, 非 spanning read 觀察
  ├── α: (3) G>A
  │      └── α-1: + (5) G>C            ← (3)→(5) direct nesting 支持
  └── beta-like program
        ├── left block:  (1) C>G + (2) G>C     ┐ 兩 block 間 = 互斥 + 甲基 coherence
        └── right block: (4) C>G（G/T/DEL jitter）+ (6) C>G   ┘ 連接, 缺直接 genetic bridge → 畫虛線
```
- α/beta-like 分叉：直接互斥支持。(3)→(5)：direct nesting。
- beta left/right：**methylation-bridged inference（虛線）**。(4)(6) 共現不可排序。

## A7 使用者推論 1-6 最終判定（verdict_02）
| 推論 | 判定 | 修正版 |
|---|---|---|
| 1 發生 LOH | **✅ 確認** | SEQC2 LOH BED + HP imbalance 一致 |
| 2 發生 subclone | **🟢 強支持、有界** | regional operational subclonal states；非完整 biological clone truth |
| 3 HP2 沒抽到/突變到不見 | **❌ 否定** | HP2 是主體；all-REF 短 linkage read 存在；缺完整 read = coverage |
| 4 (1)(2)(3)(6) 先突變 | **❌ 不成立/需重寫** | (3) 與 (1)(2)(6) 互斥，不能同 trunk；只可稱不同 early branch-defining candidates |
| 5 (5) 將群一群二分開 | **🟢 支持（provisional）** | (5)C 對 (3)A perfect nesting；direct support HKU 4/DORADO 2 |
| 6 (4) 多突變造成三群 | **❌ 否定** | 合併為 pos4-altered beta-like；G/T/DEL = homopolymer-uncertain（後接 20bp poly-T；HKU G4/T13/DEL10/REF21、DORADO G6/T15/DEL7/REF29）|

## A8 證據等級 + 誠實邊界
- **強**：LOH、somatic allele existence、局部 read linkage。
- **強支持但 operational/regional**：regional subclonal structure。
- **未證**：完整 clone identity / clone fraction / 完整 phylogeny；演化順序（僅 (3)→(5) 有 nesting）；**甲基造成突變（完全未證）**。
- **升 tier 前置（→L3+）**：(a) per-read 5mC partition permutation test；(b) normal-anchored cis/ASM control 分離 subclone-甲基 vs ASM-on-retained-LOH（**已部分跑：clean 子集 vs confounded 子集**）；(c) 第二 LOH 位點/第二樣本。
- ⚠ DORADO tagged BAM 曾被 cleanup（20260420 版）但 20260315 complete_matrix 版有 HP tag 仍在；CpG 2.2 資料方向與 user 手標相反（以實測為準）。

## A9 🔑 甲基分析有幫助到 subclone 分析嗎？（核心問題的誠實答案）

> **一句話**：**有幫助，但是「有界的輔助角色（corroborate / bridge / characterize），不是獨立判別器」，且幫助綁在 normal-anchored ASM-control 之上。** 這正是論文主軸「甲基 corroborate 非 validated discriminator」的乾淨實例。

**✅ 甲基實際幫到的 4 點（每點有 verified 證據）**：
1. **橋接基因上連不起來的 block**：beta-left (1)(2) 與 beta-right (4)(6) 相距 >36 kb，read 跨不到 → **無直接 genetic linkage**；是**甲基 coherence 把兩 block 接成同一 beta-like program**（A6 虛線）。→ 甲基補了基因連鎖的空缺。
2. **補 LongPhase-S fine-tag 分不開的地方**：fine tag（HP2 vs HP2-1）無法乾淨拆開 α/α-1/beta-like，且跨 basecaller 不穩（HKU ALT 多落 HP2-1、DORADO (3)A 落 HP2）→ 甲基提供**額外 regional molecular-state characterization**。
3. **指派缺 defining 突變的 read**：2 條 ref-genotype（C-G-G）read 帶**乾淨 β 甲基**（3.1 高/2.2 低）→ **epigenotype 追 lineage**，能標示沒蓋到突變的 read。
4. **跨 basecaller 複製 = 真訊號**：乾淨子集（2.1/2.2/3.1/3.2/5.1/5.2）HKU+DORADO 同方向過 FDR → 甲基判別非視覺印象、非單 basecaller 雜訊。

**⚠ 但「幫助」嚴格有界（必同時講，否則 overclaim）**：
- **subclone 結構是「基因」定義的**：≥3 molecular states 來自 somatic allele 互斥 linkage（00/10/01/11），**不是甲基**。→ 甲基是 corroborate/characterize，**非 primary discriminator**。
- **幫助綁在 ASM-control**：3.3/3.4/3.5/4.1 的甲基差其實是**既存 germline ASM 假象**（normal HP1-vs-HP2 已顯著）；**若不先用 normal-anchor 扣掉，naive 甲基分析會誤導**。→ 甲基「有用」的前提是 ISM 把 subclone 訊號與 ASM confound 拆開的能力（C2）。
- **甲基沒有做到的**：獨立重建 subclone tree、定 clone 數、排演化順序、證明甲基造成突變 —— 全部未證。

**論文用語**：甲基 = **lineage coherence + 區域 characterization + read-assignment**，conditioned on normal-anchored ASM-control；**不是** subclone 的獨立判別/重建工具。

---

# Part B — 方法實作構想基礎（ISM subclone 分析 pipeline 概念）

> **目的**：把 chr2 案例「手做能成立」的流程，萃取成**可系統化的 ISM subclone 分析方法構想**（對應 handoff **B2「個案→系統化+標準化判準」= 主軸升 validated 關鍵**，及 Ch3 Methods C3）。**這是構想基礎，非已實作的自動化模組**。

## B0 核心構想
> 在同一條長 read 上**聯合觀測 somatic allele + haplotag + native 5mC**，於 **LOH-constrained + somatic-haplotag-conditioned** 框架下，用**多點 read linkage** 定義 regional molecular states，並用 **normal-anchored 甲基拆解**把 tumor-acquired subclone 甲基與既存 germline ASM 分開 —— 輸出**有不確定性標記的 regional subclone characterization**（非完整 phylogeny）。

## B1 九步 pipeline（每步附案例落地的判準）
| 步 | 做什麼 | 判準/閾值（案例落地）| 護欄 |
|:--:|---|---|---|
| 1 | **萃取** read×sSNV×5mC 矩陣（pysam）| unique primary、MAPQ≥20、BQ≥10；MM/ML 解碼 5mC 映回 ref；數字先落 JSON | dedup（案例 340→280）|
| 2 | **LOH context** | SEQC2 CNV/LOH BED + HP imbalance（案例 98.3%/99.6%）；matched-normal ALT=0 | LOH 為前提；非 LOH 區另論 |
| 3 | **lineage by 互斥 linkage** | 00/10/01/11，要求跨位點 read；互斥（0 `11`）= 分叉，共現（高 `11`）= 同 state | 跨距不足 → 標「不可判」非互斥 |
| 4 | **artifact guard** | homopolymer context check（案例 pos4 後 20bp poly-T → G/T/DEL 合併為單一 altered state）| ONT homopolymer/strand-bias 必查；單 strand call 不算群 |
| 5 | **甲基判別 per allele anchor** | Mann-Whitney + **BH-FDR**；**要求 cross-basecaller 同方向 + 過 FDR** | 單 basecaller n=1 不可解；distal 未複製 → 排除 |
| 6 | 🔴 **normal-HP ASM control** | normal HP1-vs-HP2 顯著者（案例 3.3/3.4/3.5/4.1）標 confounded；normal HP 無差異者（2.1/2.2/3.1/3.2/5.x）= 乾淨 tumor-associated | **不可把所有甲基差歸 subclone** |
| 7 | **tumor/normal differential** cross-check | all-read tumor vs normal MW FDR；同方向跨 basecaller | all-read 差混 composition/LOH/ASM → 須 allele+normal-HP conditioned 再拆 |
| 8 | **建可防守樹** | parsimony；direct nesting 才排序（案例 (3)→(5)）；無直接 bridge → 虛線 methylation-bridged；root=hypothesis | 共現不可排序；VAF **不可**排 siblings（LOH 扭曲局部 VAF）|
| 9 | **誠實 tier** | regional/operational subclonal；標 L2/單位點；演化順序未證除非 nesting | 用 verdict_02 語、禁「confirmed N-subclone」 |

## B2 三條防呆（案例抓出的過度宣稱根因 → 系統化必內建）
1. **homopolymer artifact**：variant 後接 homopolymer → 多 allele call（G/T/DEL）是 ONT jitter 非多 subclone。**必查 context**（案例 pos4 20bp poly-T → 「5群」收斂為 ~3 lineage）。
2. **VAF 不可排序 siblings**：LOH 區局部 VAF 被扭曲（案例 pos5 局部 0.31 vs 全基因組 0.05）；互斥 siblings 不可由 VAF 排 parent/child → 順序只在 nested genotype 觀察到時成立。
3. **coverage 假象**：缺跨多點 all-REF read ≠ ancestral allele 消失（read-length 限制）；root 是 parsimony 推論非觀測。

## B3 ASM-control = ISM 對此案例的關鍵方法增值
- **問題**：LOH 已丟 HP1 → 無法用「丟失 allele」對照，「LOH 前既存 germline-ASM 被保留再由 somatic 切分」**無法形式排除**。
- **ISM 解法**：用 **matched-normal HP1-vs-HP2 ASM** 當 filter —— normal 已有 ASM 的 CpG（3.3/3.4/3.5/4.1）保守排除，只留 normal 無 ASM 的 CpG（3.1/3.2 等）當乾淨 tumor-acquired subclone 訊號。**這就是 C2（normal-anchored cis/ASM control）在 subclone 案例的具體展現** —— ISM 不只 flag 甲基差，還拆「subclone 形成新 epimutation」vs「既存 ASM confound」。

## B4 系統化路徑（→ handoff B2，主軸升 validated 關鍵）
- 目前：**個案手做 L2**（chr2 單位點 + BRCA2 單位點）。
- 升 validated 需：① B1 九步**自動化**成 ISM 模組（read→molecular-state assignment + ASM-control + tree-with-uncertainty）；② **標準化判準**（互斥閾值、FDR、cross-basecaller 複製、homopolymer mask、ASM-control gate）；③ 第二 LOH 位點/第二樣本複製（破單位點）；④ 第二 caller/pipeline（破 DORADO=同細胞株技術重現非獨立）。
- 對應 handoff：B2（系統化）+ B5（within-HP null/ASM-control，已部分跑）+ B7（跨樣本 ⭐4）+ B4（第二 caller）。

---

## 對應論文章節
- **Ch4 Results 4.6**：本案例 = C3 旗艦 demo（用 A1-A8 + A6 樹圖 + A5 ASM-control 表）。
- **Ch3 Methods（C3 + C2）**：B1 九步 pipeline + B3 ASM-control = subclone characterization 方法描述（標構想/proof-of-concept）。
- **Ch1 §1.4 / §1.5**：chr2 案例敘述須更新為本檔 A 部（≥3 lineage、ASM-control 子集、verdict_02 用語）—— 取代舊 2-locus（F=10.6/V=0.67）框架。
- **Ch5 Limitations**：A8 + B2/B4（單位點 L2、homopolymer/VAF/coverage 護欄、系統化=future）。
- 用語鐵則：「regional, LOH-constrained, somatic-haplotag-conditioned subclonal structure with cross-basecaller methylation coherence」。

## Provenance
- 案例數字（linkage 00/10/01/11、甲基 q、ASM-control normal HP1/HP2、homopolymer G/T/DEL、HP imbalance 98.3%/99.6%、pos3 29A/30G）全引自 `20260615_chr2_18M_subclone_independent_verdict_02.md` + `_assets/data/independent_audit.json`（機器輸出；該 session 產，本主軸 session memory §47-53 已 cross-AI 重算 byte-identical / 逐位吻合）。
- 最終裁決用語、推論 1-6、可防守樹採 verdict_02 §最終裁決/§8/§9。
- 本檔為 consolidation（不產新數字）；原驗證報告 + 教學 HTML（`docs/explain/04,05`）為來源，未修改。
