<!--
建立時間: 2026-06-15
報告類型: 單一 canonical 位點 subclone reconstruction 資料驗證（含對抗式獨立複核）
任務類型: A pilot (single-locus) + 驗證導向；scope = 單位點 chr2:18.07-18.11M, 2 basecaller (HKU_MOD baseline + DORADO cross-check)
狀態: in_progress
data_sources: docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/{hku_mod_matrix.json,dorado_matrix.json,hku_mod_analysis.txt,dorado_analysis.txt,methyl_discrimination_pos3.txt}; data/answer/SEQC2/*; data/ref/GRCh38_no_alt_analysis_set.fasta; /big8_disk/.../HCC1395_Tmode_tagged*.bam + HCC1395BL*.bam + ONT_Dorado/HCC1395{,BL}.bam
驗證方式: pysam read-level 萃取 → Read 讀回 → SEQC2 truth + 兩 basecaller 複製 → 4-agent 對抗式獨立複核（tree/ASM/artifact + synthesis）
-->

# HCC1395 chr2:18.07–18.11M Subclone Reconstruction — 資料驗證報告

> **partial-scope flag**：**單一 canonical 位點 × 單樣本 × 單 pipeline**（cross-basecaller 為主要 uplift，非獨立複製）。整體證據天花板 **L2**；推廣需第二樣本/第二 LOH 位點。

## L0 一眼結論（TL;DR）

使用者手動 subclone reconstruction 的**骨幹在實際資料中獲證並跨 basecaller 複製**；但對抗式複核抓出**三處需修正的過度宣稱**。

| Claim | 內容 | Tier | 結論 |
|---|---|:---:|---|
| **A** | 6 sSNV：5 TP + 1「FP」(vs SEQC2) | **L2** | ✅ tumor/normal allele 兩 basecaller 重現；TP/FP 標籤承襲 SEQC2（pos3 落 HC 空隙不可 truth 驗證）|
| **B** | read linkage 群組 | **L2** | ⚠️ 分叉乾淨（**0 違反**），但「**5 群」過報 → 實為 ~3 lineage**（pos4 G/T/DEL 是 poly-T homopolymer artifact）|
| **C** | 甲基把兩 lineage 分開 | **L2** | ⚠️ HKU 9 CpG 顯著；**扣 normal HP1-vs-HP2 germline-ASM 後**乾淨子集=**{2.1,2.2,3.1,3.2,5.1,5.2}**、**{3.3,3.4,3.5,4.1} 被既存 ASM confound**；clean+DORADO-FDR 最強=**{3.1,3.2}**（2026-06-15 第二 audit 取代舊「robust {…,3.3}」）|
| **D** | LOH + normal=REF + 無 REF read | **L2/部分撤回** | ✅ LOH(HP1 1.4%)+normal=REF 乾淨；❌「**無 ancestral REF read」是覆蓋假象**（C-G-G ref read 確實存在）|
| **E** | subclone tree 拓樸 + 順序 | **L2 骨幹 / L3 順序** | ✅ 骨幹（all-REF root→α/β 分叉, pos5 nest α）parsimony 強制；🔶 **突變順序 L3**（VAF 不可排序 siblings + 局部 VAF 被 LOH 扭曲）|
| **F** | pos3「FP」其實是真 somatic | **L2 最強** | ✅ 30/28 分子、mapQ60、雙股、tiled starts、兩 normal 0% → 真變異，僅落 HC 空隙 |

**最重要訊息（成立）**：LOH 區單親 haplotype 保留 + 兩 lineage 由 **somatic 變異**定義 → 甲基差異是 **somatic-genotype 內**的分群（**operational subclonal**，非 live germline ASM）。**這正是論文主軸要的 subclone-level 甲基 characterization 乾淨案例。**
**最重要 caveat**：兩個最強判別子各落 confound 區 — pos3 在 SEQC2 HC **空隙**（artifact-clean 但未 truth 驗證）、甲基 flip 的 distal CpG 在 DORADO **不複製**。視為 descriptive exception，勿推廣到 non-LOH 位點。

---

## L1 摘要

**為何**：使用者用 IGV 手動觀察推出 LOH→subclone→5 群重建樹，需對實際 BAM/VCF 核對 + DORADO 複製。
**怎麼做**：pysam 萃取每 read 在 6 sSNV 的 base(含 DEL) + HP/PS + MM/ML 解碼 5mC/5hmC（映射回參考）→ tumor+normal；對照 SEQC2 truth + HC.bed + GRCh38 ref；HKU_MOD(5mCG_5hmCG,haplotagged) baseline + DORADO(raw,5mC) cross-check；數字先落 JSON→讀回→才撰寫(§13.0)；最後 **4-agent 對抗式獨立複核**（tree-topology / ASM-confound / artifact-stats + synthesis）逐項挑戰。
**結果**：骨幹獲證、三處過度宣稱修正，整體 L2。圖 `_assets/data/verification_figure.png`。

---

## L2 各 Claim 細節（含對抗修正）

### Claim A — 6 sSNV TP/FP（L2）
5 TP(SEQC2: (1)(2)(6) HighConf Strong、(4) MedConf、(5) HighConf-原 LikelyFP) + 1「FP」((3) 落 HC.bed 空隙 18085984–18086058，未在 truth)。ref/alt 全符合。(4) ClairS 抓 T allele 標 FP（truth 為 C>G）。normal 在 somatic-defining pos2/pos3 = **0/48、0/36**（兩 basecaller 一致）。**標籤承襲 SEQC2，非本地重新判定。**

### Claim B — read linkage（L2；**「5 群」修正為 ~3 lineage**）
- **分叉乾淨**：α(pos3=A) 與 β(pos1=G/pos2=C/pos6=G) **互斥，0 違反 read**（兩 basecaller）；pos1×pos3 / pos2×pos3 / pos3×pos6 alt-alt = 0。
- **🔴 修正**：「群3/4/5 = G/T/DEL = 3 subclone」**過報**。pos4(18096269=C) **緊接 20bp poly-T homopolymer**（`AGTC·[C]·TTTTTTTTTTTTTTTTTTTT`，GRCh38 ref 親驗）= ONT 最差錯誤情境。T/DEL/G reads 在**其他每個 marker 皆相同**（pos3=G、pos5=ref、pos6=G、同甲基 signature）→ **一個 pos4-altered β sublineage + homopolymer jitter，非 3 群**。pos4 變異本身真（tumor 非C 55%、normal 乾淨），但 T/DEL/G 三分是 artifact。
- **可信 lineage ≈ 3**：α(群1)、α+pos5C(群2)、β+pos4-altered(群3，4/5 併入)。

### Claim C — 甲基判別（L2；**robust 複製限 proximal 子集**）
- HKU_MOD：α vs β 在 **9 個 CpG 全顯著**（Mann-Whitney p 1e-2→1e-7；permutation 20k 對 3.1/3.2/2.2/3.4 perm-p≤1e-4）。
- **🔴 修正（2026-06-15 第二 audit ASM 控制取代舊結論）**：跑 normal HP1-vs-HP2 甲基（本 session 獨立 code path 復現，逐位吻合）發現 **3.3/3.4/3.5/4.1 在 normal 已有顯著 germline ASM**（Δ0.57–0.82，FDR≤0.03）→ 這些 CpG 的 tumor α/β 差異**被既存 ASM confound，不可當 subclone 新生甲基**。
- **乾淨子集（normal 無 ASM）= {2.1, 2.2, 3.1, 3.2, 5.1, 5.2}**；其中 **clean + DORADO-FDR 通過 = {3.1, 3.2}（最強、最可防守的 subclone-甲基證據）**。舊「robust {…,3.3}」作廢（3.3 是 ASM-confounded）。
- **4.1 座標校正**：原報 18096041 是 A（非 CpG），真 4.1 CpG=**18096340**；校正後 4.1 也屬 ASM-confound 組。
- 0 甲基-lineage 矛盾 read；2 條 ref-genotype(C-G-G) read 帶**乾淨 β 甲基**（3.1 高/2.2 低）→ epigenotype 追 lineage（強化非弱化）。

### Claim D — LOH + normal=REF（L2，**「無 REF read」撤回**）
- ✅ **LOH**：HP1=4/280=**1.4%**；HP2 與 HP2-1 共用單一 PS=16735547（longphase-S sub-haplotype 切分，非第二親本 allele）。
- ✅ **normal=REF**：6 點 normal 全 REF（pos3 normal 100% G，兩 basecaller）。
- **🔴 撤回**：「cov≥4 reads 0 條 all-REF → HP2 突變到不見」是**覆蓋假象** — 僅 ~10(HKU)/0(DORADO) reads 覆蓋 ≥4 點；放寬至 ≥2 乾淨點有 22(HKU)/17(DORADO) all-REF；**ancestral C-G-G ref reads 確實存在**(895e6dcb/f6530558/0b2ed59d)。root 為 parsimony 推論，非觀測 read。

### Claim F — pos3「FP」實為真 somatic（L2，最強單一結果）
pos3 A：HKU 30 分子 / DORADO 28 分子、**全 mapQ60、0 低品質**、雙股(HKU 20+/10- p=0.10 ns；DORADO 13+/15- p=0.85)、ref_start tile 49.7kb 跨 29 起點(非 clipped pile)、median span 17.3kb（真長 read）、兩 normal 0/36 & 0/31。pos3 context `GAGACGG` 非 homopolymer。→ **真 tumor-specific somatic，僅落 SEQC2 HC 空隙被排除（truth-gap artifact，非真 FP）。**

---

## L3 對抗式複核 — claim E 與三大 INFERENCE

### subclone tree（骨幹 L2 / 順序 L3）
- **骨幹獲證**：reciprocal exclusivity（α reads 在 pos1/2/6 為 ref；β reads 在 pos3 為 ref；互不含對方 defining 突變）→ parsimony **強制** all-REF root + α/β sibling 分叉；pos5(C) nest 在 α 內(群2)。
- **順序 L3 + 🔴 修正**：「VAF 排序突變」**撤回** — 局部 VAF pos5=**0.31**（非 SEQC2 全基因組 0.05，LOH 扭曲局部 VAF）；α/β 是 **siblings 不可由 VAF 排 parent/child**。pos5-subclone 僅 4 reads（provisional）。

### 任何論文文字必標為「推論」(非事實) 的三點
1. **「5 群/5 subclone」count** → 改「**≥3 distinguishable lineages**」；pos4 三分若保留須標 ONT-homopolymer-uncertain。
2. **突變順序 / 「無 ancestral REF read」(two-hit order)** → 覆蓋太稀(≤10 reads 跨≥4點) + ancestral read 存在；root 是 parsimony 推論非觀測。
3. **「subclonal 而非 germline ASM」** → operational 為真（單一保留 allele → live 雙親 ASM 結構上不可能 + flip 跨 HP tag + HP-free DORADO 複製）。**🆕 ASM 控制已跑（2026-06-15 第二 audit）**：normal HP1-vs-HP2 顯示 **{3.3,3.4,3.5,4.1} 有既存 germline ASM、{2.1,2.2,3.1,3.2,5.1,5.2} 無** → 部分解 confound：乾淨 subclone-甲基訊號**收斂到無-ASM 子集（最強 3.1/3.2）**。nuance：tumor 已 LOH 丟 HP1，故 confound 機制是「這些 CpG epigenetically labile」而非「簡單繼承 germline ASM」。

### 升 tier 前置（→ L3+）
(a) per-read 5mC partition 的 binomial/permutation test；(b) ~~normal-anchored cis/ASM control~~ **✅ 已完成（第二 audit；見 `independent_audit.{json,md}` + `..._independent_verdict_02.md`）**；(c) 第二 LOH 位點或第二樣本複製（仍待）。

---

## 限制
1. 單位點 × 單樣本 × 單 pipeline；cross-basecaller 為主要 uplift（同細胞株技術重現，非獨立 biological replicate）→ 天花板 L2。
2. ~~DORADO tagged BAM 已清理~~ **校正：20260420 版被刪，但 20260315 complete_matrix 版仍存在（有 HP tag）；第二 audit 已用它做 HP-aware DORADO cross-check**。
3. tagged BAM 重複記錄 340→280 unique（dedup 後分析）。
4. 與 user 手標不符：**CpG 2.2** 資料 α=H（手標 L），方向相反。

> **🔗 2026-06-15 對齊註**：本報告（第一輪驗證）之甲基 robust set、ASM、DORADO tagged 三處，已被更嚴謹的**第二輪獨立 audit** 取代/補強，並經本 session 親驗（重跑 byte-identical + 我方 code path 復現 normal ASM 逐位吻合）。權威結論以 `InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_independent_verdict_02.md` 為準；教學頁 04 Fig 4 已 data-bound 重生對齊。

## 資產
- 腳本：`_assets/scripts/{extract_locus_matrix,analyze_matrix,methyl_discrimination,make_verification_figure}.py`
- 資料：`_assets/data/{hku_mod,dorado}_matrix.json` + `*_analysis.txt` + `methyl_discrimination_pos3.txt`
- 圖：`_assets/data/verification_figure.png`
- 對抗複核 verdicts：workflow `wf_668e592a-84a`（tree/ASM/artifact + synthesis）
