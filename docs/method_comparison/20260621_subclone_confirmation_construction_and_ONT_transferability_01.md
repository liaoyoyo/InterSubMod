<!--
建立時間: 2026-06-21
報告類型: method_comparison — subclone 如何被「確認」與「建構」+ 小區域 read-level sweep + 單細胞→ONT 移植裁決 + Tarabichi 嚴謹度學習
任務類型: D handoff — 論文 framing 知識留底（durable reference）
build_branch: research/subclonal-reconstruction-202606
status: synthesis；外部 = L3 agent-fetched（帶 PMID，投稿前過 /citation-verification）；ISM 內部數字標 tier+來源
data_sources:
  - workflow w83qhk507（6-agent：4 survey + 合成 + 對抗稽核 PASS）
  - 聚焦 agent a5c6317b（Tarabichi PMC7867630 全文 LEARN 萃取）
  - 聚焦 agent a77ae865（PCDH+Maura 全文核，進行中）
  - external_validation 庫（72 源）+ memories（見關聯）
provenance_note: ISM 內部數字溯源 memories/docs（verified）；外部 = agent WebFetch（L3）。🔴 PCDH(bioRxiv 2026.02.23.707349) + Maura 2024(PMID 39149342) 為載重反證但 full-text 403-blocked，僅 search-summary → 投稿前必過 /citation-verification（驗證進行中）。
adversarial_audit: workflow w83qhk507 evaluator = PASS（0 fabrication/0 overclaim/0 redline；transferability 邏輯 5/5 全對；2 minor 已折入）。
-->

# Subclone 如何被「確認」與「建構」+ 小區域 read-level sweep + 單細胞→ONT 移植

> **這份回答**：(1) 過去研究用哪些方法**確認** subclone 是真的、(2) 用哪些方法**建構** subclone 架構、(3) 有沒有人在**小區域用 read + 多證據**建/確認 subclone（ISM 精確 niche）、(4) **單細胞方法**深入 + 能否**移植到 bulk ONT**、(5) 可從 Tarabichi 實作指南**學什麼**。
> **可信度**：ISM 內部 = tier ⭐ + L1-L5（溯源 memories/docs）；外部 = **L3**（agent WebFetch + PMID，投稿前 `/citation-verification`）。🔴 兩個載重反證（PCDH/Maura）full-text 未取得（403），驗證進行中。

---

## §0 一頁總結

**四個裁決（先看這裡）**：

1. **「確認」subclone 的黃金標準 = single-cell 或 multi-region/longitudinal ground truth**；🔴 **單一 bulk 樣本無法「自我確認」subclone**，只能在 documented confounds 下 characterize。single-cell 是 terminal arbiter（genome-wide per-cell linkage → lineage 近 100%）。單一 bulk ONT 原生能跑的唯一確認族 = 統計顯著性（PERMANOVA+PERMDISP），但 tumor-only 無監督 = NEGATIVE/double-dip，僅條件於 a-priori 基因軸才 non-circular（~7.6%，且**只是與 unsupervised-clean 打平、優勢是合法 null 非偵測力**）。

2. **「建構」subclone 架構的五大家族都建出 ISM 不建的東西**：CCF clustering / phylogeny tree / allele-CNA / read-backed phasing / single-cell tree。ISM **消費**genetic 骨幹（haplotag+LOH+CN）但**不建 CCF/tree**；ISM 的「reconstruction」= **regional LOH-constrained same-haplotype partition（phase-block 尺度）**。

3. **🔴 小區域 read-level sweep = PARTIALLY DONE，但 ISM 三軸整合白地經全文核 SURVIVES INTACT（NARROWS 非 destroys）**：
   - **Foltz et al. 2024**（"SomaticHaplotype"，PMID 39149342，PMC11326269，**Li Ding lab**；🔴 先前誤標「Maura」）：**10X linked-read（非 ONT）骨髓瘤**，用成對 phased somatic mutation 的共享 barcode 測同-haplotype-vs-separate-subclone。**genetics-only / 零甲基 / 無結構檢定** → 佔「multi-sSNV co-occurrence as subclone 確認」軸但不碰甲基。whitespace LOW（= ISM 在其上加甲基軸的 genetic 半邊先驗）。
   - **PCDH**（Hackett…**Blundell** 2026，bioRxiv **10.64898**/2026.02.23.707349）："methylation barcodes in protocadherin cluster lineage tracing"。🔴 **更正事實**：**EM-seq + ONT 雙平台、targeted-capture**（非純 bulk-ONT WGS）；於 **PCDH 特製 hypervariable barcode locus**（非任意 somatic loci）；**single-cell WGS 驗證「未能證實」**（只確認 recapitulate driver-mutation VAF 推的 genetic clone 大小）。它證偽「read-level 甲基 lineage on long-read 尚未被嘗試」→ **此宣稱 FALSE 必丟棄**。
   - **但 PCDH 只佔 ISM 三軸的 1/3**（subclone target，且是 clonal-fraction/CH 非 somatic-driver-subclone）；四差異化**全 CONFIRMED**：(1) 甲基當 **DRIVER**（neutral epimutation clock）反 ISM 紅線；(2) **特製 barcode locus** 非任意；(3) **無 normal-cis**（philosophically 倒反 — 它要 neutral，ISM 要控掉）；(4) **無距離矩陣 PERMANOVA 結構檢定**。**ISM sanctioned 3-axis 整合 SURVIVES INTACT；唯一要收的措辭 = 別宣稱「first read-level methylation lineage on long-read」**。⚠ PCDH full-text 仍 403（abstract+snippet 級，single-cell 驗證待補全文）。

4. **單細胞→ONT 移植裁決**：genome-wide single-cell lineage **不能**移植到單一 bulk 樣本（觀測單元 = CELL vs READ，物理天花板非調參）；**能**移植的 = **phase block 內的 regional read-level 單分子 co-occurrence**，由 genetic anchor 橋到 kb-Mb haplotype 尺度 = **正是 ISM 的 regional proof-of-concept scope**。

**一句話**：ISM 的位置經這輪深掘後**更精確、也更窄**——它是「在 genetic 骨幹上、於任意 somatic locus、用 read 層次整合（無監督結構檢定 + normal-cis + subclone target + 甲基佐證）做 **regional** subclone characterization」的 **novel INTEGRATION**；不是無人有的能力，**且 2026 的 PCDH 已逼近**，故 framing 必須誠實收緊。

---

## §1 如何「確認」subclone 是真的（驗證方法論）

> 六條 escalating 確認路線；可行性由單一根本事實決定：single-cell 觀測 CELL（跨 locus linkage → lineage 近 100%，terminal arbiter）vs bulk long-read 觀測 READ（局部、無跨 locus/phase-set linkage）→ 物理上只能 regional。

| 路線 | 怎麼確認（黃金標準）| regime / 粒度 | 單一 bulk ONT 可做？ | 對映 ISM gate |
|---|---|---|---|---|
| **① single-cell ground truth** | bulk 推的 subclone 是否在 single-cell（scDNA/scWGBS）重現為離散細胞群；或物理分選 sublines（每株 = 真值）。clone-recovery F1≥0.7（>4% abundance）| single-cell / cell-level（最高解析）| ❌ **NO**（需上游細胞分選）；**terminal arbiter** | single-cell gate（bulk read-level subclone NOT-established 的唯一出路）|
| **② multi-region**（TRACERx）| 跨空間 biopsy 的 subclone 突變集是否重現（全區=trunk / 限區=branch）；供 pigeonhole 約束 | bulk 多 biopsy / mutation-cluster | ⚠ PARTIAL（方法可移植 ONT，但單樣本不可得）| **G-A**（ISM 替代 = cross-SAMPLE 6 cell-line 復現，非 multi-region）|
| **③ longitudinal**（pre/post/relapse/ctDNA）| 低-CCF 族群的時間動態（治療下擴張/MRD/relapse 再現）| bulk 多時點 | ❌ NO（單時點）| 無 ISM mapping（future work；同單樣本限制）|
| **④ 模擬/benchmark 真值** | DREAM SMC-Het（SC1 prevalence / SC2 co-clustering / SC3 tree；**演算法選擇主導、無通用贏家**）；BAMSurgeon spike-in（短讀）；cell-line 真值（SEQC2 HCC1395/COLO829）| 模擬 / cell-line | ⚠ PARTIAL（可建甲基域模擬）| 封頂 single-pipeline tier（Salcedo 2024 證連 genetic SR 都 algorithm-sensitive）|
| **⑤ 統計信心**（單樣本原生唯一可做）| posterior / bootstrap 穩定 / PAC（Şenbabaoğlu）/ M3C correlation-null / PERMANOVA+PERMDISP | 單樣本 | ✅ **YES（唯一原生）** | ISM PERMANOVA+PERMDISP；🔴 但 tumor-only 無監督 NEGATIVE（double-dip，noise≈structure 83-100%）→ 僅 a-priori 軸 non-circular |
| **⑥ 正交多組學 concordance** | SNV-tree vs CNA-tree vs methylation-tree 一致 | 多組學 | ⚠ regional（phase-block 內 chr2:18M）| ISM same-pipeline cross-check（非 orthogonal co-validation）|

**🔴 結論**：subclone 確認的黃金標準 = single-cell / multi-region / longitudinal **ground truth**，ISM 單一 bulk 樣本**皆不具備**（只有統計顯著性 + regional 正交佐證）→ 這正當化 ISM **⭐3 characterization / proof-of-concept**，**不是弱點而是領域共識的物理邊界**（Tarabichi guide 明文背書，見 §5）。

---

## §2 如何「建構」subclone 架構（五大家族）

> 每家族建出**完整架構**，皆 ISM 不建。ISM 消費 genetic 骨幹但只做 regional partition。

| 家族 | 用什麼輸入 | 輸出 | 粒度 | vs ISM |
|---|---|---|---|---|
| **① CCF clustering**（PyClone-VI/DPClust/SciClone/MOBSTER/CliP）| VAF→CCF（需 purity + allele-CN + multiplicity）；零甲基/零 phasing | CCF clusters | genome-wide per-mutation | ISM 不做 VAF→CCF |
| **② phylogeny tree**（PhyloWGS/Canopy/CITUP/LICHeE/B-SCITE）| +infinite-sites/pigeonhole | clonal **TREE** | genome-wide | **定義 "reconstruction"** = ISM scope 護欄；ISM 不建樹 |
| **③ allele-CNA**（ASCAT/Battenberg/TITAN/MEDICC2）| BAF+logR | allele-specific CN（+CN-phylogeny）| segment | ISM 上游 context；**亦 LOH-unmask confound 根源**（82-91%）→ G-B 未定 |
| **④ read-backed phasing**（longphase-S/Wakhan/Severus/MethPhaser）| long read + germline SNP | **haplotype** + somatic 歸屬 | phase-block | **建 haplotype 非 subclone**；longphase-S = ISM 骨幹 + regional 物理限制來源（跨 phase-set 無連結）|
| **⑤ single-cell tree**（SCITE/SCICoNE DNA；Epiclomal/MethylTree/Sgootr 甲基）| per-cell genome-wide profile | genome-wide lineage tree | cell-level | 真重建 + 唯一獨立真值，但**不同 regime**（觀測 CELL）|

> ⚠ prior-art 提醒：**LongPhase modcall + MethPhaser 已用甲基輔助 phasing** → ISM **不可**宣稱「甲基當 phasing marker」是首創；ISM 創新只在**整合**（見 §3）。

---

## §3 小區域 / read-level / 單分子多證據 subclone —— sweep 裁決

> **問題**：有沒有人在**小局部區域**、用**個別 reads + 多種共定位證據**（somatic SNV linkage + methylation + phasing on the SAME molecules）建/確認 subclone？這是 ISM 精確 niche。

### 🔑 總裁決：**PARTIALLY DONE**（2024-26 缺口大幅收窄，framing 必須收緊）— ⚠ 待 PCDH/Maura 全文核 finalize

**軸 A — 單分子 multi-sSNV co-occurrence 作為 subclone 確認（= chr2:18M 邏輯）→ 已嚴格 DONE，但不同 regime**
- **Foltz et al. 2024（"SomaticHaplotype"，PMID 39149342 / PMC11326269，Li Ding lab；🔴 先前誤標「Maura」，作者更正）**：**10X linked-read** 骨髓瘤 + matched normal，把 somatic mutation 指派 haplotype，分析**成對 phased somatic mutation 間共享 barcode** 測是否同分子/同 subclone（如 NRAS Q61K vs G12 無 ALT-ALT barcode → 獨立 subclone）。**full-text PMC 已核**。
- vs ISM：佔 (A) 軸但**非 ISM 全 niche** —— linked-read **非原生 ONT**；**genome-wide 非 regional**；**零甲基**；descriptive barcode count **無結構檢定**；無 read-level 甲基 cis-test。whitespace LOW（= ISM 加甲基軸的 genetic 半邊先驗）。

**軸 B — read-level 甲基 + SNV 整合（cvlr/ASMS/qFDRP/CpelNano；ROCIT/MethylBERT）→ 最近鄰，無一佔 niche**
- cvlr/ASMS/qFDRP/CpelNano：read×position 甲基聚類偵測 region 內 ASM（**皆有 randomization**：cvlr 1000-random、ASMS 5000-perm、CpelNano permutation）；但 **methylation-only、無 multi-sSNV co-occurrence、無 normal-baseline cis-test（ASM≠subclone）、target=germline ASM**。
- ROCIT：supervised/binary/無 normal/**0/3 命中 ISM 軸**（PacBio-demonstrated）。MethylBERT：supervised deconvolution。
- 🔴 cvlr/ASMS/CpelNano/MethylBERT **皆 ONT-capable 且皆有 randomization** → **禁**稱短讀或缺檢定；ISM delta = **INTEGRATION**。

**軸 C — 🔴 最強反證（headline contrast paper，Tier A，已全文核 abstract+snippet 級）：PCDH 甲基-barcode lineage tracing**
- **Hackett…Blundell 2026（bioRxiv 10.64898/2026.02.23.707349）"In vivo lineage tracing across human tissues using methylation barcodes in the protocadherin gene cluster"**：🔴 **EM-seq + ONT 雙平台、targeted-capture**（非純 bulk-ONT WGS）→ haplotype-phased → PCDH 簇 10-CpG **二元 epiallele barcode（沿單 read/read-pair）** → epiallele 矩陣**階層聚類 + phylogenetic inference 重建 subclonal phylogeny**（MPN：DNMT3A ~84% 主 clone + subclonal JAK2/CBL）。🔴 **single-cell WGS 驗證「未能證實」**（agent 只確認 recapitulate driver-mutation VAF 推的 genetic clone-fraction，n=37 windows）。**尺度近 ISM chr2:18M（小區域+單分子+read-level 甲基）。**
- **它證偽**「read-level 甲基 lineage on long-read 尚未被嘗試」→ **該宣稱 FALSE，必須丟棄**。
- **但只佔 ISM 三軸的 1/3，四差異化全 CONFIRMED（agent 全文核）**：(1) **甲基是 DRIVER**（neutral epimutation clock）= 與 ISM 紅線**相反**；(2) **特製 barcode locus**（PCDH hypervariability，不 generalize 到他 locus）非任意 somatic-SNV-spanning loci；(3) **無 normal-baseline somatic cis-test**（philosophically 倒反：它**要** neutral，ISM **要控掉**它）；(4) **無距離矩陣 PERMANOVA 結構存在檢定**（用 clustering→phylogeny + recapitulation 驗證非 hypothesis test）。

**三軸合判**：每個元件、甚至 regional read-level 甲基 subclone reconstruction 現已存在 → **PARTIALLY DONE**；但**無任何方法同時組合 ISM 全部特定元件於任意 somatic driver locus**。兩佔位者各自反轉/省略 ISM 核心（Foltz = genetics-only linked-read；PCDH = 甲基-driver/barcode/無 normal-cis/無結構檢定）→ **ISM sanctioned 3-axis 整合 SURVIVES INTACT（whitespace_threat PCDH=MEDIUM-LOW NARROWS / Foltz=LOW）**。

> 🔴 **paper 行動（裁決：minor sharpening 非 pivot）**：(1) 必引 Foltz 2024 + PCDH-2026 為 must-cite must-distinguish 鄰居（PCDH = headline contrast）；(2) reframe delta = somatic-haplotag-DRIVEN + methylation-as-corroborator + normal-anchored cis + 距離矩陣結構檢定，於**任意** driver locus，明確對比 PCDH 的「methylation-neutral-clock-at-special-barcode」；(3) **停止**主張「first read-level methylation lineage on long-read」（PCDH 在特製 locus 做過）。⚠ PCDH full-text 仍 403（single-cell 驗證待補全文確認）。

---

## §4 單細胞方法深入 + ONT 移植裁決

### PART A — single-cell 為何在 lineage 上有效
| 子類 | 觀測單元 | 為何有效 | 對 ISM |
|---|---|---|---|
| scDNA CNV tree（10x/SCICoNE/CHISEL）| per-cell genome-wide CN | 每細胞 = genome-wide 觀測，有跨 locus linkage | 上游 context（ISM 把 CN held-const）|
| scSNV tree（SCITE/COMPASS/ScisTree）| per-cell mutation 有無 | mutation 跨 locus co-occurrence | **概念橋**：ISM chr2:18M = reads 跨 ≥2 sSNV 給 **LOCAL** 單分子 co-occurrence（phase-block 內 micro-tree）|
| single-cell 甲基（Epiclomal/MethylTree~100%/EPI-Clone/Sgootr/Gaiti）| per-cell genome-wide epimutation barcode | ~100% 因每細胞攜 genome-wide barcode 物理連結 | framing tension；🔴 Sgootr 證 distance-甲基重建**非新穎** → ISM 窄化為 somatic-haplotag-conditioned + normal-cis |
| single-cell multi-omics（scNMT）| 同細胞多層 | 所有層細胞內物理連結 | DIFFERENT-REGIME aspiration（ISM read = 細胞的 tiny-bulk 類比，僅 read span 內）|

### PART B — 5 機制移植表
| # | 機制 | 移植 | 精確限制 |
|---|---|:--:|---|
| **B1** read-as-pseudo-cell（cvlr/ASMS/qFDRP）| ✅ **YES 但 LOCAL ONLY** | read 跨多 CpG/SNV → 真單分子 co-occurrence；但無跨 locus linkage、bulk 無細胞身分 → **regional cluster 永非 genome-wide lineage**（ISM regional 天花板基礎；引擎本體）|
| **B2** haplotype-block as bridge（longphase-S+LOH+CN）| ⚠ **PARTIAL** | 延伸 linkage 到 phase-block kb-Mb，**stop 在 block 邊界**（「兩 block 間無連結資訊」MethPhaser 38909018）→ 跨 phase-set lineage 在單樣本**物理不可能**；anchor 買到 **REGIONAL = 正是 ISM scope** |
| **B3** atlas/supervised deconvolution（MethylBERT/ROCIT）| ⚠ YES 但 supervised | reference-dependent；輸出 fraction/origin 非 de-novo subclone；**不分 cis-ASM vs subclone**；與 ISM 無監督 normal-anchored 正交 |
| **B4** 物理預分選 long-read 單細胞（Lee/Liu-Goretsky 2025）| ❌ **NO（sidesteps）** | 需濕實驗預分選恢復 per-clone linkage；一旦預分選**就不再是 bulk subclone reconstruction**——驗證而非跨越觀測單元差異 |
| **B5** scWGBS lineage 演算法跑 bulk reads（Epiclomal/MethylTree/Sgootr）| ❌ **NO（categorically）** | 每 row 需 genome-wide per-cell 甲基 profile；bulk read 是 LOCAL/非細胞/跨 phase-set 不可連結 → 餵 reads-as-rows 產 **region 內 cluster 非 lineage tree**（ISM regional 天花板正式陳述）|

### 🔑 CORE VERDICT
**single-cell lineage 甲基有效因觀測「細胞」（genome-wide per-cell linkage）；bulk ONT 觀測「read」（局部、無跨 locus 細胞 linkage）→ genome-wide lineage 無法移植到單一 bulk 樣本。能移植的 = phase block 內 regional read-level 單分子 co-occurrence，由 genetic anchor（haplotag+LOH+CN）橋到 kb-Mb = 正是 ISM regional proof-of-concept scope。** 天花板是**物理的**（跨 phase-set 不可連結），非調參/資料量限制。

---

## §5 ISM 可學習的嚴謹度依據 — Tarabichi 2021 practical guide（PMC7867630，全文核 2026-06-21）

> 用戶指定參考。Tarabichi 2021（Nat Methods）= 領域定義「嚴謹亞克隆重建」權威。**能且應參考**，但它是**純 DNA 零甲基** → 學的是**嚴謹度框架與 pitfall**非甲基方法。每點標 evidence-tier（[FT]=全文逐字 / [SIB]=姊妹篇 Salcedo2024）。

**可採用的具體實作（按 reviewer-防守影響排序）**：
- **FIX-1 [FT] beta-binomial over-dispersion**：guide *"incorporate a Binomial or Beta-binomial noise model"*（Gaussian 不行）→ ISM per-CpG **Fisher（binomial-family，over-dispersion 未修）改 beta-binomial**；guide = canonical citation。最高槓桿（已知缺口 + 權威依據）。
- **FIX-2 [FT] 群數紀律**：guide 警告「不當 noise model → 高/低估 cluster number」→ ISM 手動 silhouette 改報 **k 穩定性分佈（bootstrap/clusterboot）**+ 每「≥2 群」配 PERMDISP（已做）。
- **LEARN-1 [FT] uncertainty reporting**：guide *"report uncertainty"* → 輸出 **partition + 信心帶 + 替代解**，對齊 ISM keep/remove 5 類 + 低信心子集。
- **LEARN-2 [FT] NRPCC-類覆蓋下限**：guide *"NRPCC>10 ≈ 40× diploid 50% purity"* → ISM 報 per-region 有效深度/reads-per-haplotype，`C_min=3`/`±5000` frame 成甲基域解析度下限。

**🔴 最重要：guide 明文背書 ISM「單樣本 = characterization 非 confirmation」framing**：
> *"tree inference should generally be restricted to multi-sample studies… subclones identified in single samples are only **weakly informative of the underlying tree**."* + *"**illusion of clonality**"*（single-sample 即使高深度也低估 subclone 數）。

→ ISM regional partition 落在「single sample 能 identify subclone」側（防守側），**前提絕不稱輸出為 tree/phylogeny** → 正當化 **⭐3 ceiling + G-A/V-1 gate**。

**🔴 引用紅線**：(1) 不可宣稱 ISM 做 guide 的 full reconstruction（CCF+multiplicity+phylogeny）；(2) guide **零甲基** → ISM 甲基層對它 orthogonal/unevaluated，禁暗示背書甲基驅動；(3) **SMC-Het 具體數字（NRPCC≥32/31 演算法）在姊妹篇 Salcedo2024 [SIB] 非本 guide**；(4) guide **不涵蓋 ONT/long-read/mappability**。（卡 §6 已升 v2，§7 全文核已解除；⚠ supplementary PDF 404 未取得。）

---

## §6 對 ISM 的意涵 + 修正 + open questions

### ISM 位置（更精確、更窄）
ISM = 「在 genetic 骨幹上、於**任意** somatic locus、用 read 層次**整合**（無監督 read×read 距離矩陣結構檢定 + normal-baseline somatic cis-test + somatic-subclone target + 甲基佐證）做 **regional** subclone characterization」的 **novel INTEGRATION**。reconstruction 由 **somatic haplotagging 驅動**，甲基 = corroborate（germline-haplotype 層，非 driver）。⭐3 單樣本 / G-A 後 ⭐4。**不佔無人有的能力；2026 PCDH 已逼近** → framing 必須誠實收緊。

### 兩個 evaluator 修正（已折入）
- **chr2:18M tier 用詞**：非無限定「PRIMARY-verified」→ 應寫「**⭐3 ceiling，證據級 L2（single-locus×single-sample×single-pipeline，cross-basecaller 為主 uplift）**」；🔴 locus **18,086,020 落 SEQC2 HighConf 空隙 = unevaluable（無 locus-level 外部真值）**。
- **a-priori 7.6% 口徑**：非偵測力優勢 —— **與 unsupervised-clean（15.3%/6.8%）統計打平**；a-priori 優勢 = **合法 permutation null / collider-免疫 / 可解釋**，非更強偵測器。禁寫「more discriminative」。

### Open questions（決定論文強度）
- 🔴 **G-B（最關鍵）**：subclone-甲基 somatic-specific vs germline-allelic UNDETERMINED（受 LOH-unmask 82-91% 壓制）。
- 🔴 **PCDH + Foltz 全文核（2026-06-21 完成，已建卡）**：裁決 = ISM 三軸整合 **SURVIVES INTACT，minor sharpening 非 pivot**（PCDH 只佔 1/3 軸、四差異化全 CONFIRMED）。⚠ **PCDH full-text 仍 403**（abstract+snippet 級）→ 投稿前須補全文確認 single-cell 驗證與「無 PERMANOVA 結構檢定」。Foltz = PMC full-text 已核。唯一框架調整 = 丟棄「first read-level methylation lineage on long-read」措辭。
- **HD-1**（phasing 循環）/ **G-A**（COLO829 缺 + SPOF，封頂 ⭐3）/ **V-1**（cis 1/816 under-tested）/ longitudinal 確認（future work）。
- **方法借鏡待落地**：beta-binomial（Fisher 必修）+ M3C correlation-preserving null（需 binary/Bernoulli 版本）。
- **library 去重**：「Liu-Goretsky 2025 PMID 40950124」=「lee-2025-melanoma」card 同一篇（slug 沿用舊誤稱，display_name 已正確）。

---

*關聯 doc：`20260621_clone_subclone_reconstruction_landscape_and_ism_feasibility_01.md`（6 學派 + 裁決）· `20260620_somatic_locus_methylation_combination_enumeration_01.md`（組合判讀）· `20260619_subclone_analysis_interpretation_full_framework_01.md`（why-hard）。*
*關聯 memory：`project_subclone_confirmation_construction_ont_transferability`（本 doc 摘要）。外部庫：`/big7_disk/liaoyoyo2001/external_validation/`（72 源 + PCDH/Maura 待建卡）。*
