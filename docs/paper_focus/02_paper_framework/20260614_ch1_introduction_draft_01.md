<!--
建立時間: 2026-06-14
狀態: Ch1 緒論草稿 v1 (依 20260614 主軸定稿 + 五桶報告 + external_validation 庫展開; 待用戶審 + citation-verify)
報告類型: thesis_chapter_draft (Ch1 Introduction)
受眾: 廖子游 · PI · 碩論口委 (繁中可懂)
framework: SCQA(背景→缺口→方法→貢獻) + claim-tier (Tier 4)
data_sources:
  - InterSubMod/docs/paper_focus/02_paper_framework/20260614_Intro主軸定稿_thesis_statement_01.md (主軸句定稿)
  - InterSubMod/docs/paper_focus/03_references/20260613_external_vs_ism_five_bucket_comparison_01.md (五桶+趨同+獨特改進)
  - /big7_disk/liaoyoyo2001/external_validation/ (背景論文 CONTEXT 卡, 一手親讀)
  - InterSubMod/docs/method_comparison/.../01_ism_method_spec_from_source.md (ISM 6 cores)
provenance_note: 背景論文引用對應 external_validation CONTEXT 卡(標 author-year); 內部數字 grep-verified(BRCA2 Δβ=-0.122 / chr2 F=10.6 V=0.67 / 1/816 / AUC=0.5049); 未 verified 者標 {{待補}}。本檔=草稿, 投稿前過 /citation-verification。§13.0 撰寫與分析分離。
-->
<!-- provenance-verified: 背景 citation 對應庫卡 slug; 內部數字引自 ISM spec L3 + flagship HTML + knowledge/11; {{待補}}=SBS1/Anderson2001。 -->

# 第一章 緒論（Introduction）— 草稿 v1

> **待辦**：投稿/口試前 (1) 過 `/citation-verification` 補正式 .bib；(2) 填 {{待補}}（SBS1/5mC-deamination、Anderson 2001 PERMANOVA）；(3) 5 張 axis3 新卡（cited_secondary tier C）升 primary 後鎖數字；(4) 可移入 `docs/thesis/chapters/ch1_introduction.md`（避免並行 session 互撞，本草稿先放 paper_focus）。

---

## 1.0 章節導引

腫瘤是由帶不同體細胞變異的**亞克隆（subclone）**組成的演化群體；解析其組成與演化是癌症基因體學核心問題。本論文提出 **InterSubMod（ISM）**——一套在 **Oxford Nanopore（ONT）長讀** tumor/normal paired 資料上、以 **somatic haplotagging 為骨幹**、在 **read 層次**對 somatic 變異鄰域甲基化做**有界 characterization** 的軟體，並以個案示範「甲基 + 多點體細胞變異」如何共同標示亞克隆結構。本章依序鋪陳：研究背景（§1.1）、現有方法的缺口（§1.2）、ISM 的方法與主張（§1.3）、旗艦個案（§1.4）、貢獻（§1.5）與誠實範圍（§1.6）。

---

## 1.1 研究背景（Background）

### 1.1.1 亞克隆重建以基因體訊號為骨幹

亞克隆重建（subclonal reconstruction）的主流方法以體細胞**單核苷酸變異的變異等位頻率（VAF）/ 癌細胞分數（CCF）** 加上**拷貝數（CNA）**校正，推斷亞克隆組成乃至演化樹：如 PyClone / PyClone-VI（Roth 2014；Gillis & Roth 2020）、PhyloWGS（Deshwar 2015）、Canopy（Jiang 2016）、SciClone（Miller 2014）、DPClust（Nik-Zainal 2012）、CliP，以及領域操作指南 Tarabichi（2021）與 DREAM benchmark（Salcedo 2024）。這些工具**完全不使用 DNA 甲基化**（本研究源碼層核實九個工具的程式碼皆零甲基處理），確立了「**重建責任在基因體訊號**」此一領域共識——本論文據此將 reconstruction 歸於 somatic haplotagging 骨幹，甲基化則定位為被 characterize 的有界附加層〔對應五桶報告 G7〕。

### 1.1.2 長讀定序與 somatic haplotagging

ONT 長讀可同時讀出**單分子的等位連鎖**與**鹼基修飾（5mC/5hmC）**。在癌症脈絡，somatic haplotagging 工具（longphase-S，2025）以 germline 母單倍型（HP1/HP2）為基礎，進一步重建 somatic 子單倍型（HP1-1/HP2-1）與無 germline 支撐的 HP3，把每條 read 標上「屬於哪個（子）單倍型」。ISM 即**消費**此 haplotype tag 與 somatic 變異標記，作為甲基 characterization 的框架。

### 1.1.3 甲基化與體細胞變異的關係（雙向，且多為背景機制）

DNA 甲基化與體細胞變異存在**雙向**生物學關聯，兩者皆為本研究的背景動機而非 ISM 的待證主張：

- **變異 → 甲基（cis，局部）**：體細胞點變異或結構變異可改變其**鄰近**位點甲基——Do（2020）在泛癌 ASM 框架內報告「rare somatic mutations affecting CTCF/TF motif classes track with de novo ASM」，Zhang（2019）報告 somatic SV breakpoint 改變鄰近 CGI 甲基（扣除 CNA 後 1286 顯著 probe）〔`somatic-variant-local-cis-methylation`〕。**關鍵**：Do 親述此 de-novo ASM「does not make a large numerical contribution」——即**somatic-cis 甲基真實存在但稀有**。此為 ISM cis-test 量測對象的**直接前例**。
- **變異 → 甲基（trans，全基因組）**：表觀修飾酶基因（IDH1/2、TET2、DNMT3A）的體細胞突變經代謝/酶活性改變**全基因組**甲基組（Turcan 2012；Figueroa 2010）〔`somatic-mutation-methylation-modifier`〕——此為 trans 機制，**非** ISM 量測的 locus-local cis。
- **甲基 → 變異（反向機制，本研究不量測方向、僅作背景與 confound 依據）**：5-甲基胞嘧啶（5mC）在 CpG 上隨機/酵素性**去胺基化**成 thymine，產生 T:G mismatch，若複製前未修復則固定為 **C>T transition**——這是 CpG 成為突變熱點、p53 熱點落在甲基化 CpG、以及泛癌 clock-like 突變特徵 **COSMIC SBS1** 的分子基礎（Pfeifer 2006, PMID 16570852；COSMIC SBS1 / Alexandrov 2020；BER 塑形見 Nat Commun 2024）〔`methylation-to-mutation-deamination`〕。⚠ 量化（5mC 去胺基化速率約未甲基 C 的 ~4×、CpG 突變約佔癌症突變負荷 ~30%）待一手全文核。<br>🔴 **方法學關鍵**：因「甲基→突變」（SBS1）與「突變→甲基」（cis，Do 2020）**兩方向並存**，ISM 觀察到的「somatic 變異 ↔ 鄰域甲基差」共定位**不可 naive 歸因任一方向**（例如落在 CpG 上的 somatic C>T 本身可能就是 deaminated mCpG，使「甲基差」與「該突變」成為同源事件而非 cis-調控效應）。**這正是 ISM 只報「關聯強度 + normal-anchored cis/drift 分解」、不主張因果或先後順序的依據**——雙向機制並存使單向歸因不成立，誠實定位反而站得更穩。
- **germline cis baseline**：germline 變異 → 鄰近甲基（cis-meQTL）真實且普遍（GoDMC，Min 2021：45% 測試位點有遺傳效應）〔`cis-meqtl-germline-atlas`〕，且 ASM 多由序列/單倍型決定（Onuchic 2018）〔`asm-sequence-haplotype-dependent`〕——**這正是 ISM 以 matched-normal 為錨點要扣除的 baseline**。
- **技術 confound**：雜合 SNP 落在 CpG 上可造成**表面 ASM 但非調控性 cis**（Shoemaker 2010）〔`cpg-snp-pseudo-asm-confound`〕——本研究的 limitation。

癌症中 ASM 整體升高已由 Do（2020，平均 5×、淋巴瘤 8.5×、GBM 9×、hypomethylation 主導）與 Martin-Trujillo（2017，印記 DMR 之 allelic 甲基異常 82–92% 由 CN 解釋）建立；且多源（cvlr Raineri 2023：8 個印記基因 Δm P<10⁻³；ASMS Raineri 2024：16/20 ICR；O'Neill 2024：4.46M aDMR；LongHap 2026 Fig3C）一致顯示**甲基訊號在 germline-haplotype 層級最強**〔五桶報告 R2〕。

### 1.1.4 read 層次甲基分群與既有距離法

在 read 層次找甲基子群的工具與 ISM 最接近：cvlr（Raineri 2023，EM 分群）、ASMS（Raineri 2024，read-clustering ASM）共用「ONT read×CpG 二元甲基矩陣 → read 子群」底材；qFDRP（Scherer 2020）使用與 ISM 預設 **NHD 距離數學同式**的 per-pair normalized-Hamming kernel；DAMEfinder（Orjuela 2020）為 tuple-ASM 概念先例。**這些單一技術元件皆為既有先例**——本論文的可防守新穎在於其**特定組合 + 條件化**（§1.5）。

---

## 1.2 研究缺口（Gap）

綜合 §1.1，存在一個**未被現有工具覆蓋的交集**：

1. **基因體 SR 工具不碰甲基**（§1.1.1），無法利用甲基提供的額外 read 層次結構。
2. **甲基-演化/lineage 工具不在 bulk 長讀 read-level regime**：Sgootr（2023，距離式甲基 lineage tree）、Gaiti（2019）、Epiclomal（2020）、MethylTree（2025）在**單細胞** regime 重建譜系；EVOFLUx（2025）、colorectal-growth（2026）以 **bulk-clock** 模型推演化參數〔五桶報告 TENSION〕——皆非「bulk 長讀、read 層次、為體細胞變異判別」的設定。
3. **甲基-輔助 phasing 仍 germline-only**：MethPhaser（2024）**逐字**將 cancer 列為 future work；癌症 phaser（Wakhan 2025、HiCancer 2021）源碼層零甲基——**cancer-subclone 甲基-phasing 是真實白地**〔五桶報告 G4〕。
4. **somatic-cis 甲基真實但稀有且被 confound**：Do（2020）的 de-novo ASM 稀有，CpG-SNP（Shoemaker 2010）與 CNA（Martin-Trujillo 2017）皆 confound——**需 matched-normal 錨點 + cis/drift 分解才能把乾淨 somatic-cis 個案從背景中分離**，而現有工具未在 read 層次 within-sample 做此事（Do 的 de-novo ASM 雖用 within-patient matched-normal，但為 bulk site/DMR-level bisulfite，非 read-level 距離結構）。

**缺口一句話**：缺一套在**單樣本 tumor/normal 長讀**資料上、**read 層次**、**normal-anchored**、把**體細胞變異與鄰域甲基聯合 characterize**（並能分解 cis vs 背景）的工具——既補基因體 SR 不碰甲基之處，又補甲基工具不在此 regime 之處。

---

## 1.3 我們的方法與主張（ISM）

> **（主軸句定稿 v1，2026-06-14）**

**InterSubMod（ISM）** 是**首個**在 **ONT 長讀 tumor/normal paired** 資料上、整合 **① read 層次甲基結構檢定 + ② read-level normal-anchored somatic cis-test + ③ somatic-subclone 目標** 的 read-level 軟體。它以 somatic haplotagging（longphase-S）的 **haplotype tag（HP1/HP2 + 子單倍型 HP1-1/HP2-1）與 somatic 變異為框架**，在 read 層次量化每個 somatic 位點鄰域的等位基因特異性甲基化（ASM）結構，輸出三類資訊：

1. **與 somatic 變異高度關聯**、在 read×read 距離空間呈顯著結構分離的甲基位點（PERMANOVA pseudo-F，999 perm）；
2. **tumor 內 haplotype 間 / tumor-normal 間有顯著差異**的甲基位點，並以 matched-normal 為錨點**分解**此差異屬真 focal cis 抑或 copy/subclone 背景（BRCA2 promoter 為旗艦示範）；
3. 甲基差異**與哪一軸（HP / ALLELE / copy）相關**的分解。

進一步，以「**多個 somatic 變異 + 鄰域甲基結構**」的組合對單一 region 的 read 做 **subclone 結構 characterization**，並以 **chr2:18M 個案示範**甲基如何與突變共同標示同一子克隆事件、據此推導時序假說（自動化分析流程與大規模驗證為延伸）。

**定位**：在**只有單樣本**的癌症資料上，提供 read 層次、normal-anchored、比 per-position 聚合更細緻的**甲基-變異聯合 characterization**——量化體細胞變異與甲基位點的**關聯強度**（不主張因果或先後順序），補足 germline-only 甲基-phasing 與 bulk subclone 重建未覆蓋的 **somatic / single-sample / read-level niche**。

---

## 1.4 旗艦個案（詳見第四章）

- **BRCA2 promoter（HCC1395 乳癌，TSG）**：ISM 同時偵到 HP-axis 甲基差（Δβ=−0.122）與 tumor-normal 甲基差，共定位於 BRCA2 promoter；normal-anchored 分解顯示此差**主要 track somatic 子克隆 tag（HP1-1）**、focal cis 效應 marginal（d_within=−0.023, perm p=0.024）。→ 示範 ISM **不只 flag「與突變共定位的甲基差」，更量化它屬子克隆結構 vs 直接 focal cis**（誠實：方向為 hypo，非典型 TSG promoter hyper-沉默）。
- **chr2:18M（HCC1395，30 kb，57 reads×203 CpG）**：longphase-S 分出 Normal=G → HP1 母本（19 reads 全 G）→ HP1-1 子克隆（36 reads，28 A/8 G）；位點 18,086,020 的甲基分群顯著（PERMANOVA F=10.6, p=0.01）且跟隨 HP/子克隆結構（Cramér's V=0.67），同區 somatic 位點 18,096,269 的甲基則跟隨 somatic 等位本身（V=0.97）。→ 兩種互補模式構成「甲基+突變共同標示子克隆」的機制基礎，**個案示範** workflow 並推導 Normal→HP1→HP1-1 時序假說（時序為數據支撐之推論，非正交確認；18,086,020 為標記子克隆的 het，非 caller-somatic）。

---

## 1.5 貢獻（Contributions）

1. **C1｜read 層次甲基-變異聯合 characterization**（novelty=**組合**）：在 HP tag + somatic 框架下，把 read×read 距離矩陣升格為分群引擎並以 PERMANOVA 檢定 read 在距離空間的**結構分離**。單一元件（NHD=qFDRP kernel、階層分群、PERMANOVA=Anderson 2001）皆為先例，貢獻在組合與條件化。
2. **C2｜read-level normal-anchored somatic cis/drift 分解**（增量獨特）：以 matched-normal 為錨點分解體細胞位點鄰域甲基差屬真 focal cis 抑或 copy/subclone 背景。🔴 **normal-anchor 本身非獨家**（Do 2020 de-novo ASM 亦用 within-patient matched-normal）；ISM 增量為 **read-level within-sample 距離結構 + LOH-HP 軸 + 整合**。
3. **C3｜單樣本 somatic-subclone 結構 characterization（個案示範）**：以多體細胞變異 + 鄰域甲基組合對單 region 做 read→subclone 結構，個案示範時序假說推導。
4. **科學貢獻（findings）**：(a) 嚴格證偽「甲基當體細胞變異 TP/FP filter」（LOSO 100% circularity、|Δβ| AUC=0.5049、strong-ASM 在 FP 富集 5× OR=0.194）——並提出外部未報告的 anti-discriminative novel 觀察〔五桶 G6〕；(b) within-sample somatic-controlled ASM 口徑（6/6 excess-over-null）；(c) 確立 cancer-subclone 甲基-assist phasing 白地（LOH-constrained phasing same-hap 93–99% 6/6）。
5. **誠實對齊**：本研究的 somatic-cis 甲基**稀有（816 可測中 1 乾淨 exemplar）與既有文獻一致**（Do 2020：de-novo ASM「不大數值貢獻」），是領域特徵而非偵測力不足。

---

## 1.6 誠實範圍與限制（必讀；詳見第五章）

1. **個案示範（proof-of-concept）**：chr2/BRCA2 示範主軸概念，**非已驗證的自動化 subclonal reconstruction**；自動化流程 + 大規模驗證為 future work。
2. **無 lineage tree / 演化參數推論**：reconstruction = read→subclone 指派/分群結構；時序 = 個案假說（標 hypothetical / consistent-with-data）。
3. **單樣本 single-pipeline ⭐3**：TP/FP 真值取自 tumor-only caller（非獨立 orthogonal ground truth）；6 cell line = 跨模型 reproducibility，非病人 cohort 泛化。
4. **甲基為 characterization 非 utility**：不改善 variant calling（filter-DEAD）；只量關聯強度，不主張因果/順序。
5. **方向與分解誠實**：BRCA2 hypo≠典型 TSG hyper；focal cis marginal（價值在分解層）。
6. **confound 清單**：CpG-SNP 假 ASM（Shoemaker 2010）、CNA（Martin-Trujillo 2017）、trans-modifier 通道（IDH/TET2，Turcan 2012）皆須在報告註記。

---

## 引用對照（投稿前過 /citation-verification）
- 庫內 verified（slug）：longphase-s / pyclone(-vi) / phylowgs / canopy / sciclone / dpclust / clip / tarabichi / salcedo / sgootr / gaiti / epiclomal / methyltree / evoflux / colorectal-growth / methphaser / wakhan / hicancer / longhap / cvlr / asms / qfdrp / damefinder / martin-trujillo / do-2020 / **somatic-variant-local-cis-methylation(Do2020+Zhang2019)** / **somatic-mutation-methylation-modifier(Turcan2012+Figueroa2010)** / **cis-meqtl-germline-atlas(GoDMC Min2021)** / **asm-sequence-haplotype-dependent(Onuchic2018)** / **cpg-snp-pseudo-asm-confound(Shoemaker2010)**。
- 🔴 **待補**：① SBS1 / 5mC-deamination（甲基→突變方向，庫內無；ISM 不測此方向，僅動機用）② Anderson 2001（PERMANOVA 方法源）③ 5 張 axis3 新卡為 cited_secondary tier C → 正式引用前升 primary（親讀全文）。

## Provenance
- 主軸句：`20260614_Intro主軸定稿_thesis_statement_01.md`（用戶 2026-06-14 逐句確認）。
- 背景論文：external_validation CONTEXT 卡（author-year 對應 slug，2026-06-14 親讀；5 張 axis3 甲基↔變異卡為當日新增 cited_secondary）。
- 內部數字：BRCA2 Δβ=−0.122 / d_within=−0.023（ISM spec L3）；chr2 F=10.6 / V=0.67 / 18096269 V=0.97（flagship HTML grep-verified）；AUC=0.5049 / OR=0.194 / 6/6 excess / 1/816（knowledge/11_external_literature/{05,07,08,10}）。
- {{待補}}：SBS1/Anderson2001 + 5 卡升 primary。
