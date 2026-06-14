<!--
建立時間: 2026-06-14
狀態: thesis-statement 定稿 v1 (用戶逐句確認後鎖定; Intro 開頭 + 貢獻定位 + 論據 + 待補 citation)
報告類型: paper_intro_thesis_statement
受眾: 廖子游 · PI · 碩論 Ch1 Intro + Ch4/5 contribution
framework: Verdict-Pyramid + Assertion-Evidence + claim-tier (Tier 4 paper-scope)
data_sources:
  - InterSubMod/docs/paper_focus/03_references/20260613_external_vs_ism_five_bucket_comparison_01.md (五桶+趨同+獨特改進)
  - InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/01_ism_method_spec_from_source.md (ISM 6 cores)
  - InterSubMod/docs/paper_focus/05_html_staging/flagship_chr2_18086020_subclone_20260612.standalone.html (chr2 案例 grep-verified)
  - knowledge/11_external_literature/{05,06,07,08,10} (內部結果 grep-verified)
provenance_note: 主軸句逐句對應 verified 證據(內部 tier + 外部源)。BRCA2/chr2 數字 grep-verified。甲基→突變因果機制標 {{待補 citation}}(我方庫未驗, 投稿前 add+citation-verify)。本檔為 framing 文件, 不產新數字; §13.0 撰寫與分析分離。
-->
<!-- provenance-verified: BRCA2 Δβ=-0.122 / chr2 F=10.6 V=0.67 / 18096269 V=0.97 引自 flagship HTML + ISM spec; 外部 citation 引自五桶報告; 5mC-deamination 標待補。 -->

# Intro 主軸定稿 — Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing

> **這份是什麼**：用戶 2026-06-14 逐句確認後鎖定的論文主軸句（Intro 開頭）+ 貢獻定位 + 每句論據 + 待補 citation。**精修原則**：保留用戶願景，每句守住 reviewer-trap / framing-tension 紅線（依五桶報告）。

---

## L0 — 主軸句（定稿 v1，可當 Intro 開頭，逐句可防守）

> **InterSubMod (ISM)** 是**首個**在 **ONT 長讀 tumor/normal paired** 資料上、整合 **① read 層次甲基結構檢定 + ② matched-normal-anchored somatic cis-test + ③ somatic-subclone 目標** 的 read-level 軟體。它以 somatic haplotagging（longphase-S）產生的 **haplotype tag（HP1/HP2 + somatic 子單倍型 HP1-1/HP2-1）與 somatic 變異標記為框架**，在 read 層次量化每個 somatic 位點鄰域的**等位基因特異性甲基化（ASM）結構**，輸出三類資訊：
>
> 1. **與 somatic 變異高度關聯**、在 read×read 距離空間呈顯著結構分離的甲基位點（PERMANOVA）；
> 2. **tumor 內 haplotype 間 / tumor-normal 間有顯著差異**的甲基位點，並以 matched-normal 為錨點**分解**此差異屬真 focal cis 抑或 copy/subclone 背景（**BRCA2 promoter 為旗艦示範位點**）；
> 3. 甲基差異**與哪一軸（HP / ALLELE / copy）相關**的分解。
>
> 進一步，以「**多個 somatic 變異 + 鄰域甲基結構**」的組合，對單一 region 的 read 做 **subclone 結構 characterization**，並以 **chr2:18M 案例示範**甲基如何與突變**共同標示同一子克隆事件**、據此推導 *Normal→HP1→HP1-1* 的**時序假說**（自動化 subclone 分析流程與大規模驗證為 ISM 延伸）。
>
> **貢獻定位**：在**只有單樣本**的癌症資料上，提供 read 層次、normal-anchored、比 per-position 聚合**更細緻**的**甲基-變異聯合 characterization** 工具 —— 量化 somatic 變異與甲基位點的**關聯強度**（甲基→突變的因果機制為既有文獻背景〔§citation 待補〕，ISM **不宣稱**建立因果或先後順序），補足既有 germline-only 甲基-phasing 與 bulk subclone 重建工具未覆蓋的 **somatic / single-sample / read-level niche**。

---

## L1 — 用戶 4 點確認 → 精修映射（守線紀錄）

| # | 用戶本意（2026-06-14）| 精修後（定稿措辭）| 守的線 |
|---|---|---|---|
| 1 chr2/演化 | 軟體尚無法自動分析 subclone+流程，但 chr2 **案例**可用 sSNV 重建狀況、甲基確認 subclone 狀況；缺的是延伸驗證 | 「**案例示範**甲基如何與突變共同標示子克隆事件、推導 *時序假說*；自動化分析+大規模驗證=延伸」| 不輸出 lineage tree；時序=數據支撐**推論**非正交確認（gap audit）|
| 2 甲基造成突變 | 因果機制**有實際論文證實**，但我方無法確認先後順序；**可提供 sSNV-甲基關聯強度** | 「甲基→突變因果機制=**既有文獻背景**〔citation 待補〕；ISM 貢獻=量化**關聯強度**，不宣稱因果/順序」| 因果歸文獻背景，ISM 只做 association（不踩 causal overclaim）|
| 3 創新領先 scope | 接受明確敘述 | 「**首個整合** ①②③ 於 tumor/normal read-level」| 創新點=特定組合非空泛（Sgootr/qFDRP/NANOME 已佔各單件）|
| 4 BRCA2 | tumor 上 HP 差異甲基 + somatic 差異甲基**並存**於 BRCA2 promoter，與 normal 有別、HCC1395 乳癌，高度懷疑與突變有關，是最有用的例子 | 旗艦示範：HP-axis Δβ=−0.122 + tumor-normal 差並存於 BRCA2 promoter；normal-anchored **分解**顯示主要 track somatic subclone tag、focal cis marginal —— ISM **不只 flag 還量化**屬 subclone vs focal cis | 共定位/關聯=候選（defensible）；focal-cis 因果=marginal（誠實）；方向 hypo≠canonical TSG hyper（勿寫 silencing）|

---

## L2 — 三件貢獻（精確措辭 + 內部論據 tier + 外部論文 + novelty/corroboration）

| 貢獻 | 精確 claim | 內部論據（grep-verified, tier）| 屬性 | 要 cite 的論文 |
|---|---|---|:--:|---|
| **C1 read 層次甲基-變異聯合 characterization** | 在 HP tag + somatic 框架下量化每 somatic 位點鄰域 ASM 結構（read×read 距離 → 分群 → PERMANOVA）| ISM 6 cores（`StructureTest.cpp:141-205`）；BRCA2/chr2 案例 ⭐3 | **novelty=組合**（單件皆 prior art）| qFDRP(NHD kernel 先例)/cvlr/ASMS/DAMEfinder（read-level ASM 先例，須劃界）|
| **C2 read-level normal-anchored somatic cis/drift 分解** | residual=raw−normal mean + 三角 d_cis/d_drift/d_within，判別真 cis vs copy/subclone 背景 | `NormalBaseline.hpp`(residual) + `scripts/34`(三角)；chr17/TBC1D16 cis d_within=0.142 p=0.001 ⭐3 | **增量獨特**（🔴**非「獨家有 normal-anchor」**：Do2020 de-novo ASM 子分析亦用 within-patient matched-normal；ISM 增量 = **read-level within-sample 距離結構 + LOH-HP 軸 + 整合**）| **直接前例 somatic-variant-local-cis-methylation(Do2020 de-novo ASM + Zhang2019 SV)**；germline-cis baseline=cis-meqtl-germline-atlas(GoDMC)/Onuchic2018；CN 背書 Martin-Trujillo 2017 |
| **C3 somatic-subclone 結構 characterization（單樣本）** | 多 somatic 變異 + 鄰域甲基組合 → 單 region read→subclone 結構；case-demo 推時序假說 | chr2:18M case-demo ⭐3；LOH-phasing same-hap 93-99% 6/6 W=28 p=0.0078 | case-demo（自動化=future）| Sgootr/Gaiti/Epiclomal/MethylTree（距離式/單細胞甲基重建 prior art，**降 prior art + 劃 regime 界**）；MethPhaser/Wakhan/HiCancer（白地）|

**支撐性 corroboration（非 novelty，當外部效度背書）**：
- ASM 在癌症存在升高 → Do2020(5-9×, ⚠caliber 差)/cvlr/ASMS。
- 甲基=germline-haplotype 層級 → cvlr(8 gene P<10⁻³, H19 7/8)/ASMS(16/20 ICR)/O'Neill(4.46M aDMR)/longhap(Fig3C)。
- 弱-subclone 紅線 → EVOFLUx(1610/1976 neutral)/MethylTree(~100% Q 有監督, modality-bounded)。
- 甲基非 variant filter → 軸6 caller 全零甲基 + Kapoor/Soneson（filter-DEAD locked）。

---

## L3 — 兩個旗艦案例（grep-verified 事實 + 誠實分層）

### 旗艦 A — BRCA2 promoter（C1+C2 示範：偵測 + 分解）

- ✅ **實際數據**：HCC1395 BRCA2 promoter，HP-axis **Δβ=−0.122**（HP1 germline-tag vs HP1-1 somatic-subclone-tag）；tumor-normal 有差（normal-baseline residual）；共定位於 BRCA2（TSG）promoter。
- 🟡 **分解（C2 的價值）**：normal-anchored 分解顯示此甲基差**主要 track somatic subclone tag（HP1-1）**，**focal cis 效應 marginal**（d_within=−0.023, perm p=0.024, % 不 robust）。→ ISM **不只 flag「此處有與突變共定位的甲基差」，更量化它有多少屬 subclone 結構 vs 直接 focal cis**。這正是「我們最主要要提供的資訊」。
- ⚠ **誠實邊界**：方向 **hypo**（≠canonical TSG promoter hyper-silencing）→ 勿寫「promoter 高甲基沉默 BRCA2」；「與突變的關係」是**共定位/候選關聯**，focal cis 因果 marginal；單樣本 ⭐3。

### 旗艦 B — chr2:18M（C3 示範：甲基+突變共同標示 subclone）

- ✅ **實際數據**（`flagship_chr2_18086020` ISM 補跑，57 reads×203 CpG）：
  - 系譜（longphase-S）：Normal=G → **HP1 母本**（19 reads 全 G，germline）→ **HP1-1 子克隆**（36 reads = 28 A / 8 G，somatic A 富集）。
  - 18,086,020（het 標記子克隆）：甲基分群 **PERMANOVA F=10.6, p=0.01**，optimal_k=2，分群↔HP **Cramér's V=0.67**（甲基跟 haplotype/子克隆譜系）。
  - 18,096,269（somatic FP-candidate）：甲基分群↔**somatic 等位** **V=0.97**（HP=0），k=5，F=40.0 → 甲基跟 somatic 等位本身。
  - **兩種互補模式並存** = 「甲基+突變共同重建 subclone」的機制基礎。
- 🟡 **合理推論**：Normal→HP1→HP1-1 時序假說（從巢狀等位結構 + 甲基共變推導）。
- ⚠ **誠實邊界**：18,086,020 **不在 somatic VCF**（het 標記子克隆，甲基↔HP 含 germline-allelic 成分）；subclone 甲基的 **somatic-專一性**待 within-haplotype 對照（G-B 未跑）；**時序為數據支撐推論非正交確認**；自動化 subclone+流程分析與大規模驗證 = ISM 延伸（future work）；單樣本單 pipeline ⭐3。

---

## L4 — 誠實天花板（Intro/Discussion 須明標）

1. **case-demonstration**：chr2/BRCA2 是案例示範主軸概念，**非已驗證的自動化 subclone reconstruction**；大規模驗證 = future work（G-D/G-E）。
2. **無 lineage tree / 演化參數推論**：reconstruction = read→subclone 指派/分群結構；時序 = 案例假說。
3. **單樣本 single-pipeline ⭐3**：TP/FP 真值取自 tumor-only caller（非獨立 ground truth）；6 cell line = 跨模型 reproducibility 非病人 cohort 泛化。
4. **甲基為 characterization 非 utility**：不改善 variant calling（filter-DEAD AUC≈0.50）；不宣稱因果/順序。
5. **方向誠實**：BRCA2 hypo≠canonical TSG hyper；focal cis marginal（值在分解層）。

---

## L5 — Citation 清單（投稿前過 /citation-verification）

| 用途 | 論文 | 狀態 |
|---|---|---|
| read-level ASM/距離 先例 | qFDRP(Scherer 2020 NAR)/cvlr(Raineri 2023)/ASMS(Raineri 2024)/DAMEfinder(Orjuela 2020) | ✅ 庫內 verified |
| 距離式/單細胞甲基重建 prior art（劃界）| Sgootr(2023)/Gaiti(2019)/Epiclomal(2020)/MethylTree(2025) | ✅ 庫內（Sgootr 待升 primary）|
| cancer 甲基-phasing 白地 | MethPhaser(2024)/Wakhan(2025)/HiCancer(2021)/longhap(2026) | ✅ 庫內 verified |
| CN-confound 背書 | Martin-Trujillo(2017) | ✅ 庫內 |
| ASM 存在 | Do2020 | ✅ 庫內（caliber 差須標）|
| **somatic 變異→局部 cis 甲基（ISM cis-test 直接前例）** | **Do2020 de-novo ASM + Zhang2019 SV（`somatic-variant-local-cis-methylation`）** | ✅ **庫內（2026-06-14 新增）**；「rare/不大數值貢獻」**佐證我方 1/816 稀有**；tier C cited_secondary 待升 primary |
| somatic 變異→甲基（trans 背景）| Turcan2012 IDH1 + Figueroa2010（`somatic-mutation-methylation-modifier`）| ✅ 庫內（2026-06-14）；🔴 trans≠cis 須標 |
| germline-cis baseline（normal-anchor 扣的）| GoDMC Min2021（`cis-meqtl-germline-atlas`）+ Onuchic2018（`asm-sequence-haplotype-dependent`）| ✅ 庫內（2026-06-14）|
| CpG-SNP 假 ASM（limitation）| Shoemaker2010（`cpg-snp-pseudo-asm-confound`）| ✅ 庫內（2026-06-14）|
| reconstruction 骨幹 | LongPhase(Lin 2022, DOI 10.1093/bioinformatics/btac058)/longphase-S | ✅ 庫內 |
| PERMANOVA 方法源 | Anderson 2001 | 🟡 須加 |
| **甲基→突變因果機制（#2 用戶述方向，與上列反向）** | **Pfeifer 2006（PMID 16570852）+ COSMIC SBS1 / Alexandrov 2020 + Nat Commun 2024 BER**（`methylation-to-mutation-deamination`）| ✅ **2026-06-14 建卡 + 一手驗證**（Pfeifer PMID/COSMIC aetiology 親驗）；tier C cited_secondary（全文 PDF 未取、Alexandrov exact ID 待 citation-verify、量化 ~4×/~30% 待全文）。**ISM 不測方向**，作背景 + confound 依據（雙向並存 → 不可 naive 歸因）|

---

## Provenance
- 主軸句逐句對應五桶報告（`20260613_external_vs_ism_five_bucket_comparison_01.md`）+ ISM spec。
- BRCA2 Δβ=−0.122 / d_within=−0.023 引自 ISM spec L3（`knowledge/11_external_literature/07`）。
- chr2 案例（F=10.6/V=0.67/18096269 V=0.97/57×203/HP1-1 28A8G）引自 `flagship_chr2_18086020_subclone_20260612.standalone.html`（ISM 補跑 grep-verified, build 069cadb）。
- 5mC→突變 citation 標 {{待補}}（庫內 grep 無 verified 源；用戶述「有論文證實」屬大文獻 SBS1/deamination，須 add）。
- 用戶逐句確認紀錄：2026-06-14（#1 case-demo+延伸驗證 / #2 因果歸文獻+ISM 量關聯 / #3 接受 scope / #4 BRCA2 共定位+分解）。
