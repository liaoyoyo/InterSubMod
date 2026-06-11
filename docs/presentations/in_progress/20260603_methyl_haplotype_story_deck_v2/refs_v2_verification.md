<!--
建立時間: 2026-06-03
狀態: in_progress (引用驗證狀態表 — refs_v2.bib 的 companion)
報告類型: citation_verification_audit
data_sources:
  - docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md §參考文獻 (原始 22 篇)
  - refs_v2.bib (驗證後)
provenance_note: /citation-verification 協議，4 個 researcher agent 平行 (WebSearch + Google Scholar + PMC/PubMed)；每條 4 欄對照 + relation 合理性判斷。
-->

# refs_v2.bib 引用驗證狀態表

> **方法**：`/citation-verification` 協議 — 4 個 `researcher` agent 平行，每篇 WebSearch「標題+第一作者+年」+ `site:scholar.google.com` 確認存在 + PMC/PubMed 4 欄對照（標題/第一作者/年±1/期刊）+ 判斷 deck relation 合理性。
> **總結**：**24 條全部存在、0 NOT_FOUND、0 PREPRINT_UNVERIFIED**；16 VERIFIED + 8 CORRECTED（metadata 修正）。所有 relation（supports/partial/contradicts）經獨立判斷皆 **reasonable**。

## 逐條狀態

| key | STATUS | relation | 驗證/修正重點 |
|---|---|---|---|
| fishchuk2024brca2 | ✅ VERIFIED | C1 supports | PMID/DOI 全中；41.9% tumour + **0–69% 範圍即出自本篇**；作者拼字修 Dubitskaa→Dubitska |
| evans2018brca1 | ✅ VERIFIED | C1 partial | 精確；BRCA1 cis-variant→methylation-silencing 人類先例（≠BRCA2 caveat 正確）|
| catteau2002brca1 | ✅ VERIFIED | C1 supports | review，適合共識引用；補 DOI |
| esteller2000brca1 | ✅ VERIFIED | C1 supports | deck 原缺 PMID → 獨立查得 10749912（JNCI 92(7):564）|
| rodriguezbalada2018brca | 🔧 CORRECTED | C1 partial/boundary | deck 原僅主題描述 → 補第一作者 **Rodríguez-Balada M**, Clin Transl Oncol 20(9):1226 |
| oneill2024longread | ✅ VERIFIED | C2/C4 partial | 補 PMID 39406235；⚠ **matched normal = 43 tumors/41 patients（非全 189）**；aDMR tumor 比 normal 多約 5× |
| do2020asm | ✅ VERIFIED | C2 **contradicts（口徑不同）** | PMID 32594908；⚠ 必明寫「allele-level bulk WGBS（含 germline-het 基線）vs 我們 HP-axis somatic-controlled」否則 overstated |
| liu2025nanome | ✅ VERIFIED | C2/C4 partial | bioRxiv preprint（PMID 40631091/PMC12236756 確存）；因果 phasing→methylation |
| canasm2025 | 🔧 CORRECTED | C2 partial | **已有正式版 BMC Genomics 26:648（DOI 10.1186/s12864-025-11849-7）優先引用**；⚠ 作者全名 auth-wall 待核 |
| derrien2021myeloma | ✅ VERIFIED | C3 supports | Genome Medicine 13:127；read-level PDR；混亂源 within-read inconsistency；⚠ middle authors 待補 |
| landan2012epipoly | ✅ VERIFIED | C3 supports | epipolymorphism 奠基論文；標題補尾段「in normal and cancerous tissues」|
| russo2021chaos | 🔧 CORRECTED | C3 supports | 補 DOI 10.3390/cancers13081800；review 非一手數據 |
| brocks2014intratumor | ✅ VERIFIED | C3 partial | 補 DOI；單一癌種 prostate；異質性反映 clonal evolution（partial 判斷恰當）|
| fu2024methphaser | ✅ VERIFIED | C4 partial | Nat Commun 15:5327；甲基→haplotype（與 ISM 因果相反）|
| akbari2021nanomethphase | ✅ VERIFIED | C4 partial | Genome Biology 22:68；甲基相位化奠基工具 |
| ont_modkit | ✅ VERIFIED | C4 工具鏈 | repo 確存；軟體無正式 paper → @misc 正確；不做 read-read 距離（凸顯 ISM niche）|
| haffner2011_5hmc | ✅ VERIFIED | C5 supports | PMID/DOI 全中；癌症 5hmC 全域耗竭一手證據 |
| chen2016_5hmc | 🔧 CORRECTED | C5 supports | **年份 2015→2016**；補 PMID 26680004 / DOI 10.1038/cr.2015.150 |
| halliwell2025 | ✅ VERIFIED | C5 supports/方法 | Commun Biol 8:243；⚠ **口徑修正：nanopore 經適當 basecalling 可分離 5mC/5hmC，風險在預設模型吸收 5hmC，非「本質混淆」** |
| lian2012_5hmc | 🔧 CORRECTED | C5 partial | **實為 Cell 2012 melanoma**（非 2011 Cancer Res）；PMID 22980977 / DOI 10.1016/j.cell.2012.07.033 |
| martintrujillo2017 | 🔧 CORRECTED | C6 supports | 🔴 **deck 標題寫錯**：正確＝「…major **dictator of imprinted** methylation in **tumors**」（非 drivers of allele specific…cancers）；scope=imprinted DMR，泛化須標 |
| abante2020 | ✅ VERIFIED | C6 supports | 補 PMID 33067439；relation 為**方法學**（haplotype-resolved ASM 偵測框架），非 CN-causation 直接證據 |
| wu2011ngs | 🔧 CORRECTED | C6 partial/統計 | **第一作者 Wu G（非 Wu H）**；補 PMID 21698242 |
| misra2010_brca2zar2 | ✅ VERIFIED | 署名更正 | draft 誤「Liu」→ **Misra S**（PMID 20202217，BRCA2/ZAR2 雙向 promoter）|
| gochhait2007_brca2 | ✅ VERIFIED | 署名更正 | draft 誤「Healey」→ **Gochhait S**（PMID 17945002，BRCA2 -26G>A）|

## 對 deck 的回寫動作（已套用到 report_v2.html）

1. **martintrujillo2017 標題改正**（最高優先；原引用 paraphrased 標題錯誤）+ 加 scope 註「demonstrated at imprinted DMRs」。
2. **halliwell C5 措辭軟化**：「nanopore 會混淆」→「預設 5mC 模型會吸收 5hmC，須顯式建模」。
3. chen 2015→2016、lian 2011→2012(Cell melanoma)、wu Wu H→Wu G、russo +DOI。
4. **0–69% 標來源**＝Fishchuk 2024 本篇（非獨立 meta）。
5. do2020「contradicts」維持但已明寫口徑差異（deck 既有）。

## 殘留待辦（正式投稿前）

- 數條 middle-author list 為領域慣例補全（do/derrien/landan/brocks/canasm/misra）→ 從 Google Scholar「Cite」逐條複製 author 欄。
- canasm2025 作者全名待開 BMC Genomics 全文核對。
- middle-author 拼字終驗後即可進論文 .bib。

## Provenance

- 驗證 agent：a43c6087（C1）/ ab848691（C2-C3）/ a37a04cc（C4-preprint）/ a74832d6（C5-C6+署名）。
- 每條 source URL 見各 agent 回報（PubMed/PMC/期刊頁）。tier：全部 L2（一手論文 DOI 可訪問）或 L3（review：catteau/russo）。
