---
id: ism-kb-11-external-literature-goal-alignment-methods
name: "Q1 目標對齊：長/短讀方法、時間軸、paired vs tumor-only"
description: "ISM 的 ONT 長讀 read-level haplotype-linked 甲基方向，對齊文獻方法地景與 array→WGBS/EM-seq→long-read 時間軸；paired tumor-normal vs tumor-only 設計規範。內部 L1-L2 × 外部 L3。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal claims vs research story §1 + five_research_goals"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, long-read, short-read, timeline, paired, tumor-only, nanopore, methodology]
canonical_paths: [11_external_literature/01_goal_alignment_methods_timeline.md]
alias_paths: []
---

# Q1 目標對齊：長/短讀方法、時間軸、paired vs tumor-only

- **一句裁決**：方向 **CONFIRM**——(a) 長讀獨家使能 read-level haplotype-linked 甲基（array/WGBS 做不到）、(b) array→WGBS/EM-seq→long-read 時間軸、(c) paired vs TO 是公認 somatic 設計軸，**三者全被高品質 peer-reviewed 背書**。ISM 的**程式定位**穩坐活躍領域；ISM 的**特定方法**（read-read 距離 + cis-test）相對 niche 是 DIFFERENT-VIEW（未被反駁但也未對 incumbents 證優）。

## 問題與內部立場

ISM 用 ONT 長讀建 read×CpG 甲基矩陣，取得 read-level、可綁 haplotype 的甲基訊號（這是短讀/bisulfite/array 結構上做不到的）。承載[五大目標](../01_project_overview/02_five_research_goals.md)。要確認：我們的長讀 read-level 路線在領域方法地景與歷史時間軸中的位置，以及 paired(tumor+normal) vs tumor-only 是否為公認設計軸。

**內部錨點（L1-L2）**：ONT 5mCG+5hmCG modBAM；pipelines = paired_full / paired_pileup / TO；主樣本 HCC1395（有 SEQC2 ground-truth）。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑差異 | Cred |
|---|---|---|---|---|---|---|---|
| Akbari, **NanoMethPhase**, Genome Biol | 2021 | PMID 33618748 / 10.1186/s13059-021-02283-5 | 長讀 ONT；germline tool | ~10× 全基因組偵測 ASM（SNP+phasing+甲基，megabase 尺度）| **CONFIRM**（方法可行）| germline/imprinting，非 somatic；region-level 非 read-distance | high |
| O'Neill/POG, Cell Genomics | 2024 | PMID 39406235 / 10.1016/j.xgen.2024.100674 / PMC11605692 ⚠識別碼見 09 | 長讀 ONT PromethION；**paired cohort 189 腫瘤+41 normal** | long-range phasing → aDMR（含復發 RET/CDKN2A）；median 17.5× normal/26× tumor，N50 31.3kb | **CONFIRM**（最近 cohort 先例）| cohort vs 我們單樣本；region-level aDMR vs 我們 read-read+cis；他們**復發**我們 private 0/38 | high |
| Teer, Human Genomics | 2017 | PMID 28870239 / 10.1186/s40246-017-0118-2 | 短讀（突變類比）；TO vs matched | TO recall 似(43–60%) 但 **precision 20–21% vs matched 30–82%** → 不宜新發現 | **CONTEXT**（設計規範動機）| 突變非甲基；量化「為何 paired 是標準」 | high |
| Liu, **NANOME**, bioRxiv | 2025 | PMID 40631091 / 10.1101/2025.06.29.662079 | 長讀 ONT；germline imprinting | XGBoost 共識(Megalodon/Nanopolish/DeepSignal)：+11% precision、+2.4% F1、~200k 多 CpG | **DIFFERENT-VIEW**（鄰接 niche）| consensus *caller*，非 read-read 距離/somatic cis；germline | preprint |
| de Abreu, Epigenetics & Chromatin | 2025 | PMID 40855329 / 10.1186/s13072-025-00616-3 | 混合(WGBS/EM-seq/ONT/array)；3 樣本 | EM-seq 與 WGBS 最一致；**ONT 一致度較低但抓獨特 loci + 長程** | **CONTEXT**（平台定位）| benchmark 非 haplotype/somatic；佐證 ONT 獨特-locus+長程優勢 | high |
| Liu (**NanoEM**/酵素轉化), NAR | 2021 | DOI 10.1093/nar/gkab397（PMID UNVERIFIED）| 長讀 ONT+EM 轉化；乳癌株+臨床 matched normal | 保留 3–7kb 片段，偵測短讀做不到的重複/缺失區 ASM | **CONTEXT**（化學橋）| 轉化非 native modBAM（我們 native 5mC/5hmC）；乳癌域重疊 | peer-rev（識別碼未驗）|
| Capper, Nature | 2018 | PMID 29539639 / 10.1038/nature26000 | Array 450K；TO 分類器(CNS) | 甲基分類器；前瞻案例 **改診斷達 12%**；入 WHO 2021 | **CONTEXT**（array→seq 時間軸；TO 臨床規範）| array bulk-only，與 read-level 對立極 | high |
| ONT tumour-normal workflow（vendor）| 2023+ | UNVERIFIED（無 DOI）| 長讀 ONT；paired | 場域規範覆蓋 **normal 30×/tumor 60×**；paired=somatic/表觀標準 | **CONTEXT**（設計規範）| vendor 指南非 peer-review | secondary |

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| 長讀使能 read-level HP-linked 甲基，array/WGBS 做不到（L1-L2）| CONFIRM：NanoMethPhase / POG / de Abreu / Capper-對比 | 🟢 none | **Confirmed** — field-validated，well-established 方向 |
| paired vs TO 是公認 somatic 設計軸（L2）| CONFIRM：ONT workflow(30×/60×) / POG(41 normal) / Teer(precision 30-82% vs 20-21%) | 🟢 none | **Confirmed** — paired 是標準；我們 TO 的 germline-confound 是公認代價 |
| 時間軸 array→WGBS/EM-seq→long-read（L2）| CONFIRM：Capper18→Liu21/de Abreu→NanoMethPhase21/POG24/NANOME25 | 🟢 none | **Confirmed** — 可引用 |
| read-read 距離+clustering = unfilled niche（L1-L2）| DIFFERENT-VIEW：上述工具做 haplotype-aware ASM 但在 region/DMR/imprinting 經 phasing，非 read-read 距離矩陣+cis-test | 🟡 caliber-only | **部分 confirmed**：haplotype-aware ASM 擁擠；read-read+cis 機制似真未佔，但文獻 silent（非反駁）on 是否優於既有 phasing |
| somatic ASM private 0/38，單樣本 ⭐3（L2）| 對比 POG aDMR **跨 189 cohort 復發**(RET/CDKN2A)| 🔴 direction（recurrent vs private，可調和）| **Caliber-limited** — POG 證 recurrence 在 cohort 可達；我們單樣本誠實不宣稱 |

## 缺口（我們缺什麼）

1. **Cohort scale + recurrence**：POG（189 腫瘤）證**復發** aDMR（RET/CDKN2A）= 臨床 payoff；我們最強是單樣本 HCC1395 ⭐3，somatic private 0/38；cross-sample 6 樣本但單 pipeline，⭐3 封頂。
2. **獨立 caller / 多工具共識**：NANOME 顯示領域期待多工具共識(XGBoost over 3 caller) 才宣稱穩健 ASM；我們單 caller，未對 Megalodon/Nanopolish/DeepSignal benchmark。
3. **同資料 head-to-head**：宣稱「unfilled niche」但未在 HCC1395 跑 NanoMethPhase/NANOME 證 read-read+cis 加值。
4. **平台一致性驗證**：de Abreu 顯示 ONT WGBS 一致度較低；我們未對 orthogonal(WGBS/EM-seq/array) 量化 5mC truth（甲基 truth 目前內部）。
5. **5hmC 未充分利用**：領域視 native 5hmC 為 frontier；我們 5hmC marginal 0.03-0.07 大致擱置，未證 5mC/5hmC 分軌產出領域所缺的發現。

## 裁決

方向 **CONFIRM 非 refute**：三前提（長讀獨家 read-level、array→WGBS→long-read 時間軸、paired vs TO somatic 設計軸）全被 NanoMethPhase/POG/Teer/de Abreu/Capper/ONT-workflow 獨立背書，ISM **程式定位**穩坐活躍領域。ISM **特定方法**（read-read 距離矩陣 + 階層 subclonal clustering + normal-anchored cis-test + 5mC/5hmC 分軌）相對 niche 是 **DIFFERENT-VIEW**：haplotype-aware ASM 擁擠（NanoMethPhase/NANOME/POG 都做），但都在 region/DMR/imprinting 經 phasing，我們的 read-read+cis 機制似真未佔——文獻 *silent* 非反駁，niche 成立但**對 incumbents 未證**（無 head-to-head）。誠實限制：每個 ASM 發現皆單樣本 HCC1395 ⭐3；cross-sample 單 pipeline ⭐3 封頂；無 cohort-recurrence 軸（POG 的 payoff）、無獨立 caller 共識（NANOME 的門檻）、無 orthogonal 平台甲基 truth（de Abreu 的框架）。POG 對比最尖：其 aDMR cohort 復發（RET/CDKN2A），我們 private(0/38)——可調和（不同軸）但提醒：無多 pipeline/多樣本/ground-truth 我們到不了 cohort-recurrence tier。

**相關**：[跨問題依賴鏈](08_cross_question_correlation.md) ｜ [參考文獻+provenance](09_references_and_provenance.md) ｜ Q6 phasing：[06](06_methyl_assisted_phasing.md) ｜ Q7 ASM：[07](07_asm_cis_cancer_impact.md)
