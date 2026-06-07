---
id: ism-kb-11-external-literature-references-provenance
name: "參考文獻 + Provenance + 識別碼核對 + 未覆蓋工具"
description: "本地景 ~38 篇 WebFetch-verified 外部文獻完整清單(分類+真實 PMID/DOI+credibility)；跨 agent 識別碼不一致核對(POG 三 DOI / epihet 作者-PMID / HiCancer venue)；完整性 critic 點出的未覆蓋工具(WhatsHap/LongPhase/methylKit/DSS/Gaiti2019/amrfinder/MLML/COLO829 benchmark)；workflow provenance。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: reference
verified_scope: "all identifiers verified via WebFetch in wf_37b2cc97-663 (2026-06-05) unless marked UNVERIFIED; reconciliation notes per completeness critic"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-goal-alignment-methods
  - ism-kb-11-external-literature-finding-methyl-diff-loci
  - ism-kb-11-external-literature-ism-method-comparison
  - ism-kb-11-external-literature-subclone-landscape
  - ism-kb-11-external-literature-methyl-tp-fp-germline-somatic
  - ism-kb-11-external-literature-methyl-assisted-phasing
  - ism-kb-11-external-literature-asm-cis-cancer-impact
  - ism-kb-11-external-literature-cross-question-correlation
tags: [references, provenance, citation-verification, identifier-reconciliation, missing-tools, bibliography]
canonical_paths: [11_external_literature/09_references_and_provenance.md]
alias_paths: []
---

# 參考文獻 + Provenance + 識別碼核對 + 未覆蓋工具

- **一句說明**：本地景所有外部文獻經 workflow wf_37b2cc97-663 的 verify agent **實際 WebFetch** 取得真實 PMID/DOI；無法 fetch 確認者標 **UNVERIFIED**（不編造識別碼）。正式投稿引用前，仍須用 `/citation-verification` 對每篇再核一次（尤其下方 §識別碼核對 的 3 個衝突）。

## §A 完整參考文獻（依主題；標 credibility）

> cred = high(peer-reviewed 已驗) / pre(preprint) / uv(識別碼未驗) / sec(secondary/review/vendor)

### A1 長讀方法與 phasing 工具
1. Akbari V, et al. **NanoMethPhase**. *Genome Biology* 2021. PMID 33618748 / DOI 10.1186/s13059-021-02283-5. [high；Q1/Q2/Q3/Q5/Q6/Q7 方法骨幹]
2. Fu Y, et al. **MethPhaser**: methylation-based long-read haplotype phasing. *Nat Commun* 2024. PMID 38909018 / DOI 10.1038/s41467-024-49588-0. [high；Q3/Q4/Q6]
3. Liu Y, et al. **NANOME**: haplotype-aware allele-specific consensus methylation. *bioRxiv* 2025. PMID 40631091 / DOI 10.1101/2025.06.29.662079 (PMC12236756). [pre；Q1/Q2/Q3]
4. **HapBridge**: methylation-based iterative tagging of unphased reads. *bioRxiv* 2025. DOI 10.1101/2025.11.07.687303. [pre；Q6 — 我們 methyl-rescue 的 germline 機制 analog]
5. Snajder R, et al. **pycoMeth**: haplotype-aware read-level Bayesian DMR. *Genome Biology* 2023. DOI 10.1186/s13059-023-02917-w（PMID 37076875，未直接 fetch；DOI 為權威）. [high/uv-PMID；Q2/Q3 最近方法對手]
6. Raineri E, et al. **cvlr**: read clustering via multivariate Bernoulli mixture. *Bioinformatics Advances* 2023. PMID 36726731 / DOI 10.1093/bioadv/vbac101. [high；Q3 最近 read-level clustering analog]
7. Ho (HKU-BAL), et al. **LongPhase-S**: somatic haplotyping + tag→caller. *bioRxiv* 2025. DOI 10.1101/2025.11.20.689492. [pre；Q6 — HKU 協作工具，ClairS +4.5% SNV]
8. **Wakhan**: LOH/allelic-imbalance-driven whole-chromosome tumor phasing. *medRxiv* 2025. DOI 10.64898/2025.12.11.25342098. [pre；Q6 — 我們 LOH-axis 機制(非甲基)]

### A2 POG 長讀癌症 cohort（⚠識別碼見 §B1）
9. O'Neill K / Pleasance E / (Lin?), et al. **POG long-read advanced cancer cohort**. *Cell Genomics* 2024. **PMID 39406235 / DOI 10.1016/j.xgen.2024.100674 / PMC11605692**（採此為 canonical）. [high；Q1/Q2/Q6/Q7 — 最大真癌 ASM cohort，復發 RET/CDKN2A，79% aDMR in CNV/LOH]

### A3 癌症 phasing（非甲基對照）
10. Sakamoto Y, et al. Long-read NSCLC phasing. *Nat Commun* 2022. PMID 35710642 / DOI 10.1038/s41467-022-31133-6. [high；Q6/Q7 — patient-specific hap-DMR **確認我們 private**；甲基不用來 phase]
11. Mao Q, et al. **HiCancer**. *Sci Rep* 2021. PMID 33758310 / DOI 10.1038/s41598-021-86104-6（⚠venue 見 §B3）. [high；Q6 — LOH 處理原則(Hi-C)]

### A4 失序 / epipolymorphism 工具與奠基
12. Li S, et al. **methclone**. *Genome Biology* 2014. PMID 25260792 / DOI 10.1186/s13059-014-0472-5. [high；Q2/Q3/Q4 disorder 極]
13. Chen X / (Pan?), Ashoor H, et al. **epihet**. *Sci Rep* 2021. DOI 10.1038/s41598-020-79627-x / PMC7801679（⚠作者+PMID 見 §B2）. [high；Q2/Q3/Q4]
14. Lee S, et al. **Metheor** (incl. LPMD within-read distance). *PLOS Comp Biol* 2023. DOI 10.1371/journal.pcbi.1010946（PMID UNVERIFIED）. [high/uv-PMID；Q3]
15. Landau DA, et al. Locally disordered methylation in CLL (PDR). *Cancer Cell* 2014. PMID 25490447 / DOI 10.1016/j.ccell.2014.10.012. [high；Q2/Q3/Q4]
16. Landan G, et al. Epipolymorphism / stochastic DMR formation. *Nat Genet* 2012. PMID 23064413 / DOI 10.1038/ng.2442. [high；Q2]
17. Scherer M, et al. WSH score battery (coverage sensitivity). *Nucleic Acids Res* 2020. PMID 32103242 / DOI 10.1093/nar/gkaa120. [high；Q2 — **確認 O11 coverage-confound 機制**]
18. Jenkinson G, et al. **informME** (NME/dNME entropy). *BMC Bioinformatics* 2018. PMID 29514626 / DOI 10.1186/s12859-018-2086-5. [high；Q2]

### A5 Subclone / phyloepigenetics / clock
19. Hong J, Siegmund KD, et al. Infer tumor ancestry from passenger methylation. *PLOS One* 2010. PMID 20711251 / DOI 10.1371/journal.pone.0012002. [high；Q4 — phyloepigenetics 奠基]
20. Yuan K, et al. **BitPhylogeny**. *Genome Biology* 2015. PMID 25786108 / DOI 10.1186/s13059-015-0592-6. [high；Q4 — phylogeny 金標準]
21. Gabbutt C, Duran-Ferrer M, et al. **EVOFLUx** (fluctuating CpG clock). *Nature* 2025. PMID 40931062 / DOI 10.1038/s41586-025-09374-4. [high；Q4 — bulk subclonal 推論 SOTA(n=1,976)]
22. Larose Cadieux C, Castignani C, et al. (Van Loo/TRACERx) **CAMDAC**. *Nat Genet* 2025. PMC12425823 / DOI 10.1038/s41588-025-02307-x. [high；Q4 — CN-aware deconv 金標準]
23. Lee, et al. Long-read single-cell-derived melanoma subclones. *bioRxiv* 2025. DOI 10.1101/2025.08.28.672865 / PMC12424993. [pre；Q4 — 最近 ISM setup 的長讀 subclone 重建]
24. Barrett J, Mazzocco G, et al. Bayesian epiallele tumor evolution. *BMC Bioinformatics* 2017. PMID 28743252 / DOI 10.1186/s12859-017-1753-2. [high；Q3 — subclonal deconv 目標一致]
25. Zhang Z / Sun, et al. DNAm haplotype map of 11 cancers (81,567 MHB). *Cell Reports* 2025. PMID 40849904 / DOI 10.1016/j.celrep.2025.116197. [high；Q2]

### A6 ML 陷阱 + somatic variant calling
26. Kapoor S, Narayanan A. Leakage in ML (8-type taxonomy). *Patterns* 2023. PMID 37720327 / DOI 10.1016/j.patter.2023.100804. [high；Q5 — **L3.2/L1.3 命名我們 LOSO 循環**]
27. Soneson C, et al. Batch-class confounding inflates CV. *PLOS One* 2014. DOI 10.1371/journal.pone.0100335（PMID 24967636，本輪 agent 標 UNVERIFIED；前 wf 已驗）. [high；Q5]
28. Zheng Z, et al. **ClairS-TO**. *Nat Commun* 2025. DOI 10.1038/s41467-025-64547-z（PMID UNVERIFIED）. [high/uv-PMID；Q5 — **無甲基**達 TO FP 控制]
29. Zheng Z, et al. **ClairS** (LongPhase-S haplotype-consistency filter). *bioRxiv* 2023. DOI 10.1101/2023.08.17.553778. [pre；Q5]
30. Wood DE, et al. **VarNet**. *Nat Commun* 2022. PMID 35869060 / DOI 10.1038/s41467-022-31765-8. [high；Q5 — alignment 特徵非甲基]
31. Kim S, et al. **AIVariant**. *Exp Mol Med* 2023. PMID 37524869 / DOI 10.1038/s12276-023-01049-2. [high；Q5 — epigenetic 特徵在低 purity DL 有助(唯一部分反例)]
32. Krishnamachari K, et al. cfDNA ML FP filter. *Sci Rep* 2023. DOI 10.1038/s41598-023-37409-1 / PMC10300101（PMID UNVERIFIED）. [high/uv-PMID；Q5]

### A7 ASM / cis 生物與癌症影響
33. Do C, et al. ASM increased in cancers. *Genome Biology* 2020. PMID 32594908 / DOI 10.1186/s13059-020-02059-3. [high；Q7 — **+5× hypo-dominant，口徑差非矛盾**]
34. Martin-Trujillo A, et al. **CN rather than epigenetic drives allelic methylation**. *Nat Commun* 2017. PMID 28883545 / DOI 10.1038/s41467-017-00639-9. [high；Q4/Q7 — **CN-confound 最強佐證**]
35. Schalkwyk LC, et al. Allelic skewing of DNA methylation (normal). *Am J Hum Genet* 2010. PMID 20159110 / DOI 10.1016/j.ajhg.2010.01.014. [high；Q7 — normal cis-ASM baseline]
36. Do C, Shearer, Greally, Tycko. mQTL vs hap-ASM (review). *Genome Biology* 2017. PMID 28629478 / DOI 10.1186/s13059-017-1250-y. [high；Q7 — caliber framing 骨幹]
37. Hesson LB, et al. BRCA1 promoter methylation rare second hit in PDAC. *Mol Diagn Ther* 2022. DOI 10.1007/s40291-022-00614-1 / PMC9626413（PMID UNVERIFIED）. [high/uv-PMID；Q7]
38. Herman JG, Baylin SB. Gene silencing by promoter hypermethylation (VHL). *J Clin Invest* 2001 review. DOI 10.1172/JCI9462 / PMC289180（UNVERIFIED）. [sec/uv；Q7 — **canonical TSG hyper 方向**]
39. Fine-mapping regulatory variants by native CpG methylation (nanopore). *bioRxiv* 2024. PMID 39386487 / DOI 10.1101/2024.09.27.614715. [pre；Q5 — ASM genotype-driven，含 HCC1395/COLO829]
40. "Switchable / genetics-influenced ASM" paradigm. *Cell Discovery* 2017. PMID 29387450 / DOI 10.1038/celldisc.2017.38. [high；Q5]

### A8 平台 / 時間軸
41. Teer JK, et al. Tumor-only vs matched normal precision. *Human Genomics* 2017. PMID 28870239 / DOI 10.1186/s40246-017-0118-2. [high；Q1]
42. de Abreu, et al. WGBS/EM-seq/ONT/array benchmark. *Epigenetics & Chromatin* 2025. PMID 40855329 / DOI 10.1186/s13072-025-00616-3. [high；Q1]
43. Liu Y, et al. **NanoEM** (enzymatic long-read methylation). *Nucleic Acids Res* 2021. DOI 10.1093/nar/gkab397（PMID UNVERIFIED）. [high/uv-PMID；Q1]
44. Capper D, et al. Methylation-based CNS tumor classifier. *Nature* 2018. PMID 29539639 / DOI 10.1038/nature26000. [high；Q1 — array→seq 時間軸/TO 臨床規範]

## §B 識別碼核對（跨 agent 不一致 — 投稿前必核）

> §13 防假引用：以下 3 處不同 verify agent 回傳不一致；採「最一致 verified anchor」並標衝突，**不把衝突當事實並列**。

1. **POG 長讀癌症 cohort（Ref 9）**：3 個 article-number(`100674` / `100693` / `100538`) + 2 種作者署名(`O'Neill et al` / `Lin et al`) 在不同 angle 出現。**穩定 anchor = PMC11605692 + PMID 39406235**（Q2/Q6 一致，且與既有研究故事 ref #5 同）。採 **DOI 10.1016/j.xgen.2024.100674**。⚠ `100693`/`100538` 與「Lin」署名未交叉確認 → 正式引用前用 PMID 39406235 在 PubMed 核對確切 article number + 第一作者。
2. **epihet（Ref 13）**：作者 `Chen et al`(Q2/Q4) vs `Pan et al`(Q3)；PMID `33432081`(Q2) vs `33414514`(Q3)。**穩定 anchor = DOI 10.1038/s41598-020-79627-x / PMC7801679**（三 angle 一致）。既有研究故事 ref #24 用 `Chen X, Ashoor H, et al`。⚠ 投稿前核第一作者 + PMID。
3. **HiCancer（Ref 11）**：Q6 verify agent 自 flag「candidate 宣稱的 venue 有誤，實為 *Scientific Reports*」。採 **PMID 33758310 / DOI 10.1038/s41598-021-86104-6 / Sci Rep**。⚠ 勿沿用 candidate 的錯 venue。
4. **credibility 標註修正**：NanoEM/Liu2021(Ref 43)、ClairS-TO(28)、Krishnamachari(32)、pycoMeth(5)、Metheor(14) 等 PMID 未直接 fetch → 標 **uv-PMID**（DOI 已驗）；勿在投稿標「PMID 已驗」。

## §C 未覆蓋工具（completeness critic 點出 — 本輪未 WebFetch，仍須 position-against）

> 以下是 critic 認為**應對照但本輪未驗證**的工具/論文。**不在此編造其發現**；列為投稿前補查清單（用 `/citation-verification`）。

| 角度 | 應補的工具/論文 | 為何重要 |
|---|---|---|
| Q1/Q6 | **WhatsHap / WhatsHap polyphase**、**HapCUT2**、**LongPhase**(base) | 我們反覆說「via phasing」卻沒引 de-facto 標準 SNV phaser；polyphase 是「cluster reads 成 haplotype 群」最近的*演算法* analog；reviewer 必問為何 methyl-rescue 不對 polyphase benchmark |
| Q2 | **DSS / methylKit / metilene / bsseq** | canonical DMR caller；「我們找甲基差異位點」須對標，否則「ISM 比標準 DMR calling 多什麼」無防守。methylKit 是 reviewer「為何不用 methylKit windows」的第一問 |
| Q2/Q3 | **modbamtools / methylartist / NanoMethViz / pb-CpG-tools** | 標準 ONT/PacBio modBAM per-read 視覺化；我們 IGV-可見「稀有但真實」+read×CpG 矩陣應引；pb-CpG-tools 是 PacBio 側跨平台 generality |
| Q3/Q5 | **DAMEfinder / epialleleR** | allele-aware 差異甲基(DAMEfinder)；epialleleR 做 read-level epiallele/ASM 評分 = 最近的*已採用*工具 |
| Q4 | **Gaiti et al. 2019 (CLL, Nature) scWGBS clonal phylogeny**、epiCHASS/Hannon | **Q4 最大缺**：Gaiti 2019 是 canonical 單細胞甲基→clonal-phylogeny，正是我們「subclone reconstruction」動詞所需的 ground-truth/recurrence 先例 |
| Q5/Q7 | **methpipe allelicmeth / amrfinder (Fang/Smith)**、**MethCNA** | amrfinder = canonical genotype-independent ASM caller = 我們 normal-anchored cis-test 的*概念祖先*；MethCNA = CN-from-methylation orthogonal check |
| Q7 | **MLML (Qu et al.)**、**epiTOC/epiTOC2**、Knudson 表觀第二打擊、Woodhouse imprinting atlas | MLML 是 5mC/5hmC 同時估計的 canonical 方法——我們 5mC/5hmC 分軌 novelty *缺方法引用*；epiTOC2 是甲基-clock(對 EVOFLUx)；表觀第二打擊直接關 Q7 BRCA2 framing(且撐 canonical hyper 方向) |
| Q6 | **COLO829 long-read benchmark / truth-set**(Valle-Inclán 等) | COLO829 反覆被當「⭐4 ground-truth」卻無 benchmark 論文引用，「⭐4 需要什麼」無 anchor |

## §D Provenance

- **外部文獻檢索**：workflow **wf_37b2cc97-663**（2026-06-05；16 agents = 7 scout + 7 verify(researcher agentType，WebSearch/WebFetch) + cross-synthesis + completeness-critic；210 tool uses；~8 min）。每篇由 verify agent 實際 WebFetch/WebSearch 取得識別碼 + 關鍵數字；無法 fetch 確認者標 UNVERIFIED。
- **內部數字**：沿用既有 validated 來源（研究故事 §7 / master_draft / evidence_ledger / tsg snapshot @06-03），本輪**未跑新分析**。內部 L1-L2 / 外部 L3。
- **撰寫隔離（§13.0）**：本 KB 9 份 .md 的 Write 與外部 WebFetch 查證**不同 batch**——先 workflow 查證落 /tmp → Read 回 → 才撰寫。報告每個外部數字可在對應 angle 檔 + 此清單 grep 到。
- **後續**：投稿前須 (1) `/citation-verification` 核 §B 三衝突 + 所有 uv-PMID；(2) 補 §C 未覆蓋工具。

**相關**：[00 索引](00_index.md) ｜ [08 跨問題鏈](08_cross_question_correlation.md) ｜ 7 角度 Q1-Q7（見索引）
