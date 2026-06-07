---
id: ism-kb-11-external-literature-ism-method-comparison
name: "Q3 ISM 方法是否合理有效：相似/差異工具，structure vs disorder"
description: "ISM(read-read 距離+subclonal clustering+LOH/CN 耦合+normal-anchored cis-test) 對比 disorder 工具(methclone/epihet/Metheor)、phasing 工具(MethPhaser/NanoMethPhase/NANOME)、讀層 clustering(cvlr)、epiallele 演化(Mazzocco)。structure≠disorder framing 成立；niche 收窄至 cis-test+LOH 耦合+5mC/5hmC 分軌。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal vs research story §7B + §7 C4"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, ism, methclone, epihet, metheor, cvlr, pycometh, structure-vs-disorder, niche]
canonical_paths: [11_external_literature/03_ism_method_comparison.md]
alias_paths: []
---

# Q3 ISM 方法是否合理有效：相似/差異工具，structure vs disorder

- **一句裁決**：framing **CONFIRM**（structure≠disorder 經 Landau 明文背書），但 **niche 收窄**——read-level haplotype-aware DMR 已被 pycoMeth/cvlr/NanoMethPhase/Lin24 佔；**零競爭的只剩 normal-anchored somatic cis-test + LOH/CN 耦合 + 5mC/5hmC 分軌**。read-read 距離矩陣**單獨**非新（cvlr 做 clustering，Metheor 做 within-read 距離）。是**方法-niche 裁決，非生物-impact 裁決**。

## 問題與內部立場

ISM = read-read 距離矩陣(L1/NHD/Bernoulli) + 階層 subclonal clustering + LOH/CN 耦合 + normal-anchored cis-test，我們宣稱量 haplotype-linked **STRUCTURE**，相對 epipolymorphism 工具量 **DISORDER**。要確認方法合理有效、有哪些相似/目標一致但不同的工具、structure-vs-disorder 是否有效、read-read 距離+cis-test niche 是否真未填。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Li, **methclone**, Genome Biol | 2014 | PMID 25260792 / 10.1186/s13059-014-0472-5 | 短讀 bisulfite；AML dx→relapse | 組合熵差 eloci（區域 4-CpG tuple 熵）；量 epiallele clonality DISORDER | CONTEXT（disorder 極）| 短讀；熵-of-tuple 非 read-read 距離；無 haplotype/LOH | high |
| Landau, Cancer Cell | 2014 | PMID 25490447 / 10.1016/j.ccell.2014.10.012 | 短讀 RRBS+WGBS；104 CLL | **PDR**=locally disordered；67.6% 變異來自 discordant reads；明文 **stochastic 非 allelic/haplotype-linked** | **CONFIRM（驗證 structure-vs-disorder 二分）** | cohort vs 單株；他們自己說 disorder≠allelic structure → 直接撐我們兩極 framing | high |
| Pan/Chen, **epihet**, Sci Rep | 2021 | DOI 10.1038/s41598-020-79627-x / PMC7801679 ⚠作者/PMID 見 09 | 短讀 bisulfite；methclone 下游 | PDR+epipolymorphism+Shannon entropy + network | CONTEXT（disorder 極）| 短讀，無 structure/haplotype/cis | high |
| Lee, **Metheor**, PLOS Comp Biol | 2023 | DOI 10.1371/journal.pcbi.1010946（PMID UNVERIFIED）| 短讀 bisulfite；928 株 benchmark | 6+ disorder 度量 + 新 **LPMD**（固定距離 CpG-CpG discordance **同一 read 內**），300× 快 | DIFFERENT-VIEW（距離-aware 但 within-read）| ⚠ Metheor「距離」= 同分子 CpG↔CpG，非 分子↔分子；銳化我們 **read-read**(between-molecule) 距離才是未填軸 | high |
| Raineri, **cvlr**, Bioinform Adv | 2023 | PMID 36726731 / 10.1093/bioadv/vbac101 | 長讀 ONT；NA12878 | **多元 Bernoulli 混合 + EM** clustering reads 成甲基亞群（無需 prior phasing）；chr15 909 peaks | DIFFERENT-VIEW（最近 read-level clustering analog）| 最近我們 clustering 步；ONT 分子-level，**但無 read-read 距離矩陣**(model-based EM 非階層距離)、無 LOH/CN、無 normal-anchored cis；germline | high |
| Fu, **MethPhaser**, Nat Commun | 2024 | PMID 38909018 / 10.1038/s41467-024-49588-0 | 長讀 ONT；germline 株+4 血液 | haplotype-specific 甲基**延伸 SNV phasing**；N50 +78–151% @ 83.4–98.7% | DIFFERENT-VIEW（鄰目標，**反因果方向** methyl→phasing）| phasing 工具用 methyl 當 input；我們讀出 structure 當 output；germline | high |
| Akbari, **NanoMethPhase**, Genome Biol | 2021 | PMID 33618748 / 10.1186/s13059-021-02283-5 | 長讀 ONT；germline trio | SNV→WhatsHap phasing→haplotype 甲基+DMA；2,205 DMR，60% ICR；96.3% concordance | DIFFERENT-VIEW（phasing→methyl，germline）| region-level haplotype DMR 非 read-read；germline imprinting | high |
| Liu, **NANOME**, bioRxiv | 2025 | PMID 40631091 / 10.1101/2025.06.29.662079 | 長讀 ONT；imprinting 控制 | consensus methyl(XGBoost) + variant+phasing → haplotype-aware ASM；+11% precision | DIFFERENT-VIEW（最近長讀 ASM pipeline 但 imprinting）| imprinting 非 somatic；無 LOH/CN-耦合 subclonal、無 read-read、無 normal-anchored cis | preprint |
| Barrett/Mazzocco, BMC Bioinf | 2017 | PMID 28743252 / 10.1186/s12859-017-1753-2 | 短讀 RRBS；tumor-normal 多區(1 肺瘤) | **Bayesian epiallele 偵測**；推 epiallele set+比例 → 重建腫瘤演化史 | CONTEXT（最近 subclonal-deconv 目標）| 目標一致但 epiallele-freq/disorder-based，短讀，region-level；無 haplotype/LOH/read-read/cis | high |

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| structure-vs-disorder framing 有效（L2）| **CONFIRM**：Landau 明文 disorder(PDR) stochastic 非 allelic-linked；methclone/epihet/Metheor 全 read-internal discordance | 🟢 none | framing 成立——disorder 工具構造上量不到 haplotype-linked structure |
| O11 AUC 0.845→0.530 coverage 校正（L1）| CONTEXT/SILENT：無外部把 epipolymorphism 當 classifier 做 coverage 校正 | 🟡 unaddressed | 我們 coverage-confound 校正是真貢獻，未被反駁 |
| read-read 距離矩陣(L1/NHD/Bernoulli) 未填（L2）| **CONFIRM 帶 1 nuance**：cvlr 做分子-level Bernoulli *clustering*(model-based 無顯式距離矩陣)；Metheor LPMD 是 *within-read* CpG-CpG 距離；無工具建顯式 read↔read 距離矩陣+階層 clustering on ONT 5mC/5hmC | 🟡 caliber-only | niche 成立；cvlr 最近(clustering 有，距離矩陣無)|
| normal-anchored cis-test 未填（L2）| **CONFIRM**：NanoMethPhase/MethPhaser/NANOME 偵測 ASM 但**無 tumor-vs-matched-normal cis vs drift 判別** | 🟢 none | 最清楚未填 |
| LOH/CN 耦合未填（L2）| **CONFIRM**：無 surveyed 工具耦合 ASM/clustering 到 LOH/CN | 🟢 none | 未填 |
| methyl→phasing 方向他處存在（Q6）| DIFFERENT-VIEW：MethPhaser 主目標 methyl→phasing；我們 methyl-rescue 同向但 somatic/LOH(germline-het 失去處)，MethPhaser 不涵蓋 | 🟡 direction(of use) | 我們 cancer/LOH 應用是 MethPhaser 未蓋的 open extension |
| 5mC/5hmC 分軌做 ASM（Q7）| **SILENT**：無 ASM/phasing 工具分軌 | 🟡 unaddressed | 真新軸（內部 5hmC marginal 為單樣本）|

## 缺口（我們缺什麼）

1. **無 vs cvlr 的 model-based clustering benchmark**：cvlr(Bernoulli+EM) 最近且 ONT-native；未證我們距離矩陣+階層 clustering 優於/互補 EM clustering。
2. **coverage-confound 領域少報**：我們 O11 校正(0.845→0.530) 暴露 disorder 工具不常控的 confound——可發 cautionary methods note；目前領域 silent 故我們既無佐證亦無反駁。
3. **無腫瘤演化/phylogeny 輸出**：Mazzocco(Bayesian epiallele) 重建多區演化史，是我們 clustering 尚未產的 deliverable；我們做 subclone *detection* 非 *phylogeny/ordering*（與 cross-region 二次打擊順序 NEGATIVE 一致）。
4. **germline imprinting 正控**：cvlr/NanoMethPhase/NANOME 都對已知 imprinting ICR 驗證；我們對 SEQC2 variant truth 驗但**未證 ISM 復現 canonical imprinted ASM**——便宜高 credibility 的正控缺。
5. **LPMD-style within-read 距離未對照**：須在文字明確區分 Metheor LPMD 與我們 read-read 距離，否則 reviewer 易混。

## 裁決

文獻**確認我們 framing 與方法定位**，每個核心 claim 衝突度 caliber-only 或無。structure-vs-disorder 二分有實證根基——Landau(2014) 明文 locally disordered 甲基 stochastic 非 allelic/haplotype-linked，整個 disorder 工具鏈(methclone/epihet/Metheor) 算 read-*internal* discordance、零 haplotype linkage。structure 極：cvlr 做分子-level Bernoulli *clustering* 但無顯式 read-read 距離矩陣/LOH-CN/cis；MethPhaser/NanoMethPhase/NANOME 做 haplotype/ASM 但為 germline phasing/imprinting(反因果用)，無 somatic LOH/CN 耦合/normal-anchored cis；Mazzocco 做目標一致的腫瘤 subclonal deconv 但短讀 epiallele-freq(disorder 族)。特定組合——**長讀 5mC/5hmC read-read 距離矩陣 + 階層 subclonal clustering + LOH/CN 耦合 + normal-anchored cis(vs drift) test**——無任何 verified 工具佔；**normal-anchored cis-test + LOH/CN 耦合最清楚未填（零競爭）**。誠實限制：這是 *方法-niche* 裁決非 *生物-impact* 裁決——niche 未填不代表訊號可用（Q5/Q7 顯示 ASM 真實但非判別、對 FP filter anti-discriminative、單樣本 HCC1395 ⭐3 單 pipeline）。我們也缺對照論文都報的便宜正控(imprinting ICR 復現、cvlr head-to-head)，故 structure claim 目前對 variant ground-truth 已驗但未對領域標準表觀 benchmark 驗。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ 研究故事 §7B 工具地景：[../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md](../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md)
