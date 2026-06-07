---
id: ism-kb-11-external-literature-finding-methyl-diff-loci
name: "Q2 癌症甲基差異位點：做法與困難"
description: "找癌症 DMR/ASM 位點的方法地景與困難。我們的『tumor 5mC 太混亂、抓不到大量一致位點』是公認 epipolymorphism/heterogeneity；但 Lin/POG cohort 找到 ~23.6K aDMR/腫瘤，部分 REFUTE『找不到很多』。O11 coverage artifact 被 Scherer20 機制背書。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal vs research story §2/§7 C3 + O11"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, dmr, asm, epipolymorphism, heterogeneity, pdr, coverage-confound, methylation]
canonical_paths: [11_external_literature/02_finding_methyl_diff_loci.md]
alias_paths: []
---

# Q2 癌症甲基差異位點：做法與困難

- **一句裁決**：困難 **CONFIRM**（tumor 甲基失序/隨機是公認、有名字 epipolymorphism）；但「找不到很多」**部分 REFUTE**——Lin/POG 長讀 cohort 找到 ~23.6K aDMR/腫瘤（79% 落 CNV/LOH）。**我們的無法定位糾纏於單樣本+單 pipeline underpowering，非硬生物天花板。** O11 epipolymorphism 崩潰被 Scherer20（coverage-依賴）機制背書 = 真方法學 win。

## 問題與內部立場

要找癌症 DMR/ASM 位點有哪些穩健方法；我們撞到的困難（tumor 甲基太失序，抓不到大量清楚、一致於 HP/somatic 的位點，個案 IGV 可見但群體統計失敗）是否為文獻公認的 epigenetic heterogeneity / epipolymorphism。

**內部錨點（L1-L2）**：tumor 5mC 混亂；strong-ASM 稀有 **0.34%** 無方向；O11 epipolymorphism 特徵 AUC **0.845→0.530**（n_reads/coverage 校正後崩 = 表面訊號是 coverage artifact）。[src: 研究故事 §2/§7 C3]

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Landau, Cancer Cell | 2014 | PMID 25490447 / 10.1016/j.ccell.2014.10.012 | 短讀 WGBS+RRBS；104 CLL+26 normal | 提出 **PDR**；**67.6±3.2% 瘤內甲基變異來自 discordant reads**；失序 genome-wide、stochastic | **CONFIRM** | cohort/bisulfite vs 我們單株 ONT；非 haplotype-resolved | high |
| Landan, Nat Genet | 2012 | PMID 23064413 / 10.1038/ng.2442 | 短讀 bisulfite；纖維母細胞+組織 | **epipolymorphism** 奠基定義；DMR 形成 stochastic+near-deterministic | **CONFIRM**（定義）| 非 somatic-controlled；定義我們 O11 用的度量 | high |
| Scherer, NAR | 2020 | PMID 32103242 / 10.1093/nar/gkaa120 | 短讀 WGBS(sim+real)；方法 benchmark | PDR/MHL/epipolymorphism/entropy 系統比較；**epipolymorphism & entropy 需 >10×**；隨 read-length/CpG-density 變 | **CONFIRM（我們 artifact 的機制）** | 方法非癌生物；直接命名我們撞的 coverage confound | high |
| Jenkinson (**informME**), BMC Bioinf | 2018 | PMID 29514626 / 10.1186/s12859-018-2086-5 | 短讀 WGBS；3 對肺正常/癌 | Shannon-entropy(NME/dNME) Ising；癌「**全域 shift 到失序態**」；dNME>DSS | **CONFIRM** | n=3 小；region-level entropy 非 haplotype；dNME 我們沒有 | high |
| Li (**methclone**), Genome Biol | 2014 | PMID 25260792 / 10.1186/s13059-014-0472-5 | 短讀 bisulfite；白血病 dx→relapse | **epiallele clonality SHIFT**（組合熵差，eloci）| **DIFFERENT-VIEW** | 量 DISORDER/clonality-shift 非 haplotype-linked cis；Q3 對照 | high |
| Chen (**epihet**), Sci Rep | 2021 | DOI 10.1038/s41598-020-79627-x / PMC7801679 ⚠PMID 見 09 | 短讀 WGBS/RRBS；AML | R 工具算 **PDR+epipolymorphism+Shannon entropy** ITH | **DIFFERENT-VIEW** | DISORDER 量化 bisulfite，非 haplotype/cis | high |
| **Lin/POG**(長讀), Cell Genomics | 2024 | PMID 39406235 / 10.1016/j.xgen.2024.100674 ⚠見 09 | **長讀 ONT**；189 腫瘤/43 matched normal/26 癌種 | **84% CpG 可 phase**；mean **23.6K aDMR/腫瘤**(~5× normal 4.7K)；**79% aDMR 落 CNV/LOH**；復發 promoter aDMR(RET/CDKN2A/BRCA1/RAD51C)| **CONFIRM + 部分 REFUTE「找不到很多」** | cohort+ground-truth panel vs 單株；他們**確實**定位大量 aDMR → 我們缺口是 cohort-power | high |
| Akbari (NanoMethPhase), Genome Biol | 2021 | PMID 33618748 / 10.1186/s13059-021-02283-5 | **長讀 ONT**；normal 株 ~10× | 全基因組 ASM via phasing；2,205 DMR(NA19240)；recover 60% ICR | **CONTEXT**（方法骨幹）| normal/imprinting；SNP-anchored phasing（LOH 會失去的 anchor）| high |
| Snajder (**pycoMeth**), Genome Biol | 2023 | DOI 10.1186/s13059-023-02917-w（PMID 37076875*）| **長讀 ONT**；DMR 方法 | **haplotype-aware Bayesian read-level changepoint 分段 + consensus DMR**；用 bisulfite 工具忽略的 read-level 資訊 | **DIFFERENT-VIEW（最近方法對手）** | read-level+haplotype-aware 與我們重疊，但無 read-read 距離 clustering/somatic cis/5mC-5hmC 分軌 | high |
| Zhang/Sun（11 癌 DNAm haplotype map）, Cell Reports | 2025 | PMID 40849904 / 10.1016/j.celrep.2025.116197 | 短讀 WGBS(MHB)；110 腫瘤/11 癌種 | **81,567 MHB**，高癌種特異；MHB 與表現相關（獨立於 mean methylation）| **CONFIRM（haplotype-block 跨癌公認）** | cohort/bisulfite；MHB(region) 非單 read clustering | high |

*pycoMeth PMID 未直接 fetch（redirect/auth）；DOI 已 crossmark 驗證為權威識別碼。

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| tumor 5mC 混亂、抓不到大量 HP/somatic-一致位點（L1，單株）| Landau14/Landan12/Jenkinson18 CONFIRM 失序真實+genome-wide；**但 Lin24 定位 ~23.6K aDMR/腫瘤** | 🔴 direction（部分）| 困難是公認**現象**；我們的**無法定位大量**部分是 **power/scope 限制**（單株、無 cohort/ground-truth panel），非純內在 |
| strong-ASM 稀有 0.34% 無方向（L2）| 文獻 silent on 確切率；Landau14 失序方向不可知（與無偏一致）| 🟡 caliber/unaddressed | 無外部反駁；確切率新+單樣本不可佐證 |
| O11 AUC 0.845→0.530 coverage 校正後 = artifact（L1）| **Scherer20 CONFIRM**：epipolymorphism/entropy 需 >10× 且 coverage-sensitive | 🟢 none | **強 confirmed** — 正是公認失效模式；真方法學 win |
| read-read 距離+cis+5mC/5hmC 分軌 = unfilled niche（claim）| methclone/epihet 量 DISORDER（對照成立）；pycoMeth haplotype-aware read-level DMR（重疊）；Lin24/NanoMethPhase 做 phased aDMR | 🟡 caliber-only | niche **部分成立**：disorder-vs-structure 有效，但 pycoMeth+Lin24 已佔「長讀 haplotype-aware read-level DMR」；真正未填 = **somatic-controlled normal-anchored cis-test + 5mC/5hmC 分軌**，非 read-level-clustering 本身 |
| ALLELE-axis 被 germline allelic 甲基 confound；只 HP-axis somatic-controlled（L2）| NanoMethPhase/NANOME 做 SNP-anchored ASM(imprinting) = 我們 null 控制的 baseline | 🟢 none（我們 extend）| 我們 germline-het null 是這些工具不 foreground 的 rigor step |
| ASM 真實但非方向/非判別/coverage-modulated（L2，HCC1395）| Lin24: 79% aDMR 在 CNV/LOH（與 copy-component confound 一致）；無外部 TP/FP-filter claim | 🟡 unaddressed | filter-NEGATIVE 內部-only；無文獻 reopen |

## 缺口（我們缺什麼）

1. **Cohort power + ground-truth panel**：Lin24 定位 ~23.6K aDMR *因為* cohort scale + recurrence test；我們「找不到很多」糾纏單株/單 pipeline。無法分辨「腫瘤本質不可定位」vs「我們 underpowered」。
2. **CNV/LOH-aware aDMR 計算**：Lin24 headline「**79% aDMR 落 CNV/LOH**」是我們 copy-component(d_copy=−0.11) confound 的最乾淨外部框架——應改採 CNV/LOH-分層 aDMR 率，而非 flat 0.34%。
3. **dNME/差異熵框架**（informME）：我們用 read-read 距離+clustering，沒有 per-region 差異常態化甲基熵(dNME/JSD tumor-vs-normal) 通道可佐證 O11。
4. **WSH score battery**（Scherer20）：未對 canonical qFDRP/PDR/MHL/epipolymorphism battery benchmark；跑 qFDRP(最 coverage-robust) 可直接壓力測我們訊號是否在 epipolymorphism 崩潰處仍存活。
5. **vs pycoMeth head-to-head**（真方法對手）：「unfilled niche」在證明 read-read-distance clustering + normal-anchored cis-test 能給出 pycoMeth Meth_Comp 沒有的東西前不可防守。
6. **5mC/5hmC 分軌真正未被利用**：上述 ASM/DMR 工具皆不分軌；我們「ASM 5mC-driven / 5hmC marginal」是本角度最可防守的 novelty。

## 裁決

核心困難——**tumor 甲基太失序/隨機，難乾淨定位大量 HP/somatic-一致位點**——**領域堅實公認**非個別：Landau14(67.6% 變異來自 discordant reads，genome-wide，stochastic)、Landan12(epipolymorphism)、Jenkinson18(癌=全域 shift 到失序)、加 disorder 工具族(methclone/epihet/Scherer WSH battery) 全命名量化。**Scherer20 獨立確認我們 O11 崩潰的確切機制**(epipolymorphism/entropy coverage-依賴 >10×) = 真外部佐證的方法學 catch 非弱點。但文獻也**部分 refute「找不到很多」**：Lin24 長讀 cohort *確實*定位 ~23.6K aDMR/腫瘤(79% in CNV/LOH)，故我們無法定位糾纏於單株/單 pipeline underpowering 非硬生物天花板。ISM **niche 只在窄化形式存活**：read-level haplotype-aware DMR/ASM *已被* pycoMeth/NanoMethPhase/NANOME/Lin24 佔，故「read-read 距離+clustering」單獨非未填。可防守、仍開放的差異化是 (a) **somatic-controlled normal-anchored cis-test** 含 germline-het null、(b) **5mC/5hmC 分軌**——皆非上述工具所做。誠實限制：每個定位 claim 皆單樣本 HCC1395(⭐3，無 COLO829 truth)+單 ClairS-TO pipeline，尚不能證優於 pycoMeth 或匹敵 Lin24 cohort 定位力；niche 是**有真實 gap 支撐的假說，非已證的 win**。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ Q3 工具對比：[03](03_ism_method_comparison.md) ｜ 既有 O11 NEGATIVE：[../09_conclusions/00_index.md](../09_conclusions/00_index.md)
