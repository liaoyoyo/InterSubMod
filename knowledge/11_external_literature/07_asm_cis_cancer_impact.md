---
id: ism-kb-11-external-literature-asm-cis-cancer-impact
name: "Q7 ASM / cis allele-specific methylation 與癌症影響"
description: "癌症 ASM 頻率/方向(Do2020 +5x hypo-dominant)、CN confound(Martin-Trujillo CN 82-92%)、TSG silencing(BRCA2/VHL hyper)、cis-ASM、功能後果。我們 somatic-controlled HP-axis(稀有 0.34% 無方向)是口徑差非矛盾；CN-confound 被 Martin-Trujillo 最強背書；BRCA2 hypo≠canonical hyper 須標；private pattern 被 Sakamoto 確認、recurrence 被 POG 證但我們 underpowered。"
status: active
last_verified: 2026-06-08
content_nature: paper-derived
doc_type: explanation
verified_scope: "external papers via WebFetch wf_37b2cc97-663 (2026-06-05); BRCA2 cis口徑 reconciled to 06-07 fast_cnv_validation.json (BRCA2=subclone/copy-confounded〔HP1-1=somatic subclone tag 非 copy；非 CN-dosage〕, chr17/TBC1D16 lone clean cis) per convergence audit wf_9e169112-573 (2026-06-08) + 06-09 capstone reconcile (% split not-robust; ledger 94)"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
  - ism-kb-11-external-literature-paper-readiness-convergence
tags: [literature, asm, cis, copy-number-confound, tsg-silencing, brca2, do2020, martin-trujillo, oneill, private-recurrence]
canonical_paths: [11_external_literature/07_asm_cis_cancer_impact.md]
alias_paths: []
---

# Q7 ASM / cis allele-specific methylation 與癌症影響

> 🔴 **2026-06-08 新鮮度修正（本檔 06-05 建立時用 06-03 tsg 快照；06-07 cis 分析已收緊，下文 BRCA2 cis 敘述部分過時）**：
> - **BRCA2 (chr13:32,315,128) 不再是乾淨 cis 錨點** —— 06-07 重分析列為 **predominantly subclone/copy-confounded**（**HP1-1 是 longphase-S 的 somatic SUBCLONE tag〔germline-H1 + 帶 somatic ALT, HaplotagStrategy.cpp:505-516〕，非 copy；非 CN-dosage——dosage 已決定性 REFUTED**）。HP-axis Δβ=−0.122 ≈ d_copy −0.11〔subclone vs germline 同 REF〕+ d_within −0.023〔focal allele 真 cis 殘餘, perm p=0.024, n=19/19 邊際·未確立〕。誠實口徑 = 「**subclone/copy 主導；焦點突變 cis 在單樣本與 subclone 背景不可乾淨分離（CAMDAC 同此限制）**」——**勿寫精確 %（split 不 robust，分解跨位點不閉合：BRCA2 殘差 9%、chr17 達 40%）**，勿寫「copy number」，勿寫「BRCA2 真 cis-driven」，亦勿寫「純 copy-artifact」（allelic 甲基差是真的）。
> - **唯一乾淨 cis exemplar = chr17:79991120 (TBC1D16)**（within d_within=0.142 > HP-axis 0.122，perm p=0.001, CGI；單 locus 單樣本 n=16/14, nominal p → 「稀有但真實 exemplar」）。
> - 下文「**60 T3 + 4 survivor**」是 screening intermediate，06-07 Bonferroni+copy-test 後**塌成 1 個乾淨 cis**（816→10 Bonferroni）。
> - ⚠ **Martin-Trujillo (CN-DOSAGE 解釋 82-92%) 不 corroborate BRCA2 reclassification** —— 06-07 FAST 分析**決定性 REFUTED copy-DOSAGE**（MW neutral-vs-gain p=0.6183；signed ρ(\|Δβ\|,CN)=−0.083 反向）；BRCA2 的 confound 是 **subclone（細胞群差異）非 CN-dosage**。引 CN-dosage 文獻當 corroborate 是 **category error**。Martin-Trujillo 的正確角色 = 我們**控掉**的 confound 警告（dosage falsification back-stop 它），**非** BRCA2 reclassification 的佐證。正確對照 = **subclone-specific methylation**（pubmed 25066126/24356097）+ **CAMDAC**（單樣本 focal-cis vs subclone 不可分離, biorxiv 2020.11.03.366252）。
> - 完整收斂與 paper-readiness → [10_paper_readiness_convergence.md](10_paper_readiness_convergence.md)（HD-3）。本檔保留原 06-05 敘述供 audit trail，引用時以本 banner 為準。

- **一句裁決**：方法 **CONFIRM**、方向 **部分 confirm，兩個誠實 caveat**。CN-confound 是最強外部對齊（Martin-Trujillo：CN 解釋 82–92% 腫瘤 imprinted allelic 甲基 = 我們 HP-axis 設計 held-constant 的確切 artifact；ASM×CN ρ=−0.055 確認我們非 CN-driven）。cis-driven + BRCA2 hypo 與 Do2020(hypo-dominant, SNP-cis) 對齊；private/非復發 被 Sakamoto 直接確認。**caveat 1**：Do2020「ASM +5x」是 cross-individual SNP-anchored bisulfite-DMR 率，我們「稀有 0.34%」是 within-sample somatic-controlled read-率——**不同分母，口徑差非矛盾**。**caveat 2**：BRCA2 方向 **hypo**，與 canonical TSG **hyper**-silencing 相反——我們 HP-axis Δβ 是等位結構差(含 allele+copy)非 promoter-mean silencing，**不可當 TSG-silencing 證據**。recurrence 是決定性內部限制：POG(189) 偵 RET/CDKN2A 復發 aDMR，我們單 pipeline 6 樣本(主 HCC1395 ⭐3) 偵不到，private 0/38 誠實但 underpowered。

## 問題與內部立場

文獻對癌症 ASM 頻率/方向(Do2020 +5×/hypo)、CN-confound(Martin-Trujillo CN 主導)、TSG silencing(BRCA2/VHL)、cis-ASM、功能/表現後果的說法；我們更乾淨的 somatic-controlled HP-axis 結果(稀有、非方向)在哪確認/矛盾/口徑差。

**內部錨點（L1-L2，tsg snapshot @06-03）**：ZAR1L/BRCA2 ASM 真實但非方向/非判別/coverage-modulated；BRCA2 promoter chr13:32,315,128 **HP-axis Δβ=−0.122**(n=197, 含 allele+copy)/純 **within-somatic d_within=−0.023**(perm p=0.022, 小但真, **hypo**)；ALLELE −0.099；**d_copy=−0.11**；normal-anchored cis-test 初判 d_cis=−0.142 vs d_drift≈−0.022(表面值；06-07 copy-partition 後 = subclone/copy 主導，真 cis = within-somatic −0.023，乾淨 cis 改 chr17)；genome cis-scan **816 → 60 T3 cis-candidate + 4 perm-survivor**；cross-sample **6/6 excess-over-null**(mean 0.168, 3 癌種)、somatic **private 0/38**；ASM×CN HP-axis partial **ρ=−0.055**(非 CN-driven)；5mC-driven，5hmC marginal 0.03–0.07。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Do, Genome Biol | 2020 | PMID 32594908 / 10.1186/s13059-020-02059-3 | 短讀 bisulfite(capture+WGBS)；16 癌 vs 24 normal | ASM **+5×** in cancer(myeloma 5×/lymphoma 8.5×/GBM 9×)；**hypo-dominant**(49–76% 為 allele-specific loss)+16–25% focal hyper gain；15,112 ASM DMR；cis by SNP(CTCF/TF)| **DIFFERENT-VIEW(頻率)+CONFIRM(機制)**| 他們「+5×」是 cross-individual SNP-anchored bisulfite ASM 率；我們是 within-sample somatic-controlled read-率(0.34%)——**不同分母**。但 hypo-dominance+cis-by-SNP 對齊我們 hypo BRCA2+cis 倖存者 | high |
| Martin-Trujillo, Nat Commun | 2017 | PMID 28883545 / 10.1038/s41467-017-00639-9 | Array(450K+SNP6 CN)；287 癌株+TCGA | imprinted DMR **82–92% 甲基異常由 CN 解釋**非真表觀(僅 8.1–17.6% CN-independent)。「CN 是主 dictator」| **CONFIRM（我們 CN-confound+HP-axis 設計）**| 他們警告 CN 驅動 apparent allelic 甲基；我們 HP-axis 設計 held-constant CN，ASM×CN ρ=−0.055 顯非 CN-driven。正是我們方法控的 confound | high |
| O'Neill/POG, Cell Genomics | 2024 | PMID 39406235 / DOI 100538 或 100674 ⚠見 09(SSRN preprint 4743887)| 長讀 ONT；189 腫瘤/41 normal | long-range phasing 找 aDMR+ASE 含**復發 RET/CDKN2A aDMR** | **DIFFERENT-VIEW(復發)**| 最近方法同儕(長讀/paired/癌)；他們**復發**跨 189，我們 **private**(0/38)。caliber gap 決定性：cohort vs 單樣本。PMID/DOI 須標 ⚠ | high(識別碼未完全驗) |
| Sakamoto, Nat Commun | 2022 | PMID 35710642 / 10.1038/s41467-022-31133-6 | 長讀 ONT+Illumina；20 NSCLC paired | haplotype-biased 表現+DMR **patient-specific 非復發**：avg ~32 haplotype-biased 基因/案(HP1-vs-HP2 DMR)。**未**分析 MLH1/BRCA1/RAD51C(scout 宣稱有誤)| **CONFIRM（我們 private/非復發）**| 強 match：長讀 within-sample HP1-vs-HP2(≈我們 HP-axis) 找 **private case-specific** hap-DMR——正是我們發現。scout「復發 driver」framing 錯，本論文實撐我們 private | high |
| Schalkwyk, Am J Hum Genet | 2010 | PMID 20159110 / 10.1016/j.ajhg.2010.01.014 | Array(SNP)；normal 組織 | ASM **genome-wide 廣泛**(>35,000 site)、多 **cis-driven by 序列變異**、個體+組織異質 | **CONTEXT(baseline 率)+CONFIRM(cis,非方向)**| 立 normal cis-ASM baseline；我們非方向 read-level 與廣 cis-ASM 背景對齊；normal-tissue array 定 floor | high |
| Do/Greally/Tycko(review), Genome Biol | 2017 | PMID 28629478 / 10.1186/s13059-017-1250-y | Review(mQTL+hap-ASM)| 定義 **mQTL(cross-individual 相關) vs sequencing hap-ASM(within-individual 等位比較)** 區別 | **CONTEXT(caliber framing)**| 直接解釋為何我們 single-individual within-sample HP-axis 是不同 caliber 的證據 vs population mQTL/Do2020 cross-individual 率。我們「caliber 差」論點的概念骨幹 | high |
| Hesson, Mol Diagn Ther | 2022 | PMC9626413 / 10.1007/s40291-022-00614-1（PMID UNVERIFIED）| Array+pyroseq；BRCA1 PDAC | BRCA1 promoter 甲基是**罕見第二打擊**：mean 3.62%；**LOH 主導(63%)**、somatic mutation 10.5% | **CONFIRM（稀有+緩和 TSG-silencing 預期）**| 撐「稀有 ASM」+TSG silencing-by-methylation 比 LOH 少。方向：此為 **hyper** silencing(稀有) 非我們 hypo BRCA2 | high(PMID 未驗)|
| Herman & Baylin(review), J Clin Invest | 2001 | PMC289180 / 10.1172/JCI9462（UNVERIFIED）| Review/array；tumor | VHL promoter **hyper**-甲基當 allele-specific 第二打擊：>30% non-LOH RCC 沉默非突變等位 | **REFUTE-方向(vs 我們 hypo)/CONTEXT**| canonical TSG silencing 是 **hyper**；我們 BRCA2 軸是 **hypo**。我們發現須對標的 canonical 方向參考 | high(secondary/review；識別碼未驗)|
| Akbari, NanoMethPhase, Genome Biol | 2021 | PMID 33618748 / 10.1186/s13059-021-02283-5 | 長讀 ONT；normal/tool | phase 5mC+偵測 ASM @~10× | **CONTEXT(方法源流)**| 長讀 HP-axis ASM 技術骨幹；ISM extend(加 5mC/5hmC 分軌+read-read 距離+somatic-controlled HP1-vs-HP1-1+normal-anchored cis)| high |

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| ASM 真實但**稀有**(strong-ASM 0.34%) — 單樣本 ⭐3 | Do2020 ASM **+5×**；Schalkwyk >35k(normal)；Hesson TSG 甲基稀有 | 🟡 **caliber-only**（分母不同：within-sample somatic-controlled read-率 vs cross-individual SNP-anchored bisulfite DMR-率）| 非真矛盾。Do2020「increase」測根本不同軸；我們稀有在更嚴 somatic-controlled 軸。**明標；勿稱 Do2020 矛盾我們** |
| **非方向**(strong-ASM hypo 44%/hyper 56% 無偏) ⭐3 | Do2020: hypo-**dominant**(49–76%)+focal hyper gain；Herman canonical TSG: hyper | 🔴 **direction(部分)** | genome-wide aggregate 平衡，但 **BRCA2 within-somatic hypo**(−0.023, p=0.022) 與 Do2020 hypo-dominance 一致非 canonical TSG-hyper。誠實：整體平衡，唯一檢的 TSG(BRCA2) locus-level hypo |
| **CN 非驅動**(ASM×CN HP-axis ρ=−0.055；設計 held CN)⭐3 pilot | Martin-Trujillo: CN 解釋 82–92% imprinted allelic 甲基 | 🟢 none(互撐)| **最強外部佐證。** 他們警告正是我們 HP-axis 設計擊敗的 confound。**直接 CONFIRM** |
| **cis-driven ASM** 真但乾淨案稀少(copy-control 後唯一乾淨 = **chr17/TBC1D16**；BRCA2 表面 d_cis=−0.142 但 copy-partition 後僅 within-somatic −0.023 邊際；60 T3 cis-candidate→copy-control 後 chr17/18 survive)⭐3 | Do2020+Do2017+Schalkwyk: ASM cis-driven by SNP/haplotype | 🟢 none | 機制被文獻確認；乾淨 cis 在單樣本+subclone 背景稀少(CAMDAC 同限制)。caliber：我們 within-individual，文獻 population |
| somatic ASM **private**(0/38；6/6 excess，3 癌種)⭐3 單 pipeline | Sakamoto22: patient-specific hap-DMR(~32/案)；O'Neill24: **復發** RET/CDKN2A | 🔴 caliber+direction-of-recurrence | Sakamoto CONFIRM private；O'Neill 在 cohort REFUTE(189 偵復發)。**我們 6 樣本單 pipeline 偵不到 O'Neill 的復發**——關鍵誠實限制。現象-replicate(6/6 excess) 對齊所有長讀同儕 |
| **5mC-driven，5hmC marginal**(0.03–0.07)⭐3 | 文獻大多 silent(Do2020/Sakamoto/O'Neill 測 5mC；5hmC 分軌少報)| 🟡 unaddressed | 真 gap/niche——長讀 5mC/5hmC 分軌做 ASM 未探。新但單樣本低 tier |
| BRCA2 HP-axis Δβ=−0.122 / within-somatic −0.023 ⭐3 | Herman canonical: TSG silencing hyper；Hesson: BRCA1 甲基稀有 | 🔴 direction(vs canonical)+caliber | 我們 hypo **矛盾 canonical TSG-hyper-silencing**。須標：我們軸是 HP-linked allelic Δ(含 allele+copy)非 promoter-mean hyper-silencing，**非同量測**。**勿把 BRCA2 當 TSG-silencing 證據** |

## 缺口（我們缺什麼）

1. **Cohort-scale recurrence 偵測**：O'Neill(189)+Do2020 能 call 復發癌基因 aDMR(RET/CDKN2A)；我們 6 樣本單 pipeline 只能立 private+現象-replication。無法分「無復發」vs「underpowered」——最大 gap，cross-sample tier-capped ⭐3 之因。
2. **cross-individual mQTL/population ASM anchoring**：Do2017 區分 within-individual hap-ASM(我們) vs population mQTL；我們無 population mQTL 軌 anchor 我們 private loci 是 cis-genetic vs somatic-acquired(只單 normal-anchored cis-test)。
3. **功能/表現後果**：O'Neill/Sakamoto 配 ASM+allele-specific **expression**(RNA-seq)；我們**無表現/RNA 讀出**——ASM 純甲基-結構，功能後果未立。
4. **GWAS/regulatory-SNP overlap**：Do2020 map 1,842 ASM-SNP 到 GWAS peak(疾病相關 payoff)；我們未交集 cis-candidate 與 GWAS/regulatory 註解。
5. **獨立 caller/多 pipeline**：cross-sample 全 ClairS-TO(單 pipeline)→tier ⭐3。Martin-Trujillo-style 顯式 CN-分解跨全 loci(非只 pilot) 也未 genome-wide 做。
6. **Ground-truth ASM benchmark**：無 orthogonal(bisulfite/array) 驗我們長讀 ASM call 於同樣本——Akbari/NanoMethPhase 源流在但未 cross-validate。

## 裁決

方法**確認、方向部分確認帶兩誠實 caveat**。CN-confound 是最強外部對齊：Martin-Trujillo(CN 解釋 82–92% 腫瘤 imprinted allelic 甲基) 正是我們 somatic-controlled HP-axis 設計 neutralize 的 artifact，ASM×CN pilot(ρ=−0.055) 確認非 CN-driven = 真可防守 niche。cis-driven+hypo-at-BRCA2 對齊 Do2020 hypo-dominant SNP-cis ASM，private/非復發 pattern 被 Sakamoto22(patient-specific hap-DMR) 直接確認。ISM niche——長讀 read-read 距離+somatic-controlled HP1-vs-HP1-1+normal-anchored cis-test+5mC/5hmC 分軌——**成立**：無 verified 同儕合四者，5hmC 分軌文獻幾乎未 addressed。兩 caveat 須明說：(a) **頻率是 caliber-only 非矛盾**：Do2020「ASM +5×」是 cross-individual SNP-anchored bisulfite-DMR 率，我們「稀有 0.34%」是 within-sample somatic-controlled read-率——不同分母，故勿稱 Do2020 矛盾我們亦勿稱我們 refute Do2020；(b) **BRCA2 方向 hypo，違 canonical TSG hyper-silencing**(Herman/VHL、Hesson/BRCA1)：我們 HP-axis Δβ 是等位結構差(含 allele+copy)非 promoter-mean silencing，**勿把 BRCA2 當 TSG-silencing 證據**——只當過 permutation cis-test 的 cis-driven allelic-methylation 例。決定性內部限制是 **recurrence**：O'Neill24 的 189-tumor cohort 偵 RET/CDKN2A 復發 aDMR，我們單 pipeline 6 樣本(主 HCC1395 ⭐3，⭐4 需 COLO829 truth) 偵不到——private 0/38 誠實但對 cohort-scale recurrence underpowered。淨：confirmed-方法、partially-confirmed-方向、niche-成立、但因單樣本/單 pipeline 牢牢封在 characterization-tier。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ 研究故事 §5 BRCA2 錨點 + §7 C1/C2/C6：[../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md](../../docs/concepts/2026/06/20260603_研究故事與困難敘述_甲基haplotype論文_內外證據_01.md)
