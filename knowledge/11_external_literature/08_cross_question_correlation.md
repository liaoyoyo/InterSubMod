---
id: ism-kb-11-external-literature-cross-question-correlation
name: "跨問題依賴鏈 + 可信度衝突矩陣"
description: "把 7 問題串成『理想情境鏈』(找位點→助 phase→ISM 改造→解 subclone→分 TP/FP；底層 ASM/cis、長讀使能)，逐環標 REAL/BLOCKED/POSSIBLE；能力 triage、可信度×衝突矩陣、最可防守論文主張 + 3 reviewer 脆弱點防守。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "cross-synthesis of 7 angle analyses from wf_37b2cc97-663 (2026-06-05); internal claims from research story + master_draft"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-goal-alignment-methods
  - ism-kb-11-external-literature-finding-methyl-diff-loci
  - ism-kb-11-external-literature-ism-method-comparison
  - ism-kb-11-external-literature-subclone-landscape
  - ism-kb-11-external-literature-methyl-tp-fp-germline-somatic
  - ism-kb-11-external-literature-methyl-assisted-phasing
  - ism-kb-11-external-literature-asm-cis-cancer-impact
  - ism-kb-11-external-literature-references-provenance
  - ism-kb-11-external-literature-paper-readiness-convergence
tags: [literature, dependency-chain, capability-triage, credibility-matrix, reviewer-defense, thesis, synthesis]
canonical_paths: [11_external_literature/08_cross_question_correlation.md]
alias_paths: []
---

# 跨問題依賴鏈 + 可信度衝突矩陣

- **一句結論**：你的理想鏈**不是一條完整可走的路，而是一條「前段使能 REAL、中段白地 LIVE、後段死路 BLOCKED」的分段鏈**。長讀使能(L0)穩、找大量乾淨位點(L1)死、LOH 甲基 phasing(L2)是最強白地、ISM 方法(L3)distinct 但未證優、subclone(L4)只能 characterize 不能 reconstruct、分 TP/FP(L5)已死終止；底層 ASM/cis(U1)真實但不可用、不可宣稱 recurrence。

## 1. 跨問題依賴鏈（逐環裁決）

> 理想鏈：**找甲基差異位點(Q2) → 助 phasing/tagging(Q6) → ISM 改造成做這些(Q3) → 解 subclone(Q4) → 分 TP/FP(Q5)；底層 ASM/cis 生物(Q7)、長讀使能(Q1)。**

| 環 | 內容 | 裁決 | 證據 |
|---|---|---|---|
| **L0** | 長讀使能 read-level haplotype-linked 甲基(Q1)→全鏈 | 🟢 **REAL** | 內部 ONT 5mC/5hmC modBAM 矩陣可行；外部 CONFIRM(L3) NanoMethPhase(~10× ASM phasing)/POG-Lin24(84% CpG 可 phase)/de Abreu(ONT 獨特-locus)。array/WGBS 做不到 read-level phasing。**無衝突、使能前提 field-validated** |
| **L1** | 找位點：定位**大量、乾淨、可群體統計**的 HP/somatic-一致甲基位點(Q2)| 🔴 **BLOCKED（「大量乾淨群體統計」）/ 🟢 REAL（「稀有、真實、IGV 可見」）**| 內部 NEGATIVE：5mC 混亂、群體統計失敗、strong-ASM 0.34% 無方向、O11 0.845→0.530=coverage artifact。外部**分裂**：Landau/Landan/Jenkinson CONFIRM 失序真實+stochastic(67.6% 變異來自 discordant reads)；**但 Lin24(POG 189) 部分 REFUTE「找不到很多」(定位 ~23.6K aDMR/腫瘤，79% in CNV/LOH)**。Scherer20 確認 O11 coverage-confound 機制=真 win。**「單樣本群體統計找大量乾淨位點」夢死；「位點存在、稀有、cohort-scale 可定位」現實活但超我們當前 power** |
| **L2** | 甲基差異 → 助 phasing/tagging(Q6)| 🟡 **POSSIBLE（最強 LIVE 白地）— 尚非 REAL** | **單一最強 LIVE 夢**，但支持是 *可行性-only* 非 *tumor-已證*。內部：methyl-rescue pilot 存在；LOH-constrained phasing Grade A(n=7, p=0.0078; NG=2 inner 93-99% 6/6)。外部(L3)：MethPhaser/NanoMethPhase/HapBridge 證確切「甲基 tag unphase read」機制**但全 germline 株、明確留 cancer/LOH future**；LOH-as-phasing 被 ≥3 癌症論文(Wakhan/HiCancer/Sakamoto)背書但走 CNA/SV/Hi-C **非甲基**。機制已證(germline)、LOH-訊號已證(非甲基)，**但 tumor/LOH 甲基-phasing 無人做(含我們)**。真未填。POSSIBLE 非 REAL：無 switch-error/N50 head-to-head、無 LOH 內 phase-truth |
| **L3** | ISM 改造成做這些(Q3)：read-read 距離+clustering+cis-test 當方法 | 🟢 **REAL（distinct 方法）/ 🟡 POSSIBLE（*更好*的方法 — niche 成立但對 incumbents 未證）**| 定位健全、*組合*未佔。structure-vs-disorder 有實證(Landau 明文 disorder stochastic 非 haplotype-linked)。**但 niche 只在窄化形式存活**：read-level haplotype-aware DMR/ASM *已被* pycoMeth/cvlr/NanoMethPhase/NANOME/Lin24 佔。真未填 cell：**(a) normal-anchored somatic cis-test(零競爭) (b) LOH/CN 耦合(零競爭) (c) 5mC/5hmC 分軌(文獻 silent)**。read-read-距離*單獨*非未填(cvlr clustering、Metheor within-read 距離)。**無 vs pycoMeth/cvlr head-to-head**→「ISM 加值」POSSIBLE-未證 |
| **L4** | ISM 解 subclone(Q4)| 🔴 **BLOCKED（「reconstruction」）/ 🟡 POSSIBLE（「characterization」）**| 內部 mixed：HPFineNGroups ⭐4→⭐3(pipeline-dependent)；AF→NGroups +0.705~0.787(7/7) POSITIVE；**cross-region 二次打擊順序 NEGATIVE**。外部 CONFIRM-精神：read-level 甲基多樣性追蹤 clonal(Landau/Li/Hong 鐘)；長讀**可**重建 subclone 軌跡(Lee25)。**但領域「reconstruction」金標準我們都缺**：phylogeny/clock(BitPhylogeny/EVOFLUx n=1,976)、CN-aware deconv(CAMDAC/TRACERx)、單細胞 ground-truth(Lee25)。我們 *flat clustering* 非 *phylogeny*。**「reconstruct」誇大；「characterize subclonal structure」可防守**。二次打擊順序夢死(NEGATIVE，無外部 read-level 二次打擊順序 claim 可 reopen)|
| **L5** | subclone/甲基 → 分 TP/FP(Q5)| 🔴 **BLOCKED — concluded NEGATIVE，parked，不 reopen** | **最清楚死夢。** 內部 NEGATIVE(L2，難得多株)：LOSO 100% circularity；\|Δβ\| AUC=0.5049 中性；**strong-ASM 在 FP 富集 5×(OR=0.194) anti-discriminative**；complete TP/FP/FN ~280k：TP 3.95%>FP 1.07% 但 sensitivity ~4%、COLO829 TP≈FP。外部(L3)**確認負結果**：Kapoor L3.2/L1.3(LOSO 循環是教科書 ML 陷阱)、Soneson(in-dist→0/external→隨機)。**無 production caller 用甲基當 standalone FP 判別器**(ClairS-TO/VarNet/AIVariant/Krishnamachari 用 alignment/PoN/purity/CN)。AIVariant 唯一部分反例(DL image 內 epigenetic 特徵低 purity 有助)→誠實 **DIFFERENT-VIEW 非乾淨 refute**——但 methyl-as-filter 正確 locked NEGATIVE。**鏈在此終止：subclone 不流向可用 TP/FP filter** |
| **U1** | 底層：ASM/cis 生物真實(Q7)| 🟢 **REAL（存在+cis+CN-controlled）/ 🟡 部分（方向）/ 🔴 BLOCKED（recurrence、可用性）**| 存在+cis-機制確認。內部：BRCA2 HP-axis Δβ=−0.122、within-somatic −0.023(p=0.022, **hypo**)；cis-test d_cis=−0.142 vs drift≈−0.022；4 perm-survivor；6/6 excess(3 癌種)；**ASM×CN ρ=−0.055(非 CN-driven)**。外部：Martin-Trujillo(CN 解釋 82-92% imprinted allelic 甲基)=**我們 HP-axis 擊敗的確切 confound=最強佐證**；Do2020/Schalkwyk CONFIRM cis-driven；Sakamoto22 CONFIRM **private/非復發**。**兩誠實 flag**：(a) BRCA2 **hypo** 違 canonical TSG **hyper**(Herman/VHL)——我們軸是 allelic 結構 Δ 非 promoter-mean silencing，**勿當 TSG-silencing 證據**；(b) **recurrence BLOCKED**——O'Neill/Lin24(189) 偵復發 RET/CDKN2A，我們單 pipeline 6 樣本偵不到，「private 0/38」誠實但 underpowered。**ASM 真實 ≠ ASM 可用(L5 已死)** |

**夢的狀態總結**：
- 🔴 **DEAD**：甲基-as-TP/FP-filter(L5, concluded NEGATIVE)；cross-region 二次打擊**順序**(L4, NEGATIVE)；「單樣本群體統計定位大量乾淨位點」(L1)；BRCA2-as-TSG-silencing(Q7 方向)。
- 🟢 **LIVE（最強）**：tumor/LOH 甲基輔助 phasing/tagging(L2)——真白地，germline 已證可行、無人在癌症測；LOH-constrained phasing 主軸(Grade A)；ASM 存在+cis+CN-controlled(U1)。
- 🟡 **OPEN/未證**：ISM 方法*優於* pycoMeth/cvlr(L3，無 head-to-head)；subclone *characterization*(非 reconstruction)(L4)；5mC/5hmC 分軌新穎性(L3/U1，文獻 silent)。

## 2. 能力 triage

| 能力 | 價值 | 狀態 | gating 證據 |
|---|---|---|---|
| tumor/LOH 甲基輔助 phasing/tagging | **high** | **possible-未證（LIVE 白地）** | LOH-constrained Grade A(p=0.0078; NG=2 6/6)。機制 germline 已證(MethPhaser/HapBridge) 但**無工具測 tumor/LOH**。缺：switch-error/N50 head-to-head、LOH 內 phase-truth |
| normal-anchored somatic cis-test | **high** | **achievable-now（已建，單樣本）** | BRCA2 d_cis=−0.142 vs drift≈−0.022；4 perm-survivor。**零外部競爭** 合 somatic cis-test+read-level。封頂單樣本 ⭐3 |
| 5mC/5hmC 分軌做 ASM | **med** | **achievable-now（新軸但訊號弱）** | 文獻 SILENT(同儕皆 5mC-only)。內部 ASM 5mC-driven、5hmC marginal 0.03-0.07。新但 5hmC 近 noise。單樣本 |
| ASM 存在 characterization(CN-controlled/cis/跨癌 replicate)| **med-high** | **achievable-now（單 pipeline 封頂）** | 6/6 excess、3 癌種；ASM×CN ρ=−0.055(CONFIRM vs Martin-Trujillo)。封頂 ⭐3：單 pipeline。不能宣稱 recurrence(vs Lin24/O'Neill cohort)|
| subclonal-structure characterization(NGroups/AF→NGroups)| **med** | **difficult（「reconstruction」誇大）** | AF→NGroups +0.705~0.787(7/7) POSITIVE。但 HPFineNGroups ⭐4→⭐3。缺三金標準(phylogeny/CN-deconv/單細胞)。「characterize」OK，「reconstruct」不行 |
| ISM 方法優於 incumbents(vs pycoMeth/cvlr)| **med** | **possible-未證** | niche-組合未佔，但 read-read-距離*單獨*非新。**無 head-to-head**。只能稱「互補軸」非「優越工具」 |
| cross-region 二次打擊**順序** | **low**(paused)| **refuted（NEGATIVE, parked）** | 內部 NEGATIVE；無外部 read-level 二次打擊順序 claim。正確 parked |
| 甲基-as-TP/FP filter | **low**(若成立會 high)| **refuted（concluded NEGATIVE, locked）** | LOSO 100% circularity；AUC 0.5049；strong-ASM anti-discriminative(FP 5×)。外部 CONFIRM 負(Kapoor/Soneson)。**locked — 勿 reopen** |
| cohort-recurrence aDMR 偵測(RET/CDKN2A-style)| **high**(臨床 payoff)| **refuted-for-us（underpowered，非 refuted-in-field）** | Lin24/O'Neill 在 n=189 偵復發 aDMR。我們 6 樣本單 pipeline→0/38。**無法分「無復發」vs「underpowered」**。最大 scope gap |

## 3. 可信度 × 衝突矩陣

| 關鍵內部 claim | tier | 外部狀態 | 衝突度 | 註 |
|---|---|---|---|---|
| 長讀使能 read-level HP-linked 甲基 | L1-L2 | **external-confirmed** | 🟢 none | NanoMethPhase/POG-Lin24/de Abreu。無爭議使能前提 |
| tumor 5mC 混亂；單樣本群體統計抓不到大量乾淨位點 | L1-L2(1 株)| **external-different-view（部分-contradict）** | 🔴 direction(部分)| Landau/Landan/Jenkinson CONFIRM 失序；**Lin24 REFUTE「找不到很多」at cohort(~23.6K aDMR/腫瘤)**。我們限制=power 非生物 |
| O11 epipolymorphism AUC 0.845→0.530 coverage 校正 | L1 | **external-confirmed** | 🟢 none | Scherer20：epipolymorphism/entropy 需 >10×、coverage-sensitive。真 win |
| strong-ASM 稀有 0.34% 無 hypo/hyper 偏 | L2(單樣本)| **internal-only-novel** | 🟡 unaddressed(caliber)| 無外部反駁；確切率新+單樣本。Do2020「+5×」是不同分母非矛盾 |
| ASM private(0/38)；跨癌 6/6 excess | L2-L3(單 pipeline ⭐3)| **external-different-view（分裂）** | 🔴 direction-of-recurrence | Sakamoto22 CONFIRM private；**O'Neill/Lin24 證 recurrence at n=189**→我們單 pipeline 偵不到。可調和但 caliber gap 決定性 |
| ASM×CN HP-axis ρ=−0.055（非 CN-driven）| L2(HCC1395 pilot ⭐3)| **external-confirmed** | 🟢 none(互撐)| **最強佐證**：Martin-Trujillo(CN 解釋 82-92% imprinted allelic 甲基)=我們 HP-axis 擊敗的 confound |
| BRCA2 cis-driven(d_cis=−0.142 vs drift≈−0.022)| L2(單樣本 ⭐3)| **external-confirmed（機制）** | 🟢 none | Do2020/Do2017/Schalkwyk：ASM cis-driven。caliber：我們 within-individual。niche 成立 tier-capped |
| BRCA2 HP-axis Δβ hypo 方向 | L2(單樣本 ⭐3)| **external-contradicted（方向）** | 🔴 direction(vs canonical)| **須 flag**：Herman/VHL canonical TSG=hyper；我們 hypo。軸=allelic 結構 Δ(allele+copy) 非 promoter-mean silencing→**勿當 TSG-silencing 證據** |
| methyl-filter LOSO 100% circularity / AUC 0.5049 / anti-discriminative | L2(6 樣本，單 pipeline)| **external-confirmed** | 🟢 none | Kapoor L3.2/L1.3+Soneson 確認循環標誌；無 caller 用甲基當 standalone FP。負結果可信、正確 locked |
| ALLELE-axis 被 germline allelic 甲基 confound；只 HP-axis somatic-controlled | L2 | **external-confirmed** | 🟢 none | 2024 nanopore fine-mapping(ASM genotype-driven，同 HCC1395/COLO829)+Cell Discovery 2017。germline-het null 是對的 guard |
| HPFineNGroups subclone marker ⭐4→⭐3(pipeline-dependent)| L2-L3→降級 | **internal-only-novel** | 🟡 caliber-only | 概念(read-level epi-structure→clones) 精神撐(Landau/Li)；特定 marker 外部未驗+pipeline-dependent。誠實 ⭐3 成立 |
| LOH-constrained phasing Grade A(NG=2 inner 93-99% 6/6)| L2-L3(Grade A, 7 樣本)| **external-confirmed（原則）** | 🟡 caliber(germline vs LOH)| MethPhaser 證甲基*是* haplotype-phasable；Wakhan/HiCancer 證 LOH-as-phasing——**但全非甲基或 germline**。我們 methyl-under-LOH 貢獻 distinct+未被他們測 |
| read-read 距離+cis-test+5mC/5hmC 分軌 = unfilled niche | L2(claim)| **external-different-view** | 🟡 caliber-only | niche 只在窄化形式成立：cis-test+LOH-耦合零競爭；read-read-距離*單獨*非未填(cvlr/Metheor)。5hmC 分軌真 silent |
| cross-region 二次打擊順序 NEGATIVE | L2(paused)| **internal-only（無外部 claim）** | 🟡 unaddressed | 無外部 read-level 二次打擊順序 claim。正確 parked；無衝突 |

## 4. 最可防守論文主張 + reviewer 防守

**最可防守主張**：*「在 ONT 長讀癌症基因組中，allele/haplotype-specific 甲基化(ASM)真實、cis-driven、且設計上獨立於 copy number——但它是結構/連結的 characterization 訊號，不是 somatic variant filter。我們貢獻 (a) normal-anchored somatic cis-test，用 HP1-vs-HP1-1 軸 held-constant CN/ploidy/alignment 來分辨真 cis-driven ASM 與 drift；(b) 把 haplotype-linked 甲基用於 LOH 區(germline-het anchor 消失處)的 phasing/tagging——現有 germline methyl-phaser(MethPhaser/NanoMethPhase/HapBridge)明確留為 future work 的白地。」***

此主張撐過每個內部負結果：只宣稱資料支持的——定位 ISM 為 characterization+tooling(對齊 filter NEGATIVE)、foreground CN-confound 擊敗(唯一鐵證外部佐證 Martin-Trujillo)、押 LOH-phasing 白地(真未填、可行性已證) 而非過度宣稱 subclone reconstruction 或 recurrence。明確排除死夢(filter、二次打擊順序)與 underpowered claim(recurrence、「定位大量」)。

**3 大 reviewer 脆弱點 + 文獻先發**：
1. **「單樣本 HCC1395/單 pipeline ClairS-TO——ASM 率與 cis-candidate 不可復現/不復發。」**（最硬、誠實天花板）→ **先發**：明 concede ⭐3 cap，引 Lin24/O'Neill(POG 189) + Sakamoto22 正確 frame——Sakamoto *private patient-specific* hap-DMR 直接佐證我們 0/38 private，故 private-ness 是 cohort-一致*預期*非失敗；recurrence(RET/CDKN2A) 是 cohort-scale 讀出我們透明缺 power。Frame「現象 replicate(6/6 excess, 3 癌種)；recurrence 需 cohort scale(future work)」。
2. **「read-read-距離方法不新——cvlr clustering、pycoMeth haplotype-aware read-level DMR、Metheor 距離-aware discordance。」**→ **先發**：在 reviewer 之前自己窄化新穎性。引 cvlr(Raineri23)/pycoMeth(Snajder23)/Metheor(Lee23 — 註 LPMD 是 *within-read* CpG-CpG 距離非 between-molecule)，concede read-read-clustering 單獨已被佔。只押**未佔組合**：normal-anchored *somatic* cis-test + LOH/CN 耦合(零競爭) + 5mC/5hmC 分軌(文獻 silent)。承認無 head-to-head→frame「互補軸」非「優越工具」。
3. **「BRCA2 hypo-甲基矛盾 canonical TSG hyper-silencing——你的訊號是真生物嗎？」**→ **先發**：自己 flag 方向並澄清量測。引 Herman & Baylin/VHL(canonical=promoter-mean HYPER-silencing)，澄清我們 HP-axis Δβ 是 **等位結構差(含 allele+copy)非 promoter-mean silencing**——不同軸故非矛盾；再引 Do2020 hypo-dominant cancer ASM(49-76%) 顯 allele-specific *loss* 本身是公認頻繁 cancer-ASM pattern。明說 BRCA2 **不**當 TSG-silencing 證據——只當過 permutation cis-test 的 cis-driven allelic-methylation 例。

**相關**：[00 索引](00_index.md) ｜ 7 角度 [Q1](01_goal_alignment_methods_timeline.md) [Q2](02_finding_methyl_diff_loci.md) [Q3](03_ism_method_comparison.md) [Q4](04_subclone_landscape.md) [Q5](05_methyl_tp_fp_germline_somatic.md) [Q6](06_methyl_assisted_phasing.md) [Q7](07_asm_cis_cancer_impact.md) ｜ [09 參考文獻](09_references_and_provenance.md)
