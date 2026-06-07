---
id: ism-kb-11-external-literature-methyl-assisted-phasing
name: "Q6 癌症甲基輔助 phasing/tagging（somatic、LOH 脈絡）"
description: "甲基能否在癌症(somatic/LOH，germline-het anchor 消失)輔助 phasing/read-tagging。germline 機制已證(MethPhaser/NanoMethPhase/HapBridge)且明確留 cancer 為 future work；LOH-as-phasing 被 Wakhan/HiCancer/Sakamoto 非甲基背書；tag→caller 被 LongPhase-S 證(ClairS +4.5%)。tumor/LOH 甲基-phasing 沒人做過 = 最強 LIVE 白地。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal vs methyl-rescue pilot + LOH-constrained phasing Grade A + HKU tag-to-caller"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, phasing, haplotag, methphaser, hapbridge, loh, somatic, tag-to-caller, longphase-s, white-space]
canonical_paths: [11_external_literature/06_methyl_assisted_phasing.md]
alias_paths: []
---

# Q6 癌症甲基輔助 phasing/tagging（somatic、LOH 脈絡）

- **一句裁決**：方向 **LIVE 真白地、方法可行**——但文獻**尚未提供** methyl-assisted phasing 在 tumor/LOH 可行的外部佐證，*那正是白地*。三個 germline 工具(MethPhaser/NanoMethPhase/HapBridge) 證 methyl-as-phasing-linkage 及我們 pilot 的「甲基救 unphase read」確切機制，**但全 germline 株、明確把 cancer/LOH 留 future work**。LOH-as-phasing 前提被 ≥3 癌症論文(Wakhan/HiCancer/Sakamoto) 獨立背書但走 CNA/SV/Hi-C 非甲基。tag→caller 被 LongPhase-S 強驗(ClairS +4.5% SNV)，但其 tagging germline-anchored 不解我們 germline-het-absent gap。**niche 成立，升 ⭐4 需多樣本 tumor + LOH 內 phase-truth + vs germline methyl-phaser head-to-head。**

## 問題與內部立場

甲基(或其他 read-level 資訊)能否在癌症(somatic/LOH 抹去 germline-het anchor 處)輔助 haplotype phasing 與 read-tagging；MethPhaser-style methyl-phasing 在 tumor 是否可行；有無 somatic/cancer 成功案例。

**內部錨點（L1-L2）**：methyl-rescue-unphased-read pilot；**LOH-constrained phasing Grade A**(n=7, W=28, p=0.0078；NG=2 inner 93–99% 6/6) = LIVE 主軸；HKU ClairS **tag→caller**。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Fu, **MethPhaser**, Nat Commun | 2024 | 10.1038/s41467-024-49588-0（PMID 38909018）| 長 ONT；**germline 株+血液，無 tumor** | 甲基延伸 SNV phasing：N50 ×1.6–2.5(HG002)、僅 ×1.07–1.17 血液(短 read)；準度 83.4–98.7%；**cancer 是 future work** | **CONTEXT**（方法先例非癌證據）| 平台 match 但 germline-only、非 somatic-controlled、非 LOH | high |
| Akbari, **NanoMethPhase**, Genome Biol | 2021 | 10.1186/s13059-021-02283-5（PMID 33618748）| 長 ONT；**germline 株**(含 Colo829**BL**=normal)；無 tumor | phase 5mC @~10×；92.3% reads 與 trio congruent；偵測 60% ICR；討論明點「**對 LOH 癌症的效用**」但未測 | **CONTEXT** | 甲基-as-phasing 奠基；germline-only；Colo829BL=normal 非 tumor | high |
| **HapBridge**, bioRxiv | 2025 | 10.1101/2025.11.07.687303（preprint）| 長 ONT R9/R10+PacBio；**germline，無 tumor** | 迭代用 haplotype-specific 甲基 tag unphased reads 再發現新 ASM；switch error −3.07–18.72%、N50 +5.84–68.61% **vs MethPhaser**。**正是我們 pilot 的「甲基 tag unphase read」機制** | **CONTEXT（最近機制 analog）**| 機制 match 我們 pilot **但 germline-only，無 somatic/LOH**；region+迭代-site 非 read-read 距離+cis | preprint |
| Sakamoto, **肺癌 phasing**, Nat Commun | 2022 | PMID 35710642 / 10.1038/s41467-022-31133-6 | 長 ONT+Illumina；**tumor-normal paired，n=20 NSCLC** | 真 somatic phasing：834kb N50、>99% concordance；~43% 點突變指派 haplotype；190 haplotype-biased 區；**經 phasing 偵測 LOH(PTEN del)**。甲基由 Nanopolish 在 phasing **後**call——**非當 phasing 訊號** | **CONTEXT/DIFFERENT-VIEW**（SNV-only somatic phasing 有效；甲基不用來 phase）| cohort+paired+somatic-controlled(比我們強)；但甲基是下游讀出非 phasing input——與我們 methyl-assist 假說相反 | high |
| O'Neill/POG, Cell Genomics | 2024 | PMID 39406235 / 10.1016/j.xgen.2024.100674 / PMC11605692 | 長 ONT；189 腫瘤/43 matched normal | long-range WhatsHap phasing→haplotag→ASM(含 somatic-bearing alleles)；復發 RET/CDKN2A aDMR；**79% aDMR 落 het CNV+LOH**；allelic promoter 甲基多在 **trans** 與主表現等位 | **CONFIRM(現象)/CONTEXT** | 最大真癌 ASM cohort；somatic-aware；**但 ASM 經 SNV/CNV phasing 偵測，甲基不用來 PHASE**；region/aDMR-level | high |
| **Wakhan**, medRxiv | 2025 | DOI 10.64898/2025.12.11.25342098（preprint）| 長 ONT/PacBio；**tumor**(株 vs HATCHet)| 用 **tumor allelic imbalance/CNA haplotype 差異** + 一 somatic SV 只影響一 haplotype → haplotype-specific CNA，全染色體 scale phasing。直接 **LOH/imbalance-driven phasing**。無甲基 | **CONFIRM（我們 LOH-axis 機制，非甲基）**| 同我們 LOH-constrained(NG=2 inner 93–99%) 但走 CNA/SV 非 read-level 甲基；tumor 株 | preprint |
| Mao, **HiCancer**, Sci Rep | 2021 | PMID 33758310 / 10.1038/s41598-021-86104-6 ⚠venue 見 09 | 短(Hi-C)；癌株(K562/KBM-7)| phase ~100% SNP(vs HAPCUT2 70–80%)；AER 0.54–1.53%；**顯式處理 LOH**(偵測 low-het 區當單等位定序)；LOH sensitivity 94.4–100%。無甲基 | **CONFIRM（LOH 處理原則）/DIFFERENT-VIEW（短讀）**| 非甲基對照；LOH-as-single-allele 洞見先於我們；Hi-C/短讀 germline-SNP | high |
| Ho, **LongPhase-S**, bioRxiv | 2025 | 10.1101/2025.11.20.689492（preprint v1）| 長 ONT；**tumor-normal paired，6 株/8 資料** | somatic haplotyping：每 somatic read anchor 到 germline lineage(HP1-1/HP2-1)、解 germline-vs-somatic、估 tumor DNA fraction、recalib caller(ClairS SNV +4.5%/indel +7.1%；DeepSomatic +1.2%/+0.5%)。tag→caller 即 HKU 路徑。**甲基不用來 tag**(germline-anchor) | **CONFIRM（tag→caller 價值）/DIFFERENT-VIEW（無甲基）**| 直接 HKU 工具；somatic-controlled 多株(比單樣本強)；但 tagging germline-anchored，不解 germline-het-absent/LOH gap | preprint |

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| **Q6 核心**：甲基能在癌症(somatic/LOH，anchor 失去)輔助 phasing/tagging — methyl-rescue pilot | MethPhaser/NanoMethPhase/HapBridge 證機制**在 germline 有效**且明點 cancer/LOH 為 *future/未測*；**無論文展示 tumor/LOH 甲基-phasing** | 🟡 **unaddressed（白地）** | 真**開放 niche**——文獻立 germline 可行性、motivate cancer，但 tumor/LOH extension 未填。pilot 成立為 novel |
| **LOH-constrained phasing Grade A**(n=7, W=28, p=0.0078；NG=2 inner 93–99% 6/6)| Wakhan+HiCancer 獨立用 LOH/imbalance phase(非甲基)；Sakamoto 經 phasing 偵 LOH | 🟢 none/caliber | **CONFIRMED in principle。** LOH-as-phasing 被 ≥3 癌症論文背書——強化外部效度。但他們用 CNA/SV/Hi-C **非 read-level 甲基**，我們甲基貢獻 distinct(且他們未測)|
| **HKU ClairS tag→caller**(餵 phasing/haplotag tag 給 caller)| LongPhase-S 直接驗：somatic-read tagging→ClairS +4.5% SNV F1；ClairS 已整合 LongPhase phasing | 🟢 none | **CONFIRMED。** tag→caller 是公認量化 gain。我們角色是 LongPhase-S **沒提供**的 *methyl-assisted* tag 層(它 germline-anchored)|
| **ASM 真實但非方向/非判別/coverage-modulated**(單樣本 HCC1395 ⭐3)| POG cohort：aDMR 真實，79% LOH/CNV，多 **in-trans**；O11 0.845→0.530(內部)| 🟡 caliber(cohort vs 單樣本)| **CONFIRMED+extended。** POG(189) 確認 ASM 真實+結構耦合 LOH/CNV——比我們高 caliber 對齊 coverage-modulation+LOH-耦合 caveat。我們非判別與「ASM≠乾淨 phasing label」一致 |
| **甲基能解 subclone/read-level 判別**(HPFineNGroups 4★→3★；methyl→TP/FP NEGATIVE)| 無論文宣稱甲基**判別** subclone/variant；都當 phased *讀出* 或 *bridging linkage* 非 classifier | 🟢 none | **一致。** 文獻從不把甲基當 variant/subclone *判別器*；撐我們 filter NEGATIVE。甲基驗證角色是 *linkage/bridging* 正是 phasing-assist 軸(非 filter 軸)|
| **read-read 距離+cis-test+5mC/5hmC 分軌 = unfilled niche** | HapBridge(迭代 ASM-site tag)/MethPhaser(block bridge)/POG(aDMR cis-trans)/epihet(disorder) — 無人合**read-read 距離+normal-anchored cis-test+5mC/5hmC 分軌 in tumor** | 🟡 unaddressed | **NICHE 成立。** cancer/LOH 脈絡無外部佔此組合 |

## 缺口（我們缺什麼）

1. **cis vs trans 等位定向 framing**：POG headline = allelic promoter 甲基多 **in trans** 與主表現等位、79% aDMR 落 LOH/CNV。我們 cis-test 機制有，但未把 HCC1395 ASM loci 按 **cis/trans 相對表現或保留 LOH 等位** report——cohort 已驗的高 impact 讀出。
2. **CNA/allelic-imbalance 當*互補* phasing 訊號（非 confound 扣除）**：Wakhan/HiCancer *用* allelic imbalance/CNA phase 跨 LOH；我們目前把 CN 當 confound 扣(ASM×CN ρ=−0.055)。缺建設性用法——甲基*與* CNA imbalance 合用 phase germline-het-absent 區(NG=2 是種子，Wakhan 是方法模板)。
3. **methyl-rescue vs MethPhaser/HapBridge head-to-head**：HapBridge 報具體 delta vs MethPhaser(switch error/N50)；我們 pilot 無可比 switch-error/N50/read-tag-recovery 度量，無同 tumor 資料對比——領域標準 yardstick 我們缺。
4. **LOH 內真 somatic-controlled methylation-phasing benchmark**：無外部提供 LOH(germline-het absent) 內甲基-phasing ground-truth 準度。這正是白地——但也意味我們**無外部 truth set** 驗 rescue(只 SEQC2 SNV truth，非 LOH 內 phase truth)。升 tier 需 orthogonal phase-truth(trio/Hi-C/strand-seq)。
5. **patient-blood 退化訊號（MethPhaser）**：MethPhaser gain 從 ×1.6–2.5(HG002) 崩到 ×1.07–1.17 血液(短 read)——甲基-phasing uplift **read-length-sensitive**，我們未對 tumor BAM characterize。

## 裁決

方向**確認為真未填 niche、方法可行**，但文獻**尚未提供** methyl-assisted phasing 在 tumor/LOH 可行的外部佐證——*那正是白地*。三 germline 工具(MethPhaser peer-reviewed、NanoMethPhase peer-reviewed、HapBridge preprint) 證 methyl-as-phasing-linkage 及我們 pilot 確切的「甲基 tag unphase read」機制，**但全 germline 株、明確把 cancer/LOH 留 future**——無 somatic-controlled、無 LOH test。LOH-as-phasing 前提(撐我們 Grade-A 主軸)被 ≥3 癌症論文(Wakhan/HiCancer/Sakamoto)獨立背書，但走 CNA/SV/Hi-C 非 read-level 甲基，我們甲基貢獻 distinct 且他們未測。tag→caller 路徑被 LongPhase-S **強驗**(ClairS +4.5% SNV F1)——但 LongPhase-S tagging germline-anchored，**不**解我們 methyl-rescue 鎖定的 germline-het-absent gap，留出 methyl-assisted tag 層空間。POG cohort(189)**在更高 caliber 確認我們 ASM caveat**(ASM 真實、79% LOH/CNV-耦合、多 in-trans)，既驗現象又強化甲基是 *結構/linkage* 訊號非 *判別器*——與我們 filter NEGATIVE 一致。**誠實限制**：所有 ASM/phasing-assist 證據單樣本 HCC1395 ⭐3+單 pipeline；無 vs MethPhaser/HapBridge 的 switch-error/N50 head-to-head；無 LOH 內 phase-truth(只 SNV truth)；未按 cis/trans framing。niche 成立，升 ⭐3 以上需多樣本 tumor + LOH 內 orthogonal phase-truth + vs 既有 germline methyl-phaser 直接 benchmark。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ Q4 subclone：[04](04_subclone_landscape.md) ｜ Q7 ASM：[07](07_asm_cis_cancer_impact.md)
