---
id: ism-kb-11-external-literature-subclone-landscape
name: "Q4 Subclone 重建：文獻地景 vs 我們觀察"
description: "從甲基化重建 tumor subclone 的文獻(phyloepigenetics Hong2010/Landau2014、phylogeny BitPhylogeny/EVOFLUx、CN-deconv CAMDAC、long-read Lee2025)。概念 CONFIRM 但『reconstruction』過度——缺 phylogeny/CN-deconv/單細胞三大金標準；ISM 做 characterization 非 reconstruction。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal vs HPFineNGroups + LOH-constrained phasing + AF→NGroups findings"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, subclone, phyloepigenetics, evoflux, bitphylogeny, camdac, single-cell, long-read, clonal-evolution]
canonical_paths: [11_external_literature/04_subclone_landscape.md]
alias_paths: []
---

# Q4 Subclone 重建：文獻地景 vs 我們觀察

- **一句裁決**：方向 **CONFIRM 非 refute**——read-level 甲基多樣性追蹤 clonal 結構是公認(Landau/Li/Hong)，甲基 ITH mirror 基因組 subclonal 結構(CAMDAC/TRACERx)，長讀 ONT 可重建 subclone 軌跡(Lee25)。但主流 framing 是 **DIFFERENT-VIEW**：金標準是 *disorder/entropy* 或 *probabilistic phylogeny/clock*；ISM 的 read-read 距離+cis-test+5mC/5hmC 分軌是 distinct **未佔 niche，但窗口因 2024-25 論文快速收窄**。「reconstruction」**過度**——我們做 flat clustering，缺 phylogeny/CN-deconv/單細胞三大支柱，**「subclonal-structure characterization」才是可防守動詞**。

## 問題與內部立場

從甲基化偵測/重建 subclone（尤其 read-level/long-read/單細胞）的研究地景；是否確認我們觀察(HPFineNGroups、LOH-constrained phasing、AF→NGroups)，ISM 能否達成 subclone 偵測，缺哪些觀察/方法。

**內部錨點（L1-L2）**：HPFineNGroups subclone marker **⭐4→⭐3**（pipeline-dependent，raw ClairS-TO split 無法復現 89.1%）；LOH-constrained phasing **NG=2 inner 93–99% same-hap (6/6)**；Inter **AF→NGroups +0.705~0.787 (7/7)**；cross-region 二次打擊**順序 NEGATIVE**(paused)。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Hong/Siegmund（推 tumor ancestry）, PLOS One | 2010 | PMID 20711251 / 10.1371/journal.pone.0012002 | array/locus bisulfite；12 CRC TO 多區 | passenger 甲基當分子鐘；單 clonal expansion 解釋 11/12 CRC。**phyloepigenetics 奠基** | CONTEXT（AF→NGroups-as-clock 概念父）| region-level 分子鐘模擬 CRC array；非 read-level long-read；germline-uncontrolled | high |
| Landau（CLL locally disordered）, Cancer Cell | 2014 | PMID 25490447 / 10.1016/j.ccell.2014.10.012 | WGBS+RRBS；104 CLL+26 normal+14 縱向 | PDR；67.6% ITH 變異來自 discordant reads；PDR 隨遺傳 clonal 演化升(p=0.037)；高 PDR→差存活(HR 2.81)| DIFFERENT-VIEW（量 DISORDER；我們 NGroups 量 structure）| 短讀無 haplotype/LOH；cohort vs 單株；disorder≠structure | high |
| Li, methclone, Genome Biol | 2014 | PMID 25260792 / 10.1186/s13059-014-0472-5 | 短讀 bisulfite；白血病縱向 | 組合熵差「epiallele shift/million loci」clonal dynamics | DIFFERENT-VIEW（NGroups benchmark target）| 短讀無 haplotype linkage/5mC-5hmC 分軌 | high |
| Chen/Ashoor, epihet, Sci Rep | 2021 | DOI 10.1038/s41598-020-79627-x / PMC7801679 | 短讀 bisulfite；AML | PDR+epipolymorphism+Shannon entropy per locus；網路分析 | DIFFERENT-VIEW（可比的 read-level disorder 度量）| 短讀；entropy/PDR 族非 haplotype-structure；撐我們 O11 baseline | high |
| Yuan, **BitPhylogeny**, Genome Biol | 2015 | PMID 25786108 / 10.1186/s13059-015-0592-6 | 甲基+SNV 方法 | full-Bayesian 聯合推 clone 數/組成/樹(TSSB)；甲基當演化 marker | CONTEXT（canonical phylogeny 方法，**我們未實作**）| 概率 phylogeny；我們 flat clustering 無樹/clone-count posterior | high |
| Gabbutt, **EVOFLUx**, Nature | 2025 | PMID 40931062 / 10.1038/s41586-025-09374-4 | array/短讀 bisulfite；1,976 淋巴瘤 bulk TO | fluctuating-CpG 鐘從 **bulk 推**生長率/惡性年齡/epimutation rate；subclonal selection 罕見偵測 | DIFFERENT-VIEW/CONTEXT（SOTA bulk subclonal 推論，**我們缺**）| bulk array/短讀，n=1,976 vs 單株 read-level；定量演化模型我們缺 | high |
| Larose Cadieux/TRACERx, **CAMDAC**, Nat Genet | 2025 | PMC12425823 / 10.1038/s41588-025-02307-x | RRBS 短讀；217 區/59 NSCLC 多區 paired | CN+purity-aware deconv：純瘤甲基=bulk−normal 加權(CN×purity)；甲基 ITH mirror 基因組 clonal 結構 | **CONFIRM（甲基 ITH 追蹤基因組 subclonal）+ CONTEXT（CN-confound 控制法，我們缺）**| RRBS cohort 多區 paired，顯式 CN/purity 校正；我們 HP-axis held-constant CN 但無正式 deconv | high |
| Lee（長讀黑色素瘤 subclone）, bioRxiv | 2025 | DOI 10.1101/2025.08.28.672865 / PMC12424993 | **長讀 ONT(基因組+5mC)**；23 單細胞衍生 subclone(B2905 鼠) | ONT 整合 SNV+SV+CNV+差異甲基；平行演化；lineage-specific 甲基追蹤侵襲 phenotype | **CONFIRM（長讀甲基 subclone 重建可行——最近 ISM setup）**| 長讀 ONT 株，但用單細胞**衍生** subclone 當 ground-truth(clonal 培養)，非 in-silico read-level deconv；preprint | preprint |
| Fu, MethPhaser, Nat Commun | 2024 | PMID 38909018 / 10.1038/s41467-024-49588-0 | 長讀 ONT；germline 株+4 血液 | 甲基延伸 SNV phasing；N50 +78–151% @ 83.4–98.7%(germline SNV anchor 耗盡處)| CONFIRM（LOH-constrained phasing + methyl-rescue 的方法表親）| germline(HET-SNV-anchored)；**不涵蓋 somatic/LOH** = 我們的 cancer gap | high |

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| HPFineNGroups subclone marker ⭐4→⭐3 pipeline-dependent | SILENT on 此 construct；Landau/Li/epihet 顯示 read-level epi-structure 追蹤 clones 但無人驗 haplotype-tag-derived NGroups marker | 🟡 caliber-only | 概念(read-level epi-structure→clones) 精神被撐；特定 marker 外部未驗+pipeline-dependent。誠實 ⭐3 成立 |
| LOH-constrained phasing NG=2 inner 93–99% (6/6)；Grade A | MethPhaser CONFIRM haplotype-linked 甲基廣泛可 long-read phase；但僅 germline/HET | 🟡 caliber(germline vs LOH)| CONFIRMED in principle(甲基*是* haplotype-structured)；我們新穎=在 LOH(germline anchor 消失)做，MethPhaser 不蓋 |
| Inter AF→NGroups +0.705~0.787 (7/7)| Hong/Siegmund+EVOFLUx CONFIRM 甲基多樣性如鐘隨 clonal/AF 結構 scale | 🟡 caliber-only | CONFIRMED 概念(AF/clonality↔甲基多樣性是公認分子鐘關係)；機制仍需 phasing-vs-methylation 解耦 |
| cross-region 二次打擊順序 NEGATIVE (paused)| TRACERx-NSCLC 顯示甲基-基因組 CO-EVOLUTION 但**未**證 read-resolvable 二次打擊順序 | 🟡 unaddressed | 我們 NEGATIVE 與外部無正面證據一致；無衝突，正確 parked |
| read-read 距離+cis-test+5mC/5hmC 分軌 = unfilled niche | 無外部把三者合一；最近 Lee25 ONT 用單細胞衍生 clone 非 read-level deconv | 🟢 none（真空）| niche **成立但**鄰地快速填(Lee25/EVOFLUx/CAMDAC 皆 2024-25)，窗口收窄 |

## 缺口（我們缺什麼）

1. **Phylogeny/樹重建**（BitPhylogeny/EVOFLUx）：我們 flat 階層 clustering→NGroups，**無 clone-count posterior、無樹、無 branch-length/分子鐘模型**。EVOFLUx/BitPhylogeny 是「subclone reconstruction」應 benchmark 的 canonical 方法——目前我們 under-claim「reconstruction」(我們 cluster，不 phylogenize)。
2. **CN-aware purity deconv**（CAMDAC）：正是 ASM/LOH 甲基 claim 的 confound 控制法(純瘤甲基=bulk−normal 加權 CN×purity)。我們 HP-axis *held-constant* CN(ρ=−0.055，good) 但無正式 deconv；BRCA2/ZAR1L copy-component d_copy=−0.11，CAMDAC-style 校正是缺的金標準。
3. **單細胞 ground-truth**（Lee25/scWGBS phylogeny）：Lee25 從單細胞衍生 subclone 當 ground-truth 再長讀 profile；我們**無單細胞 anchor**——NGroups 是 bulk read clustering 的 pseudo-subclone 無 orthogonal 驗證。**這是任何 subclone-reconstruction claim 的最大 credibility gap。**
4. **vs 既有 disorder 度量 head-to-head**（PDR/epipolymorphism/epiallele-shift）：定位 read-read 距離=STRUCTURE vs 他們 DISORDER，但未在同 loci 跑 NGroups vs PDR/epipolymorphism/EPM。差異是 assert 非 demonstrate。
5. **Cohort scale + cross-pipeline**：每個外部 CONFIRM 皆 cohort-scale(Landau 104、EVOFLUx 1,976、TRACERx 59)；我們單株(HCC1395)+單 pipeline(ClairS-TO，⭐3 封頂)。單瘤多區取樣(Hong2010/TRACERx 的設計)我們在甲基-subclone 軸沒有。
6. **LOH 下 somatic phasing**（germline anchor 失去）：MethPhaser 解 germline methyl-phasing；**無外部解 LOH 下 methyl-assisted phasing**——真未填，且與我們 LIVE 主軸(methyl-rescue + LOH-constrained phasing) 重疊 = 最可防守 niche。

## 裁決

方向**廣泛 CONFIRM 非 refute**：領域堅實確立 (a) read-level 甲基多樣性追蹤 clonal/subclonal 結構(Landau/Li/epihet/Hong)、(b) 甲基 ITH mirror 基因組 subclonal 架構(CAMDAC/TRACERx)、(c) 長讀 ONT 甲基可重建 divergent subclone 軌跡(Lee25)；MethPhaser 獨立確認 haplotype-linked 甲基廣泛 long-read-phasable，直接撐 LOH-constrained-phasing 前提。但主流 framing 是 **DIFFERENT-VIEW**：金標準是 *disorder/entropy*(PDR/epipolymorphism/epiallele-shift) 與 *probabilistic phylogeny/clock*(BitPhylogeny/EVOFLUx)，ISM 提供 *haplotype-linked read-read 距離 + normal-anchored cis-test + 5mC/5hmC 分軌* = distinct 且(本候選集) **未佔 niche，故 niche 成立**。誠實 caveat 決定性：(i) 我們無 phylogeny、無 CN-aware deconv、無單細胞 ground-truth = 領域 legitimize「subclone reconstruction」的三支柱(BitPhylogeny/EVOFLUx、CAMDAC、Lee25/scWGBS) 都缺；(ii) 每個外部佐證 cohort-scale，我們單株/單 pipeline ⭐3，「reconstruction」目前**誇大** flat NGroups clustering 所為——「**subclonal-structure characterization**」才是可防守動詞；(iii) 鄰地快速填(三篇最強是 2024-25)，niche 窗口真實但收窄。HPFineNGroups 維持 pipeline-dependent ⭐3，二次打擊順序 NEGATIVE 與無任何外部 read-level 二次打擊順序 claim 一致——皆正確 parked。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ HPFineNGroups KB：[../07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md) ｜ Q6 phasing：[06](06_methyl_assisted_phasing.md)
