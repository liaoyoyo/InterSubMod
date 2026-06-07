---
id: ism-kb-11-external-literature-methyl-tp-fp-germline-somatic
name: "Q5 甲基區分 TP/FP 或 germline/somatic：有效嗎"
description: "用甲基/表觀特徵篩 somatic variant TP/FP 或分 germline/somatic。我們 LOSO 100% circularity / AUC 0.5049 NEGATIVE 被 ML 文獻權威背書(Kapoor leakage taxonomy / Soneson confounded-CV)；無 production caller 用甲基當 standalone FP 訊號(ClairS-TO/VarNet/AIVariant 用 alignment/PoN/purity/CN)；ALLELE-axis germline confound 被 nanopore fine-mapping 確認。"
status: active
last_verified: 2026-06-05
content_nature: paper-derived
doc_type: explanation
verified_scope: "papers verified via WebFetch in wf_37b2cc97-663 (2026-06-05); internal vs LOSO + |dBeta| AUC + strong-ASM-in-FP findings"
decision_refs: []
related_ids:
  - ism-kb-11-external-literature-cross-question-correlation
  - ism-kb-11-external-literature-references-provenance
tags: [literature, tp-fp, somatic, germline, variant-filter, circularity, data-leakage, clairs-to, negative-result]
canonical_paths: [11_external_literature/05_methyl_tp_fp_germline_somatic.md]
alias_paths: []
---

# Q5 甲基區分 TP/FP 或 germline/somatic：有效嗎

- **一句裁決**：方向 **DEAD 但被背書、非 refute**。(a) 我們撞的循環性是 *公認、領域級* ML 陷阱——Kapoor L3.2(樣本非獨立)/L1.3(select-then-train) + Soneson(in-dist→0-error/external→隨機)，LOSO 是教科書解；(b) **甲基當 standalone somatic TP/FP 判別器在文獻幾乎從未被採用**——每個 production caller(ClairS-TO/VarNet/AIVariant/Krishnamachari) 用 alignment/PoN/purity/CN/haplotype-consistency，*不用* CpG 甲基；AIVariant 是唯一部分反例(DL image 內未指明 epigenetic 特徵在低 purity 有幫助)→ 誠實 framing 是 **DIFFERENT-VIEW 非 REFUTE**；(c) ALLELE-axis germline confound 被 2024 nanopore fine-mapping(同 HCC1395/COLO829) 獨立確認。**負結果方法學可信、正確 locked。**

## 問題與內部立場

有沒有人用甲基/表觀特徵篩 somatic TP/FP 或分 germline/somatic；結果正/負；我們撞的 cross-validation 循環性/data-leakage 是否公認統計陷阱。

**內部錨點（L1-L2）**：methyl-filter LOSO **100% sample-level circularity**(in-dist ΔF1 +0.0224 vs held-out −0.0001；HCC1954 transfer −0.377)；連續 **|Δβ| AUC=0.5049**(perm p=0.496 中性)；**strong-ASM 在 FP 富集 5×**(OR=0.194, p=1.8e-28) anti-discriminative；complete TP/FP/FN 6 樣本 ~280k loci：TP sig 3.95% > FP 1.07% 但 sensitivity ~4%、COLO829 TP≈FP；ALLELE-axis 被 germline allelic 甲基 confound(TP 11.1% < het null 15.2%)，只 HP-axis somatic-controlled(TP 7.2% vs null 4.1%, OR=1.79 modest)。

## 外部文獻（WebFetch 實證）

| 文獻 | 年 | PMID/DOI | 平台·設計 | 關鍵發現 | 關係 | 口徑 | Cred |
|---|---|---|---|---|---|---|---|
| Kapoor & Narayanan, Patterns | 2023 | PMID 37720327 / 10.1016/j.patter.2023.100804 | 方法 meta(294 論文/17 領域)| 8-type leakage taxonomy；**L3.2 訓練-測試樣本非獨立** + **L1.3 在 train+test 選特徵** 直接命名我們 LOSO 循環+select-then-train；修正後 ML「優越」消失 | **CONFIRM**（我們陷阱是公認廣泛 ML 失敗）| 領域無關 taxonomy；我們 LOSO=他們修正協議 | high |
| Soneson, PLOS ONE | 2014 | 10.1371/journal.pone.0100335（PMID UNVERIFIED）| 微陣列表現；PETACC-3 模擬 | batch-class confounding 下 CV **null 資料 0% error** vs external **~50%**；ComBat 救不回 | **CONFIRM**（正是我們 in-dist +0.0224 vs held-out −0.0001）| gene-expr 非 variant-calling；機制相同 | high |
| Zheng, **ClairS-TO**, Nat Commun | 2025 | 10.1038/s41467-025-64547-z（PMID UNVERIFIED）| ONT 長讀；tumor-only | FP 控制靠 **(1)9 hard-filter (2)PoN (3)purity/CN germline-vs-somatic separator**。**無甲基**。benchmark HCC1395+COLO829 | **CONTEXT/DIFFERENT-VIEW**（我們 methyl filter 所處 TO FP-control 空間——且**無甲基就成功**）| 同平台/樣本；用 purity/CN 統計分離器，正是我們想用甲基填卻失敗的 niche | high |
| Zheng, **ClairS**, bioRxiv | 2023 | 10.1101/2023.08.17.553778（PMID UNVERIFIED）| ONT 長讀；paired | LongPhase-S haplotype-consistency flag 與 germline haplotype 不一致的 somatic call——**haplotype** 非甲基當 FP 訊號 | **DIFFERENT-VIEW**（haplotype-structure FP filter 有效，對齊我們 HP-axis>ALLELE-axis）| paired vs TO；haplotype tag vs 甲基特徵 | preprint |
| Wood, **VarNet**, Nat Commun | 2022 | PMID 35869060 / 10.1038/s41467-022-31765-8 | Illumina 短讀；paired | weakly-supervised DL FP filter(4.6M variant/356 WGS)；特徵=base/mapping quality/strand bias。**無甲基/表觀** | CONTEXT（SOTA learned FP filter 用 alignment 非甲基）| 短讀 paired vs 長讀 TO | high |
| Kim, **AIVariant**, Exp Mol Med | 2023 | PMID 37524869 / 10.1038/s12276-023-01049-2 | Illumina 短讀；paired 9 purity | DL caller 在 read image 編碼「**epigenetic features**」；PR-AUC 0.794@40%/15× vs 對手 0.048–0.524。最近「epigenetic 特徵助 somatic calling」尤其低 purity(epigenetic 未指明，非 read-level CpG 矩陣)| **DIFFERENT-VIEW**（低 purity DL setting 有幫助，對比我們中性結果）| 短讀 paired vs 長讀 TO；bundled DL 特徵 vs 我們孤立 methyl test | high |
| Krishnamachari, Sci Rep | 2023 | 10.1038/s41598-023-37409-1 / PMC10300101（PMID UNVERIFIED）| Illumina cfDNA panel；TO(70 MBC)| ML FP filter：10 tumor-free 個體偵測 **1 variant vs 無模型 521**。特徵=read/alignment/sequence，**無甲基** | CONTEXT（有效 TO ML FP filter，無 matched normal 也無甲基）| cfDNA 短讀 vs 長讀組織 | high |
| Akbari, NanoMethPhase, Genome Biol | 2021 | 10.1186/s13059-021-02283-5（PMID UNVERIFIED auth-walled）| ONT 長讀；germline phasing | phase 5mC 到 haplotype @~10×；偵測 allele/haplotype-specific 甲基 | CONTEXT/CONFIRM-method（驗證 HP-linked 5mC 可行=我們量測基質）| germline-only；無 somatic/TP-FP claim | high |
| Fine-mapping native CpG 甲基(nanopore), bioRxiv/PMC | 2024 | PMID 39386487 / 10.1101/2024.09.27.614715 | ONT 長讀；株含 **HCC1395/COLO829**+GWAS/mQTL | 6,390 haplotype-specific 甲基區綁 **germline variant**；ASM 大多 **genotype/cis-driven**(rs7247241 等 mQTL)| **CONFIRM**（直接撐我們 ALLELE-axis germline-het confound：ASM 由 genotype 決定非 somatic-specific）| germline regulatory ASM 非 somatic TP/FP；**同我們的樣本** | preprint |
| 「Switchable/genetics-influenced ASM」, Cell Discovery | 2017 | PMID 29387450 / 10.1038/celldisc.2017.38 | array/短讀；鼠 reciprocal-cross | 定義兩類 ASM：**genetics-influenced(genotype-dependent)** vs **switchable(隨機 hypo/hyper)**；比例 qualitative | **CONFIRM**（命名我們撞的兩 confound：germline genotype-driven + 隨機 switchable 模仿 somatic 訊號）| 鼠，無 somatic-variant context | high |

> 為誠實/精簡未 fetch 而排除（scout 出現但本角度邊際值低、未獨立驗證、不編造識別碼）：2025 ctDNA RF ensemble、NANOME、epialleleR、DAMEfinder、Woodhouse imprinting atlas、Koelsche 甲基分類器、Cingöz/Onuchic intermediate-methylation。

## 內部 × 外部 cross-reference

| 內部 claim（tier）| 外部關係 | 衝突度 | 裁決 |
|---|---|---|---|
| methyl-filter LOSO 100% circularity；in-dist +0.0224 vs held-out −0.0001；HCC1954 −0.377（L2）| Kapoor L3.2/L1.3；Soneson(CV→0% null vs ~50% external)| 🟢 none | **CONFIRMED 公認 ML 陷阱。** LOSO 是教科書修正；in-dist/held-out gap 是 confounded-CV 標誌。本角度最強外部佐證 |
| \|Δβ\| AUC=0.5049 perm p=0.496 中性（L2，6 樣本 ~280k）| ClairS-TO/VarNet/Krishnamachari 全**無甲基**達 FP 控制 | 🟡 caliber-only | **CONFIRMED-by-absence。** 無外部用 CpG 甲基當 standalone somatic TP/FP 判別器；我們中性結果與「甲基非採用的 FP 訊號」一致 |
| strong-ASM 在 FP 富集 5×(OR=0.194, p=1.8e-28) anti-discriminative（L2，單樣本-主導）| 無直接外部 test；Cingöz/switchable-ASM 暗示隨機 ASM noise-like | 🟡 unaddressed | **NOVEL，外部未addressed。** 無論文報甲基失序在 FP somatic call 富集；可能 switchable/隨機 ASM 聚在易錯 loci，但文獻未測——真 ISM 觀察 |
| ALLELE-axis 被 germline allelic 甲基 confound(TP 11.1% < het null 15.2%)（L2）| 2024 nanopore fine-mapping(ASM genotype-driven)；Cell Discovery 2017(genetics-influenced ASM)| 🟢 none | **CONFIRMED。** 長讀+經典遺傳皆確立 ASM 大多 germline/genotype-determined → ALT-vs-REF 甲基差不可能 somatic-specific。germline-het 負控是正確 guard |
| 只 HP-axis somatic-controlled；TP 7.2% vs null 4.1%, OR=1.79 modest（L2）| ClairS LongPhase-S(haplotype-consistency FP filter 有效)；NanoMethPhase(HP-linked 5mC 可行)| 🟡 caliber-only | **CONFIRMED-方向，幅度 modest。** haplotype-controlled 比較是對的設計；但 OR=1.79 太弱不能當可用 filter——誠實 |
| O11 epipolymorphism AUC 0.845→0.530 coverage 校正（L2）| Soneson(confound inflation)；Kapoor L2(illegitimate/proxy 特徵)| 🟢 none | **CONFIRMED 陷阱。** coverage 當 proxy 正是 L2 illegitimate feature；我們校正方法學健全 |

## 缺口（我們缺什麼）

1. **filter-negative 的 cross-pipeline 外部驗證是單 pipeline**：LOSO 用 6 株全 ClairS-TO call；Kapoor L3.2 警告即使單 pipeline held-out 仍可能殘留非獨立(shared caller artifact)。真獨立 test 需第二 caller(DeepSomatic/Mutect2) 重 call——我們負結果都沒做。
2. **無 imprinting/DMR 扣除 baseline**：宣稱 somatic ASM 前應 mask 已知 germline imprinted DMR；我們 strong-ASM-in-FP 與 ALLELE-axis 未顯式扣 imprinting 控制區。
3. **我們孤立甲基；領域 bundle 它**：AIVariant 正結果來自 DL read-image 內 epigenetic 特徵，非 standalone；我們未測甲基在 alignment/PoN 特徵上的*邊際*提升(領域唯一報它有用的 setting)。我們「甲基單獨→AUC 0.50」真實但不反駁「甲基當小 DL 特徵」。
4. **switchable-ASM 量化缺**：Cell Discovery 2017 命名隨機/switchable ASM 但未量化；我們觀察 strong-ASM-in-FP 富集但未 model 成 switchable-ASM noise——可測機制解釋目前未用。

## 裁決

方向**確認、方法 vindicated 非 refute**。(a) 我們撞的循環/leakage 是*正式命名、領域級* ML 陷阱——Kapoor L3.2(樣本非獨立)/L1.3(select-then-train) + Soneson 模擬顯示確切的 in-dist→0-error/external→隨機 標誌(我們 in-dist +0.0224、held-out −0.0001、transfer −0.377)；LOSO 是教科書解，故 methyl-as-FP-filter 的*負*結論方法學可信。(b) **甲基-as-standalone-somatic-TP/FP-判別器在文獻幾乎從未採用**：每個我們驗證的 production caller(ClairS-TO/ClairS/VarNet/AIVariant/Krishnamachari) 用 alignment/PoN/purity/CN/haplotype-consistency 控 FP，*非* CpG 甲基。AIVariant 是唯一部分反例(未指明 epigenetic 特徵在 DL image 內、低 purity 有幫助)，故誠實 framing 是 **DIFFERENT-VIEW 非 REFUTE**——甲基可能加邊際 DL lift，但無人用它當單獨判別器，與我們 AUC≈0.50 一致。(c) germline-het ALLELE-axis confound 被 2024 nanopore fine-mapping(ASM genotype/cis-driven，*在我們自己的* HCC1395/COLO829) + 2017 genetics-influenced-ASM 範式**獨立確認**。**ISM niche 在本角度成立**：read-read 距離+normal-anchored cis-test+5mC/5hmC 分軌當 *characterization/tooling* 未填——但正因甲基是 *characterization* 訊號非 *filter*，正是我們 locked 定位。**誠實限制**：filter-negative 單 pipeline(全 ClairS-TO，無第二 caller 重驗)；ASM 幅度單樣本-HCC1395-主導(⭐3，⭐4 需 COLO829 truth)；未扣 imprinted DMR；未測甲基當*邊際* DL 特徵(領域唯一報有益的配置)。

**相關**：[08 跨問題鏈](08_cross_question_correlation.md) ｜ [09 參考文獻](09_references_and_provenance.md) ｜ 既有 filter-DEAD 結論：[../09_conclusions/00_index.md](../09_conclusions/00_index.md)
