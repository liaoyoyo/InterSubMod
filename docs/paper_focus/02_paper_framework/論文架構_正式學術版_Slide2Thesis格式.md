<!--
建立時間: 2026-06-09（v2 2026-06-09 晚：併入所有裁決 — 雙主軸 rule+exception / de-confound floor / catalog R6 / characterize 非 reshapes / same-hap 3/6 修 / 🟡→🟢 T-PROV / GB venue / pre-reg+fallback）
狀態: architecture (正式學術論文架構 — Slide2Thesis Nature-journal 格式, Pandoc-ready)
報告類型: paper_focus_formal_paper_architecture
受眾: 廖子游 · PI · 投稿用
format: Slide2Thesis (github.com/ythuang0522/Slide2Thesis) Nature-journal style (sn-jnl.cls) + YAML metadata + Pandoc/Tectonic + BibTeX(sn-nature)
target_venue: Genome Biology / Bioinformatics tier（用戶 2026-06-09 定）
provenance_note: de-confound 主軸數字 🟢（catalog wf_5644ed77-082 6/6 驗 + 本 session grep）；co-validation 例外 🟡（contingent on #21/#22）；同 §13。
-->
<!-- provenance-verified: catalog 數字本 session grep（A=1/C=12868/G=291518/total 332705）；T-PROV 全 🟡→🟢；co-val/phasing 標 contingent。 -->

# 正式論文架構 v2（Slide2Thesis Nature-journal · Pandoc-ready · 雙主軸）

> **L0**：carrying claim = **雙主軸＝「規則 + 例外」**：①**規則（de-confound，🟢 floor，今天可寫）**：read-level 甲基 allelic 分群跟隨 **germline** haplotype 而非 somatic（catalog 332,705 位點：12,868 reliable 全 germline-allelic、僅 **1/332,705** 乾淨 somatic cis、甲基不能 filter）；②**例外（co-validation，🟡 contingent on #21/#22）**：chr17/TBC1D16 + phasing 軸三訊號正交佐證。
> **L1**：① **動詞用 characterize/map 不用 reshapes**（無功能資料）；② **目標 Genome Biology tier**；③ 數字 🟢 已 T-PROV 對賬 / 🟡 contingent；④ 敘述弧 = ASM 存在 → 邊界(不能 filter) → **catalog 規則(de-confound)** → 例外(chr17+phasing) → 跨樣本；⑤ pre-reg+fallback 見 `pre_registration_與fallback_carrying_claim.md`，**任一分支論文站得住（de-confound floor）**。

---

```yaml
---
title: "Read-level characterization of allele-specific methylation in long-read cancer genomes: methylation tracks germline rather than somatic haplotypes, with a rare somatic co-validation exception"
title_zh: "長讀癌症基因體中等位基因特異甲基化的單讀層級 characterization：甲基化跟隨 germline 而非 somatic 單倍型，以及一個罕見的 somatic 共驗證例外"
author:
  - {name: "Tzu-Yu Liao", affiliation: 1}
  - {name: "[PI / co-authors]", affiliation: 1}
affiliation: [{id: 1, name: "[Institution]"}]
keywords: [long-read sequencing, nanopore, allele-specific methylation, loss of heterozygosity, haplotype, somatic mutation, de-confounding, cancer epigenomics, read-level]
target_venue: "Genome Biology / Bioinformatics tier"
bibliography: references.bib
bibliographystyle: sn-nature
documentclass: sn-jnl
classoption: [pdflatex, sn-nature]
abstract_en: >
  Long-read sequencing resolves somatic mutations, haplotypes and 5mC/5hmC on single
  molecules, yet whether allele-specific methylation (ASM) in cancer reflects somatic events
  has not been characterized at read level nor de-confounded from germline-allelic and
  copy-number background. We introduce ISM, a read-level pipeline (read×CpG matrix, read–read
  distance, reliability-gated hierarchical clustering, haplotype/LOH tagging, normal-anchored
  cis-test), and catalogue 332,705 loci across six cancers. The rule: reliable methylation
  read-clustering exists at 12,868 loci but is overwhelmingly germline-allelic; genome-wide
  only 1/332,705 loci is a clean somatic cis-ASM, and methylation does not discriminate true
  from false somatic variants (four independent lines). The exception: at that single locus
  (chr17/TBC1D16) and along the LOH-constrained phasing axis, somatic, haplotype and
  methylation signals co-validate the same reshaping event. Copy-number dosage is not the
  driver of clustering magnitude. We deliver a read-level de-confounding of cancer ASM and a
  mechanism-resolved catalogue of what methylation co-validates versus cannot filter.
abstract_zh: "[中文摘要 — 投稿前定稿]"
---
```

---

## Abstract（結構式）

- **Background**：long-read 同時讀 somatic/haplotype/甲基；但癌症 ASM 是否反映 somatic、能否與 germline-allelic / copy 背景區分，未在 read-level characterize。
- **Methods**：ISM read-level pipeline + 6 樣本 332,705 位點 catalog + 驗證機制（LOSO/null/subhap-matched/ARI ruler/signed-rank/copy-test/germline-het null）。
- **Results（規則）**：reliable 甲基分群存在於 12,868 位點但壓倒性 germline-allelic；全基因組僅 **1/332,705** 乾淨 somatic cis；甲基不能判別變異（四道）。copy-dosage 非 driver【🟢 MW p=0.6183】。**（例外）**：chr17/TBC1D16 + phasing 軸三訊號 co-validate【🟡】。
- **Conclusions**：read-level de-confounding of cancer ASM + 機制解明的「能 co-validate、不能 filter」目錄。

---

## 1. Introduction
> **開場**：Long-read resolves somatic variation, haplotypes and 5mC/5hmC on one molecule, enabling a read-level test of whether cancer ASM is somatic or germline-allelic [@background].
- **1.1 重要性**：read-level 連結三訊號 + de-confound apparent ASM。🔵
- **1.2 缺口**：(i) read-level somatic-controlled de-confound 未做；(ii) 甲基-as-filter 未公開證偽；(iii) somatic-controlled ASM 樣貌未界定。🔵 [@MartinTrujillo2017]
- **1.3 貢獻**：**規則**（de-confound：甲基跟 germline 不跟 somatic、不能 filter）+ **例外**（chr17 + phasing co-validation）+ 誠實 catalog + ISM（vs 現有工具，#23 benchmark）。
- *(Fig 1：rule+exception 整合圖 — LOH→germline-allelic 主背景 + 罕見 somatic co-opt(chr17)；每軸標證據強度 🟢/🟡)*

---

## 2. Results（雙主軸：規則 floor → 例外 ceiling）

### 2.1 ASM exists and is somatically enriched
TP significant 3.95% > FP 1.07%（subhap-matched 3.77 vs 1.09 ~3.5×）【🟢 ledger:85】；但 ~4% 低 sensitivity。*(Fig 2)*

### 2.2 Boundary — methylation cannot discriminate true/false somatic variants（四道）
(a) LOSO 100% circularity【🟢 +0.02236→−0.00012; mean −0.00004, Wilcoxon p=0.125】[@Kapoor2023]；(b) |Δβ| AUC=0.505 落 null【🟢】；(c) ~4% sens + COLO829 TP≈FP【🟢 0.0089≈0.0103】；(d) regression-to-extreme（OR no-clust 8.63 / LOH 4.09 / subhap 5.84；dbeta_only TP:FP=11.9 最 FP-enriched）【🟢 T-PROV 對賬】。*(Fig 4)*

### 2.3 ⭐ RULE — read-level catalogue: methylation clustering is germline-allelic, not somatic（de-confound floor）
6 樣本 **332,705 位點 → 7 TAG**【🟢 catalog wf_5644ed77-082 6/6 驗 + grep】：reliable 甲基分群 **12,868（3.87%, TAG-C）** 但 cis-test 證**全是 germline-allelic 背景**；latent **28,254（TAG-E）** PERMANOVA 救回但 enrich≈base（無判別力）；無訊號 291,518（87.6%）；**乾淨 somatic cis-ASM 全基因組僅 1（TAG-A=chr17）**。→ **甲基 allelic 分群是 germline-driven，不是 somatic 判別器。** *(Fig 3: P5 per-TAG 長條 = 主圖)*

### 2.4 Copy-dosage is not the driver; cis vs copy partition
MW neutral-vs-gain p=0.6183；signed ρ(|Δβ|,CN)=−0.0829 → REFUTED【🟢】。816 HP-axis loci copy-partition：唯一乾淨 cis = chr17（within 0.142 ≫ d_copy 0.029, perm p=0.001）【🟢 T-PROV】；BRCA2 = copy-confounded（d_within −0.023 ≪ d_copy −0.11）【🟢】。*(Fig 6)*

### 2.5 EXCEPTION — chr17/TBC1D16 three-signal co-validation exemplar（🟡 contingent）
在唯一乾淨 cis 位點，somatic SNV / haplotype / 甲基三訊號正交佐證同一 LOH 事件。**生物動機（L3，加分）**：TBC1D16 是已知癌症甲基化基因——**hypomethylation 重活化 cryptic 轉錄本 TBC1D16-47KD → RAB5C/EGFR 驅動黑色素瘤轉移**（Vizoso/Esteller *Nat Med* 2015；*Cancer Discovery* 2015；*Clin Epigenetics* 2019）；方向 **hypo** 與本研究「hypo≠canonical hyper」一致。**誠實 caveat（必隨行）**：(a) 文獻是黑色素瘤+bulk hypomethylation，我們是乳腺 HCC1395+allele-specific cis-ASM → **跨癌種 + bulk-vs-ASM**，措辭用「與已知 TBC1D16 甲基生物學一致」非「證實」；(b) 🔴 **perm p=0.001 跨 816 cis-test，Bonferroni 後可能不顯著（#24 FDR 校準）**；(c) 單 locus 單樣本、HCC1395-private somatic 不可在第二乳腺樣本重現。⭐ **COLO829（黑色素瘤＝TBC1D16 主場）= 理想生物驗證樣本**（T-VAL-3）。*(Fig 5)*

### 2.6 LOH-constrained phasing axis（🟡 contingent on #21 R-SELFREF）
NG=2 Inner same-HP1 > Outer，**方向 7/7**（W=28, p=0.0078125）【🟢 master_draft:81】；paired 負控 NS【🟢】。**caveat**：by-construction 部分循環（R-SELFREF #21 待跑，pre-reg PR-1）；**same-hap 中位 ~92% 但僅 3/6 ≥93%**（0.840/0.939/0.759/0.429/0.932/0.920，T-PROV 修正，非 6/6）；single-pipeline ⭐3。AF→NGroups = phasing 非甲基【🟢 HD-4, r=0.656】。*(Fig 5)*

### 2.7 Cross-sample reproducibility（phenomenon-level）
6/6 excess-over-null +0.101–0.241（mean 0.168）【🟢 T-PROV *_gwasm.json】= phenomenon 復現；同-locus 0/38【🟢】但 UNDERPOWERED（E[overlap]=0.16）。[@ONeill2024]

---

## 3. Discussion
- **3.1 規則的意義（de-confound 貢獻）**：甲基 allelic 分群 read-level 主要是 germline；apparent ASM 多被 germline/copy 混淆 → 對「ASM 普遍 hypo 主導」文獻提供 somatic-controlled de-confound 視角【守 3 不對齊】[@Do2020]。
- **3.2 例外的意義 + 邊界**：chr17/TBC1D16 是罕見 somatic co-opt 甲基的 exemplar（讓規則完整），且落在已知癌症甲基驅動基因（hypo→EGFR 黑色素瘤轉移）[@Vizoso2015]；甲基不能 filter 的機制（circularity/regression-to-extreme/germline-allelic）= 防雷貢獻 [@Kapoor2023][@MartinTrujillo2017]。
- **3.2b 與 epipolymorphism 區隔**：epipolymorphism/methylation-entropy [@Landan2012] 測 read-level 異質性 —— 我們自家 O11 證實其判別力是 **n_reads-confounded artifact（AUC 0.845→0.530 after n_reads correction）**；ISM 不靠熵，而是 **reliability-gated 分群 + normal-anchored cis-test** 控掉該 confound → 把 read-level 甲基異質性從 confounded 描述升級為 de-confounded characterization。
- **3.3 Tooling 展望（future）**：甲基追蹤 haplotype 結構 → 開啟甲基輔助 phasing（umtag）【🟢 held-out 0.8852 / V10 非-copy，T-PROV 對賬】；需 yardstick/COLO829 (#23/T-VAL-3)。
- **3.4 Limitations（誠實）**：single-pipeline ⭐3；機制單樣本 HCC1395；chr17 例外單 locus 不可複製；phasing by-construction（#21 待）；6 cis untestable；matched normal 僅 HCC1395。

---

## 4. Methods（Nature 置 Discussion 後）
- **4.1 Samples**：HCC1395（SEQC2 truth）+ HCC1937/1954/H1437/H2009/COLO829；ONT 5mC/5hmC；ClairS-TO（single-pipeline）。
- **4.2 ISM pipeline**：read×CpG 矩陣 → read–read 距離(6 度量) → 階層分群 + **CramersV reliability gate**（Cochran≥5；三態 reliable/latent/none，ISM-2 audit）→ HP/LOH tag → **normal-anchored cis-test**。細節 `ISM_工具方法細節與應用意義.md`；vs 現有工具 #23 benchmark；與 epipolymorphism 區隔見 §3.2b。
- **4.2b 修飾範圍（5mC primary）**：主結論以 **5mC** 為準（well-established）；ISM 亦接受 5hmC 輸入，但 5hmC 生物意義不同（常 active gene-body），**未詳細驗證，列 future observation**；5mC 純度需 MSA Level1 dup-bug 修（#25）。
- **4.3 Catalogue**：16-欄 schema → 7-TAG（scripts 98-100；catalog_skeleton.tsv 332,705 列）。
- **4.4 驗證機制**：LOSO / 2000× shuffle null / subhap-matched / blind-ARI+imprinting 正控 / signed-rank+bootstrap+paired 負控 / copy-test/copy-partition / germline-het null。（→ `04_知識_驗證機制_研究機制_對應表.md`）
- **4.5 統計 + FDR 校準（#24，必補）**：🔴 **reliability 計數(12,868) + chr17(perm p=0.001 跨 816) 都須 null/FDR 校準**——label-shuffle null 重算 reliable rate + BH-FDR across 332,705 / chr17 Bonferroni(×816≈0.82)+genome-wide permutation；**必含 n_reads/coverage 校正**（O11 教訓：epipoly AUC 0.845→0.530 after n_reads correction）。佐證:Storey q-value/empirical-Bayes + TBC1D16 文獻 + COLO829 生物複製。effect size + CI；NEGATIVE 走 falsification + pre-registration（`pre_registration_與fallback_carrying_claim.md`）。

---

## Pre-registration & Fallback
#21 R-SELFREF + #22-reframed 跑前判準已寫死（PR-1/PR-2）；fallback tree：**de-confound floor 永遠成立**，#21/#22 任一 outcome 論文都站得住。詳 `pre_registration_與fallback_carrying_claim.md`。

## Data & Code ／ References
- Data：ONT BAM + SEQC2 truth；catalog `research/.../genome_survey_v2/catalog/`。Code：ISM(C++17+Python) + scripts 98-100。
- References：`references.bib`（→ `knowledge/11`，/citation-verification）：[@Kapoor2023][@MartinTrujillo2017][@Do2020][@ONeill2024][@MethPhaser][@Landan2012][@Derrien2021]。

---

## L1 — 投稿前 checklist
- [ ] de-confound floor 主文（§2.1-2.4 + §2.7）可即寫（全 🟢，不等 #21/#22）
- [ ] #21 R-SELFREF 回來 → 依 fallback tree 定 §2.6 phasing（positive / 降 characterization）
- [ ] #22-reframed 確認 chr17 唯一 + germline 背景量化
- [ ] same-hap 口徑改 3/6（非 6/6）— 全文同步
- [ ] #23 ISM benchmark → Methods §4.2 tooling 定位
- [ ] Fig 1–8 產出（P1-P6 ✅ 已有；Fig 1 rule+exception / Fig 8 chr17 待製）
- [ ] /citation-verification 核 PMID/DOI；FDR 正規化；COLO829 升 ⭐4（資源到位）
