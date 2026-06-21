---
title: "Ch5 限制段草稿 — LOH-unmask ASM confound 對 cis 主張的界定"
date: 2026-06-22
status: draft
task: T-C-UNMASK
audience: 論文作者（Ch5 Limitations 撰寫材料）
build_branch: docs/limitations
data_sources:
  - knowledge/11_external_literature/07_asm_cis_cancer_impact.md
  - docs/methodology/20260620_somatic_locus_methylation_combination_enumeration_01.md
  - docs/methodology/20260622_ism_method_cognition_and_open_questions_01.md
  - docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/refs_v2.bib
related_memory:
  - project_somatic_locus_methylation_combination_interpretation
  - project_ont_cnv_sv_subclone_verification_feasibility
  - project_O12_loh_methylation_scenarios
provenance_note: >
  本檔為限制段草稿，所有外部數字（82-91% 等）以引用既有結論為主，
  不在此 session 重算；數字來源見內文行內標註。投稿前外部 PMID/DOI 走 /citation-verification（見 T-E3 scaffold）。
---

# Ch5 限制 — LOH-unmask ASM confound（draft）

> **這段要解決的問題**：論文用 normal-anchored cis-test 主張「某些 somatic 變異與鄰近甲基有 cis 關聯」。在拷貝數中性區這個主張較乾淨，但在 **LOH（雜合性缺失）** 區，觀察到的「等位特異甲基（ASM）」可能**不是新生 cis 事件**，而是 LOH 把一條 allele 移除後、**揭露（unmask）了原本就存在但被雙等位平均掩蓋的 ASM**。這是論文最重要、必須在 Ch5 明確處理的 confound 之一。

---

## 5.x.1 機制：LOH 如何「揭露」既存 ASM（一段話版本）

在二倍體、拷貝數中性的位點，bulk 觀測到的甲基訊號是兩條 allele 的混合平均。即使兩條 allele 本身甲基狀態不同（即存在 germline/imprinted ASM），平均後常呈現中等、看似單一的 β 值。當該位點發生 **LOH**（一條 allele 被刪除或被另一條取代為 copy-neutral LOH/UPD），bulk 觀測就只剩存活的那一條 allele → 原本被平均掩蓋的等位差異**突然顯現**為強 allelic methylation 訊號。

> 關鍵後果：在 LOH 區看到的「強 ASM」**在時間上可能早於、且因果上獨立於** somatic 突變本身。它由「哪條 allele 還在」（拷貝數/LOH 結構）決定，而不是由突變誘發鄰近甲基改變（即不是我們想宣稱的 somatic cis）。

---

## 5.x.2 外部證據（量化 confound 的規模 — 引用既有結論）

| 證據 | 內容（方向 + 量） | 來源（投稿前過 /citation-verification）| scope 注意 |
|------|------------------|------------------------------------|-----------|
| Martin-Trujillo et al. 2017, *Nat Commun* | 腫瘤中 **imprinted DMR 的甲基改變約 82–91% 由 CN/LOH 解釋，而非真正的表觀突變（epimutation）** | `[@MartinTrujillo2017]` PMID 28883545 | 🔴 **scope 限 imprinted DMR**；泛化到全基因組 allele-specific methylation 須明寫「demonstrated at imprinted DMRs」 |
| Chase et al. 2015 | 甲基隨 acquired UPD（aUPD）的範圍尺度共變，r≈0.76 | PMID 26114957（L3，待核）| 支持「CN/LOH 結構驅動 allelic 甲基」 |
| Rosenski et al. 2025 | 提供 parental-ASM + segmental-duplication-ASM 座標，**可作 intersect blocklist** 把已知 ASM 區排除 | PMID 40069157（L3，待核）| 可操作的去 confound 工具（future work）|

> **數字校正紀錄（避免 over-claim）**：external_validation 庫的 4-lens 對抗稽核指出，Martin-Trujillo 的 82–92% 是 **CN-dependent 部分**；其互補的「真表觀（CN-independent）」部分僅約 8.1–17.6%，且同樣**限 imprinted DMR**。引用時兩個分母都要寫清楚，且**不可**把這篇拿來「corroborate」任一具體位點的 subclone/cis reclassification（那是 category error；見 `knowledge/11_external_literature/07_asm_cis_cancer_impact.md` 的 amend 註記）。

---

## 5.x.3 對本論文「cis 主張」的影響界定（誠實底線）

1. **CN-neutral 區的 cis 主張相對乾淨**：當位點落在雙 allele 存活（copy-neutral、非 LOH）區，「揭露既存 ASM」這條解釋路徑被關閉 → cis 訊號較可歸因於 somatic-linked 局部效應。**論文的 cis 主張應主要建立在 CN-neutral 區，並在文中明確 scope 化。**
2. **LOH 區的 ASM 不得單獨當作 somatic cis 證據**：在 LOH 區，`I2 = LOH-unmask 既存/imprinted ASM`（非新 somatic）與 `I3 = 真 somatic cis-ASM`、`I4 = somatic subclone` 在**單一 bulk 樣本下讀數退化、無法形式區分**（見單位點組合窮舉 D1-D7×I1-I8 的三句結論：最像 subclone 的組合 LOH×半帶 ALT×k≥2×對齊**全 undecidable**）。
3. **這直接界定 reconstruction 用詞**：本論文的 reconstruction = **regional、LOH-constrained partition**，不是完整 phylogeny/CCF tree（呼應認知文件 §5 紅線 2）。
4. **G-B（within-haplotype somatic null）未跑前禁宣 somatic-specificity**：因為 LOH-unmask confound 的 82–91%（imprinted DMR scope）尚未在我方資料上拆解，**Ch5 必須明寫「目前無法宣稱所觀測 ASM 為 somatic-specific」**（認知文件 §5 紅線 3）。

> 🔴 **紅線（與全論文一致）**：甲基在本論文的角色是 **characterize / corroborate**（刻畫讀層結構、佐證 haplotag 一致性），**不是** driver、不是 genome-wide tree reconstruction、不是 variant filter。LOH-unmask 段落不可被寫成「我們的甲基能偵測 subclone，只是被 LOH 干擾」——正確敘述是「LOH 區的甲基訊號連 cis vs unmask 都分不開，因此甲基在這裡只能 characterize，不能 attribute」。

---

## 5.x.4 緩解方向（標為 future work，不宣稱已做）

去除/降低此 confound 需要至少一項外部錨（單樣本 bulk 無法自證）：

- **normal cis-control**：用 matched normal 在同位點建立 germline/baseline ASM → 砍掉 `I2`（既存 ASM）這條路徑。這是論文宣稱「真 somatic cis」最關鍵的待補步驟。
- **multi-sSNV CCF 梯度法**（chr2:18M 原型）：用同一長 read 上多個 sSNV 的 CCF 階梯，分辨 `I3`（單事件 cis）vs `I4`（subclone）。
- **已知-ASM blocklist intersect**（Rosenski 2025 座標）：把已知 parental/SD-ASM 區先排除，降低 unmask 假陽性。
- **跨樣本 / single-cell** 終判（COLO829、single-cell methylation）。

> 落地參考（CN/LOH 工具與本機驗證現況）：SAVANA 已在 HCC1395 ONT 跑通並與 SEQC2 truth 高度一致（全基因組 LOH Jaccard 已由 het-site 汙染版 0.616 修正到 normal-het 版 0.962；見 memory `project_ont_cnv_sv_subclone_verification_feasibility`）。**但** CN 軸與甲基軸經 LOH-unmask 共同上游耦合 → **「CN 分群↔甲基分群對齊→驗 subclone」已被裁決為 REFUTED as standalone validator**；CN 只能當「固定觀察參考層」分流值得深究的區，不能當 cis 的去 confound 工具本身。SV 軸（split-read）是唯一對甲基 β 無因果耦合的非循環 anchor。

---

## 待填 / 缺資訊（投稿前補）

- `{{待填}}` 本論文資料中落在 LOH 區 vs CN-neutral 區的 cis-candidate 位點計數（須在 `feat/summary-nreadsvalid` branch 上算；本草稿不引用未驗證數字）。
- `{{待填}}` normal cis-control 跑完後，LOH 區「真 somatic cis」殘餘比例。
- PMID 26114957 / 40069157 為 L3 二手錨，**投稿前必過 /citation-verification**（見 T-E3）。
