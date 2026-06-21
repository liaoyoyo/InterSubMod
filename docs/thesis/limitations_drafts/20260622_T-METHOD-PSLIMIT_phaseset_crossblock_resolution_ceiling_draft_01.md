---
title: "方法/討論限制段草稿 — phase-set 跨 block 無 read/cell 連結（重建解析度天花板）"
date: 2026-06-22
status: draft
task: T-METHOD-PSLIMIT
audience: 論文作者（Methods 限制 + Discussion 撰寫材料）
build_branch: docs/limitations
data_sources:
  - docs/methodology/20260620_somatic_locus_methylation_combination_enumeration_01.md
  - docs/methodology/20260622_ism_method_cognition_and_open_questions_01.md
  - docs/presentations/in_progress/20260603_methyl_haplotype_story_deck_v2/refs_v2.bib
related_memory:
  - project_somatic_locus_methylation_combination_interpretation
provenance_note: >
  本檔為限制段草稿，以引用既有結論為主；不重算數字。
  外部 PMID/DOI（MethPhaser 38909018 等）投稿前過 /citation-verification（見 T-E3）。
---

# 方法/討論限制 — phase-set 跨 block 觀測單元天花板（draft）

> **這段要解決的問題**：本論文以「同一條長 read 上多個位點 + 鄰近甲基的共現」作為單分子（read-level）證據，並在 phase block 內把 read 分群、做 read×read 距離結構檢定。問題是：**這個能力的有效範圍 = 一個 phase set（phase block）之內**。跨越 phase block 邊界時，bulk 長讀資料**沒有任何 read 或 cell 把兩個 block 連起來** → genome-wide 的譜系/lineage 重建在本資料型態下**有一個物理上限**。Ch5/Methods 必須明白界定這條天花板，避免讀者誤以為論文做的是全基因組 phylogeny。

---

## 5.x.1 機制：為什麼跨 phase block 沒有連結（一段話版本）

長讀 phasing 把基因組切成一段段 **phase block**：每個 block 內部，雜合位點可被同一批跨越它們的 read 連續地指派到 H1/H2，因此 block 內的「哪些位點共在同一條 allele 上」是可解的。但**相鄰兩個 block 之間沒有任何單一 read（或單一細胞）同時跨過兩者** —— 一旦 read 長度 / 雜合位點密度耗盡，phase 就斷開，**H1 在 block A 與 H1 在 block B 之間的對應關係是未知的**（隨機方向）。

> 換句話說：**block 內**，read 提供「跨位點的物理 linkage」；**block 間**，這個 linkage 消失。bulk 資料只有 block-local 的單分子共現，沒有 block-to-block 的單分子或單細胞鏈接。

外部依據：MethPhaser（Fu et al.）明確指出**相鄰 phase block 之間沒有連結資訊**，並以甲基相關性作為「延長 phase block / 跨 block 定相」的輔助手段——但其解決的是「**把 block 接起來以延長 phasing**」，**不是**提供 block 間的 read-level lineage linkage。引用：`[@MethPhaser]`（Nat Commun 2024，PMID 38909018 待 /citation-verification 核對；參見 T-E3 條目）。

---

## 5.x.2 對本論文方法解析度的界定（誠實底線）

1. **本論文的單分子證據與 read×read 距離結構檢定，作用域 = 單一 phase set / 局部區域**。所有「同一條 read 上多 sSNV + 甲基共現」的論證（含 chr2:18M 原型）都受此邊界約束。
2. **跨 phase block 的 genome-wide lineage 重建在本資料型態下無法直接達成** —— 需要外部錨（CCF/CNA/phylogeny 推論）或 single-cell 資料把 block 串起來。這正是論文不宣稱「完整 phylogeny / CCF tree」的根本原因。
3. **這直接界定 reconstruction 用詞**（與認知文件 §5 紅線 2 一致）：本論文 reconstruction = **regional、LOH-constrained partition**（受 phase-set / 觀測單元天花板限制），**非** genome-wide tree。
4. **單位點組合窮舉的適用範圍同受此限**：D1-D7×I1-I8 的詮釋只在「同一 phase set 內」有效；跨 phase set 因無 read/cell 連結，須宣告為 **reconstruction GAP**，論文**不宣稱 ISM 已解此 gap**（只 map，不下「已解決」貢獻）。

---

## 5.x.3 觀測單元樞紐（把限制講成「為什麼難」而非「我們的缺陷」）

| 觀測單元 | 跨 locus linkage | 跨 phase block linkage | 適用 regime |
|---------|:---:|:---:|------|
| **single cell** | ✓（細胞內全基因組共存） | ✓（同一細胞跨 block） | single-cell methylation lineage（Gaiti/Epiclomal/MethylTree）|
| **long read（bulk，本論文）** | ✓（read 內，block-local） | ✗（**block 間無連結**）| **regional 單分子共現 = ISM scope** |
| short read（bulk） | 受限（read 短） | ✗ | — |

> **可移植 vs 不可移植（與 transferability 裁決一致）**：能移植到本論文 bulk ONT 的，是 **phase-block 內 regional 單分子共現**（ISM 的 scope）；**不能**移植的是 genome-wide lineage（那是觀測單元的物理天花板，不是演算法可補的洞）。

---

## 5.x.4 緩解方向（future work，不宣稱已做）

- **single-cell methylation**：唯一能在細胞層提供跨 block linkage 的資料型態（終判，但非本論文資料）。
- **外部 lineage 錨**：CCF 梯度 / CN-LOH 結構 / phylogeny 推論把 block-local partition 串成更大尺度的假設（仍是推論，非單分子證據）。
- **甲基輔助延長 phase block**（MethPhaser 類方法）：可**延長** phasing 連續性，**但仍不提供** block 間的 read-level lineage —— 引用時須區分「延長 phasing」與「提供 lineage linkage」兩件不同的事，避免 over-claim。

---

## 待填 / 缺資訊（投稿前補）

- `{{待填}}` 本論文資料中 phase block 的 N50 / median 長度與「跨 block 候選位點」比例（須在對應 branch 上算；本草稿不引用未驗證數字）。
- PMID 38909018（MethPhaser）為投稿前 /citation-verification 必核項（見 T-E3）。
- 注意 refs_v2.bib 已有一條 `fu2024methphaser`（Nat Commun 15:5327, DOI 10.1038/s41467-024-49588-0，**無 pmid 欄**）；本段引用的 PMID 38909018 須與該 bib 條目的身分核對一致後再定稿（見 T-E3 §需 web 核對清單）。
