---
title: ISM 計算複雜度定位 — 一頁速覽 CHEATSHEET（Ch2 撰寫用）
date: 2026-07-05
status: in_progress
task_type: summary-cheatsheet（本 session 計算複雜度 arc 收斂速覽；非新分析）
build_branch: research/subclonal-reconstruction-202606
verification_stamp: >
  2026-07-05 fresh 驗證（§13.7 完成 gate，實跑非憑記憶）：4 commits 在(f56d215/d2a5eba/cce3854/eaf77c7)、
  working tree 乾淨、8 新卡實存、庫 CONTEXT=82 / axis1=26、tier B×44、.bib=39 well-formed 條目、
  殘留 "Mahapatra 2024"=0、PhISCS PMID=31628256(已修)、memory 步驟1-3 已寫。
data_sources: >
  docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md,
  docs/method_comparison/20260705_ch2_related_work_computational_positioning_draft_01.md,
  docs/method_comparison/_assets/ism_positioning_refs.bib
related:
  - docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md
  - docs/method_comparison/20260705_ch2_related_work_computational_positioning_draft_01.md
---

# ISM 計算複雜度定位 — 一頁速覽

> **一句話**：ISM 的**可解階段踩在已被證明多項式可解的數學島**（Gusfield perfect-phylogeny + Pe'er 2004 IDP），**拒絕/列舉的階段全是已被證明 NP-hard 的島**。→「⭐3 characterization + 定不出來即答案」= **對齊真實複雜度分界的原則性設計**，非保守選擇。
>
> **本檔用途**：Ch2 撰寫的 cheat-sheet + 本 session（形式化 → 數學島 → 最近鄰 → citation → 建卡 → Ch2 草稿）收斂記錄。深入見 SoT `20260705_ism_computational_complexity_positioning_01.md`。

## 1. 形式化（三層誠實裁決）

ISM subclone 重建 = 布爾超立方體 $\{0,1\}^n$ 上**有向有根 perfect phylogeny**（頂點=二元 sSNV 基因型 R=0/A=1，$n\le8$；根=germline $0^n$；邊 $0\to1$ 單調；four-gamete=Gusfield 相容性；$H_*$ 隱藏節點=Steiner 點）。

| 層 | 裁決 | 要點 |
|---|---|---|
| ① hypercube + directed + rooted-at-RR | ✅ **準確** | max-parsimony over binary ≡ Hamming hypercube 樹 |
| ② "Steiner tree" | ⚠ **物件準/計算不準** | $H_*$=真 Steiner 點；但 ISM 只算**多項式 perfect-phylogeny 特例**，不解 NP-hard Steiner |
| ③ "group Steiner" | 🔴 **過度宣稱 as solved** | partial-read subcube + multiplicity 是真結構類比，但 ISM 不做 cover-minimize；白地 |

🔴 **腫瘤 casting 成 hypercube-Steiner = 原創綜合**（通用版 Mahapatra 2025〔*Acta Inf* 62(1):6〕有；腫瘤文獻用 El-Kebir 2015 spanning arborescence，ISA 下無 latent 節點故非 Steiner）。

## 2. 數學島（誰解過類似問題）

| ISM 階段 | 數學問題 | 複雜度 | 島 / Citation |
|---|---|---|---|
| 骨幹（完整 read） | directed perfect phylogeny | poly ✅ | Gusfield 1991 |
| 🎯 **partial-read 補全（X=缺失）** | Incomplete Directed Perfect Phylogeny | poly→linear ✅ | **Pe'er 2004** SIAM J Comput |
| 選丟 loci | Min Character Removal ≡ Vertex Cover | NP-hard, FPT-in-k | Day-Sankoff 1986 |
| 拆混合 bulk row | Min Conflict-Free Row Split | NP-hard + inapproximable | Hujdurović 2018 |
| 錯誤翻轉去衝突 | Min-flip | NP-hard | Chen 2006（⚠非 Kolmogorov） |
| Steiner 補全 hypercube | Steiner in {0,1}ⁿ | NP-/MAX-SNP-hard, FPT | Foulds-Graham 1982；Mahapatra 2025 |
| 允許甲基 loss（1→0） | Dollo-k | k=1 poly / k≥2 NP-complete | Bonizzoni（arXiv） |
| multiplicity/CCF 唯一性 | non-identifiable | impossibility（證明） | DeCiFer 2021 |

**🎯 兩個最強論述**：
- **partial-read 有嚴格基礎**：Pe'er 2004 IDP 使「救 40.4% 空向量」= 有理論保證的 incomplete-perfect-phylogeny 補全（可行性 L1-solid；**選哪個完成** non-identifiable，與 DeCiFer 一致）——非工程 trick。
- **繞過 NP-hard**：ISM 單分子共現**直接觀測**超立方體頂點；bulk VAF 法**推斷**=NP-hard VAFFP（El-Kebir）且需多樣本 → ISM 繞過短讀必解的反卷積。代價 $n\le8$/單樣本/hard-case 拒解=⭐3。

## 3. 定位口徑（reviewer 攻防）

**✅ 可宣稱**：① 實用性/成本（單一 bulk ONT，免單細胞/多樣本/linked-read）② 原生同分子甲基（甲基與 sSNV 同一條物理分子）③ 非循環設計（HP 來自零甲基 LongPhase-S；甲基 label-first bounded-auxiliary 不進 likelihood）。

**🔴 絕不可宣稱**：① 解析度/確認優越（單細胞 MethylTree/sgootr Q≈1、Gaiti 7.4e-9 + 多樣本 Pairtree 都贏）② 單分子共現為獨有（Foltz 2024 linked-read 已發表）③ 共現單獨=確認 subclone（necessary-not-sufficient，只互斥→分開可靠）④ 說 cvlr/ASMS/CpelNano「沒顯著檢定」（都有 permutation）⑤ ROCIT「是 ONT」（PacBio）⑥ MethylBERT「不吃 ONT」（code 吃）。

## 4. 最近鄰 4 名 + 唯一關鍵差

| 名 | 最近鄰 | 唯一關鍵差 |
|---|---|---|
| 1 | Foltz 2024（SomaticHaplotype） | 10X linked-read 非原生 ONT；零甲基；無結構檢定/樹 |
| 2 | sgootr 2023 | 單細胞；甲基 DRIVER；擁「distance-methylation reconstruction」方法類 |
| 3 | TumorLens 2026 | 偵測/註記器非重建器；窗式 DMR |
| 4 | LongPhase-S 2025 | ISM 自己骨幹工具；同實驗室 |

## 5. 產物地圖（全 committed，feature branch，未 push）

| 產物 | 路徑 / commit |
|---|---|
| 定位 doc（SoT） | `InterSubMod/docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md`（f56d215） |
| verified `.bib`（**39** 條，0 CITATION_NEEDED） | `InterSubMod/docs/method_comparison/_assets/ism_positioning_refs.bib`（d2a5eba） |
| +8 Steiner 學派卡（external，不 commit 本體） | axis1 18→26 / 庫 74→82（cce3854 = in-repo bridge 計數） |
| Ch2 相關研究草稿 | `InterSubMod/docs/method_comparison/20260705_ch2_related_work_computational_positioning_draft_01.md`（eaf77c7） |
| memory | `project_ism_computational_complexity_positioning`（已含步驟 1-3） |

**citation 關鍵修正（步驟 1）**：Mahapatra 2024→2025 / PhISCS PMID 31628257→31628256（原 SiCloneFit）/ Bonizzoni 作者 Soto Gomez·Carrieri 非 Dondi / Chen 2006 Kolmogorov 非作者 / cvlr「call variants」描述錯誤 / ROCIT=PacBio 非 ONT。

## 6. 投稿前剩餘

1. `.bib` **7 條 author=others** 補全（Foltz/Wakhan/Severus/MethPhaser/MethylTree/ROCIT/Epiclomal）+ ASMS 全作者（bioRxiv 403）。
2. 8 新卡**演算法內部細節** ⚠UNVERIFIED → 親讀全文 PDF 升 V2/V3（Pairtree scaling / Pe'er·AncesTree 複雜度歸約 / Mahapatra FPT 界）。
3. `\cite{}` 轉投稿 citation 樣式；全文再對照「絕不可宣稱」清單防 overclaim。
4. ⚠ 並行 session 正實作 group-Steiner（commit 68ced6f）——與本理論定位協調，避免 Track A 實作與定位口徑漂移。

## 關聯
- SoT：`InterSubMod/docs/method_comparison/20260705_ism_computational_complexity_positioning_01.md`
- Ch2 草稿：`InterSubMod/docs/method_comparison/20260705_ch2_related_work_computational_positioning_draft_01.md`
- CN 軸角色分工：`InterSubMod/docs/method_comparison/20260703_cn_axis_perread_feasibility_and_sSNV_methyl_role_01.md`
