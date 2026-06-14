---
title: "ISM 顯著性方法盤點對照表 + 方法學評估（輸出欄→計算→目標→是否最有效→替代）"
date: 2026-06-14
type: methodology-inventory
source: 3 並行 agent 源碼盤點 (GlobalTest / LabelTest+Δβ / PERMANOVA+PerCpgAsm) + permutation 數校正
claim_levels: 計算邏輯=L1(源碼) / 方法學評估(最有效·替代)=L2(領域知識,替代方法待 KB/文獻最終確認)
---

# ISM 顯著性方法盤點 + 方法學評估

## L0 盤點總結
ISM 顯著性分 **6 類方法**，分屬 4 個目標層（存在性檢定 / 方向描述 / 定位 / 異質性 + 2 混淆守衛）。
**大部分計算合理（Fisher-FH / PERMANOVA-Anderson / BH-FDR / NME-Epipoly 都對標標準）**，但盤出
**6 個方法 gap/缺陷**（最關鍵：🔑**Δβ 有實作但無檢定**、Dispersion p 查表近似、permutation 數不一致）。

## 1. 盤點對照表（輸出欄 → 計算 → 目標 → permutation）

| 方法 | 計算（L1 源碼）| 主要 summary 輸出欄 | 目標問題 | perm 數 / p-floor |
|------|------|------|------|:---:|
| **GlobalTest** | cluster_label × {allele/HP/HP-family/HP-fine/sample} contingency → **Fisher-Freeman-Halton MC**(p) + **Cramér's V**(effect)；chi-sq 是 placeholder | `GlobalP`(min), `CramersV`(max, Cochran-gated), `GlobalP/CramersV_HPFamily`, `GlobalP/CramersV_HPFine`, `HPFine_NGroups_CF` | cluster 是否只是重現既有標籤(allele/HP/T-N)的分群？(下游 gate) | Fisher MC |
| **LabelTest ① distance** | `delta = between_dist_mean − within_dist_mean`；delta>0 才 **permutation**；p=n_extreme/(perm+1) | `HPMergedDelta/P/Sig`, `HPFineF/P/Sig`, `AlleleDelta/P/Sig`, `SampleASM_*`, `HP_Residual_Delta/P/Sig`, `Tumor/Normal_HP_Delta`, `Unassigned*` | HP/allele/sample 是否解釋甲基**距離結構**(存在性) | **999** / 0.001 |
| **🔑 Δβ ② signed** | `mean(HP1 甲基率) − mean(HP2 甲基率)`(raw_matrix, **無 permutation/p**) | `Tumor/Normal_HP_Signed_Delta`, `HP_Signed_Residual`, `Combined_HP_Signed_Delta` | ASM **方向+大小**(hypo/hyper) | **無檢定** |
| **PERMANOVA(StructureTest)** | Anderson 2001 SS 分解(raw-accum 避 cancellation) + `pseudo-F=(SSb/dfb)/(SSw/dfw)` + **shuffle-label** permutation；NaN 先 `filter_reads_for_complete_matrix` | `ClusterPermanovaF/P/Valid`, `LabelHP_Permanova*`, `LabelAllele_Permanova*` | read 距離空間是否**整體有結構分群**(cluster/HP/allele 三套) | **99** / 0.01 |
| **Dispersion(betadisper)** | 每 read 到群心距離 proxy → one-way ANOVA F → **查表 p**(F>4→.01/F>2.5→.05/else .1) | `ClusterDispersionP/Warn`, `LabelHP/Allele_Dispersion*` | 混淆守衛：PERMANOVA 顯著是否其實是**離散度差**非位置差 | 查表(非 perm) |
| **PerCpgAsm** | per-CpG 2×2(HP×meth/unmeth, β>0.5) **Fisher exact + BH-FDR**；**NME**=Shannon(read pattern)/Hmax；**Epipoly**=1−Σpᵢ²(4-CpG window) | `PerCpgASM_Valid`, `Fisher_N_Sig/Frac_Sig/N_Tested/MaxNegLogFDR`, `NME_HP1/HP2`, `Entropy_Imbalance`, `Epipoly_HP1/HP2/Delta` | **哪些 CpG** HP 差異顯著(定位) + 各 HP epiallele **異質性/熵** | Fisher exact |
| **CramersV reliable gate** | Cochran：>20% 格期望<5 → unreliable → 下游當 0 | (gate, 影響 CramersV) | effect size 可信度守衛 | — |

## 2. 🔑 Δβ（delta β）確認 — 用戶核心問題答案
**ISM 有兩條獨立的線，必須分清**：
- **① distance line**：所有 `*P`/`*Sig` 欄（HPMergedP/Sig, AlleleP, SampleASM_P, HP_Residual_P…）= 距離 between−within delta 的 **permutation 檢定**。**這是唯一顯著性來源**。
- **② Δβ line**：`compute_signed_hp_delta`（RegionProcessor:997-1016）= **直接甲基率差** mean(HP1)−mean(HP2)，對應 4 個 `*_Signed_*` 欄。**有實作，但純描述、無 p-value/sig 欄**。
- **⚠ 命名陷阱**：`HP_Residual_Delta`(:976,距離 delta 相減) ≠ `HP_Signed_Residual`(:1031,Δβ 相減)，名字像本質不同。

→ **結論**：Δβ **已直接實作（描述方向/大小），但「Δβ 本身的顯著性」沒有檢定**（顯著性都走 distance line）。若要「這個位點 Δβ 是否統計顯著」，目前 ISM 給不出——這是 P2b 之外的另一個 gap（見 §3 #1）。

## 3. 方法學評估：是否最有效 + 替代方法（L2）

| 方法 | 是否最有效 | 替代方法（業界）| 判定 |
|------|------|------|:---:|
| GlobalTest Fisher-FH+CramersV | 合理：Fisher exact 對小樣本/稀疏 contingency 穩健 | log-linear / mutual-info / ARI(cluster-vs-label) | 🟢 OK |
| LabelTest ① distance delta | 可用但非最標準：between−within 不如 pseudo-F 標準 | **PERMANOVA**(ISM 已另有) / ANOSIM / MRPP | 🟡 **與 PERMANOVA 冗餘**(#6) |
| 🔑 Δβ ② signed(無檢定) | 描述 OK，**缺檢定** | **beta-binomial GLM (DSS) / logistic (methylKit) / M-value+limma** | 🔴 **gap #1** |
| PERMANOVA(Anderson) | ✅ distance-multivariate 金標準，實作正確 | db-RDA / MANOVA(需常態) | 🟢 OK |
| Dispersion 查表 p | 概念對(betadisper 標準用途)，**p 不嚴謹** | **PERMDISP2**(真 permutation p) | 🔴 **gap #2** |
| PerCpgAsm Fisher+BH | 合理，但 Fisher 把 read 當獨立、**忽略 over-dispersion** | **beta-binomial(DSS) / methylKit** 更 powerful | 🟡 **gap #4** |
| NME / Epipoly | ✅ 標準 epiallele 指標(Jenkinson2020 / Li2014) | PDR / MHL | 🟢 OK |

## 4. 彙整：顯著性計算判定的方法問題清單（你要的「合理沒大問題嗎」）
**整體判定：核心計算合理、無致命錯誤；但有 6 個 gap/缺陷需逐一決定修不修。**

| # | 問題 | 嚴重度 | 對應修法 |
|---|------|:---:|------|
| 1 | 🔑 **Δβ 無顯著性檢定**（signed delta 只描述）| 🔴 高 | 加 beta-binomial / logistic 對 Δβ 檢定（DSS-like）|
| 2 | **Dispersion p 查表近似**（F>4→0.01 非真分布）| 🔴 中 | 改 PERMDISP2 真 permutation p |
| 3 | **hp_residual_sig 借 tumor p**（gap2 已揭，殘差未檢定）| 🔴 高 | 對殘差做正規 permutation（P2b）|
| 4 | **PerCpgAsm Fisher over-dispersion**（read 當獨立）| 🟡 中 | beta-binomial per-CpG（memory 全碼稽核已提）|
| 5 | **permutation 數不一致**（PERMANOVA 99 / LabelTest 999，p floor 0.01 vs 0.001）| 🟡 中 | 統一或文檔明示理由 |
| 6 | **LabelTest delta + PERMANOVA 冗餘**（同測 HP/allele 兩套）| 🟢 低 | 釐清分工或擇一 |

## 5. 下一步（壓縮資訊後詳細研究）
建議從**最高槓桿 + 最上游**的 #1（Δβ 檢定）或 #3（hp_residual 殘差檢定）切入詳細研究——兩者都直接決定「ASM/somatic 訊號能否可信判定」。#2/#4 是改善統計嚴謹度。#5/#6 是釐清。
> 盤點完成 → 後續可只帶此表 + 問題清單繼續，不需重載 3 agent 的完整源碼回報。
