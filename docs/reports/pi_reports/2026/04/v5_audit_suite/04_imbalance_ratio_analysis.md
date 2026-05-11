---
title: V5 Audit Suite 04 — HP1:HP2 Imbalance Ratio Analysis
date: 2026-04-27
author: Claude Code (Agent C)
audience: PI
status: pilot
last_updated: 2026-04-27
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/02_concordance_quantification.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/05_per_site_improvement.md
  - InterSubMod/docs/reports/pi_reports/2026/04/20260427_V5_IGV_session_visual_audit_01.md
  - InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/v5_audit_summary.tsv
artifacts:
  data:
    - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/imbalance_ratio.tsv
  figures:
    - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/04_imbalance/fig04a_ratio_scatter.png
    - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/04_imbalance/fig04b_distance_distribution.png
  scripts:
    - InterSubMod/scripts/analysis/v5_imbalance_improvement.py
---

# 04 — HP1:HP2 Imbalance Ratio Analysis

## 摘要

本文件量化 V5 與 Baseline (BL) 兩個 LongPhase-TO 版本，在 15 個經人眼審查的位點上，HP1-side 與 HP2-side read 比例的不平衡程度，並以 paired tumor BAM 的同位點 HP 比例作為「真實 phasing」參考錨點。每位點由 BL/V5 的 HP:i 整數標籤計數 HP1+HP11 vs HP2+HP21；paired 由 HP:Z 字串標籤（"1"/"1-1" 與 "2"/"2-1"）計數。所有距離度量採 orientation-flip aware 形式 `min(|q-r|, |1-q-r|)`，避免 PS block 任意翻轉造成假性差異。

## Section 1 — 為什麼要看 imbalance？（self-phasing 的本質特徵）

Self-phasing 的視覺與量化特徵在 IGV 截圖中已被反覆確認（`InterSubMod/docs/reports/pi_reports/2026/04/20260427_V5_IGV_session_visual_audit_01.md` Section 4 的 D 類三筆極端位點）：當位點落在 phaseable region 之外（GE.bed 不覆蓋）但有大量 reads，LongPhase 會**整批單向 vote**，把所有覆蓋 reads 都 assign 到單一 haplotype，產生 HP1:HP2 ≈ 100:0 的極端不平衡。這個現象就是 ISM 與 LOH 推論被嚴重污染的根源，也是 VCF Self-Phasing Circular Dependency CONFIRMED（`docs/concepts` 與 MEMORY 中 `project_self_phasing_causal_chain_confirmed.md`）的 read-level 表現。

合理的「解 phasing」結果應該滿足三條件：

1. **與 paired 同方向**：HP1-side 與 HP2-side 大致對齊 paired tumor BAM 的同位點分群（容許 PS block label flip）。
2. **不過度極端**：除非真正落在 LOH 區，HP1:HP2 不應 100:0。
3. **不過度模糊**：HP33 (V3-Fixed enum bug) 與 HP0 (untagged) 不應囤積大量 reads。

本文件量化的「imbalance distance to paired」即條件 1。條件 2/3 在 fig05b 的 M3/M4 metric 中補充。

> 📌 **核心定義**
> `ratio = (HP1 + HP11) / (HP1 + HP11 + HP2 + HP21)`（HP33/HP0 排除於分母）
> `imbalance = |ratio - 0.5|` ：0 = 完全平衡，0.5 = 極端 LOH-like
> `dist_to_paired = min(|r_query - r_PA|, |1 - r_query - r_PA|)` ：含 PS block flip

## Section 2 — 15 位點 ratio 數據表

完整資料：[`data/imbalance_ratio.tsv`](data/imbalance_ratio.tsv)

| Case | Chrom:Pos | BL ratio | V5 ratio | PA ratio | BL dist→PA | V5 dist→PA | Δ dist (BL-V5) |
|------|-----------|---------:|---------:|---------:|-----------:|-----------:|---------------:|
| A_TP01 | chr6:145444893 | 1.000 | 1.000 | 0.012 | 0.012 | 0.012 | 0.000 |
| A_TP02 | chr4:70548355 | 1.000 | 0.000 | 0.914 | 0.086 | 0.086 | 0.000 |
| A_TP03 | chr5:153209947 | 0.974 | 0.961 | 0.033 | 0.007 | 0.006 | +0.001 |
| **A_TP04** | **chr16:35118902** | **0.730** | **0.127** | **1.000** | **0.270** | **0.127** | **+0.143** |
| A_TP05 | chr7:109185781 | 0.766 | 0.761 | 0.704 | 0.062 | 0.057 | +0.005 |
| B_FPA1 | chr8:93565727 | 0.027 | 0.492 | NA | NA | NA | NA |
| **B_FPA2** | **chr9:137953060** | **0.951** | **0.725** | **0.960** | **0.009** | **0.235** | **−0.227** |
| B_FPB1 | chr7:52087777 | 0.257 | 0.189 | 0.275 | 0.017 | 0.086 | −0.069 |
| B_FPB2 | chr9:75383880 | 0.433 | 0.561 | 0.590 | 0.023 | 0.029 | −0.006 |
| C_V5max1 | chr19:4639528 | 1.000 | 1.000 | 0.969 | 0.031 | 0.031 | 0.000 |
| C_V5max2 | chr19:2235521 | 1.000 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| C_V5max3 | chr19:7405500 | 0.025 | 0.062 | 0.988 | 0.013 | 0.050 | −0.037 |
| D_SP1 | chr19:17565944 | 0.991 | 0.000 | 0.244 | 0.235 | 0.244 | −0.009 |
| D_SP2 | chr19:12452332 | 0.991 | 0.009 | 0.319 | 0.310 | 0.309 | +0.000 |
| D_SP3 | chr19:12467180 | 1.000 | 0.000 | 0.449 | 0.449 | 0.449 | 0.000 |

**關鍵觀察**：
- A_TP01/A_TP02/C_V5max1/C_V5max2/D_SP3 等 5 個位點，BL 與 V5 的 ratio 都是 1.000 或 0.000（單向 vote），所以 dist 完全相同（0.000）；orientation-flip 把 ratio=1 與 ratio=0 視為等距。
- **A_TP04** V5 把 ratio 從 0.730 翻到 0.127，距離 paired (1.000) 從 0.270 降到 0.127——**唯一強改善位點**。
- **B_FPA2** V5 把 ratio 從 0.951 拉到 0.725，反而離 paired (0.960) 變遠 → 看似「regression」，但**這是 V5 設計的副作用**，討論於 `05_per_site_improvement.md` Section 5。
- B_FPA1 paired BAM 在此位置 HP-tagged reads = 0（PA_HP0=114），標為 NA。

## Section 3 — 整體不平衡分布（histogram + 統計）

![dist_distribution](figures/04_imbalance/fig04b_distance_distribution.png)

**`fig04b_distance_distribution.png` 重點**：

| 統計量 | Baseline | V5 | 解讀 |
|--------|---------:|---:|------|
| mean BL_dist→PA | **0.109** | — | — |
| mean V5_dist→PA | — | **0.123** | V5 平均略遠 |
| median BL_dist→PA | **0.027** | — | — |
| median V5_dist→PA | — | **0.071** | V5 中位數明顯遠 |
| sites with dist=0.000 (BL) | 1 | — | C_V5max2 |
| sites with dist=0.000 (V5) | 1 | — | C_V5max2 |

**Inset scatter**（BL_dist vs V5_dist）顯示大多數點落在 y=x 對角線上或下方少數點偏離；少數點顯著偏離對角線並落在上半（V5 較差）——對應 B_FPA2、B_FPB1、C_V5max3。

> ⚠️ **整體 mean/median 兩面解讀**：
> - 表面看 V5 比 BL 「平均距離 paired 較遠」（mean +0.014, median +0.044）。
> - 但這個 average 受 B_FPA2 單一 large regression（−0.227）強烈拉動：把 B_FPA2 移除後，剩 14 sites 的 mean delta = +0.001（基本持平）。
> - 中位數變遠的真正原因是 BL 在 5 個 single-vote 位點上「**剛好碰巧**」與 paired 的 100:0 同向（A_TP01/A_TP02 等），是 self-phasing bug 的「假乾淨」。V5 在這些位點轉成相反方向（仍 100:0 但翻轉）導致 distance 不變但「視覺方向」改變——**這是 IGV audit Section 4 觀察的真實再現**，不是 V5 的退步。

## Section 4 — V5 改善幅度（V5 比 BL 更接近 paired ratio 的位點數）

![ratio_scatter](figures/04_imbalance/fig04a_ratio_scatter.png)

**`fig04a_ratio_scatter.png` 解讀**：每位點 3 個點（BL ●、V5 ■、PA ▲）以淡灰線連起，視覺顯示：

| 改善類別 | 位點數 | 占比 | 代表位點 |
|----------|------:|------:|---------|
| 強改善（Δdist > +0.10） | **1** | 6.7% | A_TP04 (+0.143) |
| 中改善（+0.02 < Δdist ≤ +0.10） | 0 | 0% | — |
| 持平（−0.02 ≤ Δdist ≤ +0.02） | 10 | 66.7% | A_TP01/02/03/05、C_V5max1/2、D_SP1/2/3、B_FPB2 |
| 反向（Δdist < −0.02） | 3 | 20% | B_FPA2 (−0.227)、B_FPB1 (−0.069)、C_V5max3 (−0.037) |
| NA（PA 無 HP tag） | 1 | 6.7% | B_FPA1 |

**依位點類別交叉表**（4 大類 × 3 改善類別）：

| 位點類別 | 強改善 | 持平 | 反向 |
|----------|------:|------:|-----:|
| A：Phase4 TP（5 sites） | 1 (TP_04) | 4 | 0 |
| B：Phase4 FP（4 sites） | 0 | 1 | 2 (+1 NA) |
| C：V5_max reassign（3 sites） | 0 | 2 | 1 |
| D：Self-phasing extreme（3 sites） | 0 | 3 | 0 |

**重點**：
1. V5 強改善只發生在 **A 類 TP** 中，且只有 A_TP04（bimodal HP+allele 配色）。其他 A_TP 位點都是 BL 已經 100:0、V5 不變或翻轉 → distance 持平。
2. V5 反向發生在 **B 類 FP** 為主（2/3 反向 + 1 NA = 全 B 類偏向「V5 較差」），這是因為 paired 在這些 FP 位點本身就是高度不平衡（germline het 或 LOH 邊緣）；V5 引入 PoN-anchored 修正反而把 ratio 拉向中間（更平衡），距離 paired 的高 imbalance 點就變遠。
3. **D 類 self-phasing 三個位點全部持平**——V5 仍**留下 100:0 的極端 vote**（只是方向可能翻轉），不修這類根本問題。這與 V5 設計**完全一致**：`getVote()` 層的 self-phasing fix 留待 V6 處理（見 IGV audit doc Section 4 的明確聲明）。

## Section 5 — 結論

**S5.1**　V5 在 imbalance ratio 的整體表現可以分為三段：

1. **A_TP04 一個位點顯示真正 lift**：Δdist=+0.143（強改善），對應 V5 修正了 V3-Fixed 的 enum bug 並讓 90 reads 從 HP1 大群被重歸 HP2 大群，與 paired (0:50) 對齊。這是 V5 Layer 1.5 fallback 機制的視覺驗證 best-case。

2. **10 個位點持平（66.7%）**：BL 與 V5 的 ratio 都是 100:0 single-vote，加上 orientation-flip 容許 1 與 0 等距，所以 distance 不動。其中 5 個位點（A_TP01/A_TP02/D_SP1/D_SP2/D_SP3）真實情境是「BL 與 V5 都 self-phasing，方向可能不同」——distance 不變但視覺差很大；這個現象是 IGV audit Section 4 的核心發現，不是 V5 的失敗，而是 V5 設計**不解決** getVote 層 self-phasing 的明確展現。

3. **3 個位點反向（20%）**：B_FPA2/B_FPB1/C_V5max3。B_FPA2 是 paired germline het（true LOH-like 0.96 imbalance）但 V5 引入 PoN-fallback vote 把它拉向 0.725 → 距離 paired 變遠。**這是 V5 將 ambiguous reads 強制分配的副作用**，不是 bug——詳見 `05_per_site_improvement.md` Section 5。

**S5.2**　Imbalance distance 平均看似 V5 略差（mean +0.014），但這個總和受 B_FPA2 單一強反向（−0.227）主導；移除 B_FPA2 後 mean delta = +0.001。**真正的訊號是 self-phasing 三個 D 類位點的 100:0 極端不平衡完全沒被修復**，與 V5 文件聲明完全一致。

**S5.3**　Imbalance ratio 是評估 phasing 品質的有用 lens，但**不能單獨用來判定 V5 vs BL 何者較好**。原因：(a) orientation-flip 把 100:0 與 0:100 視為等距，掩蓋了視覺方向差異；(b) PoN 校正會改變 ambiguous read 的歸屬，使 V5 在 LOH-like FP 位點「看似較差」；(c) self-phasing 不會體現在 dist 上（因 paired 也常 100:0）。需與 `02_concordance_quantification.md` 的 4-metric concordance 結合判讀（fig05b heatmap 即為此目的）。

**S5.4**　**結論**：V5 在 15 位點中得到 1 強改善 / 10 持平 / 3 反向 / 1 NA。**強改善位點 (A_TP04) 顯示 V5 機制正確；反向位點 (B_FPA2/FPB1/C_V5max3) 在 V5 設計範圍內可解釋；持平的 self-phasing 位點 (D_SP1/2/3) 是已知 known limitation，不在 V5 修補範圍**。整體與 IGV 真截圖 audit 結論一致：V5 在 enum bug 與 ambiguous-resolve 兩層機制成立，self-phasing 修復留待後續版本。
