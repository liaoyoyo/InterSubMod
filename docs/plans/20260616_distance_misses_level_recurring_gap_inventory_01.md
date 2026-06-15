---
title: "方法學系統性盲點盤點 — distance/clustering 漏 mean-level shift（遞迴主題）"
date: 2026-06-16
type: inventory + 根因 + 修法（方法-觀察不一致的根因盤點）
trigger: 用戶觀察 within-HP 結構比方法(Stage② pattern 1.2%)多；確認方法漏什麼
data_sources: /tmp/dbeta_chr1_wh (within-HP 解構), within_hp_clean_pilot, AB_TWO_DIRECTION_ANALYSIS, structure_false_negative plan
claim_levels: 解構數據=L1(chr1 per-region) / 遞迴主題=L1(三處實證) / 修法=L2
status: 盤點完成 — 待修法決策
---

# 方法學系統性盲點：distance/clustering 漏 mean-level shift

## 1. 核心發現（方法 ≠ 觀察的根因）

ISM 的結構偵測**以 distance/clustering 為主**（HP-AUC、PERMANOVA、cluster×label、within-HP distance）。distance（BERNOULLI/NHD）= **per-CpG pattern 距離**。
**但真實甲基結構常是「水平(level)位移」**（一群 reads 整體甲基較高/較低，per-CpG pattern 相似）→ **distance 不分開 level 位移 → 系統性漏掉**。

## 2. 遞迴實證（同一盲點，三處出現）

| # | 場景 | distance/clustering 說 | level(mean β) 說 | 漏掉 |
|---|------|------|------|------|
| **A/B 兩方向** | clustering/GlobalTest | A germline Δβ 25% sig | B over-call/漏 mean-shift | A/B 互補 ~10% 不一致 |
| **最終判別 VC** | cluster×label(GlobalTest) 錨定 | Δβ 不參與 verdict | **41% 假陰性**(已用 LabelShift 修) |
| **within-HP (Stage②)** | distance clustering **1.2%** | **level 乾淨平衡 26.1%** | **level-only 25.6% 完全漏** ⭐ |

**within-HP 解構(chr1 2624)**: pattern(distance) 31(1.2%) / level(乾淨平衡) 684(26.1%) / both 12 / **level-only(漏) 672(25.6%)** / pattern-only 19 / neither 1921。
→ within-HP 結構 **96% 是 level，distance 只抓到 4%**。

## 3. 是 gap 還是真實？— 是 GAP

那 26.1% level 結構**真實**（判準：min群≥20%平衡 + gap>0.15 + var-reduction>0.5 = 非 tiny-group artifact）。distance metric(pattern-based) **設計上抓不到均勻 level 位移**。**非真實罕見，是方法盲點**。

## 4. 驗證方式（如何確認 method≠observation）

- **解構量化**(已做): level vs pattern 各抓什麼 → level-only 25.6% = 漏的。
- **個案圖**(待): 取 level-only region 的 HP1-family 甲基熱圖 → 肉眼確認有水平分群但 distance NGroups=1 + silhouette<0.5。
- **silhouette 檢視**: level-only region 的 within-HP silhouette 應 <0.5（distance 不分），但 mean β gap>0.15（level 分）。

## 5. 修法（within-HP 加 level 軸）

- **Stage② 增強**：`within_hp_clean_multigroup = pattern(distance, 現有) OR level(mean β 乾淨平衡, 新增)`。
  - level 判準：tumor HP-family per-read mean β 最佳 2-split，min群≥max(3,20%) + gap>0.15 + var-reduction>0.5。
  - 預期 within-HP 1.2%→~26% → MultiGroupNoLabel/CN-Substructure 大增 → 救回 35%→更高。
- **更廣意涵**：ISM 結構判別應**並用 level + pattern 兩軸**（distance 抓 pattern，mean β/Δβ 抓 level）。已在 VC(LabelShift) + Δβ 模組部分落實；within-HP 是最後一塊。

## 6. 決策
- D1：§1-§3 根因（distance 漏 level，遞迴盲點）+ 解構數據(level-only 25.6%)符合你的觀察嗎？
- D2：Stage② 加 level 軸（抓回 ~26%）— 實作？(改 VC 須再更新 golden + 全基因組重驗)
- D3：先產 level-only 個案圖肉眼確認再修，還是直接修？
