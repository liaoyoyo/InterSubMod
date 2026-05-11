---
title: 13 位點細緻 Nuance 審計 — V5 改善的真實樣貌與 metric 衝突
date: 2026-04-29
author: liaoyoyo2001
tags: [audit, per-site, v5, paired-concordance, nuance, methodology-warning]
status: validated_complete
audience: PI + reviewers
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/02_read_intersection_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/05_per_site_improvement.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md
---

# 13 位點細緻 Nuance 審計 — V5 改善的真實樣貌

## §0 一句話結論

> **V5 對 13 個 audit sites 中只有 1 個（A_TP04）是「強改善」**；其餘是 metric 衝突 / behavior-only / NA / regression / tie。所謂「Clean PS paired concordance +13.3pp」是 **conditional accuracy**（V5 evaluable denominator 較小），用 fixed common-read denominator 重算則 baseline 50.77% vs V5 50.70% **幾乎持平**。本檔對齊用戶細緻觀察，揭露 V5 改善的真實樣貌與 metric 解讀注意事項。

---

## §1 用戶提出的細緻 nuance 表（完整保留）

### 1.1 13 位點分類觀察

| 位點 | 類型 | V5 vs BL 觀察 | 主要問題 | 注意點 |
|------|------|--------------|---------|-------|
| **A_TP04** | TP / best improve | exact-site **+0.1427**；paired conditional **+0.9737** | **唯一強改善案例** | 可作為 V5 成功例子，但 **L4 仍略輸 BL**，代表改善主要來自 directional tag 對齊，不是所有 orientation-corrected metric 都贏 |
| A_TP05 | TP | exact-site 小幅 **+0.0048**，但 paired conditional **−0.1176** | V5 evaluable reads 大幅下降：BL 107、V5 68 | **不可稱改善**；若用 fixed common-read denominator，V5 明顯變差 |
| A_TP02 | TP / problem PS | exact-site neutral，但 paired conditional **−0.8490** | PS class = problem；V5 orientation 對 paired 很差 | 這類**不應當成 V5 失敗核心證據**，但也不能拿來支持改善 |
| B_FPA1 | FP / low germline | paired ratio **NA**；paired concordance **NA** | Paired site 沒有可用 HP1/HP2 assigned reads | **不可判斷** V5 是否更接近 paired；應從改善統計排除 |
| B_FPA2 | FP | exact-site **−0.2266**，但 paired conditional **+0.7565** | metric 方向衝突 | **方法學警訊**：conditional read concordance 說 V5 好，exact-site ratio / L4 說 V5 差；不可單獨引用 |
| **B_FPB1** | FP | exact-site **−0.0685**；paired conditional **−0.0667**；L4 也輸 | 一致 regression | **可列為 V5 真正變差位點** |
| B_FPB2 | FP | exact-site neutral/slightly worse **−0.0060**，paired conditional **+0.8571**，L4 輸 | conditional metric 與 L4 衝突 | V5 對 directional reads 對齊較好，但 orientation-corrected family consistency 不支持 |
| C_V5max1 | V5max | sanity 顯示 HP33→HP11 很多，但 exact-site / paired 都 **tie** | reassign 發生在 ±50bp window，不一定落在 exact site | **不能把「V5 有 reassign」直接等同「該 site tag 更接近 paired」** |
| C_V5max2 | V5max | 同 C_V5max1，paired/exact-site **tie** | HP33→HP11 是 **行為證據，不是 correctness 證據** | 可支持 V5 fallback 有作用，但**不支持 paired improvement** |
| **C_V5max3** | V5max | exact-site **−0.0370**；paired **−0.0299**；L4 輸 | V5 reassign 後**反而偏離 paired** | 這是 V5max 類別中**最該標註的 regression** |
| D_SP1 | self-phasing / low germline | exact-site neutral/slightly worse；paired **−0.5337** | low_germ_n，paired orientation 不穩 | 不適合作強結論；但顯示 self-phasing extreme 下 V5 不一定改善 |
| D_SP2 | self-phasing / problem PS | exact-site neutral；paired **+0.3562**；L4 小贏 | problem PS，訊號不穩 | 可列為 conditional improvement，但**可信度低於 clean PS** |
| D_SP3 | self-phasing | exact-site neutral；paired **−0.0492**；L4 小贏 | metric 輕微衝突 | 不應過度解讀，**接近持平** |

### 1.2 較不特殊但可作 baseline 參考

| 位點 | 觀察 | 注意點 |
|------|------|-------|
| A_TP01 | exact-site、paired、L4 幾乎 **tie** | V5 沒有實質影響 |
| A_TP03 | exact-site 小幅中性；paired / L4 略輸 | **不支持改善，但幅度小** |

### 1.3 跨位點共同注意點（methodology warnings）

| # | 注意點 | 影響 |
|:-:|-------|------|
| 1 | `paired_ground_truth.tsv` 的 +6.65pp / clean +13.24pp 是 **conditional accuracy**，只在 HP1/HP2/HP11/HP21 directional reads 上計算 | V5 evaluable denominator 較小 |
| 2 | 若用 **fixed common reads** 當 denominator，全體幾乎持平：BL 657/1294 = **50.77%**，V5 656/1294 = **50.70%** | conditional vs fixed 差距巨大 |
| 3 | 05 exact-site ratio 只支持 **1 個強改善、3 個 regression、10 個 neutral、1 個 NA** | 強改善案例稀少 |
| 4 | V5max 類位點要分清楚：HP33→directional 是 **sanity / behavior improvement**，不必然是 **paired correctness improvement** | 不能 over-claim |

---

## §2 V5 改善的細緻分類（更新版）

### 2.1 改善等級分類

| 等級 | 定義 | 數量 | 位點 |
|------|------|:--:|------|
| 🟢 **Strong improve** | exact-site + paired conditional + L4 都正向 | **1/15** | A_TP04 |
| 🟡 **Conditional improve** | paired conditional 正向但 L4 / exact-site 不一致 | 3/15 | B_FPA2, B_FPB2, D_SP2 |
| ⚪ **Tie / neutral** | 三 metric 都接近 0 或方向矛盾 | 5/15 | A_TP01, A_TP03, C_V5max1, C_V5max2, D_SP3 |
| 🔴 **Regression** | exact-site 與 paired 都負向 | **2/15** | B_FPB1, C_V5max3 |
| ⚫ **Cannot judge** | NA 或 low_germ_n 不穩 | 4/15 | B_FPA1, A_TP02, A_TP05, D_SP1 |

→ **strong improve 比例只有 6.7%（1/15）**，不能用「+13.3pp」單一數字 over-claim V5 全面改善。

### 2.2 改善的真實 mapping

| 既有 audit suite 結論 | 真實樣貌（13 位點 nuance 後） |
|---------------------|---------------------------|
| 「Aggregate paired concordance +6.65pp」 | conditional accuracy；fixed-denom 後幾乎持平 |
| 「Clean PS paired concordance +13.3pp」 | 同上 conditional；clean PS subset 內 V5 略勝但 evaluable reads 較少 |
| 「V5max1-3 是 V5 修復強的 sites」 | 1/3 是 tie (V5max1, V5max2)，1/3 是 regression (V5max3)；HP33→HP11 是 behavior 不是 correctness |
| 「SP1-3 是 self-phasing extreme 修復」 | 1/3 conditional improve (SP2)，2/3 不可信 (SP1 low_germ, SP3 衝突) |

---

## §3 13 位點 IGV 截圖（4 版本對比：BL / V2b / V3F / V5）

### 3.1 既有 IGV 4ver 截圖（baseline + V2b + V3F + V5 並列，HP coloring）

#### A 類：True Positive (5 sites)

| 位點 | 截圖 | 改善等級 |
|------|:----:|:------:|
| A_TP01 chr6:145444893 | ![A_TP01](../figures/igv_v5_audit/by_HP_4ver/A_TP01_chr6_145444893.png) | ⚪ tie |
| A_TP02 chr4:70548355 | ![A_TP02](../figures/igv_v5_audit/by_HP_4ver/A_TP02_chr4_70548355.png) | ⚫ problem PS, 不可判斷 |
| A_TP03 chr5:153209947 | ![A_TP03](../figures/igv_v5_audit/by_HP_4ver/A_TP03_chr5_153209947.png) | ⚪ neutral 略輸 |
| **A_TP04 chr16:35118902** | ![A_TP04](../figures/igv_v5_audit/by_HP_4ver/A_TP04_chr16_35118902.png) | 🟢 **唯一強改善** |
| A_TP05 chr7:109185781 | ![A_TP05](../figures/igv_v5_audit/by_HP_4ver/A_TP05_chr7_109185781.png) | ⚫ V5 evaluable 下降 |

**A_TP04 的解讀**：V5 改善主要來自 directional tag 對齊（exact-site +0.1427, paired +0.9737），但 L4 仍略輸 BL → 改善是「partial」而非「all-axis」。

#### B 類：False Positive (4 sites)

| 位點 | 截圖 | 改善等級 |
|------|:----:|:------:|
| B_FPA1 chr8:93565727 | ![B_FPA1](../figures/igv_v5_audit/by_HP_4ver/B_FPA1_chr8_93565727.png) | ⚫ NA (paired 無 HP1/HP2) |
| B_FPA2 chr9:137953060 | ![B_FPA2](../figures/igv_v5_audit/by_HP_4ver/B_FPA2_chr9_137953060.png) | 🟡 metric 衝突 |
| **B_FPB1 chr7:52087777** | ![B_FPB1](../figures/igv_v5_audit/by_HP_4ver/B_FPB1_chr7_52087777.png) | 🔴 **真 regression** |
| B_FPB2 chr9:75383880 | ![B_FPB2](../figures/igv_v5_audit/by_HP_4ver/B_FPB2_chr9_75383880.png) | 🟡 metric 衝突 |

#### C 類：V5max (3 sites, V5 修復幅度最大)

| 位點 | 截圖 | 改善等級 |
|------|:----:|:------:|
| C_V5max1 chr19:4639528 | ![C_V5max1](../figures/igv_v5_audit/by_HP_4ver/C_V5max1_chr19_4639528.png) | ⚪ tie (HP33→HP11 是 behavior) |
| C_V5max2 chr19:2235521 | ![C_V5max2](../figures/igv_v5_audit/by_HP_4ver/C_V5max2_chr19_2235521.png) | ⚪ tie |
| **C_V5max3 chr19:7405500** | ![C_V5max3](../figures/igv_v5_audit/by_HP_4ver/C_V5max3_chr19_7405500.png) | 🔴 **regression** (V5 偏離 paired) |

**C 類重要警語**：「V5max」這個分類名稱**容易誤導**。它原本意思是「V5 vs BL 在 sanity check 上 reassign 數量最大」，**不是「V5 結果最好的位點」**。實際上：
- V5max1, V5max2：reassign 多但 paired tie
- V5max3：reassign 多但 paired **regress**

#### D 類：Self-phasing extreme (3 sites, chr19 cluster)

| 位點 | 截圖 | 改善等級 |
|------|:----:|:------:|
| D_SP1 chr19:17565944 | ![D_SP1](../figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png) | ⚫ low_germ, 不穩 |
| **D_SP2 chr19:12452332** | ![D_SP2](../figures/igv_v5_audit/by_HP_4ver/D_SP2_chr19_12452332.png) | 🟡 conditional improve, 可信度低 |
| D_SP3 chr19:12467180 | ![D_SP3](../figures/igv_v5_audit/by_HP_4ver/D_SP3_chr19_12467180.png) | ⚪ 接近持平 |

**D 類重要警語**：SP1-3 在 0.93 baseline 下都是 self-phasing extreme（HP1:HP2 ≥ 100:1），但 paired ground truth 在這些 site 上**信號本身不穩**（low_germ_n 或 problem PS），所以 V5 在這裡的改善判斷**可信度有限**。

### 3.2.1 6-Panel IGV with Paired Ground Truth (2026-04-29 補充)

新增 6-panel IGV 截圖（含 PA_normal + PA_tumor + 8 LOH/GE BED tracks），位於 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/09_purity06/igv_with_paired_loh/`：

| Site | 6-panel 截圖 | 對應 nuance 觀察 |
|------|:----:|------|
| **A_TP04** | ![TP04](figures/09_purity06/igv_with_paired_loh/TP04_chr16_35118902.png) | 🟢 唯一強改善 — paired 與 V5_93 同方向 |
| **D_SP1** | ![SP1](figures/09_purity06/igv_with_paired_loh/SP1_chr19_17565944.png) | ⚫ paired 為 HP1 但 V5 翻 HP2（paired -0.5337）|
| **D_SP2** | ![SP2](figures/09_purity06/igv_with_paired_loh/SP2_chr19_12452332.png) | 🟡 conditional improve, paired 不穩 |
| **D_SP3** | ![SP3](figures/09_purity06/igv_with_paired_loh/SP3_chr19_12467180.png) | ⚪ 接近持平 |
| **C_V5max1** | ![V5max1](figures/09_purity06/igv_with_paired_loh/V5max1_chr19_4639528.png) | ⚪ tie (HP33 reassign 但 paired 同向) |
| **C_V5max2** | ![V5max2](figures/09_purity06/igv_with_paired_loh/V5max2_chr19_2235521.png) | ⚪ tie |
| **C_V5max3** | ![V5max3](figures/09_purity06/igv_with_paired_loh/V5max3_chr19_7405500.png) | 🔴 V5 reassign 後反 paired (regression 視覺證據) |

→ 對 SP/V5max class 的觀察：**V5 在 self-phasing extreme sites 上翻轉 PS block orientation 反而偏離 paired truth**，這是 19 號 §2.1 reggression 等級分類的視覺證據。

### 3.2 補充：0.6 vs 0.93 V5 行為對比（從 09 號 simulation）

對 9 個 sites（SP1-3, V5max1-3, TP01/02/04），新增 IGV 截圖 (4 panels: BL_93/V5_93/BL_06/V5_06):

| 位點 | 0.6 場景截圖 |
|------|:--:|
| SP1 | ![purity 06 SP1](figures/09_purity06/igv/SP1_chr19_17565944.png) |
| V5max1 | ![purity 06 V5max1](figures/09_purity06/igv/V5max1_chr19_4639528.png) |
| TP04 | ![purity 06 TP04](figures/09_purity06/igv/TP04_chr16_35118902.png) |

→ 0.6 下 baseline self-phasing 衰減（SP1 從 113:1 → 1:70）；V5 與 baseline 在 0.6 趨於一致。

---

## §4 統計層級驗證（conditional vs fixed denominator）

### 4.1 既有 +6.65pp / +13.3pp 的 conditional 限制

```
Aggregate paired concordance:
  BL:  72.20%  (eligible reads ~)
  V5:  78.85%  (eligible reads ~ 較少)
  Δ = +6.65pp  ← 但這是 conditional 不是 fixed-denom
```

→ V5 evaluable denominator 較小（部分 reads 進 HP33 ambiguous）→ conditional accuracy 自然偏高。

### 4.2 用 fixed common-read denominator 重算

```
Common reads (both BL and V5 直接 tagged):
  BL:  657 / 1294 = 50.77%
  V5:  656 / 1294 = 50.70%
  Δ = -0.07pp  ← 幾乎持平
```

→ 在「兩版本都 tag 的 reads」上，V5 vs BL 幾乎沒差別。

### 4.3 兩個 metric 的價值差異

| Metric | 計算方式 | 價值 |
|--------|--------|------|
| **Conditional accuracy** | 只看雙方都 tag 為 directional 的 reads | 反映「在 V5 願意給 directional tag 時，是否更精準」 |
| **Fixed common-read** | 雙方都 tag 的 reads（含 V5 conservative 情況）| 反映「全 read 池的整體準確率」 |
| **Strict total** | 包含未 tag 的 reads (HP:i:0) | 反映「coverage 完整性」 |

→ V5 的 conservative tagging（HP33 多）在 conditional accuracy 看起來好，但 **fixed-denom 顯示它沒有真正提升整體 accuracy**。

---

## §5 Methodology Warnings（解讀準則）

### 5.1 不可 over-claim 的 4 個訊號

| 警訊 | 表現 | 範例 |
|------|------|------|
| **#1 Conditional vs Fixed denominator 差距大** | conditional +13pp 但 fixed -0.07pp | aggregate metric |
| **#2 metric 方向衝突** | exact-site 與 paired 反向 | B_FPA2, B_FPB2 |
| **#3 evaluable denominator 大幅下降** | V5 reads 比 BL 少 30%+ | A_TP05 |
| **#4 reassign 動作不等於 correctness** | HP33→HP11 量大但 paired tie | C_V5max1, C_V5max2 |

### 5.2 該怎麼正確描述 V5

❌ **Over-claim**：「V5 在所有 sites 上都更接近 paired」
❌ **Over-claim**：「Clean PS +13.3pp 證明 V5 全面改善」
✅ **正確**：「V5 在 1/15 sites 強改善（A_TP04）；3/15 conditional improve；2/15 regression；其餘 tie 或不可判斷。整體 conditional metric 偏好 V5（+13pp），但 fixed-denom 重算後幾乎持平。」

### 5.3 V5 真實價值定位

V5 的價值**不在「全面改善 paired concordance」**，而在：

| 真實價值 | 證據 |
|---------|------|
| ✅ 修復 17.3:1 self-phasing bias (genome-wide aggregate) | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md` |
| ✅ Conservative tagging（HP33 顯式標出 ambiguous） | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` §2.2 |
| ✅ 0.6 purity 下 self-phasing 自然衰減後 V5 仍合理 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` |
| ✅ 不傷 ClairS-TO calling F1 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` |
| ⚠ Per-site paired concordance 提升幅度有限（1 強改善 / 2 regression）| **本檔** |

---

## §6 對既有 audit suite 結論的影響

### 6.1 結論不變但需 nuance 化

| 既有結論 | 是否仍成立 | Nuance 化 |
|---------|:--:|---------|
| V5 修復 17.3:1 bias | ✅ | aggregate 層級正確 |
| V5 sanity check 全 PASS | ✅ | 守恆律不變 |
| Aggregate paired +6.65pp | ⚠ | 是 conditional accuracy；fixed-denom 持平 |
| Clean PS +13.3pp | ⚠ | 同上，clean PS subset 內 conditional 改善 |
| V5max sites 是 V5 修復最強 | ❌ | 名稱誤導；實際 1 regression 2 tie |
| V5 全面更接近 paired | ❌ | 1 強改善 / 2 regression / 5 tie / 4 不可判斷 / 3 conditional |

### 6.2 該更新的章節

| 報告 | 章節 | 建議更新 |
|------|------|---------|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/05_per_site_improvement.md` | 改善排序 | 加入「強改善 1/15」明確標示 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` | +6.65pp / +13.3pp | 加 caveat：「conditional accuracy, denominator 較小」 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/08_synthesis_conclusions.md` | Q1 V5 在哪些 metric 比 BL 接近 | nuance 化呈現 |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/14_user_report_integrated.md` | §9 強數據證明 | 在「+13.3pp」標註本檔 caveat |

---

## §7 一頁速查

```
13 audit sites × V5 改善等級分類:
  🟢 Strong improve:    1 (A_TP04)
  🟡 Conditional only:  3 (B_FPA2, B_FPB2, D_SP2)
  ⚪ Tie/neutral:       5 (A_TP01, A_TP03, V5max1, V5max2, SP3)
  🔴 Regression:        2 (B_FPB1, V5max3)
  ⚫ Cannot judge:      4 (B_FPA1, A_TP02, A_TP05, D_SP1)
                       ────────
                       15 total

Methodology:
  +13.3pp = conditional accuracy (V5 evaluable denom 較小)
  Fixed common-read denom: BL 50.77% vs V5 50.70% (持平)

V5 真實價值:
  ✓ Aggregate self-phasing 修復
  ✓ Conservative tagging (HP33 顯式標)
  ⚠ Per-site paired correctness 改善有限
  ❌ 不可 claim「全面改善」
```

---

## §8 圖檔索引

### 4 版本 IGV 截圖（baseline + V2b + V3F + V5）

| 路徑 | 內容 |
|------|------|
| `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_HP_4ver/` | 13 位點 4 版本對比 |
| `InterSubMod/docs/reports/pi_reports/2026/04/figures/igv_v5_audit/by_MOD_4ver/` | 同上但 methylation coloring |

### 0.6 vs 0.93 IGV 截圖

| 路徑 | 內容 |
|------|------|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/09_purity06/igv/` | 9 位點 4 場景 (BL_93/V5_93/BL_06/V5_06) |

---

## §9 跨檔索引

| # | 文件 | 與本檔關係 |
|---|------|---------|
| 02 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/02_read_intersection_concordance.md` | L1-L4 metric 定義 |
| 05 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/05_per_site_improvement.md` | 15 位點數值排序 |
| 07 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` | +6.65pp / +13.3pp 出處 |
| 09 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` | 0.6 場景補充 |
| 16 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/16_baseline_subgenotype_clarification.md` | baseline GT2/GT3 機制 |
| 18 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/18_purity_calculator_failure_root_cause.md` | Purity calculator 修復 |
