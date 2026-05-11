---
title: V5 PPT Gallery — 所有圖示索引與美觀檢核
date: 2026-04-30
author: liaoyoyo2001
tags: [pi-presentation, gallery, figures-index, visual-review]
status: review_in_progress
audience: 用戶確認圖示美觀與正確性
purpose: 用 .md 相對路徑顯示所有 PPT 用圖，方便逐張檢視
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/22_pi_presentation_integrated_narrative.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/23_ppt_slide_by_slide_outline.md
---

# V5 PPT 用圖示 Gallery

## 設計原則（OFFICE_LIGHT 主題）

對齊 longphase-to / longphase-s 論文 slide 風格：
- 白底 + STRUCT_BLUE (#1F3A93) 主敘事
- HP1 紅 (#E91E63 IGV) / HP2 藍 (#1976D2 IGV)
- Somatic 橘 (#E67E22)
- HP33 ambiguous 紫 (#7B1FA2)

---

# § Part 1：本 session 新生成的 6 張綜合圖

## fig23a — V5 修法時間軸

**用途**：Slide 3（流程概覽）+ Slide 22（整體時間軸）

![fig23a 修法時間軸](figures/23_ppt_narrative/fig23a_version_timeline.png)

**內容**：
- 5 個版本節點（V0 → V2b → V3F → V5 → V5+）
- 每個版本有 commit hash + 日期 + 修法 issue 列表
- 紅色 = bug 階段；綠色 = fixed

**檢核**：
- [ ] 5 節點清楚 ✓
- [ ] commit hash + 日期可讀 ✓
- [ ] issue 列表簡潔 ✓
- [ ] 顏色對應狀態（紅 bug / 黃 partial / 綠 fixed）✓

---

## fig23b — Metric Map (4-Layer × 3 Versions × 1 Truth)

**用途**：Slide 1-2（核心結論支援）+ Slide 19（業界對齊）

![fig23b Metric Map](figures/23_ppt_narrative/fig23b_metric_map.png)

**內容**：
- 4 layers: Caller / Phase / Tag / Downstream
- 5 columns: Layer / Metric / Baseline / V5 / Paired Truth
- 右側 verdict 標示

**檢核**：
- [ ] 4 layer 分明 ✓
- [ ] 數值 alignment ✓
- [ ] 顏色編碼（紅 bug / 綠 ✓ / 橘 truth）✓
- [ ] verdict 結論清楚 ✓

**用途**：解釋「不同 metric 的 denominator 不同」，避免簡報時混淆 caller F1 與 tag concordance。

---

## fig23c — 研究發展樹（V5 → 五大目標）

**用途**：Slide 9（五大目標）+ Slide 21（研究發展樹）

![fig23c 研究發展樹](figures/23_ppt_narrative/fig23c_research_tree.png)

**內容**：
- 根節點：V5 / V5+ 修復
- 5 大分支：Subclonal / ASM / TP-FP / Cross-platform / Clinical
- 每分支下 3 條 sub-task
- 底部：3 層 P0/P1/P2 priority

**檢核**：
- [ ] 5 大目標彩色分明 ✓
- [ ] 連接線清楚 ✓
- [ ] sub-task 可讀 ✓
- [ ] P0/P1/P2 priority 顯眼 ✓

---

## fig23d — Bug 雙層機制 vs V5 雙層修復（左右對照）

**用途**：Slide 6（17.3:1 機制）+ Slide 13（PON-only flag）

![fig23d 雙層 bug 機制](figures/23_ppt_narrative/fig23d_dual_layer_bug.png)

**內容**：
- 左：Baseline 紅色（Phase 階段 anchor bug + Tag 階段 voting bug → 17.3:1）
- 右：V5 綠色（PON-only flag + 三層 getVote → 1.42:1）
- Phase + Tag 兩層各自說明 + 結果對比

**檢核**：
- [ ] 左右對照清楚 ✓
- [ ] 雙層機制分明 ✓
- [ ] 結果框醒目 ✓
- [ ] Baseline 紅 / V5 綠 配色一致 ✓

---

## fig23e — 數據對比視覺（4 metric ratio + distance to truth）

**用途**：Slide 6（17.3:1 量化）+ Slide 15（WG 證據）

![fig23e 數據對比](figures/23_ppt_narrative/fig23e_data_comparison.png)

**內容**：
- 左 panel：4 metric ratio bar（17.3 / 2.08 / 1.42 / 1.08）含 log scale
- 右 panel：4 場景 distance to PA + 改善 −47% / 偏離 +192% 標註

**檢核**：
- [ ] 4 個 ratio 對比清楚 ✓
- [ ] log scale 標示 ✓
- [ ] 改善 −47% 箭頭醒目 ✓
- [ ] 偏離 +192% 警示清楚 ✓

---

## fig23f — PPT 章節結構圖（6 sections × 22 slides）

**用途**：Slide 1（簡報架構）或附錄

![fig23f PPT 結構](figures/23_ppt_narrative/fig23f_ppt_structure.png)

**內容**：
- 6 大圓圈節點（顏色對應章節）
- 每節點下 slide 細項列表
- 底部主要圖示對應索引

**檢核**：
- [ ] 6 sections 顏色分明 ✓
- [ ] slide 細項可讀 ✓
- [ ] 底部圖示索引清楚 ✓

---

## fig23g — Somatic HP Bias 精確版 v2（取代舊 fig01d）

**用途**：Slide 6（17.3:1 機制詳解 — 雙層 bug + 雙層修復）

![fig23g v2 精確版 somatic bias](figures/23_ppt_narrative/fig23g_somatic_bias_precise.png)

### v2 修正的 4 個文字精確性問題

| Panel | v1 問題 | **v2 修正** |
|-------|------|------|
| [A] | 「同 tumor clone reads 共享 ALT pattern」過於籠統 | **同一 sub-clonal population reads 在多個 somatic sites 共現 ALT 等位**（精確生物學描述）+ 補充「這是 tumor sub-population 真實特性，不是 phasing 工具該模仿的訊號」 |
| [B] | `convertNonGermlineToSomatic()` 解釋抽象 | 明確說 **`origin: ORIGIN_UNDEFINED → SOMATIC`** + **`addEdge() 跳過 SOMATIC origin variants`**（具體程式碼層次）|
| [C]+[D] | 用 (HP1=5, HP2=3, HP1_1=1, HP2_1=0) 例子 → baseline 與 V5 結果**都是 hpResult=11**，無法 differentiate | 改用 **(HP1=4, HP2=6, HP1_1=2, HP2_1=0)** → baseline=11, V5=21（兩者結果**不同**才能展示 V5 修法效果） |
| [E] | distance 計算缺公式說明 | 補上 **`distance = |log10(ratio_X) − log10(ratio_PA)|`** |

### 新版的 differentiating 例子（[C][D] 核心）

```
countMap: HP1=4, HP2=6, HP1_1=2, HP2_1=0 (一條 read 同時涵蓋 germline 與 somatic)

Baseline (somatic-first):
  Priority 1: (HP1_1, HP2_1) → HP1_1=2 > HP2_1=0 → break
  → hpResult = 11 (HP1 family) ❌
  ※ 但 germline HP2=6 > HP1=4，read 真實應屬 HP2 family

V5 (germline-first):
  Layer 1: HP2=6 > HP1=4 → germlineResult = 2
  Layer 1.5: 跳過
  Layer 2: somaticTotal=2>0 → hpResult = 21 (HP2 family) ✓

⚡ Baseline=11 vs V5=21 — 不同結果！這就是 17.3:1 修復來源
```

### 5 個 panels 完整結構

- **[A] Baseline Phase**：6 variants + 3 reads（同 sub-clone）+ 紅色 graph edges 連結 4 somatic → 推進同一 phase block (HP1)
- **[B] V5 Phase**：somatic 變空心圈（excluded from anchor）+ `convertNonGermlineToSomatic()` 細節 + 只有 g1 ↔ g2 anchor
- **[C] Baseline Tag**：countMap bar (4/6/2/0/0) + Priority for-loop 3 條 + 觸發 [1] → hpResult=11
- **[D] V5 Tag**：同 countMap + 三層決策 active 標示（Layer 1 active, Layer 1.5 skipped, Layer 2 active）+ trace=21
- **[E] 結果**：Baseline 2.08:1 / V5 1.42:1 / Paired 1:1.08 三柱對比 + distance 公式

**檢核（v2 全通過）**：
- [x] 雙層 bug + 雙層修復清楚分明 ✓
- [x] 生物學精確描述（sub-clonal population）✓
- [x] convertNonGermlineToSomatic 程式碼層次 ✓
- [x] differentiating 例子展示 V5 真實效果（baseline=11 vs V5=21）✓
- [x] V5 三層 active/skipped 視覺標示 ✓
- [x] Paired truth + distance 公式對比明顯 ✓

→ **建議用此 v2 取代 fig01d 在 PPT slide 6**。

---

# § Part 2：既有 audit suite 圖（重用於 PPT）

## (1) 17.3:1 機制總圖

**Slide 6 用**

![Self-phasing 機制](figures/01_code_diff/fig01d_somatic_bias_explanation.png)

**內容**：6-panel 含 Expected / Observed / Mechanism / Aggregate / Per-site / V5-fix

---

## (2) Phase vs Tag 演算法流程

**Slide 5 用**

![Phase 演算法](figures/13_phase_vs_tag_algo/figA_phase_algorithm_flow.png)

---

## (3) getVote 三層決策樹

**Slide 11 用**

![getVote 三層](figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png)

---

## (4) 具體例子（3 reads × 5 variants）

**Slide 12 用**

![具體例子](figures/13_phase_vs_tag_algo/figC_concrete_example.png)

---

## (5) GT 分布（baseline vs V2b 對比）

**Slide 16 用**

![GT 分布](figures/12_gt_distribution/figA_gt_class_by_filter.png)

---

## (6) PON / Phase / Tag 三階段對照

**Slide 5 補充**

![PON Phase Tag](figures/01_code_diff/fig01c_pon_phase_tag_comparison.png)

---

## (7) 影響分類矩陣

**Slide 14 用**

![影響分類](figures/14_impact_classification/fig14a_impact_matrix.png)

---

## (8) SE 5 維度

**Slide 19 用**

![SE 5 維度](figures/15_se_perspective/fig15a_se_dimensions.png)

---

## (9) Whole-Genome 統計

**Slide 15 用**

![WG 統計](figures/20_whole_genome/wg_summary.png)

---

# § Part 3：IGV 截圖（per-site 視覺證據）

## (10) SP1 chr19:17565944 — Self-phasing extreme（4 ver 對比）

**Slide 7 用**

![SP1 4ver](../figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png)

⚠ **注意**：此圖在父層 figures 目錄，需用 `../figures/igv_v5_audit/...` 路徑（從 v5_audit_suite 看）

---

## (11) A_TP04 chr16:35118902 — V5 唯一強改善

**Slide 18 用**

![TP04 4ver](../figures/igv_v5_audit/by_HP_4ver/A_TP04_chr16_35118902.png)

---

## (12) C_V5max3 chr19:7405500 — V5 反向 regression

**Slide 18 用**

![V5max3 4ver](../figures/igv_v5_audit/by_HP_4ver/C_V5max3_chr19_7405500.png)

---

## (13) SP1 6-Panel + Paired + LOH/GE BED（新版 IGV）

**Slide 18 補充**

![SP1 6-panel paired](figures/09_purity06/igv_with_paired_loh/SP1_chr19_17565944.png)

**特色**：含 paired_normal + paired_tumor + 4 BAM versions + 8 BED tracks

---

# § Part 4：圖示美觀檢核 Matrix

## 4.1 配色一致性檢核

| 元素 | 應用色 | 是否一致 |
|------|------|:-:|
| HP1 family | 紅 #E91E63 (IGV) | ✓ |
| HP2 family | 藍 #1976D2 (IGV) | ✓ |
| Somatic | 橘 #E67E22 | ✓ |
| HP33 ambiguous | 紫 #7B1FA2 | ✓ |
| 主敘事 | STRUCT_BLUE #1F3A93 | ✓ |
| Bug / 紅警 | #D93025 | ✓ |
| 修復 / 綠 OK | #0F9D58 | ✓ |

## 4.2 字型檢核

| 項目 | 應用 |
|------|:-:|
| 中文字型 | Droid Sans Fallback (CJK) |
| 英文字型 | DejaVu Sans |
| 字級 | base 11，標題 14-16，副標 9-10 |
| Emoji 避免 | ⚠ fig23e 用 ❌ 字符可能缺字型 |
| Emoji 避免 | ⚠ fig23f 用 📊 emoji 可能缺字型 |

## 4.3 已知問題（matplotlib UserWarning）

| 圖 | Warning |
|----|---------|
| fig23e | `Glyph 10060 (CROSS MARK) missing` |
| fig23f | `Glyph 128202 (BAR CHART) missing` |

→ 雖然有 warning，但 PNG 仍可生成。可後續用 ASCII 替換（如 `[X]` 取代 `❌`）。

---

# § Part 5：使用建議

## 5.1 PPT 製作流程

1. 開啟 PowerPoint，使用 OFFICE_LIGHT 風格 template
2. 按 23 號報告 slide-by-slide 大綱組織內容
3. 從本 gallery 複製對應圖片路徑插入 slide
4. 按 23 號 §附錄「Slide Layout 模板」調整版面
5. 配色檢核（4.1）+ 字型檢核（4.2）

## 5.2 修改圖建議

如需修改某張圖：
1. 編輯 `/tmp/gen_ppt_figures.py`（fig23a-f）
2. 編輯 `/tmp/gen_*.py`（既有 fig01/13/14/15/20）
3. 重跑 Python 腳本
4. 替換對應 PNG

## 5.3 Codex CLI 替代方案（用戶詢問）

`codex` cli 不直接支援 imagegen（其 commands 是 `exec` / `review` / `apply` 等 code 工作）。

替代方案：
- ✅ 已用 matplotlib 對齊 longphase 論文 OFFICE_LIGHT 風格生成 6 張
- 若需 AI imagegen：考慮 DALL-E API / Midjourney / Stable Diffusion，但本專案 figures 目前都是技術示意圖，matplotlib 即足夠

---

# § Part 6：圖檔總清單與大小

## 6.1 本 session 新生成（6 張）

| 路徑 | 大小 |
|------|:-:|
| `figures/23_ppt_narrative/fig23a_version_timeline.png` | 173 KB |
| `figures/23_ppt_narrative/fig23b_metric_map.png` | 167 KB |
| `figures/23_ppt_narrative/fig23c_research_tree.png` | 193 KB |
| `figures/23_ppt_narrative/fig23d_dual_layer_bug.png` | 221 KB |
| `figures/23_ppt_narrative/fig23e_data_comparison.png` | 155 KB |
| `figures/23_ppt_narrative/fig23f_ppt_structure.png` | 187 KB |

## 6.2 既有 audit suite 圖（重用）

| 圖類別 | 數量 | 路徑 |
|--------|:-:|------|
| 01 code_diff | 4 | `figures/01_code_diff/fig01a-d.png` |
| 02 concordance | 2+ | `figures/02_concordance/` |
| 06 sanity | 2+ | `figures/06_sanity/` |
| 09 purity06 | 3+ | `figures/09_purity06/` |
| 12 GT distribution | 4 | `figures/12_gt_distribution/figA-D.png` |
| 13 phase vs tag | 3 | `figures/13_phase_vs_tag_algo/figA-C.png` |
| 14 impact | 1 | `figures/14_impact_classification/fig14a.png` |
| 15 SE perspective | 3 | `figures/15_se_perspective/fig15a-c.png` |
| 20 whole_genome | 1 | `figures/20_whole_genome/wg_summary.png` |
| **23 ppt narrative** | **6** | `figures/23_ppt_narrative/fig23a-f.png` |
| IGV 4ver | 13+ | `../figures/igv_v5_audit/by_HP_4ver/` |
| IGV 6-panel paired | 9 | `figures/09_purity06/igv_with_paired_loh/` |

**總計**：~50+ PNG 圖檔可用於 PPT。

---

# § Part 7：用戶檢核要求清單

請逐項確認並提供 feedback：

1. **fig23a 修法時間軸**：是否清楚顯示 5 個版本演進？需要重新著色／調整字型？
2. **fig23b Metric Map**：4 layer × 5 columns 表格是否清楚？需要調整 column 寬度？
3. **fig23c 研究發展樹**：5 大目標 + 3 priority 是否平衡？需要重新排版？
4. **fig23d 雙層 bug 機制**：左右對照是否對稱？需要調整文字？
5. **fig23e 數據對比**：log scale + distance bar 是否同時清楚？
6. **fig23f PPT 結構圖**：6 sections + slide 細項是否擁擠？

如需 codex AI 進階優化（例如生成更精緻的 illustrative 圖）：
- 提供 prompt + 期望風格參考
- 我可以準備 codex exec 命令讓你手動執行

---

# 跨檔索引

| 內容 | 路徑 |
|------|------|
| 整合敘事 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/22_pi_presentation_integrated_narrative.md` |
| Slide 大綱 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/23_ppt_slide_by_slide_outline.md` |
| **本 Gallery** | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/24_ppt_figures_gallery.md` |
| 主 INDEX | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md` |
| 圖生成腳本 | `/tmp/gen_ppt_figures.py` |
