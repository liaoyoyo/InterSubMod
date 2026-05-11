---
title: V5 PI 簡報 Slide-by-Slide 大綱（22 slides）
date: 2026-04-30
author: liaoyoyo2001
tags: [pi-presentation, ppt-outline, slide-design, v5]
status: draft_for_review
audience: PI presentation source
purpose: 22 張 slide 完整內容草稿（標題 / 重點 / 圖示 / 講稿）
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/22_pi_presentation_integrated_narrative.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/24_ppt_figures_gallery.md
---

# V5 PI 簡報 — Slide-by-Slide 大綱（22 slides）

## 設計原則

**對齊 longphase-to / longphase-s 論文 OFFICE_LIGHT 風格**：
- 白底 + STRUCT_BLUE (#1F3A93) 主敘事
- HP1 紅 (#E91E63 IGV) / HP2 藍 (#1976D2 IGV)
- Somatic 橘 (#E67E22)
- Ambiguous (HP33) 紫 (#7B1FA2)
- 每張 slide：標題 ≤ 15 字 + 2-3 焦點 + 結論句（非數據）+ 圖示 60-70%

**整體結構**：![Section 結構](figures/23_ppt_narrative/fig23f_ppt_structure.png)

---

# Section 1：開場（slides 1-3）

## Slide 1：Title — V5 修正 self-phasing tag bias

**標題**：LongPhase-TO V5 — 修正 17.3:1 HP tag bias 還原下游分析可信度

**副標**：從 baseline 雙層 bug 到 V5+ purity fix，22 份 audit suite 完整證明

**視覺**：
- 大標 + 副資訊（日期、作者、PI 報告）
- 兩條 HP 短帶（紅 HP1、藍 HP2 IGV 慣例）作裝飾
- 對應論文 slide 01 cover 風格

**講稿**：
> 「今天要報告 LongPhase-TO 的 self-phasing tag bias 修正。Baseline 在 somatic ALT-only metric 上 HP1:HP2 = 17.3:1，遠偏離 paired truth 的 1:1。我做了 V2b/V3F/V5/V5+ 共 5 個修法，全部對齊原 LongPhase-TO 設計理念。」

---

## Slide 2：最重要結論（4 行）

**標題**：核心結論 — 修復成功 + 邊界清楚

**4 條結論**：
```
✅ 0.93 高純度: HP1:HP2 從 2.08:1 → 1.42:1，距 paired 改善 −47%
✅ 不傷 caller: F1 = 0.7157 → 0.7154 (Δ=−0.0003 噪音)
✅ 對齊原設計: V2b 重用 LongPhase-TO Two-Pass mechanism
⚠ 0.6 低 purity: V5 conservative tagging 副作用顯著
```

**視覺**：
- 4 個圓角方塊（綠/綠/藍/黃 = 4 個 status）
- 每塊一行重點
- 對應 longphase-s slide 11 conclusion 兩條粗體 + recap 風格

**講稿**：
> 「核心結論四條：第一，0.93 高純度場景修復顯著；第二，不影響 ClairS-TO calling；第三，完全對齊論文原設計；第四，0.6 低純度仍有副作用，後續需要 purity-aware 修正。」

---

## Slide 3：流程概覽

**標題**：5 階段流程 — 從問題確認到完整 audit suite

**5 階段**：

| 階段 | 動作 | 產出 |
|------|------|------|
| 1 確認 | IGV + WG 統計 | 17.3:1 真實存在 |
| 2 找根因 | 雙層 bug 分析 | phase + tag |
| 3 修復 | V2b/V3F/V5/V5+ | 5 修法 |
| 4 驗證 | 22 份 audit suite | 47 圖 + 20 TSV |
| 5 後續 | nuance / WG / boldness | 完整 traceability |

**視覺**：
- ![時間軸](figures/23_ppt_narrative/fig23a_version_timeline.png)
- 主圖佔下半部 60%

**講稿**：
> 「整個工作流程分 5 階段。最關鍵的是用 22 份 audit suite 文件確保每個結論都有可追溯的證據鏈。」

---

# Section 2：問題定義（slides 4-7）

## Slide 4：HP tag 5 種編碼定義

**標題**：HP tag 編碼 — 5 種值 / 2 個 family

**內容**：
| HP:i: | 含義 | Family |
|:--:|------|:--:|
| 0 | untagged | — |
| 1 | germline HP1 | **HP1** |
| 11 | hp1-1 (somatic from HP1) | **HP1** |
| 2 | germline HP2 | **HP2** |
| 21 | hp2-1 (somatic from HP2) | **HP2** |
| 33 | hp3 (ambiguous) | (none) |

**視覺**：
- 表格 + 4 條彩色 HP tag 帶（紅/暗紅/藍/暗藍）
- 對應 longphase-s slide 5 objective 4 條 HP 結果風格

**講稿**：
> 「LongPhase-TO 的 haplotag 用 5 個 HP:i: 值。HP1 family = {1, 11}，HP2 family = {2, 21}，HP33 是 ambiguous somatic。」

---

## Slide 5：Phase 與 Tag 兩階段獨立

**標題**：Phase 階段 vs Tag 階段 — 職責分明

**流程圖**：
```
ClairS-TO snv.vcf + tumor BAM
        ↓
[Phase] → phased.vcf (GT/PS/GT2/PS2/GT3/PS3)
        ↓
[Haplotag] → tumor_tagged.bam (HP:i:N)
```

**重點**：
- Phase 在 VCF 上分配 GT；Tag 在 BAM 上分配 HP
- 兩階段獨立程式：PhasingProcess.cpp + HaplotagProcess.cpp

**視覺**：
- ![Phase vs Tag](figures/13_phase_vs_tag_algo/figA_phase_algorithm_flow.png) (既有)
- 對應 longphase-s slide 6 method workflow 5 階段風格

**講稿**：
> 「LongPhase-TO 分兩階段：phase 在 VCF 上，tag 在 BAM 上。獨立程式碼。後面我會證明 17.3:1 bias 是雙層問題的疊加。」

---

## Slide 6：問題量化 — 17.3:1 self-phasing bias

**標題**：HP1:HP2 = 17.3 : 1 — 嚴重偏離 1:1 預期

**3 個 metric scale**：
```
Specific subset (somatic ALT-only): HP1=614K, HP2=35.5K → 17.3 : 1 ❌
Whole-genome (all PASS sites):       HP1=1.5M, HP2=720K → 2.08 : 1 ⚠
Paired truth (PA_93):                HP1=326K, HP2=351K → 1 : 1.08 ✓
```

**視覺（推薦用 fig23g 精確版）**：
- ![17.3:1 雙層機制](figures/23_ppt_narrative/fig23g_somatic_bias_precise.png)
- 5 panels: [A] Baseline Phase + [B] V5 Phase + [C] Baseline Tag + [D] V5 Tag + [E] 結果三柱對比
- 舊版備選：`figures/01_code_diff/fig01d_somatic_bias_explanation.png`

**講稿**：
> 「17.3:1 是特定 subset 的 metric。在全 PASS sites WG metric 是 2.08:1，paired truth 是 1.08:1。本圖揭露雙層 bug：Phase 階段 somatic 進 graph anchor 強制連結同 clone reads；Tag 階段 getVote 把 somatic 列為第 1 優先。V5 雙層修復：PON-only flag 切斷 graph 連結 + getVote 三層決策。最終 baseline 2.08:1 → V5 1.42:1，距 paired truth 改善 47%。」

---

## Slide 7：極端位點 — chr19 self-phasing cluster

**標題**：Per-site 極端例子 — SP1/SP2/SP3 達 100×

**3 個位點**：
| 位點 | Baseline ratio |
|------|:-:|
| SP1 chr19:17565944 | **113 : 0** (∞) |
| SP2 chr19:12452332 | **109 : 1** |
| SP3 chr19:12467180 | **108 : 0** (∞) |

**視覺**：
- ![SP1 IGV](figures/igv_v5_audit/by_HP_4ver/D_SP1_chr19_17565944.png) (4 ver 對比)
- 對應 longphase-s slide 10 IGV 風格

**講稿**：
> 「個別位點可達 109× 以上。SP1 的 113 reads 全部 HP1，沒有 HP2。這不是統計噪音，是真實 bug。」

---

# Section 3：影響評估（slides 8-9）

## Slide 8：影響範圍 — TO pure + ISM 共用 tag

**標題**：為何優先修 — TO pure 與 ISM 都依賴 HP tag

**3 個影響**：
1. **TO pure benchmark**：所有 TP/FP characterization 基於 HP tag
2. **ISM 第一步分組**：HP1/HP2 reads 分組做 ASM
3. **量化影響**：94.6% somatic reads 標 HP1 → ISM ASM 統計力極弱

**視覺**：
- ISM 流程圖（Tagged BAM → 分組 → 比較 ASM → 推論 subclone）
- 對應 longphase-s slide 4 motivation 流程風格

**講稿**：
> 「Tag 錯 = ISM 錯。所有用到 HP1/HP2 的下游都受影響。修復 tag 是上游基礎建設工作。」

---

## Slide 9：與五大研究目標關聯

**標題**：V5 修復 = 五大目標基石

**視覺**：
- ![研究發展樹](figures/23_ppt_narrative/fig23c_research_tree.png)
- 5 大目標彩色分支（紅 / 橘 / 藍 / 紫 / 黃）

**結論**：5/5 目標都受 tag 影響 → V5 修復是 P0 任務

**講稿**：
> 「五大研究目標全部依賴 HP tag。修復 tag 後，每個目標都可以基於可信 baseline 重新展開。」

---

# Section 4：機制與修法（slides 10-14）

## Slide 10：Baseline `getVote()` 的順序錯誤

**標題**：Baseline 把 somatic 列為投票第 1 優先（順序錯）

**程式碼**：
```cpp
// Baseline (錯)：
variantKeys = {
    {HAPLOTYPE1_1, HAPLOTYPE2_1},   // ← 第 1 優先（somatic）❌
    {HAPLOTYPE3, HAPLOTYPE2_1},      //   第 2 優先
    {HAPLOTYPE1, HAPLOTYPE2}         // ← 第 3 優先（germline）❌
};
```

**問題**：read 涵蓋任何 somatic site → 第 1 優先立即決定 → germline 投票被覆蓋

**視覺**：
- 程式碼 + 流程圖
- 對應 longphase-s slide 7 filters 4 階段 priority 風格

**講稿**：
> 「commit 41ff147 直接記載 root cause：getVote 把 somatic 列為第 1 優先，導致 99.9% reads 變成 HP:i:21。」

---

## Slide 11：V5 三層決策 — germline first

**標題**：V5 修法 — Layer 1 / 1.5 / 2 三層拆分

**程式碼**：
```cpp
// Layer 1: Germline first（順序反轉）
if (germlineHP1 > 0 || germlineHP2 > 0) {
    germlineResult = (germlineHP1 >= germlineHP2) ? 1 : 2;
}
// Layer 1.5: Somatic fallback
else if (somaticHP1 > 0 || somaticHP2 > 0) {
    germlineResult = (somaticHP1 >= somaticHP2) ? 1 : 2;
}
// Layer 2: Encoding (separated)
if (somaticTotal > 0) hpResult = 11/21/33;
else hpResult = germlineResult;
```

**視覺**：
- ![三層決策](figures/13_phase_vs_tag_algo/figB_tag_getVote_decision.png) (既有)

**講稿**：
> 「V5 把 baseline 單體 if-else 拆成三層：germline 先，somatic 作 fallback，encoding 與 voting 分離。」

---

## Slide 12：具體例子 — 3 reads × 5 variants

**標題**：LOH 區具體例子 — V5 救起原本 lost 的 reads

**場景**：5 variants (2 germline + 3 somatic) × 3 reads

| Read | Baseline tag | V5 tag |
|------|:-:|:-:|
| Read A (cover g1+s1+s2+s3+g2) | HP:i:1 | **HP:i:11** |
| Read B (only s1+s2+s3) | HP:i:0 (lost) | **HP:i:11** (Layer 1.5) |
| Read C (only s1+s2) | HP:i:0 (lost) | **HP:i:11** (Layer 1.5) |

**視覺**：
- ![具體例子](figures/13_phase_vs_tag_algo/figC_concrete_example.png) (既有)

**講稿**：
> 「具體例子：三條 reads 涵蓋不同 variants。Baseline 把無 germline 的 Read B/C 完全丟失，V5 用 Layer 1.5 fallback 救起。」

---

## Slide 13：Phase 階段修法 — `--pon-only-phasing` flag

**標題**：V5 phase 修法 — 限制 anchor only PON germline

**程式碼**：
```cpp
// PhasingProcess.cpp:154-157
if (params.ponOnlyPhasing) {
    vGraph->convertNonGermlineToSomatic();
    // 把所有非 PON variants 標為 SOMATIC origin
    // → phasing graph 只用 PON germline 當 anchor
    // → 切斷 self-phasing 來源
}
```

**重要說明**：
- **`convertNonGermlineToSomatic()` 不是 V5 新增**
- 原 LongPhase-TO 已有此函式，用於 Two-Pass 高純度策略 (purity > 0.95)
- V5 把它從 conditional 升級為 user-controlled flag

**視覺**：
- 程式碼 + flow node 圖
- 對應 longphase-to slide 04 pipeline overview 風格

**講稿**：
> 「V5 重用 LongPhase-TO 既有 `convertNonGermlineToSomatic()` 函式，不是新發明。論文 Two-Pass 策略原本只在 purity > 0.95 啟用，V5 把它變成 user flag，可在任何 purity 下開啟。」

---

## Slide 14：V5 完整修法清單（5 項）

**標題**：V5 共 5 個修法 — 4 必須改 + 1 完整性補強

| # | 修法 | 性質 | 必要性 |
|:-:|------|:-:|:-:|
| 1 | getVote 順序反轉 | bug fix | 🟢 必須 |
| 2 | HP:i:33 enum 比對 | type bug | 🟢 必須 |
| 3 | UNDEFINED guards | defensive | 🟢 必須 |
| 4 | --pon-only-phasing flag | conditional | ⚠ 大膽 |
| 5 | collectPloidyRatio | completeness | 🟢 必須 |

**視覺**：
- ![影響分類矩陣](figures/14_impact_classification/fig14a_impact_matrix.png) (既有)

**講稿**：
> 「V5 共 5 項修法：4 個是必須改的 bug fix，1 個是 conditional flag 設計擴展。第 5 項 collectPloidyRatio 是本 session 修補的，解決 V5 PON-only mode purity = 0 問題。」

---

# Section 5：驗證證據（slides 15-19）

## Slide 15：證據 #1 — Whole-Genome 全 PASS sites

**標題**：WG 統計 — V5 距 paired truth 改善 −47%

**3 大數值**：
| Sample | HP1:HP2 | distance to PA |
|--------|:-:|:-:|
| Baseline @ 0.93 | **2.08:1** | 0.352 |
| **V5 @ 0.93** | **1.42:1** | **0.186** ✓ |
| Paired truth (PA_93) | **1:1.08** | 0 (基準) |

**視覺**：
- ![WG 統計](figures/20_whole_genome/wg_summary.png) (既有)
- 加上 fig23e 數據對比

**講稿**：
> 「全基因組 ~90K sites × 5 BAMs 跑了 2 小時。V5_93 比 baseline 距 PA truth 改善 47%。」

---

## Slide 16：證據 #2 + #3 — VCF GT 一致 + F1 不變

**標題**：Phase 階段無 bug + V5 不傷 caller

**兩證據**：

**(a) VCF GT 100% 一致**：
| | PASS Somatic_NoLOH | PASS Somatic_in_LOH |
|--|:-:|:-:|
| Baseline | 21,304 | 11,983 |
| V5 (V2b backbone) | 21,304 | 11,983 |
| Δ | **0** | **0** |

**(b) F1 不變**：
| 版本 | F1 | Δ |
|:-:|:-:|:-:|
| Baseline | 0.7157 | (基準) |
| V5 | 0.7154 | **−0.0003 (噪音)** |

**視覺**：
- ![GT 分布](figures/12_gt_distribution/figA_gt_class_by_filter.png) (既有)

**講稿**：
> 「Phase 階段 GT 100% 一致，差異全在 tag 階段。F1 變化僅噪音，V5 不影響 calling。」

---

## Slide 17：證據 #4 — Sanity check 5/5 PASS

**標題**：5 項硬性檢查全部通過

**5 項守恆律**：
| 檢查 | 結果 |
|------|:-:|
| 守恆律 A：Δ_HP33 + Δ_(HP11+HP21) = 0 | ✅ 15/15 |
| 守恆律 B：germline (HP1, HP2) 不變 | ✅ 15/15 |
| Layer 1.5 預期 1：新 HP11/HP21 來自 V3F HP33 | ✅ 15/15 |
| Layer 1.5 預期 2：無 germline → HP33 違規 | ✅ 0 violations |
| Untagged → directional：V5 不強行 tag | ✅ 0 violations |

**視覺**：
- 5 個綠色 ✓ 圖標 + 數字
- 對應 longphase-s slide 9 results 3 個 bar chart 風格

**講稿**：
> 「5 項硬性 sanity check 全部通過。守恆律驗證 V5 沒有破壞既有資訊，Layer 1.5 行為符合預期。」

---

## Slide 18：證據 #5 — IGV 視覺案例（強改善 + 反向 regression）

**標題**：IGV 4-version 對比 — 真實案例呈現

**兩個關鍵案例**：

**(a) A_TP04 chr16:35118902 — 唯一強改善**：
- exact-site +0.1427
- paired conditional **+0.9737**
- ![TP04](figures/igv_v5_audit/by_HP_4ver/A_TP04_chr16_35118902.png)

**(b) C_V5max3 chr19:7405500 — V5 反向 regression（誠實揭露）**：
- V5 reassign 後反偏離 paired
- ![V5max3](figures/igv_v5_audit/by_HP_4ver/C_V5max3_chr19_7405500.png)

**視覺**：兩張並排 IGV 截圖（對應 longphase-s slide 10 IGV 風格）

**講稿**：
> 「TP04 是 V5 唯一強改善案例。但 V5max3 顯示 V5 不是全面改善，有些 site 反而偏離 paired。誠實揭露讓結論可信。」

---

## Slide 19：業界對齊與矛盾盲點

**標題**：V5 設計符合業界標準 + 揭露 metric bias

**業界對齊（5 項）**：
| 業界方法 | 對應 V5 |
|---------|------|
| LongPhase v1.7 paired mode | germline-anchored phasing |
| ClairS-TO PoN-aware | --pon-only-phasing flag |
| WhatsHap haplotag | germline-first priority |
| GATK PhaseByTransmission | 主→輔 priority pattern |
| LongPhase-TO Two-Pass | reuse `convertNonGermlineToSomatic` |

**揭露 metric bias**：
- +6.65pp / +13.3pp 是 **conditional accuracy**
- WG fixed-denom 重算 ≈ 持平
- HP33 在現有 metric 被低估

**視覺**：
- ![SE 5 維度](figures/15_se_perspective/fig15a_se_dimensions.png) (既有)

**講稿**：
> 「V5 設計符合 5 項業界標準。但我們也誠實揭露：+13pp 是 conditional accuracy，WG fixed-denom 重算 V5 與 BL 在 site-level 持平。HP33 conservative 在現有 metric 被系統性低估。」

---

# Section 6：未來規劃（slides 20-22）

## Slide 20：立即接續工作

**標題**：V5 修復後 — 立即可做的 3 件事

**3 個 P0 任務**：
1. **重跑 ISM 分析**：用 V5 tagged BAM 重做 allele-specific methylation (1-2 day)
2. **TP/FP characterization**：用 V5 BAM 看 TP/FP 在 HP1/HP2 真實分布 (0.5-1 day)
3. **更新 audit suite §02-08**：加入新發現的 nuance (0.5 day)

**視覺**：
- 3 個彩色任務方塊（對應 longphase-to slide 11 conclusion 3 條結論風格）

**講稿**：
> 「V5 修復後最重要的是重新跑 ISM。之前所有 HP1/HP2 觀察都帶 17.3:1 偏差。」

---

## Slide 21：研究發展樹 — 五大目標關聯

**標題**：V5 = 五大研究目標基石

**視覺**：
- ![研究發展樹](figures/23_ppt_narrative/fig23c_research_tree.png) (新)

**講稿**：
> 「V5 修復解鎖五大研究方向。立即（P0）/ 中期（P1）/ 長期（P2）三層任務排序，從 ISM 重跑到 trio sequencing。」

---

## Slide 22：整體進展時間軸

**標題**：從 V0 → V5+ 的完整時間軸

**視覺**：
- ![時間軸](figures/23_ppt_narrative/fig23a_version_timeline.png) (新)

**結尾**：
- ✅ Audit suite 22 份報告完成
- ✅ 修復 + 驗證 + 文檔 完整 traceability
- ✅ 下一步：基於 V5 baseline 展開五大目標研究

**講稿**：
> 「總結：從 4-10 V2b 第一個 commit 到 4-29 V5+ purity fix，本 session 補上完整性。22 份 audit suite 文件提供完整 traceability，可作為後續所有 TO 研究的基線參考。」

---

# 附錄：Slide Layout 模板

每張 slide 建議 layout：

```
┌─────────────────────────────────────────────────────────┐
│  [標題 ≤ 15 字]                              [章節編號]  │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  [2-3 個重點 / 表格 / 程式碼]   [圖示 60-70%]            │
│                                                         │
│                                                         │
├─────────────────────────────────────────────────────────┤
│  [結論句（非數據）]                                      │
└─────────────────────────────────────────────────────────┘
```

**設計檢核清單**：
- [ ] 標題 ≤ 15 字
- [ ] 每頁 2-3 焦點，不超載
- [ ] 結論句是 takeaway 不是數據
- [ ] 圖佔比 60-70%
- [ ] 配色用 OFFICE_LIGHT
- [ ] 中英文混排：英文縮排 0.25", 60% 字號
- [ ] 禁止斜體（標題可粗體）

---

# 跨檔索引

| 內容 | 路徑 |
|------|------|
| 整合敘事報告 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/22_pi_presentation_integrated_narrative.md` |
| **本 slide 大綱** | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/23_ppt_slide_by_slide_outline.md` |
| Gallery + 圖示說明 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/24_ppt_figures_gallery.md` |
| 主 INDEX | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md` |
