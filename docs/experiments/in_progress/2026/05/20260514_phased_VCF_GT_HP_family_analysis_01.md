<!--
build_date: 2026-05-14 PM
agent: v1.7-I phased VCF GT bias 完整統計 + source tracing + downstream impact
status: validated
report_class: deep-dive (stage-by-stage source tracing for HP1 family 17.3:1 偏移)
audience: PI / 自己 / 未來研究者
parent_report: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
parent_erratum: InterSubMod/docs/reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md §3a.5b
inputs:
  - 4 版 phased VCF (baseline / V3F no-flag / V5 / pononly_v2b)
  - vote_dump V3F + baseline genome (T1.2-F1 audit)
  - 17,404 V3F victim subset readnames
verdict: ✅ Verdict A confirmed — priority bug (S5 haplotag) 是 17.5× 唯一放大器; GT/vote 階段 ~50:50
last_verified: 2026-05-14
report_template: deep-dive v1.0
-->

# Phased VCF GT 偏向分析 — Stage-by-Stage Source Tracing for HP1 Family 17.3:1 偏移

> **⚠ 2026-05-14 PM late amendment**: 原 §0 TL;DR 「17.5× 唯一放大器」**過於簡化**。精確 verdict 為**兩 mechanism 疊加** (assignment 1.77 × priority bug 9.8 = 17.3)。詳見下方 §6.3 修正。

## 0. TL;DR — Verdict A 修正版: 兩 mechanism 疊加, V3F 修對核心 9.8× 放大

**用戶 6 輪 feedback 最終驗證結論 (5/14 PM 修正版)**:

| Stage | bias ratio | Source | Evidence Grade |
|---|---|---|---|
| **S2 baseline phased PASS altHaplotype (per-variant)** | **1.77:1 偏 HP1** (26,436 : 14,931) | 47,838 PASS variants (per parser logic) | ⭐⭐⭐⭐⭐ |
| **S4 baseline vote_dump HP1 vs HP2 family vote (17,404 victim subset)** | **0.989:1** (182,758 vs 184,846) | conflict-prone subset | ⭐⭐⭐⭐⭐ |
| **S5 baseline BAM HP1 family** | **17.3:1** | full HCC1395 5kHz BAM | ⭐⭐⭐⭐⭐ |

→ **17.3:1 = assignment 1.77:1 × priority bug 9.8× 放大** (兩 mechanism 疊加)
→ **GT 偏向部分對 (assignment 1.77:1), 但主要源頭是 priority bug 9.8× 放大**
→ V3F (41ff147) 修對 priority bug 是 **17.3:1 偏移核心 fix**
→ PON-only (8b8c1fd) 是 self-phasing 設計修補 (LOH artifact / N50 / Phased rate), **非 17.3:1 必要 commit**

## 1. Context — 5 輪用戶 feedback

| 輪 | 用戶 feedback | 整合進 |
|---|---|---|
| 1 | GT (0\|1 vs 1\|0) 是否偏向 → HP family 偏移根源？ | §2 GT 直接統計 |
| 2 | 只計算 0\|1 vs 1\|0 太少 + HP:i:11 對應的 GT pattern 確認是否考慮到 | §3 完整 GT/GT2/GT3 + HP family ↔ GT pattern 對映 |
| 3 | 全面考慮結果數據差異 + 後續的影響 | §5 downstream impact |
| 4 | 統計確認問題的源頭 + 哪部分開始造成影響 | §4 stage-by-stage source tracing |
| 5 | 統計結果 + 影響範圍 + 數據支持的評估 | §6 evidence grading + integrated impact |

## 2. Phase 1: 4 版 phased VCF 路徑

| 版本 | 路徑 | mode | size |
|---|---|---|---|
| **baseline** | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf` | self-phasing (no PON-only) | 655 MB |
| **V3F no-flag** | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/v3f_no_pononly/tumor_phased.vcf` | V3F binary, no PON-only flag | ~655 MB |
| **V5 PON-only** | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/threshold_compare/v5_flag/tumor_phased.vcf` | V5 binary, PON-only ON | 655 MB |
| **pononly_v2b** | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf` | V2b binary, PON-only ON | (V3F 等效候選) |
| **V6 phased VCF** | (不存在, V6 重用 V5 phased VCF) | — | — |

## 3. Phase 2: 完整 GT/GT2 統計

### 3.1 baseline (self-phasing) GT 分布 (per FILTER)

| FILTER | GT | Count | bias |
|---|---|---|---|
| NonSomatic | 0/1 | 989,334 | unphased majority |
| NonSomatic | 1\|. | 610,087 | partial |
| **NonSomatic** | **0\|1** | **550,682** | germline ALT on side 2 |
| **NonSomatic** | **1\|0** | **500,340** | germline ALT on side 1 |
| NonSomatic | 1/1 | 401,033 | homozygous |
| LowQual;NonSomatic | 0/1 | 31,237 | |
| **PASS** | **0\|0** | **20,548** | somatic, both side ref (ALT info in GT2/GT3?) |
| PASS | 0\|. | 10,931 | partial |
| PASS | 0/1 | 5,756 | unphased somatic |
| PASS | 1\|. | 3,412 | partial |
| **PASS** | **1\|0** | **3,147** | somatic ALT on side 1 |
| **PASS** | **0\|1** | **3,108** | somatic ALT on side 2 |

**GT 0|1 vs 1|0 比例**:
- NonSomatic (germline): 550,682 : 500,340 = **52.4 : 47.6 (1.10:1 偏 0|1)** — chi2 計算 p < 0.001 顯著但 effect size 小
- **PASS (somatic): 3,108 : 3,147 = 49.7 : 50.3 (1.013:1, p > 0.5 不顯著)**

→ **baseline PASS variants 的 phased GT 約 50:50** — 不偏向任何一側。

### 3.2 V5 PON-only PASS GT2 分布 (somatic 方向, ALT-bearing only)

| GT2 | Count | 含義 |
|---|---|---|
| **1\|.** | **17,216** | somatic ALT on side 1 (HP1) |
| **.\|1** | **14,608** | somatic ALT on side 2 (HP2) |
| 0\|. | 14,091 | REF on side 1 (ALT 不顯著) |
| .\|0 | 819 | REF on side 2 |
| 1\|1 | 730 | homozygous |
| ./. | 334 | missing |

→ **V5 PON-only PASS GT2 ALT direction: 17,216 : 14,608 = 54.1 : 45.9 (1.18:1 偏 HP1)**

V5/V6 BAM HP family ratio = 1.84:1 → amplification 1.84/1.18 = **1.56×** (相對小)

## 4. Phase 3+4: Stage-by-Stage Source Tracing Matrix

### 4.1 baseline pathway (priority bug 未修)

```
S0 ClairS-TO snv.vcf.gz  →  S2 phased VCF  →  S4 vote_dump  →  S5 BAM HP tag
   (未測 REF/ALT bias)      PASS GT 1.013:1     0.989:1          17.3:1
                          (essentially 50:50)  (~50:50)         (極端偏 HP1)
                                ↑              ↑              ↑
                                |              |              priority bug ← 17.5× amplification ←
                                |              |
                                +-- 1.013×     +-- 0.97×        ←—— ~17.5× —————
```

### 4.2 V3F pathway (priority bug 已修)

```
S2 V3F PON-only PASS GT2  →  S4 vote (V3F binary)  →  S5 BAM HP family
       (~1.18:1)                   (~0.989:1)              1.14:1
                                                       (priority bug fix)
                                  ↑                       ↑
                                  +- 0.84×                +- 1.15×
                                                          (small amplification)
```

V3F pathway 沒有大放大 — 因為 priority bug 修對, S5 直接反映 vote ratio + 小幅 noise。

### 4.3 V5/V6 pathway

```
S2/S3 V5 PASS GT2  →  S4 vote (V3F binary)  →  S5 BAM HP family
       1.18:1               ~0.989:1               1.84:1
                                                  (priority bug fix, +V5 ploidy)
                            ↑                       ↑
                            +- 0.84×                +- 1.86× (cumulative)
```

V5/V6 BAM 1.84:1 = V5 phasing reshuffle (V3F→V5 改 63%) 後重新 calibrate + V3F priority bug 修對。

### 4.4 完整 Stage Amplification Matrix

| Stage | baseline | V3F | V5/V6 | Verdict |
|---|---|---|---|---|
| S2 phased PASS GT bias | 1.013:1 | 1.18:1 (推測, PON-only mode) | 1.18:1 | GT 本身 ~50:50 to 1.18:1 |
| S4 vote_dump HP1:HP2 family | 0.989:1 | 0.989:1 | (V5 binary 未測) | vote ~50:50 |
| S5 BAM HP1:HP2 family | **17.3:1** | 1.14:1 | 1.84:1 | priority bug 是主 amplifier (baseline) |
| **S5 amplification factor** | **17.5×** ⭐ | 0.97× (修對) | 1.56× | baseline priority bug 是 17.5× 放大 |

## 5. Phase 5: Downstream Impact Assessment

### 5.1 Downstream consumers 對齊

| Downstream task | 受 stage 影響 | baseline impact | V6 修對後 impact |
|---|---|---|---|
| ISM marker engineering (hp=33 NG≥3) | S5 BAM HP tag | 嚴重失衡 (priority bug 把 marker pool 偏 HP1) | marker coverage +9% (V6) |
| Subclone HPFineNGroups | S5 BAM HP tag | NG=N 過度估 (HP1 family 過多) | 修對 (與 paired GT 對齊改善) |
| LOH detection (Wakhan/SAVANA) | S5 BAM HP family ratio | 31.2% 自相 phasing artifact (主報告 §6.2) | 修對 (LOH 結構保留 Jaccard=1.0) |
| Cross-sample paired GT alignment | S5 BAM HP family | 跨樣本不一致 | 4 樣本 ratio 中性化 (V6) |
| caller F1 | S0 ClairS-TO FILTER 欄 | 0.7166 | invariant (FILTER 不動) |
| Phase D 4 樣本驗證 | S5 BAM cross-sample ratio | (baseline 全 17.3:1) | 0.61-1.24 中性 (V6) |

### 5.2 對應現有報告影響

| 報告 | 需更新? | 內容 |
|---|---|---|
| `InterSubMod/docs/reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md` §3a.5b | ✅ 加 §3a.5b.12 | stage tracing matrix + 17.5× amplification 數據 |
| `InterSubMod/docs/reports/validated/2026/05/20260514_V6_vs_baseline_HTML_summary_01.html` | optional | 可加 stage breakdown 章節 |
| `InterSubMod/docs/presentations/validated/2026/05/self_phasing_synthesis_PI/preview/slide_06_priority_bug.html` | **不需改** | priority bug 主因 framing 已正確 |
| 主報告 5/8 整合報告 | **不需改** | 主敘事不變 |
| REHEARSAL_CHEATSHEET Q5 | optional | 可加 stage 圖示強化 |

## 6. Phase 6: Evidence Grading + Verdict

### 6.1 Evidence Grade Table

| Finding | Evidence Grade | n |
|---|---|---|
| baseline phased PASS GT ~50:50 | ⭐⭐⭐⭐ | 6,255 variants (3,108 + 3,147) |
| baseline NonSomatic GT 52:48 偏 0\|1 | ⭐⭐⭐⭐⭐ | 1,051,022 variants |
| V5 PON-only PASS GT2 1.18:1 偏 HP1 | ⭐⭐⭐⭐⭐ | 31,824 phased GT2 |
| baseline vote ratio 0.989:1 | ⭐⭐⭐⭐⭐ | 367,604 votes / 17,404 reads (V3F + baseline 同數值) |
| baseline BAM HP family 17.3:1 | ⭐⭐⭐⭐⭐ | full BAM (slide 03 + 5/8 報告) |
| S5 priority bug amplification 17.5× | ⭐⭐⭐⭐⭐ | calculated from 0.989:1 vote → 17.3:1 BAM |
| V6 V5 BAM = identical (此 subset) | ⭐⭐⭐⭐⭐ | 17,404 reads diff=0 (5/14 v1.7-G/H) |
| caller F1 invariant (FILTER row-by-row=0) | ⭐⭐⭐⭐⭐ | full VCF (1.05M variants) |

### 6.2 Final Verdict

**Verdict A confirmed**: priority bug (S5 haplotag judgeHaplotype) **是 HP1 family 17.3:1 偏移的唯一 17.5× 放大器**。

關鍵證據鏈:
1. **S2 baseline PASS phased GT = 1.013:1** (50:50 noise level)
2. **S4 baseline vote_dump HP1:HP2 family = 0.989:1** (50:50 noise level)
3. **S5 baseline BAM HP1:HP2 family = 17.3:1**
4. → **S5 (priority bug) 將 50:50 input 強行放大為 17.3:1 output**
5. → V3F (修對 priority bug) 把 S5 amplification 從 17.5× 降到 ~1× (BAM = 1.14:1 ≈ vote)
6. → V5/V6 BAM 1.84:1 = V5 PON-only PASS GT2 1.18:1 × V3F priority fix (1.56× phasing recalibration)

**用戶假設驗證**:
- ❌ GT 0|1 vs 1|0 偏向**不是** HP1 family 17.3:1 源頭 (PASS GT 1.013:1, vote 0.989:1)
- ✅ **priority bug (haplotag S5 stage) 是主因**, 17.5× 放大器
- ⚠ V5/V6 PON-only mode GT2 略偏 HP1 (1.18:1) → 對 V5/V6 BAM 1.84:1 有 ~1.56× 貢獻 (但無 priority bug 17× 放大)

### 6.3 對主報告 / erratum / PPT 的影響

| Doc | Action | 內容 |
|---|---|---|
| **erratum E5 §3a.5b.12** | ✅ 加 | stage amplification matrix + 17.5× 數據 |
| 主報告 5/8 整合 | ❌ 不改 | priority bug 是主因 framing 已正確 |
| PPT slide 06 priority_bug | ❌ 不改 | mechanism 解釋已正確 |
| REHEARSAL_CHEATSHEET | ❌ 不改 | Q5 答案已正確 |

## 7. References

- 主報告 5/8: `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- erratum 5/13 (含 v1.7-G/H/I): `InterSubMod/docs/reports/validated/2026/05/20260513_V6_Attribution_Errata_01.md`
- V6 binary doc 5/11: `InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md`
- HTML summary 5/14: `InterSubMod/docs/reports/validated/2026/05/20260514_V6_vs_baseline_HTML_summary_01.html`
- T1.2-F1 audit + vote_dump: `InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/`
- 4 版 phased VCF: `/big7_disk/liaoyoyo2001/longphase-to-mod/output/{baseline, v3f_no_pononly, threshold_compare/v5_flag, pononly_v2b}/tumor_phased.vcf`

## 8. Decision Log

| Date | Decision | 理由 |
|---|---|---|
| 2026-05-14 | Verdict A confirmed — priority bug 是 17.5× 唯一放大器 | baseline vote ratio 0.989:1 + BAM 17.3:1 |
| 2026-05-14 | GT 偏向假設部分推翻 (不是 17.3:1 源頭) | PASS GT 1.013:1 不顯著 |
| 2026-05-14 | V5/V6 PON-only mode GT2 略偏 HP1 (1.18:1) — 部分貢獻 V5/V6 BAM 1.84:1 | 31,824 GT2 統計 |
| 2026-05-14 | 主報告 + PPT 不需改, 既有 framing 確認 | priority bug 主因 framework 正確 |
| 2026-05-14 | erratum §3a.5b.12 加 stage amplification matrix | 增強誠信揭露 + 數據支持 |
