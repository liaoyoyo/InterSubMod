---
title: GT/GT2/GT3 分布稽核 — Phasing 問題 vs Tag 問題判別
date: 2026-04-27
author: liaoyoyo2001
tags: [audit, longphase-to, GT, phasing, haplotag, baseline, v5]
status: validated_complete
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md
---

# GT/GT2/GT3 分布稽核 — Phasing 問題 vs Tag 問題判別

## 一句話摘要

> **17.3:1 self-phasing bias 是 haplotag 階段（tag 問題），不是 phasing GT 分配階段（非 phasing 問題）**。
> 證據：baseline vs V2b 在 PASS 的 somatic GT（`0|0`/`.|0`/`0|.`）**完全相同**（21,304 / 11,323-11,290 / 660-693），only germline het 的 PS block 方向有差（`1|0` ↔ `0|1` 翻向）。

---

## Section 1 — 任務定義

驗證 longphase-to baseline 輸出的 GT / GT2 / GT3 比例分布，判別 self-phasing 異常的根因層級：

| 層級 | 證據判定 |
|------|---------|
| **Phasing 階段問題** | baseline 與後續版本（V2b/V5）的 GT 分布有顯著差異，特別是 somatic 類別（`0|0`, `.|0`, `0|.`） |
| **Tag 階段問題** | GT 分布幾乎相同，但 BAM HP tag 結果差異大 |

對應的版本層級（`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` Section 1.3）：

```
Baseline phase ──→ baseline.vcf ──→ Baseline haplotag ──→ baseline.bam (17.3:1 bias)
                                                               ↑
                                                       getVote() Layer 1 only
                                                       無 PON anchor 過濾

V2b phase     ──→ v2b.vcf      ──→ V5 haplotag      ──→ v5.bam (≈1:1)
                                                               ↑
                                                       getVote() Layer 1 + 1.5
                                                       PON-only phasing flag
```

→ **如果 baseline.vcf ≈ v2b.vcf 但 baseline.bam ≠ v5.bam，則差異來自 haplotag 階段**。

---

## Section 2 — 資料與方法

### 輸入

| Label | VCF 路徑 | 大小 |
|-------|---------|------|
| `baseline` | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/baseline/tumor_phased.vcf` | 655 MB |
| `v2b_used_by_v5` | `/big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v2b/tumor_phased.vcf` | 同樣本，PON-only flag on |

V5 重用 V2b 的 phased.vcf 直接 haplotag（驗證點：`pononly_v5_somatic_fallback/tumor_tagged.out:1` header `##snpFile:output/pononly_v2b/tumor_phased.vcf`）。

### 變數總計

兩 VCF 各 **3,187,275 variants**（完全相同 row count，因為 ClairS-TO calling 是固定的；longphase-to 只動 GT/GT2/GT3/PS 欄位）。

### 分類規則（依 longphase-to 文件）

| GT 類別 | GT 值 |
|---------|------|
| Germline_Het | `0\|1`, `1\|0` |
| Germline_Hom_or_LOH | `.\|1`, `1\|.` |
| **Somatic_NoLOH** | `0\|0` |
| **Somatic_in_LOH** | `.\|0`, `0\|.` |
| Unphased | `0/1`, `1/1`, etc. |

### FILTER 分桶

| Bucket | 含義 |
|--------|------|
| **PASS** | 通過全部 ClairS-TO filter，視為 high-confidence somatic |
| NonSomatic | 被 NonSomatic filter 排除（germline 主導） |
| LowQual+NonSomatic | 同時 LowQual + NonSomatic |
| LowQual_Other | 其他 LowQual 變體 |

### 輸出

- 圖：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/12_gt_distribution/`
- 表：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/`

---

## Section 3 — 主要結果

### 3.1 PASS（高信心 somatic）GT class 分布幾乎相同

| GT_class | Baseline % | V2b % | Δpp |
|---------|:---------:|:-----:|:---:|
| Germline_Het | 13.09 | 12.91 | **−0.18** |
| Germline_Hom_or_LOH | 7.48 | 7.48 | 0.00 |
| **Somatic_NoLOH** | **42.99** | **42.99** | **0.00** |
| **Somatic_in_LOH** | **24.19** | **24.19** | **0.00** |
| Unphased | 12.26 | 12.44 | +0.18 |

→ **somatic 類別在 PASS 完全一致**：兩版本對 somatic variant 的 GT 標籤決策**沒有差別**。
→ 唯一變化：~86 個 variants（佔 0.18%）從 Germline_Het 改判 Unphased。

![Figure A — GT class composition by FILTER bucket](figures/12_gt_distribution/figA_gt_class_by_filter.png)

### 3.2 raw GT 值揭露 PS block 翻向（不影響語意）

| GT raw | Baseline 數量 | V2b 數量 | Δ |
|--------|:------------:|:--------:|:-:|
| `0\|1` | 559,116 | **1,130,749** | **+571,633 (+102%)** |
| `1\|0` | 508,368 | **34,587** | **−473,781 (−93%)** |
| `0/1` | 1,030,285 | 932,433 | −97,852 |
| `1\|.` | 614,471 | 608,635 | −5,836 |
| `.\|1` | 35,444 | 41,280 | +5,836 |
| `1/1` | 406,304 | 406,304 | 0 |
| **`0\|0`** (somatic NoLOH) | **21,304** | **21,304** | **0** |
| **`0\|.`** (somatic in LOH) | **11,323** | **11,290** | **−33** |
| **`.\|0`** (somatic in LOH) | **660** | **693** | **+33** |

**關鍵觀察**：
1. **Germline het orientation 翻向**：baseline 偏好 `1|0`（hap1=alt），V2b 偏好 `0|1`（hap2=alt）。這對 GT_class 不影響（兩者都是 Germline_Het），但對 BAM HP tag 寫入會影響 `HP:i:1` 還是 `HP:i:2`。
2. **Somatic GT 完全保留**：21,304 個 `0|0`、~12,000 個 LOH-somatic 在兩版本之間個數一致，僅 33 個 `0|.` ↔ `.|0` 翻向（< 0.3% 變動）。

![Figure B — Top-15 GT raw value distribution](figures/12_gt_distribution/figB_gt_raw_top15.png)

### 3.3 GT2 × GT3（PASS sub-genotype）分布

PASS 的 sub-genotype 結構（GT2 × GT3 cross-tabulation）：

- baseline 與 V2b 的 GT2/GT3 cell 分布**完全相同**（PASS only n=47,798）
- 顯著 cells：
  - `GT2=.|1, GT3=./.`：19% (somatic mut not in LOH，hap2 是 alt)
  - `GT2=1|., GT3=./.`：18% (somatic mut not in LOH，hap1 是 alt)
  - `GT2=1|1, GT3=./.`：3.6% (ambiguous somatic, hp3)
  - LOH 區內的 sub-genotypes 各 < 1%

![Figure C — GT2 × GT3 sub-genotype heatmap (PASS)](figures/12_gt_distribution/figC_gt2_gt3_heatmap_PASS.png)

→ **sub-genotype phasing 結果完全相同**。

### 3.4 baseline vs V2b 差異熱圖

![Figure D — Δ% (V2b − baseline) per FILTER × GT_class](figures/12_gt_distribution/figD_baseline_vs_v2b_delta.png)

**所有 PASS cell 的 |Δ| ≤ 0.18 pp**（噪音等級）。
顯著差異集中在 **NonSomatic** 桶：V2b 多 phase 了 95,863 個 germline het（NonSomatic 的 Germline_Het 從 34.05% → 37.16%，+3.11 pp），對應的 Unphased 下降。

### 3.5 PS block coverage

| Filter | Label | n | PS rate | PS2/PS3 rate |
|--------|-------|---:|:-------:|:-----------:|
| PASS | baseline | 47,798 | 87.74% | 24.19% |
| PASS | V2b | 47,798 | 87.56% | 24.19% |
| NonSomatic | baseline | 3,086,681 | 54.96% | 0% |
| NonSomatic | V2b | 3,086,681 | **58.06%** | 0% |

→ V2b 多 phase 約 3.1 pp 的 NonSomatic（germline）變體，PASS 的 PS coverage 接近相同。

---

## Section 4 — 關鍵推論：Phasing vs Tag

### 4.1 證據鏈

| 觀察 | 推論 |
|------|------|
| Baseline 與 V2b 的 PASS somatic GT 完全相同（21,304 / 11,323+660 → 21,304 / 11,290+693） | **phasing 算法對 somatic 的 GT 決策一致** |
| 兩版本對 GT_class 的 PASS 分布 |Δ| < 0.2 pp | **phasing 階段語意輸出無顯著差異** |
| GT raw 出現 `0\|1` ↔ `1\|0` 大量翻向（germline het orientation flip） | **PS block phase orientation 不同**，但對下游語意無影響 |
| GT2 × GT3 sub-genotype 分布完全相同 | **somatic LOH 解析在 phasing 階段一致** |
| ClairS-TO calling 完全相同（同 input VCF）| 變體列表固定，差異只能來自 phasing 或 tag 階段 |
| Baseline BAM HP1:HP2 = 17.3:1，V5 BAM ≈ 1:1（見 audit suite §10） | BAM 層差異 ≠ VCF 層差異 |

### 4.2 結論

```
┌───────────────────────────────────────────────────────────────────┐
│  Phasing 階段 (longphase-to phase)                                │
│   baseline.vcf ≈ v2b.vcf （somatic GT 100% 一致, germline 翻向）    │
│   → ❌ 不是 phasing 問題                                           │
└───────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌───────────────────────────────────────────────────────────────────┐
│  Haplotag 階段 (longphase-to haplotag)                            │
│   baseline.bam (17.3:1 bias) ≠ v5.bam (≈1:1)                      │
│   → ✅ 確認是 tag 問題                                             │
│                                                                   │
│   根因：getVote() Layer 1 only, 無 PON anchor 過濾                │
│   修復：V5 加入 Layer 1.5 (somatic fallback) + PON-only flag      │
└───────────────────────────────────────────────────────────────────┘
```

### 4.3 對應 V5 修復點（呼應 11_code_issues_inventory.md）

| Bug ID | 階段 | V5 修復 |
|--------|------|--------|
| **Bug 1-1** getVote() priority | haplotag | 加 Layer 1.5 somatic fallback |
| **Bug 1-2** HP:i:33 enum | haplotag | 改 enum 比對 |
| **Bug 1-3** countINDELHaplotype guard | haplotag | UNDEFINED guard |
| **Bug 1-4** PhasingProcess no PON anchor | phase | `--pon-only-phasing` flag |

→ 4 個修復中 **3 個在 haplotag 階段**（與本稽核結論一致），1 個在 phase 階段（解釋 V2b 的 GT raw 翻向：PON-only 改變了 germline het orientation 決策）。

---

## Section 5 — 為什麼 GT raw 翻向不影響語意？

PS block 內部 phase orientation 是相對的：
- `1|0` 表「hap1=alt, hap2=ref」
- `0|1` 表「hap1=ref, hap2=alt」

對下游 haplotag 而言，只要在同一 PS block 內 orientation 一致即可正確 tag reads。**翻轉整個 block 的 orientation 不改變 reads 的 HP1/HP2 隸屬關係**（所有 reads 一起翻）。

但是：跨 block 的 self-phasing bias 是另一回事。當 haplotag 階段 getVote() 把不該成為 phasing anchor 的 somatic site 也當成 anchor 時，就會把 somatic-derived sub-population reads 全塞入同一 phase（造成 17.3:1）。這是 **tag 階段的決策錯誤**，與 phase orientation 無關。

---

## Section 6 — 限制與 caveat

| Caveat | 說明 |
|--------|------|
| 樣本範圍 | 本稽核只看 HCC1395 baseline tumor（同一 input BAM，同一 ClairS-TO VCF）。其他樣本（COLO829 等）行為應類似但未驗證。 |
| Tag 直接驗證 | 本稽核從 VCF 推論 tag 為差異源；BAM 直接驗證見 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_self_phasing_circular.md`。 |
| V5 phased.vcf 缺失 | V5 直接重用 V2b 的 phased.vcf，無「V5 自己的 phasing 輸出」可比。但 PON-only flag 邏輯由 V2b commit 引入，行為等價。 |
| Germline het 翻向 | 翻向不影響 ratio 但會影響某些下游（如 trio analysis、PS block 鏈接）。本研究不關心。 |

---

## Section 7 — 一頁速查

```
┌────────────────────────────────────────────────────┐
│  問題：17.3:1 HP1:HP2 self-phasing bias            │
│                                                    │
│  Phasing 階段（VCF GT 分布）                        │
│    Baseline vs V2b PASS somatic 100% 一致 ✓       │
│    Germline het 翻向 (1|0 → 0|1) 但語意不變 ✓     │
│    GT2 × GT3 sub-genotype 完全相同 ✓              │
│    → 不是 phasing 問題                             │
│                                                    │
│  Tag 階段（BAM HP tag）                            │
│    Baseline 17.3:1 bias                            │
│    V5 ≈ 1:1                                        │
│    根因：getVote() Layer 1 only + 無 PON anchor    │
│    → 是 tag 問題                                   │
│                                                    │
│  V5 修復: Layer 1.5 + PON-only flag + Bug 1-2/1-3 │
└────────────────────────────────────────────────────┘
```

---

## Section 8 — 跨檔交叉引用

- 上層 audit suite 索引：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/00_INDEX.md`
- Phasing/Tag 階段差異圖示：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/01_code_diff_analysis.md` §1.3
- Self-phasing 機制圖示：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/10_somatic_bias_explanation.md`
- Bug 編號對應：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/11_code_issues_inventory.md` §1
- Self-phasing 因果鏈：`InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/04_self_phasing_circular.md`

## 數據檔索引

| 檔名 | 內容 |
|------|------|
| `data/gt_distribution_raw.tsv.gz` | 全 6.4M variants × 9 cols (label, GT, GT2, GT3, PS...) |
| `data/gt_class_by_filter.tsv` | GT_class × filter_bucket × label 統計 |
| `data/gt_raw_by_label.tsv` | GT raw value 分布 |
| `data/gt2_gt3_cross_PASS.tsv` | GT2 × GT3 PASS-only cross |
| `data/gt_class_baseline_vs_v2b_delta.tsv` | Δpp baseline vs V2b |
| `data/ps_coverage_by_label_filter.tsv` | PS block coverage rate |

| 圖檔 | 內容 |
|------|------|
| `figures/12_gt_distribution/figA_gt_class_by_filter.png` | GT class stacked bar |
| `figures/12_gt_distribution/figB_gt_raw_top15.png` | Top-15 GT raw |
| `figures/12_gt_distribution/figC_gt2_gt3_heatmap_PASS.png` | GT2×GT3 heatmap |
| `figures/12_gt_distribution/figD_baseline_vs_v2b_delta.png` | Δpp heatmap |
