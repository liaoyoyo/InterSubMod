---
title: 全基因組 (Whole-Genome) Per-Site HP + Paired Concordance Audit
date: 2026-04-29
author: liaoyoyo2001
tags: [audit, whole-genome, paired-concordance, v5, baseline, scale-up]
status: validated_complete
audience: PI + reviewers
related_docs:
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md
  - InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/19_per_site_nuance_audit.md
---

# 全基因組 Per-Site HP + Paired Concordance Audit

## §0 一句話結論

> **全基因組 (~90K PASS sites × 5 BAMs) 數據揭露：V5 在 aggregate HP1:HP2 ratio 上比 baseline 接近 paired truth (V5_93=1.42:1 vs BL_93=2.08:1，PA_93=1:1.08)，但在 site-by-site paired concordance 上 V5 vs BL 幾乎相同 (48% vs 49%，皆接近 50% random level)，再次證實 19 號 nuance audit 的結論「V5 改善是 aggregate-level 不是 site-by-site」**。

---

## §1 分析設計

### 1.1 Why Whole-Genome？

既有 audit suite 多用 15 audit sites + chr19-22 sample，**全 PASS sites 從未跑過完整對比**。本檔將分析尺度提升到 whole-genome 層級。

### 1.2 5 BAMs × 全 PASS sites

| Sample | BAM | n PASS sites |
|--------|-----|:-:|
| **PA_93** | `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam` (paired tagged) | 47,798 |
| **BL_93** | `/big7_disk/.../baseline/tumor_tagged.bam` | 47,798 |
| **V5_93** | `/big7_disk/.../pononly_v5_somatic_fallback/tumor_tagged.bam` | 47,798 |
| **BL_06** | `/big7_disk/.../baseline_06/tumor_tagged.bam` | 42,770 |
| **V5_06** | `/big7_disk/.../v5_06/tumor_tagged.bam` | 42,770 |

ClairS-TO PASS sites 從各自 sample 的 VCF 取（whole genome chr1-22+X+Y）。

### 1.3 多執行緒設計

```
Per BAM:
  Pool(24 workers) × 24 chromosomes
  Each worker: own pysam.AlignmentFile handle
5 BAMs sequential (avoid IO contention on shared NFS/local disks)
```

實際耗時：
- BAM index loading: ~30s/BAM
- pysam.pileup per site: ~10-50ms
- **Total: 7,149s ≈ 2 小時**（5 BAMs × ~1,300s avg）

---

## §2 Whole-Genome Aggregate Results

### 2.1 HP1:HP2 ratio (somatic ALT-only reads)

| Sample | n_sites | total_alt_reads | HP1 | HP2 | HP33 | **HP1:HP2 ratio** |
|--------|:-:|:-:|:-:|:-:|:-:|:-:|
| **PA_93** | 47,798 | 2,344,492 | 325,505 | 350,723 | 0 | **1:1.08** ✓ ≈ 1:1 |
| **BL_93** | 47,798 | 2,344,492 | 1,497,129 | 719,367 | 19,774 | **2.08:1** |
| **V5_93** | 47,798 | 2,344,492 | 1,218,134 | 856,769 | 124,834 | **1.42:1** ← 更接近 PA |
| **BL_06** | 42,770 | 702,113 | 273,269 | 393,102 | 11,642 | **1:1.44** |
| **V5_06** | 42,770 | 702,113 | 166,151 | 412,931 | 76,257 | **1:2.49** |

![WG aggregate summary](figures/20_whole_genome/wg_summary.png)

### 2.2 重要對比 — V5 vs BL 距離 PA 的程度

定義 distance = |log10(ratio)|（PA = 0 為原點）：

| Sample | log10(ratio) | distance to PA |
|--------|:-:|:-:|
| **PA_93** | -0.034 | 0 (基準) |
| BL_93 | +0.318 | **0.352** |
| **V5_93** | +0.152 | **0.186** ← 比 BL 接近 PA 0.166 (47% closer) |
| BL_06 | -0.158 | 0.124 |
| V5_06 | -0.396 | 0.362 ← 比 BL_06 偏離 PA |

→ **V5 @ 0.93 比 BL @ 0.93 接近 PA truth**（aggregate level）
→ **V5 @ 0.6 反而比 BL @ 0.6 更偏離 PA**（HP33 conservative tagging 副作用）

### 2.3 HP family composition

| Sample | HP1% | HP2% | HP33% | untagged% |
|--------|:-:|:-:|:-:|:-:|
| **PA_93** | 14% | 15% | **0%** | **71%** ⚠ |
| BL_93 | 64% | 31% | 0.8% | 5% |
| V5_93 | 52% | 37% | 5% | 6% |
| BL_06 | 39% | 56% | 1.7% | 3% |
| V5_06 | 24% | 59% | **11%** | 7% |

⚠ **重要發現**：**PA_93 71% reads 是 untagged**（HP:i:0）。原因：paired phasing 只用 normal germline anchor，多數 PASS somatic 周圍 PS block 太弱/短，無法 tag。**這是 paired ground truth 的根本限制**。

→ V5 / BL 在 **aggregate 層級** 都比 PA tag 多很多 reads（24-64% HP1 vs 14% in PA），但這不必然代表更精準。

---

## §3 Site-Level Paired Concordance

### 3.1 對比方法

對每個 PASS site，定義 dominant HP family（HP1/HP2/HP33/tie），然後計算 BL_93/V5_93/BL_06/V5_06 的 dominant 與 PA_93 dominant 的匹配率。

只在 PA_93 dominant ∈ {HP1, HP2} 的 sites 上計算（因 HP33 在 PA_93 是 0，無法當 ground truth）。

### 3.2 結果

| Comparison | n_sites_evaluable | match | match% | opposite | opposite% | hp33_conservative | hp33% |
|------------|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| **BL_93 vs PA_93** | 12,657 | 6,148 | **48.57%** | 6,337 | 50.07% | 172 | 1.36% |
| **V5_93 vs PA_93** | 12,609 | 6,060 | **48.06%** | 6,065 | 48.10% | 484 | 3.84% |
| **BL_06 vs PA_93** | 7,328 | 3,797 | **51.81%** | 3,407 | 46.49% | 124 | 1.69% |
| **V5_06 vs PA_93** | 7,240 | 3,566 | **49.25%** | 3,284 | 45.36% | 390 | 5.39% |

### 3.3 結論：site-by-site paired concordance ≈ random

**全部 4 個 V5/BL 場景 vs PA_93 的 match% 落在 48-52%**，接近 random level (50%)。

| 對比 | V5 vs BL 改善 |
|------|:----:|
| @ 0.93: BL 48.57% vs V5 48.06% | **−0.51pp** ❌ V5 略差 |
| @ 0.6: BL 51.81% vs V5 49.25% | **−2.56pp** ❌ V5 略差 |

→ **V5 在 site-by-site paired concordance 上沒有改善**（甚至略差）。

### 3.4 為什麼 audit suite §07 報的 +6.65pp / +13.3pp 不見了？

| 報告 | metric | denominator | match% |
|------|------|------|:-:|
| §07 paired_ground_truth | conditional concordance | 雙方都 directional tag | +6.65pp / +13.3pp |
| **本檔 (§3)** | **site dominant**（含 HP33）| 全 evaluable PASS sites | **−0.5 ~ −2.5pp** |

→ **+6.65pp / +13.3pp 是 conditional accuracy**（V5 evaluable denominator 較小）。**全 PASS sites scale 看，V5 paired concordance 沒改善**。這完全證實 19 號 nuance audit §4.2 的「fixed-denom 重算後幾乎持平」結論。

---

## §4 V5 改善的真實樣貌（重新定位）

### 4.1 V5 改善的 metric

| Metric | V5 vs BL @ 0.93 | V5 vs BL @ 0.6 |
|--------|:-:|:-:|
| **Aggregate HP1:HP2 ratio distance to PA** | ✅ **-47%** (0.352 → 0.186) | ❌ +192% (0.124 → 0.362) |
| **Aggregate HP33% (conservative tagging)** | ↑ (0.8% → 5%) | ↑ (1.7% → 11%) |
| **Site-level paired concordance match%** | ❌ −0.51pp | ❌ −2.56pp |
| Untagged% | similar | similar |

### 4.2 重新定位

V5 的真實價值：
- ✅ **修復 aggregate self-phasing imbalance**（17.3:1 → ~1.5:1，全 ALT reads scale；本檔 BL_93 2.08:1 → V5_93 1.42:1）
- ✅ **Conservative tagging**（HP33% 顯著上升）
- ❌ **Site-level paired correctness 未改善**（與 19 號 nuance audit 一致）
- ⚠ **0.6 下 V5 反而偏離 PA**（HP33 conservative 過度推到 HP2，稀釋 HP1 votes）

### 4.3 對既有 audit suite 結論的進一步精確化

| 既有結論 | WG-level 數據 | Nuance |
|---------|---------|--------|
| V5 修復 17.3:1 bias | ✅ 全 ALT reads scale: 2.08:1 → 1.42:1 (WG)；614K:35.5K (legacy subset) | 17.3:1 是 high-AF subset metric |
| Aggregate paired +6.65pp | ⚠ Conditional only；site-level **−0.51pp** | conditional 偏 V5（HP33 排除）|
| Clean PS +13.3pp | ⚠ 同上，clean PS subset 內 conditional 改善 | site-level 不改善 |
| V5 conservative tagging | ✅✅ 全 WG 確認（HP33% 0.8% → 5%）| 0.6 下更顯著（11%）|

---

## §5 Paired Ground Truth 的根本限制

### 5.1 71% Untagged 的問題

PA_93 在 47,798 PASS sites 上有 1,668,264 untagged reads（71% of total ALT reads）：
- Paired phasing **只用 normal germline anchor**
- 多數 PASS somatic 在 normal 的 PS block 弱/短 → reads 無 HP tag

→ PA_93 不是 universal ground truth；它是「**germline-phasing-anchored ground truth**」。

### 5.2 這個限制的影響

- 49.25% match% **不代表 V5 49% 對 51% 錯**
- 真實情況：V5 在 PA_93 給出 directional tag 的 12,657 sites 上 ~50% 對 ~50% 錯
- 在 PA_93 untagged 的 sites 上，V5 vs BL **無法判斷對錯**

→ 「Site-level paired concordance」是 **conditional metric**，但 evaluable subset 已經是 PA_93 較有信心的 sites（germline-rich PS blocks），V5 在這些 sites 上沒有勝出，意味著 V5 的 phase orientation 選擇與 PA 並不一致。

### 5.3 對未來研究的影響

需要 **更強的 ground truth**：
- Trio sequencing（父母 + 子女）建立 read-level haplotype phasing
- 或 long-read assembly + variant calling 的 phased VCF
- 目前可用的 paired-tagged tumor BAM 僅是 weak truth

---

## §6 一頁速查

```
WG-level metrics (5 BAMs × ~90K sites total):

[Aggregate HP1:HP2 ratio (距離 PA = 1.08:1)]
  PA_93:  1:1.08   ←  truth
  BL_93:  2.08:1   |Δ| = 0.352
  V5_93:  1.42:1   |Δ| = 0.186  ← V5 改善 -47%  ✓
  BL_06:  1:1.44   |Δ| = 0.124
  V5_06:  1:2.49   |Δ| = 0.362  ← V5 偏離 +192% ⚠

[Site-level paired concordance match%]
  BL_93 vs PA_93:  48.57%  ←  random level
  V5_93 vs PA_93:  48.06%  ← V5 -0.51pp 略差
  BL_06 vs PA_93:  51.81%
  V5_06 vs PA_93:  49.25%  ← V5 -2.56pp 略差

[HP33 Conservative tagging%]
  BL_93: 0.8%, V5_93: 5%  (V5 顯著保守)
  BL_06: 1.7%, V5_06: 11% (V5 更保守)

[Untagged in PA_93]
  71% reads 無 paired tag  ← 限制 ground truth 可用性
```

---

## §7 結論與建議

### 7.1 V5 改善的真實定位

✅ **Aggregate self-phasing imbalance 修復** — 全 WG ALT reads HP1:HP2 從 2.08:1 → 1.42:1
✅ **Conservative tagging** — HP33% 從 0.8% → 5%（WG 級確認）
❌ **Site-by-site paired correctness 未改善** — match% 48% vs 49% (≈ random)
⚠ **0.6 下 V5 反而偏離 PA** — 因 HP33 過度推到 HP2

### 7.2 audit suite 內部結論一致性

| 報告 | 結論 | 與本檔 (WG) 一致？ |
|------|------|:----:|
| 07 (paired) | conditional +6.65/+13.3pp | ⚠ 是 conditional 的，site-level 不一致 |
| 09 (0.6 simulation) | V5 conservative tagging | ✅ 完全一致 |
| 14 (user report) | V5 修復 self-phasing | ✅ aggregate 一致 |
| 17 (design check) | V5 對齊 LongPhase-TO 設計 | ✅ 一致 |
| 19 (nuance) | 1 strong / 2 regression / 5 tie / 4 NA | ✅ **WG 級完全證實** |
| **20 (本檔)** | aggregate 改 / site-level 沒改 | — |

### 7.3 建議使用 V5 的場景

| 用途 | 推薦? |
|------|:-:|
| 全基因組 aggregate self-phasing 統計 | ✅ V5 |
| 個別位點 paired truth 對比 | ⚠ 兩者皆無顯著差別 |
| Conservative downstream（需要 HP33 顯式標）| ✅ V5 |
| 低 purity 樣本（< 0.7） | ⚠ V5 conservative 可能過度，paired 距離增加 |

---

## §8 數據與圖檔索引

### 數據

| 路徑 | 內容 | 大小 |
|------|------|:-:|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/wg_per_site_hp.tsv.gz` | 228,934 rows × 9 cols (long format per-site) | 1.7 MB |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/wg_summary.tsv` | 5 rows (per-sample aggregate) | 838 B |
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/data/wg_paired_concordance.tsv` | 4 rows (BL_93/V5_93/BL_06/V5_06 vs PA_93) | 454 B |

### 圖

| 路徑 | 內容 |
|------|------|
| `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/figures/20_whole_genome/wg_summary.png` | 2 panels: HP1:HP2 ratio + HP family composition |

### 配套腳本

| 路徑 | 用途 |
|------|------|
| `InterSubMod/scripts/analysis/v5_whole_genome_per_site_audit.py` | Multi-process WG audit (24 workers/BAM, 5 BAMs sequential) |

---

## §9 跨檔索引

| # | 文件 | 與本檔關係 |
|---|------|---------|
| 07 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` | conditional accuracy 來源（+6.65/+13.3pp）|
| 09 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/09_purity06_simulation.md` | 0.6 場景 + chr19-22 範圍 |
| 19 | `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/19_per_site_nuance_audit.md` | site-by-site nuance（13 audit sites）|
| **20** | **本檔** | WG-level (~90K sites) 完整證實 |
