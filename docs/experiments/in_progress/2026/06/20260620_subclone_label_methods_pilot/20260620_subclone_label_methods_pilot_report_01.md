---
title: 三個甲基↔標籤關聯偵測方法 — BRCA2 pilot 結果報告
date: 2026-06-20
status: in_progress (BRCA2 pilot, 單樣本 HCC1395)
report_type: results_report (method comparison pilot)
tier: L2 (單樣本 × BRCA2 5 region × 單 pipeline; 偵測非驗證)
scope_flag: ⚠ PARTIAL — 單樣本 HCC1395 + BRCA2 5 region (想法1/2); 想法3 含全基因組 29,754 region 分佈(從既有 CSV)。跨樣本(COLO829)未做。
audience: 廖子游 · subclone 主軸
data_sources: _assets/data/idea1_summary.json, idea1_crosscheck_brca2.json, idea1_per_cpg_axis_attribution.tsv, idea3_fullscope_venn.json, idea3_pilot_dbeta.tsv, idea2_summary.json, idea2_concordance.json, idea2_dss_out_22305.tsv
provenance_note: 所有數字由 _assets/scripts/{run_idea1,run_idea3,run_idea2,run_idea2_compare}.py + audit_B_dss.R 跑出落 JSON/TSV，本報告 Read 回後撰寫(§13.0 撰寫與分析分批)。
---
<!-- provenance-verified: 全部 metric 引自 _assets/data/*.json (本 session 跑出讀回); modkit 07 區域值引既有 benchmark/07_benchmark_results.md。 -->

# 三個甲基↔標籤關聯方法 — BRCA2 pilot 結果

> ⚠ **PARTIAL / L2**：單樣本 HCC1395、BRCA2 5 region（想法1/2）；想法3 另含全基因組 29,754 region 分佈（讀既有 CSV）。**characterization 非 TP/FP 判別**；跨樣本未做。

> 🔴 **2026-06-21 散度 caveat（補，全 scope 已定論）**：本報告（及全 scope Venn）用的 label-PERMANOVA **未控 dispersion** —— 20260420 run 的 PERMDISP 欄全 = 1.0 stub（舊 binary，早於 C++ `check_dispersion`）。後處理 analytic betadisper 全 scope（27,304 label-sig region）：**「標籤顯著」中 ~34-36% 是 dispersion-confounded（過於樂觀）/ ~64% 是真 location 差異**，**HP 軸 34.4% / allele 軸 36.1% 相近**（⚠ BRCA2 pilot 看似「allele 乾淨」是 n=3 巧合，全 scope 推翻）。→ over-optimism 約 **1/3**（非 L3 估的 72%）；任何 label-PERMANOVA 顯著數字須附 location-clean 比例。詳 `docs/plans/20260620_subclone_situation_verification_reasoning_01.md §11.6`。⚠ BERNOULLI 非歐 → PERMDISP 本身近似。

## L0 — TL;DR（3 層）

1. **核心結論**：三法是「**甲基↔標籤關聯**」的 granularity ladder，pilot 證實它們**偵測不同層次的訊號、非同義**——per-CpG（想法1/2）抓局部訊號、region-mean Δβ（想法3）稀釋掉局部訊號、read-distance（PERMANOVA）抓 haplotype-coherence 且近飽和。
2. **為何**：BRCA2 ASM 訊號**局部集中在少數 CpG**（22305: 18/229 subclone-sig CpG）→ 想法1 抓到（方向 −0.106 符真值 −0.122）；想法3 把整窗 229 CpG 平均 → 稀釋成 Δβ=0.022（n.s.）**漏掉**；距離法 PERMANOVA 顯著但「到處都顯著」（91.8%）。
3. **可行性裁決**：三法**都能用既有 per-region 矩陣 Python 跑，零重跑**；想法1/2 證實可行且互補；想法3 證實是距離法的**保守子集**（差集近 0）。modkit/DSS 與我方在強訊號位點高度一致（Jaccard 0.76），但**二值化/無重複各有失真**。

---

## §1 Pilot 設計與資料（L1）

- **資料**：20260420 HCC1395 paired_full run，BRCA2 5 個 TP region（canonical=**22305** chr13:32,315,128 G>A）。
- **read 組成（22305）**：101 reads = 56 tumor〔HP1 germline 18 + HP1-1 carrier 26 + HP2 12〕+ 45 normal。
- 🔴 **normal reads 全 hp=0（unphased）** → normal-anchored ASM-control（chr2 §A5 式）**無法從 ISM 標準輸出做**（5/5 region 皆 `normal_control_evaluable=False`）。要做需 bespoke normal phasing（`research/tsg_promoter_asm_reviewer` step3/4 有）。
- 全程 read-set 一致（同矩陣），隔離「統計量本身」差異。

---

## §2 想法1 — per-CpG 三軸歸屬（L2；方法可行 ✓）

每 region 窗內逐 CpG × {HP-family / HP-fine / subclone / allele} MWU/Kruskal + 各軸獨立 BH-FDR。

| RegionID | SNV | CpG | HP-fam sig | subclone sig | allele sig | normal-ctrl |
|---|---|--:|--:|--:|--:|---|
| **22305** | chr13:32315128 | 229 | 7 | **18** | 4 | ✗ unphased |
| 22306 | chr13:32317522 | 183 | 4 | **65** | 28 | ✗ |
| 22307 | chr13:32324831 | 91 | 0 | 0 | 0 | ✗ |
| 22308 | chr13:32339132 | 78 | 0 | 1 | 1 | ✗ |
| 22309 | chr13:32345630 | 86 | 0 | 0 | 0 | ✗ |

- **真值 cross-check（22305）✓**：HP-family sig CpG（n=7）平均 Δβ = **−0.106**（與 bespoke 全域真值 **Δβ=−0.122 同方向同量級**）；promoter proximal（−800..+200bp）平均 Δβ=−0.071。
- **dominant_axis（22305）**：none 207 / subclone 16 / HP-family 4 / allele 2 → **多數 CpG 無訊號（合法）**，有訊號者 subclone 軸主導。
- ⚠ **allele 軸未控 germline-het baseline null**（deferred）→ allele 結果暫標「未控 baseline confound」。
- **輸出表**：`idea1_per_cpg_axis_attribution.tsv`（667 CpG×行，含各軸 p/q/delta + dominant_axis + dist_to_snv）。

**裁決**：方法**可行**——per-CpG 能定位「哪些 CpG、與哪軸相關」，22305 重現 BRCA2 ASM。但僅 22305/22306 有訊號（BRCA2 5 region 中 2 個），符合「BRCA2 是手挑 showcase」。

---

## §3 想法2 — modkit/DSS 逐 CpG vs 我方 read-distance（L2）

subclone 軸（HP1 germline vs HP1-1 carrier），四法吃同一矩陣（modkit-style 二值化 Fisher / DSS beta-binomial / 我方連續 MWU / read-distance PERMANOVA）。

| Region | modkit-Fisher | DSS | 我方 MWU | all-3 重疊 | Jaccard(modkit,我方) | Jaccard(DSS,我方) | read-dist PERMANOVA |
|---|--:|--:|--:|--:|--:|--:|--:|
| **22305** | 19 | **52** | 18 | 16 | **0.762** | 0.296 | p=0.01 ✓結構顯著 |
| 22306 | 20 | 41 | **65** | 17 | 0.308 | 0.293 | p=0.34 |

**發現**：
1. **強訊號位點（22305）modkit ≈ 我方**（Jaccard 0.76，19 vs 18）→ 二值化率差與連續 MWU 在乾淨二態訊號上相當。
2. **漸進訊號（22306）我方連續 MWU 多抓 3×**（65 vs 20，our-only 38）→ **二值化（modkit @ β>0.5）丟掉漸進甲基訊號**，read-level 連續保留。
3. **DSS 過度敏感**（22305: 52 vs 18-19；dss-only 33）→ **無生物重複**下 smoothing=TRUE 跨 CpG 借力 = anti-conservative（audit_B_dss.R 自承 no-replicate）。
4. **真 modkit 工具（既有 07，pooled 全 reads）**：BRCA2 somatic 軸 effect_size=**−0.159** vs ISM Δβ=**−0.122**（同向同量級，跨工具交叉驗證）。
5. **位點級 2×2**：22305 結構顯著(0.01)×有 sig CpG = 一致；22306 HP-family PERMANOVA n.s.(0.34) 但 subclone per-CpG 有訊號（⚠ 不同軸，非乾淨 2×2，因 significance.json 無 subclone-軸 PERMANOVA）。

**裁決（優缺點）**：
- **modkit/DSS 優**：CpG 解析度邊際 call、成熟工具、原生 effect-size+CI。
- **modkit/DSS 缺**：二值化丟漸進訊號（modkit）；需重複才有效推論（DSS 無重複過 liberal）；**皆丟 read 內 haplotype-coherence**（pooled 邊際率差）；皆無 cis-control。
- **我方 read-distance 優**：保留 read-level 連續 + haplotype-coherence + normal-anchored cis-test（07 證 ISM 獨有）。**缺**：給 region 級結構非 per-CpG（per-CpG 需想法1 補）。
- → **不是「誰更好」，是分工**：per-CpG（想法1/modkit/DSS）定位 CpG；read-distance 抓 coherence + cis 隔離。

---

## §4 想法3 — per-label Δβ + 集合 Venn（L2 pilot + 全 scope 分佈）

### 4.1 pilot（per-read region-mean Δβ，permutation p）
| Region | HP Δβ(p) | subclone Δβ(p) | allele Δβ(p) |
|---|---|---|---|
| **22305** | −0.027 (0.42) | 0.022 (0.53) | −0.042 (0.19) |
| 22306 | 0.070 (0.19) | **0.125 (0.008)** | −0.072 (0.10) |

🔴 **關鍵**：22305 想法3 三軸**全不顯著**——region-mean 把局部 ASM 訊號（想法1 抓到 18 個 subclone-sig CpG）**稀釋掉**。→ 回答你「偵測更多嗎」：**想法3 偵測更少、更保守**（region-mean 對「局部 pattern」天生盲）。

⚠ **cross-check 揪出口徑差**：Python per-read 想法3 subclone Δβ=0.022 vs CSV 既有 `HPFineD_HP1_HP1S`=**0.236**（差 ~10×）→ **CSV 既有欄 ≠ 想法3 字面定義**（CSV per-CpG-aware，想法3 region-mean）。**full-scope 用 CSV 為 proxy；exact 想法3 更保守**。

### 4.2 全基因組 Venn（29,754 region，從既有 CSV，proxy）
**set B 兩定義並排**（你問「是否完全包含距離切顯著」）：
- **set B = label-PERMANOVA**（HP 或 allele）= **27,304（91.8%）** ← near-saturated（結構到處有）
- set B = GlobalP/對齊 = 2,520

| τ | set A(Δβ) | A∩B | **A−B** | B−A | Jaccard | 對 GlobalP: A−B / B−A |
|--:|--:|--:|--:|--:|--:|--|
| 0.10 | 7,333 | 7,331 | **2** | 19,973 | 0.27 | 6,373 / 1,560 |
| 0.15 | 3,910 | 3,908 | **2** | 23,396 | 0.14 | 3,253 / 1,863 |
| 0.20 | 2,367 | 2,366 | **1** | 24,938 | 0.09 | 1,901 / 2,054 |
| 0.25 | 1,595 | 1,595 | **0** | 25,709 | 0.06 | 1,247 / 2,172 |

🔴 **答案**：**想法3（Δβ法）幾乎完全被距離法（label-PERMANOVA）包含**（A−B = 0~2，τ=0.25 時完全子集），距離法多偵測 ~7×。**但對 GlobalP/對齊集則有真實雙向差集**（各 ~2,000-6,000）→ **答案取決於「距離法」用哪個定義**（label-PERMANOVA 近飽和 vs GlobalP/對齊嚴格）。

### 4.3 差集案例（τ=0.15）
- **A−B（Δβ 抓到、距離沒切）= 僅 2 個**：`chr8:118671269`、`chr5:175081627`，**都 Potential_LOH=True + Fisher_N_Sig=0**（無 per-CpG driver）→ **LOH 區純整體位移**（level-shift），無局部 pattern 故 PERMANOVA 不顯著。
- **B−A（距離切出、Δβ 沒過）**：top 由高 PERMANOVA F（47-140）但 HPMergedDelta 小（0.08-0.14）、Fisher_N_Sig 高（6-65）→ **局部 per-CpG 結構**，region-mean Δβ 抓不到（正是想法3 盲點）。

### 4.4 LEVEL-shift confound
HP-Δβ-sig 728 個中 **134（18.4%）無 per-CpG driver**（Fisher_N_Sig=0）= 純整體位移非 per-CpG pattern。

**裁決**：想法3 = 距離法的**保守子集**；偵測「整體 level-shift」、漏「局部 pattern」；18% 命中是純位移。**與 tumor-only 非監督 NEGATIVE 區隔**（這是 a-priori label 軸 Δβ，非無監督切群）。

---

## §5 三想法統一比較

| 軸 | 想法1 (per-CpG) | 想法2 (modkit/DSS) | 想法3 (region-mean Δβ) | read-distance (我方) |
|---|---|---|---|---|
| 粒度 | per-CpG | per-CpG (邊際率差) | per-region (1 數) | per-region 結構 |
| 22305 抓到? | ✓ 18 subclone CpG | ✓ 19 (modkit) | ✗ 稀釋(n.s.) | ✓ PERMANOVA 0.01 |
| 對漸進訊號 | ✓ 連續 | modkit ✗二值化 / DSS borrow | ✗ 平均洗掉 | ✓ 連續距離 |
| haplotype-coherence | ✗ 邊際 | ✗ 邊際 | ✗ 平均 | ✓ read×read |
| 全 scope 關係 | — | — | **⊂ 距離法**(A−B≈0) | 近飽和 91.8% |

**統一收斂**：
- **(a) 位點層**（想法3 vs 距離）：想法3 ⊂ 距離法（差集近 0）。
- **(b) CpG 層**（想法1 vs 想法2）：強訊號高度一致（Jaccard 0.76）；漸進訊號連續法多抓 3×；DSS 無重複過 liberal。
- **(c) 軸一致**：subclone 軸是 BRCA2 最強訊號軸（想法1 18 sig、modkit −0.159、想法3 22306 sig）；HP/allele 軸較弱。

---

## §6 誠實邊界 + tier + 待補

- **tier = L2**（單樣本 × BRCA2 5 region × 單 pipeline；BRCA2 是 genome-wide TSG-enrichment=0 的手挑 showcase，不可外推）。
- **characterization 非判別**：全程未宣告 TP/FP 判別力 / filter（ASM B-discrimination 已 NEGATIVE）。
- **待補（deferred）**：
  1. **normal-anchored ASM-control**（normal unphased → 需 bespoke phasing；想法1 #4 的 clean/confounded 拆解未能跑）。
  2. **allele germline-het baseline null**（想法1/3 的 allele 軸未控 cis-ASM baseline）。
  3. **exact 全-scope 想法3**（per-read region-mean，非 CSV proxy；CSV 口徑差已揪出）。
  4. **跨樣本**（COLO829）破單樣本 → ⭐4。
  5. **真 modkit/DSS pooled 全-reads** 全 panel head-to-head（pilot 用同矩陣隔離統計量；真工具 read-set 不同）。

---

## Provenance
- 數字 SoT：`_assets/data/{idea1_summary,idea1_crosscheck_brca2,idea3_fullscope_venn,idea2_summary,idea2_concordance}.json` + `idea1_per_cpg_axis_attribution.tsv` + `idea2_dss_out_*.tsv`。
- 腳本：`_assets/scripts/{lib_region,run_idea1,run_idea3,run_idea2,run_idea2_compare}.py` + DSS via `benchmark/runs/audit_B_dss.R`（env 自帶 Rscript）。
- 真 modkit 07：`docs/method_comparison/20260609_ism_vs_external_methylation_tools/benchmark/07_benchmark_results.md`（BRCA2 −0.159 vs ISM −0.122）。
- 推論層 SoT：`InterSubMod/docs/plans/20260620_subclone_situation_verification_reasoning_01.md`。
