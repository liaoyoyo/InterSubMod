---
title: "U1 — MAX_DIST vs SKIP nan-distance-strategy 對 subclone 分群品質比對"
date: 2026-06-14
type: experiment / methodology
status: chr1 pilot + 全基因組 confirmed (單樣本 HCC1395)
scope_flag: "⚠ 單樣本 HCC1395 paired; chr1(2624)+全基因組(30490 位點) 模式一致; 多樣本(COLO829 等)待擴"
binary: fix/clustering-nan-skip-segfault @ 517ed90 (SKIP segfault 修復後)
config: HCC1395 paired tumor/normal, window ±5000, BERNOULLI, -j16, 99-perm PERMANOVA
data_sources: comparison_summary.json, interpretation.json, vc_crosstab.tsv, region_full_compare.tsv
---

# U1：MAX_DIST vs SKIP 對 subclone 分群品質

> ⚠ **PARTIAL pilot**：chr1 TP SNV 2624 位點（單染色體 subset）；apples-to-apples 同配置只差 `--nan-distance-strategy`。

## L0 結論（謹慎）
SKIP 去除 MAX_DIST 的 **1.0 假距離污染**後，甲基分群訊號**大幅浮現**（Significant **91→1007**、Strong **389→1767**、CramersV_mean **0.017→0.230**），且 confound 已排除（**非小樣本、非 permutation 邊緣 → 真訊號**）。**但浮現的主導是 HP-axis germline haplotype 分群（DominantLabel hp 83–90%），非 somatic subclone** — 與主軸「甲基 = germline-haplotype 層級」一致。**對 somatic subclone 偵測力的淨提升尚需 normal 對照確認（記錄缺口 gap2）。**

## 全維度比對（chr1 2624 位點，數字 ← comparison_summary.json）

| 維度 | MAX_DIST | SKIP | 倍率 |
|------|---------:|-----:|:----:|
| **VerificationClass** Strong | 389 | **1767** | 4.5× |
| ‧ Weak | 1961 | 562 | 0.29× |
| ‧ Noise | 262 | 194 | 0.74× |
| ‧ Subclone | 12 | **101** | 8.4× |
| **Significant** | 91 | **1007** | 11× |
| PassedGating / PermanovaValid | 529 | **1918** | 3.6× |
| CramersV_mean | 0.0173 | **0.2304** | 13× |
| HeuristicScore_mean | 1.75 | **12.99** | 7.4× |
| NumReads / NumCpGs median | 118 / 78 | 118 / 78 | 不變（fetch 不受 strategy 影響）|

**Per-region 變化**：VerificationClass 改變 **1751** region（67%）、Significant 改變 **964**、任一改變 **1911**（73%）。

## VerificationClass 流向（row=MAX_DIST → col=SKIP，← vc_crosstab.tsv）
主流向 = **Weak(MAX_DIST) → Strong(SKIP) = 1403**（MAX_DIST 的 1.0 污染把真 Strong 壓成 Weak；SKIP 還原）。
其餘：Weak→Weak 433、Strong→Strong 332、Weak→Subclone 39、Weak→Noise 86。

## Confound 判別（← interpretation.json confound_checks）— 為何是真訊號非 artifact
- **非小樣本 artifact**：SKIP 新增 940 個 Significant region 的 NumReads median **116**（≈全體 118）→ 剔除 NaN read 後樣本仍正常大小。
- **非 permutation 邊緣**：SKIP 新 sig 的 PermanovaP median = **0.01**（= 99-perm 下限 1/100，observed F 超過所有 shuffle），且 CramersV_skip median = **0.573**（強關聯）→ 是強分群非弱證據。
- **結論**：差異是 SKIP 去污染後**真實甲基分群訊號浮現**，MAX_DIST 的 1.0 假距離稀釋了 PERMANOVA F 統計量。

## 🔑 本質判別（← interpretation.json）— 浮現的是 germline-HP，非 subclone
| 軸 | MAX_DIST sig (91) | SKIP sig (1007) |
|----|------------------:|----------------:|
| DominantLabel = **hp** | 82 (90%) | 835 (83%) |
| DominantLabel = allele | 9 | 172 |
| HPMergedSig（germline HP 分群）| 81 | 897 |
| AlleleSig | 60 | 637 |

→ 兩 strategy 的顯著分群**主導都是 HP-axis（germline haplotype）**，SKIP 只是讓這個 germline 甲基訊號從污染中浮現。**這不等於「更多 somatic subclone」** — 甲基的 germline-HP 結構在 normal 也存在。

## 🔴 記錄缺口（用戶訴求「足夠比對的資訊」— 2 項待補才能完整回答）
1. **gap1 — n_reads_valid 不在 summary**：SKIP 每個 region 實際剔除多少 read 無法從 `significance_summary.csv` 取得（只有全域 runlog：invalid pairs 3,489,658 / valid pair ratio 82% / avg common CpG 41.44）。per-region 剔除量缺 → 個別 region 的剔除強度無法追溯。**補法**：cpp-change 加 `NReadsValid` 欄到 summary（小改），或重跑存 distance matrix 另算。
2. **gap2 — germline vs somatic-subclone 分離**：要判別 1007 個 sig 中「真 tumor-specific subclone」比例，需 **normal-only HP 分群對照**（normal 也分同 HP = germline；僅 tumor 分 = subclone）。本 run 為 paired 合一，未產 normal-only 對照。

## ✅ 全基因組驗證（30490 TP 位點，數字 ← wg/comparison_summary.json）
chr1 pilot 模式**完全複現於全基因組**（單樣本 HCC1395），確認非 chr1 特異：

| 維度 | MAX_DIST | SKIP | chr1 對照(倍率) |
|------|---------:|-----:|:---:|
| Significant | 1276 | **11912** | 9.3× (chr1 11×) |
| Strong / Weak / Noise / Subclone | 5133/22305/2856/196 | **21015**/6122/2058/**1295** | Strong 4.1× (chr1 4.5×) |
| CramersV_mean | 0.0213 | **0.2383** | 11× (chr1 13×) |
| **sig_hp_fraction** | 0.851 | 0.853 | **85% germline-HP**(chr1 0.90/0.83) |
| VC changed / any changed | — | — | 19961(65%) / 21901(72%) |

**SKIP-new-significant confound 全基因組同樣 clear**：10946 個（lost_by_skip 僅 310），NumReads median **123**（非小樣本）、CramersV median **0.578**、PermanovaP median **0.01**（99-perm 下限）、PermanovaP<0.01 = **0**。→ 與 chr1 逐項一致 = 真訊號、germline-HP 主導，跨全染色體 robust。

## 🖼 案例熱圖視覺確認（figures/case_*.png，4 個明顯改善 region）
chr1 挑 4 個 Weak→Strong / CramersV 0→1.0 / DominantLabel=hp 案例（n_reads 83–98、NumCpGs 42–89、SKIP NaN 剔除 16–38%），畫距離熱圖（read×read，HP 排序）+ 甲基圖：
- **MAX_DIST 距離**：黃色 1.0 假距離污染**淹沒** HP block 結構。
- **SKIP 距離**：灰=NaN 剔除後，**HP1/HP2 對角 block 深藍（群內低距離）+ 跨群綠（高距離）清晰浮現**（黃虛線=HP 邊界）。
- **甲基圖**：HP1 vs HP2 群甲基模式分層（藍 低/紅 高）→ 證實是 **germline-haplotype 甲基差異**（read 按 hp tag 即分開）。
- **視覺結論**：確認「結果改善」= SKIP 還原被污染掩蓋的 **germline-HP 甲基分群**；非 somatic subclone（與全基因組 85% hp-dominant 一致）。
- 可重現：`plot_cases.py`（建小 VCF 重跑開 distance/methylation 輸出 → matplotlib）。

## 結論與下一步
- **技術層**：SKIP 是更正確的 NaN 處理（不用 1.0 假距離污染）；MAX_DIST 系統性低估甲基分群（之前「Noise 9.4% 運作良好」其實是污染下的保守表象）。
- **對 subclone 重建**：SKIP 浮現的主要是 germline-HP 訊號，對 somatic subclone 的淨增益**尚未確立** → 需補 gap2（normal 對照）。
- **建議優先序**：① 補 gap1（`NReadsValid` 欄，cpp-change 小改）讓記錄完整可追溯；② 補 gap2（normal-only HP 對照）才能回答「SKIP 對 somatic subclone 是否更好」。
- ⚠ 本結論為 **chr1 partial**；若推 subclone 結論需全基因組 + COLO829 等多樣本（主軸 6 樣本資產）。

## 產物
- `comparison_summary.json` — 總體分布 + change counts
- `interpretation.json` — 軸判別 + confound checks + record gaps
- `vc_crosstab.tsv` — VerificationClass 流向
- `region_full_compare.tsv`（2624 region）/ `region_diff_changed.tsv`（1911 變化 region）— per-region 可追溯
- `compare_u1.py` — 可重現比對 script
