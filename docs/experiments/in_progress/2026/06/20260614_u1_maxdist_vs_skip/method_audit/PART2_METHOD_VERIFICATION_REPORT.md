---
title: "Part 2 方法驗證報告 — 分群/標籤/關聯判定方法細節 + 邊緣案例 + 個案肉眼驗證"
date: 2026-06-15
type: method-verification (2A 方法audit + 2B 邊緣案例 + 2C 個案圖)
scope: B comprehensive — chr1(2624) 經驗 + 全基因組(30490)交叉; 單樣本 HCC1395 paired (partial flag)
binary: develop a8d8d61 (SKIP 預設 + Δβ 模組 stage1+2 + HP-AUC)
data_sources: dbeta_stage1_verify_wg.json, dbeta_stage2_verify_wg.json, dbeta_stage3_crosscheck.json, inspection_cases.json, /tmp/w4ukyquky.output(2A audit), /tmp/dbeta_chr1_s2/significance_summary.csv
claim_levels: 方法audit/edge=L1(源碼+6 fresh-context agent) / 經驗數字=L1(本run CSV) / 建議=L2
related: SIGNIFICANCE_INVENTORY / METHODS_BY_GOAL / TWO_DIRECTIONS_FRAMEWORK / COMPUTE_FLOW_ARCHITECTURE
---

# Part 2 方法驗證報告

> **L0 總 verdict**：**無翻轉結論的 bug**；分群(UPGMA)/PERMANOVA/Fisher/BH-FDR/NME-Epipoly/Δβ 核心計算**全部標準-正確**；
> 所有邊緣案例(tiny-N / all-same / single-HP / NaN-heavy / empty)**全 graceful**(回 NA/skip，無 crash/garbage)。
> 6 個 audit 全 SOUND 或 MINOR_GAP(無 NEEDS_FIX)。發現 **6 個可精修點**(rigor/hygiene 層，非錯誤)+ **1 個 subclone 統計脆弱性**(本 session 新發現)。

---

## 2A. 方法細節 audit（6 fresh-context agent 源碼層，read-only）

| 元件 | verdict | 嚴重度 | 標準（精確） | 關鍵 issue |
|------|:---:|:---:|------|------|
| **分群 UPGMA** | SOUND | low | size-weighted Lance-Williams average linkage（教科書正確，`HC.cpp:131`）；height=dist/2 單調 | build_complete max_dist 初值 0.0（COMPLETE 未用，dead）；tie-break 依 read 順序；Ward 是近似（未用）|
| **tree-cut + k 選取** | MINOR_GAP | **med** | k=**silhouette-optimal** in [2, min(6,n_valid/2)]（非固定 cut height、非 k=N）+ size-rebalance bump | ① tie-break 偏小 k（near-degenerate 距離→欠分辨子結構）② **cluster-PERMANOVA 缺 min_reads_per_group(3) gate**（label PERMANOVA 有）→ singleton cluster 可通過 |
| **標籤指派** | MINOR_GAP | low | HP=BAM `HP` tag（longphase-s 字串/-to 整數 1/2/11/21/33）；allele=SNV base vs VCF；sample=來源 BAM | 🔴 **stale 危險註解**：`RegionProcessor:810`+`ReadAggregator:46-48` 寫 normal hp="0"/REF **與碼相反**（normal 保留 germline HP，2509/2624 Normal_HP_Valid=true，germline baseline 依賴它）→ 未來「修正」會弄壞 |
| **gap#2 Dispersion p** | MINOR_GAP | low | F→p 查表（F>4→.01/>2.5→.05/else .1）非真 p | **production dead code**（`enable_dispersion=false@:1926`，全 2624 DispersionP=1.0）→ 今天誤導 0 下游；若重啟用是 foot-gun |
| **gap#4 PerCpgAsm Fisher** | MINOR_GAP | **med** | per-CpG 2×2 Fisher exact + BH-FDR（皆正確）+ NME/Epipoly（標準）| Fisher 把 read 當獨立→**over-dispersion 膨脹**（MaxNegLogFDR≤12.9、median Frac_Sig 0.19 過強）；**characterization-only 非 filter**→不直接撐判別 claim |
| **gap#5 perm 數** | (含上) | low | cluster PERMANOVA 99（floor 0.01）vs LabelTest 999（floor 0.001）| 1916/1918 cluster 釘在 0.01；benign(GlobalP 用 Fisher 非 PERMANOVA，無 live 比較混兩 floor)；建議統一 999 |
| **gap#6 冗餘** | SOUND | low | distance-delta vs PERMANOVA-on-HP-label：同假設不同統計量 | **非自由冗餘**（89.5% 一致、267/2543 不一致）；🔴 **distance-delta 不可移除**（SubcloneAnalyzer:80/109 + quality_score 硬依賴）；可移除的是 **label_hp_permanova**（0 個 C++ 消費者）|

**🔑 精簡冗餘修正（重要）**：先前假設要移除 distance 版（1b/2b）是**錯的** — `hp_merged_delta` 是 subclone 指派 + quality_score 的硬輸入，移除會破壞 pipeline。真正低風險可移除的是 **`label_hp_permanova`**（無 C++ 消費者，僅 CSV+Python script 用）。

---

## 2B. 邊緣案例驗證

**全 graceful（6 agent + 經驗確認，無 crash/garbage）**：tiny-N（gate 擋/早返）、single read（單葉樹）、all-same（silhouette 退化→fallback k=2）、NaN-heavy（SKIP peel-off 到完整子矩陣，2026-06-14 segfault fix）、single-HP tumor（812/2624，fallback normal-delta）、empty（早返預設）。

**🔴 本 session 新發現 — subclone Δβ 的 tiny-group 統計脆弱性**（最高槓桿精修點）：
- 獨立從 methylation.csv 重算所有 Δβ → 與 module **完全一致 ✓**（第三次交叉驗證）
- 但最強的兩個 subclone「sig」案例（chr1:185882341 subHP1=−0.864 / chr1:172309527 subHP1=−0.808）都是 **germline 僅 1 read** vs 數十 carrier → 單 read 主導
- 全 chr1 量化：subclone sig 中 **63-65% 穩健**（min-group≥5），但 **HP1 26%/HP2 23% min-group<3、15%/13% 由 1 read 主導**
- 與 2A 的 M2（cluster-PERMANOVA 缺 min_reads_per_group gate）**同類**：兩處都缺最小群守衛
- **→ 精修**：subclone Δβ（`compute_group_dbeta_test`）加 **min-group 守衛**（建議 ≥3，sig 才算）+ 輸出 fine 群大小（germline/carrier count）透明化

---

## 2C. 個案圖肉眼驗證（14 case / 6+1 層，`figures/case_*.png`）

> 每圖：左=甲基 β read×CpG / 右=read×read 距離（BERNOULLI）；側欄 HP色+T/N(橙=tumor/綠=normal)；標題含全方法數值。
> 選取 `select_inspection_cases.py`（4-層 task-gate）+ 穩健 subclone 補選。值全 grep 自本 run CSV。

| 層 | 案例 | 肉眼看點 |
|----|------|------|
| canonical germline ASM | 209109399（germΔβ=0.386 normal HP1/HP2=28/23）、103580838 | normal HP1 vs HP2 甲基均移；clustering/CramersV 是否同意 |
| **subclone 脆弱** | **185882341（subHP1=−0.864 但 germline=1）**、172309527（CramersV=0）| 「強」由單 read 主導 = artifact（對照穩健組）|
| **subclone 穩健** | **95863662（germ31/car62, CramV=0.74）**、**58061351（germ32/car13, sig 但 CramV=0）**| 真 subclone 視覺明顯；58061351 = **label-first 抓到 structure-first clustering 漏的**（mean-shift 非雙峰）→ A/B 兩方向互補實證 |
| somatic residual | 78328503（som=0.953 但 tumor HP2=2）、60459214 | 殘差「sig」但 tumor 一側 reads 極少 = validity 邊緣 |
| edge single-HP | 40805757（tumor 全 HP1=142）、40867755 | somatic haplotag 主導；normal 仍雙色 |
| edge NaN-heavy | 16637933（408→309 valid, VC=Noise）、16725885 | SKIP 剔除後分群是否仍合理 |
| edge tiny-N | 72933823、75106326（N=60 valid 31/40）| 小樣本分群退化行為 |

**outlier_cramersv1_lowN = 0 案例**（SKIP run 無「CramersV=1.0 + 低 N overfit」）= 好訊號（SKIP 去污染後無此類 overfit 假象）。

---

## 3. 收斂：可精修點排序（皆 rigor/hygiene，非翻轉錯誤）

| 優先 | 精修 | 類型 | 槓桿 |
|:---:|------|------|------|
| **1** 🔴 | subclone Δβ + cluster-PERMANOVA 加 **min-group 守衛(≥3)** + 輸出 fine 群大小 | cpp-change | 直接影響 subclone 可信判定（你的核心）|
| 2 🟡 | gap#4 Fisher 標 anti-conservative；ASM claim 勿單靠 Fisher 升 tier-4 | 文檔+紀律 | 防過度宣稱 |
| 3 🟡 | 修 stale 危險註解（normal hp="0"/REF）+ 加未知 HP 值 counter | cpp-change(註解+log) | 防未來誤修弄壞 germline baseline |
| 4 🟡 | gap#5 perm 數統一 999（或文檔明示 cluster 0.01 floor 不可 min-p 混比）| cpp-change(1 行) | 跨方法可比性 |
| 5 🟢 | gap#2 dispersion：struct 預設改 false（消 foot-gun）或加 placeholder 註解 | cpp-change | 移除 latent 風險 |
| 6 🟢 | 精簡冗餘：移除 **label_hp_permanova**（非 distance-delta）+ 更新 Python script | cpp+py | 真正可省的冗餘 |

> **無 NEEDS_FIX、無結論翻轉**。核心方法可信；精修都是「提升嚴謹度 / 防未來 footgun / 透明化」。subclone min-group 守衛（#1）是唯一直接影響研究判定的，建議優先。
