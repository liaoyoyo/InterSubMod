---
title: ISM 方法認知校正 + 開放問題（Q1-Q4）+ 邊界紅線
date: 2026-06-22
status: in_progress
audience: 論文作者（自我認知整理）+ 跨 agent
data_sources:
  - src/core/DistanceMatrix.cpp
  - src/core/SignificanceAnalyzer.cpp
  - src/core/HierarchicalClustering.cpp
  - include/core/Config.hpp
  - include/core/SignificanceAnalyzer.hpp
  - include/core/DistanceMatrix.hpp
related_memory:
  - project_ism_method_soundness_validation
  - project_subcluster_cluster_count_determination
  - project_subclone_snv_difficulty_methylation_framework
  - project_zar1l_brca2_asm_verification
build_branch: research/subclonal-reconstruction-202606
---

# ISM 方法認知校正 + 開放問題 + 邊界紅線

> **論文**：Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing
> **定位**：把作者對 ISM 方法的完整認知整理校正（對照第一手原始碼），並列出論文 motivation 的開放計算問題（Q1-Q4）與敘述紅線。
> **方法層校正全部第一手 grep 自原始碼（file:line），研究數字以路徑交叉引用（不在此重述數字，§13 反捏造）。**

---

## 0. 結論先講：你的認知骨幹正確，3 處要校正（其中 1 處是你原本就對）

| # | 你的認知 | 第一手原始碼 | 判定 |
|---|---|---|---|
| **C1 距離** | read 間用 **Bernoulli** 距離 | canonical pipeline **確實是 BERNOULLI**（`Config.hpp:40` 預設 `{BERNOULLI}`；`run_vcf_all_snv.sh:65`/`run_batch_vcf_analysis.sh:24` `METRICS="BERNOULLI"`；regression「主研究基準=BERNOULLI」）| ✅ **你是對的**（細節見 1.2-C1） |
| **C2 檢定分工** | **Fisher 測「結構是否存在」** | ❌ Fisher 測「cluster×標籤**關聯**」；**PERMANOVA 才測「結構/分離是否存在」** | ⚠ **需校正**（見 1.2-C2） |
| **C3 切 k** | silhouette 有閾值（≥0.1 才算可分）、可能多 k / 非等高切 | ❌ **無硬閾值**、**永遠回單一 best_k**、非等高、bootstrap 預設關 | ⚠ **需校正**（見 1.2-C3） |

> 重要副產品：先前一個 Explore agent 曾誤報「距離預設是 NHD 不是 Bernoulli」——那是把 `DistanceMatrix.hpp:22` 的 **struct 欄位 fallback 預設（NHD）** 誤當成 **pipeline 實際設定**。實際 run 一律帶 `--distance-metric BERNOULLI`。**這就是「不憑二手回報、要回原始碼」的價值。**

---

## 1. 校正後的 ISM 方法流程（程式碼實證）

### 1.1 骨幹 6 步操作流程（逐步對照）

| 步 | 你的認知 | 原始碼實際 | 錨點 | 判定 |
|---|---|---|---|---|
| 1 | 依 SNV 定位、取經過 SNV 的 read | ✓ 如述（tumor 在 SNV point 檢查 alt-support；normal 取 window） | `RegionProcessor.cpp` ReadParser | ✅ |
| 2 | ±5000bp window 取甲基 + 壓縮 | window ✓（±5000）；但 **canonical Bernoulli 用「原始機率值」非 binary 壓縮**（binary 只用於 NHD/Jaccard，非 canonical）| `Config.hpp` window；`DistanceMatrix.cpp:341-345` | ⚠ 「壓縮」用詞偏 |
| 3 | read 間 Bernoulli 距離矩陣 | ✓ BERNOULLI（期望不一致 + 信賴加權；見 1.2-C1）；C_min=3 共同位點；NaN 策略**預設 SKIP**（非 MAX_DIST，2026-06-14 起，去 1.0 污染）| `DistanceMatrix.cpp:254-302`；`Config.hpp:39-40` | ✅ |
| 4 | UPGMA 聚類 | ✓ UPGMA 預設（亦支援 Ward/Single/Complete）| `HierarchicalClustering.*` | ✅ |
| 5 | best silhouette 切 k（無法分群 / k≥2）| ⚠ **無閾值**、永遠回單一 best_k（見 1.2-C3）；「無法分群」**不是在這層判**，是下游檢定判 | `HierarchicalClustering.cpp:488-549` | ⚠ |
| 6 | Fisher 顯著 + Cramér's V 強度 | ✓ 兩者都有、且是 **4 層**（ALT/REF、HP-pure、HP-family、HP-fine 4 群）；但角色 ≠「結構存在性」（見 1.2-C2）| `GlobalTest.*`、`SignificanceAnalyzer.cpp` | ⚠ 角色 |

### 1.2 三處校正詳述

**C1 — 距離是 Bernoulli（你對），但不是「binary 壓縮」**
canonical 路徑用**原始甲基機率**（`raw_matrix`），公式（`DistanceMatrix.cpp:254-302`）：
- 期望不一致 `delta(p,q) = p(1−q) + (1−p)q`
- 信賴權重 `weight(p) = 2·|p−0.5|`（p≈0 或 1 高信賴、p≈0.5 權重→0）
- 距離 `= Σ(w_i·w_j·delta) / Σ(w_i·w_j)`，範圍 [0,1]；共同位點 < C_min(=3) 或總權重≈0 → 回 −1（無效）

→ 論文寫法：ISM 距離 = **信賴加權的期望甲基狀態不一致率**（比 binary Hamming 多用了甲基「確定度」資訊）。**勿寫「binary 壓縮」**（那是 NHD/Jaccard 路徑）。NaN-pair 預設 **SKIP**（不灌 1.0），是現行主基準。

**C2 — Fisher 測「關聯」、PERMANOVA 測「結構/分離」（你把兩者角色對調了）**
原始碼 5 phase（`SignificanceAnalyzer.cpp:79-340`）精確分工：

| 檢定 | 程式變數 | 測什麼 | 是否進主 verdict |
|---|---|---|---|
| **Fisher / Cramér's V** | `global_alt` / `global_hp_family` | 無監督 cluster **是否對齊**生物標籤（cluster-first 關聯）| ✅ 是 gate + `cluster_significant` |
| **cluster-PERMANOVA** | `result.permanova`（`cpp:123`，跑在 **cluster_labels** 上）| 無監督群在距離空間**是否真的分離**（結構存在性）| ⚠ 餵信心/ambiguity，**非**主 verdict 條件 |
| **label-PERMANOVA** | `label_hp_permanova` / `label_allele_permanova`（`cpp:301-304`）| **標籤定義的組**在距離空間是否分離（label-first）| ✅ 是 `label_sig` |
| **PERMDISP** | `check_dispersion`（`cpp:126,264`）| 群內離散度是否相同（防 location-vs-dispersion 混淆）| ⚠ warning 旗，不自動否決 |

→ 所以你問「PERMANOVA 有沒有用」：**它就是結構存在性的核心檢定**（兩個位置都用 PERMANOVA：一個測無監督群、一個測標籤組）。Fisher 不是測結構存在，是測「cluster 與標籤的關聯」。

**C3 — silhouette 無閾值、永遠回單一 best_k**
`find_optimal_clusters`（`HierarchicalClustering.cpp:488-549`）：`best_score` 從 −2.0 起，迴圈 k=min_k..max_k 取平均 silhouette **最高**者，回 `{best_k, best_labels}` 單一解。
- **沒有「silhouette≥0.1 才算可分」的閾值** → 你看到 0.1~0.8 全有，是因為**系統一定挑同批最高的 k，即使絕對值很低**。
- **永遠會切出 best_k≥min_k**（通常 ≥2）→ 「這個 region 其實沒有結構」**不是 silhouette 判的**，是下游 cluster-PERMANOVA（結構存在）+ Fisher（標籤關聯）判的。
- 回單一 k（+ 離群 retry k+1），**非多 k、非每子群不同高度切**。
- bootstrap 穩定性**已實作但預設關**（`SignificanceAnalyzer.hpp:37 enable_bootstrap=false`；`RegionProcessor.cpp:1717` 亦設 false）。

> **這直接決定了「切幾群」的方法論方向（已拍板）**：因為無監督 silhouette **一定會切**、不能當「有沒有 / 幾群」的真值 → **改以「切法是否對齊獨立 a-priori 軸（CramérV≥0.7 + Fisher）」來確認 ≥2 群**，不追無監督「抽象幾群」（memory `project_subcluster_cluster_count_determination`：無監督抽象群數無乾淨解，模擬證 gap 會把 1-clone 判成 k=3-4）。

---

## 2. 最終 4-state verdict 邏輯（第一手，`SignificanceAnalyzer.cpp:307-340`）

```
cluster_significant = passed_gate AND (global_alt.fisher_p ≤ 0.05 OR global_hp_family.fisher_p ≤ 0.05)   # Fisher 關聯
label_sig          = label_significant                                                                    # label-PERMANOVA

 label_sig &  cluster_significant → "Strong"    （雙路一致；HP 極端比例另以 LOH_Subtype track）
!label_sig &  cluster_significant → "Subclone"  （有距離結構、但與已知標籤無關 = 可能新 subclone 軸）
 label_sig & !cluster_significant → "Weak"      （標籤有訊號、但 cluster 弱）
其餘                              → "Noise"
```

verdict 字串確切是 `Strong` / `Subclone` / `Weak` / `Noise`（非其他命名）。
**關鍵洞見**：主 verdict 由 **Fisher（cluster×標籤）+ label-PERMANOVA** 決定；cluster-PERMANOVA（無監督結構）只調信心。這也是「PERMANOVA 大-N 偏 liberal（valid 率偏高）」caveat 的位置 → 待加 effect-size / TP-vs-FP 門檻（計劃 WS-A3 `T-METHOD-EFFECTSIZE`）。

---

## 3. 對「方法穩定性自評 / 待驗證」的逐項回應

| 你的疑慮 | 回應 | 對應任務 |
|---|---|---|
| silhouette 切線是否該含 0.1 以上 | 程式**本來就不設閾值**、一定挑最高 → 問題不在「閾值多少」，而在「無監督群數不可信」→ 改 a-priori 軸驗證 | WS-A1（`T-SL-C1` reframe）|
| k 該用單個還是多個、怎麼定義紀錄 | 程式回單一 best_k；論文**不主張無監督 k 真值**，主張「切法對齊 a-priori 軸時 ≥2 群可確認」 | WS-A1 |
| 結構是否穩定 | cluster-PERMANOVA（結構存在）+ PERMDISP（防離散度假象）已實作；bootstrap optional | C2；WS-A1 |
| 非等高切（每子群不同高度）| 現行 `cut_by_num_clusters` 是等 k 切、非等高；非等高切是另一種演算法（暫不追，記入 WS-D 邊界）| WS-D |
| Fisher 判別標準不穩 | cluster×標籤 Fisher 是精確檢定（OK）；**真正要修的是 per-CpG 甲基層 Fisher 的 over-dispersion → beta-binomial**（兩個 Fisher 不同層，勿混）| WS-A2（`T-V3`，走 `/cpp-change`）|
| PERMANOVA 是否有用 | **是核心結構檢定**（C2）；`/V-1` 判 SOUND_WITH_CAVEATS；唯一真 caveat = 大-N 過敏 → effect-size 門檻 | WS-A3（`T-METHOD-EFFECTSIZE`）|

---

## 4. Q1-Q4 論證狀態（論文 motivation + 甲基價值論證）

| Q | 內容 | 狀態 | 證據 / 待補 |
|---|---|---|---|
| **Q1** | ~30,490 sSNV 全基因組兩兩間距 × ONT 讀長 → 單讀可跨幾個 sSNV、理論覆蓋 | 🔴 **全新計算（無任何產物）** | 待寫 `scripts/analysis/ssnv_spacing_coverage.py`（WS-B1 `T-Q1-COVERAGE`）。**這是「為何只用 sSNV 難」的量化基石** |
| **Q2** | 為何直接用 read 的 sSNV 建 subclone 困難 | ✅ 已有完整框架可引 | `InterSubMod/docs/method_comparison/20260619_subclone_analysis_interpretation_full_framework_01.md` §5 + `20260621_clone_subclone_reconstruction_landscape_and_ism_feasibility_01.md` §1（multiplicity 歧義 / 低-AF floor / 非唯一性 / bulk 無跨-locus 連結）。投稿前核對 Tarabichi 原句 |
| **Q3** | 甲基對「找突變相關 cis」有意義 | ⚠ pilot 證據有、須慎框架 | BRCA2：`InterSubMod/docs/experiments/in_progress/2026/06/20260620_per_cpg_multiaxis_attribution_brca2_pilot_01.md`、memory `project_zar1l_brca2_asm_verification`。**紅線**：normal unphased 受限、乾淨 cis 稀有 → 禁寫「甲基驅動偵測 / 當 filter」 |
| **Q4** | 甲基能找潛在無法分析狀況、當獨立資訊輔助 | ⚠ 證據強、邊界嚴 | HP 軸非循環（haplotag 零甲基來源）+ latent cluster（有標籤可見、無自動分群）；memory `project_ism_method_soundness_validation`、`project_asm_locus_display_and_cramersv_reliability_gate`。**紅線見 §5** |

---

## 5. 邊界紅線（敘述誠實底線 — 投稿前必守）

1. 甲基角色 = **characterize / corroborate**（刻畫讀層結構、佐證 haplotag 一致性），**禁** driver / genome-wide tree reconstruction / subclone filter。
2. **reconstruction** 在本論文 = regional、LOH-constrained partition（受 phase-set/觀測單元天花板限制），**非**完整 phylogeny-CCF tree。
3. G-B（within-hap somatic null）未跑前 → **禁宣 somatic-specificity**（LOH-unmask confound 82-91% 由 CN 解釋，未拆）。
4. 對手定位 = 「無監督 read×read 距離結構 PERMANOVA + normal-baseline cis-test + somatic-subclone 目標」；**禁**用「對手二代定序 / 對手缺顯著性檢定」當差異（cvlr/ASMS/MethylBERT 都 ONT-capable 且有 randomization 檢定）。

---

## 6. 接下來（對應計劃 WS-A 起點）

method-verify-first（已拍板）→ 起手 = WS-A：
- A1 `T-SL-C1` reframe 成 a-priori 軸驗證（不追無監督群數）
- A2 `T-V3` per-CpG Fisher → beta-binomial（`/cpp-change`）
- A3 `T-METHOD-EFFECTSIZE` PERMANOVA 大-N guard
- A4 `T-SL-C3` 分類註解 emit（對齊軸 / HP 組合 / LOH / CNV / 一標籤多群 / 一群多標籤）
- 平行 WS-B1 `T-Q1-COVERAGE`（Q1 計算，論文 motivation 基石）

> 完整工作流與任務對應見計劃檔；任務節點狀態見 `state/tasks/tasks_board.html`。
