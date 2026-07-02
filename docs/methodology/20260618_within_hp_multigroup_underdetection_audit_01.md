<!--
建立時間: 2026-06-18
問題類型: C++ 方法 | 統計方法（within-HP 多群偵測 gate）
影響 track: paired（_wg_d1_unified）+ TO
狀態: pending_decision
build_branch: research/subclonal-reconstruction-202606 (worktree ism-review-infra)
data_sources: _assets/20260616_readset_provenance/within_hp_rederive.json, verify_categories.json
-->

# within-HP 多群偵測（一標籤多群 / subclone 方向）漏判 — 審查文件

## 問題描述

**「一標籤多群」**（單一 germline HP 內甲基再分 ≥2 群 ＝ within-HP substructure ＝ subclone 最相關方向）此現象**出現了，但被系統性漏判**：主聚類已把單一 HP 的 reads 拆成多群（read-level `max_label_span≥2`），但 `WithinHP_CleanMultigroup=False`（沒被採計為多群），位點被歸成「對齊-單群」（乾淨 1:1）。

**程式碼位置**：
- `compute_within_hp_substructure`（`src/core/RegionProcessor.cpp:2199`）：`find_optimal_clusters` 已找到 k≥2 分群，但**只在 `silhouette≥0.5 且 min_sz≥max(3, n/5)` 才採計** `out_ngroups=k`；否則回報單群。
- `within_hp_level_clean`（`RegionProcessor.cpp:2301-2342`）：LEVEL 軸，`best_gap>0.15 且 varexpl>0.5` 才算。
- `WithinHP_CleanMultigroup = (NG1≥2 或 NG2≥2 或 level_bimodal)`（:1189）。

## 量化影響（來源：within_hp_rederive_v2.json，從 methylation.csv 忠實重建 BERNOULLI + 同 gate；N 擴大版）

**重建忠實度（控制組 B：PATTERN, CleanMG=True, NG≥2, N=50）**：50/50 重建 silhouette≥0.5（49/50 過完整 gate，median 0.54）→ **重建方法忠實，量化可信**。

**漏判量化（Group A：113 個 對齊-單群, CleanMG=False, seed 11；擴大自 pilot N=38）**：

| 重建後狀況 | 數量 | 比例 | 被什麼擋 |
|---|--:|--:|---|
| 強但不平衡（sil≥0.5, minor 子群 min_frac 2-4%）| **30** | **26.5%** | **balance 要求**（min_sz≥n/5＝20%）→ 擋掉 minor/低-CCF subclone |
| 弱分群（sil 0.3–0.5）| **64** | **56.6%** | **silhouette≥0.5 閾值** → 擋掉弱-中分群 |
| 真噪音（sil<0.3）| 19 | 16.8% | 正確不採計 |

- **94/113（83.2%）有 within-HP 分群（sil≥0.3）但沒被 flag**；median best_sil = 0.413（剛好在 0.5 閾值下）。
- `rederive_gate_pass=0` → **無自我矛盾**（C++ gate 一致），是**刻意保守的閾值**造成漏判，非 bug。
- **minor-subclone（被 balance 擋）在大樣本下更高（26.5% vs pilot 15.8%）** → 強化方案 B。min_frac 2–4% 的小子群 ＝ 正是低-CCF subclone。
- ⚠ **單樣本（_wg_d1_unified, HCC1395）+ 同 binary 不同 config 冗餘**；真正跨**生物樣本**確認 deferred 到 /cpp-change 驗證（改後 binary 跑 COLO829/H1437 標準 3 樣本回歸）。這次 N=113 數字＝改前 baseline。

## 根因（method/threshold，非 code bug）

1. **silhouette≥0.5 是「粗略 null 替代品」**：silhouette 在隨機 reads 上也會過度偵測（random 仍分出緊湊群），故用高閾值 0.5 壓 FP。代價＝**弱但真的 subclone（sil 0.3–0.5, 58%）被一併擋掉**。
2. **balance（min_sz≥max(3, n/5)＝20%）反 subclone**：subclone 常是 minor（低細胞比例 → 少 carrier reads → 不平衡）。此要求**系統性排除最有意義的低-CCF subclone**（本例 16% 是這類）。
3. **二元 CleanMG 一刀切**：把連續的「分群強度」壓成 True/False，邊界（sil≈0.5）位點全歸 False，丟失 gradient。
4. **方法兩難（為何不能直接降閾值）**：直接降 silhouette 會放進噪音；而「用 PERMANOVA 對衍生 cluster 做顯著性」＝ **double-dip 循環**（cluster 由距離生、又用同距離檢定 → 恆顯著；見 memory `project_tumor_only_axis_negative_subclone_classification`）。**正解需 data-permutation null（gap statistic / clusterboot Jaccard）**，非 label-permutation（見 `project_distance_matrix_cluster_validation_methods`）。

> 對齊本身（cluster×label 列聯表）方法已驗證可重現（mapping 60/60＝100%，見 verify_categories.json）；本問題**只在 within-HP 多群這個子軸**，且方向是**假陰性（漏判 subclone）**，與既有 FN 審查一致（memory `project_ism_verdict_false_negative_audit_2026_06_16`）。

## 修改選項

### 方案 A：不修改（接受現狀 + 補文件）
保留 sil≥0.5 + balance；在報告標註「within-HP 多群為高特異性、低敏感度偵測，漏判 ~74% 弱/minor subclone」。

### 方案 B：放寬 balance（救 minor subclone）⭐ 低風險高價值
`min_sz≥max(3, n/5)` → `min_sz≥3`（絕對下限，移除 20% 相對要求）。只影響「sil≥0.5 已乾淨、僅不平衡」的 minor subclone（**本例 30/113＝26.5%**）。silhouette 仍守乾淨度。
- 位置：`RegionProcessor.cpp:2199`（單行條件）。

### 方案 C：降 silhouette 閾值 + data-permutation null（救弱分群）中-高風險
sil≥0.5 → sil≥0.35，但**加 gap-statistic null**（permute 每 CpG 的甲基值、重聚類、比 sil 對 null 分布）守 FP。救弱-中分群（**64/113＝56.6%**）但需新增 null 計算（~999 permutation/HP）。
- 位置：`RegionProcessor.cpp:2142-2199` 重構 + 新 `within_hp_gap_null()`。

### 方案 D：分級輸出（不動邊界，暴露 gradient）⭐ 與 B/C 正交可疊加
新增欄位 `WithinHP1_BestSil / WithinHP1_MinFrac / WithinHP_GapP`（連續值），不只 binary CleanMG。下游（VC / report）自選閾值。
- 位置：`RegionProcessor.cpp` emit 段（+3 欄）+ 表頭。

## SWOT

**方案 B（放寬 balance）**
|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** 單行改動、silhouette 仍守乾淨度、直接救低-CCF subclone（論文主軸最在意）| **W** minor 子群統計力弱（n 小），Δβ/下游檢定可能不顯著 |
| External | **O** 對齊 subclone reconstruction 目標、低-CCF 正是難點 | **T** 若無下游 null 守，minor 群可能是 mapping/strand artifact |

**方案 C（降閾值 + gap null）**
|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** 統計嚴謹（data-permutation 非 double-dip）、救最大宗（58%）| **W** 計算成本（999 perm/HP × 2 HP × 30k 位點）、需新 module + 測試 |
| External | **O** 對齊業界 chance-partition（Şenbabaoğlu2014）、可寫方法貢獻 | **T** gap-stat 在小 n 不穩、閾值仍需校準恐引入新 FP |

**方案 D（分級輸出）**
|  | Helpful | Harmful |
|---|---|---|
| Internal | **S** 不動既有邊界（零回退風險）、暴露 gradient 供研究 | **W** 不「自動」救漏判，需下游再決策 |
| External | **O** 與 reclassify-v2 / verify-workstation 判讀整合、可人工複核 | **T** 多欄位增 CSV 寬度 + 下游需更新 |

> W+T vs S+O：B = 2 vs 3（推進）；C = 2 vs 3（推進但需 Step2 補 null 校準）；D = 1 vs 3（推進）。

## 驗收標準
- [ ] test-quick 通過 + 編譯通過（C++ Hard Gate）
- [ ] 方案 B：救回的 minor-subclone 子群在 within_hp_rederive 重算可重現（sil≥0.5）
- [ ] 方案 C：gap-null 在控制組（隨機 reads）FP < 5%；PATTERN 控制組仍 100% 偵測
- [ ] 跨 ≥3 樣本確認 CleanMG 數量變化方向一致、無 noise 暴衝
- [ ] 不回退既有對齊（mapping validation 仍 100%）

## 用戶決策
**選擇**：[x] **B + C + D**（綜合修正：放寬 balance + 降 sil 閾值＋gap-null + 分級輸出）
**日期**：2026-06-18
**理由**：救回低-CCF subclone（B）+ 弱分群最大宗（C）+ 暴露 gradient 供研究/人工複核（D）。
**前置條件**：**動碼前先擴大量化**（每組 ~100 位點、跨樣本確認漏判率/控制組重現方向一致），數字穩固後才進 /cpp-change。
**狀態**：✅ 已實作並驗證（見下節）。

## 實作 + 驗證（2026-06-18 — cpp-change 完成，branch feat/summary-nreadsvalid）

### 改動（清楚敘述）— 3 commit
- **c6b5e6c**（接手）：先前未提交的 tumor-only PERMANOVA 軸 + 5-component StrengthScore（用戶要保留作分析標記）。
- **081a556 B+D**：
  - **B** `RegionProcessor.cpp:2199`：`min_sz >= max(3, n/5)` → `min_sz >= 3`（移除 20% balance，救 minor subclone；silhouette≥0.5 仍守乾淨度）。
  - **D**：`compute_within_hp_substructure` +`out_min_frac`；struct +`within_hp_best_sil/min_frac`；emit +`WithinHP_BestSil/MinFrac`（分級暴露 gradient）。
- **49a806a C′（取代原 C gap-null）**：🔴 **原方案 C（降 sil + gap-null）放棄** — memory 證 gap-null 對 tumor-only 軸 NEGATIVE（read-內相關）。改 **C′ = `compute_within_hp_subclone_permanova`**：within-HP 子距離矩陣 + germline-tag(0)/carrier-tag(1) **a-priori 標籤** → `StructureTest::run_permanova`+`check_dispersion`（重用既有 module，無 double-dip）；emit +5 欄。對齊用戶定向「結構偵測(silhouette) + a-priori 標記確認(PERMANOVA)」。

### 驗證（更正確達成目標，非只數字變）
- **B 救援真實**（chr1 HCC1395）：PATTERN 24→51，救回 27 minor subclone，**全部 carrier 支持**（HP1S+HP2S≥3=有 somatic reads），BestSil 全≥0.51。5/27 邊際(min_sz 3-4)但 carrier 支持 + D 暴露低 min_frac。
- **C′ 判別（解 double-dip）**：C′ sig 在 **Noise VC 僅 12.1%** vs 有結構 28.1%；對比 `tumor_intrinsic`（無監督）Noise **97%** → **a-priori 標籤確實解掉 double-dip**。C′ 覆蓋 61%、與 B 互補（交集僅 99/711/442）。
- **不回退對齊**：regression 雙守護 PASS（SKIP 主基準 byte-identical 2624 位點；對齊欄 GlobalP/CramersV/PassedGating 0 變動）；ctest 221/221。MAX_DIST golden 因 B 故意改 within-HP 偵測而 --update-golden（非回退）。
- ✅ **strength_tumor 修正完成（44f434b + 註解修正 01bc9d3）**：見下節「strength_tumor 修正 + 對抗驗證」。
- ✅ **全基因組 HCC1395 scale 驗證完成**：30,490 TP SNV，對齊欄（GlobalP/CramersV/PassedGating）0 變動；跨生物樣本(COLO829/H1437) deferred（原始 BAM/VCF 資料不可及）。

### strength_tumor 修正 + 對抗驗證（2026-06-18 — commit 44f434b / 註解 01bc9d3）

**問題（task#7 發現）**：committed `tumor_intrinsic`（無監督 tumor-only cluster PERMANOVA，含 PERMDISP gate）重驗仍 double-dip（chr1 Noise VC 96.3%、WG 91.5%）；而 5-元件 `StrengthScore` 的 `strength_tumor` 元件原本 `= tumor_intrinsic ? clamp01(tumor_only_silhouette/0.5) : 0` → **把 double-dip 灌進評分**。

**改動（清楚敘述）** `RegionProcessor.cpp:1257-1258`：
- `strength_tumor = within_hp_subclone_sig ? clamp01(within_hp_best_sil/0.5) : 0`（gate 改 C′ a-priori，magnitude 用 within-HP best silhouette）。
- `tumor_only` emit 欄保留為 raw 診斷（observability），只剔出 StrengthScore verdict。
- **不進 VerificationClass**：reclass（:1330-1356）只 branch on Δβ/HP-AUC/cluster/LOH/within_hp_clean_multigroup/PERMANOVA，**不讀 strength**（C5 對抗確認單向資料流）。

**驗證（更正確達成目標，非只數字變；WG 30,490 重derive + chr1 2,624 fresh ground-truth）**：
- **double-dip 移除（判別方向翻轉）**：`strength_tumor>0` 比例 — OLD(`tumor_intrinsic`) noise **91.5% > structure 70.2%**（反判別）→ NEW(C′) noise **7.2% < structure 17.2%**（all-VC；valid-coverage denominator 12.7% vs 28.4%）。noise mean-β 0.879→0.064。chr1 fresh gate-consistent 100%（struct+noise 2459/2459）；392 個 noise VC StrengthGrade 去污染降級。
- **口徑（守 memory 非循環紅線）**：成果＝「double-dip removed / **~2.4× enriched on legacy VC**」，**非**「more discriminative / 0% noise」。C′ 是 **低 recall gate**（只 ~17% 真結構位點開火、~40-44% 位點因 germline/carrier 各 <3 reads 而 abstain），判別是 **noise-vs-structure RATIO** 非高 sensitivity。
- **個案親驗（meth 層）**：去污染 48 例（Noise+TI=T+C′=F）**0/48 帶 a-priori Δβ 標記**（⚠ 此乾淨僅 scope 到 VC=Noise；全 VC 去污染 1630 例中 ~502 非-Noise 帶標記，故「去污染不丟訊號」必標 VC=Noise 限定）。最極端 `Noise_Uniform` 但 strength_tumor=1.0 案（chr1:220099749）：全域均值 germline 0.869≈carrier 0.891（故 Noise_Uniform）但**逐 CpG 9 個 |Δβ|>0.3、max 0.925** → C′ 抓全域均值漏掉的「a-priori 條件軸 per-CpG 子結構」＝**feature 非 bug**。正控 1760（Strong, C′F=13.7）平衡 n（20/18）梯度合理。
- **對抗驗證（6-agent Workflow，全 CONFIRMED / 0 REFUTED / 0 blocking）**：C1 C′ 標籤源碼證 a-priori（`hp_tag=="1-1"`，非 distance-derived）→ 無 double-dip；C2 翻轉三路交叉復現（WG OLD 欄 / chr1 native NEW 欄 / chr1 公式重derive 2624/2624 zero-mismatch）；C3 去污染 0/48（scope=Noise）；C4 Noise+sig=feature（rare ~0.15%）；C5 strength 不回流 VerificationClass；C6 與 memory 框架一致。唯一具體缺陷＝stale 註解數字，已修（01bc9d3）。

### 與背景一致
C′ 用 a-priori 標籤 PERMANOVA（守 double-dip 紅線，對齊 D1/D2 既有 PERMANOVA 路線）；不動「甲基非重建驅動 / subclone ⭐3」框架；`strength_tumor`/C′ 為「**mark for confirmation**」非 hard subclone verdict（「subclone」＝ within-HP germline-tag vs somatic-carrier-tag 甲基分離，非 multi-subclone 重建）。對抗驗證 6/6 CONFIRMED 證實口徑與 [[project_tumor_only_axis_negative_subclone_classification]]（unsupervised NEGATIVE / a-priori 軸合法）一致。

---
*數據：`_assets/20260616_readset_provenance/{within_hp_rederive,verify_categories}.json`（seed 11/7，控制組 20/20 重現）。相關 memory：`project_ism_verdict_false_negative_audit_2026_06_16` / `project_tumor_only_axis_negative_subclone_classification` / `project_distance_matrix_cluster_validation_methods`。*
