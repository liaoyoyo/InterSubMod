---
title: "顯著判定過嚴 + within-HP 是否算顯著 + 輸入合規性 — 設計分析"
date: 2026-06-16
type: methodology design + 根因確認（significance verdict 層 + input-compliance 護欄）
trigger: 用戶問 ① within-HP clean multigroup 該不該算「顯著特徵」② 救回的 PERMANOVA 顯著是否真結構/該否算顯著判定 ③ SEQC2 不可當輸入(只驗證)，LOH 用 HP-tag 推估、CN 用 coverage-peak 推估
data_sources: output/_wg_stage2_validate/significance_summary.csv (HCC1395 全基因組 30490 TP SNV, SKIP, commit b6dde2b), SignificanceAnalyzer.cpp:270-343, RegionProcessor.cpp:1385
claim_levels: 架構=L1(親讀 code) / 量化=L1(第一手讀回 WG CSV) / 修法=L2(設計)
status: 設計分析 — 待用戶決策（未動 C++）
related: 20260616_strength_score_reclassify_cpp_change_spec_01 / 20260616_structure_label_verdict_false_negative_audit_01 / memory project_ism_verdict_false_negative_audit_2026_06_16
---

# 顯著判定過嚴 + within-HP 顯著性 + 輸入合規

## 0. 用戶三問（本檔回答）
1. within-HP clean multigroup 該不該被算成「有顯著的特徵」？
2. 救回的 PERMANOVA 顯著是不是真結構 → 該不該也算「顯著判定」？是否要修「判斷顯著太嚴格」？
3. 合規護欄：**SEQC2 只能驗證後觀察、不可當輸入**；LOH 用 HP-tag 推估、CN 用 read-coverage peak 推估。

## 1. 現行兩個獨立判定（親讀 code，先釐清）

ISM 其實有**兩個不同的「判定」**，常被混為一談：

| 判定 | 來源 | 公式 | 性質 |
|------|------|------|------|
| **`Significant`（布林欄）** | `RegionProcessor.cpp:1385` | `passed_gating && global_p<=0.05 && CramersV>=0.1 && num_reads>=20` | 離散 cluster×label Fisher 的**嚴格 gate** |
| **`VerificationClass`（legacy）** | `SignificanceAnalyzer.cpp:326-340` | 2×2 of `cluster_significant`(=passed_gate&&FFH) × `label_sig` | Strong/Subclone/Weak/Noise |

🔴 **關鍵事實（FN1）**：`label_hp_permanova` / `label_allele_permanova`（`SignificanceAnalyzer.cpp:301-304`）**有算**（含 permutation p + dispersion），但**完全不進** `Significant` 也不進 legacy verdict。Δβ-perm 同樣被排除（Stage④ 已把它接進 VC，但仍未進 `Significant`）。

## 2. 量化坐實（全基因組 30490 TP SNV，第一手讀回）

| 指標 | 值 | 意涵 |
|------|----|------|
| legacy Noise+Weak 池 | 8180 | 「被判沒結構」的池 |
| 其中 `Significant`(嚴格 gate)=true | **65 (0.8%)** | gate 幾乎全擋 |
| 其中**有 valid 顯著 label-PERMANOVA** (p<0.05) | **7476 (91.4%)** | **真 permutation 檢定顯著卻被丟** |
| — HP 軸 / allele 軸 | 6484 / 6658 | allele 略多（貼 somatic/cis，符合你的觀察 3） |
| — dispersion **未** warned（乾淨 location-shift） | **7476（全部）** | 🔴 見 §4 |
| — dispersion warned（離散度混淆風險） | **0** | 🔴 PERMDISP 疑未啟用 |
| within-HP clean 但**無**顯著 PERMANOVA（純啟發幾何） | 124 | 見 §3 |

## 3. 回答①：within-HP clean multigroup 算不算「顯著」？— **分兩層，不可混**

within-HP clean multigroup 是 **heuristic 幾何判準**（silhouette>=0.5 + 平衡，或 mean-β gap>0.15 + var-reduction>0.5），**不是統計顯著性檢定**（無 permutation null、無 p-value）。

- **作為「結構存在 / 該保留」（VerificationClass 層）**：✅ **合適，已做**（Stage② 接進 MultiGroupNoLabel）。它回答的是「這 HP 內部有沒有甲基亞結構值得保留」。
- **作為「統計顯著」（`Significant` 布林）**：❌ **不該直接接**。幾何乾淨 ≠ 統計顯著（會犯 L2/over-claim）。要讓它「顯著」需補 **permutation null**（HP-family 內打亂 read，重算 silhouette/gap → 看觀測值是否超過 null 分佈 → p）。
- 實證支撐：124 個 region 是「within-HP 幾何乾淨但無顯著 PERMANOVA」→ 正是「幾何有、統計未證」的那群，**不應自動稱顯著**。其餘 within-HP clean 多與 PERMANOVA 顯著重疊（彼此佐證）。

> **一句話**：within-HP = 結構**存在**的證據（VC 該收），但**顯著性**要靠檢定（PERMANOVA / 或補 within-HP permutation null）才能宣稱。

## 4. 回答②：救回的 PERMANOVA 顯著是真結構嗎？該算顯著判定嗎？

- **是真的「顯著」嗎**：label-PERMANOVA 是 **permutation 檢定**，p<0.05 是正規統計顯著 → **這 7476 個確實是統計顯著、卻被 verdict 丟掉** = 「判斷顯著太嚴格」的**設計缺陷確認**。把 valid 顯著 PERMANOVA 認列為一種顯著，**方法學上正當**。
- 🔴 **但有一個必須先解的 blocker — dispersion 混淆**：PERMANOVA 在「群間離散度不同」時也會顯著（**非位置差異**，Anderson 2006）。要分辨「真 location-shift 結構」vs「離散度異質假象」需 **PERMDISP（dispersion 檢定）**。本資料 **DispersionWarn 全 0 / 7476 全 clean** → 高度懷疑 **dispersion 檢定未真正啟用/未生效**（不可能 7476 個全無離散度差異）。
  - **若直接把 PERMANOVA 顯著接進 `Significant`**：等於把「過嚴」換成「過鬆」——可能把離散度假象也標成顯著。
  - **必修前置**：先讓 PERMDISP 真正運作（查 `config_.enable_permanova` / dispersion 計算路徑為何全 false），用它**過濾** dispersion-driven 假顯著，剩下的 location-shift 才認列。
- **真結構是不是顯著判定**：可以，但 verdict 要**多軸並陳**而非單一布林——建議 `Significant` 升級成「顯著來源旗標集」：`{cluster_FFH, label_PERMANOVA(HP/allele, PERMDISP-passed), Δβ_perm}` 任一 valid 顯著 → significant，並標**哪個軸**顯著（HP=germline-ASM 嫌疑、allele=somatic/cis 嫌疑、cluster=共識）。

## 5. 回答③：輸入合規護欄（重要 — 重塑 Stage③）

**原則**：分類 / 顯著判定的**輸入只能用資料本身可導出的量**（BAM 衍生：甲基矩陣、HP tag、read coverage）。**SEQC2 CN/LOH BED 是外部真值 → 只能「判定後」做驗證/觀察，不可當分類輸入**。理由：(a) circular/leakage（用答案分類）；(b) 不可泛化（SEQC2 僅 HCC1395 有，方法跑不了其他樣本）。

| 量 | 合規自含估法（輸入用） | SEQC2 角色 |
|----|----------------------|-----------|
| **LOH** | **HP-tag 比例失衡**：`potential_loh = hp_ratio<0.1 \|\| >0.9`（`SignificanceAnalyzer.cpp:322`，**已存在**）| 驗證 ISM 自估 LOH 的 recall（post-hoc） |
| **CN** | **read coverage / depth peak** 分析（diploid 基線 expected_coverage 偏離）| 驗證自估 CN 落點（post-hoc） |
| 結構/顯著 | 甲基 read×read 距離 + HP/allele 標籤 + PERMANOVA/Δβ | — |

🔴 **Stage③ 撤回重設計**：原 spec §4「載 SEQC2 BED → `SEQC2_CN_Type` 欄 → LOH-aware 分類」**違反合規護欄，作廢**。改為：
- **LOH-aware** = 用既有 `potential_loh`(HP-ratio)（已合規，無需新輸入）。
- **CN-aware（CN-Substructure 細類）** = 新增**coverage-peak 自估 CN**（self-contained），非 SEQC2。
- **SEQC2** = 獨立 **validation 腳本**（post-hoc concordance：ISM 自估 LOH/CN vs SEQC2 真值），**不進 C++ 分類輸入**。

## 6. 決策點（待用戶定，未動 C++）
- **D1（顯著架構）**：是否把 `Significant` 從「單一嚴格 gate」改成「多軸顯著旗標集」（cluster_FFH ∪ label-PERMANOVA ∪ Δβ-perm，**經 PERMDISP 過濾**），並標顯著軸？
- **D2（前置 blocker）**：先修 **PERMDISP dispersion 檢定**（查為何 DispersionWarn 全 0）——這是把 7476 認列為顯著前的**必要關卡**（否則過嚴變過鬆）。建議 **D2 先於 D1**。
- **D3（within-HP 顯著性）**：within-HP clean 維持只進 VC（結構保留）；要不要再補 **within-HP permutation null** 才讓它能進 `Significant`？（124 個純幾何待證）
- **D4（Stage③ 重設計）**：確認 Stage③ 改為 coverage-peak 自估 CN + HP-ratio LOH + SEQC2 降為驗證腳本？

## 7. 用戶決策（2026-06-16）+ 實作 log

**決策**：D1 顯著架構=**先修 PERMDISP 再認列**（D2→D1）；within-HP=**補 permutation null 並輸出記錄，是否進 Significant 後憑紀錄再分析**；Stage③=**確認重設計**（SEQC2 降驗證、自估 LOH/CN）。

**D2 PERMDISP 修復（2026-06-16, 完成）— 改動合理正確確認**：
- **根因（親讀 code 坐實）**：PERMDISP 沒生效有兩因 — ① `RegionProcessor.cpp:2348` 明寫 `enable_dispersion=false`（pipeline 關掉，無註解理由）；② 即使開，`StructureTest.cpp:292` 的 dispersion p 是「F 統計量 3 段查表」(F>4→0.01 / >2.5→0.05 / else 0.1)，非真 F 分佈 p → warning⟺F>4 太粗。`check_dispersion` 核心（距質心距離 ANOVA = 標準 PERMDISP/Anderson 2006）方向正確。
- **改動**：① `RegionProcessor.cpp:2348` `enable_dispersion=true`；② `StructureTest.cpp` 加 `betai/betacf/f_dist_sf`（Numerical Recipes 正則化不完全 beta → F 分佈上尾 p），check_dispersion 改用 analytic p（df1=k群-1, df2=N_obs-k群）取代 3 段查表。**純加 dispersion 欄，不進 verdict（D1 才接）**。
- **驗證（合理正確）**：① build OK、**ctest 221/221** ② **regression 雙守護 byte-identical PASS（未 --update-golden）** → 證 dispersion 對 9 golden 欄（含 Significant/VC）**純加無耦合** ③ **analytic F-p vs scipy.stats.f.sf 8 組 (F,df1,df2) max abs diff 1.1e-13** → 數學正確 ④ 全基因組 dispersion 分佈（30490 TP SNV, SKIP, 891s）：DispersionWarn 由**全 0 → 29695/30490(97.4%)** 有警告（dispersion p 範圍 **5.4e-43~0.9999**，非 bucket → analytic p live）；7476 顯著 PERMANOVA 中 **6507(87.0%) dispersion-WARNED（混淆，dispersion p 常 1e-10 級，非大-N 邊際）、僅 969(13.0%) clean location-shift**（HP 6484/allele 6658 重疊）。

> 🔴 **D2 關鍵結論**：真正可認列的「顯著 location 結構」是 **969，非 7476**。若略過 D2 直接認列 7476，**6507(87%) 離散度異質假象會被誤標顯著** = 過嚴換過鬆。D2→D1 排序坐實必要。969 為 D1 認列保守候選；analytic-F 為 parametric（大-N 敏感），但本資料 dispersion p 極端(5e-43)故結論 robust（permutation PERMDISP 可作未來 refine）。
- **git**：D2 自成一 commit（StructureTest.cpp + RegionProcessor.cpp；regression-safe 無需更新 golden）。
