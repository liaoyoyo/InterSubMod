---
title: "強度分數 + reclassify-v2 判別 — cpp-change 完整 spec（定稿）"
date: 2026-06-16
type: cpp-change spec (定稿，待實作)
trigger: 原最終判別不足(HeuristicScore 飽和 / VC 粗+標籤錨定 / 漏 Δβ-HP-AUC)；pilot 驗證修正(46.6% 救回/0 FP)
data_sources: reclassify_v2_pilot.json, within_hp_clean_pilot.json, /tmp/dbeta_wg_abc(強度分數驗證)
claim_levels: 公式/邏輯=L2(設計) / pilot 數據=L1
status: SPEC 定稿 — 待用戶 GO 進 cpp-change
related: 20260616_cn_loh_aware_structure_classification_design_01 / 20260615_structure_false_negative_classification_plan_01
---

# 強度分數 + reclassify-v2 — cpp-change spec

## 0. 為何改（原方法不足，已確認）
| 原 | 問題 | 證據 |
|---|------|------|
| HeuristicScore | 飽和(median20/max22) | 無梯度無法排序 |
| VerificationClass 4類 | 粗+標籤錨定 | 46.6% 假陰性可救(0 FP) |
| 判別只認離散 cluster×標籤 | 漏 Δβ/HP-AUC | 假陰性圖 6202793 HP-AUC0.969 判Weak |

## 1. 強度分數（取代 HeuristicScore）— 已驗證

```
struct = clamp((HP_AUC_All - 0.5)/0.5, 0, 1)        # 距離結構
eff    = clamp(|GermlineAsmDbeta|/0.3, 0, 1)         # Δβ effect size
assoc  = CramersV                                    # 標籤關聯(reliable-gated)
sig    = clamp(-log10(max(GlobalP,1e-6))/4, 0, 1)    # 顯著性(GlobalP=0→1)
StrengthScore = 0.30*struct + 0.25*eff + 0.25*assoc + 0.20*sig   # ∈[0,1]
```
- **移除 ClusterPermanovaF**（過敏/循環，always-high 不分辨）
- **5 級 StrengthGrade**：A≥0.65 / B 0.5-0.65 / C 0.35-0.5 / D 0.2-0.35 / E<0.2
- **驗證**(全基因組)：median 0.344/max 0.975；Top500=498 Strong(concord)；Bottom500=358 Noise+142 Weak；抓 31 個 VC 漏的高分真強。

## 2. 新 VerificationClass（reclassify-v2 邏輯）— 已驗證 46.6% 救回/0 FP

```
# VALID = 任一結構證據(預設保留, 用戶哲學)
valid = Δβ_sig                              # any Δβ module sig (level-shift, 主因)
     OR HP_AUC_All>=0.7                       # 連續距離結構
     OR cluster_matches_label                 # 原 GlobalTest cluster_sig
     OR within_HP_clean_multigroup            # §3 新指標
     OR (SEQC2_LOH AND structure)             # §4

# 細類
Strong            : cluster_matches_label (+Δβ 最強)
LabelShift        : Δβ_sig 但無 distance 結構 (2131, 主救回)
CN-Substructure   : within-HP多群 AND SEQC2-gain (627)
MultiGroupNoLabel : within-HP乾淨多群 不對標籤 (323)
LOH-Structure     : SEQC2-loh + 結構 (609)
StructureNoLabel  : HP-AUC>=0.7 但無標籤對應 (120)

# NOISE = 無上述任一 + 標籤交錯, 三子型(用 Epipoly/PairwiseMeanDist; NME 飽和不用)
Noise_Uniform     : PairwiseMeanDist<0.15 或 Epipoly<0.2 (均勻)
Noise_Chaotic     : Epipoly>0.5 或 PairwiseMeanDist>0.35 (混亂高熵無群)
Noise_Uncorrelated: 其餘
```

## 3. 新指標：within-HP clean multigroup（C++ 待加）
- **pilot 用 1-D mean β**（loose 40.7%→clean 14.2% 平衡+gap）；🔴 **C++ 應改 pattern-based**：對單一 HP-family 的 read **per-CpG 距離 clustering**（複用 HierarchicalClustering），silhouette>=0.5 + 兩群各>=max(3,20%) → clean multigroup。比 1-D meanβ 嚴謹（抓 pattern 異質非僅 level）。
- 新欄：`WithinHP1_NGroups` / `WithinHP2_NGroups` / `WithinHP_CleanMultigroup`(bool)。

## 4. 新輸入：SEQC2 CN/LOH BED（C++ 待加）
- 載入 `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`（chr/start/end/{gain,loh,loss}）。
- 每 region overlap → `SEQC2_CN_Type` 欄。**LOH-aware 用此真值**（ISM Potential_LOH recall 僅 13%）。
- ⚠ 路徑為 CLI 參數（可選；無 BED 時 LOH-Structure 退回 ISM Potential_LOH）。

## 5. 新輸出欄
`StrengthScore, StrengthGrade, VerificationClass(v2 細類), WithinHP1_NGroups, WithinHP2_NGroups, SEQC2_CN_Type`
（保留舊 VerificationClass 為 `VerificationClass_Legacy` 以利對照/回溯）

## 6. cpp-change 實作步驟（Hard Gate，主 agent）
1. `compute_strength_score`(SignificanceResult) → StrengthScore + Grade（取代 compute_heuristic_score 或並存）
2. `compute_within_hp_clustering`（per HP-family 距離 clustering，複用 HierarchicalClustering + StructureTest）
3. `LohCnAnnotator` 載 SEQC2 BED（複用 LohBedAnnotator 模式）
4. 改 SignificanceAnalyzer.cpp:326-339 判別邏輯為 §2 v2（保留 legacy 欄）
5. header + emit 新欄；hpp RegionResult 欄位
6. build + ctest + **regression（VerificationClass/Significant 會變 → --update-golden 並文件記錄理由）** + 全基因組重驗新分佈
7. evidence_ledger + memory

## 7.5 實作驗證 log（分階段）

**Stage ① 強度分數（2026-06-16, 完成）— 改動合理正確確認**：
- **改動**：RegionProcessor `process_single_region` HP block 後加 strength_score 計算（hpp +2 欄 strength_score/strength_grade，emit +2 欄 StrengthScore/StrengthGrade）。保留 HeuristicScore（標 LEGACY）。
- **🐛 修正的 bug**：初版把計算放 `perform_clustering_and_significance`(:872 呼叫)，但 Δβ 在 process_single_region :897 之後才算 → germline_asm_dbeta=NaN → eff 項漏（C++ vs Python diff≈0.25×eff）。**移到 HP block 後**修正。
- **驗證（合理正確）**：① **C++ StrengthScore = Python 獨立重算 0/2624 不一致**（公式實作正確）② ctest **221/221** ③ regression 雙守護 **PASS**（StrengthScore/Grade 非 golden 9 欄，純加欄不動既有結果）④ StrengthGrade chr1 分佈 A53/B370/C849/D621/E731 合理（梯度，非飽和如 HeuristicScore）。
- **git**：feat/summary-nreadsvalid 累積（不污染 develop）；commit 後 feat 領先 develop +1 C++。

## 7. 風險 / 決策
- 🔴 改 VerificationClass = golden 欄大量變動 → **必更新 golden + 全基因組重驗**（pilot 已證 0 FP，但須 fresh 驗）
- within-HP pattern clustering 增 compute（per-region 多一次 HP-subset clustering）→ 須計時（預期 <2x，perm 非瓶頸）
- SEQC2 BED 路徑硬編碼風險 → 設 CLI 參數 + 缺檔 graceful
- 細類粒度（6+3）若太細 → 可 config 合併
