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

**Stage ④ 新 VerificationClass 判別（2026-06-16, 簡版完成）— 改動合理正確確認**：
- **改動**：`process_single_region` try-block 末（Δβ+Per-CpG 後，所有證據 final）override verification_class，原值存 verification_class_legacy（+1 欄）。判別=用戶哲學「任一結構證據→保留」：valid = Δβ_sig OR HP_AUC>=0.7 OR cluster_match(legacy Strong/Subclone) OR (Potential_LOH+結構)；細類 Strong/LabelShift/StructureNoLabel/LOH-Structure；noise 子型 Uniform/Chaotic/Uncorrelated（Epipoly+PairwiseMeanDist；NME 飽和不用）。**within-HP multigroup + SEQC2-LOH 簡版略過（Stage②③ 補）**。
- **為何 override 而非改 SignificanceAnalyzer**：legacy VC 在 SignificanceAnalyzer(:872 呼叫) 算，早於 Δβ(:897 後)→Δβ 不可用，故末段 override（同 strength 修法）。
- **驗證（合理正確）**：① regression diff **只第 9 欄 VerificationClass 變**(前 8 欄 RegionID/NumReads/NumCpGs/GlobalP/CramersV/PassedGating/ClusterPermanovaValid/Significant byte-identical；Weak→LabelShift/Noise→Noise_Uniform 皆預期) ② **--update-golden**(新 VC 為基準, SKIP+MAX_DIST 各 2624 @fd856e7) + 重跑雙守護 **PASS** ③ ctest **221/221** ④ 全基因組：**救回 3321/8180=41%**(簡版; 全版 46.6%), **legacy Strong→非valid FP=0**(純加識別無副作用) ⑤ 新分佈 Strong22310/LabelShift2541/StructureNoLabel404/LOH-Structure376/Noise{Uniform1559/Chaotic829/Uncorrelated2471}。
- **git**：feat 累積；commit 含更新 golden(text TSV 非 binary)。

**Stage ② within-HP substructure（pattern + level 軸）（2026-06-16, 完成）— 改動合理正確確認**：
- **改動**：`process_single_region` HP block 後加 within-HP 區塊（hpp +4 欄 within_hp1_ngroups/within_hp2_ngroups/within_hp_level_bimodal/within_hp_clean_multigroup + 2 method decl）。
  - **pattern 軸** `compute_within_hp_substructure`：對單一 germline-HP（tumor reads）的 per-CpG 距離子矩陣（複用 all_dist，NaN-peel）跑 build_tree + TreeCutter + silhouette；clean = silhouette>=0.5 + 平衡(最小群>=max(3,20%))。
  - **level 軸** `within_hp_level_clean`：distance 只抓 PATTERN，**漏 mean-β LEVEL 位移**（chr1 解構：distance 1.2% vs level 26.1%，level-only 25.6% 完全漏 — 遞迴盲點，見 `20260616_distance_misses_level_recurring_gap_inventory_01.md`）。level 軸 = per-read mean β 最佳 2-split，min群>=max(3,20%)+gap>0.15+var-reduction>0.5。
  - 合併 `within_hp_clean_multigroup = pattern OR level`，接進 Stage④ VC 的 valid 證據（MultiGroupNoLabel 細類）。
- **🐛 修正的 bug（MAX_DIST segfault — 根因優先 4-phase）**：Stage② pattern clustering 在 **MAX_DIST nan-distance-strategy 下 segfault**（SKIP 正常）。gdb 單/多執行緒重現 → 崩在 jemalloc `free()`（`BamReader::~BamReader`），faulting ptr=`0x4038000000000000`（= IEEE-754 double **24.0**）= **heap 被一個 double OOB 寫壞**。根因：`TreeCutter::find_optimal_clusters` 回 `best_k` 但其 labels **不保證連續 [0,best_k)** — `cut_by_num_clusters` 切到「>=k 群」高度，**MAX_DIST sentinel 距離造成大量 tied 合併高度 → 過切成 >k 群** → labels 值 >= best_k。我的 silhouette 迴圈用 `std::vector<double> other_sum(k)`／`sizes(k)` 以 raw label 直接索引 → **越界寫 double**。主 clustering 路徑免疫（用 `std::map` 不用固定大小 vector）。SKIP peel 掉 NaN → 子矩陣連續 → 從未觸發（純資料相依，故只 MAX_DIST 炸）。
- **修法（localized at consumer）**：取 labels 實際最大值 → remap 到連續 `[0,K)` → 全程用實際群數 K，絕不以 raw label 當索引界。只改 within-HP consumer，**不動 find_optimal_clusters / 主路徑**（regression-safe）。`out_ngroups=K`（實際群數）語意更正確；MAX_DIST tie 過切成微群 → 失平衡 → 不誤判 multigroup（保守）。
- **驗證（合理正確）**：① MAX_DIST 隔離跑 **EXIT 0**、significance_summary.csv **2624** 列（segfault 消失）② **ctest 221/221** ③ regression（未更新 golden）：**SKIP PASS 不變**（修法對連續 label 為 no-op）；MAX_DIST diff = **前 8 欄 byte-identical**、僅第 9 欄 VC：**607 → MultiGroupNoLabel**（256 Noise_Chaotic + 252 LabelShift + 46 Noise_Uniform + 35 StructureNoLabel + 18 Noise_Uncorrelated；其中 **320 從 Noise 救回**）④ **--update-golden**（SKIP+MAX_DIST 各 2624）+ 重跑 **雙守護 PASS** ⑤ chr1 within-HP CleanMultigroup **705/2624=26.9%**（level 26.2% 對齊 Python 解構 26.1%）⑥ 全基因組（SKIP, 修復 binary, 30490 TP SNV, 898s, 無 crash）：**WithinHP_CleanMultigroup 8077/30490=26.5%**（level 軸 7910/25.9% 為主力；對齊 chr1 26.9%/26.2% → 跨基因組成立）；**救回 4407/8180=53.9%**（分母 8180 = legacy Weak 6122+Noise 2058，與 Stage④ 同池；Stage④ pattern-only=41% → **level 軸再加 ~13pt 到 54%**，對齊 chr1 35%→50%）；**legacy Strong→非valid FP=0**（純加無副作用）；新 VC 分佈 Strong22310/MultiGroupNoLabel1899/LabelShift1809/LOH-Structure376/StructureNoLabel323/Noise{Uncorrelated1943/Uniform1275/Chaotic555}。
- **git**：feat/summary-nreadsvalid 累積；commit 含 cpp + 更新 golden(SKIP+MAX_DIST text TSV) + 本 spec。

## 7. 風險 / 決策
- 🔴 改 VerificationClass = golden 欄大量變動 → **必更新 golden + 全基因組重驗**（pilot 已證 0 FP，但須 fresh 驗）
- within-HP pattern clustering 增 compute（per-region 多一次 HP-subset clustering）→ 須計時（預期 <2x，perm 非瓶頸）
- SEQC2 BED 路徑硬編碼風險 → 設 CLI 參數 + 缺檔 graceful
- 細類粒度（6+3）若太細 → 可 config 合併
