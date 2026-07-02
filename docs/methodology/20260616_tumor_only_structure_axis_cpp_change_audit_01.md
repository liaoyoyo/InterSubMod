<!--
建立時間: 2026-06-16
問題類型: C++ 方法 + 評分設計
影響 track: paired (tumor+normal)
狀態: pending_decision
build_branch: feat/summary-nreadsvalid @ 990175d
data_sources: docs/methodology/_assets/20260616_readset_provenance/readset_stats.json, docs/methodology/_assets/20260616_readset_provenance/tumor_intrinsic_ref.json, docs/methodology/_assets/20260616_readset_provenance/tumor_cluster_pilot_agg.json
-->

# tumor-only 結構軸 + StrengthScore tumor 權重 — C++ 變更審查

## §1 問題描述（IDENTIFY）

**問題 1 — 判定軸全在混池上跑，無 tumor-only 結構軸。**
主距離矩陣 / clustering / GlobalTest / PERMANOVA 用 `read_list = get_reads()` = tumor+normal 全部 reads（`RegionProcessor.cpp:848,866-872`）。只有 within-HP（`:1169`）、HP_AUC_Tumor（`:2237`）、Tumor_HP_Delta（`:970`）是純 tumor。**判定的「有結構」無法自證是 tumor-intrinsic 還是 tumor-vs-normal 混池假象（sampleOrphan）。**

**問題 2 — StrengthScore 0% tumor 權重。**
`RegionProcessor.cpp:1196-1204`:
```
struct_s = (hp_auc_all−0.5)/0.5     # all-pool（normal-dominated）
eff_s    = |germline_asm_dbeta|/0.3 # GERMLINE!
assoc_s  = cramers_v / sig_s = −log10(global_p)/4   # all-pool
strength = 0.30·struct + 0.25·eff + 0.25·assoc + 0.20·sig
```
強度完全由 all-pool + germline 訊號構成，tumor 特異結構權重 = 0。

## §2 量化影響（QUANTIFY；本 session 全基因組 HCC1395 30,490 TP）

| 證據 | 數字 | 來源 |
|---|---|---|
| Significant 只 48.0% 有純 tumor 軸佐證；28.6% 純混池無佐證 | 13,605 / 8,103 | readset_stats.json |
| germline-HP 軸 normal 比 tumor 強 | HP_AUC median tumor 0.641 < normal 0.788；70.7% normal 較強；tumor-intrinsic 僅 1.6% | tumor_intrinsic_ref.json |
| tumor-only clustering **可行** 且結構**不崩** | 100% 位點 ≥6 tumor reads；抽掉 normal collapse 僅 2%；tum≥pool 61% | tumor_cluster_pilot_agg.json |
| 但 raw silhouette **過度偵測**（Noise 也高）| Noise_Uniform sil 0.613 > Significant 0.475 | tumor_cluster_pilot_agg.json |

**結論**：tumor-only clustering 可行、結構 tumor-intrinsic（值得做），但需 **PERMANOVA permutation null** 才能 gate 真結構（crude silhouette 不夠）。這只有 C++ 忠實複製既有 PERMANOVA 才做得對。

## §3 修改選項（OPTIONS）

### 方案 A：不修改（接受現狀）
判定維持混池，StrengthScore 維持 0% tumor 權重，補觀察文件。

### 方案 B：加 tumor-only PERMANOVA 軸 + TumorIntrinsic flag + StrengthScore tumor 項（**推薦**）
**新函式** `compute_tumor_only_cluster_structure`（模板 = `compute_within_hp_substructure` `:2079`）：
1. 抽 tumor 子距離矩陣（`tumor_indices` 已存在於 `:953` / extract_and_test_hp）
2. `extract_complete_submatrix_indices` NaN peel-off（已存在）
3. `HierarchicalClustering::build_tree` + `TreeCutter::find_optimal_clusters` + label remap（segfault fix 已內建 `:2114`）
4. **`StructureTest::run_permanova(tumor_dist, cluster_labels)`**（public，`:70`，複用）→ tumor cluster PERMANOVA F/p/valid
5. **dispersion check**（`check_dispersion`，D2 已加 analytic-F）→ tumor_dispersion_warn
6. silhouette（複用 `:2126` 模式）

**新欄位**：`TumorOnlyClusterK / TumorOnlySilhouette / TumorOnlyPermanovaF / TumorOnlyPermanovaP / TumorOnlyPermanovaValid / TumorOnlyDispersionWarn / TumorIntrinsic`
**split 維度**：`TumorIntrinsic = TumorOnlyPermanovaValid && p<0.05 && !TumorOnlyDispersionWarn`（tumor 自己就有 clean location 結構）
**StrengthScore 改**（見下 D2 子選項）：加 `tumor_s` 項。
修改位置：`RegionProcessor.cpp`（call site ~`:1166` 區、StrengthScore `:1196`、CSV `:1432`、header 欄位）。**不改 Significant 定義**（純加 flag + 強度）。

### 方案 C：方案 B + 重定義 Significant 納入 TumorIntrinsic
在 B 基礎上，把 tumor-intrinsic 變成 valid 軸之一 / 或反過來把「純混池無 tumor 佐證」降級。**較大變動、影響回歸基準大**。

---

## §3.1 SWOT

### 方案 B
|  | Helpful | Harmful |
|---|---|---|
| **Internal** | **S** 全複用既有 module（run_permanova/build_tree/find_optimal_clusters/segfault-fix），零重寫；不改 Significant 故回歸風險低；補上唯一缺的 null-gated tumor 軸 | **W** 增加 7 欄 + 一個 clustering pass（per-region ~tumor reads，輕量）；StrengthScore 權重需調參 |
| **External** | **O** 直接服務 subclone 論文（tumor-intrinsic = somatic 候選的必要條件）；對齊 sampleOrphan confound 修正方向 | **T** 單樣本，tumor-intrinsic≠somatic（仍含 germline within-tumor + cis-ASM），需 cross-sample + cis-control 後續 |

### 方案 C
|  | Helpful | Harmful |
|---|---|---|
| **Internal** | **S** 判定直接反映 tumor-intrinsic | **W** 改 Significant → 回歸基準大變、雙守護 golden 需大改、與既有 D1 統一判定耦合風險 |
| **External** | **O** 一步到位的 somatic-aware 判定 | **T** 過早把未跨樣本驗證的 tumor-intrinsic 寫進判定，違 single-sample 紀律；不可逆性高 |

→ 方案 C 的 W+T > S+O（不可逆 + 回歸大 + 紀律風險）→ **降為條件性方案**（待 B 落地 + cross-sample 後再議）。

## §4 驗收標準
- [ ] ctest 全綠（221/221）
- [ ] 雙守護 regression：新欄位為**附加**（既有 9 欄 col1-9 byte-identical）→ 用 `--update-golden` 重建含新欄基準後雙守護 PASS
- [ ] tumor-only PERMANOVA p 對 scipy 交叉驗證（複用 D2 的 f_dist_sf 驗證法）
- [ ] 全基因組重跑：TumorIntrinsic 比例合理（與 pilot 量級對照，但 PERMANOVA-gated 應 < pilot silhouette 的 96%）
- [ ] StrengthScore 改後分佈合理（連續、未飽和、tumor 項真有貢獻）
- [ ] FP 安全：legacy Strong 不因加 flag 被降級（不改 Significant 故應 =0）

## §5 用戶決策（DECIDED 2026-06-16）
**D1 — 範圍**：**[X] B（flag + 強度，不改 Significant）**
**D2 — StrengthScore 公式**：**[X] B2（加 tumor 項 + somatic-residual 取代 germline）**
**新 StrengthScore 公式（D3 修正 — 等分 + 元件輸出，2026-06-16）**：
> 敏感度分析（weight_sensitivity.json）：equal vs 判斷權重 Spearman ρ=**0.9976**、top-1000 重疊 94.6% → **權重幾乎不影響排序**。原因：`sig` 元件飽和（median=1.0 不鑑別）。無 ground-truth label（F1 DEAD/SEQC2 不可輸入）→ 無法校準「最佳」權重，robustness 是唯一可做的驗證，而它說等分就夠。
```
# 移除飽和 sig，5 元件等權（各 0.20），每元件 sub-score 輸出成欄
strength = ( struct(hp_auc_all) + tumor(silhouette, tumor_intrinsic-gated)
           + somatic(|somatic_residual_dbeta|/0.3) + assoc(cramers_v)
           + germline(|germline_asm_dbeta|/0.3) ) / 5
輸出欄：StrengthStruct/StrengthTumor/StrengthSomatic/StrengthAssoc/StrengthGermline + StrengthScore
```
**日期**：2026-06-16
**理由**：tumor-only PERMANOVA 補唯一缺的 null-gated tumor 軸；不改 Significant 降回歸風險；somatic 取代 germline 對齊 subclone；**等分（實證權重無差）+ 元件輸出使強度可觀察驗證、可下游重加權**。C 條件性方案待 cross-sample。加 tumor/somatic 後重跑敏感度。
→ 進 `/cpp-change` 6 步實作。
