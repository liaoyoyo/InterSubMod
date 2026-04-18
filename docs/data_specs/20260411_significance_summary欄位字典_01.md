<!--
建立時間: 2026-04-11 21:00
目標: significance_summary.csv 全部 59 欄完整欄位字典
處理範圍: C++ 核心輸出 significance_summary.csv 的每個欄位定義、資料型態、值域、生物學意義、計算來源
關聯檔案:
  - include/core/RegionProcessor.hpp (struct RegionResult)
  - src/core/RegionProcessor.cpp (write_significance_summary)
  - src/core/SignificanceAnalyzer.cpp (classification logic)
  - src/core/GlobalTest.cpp (gating logic)
  - docs/data_specs/20260411_master_dataset欄位字典_01.md
-->

# significance_summary.csv 欄位字典

> **版本**: v1.0 | **日期**: 2026-04-11 | **欄位數**: 59

**來源**: `build/bin/inter_sub_mod` C++ 核心處理，每次 run 產出一份，位於 `{output_dir}/significance_summary.csv`。每行代表一個 SNV region 的甲基化異質性分析結果。

**產出函式**: `RegionProcessor::write_significance_summary()` (`src/core/RegionProcessor.cpp:806`)

---

## 欄位分群總覽

| 群組 | 欄位編號 | 說明 |
|------|---------|------|
| [A. 區域識別](#a-區域識別) | #1-5 | SNV 位點座標與變異 |
| [B. 基本統計](#b-基本統計) | #6-7 | Read 深度與 CpG 數量 |
| [C. 全域檢驗](#c-全域檢驗) | #8-11 | Fisher 檢定 + Cramér's V + 啟發式評分 + 門檻過濾 |
| [D. 距離度量](#d-距離度量) | #12-13 | Read 間成對距離統計 |
| [E. 分群 PERMANOVA](#e-分群-permanova) | #14-18 | 無監督分群的結構檢驗 |
| [F. Stage 1: HP Merged 檢驗](#f-stage-1-hp-merged-檢驗) | #19-23 | HP1-family vs HP2-family 合併檢驗 |
| [G. Stage 2: HP Fine-Grained 檢驗](#g-stage-2-hp-fine-grained-檢驗) | #24-27 | 4 群 HP 精細檢驗 |
| [H. Stage 3: Allele 檢驗](#h-stage-3-allele-檢驗) | #28-30 | ALT vs REF 等位基因檢驗 |
| [I. Label×HP/Allele PERMANOVA](#i-labelhpallele-permanova) | #31-40 | 標籤分群結構檢驗（HP + Allele） |
| [J. Stage 4: Unassigned Affinity](#j-stage-4-unassigned-affinity) | #41-45 | HP3/HP0 未分配 read 歸屬 |
| [K. 判定欄位](#k-判定欄位) | #46-48 | 主導標籤 + 穩定性 + 驗證分類 |
| [L. 品質評估](#l-品質評估) | #49-55 | HP 比例 + LOH + 覆蓋度 + 品質分數 |
| [M. 局部檢驗與最終判定](#m-局部檢驗與最終判定) | #56-59 | 局部最佳分群 + 顯著性 + 過濾建議 |

---

## A. 區域識別

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 1 | **RegionID** | int | 0-indexed | 該次 run 中的 region 流水號 | `RegionResult::region_id` |
| 2 | **Chr** | string | chr1-chr22, chrX, chrY | 染色體名稱 | `chrom_index_.get_name(snv.chr_id)` |
| 3 | **Pos** | int | 1-based | SNV 在染色體上的位置（VCF POS） | `SNV::pos` |
| 4 | **Ref** | char | A/C/G/T | 參考等位基因 | `SNV::ref_base` |
| 5 | **Alt** | char | A/C/G/T | 替代等位基因 | `SNV::alt_base` |

---

## B. 基本統計

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 6 | **NumReads** | int | ≥1 | 涵蓋此 region (±window_size bp) 的有效 reads 數量，經 MAPQ≥20 + 長度≥500bp 過濾 | `RegionResult::num_reads` |
| 7 | **NumCpGs** | int | ≥0 | 分析視窗內的 CpG 位點數量，定義了甲基化矩陣的列數 | `RegionResult::num_cpgs` |

---

## C. 全域檢驗

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 8 | **GlobalP** | float | [0, 1] | 全域 Fisher Exact Test (FFH) p-value。測試甲基化狀態與 read 分群之間是否有顯著關聯。以 ALT allele 為分群依據。0 = 極顯著 | `global_alt.fisher_ffh.p_value` |
| 9 | **CramersV** | float | [0, 1] | Cramér's V 效應量。衡量甲基化模式與 read 標籤的關聯強度。**研究發現**: 93% 的值為零，因 2×2 列聯表結構限制（R1 結論） | `global_alt.cramers_v` |
| 10 | **HeuristicScore** | float | [0, ~20] | 綜合啟發式評分 = `-log10(best_p) + best_v × 2`。best_p 取 ALT/HP/HP-family 最小 p。有 dispersion warning 時 ×0.5，PERMANOVA 不顯著時 ×0.7 | `compute_heuristic_score()` |
| 11 | **PassedGating** | bool | true/false | 門檻過濾結果。通過條件：(ALT 或 HP-family) 的 p ≤ gating_threshold **且** V ≥ min_cramers_v。任一維度通過即為 true | `passed_gate` via `GlobalTest::apply_gating()` |

---

## D. 距離度量

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 12 | **PairwiseMeanDist** | float | [0, 1] | 所有有效 read pairs 的甲基化距離平均值（使用指定的 distance metric，預設 Bernoulli） | `summarize_pairwise_distances()` mean |
| 13 | **PairwiseMedianDist** | float | [0, 1] | 所有有效 read pairs 的甲基化距離中位數。比 mean 更穩健，不受極端值影響 | `summarize_pairwise_distances()` median |

---

## E. 分群 PERMANOVA

Cluster-First 路徑：先做無監督階層式分群，再以 PERMANOVA 檢驗分群結構是否具統計意義。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 14 | **ClusterPermanovaF** | float | ≥0 | PERMANOVA pseudo-F 統計量。越大表示群間差異相對群內差異越大 | `StructureTest` cluster-based |
| 15 | **ClusterPermanovaP** | float | [0, 1] | PERMANOVA 置換檢定 p-value（999 次置換） | `StructureTest` cluster-based |
| 16 | **ClusterPermanovaValid** | bool | true/false | PERMANOVA 是否有效執行（需 ≥2 群且每群 ≥2 reads） | `StructureTest` validity |
| 17 | **ClusterDispersionP** | float | [0, 1] | 群間離散度檢定 p-value（PERMDISP2）。顯著時表示群間變異度不均，PERMANOVA F 可能膨脹 | `StructureTest` dispersion |
| 18 | **ClusterDispersionWarn** | bool | true/false | 離散度警告（DispersionP ≤ 0.05）。若為 true，PERMANOVA 結果需謹慎解讀 | `StructureTest` |

---

## F. Stage 1: HP Merged 檢驗

多階段 HP 驗證的第一階段：將同一 haplotype 的 germline + somatic reads 合併（HP1+HP1-1 vs HP2+HP2-1），檢驗兩個 haplotype family 之間的甲基化差異。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 19 | **HPMergedDelta** | float | [0, ~1] | HP1-family 與 HP2-family 之間的平均甲基化距離差 (d_between - d_within)。正值=群間距離大於群內 | `hp_multistage.merged_delta` |
| 20 | **HPMergedP** | float | [0, 1] | 置換檢定 p-value（999 次）。檢驗 HP 合併分群是否具顯著甲基化差異 | `hp_multistage.merged_p` |
| 21 | **HPMergedSig** | bool | true/false | p ≤ 0.05 且 delta > 0 | `hp_multistage.merged_sig` |
| 22 | **HP1FamilyN** | int | ≥0 | HP1-family 的 read 數量（HP1 + HP1-somatic）。**研究發現**: TO 模式中因 self-phasing bias 導致 94.6% 偏向 HP1 | `hp_multistage.n_hp1_family` |
| 23 | **HP2FamilyN** | int | ≥0 | HP2-family 的 read 數量（HP2 + HP2-somatic） | `hp_multistage.n_hp2_family` |

---

## G. Stage 2: HP Fine-Grained 檢驗

第二階段：將 reads 分為最多 4 群（HP1-germline, HP1-somatic, HP2-germline, HP2-somatic），以 PERMANOVA 檢驗群間差異。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 24 | **HPFineF** | float | ≥0 | Fine-grained PERMANOVA pseudo-F 統計量 | `hp_multistage.fine_f` |
| 25 | **HPFineP** | float | [0, 1] | Fine-grained PERMANOVA p-value | `hp_multistage.fine_p` |
| 26 | **HPFineSig** | bool | true/false | p ≤ 0.05 | `hp_multistage.fine_sig` |
| 27 | **HPFineNGroups** | int | 0-4 | 有效群組數（每群需 ≥3 reads）。**研究發現**: N≥4 且 NumReads≥80 → TP rate 89.1%，是 somatic heterogeneity marker，但不能用於 variant filter | `hp_multistage.fine_n_groups` |

---

## H. Stage 3: Allele 檢驗

第三階段：以 ALT/REF allele 支持狀態分群，檢驗等位基因間的甲基化差異。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 28 | **AlleleDelta** | float | [0, ~1] | ALT-supporting vs REF-supporting reads 之間的平均甲基化距離差。**研究發現**: 與 caller AF 高度相關，是 AF confound 而非獨立信號（O12） | `hp_multistage.allele_delta` |
| 29 | **AlleleP** | float | [0, 1] | Allele 分群置換檢定 p-value | `hp_multistage.allele_p` |
| 30 | **AlleleSig** | bool | true/false | p ≤ 0.05 且 delta > 0 | `hp_multistage.allele_sig` |

---

## I. Label×HP/Allele PERMANOVA

Label-First 路徑：使用已知標籤（HP tag、Allele）作為分群依據，以 PERMANOVA 檢驗標籤是否解釋甲基化變異。

### HP 維度

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 31 | **LabelHPPermanovaF** | float | ≥0 | 以 HP tag 為分群的 PERMANOVA pseudo-F | `label_hp_permanova` |
| 32 | **LabelHPPermanovaP** | float | [0, 1] | HP 標籤 PERMANOVA p-value | `label_hp_permanova` |
| 33 | **LabelHPPermanovaValid** | bool | true/false | HP PERMANOVA 是否有效（≥2 群且每群 ≥2 reads） | `label_hp_permanova` |
| 34 | **LabelHPDispersionP** | float | [0, 1] | HP 群間離散度 p-value | `label_hp_dispersion` |
| 35 | **LabelHPDispersionWarn** | bool | true/false | HP 離散度警告 | `label_hp_dispersion_warning` |

### Allele 維度

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 36 | **LabelAllelePermanovaF** | float | ≥0 | 以 Allele 為分群的 PERMANOVA pseudo-F | `label_allele_permanova` |
| 37 | **LabelAllelePermanovaP** | float | [0, 1] | Allele PERMANOVA p-value | `label_allele_permanova` |
| 38 | **LabelAllelePermanovaValid** | bool | true/false | Allele PERMANOVA 是否有效 | `label_allele_permanova` |
| 39 | **LabelAlleleDispersionP** | float | [0, 1] | Allele 離散度 p-value | `label_allele_dispersion` |
| 40 | **LabelAlleleDispersionWarn** | bool | true/false | Allele 離散度警告 | `label_allele_dispersion_warning` |

---

## J. Stage 4: Unassigned Affinity

第四階段：評估未分配 reads (HP3 = multi-mapped, HP0 = unphased) 與哪個 haplotype 的甲基化模式更接近。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 41 | **UnassignedAffinity** | float | [0, 1] | 未分配 reads 對 HP1/HP2 的 affinity score。越高表示偏好某一 haplotype | `unassigned_affinity` |
| 42 | **UnassignedAffinityP** | float | [0, 1] | Affinity 置換檢定 p-value | `unassigned_affinity_p` |
| 43 | **UnassignedDir** | string | HP1/HP2/None | Affinity 方向：未分配 reads 更接近哪個 haplotype | `unassigned_dir` |
| 44 | **NHP3** | int | ≥0 | HP3 reads 數量（multi-mapped，同時屬於多個 haplotype） | `unassigned_n_hp3` |
| 45 | **NHP0** | int | ≥0 | HP0 reads 數量（unphased，無 haplotype 資訊） | `unassigned_n_hp0` |

---

## K. 判定欄位

綜合多階段檢驗結果的分類判定。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 46 | **DominantLabel** | string | hp / allele / none | 主導差異來源。hp=haplotype 間差異最大（ASM），allele=等位基因間差異最大，none=無顯著差異 | `dominant_label` |
| 47 | **Stability** | float | [0, 1] | 分群穩定性分數（subsampling 驗證）。目前版本固定為 0.0（未啟用） | `cluster_stability` |
| 48 | **VerificationClass** | string | Strong / Subclone / Weak / Noise | 雙向驗證分類，見下方決策邏輯 | `SignificanceAnalyzer` Phase 5 |

### VerificationClass 決策邏輯

```
if label_sig AND cluster_significant:
    → "Strong"      # 雙路徑一致：Label-First + Cluster-First 皆顯著
elif NOT label_sig AND cluster_significant:
    → "Subclone"    # 有分群結構但無標籤相關：可能是亞克隆差異
elif label_sig AND NOT cluster_significant:
    → "Weak"        # 有標籤信號但分群結構弱
else:
    → "Noise"       # 兩路徑皆不顯著
```

其中 `cluster_significant` = `passed_gate && (global_alt_p ≤ 0.05 || global_hp_family_p ≤ 0.05)`

---

## L. 品質評估

Multi-Layer Validation Quality Assessment（Phase 5 新增）。

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 49 | **HP_Ratio** | float | [0, 1] | HP1-family / (HP1-family + HP2-family)，含 Laplace 平滑。0.5=平衡，接近 0 或 1 表示 LOH-like。**注意**: TO 模式因 self-phasing 導致系統性偏移 | `compute_hp_ratio()` |
| 50 | **Potential_LOH** | bool | true/false | HP_Ratio < 0.1 或 > 0.9。指示可能的 LOH（Loss of Heterozygosity）區域。**TO 模式警告**: 因 self-phasing artifact，62% 的 LOH 判定會消失（因果鏈已驗證） | `is_potential_loh()` |
| 51 | **Coverage_Multiple** | float | >0 | NumReads / diploid_coverage。diploid_coverage 由 KDE mode 估計（所有 region 的 NumReads 分布取眾數），fallback=75.0 | `compute_coverage_multiple()` |
| 52 | **Coverage_Category** | string | 見下表 | 基於 Coverage_Multiple 的分類 | `determine_coverage_category()` |
| 53 | **LOH_Subtype** | string | 見下表 | 結合 Potential_LOH 與 VerificationClass 的 LOH 亞型 | `determine_loh_subtype()` |
| 54 | **Quality_Score** | float | [0, 100] | 綜合品質分數，從 100 開始扣分/加分。**模式差異**: Paired 使用全部權重，TO 停用 LOH penalty (=0) 和 verify bonus (=0) | `compute_quality_score()` |
| 55 | **Quality_Tier** | string | High / Medium / Low | 品質等級：High (≥70), Medium (40-69), Low (<40) | `determine_quality_tier()` |

### Coverage_Category 對應表

| Coverage_Multiple | Category | 生物學意義 |
|-------------------|----------|-----------|
| < 0.5 | CNV_Loss | 拷貝數缺失，可能整段刪除 |
| 0.5-0.8 | Low | 低覆蓋，統計力不足 |
| 0.8-1.2 | Normal | 正常二倍體覆蓋 |
| 1.2-1.5 | Elevated | 略高，可能有低倍數增益 |
| 1.5-2.0 | CNV_Gain | 拷貝數增加 |
| > 2.0 | High_Copy | 高拷貝數，可能是 amplification |

### LOH_Subtype 對應表

| Potential_LOH | VerificationClass | LOH_Subtype | 意義 |
|---------------|-------------------|-------------|------|
| false | any | None | 非 LOH |
| true | Noise | LOH_Noise | LOH 區域但無信號 |
| true | Weak | LOH_Weak | LOH 區域 + 弱信號 |
| true | Strong | LOH_Strong | LOH 區域 + 強信號 |
| true | Subclone | LOH_Subclone | LOH 區域 + 亞克隆結構 |

### Quality_Score 權重表

| 項目 | 條件 | Paired 權重 | TO 權重 | 說明 |
|------|------|-------------|---------|------|
| read_penalty_severe | NumReads < 30 | -20 | -20 | 嚴重深度不足 |
| read_penalty_moderate | NumReads < 50 | -10 | -10 | 深度偏低 |
| cpg_penalty_severe | NumCpGs < 20 | -15 | -15 | CpG 不足 |
| cpg_penalty_moderate | NumCpGs < 30 | -10 | -10 | CpG 偏少 |
| cov_penalty_cnv_loss | Coverage_Multiple < 0.5 | -30 | -30 | CNV 缺失 |
| cov_penalty_low | Coverage_Multiple < 0.8 | -15 | -15 | 低覆蓋 |
| cov_penalty_high_copy | Coverage_Multiple > 2.0 | -20 | -20 | 高拷貝 |
| **loh_penalty** | Potential_LOH = true | **-25** | **0** | TO 停用因 self-phasing artifact |
| **verify_bonus** | HPMergedSig ∧ AlleleSig | **+10** | **0** | TO 停用因 FP > TP 反效果 |
| effect_bonus_strong | CramersV ≥ 0.3 | +15 | +15 | 強效應量 |
| effect_bonus_moderate | CramersV ≥ 0.2 | +5 | +5 | 中效應量 |

---

## M. 局部檢驗與最終判定

| # | 欄位名 | 型態 | 值域 | 說明 | C++ 來源 |
|---|--------|------|------|------|---------|
| 56 | **LocalBestCluster** | int | -1, 0, 1, ... | 局部檢驗中最顯著的分群 ID。-1=無有效分群 | `local_best_cluster` |
| 57 | **LocalBestP** | float | [0, 1] | 最佳局部分群的 p-value | `local_best_p_value` |
| 58 | **Significant** | bool | true/false | **最終顯著性判定**（寫入時計算）：`PassedGating && GlobalP ≤ 0.05 && CramersV ≥ 0.1 && NumReads ≥ 20` | Write-time computation |
| 59 | **SuggestFilter** | bool | true/false | 建議過濾標記（F1 優化用）：`label_delta > 0.3`。**注意**: label_delta 本身不在 CSV 中，它是 deprecated 的 Label-First delta | Write-time computation |

---

## 附註

### 生產檔案 vs 程式碼 header 差異

目前 C++ source (`src/core/RegionProcessor.cpp:818-842`) 的 header 已更新為包含額外的 multi-layer HP 欄位（`GlobalP_HPFamily`, `CramersV_HPFamily`, `GlobalP_HPFine`, `CramersV_HPFine`, `HPFine_NGroups_CF`, `HPFineN_HP1`, `HPFineN_HP1S`, `HPFineN_HP2`, `HPFineN_HP2S`, `HPFineD_HP1_HP1S` 等 6 個 pairwise distance），但 **當前生產 canonical run 的檔案仍為 59 欄**。未來重新 build + 重跑後會產出擴展版本。

### 已知限制（研究結論）

1. **CramersV**: 93% 值為零 — 2×2 列聯表結構限制（R1 結論）
2. **HPMergedDelta**: 多群場景下可能反向（R2 結論）
3. **AlleleDelta**: 與 caller AF 高度共線，是 confound 而非獨立信號（O12 結論）
4. **HP_Ratio / Potential_LOH**: TO 模式受 self-phasing 嚴重影響，62% 為 artifact（因果鏈驗證）
5. **Quality_Score**: TO 模式 AUC=0.497（等同隨機），Paired 模式有效（O2 結論）
6. **HPFineNGroups**: N≥4 是 somatic heterogeneity marker，但不能作為 variant filter（R5/Beyond-AUC 結論）
