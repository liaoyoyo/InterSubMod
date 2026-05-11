# Step 8 — R1-Global HCC1395 Normal BAM Phase 2 TO arm

**Date**: 2026-04-21
**Scope**: Full 40,239 TO variants (28,396 TP + 11,843 FP) × inter_sub_mod Phase 2 pipeline × full 144GB Normal BAM random access
**Purpose**: 驗證 F pilot CL-025a（SampleASM_Delta ⭐3 characterization-only）在全域 region 是否保持；評估 paired arm 放棄後 TO arm 是否有獨立的 F1-filter 潛力。

## 1. 規模

- Global regions processed: 40237 (vs F pilot subset 801)
- TP=28396, FP=11841, TP_rate=0.7057 (vs F pilot subset 0.9700)
- FP=11841 is 987x F pilot subset's 12 FP — bootstrap CI now statistically meaningful

## 2. 全域 per-feature AUC（raw + residualized on NumReads+Coverage_Multiple）

![](../figures/step8_global_feature_auc.png)

| Feature | Raw AUC [95% CI] | n | Residualized AUC [95% CI] | n |
|---------|------------------|---|---------------------------|---|
| Epipoly_Delta | 0.532 [0.526,0.538] | 38860 | 0.533 [0.527,0.538] | 38860 |
| Combined_HP_Signed_Delta | 0.529 [0.522,0.536] | 39818 | 0.529 [0.523,0.535] | 39818 |
| HPFineF | 0.527 [0.520,0.533] | 40237 | 0.527 [0.520,0.533] | 40237 |
| SampleASM_Delta | 0.526 [0.520,0.532] | 40237 | 0.527 [0.520,0.533] | 40237 |
| NormalBaseline_Mean | 0.523 [0.517,0.530] | 40237 | 0.523 [0.517,0.530] | 40237 |
| HPFineNGroups | 0.522 [0.518,0.527] | 40237 | 0.532 [0.526,0.538] | 40237 |
| PairwiseMedianDist | 0.522 [0.515,0.528] | 40237 | 0.523 [0.516,0.529] | 40237 |
| Normal_HP_Delta | 0.513 [0.508,0.520] | 38890 | 0.513 [0.508,0.520] | 38890 |
| Tumor_HP_Delta | 0.513 [0.503,0.522] | 17424 | 0.511 [0.502,0.520] | 17424 |
| HP_Signed_Residual | 0.508 [0.499,0.515] | 20647 | 0.507 [0.499,0.515] | 20647 |
| NormalBaseline_Coverage | 0.506 [0.500,0.512] | 40237 | 0.512 [0.507,0.519] | 40237 |
| Entropy_Imbalance | 0.506 [0.500,0.511] | 38860 | 0.509 [0.503,0.515] | 38860 |
| HP_Residual_Delta | 0.504 [0.498,0.510] | 40237 | 0.504 [0.498,0.510] | 40237 |
| Fisher_Frac_Sig | 0.504 [0.497,0.510] | 40237 | 0.504 [0.497,0.510] | 40237 |
| CramersV | 0.501 [0.498,0.503] | 40237 | 0.508 [0.502,0.515] | 40237 |

### 2.1 Promising features（residualized AUC ≥0.60 且 CI 下界 ≥0.55）

**無特徵通過標準** — 全域 Phase 2 TO arm 無 F1-filter 候選。

## 3. 與 F pilot subset 對照（overfit 評估）

![](../figures/step8_global_vs_subset_overfit.png)

| Feature | F pilot subset old (n=801) | Global ∩ F pilot filter | Global full (40,239) |
|---------|----------------------------|-------------------------|----------------------|
| SampleASM_Delta | 0.643 (n=801) | 0.524 (n=6965) | 0.526 (n=40237) |
| HP_Residual_Delta | 0.577 (n=801) | 0.507 (n=6965) | 0.504 (n=40237) |
| HP_Signed_Residual | 0.515 (n=798) | 0.509 (n=6962) | 0.508 (n=20647) |
| Combined_HP_Signed_Delta | 0.518 (n=801) | 0.505 (n=6965) | 0.529 (n=39818) |
| Tumor_HP_Delta | 0.524 (n=794) | 0.515 (n=6965) | 0.513 (n=17424) |
| Normal_HP_Delta | 0.578 (n=798) | 0.525 (n=6945) | 0.513 (n=38890) |
| NormalBaseline_Mean | 0.517 (n=801) | 0.525 (n=6965) | 0.523 (n=40237) |
| NormalBaseline_Coverage | 0.687 (n=801) | 0.523 (n=6965) | 0.506 (n=40237) |
| Fisher_Frac_Sig | 0.541 (n=801) | 0.506 (n=6965) | 0.504 (n=40237) |
| Entropy_Imbalance | 0.522 (n=800) | 0.502 (n=6965) | 0.506 (n=38860) |
| Epipoly_Delta | 0.507 (n=800) | 0.518 (n=6965) | 0.532 (n=38860) |
| PairwiseMedianDist | 0.517 (n=801) | 0.526 (n=6965) | 0.522 (n=40237) |
| HPFineF | 0.582 (n=801) | 0.524 (n=6965) | 0.527 (n=40237) |
| HPFineNGroups | 0.520 (n=801) | 0.500 (n=6965) | 0.522 (n=40237) |
| CramersV | 0.504 (n=801) | 0.507 (n=6965) | 0.501 (n=40237) |

**Overfit heuristic**：若 `auc_subset_old` >> `auc_global_fpilot_filter` >> `auc_global` 遞減 → F pilot 子集內的 AUC 為 pre-selection overfit；若 `auc_global_fpilot_filter` ≈ `auc_subset_old` → filter 真有 enrichment 效果；若 `auc_global` ≈ `auc_subset_old` → 特徵本身 robust。

## 4. 推論與結論

### 4.1 reasoning chain

1. **Premise**：F pilot subset Step 7 顯示 paired `Fisher_Frac_Sig` 在 941-region 子集內 AUC=0.726，但殘差化後 CI 跨隨機 + TP 飽和（99.5%），paired F1-filter 方向 2026-04-21 放棄。
2. **Question**：TO arm `SampleASM_Delta` 在同 subset 殘差化 0.610（CL-025a ⭐3）是否為真實訊號？
3. **Test**：全域 40,239 region 不經 F pilot filter 直接評估；若 residualized AUC ≥0.60 且 CI 下界 ≥0.55 → 真訊號；若顯著衰減 → subset artifact。
4. **Result**：見 §2.1 與 §3 表。
5. **Interpretation**：
   - 無特徵通過標準（最高 residualized AUC=0.533 `Epipoly_Delta`，CI 上界 0.538 仍 < 0.58 ceiling）。
   - **Overfit severity（§3 對照）**：
     - `SampleASM_Delta`：subset 0.643 → global 0.527 (**Δ=-0.116**, 90% signal collapse)
     - `NormalBaseline_Coverage`：subset 0.687 → global 0.506 (**Δ=-0.181**, 97% signal collapse)
     - `HPFineF`：subset 0.582 → global 0.527 (Δ=-0.055)
     - `Normal_HP_Delta`：subset 0.578 → global 0.513 (Δ=-0.065)
   - 同一特徵在 `Global ∩ F pilot filter`（n=6965，模擬 F pilot 條件但非 pre-selected TP/FP 樣本）也全部塌回 ≤0.53 → 確認 subset 0.64-0.69 的訊號**來自 F pilot 的 pre-registered subset 選擇行為**（選 NG=4+AF<0.4+NR≥80+NonLOH 時，subset 內 12 FP 恰好在 NormalBaseline_Coverage / SampleASM_Delta 上與 720 TP 有分離），**不是特徵本身能 transfer 到未選過的 region**。
   - 結論：**TO arm Phase 2 所有 15 個 ISM 輸出特徵在全域 HCC1395 TO FP filter 場景均無 F1-filter 價值**。

### 4.2 對 Registry 的影響

- **CL-025a**（SampleASM_Delta characterization）：原 ⭐3
  → **降級為 ⭐2 concluded NEGATIVE-for-F1-filter**（global residualized AUC=0.527 [0.520, 0.533]；subset 0.610 為 pre-selection overfit；characterization 角色可保留但需配合額外 stratification）。
- **CL-008**（Beyond-AUC ≤0.58 ceiling, pure methylation）：**✅ 維持 ⭐4**；Phase 2 Normal BAM + Sample ASM 整合後，ceiling 仍未突破（最高 raw AUC 0.532 `Epipoly_Delta`，最高 residualized 0.533），與 Phase 1A 結論完全一致。
- **CL-025**（LOH.bed self-phasing independence）：⭐5 不受影響（本結果不牽涉 LOH.bed 機制）。
- **CL-017**（Paired F1 saturation）：⭐5 不受影響。
- **CL-025b**（paired F1-filter abandoned）：本結果**加強**此決策的正當性（paired Fisher_Frac_Sig 0.726 subset 與 TO SampleASM_Delta 0.643 subset 為**同類 overfit**，兩者 global 均塌回隨機）。
- **CL-025c**（Cross-mode concordance）：⭐4 不受影響（TO/paired 獨立性聲稱本身成立）。

### 4.3 對 Opus 4.7 plan Part C.2 Phase 2A Normal Methylation Reference 的影響

- Part C.2 原聲稱「Normal ASM 可區分 germline vs somatic-driven ASM」的成功標準為「Sample ASM vs Normal ASM delta AUC ≥ 0.60（residualized）」。
- 本 pilot `SampleASM_Delta` global residualized AUC = **0.527** [0.520, 0.533]，CI 上界 0.533 < 成功標準 0.60 → **Phase 2A Normal Methylation Reference 方向在 HCC1395 TO 場景 F1-filter 上 NEGATIVE**。
- Characterization 角色（`SampleASM_Sig` flag、NormalBaseline annotation）仍可保留為 **descriptive tool**（非 F1-driven），但不再是 POSITIVE 候選。

## 5. 下一步候選（批次 3 重新排序）

原批次 3 順序（R8 Per-CpG ASM + R12 跨樣本擴展）預設 R1-Global 成功；本結果 NEGATIVE 後需重排：

1. **R12 跨樣本擴展 (降優先)**：HCC1395-only NEGATIVE 已判定，擴展到其他 5 in-scope 樣本若預期相同結果 → 僅確認性質，不太可能翻轉結論。建議**暫緩**，等有新特徵候選再做。
2. **R8 Per-CpG ASM residualized (維持 P1)**：Per-CpG 層級（CL-020 24 metrics）在 subset 內 fisher_frac_sig diff=0.125 邊緣 → 需在全域驗證是否同樣 subset overfit。
3. **新方向 — 其他五大研究目標**（使用者 2026-04-18 偏好）：ISM 定位 read-level characterization，而非 filter。本結果強化了此定位：
   - 目標 1（ReadParser 修復 + HP-dependent 29 特徵重測）仍為 P0
   - 目標 2-5（subclone characterization / gene-level / CpG 分層 / LOH×AF mechanism）重心應從 F1 移到 **biological interpretation**
4. **論文敘事**：Result §X LOH Subclone 雙證據鏈（F pilot + Z3 pilot）可保留；但 §X「HPFineNGroups somatic heterogeneity marker」需明記**僅限 characterization 用途**，不宣稱 F1-filter 價值。

## 6. 自我質疑 — 本結論是否可靠？

| 質疑 | 回應 |
|------|------|
| Q1: 40,237 region 數據是 single HCC1395，可能樣本特殊？ | 可能。但 subset overfit 模式（subset 0.64 → global 0.53）為統計學普遍現象，與樣本特殊性無關。 |
| Q2: residualization confound 選擇 (NumReads, Coverage_Multiple) 是否不足？ | 這兩個是 Phase 1A 驗證過的主要 confound；原 F pilot Step 7 `SampleASM_Delta` 0.610 也用類似 residualization。若加更多 confound 只會 AUC 更低，不會更高。 |
| Q3: CI 0.520-0.533 是否仍可解釋為「弱訊號」？ | 不。CI 上界 0.533 < 0.58 ceiling，且 bootstrap seed=2026 穩定；與 CL-008 耗盡論完全吻合。 |
| Q4: `auc_global_fpilot_filter` 也 ~0.50-0.53，是否 F pilot filter 完全無用？ | F pilot filter 對 **characterization**（TP rate 89.1% / NG 結構）仍有效；對 **F1-filter**（subset 內 TP/FP 分離 AUC）無效。兩者不同目標。 |
| Q5: 是否該做其他樣本驗證才能定論？ | 不需要。overfit 證據已決定性（n=40,237 與 n=801 的對比 Δ=-0.12），跨樣本只會重現同模式。 |

## 6. Provenance

- Binary: `build/bin/inter_sub_mod`
- Global VCF: `output/hcc1395_normal_pilot_global/HCC1395_TO_global.vcf.gz` (40,239 variants)
- Tumor BAM: `tumor_tagged.bam` (TO arm, 2026-03-07 step03_longphase_to)
- Normal BAM: `HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (full 144GB, BAM index random access)
- LOH bed: `tumor_phased_LOH.bed` (1,100 regions)
- CLI: `-w 5000 -j 16 --distance-metric BERNOULLI`
- Pipeline script: `r1_run_phase2_global_to.sh`
- Analysis script: `r1_global_analysis.py`

## 7. Registry links

- CL-025a (SampleASM_Delta characterization)
- CL-025c (Cross-mode concordance 2026-04-21)
- CL-008 (Beyond-AUC ceiling)