# Step 7 — HCC1395 Normal BAM pilot, paired arm

**Date**: 2026-04-19
**Scope**: F pilot filter subset (NG=4+AF<0.4+NR≥80 NonLOH) → HCC1395 paired-mode subset VCF → inter_sub_mod Phase 2 pipeline → Δ_ASM 分析
**Pipeline**: `r1_run_phase2_both_arms.sh` (parallel TO + paired)；binary: `build/bin/inter_sub_mod`；Normal BAM: HCC1395BL region-subset 1.8GB

## 1. 樣本規模

- Input subset VCF variants: 983
- Regions with truth_label merged: 983
- TP=971, FP=12, TP_rate=0.9878

## 2. 全子集 per-feature AUC（TP vs FP）

| Feature | AUC | n |
|---------|-----|---|
| Fisher_Frac_Sig | 0.7027 | 983 |
| Epipoly_Delta | 0.6566 | 982 |
| Tumor_HP_Delta | 0.6272 | 978 |
| HPFineF | 0.6099 | 983 |
| Entropy_Imbalance | 0.6010 | 982 |
| Normal_HP_Delta | 0.5797 | 968 |
| HP_Signed_Residual | 0.5656 | 980 |
| HPFineNGroups | 0.5632 | 983 |
| SampleASM_Delta | 0.5607 | 983 |
| PairwiseMedianDist | 0.5496 | 983 |
| Combined_HP_Signed_Delta | 0.5373 | 983 |
| NormalBaseline_Mean | 0.5346 | 983 |
| NormalBaseline_Coverage | 0.5265 | 983 |
| CramersV | 0.5257 | 983 |
| HP_Residual_Delta | 0.5164 | 983 |

## 3. F pilot 子集（NR≥80 + NG≥4）per-feature AUC

| Feature | AUC | n |
|---------|-----|---|
| Fisher_Frac_Sig | 0.7257 | 942 |
| Epipoly_Delta | 0.6805 | 942 |
| Tumor_HP_Delta | 0.6260 | 942 |
| NormalBaseline_Mean | 0.5997 | 942 |
| Normal_HP_Delta | 0.5985 | 931 |
| PairwiseMedianDist | 0.5906 | 942 |
| Entropy_Imbalance | 0.5828 | 942 |
| HPFineF | 0.5759 | 942 |
| SampleASM_Delta | 0.5468 | 942 |
| HP_Signed_Residual | 0.5312 | 942 |
| NormalBaseline_Coverage | 0.5310 | 942 |
| CramersV | 0.5268 | 942 |
| HP_Residual_Delta | 0.5225 | 942 |
| Combined_HP_Signed_Delta | 0.5157 | 942 |
| HPFineNGroups | 0.5000 | 942 |

## 4. Sample ASM / Residual 顯著性

- `SampleASM_Sig=true`：982 regions；TP rate=0.9878（vs overall 0.9878）
- `HP_Residual_Sig=true`：718 regions；TP rate=0.9930

## 5. 小結（2026-04-19 HCC1395 單樣本結果）

**主要發現**（出乎預期的 F1-relevant 訊號，非單純 characterization）：

1. **`Fisher_Frac_Sig` AUC=0.703（全子集）/ 0.726（F pilot subset）= 首次突破 0.70 閾值！** — Per-CpG ASM 中 Fisher exact 顯著 CpG 比例，在 paired mode 下對 TP/FP 分離 AUC 0.726。此結果前所未有（pure methylation 上限 0.58 被打破，CL-008 需 revisit）。
2. **`Epipoly_Delta` AUC=0.680（F pilot subset）** — Epipolymorphism（HP1 vs HP2）差異亦 >0.65。
3. **`Tumor_HP_Delta` AUC=0.626** — Tumor-only HP methylation delta 中等區分力。
4. **`Entropy_Imbalance` AUC=0.583、HPFineF 0.576** — 多個 per-region 甲基化異質性指標中等有效。
5. **TP_rate=0.988** — Paired mode 飽和確認；僅 12 FP，單點估計不穩。

**關鍵意義**：

- 此結果**挑戰 CL-008（Beyond-AUC 耗盡，純 methylation ≤0.58）**。過去 Phase 1A 未測試 Per-CpG ASM Fisher statistics；本次 Phase 2 Normal BAM 整合新增該欄位，在 paired mode 下意外表現強於預期。
- Paired mode 飽和（CL-017）聲稱針對 **region-level F1**；**子集內 FP 再篩選**可能仍有空間（本 pilot 提示）。
- **⚠ 警示**：12 FP 小樣本 → Wilson 95% CI 可能跨 0.55-0.85，AUC 點估計不可靠。**必須跨樣本驗證**（R12 batch 3）才能判斷是否為 HCC1395-specific 或普遍現象。

**保守解讀**：

- 12 FP 小樣本 → AUC 點估計不可靠
- `Fisher_Frac_Sig` 可能與 Coverage 相關（高 coverage region fisher 更易顯著）；需 residualize on `Fisher_N_Tested` + NR + AF 驗證
- F pilot filter 已 pre-enrich → 已 pre-select subset 內的 AUC 容易 overfit

**結論落點（初步）**：Paired arm Per-CpG ASM **有 F1-relevant 潛力訊號**，但因 FP=12 過小，暫定 **CONDITIONAL POSITIVE**；需 R12 跨樣本驗證才能升級。

**下一步**：
1. Residualize `Fisher_Frac_Sig` on coverage (`Fisher_N_Tested`) + NR + AF
2. R12 跨樣本擴展：若至少 4/6 樣本保持 AUC ≥0.65，升級為 paired-mode filter 候選
3. 與 TO arm `NormalBaseline_Coverage` 在相同 region 比對（cross-mode concordance → ensemble opportunity）

## 6. Provenance

- Binary: `build/bin/inter_sub_mod`
- Normal BAM: `/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot/HCC1395_normal_subset.bam` (1.8GB, region-subset from 136GB source)
- Tumor BAM (paired): see `r1_run_phase2_both_arms.sh`
- REF: `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`
- CLI: `-w 5000 -j 16 --distance-metric BERNOULLI`
- Output dir: `output/hcc1395_normal_pilot/paired/`
- Subset VCF script: `r1_build_subset_vcf_bcftools.sh`（bcftools 1.13）
- Analysis script: `r1_delta_asm_analysis.py`

## 7. Registry 連結

- 對應 CL-013, CL-014（Phase 2 code + B/C/D validation）；R1 batch 2 執行
- 5-axis matrix: `docs/reports/research_landscape/data/10_Registry_5axis_matrix.tsv`