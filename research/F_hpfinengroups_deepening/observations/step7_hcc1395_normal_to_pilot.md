# Step 7 — HCC1395 Normal BAM pilot, TO arm

**Date**: 2026-04-19
**Scope**: F pilot filter subset (NG=4+AF<0.4+NR≥80 NonLOH) → HCC1395 TO-mode subset VCF → inter_sub_mod Phase 2 pipeline → Δ_ASM 分析
**Pipeline**: `r1_run_phase2_both_arms.sh` (parallel TO + paired)；binary: `build/bin/inter_sub_mod`；Normal BAM: HCC1395BL region-subset 1.8GB

## 1. 樣本規模

- Input subset VCF variants: 801
- Regions with truth_label merged: 801
- TP=720, FP=81, TP_rate=0.8989

## 2. 全子集 per-feature AUC（TP vs FP）

| Feature | AUC | n |
|---------|-----|---|
| NormalBaseline_Coverage | 0.6867 | 801 |
| SampleASM_Delta | 0.6434 | 801 |
| HPFineF | 0.5820 | 801 |
| Normal_HP_Delta | 0.5783 | 798 |
| HP_Residual_Delta | 0.5768 | 801 |
| Fisher_Frac_Sig | 0.5406 | 801 |
| Tumor_HP_Delta | 0.5237 | 794 |
| Entropy_Imbalance | 0.5220 | 800 |
| HPFineNGroups | 0.5197 | 801 |
| Combined_HP_Signed_Delta | 0.5183 | 801 |
| NormalBaseline_Mean | 0.5175 | 801 |
| PairwiseMedianDist | 0.5167 | 801 |
| HP_Signed_Residual | 0.5150 | 798 |
| Epipoly_Delta | 0.5071 | 800 |
| CramersV | 0.5036 | 801 |

## 3. F pilot 子集（NR≥80 + NG≥4）per-feature AUC

| Feature | AUC | n |
|---------|-----|---|
| NormalBaseline_Coverage | 0.6841 | 743 |
| SampleASM_Delta | 0.6446 | 743 |
| Normal_HP_Delta | 0.5870 | 742 |
| HP_Residual_Delta | 0.5785 | 743 |
| HPFineF | 0.5746 | 743 |
| Fisher_Frac_Sig | 0.5398 | 743 |
| PairwiseMedianDist | 0.5332 | 743 |
| Epipoly_Delta | 0.5302 | 743 |
| NormalBaseline_Mean | 0.5263 | 743 |
| Entropy_Imbalance | 0.5214 | 743 |
| Combined_HP_Signed_Delta | 0.5178 | 743 |
| Tumor_HP_Delta | 0.5134 | 743 |
| HP_Signed_Residual | 0.5121 | 743 |
| CramersV | 0.5042 | 743 |
| HPFineNGroups | 0.5000 | 743 |

## 4. Sample ASM / Residual 顯著性

- `SampleASM_Sig=true`：795 regions；TP rate=0.9006（vs overall 0.8989）
- `HP_Residual_Sig=true`：581 regions；TP rate=0.8881

## 5. 小結（2026-04-19 HCC1395 單樣本結果）

**主要發現**：

1. **`NormalBaseline_Coverage` AUC=0.687（全子集）/ 0.684（F pilot subset）= 最強訊號** — Normal BAM 在 tumor region 的覆蓋度差異，是 region-level TP/FP 最強的新特徵；**源自 Phase 2 Normal BAM 整合**，原 TO-only pipeline 無此維度。
2. **`SampleASM_Delta` AUC=0.643** — Sample ASM (tumor ASM − normal ASM) 在 F pilot 已預篩子集內仍展現 +0.143 over 隨機的區分力。
3. **`HP_Residual_Delta` AUC=0.577 / `Normal_HP_Delta` 0.578** — 均顯示 normal BAM 資訊對 TP/FP 分離有獨立貢獻。
4. **傳統 ISM 特徵在此子集飽和**（HPFineNGroups AUC=0.500，CramersV=0.504）— F pilot filter 已 pre-select TP-rich subset，傳統特徵無殘餘區分力，符合 Beyond-AUC 耗盡（CL-008）預測。

**保守解讀（single sample + 81 FP）**：

- TP_rate=0.899（F pilot filter 已 pre-enrich）；FP 數 81 較小，AUC 需謹慎。
- `NormalBaseline_Coverage` 可能包含 coverage drift confound（non-random tumor/normal sampling）；需 residualize on NR + AF + Coverage_Multiple 做正式驗證。
- `SampleASM_Delta` 與 `NormalBaseline_Coverage` 相關性未檢查；若高相關則 information gain 重疊。

**結論落點（初步）**：TO arm Normal BAM **具 filter 潛力**（AUC ≥0.60 達成）；但需 R12 batch 3 跨樣本驗證（至少 3/6 樣本保持 AUC ≥0.60）才能提升 stability rating。

**下一步**：
1. Residualize `NormalBaseline_Coverage` on NR + Coverage_Multiple 後重算 AUC（排除 coverage drift confound）
2. R12 跨樣本擴展（HCC1954、COLO829 除外的 4 in-scope 樣本）
3. 與 F pilot filter 聯動：確認 `NormalBaseline_Coverage` + F pilot filter 的 ensemble F1 是否 > 單獨 F pilot filter

## 6. Provenance

- Binary: `build/bin/inter_sub_mod`
- Normal BAM: `/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot/HCC1395_normal_subset.bam` (1.8GB, region-subset from 136GB source)
- Tumor BAM (TO): see `r1_run_phase2_both_arms.sh`
- REF: `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`
- CLI: `-w 5000 -j 16 --distance-metric BERNOULLI`
- Output dir: `output/hcc1395_normal_pilot/TO/`
- Subset VCF script: `r1_build_subset_vcf_bcftools.sh`（bcftools 1.13）
- Analysis script: `r1_delta_asm_analysis.py`

## 7. Registry 連結

- 對應 CL-013, CL-014（Phase 2 code + B/C/D validation）；R1 batch 2 執行
- 5-axis matrix: `docs/reports/research_landscape/data/10_Registry_5axis_matrix.tsv`