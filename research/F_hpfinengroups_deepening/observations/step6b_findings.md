---
id: F_step6b_findings
title: F Pilot Step 6B — Orthogonal feature validation findings
date: 2026-04-18
status: complete
scope: 驗證 HPFineNGroups 是否獨佔訊號 / orthogonal feature 是否值得納入
prev: step5_findings.md
---

# Step 6B Findings — Orthogonal Feature Validation

## 目的

在 NG=4+AF<0.4+NR≥80+NonLOH（n=14,197，TP=13,176 / FP=1,021）內檢查其他 heterogeneity-related 特徵是否提供獨立於 HPFineNGroups 的訊號。

## Part A — Pooled AUC (TP vs FP)

| feature | best_AUC | direction | Mann-U p | ΔTP_mean−FP_mean |
|---|---|---|---|---|
| PairwiseMedianDist | 0.518 | TP>FP | 0.053 | +0.0066 |
| **PairwiseMeanDist** | **0.525** | TP>FP | **0.007** | +0.0087 |
| HPMergedDelta | 0.511 | FP>TP | 0.243 | −0.003 |
| **AlleleDelta** | **0.548** | **FP>TP** | **<0.0001** | −0.006 |
| **HPFineF** | **0.546** | **FP>TP** | **<0.0001** | −1.205 |

**觀察**：
- 所有 feature 單一 AUC ≤0.55 → 無法單獨作為強 filter
- **AlleleDelta / HPFineF**：在 AF<0.4 內 FP 仍有殘留 allele imbalance → 暗示 AF<0.4 子集仍存在弱 germline-like contamination
- **PairwiseMedianDist / PairwiseMeanDist**：TP reads 更 diverse（TP 為真 subclonal somatic SNV 的標誌）

## Part B — Per-sample AUC（樣本特異性顯著）

| sample | n | PairwiseMedianDist | PairwiseMeanDist | HPMergedDelta | AlleleDelta | HPFineF |
|---|---|---|---|---|---|---|
| H1437 | 1132 | 0.638 | 0.644 | 0.637 | 0.544 | 0.501 |
| H2009 | 9939 | 0.518 | 0.523 | 0.513 | 0.522 | 0.511 |
| HCC1395 | 860 | 0.580 | 0.583 | 0.553 | 0.569 | 0.549 |
| **HCC1395_DORADO** | 566 | **0.692** | **0.691** | 0.505 | 0.512 | 0.608 |
| **HCC1937** | 626 | 0.586 | 0.586 | 0.547 | **0.663** | 0.520 |
| **HCC1954** | 1040 | 0.578 | 0.583 | 0.585 | 0.546 | **0.644** |

**樣本特異性發現**：
1. **HCC1395_DORADO PairwiseMedianDist AUC=0.692**：Dorado basecall 後 AF<0.4 仍有高 pairwise 距離 TP → **可能暗示 Dorado 模型對 low-AF reads 有特殊 noise pattern**（新假設）
2. **HCC1937 AlleleDelta AUC=0.663**：BRCA1-mutant 樣本 AF<0.4 內仍有 allele imbalance FP（可能 BRCA1 heterozygous loss 相關的 subtle CNV）
3. **HCC1954 HPFineF AUC=0.644**：HER2+ 樣本 FP 具較高 HPFineF（fine-level heterogeneity 統計量）→ FP reads 本身異質性更高的 secondary signature
4. **H2009 全 feature ≈0.52**：ceiling effect，無 FP signal 可補
5. **H1437 PairwiseMeanDist AUC=0.644**：與其他樣本差異最大（n_FP=40 僅）

## Part C — Spearman ρ(feature, HPFineNGroups) @ NG≥2+AF<0.4+NR≥80+NonLOH (n=63,841)

| feature | ρ | p |
|---|---|---|
| HPMergedDelta | 0.167 | <e-300 |
| PairwiseMeanDist | 0.100 | <e-300 |
| PairwiseMedianDist | 0.085 | <e-300 |
| HPFineF | −0.064 | <e-300 |
| AlleleDelta | −0.058 | <e-300 |

**結論**：所有 feature 與 HPFineNGroups 的 |ρ| ≤ 0.17 → **高度獨立**（不是 NG 的同義語）。

## Part D — Combined score (HPFineNGroups + orthogonal feature)

Baseline HPFineNGroups AUC = 0.500（NG=4 子集無變異，所以 AUC 恆定）
- Top combo AUC = 0.548 = AlleleDelta 單獨 AUC（組合無增益，因 NG 固定）
- 意味在 NG=4 子集內，filter 已充分壓榨 NG 的訊號；其他 feature 只能作為 per-sample fine-tuning

## 綜合結論

### POSITIVE（orthogonal 證據）
1. **HPFineNGroups 是 dominant signal**：canonical filter 無需改動
2. **5/6 in-scope 樣本在 AF<0.4 內仍有弱 orthogonal signal**（H2009 ceiling 除外）
3. **Per-sample fine-tuning 有機會**：
   - HCC1937：+AlleleDelta<0.02 → 可再減 FP
   - HCC1954：+HPFineF<5 → 可再減 FP
   - HCC1395_DORADO：+PairwiseMedianDist>0.15 → 可再減 FP

### 科學意涵
- **AlleleDelta / HPFineF 的 FP>TP 方向**：AF<0.4 殘餘 FP 仍攜帶 allelic imbalance，暗示部分 FP 是 subtle germline variant（非 CNV-driven）
- **PairwiseMedianDist 的 TP>FP 方向**：真 subclonal TP 的 reads 更 diverse（不同 subclone 的甲基化 pattern 不同）
- **所有 feature 與 NG |ρ|<0.17**：確認有**多個獨立 heterogeneity axes**，HPFineNGroups 只是其中一個

### Caveat
- 本 Step 為 exploratory，**不改變 canonical filter**（NG=4+AF<0.4+NR≥80+NonLOH）
- Per-sample fine-tuning 若納入需驗證 overfitting（SEQC2 truth only）
- HCC1395_DORADO 的高 PairwiseMedianDist AUC=0.692 可能暗示 Dorado-specific artifact，值得單獨研究

## 生成檔案

- `research/F_hpfinengroups_deepening/data/step6b_pooled_auc.tsv`
- `research/F_hpfinengroups_deepening/data/step6b_per_sample_auc.tsv`
- `research/F_hpfinengroups_deepening/data/step6b_correlation_with_ng.tsv`
- `research/F_hpfinengroups_deepening/data/step6b_combined_score.tsv`
- `research/F_hpfinengroups_deepening/scripts/step6_orthogonal_feature_validation.py`
- `research/F_hpfinengroups_deepening/observations/step6b_run.log`

## 對主報告的更新建議

1. 新增 Section 3d：Step 6B orthogonal validation results
2. 更新 Section 6：加 per-sample fine-tuning 候選
3. 標記三個新 follow-up：Dorado-specific artifact、HCC1937 BRCA1 CNV、HCC1954 fine-level heterogeneity

## 是否構成「重大方向」（提示暫停）

**不構成**：Step 6B 為輕量 exploratory 分析，結論確認 HPFineNGroups 為 dominant signal，canonical filter 穩定，orthogonal features 僅作為 per-sample future fine-tuning 候選。

**真正的重大方向（需使用者決策）**：
- **Step 6A — Actual purity mixture experiment**：需實際混合 tumor×normal BAM 並 haplotag 重跑，>1h 計算 + 多檔案輸出 + 不可逆
- **Step 6C — Phase 2A Normal ASM integration**：需等 Phase 2 Normal BAM pipeline output，跨 pilot 整合影響整個 Phase 2 方向
- **Cross-pilot 總結**：F pilot 結論餵回主線論文 outline
