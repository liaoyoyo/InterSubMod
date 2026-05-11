---
id: ism-kb-09-conclusions-negative-17-to-feature-deep-study
name: "N17: TO Feature Deep Study Q1-Q6 — RF ceiling 0.69-0.77"
description: "2026-04-05 6 問 32 特徵 5 分層；HP2FamilyN AUC 0.72 去混淆後→0.54（73% artifact）；RF CV ceiling LOH 0.71-0.74/LOH_Weak 0.75-0.77；甲基化距離全 AUC≈0.50；AlleleSig Non-LOH 反轉；Baseline > v2b in LOH。"
status: archived
last_verified: 2026-04-05
content_nature: historical-note
doc_type: explanation
verified_scope: "ADR structured from MEMORY project_to_feature_deep_study.md"
decision_refs: []
related_ids:
  - ism-kb-09-conclusions-concluded-negative
tags: [negative, no-go, to-feature, hp2familyn, rf-ceiling, circular-artifact]
canonical_paths: [09_conclusions/negative/N17_to_feature_deep_study.md]
alias_paths: []
---

# N17: TO Feature Discrimination Deep Study (Q1-Q6)

## Status
`concluded-negative`（HP2FamilyN 硬閾值 filter STOP；TO ISM 重新定位為輔助特徵）。

## Context
矯正 ClairS-TO VCF 後 TO 特徵區分力需重新量化。Q1-Q6 六問深度比較 longphase-to baseline vs PON-only v2b 兩版 haplotagging。

## Decision-tested
1. HP2FamilyN 等高 AUC 特徵是否循環依賴 artifact
2. RF multi-feature ensemble 的 CV AUC ceiling
3. Baseline vs v2b haplotagging 哪版 discrimination 更好
4. 甲基化距離特徵是否有區分力

## Method
- 32 特徵 × 5 分層 × 2 haplotag 版本
- HP2FamilyN 去混淆（剔除循環 features）
- 16 特徵 RF 5-fold CV

## Result
1. **HP2FamilyN AUC 0.72 → 去混淆後 0.54**（73% 是循環依賴 artifact）
2. **RF 16 特徵 CV AUC ceiling**：All 0.69, LOH 0.71-0.74, LOH_Weak 0.75-0.77, Non-LOH 0.64
3. **Baseline > v2b in LOH**：v2b 擴展 LOH coverage (+14pp) 但稀釋 discrimination
4. **甲基化距離特徵全 AUC ≈ 0.50**（PairwiseMeanDist, CramersV）
5. **AlleleSig Non-LOH 反轉**（FP>TP，germline ASM）
6. **AlleleDelta 僅在 minor HP=0-2 有效** (AUC 0.58-0.62)，其餘 <0.50

### Baseline vs v2b 關鍵對比
| 指標 | v2b | Baseline |
|------|-----|----------|
| LOH% | 72.6 | 58.7 |
| RF LOH AUC | 0.712 | **0.736** |
| HPFineP AUC (LOH) | 0.615 | **0.661** |
| HPFineSig Δ (LOH_Strong) | +0.077 | **+0.177** |

## Consequences
- **根因**：HP2FamilyN 本身是 HP tag count 的衍生，TP 本身含 somatic HP 較多，造成循環；v2b 改善 haplotagging 但擴入邊界區域，discrimination 反而稀釋。
- **替代方向**：
  - HP2FamilyN 硬閾值 filter：**STOP**（循環 + 方向錯誤）
  - TO ISM 定位：輔助特徵（RF AUC 0.69）而非獨立 filter
  - 研究方向：FN rescue（後 N09 否決）/ ISM 作 caller 輔助特徵 / 低純度樣本驗證
- **決策**：2026-04-05 TO ISM 重新定位為輔助特徵。

## References
- MEMORY：`project_to_feature_deep_study.md`
- Report：`20260405_TO_ClairSTO特徵區分力深度研究_01.md`
- Figures：`output/clairsto_correction_analysis/feature_study/Q1-Q6/`
