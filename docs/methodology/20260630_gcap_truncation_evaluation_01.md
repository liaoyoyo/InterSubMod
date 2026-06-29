---
title: GCAP=8 genotype 截斷評估（Part C）
date: 2026-06-30
status: evaluation
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json, /big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone/{sample}/topology_per_region.json
---

# GCAP=8 genotype 截斷評估

## TL;DR
genotype 向量截斷在 8 位點（`GCAP=8`，套用於上游 `sm_region_integration` 的 populations 建構）。**評估結論：不提高 GCAP**，改以前端誠實標 `truncated` 旗標。理由：影響高度偏斜（5/7 樣本 <1.1%；唯 H2009 17.7% 但仍在生成中），而提高需重跑整條管線（re-genotyping → topology → backbone → scoring → enumerate → build × 7 樣本），會作廢凍結 HCC1395 參考 + 剛驗證的 R6/enumerate 全部工作；成本遠大於收益。

## 量化（detail.truncated 旗標，全 7 樣本實測）

| 樣本 | 總區 | 截斷區 | 占比 | max n_sSNV | 截斷區 determinacy 主成分 |
|---|---|---|---|---|---|
| HCC1395 | 3885 | 42 | 1.1% | 150 | B_pairwise 25 / incompatible 8 |
| COLO829 | 3786 | 5 | 0.1% | 26 | B_pairwise 4 |
| H1437 | 4768 | 220 | 4.6% | 1007 | B_pairwise 105 / incompatible 50 / A_determined 58 |
| **H2009** | 4243 | **752** | **17.7%** | 1119 | A_determined 266 / **incompatible 250** / B_pairwise 190 |
| HCC1395_DORADO | 2379 | 9 | 0.4% | 63 | incompatible 4 / B_pairwise 4 |
| HCC1937 | 1636 | 2 | 0.1% | 13 | — |
| HCC1954 | 1979 | 5 | 0.3% | 190 | B_pairwise 3 |
| **合計** | 22676 | **1035** | **4.6%** | — | — |

🔴 **截斷區的 incompatible 可能為假象**：四配子違反對若落在第 8 位之後，stored populations 查不到真違反 → 被分類 incompatible 實為截斷 artifact（`backbone_resolution.py` 的 `cap_caveat` 已標）。H2009 250 個截斷 incompatible 最須警惕。

## 截斷位置與提高成本

- **位置**：`sm_region_integration.py` 建 populations 時每區最多 8 sSNV 位點（`topology_analysis.py` provenance `genotype_cap=8`）。populations 依賴 per-read genotyping（「groups」，源自 BAM pileup，見 `sm_linkage_genomewide.py`）。
- **2^k 不爆炸**：perfect-phylogeny 的 2^k 是理論空間,但 populations 是**觀察到的 read 基因型**（數量受 read 數限制，非 2^k）→ solve_topology 多算幾位 = compute 便宜。
- **真成本 = 管線重 derivation**：提高 GCAP 須 re-genotype（BAM-based）+ 重跑 topology/backbone/scoring/enumerate/build × 7 樣本 → 會改動所有 topology_per_region.json，作廢凍結 HCC1395 參考（R1/R6/enumerate 全部 MATCH 基準）+ 已驗證成果。

## 決定與處置
1. **不提高 GCAP**（成本 >> 收益；5/7 樣本影響 <1.1%；H2009 17.7% 但為並行 session 仍在生成的不完整樣本）。
2. **誠實標 truncated 旗標**（§13 路徑）：工作站 detail kv 加紅 badge「⚠ 截斷 n_sSNV>8(偵測不完整)」；`enumerate_candidate_trees.json` 每區帶 `truncated` 欄;此區 ambig/四配子/機率偵測標記不完整。
3. **Reopen 條件**：H2009 定稿後，若「截斷假成環」影響特定下游分析 → 僅對受影響該批區 re-genotype（cap=16）+ 局部重跑，不全管線重derivation。

## 驗證
- 量化數字皆 grep 自各樣本 `topology_per_region.json` 的 `detail[].truncated` 旗標（可重算）。
- 前端 badge：H2009 chr1:194143326-196698970（402 sSNV）driver 驗證 `has_badge=True`、0 pageerror。
