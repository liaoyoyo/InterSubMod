---
title: 全 7 樣本 subcluster 切割總帳（canonical paired_full longphase-s, tumor-only）
date: 2026-06-20
status: in_progress
task_type: B_validation
scope: 全基因組 + 全 7 樣本（comprehensive）
build_branch: feat/summary-nreadsvalid @ 2d5692f
binary: .claude/worktrees/ism-review-infra/build/bin/inter_sub_mod
data_sources: results/result_*.json, results/significance_summary_*.csv, results/summary_all_samples.json
---

# 全 7 樣本 subcluster 切割總帳

## 0. 一句話結論

把 HCC1395 的「切割總帳」(A=Fisher+V cansplit vs C=PERMANOVA clean-location 的 Venn) **逐字相同流程**跑遍全 7 樣本 canonical paired_full longphase-s tagged BAM（tumor-only）。
**A 與 C 的 Jaccard 在 7 樣本全部偏低（0.091–0.161）** — 無監督切群 (cansplit) 與標籤對齊 PERMANOVA (clean-location) 在每個樣本、跨 3 癌種都擷取到**大致不相交**的位點集合，穩健重現 HCC1395 的單樣本發現。

## 1. Provenance（可重現）

| 項目 | 值 |
|------|----|
| Binary | `inter_sub_mod` @ `feat/summary-nreadsvalid` commit 2d5692f（含 reclassify-v2 / within-HP B,D,C′ / strength_tumor）|
| 輸入 BAM | canonical **paired_full longphase-s tagged**（`{sample}_tagged.bam`，HP:Z:X-Y 格式，含 MM/ML 甲基化 tag）|
| 模式 | **tumor-only**（無 normal BAM）|
| ISM 參數 | `-w 5000 --distance-metric BERNOULLI --min-common-coverage 3 --nan-distance-strategy SKIP --methyl-low 0.2 --methyl-high 0.8 -j 48` |
| Reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta` |
| 輸入 VCF | 各樣本 canonical paired_full `filtered_snv_tp.vcf.gz`（ClairS somatic TP）|

**流程鏈**（單次 ISM + 逐字後處理）：
1. ISM scan（保留 per-region distance matrix）→ `significance_summary.csv`（117 欄）
2. `compute_split_accounting.py`：逐字複製 `wg_contingency.py` 的 `peel`/`cluster_bestk`（UPGMA + silhouette，tumor reads × LABMAP）→ 每 region best_k → **A=cansplit (best_k≥2)**
3. 讀 significance_summary 的 `LabelHP/AllelePermanovaP<0.05 & !DispersionWarn` → **C=clean PERMANOVA**
4. Q1 Venn（逐字 `permanova_clean_and_4group.py`）

> **優化說明**：原始 HCC1395 流程跑 2 次 ISM（sig scan + wg_contingency per-chr）。因 per-region 距離矩陣與 binary/參數/region 一一對應為位元相同，本流程改**單次 ISM scan 同時取得 significance + 距離矩陣**，計算減半、數字一致。

## 2. 邏輯驗證（信任基礎）

用 HCC1395 **pileup-track** 既有 `records_wg2.json` + `significance_summary.csv` 重跑本 driver 的 Q1+分帶邏輯，對比已 commit 的參考數字：

| 指標 | driver | 參考 | |
|------|--------|------|---|
| A=cansplit | 5997 (19.7%) | 5997 (19.7%) | ✅ 位元一致 |
| C=clean PERMANOVA | 8833 (29.0%) | 8833 (29.0%) | ✅ |
| A∩C / A−C / C−A | 1560 / 4437 / 7273 | 1560 / 4437 / 7273 | ✅ |
| Jaccard | 0.118 | 0.118 | ✅ |
| 分帶 insuff/cant/cansplit/degen/vweak | 412/24081/5997/536/1249 | 全中 | ✅ |
| clear/mid/weak | 2942/767/503 | 2941/765/506 | ±1-3（總和 4212 同；分帶腳本未 commit，屬定義重建容差）|

## 3. 主表 — A/C clean-location Venn（7 樣本，全基因組）

| 樣本 | N | A=cansplit | C=clean PERMANOVA | A∩C | A−C | C−A | Jaccard |
|------|---|-----------|-------------------|-----|-----|-----|---------|
| HCC1395 | 29754 | 5816 (19.5%) | 8247 (27.7%) | 1359 (4.6%) | 4457 (15.0%) | 6888 (23.1%) | **0.107** |
| HCC1395_DORADO | 29889 | 7111 (23.8%) | 6662 (22.3%) | 1419 (4.7%) | 5692 (19.0%) | 5243 (17.5%) | **0.115** |
| COLO829 | 35185 | 7799 (22.2%) | 5522 (15.7%) | 1487 (4.2%) | 6312 (17.9%) | 4035 (11.5%) | **0.126** |
| H2009 | 132909 | 54939 (41.3%) | 29069 (21.9%) | 11622 (8.7%) | 43317 (32.6%) | 17447 (13.1%) | **0.161** |
| H1437 | 67468 | 16582 (24.6%) | 15181 (22.5%) | 3114 (4.6%) | 13468 (20.0%) | 12067 (17.9%) | **0.109** |
| HCC1937 | 12393 | 1875 (15.1%) | 2800 (22.6%) | 388 (3.1%) | 1487 (12.0%) | 2412 (19.5%) | **0.091** |
| HCC1954 | 17909 | 6493 (36.3%) | 5423 (30.3%) | 1646 (9.2%) | 4847 (27.1%) | 3777 (21.1%) | **0.160** |

## 4. 切割品質總帳（分帶為依 build_spectrum_html.py 定義重建，邊界±容差）

| 樣本 | N | insuff(n<6) | cant(切不出) | cansplit | clear(V≥.7+sig) | mid(.5-.7) | weak(.3-.5) | vweak(<.3) | degen |
|------|---|-------------|--------------|----------|-----------------|------------|-------------|------------|-------|
| HCC1395 | 29754 | 700 (2.4%) | 23238 (78.1%) | 5816 (19.5%) | 2897 (9.7%) | 711 | 465 | 1192 | 551 |
| HCC1395_DORADO | 29889 | 3913 (13.1%) | 18865 (63.1%) | 7111 (23.8%) | 3572 (12.0%) | 754 | 538 | 1249 | 998 |
| COLO829 | 35185 | 2142 (6.1%) | 25244 (71.7%) | 7799 (22.2%) | 859 (2.4%) | 706 | 1645 | 3540 | 1049 |
| H2009 | 132909 | 3080 (2.3%) | 74890 (56.3%) | 54939 (41.3%) | 1466 (1.1%) | 1481 | 3805 | 41594 | 6593 |
| H1437 | 67468 | 21099 (31.3%) | 29787 (44.1%) | 16582 (24.6%) | 3382 (5.0%) | 2421 | 2700 | 4851 | 3228 |
| HCC1937 | 12393 | 120 (1.0%) | 10398 (83.9%) | 1875 (15.1%) | 515 (4.2%) | 173 | 206 | 787 | 194 |
| HCC1954 | 17909 | 451 (2.5%) | 10965 (61.2%) | 6493 (36.3%) | 925 (5.2%) | 863 | 1453 | 3175 | 77 |

## 4b. TP vs FP 比例對照（FP=filtered_snv_fp.vcf.gz，同流程 tumor-only）

| 樣本 | TP N | TP A% | TP C% | FP N | FP A% | FP C% | ΔA%(TP−FP) | ΔC%(TP−FP) |
|------|------|-------|-------|------|-------|-------|------------|------------|
| HCC1395 | 29754 | 19.5 | 27.7 | 627 | 38.8 | 28.1 | −19.3 | −0.4 |
| HCC1395_DORADO | 29889 | 23.8 | 22.3 | 240 | 26.2 | 12.5 | −2.4 | +9.8 |
| COLO829 | 35185 | 22.2 | 15.7 | 2273 | 21.9 | 17.7 | +0.3 | −2.0 |
| H2009 | 132909 | 41.3 | 21.9 | 86 ⚠低N | 27.9 | 9.3 | +13.4 | +12.6 |
| H1437 | 67468 | 24.6 | 22.5 | 8 ⚠低N | 12.5 | 0.0 | +12.1 | +22.5 |
| HCC1937 | 12393 | 15.1 | 22.6 | 195 | 15.9 | 3.6 | −0.8 | +19.0 |
| HCC1954 | 17909 | 36.3 | 30.3 | 29 ⚠低N | 13.8 | 34.5 | +22.5 | −4.2 |

- **Pooled**：TP N=325507 A=30.9% C=22.4%；FP N=3458 A=25.0% C=18.3%（pooled 受最大樣本主導，僅參考）。
- 🔴 **cansplit 不區分 TP/FP**：可靠的大-N FP（HCC1395 FP=627、COLO829 FP=2273、DORADO FP=240、HCC1937 FP=195）的 A% 與 TP **相當或更高**（HCC1395 FP A=38.8% > TP 19.5%）→ 切群能力非 TP/FP filter，重現「ASM/結構非 usable filter」。
- ⚠ H1437(8)/HCC1954(29)/H2009(86) FP N<100，個別比例為噪音不可解讀。

## 5. 見樹也見林（四層觀察）

- **Aggregate**：Jaccard(A,C) 全 7 樣本 0.091–0.161，**一律偏低**。「切得出群」與「標籤對齊乾淨 PERMANOVA」是大致不相交的兩件事 — 跨 3 癌種穩健。
- **Canonical**：HCC1395 paired_full（A 19.5% / C 27.7% / J 0.107）≈ pileup-track 參考（19.7% / 29.0% / 0.118）。**換 BAM track 幾乎不動主數字**，互證穩定。
- **Outlier**：
  - **H2009**：cansplit 高達 41.3%，但其中 **vweak(V<0.3) 41594 = 76%** → 切得多但多數無意義（高覆蓋造成過度切分）。
  - **H1437**：insuff 31.3%（最高）→ 大量 region tumor reads <6，coverage 偏低。
  - **COLO829**：cansplit 22.2% 但 clear 僅 2.4%（最低）→ 切得出但對齊獨立軸的極少。
- **Well-explained**：低 Jaccard 重現 memory `project_apriori_subclone_classification_model`「within_clean≠subclone, Jaccard 0.123」與「切不出≠沒訊號」框架 — **無監督切群數 ≠ 可確認 subclone**。

## 6. Caveats（口徑）

1. **tumor-only**：normal BAM 各樣本來源不對稱（COLO829 normal=PAO33946 等），且 cansplit 分群本就只用 tumor reads（`is_tumor` filter），故統一 tumor-only 使 7 樣本完全可比。對 A 無影響；對 C（LabelHP/Allele PERMANOVA）在 tumor-only 下可算且自洽。
2. **track 差異**：你原始參考表（5997/8833/0.118）為 **pileup-track**；本表 HCC1395 為 **paired_full longphase-s**，故略有差異（見 §5 Canonical，差異極小）。
3. **分帶重建**：clear/mid/weak 的精確計算腳本未 commit，本表依 build_spectrum_html.py 顯示定義（V≥0.7+sig 等）重建，邊界 ±1-3 容差；insuff/cant/cansplit/degen/vweak 為確定性定義無容差。
4. **unmatched significance**：records_wg2 與 significance_summary 的 Pos 鍵少量不匹配（0–216，最高 COLO829 0.6%），對 C 影響可忽略。
5. **單 binary / 單 pipeline**：tier 上限 ⭐3（cross-platform 復現需 SEQC2/外部工具，見既有框架）。

## 7. 檔案

- `results/SUMMARY_all_samples.md` — 自動彙整對照表
- `results/result_<S>.json` — 每樣本切割總帳真值
- `results/significance_summary_<S>.csv` — 每樣本完整 117 欄 ISM 輸出
- `results/records_wg2_<S>.json` — 每樣本 per-region best_k + contingency
- `scripts`：`InterSubMod/scripts/analysis/subcluster_split_allsamples/{compute_split_accounting,validate_against_reference,aggregate_results}.py + run_all_samples.sh`
