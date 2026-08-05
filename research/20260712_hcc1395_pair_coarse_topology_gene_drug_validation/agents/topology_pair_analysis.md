<!--
建立時間: 2026-07-12T00:39:54+08:00
目標: 計算七個 dataset rows 的五類粗拓撲，並驗證 HCC1395 兩列資料的 region interval concordance
處理範圍: chr1-22; historical layered-v2 engineering snapshot; HCC1395 vs HCC1395_DORADO
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json
  - /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395/layered_reconstruction_HCC1395.json
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395_DORADO/layered_reconstruction_HCC1395_DORADO.json
-->

# HCC1395 兩列資料的區域粗拓撲與區間一致性驗證

> **TL;DR（Task B，服務 G4/G5）**：Exact-coordinate 可配對 6,252 個區域，其中 5,720 個兩邊皆 complete；五類 raw agreement=69.39%（95% CI 68.18%–70.57%）、κ=0.497。這是高於 chromosome-preserving null 的可重現訊號，但 aggregate 組成差異與 exact tree-set digest 一致率顯示兩列資料並非可互換；支持「區域級方向性 reproducibility」，不足以證明方法已生物學驗證。

## 1. 任務分類、研究問題與邊界

- Task type：**B Comprehensive validation**（全 7 dataset rows；不得以 HCC pair 子集代替全樣本 census）。
- 研究問題：同一 biological cell line 的 HCC1395 5kHz 與 DORADO，是否在相同／相近 genomic interval 得到一致的五類 mutation-state 粗拓撲？
- 成功條件：分類守恆、配對鍵唯一、exact 與 interval sensitivity 方向一致、觀察 agreement 顯著高於染色體內 permutation null。
- 失敗條件：五類無法回加 complete、配對多對一、結論只由 dominant Topo>1 驅動，或 exact tree signatures 幾乎不重現。
- Claim ceiling：tree node 是 mutation state；`H_*` 是 Steiner／partial-supported completion state，**不是已確認 clone**。本結果不驗證細胞族群真值，也不建立癌症基因或藥物因果。

## 2. 可稽核 operational classifier

Primary grain 是 **complete primary region 的 ordered HP forest**。每個 HP component 保留 HP1/HP2 順序；region 五類用 component feature 的 OR 合併：

1. `Topo>1 未定`：ordered regional topology 超過一種；即使候選都屬相同粗類，仍不硬選。
2. 其餘 Topo=1：在每個 HP rooted tree 中，`max_outdegree≥2` 表示姐妹分枝；`max_depth≥2` 表示直系鏈。
3. 兩者皆否＝無 HP 內關係；只有分枝＝姐妹 only；只有鏈＝直系 only；兩者皆是＝姐妹＋直系。
4. HP1/HP2 是兩個 forest components，跨 HP 節點**永遠不判姐妹或直系**。
5. 此 primary graph classifier 包含 `H_*`；另外的 observed-state sensitivity 只可把非 ROOT、非 `H_*` endpoint 當直接觀測。本分析不把 hidden node 改名成 clone。

### Golden logic cases

| 類別 | 邏輯 | HCC1395 實際 region | ordered HP coarse | rooted shape set |
|---|---|---|---|---|
| 無 HP 內關係 | ROOT→A（或 HP1/HP2 各自單點） | `chr1:47652992-47686811` | `HP1={no_within_hp_relation}|HP2={no_within_hp_relation}` | `HP1={TS-aa75362b1c}|HP2={TS-aa75362b1c}` |
| 姐妹 only | ROOT→A 且 ROOT→B | `chr1:16637933-16647080` | `HP2={sister_only}` | `HP2={TS-ca3708ca57}` |
| 直系 only | ROOT→A→B | `chr1:94893-113820` | `HP2={direct_only}` | `HP2={TS-a99c6ae6b3}` |
| 姐妹＋直系 | ROOT→A，A→B 且 A→C | `chr1:16583807-16585638` | `HP2={sister_and_direct}` | `HP2={TS-f3932cb25b}` |
| Topo>1 未定 | 候選 rooted shapes >1 | `chr1:1004726-1049980` | `HP1={direct_only,sister_only}` | `HP1={TS-a99c6ae6b3,TS-ca3708ca57}` |

## 3. 全 7 dataset rows 五類組成

五類分母是 complete；Incomplete 另列，並以 `五類=complete`、`complete+incomplete=primary` 雙重守恆。

| Dataset | Primary | Complete | Incomplete | 無 HP 內 | 姐妹 only | 直系 only | 姐妹＋直系 | Topo>1 未定 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 7,590 | 6,940 | 650 | 910 (13.11%) | 334 (4.81%) | 2,053 (29.58%) | 141 (2.03%) | 3,502 (50.46%) |
| HCC1395_DORADO | 7,268 | 6,750 | 518 | 1,337 (19.81%) | 140 (2.07%) | 933 (13.82%) | 34 (0.50%) | 4,306 (63.79%) |
| COLO829 | 7,659 | 6,949 | 710 | 1,259 (18.12%) | 17 (0.24%) | 2,592 (37.30%) | 8 (0.12%) | 3,073 (44.22%) |
| H1437 | 8,630 | 6,984 | 1,646 | 877 (12.56%) | 91 (1.30%) | 3,520 (50.40%) | 96 (1.37%) | 2,400 (34.36%) |
| H2009 | 9,581 | 5,882 | 3,699 | 406 (6.90%) | 112 (1.90%) | 2,745 (46.67%) | 265 (4.51%) | 2,354 (40.02%) |
| HCC1937 | 2,674 | 2,557 | 117 | 344 (13.45%) | 318 (12.44%) | 1,209 (47.28%) | 70 (2.74%) | 616 (24.09%) |
| HCC1954 | 3,975 | 3,823 | 152 | 894 (23.38%) | 362 (9.47%) | 847 (22.16%) | 62 (1.62%) | 1,658 (43.37%) |

HCC1395 兩列的 marginal composition 並不相同：Topo>1 未定為 50.46% vs 63.79%，直系 only 為 29.58% vs 13.82%；五類 total-variation distance=20.03%。所以不能只因為是同一 biological sample 就宣稱分布一致。

## 4. HCC1395 region interval concordance

配對規則：同染色體、一對一；exact 要求 `(chrom,start,end)` 全同。Sensitivity 用 reciprocal overlap 0.80/0.50，或 start/end 雙端 anchor 各在 1 kb/5 kb 內。非 exact 情境用 maximum-cardinality、再 maximum-quality assignment；matched A/B keys 均唯一。

| Scenario | Matched all | Complete both | A/B complete coverage | Raw agreement (95% CI) | κ | Symmetric BA | Macro category Jaccard | Non-dominant-only agree | Remove dominant concordant cell | Permutation null mean (95%) | p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| exact_coordinate | 6,252 | 5,720 | 82.4%/84.7% | 69.39% (68.18%–70.57%) | 0.497 | 0.671 | 0.423 | 87.18% (n=1,677) | 45.50% (n=3,213) | 39.51% (38.51%–40.51%) | 0.00020 |
| reciprocal_overlap_0.80 | 6,369 | 5,793 | 83.5%/85.8% | 69.34% (68.14%–70.52%) | 0.495 | 0.665 | 0.420 | 86.82% (n=1,685) | 45.17% (n=3,239) | 39.72% (38.72%–40.69%) | 0.00020 |
| reciprocal_overlap_0.50 | 6,594 | 5,943 | 85.6%/88.0% | 69.43% (68.24%–70.58%) | 0.492 | 0.661 | 0.416 | 86.59% (n=1,693) | 44.65% (n=3,283) | 40.24% (39.27%–41.22%) | 0.00020 |
| endpoint_anchor_1kb | 6,271 | 5,735 | 82.6%/85.0% | 69.29% (68.09%–70.47%) | 0.496 | 0.666 | 0.421 | 86.98% (n=1,682) | 45.38% (n=3,224) | 39.49% (38.47%–40.47%) | 0.00020 |
| endpoint_anchor_5kb | 6,305 | 5,763 | 83.0%/85.4% | 69.22% (68.01%–70.40%) | 0.495 | 0.664 | 0.419 | 86.83% (n=1,686) | 45.21% (n=3,238) | 39.51% (38.50%–40.53%) | 0.00020 |

### Exact-coordinate 的解讀

- Raw agreement=69.39%、κ=0.497，高於 chromosome-preserving null 39.51%（p=0.00020）：不是只靠染色體別 class composition 的巧合。
- 但 Topo>1 是 pooled dominant class。若只看兩邊都非 dominant 的子集，agreement=87.18%；這個數字排除了 Topo↔determinate 的困難案例，不能單獨當性能。若只移除 `Topo>1/Topo>1` 的同格貢獻，剩餘 agreement 降為 45.50%，說明 headline agreement 明顯受 dominant concordance 驅動。
- Binary resolved/unresolved 2×2（A rows × B columns）：resolved/resolved=1,677、resolved/unresolved=1,096、unresolved/resolved=440、unresolved/unresolved=2,507；agreement=73.15%、κ=0.459。兩邊都 resolved 後 87.18% 是**條件式一致**，不可作全體 accuracy。
- 五類 per-category Jaccard：無 HP 內關係=0.498；姐妹 only=0.361；直系 only=0.400；姐妹＋直系=0.237；Topo>1 未定=0.620。
- Balanced accuracy 以任一 dataset 當 reference 都不等於 truth。A→B=0.561、B→A=0.780，不對稱反映 marginal class shift。
- 染色體分層 exact agreement 範圍為 60.14%（chr8, n=291）到 75.86%（chr21, n=87）；22 條 autosomes 全有 exact matched complete regions。Permutation 已在染色體內打亂，避免把 chromosome composition 當一致性。

### Ordered HP 與 exact tree-set signature

- 兩邊皆 Topo=1 的 exact-coordinate regions：ordered HP coarse signature agreement=59.93%；允許全區 HP1↔HP2 phase flip 後=80.26%。Primary 輸出仍保留 ordered pair，swap 只作跨資料集 sensitivity。
- Rooted unlabeled candidate-shape set：ordered=30.65%、phase-swap tolerant=47.22%。
- 原始 layered exact candidate-tree-set digest：ordered=19.74%、phase-swap tolerant=35.21%（n=5,720）。這是比五類更嚴格的 reproducibility ceiling。

## 5. 科學判定

1. **可支持**：同一 biological sample 的兩列資料，在大量 exact／高 overlap intervals 上有高於 permutation null 的粗拓撲 concordance，方法抓到一部分可重現的區域結構訊號。
2. **不可宣稱完全一致**：marginal 五類 TV distance、Topo>1 比例、balanced accuracy 不對稱與 strict tree-set signature 都顯示 dataset/basecalling/read-support 依賴。
3. **不可由此證明方法有效**：兩列不是獨立 biological replicates，也沒有 clone ground truth；同一樣本 cross-dataset agreement 主要是 reproducibility evidence，不是 accuracy/validity proof。
4. **基因／藥物資料的角色**：癌症基因或藥物 annotation 可檢查已知 biology 的 enrichment/face validity，但不是 topology ground truth；除非有獨立克隆、單細胞、longitudinal 或功能性證據，不能把 annotation overlap 寫成方法被證明。
5. **目前狀態**：`SHARE WITH CAVEATS / SCIENTIFIC NO-GO for proof-of-effectiveness`。clean layered-v3 全 7 樣本尚未 closeout，這裡仍是 historical layered-v2 engineering snapshot。

## 6. Step → Verify 與 QA

1. 讀取 exact catalog、C/T report、regional composition、region TSV → 驗證：7 dataset rows 齊全、region key 唯一。
2. 依 rooted shape metadata 重算五類 → 驗證：每列五類回加 complete；complete+incomplete 回加 primary。
3. 建立 exact/overlap/anchor 一對一配對 → 驗證：每 scenario A/B matched keys 各自唯一。
4. 重算 agreement/κ/BA/Jaccard/null → 驗證：seed=20260712、5,000 permutations、chromosome-preserving。
5. 追到 layered raw exact tree digest → 驗證：不以 coarse shape 猜 exact tree equality。

- Checks：**64/64 PASS**。
- 主要資料表：`/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_dataset_summary.tsv`、`/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv`、`/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_match_metrics.tsv`、`/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv`。

## 7. 完整執行命令、輸入、輸出與實際片段

### 執行命令

```bash
python3 /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/scripts/topology_pair_analysis.py --exact-catalog /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json --ct-report /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json --regional-composition /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json --region-tsv /big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv --output-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data --output-report /big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/topology_pair_analysis.md
```

### 輸入

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json`（SHA-256 `e0f29a01f3ddffd4a76083f922eed713f0fd1fbf90b2aaf4df497548c4f82c58`）
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json`（SHA-256 `df60db7e26bab2f47620caa0a3cb6583d6d5488ec8ae0f574bf2fab4acf9a759`）
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json`（SHA-256 `b71c663cf7b31c4e0fd18ad6fef9f878cec7ae6763f0438a8cd14de914150b67`）
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_regions.tsv`（SHA-256 `a6649bf12dba686fad3d74d33243e16f4ec2b74361cac5a6b6d09eb09b794996`）
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395/layered_reconstruction_HCC1395.json`（SHA-256 `6b0398c4f32cdc1f8380e675440c2bffb5a83e159768b9111ce9596cf3280b60`）
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/samples/HCC1395_DORADO/layered_reconstruction_HCC1395_DORADO.json`（SHA-256 `16f38a178952f1ae5ef1c5ac3d2bd9fd1db93c162bb945b2a4b2946d37aa7e7f`）

### 輸出

- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/topology_pair_analysis.json`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_dataset_summary.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_all_regions.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_match_metrics.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_confusion.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_resolution_binary.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_interval_by_chrom.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/coarse_topology_golden_cases.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/topology_pair_validation_checks.tsv`
- `/big7_disk/liaoyoyo2001/InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/topology_pair_analysis.md`

### 實際輸出片段

```text
ACTUAL SUMMARY -> exact matched=6,252; complete_both=5,720; agreement=0.693881; kappa=0.497329; null_p=0.000200
VALIDATION -> 64/64 PASS
```
