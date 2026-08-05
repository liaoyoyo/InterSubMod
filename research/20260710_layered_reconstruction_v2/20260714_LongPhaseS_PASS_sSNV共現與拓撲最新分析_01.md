<!--
建立時間: 2026-07-14T05:56:30+08:00
目標: 以最新 LongPhase-S recalibrated FILTER=PASS canonical run 重算 7-dataset sSNV 共現、候選樹與 topology 狀態，並與 ClairS PASS sensitivity 對照
處理範圍: 7 datasets / 6 biological samples / chr1-22；canonical v5 與 sensitivity v6
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json；InterSubMod/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json
-->

# LongPhase-S PASS sSNV 共現與拓撲最新分析

## 結論

本次 **Task Type B comprehensive validation** 已完成。最新 canonical 輸入明確鎖定同一次 raw-all LongPhase-S production run 的 recalibrated `FILTER=PASS` VCF；7/7 dataset 的 sSNV 共現、partial/full read-state、候選樹、拓撲、hidden state 與全位點 ledger 都有實體輸出，canonical 與 ClairS-PASS sensitivity 都是 7/7 PASS。

canonical aggregate 為 582,820 個 LongPhase-S PASS records、469,849 個 chr1-22 biallelic sSNV、194,149 個 retained sSNV、51,815 個 regional groups。50,215 個有 primary HP1/HP2 lineage 的 regions 中，42,240 個 candidate-complete、7,975 個 capped/incomplete；三種合法 C/Topo 狀態分別為 11,582、10,737、19,921，數學不可能的 `C=1/Topo>1` 為 0。

## 最新輸入權威

Canonical run root：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`

7/7 `selected_tree_role` verifier checks 均為：

```text
expected      = longphase_s_recalibrated_FILTER_PASS
selected_role = longphase_recalibrated_pass_vcf
pass          = true
```

實際 tree VCF pattern：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/<sample>/<sample>.longphase_s.recalibrated.pass.vcf.gz`

這不是原始 ClairS PASS，也不是 2026-07-12 尚未完成的 layered v3 root。ClairS PASS 只存在於下列 sensitivity root：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6`

## sSNV 共現資料是否完整

共現 evidence 已寫入每個 sample 的五片 `mlhp_part_*.json` 與 `layered_reconstruction_<sample>.json`：

- full-span observed populations：`obs_populations` / `n_full_pops`。
- overlapping partial-read states：`obs_subreads` / `n_partial`。
- primary HP1/HP2 units：72,994。
- 同時有 full + partial evidence：28,322 units。
- partial-only，但可由重疊 reads 約束：44,672 units。
- full-only：0 units。
- read-tag alignment-group exposures：11,513,224。
- exact sidecar matches：11,513,224。
- missing、conflict、extra、malformed、allele conflict：全部 0。
- mixed-PS regions：5,623；PS 只作 phase-block QC，不作 topology edge。

因此正式用語仍是 **partial-read-constrained regional mutation-state candidate tree**。partial-only 不代表缺資料或錯誤，也不等於有單一 read 跨完整 region；它表示候選集合由互相重疊的 read-state constraints 建立。

## Canonical C/Topo 結果

定義：

- `C_region = product(n_trees)`，跨該 region 的 analysis-complete primary HP1/HP2 units。
- `Topo_region = product(n_distinct_shapes_exact)`，跨相同 units。
- complete：所有 primary units 都 non-capped、candidate-complete、`verification_status=full_pass`，且完整執行 V1-V7。
- incomplete：至少一個 primary unit capped 或 candidate set 不完整；不得以 stored prefix 宣稱唯一。

| Dataset | W_tree | W_primary | Complete | Incomplete | C=1/T=1 | C>1/T=1 | C>1/T>1 | Impossible |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829 | 8,007 | 7,878 | 7,128 | 750 | 1,471 | 2,459 | 3,198 | 0 |
| H1437 | 9,238 | 9,005 | 7,236 | 1,769 | 1,803 | 2,695 | 2,738 | 0 |
| H2009 | 9,674 | 9,538 | 5,813 | 3,725 | 1,145 | 2,227 | 2,441 | 0 |
| HCC1395 | 8,222 | 7,932 | 7,151 | 781 | 2,019 | 1,406 | 3,726 | 0 |
| HCC1395_DORADO | 8,385 | 7,665 | 7,060 | 605 | 1,847 | 621 | 4,592 | 0 |
| HCC1937 | 3,612 | 3,585 | 3,470 | 115 | 1,609 | 662 | 1,199 | 0 |
| HCC1954 | 4,677 | 4,612 | 4,382 | 230 | 1,688 | 667 | 2,027 | 0 |
| **合計** | **51,815** | **50,215** | **42,240** | **7,975** | **11,582** | **10,737** | **19,921** | **0** |

Hidden-state 分布：

- complete 且 hidden=0：7,770 regions。
- complete 且 hidden>0：34,470 regions。
- incomplete：7,975 regions。

Candidate-combination 分布：`C=1` 11,582；`C=2` 9,726；`C=3` 6,801；`C=4` 782；`C=5` 743；`C=6` 2,173；`C>6` 10,433；incomplete 7,975。

## 驗證狀態

- canonical root `_SUCCESS`：存在。
- canonical verifier：`all_pass=true`、`n_pass=7`、`n_fail=0`、`error_codes=[]`。
- 每個 sample 的 `layered_reconstruction`、`layered_region_view`、site ledger 與五片 MLHP 都存在。
- 所有 eligible non-capped units：V1-V7 PASS，`n_verification_fail=0`、`n_eligible_skipped_V4V5=0`。
- 每個 sample 的 group count = region count；全體 51,815 = 51,815。
- 每個 complete region 均滿足 `C >= Topo >= 1`。
- `complete + incomplete = W_primary`，三類合法 topology + incomplete = W_primary。
- impossible topology：0。

H2009 有 4 個 canonical regions 同時包含 recurrence-required 與 capped primary units。legacy 單一 `region_determinacy` label 依 precedence 顯示 `has_recurrence`，但本次 candidate-level 彙總仍正確歸入 incomplete；這是 label precedence，不是遺失候選或誤判 complete。sensitivity 中同類情況有 2 個。

## 與 ClairS PASS sensitivity 的差異

| 指標 | LongPhase-S PASS canonical | ClairS PASS sensitivity |
|---|---:|---:|
| Tree-input records | 582,820 | 568,080 |
| chr1-22 biallelic sSNV | 469,849 | 455,210 |
| Retained sSNV | 194,149 | 182,400 |
| W_tree | 51,815 | 48,960 |
| W_primary | 50,215 | 47,407 |
| Complete | 42,240 | 39,883 |
| Incomplete | 7,975 | 7,524 |
| C=1/T=1 | 11,582 | 10,812 |
| C>1/T=1 | 10,737 | 11,194 |
| C>1/T>1 | 19,921 | 17,877 |
| Impossible | 0 | 0 |

跨 backbone 判定仍為 `backbone_sensitive`：最大 aggregate metric delta 2.056768 pp；最低 retained-position Jaccard 0.577257；最低 primary-unit-key Jaccard 0.474027；最低 shared topology digest concordance 0.936110。故 main result 必須使用 LongPhase-S PASS，ClairS PASS 僅作 sensitivity。

2026-07-11 教授版 historical 報告中的 `W_tree=48,959`、`W_primary=47,377`、complete 39,885 與三類 10,832/11,144/17,909，是 truth-BED-conditioned historical engineering snapshot，不能再代表目前 canonical。最新可引用數字以本報告與 machine-readable JSON 為準。

## 執行命令

輸入：canonical v5、sensitivity v6、backbone comparison v6 與既有 C/Topo 定義函數。

```bash
env PYTHONDONTWRITEBYTECODE=1 \
  python3 research/20260710_layered_reconstruction_v2/scripts/summarize_current_layered_topology.py \
  --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --sensitivity-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6 \
  --backbone-comparison research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json \
  --topology-module research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py \
  --output-json research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json \
  --output-tsv research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.tsv
```

實際輸出：

```text
ALL_PASS=true
CANONICAL_W_TREE=51815
CANONICAL_W_PRIMARY=50215
CANONICAL_COMPLETE=42240
CANONICAL_INCOMPLETE=7975
```

Machine-readable outputs：

- `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json`，SHA-256 `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`。
- `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.tsv`，SHA-256 `44af3dda7059b5d3601b5f46e3176644a83d1dba44cf292ad3d0e8841904fdb1`。

## Claim ceiling

以上結果確認的是 bulk ONT read-supported、regional mutation-state 候選集合與模型內 topology 唯一性。它不能單獨確認 cellular clone 數、clone fraction 或真實祖先關係；後者仍需要 single-cell、multi-region/longitudinal、colony 或合成 truth。這個限制不影響「流程已完整執行」的工程結論，但限制生物演化主張。
