<!--
建立時間: 2026-07-14 Asia/Taipei
目標: 唯讀稽核 layered_workstation 是否綁定最新 canonical v5，並驗證資料 provenance、grain、完整性、唯一性、一致性、時效性與展示語意
處理範圍: Task B — 全 7 datasets / 6 biological samples / chr1-22；canonical v5 為唯一主結果；ClairS PASS v6 僅作 sensitivity
關聯檔案: InterSubMod/docs/CURRENT_FOCUS.md; InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md; InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md; InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json; InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py; InterSubMod/docs/methodology/_assets/layered_workstation/index.html
任務類型: B — Comprehensive validation；不可 subset
服務目標: G3 read-level HP/auxiliary evidence 邊界；G4 多樣本一致與 reproducibility；G5 外部可驗證 provenance
狀態: canonical_v5_data_PASS; current_workstation_binding_NO-GO
限制: 本輪只寫本 audit note；未修改 generator、HTML 或 canonical data；未執行視覺稽核；v6 不得升格為主結果
-->

# Layered workstation canonical v5 資料、provenance 與 grain 稽核

用 Evidence-to-Claim + grain-first data-quality audit：先鎖定每個 artifact 的一列代表什麼，再驗證 completeness、uniqueness、validity、consistency、integrity、timeliness，最後判斷 HTML 可否宣稱 canonical。

> **TL;DR — canonical v5 原始資料本身通過全 7-dataset 稽核：7/7 `_SUCCESS`、63/63 manifest-declared outputs hash/size 相符、51,815 MLHP groups = 51,815 region rows、194,149 retained sites 集合逐位點相符、所有主鍵重複與跨層 mismatch 都是 0；但目前 workstation 是 2026-07-10 的 unversioned ClairS-PASS-like 舊快照，並以 `CURRENT / CANONICAL` 顯示，來源、分母與 capped/recurrence 語意都不符合 v5，因此目前展示 **NO-GO**（影響：高；信心：高）。** `[O-L1]`

## 0. 最終裁決

| 稽核面 | 裁決 | 影響 | 信心 | 直接理由 |
|---|---|---:|---:|---|
| canonical v5 run root | **PASS** | 高 | 高 | `_SUCCESS` 存在、run state=`SUCCEEDED`、verification=`7/7 PASS`。 |
| 7-sample region-view 可用性 | **PASS** | 高 | 高 | 7/7 `layered_region_view_<sample>.json` 存在、schema `2.0`、JSON 可解析；**不需要另做 region-view converter**。 |
| grain / key / conservation | **PASS** | 高 | 高 | group、region、detail、site 主鍵 duplicate=0；MLHP↔region↔ledger retained-site 集合完全一致。 |
| topology summary | **PASS** | 高 | 高 | 以既有 summarizer 重算後，canonical、sensitivity 與 comparison aggregate 和 current machine JSON 完全相等。 |
| 現行 workstation 的 canonical 綁定 | **FAIL / P0** | 高 | 高 | builder 硬編 7/10 舊路徑；現行 7 個來源無 schema/provenance；首頁數字不是 v5。 |
| 現行 landing-table 語意 | **FAIL / P0** | 高 | 高 | 用 legacy single-label `region_determinacy` 代替 candidate-complete C/Topo；`Has capped` 在 H2009 少算 4 個 recurrence+capped region。 |
| per-sample HTML freshness gate | **FAIL / P0** | 高 | 高 | `page_ready` 只檢查 HTML 是否存在，沒有把 HTML 綁回 source hash；`--index-only` 可產生新 index + 舊 sample pages。 |
| v6 sensitivity 角色 | **PASS as sensitivity only** | 中 | 高 | verdict=`backbone_sensitive`；不得覆蓋 v5 主結果。 |

**發布判斷：**在 workstation 重新綁定 v5、重建 7 張 sample pages、改正 W/C/Topo/capped 語意並加入 hash freshness gate 前，不得再把目前 `layered_workstation/index.html` 稱為 current canonical interface。`[I-L2]`

### 0.1 證據標記

| Tag | 定義 |
|---|---|
| `[F-L1]` | 直接由程式碼或 machine artifact 讀出的事實。 |
| `[O-L1]` | 本稽核實際執行、重新解析或重算得到的觀察。 |
| `[I-L2]` | 至少兩項 L1 證據支持的工程／方法判斷。 |
| `[U-L5]` | 本資料無法回答或本輪未執行的事項。 |

## 1. 任務 gate、scope 與 canonical authority

### 1.1 Task B 固定邊界

- Task type：**B Comprehensive validation**。
- Dataset scope：**7 datasets / 6 biological samples**：COLO829、H1437、H2009、HCC1395、HCC1395_DORADO、HCC1937、HCC1954。
- Primary genomic scope：**chr1–22**；不是含 chrX/chrY 的 whole-genome topology。chrX/chrY 只保留在 site ledger 的 out-of-scope census。
- Main backbone：`longphase_s_recalibrated_filter_pass`。
- Main claim ceiling：`regional mutation-state candidate trees; not confirmed cell clones`。
- v6 `clairs_FILTER_PASS_sensitivity`：只可放在明顯標示的 sensitivity 區塊，不能與 main cards/table 混算。

### 1.2 唯一 canonical root

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260713_layered_reconstruction_v3_raw_all_lps_pass_v5
```

Root authority 實讀值：

```json
{
  "run_id": "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5",
  "state": "SUCCEEDED",
  "dataset_count": 7,
  "biological_sample_count": 6,
  "mode": "full",
  "scope": "chr1-22",
  "verification_all_pass": true,
  "verification_n_pass": 7,
  "verification_n_fail": 0,
  "verified_at_utc": "2026-07-12T21:09:40.762388Z"
}
```

`_SUCCESS` 建立時間為 `2026-07-12T21:09:40.979899Z`（Asia/Taipei `2026-07-13 05:09:40 +08:00`）；`_SUCCESS` SHA-256=`d382e0f2e5ef5ea6ebe9145688bc8fa1fe2f0e040fb930228ee050472ae6313a`，verification SHA-256=`9cdce65e7cb7488c1ae7b51bf0f5c87cbafee53605442503887dc7cf81bf983e`。`[F-L1]`

### 1.3 本次讀取的 SoT snapshot

| Input | SHA-256 |
|---|---|
| `InterSubMod/docs/CURRENT_FOCUS.md` | `a9bdd0dc52a34ce2dd263c370c448b8283642b758de29c17c49dd2b31b54e218` |
| `InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md` | `768efa09d4e7632e7af530a4195bf7c76da8ca12fc47a226211eff9cfda5acd2` |
| `InterSubMod/research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md` | `6849f8161e8b8daffe9743e7e756deee6d496ea02af2f0f9edff5a26fdb37042` |
| `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` | `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7` |
| `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py` | `5254896fc6fe6edcd18315014588c31ab5d17e956d039a3336ccc4a65061b24a` |
| `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` | `60a2c6df6866c55a2b99c73db647dbff3b8ab047e069ab1a95c6dea1aa273070` |
| topology summarizer | `5a75056d8ab3a04f46e64d701b0d50289011d2ea5311fbbd73d0ebd21c55db69` |
| topology module | `250a41d22287b6605361b3e677da4ec766e6143fddbc18a62ba4296318fab0e7` |

## 2. 最新 per-sample canonical paths、freshness 與 parse 結果

以下相對路徑都以 §1.2 canonical root 為基準。每個 dataset 實查 **10/10** 必要產物：5 個 `mlhp_part_*.json`、layered JSON、region-view JSON、site-ledger TSV.gz、site-ledger summary JSON、output manifest JSON。

| Sample | Latest region-view path | mtime +08 | SHA-256 | 10/10 / parse |
|---|---|---|---|---|
| COLO829 | `samples/COLO829/layered_region_view_COLO829.json` | 2026-07-13 02:17:38 | `438950cf287e18802065bd67ce12e0ea5b6fc3b77eaac96e9fbd71543b54cd32` | PASS |
| H1437 | `samples/H1437/layered_region_view_H1437.json` | 2026-07-13 03:13:51 | `d6be2187ad209512b75d992d3e645a68391bee8150a7847f804637f1cc8d927e` | PASS |
| H2009 | `samples/H2009/layered_region_view_H2009.json` | 2026-07-13 03:48:38 | `77fdad8d4466a9fdd5c3e07f06ecbe6460e7dc169e6c4c9d8389b53360d54629` | PASS |
| HCC1395 | `samples/HCC1395/layered_region_view_HCC1395.json` | 2026-07-13 03:03:06 | `9e28d98087d73cf4ebbe2185ab8cb4159bbc382799291c8c4025adc55f34c9d8` | PASS |
| HCC1395_DORADO | `samples/HCC1395_DORADO/layered_region_view_HCC1395_DORADO.json` | 2026-07-13 04:51:26 | `0a4e9e2eb357d047f874faf28a25ccb39c91ee37283c617187a3699801b8b1b2` | PASS |
| HCC1937 | `samples/HCC1937/layered_region_view_HCC1937.json` | 2026-07-13 04:37:55 | `e238c6f28b55901f43b244d686a97c1be71132cd725c4dae9544d2bbdf0dc149` | PASS |
| HCC1954 | `samples/HCC1954/layered_region_view_HCC1954.json` | 2026-07-13 04:25:55 | `e5a42c012c1ea36c1bfb2b77fe9ff2b025007d8ecbbf2b1ffbf3c6fcfbac0e28` | PASS |

對應 layered reconstruction：

| Sample | Latest layered path | SHA-256 |
|---|---|---|
| COLO829 | `samples/COLO829/layered_reconstruction_COLO829.json` | `b1bc974f0d48f9d4c59c689ad01056c0da916a8beba0774611668391d4dc601d` |
| H1437 | `samples/H1437/layered_reconstruction_H1437.json` | `f248c6e933e27173681a87b2dc952d6db9a2d4cb3d6e3ef1d4dfe427be434407` |
| H2009 | `samples/H2009/layered_reconstruction_H2009.json` | `2cb00cdd45b7db326ad232b26e8f12fb8620ae39be64a8667c8e861f6aa7ca73` |
| HCC1395 | `samples/HCC1395/layered_reconstruction_HCC1395.json` | `666d2d8a54322e0b7e1fa174c7b6ff2878891c4c88243669c05b63695ef82099` |
| HCC1395_DORADO | `samples/HCC1395_DORADO/layered_reconstruction_HCC1395_DORADO.json` | `f75588601644a8dd5cb4ec1228b912eca711f397246f19f56ad1c1ae71ab15f5` |
| HCC1937 | `samples/HCC1937/layered_reconstruction_HCC1937.json` | `35a514b2c57c424e848276aebc4810a2c24e0e020bc03aa5860a73a8ec9b8b7b` |
| HCC1954 | `samples/HCC1954/layered_reconstruction_HCC1954.json` | `718d8914773f433962a96fa9b8171911558d98ef033b6a41d9346d52e805fb39` |

`[O-L1]` 實際解析 63 個 sample-level JSON documents、7 個 gzipped site ledgers，以及 root verification／success／state／current topology JSON；零 parse error。v5 的 region-view schema 已存在且符合現有 workstation 所需的核心 `census` / `regions` 結構，所以**不需要先建 converter**；需要的是正確 source binding 與 metric semantic refactor。

## 3. Grain 先行：每個數字的一列到底是什麼

### 3.1 全 cohort funnel

```text
638,259 raw ClairS ledger records
 ├─ 55,439 excluded_by_longphase_filter
 └─ 582,820 LongPhase-S recalibrated FILTER=PASS tree-input records
     ├─ 112,971 chrX/chrY out-of-scope records
     └─ 469,849 chr1–22 biallelic sSNV
         ├─ 50,432 positional singletons
         ├─ 225,268 MAX_SNV-excluded sites
         └─ 194,149 retained sites
             └─ 51,815 retained MLHP groups = 51,815 W_tree regions
                 ├─ 1,600 no primary HP1/HP2 mutation lineage
                 └─ 50,215 W_primary regions
                     ├─ 42,240 candidate-complete regions
                     └─ 7,975 incomplete regions
```

所有分支逐樣本與 aggregate 都守恆，`read_unsupported=0`。`[O-L1]`

### 3.2 Artifact grain 與合法 key

| Artifact / metric | Grain | 合法 key | Count | 不可混用的分母 |
|---|---|---|---:|---|
| site ledger | 1 raw caller record | `(sample, chrom, pos, ref, alt)` | 638,259 | 不能當 LPS PASS、retained 或 region 數。 |
| LPS tree input | 1 recalibrated PASS VCF record | VCF record key | 582,820 | 含 chrX/Y；不能當 chr1–22 scope。 |
| autosomal sSNV | 1 chr1–22 biallelic sSNV | VCF record key | 469,849 | 尚未經 singleton/MAX_SNV 分支。 |
| retained site | 1 進入 read census 的 sSNV | `(sample, chrom, pos)` | 194,149 | 一個 region 可含 2–8 sites。 |
| MLHP group / `W_tree` | 1 retained multilocus group | `(sample, region)` | 51,815 | 含 1,600 no-primary regions。 |
| `W_primary` | 1 至少有 1 個 mutation-bearing HP1/HP2 primary unit 的 region | `(sample, region)` | 50,215 | C/Topo 的 region denominator。 |
| layered detail unit | 1 emitted family/control/auxiliary unit | `(sample, region, family)` | 118,234 | 包含 primary、reference-only、unphased、H3 auxiliary。 |
| primary reconstruction unit | 1 mutation-bearing HP1/HP2 unit | `(sample, region, family∈{1,2})` | 72,994 | unit-level，不可與 W_primary 互換。 |
| `C_region` | complete region 內各 primary unit `n_trees` 的乘積 | `(sample, region)` | 42,240 有值 | incomplete 必須是 NA，不是 0。 |
| `Topo_region` | complete region 內各 primary unit exact unlabeled-shape count 的乘積 | `(sample, region)` | 42,240 有值 | 只能在 candidate-complete unit 上計算。 |

### 3.3 L0 unit role census

| Unit role | Count | 可否承擔 primary topology claim |
|---|---:|---|
| `primary_mutation_lineage` | 72,994 | 是；HP1/HP2 primary。 |
| `reference_only_control` | 10,680 | 否；不可算入 HP multiplicity。 |
| `unphased_auxiliary` | 24,594 | 否；auxiliary。 |
| `unresolved_H3_auxiliary` | 9,966 | 否；auxiliary。 |
| H4 auxiliary | 0 regions | 本 run 無 H4。 |
| **Total detail units** | **118,234** | — |

因此 v5 的 exact 定義是：

```json
"hp_multiplicity_definition": "mutation-bearing HP1/HP2 primary lineages"
```

aggregate HP multiplicity：`0=1,600`、`1=27,436`、`2=22,779`，總和 51,815。現行首頁「可能包含 root-only control、不等同 both-mutation-bearing」是舊定義，對 v5 已經錯誤。`[F-L1]`

## 4. 五大資料品質維度與逐樣本證據

### 4.1 DQ scorecard

| Dimension | Status | 實測證據 | Severity / action |
|---|---|---|---|
| Availability / completeness | PASS | 7/7 samples；每樣本 10/10 required artifacts；chr1–22 每樣本均 22 contigs。 | 無 data blocker。 |
| Uniqueness | PASS | duplicate region keys=0；detail keys=0；MLHP group keys=0；site keys=0。 | 無。 |
| Validity | PASS | region/layered/site-summary schema均 `2.0`；sample IDs 正確；region string、start/end、positions、n_sSNV 全一致。 | 無。 |
| Consistency | PASS | MLHP group set = region set；retained ledger site set = region positions；group→positions mapping完全一致。 | 無。 |
| Integrity | PASS with metadata warning | 63/63 manifest-declared outputs實際 size/hash相符；read-tag exact join 11,513,224/11,513,224。 | site-index provenance 見 §9.2。 |
| Timeliness — v5 data | PASS | latest sample artifact 2026-07-13 04:52；root verified 05:09。 | 無。 |
| Timeliness — workstation | **FAIL** | index generated 2026-07-10 23:53；source snapshot 15:31；未綁 v5。 | P0 rebuild。 |
| Semantic validity — workstation | **FAIL** | legacy region label、錯誤 HP copy、錯誤 avg denominator、無 W_primary/C/Topo。 | P0/P1 refactor。 |

### 4.2 逐樣本 row/key/completeness/uniqueness/consistency

`Dup` 欄依序為 region/detail/MLHP/site；`Mismatch` 是 MLHP↔region position set、ledger retained↔region site set 與 manifest output mismatch 的合計，全部為 0。

| Sample | Region rows | Detail rows | MLHP groups | Ledger rows | chr count | Dup R/D/G/S | Mismatch | Machine topology match |
|---|---:|---:|---:|---:|---:|---|---:|---|
| COLO829 | 8,007 | 17,613 | 8,007 | 43,014 | 22 | 0/0/0/0 | 0 | true |
| H1437 | 9,238 | 19,668 | 9,238 | 83,261 | 22 | 0/0/0/0 | 0 | true |
| H2009 | 9,674 | 21,973 | 9,674 | 168,638 | 22 | 0/0/0/0 | 0 | true |
| HCC1395 | 8,222 | 20,119 | 8,222 | 134,122 | 22 | 0/0/0/0 | 0 | true |
| HCC1395_DORADO | 8,385 | 20,267 | 8,385 | 123,240 | 22 | 0/0/0/0 | 0 | true |
| HCC1937 | 3,612 | 7,174 | 3,612 | 59,161 | 22 | 0/0/0/0 | 0 | true |
| HCC1954 | 4,677 | 11,420 | 4,677 | 26,823 | 22 | 0/0/0/0 | 0 | true |
| **Total** | **51,815** | **118,234** | **51,815** | **638,259** | **7×22** | **0/0/0/0** | **0** | **7/7** |

### 4.3 Manifest 與 read-tag integrity

- 7 個 output manifests 各宣告 9 outputs，共 **63/63**：actual path exists、size 相符、SHA-256 相符。
- read-tag `alignment_group_exposures=11,513,224`，`sidecar_exact_matches=11,513,224`。
- `sidecar_missing/conflicts/extra/malformed=0`，`alignment_identity_allele_conflicts=0`。
- post-run BAM/sidecar identity check 7/7 通過；所有 sample `check_exact_sidecar_join=true`。
- mixed-PS regions=5,623；**PS 只作 QC context，不是 topology edge 或 lineage label**。

## 5. v5 canonical 數字：workstation 應展示的主表

| Sample | LPS PASS tree input | chr1–22 sSNV | Retained sSNV | W_tree | W_primary | Complete | Incomplete | C1/T1 | C>1/T1 | C>1/T>1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829 | 39,458 | 37,788 | 27,756 | 8,007 | 7,878 | 7,128 | 750 | 1,471 | 2,459 | 3,198 |
| H1437 | 79,655 | 77,080 | 37,543 | 9,238 | 9,005 | 7,236 | 1,769 | 1,803 | 2,695 | 2,738 |
| H2009 | 161,595 | 154,465 | 51,164 | 9,674 | 9,538 | 5,813 | 3,725 | 1,145 | 2,227 | 2,441 |
| HCC1395 | 113,061 | 79,687 | 27,102 | 8,222 | 7,932 | 7,151 | 781 | 2,019 | 1,406 | 3,726 |
| HCC1395_DORADO | 113,145 | 79,739 | 27,615 | 8,385 | 7,665 | 7,060 | 605 | 1,847 | 621 | 4,592 |
| HCC1937 | 52,115 | 18,690 | 9,331 | 3,612 | 3,585 | 3,470 | 115 | 1,609 | 662 | 1,199 |
| HCC1954 | 23,791 | 22,400 | 13,638 | 4,677 | 4,612 | 4,382 | 230 | 1,688 | 667 | 2,027 |
| **Total** | **582,820** | **469,849** | **194,149** | **51,815** | **50,215** | **42,240** | **7,975** | **11,582** | **10,737** | **19,921** |

Invariants：

```text
W_tree = W_primary + no_primary = 50,215 + 1,600 = 51,815
W_primary = complete + incomplete = 42,240 + 7,975 = 50,215
W_primary = C1/T1 + C>1/T1 + C>1/T>1 + incomplete + impossible
          = 11,582 + 10,737 + 19,921 + 7,975 + 0
```

### 5.1 C、Topo 與 hidden-node 主語意

- `C_region`：candidate-complete primary HP1/HP2 units 的 exact tree counts 乘積。
- `Topo_region`：相同 units 的 exact unlabeled topology-shape counts 乘積。
- complete：所有 primary units 都是 non-capped、`verification_status=full_pass`、`verification_complete=true`、`analysis_candidate_set_complete=true`，且生成 tree count 與 stored analytical count 一致。
- incomplete：任一 primary unit 不完整時，C 與 Topo 都應顯示 **NA / incomplete**，不能顯示 0 或把 provisional stored trees 當全集。
- impossible `C=1, Topo>1`：0 regions，符合 `C_region ≥ Topo_region ≥ 1`。

C bins：

| C bin | Region count |
|---|---:|
| C=1 | 11,582 |
| C=2 | 9,726 |
| C=3 | 6,801 |
| C=4 | 782 |
| C=5 | 743 |
| C=6 | 2,173 |
| C>6 | 10,433 |
| incomplete | 7,975 |

Joint exact combinations total=`2,414,607`；joint topology combinations total=`115,088`；max C=`15,625`；max Topo=`81`。Hidden node classes：0 hidden=`7,770`、positive hidden=`34,470`、incomplete=`7,975`。

### 5.2 Full/partial read evidence

Primary reconstruction units=72,994：

| Coverage evidence | Units | Meaning |
|---|---:|---|
| full + partial | 28,322 | 同 unit 同時有 full-pattern 與 partial reads。 |
| partial only | 44,672 | 沒有 full-pattern population；candidate 是由 partial constraints 支持。 |
| full only | 0 | 本 run 無此類。 |

這是 unit-level evidence coverage，不是 region-level biological clone fraction。首頁若要展示，分母必須寫 `72,994 primary units`。

## 6. L0→L3、HP、CN、recurrence 與 capped 的正確界線

### 6.1 實際 contract excerpt

```json
{
  "schema_version": "2.0",
  "U1_backbone_source": "longphase_s_recalibrated_filter_pass",
  "U1_input_contract": "LongPhase-S _sc.vcf FILTER=PASS biallelic sSNV",
  "analysis_scope": "chr1-22 primary; chrX/chrY out-of-scope census only",
  "claim_scope": "regional mutation-state trees; not confirmed cell clones",
  "hp_multiplicity_definition": "mutation-bearing HP1/HP2 primary lineages",
  "L3_methyl": {
    "role": "bounded_auxiliary",
    "status": "not_evaluated",
    "allowed_uses": ["negative_screen", "residual_flag"],
    "prohibited_uses": ["tree_ranking", "lineage_confirmation", "clone_confirmation"]
  }
}
```

### 6.2 每一層可說與不可說

| Layer | v5 中的角色 | 可說 | 不可說 |
|---|---|---|---|
| L0 | HP family / role partition | HP1/HP2 primary 分開；H3、unphased、reference-only 另列。 | 不能把 HP3、none 或 root-only 算成 primary multiplicity。 |
| L1 | read-level sSNV candidate enumeration | 哪些 minimal candidate trees/shape 與觀測相容；是否 complete/capped。 | 不能等同細胞 clone 或真實 ancestry。 |
| L2 | post-tree CN context | 對 recurrence-like outcome 作 CN confound annotation。 | 不能回頭改寫 observed read state，也不能把 CN annotation 當 clone confirmation。 |
| L3 | bounded methyl auxiliary | 本 run 只標 `not_evaluated`；未來最多 negative/residual flag。 | 禁止 rank、select、confirm tree/lineage/clone。 |

CN availability：H1437、H2009、HCC1395、HCC1395_DORADO、HCC1954 為 measured；COLO829、HCC1937 為 unavailable。L3 是 **0/7 evaluated**。

### 6.3 Recurrence 與 capped 不是互斥集合

v5 legacy `region_determinacy` 是**單一優先標籤**，不是 multi-label truth。Aggregate legacy labels：

```text
all_determined=11,581
has_ambiguous=30,505
has_capped=7,971
has_recurrence=158
no_primary_lineage=1,600
```

H2009 有 4 個 region 同時含 recurrence-required primary unit 與 capped primary unit，但 legacy region label 只顯示 `has_recurrence`：

```json
[
  {"region":"chr8:79992384-80149786", "legacy":"has_recurrence", "HP1":"capped/incomplete", "HP2":"recurrence_required/CN-gain_confounded"},
  {"region":"chr13:93837736-93888639", "legacy":"has_recurrence", "HP1":"recurrence_required/CN-gain_confounded", "HP2":"capped/incomplete"},
  {"region":"chr15:31733893-31800487", "legacy":"has_recurrence", "HP1":"capped/incomplete", "HP2":"recurrence_required/LOH_unresolved"},
  {"region":"chr9:275701-337149", "legacy":"has_recurrence", "HP1":"capped/incomplete", "HP2":"recurrence_required/CN-gain_confounded"}
]
```

因此：

```text
legacy has_capped labels = 7,971
candidate-level incomplete regions = 7,975
discordance = 4（全部在 H2009）
```

首頁欄名 `Has capped` 卻直接讀 `region_determinacy.has_capped`，會漏掉這 4 個 region；這不是 rounding，而是 representation bug。應以 primary-unit flags 重算 `any_capped`，或直接用 candidate-complete `incomplete` 作主指標，legacy label 只放 diagnostics。

### 6.4 V1–V7 顯示必須保留 N/A

所有 118,234 detail units 中：

- `full_pass=109,435`；
- `not_applicable_capped=8,799`；
- capped units 的 V4/V5 各有 8,799 次 `skipped/not_applicable`。

所以首頁的單一 badge `V1–V7 PASS` 太寬。正確文字應是：**all non-capped eligible units full-pass；capped units V4/V5 not applicable and remain incomplete**。不能讓使用者誤讀為 capped candidate set 也已被 V4/V5 完整驗證。

## 7. 現行 7/10 workstation 與 canonical v5 的數值差異

### 7.1 現行 HTML 的真正 snapshot

`InterSubMod/docs/methodology/_assets/layered_workstation/index.html` 實際顯示：

```text
Data version / backbone = ClairS PASS / LongPhase-S _sc.vcf(移除 is_somatic 粗重檢)
Data snapshot = 2026-07-10T15:31+08:00
Page generated = 2026-07-10T23:53+08:00
Sample workstations = 7/7
Census inputs available = 7/7
L3 = PENDING 0/7
```

現行 index 生成時間比 v5 `_SUCCESS` 早 **53 小時 16 分**；其 source snapshot 比 v5 verification 早 **61 小時 38 分**。

### 7.2 Per-sample old → v5 main

| Sample | Current HTML `Somatic sSNV` → v5 LPS PASS tree input | Delta | Current W → v5 W_tree | Delta |
|---|---:|---:|---:|---:|
| HCC1395 | 113,997 → 113,061 | −936 | 7,927 → 8,222 | +295 |
| COLO829 | 38,196 → 39,458 | +1,262 | 7,774 → 8,007 | +233 |
| H1437 | 75,578 → 79,655 | +4,077 | 8,865 → 9,238 | +373 |
| H2009 | 157,405 → 161,595 | +4,190 | 9,717 → 9,674 | −43 |
| HCC1395_DORADO | 112,387 → 113,145 | +758 | 7,958 → 8,385 | +427 |
| HCC1937 | 49,548 → 52,115 | +2,567 | 2,695 → 3,612 | +917 |
| HCC1954 | 20,969 → 23,791 | +2,822 | 4,023 → 4,677 | +654 |
| **Total** | **568,080 → 582,820** | **+14,740** | **48,959 → 51,815** | **+2,856** |

目前 HTML 的 568,080 恰好等於 current v6 sensitivity 的 tree-input aggregate，但其 W=48,959 又比 final v6 的 48,960 少 1。結論是：它是**歷史 ClairS-PASS-like snapshot**，既不能稱 v5 main，也不應直接重貼標籤稱 final v6 sensitivity。

### 7.3 v6 僅作 sensitivity

| Metric | v5 canonical main | v6 ClairS PASS sensitivity |
|---|---:|---:|
| Tree input | 582,820 | 568,080 |
| chr1–22 biallelic sSNV | 469,849 | 455,210 |
| Retained sSNV | 194,149 | 182,400 |
| W_tree | 51,815 | 48,960 |
| W_primary | 50,215 | 47,407 |
| Complete | 42,240 | 39,883 |
| Incomplete | 7,975 | 7,524 |

Sensitivity verdict=`backbone_sensitive`；min retained-site Jaccard=`0.577257`、min primary-unit-key Jaccard=`0.474027`、shared topology-digest concordance=`0.936110`。因此 v6 可以說明共享 unit 上 topology 高度一致，但不能替代 v5，也不能把 marginal percentage similarity寫成 backbone invariance。

## 8. 現行 builder / HTML 的具體失效點

### 8.1 P0 — hardcoded source 指向舊、無 provenance 的資料

`build_layered_per_sample.py:31-48`：HCC1395 指向 repo 內 `20260618_subcluster_pilot`；其餘六個指向 `/big7.../multisample_subclone/<sample>/...`。實查 7/7 都存在，但：

- `schema_version` 全部缺失；
- `provenance` 全部缺失；
- `analysis_scope`、`claim_scope`、`hp_multiplicity_definition` 全部缺失；
- SHA 與 v5 region views 全部不同；
- backbone 文字為歷史 `ClairS PASS / LongPhase-S _sc.vcf...`。

然而現行 sample pages 顯示 `CURRENT / CANONICAL`。這是 provenance mislabel，不只是 freshness warning。

### 8.2 P0 — `page_ready` 只代表檔案存在

`build_layered_per_sample.py:86-88, 209-223, 428-439`：

- `page_ready = os.path.exists(output)`；
- cohort `complete` 只看 page exists；
- 沒有比對 HTML 內的 source SHA、run ID、generated-from timestamp；
- `--index-only` 明確不重建 sample HTML。

所以只改 index source 後執行 footer 建議的 `--index-only`，會得到 **v5 index + 7/10 old sample pages**，卻仍顯示 7/7 complete。這必須 fail closed：每張 page 應嵌入 source SHA，index 只在 `embedded_source_sha == expected_v5_region_sha` 時標 complete。

### 8.3 P0 — legacy region labels 被當主 topology result

`build_layered_per_sample.py:112-128` 直接讀：

- `n_regions`；
- `region_determinacy.all_determined/has_ambiguous/has_capped/has_recurrence`；
- `hp_multiplicity["2"]`。

這些可以作 legacy census diagnostics，但不能取代 v5 的 `W_primary`、candidate complete/incomplete、C/Topo classes。尤其 `has_capped` 是單一 label count，不是 any-capped boolean。

### 8.4 P1 — `Avg trees / lineage` 分母不一致

目前公式：

```python
sum_ntrees_noncapped / n_lineage_units
```

分子排除 capped，分母卻包含 capped primary units，得到的不是 non-capped 平均，也不是 joint `C_region`。此欄應移除；若保留，必須把分母改為 analysis-complete non-capped primary units，並明示是 **per-unit arithmetic mean**，不可與 region-level joint C 混稱。

### 8.5 P1 — L3 status 讀錯 schema 位置

目前 `l3_ready = bool(census.get("L3"))`；v5 contract 位於 top-level `L3_methyl.status="not_evaluated"`。目前顯示 pending 只是因缺 key 而碰巧正確，未來即使加入明確 L3 payload也可能判斷錯誤。應讀 exact status enum，而不是對任意 dict 做 truthiness。

### 8.6 P1 — `Data snapshot=max(mtime)` 不是 atomic provenance

`cohort_metadata()` 取 7 個 source mtimes 的最大值。這不能證明 7 個檔案屬於同一 run，也不能證明 verification 完成。主頁應展示：

- run ID；
- main backbone role；
- scope 7/6/chr1–22；
- `_SUCCESS.created_at_utc`；
- verification SHA；
- current topology summary SHA；
- 每頁 region-view SHA。

### 8.7 P1 — current table 的 copy 對 v5 已失真

- `Somatic sSNV` 未說明是 raw ledger、LPS PASS 還是 chr1–22；三者分別為 638,259 / 582,820 / 469,849。
- `Eligible regions` 未區分 W_tree 與 W_primary。
- HP multiplicity copy 說可能含 root-only，與 v5 exact definition 相反。
- `All determined` 不等於 `C1/T1` 的所有情況；H2009 有 1 個 recurrence-required unit 形成 C1/T1，造成 legacy all_determined 11,581 vs canonical C1/T1 11,582。
- `Has capped` 漏 recurrence+capped overlap。
- `V1–V7 PASS` 未揭露 capped V4/V5 N/A。

## 9. Canonical data 仍有三個低於 workstation P0 的 provenance debt

這兩項不改變 v5 topology 數字，也不阻擋 workstation 讀 region-view，但在 G5 外部交付前應修。

### 9.1 P2 — site-ledger summary 內嵌 staging path 已失效

7/7 `ssnv_site_ledger_<sample>.summary.json` 的：

```json
"output_tsv_gz": ".../<sample>/.site_ledger.stage.<pid>/ssnv_site_ledger_<sample>.tsv.gz",
"output_index": ".../<sample>/.site_ledger.stage.<pid>/ssnv_site_ledger_<sample>.tsv.gz.tbi"
```

都指向已不存在的 staging path。實際 final TSV.gz 與 `.tbi` 都存在，而且 `tabix -l` 7/7 可讀 chr1–22、chrX、chrY；manifest 也正確綁 final TSV.gz。最小修正是 producer 在 atomic publish 後，把 summary 寫成 logical final path，或用 relative path。

### 9.2 P2 — final `.tbi` 可用但未納入 output manifest

7 個 final `.tbi` 都存在且可查詢，但 7 個 output manifests 的 9 outputs 不含 index；因此本次「63/63 hash/size PASS」不包含 `.tbi`。G5 handoff 應把 `.tbi` 納入 manifest role/hash/size。

### 9.3 P3 — verification summary 的 marker 欄位是 pre-marker 語意

`verification_summary.json.success_marker_created=false`，但 `_SUCCESS` 實際存在且時間晚約 0.22 秒，run state 也是 `SUCCEEDED`。這是 lifecycle 順序造成的欄名歧義，不是 run failure；外部介面應以實際 `_SUCCESS` readback 為準，或另寫 post-success receipt，避免把 `false` 誤顯示為失敗。

## 10. 建議的 workstation canonical contract（不在本輪實作）

### 10.1 Source gate

1. 唯一 main root 固定為 §1.2 v5；不可從 `multisample_subclone` 或 pilot path 回退。
2. 開始 build 前要求 `_SUCCESS`、run state `SUCCEEDED`、verification `all_pass=true,n_pass=7,n_fail=0`。
3. 對 7 個 region-view 與 layered JSON 重算 SHA，和 machine topology summary 記錄一致才繼續。
4. 每張 HTML 嵌入 `run_id/source_path/source_sha/generated_at`；index 驗證 7/7 embedded SHA 後才顯示 complete。
5. v6 放獨立 sensitivity tab/card，永遠帶 `NOT MAIN` badge。

### 10.2 Landing page 的最小正確資訊架構

1. **Authority strip**：run ID、LPS recalibrated PASS、7/6、chr1–22、verified time、hash。
2. **Funnel**：638,259 → 582,820 → 469,849 → 194,149 → 51,815。
3. **W definitions**：W_tree 51,815；W_primary 50,215；no-primary 1,600。
4. **Candidate topology**：C1/T1、C>1/T1、C>1/T>1、incomplete，分母 W_primary。
5. **Evidence coverage**：full+partial / partial-only，分母 72,994 primary units。
6. **Per-sample table**：使用 §5 數字；欄名內直接寫 grain/denominator。
7. **Auxiliary status**：CN 5/7 measured；L3 0/7 evaluated；PS QC-only。
8. **Diagnostics drawer**：legacy region labels、recurrence、any-capped overlap、V1–V7 N/A 明細。

### 10.3 No-converter 結論

7/7 v5 region-view 都已存在、schema一致、keys完整，且現行 landing builder 所讀的基本 census keys 都可取得。**不需要另建 data converter 或新 region-view generator。**應直接讓 builder 讀 v5 paths；但不能只換 `SAMPLES`：還必須一起修正 W/C/Topo、L3 status、capped overlap、avg denominator、hash freshness gate，並完整重建 7 張 sample pages。

## 11. 實際執行命令、輸入、輸出與驗證片段

### 11.1 Canonical/sensitivity 第二路重算

輸入：v5 root、v6 sensitivity root、backbone comparison JSON、current topology module。

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260710_layered_reconstruction_v2/scripts/summarize_current_layered_topology.py \
  --canonical-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5 \
  --sensitivity-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6 \
  --backbone-comparison research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json \
  --topology-module research/20260710_layered_reconstruction_v2/scripts/build_layered_observation_report.py \
  --output-json /tmp/layered_workstation_v5_data_audit_recomputed.json \
  --output-tsv /tmp/layered_workstation_v5_data_audit_recomputed.tsv
```

實際輸出片段／exit：

```text
ALL_PASS=true
CANONICAL_W_TREE=51815
CANONICAL_W_PRIMARY=50215
CANONICAL_COMPLETE=42240
CANONICAL_INCOMPLETE=7975
CANONICAL_TOPOLOGY={"exact_and_topology_unique":11582,
  "topology_unique_exact_multiple":10737,
  "topology_multiple_exact_multiple":19921,
  "incomplete":7975,
  "impossible_exact_unique_topology_multiple":0}
exit=0
```

讀回比較：

```text
canonical_equal True
sensitivity_equal True
comparison_aggregate_equal True
```

### 11.2 Full path / JSON / grain / key / set / manifest audit

輸入：v5 `samples/*` 下的 7×10 artifacts；current topology JSON；verification summary。

核心檢查的可重現邏輯：

```python
region_key = (sample, region["region"])
detail_key = (sample, unit["region"], str(unit["family"]))
group_key = (sample, f'{group["chrom"]}:{group["start"]}-{group["end"]}')
site_key = (sample, chrom, pos, ref, alt)

assert len(regions) == census["n_regions"] == funnel["n_groups_retained"]
assert set(region_keys) == set(mlhp_group_keys)
assert retained_ledger_position_set == region_position_set
assert retained_ledger_group_to_positions == region_to_positions
assert len(detail) == layered["n_detail_units"] == verification_metrics["n_detail_units"]
assert manifest_sha256_and_size_mismatches == 0
assert duplicate_region == duplicate_detail == duplicate_group == duplicate_site == 0
assert read_tag_exposures == read_tag_exact_matches
```

實際輸出片段：

```text
SAMPLE             paths region detail mlhp ledger_rows chroms dupR dupD dupG dupS mismatch manifest topo_match
COLO829             10/10   8007  17613 8007       43014     22    0    0    0    0        0        0 true
H1437               10/10   9238  19668 9238       83261     22    0    0    0    0        0        0 true
H2009               10/10   9674  21973 9674      168638     22    0    0    0    0        0        0 true
HCC1395             10/10   8222  20119 8222      134122     22    0    0    0    0        0        0 true
HCC1395_DORADO      10/10   8385  20267 8385      123240     22    0    0    0    0        0        0 true
HCC1937             10/10   3612   7174 3612       59161     22    0    0    0    0        0        0 true
HCC1954             10/10   4677  11420 4677       26823     22    0    0    0    0        0        0 true
PROBLEMS []
```

### 11.3 Site-ledger index readback

```bash
ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/samples
for s in COLO829 H1437 H2009 HCC1395 HCC1395_DORADO HCC1937 HCC1954; do
  printf '%s\t' "$s"
  tabix -l "$ROOT/$s/ssnv_site_ledger_$s.tsv.gz" | paste -sd, -
done
```

實際：7/7 exit 0，每個都列出 `chr1,...,chr22,chrX,chrY`。

### 11.4 本輪唯一 repo 輸出

```text
InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/10_layered_workstation_v5_data_audit.md
```

未修改：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py`
- `InterSubMod/docs/methodology/_assets/layered_workstation/index.html`
- canonical v5 / sensitivity v6 data roots

## 12. Claim-to-evidence mapping

| Claim | Direct evidence | Verdict |
|---|---|---|
| v5 是唯一 main | CURRENT_FOCUS、latest analysis、00_INDEX、v5 `_SUCCESS` | PASS |
| 7/7 region views 可直接供 workstation | 7 paths、schema 2.0、parse、key audit | PASS；no converter |
| v5 數字可重現 | existing summarizer re-run、current JSON exact equality | PASS |
| current HTML 是 v5 | source paths、mtime、old totals、missing provenance | **FAIL** |
| current `Has capped` 是 any-capped | H2009 4 recurrence+capped excerpts | **FAIL** |
| current table 分母足以解釋 topology | 缺 W_primary、C/Topo、complete contract | **FAIL** |
| L3 已可展示結果 | 7/7 `status=not_evaluated` | **FAIL；只能顯示 unavailable/not evaluated** |
| regional candidates 是 confirmed clones | claim_scope 明確禁止 | **FAIL / prohibited claim** |
| v6 可取代 v5 | sensitivity verdict=`backbone_sensitive` | **FAIL** |

**最終結論：**資料修復不是當前 blocker；展示層的 source authority、grain、denominator 與 semantic binding 才是 blocker。先以 v5 直接重建，再做桌面／手機視覺稽核；任何視覺 polish 都不能先於這個 canonical gate。
