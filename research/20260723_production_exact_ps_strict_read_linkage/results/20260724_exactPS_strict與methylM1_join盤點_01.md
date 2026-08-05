<!--
建立時間: 2026-07-24
目標: 核對 latest production exact-PS strict read-linkage 的 W/k 正式定義、全 7 datasets 與 HCC1395 數量，並鎖定與 methyl M1 site table 的可稽核 join contract
處理範圍: 7 technical datasets × chr1–22；primary threshold=3；唯讀重算
關聯檔案:
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/READY.json
  - InterSubMod/docs/CURRENT_FOCUS.md
服務目標: G3（read-level epigenetic integration）、G4（多樣本一致性與 reproducibility）
Task type: B — Comprehensive validation
-->

# Latest exact-PS strict read-linkage 與 methyl M1 join 盤點

用 PREP：**權威輸出已通過 7 datasets × chr1–22 完整性檢查；正式 `W` 是 exact-PS×HP 容器內的無向 read-linkage component，不是演化樹。全資料有 85,621 個 `k≥2` W、443,349 個 W node memberships；HCC1395 有 11,462 個 W、34,267 個 W memberships。M1 與 strict graph 可用 frozen `(dataset, chrom, pos)` 做 site-level left join；v8 的 7/7 formal pairs 也都在 threshold=3 共享同一 W，但這不等於 ancestry。**（影響：高；信心：高）

## 1. 正式定義

### 1.1 容器、邊、W 與 k

- 容器  
  `H = (dataset, chromosome, exact nonmissing PS, primary HP family)`。Primary HP family 只含 HP1／HP2；HP3、HP4、unphased 與 missing PS 不進 primary graph。
- 節點  
  在同一 `H` 中，至少被一條 canonical molecule 以 fixed `R/A` 觀察到的 sSNV genomic locus membership。
- 合格邊  
  同一 canonical molecule fixed-observes 兩個 endpoint；primary threshold 要求至少 3 個 distinct canonical molecules。邊是 **undirected read-linkage**。
- 正式區域  
  `W_i` 是同一 `H` 中由合格 endpoint edges 形成的 maximal connected component，且 `k_i=|W_i|≥2`。
- `k_sSNV`／`k`  
  是 `W_i` 內的 sSNV node membership 數；不是 bp span、不是 read 群數 `C`，也不是 clone 數。
- 顯示區間  
  `start=min(pos in W)`、`end=max(pos in W)`、`span=end−start` 只是 coordinate envelope；區間內沒有成為 W node 的其他位點不會自動加入。
- 大區塊  
  `k>12` 的 `W_i` 可切成 computational blocks `B_ij`；`B` 不是新的生物區域，必須用 source component ID 回聚到原始 `W_i`。

### 1.2 `k=1` 的處理

`k=1` component 被記為 `ABSTAIN_SINGLETON_UNLINKED`：

- 表示該 `(locus, exact PS, HP)` membership 在 threshold=3 沒有合格 edge。
- 保留在 component audit registry。
- 不進 W、k≤12 partition、Steiner tree、topology 或 clone/subclone 分母。

同一 physical locus 可在不同 PS／HP 容器中有多個 memberships，因此可能在一個容器是 singleton、另一個容器又屬於 W。排除必須在 membership 層做，不能把該 genomic locus 全域刪除。

### 1.3 claim ceiling

目前 production strict 結果完成的是：

1. exact-PS×HP read-linkage container；
2. threshold-qualified undirected endpoint edges；
3. singleton 與 `k≥2` W partition。

目前 **不能**由此直接主張：

- R→A、A→R、父→子或 mutation ancestry；
- 唯一真實演化樹；
- cellular clone 數；
- exact cellular parent-child；
- 跨 HP 生物 fusion。

`RA/AR` 只表示 genomic-order 兩個 endpoint 的觀察狀態。權威 `topology_completion_status.tsv` 對 7/7 datasets 的 strict directed topology 均為 `NOT_RUN` 或 `NOT_RUN_LATEST_V4_PARTITION_ONLY`，VAF/read-likelihood ranking 也為 `NOT_PRODUCTION_VALIDATED`。

## 2. 全資料與 HCC1395 數量

### 2.1 分母守恆

| Scope | Candidate sSNV `S` | Active unique physical loci | Active node memberships | `k=1` memberships | `k≥2` W | W node memberships | Unique physical loci in ≥1 W |
|---|---:|---:|---:|---:|---:|---:|---:|
| 7 datasets 合計 | 469,849 | 342,374 | 613,480 | 170,131 | 85,621 | 443,349 | 263,127 |
| HCC1395 | 79,687 | 36,384 | 62,651 | 28,384 | 11,462 | 34,267 | 22,689 |

解讀：

- `W node memberships = Σ_W k`，是容器 membership 數；全資料 443,349，HCC1395 34,267。
- `Unique physical loci in ≥1 W` 是 dataset 內以 genomic coordinate 去重後的位點數；全資料逐 dataset 加總 263,127，HCC1395 22,689。
- 上述兩者不可混用。全資料也不是跨不同 datasets 再做 biological-sample 去重。
- 守恆：
  - 全資料：`613,480 = 170,131 + 443,349`；`255,752 components = 170,131 k1 + 85,621 W`。
  - HCC1395：`62,651 = 28,384 + 34,267`；`39,846 components = 28,384 k1 + 11,462 W`。
- W memberships 佔 active memberships：全資料 72.2679%，HCC1395 54.6951%。
- Unique W loci 佔 S：全資料 56.0025%，HCC1395 28.4726%。

### 2.2 `k≥2` 分布

| Scope | k band | W 數 | W 比例 | Node memberships | Membership 比例 |
|---|---|---:|---:|---:|---:|
| ALL | 2 | 40,599 | 47.4171% | 81,198 | 18.3147% |
| ALL | 3 | 16,417 | 19.1740% | 49,251 | 11.1089% |
| ALL | 4 | 8,235 | 9.6180% | 32,940 | 7.4298% |
| ALL | 5 | 4,758 | 5.5570% | 23,790 | 5.3660% |
| ALL | 6–8 | 6,799 | 7.9408% | 46,047 | 10.3862% |
| ALL | 9–12 | 3,490 | 4.0761% | 35,737 | 8.0607% |
| ALL | >12 | 5,323 | 6.2169% | 174,386 | 39.3338% |
| **ALL total** | **≥2** | **85,621** | **100%** | **443,349** | **100%** |
| HCC1395 | 2 | 6,840 | 59.6754% | 13,680 | 39.9218% |
| HCC1395 | 3 | 2,578 | 22.4917% | 7,734 | 22.5698% |
| HCC1395 | 4 | 1,041 | 9.0822% | 4,164 | 12.1516% |
| HCC1395 | 5 | 455 | 3.9696% | 2,275 | 6.6390% |
| HCC1395 | 6–8 | 373 | 3.2542% | 2,453 | 7.1585% |
| HCC1395 | 9–12 | 85 | 0.7416% | 855 | 2.4951% |
| HCC1395 | >12 | 90 | 0.7852% | 3,106 | 9.0641% |
| **HCC1395 total** | **≥2** | **11,462** | **100%** | **34,267** | **100%** |

HCC1395 的 W：mean k=2.9896、median=2、p90=4、p95=5、p99=10、max=153。全資料最大觀察 k=567（H2009）。大 W 雖少，承載的 membership mass 可很高，因此報告 W 數與 node memberships 應同時呈現。

## 3. 與 methyl M1 site table 的 join contract

### 3.1 Frozen versions

| 角色 | Authority | Schema／狀態 |
|---|---|---|
| 全資料 strict bundle | `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/READY.json` | `intersubmod.strict_region_all7_bundle_ready` v1.0.0；`all_pass=true`；7×22 |
| strict report data | 同 bundle 的 `data/receipt.json` | `intersubmod.strict_region_all7_report_data_receipt` v1.1.0；25/25 checks true |
| raw strict chromosome | 各 chromosome `receipt.json` | `intersubmod.strict_exact_ps_hp_regions` v1.1.0；primary threshold=3 |
| methyl M1 site table | `.../all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz` | run manifest `intersubmod.all_ssnv_focal_alt_multigroup_run_manifest` v1.2.0；`EXECUTION_PASS` |

M1 table SHA-256：

`a8871af3a8c3955bf31aec5eeef0c93aca0683f52cf6d6f1e06fbbb713324f74`

M1 stable assignments SHA-256：

`82af076d8ce70c66c7f75c4ed4bdeacb7d4c5c43d0905859303b121df483a4ba`

### 3.2 Site-level join

M1 site table有：

`(dataset, chrom, pos, ref, alt)`

strict membership table有：

`(dataset, chrom, pos1, phase_set, linkage_basis, threshold, component_id, tree_eligible)`

Frozen M1 table的 469,849 rows 經重算：

- `(dataset, chrom, pos)`：469,849 unique，0 duplicate keys。
- `(dataset, chrom, pos, ref, alt)`：469,849 unique，0 duplicate keys。
- HCC1395：79,687 rows／79,687 unique coordinates，0 duplicate keys。
- 七個 dataset 的 M1 row counts 與 strict `candidate_loci_S` 逐 dataset 完全相等。

所以在這兩個 **版本鎖定** 的 authority 間，physical-site join 可用：

```text
m1.dataset = membership.dataset
AND m1.chrom = membership.chrom
AND m1.pos = membership.pos1
AND membership.threshold = 3
```

建議由 M1 做 `LEFT JOIN`，保留三種可分辨結果：

1. 沒有 strict row：沒有 active exact-PS×primary-HP fixed R/A membership；
2. strict `k=1` row：有 active membership，但 threshold=3 無合格 read-linkage edge；
3. strict W row：位點 membership 屬於 `k≥2` read-linked W。

### 3.3 一對多與 component join

site-level join **不是保證一對一**。同一 physical site 可能同時出現在 HP1、HP2 或多個 exact PS memberships；join 後要保留：

`phase_set, linkage_basis, component_id, component_class, tree_eligible`

若要取得 W 的 `k/start/end/span`，再以：

```text
(dataset, chrom, threshold=3, component_id)
```

把 membership join 到 components。不要用 `component_id` 當 M1 site key。

若將來更換 M1／VCF authority，strict membership 本身沒有 `ref/alt`，必須先重新驗證 `(dataset,chrom,pos)` 的 biallelic 唯一性，或從同一 frozen VCF 補入 `ref/alt`，不可默認座標永遠唯一。

## 4. v8 formal pairs 對 latest W 的核對

### 4.1 Authority 與分類準則

v8 pair authority：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/methyl_ssnv_pair_results.tsv.gz`

SHA-256：

`62f72f9535e88cb1cbf16872e566954b41de00bf465d040b983ed8ca9139bef1`

- run receipt：`intersubmod.methyl_ssnv_cooccurrence_validation` v2.0.0；task B、full scope、`pass=true`。
- release receipt：`intersubmod.cooccurrence_release_receipt` v1.0.0；`RELEASE_RECONCILIATION_PASS`。
- release 獨立重算：407,738 focal-partner pairs 中，`endpoint_a_formal_pair_by_confirmed=true` 共 7 筆。

分類規則：

- **同 W**：至少存在一組 threshold=3 memberships，使兩端的  
  `(dataset, chrom, exact phase_set, linkage_basis, component_id)` 完全相同，且兩端都是 `READ_LINKED_MULTISITE_REGION`、`tree_eligible=true`。
- **同 exact 容器不同 W**：兩端有相同 `(dataset, chrom, phase_set, linkage_basis)`，但 component IDs 不同。
- **無共同 exact membership**：找不到兩端共同的 exact PS×HP 容器。

### 4.2 七個 pair 的逐列結果

| Dataset | Focal ↔ partner | 共同 exact PS×HP | W k | Direct edge support | 分類 |
|---|---|---|---:|---:|---|
| H2009 | chr12:4413414 C>A ↔ 4414974 G>C | PS=2875649, HP2 | 10 | 87 | 同 W |
| H2009 | chr13:28815939 G>A ↔ 28817498 G>A | PS=28259005, HP1 | 7 | 74 | 同 W |
| H2009 | chr18:567920 T>C ↔ 563687 A>G | PS=307022, HP1 | 4 | 84 | 同 W |
| H2009 | chr4:2307521 T>C ↔ 2304921 G>T | PS=2165593, HP1 | 8 | 31 | 同 W |
| H2009 | 同上 | PS=2165593, HP2 | 8 | 68 | 同 W |
| H2009 | chr5:44615693 G>T ↔ 44612223 G>C | PS=44391260, HP1 | 62 | 155 | 同 W |
| H2009 | 同上 | PS=44391260, HP2 | 62 | 58 | 同 W |
| HCC1395 | chr3:127986757 G>A ↔ 127981978 C>G | PS=126975958, HP2 | 5 | 55 | 同 W |
| HCC1954 | chr8:100071382 A>G ↔ 100070832 C>T | PS=99963594, HP1 | 2 | 115 | 同 W |
| HCC1954 | 同上 | PS=99963594, HP2 | 2 | 63 | 同 W |

Pair-level 彙總：

```text
pairs=7
SAME_W=7
SAME_EXACT_CONTAINER_DIFFERENT_W=0
NO_COMMON_EXACT_MEMBERSHIP=0
shared_W_containers=10
```

這 7 個 pair 共命中 10 個共同 W containers；H2009 chr4、chr5 與 HCC1954 chr8 各同時命中 HP1、HP2。進一步逐列查 `endpoint_edges.tsv.gz`，10/10 containers 都有 focal-partner 的 direct edge，且 `passes_primary_threshold=true`；因此本次不是只靠中間節點的 transitive connectivity。

代表 component IDs：

```text
H2009:chr12:PS_HP2:PS0149795ba64f:E3:18:4368034-4414974:cd2ef3e4fd80fd39
H2009:chr13:PS_HP1:PSd7af88270806:E3:5:28782637-28872719:d6abc7f14f126b55
H2009:chr18:PS_HP1:PS4a242b301b8b:E3:4:563687-632412:a1cdb75c59e33390
H2009:chr4:PS_HP1:PS60aa69a23128:E3:3:2264402-2360811:6e956b5a64b73882
H2009:chr4:PS_HP2:PS60aa69a23128:E3:3:2264402-2360811:5f5dddb51221fc17
H2009:chr5:PS_HP1:PSb772b178dd9c:E3:2:44349008-44785725:08dba25256892d33
H2009:chr5:PS_HP2:PSb772b178dd9c:E3:2:44349008-44785725:b67da36cd479a6e3
HCC1395:chr3:PS_HP2:PScd848fb0dfa8:E3:9:127912458-127986757:485700e0e87ff350
HCC1954:chr8:PS_HP1:PS9e6111b08bce:E3:3:100070832-100071382:3dc98183b9751fcc
HCC1954:chr8:PS_HP2:PS9e6111b08bce:E3:3:100070832-100071382:a44eb407d5e52257
```

### 4.3 解讀上限

可以說：

> v8 篩出的 7 個正式 methyl-group–partner-sSNV association pairs，全部與 latest exact-PS×HP strict read-linkage W 相容；而且在共同容器內有 threshold-qualified direct read linkage。

不能說：

- 7 pairs 已驗證父子順序或 cellular ancestry；
- 7/7 是全體 pair 的準確率或 sensitivity；
- 這是完全獨立資料驗證，因為兩個分析產物仍共享上游 reads／variant set。

7/7 pair 的 `endpoint_b_status` 均為 `NOT_IDENTIFIABLE_*`，沒有唯一 mutation order。Direct edge 仍是無向 read-linkage；RR／RA／AR／AA 是 endpoint state counts，不是演化箭頭。

## 5. 實際列證據

### 5.1 HCC1395 chr22 一個 `k=5` W

Primary-threshold component：

```tsv
dataset  chrom  linkage_basis  phase_set  threshold  start1    end1      k  retained_edge_count
HCC1395  chr22  PS_HP1         10686903   3          10714123  10739431  5  10
```

其中一條 retained endpoint edge：

```tsv
dataset  chrom  linkage_basis  hp_family  phase_set  pos_i1    pos_j1    gap_bp  support_total  RR  RA  AR  AA  passes_primary_threshold
HCC1395  chr22  PS_HP1         1          10686903   10714123  10715921  1798    38             28  6   4   0   true
```

這列只證明兩 endpoint 在同一 exact PS×HP 容器具有 ≥3 canonical-molecule read support；`RA=6`、`AR=4` 都不是方向或 ancestry。

### 5.2 methyl 圖例位點 `HCC1395 chr22:47466517`

Latest strict threshold=3 membership：

```tsv
dataset  chrom  linkage_basis  phase_set  inference_role              component_class               tree_eligible  threshold  pos1      component_id
HCC1395  chr22  PS_HP1         46438690   ABSTAIN_SINGLETON_UNLINKED  UNLINKED_SINGLETON_COMPONENT  false          3          47466517  HCC1395:chr22:PS_HP1:PSb182b668b57d:E3:10:47466517-47466517:c38dc4ad68be398f
```

Component：

```tsv
dataset  chrom  linkage_basis  phase_set  threshold  start1    end1      k  solver_route
HCC1395  chr22  PS_HP1         46438690   3          47466517  47466517  1  EXCLUDE_SINGLETON_NO_READ_LINKAGE
```

因此該 methyl multigroup 圖例可用 coordinate join 找到 strict membership，但在 latest exact-PS strict graph 中是 `k=1` abstention，不是 W，也不能拿來聲稱 multi-sSNV topology。

## 6. 重算命令與實際輸出

### 6.1 輸入

- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/dataset_summary.tsv`
- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/exact_k_distribution.tsv`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2/chromosomes/chr22/*.tsv.gz`

### 6.2 M1 key 唯一性

執行命令核心：

```bash
gzip -dc /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz |
awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
{kc=$c["dataset"] FS $c["chrom"] FS $c["pos"];
 ka=kc FS $c["ref"] FS $c["alt"];
 coord[kc]++; allele[ka]++; n++}
END {print n, length(coord), length(allele)}'
```

實際輸出：

```text
469849 469849 469849
```

完整 duplicate audit 的實際結果：

```text
scope    rows    unique_coord  coord_duplicate_keys  unique_allele  allele_duplicate_keys
ALL      469849  469849        0                     469849         0
HCC1395  79687   79687         0                     79687          0
```

### 6.3 Summary 守恆

從 `dataset_summary.tsv` 逐 dataset 加總的實際輸出：

```text
datasets  S       active_unique_loci  unique_loci_in_any_W  active_node_memberships  singleton_memberships  W_memberships  all_components  k1_components  W_k_ge_2
7         469849  342374              263127                613480                   170131                 443349         255752          170131          85621
```

從 `exact_k_distribution.tsv` 以 `ΣW` 與 `Σ(k×W)` 重算：

```text
scope     W_total  node_memberships
ALL       85621    443349
HCC1395   11462    34267
```

### 6.4 v8 pair 對 strict memberships

先由 frozen pair TSV 抽出 formal positives：

```bash
gzip -dc /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/methyl_ssnv_cooccurrence_v8_m2v5_raw_identity_contract_source_locked_command_parity/methyl_ssnv_pair_results.tsv.gz |
awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
$c["endpoint_a_formal_pair_by_confirmed"]=="true" {
  print $c["sample"],$c["focal_chrom"],$c["focal_pos"],
        $c["partner_chrom"],$c["partner_pos"]
}'
```

每一 pair 再讀對應 chromosome authority：

```bash
gzip -dc DATASET.CHROM.site_component_membership.tsv.gz |
awk -F '\t' -v p1=FOCAL_POS -v p2=PARTNER_POS \
  '$9==3 && ($11==p1 || $11==p2) {
    print $1,$2,$3,$4,$7,$8,$9,$11,$12
  }'
```

最後以相同 `linkage_basis + phase_set + component_id` 配對兩端，實際輸出：

```text
pairs=7
SAME_W=7
SAME_EXACT_CONTAINER_DIFFERENT_W=0
NO_COMMON_EXACT_MEMBERSHIP=0
shared_W_containers=10
```

以相同 endpoint positions 查各 `endpoint_edges.tsv.gz`，實際為 10/10 `passes_primary_threshold=true`；support totals 依序為 `87, 74, 84, 31, 68, 155, 58, 55, 115, 63`。

### 6.5 輸出

本盤點：

`InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/results/20260724_exactPS_strict與methylM1_join盤點_01.md`

Production TSV／JSON 均未修改。
