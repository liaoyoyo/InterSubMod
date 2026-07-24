<!--
建立時間：2026-07-24
目標：對 20260724 exact-PS C++ 全域最佳 AF 樹做精確 topology-signature census
處理範圍：唯讀重算既有 MLHP 與 canonical topology JSONL；不修改 canonical outputs
關聯檔案：cpp/exact_topology_signature_census.cpp
-->

# exact-PS C++ 最佳樹 topology-signature census

此 companion tool 使用與 2026-07-24 canonical run 相同的：

- MLHP input；
- frozen C++ parser、obligation B&B 與 exact rational AF score；
- `max_family_size=100000`、`max_search_nodes=1000`。

它只對 canonical `ranked_complete` units 重算，展開每個 global-best
vertex set 的 per-child best-parent Cartesian product。每棵最佳樹先轉成
root-preserving、unlabeled rooted-tree canonical parenthesis signature，
再分類為：

- `UNIQUE_TREE`：最佳完整樹只有一棵；
- `TIED_SAME_TOPOLOGY`：最佳樹多棵，但 rooted-unlabeled signature 只有一種；
- `TIED_CROSS_TOPOLOGY`：最佳樹多棵，且 signature 超過一種。

另外保存四類 coarse geometry：
`Single-only`、`Sister-only`、`Direct-only`、`Sister+direct`。

## 最終結果

全七 `71,955` 個 `ranked_complete` units、`680,527` 棵 global-best
trees 已精確展開並全部通過 canonical reproduction：

- `UNIQUE_TREE = 39,648/71,955 = 55.1011%`
- `TIED_SAME_TOPOLOGY = 23,858/71,955 = 33.1568%`
- `TIED_CROSS_TOPOLOGY = 8,449/71,955 = 11.7421%`
- 單一 exact rooted-unlabeled topology：
  `63,506/71,955 = 88.2579%`
- 單一四類 coarse class：
  `63,511/71,955 = 88.2649%`

因此舊 legacy grouping 的 `92.18%` 不能直接當成新版 exact-PS×HP
結果；新版同口徑結果是 `88.26%`。

工具逐 unit 強制驗證：

1. minimum-family SHA-256 與 canonical 相同；
2. best AF score 相同；
3. best vertex-set count 相同；
4. factorized tie count、實際展開樹數與 canonical
   `best_tree_tie_count` 三者相同；
5. canonical deterministic representative shape 包含在精確 census 中。

## 分層工作站配套資料

### 完整 minimum candidate factorization

Region detail 不再只依賴 deterministic representative。配套 C++／Python
producer 對 71,955 ranked units 保存：

- 每個 minimum vertex set；
- 每個 child 的所有 exact best parents；
- 完整 minimum candidate edge incidence；
- AF-global-best set／tree incidence；
- exact rooted-shape census。

正式 output：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260725_exact_ps_candidate_factorization/all7_v2/
```

Receipt schema `1.1.0`，SHA-256
`54a8a0a00c1cdb1a40fe96b3a528e9142d21c76af55dc8eae553a6f1432b8164`；
71,955 ranked units、972,592 minimum trees、680,527 global-best trees
全部守恆。它仍繼承上游
`all_mutation_bearing_families_complete=false`，不得替 ABSTAIN 補造候選。

### 跨樣本 profile similarity

正式 sidecar：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260725_exact_ps_cohort_similarity/all7_v1/cohort_similarity.v1.json
```

主頁分別比較 state、determinacy、topology/coarse、active-k 與對稱 HP
balance；各維度在自身分母內正規化，再以
`1 − Jensen–Shannon distance` 表示 composition similarity。

HCC1395 technical pair 的五維 overview 為 `0.926381`（21 pairs rank 1），
strict selected labeled tree 相同為 `443/564 = 78.5%`。前者是 cohort
composition，後者有共同座標與 determinacy gate；兩者都不證明相同
cellular clone 或 biological ancestry。

讀者介面 aliases 為 `HCC1395_HKU` 與 `HCC1395_NYGC`；內部 authority
keys／檔名仍為 `HCC1395` 與 `HCC1395_DORADO`。

### Exact-locus 癌症基因／藥物 annotation

正式 write-once sidecar：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260725_exact_ps_gene_drug_annotation/all7_v2/
```

- schema／receipt `1.1.0`；
- annotation SHA-256
  `148f74e9f958cccffe14a2297e9453e28b2346a424adcb30e350f79e91ae4f50`；
- receipt SHA-256
  `9639b36df7a2412afc2aca97eb1f20819e7a41394394ef40b6f3108e809dd99e`；
- GENCODE v46 GRCh38 gene body 是座標 authority，COSMIC CGC v104
  提供 cancer-gene membership，DGIdb 僅取
  `approved=TRUE AND anti_neoplastic=TRUE`；
- 721 個 autosomal CGC GENCODE genes，279 個有上述 DGIdb
  association；
- DGIdb 224 筆缺 HGNC 的候選列全部未通過 GENCODE exact-symbol gate，
  production CGC overlap key=0。

all7 以 actual sSNV locus 而非 region span 做 same-gene join：
`3,554/98,955` final groups 命中 CGC，`1,252/98,955` 同時命中
CGC 與 approved-antineoplastic association。HCC1395
`chr10:87818272-87928739` 的新 exact unit 只有 87,818,272 /
87,840,023，兩點皆在 PTEN gene body 之外，因此是明確
`NO_HIT_EVALUATED` negative control。

此 sidecar 的 claim ceiling 是 gene-level annotation context；不可轉譯為
topology validation、driver status、藥物反應或 clinical actionability。

## 已執行命令與輸出

Frozen canonical inputs：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/
samples/{SAMPLE}/{SAMPLE}.exact_ps_mlhp.json

/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/
samples/{SAMPLE}/{SAMPLE}.topology.jsonl
```

Build：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
  -I/big7_disk/liaoyoyo2001/LongLineage/include \
  InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/cpp/exact_topology_signature_census.cpp \
  /big7_disk/liaoyoyo2001/LongLineage/src/solver/obligation_bnb.cpp \
  /big7_disk/liaoyoyo2001/LongLineage/src/solver/parent_mapping.cpp \
  -ljansson -lcrypto \
  -o /bip7_disk/liaoyoyo2001/build-exact-topology-signature-20260724/exact_topology_signature_census
```

全七每個樣本的執行參數：

```bash
/bip7_disk/liaoyoyo2001/build-exact-topology-signature-20260724/exact_topology_signature_census \
  --input .../{SAMPLE}/{SAMPLE}.exact_ps_mlhp.json \
  --canonical .../{SAMPLE}/{SAMPLE}.topology.jsonl \
  --output .../all7_v1/{SAMPLE}.census.jsonl \
  --max-family-size 100000 \
  --max-search-nodes 1000
```

實際 output root：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260724_exact_ps_cpp_topology_signature_census/all7_v1/
```

正式 summary 與 receipt：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.json

/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.tsv

/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260724_exact_ps_cpp_topology_signature_census/all7_v1/receipt.v2.json
```

彙總：

```bash
python3 \
  InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/summarize_topology_signature_census.py \
  --input-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1 \
  --output-json /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.json \
  --output-tsv /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.tsv
```

## 配套報告與完整分母

- [AF 最佳樹 topology-signature 精確統計](./20260724_exactPS最佳樹拓撲Signature精確統計_01.md)
- [分層拓撲工作站重建與 Chromium 視覺稽核](./20260724_exactPS分層拓撲工作站重建與視覺稽核_01.md)
- [Claude Code Opus 5 獨立審查](./20260725_ClaudeCode分層工作站獨立審查_01.md)
- [癌症基因 × 藥物顯示 pre-decision audit](./pre-decision-audit.md)
- [exact-PS×HP 的 k、HP 與分母守恆重算](./data/20260724_exactPS_k_HP與分母重算紀錄_01.md)
- [Source component k（區域數）](./data/20260724_exactPS_source_component_k_by_sample_01.tsv)
- [Source component k（sSNV memberships）](./data/20260724_exactPS_source_component_k_memberships_by_sample_01.tsv)
- [Final topology group k（區域數）](./data/20260724_exactPS_final_group_k_by_sample_01.tsv)
- [Final topology group k（sSNV memberships）](./data/20260724_exactPS_final_group_k_memberships_by_sample_01.tsv)
- [HP1／HP2 分流](./data/20260724_exactPS_hp_split_by_sample_01.tsv)
- [分區 funnel](./data/20260724_exactPS_funnel_by_sample_01.tsv)
