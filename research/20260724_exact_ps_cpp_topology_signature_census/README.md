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
- [exact-PS×HP 的 k、HP 與分母守恆重算](./data/20260724_exactPS_k_HP與分母重算紀錄_01.md)
- [Source component k（區域數）](./data/20260724_exactPS_source_component_k_by_sample_01.tsv)
- [Source component k（sSNV memberships）](./data/20260724_exactPS_source_component_k_memberships_by_sample_01.tsv)
- [Final topology group k（區域數）](./data/20260724_exactPS_final_group_k_by_sample_01.tsv)
- [Final topology group k（sSNV memberships）](./data/20260724_exactPS_final_group_k_memberships_by_sample_01.tsv)
- [HP1／HP2 分流](./data/20260724_exactPS_hp_split_by_sample_01.tsv)
- [分區 funnel](./data/20260724_exactPS_funnel_by_sample_01.tsv)
