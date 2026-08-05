<!--
建立時間: 2026-07-24
目標: 記錄 HCC1395 exact parity、C++ 加速與七樣本 chr1-22 全流程結果
處理範圍: Task Type B / comprehensive validation / 7/22 HCC1395 oracle + 7/23 all-seven strict production inputs
關聯檔案: InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/implementation-notes.md
-->

# C++ exact-PS 拓撲與 read-AF 全七樣本驗證

## Assertion–Evidence 結論

**HCC1395 gate 通過，但全 cohort 尚未 topology-complete。**

- HCC1395 舊 Python complete units：11,122/11,122 逐 unit exact parity，
  active k、最少 hidden-node 數、完整 minimum vertex-set family digest、
  family 數與總 parent-tree 數皆 0 mismatch。
- HCC1395 read-AF oracle：4,601 units、exact `Fraction` 分數與最佳樹並列數
  0 mismatch。
- HCC1395 controlled wall time：Python 101.54 s，C++ 6.73 s，工程吞吐
  約 15.09 倍。
- 七樣本 chr1–22 已全部執行且 technical PASS；但 85,941 個
  mutation-bearing units 中有 10,717 個命中 search guard 而 fail-closed
  `ABSTAIN_RESOURCE_LIMIT`，因此不可宣稱每個區域都已解出完整候選樹族。

## AF 分數為何可以不展開每棵樹

固定一個可行 vertex set \(V\)。對每個非根 child \(c\)，其合法父節點集合為
\(P(c)\)，而 edge score 為：

\[
s(p,c)=\sum_{i\in p}\left(AF_i-AF_{j(c)}\right),
\qquad
AF_j=\frac{ALT_j}{REF_j+ALT_j}.
\]

一棵 parent tree 的總分為：

\[
S(T)=\sum_{c\ne r}s(p_c,c).
\]

在目前的 upward Hamming-1、單一根、recurrence-allowed 模型下，每個 child
選哪個合法 parent 不會改變其他 child 的合法 parent 集合。因此：

\[
\max_T S(T)
=
\sum_{c\ne r}\max_{p\in P(c)}s(p,c),
\]

\[
\#\operatorname*{argmax}_T S(T)
=
\prod_{c\ne r}
\left|\operatorname*{argmax}_{p\in P(c)}s(p,c)\right|.
\]

若多個 minimum vertex sets 得到相同全域最佳分數，再把各 vertex set 的
最佳樹數相加。C++ 因此可直接得到 exact 最佳分數、最佳樹並列數與唯一性，
只保存一棵 deterministic representative，不必 materialize 全部 parent trees。

這個分解依賴目前模型；加入全域 infinite-sites、homoplasy penalty、
最大 out-degree、CN mixture 或其他跨 edge coupling 後，不能直接沿用。

## HCC1395 7/22 oracle gate

### 輸入

- MLHP：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json`
- SHA256：
  `68bddbc6e6080194f031fb6d58b131933036614d49f5063a63dcdc7703d15ee6`
- Historical layered Python output：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/layered_reconstruction_HCC1395_chr1_22_display_all.json`
- SHA256：
  `c5fb1457f524c29ff1171358873511465911dd31efe50c8830083f77dd7e85cf`

### C++ 執行命令

```bash
taskset -c 40 \
  /bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af \
  --input /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_exact_ps_topology_hcc1395_v1/full_chr1_22_v1/exact_ps_mlhp_HCC1395_chr1_22.json \
  --output /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/hcc1395_7_22_pilot_gate_guard1000_v3/HCC1395.topology.jsonl \
  --receipt /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/hcc1395_7_22_pilot_gate_guard1000_v3/HCC1395.topology.receipt.json \
  --max-family-size 100000 \
  --max-search-nodes 1000
```

### 實際輸出片段

```text
sample=HCC1395 units=11542 objective_certified=11268
mutation_family_complete=9326/9600 ranked=9084
elapsed_seconds=5.983738
```

### Gate 結果

| 檢查 | 結果 |
|---|---:|
| Contract checks | 23/23 PASS |
| 舊 Python complete units | 11,122/11,122，0 mismatch |
| read-AF oracle units | 4,601，0 mismatch |
| 舊 Python capped units | 420 |
| C++ 對舊 capped units補解 | 146 |
| C++ 對舊 capped units保守棄權 | 274 |
| invalid outcomes | 0 |

Gate receipt：
`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/hcc1395_7_22_pilot_gate_guard1000_v3/HCC1395.cpp_gate.receipt.json`

## 七樣本 chr1–22 全流程

七樣本 production run 使用 7/23 strict exact-PS membership；它與 7/22
HCC1395 historical pilot 的 singleton/filter denominator 不同，因此只以前者作
cohort 統計，不混用兩套分母。

### 輸入 manifest

`InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/input_contract/all7_exact_ps_inputs.local.json`

### 執行命令

```bash
PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -v \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/scripts/run_exact_ps_cpp_topology_all7.py \
  --output-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1 \
  --topology-binary /bip7_disk/liaoyoyo2001/build-exact-af-20260724/bin/exact_ps_topology_af \
  --max-family-size 100000 \
  --max-search-nodes 1000 \
  --workers 4 \
  --resume
```

### 每樣本結果

`complete %` 的分母是 mutation-bearing units；`unique/tied` 只計 AF 可排序者。

| 樣本 | groups | mutation | complete % | ABSTAIN | ranked | unique | tied | C++ s | adapter s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 11,590 | 9,624 | 97.23% | 267 | 9,130 | 7,047 | 2,083 | 6.72 | 96.00 |
| HCC1395_DORADO | 6,865 | 5,485 | 98.67% | 73 | 5,308 | 4,165 | 1,143 | 2.55 | 89.57 |
| COLO829 | 13,919 | 11,395 | 98.72% | 146 | 10,757 | 5,119 | 5,638 | 5.03 | 33.63 |
| H1437 | 17,598 | 15,661 | 89.87% | 1,586 | 13,740 | 6,369 | 7,371 | 25.25 | 96.41 |
| H2009 | 36,042 | 33,548 | 74.79% | 8,457 | 23,128 | 8,161 | 14,967 | 116.37 | 280.57 |
| HCC1937 | 5,195 | 4,400 | 97.61% | 105 | 4,245 | 3,672 | 573 | 2.93 | 54.24 |
| HCC1954 | 7,746 | 5,828 | 98.58% | 83 | 5,647 | 5,115 | 532 | 2.91 | 45.82 |
| **合計** | **98,955** | **85,941** | **87.53%** | **10,717** | **71,955** | **39,648** | **32,307** | **161.76** | **696.23** |

守恆：

- mutation family：75,224 complete + 10,717 abstain = 85,941。
- complete mutation families：
  71,955 ranked + 3,224 zero denominator + 45 recurrence-required = 75,224。
- ranked：
  39,648 unique + 32,307 tied = 71,955。

### 唯一 AF winner 的 graph geometry

只有 `best_tree_unique=true` 才把 morphology 當作最終分類；並列者不使用
deterministic representative 冒充已解拓撲。

| 分類 | 數量 | 占 unique winners |
|---|---:|---:|
| Single-only | 22,135 | 55.83% |
| Direct-only | 13,176 | 33.23% |
| Sister-only | 2,219 | 5.60% |
| Sister＋direct | 2,118 | 5.34% |
| AF 並列、shape 未解 | 32,307 | 不納入分母 |

這四類是 graph geometry，不直接等同於 cell clone/subclone；read 是 molecule，
而 CNA、LOH、allele-specific amplification、purity 與 sampling 仍可能改變
cell-level 解釋。

## 實際時間與資源

- 成功 resume wall time：1,191.49 s = 19 m 51.49 s。
- 初次啟動在 legacy `bytes`/current `size_bytes` receipt 欄位不相容後
  fail-closed，保留 45.29 s 診斷紀錄；修正後 resume。
- 從首次空目錄啟動到完成、包含該診斷失敗：20 m 36.78 s。
- 七樣本 C++ topology stage 合計：161.76 s。
- MLHP adapter stage 合計：696.23 s。
- 成功 run 最大 RSS：1,685,008 KiB。
- Summary verification：10.06 s，最大 RSS 23,124 KiB。

## Guard 結果與限制

`k<=12` 只把單一 unit 的 vertex universe 限在最多 \(2^{12}=4096\)；
它不會限制同分 optimal vertex sets 的數量。guard=1,000 的實際 abstain
分布為：

| active k | units | ABSTAIN | ABSTAIN % |
|---:|---:|---:|---:|
| 0–4 | 76,932 | 0 | 0% |
| 5 | 5,398 | 780 | 14.45% |
| 6 | 3,771 | 3,440 | 91.22% |
| 7 | 2,548 | 2,514 | 98.67% |
| 8 | 1,596 | 1,593 | 99.81% |
| 9 | 925 | 925 | 100% |
| 10 | 590 | 590 | 100% |
| 11 | 476 | 476 | 100% |
| 12 | 399 | 399 | 100% |

因此本輪證明：

1. C++ exact 語意與舊 Python 在可核對範圍一致。
2. Compact AF ranking 可以正確避免 parent-tree 全展開。
3. 全七樣本技術流程可在約 21 分鐘內 fail-closed 跑完。
4. 目前 guard 與 exact vertex-set solver 尚不足以把 k≥6 長尾全部解完；
   尤其 H2009 不能列為 topology-complete。

## 驗證命令與結果

```bash
PYTHONDONTWRITEBYTECODE=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q -p no:cacheprovider \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/tests \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/oracles/test_validate_hcc1395_cpp_topology_af_gate.py \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_partition_to_mlhp.py
```

實際輸出：

```text
28 passed in 20.55s
```

另外，cohort receipt SHA256、所有 source identities、列數與八項算術守恆
皆 PASS。

## 交付位置

- Cohort receipt：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/cohort_receipt.json`
- 完整 JSON summary：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/summary/all7_summary.json`
- 每樣本 TSV：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/summary/all7_summary.tsv`
- C++ source：
  `InterSubMod/research/20260724_exact_ps_cpp_topology_af_all_samples/cpp/exact_ps_topology_af.cpp`
- C++ binary SHA256：
  `ba13ccf23d091854c191f81dd97fa891368d11179df9a69e915f0340b7233b2e`

## 尚未完成的功能

Production compact path 已實作最佳分數、最佳樹並列數、唯一性與一棵代表樹；
**尚未實作「按需輸出全部完整 tied trees」的 CLI 選項**。這不影響本輪
compact ranking 的 exactness，但若要逐棵檢視所有並列 shape，需另加受 output
budget 控制的 enumerator。
