<!--
建立時間: 2026-07-24
目標: 驗證 read-pattern likelihood、read-AF/VAF 與 directed mutation-state tree 路徑的可辨識範圍
處理範圍: Task Type B comprehensive validation audit；7 datasets 歷史結果、HCC1395 exact-PS pilot、現行 C++/Python 原始碼
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py
-->

# Likelihood 與 VAF 樹路徑可辨識性驗證

> **TL;DR**：read-pattern likelihood 可以比較不同 mutation-state vertex sets，但不能區分同一 vertex set 內的 parent-edge 路徑。舊 Python 的確曾用 HP-specific read-AF 排序完整候選樹；該方法是祖先順序 heuristic，不是 CCF 或 calibrated likelihood。兩條歷史數據都顯示 k 增加後唯一率下降。現行 C++ 只完成 strict read-linkage 與 exact-PS bounded partition，尚未實作完整候選、BQ likelihood 與 directed topology。

> **Scope / evidence ribbon**：歷史 M0、legacy read-AF 與 current exact-PS 分母不同，只能分開說明；正式新版 strict exact-PS 全 7 datasets likelihood/directed-topology receipt 目前為 0/7。

## 1. 任務與假說

- Task Type：B — comprehensive validation。
- 服務目標：G1、G3、G4、G5。
- H1：BQ-aware read-pattern likelihood 可以選出不同的 candidate vertex set。
- H2：同一 vertex set 的不同 parent-edge assignment 可由同一批 snapshot reads 或其 marginal VAF 分出先後。
- H3：k ≥ 3 因資訊較多，唯一拓撲比例應高於 k=2。
- H4：現行 C++ 已可直接執行 candidate → likelihood → directed topology。

判定：H1 在 candidate family 完整且資料可辨識時成立；H2、H3、H4 不成立。

## 2. Likelihood 真正評分的對象

對候選狀態集合 \(V\)，現行 M2 估計 mixture weights \(\pi\)：

\[
\ell(V)=\max_{\pi}\sum_r n_r\log\left(\sum_{v\in V}\pi_vP(r\mid v,\mathrm{BQ})\right)
\]

- R/A 依每個 call 的 BQ 建立 emission probability。
- X 直接 marginalize，不任意補成 R 或 A。
- 不同 vertex set \(V\) 可以有不同 likelihood。
- parent-edge set \(E\) 沒有進入公式；因此相同 \(V\) 的 \(E_1\) 與 \(E_2\) 必然同分。
- 只有「唯一 winning vertex set，且該 vertex set 只有一個合法 parent assignment」時，才能稱 exact topology unique。
- candidate enumeration incomplete、optimizer 無法取得 KKT certificate 或 winner tie 時必須 abstain。

最小反例：若狀態為 `00, 10, 01, 11`，`11` 可由 `10` 或 `01` 流入。兩棵樹具有相同 states、相同 mixture 與相同 site VAF，因此 snapshot read likelihood 與 marginal VAF 均不能決定真正 parent。

原始碼證據：

- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:17-26`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:79-92`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:838-847`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:1275-1324`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:1414-1436`
- `InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/build_m2_patterns_and_rank.py:2898-2916`

## 3. VAF／read-AF 能否判斷突變先後

在以下前提成立時，祖先突變的細胞盛行率理論上不低於後代，因此 `AF_ancestor ≥ AF_descendant` 可作弱排序證據：

1. infinite-sites／無 recurrence；
2. diploid copy-neutral；
3. 無 LOH、WGD、allele-specific CN 或 mutation multiplicity 差異；
4. purity 與取樣偏差已控制；
5. coverage 足夠且 mapping／allelic bias 可忽略。

但 raw VAF/read-AF 會受 CN、LOH、purity、multiplicity、深度與 allele bias 影響；相同 AF 也不代表 parent-child。結論是：

- 可作 prespecified compatibility check、排序 sensitivity 或診斷 sidecar。
- 若有 CN/purity/multiplicity，應改用帶不確定區間的 CCF model。
- 不宜把同一批 reads 算出的 VAF score 再加進 M2 likelihood，否則重複使用 allele-count 資訊。
- 即使 VAF 完全準確，通常也只能限制 prevalence order，不能單獨識別同一 vertex set 內的 edge。

## 4. 舊 Python 實際做法

舊版不是現行 M2 likelihood；它用每個 HP family 的：

\[
\mathrm{readAF}_i=\frac{ALT_i}{REF_i+ALT_i}
\]

對候選樹每條 `parent → child` 新增 mutation \(j\)，累加 parent 已有 mutation \(i\) 的：

\[
\mathrm{Score}(T)=\sum_{p\to c}\sum_{i\in mutations(p)}
(\mathrm{readAF}_i-\mathrm{readAF}_j)
\]

hidden `H_*` state 會先還原成 mutation set，因此也參與 edge score。它確實可以替完整候選樹排序，但：

- 是 heuristic，不是 probability、CCF 或 calibrated likelihood；
- 不同候選的 ancestry comparison 數可能不同；
- softmax 只是 score 轉換，不能稱 posterior；
- 使用相同 read-AF 再檢查方向一致性，不是獨立驗證。

原始碼證據：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/read_af_tree_ordering_multisample.py:13-56`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/read_af_tree_ordering_multisample.py:172-187`
- `InterSubMod/research/20260711_read_group_C_tree_T_topology_report/scripts/build_vaf_top_tie_census.py:100-155`
- `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py:416-537`
- `InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/scripts/build_current_v5_read_af_topology.py:854-879`

## 5. 數據驗證：k ≥ 3 並沒有更容易唯一

### 5.1 Historical M0 likelihood

單位為 legacy eligible HP lineage；使用舊 mixed-PS／50-kb regions 與 thresholded alignment exposure，不是 current strict W。只看真正有多候選的 T>1 units：

| k | T>1 分母 | exact single-edge winner | vertex-set tie | edge unresolved | optimizer abstain |
|---:|---:|---:|---:|---:|---:|
| 2 | 10,766 | 3,904（36.26%） | 6,775（62.93%） | 67（0.62%） | 20（0.19%） |
| ≥3 | 27,982 | 4,833（17.27%） | 20,495（73.24%） | 48（0.17%） | 2,606（9.31%） |

exact single-edge winner 由 k=3 到 k=8 的比例依序為：

- k=3：2,216 / 9,144 = 24.23%
- k=4：1,155 / 6,377 = 18.11%
- k=5：552 / 3,426 = 16.11%
- k=6：326 / 2,202 = 14.80%
- k=7：180 / 1,357 = 13.26%
- k=8：404 / 5,476 = 7.38%

輸入：

`InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz`

### 5.2 Historical exact read-AF ordering

單位為 legacy layered-v2 可評估 ambiguous region；同樣不是 current strict exact-PS W：

| k | unique first / evaluable | 唯一率 |
|---:|---:|---:|
| 2 | 8,233 / 10,328 | 79.72% |
| 3 | 6,117 / 8,379 | 73.00% |
| 4 | 3,077 / 4,885 | 62.99% |
| 5 | 960 / 1,967 | 48.81% |
| 6 | 341 / 953 | 35.78% |
| 7 | 116 / 416 | 27.88% |
| 8 | 248 / 1,255 | 19.76% |
| 合計 | 19,092 / 28,183 | 67.74% |

合計另有 tied-first same topology 6,886、tied-first different topology 2,205。這裡的 topology 是移除 mutation labels 後的 rooted shape，不是 exact labeled parent-child tree。

輸入：

`InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/vaf_top_tie_census.json`

### 5.3 Current exact-PS 狀態

- 全 7 datasets strict read-linkage W 已建立，共 85,621；k=2 為 40,599，k≥3 為 45,022。
- 但 7 datasets 的 production likelihood/directed topology 均未執行完成；winner/tie/edge 結果是 N/A，不是 0。
- HCC1395 partial pilot 的 9,600 mutation-bearing route-blocks中，complete 9,180；`T=1|Topo=1` 4,579、`T>1|Topo=1` 2,208、`T>1|Topo>1` 2,393、incomplete 420。
- 該 HCC1395 pilot 明示 `VAF ranking = 未執行`，不能當新版排序結果。

來源：

- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/exact_k_distribution.tsv`
- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/data/topology_completion_status.tsv`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/20260723_HCC1395_exactPS拓撲重建觀察_01.md:35-69`

## 6. 為何 k 越大反而更難唯一

k 增加同時帶來更多 read information，也使：

- hypercube vertex universe 以 \(2^k\) 成長；
- partial-pattern compatible states 增加；
- minimum-extra vertex sets 可能成為大型 family；
- 每個 vertex set 內可能有多個合法 parent assignments；
- hidden states 與不同 edge paths 更常 observationally equivalent。

所以不是「k 大必然更好」，而是新增 informative full/partial patterns 的速度必須超過候選空間成長，才可能提高辨識度。現有歷史數據顯示沒有發生這件事。

## 7. C++ 能力判定

現有 C++：

- 已有 `dataset × chromosome × exact PS × HP` strict endpoint graph；
- 已有 read-linked components 與 k≤12 bounded partition；
- `Observation` 只有 node 與 call code，沒有 per-call BQ；
- 沒有 symbolic group-Steiner candidate solver、minimum-hidden certificate、BQ mixture/KKT、winner tie、parent-edge/topology receipt。

證據：

- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/cpp/strict_endpoint_graph.hpp:30-135`
- `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/cpp/strict_endpoint_graph.hpp:234-360`
- `InterSubMod/CMakeLists.txt:133-138`

因此技術上可以移植，但現況不能直接產出可信的 full result。

建議順序：

1. 先凍結 Python oracle 語意與 deterministic golden corpus  
   → 驗證：candidate family、minimum hidden、winner/tie、edge count、abstain 完全可重現。
2. 先移植 C++ candidate backend，Python 保留 likelihood/ranking  
   → 驗證：tractable fixtures 與 hard-tail cases 的 objective、family count、candidate digest 完全相同。
3. 再移植 BQ mixture/KKT  
   → 驗證：LL、fitted observable probabilities、KKT gap 與 winner partition 在容許誤差內一致。
4. 最後移植 topology/receipt，跑 HCC1395_DORADO chr6 與 H2009 chr2  
   → 驗證：terminal receipt、完整/abstain 分母、runtime/RSS 全部 PASS，再進 7×22。

這是 Python oracle-first、C++ backend-first；不是要求全量永遠由 Python 執行。

## 8. 本次唯讀驗證命令與輸出片段

執行根目錄：

`/big7_disk/liaoyoyo2001/InterSubMod`

```bash
python - <<'PY'
import gzip, csv, collections
p = "research/20260716_read_linked_hypercube_exact_likelihood_validation/results/m0_full_v4/m0_hp_lineage_likelihood_census.tsv.gz"
with gzip.open(p, "rt") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
for name, subset in [
    ("k=2", [r for r in rows if int(r["k"]) == 2 and int(r["raw_tree_candidates_T"]) > 1]),
    ("k>=3", [r for r in rows if int(r["k"]) >= 3 and int(r["raw_tree_candidates_T"]) > 1]),
]:
    print(name, len(subset), collections.Counter(r["selection_status"] for r in subset))
PY
```

實際輸出摘要：

```text
k=2 10766:
  unique_vertex_single_edge=3904, tied_vertex_sets=6775,
  edge_unresolved=67, optimizer_abstain=20
k>=3 27982:
  unique_vertex_single_edge=4833, tied_vertex_sets=20495,
  edge_unresolved=48, optimizer_abstain=2606
```

C++ 稽核的既有測試結果：

```text
strict endpoint graph golden diff: PASS
exact-PS C++/Python partition focused tests: 23 PASS
Python candidate + BQ likelihood oracle tests: 61 PASS
all-7 production strict directed topology: 0/7
```

## 9. 最終 claim ceiling

目前可說：

> 我們能枚舉與 read constraints 相容的 minimum-extra mutation-state candidate sets，並用 BQ-aware read-pattern likelihood 判斷不同 state sets 是否可由資料區分。舊 read-AF 可提供條件式祖先順序 heuristic。

目前不可說：

> 已由 bulk snapshot reads／VAF 唯一重建所有 hidden states 的真實 parent-child tree、真實 clone 數、subclone lineage 或全樣本完整演化史。

