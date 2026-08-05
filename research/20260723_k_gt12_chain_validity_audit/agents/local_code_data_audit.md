<!--
建立時間: 2026-07-23 12:32 +08:00
目標: 只讀核對 exact-PS strict read-linkage、k>12 partition 與 A-B/B-C chaining 的程式語意及實際 HCC1395 數據
處理範圍: 本地程式、HCC1395 chr1-22 receipts、七 technical datasets 即時完成狀態；不修改程式、不宣稱生物 truth
關聯檔案:
  - InterSubMod/scripts/build_strict_ps_hp_regions.py
  - InterSubMod/tools/strict_endpoint_graph.py
  - InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/scripts/exact_ps_k12_partition.py
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md
-->

# k>12 與 A–B／B–C chaining：本地程式及數據稽核

> **任務類型**：B（comprehensive validation；但七資料集 production 仍在執行，跨七資料集部分帶 `PARTIAL` 標記）  
> **服務目標**：G1、G4、G5  
> **框架**：定義 → 程式實作 → receipt 重算 → claim ceiling  
> **TL;DR**：`k>12` 應同時保留成一個完整 `source region W`，並在現行運算路徑切成 `k≤12 computational blocks B`；不能把 block 當成新的生物區域。A–B、B–C 的 transitive chaining 對「定義 read-linked connected component」合理，但只代表連通，不代表 A–C 直接共現、同一 molecule、同一 cell 或同一 clone。現行 block DP 也未強制每個 child block 在 threshold-3 graph 內仍連通；HCC1395 實測有 4 個 `k≥2` child blocks 不連通，另有 134 個 `k≥2` child blocks 沒有 retained constraint weight，故下游必須另設 tree-eligibility gate。

## 1. 本次問題、假設與驗收

### 1.1 研究問題

1. 程式與數據中的 `k` 各指什麼？
2. `k>12` 是切掉原始 region，還是只切 computational blocks？
3. A–B、B–C 沒有 A–C direct read 時，串成同一個 W 是否符合程式與理論語意？
4. partition 後的 block 是否仍保留 graph connectivity 與 read constraints？
5. HCC1395 及七 technical datasets 目前實際完成到哪裡？

### 1.2 關鍵假設

- `W` 僅宣稱 read-linkage graph 的 connected analysis domain，不宣稱 biological clone。
- `B` 是計算上界所產生的 child block，不得改寫成新 region。
- `support_total≥3` 是工程 threshold；尚未證明是最優統計閾值。
- `RR`、`RA`、`AR`、`AA` 目前全都計入 edge 的 `support_total`。

### 1.3 成功／失敗條件

- 成功：程式碼、receipt 與獨立重算對 `W/B/k/constraint disposition` 一致。
- 失敗：將 block 數當 region 數；將 graph transitivity 當 allele/clone transitivity；將進行中的七資料集輸出當作 7/7 完成。

## 2. `k` 的四種相關定義

| 名稱 | 正確定義 | 不代表什麼 |
|---|---|---|
| `k_W` 或 source `k` | `k_W=|W|`：一個 threshold-qualified maximal connected component 中的 sSNV site 數 | read 數、clone 數、cell 數、候選樹數 |
| `k_B` 或 block `k` | `k_B=|B|`：一個 computational block 中的 site 數；現行強制 `k_B≤12` | 新的 biological region |
| nominal `k` | Boolean state vector 的總 bit 數；完整 state space 為 `2^k` | 實際 solver 一定使用的有效維度 |
| active `m`（active-k） | observed-alt contract 下，至少在一個 pattern 明示 ALT 的 bit 數；部分 solver 的顯式 state space 為 `2^m` | source region 的總 site 數 |

另有兩個容易混淆的數：

- 一條 pattern 實際 fixed-observe 的 endpoint 數是 `|sites(c)|`，不是 source `k`。
- `span_sites = hi−lo+1` 是 constraint 左右端點在 **W 的排序索引**間跨過的 site 數；即使 read 只 fixed-call 兩個 endpoint，中間若夾了很多 W sites，`span_sites` 仍可大於 12。

所以 `k>12` 應稱為「large multi-site read-linked component」。若要稱「基因體密集區域」，還必須另外用 `k/span_bp` 或相鄰 gap 分布定義；`k>12` 本身不是 genomic density。

## 3. `k>12`：保留 W，切 B

### 3.1 Canonical 語意

```text
source region W_i（完整保留；component_id 不變）
       │
       ├─ k_W ≤ 12 → 一個 B_i1
       │
       └─ k_W > 12 → ordered-hypergraph DP → 多個 k_B≤12 的 B_ij
```

因此答案不是「切開」或「當成特殊密集區域」二選一，而是兩層並存：

1. **研究與報告層**：仍是一個完整 W，標記 `large_region / k>12`。
2. **現行 solver 運算層**：切成數個 bounded blocks。
3. **生物解讀層**：只能報 local block topology；尚不能把多個 local trees 無損接回唯一 global topology。

### 3.2 現行 DP 的目標

`read_support_partition.py` 把 W 中依 genomic position 排序的 sites 切成不重疊、連續、完整覆蓋 W 的 blocks，依序最佳化：

1. retained molecule weight 最大；
2. retained pattern count 最大；
3. block 數最少；
4. cut 所落 genomic gaps 總和最大；
5. cut-index tuple lexicographically smaller。

這是 **constraint-preserving bounded partition**，不是 biological boundary detection，也沒有把「block graph 必須連通」寫進最佳化限制。

### 3.3 HCC1395 實際切分

| 單位 | 數量 |
|---|---:|
| source W | 11,462 |
| source W memberships | 34,267 |
| `k_W>12` source W | 90 |
| `max(k_W)` | 153 |
| 小 W（`2≤k_W≤12`） | 11,372 |
| 全部 computational B | 11,712 |
| 由 90 個 large W 產生的 B | 340 |
| large-W child B 中 retained weight >0 | 194 |
| large-W child B 中 retained weight =0 | 146 |
| 其中 `k_B=1` 且 retained weight=0 | 12 |
| 其中 `k_B≥2` 且 retained weight=0 | 134 |

因此 11,712 不能作為 region 分母，更不能直接當 topology 分母。

## 4. A–B、B–C chaining 的精確語意

### 4.1 程式實作

在同一個：

```text
dataset × chromosome × exact nonmissing PS × primary HP family
```

容器內：

1. 一個 canonical molecule 只有在兩 endpoint 都有 fixed `R/A` call 時，才對該 pair 投一票。
2. 每個 molecule 對同一 pair 最多一票。
3. edge 的 `support_total = RR+RA+AR+AA`。
4. primary edge 要求 `support_total≥3`。
5. `connected_components()` 對合格 edge 做 union-find。

所以：

```text
A—B  support≥3
B—C  support≥3
A—C  support=0
```

仍會得到 `{A,B,C}` 一個 W。A–B 與 B–C 可由完全不同的 molecule sets 支持。

### 4.2 哪個推論合理、哪個不合理

| 主張 | 判定 |
|---|---|
| A、B、C 屬於同一個 maximal read-linked connected component | 合理，也是現行程式的定義 |
| A、B、C 可放進同一個候選 mutation-state 問題 | 合理，但需保留 partial observations 與多解 |
| 已直接觀察 A–C 共現 | 不成立 |
| 有一條 molecule 橫跨 A、B、C | 不成立 |
| A、B、C 一定屬於同一 cell／clone | 不成立 |
| pairwise allele pattern 可經圖的傳遞性直接當成完整 haplotype | 不成立 |

最安全口徑：

> A–B 與 B–C 的 threshold-qualified direct endpoint evidence，允許把三點納入同一個 read-linked connected analysis domain；但完整 R/A state、A–C 共現與演化順序仍需由 constraint model 列舉，且可能不可辨識。

### 4.3 block 層的新發現

我用 strict threshold-3 endpoint edges，對全部 output blocks 做 induced-subgraph connectivity 重算：

| block 類別 | 數量 |
|---|---:|
| 全部 B | 11,712 |
| `k_B=1` | 12 |
| `k_B≥2` | 11,700 |
| `k_B≥2` 且 induced graph connected | 11,696 |
| `k_B≥2` 但 induced graph disconnected | **4** |
| disconnected 中完全沒有 threshold-3 internal edge | 1 |

4 個例外全部來自 large W 的 child blocks；其 retained weight 都是 0。這證明：

> W 的 graph connectivity 不會自動被「按 site order 切 contiguous blocks」保留。

目前應補的 gate：

```text
block_tree_eligible =
    k_B >= 2
    AND induced_threshold3_graph_connected
    AND retained_pattern_count > 0
    AND downstream_MINREAD_pass
```

在 gate 未正式落地前，不能把所有 11,712 blocks 都稱為 tree-ready。

### 4.4 四個 disconnected `k_B≥2` blocks 的精確證據

下表以每個 block 的 `linkage_basis × phase_set × positions` 回查 strict endpoint-edge TSV，只保留 `support_total≥3` 的 internal edges，再對 block induced subgraph 重新求 connected components。

| chrom | source `k_W` | block ID | `k_B` | retained weight／patterns | induced threshold-3 edges | induced component sizes |
|---|---:|---|---:|---:|---:|---|
| chr6 | 74 | `U12b049833dea3ecf49966f5e:B0007` | 8 | 0／0 | 7 | 5, 2, 1 |
| chr6 | 16 | `U500892875397fd2b1986e107:B0001` | 4 | 0／0 | 3 | 3, 1 |
| chr6 | 73 | `Uea95c3bd7e2a91b6fdd3fb32:B0007` | 5 | 0／0 | 3 | 3, 2 |
| chr8 | 14 | `Uf209e406cd66a00dfb91598b:B0001` | 2 | 0／0 | 0 | 1, 1 |

對應 source W：

```text
chr6  HCC1395:chr6:PS_HP2:PSb182726712d9:E3:7:837998-898494:e70a52c058e30cb4
chr6  HCC1395:chr6:PS_HP2:PS273c703a9279:E3:14:3062810-3111722:24603196792e782a
chr6  HCC1395:chr6:PS_HP1:PSb7b164e73255:E3:24:223039-329869:0b5fccd349a1fcb0
chr8  HCC1395:chr8:PS_HP2:PS61df106b75a8:E3:6:121832214-121895076:1440c1a7cd908ca7
```

其 internal-edge 明細：

```text
U12b...:B0007
  887965-888163 support=3 RR/RA/AR/AA=0/0/0/3
  889966-889985 support=3 RR/RA/AR/AA=0/0/0/3
  889966-890000 support=3 RR/RA/AR/AA=0/0/0/3
  889985-890000 support=3 RR/RA/AR/AA=0/0/0/3
  889985-890942 support=3 RR/RA/AR/AA=0/0/0/3
  889985-898494 support=3 RR/RA/AR/AA=0/0/0/3
  890942-898494 support=3 RR/RA/AR/AA=0/0/0/3

U500...:B0001
  3062810-3064749 support=3 RR/RA/AR/AA=0/0/0/3
  3062810-3066582 support=3 RR/RA/AR/AA=0/0/1/2
  3064749-3066582 support=3 RR/RA/AR/AA=0/0/1/2

Uea95...:B0007
  319519-320201 support=3 RR/RA/AR/AA=0/0/0/3
  319519-321990 support=3 RR/RA/AR/AA=1/0/0/2
  325357-329869 support=3 RR/RA/AR/AA=2/0/0/1

Uf209...:B0001
  no internal edge with support_total>=3
```

這裡「有 internal pair edges，但 retained full-pattern weight=0」並不矛盾：edge graph 把一條 molecule 的所有 fixed endpoint pairs 分開累計；partition constraint 則保留該 molecule 在 source W 內形成的**完整 aggregated pattern**。如果該完整 pattern 還含 block 外 site，它會被記為 cut／unavoidable，而不會因 block 內存在一對 endpoints 就拆成新的 retained constraint。

### 4.5 134 個 zero-weight `k_B≥2` blocks 的來源摘要

獨立重算確認：

- 數量：134。
- 來源：58 個不同 source W。
- 58/58 source W 全部為 `k_W>12`；小型 W 沒有此問題。
- 134/328（40.8537%）個 large-W child blocks 在 `k_B≥2` 下 retained weight=0。
- 這不等於 134 個 block 都不連通：只有前述 4 個不連通；其餘可有 threshold-3 pair edges，但沒有完整 retained pattern。

依染色體：

| chromosome | zero-weight `k_B≥2` blocks | 比例 |
|---|---:|---:|
| chr2 | 1 | 0.7463% |
| chr4 | 1 | 0.7463% |
| chr6 | 56 | 41.7910% |
| chr8 | 26 | 19.4030% |
| chr14 | 2 | 1.4925% |
| chr16 | 44 | 32.8358% |
| chr17 | 1 | 0.7463% |
| chr22 | 3 | 2.2388% |
| **total** | **134** | **100%** |

chr6＋chr16 合計 100/134（74.6269%）；加上 chr8 後為 126/134（94.0299%）。因此 zero-weight blocks 高度集中，需和 HCC1395 upstream chr6／chr16 QC 異常一起看，不能先解讀成一般生物現象。

依 child-block `k_B`：

| `k_B` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| count | 5 | 11 | 7 | 7 | 9 | 9 | 12 | 13 | 9 | 19 | 33 |

zero-weight 數最多的 source W：

| chrom | unit ID | `k_W` | zero-weight child B |
|---|---|---:|---:|
| chr6 | `U12b049833dea3ecf49966f5e` | 74 | 6 |
| chr6 | `U2e40e939dfd78acd6d03664a` | 80 | 6 |
| chr6 | `U155b145c7aa60f0f34b5b0f7` | 79 | 5 |
| chr6 | `Ua138bc87ea0a8d47c2bce00e` | 93 | 5 |
| chr16 | `U17d37148e1acadc60d0509b8` | 65 | 5 |
| chr16 | `U8195381746c8cb28eee31488` | 58 | 5 |
| chr16 | `Uae1dfc65659d24af68a00c25` | 77 | 5 |

完整 134-block 名單可由 §8.5 的命令直接輸出；本表用 source unit、chromosome 與 block-k 三條守恆路徑摘要其來源。

## 5. RR-only edge 是 chaining 的重要未知

現行 primary graph 使用 `support_total=RR+RA+AR+AA`。HCC1395 76,202 條 qualifying edges 中：

| edge 類別 | 數量 | 比例 |
|---|---:|---:|
| RR-only | 18,116 | 23.7737% |
| 至少一筆 RA／AR／AA | 58,086 | 76.2263% |

RR-only edge 證明同一 molecule 對兩位點都可 callable 且為 REF，但不證明兩個 somatic ALT 共現。這不妨礙它作為「可共同觀測」訊號，但是否足以定義 mutation-evolution W，仍是方法學待確認。

獨立反事實重算（**exploratory，不是 canonical 結論**）：

| counterfactual edge rule | 原 W 保持連通 | W partition 改變 | 完全失去 multisite ALT-informative linkage |
|---|---:|---:|---:|
| total support≥3，且至少 1 個 ALT-informative molecule | 9,025/11,462（78.7384%） | 2,437（21.2616%） | 1,814（15.8262%） |
| ALT-informative support 自身≥3 | 7,929/11,462（69.1764%） | 3,533（30.8236%） | 2,392（20.8690%） |

這不是說 RR 必須移除，而是說 W 對 edge-state policy 並不完全穩健；正式論文至少要報 RR-inclusive／ALT-informative sensitivity，或把 RR-only W 降為 callability-linked 類別。

## 6. Constraint retained／cut／unavoidable 分母重算

### 6.1 兩種分母必須分開

1. **row denominator**：57,629 個 aggregated unique pattern constraints；每列 `pattern_count=1`。
2. **molecule-weight denominator**：285,596 個 molecule–unit constraint incidences；不是全基因組 unique molecule 數。

### 6.2 獨立重算

| disposition | constraint rows | row % | molecule weight | weight % |
|---|---:|---:|---:|---:|
| retained | 53,918 | 93.5605% | 281,685 | 98.6306% |
| cut | 1,078 | 1.8706% | 1,242 | 0.4349% |
| unavoidable | 2,633 | 4.5689% | 2,669 | 0.9345% |
| **total** | **57,629** | **100%** | **285,596** | **100%** |

定義：

- `retained`：constraint 的所有 sites 落在同一 B。
- `cut`：constraint 理論上可放進某個 `≤12` block（`span_sites≤12`），但被選定 cut 分開。
- `unavoidable_span_gt_max_block_size`：`span_sites>12`，任何不重疊 contiguous `≤12` partition 都不可能完整保留。

注意：unavoidable 不等於 sequencing error，也不等於 biological conflict；它是目前 block-size 與 contiguous-disjoint partition contract 造成的結構性不可保留。

## 7. HCC1395 分布與七資料集即時狀態

### 7.1 HCC1395 strict W（chr1–22，all-pass）

| 指標 | 數值 |
|---|---:|
| candidate loci S | 79,687 |
| active unique loci | 36,384 |
| PS×HP node memberships | 62,651 |
| all components（含 singleton） | 39,846 |
| singleton abstain | 28,384 |
| tree-eligible W | 11,462 |
| `k_W>12` | 90（0.7852% of W） |
| max `k_W` | 153 |

`k_W=2/3/4/5` 分別為 6,840／2,578／1,041／455，占 W 的大多數。大型 W 是少數，但可能造成最多 computational ambiguity，故不應因比例小就忽略。

### 7.2 七 technical datasets（2026-07-23 12:32 +08:00 snapshot）

| dataset | extraction | strict chr1–22 summary | W | `k>12` | max k |
|---|---:|---|---:|---:|---:|
| HCC1395 | 22/22（既有正式輸出） | PASS | 11,462 | 90 | 153 |
| HCC1395_DORADO | 22/22 | PASS | 6,828 | 34 | 106 |
| COLO829 | 22/22 | PASS | 13,933 | 30 | 40 |
| H1437 | 22/22 | PASS | 16,326 | 904 | 138 |
| H2009 | 5/22 | 尚無 summary | — | — | — |
| HCC1937 | 0/22 | 尚無 summary | — | — | — |
| HCC1954 | 0/22 | 尚無 summary | — | — | — |

**限制**：

- 目前只能說 strict region summary 完成 4/7。
- 只有 HCC1395 已完成本稽核中的 chr1–22 `k≤12` Python/C++ partition。
- H2009 extraction 當時仍在背景執行；表格是時間點快照，不能視為最終狀態。
- 不可把 4/7 strict summary 寫成「七樣本完整驗證完成」。

## 8. 執行命令、輸入、輸出與實際片段

### 8.1 Targeted tests

輸入：

- `InterSubMod/tests/test_strict_endpoint_graph.py`
- `InterSubMod/tests/test_build_strict_ps_hp_regions.py`
- `InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_k12_partition.py`

命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  tests/test_strict_endpoint_graph.py \
  tests/test_build_strict_ps_hp_regions.py \
  research/20260722_exact_ps_k12_hcc1395_pilot/tests/test_exact_ps_k12_partition.py
```

輸出：terminal（無新結果檔）。

實際片段：

```text
............................. [100%]
29 passed in 2.08s
```

其中 `test_ab_and_bc_transitively_connect_one_component` 明確測試 A–B、B–C 合成一個 component，並驗證 `(A,C)` 不會被虛構成 direct support。

### 8.2 Constraint disposition 重算

輸入：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2/
chromosomes/chr{1..22}/python_partition/dispositions.tsv.gz
```

命令：

```bash
ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/\
20260723_production_exact_ps_strict_read_linkage/hcc1395_partition_v2/chromosomes
find "$ROOT" -path '*/python_partition/dispositions.tsv.gz' -type f -print0 |
  sort -zV | xargs -0 zcat |
  awk -F '\t' '
    NR==1 {next}
    $1=="dataset" {next}
    {
      rows++; rc[$11]++; wt[$11]+=$9; pc[$11]+=$10; ids[$4]++
    }
    END {
      for (d in rc)
        printf "disposition=%s\trows=%d\tmolecule_weight=%.0f\tpattern_count=%.0f\n",
               d,rc[d],wt[d],pc[d]
      printf "TOTAL\trows=%d\tmolecule_weight=%.0f\tpattern_count=%.0f\n",
             rows,
             wt["retained"]+wt["cut"]+wt["unavoidable_span_gt_max_block_size"],
             pc["retained"]+pc["cut"]+pc["unavoidable_span_gt_max_block_size"]
      dup=0; for(i in ids) if(ids[i]!=1)dup++
      printf "constraint_id_nonunique=%d\n",dup
    }'
```

輸出：terminal（本報告記錄結果）。

實際片段：

```text
TOTAL rows=57629 molecule_weight=285596 pattern_count=57629
constraint_id_nonunique=0
retained rows=53918 molecule_weight=281685
cut rows=1078 molecule_weight=1242
unavoidable_span_gt_max_block_size rows=2633 molecule_weight=2669
```

### 8.3 Receipt 與 parity

輸入：

- HCC1395 22 個 `python_partition/receipt.json`
- HCC1395 22 個 `comparison.json`
- `hcc1395_cpp_parity_v2/genome_parity.json`

命令：

```bash
find "$ROOT" -path '*/python_partition/receipt.json' -type f -print0 |
  sort -zV | xargs -0 jq -r \
  '[.scope.chrom,.all_pass,.counts.constraints,.counts.retained_constraints,
    .counts.cut_constraints,.counts.unavoidable_constraints,.counts.blocks,
    .constraint_mass.total,.constraint_mass.retained,.constraint_mass.cut,
    .constraint_mass.unavoidable] | @tsv'
```

實際摘要：

```text
receipts=22 pass=22 constraints=57629 retained_rows=53918
cut_rows=1078 unavoidable_rows=2633 blocks=11712
total_weight=285596 retained_weight=281685 cut_weight=1242
unavoidable_weight=2669
partition_comparisons=22 all_pass=22 total_mismatches=0
```

### 8.4 凍結 identity

```text
7260a763...29aedd  scripts/build_strict_ps_hp_regions.py
df3d6f37...20e37e  tools/strict_endpoint_graph.py
b71215a8...50e39d  exact_ps_k12_partition.py
f5e77ef2...b340aa  read_support_partition.py
b82e7117...9362d  HCC1395 strict summary.json
2180fe2a...5f4    HCC1395 genome_parity.json
f0b3bb54...a546   ordered digest of 22 Python partition receipts
```

### 8.5 Disconnected 與 zero-weight block 再現命令

輸入：

- `hcc1395_strict_regions_v2/chromosomes/chr*/HCC1395.chr*.endpoint_edges.tsv.gz`
- `hcc1395_partition_v2/chromosomes/chr*/python_partition/{units,blocks}.tsv.gz`

執行命令：

```bash
python - <<'PY'
import csv, gzip
from collections import Counter, defaultdict
from pathlib import Path

strict = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage/"
    "hcc1395_strict_regions_v2/chromosomes"
)
part = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260723_production_exact_ps_strict_read_linkage/"
    "hcc1395_partition_v2/chromosomes"
)
disconnected, zeros = [], []

for ci in range(1, 23):
    chrom = f"chr{ci}"
    p = part / chrom / "python_partition"
    with gzip.open(p / "units.tsv.gz", "rt") as handle:
        units = {
            row["unit_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    edges = defaultdict(list)
    edge_path = strict / chrom / f"HCC1395.{chrom}.endpoint_edges.tsv.gz"
    with gzip.open(edge_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["support_total"]) >= 3:
                edges[(row["linkage_basis"], row["phase_set"])].append(
                    (
                        int(row["pos_i1"]),
                        int(row["pos_j1"]),
                        int(row["support_total"]),
                        row["support_RR"],
                        row["support_RA"],
                        row["support_AR"],
                        row["support_AA"],
                    )
                )
    with gzip.open(p / "blocks.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            positions = tuple(map(int, row["positions1"].split(",")))
            unit = units[row["unit_id"]]
            if len(positions) >= 2 and int(row["retained_molecule_weight"]) == 0:
                zeros.append(
                    (
                        chrom,
                        unit["component_id"],
                        row["unit_id"],
                        int(unit["k"]),
                        row["block_id"],
                        len(positions),
                        row["positions1"],
                    )
                )
            if len(positions) < 2:
                continue
            position_set = set(positions)
            parent = {position: position for position in positions}

            def find(value):
                while parent[value] != value:
                    parent[value] = parent[parent[value]]
                    value = parent[value]
                return value

            def union(left, right):
                left, right = find(left), find(right)
                if left != right:
                    parent[right] = left

            internal = []
            for edge in edges[(row["linkage_basis"], row["phase_set"])]:
                left, right = edge[:2]
                if left in position_set and right in position_set:
                    union(left, right)
                    internal.append(edge)
            groups = defaultdict(list)
            for position in positions:
                groups[find(position)].append(position)
            if len(groups) > 1:
                disconnected.append(
                    (
                        chrom,
                        unit["component_id"],
                        row["unit_id"],
                        int(unit["k"]),
                        row["block_id"],
                        len(positions),
                        int(row["retained_molecule_weight"]),
                        int(row["retained_pattern_count"]),
                        len(internal),
                        tuple(sorted(map(len, groups.values()), reverse=True)),
                        row["positions1"],
                        internal,
                    )
                )

print("DISCONNECTED")
for row in disconnected:
    print("\t".join(map(str, row[:-1])))
    for edge in row[-1]:
        print("  EDGE", *edge, sep="\t")

print("ZERO_KGE2")
for row in zeros:
    print("\t".join(map(str, row)))
print("zero_count", len(zeros))
print("unique_source_units", len({row[2] for row in zeros}))
print("all_source_k_gt12", all(row[3] > 12 for row in zeros))
print(
    "by_chrom",
    sorted(
        Counter(row[0] for row in zeros).items(),
        key=lambda item: int(item[0][3:]),
    ),
)
print("by_block_k", sorted(Counter(row[5] for row in zeros).items()))
PY
```

輸出：terminal；不寫入或修改 production。

實際結尾：

```text
zero_count 134
unique_source_units 58
all_source_k_gt12 True
by_chrom [('chr2', 1), ('chr4', 1), ('chr6', 56), ('chr8', 26),
          ('chr14', 2), ('chr16', 44), ('chr17', 1), ('chr22', 3)]
by_block_k [(2, 5), (3, 11), (4, 7), (5, 7), (6, 9), (7, 9),
            (8, 12), (9, 13), (10, 9), (11, 19), (12, 33)]
```

## 9. 結論與下一步 gate

### 已確認

1. `W` 與 `B` 是不同單位；`k>12` 不刪除 W，只切 B。
2. A–B、B–C transitivity 是 connected-component 定義，程式與 tests 一致。
3. HCC1395 source/partition 守恆與 Python/C++ parity 通過。
4. constraint disposition 的 row 與 weight 兩種分母均獨立重算守恆。

### 尚未解決／應修正

1. **block connectivity gate**：現行 DP 未強制 child block read-linked connected；實測 4 個例外。
2. **zero-evidence block gate**：134 個 `k≥2` blocks retained weight=0，不得直接建樹。
3. **RR-only sensitivity**：23.7737% primary edges 為 RR-only；W 對 edge-state policy 有明顯敏感度。
4. **large-W global inference**：local blocks 尚無經證明的無損 global stitching。
5. **全七資料集**：只完成 4/7 strict summary，其他仍 in progress。
6. **k 上限定位**：12 是現行工程上限，不是生物常數，也不是理論上所有 `k>12` 都必須切的定理；未來可依 active-k 與 certified solver 路由，但尚未 production 化。

## 10. Claim ceiling

本稽核可以支持：

> exact PS×HP 容器內，threshold-qualified direct endpoint edges 的 connected components 定義 W；large W 經 constraint-preserving DP 切成 k≤12 computational blocks，且保留 source component identity。

本稽核不能支持：

- W 或 B 就是真實 clone/subclone。
- A–B、B–C 代表 A–C direct co-occurrence。
- 所有 output blocks 都可建樹。
- local block trees 可唯一、無損地接成 global tree。
- 七資料集完整驗證已完成。

---

`PARTIAL scope footer`：七資料集 production 尚在執行；跨資料集數字只反映 2026-07-23 12:32 +08:00 的本地檔案快照。HCC1395 chr1–22 的 strict-region 與 partition 重算為完整 scope。
