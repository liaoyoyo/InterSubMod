<!--
建立時間: 2026-07-11T17:43:38+08:00
目標: 校正 C/T/Topo 定義，以 HCC1395 教學後提供全 7 dataset rows 完整工程 census
處理範圍: chr1–22；historical layered-v2；PARTIAL（clean layered-v3 7/7 pending）
服務目標: G3 / G4 / G5
Git: branch=research/subclonal-reconstruction-202606; commit=6067568637088838a9f518955e41d222f057e4f1; dirty=True
關聯檔案: InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json
-->

# HCC1395 read 群 C、Steiner 樹 T 與拓撲形狀：全樣本數據報告

用 SCQA + Feynman：先直答定義與數量，再用兩個 HCC1395 例子拆解，最後給全樣本表、方法與限制。

## 🔴 PARTIAL / SCIENTIFIC NO-GO：完整 7-row 歷史工程 census，不是 clean-v3 終版

截至 `2026-07-11T17:06:40+08:00`，normalized raw-all producer 為 **4/7 PASS**；目前 status ledger 的非終態樣本為 **H2009**，尚未開始為 **HCC1937, HCC1954**。clean layered-v3 的 7/7 aggregate 驗證仍未產生。

因此本 HTML 可用於 **C/T/Topo 定義、HCC1395 教學、歷史 7-row 工程數量與拓撲形狀 census**；不可升級成 clean-v3 biological validation 或真實 clone phylogeny 結論。

## TL;DR

- `C` 已修正為 **reads 支持的 mutation-bearing full genotype 群數**，不再表示候選樹組合數。
- `T` 表示 exact rooted directed Steiner tree；雙 HP 是 ordered forest，joint exact-T 候選數為 Cartesian product。
- HCC1395 有 **7,590 primary regions / 10,355 primary HP units**；complete units 中共枚舉 **53,467 exact T**，看到 **32 種 unit tree shapes**。
- 全 7 dataset rows 有 **421,687 exact unit-tree candidates、46 種 rooted shapes**；HCC1395 的 ordered regional forest 另有 **311 種 signatures**。
- HCC1395 的 5,555 ambiguous units 中，5,473 有完整逐位點 read-AF；**3,747/5,473 = 68.46%** 達 top weight≥0.60。這是「VAF-supported most-likely candidate」，不是獨立確認的真樹。

## 先把 C、T、Topo 三件事分清楚

`C_(r,h)`：region `r`、primary HP family `h` 中，`populations_by_hp[h]` 內 **含 A、完整跨越 k 個位點、且 read count≥MINREAD=3** 的 R/A genotype keys 數。`ROOT=RR…R`、partial R/A/X constraints、以及額外 `H_` 節點都不算 C。

對一棵具體 HP tree：

```text
|V(T)| = 1 ROOT + C + H_T
```

`T`：一棵 exact rooted directed Steiner arborescence candidate。`N_T,unit` 是一個 region×HP unit 的 exact T 數；雙 HP region 是保留 HP 身分的 ordered forest：

```text
F_region = (T_HP1, T_HP2)
N_T,region = N_T,HP1 × N_T,HP2
```

`τ(T)`：移除 mutation label 與 sibling 順序、但保留 ROOT 與方向後的形狀。`Topo_1 / Topo_n` 是「候選集合內有 1 / >1 種形狀」；`TopoShape-ID` 則是全資料實際看到哪一種 canonical shape，兩者不可混稱。邏輯上 `T=1 ⇒ Topo=1`，所以 `T=1 / Topo>1` 不可能。

## HCC1395 教學例 1：C=2、T=1、Topo=1

| R/A group | 支持 reads | 計入 C |
|---|---|---|
| AR | 5 | 是（含 A） |
| AA | 3 | 是（含 A） |

| Region | HP | k | C | T | Topo | H | Edges |
|---|---|---|---|---|---|---|---|
| chr1:3294434-3310766 | HP2 | 2 | 2 | 1 | 1 | 0 | ROOT→AR; AR→AA |

## 為何逐位點 VAF 可以在多棵 T 中找 top candidate？

在無 copy-number／purity／multiplicity confounding、沒有回突變的理想情況，descendant mutation 所在細胞集合應是 ancestor mutation 的子集合，因此通常預期 `VAF_ancestor ≥ VAF_descendant`。本流程用同一 HP family 內 BAM reads 重算：

```text
r_j = ALT_reads_j / (REF_reads_j + ALT_reads_j)
Score(T) = Σ(parent→child) Σ(ancestor i in parent) (r_i − r_new-child)
weight(T) ∝ exp(Score(T) / 0.05)
```

教學例 2 是 HCC1395 `chr4:40979546-40998095`、HP2、CN=neutral：兩棵 exact T 的 topology 都是 chain，但一棵假設位點 1 先發生，另一棵假設位點 2 先發生。位點 read-AF 為 **0.9844** 與 **0.5098**，因此 `ROOT→H_AR→AA` 得到正的 VAF 差，relative weight=0.999999994。

這能在**固定 heuristic 與候選集合內**指出最有支持的相對 ordering；但 softmax weight 未校準成 Bayesian posterior，而且 raw read-AF 不是 purity/CN-corrected CCF，所以安全用語是 **VAF-supported most-likely candidate**。

| 位點 | 座標 | REF reads | ALT reads | read-AF |
|---|---|---|---|---|
| 位點 1 | 40,979,546 | 1 | 63 | 98.44% |
| 位點 2 | 40,998,095 | 25 | 26 | 50.98% |

| Candidate | Exact T | Shape | Score | Weight | 判定 |
|---|---|---|---|---|---|
| T_1 | ROOT→H_AR; H_AR→AA | ((())) | 0.474571 | 100.00% | VAF-supported top |
| T_2 | ROOT→H_RA; H_RA→AA | ((())) | -0.474571 | 0.00% | 否 |

## HCC1395 完整數據

| 層級 | 數量 | 意義 |
|---|---|---|
| S：輸入 sSNV（全 contigs） | 113,997 | VCF 中全部 sSNV |
| S：chr1–22 | 80,234 | 本次 W/C/T/Topo 分析母體 |
| W：k=1 | 8,447 | 單點區域，不做多位點樹 |
| W：k>1 pre-cap | 7,931 | 相鄰 gap≤50 kb 的多位點區域 |
| W：retained | 7,928 | 通過讀段支援並保留 |
| W：tree view | 7,927 | 至少產生一個 lineage unit |
| W：primary HP | 7,590 | HP1 xor HP2 或 HP1 and HP2 |
| primary HP units | 10,355 | region×HP，雙 HP 區域算兩個 unit |

### HP × HP3

| Primary HP class | not H3 | with H3 | 合計 |
|---|---|---|---|
| HP1 xor HP2 | 4,146 | 679 | 4,825 |
| HP1 and HP2 | 1,859 | 906 | 2,765 |
| 無 primary HP1/HP2 | 231 | 106 | 337 |

### C：tree-aligned per-HP unit

| C | Primary HP units |
|---|---|
| 0 | 6,455 |
| 1 | 2,703 |
| 2 | 1,132 |
| 3 | 61 |
| 4 | 4 |
| 5 | 0 |
| 6 | 0 |
| >6 | 0 |

### C：pooled region raw census

| C | Regions |
|---|---|
| 0 | 4,438 |
| 1 | 1,628 |
| 2 | 1,716 |
| 3 | 138 |
| 4 | 6 |
| 5 | 2 |
| 6 | 0 |
| >6 | 0 |

Primary regions 內 pooled C 與 primary HP1/2 C_sum 不同：**98/7,590**；全部 retained regions 對 HP1+HP2 為 **291/7,928**；若把 HP3/none 等全部 family 都納入，則為 **145/7,928**。差異來自跨 family pooling 的 collapse／threshold crossing，因此 pooled C 不可直接對 per-HP T。

### 雙 HP ordered C pair

| (C_HP1,C_HP2) | Regions |
|---|---|
| HP1=0, HP2=0 | 1,885 |
| HP1=0, HP2=1 | 127 |
| HP1=0, HP2=2 | 5 |
| HP1=1, HP2=0 | 103 |
| HP1=1, HP2=1 | 562 |
| HP1=1, HP2=2 | 31 |
| HP1=2, HP2=0 | 12 |
| HP1=2, HP2=1 | 35 |
| HP1=2, HP2=2 | 2 |
| HP1=2, HP2=3 | 1 |
| HP1=3, HP2=0 | 1 |
| HP1=3, HP2=1 | 1 |

### T / Topo identifiability

| Class | Regions |
|---|---|
| T=1；Topo=1 | 2,032 |
| T>1；Topo=1 | 1,406 |
| T>1；Topo>1 | 3,502 |
| 候選集不完整 | 650 |

| Class | H=0 | H>0 | Total |
|---|---|---|---|
| T=1；Topo=1 | 1,443 | 589 | 2,032 |
| T>1；Topo=1 | 16 | 1,390 | 1,406 |
| T>1；Topo>1 | 1 | 3,501 | 3,502 |

### 單/雙 HP × C × H × T/Topo 完整巢狀交叉

此表是最初問題的直接答案：單 HP 顯示 `HP1=n` 或 `HP2=n`；雙 HP 顯示 ordered `HP1=n, HP2=m`，不以 C_sum 取代。雙 HP 的 `H=0` 表示兩棵都無 extra state；`H>0` 表示至少一棵有 extra state。

| HP scope | C / ordered pair | H | T / Topo | Regions |
|---|---|---|---|---|
| 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 308 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 140 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,162 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 275 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 66 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 30 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 30 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 3 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 52 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 26 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 25 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 425 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 9 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 83 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 38 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 7 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 26 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T=1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 3 |
| 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 10 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo>1 | 2 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 29 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T=1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 5 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=2, HP2=3 | H=0 | T=1；Topo=1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=3, HP2=0 | H>0 | T>1；Topo>1 | 1 |
| 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H=0 | T=1；Topo=1 | 1 |
| 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 59 |
| 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 15 |
| 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 984 |
| 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 169 |
| 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 26 |
| 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 483 |
| 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 152 |
| 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 15 |
| 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 460 |
| 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 15 |
| 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 64 |
| 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 7 |
| 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 3 |
| 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 13 |
| 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T=1；Topo=1 | 2 |
| 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 6 |
| 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 4 |
| 單 HP（HP1 xor HP2） | HP1=4 | H=0 | T>1；Topo>1 | 1 |
| 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 53 |
| 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 21 |
| 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 962 |
| 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 174 |
| 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 25 |
| 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T=1；Topo=1 | 2 |
| 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 443 |
| 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 131 |
| 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 5 |
| 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 418 |
| 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 8 |
| 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 64 |
| 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo>1 | 4 |
| 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 1 |
| 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 17 |
| 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 10 |
| 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo=1 | 6 |
| 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T=1；Topo=1 | 1 |
| 單 HP（HP1 xor HP2） | HP2=4 | H>0 | T>1；Topo=1 | 1 |
| 單 HP（HP1 xor HP2） | HP2=4 | H>0 | T>1；Topo>1 | 1 |

HCC1395 complete primary regions=6,940；ordered topology alternatives 加總=17,586；joint exact-T candidate combinations 加總=195,473；candidate catalog 枚舉到 ordered regional signatures=311。

## 拓撲形狀怎麼看？

Canonical rule 是 `τ(v) = "(" + sort(τ(children)) + ")"`：

```text
(())        ROOT──●             單一 mutation-state node
((()))      ROOT──●──●          linear chain
(()())      ROOT──┬─●           root star（兩個 leaves）
                  └─●
((()()))    ROOT──●──┬─●        internal branch
                     └─●
```

完整 catalog 以 **unit support** 排名，因為 candidate occurrences 會讓高歧義 unit 被重複加權；同一 unit 可支持多種 shape，所以 unit support 跨 shape 不互斥。`Shape-only units` 表示該 unit 的完整候選集合只有這一種 topology shape。相同 graph shape 可以跨不同 k 出現，只代表 branching skeleton 相同，不代表相同突變或相同生物歷史。Incomplete/capped units 未列入 exact catalog，因此 46 種不是理論上限。

### HCC1395 unit support 前 10 種 shape

| TopoShape-ID | Signature | 形狀族 | HCC units | Global units | Global candidate T |
|---|---|---|---|---|---|
| Topo_01 | ((())) | 鏈狀 | 3,962 | 22,833 | 41,028 |
| Topo_02 | (()) | 單節點 | 2,885 | 19,047 | 19,047 |
| Topo_03 | (()()) | 根部分岔 | 2,329 | 11,268 | 11,268 |
| Topo_04 | (((()))) | 鏈狀 | 1,448 | 9,980 | 49,242 |
| Topo_05 | ((()())) | 內部分岔 | 1,327 | 7,656 | 14,063 |
| Topo_06 | ((())()) | 根部分岔 | 1,327 | 6,999 | 20,752 |
| Topo_08 | (((())())) | 內部分岔 | 650 | 4,706 | 39,289 |
| Topo_09 | (((()()))) | 內部分岔 | 582 | 4,658 | 27,628 |
| Topo_10 | (((()))()) | 根部分岔 | 556 | 3,323 | 30,503 |
| Topo_07 | ((((())))) | 鏈狀 | 509 | 5,234 | 96,933 |

## 全部 7 dataset rows 數據

### S / W

| Dataset | S all | S chr1–22 | W k=1 | W k>1 | W retained | W tree | W primary | Primary units |
|---|---|---|---|---|---|---|---|---|
| HCC1395 | 113,997 | 80,234 | 8,447 | 7,931 | 7,928 | 7,927 | 7,590 | 10,355 |
| HCC1395_DORADO | 112,387 | 79,120 | 8,406 | 7,958 | 7,958 | 7,958 | 7,268 | 10,025 |
| COLO829 | 38,196 | 36,585 | 7,841 | 7,774 | 7,774 | 7,774 | 7,659 | 11,772 |
| H1437 | 75,578 | 73,243 | 6,795 | 8,865 | 8,865 | 8,865 | 8,630 | 12,821 |
| H2009 | 157,405 | 150,370 | 3,017 | 9,717 | 9,717 | 9,717 | 9,581 | 14,247 |
| HCC1937 | 49,548 | 15,915 | 7,705 | 2,695 | 2,695 | 2,695 | 2,674 | 3,364 |
| HCC1954 | 20,969 | 19,743 | 7,767 | 4,023 | 4,023 | 4,023 | 3,975 | 5,960 |

### HP × HP3

| Dataset | xor/not H3 | xor/H3 | and/not H3 | and/H3 | none/not H3 | none/H3 |
|---|---|---|---|---|---|---|
| HCC1395 | 4,146 | 679 | 1,859 | 906 | 231 | 106 |
| HCC1395_DORADO | 3,500 | 1,011 | 1,902 | 855 | 182 | 508 |
| COLO829 | 3,330 | 216 | 3,696 | 417 | 3 | 112 |
| H1437 | 4,191 | 248 | 3,608 | 583 | 16 | 219 |
| H2009 | 4,288 | 627 | 3,379 | 1,287 | 19 | 117 |
| HCC1937 | 1,878 | 106 | 560 | 130 | 5 | 16 |
| HCC1954 | 1,834 | 156 | 1,771 | 214 | 3 | 45 |

### Pooled region C

| Dataset | C=0 | C=1 | C=2 | C=3 | C=4 | C=5 | C=6 | C=>6 |
|---|---|---|---|---|---|---|---|---|
| HCC1395 | 4,438 | 1,628 | 1,716 | 138 | 6 | 2 | 0 | 0 |
| HCC1395_DORADO | 6,394 | 766 | 755 | 42 | 1 | 0 | 0 | 0 |
| COLO829 | 4,277 | 2,348 | 1,141 | 8 | 0 | 0 | 0 | 0 |
| H1437 | 3,979 | 2,594 | 2,097 | 187 | 8 | 0 | 0 | 0 |
| H2009 | 5,562 | 1,764 | 2,035 | 306 | 41 | 9 | 0 | 0 |
| HCC1937 | 530 | 828 | 1,236 | 95 | 5 | 0 | 1 | 0 |
| HCC1954 | 2,023 | 839 | 1,090 | 66 | 5 | 0 | 0 | 0 |

### Per-primary-HP unit C

| Dataset | C=0 | C=1 | C=2 | C=3 | C=4 | C=5 | C=6 | C=>6 |
|---|---|---|---|---|---|---|---|---|
| HCC1395 | 6,455 | 2,703 | 1,132 | 61 | 4 | 0 | 0 | 0 |
| HCC1395_DORADO | 8,550 | 1,020 | 440 | 15 | 0 | 0 | 0 | 0 |
| COLO829 | 7,408 | 4,192 | 172 | 0 | 0 | 0 | 0 | 0 |
| H1437 | 6,817 | 4,957 | 1,006 | 40 | 1 | 0 | 0 | 0 |
| H2009 | 9,006 | 3,770 | 1,392 | 72 | 7 | 0 | 0 | 0 |
| HCC1937 | 870 | 1,496 | 939 | 56 | 2 | 0 | 1 | 0 |
| HCC1954 | 3,428 | 1,925 | 590 | 16 | 1 | 0 | 0 | 0 |

### T / Topo / unit shapes

| Dataset | T1/Topo1 (regions) | Tn/Topo1 (regions) | Tn/Topon (regions) | Incomplete (regions) | H=0 (regions) | H>0 (regions) | Complete units | Incomplete units | Exact T (unit candidates) | Unit-shape incidences | Shape types |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HCC1395 | 2,032 | 1,406 | 3,502 | 650 | 1,460 | 5,480 | 9,702 | 653 | 53,467 | 18,024 | 32 |
| HCC1395_DORADO | 1,814 | 630 | 4,306 | 518 | 593 | 6,157 | 9,506 | 519 | 77,764 | 20,090 | 23 |
| COLO829 | 1,406 | 2,470 | 3,073 | 710 | 822 | 6,127 | 11,061 | 711 | 83,813 | 20,066 | 20 |
| H1437 | 1,657 | 2,927 | 2,400 | 1,646 | 1,400 | 5,584 | 11,168 | 1,653 | 70,796 | 17,933 | 27 |
| H2009 | 1,170 | 2,358 | 2,354 | 3,699 | 994 | 4,888 | 10,514 | 3,733 | 104,623 | 18,913 | 39 |
| HCC1937 | 1,217 | 724 | 616 | 117 | 1,073 | 1,484 | 3,243 | 121 | 9,389 | 4,235 | 24 |
| HCC1954 | 1,536 | 629 | 1,658 | 152 | 983 | 2,840 | 5,808 | 152 | 21,835 | 9,144 | 23 |

### 全 7 dataset 單/雙 HP × C × H × T/Topo 長表

| Dataset | HP scope | C / ordered pair | H | T / Topo | Regions |
|---|---|---|---|---|---|
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 308 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 140 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,162 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 275 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 66 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 30 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 30 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 3 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 52 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 26 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 25 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 425 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 9 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 83 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 38 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 7 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 26 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T=1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 3 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 10 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo>1 | 2 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 29 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T=1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 5 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=2, HP2=3 | H=0 | T=1；Topo=1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=3, HP2=0 | H>0 | T>1；Topo>1 | 1 |
| HCC1395 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H=0 | T=1；Topo=1 | 1 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 59 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 15 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 984 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 169 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 26 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 483 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 152 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 15 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 460 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 15 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 64 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 7 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 3 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 13 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T=1；Topo=1 | 2 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 6 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 4 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP1=4 | H=0 | T>1；Topo>1 | 1 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 53 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 21 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 962 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 174 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 25 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T=1；Topo=1 | 2 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 443 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 131 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 5 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 418 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 8 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 64 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo>1 | 4 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 1 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 17 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 10 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo=1 | 6 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T=1；Topo=1 | 1 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=4 | H>0 | T>1；Topo=1 | 1 |
| HCC1395 | 單 HP（HP1 xor HP2） | HP2=4 | H>0 | T>1；Topo>1 | 1 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 657 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 90 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,484 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 233 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 31 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 6 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 3 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 3 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 27 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 11 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 5 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 166 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 2 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 12 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 3 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 2 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 6 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 2 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo=1 | 1 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 9 |
| HCC1395_DORADO | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 2 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 252 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 21 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 1,304 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 126 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 16 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 186 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 44 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 9 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 190 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 1 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 24 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 1 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 6 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 1 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 2 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 244 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 23 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 1,425 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 137 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 14 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 236 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 37 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 8 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 177 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 3 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 17 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 1 |
| HCC1395_DORADO | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 6 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 292 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 270 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,636 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 463 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 119 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 86 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 34 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 121 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 85 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 49 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | not_evaluated | 候選集不完整 | 2 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 656 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 1 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 279 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 6 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 4 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 2 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 3 |
| COLO829 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 2 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 26 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 7 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 653 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 108 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 15 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 848 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 36 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 7 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 68 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 1 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 15 |
| COLO829 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 1 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 24 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 16 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 633 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 121 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 15 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 847 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 26 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 6 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 59 |
| COLO829 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 14 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 75 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 231 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,199 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 861 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 74 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 88 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 64 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 17 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 3 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 4 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 59 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 85 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 73 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | not_evaluated | 候選集不完整 | 23 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 691 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 12 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 424 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 34 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 13 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 33 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T=1；Topo=1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 31 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | not_evaluated | 候選集不完整 | 5 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=1, HP2=3 | H>0 | T>1；Topo>1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 7 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo=1 | 4 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo>1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 35 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T=1；Topo=1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 28 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo>1 | 2 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | not_evaluated | 候選集不完整 | 5 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 2 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=2, HP2=3 | H>0 | T=1；Topo=1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=3, HP2=0 | H>0 | T>1；Topo>1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H=0 | T=1；Topo=1 | 1 |
| H1437 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H>0 | T>1；Topo=1 | 1 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 9 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 17 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 435 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 374 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 1 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T=1；Topo=1 | 2 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 941 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 91 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 12 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 286 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 7 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 116 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 5 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 5 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 8 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 6 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 4 |
| H1437 | 單 HP（HP1 xor HP2） | HP1=3 | not_evaluated | 候選集不完整 | 2 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 5 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 19 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 404 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 319 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 5 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T=1；Topo=1 | 2 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 836 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 85 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 7 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 318 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 6 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 88 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo>1 | 5 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 3 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 12 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 1 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo=1 | 2 |
| H1437 | 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T=1；Topo=1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 54 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 170 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 1,290 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 1,523 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 28 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 68 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 69 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 47 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 3 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 7 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo>1 | 2 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | not_evaluated | 候選集不完整 | 9 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 22 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 70 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 47 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | not_evaluated | 候選集不完整 | 42 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 303 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 6 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 429 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 25 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 107 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 34 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T=1；Topo=1 | 2 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 72 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 3 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | not_evaluated | 候選集不完整 | 31 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=3 | H=0 | T=1；Topo=1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=3 | H>0 | T>1；Topo=1 | 3 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=3 | not_evaluated | 候選集不完整 | 2 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=4 | H>0 | T>1；Topo=1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=1, HP2=4 | not_evaluated | 候選集不完整 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 11 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo=1 | 8 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo>1 | 3 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | not_evaluated | 候選集不完整 | 11 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 25 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T=1；Topo=1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 59 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo>1 | 5 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | not_evaluated | 候選集不完整 | 23 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 5 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H>0 | T>1；Topo=1 | 16 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | not_evaluated | 候選集不完整 | 5 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=2, HP2=3 | H>0 | T>1；Topo=1 | 2 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=0 | H>0 | T=1；Topo=1 | 2 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=0 | H>0 | T>1；Topo=1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H>0 | T>1；Topo=1 | 6 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | H>0 | T>1；Topo>1 | 1 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | not_evaluated | 候選集不完整 | 6 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=3, HP2=2 | H>0 | T>1；Topo=1 | 4 |
| H2009 | 雙 HP（HP1 and HP2） | HP1=4, HP2=1 | H>0 | T>1；Topo=1 | 1 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 11 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 23 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 381 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 819 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 5 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T=1；Topo=1 | 4 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 516 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 64 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 68 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 260 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 19 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 145 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 7 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 47 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 5 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 10 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 6 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo>1 | 1 |
| H2009 | 單 HP（HP1 xor HP2） | HP1=3 | not_evaluated | 候選集不完整 | 1 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 11 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 22 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 391 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 824 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 4 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T=1；Topo=1 | 3 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 556 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 65 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 75 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 324 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 20 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 149 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 54 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 7 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 9 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo=1 | 1 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=3 | not_evaluated | 候選集不完整 | 4 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T>1；Topo=1 | 2 |
| H2009 | 單 HP（HP1 xor HP2） | HP2=4 | H>0 | T>1；Topo=1 | 2 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 45 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 30 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 121 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 12 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 31 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 14 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 8 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | not_evaluated | 候選集不完整 | 2 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 8 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 34 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 7 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 9 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | not_evaluated | 候選集不完整 | 1 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 240 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 6 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 44 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 17 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | not_evaluated | 候選集不完整 | 5 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 10 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 3 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 1 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 6 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo=1 | 1 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | not_evaluated | 候選集不完整 | 1 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 20 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 4 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo>1 | 2 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 6 |
| HCC1937 | 雙 HP（HP1 and HP2） | HP1=3, HP2=1 | not_evaluated | 候選集不完整 | 1 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 10 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 6 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 130 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 11 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 5 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 258 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 89 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 23 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 370 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 7 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 38 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 6 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 5 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 13 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 8 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 5 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=3 | not_evaluated | 候選集不完整 | 2 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=4 | not_evaluated | 候選集不完整 | 1 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP1=6 | H=0 | T>1；Topo>1 | 1 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 8 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 2 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 149 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 15 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 2 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 241 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 78 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 29 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 375 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 7 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 49 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo>1 | 4 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 9 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 14 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 8 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo=1 | 4 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=3 | H>0 | T>1；Topo>1 | 1 |
| HCC1937 | 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T>1；Topo=1 | 1 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T=1；Topo=1 | 346 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo=1 | 60 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | H>0 | T>1；Topo>1 | 683 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=0 | not_evaluated | 候選集不完整 | 82 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T=1；Topo=1 | 81 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo=1 | 20 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=1 | H>0 | T>1；Topo>1 | 19 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=0, HP2=2 | H>0 | T=1；Topo=1 | 6 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T=1；Topo=1 | 66 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo=1 | 19 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | H>0 | T>1；Topo>1 | 15 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=0 | not_evaluated | 候選集不完整 | 2 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H=0 | T=1；Topo=1 | 443 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T=1；Topo=1 | 3 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo=1 | 61 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=1 | H>0 | T>1；Topo>1 | 14 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H=0 | T=1；Topo=1 | 24 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo=1 | 1 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=2 | H>0 | T>1；Topo>1 | 2 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=1, HP2=3 | H=0 | T=1；Topo=1 | 1 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T=1；Topo=1 | 7 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=0 | H>0 | T>1；Topo=1 | 2 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H=0 | T=1；Topo=1 | 20 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo=1 | 4 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=1 | H>0 | T>1；Topo>1 | 1 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | H=0 | T=1；Topo=1 | 2 |
| HCC1954 | 雙 HP（HP1 and HP2） | HP1=2, HP2=2 | not_evaluated | 候選集不完整 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T=1；Topo=1 | 20 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo=1 | 4 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=0 | H>0 | T>1；Topo>1 | 395 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=0 | not_evaluated | 候選集不完整 | 20 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=1 | H=0 | T=1；Topo=1 | 6 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T=1；Topo=1 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo=1 | 207 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=1 | H>0 | T>1；Topo>1 | 83 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=1 | not_evaluated | 候選集不完整 | 9 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=2 | H=0 | T=1；Topo=1 | 228 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T=1；Topo=1 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo=1 | 23 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=2 | H>0 | T>1；Topo>1 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=2 | not_evaluated | 候選集不完整 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T=1；Topo=1 | 5 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=3 | H=0 | T>1；Topo=1 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP1=3 | H>0 | T>1；Topo=1 | 1 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T=1；Topo=1 | 19 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo=1 | 4 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=0 | H>0 | T>1；Topo>1 | 364 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=0 | not_evaluated | 候選集不完整 | 23 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=1 | H=0 | T=1；Topo=1 | 4 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T=1；Topo=1 | 3 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo=1 | 205 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=1 | H>0 | T>1；Topo>1 | 79 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=1 | not_evaluated | 候選集不完整 | 11 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=2 | H=0 | T=1；Topo=1 | 240 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T=1；Topo=1 | 4 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo=1 | 14 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=2 | H>0 | T>1；Topo>1 | 2 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=2 | not_evaluated | 候選集不完整 | 3 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T=1；Topo=1 | 5 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=3 | H=0 | T>1；Topo=1 | 3 |
| HCC1954 | 單 HP（HP1 xor HP2） | HP2=4 | H=0 | T=1；Topo=1 | 1 |

### Ordered regional forest

| Dataset | Primary R | Complete | Incomplete | Topo alternatives | Joint exact T | Ordered signatures |
|---|---|---|---|---|---|---|
| HCC1395 | 7,590 | 6,940 | 650 | 17,586 | 195,473 | 311 |
| HCC1395_DORADO | 7,268 | 6,750 | 518 | 20,778 | 417,964 | 296 |
| COLO829 | 7,659 | 6,949 | 710 | 21,552 | 663,218 | 297 |
| H1437 | 8,630 | 6,984 | 1,646 | 16,305 | 357,034 | 322 |
| H2009 | 9,581 | 5,882 | 3,699 | 18,043 | 724,719 | 371 |
| HCC1937 | 2,674 | 2,557 | 117 | 3,645 | 12,676 | 181 |
| HCC1954 | 3,975 | 3,823 | 152 | 8,067 | 72,458 | 297 |

### VAF-supported top candidate

| Dataset | Ambiguous | Prepared | Missing AF | Top≥.60 | Reach/prepared | Direction OK | Neutral prep | Neutral top | Neutral reach |
|---|---|---|---|---|---|---|---|---|---|
| HCC1395 | 5,555 | 5,473 | 82 | 3,747 | 68.46% | 3,742 | 198 | 106 | 53.54% |
| HCC1395_DORADO | 5,576 | 4,973 | 603 | 3,162 | 63.58% | 3,160 | 174 | 97 | 55.75% |
| COLO829 | 6,926 | 6,843 | 83 | 2,071 | 30.26% | 2,071 | 0 | 0 | — |
| H1437 | 7,043 | 6,991 | 52 | 2,804 | 40.11% | 2,804 | 105 | 33 | 31.43% |
| H2009 | 7,430 | 7,264 | 166 | 2,778 | 38.24% | 2,775 | 2,758 | 801 | 29.04% |
| HCC1937 | 1,380 | 1,377 | 3 | 898 | 65.21% | 897 | 0 | 0 | — |
| HCC1954 | 2,552 | 2,536 | 16 | 1,765 | 69.60% | 1,762 | 213 | 145 | 68.08% |

## 完整 46 種 TopoShape catalog

| TopoShape-ID | Stable ID | Signature | 形狀族 | Nodes incl. ROOT | Depth | Root degree | Leaves | Global unit support | Shape-only units | Candidate T | HCC1395 | DORADO | COLO829 | H1437 | H2009 | HCC1937 | HCC1954 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Topo_01 | TS-a99c6ae6b3 | ((())) | 鏈狀 | 3 | 2 | 1 | 1 | 22,833 | 13,059 | 41,028 | 3,962 | 3,749 | 4,316 | 4,163 | 3,272 | 1,411 | 1,960 |
| Topo_02 | TS-aa75362b1c | (()) | 單節點 | 2 | 1 | 1 | 1 | 19,047 | 19,047 | 19,047 | 2,885 | 3,402 | 3,970 | 3,192 | 2,091 | 939 | 2,568 |
| Topo_03 | TS-ca3708ca57 | (()()) | 根部分岔 | 3 | 1 | 2 | 2 | 11,268 | 1,494 | 11,268 | 2,329 | 2,959 | 1,742 | 1,073 | 972 | 691 | 1,502 |
| Topo_04 | TS-23ba77f4e8 | (((()))) | 鏈狀 | 4 | 3 | 1 | 1 | 9,980 | 2,971 | 49,242 | 1,448 | 1,422 | 1,790 | 2,185 | 2,356 | 293 | 486 |
| Topo_05 | TS-30a85152fe | ((()())) | 內部分岔 | 4 | 2 | 1 | 2 | 7,656 | 142 | 14,063 | 1,327 | 1,443 | 1,401 | 1,326 | 1,420 | 202 | 537 |
| Topo_06 | TS-8a90a71679 | ((())()) | 根部分岔 | 4 | 2 | 2 | 2 | 6,999 | 325 | 20,752 | 1,327 | 1,441 | 1,253 | 997 | 1,170 | 231 | 580 |
| Topo_07 | TS-85f04c9db4 | ((((())))) | 鏈狀 | 5 | 4 | 1 | 1 | 5,234 | 963 | 96,933 | 509 | 560 | 929 | 1,217 | 1,814 | 79 | 126 |
| Topo_08 | TS-f3932cb25b | (((())())) | 內部分岔 | 5 | 3 | 1 | 2 | 4,706 | 56 | 39,289 | 650 | 619 | 863 | 943 | 1,402 | 58 | 171 |
| Topo_09 | TS-c17bab058f | (((()()))) | 內部分岔 | 5 | 3 | 1 | 2 | 4,658 | 35 | 27,628 | 582 | 602 | 873 | 993 | 1,408 | 49 | 151 |
| Topo_10 | TS-a5d54651bb | (((()))()) | 根部分岔 | 5 | 3 | 2 | 2 | 3,323 | 39 | 30,503 | 556 | 584 | 651 | 530 | 793 | 43 | 166 |
| Topo_11 | TS-3ad54669f9 | ((())(())) | 根部分岔 | 5 | 2 | 2 | 2 | 2,577 | 22 | 12,648 | 446 | 537 | 551 | 331 | 547 | 29 | 136 |
| Topo_12 | TS-3ffd47e1d7 | (()()()) | 根部分岔 | 4 | 1 | 3 | 3 | 2,398 | 17 | 2,398 | 507 | 874 | 341 | 156 | 173 | 72 | 275 |
| Topo_13 | TS-86db631b4d | ((()())()) | 根與內部皆分岔 | 5 | 2 | 2 | 3 | 2,322 | 1 | 10,849 | 474 | 547 | 449 | 267 | 400 | 37 | 148 |
| Topo_14 | TS-7489157280 | ((()()())) | 內部分岔 | 5 | 2 | 1 | 3 | 2,167 | 0 | 4,304 | 415 | 532 | 445 | 247 | 371 | 29 | 128 |
| Topo_15 | TS-a9ec036834 | ((())()()) | 根部分岔 | 5 | 2 | 3 | 3 | 1,980 | 2 | 8,984 | 427 | 527 | 403 | 165 | 287 | 31 | 140 |
| Topo_16 | TS-5a8f3dfd2d | (()()()()) | 根部分岔 | 5 | 1 | 4 | 4 | 603 | 0 | 603 | 146 | 282 | 79 | 9 | 28 | 5 | 54 |
| Topo_17 | TS-aaa42b13e4 | (((((()))))) | 鏈狀 | 6 | 5 | 1 | 1 | 370 | 336 | 30,115 | 5 | 4 | 7 | 73 | 245 | 28 | 8 |
| Topo_18 | TS-ec7b3dec5d | ((((())()))) | 內部分岔 | 6 | 4 | 1 | 2 | 68 | 20 | 387 | 2 | 0 | 1 | 18 | 43 | 2 | 2 |
| Topo_19 | TS-569b3dc562 | ((((()())))) | 內部分岔 | 6 | 4 | 1 | 2 | 57 | 20 | 471 | 2 | 0 | 1 | 25 | 26 | 1 | 2 |
| Topo_20 | TS-222a69fe42 | ((((()))())) | 內部分岔 | 6 | 4 | 1 | 2 | 44 | 8 | 282 | 2 | 1 | 0 | 8 | 32 | 0 | 1 |
| Topo_21 | TS-d600f311f1 | ((((())))()) | 根部分岔 | 6 | 4 | 2 | 2 | 27 | 7 | 401 | 0 | 1 | 0 | 8 | 17 | 0 | 1 |
| Topo_22 | TS-0d44385f05 | (((()))(())) | 根部分岔 | 6 | 3 | 2 | 2 | 17 | 8 | 149 | 3 | 0 | 1 | 1 | 11 | 0 | 1 |
| Topo_23 | TS-0439f089d9 | (((())(()))) | 內部分岔 | 6 | 3 | 1 | 2 | 12 | 6 | 45 | 1 | 0 | 0 | 0 | 9 | 1 | 1 |
| Topo_24 | TS-23650212ff | (((())())()) | 根與內部皆分岔 | 6 | 3 | 2 | 3 | 9 | 2 | 17 | 5 | 1 | 0 | 1 | 2 | 0 | 0 |
| Topo_25 | TS-85b33acb24 | (((()())())) | 內部分岔 | 6 | 3 | 1 | 3 | 7 | 0 | 10 | 1 | 1 | 0 | 1 | 4 | 0 | 0 |
| Topo_26 | TS-2a6f3a1b20 | (((()()()))) | 內部分岔 | 6 | 3 | 1 | 3 | 6 | 0 | 12 | 0 | 0 | 0 | 2 | 4 | 0 | 0 |
| Topo_27 | TS-a08d0d50d7 | ((()())(())) | 根與內部皆分岔 | 6 | 2 | 2 | 3 | 4 | 0 | 9 | 3 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_28 | TS-47a04b35da | (((((())())))) | 內部分岔 | 7 | 5 | 1 | 2 | 4 | 4 | 42 | 0 | 0 | 0 | 0 | 3 | 1 | 0 |
| Topo_29 | TS-09d313affb | (((())()())) | 內部分岔 | 6 | 3 | 1 | 3 | 3 | 0 | 5 | 1 | 1 | 0 | 1 | 0 | 0 | 0 |
| Topo_30 | TS-6c72c312e1 | (((()))()()) | 根部分岔 | 6 | 3 | 3 | 3 | 3 | 0 | 16 | 1 | 0 | 0 | 1 | 1 | 0 | 0 |
| Topo_31 | TS-a209143e81 | ((())(())()) | 根部分岔 | 6 | 2 | 3 | 3 | 3 | 0 | 9 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |
| Topo_32 | TS-a41e0aa077 | (((((()))))()) | 根部分岔 | 7 | 5 | 2 | 2 | 3 | 3 | 64 | 1 | 1 | 0 | 0 | 1 | 0 | 0 |
| Topo_33 | TS-4d9438c9f0 | (((()()))()) | 根與內部皆分岔 | 6 | 3 | 2 | 3 | 2 | 0 | 4 | 0 | 0 | 0 | 0 | 2 | 0 | 0 |
| Topo_34 | TS-a376051cd0 | ((())()()()) | 根部分岔 | 6 | 2 | 4 | 4 | 2 | 1 | 3 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
| Topo_35 | TS-7a25448194 | ((((((())))))) | 鏈狀 | 7 | 6 | 1 | 1 | 2 | 2 | 30 | 1 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_36 | TS-dc4dc2ebd3 | ((()())()()) | 根與內部皆分岔 | 6 | 2 | 3 | 4 | 1 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| Topo_37 | TS-7e9b0f7389 | (((((()()))))) | 內部分岔 | 7 | 5 | 1 | 2 | 1 | 1 | 24 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_38 | TS-a7b99408f6 | (((((()))()))) | 內部分岔 | 7 | 5 | 1 | 2 | 1 | 1 | 4 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| Topo_39 | TS-49d6f77cb0 | (((((())))())) | 內部分岔 | 7 | 5 | 1 | 2 | 1 | 1 | 4 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_40 | TS-33ad7c612f | ((((()())()))) | 內部分岔 | 7 | 4 | 1 | 3 | 1 | 0 | 2 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_41 | TS-218a883fa4 | ((((())(())))) | 內部分岔 | 7 | 4 | 1 | 2 | 1 | 0 | 4 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_42 | TS-d6066189ed | ((((())())())) | 內部分岔 | 7 | 4 | 1 | 3 | 1 | 0 | 2 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_43 | TS-17a065c268 | ((((()))(()))) | 內部分岔 | 7 | 4 | 1 | 2 | 1 | 0 | 6 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_44 | TS-aeebb18f28 | ((((())))(())) | 根部分岔 | 7 | 4 | 2 | 2 | 1 | 1 | 24 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| Topo_45 | TS-9310d14e5e | (((())())(())) | 根與內部皆分岔 | 7 | 3 | 2 | 3 | 1 | 0 | 4 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| Topo_46 | TS-c90e49d35c | (((()))(()())) | 根與內部皆分岔 | 7 | 3 | 2 | 3 | 1 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |

## 驗證

| Dataset | Checks | Passed | Failed | Status |
|---|---|---|---|---|
| HCC1395 | 26 | 26 | 0 | PASS |
| HCC1395_DORADO | 26 | 26 | 0 | PASS |
| COLO829 | 26 | 26 | 0 | PASS |
| H1437 | 26 | 26 | 0 | PASS |
| H2009 | 26 | 26 | 0 | PASS |
| HCC1937 | 26 | 26 | 0 | PASS |
| HCC1954 | 26 | 26 | 0 | PASS |

## 科學界線與限制

1. 這批數字來自 historical layered-v2 engineering snapshot；latest clean-v3 仍缺 7/7 aggregate validation。
2. `C=0` 只表示沒有達 MINREAD 的 full ALT genotype group；仍可能有 ALT partial-read subcube constraint 並產生 tree。
3. `H_` 目前同時包含真正未觀測 intermediate 與 partial-supported completion state，不能直接解讀為 hidden clone。
4. VAF ranking 是 unit-level、同資料 heuristic；winner direction consistency 不是獨立驗證。
5. Copy number、LOH、purity、mutation multiplicity 與 read-depth uncertainty 都可能改變 VAF ordering。HCC1395 的 top units 多數為 CN-altered，因此報告另列 CN-neutral 子集。
6. HCC1395 與 HCC1395_DORADO 是同一 biological sample 的兩個 dataset rows；7 rows 只代表 6 個 biological samples。

## Provenance

| Artifact | Path | SHA-256 | Bytes |
|---|---|---|---|
| corrected report data | research/20260711_read_group_C_tree_T_topology_report/data/c_t_topology_report_data.json | df60db7e26bab2f47620caa0a3cb6583d6d5488ec8ae0f574bf2fab4acf9a759 | 275,836 |
| validation checks | research/20260711_read_group_C_tree_T_topology_report/data/validation_checks.tsv | 2baf2ad1a65aa0ed70eff07289e22e6e5a7ead447f2fae64234c8600446a131d | 7,482 |
| nested HP/C/T/Topo/H cross | research/20260711_read_group_C_tree_T_topology_report/data/primary_region_HP_C_T_Topo_H_cross.tsv | 81a1883837b1e048a81d6817069efe3d79014640e28070b267004716dd68055e | 27,773 |
| exact topology catalog | research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json | e0f29a01f3ddffd4a76083f922eed713f0fd1fbf90b2aaf4df497548c4f82c58 | 16,229,677 |
| regional topology composition | research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json | b71c663cf7b31c4e0fd18ad6fef9f878cec7ae6763f0438a8cd14de914150b67 | 64,578 |
| read-AF ordering | research/20260711_read_group_C_tree_T_topology_report/data/read_af_tree_ordering_historical.json | 602be0684b797a4d77ccda4af5e008d5e3610b9b1a05769fb676ad20a5bad516 | 85,857 |
| HCC1395 VAF example | research/20260711_read_group_C_tree_T_topology_report/data/HCC1395_VAF_multi_T_example.json | 694f615f256c5fac60d62ca82f2547d73944063ab2b5ba7df60a0a65b255b207 | 1,845 |

## 可重現命令與版本

Repo：branch `research/subclonal-reconstruction-202606`；commit `6067568637088838a9f518955e41d222f057e4f1`；worktree dirty=`True`。本任務只新增自己的 report workspace 與 in-progress report，不覆寫既有使用者變更。

```bash
python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_exact_topology_catalog.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 \
  --output research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json

python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_regional_topology_composition.py \
  --topology research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json \
  --output research/20260711_read_group_C_tree_T_topology_report/data/regional_topology_composition.json

python3 research/20260711_read_group_C_tree_T_topology_report/scripts/build_report_dataset.py [...]
```

Historical engineering snapshot 內部守恆：`PASS`，failures=`0`；clean-v3 scientific validation 仍 pending。

---

PARTIAL — historical 7-row engineering snapshot；clean layered-v3 7/7 pending。本文件不可作 final biological validation evidence。
