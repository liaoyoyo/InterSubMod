<!--
建立時間：2026-07-24
目標：精確區分 AF 最佳樹的唯一樹、同拓撲並列與跨拓撲並列
處理範圍：7 datasets、chr1-22、20260724 exact-PS×HP canonical ranked units
關聯檔案：InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/README.md
-->

# exact-PS×HP AF 最佳樹 topology-signature 精確統計

## 結論

既有 C++ JSONL 只有完整最佳樹數、minimum-family digest 與一棵
deterministic representative，不能直接回答同分樹是否屬於相同拓撲。
本分析以 frozen source、原始 MLHP input 與相同 guard 重算全部
`71,955` 個 `ranked_complete` units，精確展開 `680,527` 棵 global-best
trees。

新版結果為：

- `UNIQUE_TREE`：`39,648/71,955 = 55.1011%`。
- `TIED_SAME_TOPOLOGY`：`23,858/71,955 = 33.1568%`。
- `TIED_CROSS_TOPOLOGY`：`8,449/71,955 = 11.7421%`。
- 因此單一 exact rooted-unlabeled topology：
  `63,506/71,955 = 88.2579%`。

這是新版 exact-PS×HP 分母；不可把舊 legacy grouping 的 `92.18%`
直接移植到新版。

## 每組資料的精確拓撲解析度

百分比的共同分母都是該資料集的 `ranked_complete` units。

| Dataset | Ranked | Unique tree | Tied, same topology | Tied, cross topology | 單一 exact topology |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 9,130 | 7,047 (77.19%) | 1,861 (20.38%) | 222 (2.43%) | 8,908 (97.57%) |
| HCC1395_DORADO | 5,308 | 4,165 (78.47%) | 1,100 (20.72%) | 43 (0.81%) | 5,265 (99.19%) |
| COLO829 | 10,757 | 5,119 (47.59%) | 4,454 (41.41%) | 1,184 (11.01%) | 9,573 (88.99%) |
| H1437 | 13,740 | 6,369 (46.35%) | 5,477 (39.86%) | 1,894 (13.78%) | 11,846 (86.22%) |
| H2009 | 23,128 | 8,161 (35.29%) | 9,927 (42.92%) | 5,040 (21.79%) | 18,088 (78.21%) |
| HCC1937 | 4,245 | 3,672 (86.50%) | 547 (12.89%) | 26 (0.61%) | 4,219 (99.39%) |
| HCC1954 | 5,647 | 5,115 (90.58%) | 492 (8.71%) | 40 (0.71%) | 5,607 (99.29%) |
| **ALL7** | **71,955** | **39,648 (55.10%)** | **23,858 (33.16%)** | **8,449 (11.74%)** | **63,506 (88.26%)** |

## 四類 coarse graph geometry

定義：

- `Single-only`：沒有 branching，也沒有 root 之下兩層以上路徑。
- `Sister-only`：有 branching，但沒有兩層以上路徑。
- `Direct-only`：有兩層以上路徑，但沒有 branching。
- `Sister+direct`：同時有 branching 與兩層以上路徑。

這四類是 mathematical graph geometry，不等同於 cellular
single-clone、multi-clone 或 confirmed subclone。

| Dataset | Single-only | Sister-only | Direct-only | Sister+direct | 跨 coarse class 未解 |
|---|---:|---:|---:|---:|---:|
| HCC1395 | 3,146 (34.46%) | 466 (5.10%) | 4,766 (52.20%) | 533 (5.84%) | 219 (2.40%) |
| HCC1395_DORADO | 1,913 (36.04%) | 412 (7.76%) | 2,778 (52.34%) | 162 (3.05%) | 43 (0.81%) |
| COLO829 | 4,599 (42.75%) | 27 (0.25%) | 4,919 (45.73%) | 28 (0.26%) | 1,184 (11.01%) |
| H1437 | 4,094 (29.80%) | 129 (0.94%) | 7,305 (53.17%) | 318 (2.31%) | 1,894 (13.78%) |
| H2009 | 4,516 (19.53%) | 213 (0.92%) | 12,084 (52.25%) | 1,276 (5.52%) | 5,039 (21.79%) |
| HCC1937 | 1,223 (28.81%) | 449 (10.58%) | 2,261 (53.26%) | 287 (6.76%) | 25 (0.59%) |
| HCC1954 | 2,644 (46.82%) | 523 (9.26%) | 2,154 (38.14%) | 286 (5.06%) | 40 (0.71%) |
| **ALL7** | **22,135 (30.76%)** | **2,219 (3.08%)** | **36,267 (50.40%)** | **2,890 (4.02%)** | **8,444 (11.74%)** |

`63,511/71,955 = 88.2649%` 可唯一指派到一個 coarse class。
其中只有 5 units 是「exact rooted shape 超過一種，但四類 coarse
class 相同」；其餘 `8,444` 個跨 exact shape units 也跨 coarse class。

## Signature 與計算方法

每棵樹由 root 開始遞迴轉成括號字串；同一 parent 的 child signatures
排序後串接。因此 signature：

- 保留 root 與 parent-child geometry；
- 忽略 mutation/site label；
- 忽略 sibling 順序。

所以 `00→10→11` 和 `00→01→11` 都是同一條三節點 chain signature
`((()))`，雖然 acquisition label 不同。

對每個 minimum vertex set：

1. 以 exact rational AF 計算每個 child 的所有 Hamming-1 parents。
2. 保存該 child 的所有最佳 parents。
3. 展開各 child 最佳 parent choices 的 Cartesian product。
4. 收集每棵 global-best tree 的 canonical signature 與 coarse class。
5. 依 signature 數量分類 resolution。

本批每 unit 最多只有 276 棵 global-best trees，故直接精確展開可行，
不需近似或抽樣。

## 驗證

- 1-unit synthetic oracle：3 棵同分最佳樹，正確得到
  `Direct-only=2`、`Sister-only=1`、兩種 signatures。
- 全七逐列檢查：
  - minimum-family SHA-256 與 canonical 一致；
  - best AF score 一致；
  - best vertex-set count 一致；
  - factorized tie count = 實際展開樹數 =
    canonical `best_tree_tie_count`；
  - representative signature 必須存在於 exact census。
- `71,955/71,955` rows 全通過；7/7 processes exit 0。
- 每樣本 runtime 1.63–33.49 秒；最高 RSS 365,276 kB。

正式 machine-readable 證據：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.json`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/summary.tsv`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260724_exact_ps_cpp_topology_signature_census/all7_v1/receipt.v2.json`

## Claim ceiling

本分析證明的是「AF 最佳 mathematical tree 的 geometry 是否唯一」。
它不證明：

- 一個 tree node 就是一個細胞 clone；
- Sister 就是多 clone；
- Direct 就是已確認的 subclone ancestry；
- AF 最佳樹是 calibrated biological posterior；
- CNA、LOH、allele-specific amplification 已被排除。
