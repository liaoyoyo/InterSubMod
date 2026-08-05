<!--
建立時間: 2026-07-12
目標: 將 separate PyClone-VI mutation clusters 接回固定 5,720 regions，量化區域子集合再現與 read-edge CP 相容性
處理範圍: HCC1395 vs HCC1395_DORADO；chr1-22；fixed exact-coordinate complete-both regions
關聯檔案: InterSubMod/research/20260712_vaf_ccf_subclone_crosssoftware_validation/results/clone_region_bridge_v1/
-->

# Clone ↔ region bridge（v1）

> **PASS；但只支持 conditional cluster-to-region reproducibility，不是 tree truth。** Primary 使用兩側獨立 separate fits；joint labels 未列為獨立再現證據。

## Join coverage

- Fixed regions：5,720。
- Shared exact regional alleles：15,713；raw VAF joined 15,713；兩側 PyClone joined 14,369。
- 兩側 assignment ≥0.8：10,965 regional alleles。
- Region×mutation key unique、PyClone keys unique：PASS。Raw VAF 的 153,569 source rows 中，73,695 列是 payload 完全相同、僅 source label 不同的重複證據列；依明示規則收斂為 79,874 unique mutations，before/after rows 守恆 PASS。

## Global vs regional mutation-cluster concordance

| Population | n | ARI | NMI | Subclonal Jaccard | κ | CP Spearman |
|---|---:|---:|---:|---:|---:|---:|
| Global separate-fit universe | 30,262 | 0.539 | 0.339 | 0.381 | 0.544 | 0.547 |
| Fixed-regional joined subset | 14,369 | 0.375 | 0.195 | 0.239 | 0.379 | 0.380 |
| Regional both assignment≥0.8 | 10,965 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

`assignment≥0.8` 的完美值不可單獨解讀成 clone 再現：10,965 個 regional high-confidence mutations 的 subclonal union 只有 1（0.0091%），幾乎全被兩側主 clonal cluster 支配，屬 **vacuous selection-induced concordance**；不可稱「高信心區域 clone 高度一致」。

## k>1 multi-locus partition pattern

- All-joined k>1：fixed stratum 5,720；evaluable 5,189；partition exact 5,028（96.90% of evaluable）。
- 其中 both-single-cluster 5,007（96.49%）；真正兩側皆 multi-cluster 且 partition exact 只有 21（0.40%）。因此 96.90% 不可稱為「多 clone 結構 96.90% 一致」。
- 限定兩側都真的含 >1 cluster 的 informative denominator 是 34：exact 21、different 13，即 61.76%（21/34）。這才是較接近「區域內多 clone partition 是否一致」的條件式數字，但 n=34 很小。
- 非 exact 分為 one-single/one-multi 148 與 both-multi-different 13；完整列在 `region_cluster_patterns.tsv.gz`。

## Read-directed edge CP compatibility

- Edge contract已從 source驗證：F=低座標祖先、R=高座標祖先、P=平行；只有 singleton F/R 進 directed endpoint。
- CP判斷：同群為 uninformative；不同群用 CP±1.96×std 與 tolerance 0.02，清楚滿足才 compatible、清楚反向才 conflict、區間重疊則 uninformative。
- HCC1395：directed+joined 1,117/8,096；compatible 19、conflict 0、same-cluster uninformative 1,098。DORADO：directed+joined 514/8,096；compatible 12、conflict 0、same-cluster uninformative 502。
- assignment≥0.8 後，HCC1395 的 975 與 DORADO 的 332 條 evaluable edges 全部是 same-cluster；compatible/conflict 皆為 0。因此 high-confidence endpoint **沒有不同 cluster 的 CP ordering 資訊**，0 conflict 不能解讀成已驗證 ancestry。
- 兩側都有 determinate direction 的只有 598/8,096 pairs；其中 same direction 596/598，opposite 2。這是高條件化分母，不能外推成全 5,720 regions 的方向一致率。
- 這是 read-pattern edge 對 conditional PyClone CCF 的診斷；不是對真實 ancestry 的證明。VAF-selected edge 未拿來當獨立驗證，以避免同 VAF 循環。

## Claim ceiling

本橋接可回答「固定區域中的哪些 shared mutations 被兩側 separate PyClone fit 接上、區域內 co-clustering partition 是否再現、read-compatible directed constraint 是否與 CCF ordering 明顯衝突」。它不能把 PyClone clusters 升格為 clone truth，也不能證明每區唯一演化樹。Historical layered-v2 region artifacts仍受 upstream engineering-snapshot ceiling。
