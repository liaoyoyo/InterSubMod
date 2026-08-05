<!--
建立時間: 2026-07-27
更新時間: 2026-07-27（read 聚集圖 v06）
目標: 回答同一 exact raw HP 下多 sSNV R/A pattern 是否具有可重現的 read-level 區域甲基差異
處理範圍: Task B；7 technical datasets / 6 biological samples / chr1-22 / exact PS × exact raw HP
關聯檔案:
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/00_INDEX.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/analysis-contract.md
  - InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/implementation-notes.md
  - InterSubMod/output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/
-->

用 SCQA + Evidence Chain：先直接回答「有沒有」，再交代母體、正式檢定、位點證據、拓撲關係與不可越過的解讀邊界。

# 多 sSNV pattern × read-level 甲基關聯全面驗證

> **TL;DR：找到 3 個較長 R/A signature 的 robust regional methylation
> association，但沒有任何 exact raw-HP full-four RR/RA/AR/AA 單元，也沒有
> robust 的 exact 二位點 RA 對照；因此可說「某些多點 pattern 的 reads
> 具有明顯甲基差異」，不可說已辨識 clone、mutation ancestry 或 lineage
> 方向。**（影響：高；統計關聯信心：高；生物 lineage 信心：不足）

互動交付：

- [完整互動 HTML（含 read 聚集圖）](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html)
- [HTML 對應資料 JSON](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_data_06.json)
- [瀏覽器 QA receipt](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/browser_qa.json)
- [Builder-ready sidecar receipt](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/sidecar/pattern_methyl_sidecar.v1.receipt.json)

![完整驗證報告桌機總覽](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/01_desktop_overview.png)

## 1. 直接回答使用者問題

### 1.1 同一 HP 內有沒有 RR、RA、AR、AA 或更長組合？

有。所有正式比較都固定在同一
`dataset × chromosome × exact PS × exact raw HP × region`，沒有把
`1`、`1-1`、`2`、`2-1` 合併。正式 `n≥5` 母體有 1,045 個單元：

- 883 個二位點單元。
- 162 個 `k≥3` 多點單元。
- `k=2/3/4/5/6/9` 分別為 `883/120/31/6/4/1`。
- 正式單元只落在 raw HP `1-1` 493 個與 `2-1` 552 個。
- **同一 exact raw HP 下同時具 RR、RA、AR、AA 的 full-four 單元為 0。**

因此，資料確實有 RA、AR、AA、RR 與更長組合，但沒有達到預註冊正式門檻的
full-four 單元可一次比較四態。

### 1.2 有沒有「某 pattern 彼此相似，但和其他 pattern 不同」？

**有，而且最清楚的例子是 H1437 chr22 的較長 signature。**

同一 raw HP `1-1`、同一 PS `10691640` 內：

- `AARR` 與 `RRAR` 最相似：profile Pearson `r=0.815`、RMSE `0.205`、
  Bernoulli distance contrast `0.0179`。
- `RRAR` 與 `RRRA` 最不同：`r=0.590`、RMSE `0.387`、
  distance contrast `0.1398`、standardized effect `1.508`。
- `AARR` 與 `RRRA` 居中：`r=0.790`、RMSE `0.281`、
  distance contrast `0.0649`。

這正是「一組 read pattern 的甲基輪廓相近，另一組明顯不同」的案例；但
三組 pair 的 Hamming distance 都大於 1，所以只可畫無方向 **pair band**，
不可把差異畫成突變先後箭頭。

### 1.3 RA 本身有沒有正式 robust 結果？

以 exact 二位點、可評估單元計算，RA 對照的 membership 為：

| 對照 | 可評估單元 | Robust |
|---|---:|---:|
| RA + AA | 266 | 0 |
| RA + AR | 186 | 0 |
| RA + RR | 3 | 0 |

這些 membership 可以重疊，因為一個單元可含兩個以上 complete states。
**沒有任何 exact 二位點 RA 對照通過 robust claim。**

最強二位點訊號是 HCC1395_DORADO chr2 的 `AR` 對 `RA`：
`R²=0.7914`、distance contrast `0.3431`、standardized effect `1.697`；
但 secondary family 只有 199 permutations，`p=0.005` 位於解析度下限，
BY `q=0.0642`、Holm `p=1`，所以正式分類仍是
`EVALUABLE_NO_ROBUST_ASSOCIATION`。

![最強 secondary AR 對 RA 案例](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/03_secondary_exploratory_detail.png)

## 2. 分析母體與資料品質

### 2.1 全範圍 census

| 項目 | 數量 |
|---|---:|
| Technical datasets | 7 |
| Biological samples | 6 |
| Chromosomes | chr1–chr22 |
| Sparse source rows | 21,601,606 |
| Candidate read projections | 8,893,098 |
| Exact raw-HP pattern strata | 154,132 |
| Formal `n≥5` units | 1,045 |
| Formal `n≥8` units | 957 |
| Formal `n≥10` units | 910 |
| Formal active markers | 2,313 |

HCC1395 與 HCC1395_DORADO 是同一 biological sample 的 technical
datasets，跨樣本分母不把它們算成兩個生物樣本。

各 technical dataset 正式單元：

| Dataset | Formal units |
|---|---:|
| COLO829 | 0 |
| H1437 | 129 |
| H2009 | 436 |
| HCC1395 | 80 |
| HCC1395_DORADO | 122 |
| HCC1937 | 236 |
| HCC1954 | 42 |

COLO829 的 0 表示沒有單元通過本次完整 R/A、read 數與 state support gate，
不是證明 COLO829 全基因組沒有任何甲基異質性。

### 2.2 Bernoulli 舊圖計算的 parity

本輪沿用舊 `BERNOULLI/` 圖的 read-pair Bernoulli distance 定義，但重新從
hash-bound long-form CpG 資料計算，不直接拼貼歷史 matrix：

- 2,313 / 2,313 formal markers PASS。
- 277,560 個抽查 pair cells。
- invalid mask mismatch = 0。
- 最大絕對誤差 `1.0538×10⁻⁴`，小於凍結門檻 `5×10⁻⁴`。

所以新 HTML 的 heatmap、state-pair distance 與舊 Bernoulli 圖採同一數學
語意，同時修正了舊圖 raw HP 標籤可能過期的問題。

## 3. 正式統計流程

### 3.1 候選與甲基完全分離

候選只由 sSNV、同一 molecule 的 R/A/X calls、exact PS/raw HP 與 coverage
決定；沒有先看甲基 effect 再挑位點。含 X 的 pattern 保留為 partial
subcube evidence，不填補成完整 vertex，也不進正式 state test。

正式 gate：

1. active `k≥2`。
2. complete R/A `N≥40`。
3. 至少兩個 complete states 各 `n≥5`。
4. unit 至少 10 個共同 CpG；state/read CpG coverage 均至少 80%。
5. read-pair 至少 3 個共同 CpG。
6. restricted permutation 只在 exact `read_group × strand` strata 內交換；
   無法交換即 `NOT_EVALUABLE`，不退回 unrestricted test。

### 3.2 多重檢定與 robust gate

- Confirmatory family：full-four 或 `k≥3`，共 162 個；初跑 999
  permutations。
- Secondary family：其餘二位點，共 883 個；199 permutations，只作完整
  census，不升格 robust claim。
- 只有三個預先通過 effect、dispersion、geometry、`n≥8`、等 N 與 CpG
  rarefaction gate 的 confirmatory 單元提高到 49,999 permutations。
- BY 與 Holm 都在各自 frozen family 內重算。

Robust association 必須同時通過：
`q_BY≤0.05`、`R²≥0.10`、distance contrast `≥0.10`、
standardized effect `≥0.5`、PERMDISP `p≥0.05`、
geometry max SMD `<0.5`、所有 state `n≥8`、
equal-N 與 rare-CpG retention `≥0.5`，以及 distal sensitivity。

## 4. 三個 robust regional methylation associations

| Dataset / locus / exact HP | States（analysis n） | CpGs | R² | BY q | Best pair | Contrast / std. effect | Topology relation |
|---|---|---:|---:|---:|---|---|---|
| H1437 chr22, PS 10691640, HP 1-1 | AARR 19 / RRAR 9 / RRRA 13 | 131 | 0.3338 | 0.00498 | RRAR / RRRA | 0.1398 / 1.508 | Hamming 2 pair band |
| H2009 chr3, PS 71518648, HP 1-1 | AAA 17 / RAA 27 | 44 | 0.3072 | 0.00498 | AAA / RAA | 0.1217 / 1.001 | Hamming 1，非 unanimous edge |
| HCC1937 chr10, PS 25994308, HP 2-1 | AAR 30 / RRA 18 | 169 | 0.4064 | 0.00498 | AAR / RRA | 0.1124 / 1.560 | Hamming 3 pair band |

三者 raw `p=2×10⁻⁵`、Holm `p=0.00272`，且 PERMDISP
`p=0.314/0.388/0.731`。equal-N retention
`1.105/1.073/1.069`、rare-CpG retention
`1.033/1.011/0.982`、distal retention
`1.057/0.998/1.094`，沒有因平衡 read 數、稀有 CpG 或移除近 marker CpG
而消失。

### 4.1 H1437 chr22：最符合「部分相似、部分不同」

- Active markers：`10,695,463 / 10,705,400 / 10,710,397 / 10,714,112`。
- 全 131 common CpG 的 state mean profile：
  `AARR 0.652`、`RRAR 0.716`、`RRRA 0.493`。
- Median：
  `AARR 0.789`、`RRAR 0.888`、`RRRA 0.509`。
- `AARR` 與 `RRAR` 的輪廓最接近；`RRAR` 與 `RRRA` 的
  between-minus-within Bernoulli separation 最大。

### 4.2 H2009 chr3：RAA 與 AAA 的單一 bit 差異

- Active markers：`72,043,068 / 72,043,379 / 72,044,016`。
- 輸入 `AAA 24 / RAA 35`；共同 basis 收斂後為 `AAA 17 / RAA 27`。
- 全 44 CpG mean / median：
  `AAA 0.399 / 0.252`，`RAA 0.586 / 0.689`。
- 這是 Hamming-1 關係，但 frozen candidate topology 對此 edge
  **不是 unanimous**，所以圖上不加 edge halo。

### 4.3 HCC1937 chr10：AAR 與 RRA 的多 bit 分離

- Active markers：`26,248,503 / 26,248,504 / 26,255,345`。
- 輸入 `AAR 30 / RRA 21`；分析為 `AAR 30 / RRA 18`。
- 全 169 CpG mean / median：
  `AAR 0.635 / 0.801`，`RRA 0.558 / 0.707`。
- Hamming 3，只能作 pair band，不可解讀成三步演化路徑。

### 4.4 舊 BERNOULLI 樣式的 read 聚集關係

互動 HTML 已在每個可用位點最前方加入：

1. label-independent UPGMA average-linkage dendrogram；
2. 同一排序的對稱 read × read Bernoulli distance heatmap；
3. 排序完成後才疊加的完整 R/A pattern、strand 與匿名 read-group strips；
4. 匿名 `r001` 型 read pair hover 與鍵盤方向鍵讀值。

聚集排序只使用 read distance，**不使用 pattern、strand 或 RG 標籤**。
因此色條若形成連續區塊，是排序後才看見的關係，不是先按 RA/AR 分組畫圖。
exact-distance ties 使用 source leaf ordinal 作 deterministic tie-break，因此
方法是 label-independent，但不宣稱任意 row permutation 下完全 invariant。
淺色代表較相似、深色代表較不同；圖不設定 cluster cut，也不推論 clone、
祖源或演化方向。

557/1,045 個 formal units 有合法原始 read matrix，共 29,316 reads、
1,648,104 個 matrix cells。其餘 488 個明示不可用：

- 234 個沒有通過正式可評估條件；
- 253 個可評估但未達 source detail trigger
  （primary `p>0.05` 且 `R²<0.05`），因此原始矩陣沒有被保存；
- 1 個 H2009 chr5 單元有 169 reads，超過既定 `N≤160` matrix 保存上限。

四個代表案例的圖形描述量如下：

| Case | Reads | Within-pattern mean | Between-pattern mean | Difference | 排序後 pattern blocks |
|---|---:|---:|---:|---:|---:|
| H1437 chr22，AARR/RRAR/RRRA | 41 | 0.236 | 0.302 | 0.066 | 11 |
| H2009 chr3，AAA/RAA | 44 | 0.300 | 0.421 | 0.122 | 4 |
| HCC1937 chr10，AAR/RRA | 48 | 0.189 | 0.301 | 0.112 | 2 |
| HCC1395_DORADO chr2，AR/RA | 42 | 0.126 | 0.470 | 0.343 | 2 |

這也顯示「視覺分塊」與「正式 robust gate」不能混為一談：H1437 是正式
robust association，但 read order 仍有 11 個 pattern blocks；反之
HCC1395_DORADO `AR/RA` 形成兩個清楚區塊，`R²=0.791`，但
secondary BY `q=0.0642`，仍標為 **evaluable、not robust**。

![H1437 robust 案例的 read 聚集圖與匿名 hover](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/06_desktop_read_cluster_hover.png)

![HCC1395_DORADO AR/RA exploratory 案例的兩個 read blocks](../../output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/qa_v06/03_secondary_exploratory_detail.png)

## 5. 全部 1,045 單元的結果分類

| Assessment | Units | 正確解讀 |
|---|---:|---|
| ROBUST_ASSOCIATION | 3 | 通過完整預註冊 gate 的區域甲基關聯 |
| EVALUABLE_NO_ROBUST_ASSOCIATION | 627 | 可評估，但未通過 robust claim |
| CONFOUNDED | 181 | 至少一個 dispersion/geometry/robustness gate 失敗；**不是已證明混雜造成結果** |
| NOT_EVALUABLE | 234 | 正式 support、共同 CpG 或 exchangeability 不足 |

234 個不可評估單元的原因：

- post-CpG-filter formal state support lost：169。
- post-exchangeability formal state support lost：43。
- insufficient common CpG：14。
- post-join complete support不足：8。

沒有產生 `TAG_DEPENDENT`。這不代表 tag effect 已被證明為零，只代表沒有單元
符合該 frozen reader-facing 類別。

## 6. 指定位點回查

### H2009 chr5:18,096,980，HP 2-1

找到 active marker-containing stratum，但 complete `N=26`：
`AAAA 1 / AAAR 4 / AARA 19 / AARR 2`。未達 `N≥40` 與每態 `n≥5`，
所以不在正式 1,045-unit evidence universe，不能從本輪宣稱陽性或陰性。

### HCC1395 chr22:46,257,699，HP 1-1

找到 `AARR 18 / RRAA 3`，complete `N=21`。同樣未達正式門檻。
HTML 將兩個 sentinel 明確標為 `NOT_IN_FORMAL_EVIDENCE`。

## 7. 多點關係與跨樣本解讀

本輪把每個 state pattern 標到既有 frozen topology，但 methylation 只作
association overlay：

- Hamming 1 且 global-best candidate family unanimous 支持時，才可加
  edge halo。
- Hamming >1 一律畫 pair band。
- 不改 selected tree、AF、edge incidence、support 或方向。
- 不使用箭頭推斷 parent/child。

Sidecar 包含 1,045 個 exact raw-HP bundles、513 個合法 edge halos、
288 個 Hamming>1 pair bands、24 個大型矩陣；7 個 topology authority
hash 均綁定，topology counts/AF/selected tree 未更動。

三個 robust 結果分別出現在 H1437、H2009、HCC1937，顯示這種關聯不是只在
單一 technical dataset 才能觀察；但是它們位於不同染色體、PS、markers與
patterns，**不是同一 locus/signature 的跨樣本重現**。預註冊 H4 要求至少
3/6 biological samples 的同方向可比 finding，因此 H4 本輪不成立。

## 8. 假說判定與證據層級

| 假說 | Verdict | 證據層級 |
|---|---|---|
| H1：exact raw-HP full-four 可解釋甲基結構 | Not supported；formal full-four = 0 | L2 統計母體 |
| H2：至少一個 state 與其他 state 呈穩健分離 | Supported，但證據來自較長 signatures；exact 二位點 RA 無 robust | L2 統計關聯 |
| H3：`k≥3` 保留 pairwise 以外的多點關聯資訊 | Supported at association level；3/162 confirmatory units robust | L2 統計關聯 |
| H4：sample-level finding 跨生物樣本同方向重現 | Not supported | L3 跨樣本不足 |
| Clone、ancestry、causal lineage | 未檢定且禁止宣稱 | L5／無直接證據 |

## 9. 驗證與 QA

- 研究 Python tests：88/88 PASS。
- Combined-run merge tests：21/21 PASS。
- C++ `run_tests`：258/258 PASS。
- Report-specific tests：18/18 PASS。
- Artifact catalog：2,313/2,313 PASS。
- Bernoulli parity：2,313/2,313 PASS。
- Sidecar receipt：PASS；所有 identity、state、pair、topology relation
  fail-closed checks 通過。
- 瀏覽器 QA：初始 1,045 rows、robust filter 3 rows、COLO829 0 rows、
  全域/篩選指標無歧義、桌機/390 px mobile/print 無 global overflow、
  console/page errors = 0；read-cluster matrix byte count、實際 canvas pixels、
  匿名 hover、鍵盤讀值與手機局部橫向捲動均 PASS。
- HTML 與外部 JSON 完全一致。
- 獨立 final closeout audit：P0=0、P1=0、P2=0。

Final hashes：

- HTML SHA-256：
  `cdc2e55165f05c430177d02993ee52535919ba4bda99a7036c23e31961361b1d`
- Data JSON SHA-256：
  `98a980c9fa28ed9859fe9656f3ab089ccfddd7eb907bd3580947ae375d56e3fc`
- Sidecar SHA-256：
  `45488a5fd2bdbbe5155a38aeaea3faa19ed02fa49a15be3fdb6139b85bf243f1`

## 10. 實際執行與輸出

主要輸入：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/census/pattern_census.receipt.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/formal/artifact_catalog.v1.tsv
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/formal/bernoulli_parity.v1.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/pattern_methyl_evidence.combined.v1.tsv.gz
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/pattern_methyl_details.combined.v1.jsonl.gz
```

正式統計參數、三個 family run 的 input hashes、merge 規則與輸出 hashes
完整記錄在：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/analysis_summary.combined.json
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/merge_receipt.json
```

最後報告重建命令：

```bash
PYTHONDONTWRITEBYTECODE=1 python \
  research/20260727_multisite_pattern_methyl_topology_validation/scripts/build_pattern_methyl_report.py \
  --census-receipt /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/census/pattern_census.receipt.json \
  --pattern-counts /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/census/pattern_counts.tsv \
  --artifact-catalog-summary /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/formal/artifact_catalog.v1.json \
  --artifact-catalog-tsv /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/formal/artifact_catalog.v1.tsv \
  --bernoulli-parity-summary /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/formal/bernoulli_parity.v1.json \
  --evidence /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/pattern_methyl_evidence.combined.v1.tsv.gz \
  --details /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/pattern_methyl_details.combined.v1.jsonl.gz \
  --analysis-summary /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/analysis/combined/analysis_summary.combined.json \
  --output-html /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_06.html \
  --output-data-json /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/20260727_multisite_pattern_methyl_association_data_06.json
```

最後輸出：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/report/
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260727_multisite_pattern_methyl_topology_validation/all7_v1/sidecar/
InterSubMod/research/20260727_multisite_pattern_methyl_topology_validation/20260727_多sSNV_pattern與甲基關聯全面驗證_01.md
```

## 11. 結論與下一個合理實驗

本輪已清楚找出三個「同一 exact raw HP、不同較長 R/A signature 的 reads
具有區域甲基差異」位點，並用 Bernoulli-style profile、pair distance matrix
與無方向 multi-point topology relation 呈現。最強 exact 二位點 `AR/RA`
訊號雖大，仍未過 secondary BY 門檻；full-four 與跨樣本同 locus replication
也沒有成立。

若要升級成 lineage 或 clone claim，下一步不是放寬本輪門檻，而是對這三個
位點做 marker-held-out 或 orthogonal phasing、第二個獨立生物樣本的同 locus
驗證，以及不以同一 sSNV labels 定義又驗證群組的外部證據。

---

Claim ceiling：本報告只支持 pattern-conditioned regional methylation
association；不支持 cellular subclone、mutation ancestry、因果方向或
methylation-rescored topology。
