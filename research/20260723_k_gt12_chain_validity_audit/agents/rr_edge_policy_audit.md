<!--
建立時間: 2026-07-23
目標: 獨立只讀驗證 strict graph 的 RR-only edge 語意，並重算 HCC1395 在 ALT-informative edge 規則下的 W connectivity 敏感度
處理範圍: HCC1395 chr1-22；exact nonmissing PS × HP1/HP2；primary threshold=3
關聯檔案:
  - InterSubMod/tools/strict_endpoint_graph.py
  - InterSubMod/scripts/build_strict_ps_hp_regions.py
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage生產流程驗證_01.md
  - InterSubMod/research/20260723_k_gt12_chain_validity_audit/agents/rr_edge_policy_recompute.py
狀態: VALIDATED READ-ONLY SENSITIVITY；未修改 production
-->

# RR-only edge 語意與 HCC1395 connectivity 敏感度獨立稽核

**TL;DR：production 的 RR-only edge 在「共同可判讀的 mutation-state 區域」語意下合理，但不能稱為 somatic ALT 共現；HCC1395 有 18,116／76,202（23.77%）條 primary edge 是 RR-only，移除後會影響 2,437／11,462（21.26%）個 source W，因此論文應保留 RR 作 primary callability/state graph，同時強制報告 ALT-informative sensitivity，不宜把兩種語意混成同一個「突變共現區域」。（影響：高；信心：高）**

用 SCQA：

- **Situation**：現行 strict graph 以同一 canonical molecule 對兩 endpoint 的 fixed `R/A` call 建 edge。
- **Complication**：`RR`、`RA`、`AR`、`AA` 都會累加到 `support_total`；因此三條全為 `RR` 的 edge 也可連接 W。
- **Question**：RR-only 是否應保留？若只允許有 ALT 證據的 edge，HCC1395 的 W 會改變多少？
- **Answer**：RR 是合法的共同 callability／負向 state evidence，但不是 ALT 共現。它對 connectivity 的影響不可忽略，應採「primary 全四態 graph＋ALT-informed sensitivity／標記」的雙層敘述。

本任務分類為 **(B) Comprehensive validation**，範圍是 HCC1395 chr1–22 全部現有 strict edge 與 membership artifact；服務 **G1、G4、G5**。

## 1. 驗證結論

### 1.1 程式契約

現行程式的 edge 語意可明確寫成：

\[
s_{ij}=n_{RR}+n_{RA}+n_{AR}+n_{AA}
\]

\[
e_{ij}^{primary}=1[s_{ij}\ge 3]
\]

其中每一票來自一個 distinct canonical molecule，且兩個 endpoint 都必須是 fixed `R` 或 `A`。程式並沒有要求 `RA+AR+AA>0`：

- `InterSubMod/tools/strict_endpoint_graph.py:97-111` 強制 `total = RR+RA+AR+AA`。
- `InterSubMod/tools/strict_endpoint_graph.py:162-171` 對同一 molecule 的所有 fixed endpoint pair 累加四態之一。
- `InterSubMod/tools/strict_endpoint_graph.py:274-295` 只以 `row.total >= threshold` union endpoints。
- `InterSubMod/scripts/build_strict_ps_hp_regions.py:370-393` 也以 `row.total >= primary_threshold` 寫入 primary flag。

因此：

| edge 類型 | 正確定義 | 可以支持 | 不可以直接支持 |
|---|---|---|---|
| `RR-only` | qualifying edge 且 `RA+AR+AA=0` | 兩個候選位點在同一批 molecule 上共同 callable，且觀察到 reference/reference state | somatic ALT 共現、同一 clone、祖先順序 |
| `ALT-any` | primary qualifying 且 `RA+AR+AA≥1` | 至少一個 endpoint 有 ALT 的 pair-state evidence | 兩個 ALT 一定共現；其中 `RA/AR` 其實是一 ALT、一 REF |
| `ALT-support≥3` | `RA+AR+AA≥3` | 至少三個 molecules 帶有一個或兩個 ALT endpoint | 唯一 clone／真實演化邊 |
| `AA` | molecule 在兩 endpoint 都為 ALT | 直接的 ALT–ALT 同分子共現 | 同一細胞 clone（bulk 中仍缺 cellular identity） |

`RA/AR/AA` 合稱 ALT-informative 可以用於 sensitivity，但若口試講「兩個突變共現」，只有 `AA` 符合該字面；`RA/AR` 是 mutation-state constraint，不是雙 ALT 共現。

### 1.2 HCC1395 edge 守恆

輸入粒度是一列 `(chromosome, exact PS, HP, site_i, site_j)` edge；support 粒度是 molecule–endpoint-pair vote。22／22 chromosome receipts 均 `all_pass=true`，本次實讀的 edge／membership path 與 SHA-256 mismatch 均為 0。全 117,760 列 observed edges 的逐列守恆錯誤為 0：

\[
1{,}133{,}489
=617{,}321_{RR}+97{,}893_{RA}+95{,}044_{AR}+323{,}231_{AA}
\]

primary 76,202 edges 的 support mass 亦守恆：

\[
1{,}077{,}043
=592{,}507_{RR}+93{,}350_{RA}+91{,}057_{AR}+300{,}129_{AA}
\]

注意：這是 **edge-support mass**，同一 molecule 可對多個 endpoint pair 投票，所以不可當成 unique read 數或 biological prevalence。

primary edge 的 evidence class：

| 類別 | edge 數 | primary edge 比例 |
|---|---:|---:|
| RR-only | 18,116 | 23.7737% |
| 至少 1 個 ALT-informative molecule | 58,086 | 76.2263% |
| ALT-informative support ≥3 | 48,292 | 63.3737% |
| **primary total support ≥3** | **76,202** | **100%** |

## 2. Connectivity 反事實重建

三個 policy 都固定使用相同的 62,651 個 exact-container node memberships；只改 edge 納入規則：

1. `primary_total_ge3`：`RR+RA+AR+AA≥3`。
2. `alt_any_with_total_ge3`：primary 合格，且 `RA+AR+AA≥1`。
3. `alt_support_ge3`：`RA+AR+AA≥3`。

重建 primary graph 完整復現正式 summary：7,118 containers、76,202 edges、28,384 singletons、11,462 tree-eligible W、34,267 linked memberships、90 個 k>12、max k=153。這表示 sensitivity 的 container/node grain 與正式輸出一致。

| graph policy | retained edges | singleton | tree-eligible W | W 相對 primary | linked memberships | k>12 | max k |
|---|---:|---:|---:|---:|---:|---:|---:|
| primary：total≥3 | 76,202 | 28,384 | 11,462 | baseline | 34,267 | 90 | 153 |
| ALT-any：total≥3 且 ALT≥1 | 58,086 | 33,604 | 9,734 | −1,728（−15.08%） | 29,047（−15.23%） | 80 | 153 |
| ALT support≥3 | 48,292 | 35,783 | 9,285 | −2,177（−18.99%） | 26,868（−21.59%） | 69 | 148 |

W 總數下降不是完整影響指標，因為移除 edge 後一個 source W 可能分裂成多個 W。以原 11,462 個 source W 逐一對應：

| 影響類別 | ALT-any | ALT support≥3 |
|---|---:|---:|
| 完全不變 | 9,025 | 7,929 |
| **membership／partition 有變** | **2,437（21.26%）** | **3,533（30.82%）** |
| 全部退化為 singleton | 1,814 | 2,392 |
| 留一個較小 W＋部分 singleton | 538 | 943 |
| 分裂成兩個以上 W | 85 | 198 |
| 原 linked memberships 退為 singleton | 5,220 | 7,399 |

### 2.1 k>12 source W 的影響

RR edge 對大型 component 也不是可忽略的小尾端：

| 原 primary k>12 source W（n=90） | ALT-any | ALT support≥3 |
|---|---:|---:|
| 完全不變 | 62 | 32 |
| **受影響** | **28（31.11%）** | **58（64.44%）** |
| 全部退為 singleton | 4 | 8 |
| 留一個較小 W＋部分 singleton | 14 | 33 |
| 分裂成兩個以上 W | 10 | 17 |
| 至少保留一個 k>12 子 component 的 source W | 77 | 67 |
| alternative graph 的 k>12 component 數 | 80 | 69 |

這證明 `k>12` 不只受 computational cap 影響，也對 edge evidence policy 敏感。論文若只呈現 primary 的 90 個而沒有 ALT-informed sensitivity，容易讓讀者誤以為 90 個大型區域都含穩健 somatic ALT linkage。

## 3. 方法判斷：保留 RR，但必須分層

### 3.1 為何 primary 可以保留 RR

本方法後續要比較完整與 partial mutation-state patterns。`RR` 並非「沒有資訊」；它是兩個候選位點在同一 molecule 上都可判讀，且兩者皆為 reference 的聯合負向證據。若直接刪除，會：

- 對低 VAF／低 purity 的 mutation 系統性降低連通；
- 只保留看到 ALT 的 molecule，形成 outcome-conditioned selection；
- 丟掉 likelihood 判別 candidate state 所需的 negative evidence。

所以若 `W` 的名稱是 **read-linked jointly callable mutation-state region**，保留 RR 合理。

### 3.2 為何不能只報 primary W

RR-only edge 完全不含 somatic ALT，不能拿來支撐「突變共現」。而且 23.77% primary edges 為 RR-only；移除後 21.26% source W membership／partition 改變，已超過可忽略範圍。故 primary 結果必須伴隨 ALT-informed sensitivity 或 edge-class annotation。

### 3.3 建議的論文與實作分層

不修改目前 production 的前提下，建議固定三層名稱：

| 層 | 建議名稱 | 用途 |
|---|---|---|
| `W_state` | read-linked mutation-state region | primary；四態都可建 edge，供 joint state likelihood |
| `W_ALT-any` | ALT-informed sensitivity region | 至少一個 ALT molecule；辨識是否只靠 RR callability 連接 |
| `W_ALT3` | conservative ALT-supported sensitivity region | 至少三個 ALT-informative molecules；檢查低支持 edge 敏感度 |

每個 primary W 至少新增／報告：

- `n_primary_edges`
- `n_RR_only_edges`
- `n_ALT_any_edges`
- `n_ALT_support_ge3_edges`
- `RR_bridge_critical`：移除 RR-only 後 membership 或 partition 是否改變
- `ALT_connectivity_class`：unchanged／partial／split／fully-unlinked

對口試最安全的說法：

> 我們先用同一 molecule 在兩個 endpoint 都有固定 R/A call 的證據，建立共同可判讀的 mutation-state graph；RR、RA、AR、AA 都是 state evidence。這個 graph 不是 ALT 共現圖，因此我們另外移除 RR-only edge，並用 ALT support≥1 與 ≥3 重建 sensitivity graph。HCC1395 有約 21% 的 source W 對 RR-only edge 敏感，所以後續 topology 必須標示 evidence class，不能把所有 W 都稱為 somatic mutation co-occurrence region。

不建議說：

> 三條 read 穿過兩個突變，就證明兩個突變在同一 clone。

建議改為：

> 至少三個 distinct molecules 在兩個候選 sSNV endpoint 都有固定 allele call，支持兩位點能在同一 read-level state model 中共同分析；是否含 ALT、是否雙 ALT 共現，以及是否屬於同一 cellular clone，仍是不同層次的證據。

## 4. Step → Verify 與可重現命令

### 4.1 程式與輸入

- 輸入：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2/`
- 重算程式：
  `InterSubMod/research/20260723_k_gt12_chain_validity_audit/agents/rr_edge_policy_recompute.py`
- production core：
  `InterSubMod/tools/strict_endpoint_graph.py`
- production builder：
  `InterSubMod/scripts/build_strict_ps_hp_regions.py`

1. 核對 edge 守恆與 primary flag  
   → 驗證：117,760 edge rows 中 `total != RR+RA+AR+AA` 為 0；`passes_primary_threshold != (total≥3)` 為 0。
2. 用相同 nodes 重建 primary components  
   → 驗證：76,202 edges、28,384 singleton、11,462 W、34,267 linked memberships，全部與正式 summary 相同。
3. 移除 RR-only 或提高 ALT support 後重建  
   → 驗證：逐 container 執行 union-find，alternative component 必為 primary component 的子集，不跨 chromosome／PS／HP。
4. 逐 source W 比對 membership set  
   → 驗證：分類為 unchanged／fully lost／partial／split，四類加總等於 11,462。

執行命令：

```bash
python3 \
  research/20260723_k_gt12_chain_validity_audit/agents/rr_edge_policy_recompute.py \
  --root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260723_production_exact_ps_strict_read_linkage/hcc1395_strict_regions_v2
```

實際輸出重點：

```text
edge_state_conservation.pass = true
input_receipt_checks.pass = true (22 receipts; path/SHA mismatch=0)
primary_RR_only_qualifying_edges = 18116
primary_ALT_informative_qualifying_edges = 58086
primary_total_ge3: edges=76202, W=11462, singleton=28384
alt_any_with_total_ge3: edges=58086, W=9734, singleton=33604
alt_support_ge3: edges=48292, W=9285, singleton=35783
```

production semantics tests：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -m pytest -q \
  tests/test_strict_endpoint_graph.py \
  tests/test_build_strict_ps_hp_regions.py
```

實際輸出：

```text
15 passed in 0.18s
```

## 5. 限制與下一個 gate

- 本次只重算 HCC1395；不能直接外推其他 technical datasets。
- `ALT≥1` 與 `ALT≥3` 是 policy sensitivity，不是經 error model／VAF／depth 校準後的最佳統計 threshold。
- `RA/AR/AA` 的 edge-class 只說至少一個 ALT；若主張「兩個 ALT 共現」，必須另報 `AA` support。
- 同一 molecule 可對多個 pair 投票，edge-support mass 不是獨立樣本數。
- bulk tumor 中即使 `AA` 也只證明同分子共現，不足以單獨證明 cellular clone identity。
- 本次沒有修改 production；若未來要把 `W_ALT-any` 升為新 primary，必須先完成 7 technical datasets × chr1–22 的全量 sensitivity、coverage/VAF/CN 分層與 preregistered decision gate。

**建議 verdict：保留現行 RR-inclusive primary graph，但把 region 正式命名為 read-linked mutation-state／callability region；將 ALT-any 與 ALT≥3 納入論文必要 sensitivity，而不是把 RR-inclusive W 直接稱為 somatic mutation co-occurrence region。**
