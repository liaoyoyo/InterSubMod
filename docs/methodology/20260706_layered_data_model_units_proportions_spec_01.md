<!--
建立: 2026-07-06
類型: 資料模型定義 spec(單位階層 + 比例分母字典 + 關係 + 最終格式)
用途: 分層 per-HP-家族樹重建的「數據檢視框架」SoT — HTML/報告所有數字的分子/分母定義依此
狀態: v1 草案(待使用者確認);數字來源逐項標;⚠標者待 Bash 恢復後 provenance 重算
build_branch: research/subclonal-reconstruction-202606
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/layered_reconstruction_HCC1395.json, sm_region_integration.json, sm_linkage_genomewide.json
-->

# 分層樹重建 — 資料模型定義（單位 · 階層 · 比例分母 · 關係 · 格式）

> **目的**：把「從 sSNV 比例 → 區域比例 → 樹/read 層次」的每一個數字，定義清楚**單位、分子、分母、關係**，讓每次檢視都能確認完整狀況、比例計算方式唯一、關係與格式明確。
> **鐵則（§13.0）**：每個比例都標「分子 ÷ 分母（分母是什麼層的單位）」；分母不同 → 數字不可直接比。

---

## §1 單位階層（6 層，細 → 粗 → 樹內）

| Lv | 單位 | 定義 | HCC1395 計數 | 來源 |
|----|------|------|------|------|
| **U1** | **sSNV** | somatic SNV（census `somatic==True`）= 重建骨幹 | **23,810**（✅ 實測 census somatic==True）| sm_linkage_genomewide.json census |
| U1-tot | census 位點（總）| 所有 union pileup 位點（somatic + germline-het + 其他）| 35,332 | 同上（🔴 handoff/master_spec 的「35,332 sSNV」實為**總 census 非 somatic**，待更正）|
| **U2** | **read / molecule** | ONT 單分子 read（覆蓋 ≥1 個目標 sSNV）| **read×region 覆蓋次 = 1,533,402**（⚠非 unique molecule；一 read 跨多 region 重複計）| layered json `L0.reads_by_family_total` 加總 |
| **U3** | **region**（🔴 **主分母**）| maximal linkage region（somatic-sSNV chain）；本 pipeline 計算單位 = multilocus 分析群（最密 ≤8 sSNV 窗）| **7,100**（multilocus 分析群，主分母）；7,143（integration linkage 全 span，同批區，見 §3）| sm_multilocus |
| **U3-attr** | **HP-multiplicity**（region 屬性）| 一個 region 有幾個 germline(1/2) lineage 有樹 | 0:104 · **1(single-HP):3004** · **2(multi-HP):3992** | layered detail（使用者定案 2026-07-06：先呈現多-HP 再定義）|
| **U4** | **lineage-unit**（次分母）| region × HP-family（1/2/3）= **一次重建的單位（=一組枚舉樹）** | **12,475** lineage（germline1/2:10,988 + somatic3:1,487）+ 4,598 none = 17,073 units | layered json `L1` |
| **U5** | **tree** | 一個 lineage-unit 枚舉出的最小樹之一（全等機率集）| Σ=45,011（non-capped）；**avg 3.78 樹/unit**；分布 1:7132·2-5:3899·6-20:1059·>20:385·max:125 | detail[].n_trees |
| **U6** | **node/隱藏祖先** | 樹節點：觀測 genotype ∪ 隱藏祖先（H_*=推斷未觀測 clone）| Σn_hidden=**17,319**；分布 0:4364·1:3934·2:1982·3:1073·4:549·tail(5-17,capped) | detail[].n_hidden |

> 🔴 **主分母 = U3 region（7,100）**（使用者定案 2026-07-06）。lineage-unit（U4）為次分母（樹的自然單位）。region-level 判定需先看 **HP-multiplicity**（多少區有多-HP 樹）再定 determinacy（§2）。

---

## §2 比例字典（每個標「分子 ÷ 分母」，主分母 = region U3）

### §2A region 層（主）— 先 HP-multiplicity，再 determinacy
| # | 比例名 | 分子 | 分母 | 值 |
|---|--------|------|------|-----|
| P1 | somatic sSNV 密度 | 23,810 somatic sSNV | 7,100 region | ≈3.4 sSNV/區 |
| **P2** | **multi-HP 區（germline 1&2 都有樹）** | 3,992 | 7,100 region | **56.2%** |
| P2b | single-HP 區（只 1 個 germline lineage）| 3,004 | 7,100 region | 42.3% |
| P2c | 無 germline lineage（只 3/none）| 104 | 7,100 region | 1.5% |
| **P3** | **region all-determined**（所有 germline lineage 都唯一樹）| 2,317 | 7,100 region | **32.6%** |
| P4 | region has-ambiguous（≥1 germline lineage 多樹）| 4,093 | 7,100 region | 57.6% |
| P5 | region has-capped | 543 | 7,100 region | 7.6% |
| P6 | region has-recurrence | 43 | 7,100 region | 0.6% |

> 🔴 **region-level determinacy 定義（✅ 使用者確認 2026-07-06）**：region「確定」⟺ 其**所有** germline(1/2) lineage 都 determined（多-HP 區要兩家族都唯一樹）。故 region all-determined 32.6% **嚴於** lineage-unit 55.2%（多-HP 區需雙確定）。多-HP 區的 (fam1,fam2) 組合見 §2C。

### §2B lineage-unit 層（次）
| # | 比例名 | 分子 | 分母 | 值 |
|---|--------|------|------|-----|
| P7 | determined（lineage）| 6,883 | 12,475 lineage-unit | 55.2% |
| P7b | ├ germline1/2 | 5,772 | 10,988 | 52.5% |
| P7c | └ somatic3 | 1,111 | 1,487 | 74.7% |
| P8 | ambiguous（lineage）| 4,976 | 12,475 | 39.9% |
| P9 | capped | 573 | 12,475 | 4.6% |
| P10 | recurrence | 43 | 12,475 | 0.34% |
| P11 | CN artifact（recurrence 內）| 29 | 43 recurrence | 67.4% |

### §2C 多-HP 區 determinacy 交叉（3,992 區；D=determined A=ambiguous C=capped R=recurrence）
| (fam1,fam2) | 數 | | (fam1,fam2) | 數 |
|---|---|---|---|---|
| D,D | 1,313 | | D,C / C,D | 99 / 91 |
| D,A / A,D | 980 / 961 | | C,A / A,C | 74 / 56 |
| A,A | 397 | | 其他(含R/C,C) | 21 |

> **絕不混分母**：P3（32.6%，分母 region）≠ P7（55.2%，分母 lineage-unit）≠ 舊「A_determined 28.1%」（分母舊 pooled 6,288，已淘汰）。三者單位不同、不可比。

---

## §3 分母漂移對照（釐清 5 個「區數」+ 7143 vs 7100 精確歸屬）

| 數字 | 是什麼 | 該當哪種分母 |
|------|--------|------------|
| **35,332** | somatic sSNV total | U1 sSNV 層 |
| **7,143** | integration linkage region（全 sSNV span，可達 94 sSNV/區）| U3 **linkage 完整版**（生物 region）|
| **7,100** | multilocus 分析群（最密 ≤8 sSNV 子窗，樹在此建）| U3 **主分母**（reconstruction 操作單位）|
| ~~6,288~~ | 舊 pooled「有向量」區 | ⚠**已淘汰**（單樹+混家族）|
| **12,475** | lineage-unit（region×家族 1/2/3）| U4 次分母 |

### 🔴 7,143 vs 7,100 歸屬（位置重疊實測，非缺區）
- **位置重疊 99.9%**：layered 7,100 有 **7,091** 落在某 integration 區內（未重疊僅 9）；integration 7,143 有 **7,141** 落在 layered（未重疊僅 2）。
- **同一批區、不同邊界定義**：integration 有 **178 個 n_sSNV>8 的大區** → multilocus 取「最密 8 sSNV 子窗」→ region-id 字串不同、基因座相同（故字串比對 6,836 重疊看似差很多，實為邊界標籤差異）。
- **哪個更完整**：7,143 = linkage 全 span（生物完整）；7,100 = 樹實建的操作窗。**reconstruction 分母固定用 7,100**；敘述提「linkage 全 span 7,143（含 178 大區）」並標二者為同批區、位置 99.9% 重疊。真正未配對 9+2 區 = 邊界邊角，逐一列於附錄。

> **去向全交代**：7,100 = 7,091(重疊 integration) + 9(邊界未配對)；7,143 = 7,141(重疊 layered) + 2(邊界未配對)；178 大區 = 同區不同窗。

---

## §4 關係表（層與層的映射 + 基數）

| 關係 | 基數 | 說明 |
|------|------|------|
| sSNV → region | many : 1 | 多 sSNV 組成一 region（≥2）|
| read → HP-family | 1 : 1 | 每 read 依 germline tag 歸一家族（1/2/3/none）|
| read → region | many : many | 一 read 可跨多 region；一 region 多 read |
| **region → lineage-unit** | 1 : 1..3 | 每個有 mutation 的 germline 家族 → 一 lineage-unit |
| **lineage-unit → tree** | 1 : 1..N | determined→1；ambiguous→N（枚舉全集）；capped→部分 |
| tree → node | 1 : many | 節點 = 觀測 genotype ∪ 隱藏祖先 H_* |
| node → sSNV | 1 : subset | 節點 = 一組已獲得的 sSNV（altset）|

**判斷軌跡（每 unit）**：`L0 家族 → L1 sSNV算法(determined/ambiguous/capped/recurrence + n_trees) → L2 CN(recurrence 才送) → L3 甲基(事後,不 rank)`。

---

## §5 最終資料格式（json schema）

**top-level**（`layered_reconstruction_{SAMPLE}.json`）:
```
L0_hp_family : {regions_total, regions_with_germ_family{1,2,3,none}, regions_mixed_germline, reads_by_family_total{...}}
L1_ssnv_algorithm : {n_lineage_units, determinacy_lineage{...}, determinacy_germline_1_2{...},
                     determinacy_somatic_3{...}, determinacy_unphased_none{...}, proportion_lineage{...},
                     n_verify_fail, all_V1V7_pass}
L2_cn : {n_recurrence_sent_to_cn, cn_split{artifact,LOH_unresolved,candidate}}
L3_methyl : {status}
detail : [ unit... ]
```
**每 unit（U4）**:
```
region, chrom, start, end, family(1/2/3/none), is_lineage, fam_label,
n_sSNV, cn, n_reads, n_full_pops, n_partial,
L1_class, L1_base_class, n_trees, n_hidden, capped, cap_reason,
trees:[{edges:[[parent,child]...], recurrence:[bit...], n_hidden}]  (≤6 範例),
n_trees_stored, L2_cn_verdict, L2_m_channel, verify{V1..V7}, verify_pass, max_vaf, trace:[...]
```

---

## §6 檢視表格範本（HTML/報告每次都呈現這 4 張，確保「完整狀況」）

1. **層級計數表**：U1–U6 各層計數 + 來源（§1）。
2. **比例字典表**：每比例 分子÷分母(單位) = 值（§2）→ 滑鼠 hover 顯示分母定義。
3. **determinacy × lineage 交叉表**：列=determined/ambiguous/capped/recurrence，欄=germline1/germline2/somatic3/none，格=count(比例)。
4. **per-region 逐 lineage 展開**：region → 各 lineage-unit 的 L1 class + 枚舉樹 carousel + L2 CN + 判斷軌跡。

> ⚠ **待 Bash 恢復**：U5/U6（Σtree、Σnode）+ P8 需從 detail 重算；並跑 `number_provenance` 產「每格→json:key」溯源表附於報告末（§13-C）。
