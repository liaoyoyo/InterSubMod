<!--
建立時間: 2026-06-26
類型: 數據紀錄 — HCC1395 全基因組 sSNV 連鎖區域統計畫像（明確驗證紀錄）
狀態: in_progress（單樣本 ⭐3, Tier-R）
build_branch: feat/summary-nreadsvalid
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/sm_region_stats.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_region_integration.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_configuration_census.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_completeness_ledger.json
-->

# HCC1395 全基因組 sSNV 連鎖區域 — 統計畫像（明確紀錄）

> 單位 = **最大可關聯區域**（≤50kb somatic-sSNV chain）。HCC1395 ⭐3 單樣本，Tier-R。
> 全部數字落檔可重算（sum-check 已過）；每個比例分母明確標示。

## §1 區域總覽（region inventory）

| 指標 | 值 |
|---|---|
| **含 ≥2 somatic sSNV 的區域數** | **7,143** |
| 這些區域內的 somatic sSNV 總數 | 23,544 |
| 每區域 sSNV 數 中位數 / 最大 | 2 / 150 |
| 完整性帳本（全 35,332 union sSNV）| linked 21,554 / underpowered 5,458 / isolated_singleton 8,320（sum=35,332 ✓）|

> 注：「區域」= 至少 2 個 somatic sSNV 在 ≤50kb 內相連。單一 somatic（無 50kb 內鄰居）= isolated，不成區域。

## §2 每區域 sSNV 數量分布（n=7,143）

| sSNV/區域 | 區域數 | 比例 |
|---|---|---|
| 2 | 3,707 | **51.9%** |
| 3 | 1,581 | 22.1% |
| 4 | 789 | 11.0% |
| 5 | 426 | 6.0% |
| 6–7 | 370 | 5.2% |
| 8–10 | 180 | 2.5% |
| 11–20 | 66 | 0.9% |
| 21–50 | 15 | 0.2% |
| 51+ | 9 | 0.1% |

🔑 **74% 區域只有 2–3 個 sSNV**；多位點密集區（≥8 sSNV，共 270 區域 / 3.8%）稀少，且部分為 CN-gain/segdup 偽影簇（見 §6）。

## §3 區域跨度（span）分布

| span | 區域數 | 比例 |
|---|---|---|
| <1kb | 319 | 4.5% |
| 1–5kb | 529 | 7.4% |
| 5–10kb | 609 | 8.5% |
| **10–50kb** | **4,005** | **56.1%** |
| 50–200kb | 1,630 | 22.8% |
| ≥200kb | 51 | 0.7% |

## §4 連接關係 census（每對 2×2 cell-pattern；powered+somatic）

cell: RR/RA/AR/AA（a/b 的 REF/ALT）。同HP=克隆相關，異HP=allelic。

| 連接關係 | cell 模式 | 同HP 對 | 異HP 對 | 意義 |
|---|---|---|---|---|
| 互斥 mutual_excl | RA+AR, AA=0 | 2,825 | 4,080 | sibling(同) / allelic(異) |
| nested a⊂b | RA+AA, AR=0 | 5,563 | 1,211 | b 祖先→a 後代 |
| nested b⊂a | AR+AA, RA=0 | 5,763 | 576 | a 祖先→b 後代 |
| 共連 co_linked | AA only | 10,254 | 1,496 | 同 lineage |
| independent | 全有 | 5,547 | 734 | 無乾淨結構 |

- **nested 同HP 11,326 >> 異HP 1,787**（巢狀=克隆階層）；**互斥異HP 4,080 > 同HP 2,825**（互斥多為 allelic 非 subclone）。
- per-sSNV 去重（非 per-pair 灌水）：有確認同HP連結的 distinct sSNV = **14,743**（nested 7,947 / sibling 5,632 / co_linked 5,022）。

## §5 區域結構（樹形）分布 + 深度（n=7,143）

| 樹形 | 區域數 | 比例 | 乾淨(LOH+neutral) |
|---|---|---|---|
| **full_tree**（分支+深度）| 677 | 9.5% | **205** |
| linear_nested（祖先→後代鏈）| 1,908 | 26.7% | 763 |
| sibling_only（平行分支）| 1,235 | 17.3% | 408 |
| co_linked_lineage（單 lineage）| 858 | 12.0% | 357 |
| no_confirmed_structure | 2,443 | 34.2% | 702 |
| inconsistent（cycle）| 22 | 0.3% | 5 |

**有確認克隆結構**（前 4 類）= **4,678 區域（65.5%）**；其中 **分支樹（full_tree）677（9.5%）**。
**樹深度**：depth 1 = 4,545、depth 2 = 2,402、depth 3 = 169、depth 4 = 17、depth 5 = 4、depth 9 = 1。→ 絕大多數為淺樹（深度 1–2）；深度 ≥3（191 區域）多在 CN-gain 密集區（偽影嫌疑）。

## §6 CN 分層（誠實關鍵）

| CN 狀態 | 區域數 | 比例 |
|---|---|---|
| **gain（混淆）** | **4,580** | **64.1%** |
| loh（乾淨）| 2,188 | 30.6% |
| neutral（乾淨）| 252 | 3.5% |
| loss | 123 | 1.7% |

🔴 **64% 區域在 CN-gain**（multiplicity/amplicon/segdup 混淆）→ 論文級可信集 = **LOH+neutral = 2,440 區域（34.1%）**，其中 full_tree 205。

## §7 驗證（明確）
- **sum-check**：sSNV-per-region 直方圖加總 = 7,143；樹形分布加總 = 7,143；帳本三桶加總 = 35,332（= union）。全部 ✓。
- **逐區域可查**：`_assets/20260618_subcluster_pilot/lists/regions.tsv`（每區域一列）+ `sm_region_integration.json`（含 nested_edges / sibling_pairs / populations）。
- **可重現**：`README_sm_linkage_pipeline.md`（pipeline DAG + 重現指令 + 數字→檔案對照）。

## §7b 分支型態 + 相鄰區域一致性（data_source: sm_branch_analysis.json）

**分支型態 taxonomy（4,678 有結構區域）**：
| 型態 | 區域數 | 比例 |
|---|---|---|
| **linear_chain（線性祖先→後代鏈，無分支）** | **1,908** | **40.8%** |
| sibling_only（單一 2-way 分支）| 1,003 | 21.4% |
| co_linked_single_lineage（單 lineage）| 858 | 18.3% |
| full_tree（分支+深度）| 677 | 14.5% |
| sibling_multi（≥3 分支）| 232 | 5.0% |

- **最常見結構 = 線性鏈（40.8%）**，多數區域是 ancestor→descendant 無分支；**有分支（sibling≥1）= 1,923 區域（26.9%）**，其中 1,235 純平行分支、688 分支+深度。
- **分支型態 = 簡單二分支主導**：sibling 數 1=1,310 / 2=344 / 3=126 / 4=56 / 5=28 → 絕大多數只 1 個 sibling 對（2-way split），深/多分支稀少。

**相鄰區域 shape 一致性（7,121 相鄰對，中位 gap 213kb）**：
- observed same-shape = **0.279** vs shuffle-null **0.252±0.005**，**z=5.7** → **統計顯著但效應小（僅 +2.7pp）** = 弱空間自相關。相鄰區域的樹形主要由「該區域剛好有哪些 sSNV」決定（局部），非強烈受全基因組克隆結構壓印。

**CN-clean VAF 峰（subclone CCF 一致性，6,138 clean VAF 值）**：
- 主峰 VAF≈**0.95(987) / 0.90(721) / 0.45(504) / 0.40(464) / 0.05(624)** → **有少數一致 CCF 層級**（高≈truncal、中≈0.4-0.45 主 subclone、低≈0.05 罕見/偽影）= 跨區域 subclone CCF 有一致性（即使各區樹形不同，因每區只取樣到部分 subclone）。

🔑 **結論**：HCC1395 的局部克隆結構以**線性鏈**為主、分支多為**簡單二分支**；相鄰區域樹形**弱一致**（z=5.7 但小），但 CN-clean 的 **VAF/CCF 層級有一致峰** → 真克隆的固定 subclone 集存在、但被各區域局部取樣稀釋成不同樹形。

## §8 🔴 限制
⭐3 單樣本；Tier-R only；64% 在 CN-gain；偽影未清（缺 mappability）；Fisher-sig 不分 subclone/allelic（HP-gate 才分）；regional（≤read-span）非 genome-wide tree；分子證據非 single-cell confirmation。
