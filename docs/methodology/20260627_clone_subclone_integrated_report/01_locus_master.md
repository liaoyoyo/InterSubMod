<!--
建立時間: 2026-06-27
類型: L0 per-locus 主紀錄 — clone/subclone 整合報告（基礎層）
狀態: in_progress（HCC1395 ⭐3, Tier-R）
data_sources: data/sm_locus_master.json, data/sm_locus_master_summary.json
-->

# L0 — Per-locus 主紀錄（所有後續層的 join 基底）

> 目的：把每個 somatic SNV 的所有狀態 join 成**單一 SoT 表**，後續 L1-L4（HP/CCF/PS/甲基）皆以此為基底逐層加註。
> 圖：`figures/01_locus_master_overview.png`。資料：`data/sm_locus_master.{tsv,json}`（35,332 列，可逐位點 grep）。

## §1 欄位定義（format / 意義）

| 欄位 | 定義 |
|---|---|
| `locus` | `chrom:pos`（1-based）|
| `src` | TP / FP（ClairS-TO filter 標籤；FP = 未過 TP filter，仍可 clonally-informative）|
| `normal_ref/alt` | normal BAM pileup 在該位的 REF/ALT 讀數（MAPQ≥20）|
| `somatic` | normal_ref≥4 且 normal-VAF<0.05（容忍噪聲）= ONT somatic-like 候選（**非** SEQC2 金標準）|
| `tumor_ref/alt`, `vaf` | tumor 讀數 + ALT 分率 |
| `hp` | 該 sSNV ALT reads 的主要 longphase-S HP tag：1/2=germline、1-1/2-1=germline 背景帶 phased somatic、3=singleton phase block |
| `region` | 所屬「最大可關聯區域」（≤50kb chain），或 isolated |
| `tree_role` | ancestor / descendant / intermediate / sibling / region_member（co_linked 併或區內未連）/ isolated |
| `cn` | SEQC2 CN 狀態：gain / loh / neutral / loss |
| `seqc2_class` | highconf / superset_only / in_HC_not_called（γ 類）/ outside_HC |
| `n_partners`, `n_links` | ≤50kb 內 sSNV 鄰居數 / powered 同讀連結數 |

## §2 數據觀察（verified）

**規模 + sum-check**：35,332 loci（= union TP 30,490 + FP 4,842），**row 數 == union ✓**。

**tree_role（somatic sSNV，總 23,810）**：
| role | 數 | 意義 |
|---|---|---|
| region_member | 13,082 | co_linked 併入節點或區內未成 distinct 樹節點 |
| sibling | 3,994 | 互斥分支 |
| descendant | 3,156 | nested 後代 |
| ancestor | 3,014 | nested 祖先 |
| intermediate | 298 | 既祖先又後代（中間層）|
| isolated | 266 | 無 50kb 鄰居（單分子不可連）|

**HP 分布（linked somatic）**：**HP1-1 10,006 ≈ HP2-1 9,491**（somatic 在兩 germline 單倍型上**大致平衡**）+ HP3 1,317（longphase 無法 phase 到 germline 單倍型）+ HP1 1。
→ 觀察：somatic 突變並非偏一條單倍型；HP3（1,317）= 無 germline 錨、無法用 HP gate 判 subclone 的一群。

**SEQC2 truth**：highconf 30,492 / in_HC_not_called 4,834（γ 類）/ superset_only 6。

## §3 Canonical 驗證（chr17:48360161）
| sSNV | src | vaf | hp | tree_role | cn | seqc2 |
|---|---|---|---|---|---|---|
| γ 48357368 | FP | 0.176 | 1-1 | **sibling** | loh | in_HC_not_called |
| β1 48362515 | TP | 0.474 | 1-1 | region_member（併 β2）| loh | highconf |
| α 48365089 | TP | 0.823 | 1-1 | **ancestor** | loh | highconf |
| β2 48365161 | TP | 0.484 | 1-1 | **descendant** | loh | highconf |
→ 與先前 read-level worked-example 結構一致（α 祖先→β2 後代、γ sibling、β1 併 β2）；四者同 hp=1-1（同 germline 單倍型，**必要非充分於同克隆** — 見 `00b` §3）、同 loh。⚠ ancestor/descendant 標籤此處由 read-nesting 給定；VAF 方向驗證見 L2（且 loh 下 VAF 有 multiplicity 歧義）。

## §4 限制（含修對抗稽核）
- ⭐3 單樣本；`somatic` = ONT normal-VAF 過濾**非金標準**（偽影未清，見骨幹對抗稽核）。
- 🔴 **denominator 透明化**：35,332 loci 中 somatic=True **23,810（67%）**、somatic=**None 10,684（30%，normal 無 pileup 證據無法評估）**、somatic=False 838（2.4%）。本表觀察聚焦 23,810 somatic 子集；30% 缺 normal 證據的 sSNV 不納入 subclone 推論。
- 「HP 分布（linked somatic）」的 **「linked」操作定義 = somatic=True 且 n_links>0**（同讀 powered 連結 ≥1）；不加此濾則 HP1-1/2-1 = 11,403/10,926。
- tree_role 的 region_member 含 co_linked-併節點（細節在 `_assets/20260618_subcluster_pilot/sm_region_integration.json` — 注意該檔在 06-18 pilot 夾、非本報告 data/，為骨幹來源）。

## §5 此層輸出供後續
- L1 HP：用 `hp` 欄量化 subclone vs allelic 鑑別。
- L2 CCF：用 `vaf`+`cn` 換 CCF tier。
- L4 甲基：用 `region`+`tree_role` 定 genotype-anchored 群做 corroboration。
