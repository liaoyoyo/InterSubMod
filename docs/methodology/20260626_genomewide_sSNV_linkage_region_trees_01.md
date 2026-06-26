<!--
建立時間: 2026-06-26
類型: methodology — 全基因組 sSNV 單分子連鎖 → 每區域克隆樹重建（HCC1395 ⭐3）
狀態: in_progress（單樣本 ⭐3, Tier-R, 已過對抗稽核並套 5 修正）
build_branch: feat/summary-nreadsvalid
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/sm_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_region_integration.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_configuration_census.json, docs/methodology/_assets/20260618_subcluster_pilot/sm_completeness_ledger.json
-->

# 全基因組 sSNV 單分子連鎖 → 每區域克隆樹重建 — 完整方法與結果

> **敘述框架**：Verdict-Pyramid（結論金字塔）。HCC1395 單樣本 ⭐3；Tier-R（same-read ≤50kb，same-PS deferred）。
> 數字全來自 data_sources JSON（§13-A），每個可 grep 重算（見 `_assets/.../README_sm_linkage_pipeline.md`）。
> **已過 fresh-context 對抗稽核**（`_assets/.../sm_adversarial_audit.md`，NEEDS_WORK → 套 5 修正）。

## §0 一句話結論

從 chr17 單一 worked example（n=1）擴展到**全基因組 exhaustive**：union TP∪FP 35,332 sSNV → 單分子 2×2 共現 → 每「最大可關聯區域」重建局部克隆樹。**7,143 個區域中 677（9.5%）有完整分支樹（sibling+nested）、205 在乾淨 CN 區**。方法可行、chr17 被自動重建；但**~69% 訊號在 CN-gain 混淆區、Fisher-sig 不分 subclone/allelic** — 乾淨可信集 = LOH+neutral，且為分子證據非 single-cell confirmation（⭐3）。

## §1 方法三層遞進（2-locus → multi-locus → region）

| 層 | 單位 | 做什麼 | 輸出 |
|---|---|---|---|
| **L1 pairwise** | 2 sSNV | per-read 2×2 共現（RR/RA/AR/AA）→ classify | configuration census |
| **L2 multi-locus** | ≥3 sSNV | 覆蓋全部 sSNV 的 read 基因型向量 → 不同向量 = 局部 population | populations/region |
| **L3 region** | 最大可關聯區域 | 整合 nested(祖先→後代)+sibling(分支)+co_linked(併lineage) → 樹 + 樹形分類 | per-region tree |

**高效**：group-based single-pass（每 sSNV 只 pileup 一次），5 平行 worker，全基因組 ~40min（linkage）+ ~15min（multi-locus）。[L1]

## §2 完整 sSNV 宇宙 + 完整性帳本（no-miss, Tier-R）[L1]

- **宇宙 = 35,332 sSNV（TP 30,490 + FP 4,842）**；union 直接修補 chr17 γ 類漏失（γ 在 FP 不在 TP）。
- **完整性帳本**（每 sSNV 歸恰好一桶，sum-check ✓ = 35,332）：**linked 21,554（61%）/ underpowered 5,458（15%）/ isolated_singleton 8,320（24%）**。
- 🔴 「不漏」限定 **Tier-R**（same-read ≤50kb）；same-PS（phase-set >50kb）deferred → isolated 是上界非真值。

## §3 配置 census（2×2 cell-pattern 實際數量）[L1]

powered（coread≥6）+ 兩端 somatic，strict cell-occupancy（cell: RR/RA/AR/AA）：

| 配置 | cell 模式 | 同HP 對 | 異HP 對 | distinct sSNV 同/異 |
|---|---|---|---|---|
| 互斥 mutual_excl | RA+AR, AA=0 | 2,825 | 4,080 | 4,174 / 5,982 |
| nested a⊂b（b祖先→a後代）| RA+AA, AR=0 | 5,563 | 1,211 | 4,915 / 510 |
| nested b⊂a（a祖先→b後代）| AR+AA, RA=0 | 5,763 | 576 | 4,979 / 479 |
| 共連 co_linked | AA only | 10,254 | 1,496 | 5,256 / 700 |
| independent | 全有 | 5,547 | 734 | 3,905 / 1,101 |

- **nested 同HP 11,326 >> 異HP 1,787**（巢狀=克隆階層，幾乎同單倍型 ✓）；**互斥異HP 4,080 > 同HP 2,825**（互斥多為 allelic 非 subclone，已隔開）。
- per-sSNV distinct（修對抗稽核 F1 per-pair 灌水）：any confirmed same-HP relation = **14,743 distinct sSNV**（nested 7,947 / sibling 5,632 / co_linked 5,022）。

## §4 每區域克隆樹分布（最終整合，單位=最大可關聯區域 n=7,143）[L1]

| 樹形 | 數量 | 比例 | 乾淨(LOH+neutral) |
|---|---|---|---|
| **full_tree**（分支+深度，chr17式）| 677 | 9.5% | **205** |
| linear_nested（祖先→後代鏈）| 1,908 | 26.7% | 763 |
| sibling_only（平行分支）| 1,235 | 17.3% | 408 |
| co_linked_lineage（單 lineage）| 858 | 12.0% | 357 |
| no_confirmed_structure | 2,443 | 34.2% | 702 |
| inconsistent（cycle）| 22 | 0.3% | 5 |

→ **3,820 區域（53%）有確認的克隆分支/階層結構**；其中 **full_tree 677（9.5%）最豐富**。乾淨 CN 區（LOH+neutral）的 full_tree = **205** = 論文級可信子集。

## §5 chr17:48360161 canonical（被 pipeline 自動重建）[L1]

shape=**full_tree**、CN=**loh**、4 sSNV、3 populations：`ARRR`（γ subclone）7 reads / `RRAR`（L1 α-only）10 / `RAAA`（L2 α+β）15。樹：germline →（γ sibling ∥ α 祖先 → β 後代），VAF α 0.82 > β 0.48 > γ 0.18，0 violations。詳見 `20260625_chr17_subclone_worked_example_explained_01.md`。

## §6 驗證 [L1/L2]

- **null**：Fisher exact（= 條件 2×2 檢定）；**sparse 類 2.2% 顯著 ≈ 機率底線**（power-gate 不從噪聲造結構）vs co_linked 92% / nested 71–81%。
- **Fisher vs free-margin Monte-Carlo**：300 取樣中 agree 258/264（97.7%），median |Δp|=0.009 → Fisher 結論穩健。
- 🔴 **Fisher-sig 不分 subclone/allelic**（diffHP negative control ~50% 也顯著）→ 鑑別**全靠 HP-gate**。
- **SEQC2**：宇宙 86% 為 HighConf-concordant。
- **chr17 + 完整 sum-check** 重現。**對抗稽核** fresh-context evaluator NEEDS_WORK → 5 修正已套。

## §7 🔴 限制（誠實，必守）

1. **CN-gain 混淆主導**：~69% linked-somatic 在 CN-gain（multiplicity/amplicon/segdup）→ 乾淨集 = LOH+neutral。[L1]
2. **偽影未清**：dense-uniform-cluster（chr8:81–83M 97+82+65 sSNV、chr9:41.8M、chr14:16.09M ~3% 均勻 VAF）= mapping/segdup 嫌疑，缺 mappability track 未移除。[L1]
3. **γ 類非確證**：3,204 FP-source somatic-linked = 候選需正交驗證，**非「SEQC2 漏掉的確證 somatic」**；「3204≈3200 收斂」為套套邏輯。[L2]
4. **Tier-R only**：same-PS 連鎖 deferred。
5. **⭐3 單樣本**：regional（≤read-span）**非 genome-wide clone tree**；分子證據**非 single-cell confirmation**；對外勿稱「甲基偵測 subclone」「genome-wide tree」「對手缺檢定」。

## §8 檔案與重現（可驗證 + 可查詢）

- **manifest**：`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/README_sm_linkage_pipeline.md`（pipeline DAG + 重現指令 + 數字驗證指南）。
- **headline SoT**：`sm_summary.json`（小檔可獨立 recompute）+ `sm_region_integration.json`。
- **查詢 TSV**：`lists/per_sSNV_census.tsv`（35,332 sSNV）/ `lists/{config}_{HP}.tsv`（逐對）/ `lists/regions.tsv`（逐區域）。
- **對抗稽核**：`sm_adversarial_audit.md`。大 JSON 不入版控（manifest 可重生）；腳本+小summary+TSV+文件入版控。
