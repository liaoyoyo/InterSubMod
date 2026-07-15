<!--
建立時間: 2026-07-15 13:00
目標: current canonical v5 layered workstation 加入全基因組拓撲／形態圖層與多選互動前審計
處理範圍: chr1-22、7 datasets、current-v5 read-AF 描述性排序、HTML 互動與視覺 QA
cycle_id: cycle_20260715-1300-layered-workstation-genome-topology
topic: layered_workstation_genome_topology_multiselect
status: verdict_GO
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md
  - InterSubMod/research/20260712_vaf_selected_shape_four_class_census/20260712_VAF最終單一Shape四類全樣本重算_01.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation_v5.py
-->

# Pre-Decision Audit：全基因組拓撲圖層與多選互動

> **Verdict**：描述性 current-v5 重算與介面改版 **GO**；「read-AF 確認真樹／confirmed clone」維持 **NO-GO**。
> **任務類型**：B Comprehensive validation；全 chr1–22 × 7 datasets，不可用 historical subset 代替。

## §0 🟤 Cynefin Domain Gate

- **Domain**：Complicated。
- **Test**：相同 exact candidate enumeration、family-specific `ALT/(REF+ALT)` 與 AHU shape canonicalization 已可重複得到決定性結果；資料版本接合與語意呈現仍需專家稽核。
- **處置**：先以 HCC1395 current-v5 做 probe，再擴全 7 datasets；介面採既有產生器原生 SVG/JS。

## §1 🟢 Observation Completeness Checklist

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| current-v5 7/7 canonical run 與 chr1–22 scope 已通過 | ✓ | L1 ⭐⭐⭐⭐⭐ | `research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` |
| HCC1395 W_primary=7,932，current structural 四類可完全回對 | ✓ | L1 ⭐⭐⭐⭐⭐ | 2026-07-15 current probe；本 cycle 後續保存 durable artifact |
| HCC1395 5,029 ambiguous complete regions 可評估；4,574 unique-first、366 tie-same-shape、89 tie-multi-shape | ✓ | L2 ⭐⭐⭐⭐ | 2026-07-15 current probe；exact Fraction tie |
| current HTML 只保存每 unit 前 32 棵 display trees，會漏 candidate_index>32 的 top tree | ✗ | L1 ⭐⭐⭐⭐⭐ | `layered_tree_reconstruction.py:78-79,179-205`；HCC probe 241 units 受影響 |
| current 7 datasets 的 read-AF／形態完整重算 | ✓ | L1 ⭐⭐⭐⭐⭐ | `data/current_v5_read_af_topology/current_v5_read_af_topology.index.json`；7/7 checks pass |
| 桌面與手機多選互動／視覺稽核 | ✓ | L1 ⭐⭐⭐⭐⭐ | `qa/full/validation_receipt.json`；index 2 + sample 14 = 16/16 runs、219/219 checks |

## §2 🔵 Credibility Score

| 維度 | 評分 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | exact Fraction ordering、rooted-unordered shape 與 set union 互動皆定義明確 |
| 觀察支撐 | 10 | current HCC1395 probe 通過；全 7 datasets 待本 cycle 重算 |
| 機制清晰度 | 20 | source counts → exhaustive candidates → exact tie → region class → ideogram mark |
| 反例風險 | 10 | 已知 CN/purity/coverage 與 biological overclaim 風險；以 claim ceiling 與 N/A fail-closed 管理 |
| 所需資源 | 10 | 7-sample solver 約數分鐘，完整 HTML/Playwright 約 1–6h |
| **TOTAL** | **70 / 100** | **GO** |

**Falsifier observable**：若接合或分類錯誤，會看到任一樣本 structural bins 不回到 current `W_primary`、候選／shape mismatch >0、形態 bins 不守恆、0/0 被標成 0%、或多選 A+B 後顯示不是聯集。

**Reality-test 三個反例觀察**：

1. read-AF exact first 可能落在原 display index >32；新 artifact 必須帶 top tree edges。
2. primary-family 任一 locus 無有效 REF+ALT 分母時必須 N/A，不可當 0% 排序。
3. 多個 exact-top shapes 或 incomplete region 不得被標為唯一 clone/subclone。

## §3 🟡 Assumption Map

| Assumption | Importance | Known | Quadrant／處置 |
|---|---|---|---|
| current run 與 HTML region key 可 exact join | HIGH | unknown | (2) 先驗：7/7 join conservation |
| read-AF score 是描述性順位、非 calibrated probability | HIGH | known | (1) UI claim ceiling |
| Topo=1 代表 mutation-state shape 唯一、非 biological truth | HIGH | known | (1) UI claim ceiling |
| 形態可由 top shape 的 max depth/outdegree 分四類 | HIGH | known | (1) historical deterministic method + current 重算 |
| category Set 空集合語意為全部 | HIGH | known | (1) protected user requirement + Playwright |
| gene/cytoband/hotspot 可由座標自行推論 | LOW | false | (4) 不顯示未提供註解 |

## §4 🟠 Quick Pilot（已執行）

1. 讀 current manifest、machine summary 與 historical exact-ranking code。
2. 以 HCC1395 current-v5 exhaustive `tree_cap=0` 重枚舉；實測 33.7 秒、RSS 約 290 MB。
3. Check：W_primary/structural bins 完全回對、candidate/shape mismatch=0、Fraction-vs-float error ≤1e-12。

**Checkpoint**：全部通過，允許擴成 7 datasets；若任何樣本 mismatch >0，該樣本 fail closed，不產生「第一順位」分類。

## §5 🔴 Gap Diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| current-only durable ranking artifact（含 top edges） | HIGH | 1–2h | P0 |
| 7-sample morphology conservation checks | HIGH | <1h | P0 |
| Chromium desktop/mobile interaction receipts | HIGH | 1–2h | P0 |
| CN/purity-calibrated clone probability | HIGH for biological truth；不影響描述性 UI | >6h | P2 / out of scope |

## §6 🟣 Evidence Conflict Scan

| Prior conclusion | Tier | Relation | Source |
|---|---|---|---|
| read-AF ordering 是 heuristic，不是 CCF／獨立真樹確認 | L2 | constraint；禁止「confirmed／機率」措辭 | `research/20260711_read_group_C_tree_T_topology_report/pre-decision-audit.md:20-33` |
| confirmed clone/subclone／true ancestry 維持 NO-GO | L2 | constraint；只顯示相容形態 | `research/20260712_vaf_selected_shape_four_class_census/20260712_VAF最終單一Shape四類全樣本重算_01.md:15-24` |
| historical layered-v2 不代表 current canonical v5 | L1 | conflict；禁止沿用舊 TSV | current probe + current summary schema |
| repo 根目錄無 `MEMORY.md` | — | governance gap；已改查 topic audit 與 validated reports | 2026-07-15 `rg --files -g MEMORY.md` 無結果 |
| LOH methylation hypothesis negative | validated | unrelated | `docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md` |

**Conflict count**：3 個相關限制；已將反例風險維持 10，不升 20。

## Red-team Gate

- **Failure mode 1**：把舊版 VAF ranking join 到 current regions，導致錯位。防線：current input-only 重算 + 7/7 exact key conservation。
- **Failure mode 2**：只重排 stored top-32，漏掉真正第一順位。防線：artifact 輸出 exhaustive top edges；HCC 已觀察 241 units 會受影響。
- **Failure mode 3**：使用「確認 clone／最可能機率」造成過度推論。防線：標 `read-AF 唯一第一順位（描述性）`、`clone/subclone-compatible`，並保留 unresolved。
- **Falsifier 可否證**：可；各 conservation/mismatch/interaction 指標皆能直接 fail build 或 test。

## §7 ⚫ Decision Path

- **TOTAL**：70 / 100
- **Verdict**：**GO for current-v5 descriptive reconstruction + UI**。
- **Claim ceiling**：regional mutation-state candidate topology；不是 confirmed clone、true ancestry、CCF 或 calibrated likelihood。
- **Next action**：建立 implementation notes → 7-sample durable artifact → HTML rebuild → Playwright QA。
- **Decision lock**：Y；只有新 current data、新校準方法或 upstream correction 才可提升 biological claim。

## Provenance Footer

- **Commit hash**：`b7826318a1813878a9156bf63401e2b45fbf74eb`
- **Build time**：2026-07-15 13:00:30 +0800
- **Skill**：`/pre-decision-audit` v0.1
- **Cycle**：`InterSubMod/state/cycles/cycle_20260715-1300-layered-workstation-genome-topology/`
- **Next**：`InterSubMod/research/20260715_layered_workstation_genome_topology_multiselect/implementation-notes.md`
