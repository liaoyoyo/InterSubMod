<!--
建立時間: 2026-07-25 07:00
目標: 在 ALT/REF 甲基證據整合進 exact-PS layered workstation 前，先固定 grain、authority 與 claim ceiling
處理範圍: 7 個 formal G1 pairs、7 個 focal ALT/REF controls、10 個 exact-PS×HP lanes、目前 all7_v2 candidate factorization
關聯檔案:
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/results/validation_receipt.json
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/formal_pair_alt_ref_topology_join.tsv
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/focal_alt_ref_control_join.tsv
  - InterSubMod/research/20260725_methyl_alt_ref_topology_overlay_validation/data/strict_hp_lane_cpp_topology_join.tsv
  - InterSubMod/research/20260724_exact_ps_cpp_topology_signature_census/scripts/build_exact_ps_layered_workstation.py
-->

# ALT/REF 甲基證據整合工作站：Pre-decision audit

> Task type B — 對預先定義的 7 個 formal G1 positives 做完整資料與介面驗證；不是 98,955 groups 的全面甲基陰性／陽性普查。服務 G3、G4。

## §0 Cynefin domain

**Complicated。** 現有工作站已有 receipt-bound annotation、region detail 與候選樹 overlay 契約；本次可沿用相同 fail-closed pattern，但 pair、site-control、HP lane 三種 grain 必須分開。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| formal G1 positive pair = 7 | ✓ | L1 | `data/formal_pair_alt_ref_topology_join.tsv` |
| focal ALT/REF joint testable = 3、axis-aligned = 1 | ✓ | L1 | `data/focal_alt_ref_control_join.tsv` |
| exact-PS×HP lanes = 10；signal support = 638、RR-only background support = 152、total = 790 | ✓ | L1 | `data/strict_hp_lane_cpp_topology_join.tsv` |
| resolved focal→partner relation = 3；三者皆跨所有 best trees 一致 | ✓ | L1 | validation receipt + candidate factorization |
| 未進 7-pair overlay 的工作站 group 是否為甲基陰性 | ✗ | L1 | 本研究只保存 formal-positive targeted subset |
| HCC1395_DORADO methyl association replication | ✗ | L1 | 只有 genetic candidate order 重現；methyl multigroup 未重現 |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20/20 | 三層 grain 與 read population 已有明確定義 |
| 觀察支撐 | 10/20 | formal subset 涵蓋三個 dataset key，但屬陽性篩選集合 |
| 機制清晰度 | 20/20 | pair → site control → PS×HP lane → model candidate 可分層追溯 |
| 反例風險 | 10/20 | selection bias、pair-to-lane expansion 與 resource ABSTAIN 仍可能被誤讀 |
| 所需資源 | 10/20 | 全量 7-page rebuild 與 Chromium audit 預估超過 1 小時 |
| **TOTAL** | **70/100** | **GO，附 fail-closed gates** |

**Falsifier observable：** 任一 pair 找不到 focal control、任一 lane 找不到目前 exact `sample×region_id`、join 出現重複膨脹、H2009 chr4 不再是唯一 `ALIGNED`、或目前 all7_v2 重算出的 pair relation/count 與 TSV 不一致，即停止 HTML 整合。

## §3 Assumption map

| 重要度 | 已知 | 未知 |
|---|---|---|
| 高 | 7 pair、10 lane、1 aligned、3 resolved 的預期守恆 | all7_v2 是否逐列重現舊 all7_v1 join；必先驗 |
| 低 | 現有 UI 可新增 mode、filter、detail panel | 甲基 overlay 是否值得未來擴成全 cohort 普查；本輪不回答 |

## §4 Quick pilot

1. 以 `sample×cpp_region_id` 對目前 topology/candidate sidecar。
2. 只先生成 H2009 頁，驗證 chr4:2,307,521 顯示 `ALIGNED` 但 candidate `ABSTAIN_RESOURCE_LIMIT`。
3. 驗證 H2009 chr18 顯示 full tree tied、局部 focal→partner = 6/6。

Checkpoint：上述三項皆通過才重建全部 7 頁；任一失敗改為 PROBE。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| validation report 尚綁 all7_v1 | 可能與工作站 all7_v2 authority 漂移 | <1h | P0 |
| 工作站沒有 pair/site/lane 分層 payload | 容易把 7/7 誤讀成 focal ALT/REF 7/7 | 1–3h | P0 |
| Chromium 尚無 methyl-specific regression | 可能有篩選或 mobile overflow 回歸 | 1–3h | P1 |

## §6 Evidence conflict scan

- `MEMORY.md`：repository root 無此檔，不能完成該入口的 conflict scan。
- 找到 `docs/reports/validated/2026/04/20260401_LOH_weekly_review/06_methylation_hypothesis_negative.md`；它限制的是既有 LOH／甲基假說，不否定本次 read-level annotation，但支持「不可升級為因果或 clone identity」。
- 目前 topic 報告明示 HCC1395_DORADO methyl multigroup 未重現；這是必須保留的反例。

## §7 Decision

**Verdict: GO（70/100）。**

Decision lock：

1. 甲基層只作 evidence overlay，不改變 topology candidate、read-AF 排序或 selected tree。
2. `7/7 formal G1` 永遠與 `focal ALT/REF aligned 1/7` 分開顯示。
3. `NOT_IN_FORMAL_OVERLAY` 代表未進本次 targeted subset，不是甲基陰性。
4. 只有在 focal/partner 都是 active loci，且目前 all7_v2 的所有 global-best trees 關係一致時，才把 focal→partner 標到 candidate view。

### Red-team gate

- Failure mode 1：7 個 formal positives 是從 407,738 evaluated rows 篩出的陽性集合；把 7/7 當準確率會造成 selection bias。
- Failure mode 2：一個 pair 可展開為多個 HP lanes；直接以 pair row 塗滿 lanes 會混淆 site-level control 與 HP-specific topology。
- Failure mode 3：來源報告原先綁 all7_v1；未更新到 all7_v2 就嵌入會形成跨 authority drift。

三項皆可由 current-v2 identity、grain-specific payload 與 UI wording 驗證，因此維持 GO。
