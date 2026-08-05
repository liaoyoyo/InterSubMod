<!--
建立時間: 2026-07-13 00:00 +0800
目標: 在全 7 dataset 四／五類 compositional adjustment 實作前，鎖定 estimand、反例與 claim ceiling
處理範圍: chr1-22；7 dataset rows；6 biological IDs；historical layered-v2 engineering snapshot
關聯檔案:
  - InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_summary.tsv
  - InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_regions.tsv
  - InterSubMod/research/20260712_vaf_selected_shape_four_class_census/data/20260712_vaf_final_single_shape_four_class_by_source.tsv
-->

# Pre-decision audit：跨樣本結構 composition／bulk-sampling adjustment

> Verdict：**GO_DESCRIPTIVE_WITH_FAIL_CLOSED_CLAIMS（70/100）**。允許完整統計重算；禁止把跨 cell-line 差異當方法失效，也禁止把技術 pair 接近升格為 biological clone-tree truth。

## §0 Cynefin domain

- Domain：**Complicated（統計估計）＋ Complex（生物詮釋）**。
- 理由：Dirichlet posterior、compositional distance、rarefaction 與 cluster bootstrap 是可重現標準方法；但 source selection、unresolved closure 與不同 cell-line biology 會交互影響，不能套單一「越相同越好」門檻。
- 行動：完整跑預先列出的 sensitivity grid；不設定強迫六個 biological IDs 相同的 null。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---|---|---|
| summary 有 7 dataset rows＋aggregate，四類 resolved counts 與 complete/unresolved 守恆 | ✓ | L2 engineering | `...four_class_summary.tsv` |
| region-level resolved 表有 37,039 unique sample×region keys | ✓ | L2 engineering | `...single_shape_regions.tsv` |
| complete universe 可由 upstream coarse/VAF region 表重建為 39,885，補回 2,846 unresolved | ✓ | L2 engineering | `coarse_topology_all_regions.tsv` + `vaf_top_tie_regions.tsv` |
| 37,039 region strings 全可解析為 chr1-22:start-end；每 dataset 皆有 22 autosomes | ✓ | L2 engineering | 本 audit read-only schema probe |
| 5 Mb blocks：全 7 共同 490、union 569，可做 common-block bootstrap | ✓ | L2 engineering | 本 audit read-only schema probe |
| HCC1395 與 DORADO 是同 biological ID 的 technical pair；全表為 7 rows／6 IDs | ✓ | L2 engineering | region `biological_id` |
| fresh layered-v3 scientific release | □ | L5 pending | `docs/CURRENT_FOCUS.md`：尚未完成 7/7 downstream gate |

## §2 Credibility score

| Dimension | Score | Reason |
|---|---:|---|
| 理論基礎 | 20/20 | multinomial/Dirichlet、JSD/Aitchison、rarefaction、cluster bootstrap 皆為標準估計框架 |
| 觀察支撐 | 10/20 | 全 7 rows 齊，但來源是 historical layered-v2 engineering snapshot |
| 機制清晰度 | 10/20 | 已知 unresolved closure、sample size 與 selection-source mixture；仍無 biological truth |
| 反例風險 | 10/20 | O13 shared-read confound、source-selection collider、跨 cell-line 真異質性仍在 |
| 資源可行性 | 20/20 | 37k rows；純 Python；預期 <1 hr |
| **Total** | **70/100** | **GO for descriptive validation** |

Falsifier observable：若 HCC technical-pair 的接近只是 resolved-only closure、樣本量或 source-mixture cancellation，加入 unresolved、equal-n rarefaction、common-block bootstrap或 source standardization 後，其距離／跨生物 pair rank會明顯改變。此時只能報 adjustment-sensitive，不可報穩健再現。

## §3 Assumption map

| Importance | Known | Unknown |
|---|---|---|
| High | 7 rows/6 IDs；四類互斥；complete與resolved分母 | source status 是否是 confounder或 outcome-dependent selector；resolved closure 對距離的偏差量 |
| Low | 所有 row 的 n 皆足以估計比例 | 5 Mb 是否為最佳 correlation block；因此另做 10 Mb sensitivity |

MUST validate first：

1. 五類 complete（含 unresolved）為 primary，四類 resolved-only只作 secondary。
2. `structural_topo1`／`vaf_resolved_topogt1` 分層後，以 pooled global resolved source weights direct-standardize，並明示可能 over-adjust。
3. 技術 pair 的 percentile 只相對「biological_id 不同」的 source-level pairs；另輸出 6 biological-ID composition，不把 7 rows 當 7 biological replicates。

## §4 Quick pilot

1. 讀 summary／region／by-source schema → 驗證：類別、分母、座標與 source totals 守恆。
2. 只讀算 HCC raw TV → 驗證：resolved-only約 8.60%；complete五類約 14.55%。
3. 用 pooled global resolved source weights（structural=21,976/37,039）direct-standardize → 驗證：HCC resolved TV約 16.75%。

Checkpoint：三項已符合；進入全 7 rows 正式執行。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---|
| fresh layered-v3 7/7 aggregate | historical-v2不能作正式 biological Results | upstream cycle | P0 claim ceiling |
| single-cell／multi-region truth | 無法驗證 clone tree truth | 外部資料 | P0 scientific ceiling |
| read-overlap／region-construction完整抽樣模型 | bootstrap仍是近似 | 中 | P1 future |
| 最佳 spatial block length | CI 可能隨 block size | 低 | P1，做5/10 Mb sensitivity |

## §6 Evidence conflict scan

| Prior conclusion | Relation | Source |
|---|---|---|
| O13 cross-region methylation correlation NEGATIVE；shared-read-count confound | dependent warning；本輪不重開甲基跨區 correlation，也不稱 clone | `docs/reports/validated/2026/04/20260408_TO_LOH內外ISM特徵區分力完整推論鏈報告_01.md:426` |
| Cross-region subclone consistency 已標 NEGATIVE | conflict if overclaimed；本輪只估 pattern composition | `docs/experiments/INDEX.md` 附錄 |
| region tree ≠ biological clone；historical v2 publication NO-GO | binding ceiling | `docs/CURRENT_FOCUS.md`、20260712 四類報告 |
| repo `MEMORY.md` | unavailable | 啟動時檔案不存在；以 CURRENT_FOCUS／INDEX／validated report 代替 |

## §7 Decision path

- Verdict：**GO_DESCRIPTIVE_WITH_FAIL_CLOSED_CLAIMS**。
- Decision lock：本輪不可用 `same distribution`、`same evolution`、`validated clone tree`；只能用 technical-pair compositional proximity／sensitivity。
- Primary estimand：五類 complete composition（Single、Sister、Direct、Sister+direct、Unresolved）。
- Secondary estimand：resolved-only四類 composition。
- Required sensitivity：Jeffreys Dirichlet intervals、JSD＋Aitchison、equal-n rarefaction、5/10 Mb common-block bootstrap、biological-ID-aware EB shrinkage、source direct standardization。
- Red team failure modes：closure bias、source-mixture cancellation、spatial pseudoreplication、technical pair double-counting、不同 cell-line biology被錯當應相同。

