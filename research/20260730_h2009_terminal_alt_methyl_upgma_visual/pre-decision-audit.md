<!--
建立時間: 2026-07-30
目標: 在建置 H2009 chr18:563687 ALT-read 甲基距離熱圖與 UPGMA HTML 前，鎖定資料範圍、可否證條件與 claim ceiling
處理範圍: 單一 H2009 位點、86 個 after-peel focal-ALT reads；Task F 示意，不作全樣本 prevalence 驗證
關聯檔案: InterSubMod/research/20260730_h2009_terminal_alt_methyl_upgma_visual/scripts/build_h2009_terminal_alt_methyl_visual.py
-->

# Pre-decision audit：H2009 terminal-ALT methyl UPGMA visual

> **DEMO / PARTIAL SCOPE — Task F；服務 G3。** 本工作只把既有 M1 與 exact-PS 拓撲證據轉成可核對的 HTML，不新增 clone/subclone 判定規則，也不寫入 validation ledger。

## §0 Cynefin domain

**Complicated / best-practice。** 相同資料矩陣、average-linkage UPGMA 與 portable artifact 流程已有可重複成功案例；本次風險主要是資料範圍與語意標示，不是未知演算法行為。

## §1 Observation completeness

| Observation | 狀態 | Tier | Source |
|---|---:|---:|---|
| H2009 chr18:563687 有 87 ALT、peel 1、保留 86 | ✓ | L1 | `all_ssnv_site_results.tsv.gz` exact row |
| 86 reads 的 M1 labels 為 81／5 | ✓ | L1 | `all_ssnv_stable_multigroup_read_assignments.jsonl.gz` exact row |
| root split between/within=2.8944385871，高於 null95=1.3814517633 | ✓ | L1 | assignment `coarse_split_trace[0]` |
| 原始 region tree 是 124 ALT+REF reads，不能直接當 M1 UPGMA | ✓ | L1 | region metadata、matrix shape 與 assignment read-ID join |
| 6/6 VAF-global-best mutation trees 共用 terminal edge 14→15（新增 563687） | ✓ | L1 | exact-PS candidate factorization |

## §2 Credibility score

| Dimension | Score | Reason |
|---|---:|---|
| 理論基礎 | 20/20 | read-distance heatmap 與 average-linkage UPGMA 是既有 M1 距離的直接視覺化 |
| 觀察支撐 | 20/20 | 單一指定案例的 primary matrices、read IDs、labels 與 split trace 完整 |
| 機制清晰度 | 20/20 | assignment IDs → 86×86 subset → symmetrize → UPGMA → same-order heatmap |
| 反例風險 | 10/20 | 單位點容易被誤解為 clone truth；以 claim ceiling 與 124/86 guard 降低風險 |
| 所需資源 | 20/20 | 小矩陣、單次 HTML packaging，預期低於 1 小時 |
| **TOTAL** | **90/100** | **GO** |

**Falsifier observable：** 任一情況成立即停止建置：86 IDs 無法全數 join、矩陣不對稱/對角非零、重算 ratio 不等於 formal trace、UPGMA root split 不形成 81/5 連續兩塊、VAF terminal edge 不是 6/6、或 packaged HTML 發生外部網路請求。

## §3 Assumption map

| Importance | Known | Unknown |
|---|---|---|
| High | 86 after-peel ALT read IDs；formal labels；primary matrix SHA；average linkage | HTML renderer 對超寬圖在窄螢幕的實際縮放 |
| Low | 色彩與 aliases 只影響閱讀，不影響統計 | 個別 branch rotation 的視覺方向 |

高重要性未知項以 desktop/mobile browser QA 驗證；branch rotation 明示不改變距離或樹結構。

## §4 Quick pilot

1. 重建 86×86 subset → checkpoint：shape=86×86、diag=0、symmetric。
2. 重算 average-linkage → checkpoint：85 merges、root height=0.5044498568、group transition=1。
3. 打包 HTML → checkpoint：artifact validation/package verification PASS，0 external requests。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---|---:|---:|
| single-cell/colony/spatial orthogonal truth | 不能把 MG 稱為 cellular clone | 高 | P0 claim ceiling |
| matched-normal、CN/purity/multiplicity/CCF | 不能判斷甲基群的腫瘤特異性與細胞比例 | 中高 | P1 後續研究 |
| 跨樣本 prevalence | 不能由此位點外推 7 samples | 高 | P1，非本 Task F |

## §6 Evidence conflict / red-team scan

- Failure mode 1：直接使用 region 現成 124-read `tree.nwk`，會把 REF reads 混入 ALT-only M1。對策：只允許 assignment 的 86 IDs 重算。
- Failure mode 2：將 UPGMA 稱為 mutation phylogeny 或父子演化樹。對策：全報告固定稱「甲基距離 UPGMA」，並另列 mutation topology。
- Failure mode 3：熱圖與樹使用不同 leaf order，製造假 block。對策：單一 `leaf_order` 同時索引樹、rows 與 columns。
- Failure mode 4：81/5 視為兩個 cellular clones。對策：標示為 methyl groups；strict confirm=`NOT_RUN`。

未發現與「把既有距離資料正確視覺化」本身衝突的 concluded/NO-GO；既有負面結論主要限制生物學升級，不阻止 Task F 圖解。

## §7 Verdict

**GO；decision lock=N。** 可建置單一位點 HTML，但必須保留 DEMO/PARTIAL ribbon、完整 86-read 統計、124-vs-86 guard、同 leaf order、以及「不是 clone truth / 不是 mutation ancestry」限制。
