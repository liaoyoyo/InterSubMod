<!--
建立時間: 2026-07-13 00:00 +0800
目標: 記錄跨樣本 composition adjustment 實作中的決策、偏離、折衷與未決項
處理範圍: 僅本 research topic 新目錄
關聯檔案:
  - InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/pre-decision-audit.md
-->

# Implementation notes

## 設計決定

1. Primary 使用 complete 五類，避免 resolved-only closure 忽略 unresolved。
2. Secondary 使用 resolved-only 四類，與 20260712 既有報告對齊。
3. HCC1395 technical pair 與 6 biological IDs 分開報告；cross-biological reference pairs 必須 `biological_id_a != biological_id_b`。
4. source standardization 使用全 7 resolved regions 的 pooled source mix：structural `21,976/37,039`、VAF-resolved `15,063/37,039`；另做 HCC-pair pooled／equal 0.5 sensitivity。
5. 固定 seed；所有 stochastic outputs以 stable sort＋固定浮點格式輸出。
6. Region annotation 將 `C_read_groups`、`T_exact_candidate_forest_count`、`Topo_shape_count` 與 `external_cluster_count` 分欄，禁止把 PyClone fit-local cluster 稱為 `C` 或真實 clone 數。
7. 每區保留 `tree_selection_source`：`structural_topo1` 是結構唯一；`vaf_resolved_topogt1` 是以每位點 VAF 排序後的最可能推測；unresolved／incomplete 不強迫選樹。
8. Region possible state 的 primary endpoint 要求所有 exact mutation 完整 join；`CP>=0.90` 只稱 clonal-like point estimate，`assignment>=0.80` 僅作 sensitivity。

## 偏離

- fresh layered-v3 尚未 root `_SUCCESS`，因此 region annotation 先以 historical layered-v2 實作並加 claim ceiling；不可升為 scientific release。

## 折衷

1. region 表只有 resolved rows；complete五類 block bootstrap需從 upstream complete-region表補回 unresolved座標，且必做 key守恆。
2. source是方法選擇狀態、可能為 outcome-dependent selector；direct standardization是 sensitivity，不是「校正後真值」。
3. 5 Mb common blocks作主 spatial bootstrap，10 Mb作敏感度；不宣稱5 Mb為生物最佳 block length。

## 未決

- fresh layered-v3 7/7完成後必須重跑，才能升正式 scientific Results。

## 執行後結果

1. 81/81 hard checks PASS；兩次完整5,000-replicate run的18/18輸出SHA-256一致。
2. Complete五類primary：HCC JSD=0.1984、relative rank=9；posterior／rarefaction／5 Mb／10 Mb／EB median rank皆9。
3. Resolved-only raw TV=8.60%，但structural conditional TV=28.24%、VAF conditional TV=8.17%；global pooled source weights direct-standardized TV=16.75%，確認mixture cancellation sensitivity。
4. Verdict：`MODERATE_NOT_EXCEPTIONAL_TECHNICAL_PROXIMITY`；不能稱高度一致或clone-tree驗證。
5. Region annotation：47,377/47,377 rows、32/32 checks PASS；兩次重跑四個 canonical outputs 全部 byte-identical；HCC pair 表另完整 join 5,720/5,720 GENCODE／CGC／DGIdb／CLP region flags。
6. HCC fixed 5,720 regions：兩側可評估 4,438，same possible state 4,296，但其中 4,267 是兩側皆 single modeled clonal-like；either possible-subclone 172、both 40，Jaccard=23.26%。
7. 真正 informative 的 both-multicluster partition 只有 34 區，label-invariant exact 21/34=61.76%；不可用 5,028/5,189 的表面 partition exact 當 headline，因 5,007 區都是 both-single-cluster。
8. `global_vs_regional_cluster_concordance.tsv` 發現既有 k-bin loop stale-variable bug；本輪報告排除該表，不以其 high-confidence k strata 下結論。
9. 紅隊修正 pair 明細：subclonal union=0 時 Jaccard 改為 `NA + both_absent`；partition exact 只在 both-multicluster 34 區定義，其他 rows 以 informative/vacuity fields fail closed。
10. 紅隊拆開三個 estimands：aggregate final-shape、matched pre-VAF coarse category、external PyClone state/partition；補列 matched final-shape agreement 4,243/5,720=74.18% 與 distance 0.1967，但不為其補造 permutation null。
11. Release gate 改為雙層：internal engineering `PASS_WITH_CAVEATS`；fresh layered-v3 root `_SUCCESS` 與 run-level verification 前，scientific/external `NO-GO`。
