<!--
建立時間: 2026-08-13
目標: 在建立多 BAM／多樣本一頁式 HTML dashboard 前，確認 scope、資料可比性、claim ceiling 與交付路徑
處理範圍: 7 份 topology datasets、HCC1395 v1/v3 extensions、portable canonical dashboard artifact；不重跑 raw BAM/C++/全樣本 ISM
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_HCC1395_drilldown完整驗證與多樣本改進_01.md
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/implementation-notes.md
-->

# Pre-decision audit：多 BAM／多樣本一頁式分析總覽

## Verdict

**GO for a source-backed portable snapshot dashboard；NO-GO for production/live or cohort methylation claims。**

本輪可以用已驗證的 7 份 topology 摘要與 HCC1395-only extensions 建立「總覽→比較→單樣本 detail→方法/來源」的一頁式 dashboard。它是整合介面與資料 contract 的候選版，不代表已重跑多 BAM、不代表 7 樣本 ISM/lineage 完成，也不是 live refresh 產品。

## 任務分類與服務目標

- **Task Type**：B — Comprehensive validation；7 datasets 全納入，不挑 subset。
- **服務 G3**：把 read-level epigenetic 可用性、缺失與 axis validity 明示。
- **服務 G4**：用 biological sample / technical replicate 分層，支援逐樣本比較而不 pooled。
- **服務 G5**：每個 card/chart/table 有來源、公式、分母與 portable verification。

## 啟動研究 5 問

| 問題 | 判定 | 實作約束 |
|---|---|---|
| Thread D read-level epigenetic？ | **是** | HCC1395 甲基只作 extension panel；其他樣本顯示 ABSENT/PARTIAL |
| Thread B 撤回範圍？ | **否** | 不把 legacy A/B 或舊 raw linkage 搬入 cohort 指標 |
| KDE-corrected？ | **未知** | dashboard 明示 provenance unknown；不可默認已校正 |
| 需要 VCF caller AF？ | **否（本頁）** | topology 指標使用保存的 exact-PS family-specific AF basis；不聲稱 caller performance |
| 長計算/C++/搬移/NO-GO gate？ | **否** | 僅聚合已審核 CSV/JSON；不啟動 raw BAM 或全資產重建 |

## 決策選項

| 選項 | 收益 | 風險 | 決定 |
|---|---|---|---|
| A. Canonical artifact → self-contained HTML | 可追溯、可驗證、sample filter、semantic fallback | snapshot 非 live；shared renderer 可能有 overflow | **採用** |
| B. 直接改 production generator 加 cohort route | 最接近長期產品 | manifest/refresh/樣本 contract 未凍結，變動面太大 | 後續 cycle |
| C. 手刻 bespoke standalone | 視覺自由 | 產生第二套 runtime，來源與數字易 drift | 不採用 |

## 關鍵假設

1. **整體母體**是 7 datasets = 6 biological samples + 1 HCC1395_DORADO technical replicate；不得把 dataset n=7 寫成 biological n=7。
2. **跨樣本可比層**只有 topology/MLHP contract；HCC1395 v1/v3 的 ISM、lineage、panel、IGV 是 sample-specific extension。
3. **All samples** 的 rate 預設用 6 biological samples 的 macro median，不用 pooled loci；volume 可用 region records 加總，但標明是 records、不是 unique biological events。
4. **樣本 selector** 必須同時作用於 cards、charts、funnel 與 details；選到沒有 extension 的樣本要顯示空狀態，不可沿用 HCC1395。
5. **HTML 是 snapshot**，不是 BAM watcher、排程器或 production refresh；產生時間與來源檔必須可見。

## Pre-mortem 與防線

| 可能失敗 | 最早可見訊號 | 防線 / 驗證 |
|---|---|---|
| 技術 replicate 被當獨立 biological n | KPI 顯示 biological=7 | reconciliation assert：datasets=7、biological=6、technical=1 |
| selected sample 仍看到 HCC methyl | 非 HCC filter 下 methyl chart 有 rows | 每個 dataset 使用共同 `sample_filter`；negative browser/filter check |
| pooled 大樣本主導 rate | H2009 改變 cohort headline 過多 | All-sample rate 使用 biological macro median；逐樣本 numerator/denominator保留 |
| missing 被畫成 0 | 無 extension 樣本出現 0% methyl | 狀態欄使用 ABSENT/PARTIAL；無資料不進數值 chart |
| raw p 被讀成驗證 | axis chart headline 寫 significance success | axis 放 diagnostic/detail；cluster invalid、raw p/FDR/causal caveat就近顯示 |
| shared renderer overflow | verifier desktop/narrow scrollWidth > clientWidth | canonical deliver 必須 `verification=passed`；失敗則不交付 HTML |
| custom page 與 artifact drift | HTML 數字和 artifact 不一致 | HTML 只能由 canonical `multi_bam_dashboard_artifact.json` 產生 |

## Step → Verify

1. 建 metric/data contract → 驗證：每個 metric 有 grain、公式、分母、來源與跨樣本 eligibility。
2. 生成 canonical dashboard artifact → 驗證：JSON parse；sourceId 全解析；percent 全為 0–1；status=`partial`。
3. 打包 self-contained HTML → 驗證：canonical deliver receipt `stages.verification=passed`、無外部 request／browser error／horizontal overflow。
4. Reconcile 預設與 sample selector → 驗證：All=7/6/1；H2009=36,042 regions/64.17% tree/35.29% unique-among-tree；HCC extension 只在 HCC1395/All 顯示。
5. 交付與限制 → 驗證：頁面首屏見 PARTIAL、claim ceiling、freshness；detail/來源在後段可展開或查詢。

## Exit criteria

- **可交付**：portable HTML verifier PASS；cards/charts/tables 對回 reviewed CSV/JSON；selected sample 不污染；明示 snapshot/PARTIAL。
- **不可交付**：任何來源無法對回、filter contamination、renderer overflow、或將 HCC extension 外推 cohort。
