<!--
建立時間: 2026-08-01
更新時間: 2026-08-13
狀態: active index; publication/release gates remain fail-closed
目標: 提供 InterSubMod 跨研究者／跨 AI 交接文件入口
處理範圍: docs/handoff 下的 active handoff bundles
關聯檔案:
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/README.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md
驗證方式: 入口相對連結存在
-->

# InterSubMod Handoff Index

## Active

- [**2026-08-13 完整研究資料與軟體交接 snapshot**](20260813_完整研究資料與軟體交接_01/00_INDEX.md)
  狀態：`RESEARCH_HANDOFF_SNAPSHOT`；公開 claim、資料路徑、軟體 I/O、bip7/bip8 preflight 與 release gates 的第一入口。
  邊界：目前不得稱 release-ready；科學數值 authority **仍是 2026-08-01 frozen bundle**，08-13 只負責治理、索引與驗證鏈。
- [Exact-PS×HP、read-AF、CNV 與區域樹研究交接](20260801_exactPS_readAF_CNV_AI交接_01/README.md)  
  狀態：`FROZEN_TECHNICAL_RESULTS_WITH_EXPLICIT_ABSTAIN`；仍是科學數值 authority。
  適用：新研究者、新 AI、論文更新、read-AF/CNV 方法審查與單一樹輸出決策。
- [Exact-PS 全 7 資料集 HTML observation report](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report.standalone.html)  
  狀態：`VALIDATED_DERIVED_OBSERVATION`；responsive、no-JS、A4 與 receipt QA 均 PASS。  
  重生方式：[Python builder／finalizer 說明](../../research/20260801_exact_ps_observation_report/00_INDEX.md)。

## Historical / superseded

> 以下只是有日期的歷史快照，不得覆蓋本頁列出的 active handoff 導覽。科學數值的唯一權威仍是 2026-08-01 frozen authority bundle。Markdown / HTML 已加上 `HISTORICAL / SUPERSEDED` 轉址；原始 JSON receipt 為保留歷史證據而不改寫。

- [2026-08-05 系統交接與驗收 bundle](20260805_系統交接與驗收_01/README.md)：歷史資料層／程式層／輸出層與當日驗收紀錄。
  - [當日交接 HTML](20260805_系統交接與驗收_01/00_交接總覽與驗收.standalone.html)
  - [當日 repo 整理 HTML](20260805_系統交接與驗收_01/05_repo整理與可觀察性盤點.standalone.html)
- [2026-08-06 LongLineage 充分性稽核與路線裁決](20260806_LongLineage充分性稽核與路線裁決_01.md)：歷史 candidate 稽核，不代表現行 public-preview gate PASS。
- [2026-08-06 兩 repo 端到端串接盤點](20260806_兩repo端到端串接可行性盤點_01.md)：歷史串接判斷；[當日 HTML](20260806_兩repo端到端串接可行性盤點_01.standalone.html)同樣已降級。

## 使用原則

1. 每個 handoff bundle 的 machine-readable authority manifest 是數字與路徑入口。
2. 舊文件可供理解歷史，但不得覆蓋 active bundle 的 claim boundary 與分母。
3. 若 artifact hash、scope 或 denominator 不一致，停止引用並重新稽核。
