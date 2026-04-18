<!--
建立時間: 2026-04-11 22:30
目標: docs/data_specs/ 目錄索引
處理範圍: 資料規格文件的導覽
-->

# 資料規格文件索引

> InterSubMod 資料流程、欄位定義、輸出結構的集中式規格文件。

## 文件清單

| 文件 | 說明 | 關鍵內容 |
|------|------|---------|
| [資料盤點快照](20260411_資料盤點快照_01.md) | 儲存拓撲、輸出目錄樹、路徑驗證結果 | 儲存配置速查、44 broken paths 分類 |
| [完整資料流程圖](20260411_完整資料流程圖_01.md) | 端對端資料流 (Mermaid) + 輸入/輸出規格 | 5 層流程圖、C++ flags、腳本對應表 |
| [significance_summary 欄位字典](20260411_significance_summary欄位字典_01.md) | C++ 輸出 59 欄完整定義 | 12 群組、值域、計算公式、已知限制 |
| [master dataset 欄位字典](20260411_master_dataset欄位字典_01.md) | 116 欄完整定義 (59 原始 + 57 聚合) | 模式差異、大寫/小寫對照、caller 欄位 |
| [工作區命名與目錄結構](20260411_工作區命名與目錄結構_01.md) | Canonical run + Workspace 命名規則與完整清單 | 56 workspace 狀態標記、後綴規則 |
| [Output 資料信任度與生命週期](20260414_output資料信任度與生命週期_01.md) | 信任度 4 級、PRE-FIX 範圍、haplotag 重跑 checklist | CURRENT/PRE-FIX/SUPERSEDED/DEPRECATED 定義 |

## 輔助資料

| 檔案 | 說明 |
|------|------|
| [path_inventory.tsv](20260411_path_inventory.tsv) | 路徑驗證原始結果（589 行 TSV，可排序篩選） |

## 快速查閱指引

| 我想知道... | 去哪裡找 |
|-------------|---------|
| 某個 CSV 欄位什麼意思 | significance_summary 或 master dataset 欄位字典 |
| 資料從哪裡來、經過什麼處理 | 完整資料流程圖 |
| 某個 workspace 的狀態是否有效 | 工作區命名與目錄結構 §3 |
| 某個路徑是否存在 | path_inventory.tsv 或執行 `verify_data_paths.py` |
| 磁碟/目錄配置 | 資料盤點快照 §1-2 |
