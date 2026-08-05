<!--
建立時間: 2026-07-16 08:01 +08:00
目標: 記錄 schema 1.1 HCC1954 chr22 診斷重跑被科學契約更新取代後的終止原因
處理範圍: 此目錄內未完成的 extraction artifacts；不可作數據來源
關聯檔案: InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/implementation-notes.md
-->

# ABORTED — 不可使用此目錄內的資料

- 原執行 scope：HCC1954 chr22，extractor schema 1.1，HP-family components。
- 終止時間：2026-07-16 08:01 +08:00。
- 終止原因：獨立 red-team 證實 schema 1.1 會把不同 phase-set 的同號 HP 合併；正式 primary 必須改為 `HP family × known PS × read-linked component`。
- 執行結果：收到 SIGINT，exit code 130；wall time 38m21.64s；peak RSS 93,436 KB。
- 完整性：`molecule_sparse_calls.tsv.gz` 為中途檔，缺 BGZF/GZIP EOF marker；沒有 `receipt.json`；不得讀取、修復、續用或納入任何比例。
- 後續：保留本目錄作可稽核失敗紀錄；schema 1.2 PS-aware pilot 必須使用全新輸出目錄。
