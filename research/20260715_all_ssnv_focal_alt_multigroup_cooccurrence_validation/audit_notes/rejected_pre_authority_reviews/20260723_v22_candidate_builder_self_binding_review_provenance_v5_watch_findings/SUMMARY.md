<!--
建立時間: 2026-07-23
目標: 封存 recovery v22 在 authority 建立前被獨立審查拒絕的完整證據
處理範圍: 七個凍結來源、三方審查 transport、未使用金鑰與 15 個正式輸出槽位
關聯檔案: rejection_evidence.json
-->

# Recovery v22 pre-authority rejection

v22 依 strictest-review-wins 規則於 authority 建立前拒絕。Mendel 與 Nash 都指出
ceremony builder 未證明實際執行程式碼等於凍結來源；兩者也指出內部 reviewer UUID
只能作為主控 session 記錄的 transport attribution，不能作為密碼學 authorship 證明。
Nash 另發現 v21 原始 v5 中間輸出只做一次性 absence check，未納入持續 watcher；
Mendel 另指出 v21 封存欄位需要 exact-key 與直接 cross-link。

外部 Claude Opus 回報 APPROVE，但不凌駕兩份 REQUEST_CHANGES。v22 authority、正式 review、
V/R/C 輸出均未建立；兩組 v22 私鑰保持未使用並標記永不重用。科學資料、469,849 位點、
read-tag join 與結論上限均未改變。
