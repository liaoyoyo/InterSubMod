<!--
建立時間: 2026-07-23
目標: 封存 recovery v24 在 authority 建立前被獨立審查拒絕的完整證據
處理範圍: 七個凍結來源、三方審查 transport、未使用金鑰與 15 個正式輸出槽位
關聯檔案: rejection_evidence.json
-->

# Recovery v24 pre-authority rejection

v24 依 strictest-review-wins 規則於 authority 建立前拒絕。Mendel 與 Nash
一致指出 regression 未實際執行 Node 與 Chromium，亦缺少 missing/substituted
runtime 的 fail-closed 測試。Nash 另指出 occupied-state regression 只覆蓋
334/349 slots，scope 漏列 v23，且 authority 宣稱 22 patterns 而實際為 23。

外部 Claude Opus 回報 APPROVE，但不凌駕兩份 REQUEST_CHANGES。v24
authority、正式 review 與 V/R/C 輸出均未建立；authority-v24 與 terminal-v14
私鑰保持未使用並標記永不重用。科學資料、469,849 位點、read-tag join 與
結論上限均未改變。
