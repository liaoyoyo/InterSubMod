<!--
建立時間: 2026-07-23
目標: 封存 recovery v23 在 authority 建立前被獨立審查拒絕的完整證據
處理範圍: 七個凍結來源、三方審查 transport、未使用金鑰與 15 個正式輸出槽位
關聯檔案: rejection_evidence.json
-->

# Recovery v23 pre-authority rejection

v23 依 strictest-review-wins 規則於 authority 建立前拒絕。Mendel 與 Nash
獨立指出 continuation 將必需的 Node runtime 指向不存在的 `v23.22.1`；
固定 SHA-256 實際對應仍存在的 `v22.22.1`。Mendel 另指出 558 項測試與
唯讀 probe 沒有實際驗證所有必需執行檔的存在性與 SHA-256。

外部 Claude Opus 回報 APPROVE，但不凌駕兩份 REQUEST_CHANGES。v23
authority、正式 review 與 V/R/C 輸出均未建立；兩組 v23 私鑰保持未使用並
標記永不重用。科學資料、469,849 位點、read-tag join 與結論上限均未改變。
