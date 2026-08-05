<!--
建立時間: 2026-07-22
目標: 保存 portable HTML 產生過程的 fail-closed browser QA 截圖
處理範圍: report renderer failures only；不代表 HCC1395 scientific-data failure
關聯檔案: InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/report_delivery_receipt.json
-->

# Portable report QA failure evidence

這些 screenshot 是成功發布前保留的 renderer failure evidence：

1. 初版寬表與 shared reader top-bar 都觸發 horizontal-overflow gate。
2. 將 validation／definition／outlier tables 拆窄後，剩餘 overflow 固定為 shared top bar 在 Linux non-overlay scrollbar 下的 8 px `100vw` 差異。
3. 最終報告重用 20260718 scrollbar-safe runtime wrapper，並由 canonical delivery lifecycle 通過 1440/390 viewport 與 source-dialog QA。

最終成功狀態只以 `../report_delivery_receipt.json` 為 authority；本目錄不可用來否定或替代最終 receipt。
