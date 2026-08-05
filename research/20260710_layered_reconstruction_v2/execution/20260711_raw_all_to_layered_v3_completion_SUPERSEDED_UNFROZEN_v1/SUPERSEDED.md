<!--
建立時間: 2026-07-11 08:41 Asia/Taipei
目標: 記錄未進入科學執行的 supervisor 為何停止
處理範圍: execution supervisor only
關聯檔案: execution.log, _FAILED, entrypoint_code.sha256
-->

# Superseded Before Scientific Execution

此 supervisor 只執行到等待 producer terminal state，未建立 production、canonical 或 sensitivity run root。

停止原因：啟動時只記錄 live entrypoint hashes，等待後未重新檢查，存在 TOCTOU source drift 風險。`SIGINT` exit 130 是主動 fail-closed；後繼版本必須在 producer 完成後、每個 manifest/run/summary 步驟前重驗 source 與 input authority hashes。
