<!--
建立時間: 2026-07-17
目標: 保存第一次 all-sSNV source-authority 簽署失敗的完整稽核痕跡
處理範圍: Task Type B source authority v1 manifest 與 approval
關聯檔案: InterSubMod/docs/provenance/source_authorities/ARCHIVED.md
-->

# All-sSNV source authority v1 失敗紀錄

## 狀態

`NOT_SIGNED / NOT_AUTHORIZED / NOT_VALID_FOR_EXECUTION`

原始 JSON 保持內容不變，但檔名已加上不可忽略的狀態前綴：

- `UNSIGNED_NOT_AUTHORIZED.source_authority.v1.json`
- `UNSIGNED_NOT_AUTHORIZED.source_authority.v1.approval.json`

其中 approval JSON 是簽署前的候選 payload；即使內部曾帶有
`approval_status` 欄位，也因從未產生 detached signature 而不具任何授權效力。

## 原因

等待式 signer 繼承 Conda 舊版 `openssl`，執行 Ed25519
`pkeyutl -sign -rawin` 時回報 `Option unknown option -rawin`。因此：

- detached signature 從未產生；
- v1 私鑰維持 mode `0400`，未被誤標為已退休；
- validator 持續 fail closed；
- v1 manifest/approval 不得作為任何正式執行或科學結論的依據。

## 接替

正式鏈改用 v2 authority 與固定 `/usr/bin/openssl` 的一次性 signer。v2 必須重新通過完整測試、來源 digest 審查、兩位獨立 reviewer 核准與 detached-signature 驗證。
