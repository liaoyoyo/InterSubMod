<!--
建立時間: 2026-07-17
目標: 指向已封存的失敗 source-authority 嘗試
處理範圍: all-sSNV Task Type B source authority v1
關聯檔案: InterSubMod/docs/archive/2026/07/20260717_failed_all_ssnv_source_authority_v1/SUMMARY.md
-->

# 封存提示

未簽署的 all-sSNV source-authority v1 已封存至：

`InterSubMod/docs/archive/2026/07/20260717_failed_all_ssnv_source_authority_v1/`

該版本狀態為 `NOT_SIGNED / NOT_AUTHORIZED`，不得作為正式執行依據。正式來源授權使用同目錄下的 v2 manifest、approval 與 detached signature。

封存 JSON 的實際檔名均以 `UNSIGNED_NOT_AUTHORIZED.` 開頭；即使不開啟
`SUMMARY.md`，也不能把它們誤認為有效 authority。

source-authority v2 雖有有效 detached signature，但 live validator 因
`/dev/stdin` pipe 不相容而拒絕，亦已封存至：

`InterSubMod/docs/archive/2026/07/20260717_signed_but_validator_incompatible_all_ssnv_source_authority_v2/`

其狀態為 `SIGNATURE_VALID / NOT_AUTHORIZED_FOR_EXECUTION`；正式執行只接受後續 v3 authority。

2026-07-18 第一次 v5 ceremony 已組裝 authority/approval，但 source signer
在產生 detached signature 前失敗，且 one-time passphrase 隨 process
退出而遺失。未簽署 JSON 與舊 key pair 已封存至：

`InterSubMod/docs/archive/2026/07/20260718_unsigned_source_authority_v5_signer_failure_01/`

其狀態為 `UNSIGNED / NOT_AUTHORIZED / NEVER_USED_FOR_FORMAL_EXECUTION`。
後續正式 v5 必須綁定新的 key、source-set digest、canonical tests 與兩位
新的 reviewer approval；不得重用本次 approval。
